//! Small-indel site model and length-changing read rewrite for the ONT denoiser.
//! POS, reference end, and coordinate-sort order are therefore invariant, so the output never needs
//! re-sorting or re-indexing.

use std::collections::{HashMap, HashSet};

use rust_htslib::bam::record::{Cigar, CigarString};

// ===========================================================================
/// An indel allele at a normalized site, indexed conceptually by its net length change: `Ref` = 0.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Allele {
    Ref,
    Ins(Vec<u8>),
    Del(u32),
}

impl Allele {
    /// Net length change this allele applies to the read relative to the reference.
    pub fn net_len(&self) -> i32 {
        match self {
            Allele::Ref => 0,
            Allele::Ins(s) => s.len() as i32,
            Allele::Del(n) => -(*n as i32),
        }
    }

    /// True when both alleles are the same variety (both insertions, both deletions, or both reference).
    pub fn kind_matches(&self, other: &Allele) -> bool {
        matches!(
            (self, other),
            (Allele::Ref, Allele::Ref)
                | (Allele::Ins(_), Allele::Ins(_))
                | (Allele::Del(_), Allele::Del(_))
        )
    }
}

/// Minimum plain match-context flank required on each side of an indel edit.
pub const INDEL_FLANK: usize = 5;

/// Tunable parameters for indel denoising.
#[derive(Debug, Clone)]
pub struct IndelOpts {
    pub enabled: bool,
    /// OFF by default (CLI `--indel-hp-length-model`) until an eval shows it helps.
    pub hp_length_model: bool,
    /// Exclusive upper bound: with 5, lengths 1..=4 are correctable.
    pub max_len: u32,
    pub vaf: f64,
    pub err0: f64,
    pub err_scale: f64,
    pub err_cap: f64,
    pub floor_mult: f64,
    pub delta: f64,
    /// An allele carried by at least this fraction of reads at a site is never corrected away.
    pub protect_vaf: f64,
    /// Minimum VAF for an indel alt to be kept.
    pub vaf_abs_min: f64,
    pub excess_alpha: f64,
}

impl Default for IndelOpts {
    fn default() -> Self {
        IndelOpts {
            enabled: false,
            hp_length_model: false,
            max_len: 5,
            vaf: 0.05,
            err0: 0.01,
            err_scale: 1.5,
            err_cap: 0.4,
            floor_mult: 3.0,
            delta: 0.3,
            protect_vaf: 0.2,
            vaf_abs_min: 0.05,
            excess_alpha: 1e-3,
        }
    }
}

/// Shift an indel maximally left so the same biological event, however the aligner anchored it.
pub fn normalize_left(refseq: &[u8], anchor: u32, allele: &Allele) -> (u32, Allele) {
    match allele {
        Allele::Ref => (anchor, Allele::Ref),

        // Deleted bases are refseq[a+1 ..= a+n].
        Allele::Del(n) => {
            let n = *n as usize;
            let mut a = anchor as usize;
            while a > 0 {
                match (refseq.get(a), refseq.get(a + n)) {
                    (Some(entering), Some(leaving)) if entering == leaving => a -= 1,
                    _ => break,
                }
            }
            (a as u32, Allele::Del(n as u32))
        }

        // `s` is inserted after refseq[a].
        Allele::Ins(s) => {
            let mut a = anchor as usize;
            let mut s = s.clone();
            while a > 0 && !s.is_empty() {
                match refseq.get(a) {
                    Some(&b) if b == *s.last().unwrap() => {
                        s.rotate_right(1);
                        a -= 1;
                    }
                    _ => break,
                }
            }
            (a as u32, Allele::Ins(s))
        }
    }
}

/// Smallest prefix `p` such that `unit` is `p` repeated (e.g. "AA"->"A", "ATAT"->"AT").
fn primitive_period(unit: &[u8]) -> &[u8] {
    let n = unit.len();
    for p in 1..=n {
        if n % p == 0 && unit.chunks(p).all(|c| c == &unit[..p]) {
            return &unit[..p];
        }
    }
    unit
}

/// Reference repeat context `L` at a normalized site.
pub fn repeat_context(refseq: &[u8], norm_pos: u32, allele: &Allele) -> u32 {
    let unit: Vec<u8> = match allele {
        Allele::Ref => return 1,
        Allele::Ins(s) => s.clone(),
        Allele::Del(n) => {
            let a = norm_pos as usize + 1;
            match refseq.get(a..a + *n as usize) {
                Some(w) => w.to_vec(),
                None => return 1,
            }
        }
    };
    if unit.is_empty() {
        return 1;
    }
    let unit = primitive_period(&unit);
    let u = unit.len();
    let start = norm_pos as usize + 1;

    // Copies running rightward from the event position...
    let mut copies = 0u32;
    let mut p = start;
    while refseq.get(p..p + u).map_or(false, |w| w == unit) {
        copies += 1;
        p += u;
    }
    // ...and leftward.
    let mut p = start;
    while p >= u && refseq.get(p - u..p).map_or(false, |w| w == unit) {
        copies += 1;
        p -= u;
    }
    copies.max(1)
}

/// Returns `anchor + 1` when there is no run, i.e. when that base does not repeat.
pub fn hp_tract_end(refseq: &[u8], anchor: u32) -> u32 {
    let start = anchor as usize + 1;
    let base = match refseq.get(start) {
        Some(&b) => b,
        None => return start as u32,
    };
    let mut p = start + 1;
    while refseq.get(p) == Some(&base) {
        p += 1;
    }
    // A lone base is not a run: fall back to `anchor + 1`.
    if p - start >= 2 {
        p as u32
    } else {
        start as u32
    }
}

/// If so, returns the run's exclusive reference end `e` and its base `b`.
pub fn hp_site(refseq: &[u8], anchor: u32) -> Option<(u32, u8)> {
    let start = anchor as usize + 1;
    let e = hp_tract_end(refseq, anchor);
    if (e as usize) < start + 2 {
        return None;
    }
    let b = *refseq.get(start)?;
    if anchor > 0 && refseq.get(anchor as usize) == Some(&b) {
        return None;
    }
    Some((e, b))
}

pub fn hp_length_allele(k: i32, b: u8) -> Allele {
    match k {
        0 => Allele::Ref,
        k if k > 0 => Allele::Ins(vec![b; k as usize]),
        k => Allele::Del((-k) as u32),
    }
}

/// Per-junction indel error probability, growing with repeat context and capped.
pub fn error_rate(l: u32, o: &IndelOpts) -> f64 {
    let exp = (l as i32 - 1).max(0);
    (o.err0 * o.err_scale.powi(exp)).min(o.err_cap)
}

/// VAF floor for an allele in context `L`: an absolute minimum.
#[cfg_attr(not(test), allow(dead_code))]
pub fn vaf_floor(l: u32, o: &IndelOpts) -> f64 {
    o.vaf.max(o.floor_mult * error_rate(l, o))
}

fn indel_excess_pvalue(n: u32, depth: u32, eps: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    let p = eps.clamp(1e-9, 1.0 - 1e-9);
    match statrs::distribution::Binomial::new(p, depth as u64) {
        Ok(b) => statrs::distribution::DiscreteCDF::sf(&b, (n - 1) as u64),
        Err(_) => 1.0,
    }
}

/// Accumulated observations at one left-normalized indel site.
#[derive(Debug, Default, Clone)]
pub struct IndelSite {
    /// Non-reference events only: allele -> (forward count, reverse count).
    pub obs: HashMap<Allele, (u32, u32)>,
    /// Total reads spanning this site (the VAF denominator).
    pub depth: u32,
    /// Column strand totals, used as the null for the strand-bias test.
    pub col_fwd: u32,
    pub col_rev: u32,
    /// Hashes of the qnames of reads that have already cast an ALT vote at this site.
    voted: std::collections::HashSet<u64>,
    /// At an HP site this, and NOT `depth - alt_total`, is the reference support.
    hp_ref: u32,
}

impl IndelSite {
    fn add(&mut self, allele: Allele, reverse: bool) {
        let e = self.obs.entry(allele).or_insert((0, 0));
        if reverse {
            e.1 += 1;
        } else {
            e.0 += 1;
        }
    }

    /// Like `add`, but caps a single read to at most one ALT vote at this site.
    pub fn add_from_read(&mut self, qname: &[u8], allele: Allele, reverse: bool) -> bool {
        if !self.mark_voted(qname) {
            return false;
        }
        self.add(allele, reverse);
        true
    }

    /// Shares the one-vote-per-read cap with `add_from_read`, so a read counts once here or once there.
    pub fn add_hp_ref(&mut self, qname: &[u8]) -> bool {
        if !self.mark_voted(qname) {
            return false;
        }
        self.hp_ref += 1;
        true
    }

    /// Claim this read's single vote at this site.
    fn mark_voted(&mut self, qname: &[u8]) -> bool {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        qname.hash(&mut h);
        self.voted.insert(h.finish())
    }

    fn alt_total(&self) -> u32 {
        self.obs.values().map(|&(f, r)| f + r).sum()
    }
}

/// The decision for one site: which alleles are kept, their renormalized frequencies.
#[derive(Debug, Clone)]
pub struct IndelSiteModel {
    pub kept: Vec<(Allele, f64)>,
    pub eps: f64,
    /// The repeat context length that produced `eps`.
    pub context_len: u32,
    pub strand_biased: usize,
    /// The RAW observed VAF of every allele at the site (kept or not), including `Allele::Ref`.
    /// Unlike `kept`, this is never renormalized, so it is the correct source for the `protect_vaf`
    /// substantiality check in `assign_allele` -- candidacy rejection must not affect protection.
    pub obs_vaf: Vec<(Allele, f64)>,
    /// Alleles rejected specifically by a strand gate.
    /// Sorted for the same determinism reason as `kept` and `obs_vaf`.
    pub strand_rejected: Vec<Allele>,
    /// `None` otherwise -- both when the model is off and when the site is not a homopolymer.
    pub hp_run: Option<(u32, u8)>,
}

/// Decide which alleles survive at a site and compute their frequencies.
pub fn fit_indel_site(
    site: &IndelSite,
    refseq: &[u8],
    norm_pos: u32,
    min_strand: u32,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
    n_sites: usize,
    o: &IndelOpts,
) -> IndelSiteModel {
    let depth = site.depth.max(1);
    let col_total = site.col_fwd + site.col_rev;
    let expected_fwd_frac = if col_total == 0 {
        0.5
    } else {
        site.col_fwd as f64 / col_total as f64
    };

    // One eps per site.
    let context_len = site
        .obs
        .keys()
        .map(|a| repeat_context(refseq, norm_pos, a))
        .max()
        .unwrap_or(1);
    let eps = error_rate(context_len, o);

    debug_assert!(
        site.alt_total() <= site.depth,
        "indel site accounting invariant violated: alt_total ({}) exceeds depth ({})",
        site.alt_total(),
        site.depth
    );
    // The slack is the uninformative reads.
    debug_assert!(
        site.hp_ref + site.alt_total() <= site.depth,
        "HP site accounting invariant violated: hp_ref ({}) + alt_total ({}) exceeds depth ({})",
        site.hp_ref,
        site.alt_total(),
        site.depth
    );

    // Every downstream consumer keys off `hp_run.is_some()`.
    let hp_run = if o.hp_length_model { hp_site(refseq, norm_pos) } else { None };

    let mut strand_biased = 0usize;
    // Raw counts over kept alleles.
    let mut raw: Vec<(Allele, f64)> = Vec::new();
    // Everywhere else REF still gets whatever depth the events did not claim.
    let ref_support = if hp_run.is_some() {
        site.hp_ref
    } else {
        depth.saturating_sub(site.alt_total())
    };
    raw.push((Allele::Ref, ref_support as f64));

    // RAW observed VAF for every allele, whether kept or not.
    let mut obs_vaf: Vec<(Allele, f64)> = Vec::new();
    obs_vaf.push((Allele::Ref, (ref_support as f64 / depth as f64).min(1.0)));
    let mut strand_rejected: Vec<Allele> = Vec::new();

    for (allele, &(fwd, rev)) in site.obs.iter() {
        let n = fwd + rev;
        obs_vaf.push((allele.clone(), (n as f64 / depth as f64).min(1.0)));
        if n == 0 {
            continue;
        }
        if fwd < min_strand || rev < min_strand {
            strand_rejected.push(allele.clone());
            continue;
        }
        let l = repeat_context(refseq, norm_pos, allele);
        let vaf = (n as f64 / depth as f64).min(1.0);
        let eps_l = error_rate(l, o);
        if vaf < o.vaf_abs_min.max(1.5 * eps_l)
            || indel_excess_pvalue(n, depth, eps_l) > o.excess_alpha / n_sites.max(1) as f64
        {
            continue;
        }
        // Near-homoplasmic exemption: a single-strand artifact is absent from the other strand's reads.
        if strand_bias_p > 0.0
            && vaf < homoplasmic_vaf
            && crate::call::strand_bias_pvalue(fwd as usize, rev as usize, expected_fwd_frac)
                < strand_bias_p
            && (fwd as f64 / (fwd + rev) as f64 - expected_fwd_frac).abs()
                >= crate::STRAND_BIAS_MIN_EFFECT
        {
            strand_biased += 1;
            strand_rejected.push(allele.clone());
            continue;
        }
        raw.push((allele.clone(), n as f64));
    }

    // Renormalize over the kept set, so mass on rejected alleles does not leak.
    let sum: f64 = raw.iter().map(|(_, c)| *c).sum();
    let mut kept: Vec<(Allele, f64)> = if sum > 0.0 {
        raw.into_iter().map(|(a, c)| (a, c / sum)).collect()
    } else {
        vec![(Allele::Ref, 1.0)]
    };

    // Deterministic order: `site.obs` is a HashMap with a randomly-seeded hasher.
    kept.sort_by(|(a, _), (b, _)| a.cmp(b));
    obs_vaf.sort_by(|(a, _), (b, _)| a.cmp(b));
    strand_rejected.sort();

    IndelSiteModel { kept, eps, context_len, strand_biased, obs_vaf, strand_rejected, hp_run }
}

/// Outcome of MAP assignment for one read at one site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Assignment {
    /// The read already carries the MAP allele, or the site is degenerate.
    Keep,
    /// Reassign the read to this allele.
    Move(Allele),
    /// The MAP allele is only reachable through a transition the rewrite walk does not support.
    Refused,
}

/// Maximum-a-posteriori allele for a read observing `observed`.
pub fn assign_allele(observed: &Allele, m: &IndelSiteModel, o: &IndelOpts) -> Assignment {
    let eps = m.eps.clamp(0.0, 0.5);

    // Error mass spreads over the alleles that could have been observed given truth k, decaying in net-length distance.
    let mut best: Option<(&Allele, f64)> = None;
    for (k, f) in &m.kept {
        let like = if k == observed {
            1.0 - eps
        } else {
            // Z normalizes over the kept set PLUS the observed allele.
            let mut z: f64 = m
                .kept
                .iter()
                .filter(|(a, _)| a != k)
                .map(|(a, _)| o.delta.powi((a.net_len() - k.net_len()).abs()))
                .sum();
            if !m.kept.iter().any(|(a, _)| a == observed) {
                z += o.delta.powi((observed.net_len() - k.net_len()).abs());
            }
            if z <= 0.0 {
                continue;
            }
            eps * o.delta.powi((observed.net_len() - k.net_len()).abs()) / z
        };
        let val = f * like;
        if best.map_or(true, |(_, b)| val > b) {
            best = Some((k, val));
        }
    }

    let target = match best {
        Some((t, _)) => t,
        None => return Assignment::Keep,
    };
    if target == observed {
        return Assignment::Keep;
    }
    // Restriction 1: the rewrite walk cannot turn an insertion into a deletion or vice versa.
    if m.hp_run.is_none() {
        match (observed, target) {
            (Allele::Ref, _) | (_, Allele::Ref) => {}
            (a, b) if !a.kind_matches(b) => return Assignment::Refused,
            _ => {}
        }
        if let (Allele::Del(n), Allele::Del(mm)) = (observed, target) {
            if mm > n {
                return Assignment::Refused;
            }
        }
    }
    // Protection: an allele carried by a substantial fraction of reads is never corrected away, even when it lost candidacy.
    let strand_rejected = m.strand_rejected.iter().any(|a| a == observed);
    let observed_vaf = m
        .obs_vaf
        .iter()
        .find(|(a, _)| a == observed)
        .map_or(0.0, |(_, v)| *v);
    if !strand_rejected && observed_vaf >= o.protect_vaf {
        return Assignment::Keep;
    }
    Assignment::Move(target.clone())
}

/// One indel carried by a read, in normalized site coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadEvent {
    pub norm_pos: u32,
    /// The reference position of the last aligned base BEFORE the event.
    pub anchor: u32,
    /// `None` means the indel is at or above `max_len`: untouchable.
    pub allele: Option<Allele>,
}

/// Length of the contiguous `M`/`=`/`X` run starting AT `target`.
pub fn match_run_len(pos: i64, cigar: &CigarString, target: u32) -> u32 {
    let mut rp = pos.max(0) as u32;
    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                if target >= rp && target < rp + n {
                    return rp + n - target;
                }
                rp += n;
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => {
                if target >= rp && target < rp + n {
                    return 0;
                }
                rp += n;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    0
}

/// Enumerate a read's own indel events straight from its CIGAR, left-normalized.
pub fn read_events(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    refseq: &[u8],
    max_len: u32,
) -> Vec<ReadEvent> {
    let mut out = Vec::new();
    let mut rp = pos.max(0) as u32; // next reference position
    let mut qp = 0usize; // next query position
    // Reference position of the last ALIGNED (M/=/X) base seen so far.
    let mut last_aligned: Option<u32> = None;

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                rp += n;
                qp += n as usize;
                if n > 0 {
                    last_aligned = Some(rp - 1);
                }
            }
            Cigar::Ins(n) => {
                // No aligned base has been seen yet: there is no anchor to report, so skip the event entirely.
                if let Some(anchor) = last_aligned {
                    if n < max_len {
                        let s = seq[qp..qp + n as usize].to_vec();
                        let (np, a) = normalize_left(refseq, anchor, &Allele::Ins(s));
                        out.push(ReadEvent { norm_pos: np, anchor, allele: Some(a) });
                    } else {
                        out.push(ReadEvent { norm_pos: anchor, anchor, allele: None });
                    }
                }
                qp += n as usize;
            }
            Cigar::Del(n) => {
                if let Some(anchor) = last_aligned {
                    if n < max_len {
                        let (np, a) = normalize_left(refseq, anchor, &Allele::Del(n));
                        out.push(ReadEvent { norm_pos: np, anchor, allele: Some(a) });
                    } else {
                        out.push(ReadEvent { norm_pos: anchor, anchor, allele: None });
                    }
                }
                rp += n;
            }
            Cigar::SoftClip(n) => qp += n as usize,
            Cigar::RefSkip(n) => rp += n,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    out
}

/// Size of the overlap between the half-open reference intervals `[s1, e1)` and `[s2, e2)`.
fn overlap(s1: u32, e1: u32, s2: u32, e2: u32) -> u32 {
    e1.min(e2).saturating_sub(s1.max(s2))
}

/// The value is the net length change the read applies to that run, `k = -`.
pub fn hp_allele(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    refseq: &[u8],
    a: u32,
    e: u32,
) -> Option<i32> {
    if e < a + 2 {
        return None;
    }
    let b = *refseq.get(a as usize + 1)?;
    if (pos.max(0) as i64) > a as i64 {
        return None;
    }
    // The canonical rewrite emits the anchor base first, so the read must actually have one there.
    if match_run_len(pos, cigar, a) == 0 {
        return None;
    }

    let mut rp = pos.max(0) as u32; // next reference position
    let mut qp = 0usize; // next query position
    let mut last_aligned: Option<u32> = None;
    let mut k: i32 = 0;

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                rp += n;
                qp += n as usize;
                if n > 0 {
                    last_aligned = Some(rp - 1);
                }
            }
            Cigar::Ins(n) => {
                // An insertion belongs to this run when it is anchored on the anchor base itself or on any base of the run.
                let inside = last_aligned.is_some_and(|l| l >= a && l < e);
                if inside {
                    if seq.get(qp..qp + n as usize)?.iter().any(|&c| c != b) {
                        return None;
                    }
                    k += n as i32;
                }
                qp += n as usize;
            }
            Cigar::Del(n) => {
                k -= overlap(rp, rp + n, a + 1, e) as i32;
                rp += n;
            }
            Cigar::RefSkip(n) => {
                if overlap(rp, rp + n, a + 1, e) > 0 {
                    return None;
                }
                rp += n;
            }
            Cigar::SoftClip(n) => qp += n as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    // `rp` is now the read's exclusive reference end.
    if rp < e + INDEL_FLANK as u32 {
        return None;
    }
    Some(k)
}

#[cfg(test)]
mod indel_tests {
    use super::*;

    #[test]
    fn net_len_and_kind() {
        assert_eq!(Allele::Ref.net_len(), 0);
        assert_eq!(Allele::Ins(b"AC".to_vec()).net_len(), 2);
        assert_eq!(Allele::Del(3).net_len(), -3);
        assert!(Allele::Ins(b"A".to_vec()).kind_matches(&Allele::Ins(b"GG".to_vec())));
        assert!(!Allele::Ins(b"A".to_vec()).kind_matches(&Allele::Del(1)));
        assert!(!Allele::Ref.kind_matches(&Allele::Del(1)));
    }

    #[test]
    fn indel_opts_default_is_disabled() {
        let o = IndelOpts::default();
        assert!(!o.enabled, "indel correction must be opt-in");
        assert!(!o.hp_length_model, "the HP length model must be opt-in");
        assert_eq!(o.max_len, 5);
    }

    // ref:  0123456789.
    const HP_REF: &[u8] = b"CAAAAAAAG";

    #[test]
    fn deletion_in_homopolymer_normalizes_to_one_anchor() {
        // A 1bp deletion of an A can be anchored anywhere in the run (1..=6); every placement must collapse to anchor 0.
        for anchor in 0..=6u32 {
            let (pos, allele) = normalize_left(HP_REF, anchor, &Allele::Del(1));
            assert_eq!(pos, 0, "anchor {anchor} did not normalize to 0");
            assert_eq!(allele, Allele::Del(1));
        }
    }

    #[test]
    fn insertion_in_homopolymer_normalizes_to_one_anchor() {
        // Inserting an extra A is legal after any A in the run; all collapse to 0.
        for anchor in 0..=7u32 {
            let (pos, allele) = normalize_left(HP_REF, anchor, &Allele::Ins(b"A".to_vec()));
            assert_eq!(pos, 0, "anchor {anchor} did not normalize to 0");
            assert_eq!(allele, Allele::Ins(b"A".to_vec()));
        }
    }

    #[test]
    fn insertion_rotates_when_shifting_left() {
        // ref: C A T A T A T G; inserting "TA" after index 1 is the same event as inserting "AT" after index 0.
        let r = b"CATATATG";
        let (pos, allele) = normalize_left(r, 1, &Allele::Ins(b"TA".to_vec()));
        assert_eq!(pos, 0);
        assert_eq!(allele, Allele::Ins(b"AT".to_vec()));
    }

    #[test]
    fn unique_sequence_event_does_not_move() {
        // ref: A C G T A C G T ; deleting the G at index 2 is already left-aligned.
        let r = b"ACGTACGT";
        let (pos, allele) = normalize_left(r, 1, &Allele::Del(1));
        assert_eq!(pos, 1);
        assert_eq!(allele, Allele::Del(1));
    }

    #[test]
    fn normalization_stops_at_contig_start() {
        // Whole reference is one homopolymer: shifting can never pass index 0.
        let r = b"AAAAA";
        let (pos, _) = normalize_left(r, 3, &Allele::Del(1));
        assert_eq!(pos, 0);
    }

    #[test]
    fn ref_allele_never_moves() {
        let (pos, allele) = normalize_left(HP_REF, 4, &Allele::Ref);
        assert_eq!(pos, 4);
        assert_eq!(allele, Allele::Ref);
    }

    #[test]
    fn repeat_context_counts_homopolymer_run() {
        // HP_REF = C AAAAAAA G -> a 1bp A event normalized to anchor 0 sits in a 7-run.
        assert_eq!(repeat_context(HP_REF, 0, &Allele::Del(1)), 7);
        assert_eq!(repeat_context(HP_REF, 0, &Allele::Ins(b"A".to_vec())), 7);
    }

    #[test]
    fn hp_tract_end_spans_the_whole_run() {
        // Anchor 7 is the base before the A-run; the run is 8..=15, so its exclusive end is 16.
        assert_eq!(hp_tract_end(REF29, 7), 16);
        // From inside the run the remaining tract still ends at 16.
        assert_eq!(hp_tract_end(REF29, 11), 16);
    }

    #[test]
    fn hp_tract_end_is_anchor_plus_one_without_a_run() {
        // refseq[21] == 'G', refseq[22] == 'T': a lone base, not a run.
        assert_eq!(hp_tract_end(REF29, 20), 21);
        // Past the contig end there is nothing to span either.
        assert_eq!(hp_tract_end(REF29, 28), 29);
    }

    #[test]
    fn repeat_context_counts_tandem_copies() {
        // ref: C ATATATAT G -> a 2bp "AT" event at anchor 0 sits in 4 copies.
        let r = b"CATATATATG";
        assert_eq!(repeat_context(r, 0, &Allele::Ins(b"AT".to_vec())), 4);
        assert_eq!(repeat_context(r, 0, &Allele::Del(2)), 4);
    }

    #[test]
    fn repeat_context_is_one_in_unique_sequence() {
        let r = b"ACGTACGT";
        // Deleting the single G at index 2 (anchor 1): "G" occurs once here.
        assert_eq!(repeat_context(r, 1, &Allele::Del(1)), 1);
        assert_eq!(repeat_context(r, 1, &Allele::Ref), 1);
    }

    #[test]
    fn repeat_context_uses_primitive_period_for_multibase_hp_events() {
        let r = b"CAAAAAAAAG"; // 8 A
        for a in [Allele::Del(1), Allele::Del(2), Allele::Del(3),
                  Allele::Ins(b"A".to_vec()), Allele::Ins(b"AA".to_vec())] {
            assert_eq!(repeat_context(r, 0, &a), 8, "{a:?}");
        }
    }

    #[test]
    fn repeat_context_reduces_tandem_multiples_to_their_period() {
        let r = b"CATATATATG"; // (AT)x4
        assert_eq!(repeat_context(r, 0, &Allele::Ins(b"ATAT".to_vec())), 4);
        assert_eq!(repeat_context(r, 0, &Allele::Del(4)), 4);
    }

    #[test]
    fn error_rate_grows_with_context_and_is_capped() {
        let o = IndelOpts::default(); // err0 0.01, scale 1.5, cap 0.4
        assert!((error_rate(1, &o) - 0.01).abs() < 1e-12);
        assert!((error_rate(2, &o) - 0.015).abs() < 1e-12);
        assert!(error_rate(3, &o) > error_rate(2, &o));
        // Monotonic non-decreasing, never above the cap.
        for l in 1..40u32 {
            assert!(error_rate(l + 1, &o) >= error_rate(l, &o), "not monotonic at {l}");
            assert!(error_rate(l, &o) <= o.err_cap + 1e-12, "cap exceeded at {l}");
        }
        // Deep enough runs saturate exactly at the cap.
        assert!((error_rate(30, &o) - o.err_cap).abs() < 1e-12);
    }

    #[test]
    fn vaf_floor_is_the_larger_of_absolute_and_error_scaled() {
        let o = IndelOpts::default(); // vaf 0.05, floor_mult 3.0
        // L=1: 3 * 0.01 = 0.03 < 0.05, so the absolute floor wins.
        assert!((vaf_floor(1, &o) - 0.05).abs() < 1e-12);
        // Deep homopolymer: 3 * 0.4 = 1.2 dominates.
        assert!(vaf_floor(30, &o) > 1.0);
        for l in 1..30u32 {
            assert!(vaf_floor(l + 1, &o) >= vaf_floor(l, &o), "not monotonic at {l}");
        }
    }

    // Build a site: `events` is (allele, fwd_count, rev_count).
    fn site(depth: u32, events: &[(Allele, u32, u32)]) -> IndelSite {
        let mut s = IndelSite { depth, ..Default::default() };
        let mut used_f = 0;
        let mut used_r = 0;
        for (a, f, r) in events {
            for _ in 0..*f { s.add(a.clone(), false); }
            for _ in 0..*r { s.add(a.clone(), true); }
            used_f += f;
            used_r += r;
        }
        // Remaining depth votes REF, strand-balanced.
        let rest = depth - (used_f + used_r);
        s.col_fwd = used_f + rest / 2;
        s.col_rev = used_r + (rest - rest / 2);
        s
    }

    const UNIQ_REF: &[u8] = b"ACGTACGTACGTACGT";
    const SB_P: f64 = 0.01;
    const HOM_VAF: f64 = 0.7;

    fn kept_of<'a>(m: &'a IndelSiteModel, a: &Allele) -> Option<&'a (Allele, f64)> {
        m.kept.iter().find(|(k, _)| k == a)
    }

    #[test]
    fn balanced_indel_in_unique_sequence_is_kept() {
        // 20% 1bp deletion, strand-balanced, in unique sequence (L=1, floor 0.05).
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 10, 10)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_some(), "balanced 20% del in unique seq must be kept");
        assert!(kept_of(&m, &Allele::Ref).is_some(), "REF is always kept");
    }

    #[test]
    fn same_indel_inside_long_homopolymer_is_rejected() {
        // Identical evidence, but now in a 7-mer A-run: floor = 3 * eps(7) which is far above 0.20.
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 10, 10)]);
        let o = IndelOpts::default();
        assert!(vaf_floor(7, &o) > 0.20, "test premise: floor must exceed the 20% VAF");
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none(), "20% del in a 7-mer run must be rejected");
    }

    #[test]
    fn strand_skewed_indel_is_rejected_and_counted() {
        // 20 alt reads, 18/2 split, against a balanced column -> strand-biased.
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 18, 2)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none());
        assert_eq!(m.strand_biased, 1);
    }

    #[test]
    fn near_homoplasmic_indel_is_exempt_from_strand_bias() {
        // 80% of the column, badly skewed -- but a single-strand artifact cannot reach 80%.
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 75, 5)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_some(), "near-homoplasmic allele must survive");
        assert_eq!(m.strand_biased, 0);
    }

    #[test]
    fn deep_indel_with_small_strand_skew_is_kept() {
        let del = Allele::Del(1);
        let mut s = site(5240, &[(del.clone(), 328, 196)]);
        s.col_fwd = 2693;
        s.col_rev = 2547;
        assert!(
            crate::call::strand_bias_pvalue(328, 196, 2693.0 / 5240.0) < SB_P,
            "test premise: the p-value alone would reject this allele"
        );
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_some(), "a 0.11 strand effect must not erase a real indel");
        assert_eq!(m.strand_biased, 0);
    }

    #[test]
    fn deep_single_strand_indel_artifact_is_still_rejected() {
        let del = Allele::Del(1);
        let mut s = site(4640, &[(del.clone(), 449, 15)]);
        s.col_fwd = 2171;
        s.col_rev = 2469;
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none(), "a 0.50 strand effect is still an artifact");
        assert_eq!(m.strand_biased, 1);
    }

    #[test]
    fn low_count_indel_fails_the_per_strand_gate() {
        // Only one read on the reverse strand: fails min_strand=2 regardless of VAF.
        let ins = Allele::Ins(b"T".to_vec());
        let s = site(100, &[(ins.clone(), 19, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &ins).is_none());
    }

    const C5_REF: &[u8] = b"TACCCCCTCTA";

    #[test]
    fn ten_percent_deletion_in_c5_run_at_depth_431_survives() {
        let del = Allele::Del(1);
        assert_eq!(repeat_context(C5_REF, 1, &del), 5, "test premise: L=5");
        let o = IndelOpts::default();
        assert!(vaf_floor(5, &o) > 0.15, "test premise: the old floor rejects 10%");
        let s = site(431, &[(del.clone(), 22, 21)]);
        let m = fit_indel_site(&s, C5_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_some(), "10% at depth 431 in a 5-C run must survive");
    }

    #[test]
    fn six_percent_deletion_in_c5_run_at_depth_60_is_rejected() {
        let del = Allele::Del(1);
        let s = site(60, &[(del.clone(), 2, 2)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, C5_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none(), "6% at depth 60 in a 5-C run must be rejected");
    }

    #[test]
    fn ten_percent_deletion_in_c5_run_at_depth_60_is_not_in_excess() {
        let del = Allele::Del(1);
        let s = site(60, &[(del.clone(), 3, 3)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, C5_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none(), "a non-significant excess must be rejected");
    }

    #[test]
    fn six_percent_deletion_in_unique_sequence_at_depth_60_is_rejected() {
        let del = Allele::Del(1);
        let s = site(60, &[(del.clone(), 2, 2)]);
        let o = IndelOpts::default();
        assert!(4.0 / 60.0 >= vaf_floor(1, &o), "test premise: the old floor kept it");
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none(), "6% at depth 60 must be rejected");
    }

    #[test]
    fn excess_test_is_bonferroni_corrected_over_n_sites() {
        let del = Allele::Del(1);
        let o = IndelOpts::default();
        let s = site(431, &[(del.clone(), 22, 21)]);
        let m = fit_indel_site(&s, C5_REF, 1, 2, SB_P, HOM_VAF, 16569, &o);
        assert!(kept_of(&m, &del).is_none(), "43/431 must not survive 16569 tests");
        let s = site(5334, &[(del.clone(), 210, 209)]);
        let m = fit_indel_site(&s, C5_REF, 1, 2, SB_P, HOM_VAF, 16569, &o);
        assert!(kept_of(&m, &del).is_some(), "419/5334 must survive 16569 tests");
    }

    #[test]
    fn kept_and_obs_vaf_are_returned_in_deterministic_sorted_order() {
        // Three alleles survive candidacy in unique sequence (L=1, floor 0.05): REF, an insertion, and a deletion.
        let ins = Allele::Ins(b"A".to_vec());
        let del = Allele::Del(1);
        let s = site(100, &[(ins.clone(), 10, 10), (del.clone(), 10, 10)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert_eq!(m.kept.len(), 3, "REF, Ins, and Del must all survive candidacy here");
        // Assert the EXACT expected sequence rather than comparing the list to a sorted copy of itself.
        let alleles: Vec<&Allele> = m.kept.iter().map(|(a, _)| a).collect();
        assert_eq!(
            alleles,
            vec![&Allele::Ref, &Allele::Ins(b"A".to_vec()), &Allele::Del(1)],
            "kept must be in deterministic sorted order"
        );
        let obs_alleles: Vec<&Allele> = m.obs_vaf.iter().map(|(a, _)| a).collect();
        assert_eq!(
            obs_alleles,
            vec![&Allele::Ref, &Allele::Ins(b"A".to_vec()), &Allele::Del(1)],
            "obs_vaf must be in deterministic sorted order"
        );
    }

    // ---- assignment (the three rows of the design's table) ----

    #[test]
    fn heteroplasmic_site_preserves_both_alleles() {
        // 40% INS / 60% REF in a long homopolymer (large eps).
        let ins = Allele::Ins(b"A".to_vec());
        let s = site(100, &[(ins.clone(), 20, 20)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0; // isolate assignment from the gate for this test
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &ins).is_some() && kept_of(&m, &Allele::Ref).is_some());
        assert_eq!(assign_allele(&ins, &m, &o), Assignment::Keep, "40% allele must not flip");
        assert_eq!(assign_allele(&Allele::Ref, &m, &o), Assignment::Keep, "60% allele must not flip");
    }

    #[test]
    fn near_homoplasmic_site_repairs_reads_that_dropped_the_indel() {
        // 90% INS / 10% REF in a long homopolymer: REF reads are dropped-indel errors.
        let ins = Allele::Ins(b"A".to_vec());
        let s = site(100, &[(ins.clone(), 45, 45)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(m.eps > 0.1, "test premise: long homopolymer must have a large eps");
        assert_eq!(assign_allele(&Allele::Ref, &m, &o), Assignment::Move(ins.clone()));
        assert_eq!(assign_allele(&ins, &m, &o), Assignment::Keep);
    }

    #[test]
    fn failed_allele_snaps_to_the_kept_allele_of_the_same_kind() {
        // Site keeps a 1bp deletion at 80%; a stray 2bp deletion fails candidacy and must snap onto the kept 1bp allele.
        let d1 = Allele::Del(1);
        let d2 = Allele::Del(2);
        let s = site(100, &[(d1.clone(), 40, 40), (d2.clone(), 1, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &d1).is_some());
        assert!(kept_of(&m, &d2).is_none());
        assert_eq!(assign_allele(&d2, &m, &o), Assignment::Move(d1.clone()));
    }

    #[test]
    fn cross_kind_reassignment_is_refused() {
        // Site keeps only an insertion; a failed deletion would MAP across kinds.
        let ins = Allele::Ins(b"T".to_vec());
        let del = Allele::Del(1);
        let s = site(100, &[(ins.clone(), 45, 45), (del.clone(), 1, 1)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none());
        // The MAP target here is Ins, which is a cross-kind transition the rewrite walk does not support.
        assert_eq!(assign_allele(&del, &m, &o), Assignment::Refused);
    }

    #[test]
    fn growing_a_deletion_is_refused() {
        // Site keeps a 3bp deletion; a read carrying a 1bp deletion would have to grow.
        let d1 = Allele::Del(1);
        let d3 = Allele::Del(3);
        let s = site(100, &[(d3.clone(), 45, 45), (d1.clone(), 1, 1)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        // d3 dominates the site and is the MAP target.
        assert_eq!(assign_allele(&d1, &m, &o), Assignment::Refused);
    }

    // ---- Plan 2026-09-26 T8c: at an HP site every kept length is reachable ----

    #[test]
    fn hp_length_model_moves_an_insertion_onto_a_kept_deletion() {
        // A site that keeps Del(1) at ~0.9 with a lone Ins("A") read.
        let ins = Allele::Ins(b"A".to_vec());
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 45, 45), (ins.clone(), 1, 1)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;

        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(m.hp_run.is_none(), "test premise: the model is off by default");
        assert!(kept_of(&m, &del).is_some_and(|(_, f)| (*f - 0.9).abs() < 0.03));
        assert_eq!(assign_allele(&ins, &m, &o), Assignment::Refused);

        o.hp_length_model = true;
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert_eq!(m.hp_run, Some((8, b'A')), "anchor 0 is the A-run's HP site");
        assert_eq!(assign_allele(&ins, &m, &o), Assignment::Move(del));
    }

    #[test]
    fn hp_length_model_can_grow_a_deletion() {
        // The second half of "any kept length is reachable": Del(1) -> Del(3) is the transition `rewrite_read` cannot make.
        let d1 = Allele::Del(1);
        let d3 = Allele::Del(3);
        let s = site(100, &[(d3.clone(), 45, 45), (d1.clone(), 1, 1)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;

        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert_eq!(assign_allele(&d1, &m, &o), Assignment::Refused);

        o.hp_length_model = true;
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert_eq!(assign_allele(&d1, &m, &o), Assignment::Move(d3));
    }

    #[test]
    fn hp_length_model_leaves_protect_vaf_and_the_strand_gates_alone() {
        // The HP length model widens which TRANSITIONS are reachable and nothing else.
        let ins = Allele::Ins(b"A".to_vec());
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 30, 30), (ins.clone(), 15, 15)]);
        let mut o = IndelOpts::default();
        o.hp_length_model = true;
        let m = fit_indel_site(&s, DEEP_HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(m.hp_run.is_some(), "test premise: this must be an HP site");
        assert!(kept_of(&m, &ins).is_none(), "test premise: the insertion loses candidacy");
        assert!(
            !m.strand_rejected.contains(&ins),
            "test premise: it must be rejected by the floor, not a strand gate"
        );
        assert_eq!(
            assign_allele(&ins, &m, &o),
            Assignment::Keep,
            "a 30% allele must stay protected at an HP site too"
        );

        // Same allele, now single-stranded: rejected by `min_strand`.
        let s = site(100, &[(del.clone(), 30, 30), (ins.clone(), 30, 0)]);
        let m = fit_indel_site(&s, DEEP_HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(m.strand_rejected.contains(&ins), "test premise: strand-gate rejection");
        assert!(
            matches!(assign_allele(&ins, &m, &o), Assignment::Move(_)),
            "a strand-rejected allele must not be shielded by protect_vaf"
        );
    }

    #[test]
    fn lone_spurious_indel_at_a_ref_only_site_is_still_corrected() {
        // THE most common real case: a single read carries a random indel error.
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 1, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert_eq!(m.kept.len(), 1, "only REF should survive");
        // A read already carrying the sole kept allele is untouched...
        assert_eq!(assign_allele(&Allele::Ref, &m, &o), Assignment::Keep);
        // ...but the stray deletion is reverted.
        assert_eq!(assign_allele(&del, &m, &o), Assignment::Move(Allele::Ref));
        // Same for a stray insertion nobody else carries.
        assert_eq!(
            assign_allele(&Allele::Ins(b"T".to_vec()), &m, &o),
            Assignment::Move(Allele::Ref)
        );
    }

    // ref: C AAAAAAAAAAAA G -> a 1bp event here sits in L=12.
    const DEEP_HP_REF: &[u8] = b"CAAAAAAAAAAAAG";

    #[test]
    fn substantial_allele_erased_from_candidacy_is_still_protected_at_assignment() {
        // A genuine 50/50 Del(1) heteroplasmy in a 12-mer homopolymer.
        let del = Allele::Del(1);
        let o = IndelOpts::default();
        assert!(vaf_floor(12, &o) > 1.0, "test premise: floor must exceed 1.0 at L=12");
        let s = site(100, &[(del.clone(), 25, 25)]);
        let m = fit_indel_site(&s, DEEP_HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none(), "candidacy must genuinely reject the allele");
        assert_eq!(
            assign_allele(&del, &m, &o),
            Assignment::Keep,
            "a substantial allele must be protected from erasure, not destroyed"
        );
    }

    #[test]
    fn low_frequency_indel_below_protect_vaf_is_still_reverted() {
        // Same deep-homopolymer context, but now the indel is a rare ~5% event, well below the default `protect_vaf` of 0.2.
        let del = Allele::Del(1);
        let o = IndelOpts::default();
        let s = site(100, &[(del.clone(), 3, 2)]); // 5/100 = 0.05 < protect_vaf 0.2
        let m = fit_indel_site(&s, DEEP_HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none(), "candidacy must reject the allele");
        assert_eq!(
            assign_allele(&del, &m, &o),
            Assignment::Move(Allele::Ref),
            "a low-frequency indel below protect_vaf must still be reverted"
        );
    }

    #[test]
    fn strand_rejected_allele_is_not_shielded_by_protect_vaf() {
        // 30% VAF in unique sequence, but entirely on the forward strand: a classic single-strand artifact.
        let ins = Allele::Ins(b"T".to_vec());
        let s = site(100, &[(ins.clone(), 30, 0)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &ins).is_none(), "candidacy must reject the strand-skewed allele");
        assert!(
            m.strand_rejected.contains(&ins),
            "allele failing min_strand must be recorded in strand_rejected"
        );
        assert_eq!(
            assign_allele(&ins, &m, &o),
            Assignment::Move(Allele::Ref),
            "protect_vaf must not shield an allele rejected by the strand-count gate"
        );
    }

    #[test]
    fn floor_rejected_allele_in_deep_homopolymer_is_still_shielded_by_protect_vaf() {
        // Same 30% VAF, but now strand-balanced in a deep homopolymer where `vaf_floor` exceeds 1.0.
        let del = Allele::Del(1);
        let o = IndelOpts::default();
        assert!(vaf_floor(12, &o) > 1.0, "test premise: floor must exceed 1.0 at L=12");
        let s = site(100, &[(del.clone(), 15, 15)]);
        let m = fit_indel_site(&s, DEEP_HP_REF, 0, 2, SB_P, HOM_VAF, 1, &o);
        assert!(kept_of(&m, &del).is_none(), "candidacy must reject the allele on the floor gate");
        assert!(
            !m.strand_rejected.contains(&del),
            "an allele rejected only by the context floor must not be in strand_rejected"
        );
        assert_eq!(
            assign_allele(&del, &m, &o),
            Assignment::Keep,
            "protect_vaf must still shield an allele rejected only by the floor"
        );
    }

    // Reference shared by the read_events tests.
    const REF29: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC";

    fn cs(ops: Vec<Cigar>) -> CigarString {
        CigarString(ops)
    }

    #[test]
    fn match_run_len_reports_remaining_run_inside_a_match() {
        // 10M starting at pos 0: from target 3, 7 match bases remain (3..10).
        let c = cs(vec![Cigar::Match(10)]);
        assert_eq!(match_run_len(0, &c, 3), 7);
        // At the very first base of the run, the whole run remains.
        assert_eq!(match_run_len(0, &c, 0), 10);
        // At the last base, exactly 1 remains.
        assert_eq!(match_run_len(0, &c, 9), 1);
    }

    #[test]
    fn match_run_len_is_zero_inside_a_deletion_or_refskip() {
        // 5M 3D 5M: targets 5..8 fall inside the deletion, and have no aligned query base at all.
        let c = cs(vec![Cigar::Match(5), Cigar::Del(3), Cigar::Match(5)]);
        for t in 5..8 {
            assert_eq!(match_run_len(0, &c, t), 0, "target {t} is inside the deletion");
        }
        // The match run resumes after the deletion.
        assert_eq!(match_run_len(0, &c, 8), 5);

        let n = cs(vec![Cigar::Match(5), Cigar::RefSkip(3), Cigar::Match(5)]);
        for t in 5..8 {
            assert_eq!(match_run_len(0, &n, t), 0, "target {t} is inside the refskip");
        }
    }

    #[test]
    fn match_run_len_is_zero_before_start_or_past_end() {
        let c = cs(vec![Cigar::Match(10)]);
        assert_eq!(match_run_len(0, &c, 10), 0, "target at end_pos is past the read");
        assert_eq!(match_run_len(5, &c, 3), 0, "target before the read's own POS");
    }

    #[test]
    fn match_run_len_skips_insertions_without_advancing_reference() {
        // An insertion consumes no reference, so the two Match ops it splits are still contiguous in REFERENCE space.
        let c = cs(vec![Cigar::Match(5), Cigar::Ins(2), Cigar::Match(5)]);
        assert_eq!(match_run_len(0, &c, 5), 5, "second run's own length at its start");
        assert_eq!(match_run_len(0, &c, 4), 1, "first run's own remaining length");
    }

    #[test]
    fn read_with_no_indels_yields_no_events() {
        let evs = read_events(0, &cs(vec![Cigar::Match(29)]), REF29, REF29, 5);
        assert!(evs.is_empty());
    }

    #[test]
    fn insertion_event_is_normalized_out_of_the_homopolymer() {
        // Read carries an extra A, anchored by the aligner at reference index 15 (the last A of the run).
        let seq = b"ACGTACGTAAAAAAAAACGTACGTACGTAC"; // 30bp: 16 + 1 inserted + 13
        let cig = cs(vec![Cigar::Match(16), Cigar::Ins(1), Cigar::Match(13)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        // The contrast this test exists for: the aligner's own CIGAR anchor (15) and the left-normalized site genuinely differ.
        assert_eq!(evs[0].anchor, 15, "the aligner's own CIGAR anchor must be reported");
        assert_eq!(evs[0].norm_pos, 7, "left-normalization must still pull the site to 7");
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"A".to_vec())));
    }

    #[test]
    fn deletion_event_is_normalized_out_of_the_homopolymer() {
        // Read is missing one A from the run; the aligner anchored the deletion at reference index 14.
        let seq = b"ACGTACGTAAAAAAACGTACGTACGTAC"; // 28bp: 29 ref bases minus one A
        let cig = cs(vec![Cigar::Match(15), Cigar::Del(1), Cigar::Match(13)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].anchor, 14, "the aligner's own CIGAR anchor must be reported");
        assert_eq!(evs[0].norm_pos, 7, "left-normalization must still pull the site to 7");
        assert_eq!(evs[0].allele, Some(Allele::Del(1)));
    }

    #[test]
    fn unique_sequence_insertion_keeps_its_anchor() {
        // Insert "G" after reference index 18 (a T).
        let seq = b"ACGTACGTAAAAAAAACGTGACGTACGTAC"; // 30bp
        let cig = cs(vec![Cigar::Match(19), Cigar::Ins(1), Cigar::Match(10)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        // In unique sequence anchor and normalized site coincide -- unlike the homopolymer fixtures above.
        assert_eq!(evs[0].anchor, 18);
        assert_eq!(evs[0].norm_pos, 18);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"G".to_vec())));
    }

    #[test]
    fn oversized_indel_is_reported_as_untouchable() {
        // A 6bp insertion with max_len 5: recorded with allele None.
        let seq = b"ACGTACGTAAAAAAAACGTGGGGGGACGTACGTAC"; // 35bp
        let cig = cs(vec![Cigar::Match(19), Cigar::Ins(6), Cigar::Match(10)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].allele, None);
        // The oversized case reports the same raw anchor in both fields.
        assert_eq!(evs[0].anchor, 18);
        assert_eq!(evs[0].norm_pos, 18);
    }

    #[test]
    fn multiple_events_in_one_read_are_all_reported() {
        // An extra A in the run (normalizes to 7) plus a 1bp deletion at index 19.
        let seq = b"ACGTACGTAAAAAAAAACGTCGTACGTAC"; // 29bp
        let cig = cs(vec![
            Cigar::Match(16),
            Cigar::Ins(1),
            Cigar::Match(4),
            Cigar::Del(1),
            Cigar::Match(8),
        ]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 2);
        assert_eq!(evs[0].norm_pos, 7);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"A".to_vec())));
        assert_eq!(evs[1].norm_pos, 19);
        assert_eq!(evs[1].allele, Some(Allele::Del(1)));
    }

    #[test]
    fn soft_clips_do_not_shift_reference_coordinates() {
        // Leading soft clip consumes query but no reference; the event must still land at reference index 18.
        let seq = b"TTTACGTACGTAAAAAAAACGTGACGTACGTAC"; // 3 clipped + 30
        let cig = cs(vec![
            Cigar::SoftClip(3),
            Cigar::Match(19),
            Cigar::Ins(1),
            Cigar::Match(10),
        ]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].norm_pos, 18);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"G".to_vec())));
    }

    #[test]
    fn read_aligned_at_nonzero_pos_reports_absolute_coordinates() {
        // Read starts at reference index 16; insert "G" after index 18.
        let seq = b"CGTGACGTACGTAC"; // ref[16..19] + "G" + ref[19..29]
        let cig = cs(vec![Cigar::Match(3), Cigar::Ins(1), Cigar::Match(10)]);
        let evs = read_events(16, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].norm_pos, 18);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"G".to_vec())));
    }

    // ---- an insertion right after a D/N anchors at the

    #[test]
    fn insertion_after_a_deletion_anchors_at_the_last_aligned_base_not_the_deleted_one() {
        // `20M 2D 1I 20M` at POS 100.
        let refseq = vec![b'A'; 150];
        let seq = vec![b'A'; 41]; // 20 (M) + 1 (I) + 20 (M); the D consumes no query
        let cig = cs(vec![
            Cigar::Match(20),
            Cigar::Del(2),
            Cigar::Ins(1),
            Cigar::Match(20),
        ]);
        let evs = read_events(100, &cig, &seq, &refseq, 5);
        assert_eq!(evs.len(), 2, "both the deletion and the insertion must be reported");
        assert_eq!(evs[0].allele, Some(Allele::Del(2)));
        assert_eq!(evs[0].anchor, 119, "the deletion's own anchor is the last aligned base");
        assert!(matches!(evs[1].allele, Some(Allele::Ins(_))));
        assert_eq!(
            evs[1].anchor, 119,
            "the insertion must anchor at the last ALIGNED base (119), not the last \
             DELETED one (121)"
        );
    }

    // ---- Plan 2026-09-26 T8b: HP site recognition ----

    #[test]
    fn hp_site_is_the_runs_leftmost_anchor_only() {
        // REF29's A-run is 8..=15.
        assert_eq!(hp_site(REF29, 7), Some((16, b'A')));
        for anchor in 8..=14u32 {
            assert_eq!(hp_site(REF29, anchor), None, "anchor {anchor} is not normalized");
        }
    }

    #[test]
    fn hp_site_rejects_unique_sequence_and_lone_bases() {
        // A single following base is not a run, and the last base of the run has only unique sequence after it.
        assert_eq!(hp_site(REF29, 15), None, "REF29[16] = C is a lone base");
        assert_eq!(hp_site(REF29, 20), None);
        assert_eq!(hp_site(REF29, 28), None, "nothing follows the last base");
    }

    #[test]
    fn hp_length_allele_maps_a_length_change_to_an_allele() {
        assert_eq!(hp_length_allele(0, b'A'), Allele::Ref);
        assert_eq!(hp_length_allele(2, b'A'), Allele::Ins(b"AA".to_vec()));
        assert_eq!(hp_length_allele(-3, b'A'), Allele::Del(3));
    }

    // ---- Plan 2026-09-26 T8a: homopolymer-length allele ----
    const HP_A: u32 = 7;
    const HP_E: u32 = 16;

    #[test]
    fn hp_allele_of_a_pure_match_read_is_zero() {
        let k = hp_allele(0, &cs(vec![Cigar::Match(29)]), REF29, REF29, HP_A, HP_E);
        assert_eq!(k, Some(0), "a read matching the whole run carries the reference length");
    }

    #[test]
    fn hp_allele_counts_an_insertion_of_the_run_base() {
        // 16M 1I 13M: the extra A is anchored on REF29[15], the last base of the run.
        let mut seq = REF29[..16].to_vec();
        seq.push(b'A');
        seq.extend_from_slice(&REF29[16..]);
        let cig = cs(vec![Cigar::Match(16), Cigar::Ins(1), Cigar::Match(13)]);
        assert_eq!(hp_allele(0, &cig, &seq, REF29, HP_A, HP_E), Some(1));
    }

    #[test]
    fn hp_allele_counts_a_deletion_inside_the_run() {
        // 10M 2D 17M: REF29[10] and REF29[11] are deleted, both inside the run.
        let mut seq = REF29[..10].to_vec();
        seq.extend_from_slice(&REF29[12..]);
        let cig = cs(vec![Cigar::Match(10), Cigar::Del(2), Cigar::Match(17)]);
        assert_eq!(hp_allele(0, &cig, &seq, REF29, HP_A, HP_E), Some(-2));
    }

    #[test]
    fn hp_allele_of_an_offsetting_insertion_and_deletion_is_zero() {
        // 10M 1I 2M 1D 16M: one A gained and one A lost, both inside the run, so the run's LENGTH is unchanged.
        let mut seq = REF29[..10].to_vec();
        seq.push(b'A');
        seq.extend_from_slice(&REF29[10..12]);
        seq.extend_from_slice(&REF29[13..]);
        let cig = cs(vec![
            Cigar::Match(10),
            Cigar::Ins(1),
            Cigar::Match(2),
            Cigar::Del(1),
            Cigar::Match(16),
        ]);
        assert_eq!(hp_allele(0, &cig, &seq, REF29, HP_A, HP_E), Some(0));
    }

    #[test]
    fn hp_allele_is_undefined_when_the_inserted_base_is_not_the_run_base() {
        // A G spliced into an A-run is not a length change at all -- it breaks the run.
        let mut seq = REF29[..16].to_vec();
        seq.push(b'G');
        seq.extend_from_slice(&REF29[16..]);
        let cig = cs(vec![Cigar::Match(16), Cigar::Ins(1), Cigar::Match(13)]);
        assert_eq!(hp_allele(0, &cig, &seq, REF29, HP_A, HP_E), None);
    }

    #[test]
    fn hp_allele_is_undefined_for_a_read_ending_inside_the_run() {
        // Ends at reference 13, five bases into the run: it cannot testify to the run's length.
        let cig = cs(vec![Cigar::Match(13)]);
        assert_eq!(hp_allele(0, &cig, &REF29[..13], REF29, HP_A, HP_E), None);
    }

    #[test]
    fn hp_allele_is_undefined_without_the_right_hand_flank() {
        // Reaches the end of the run (16) but not `INDEL_FLANK` past it, so it has no match context to anchor a rewrite against.
        let need = HP_E as usize + INDEL_FLANK;
        let cig = cs(vec![Cigar::Match(need as u32 - 1)]);
        assert_eq!(hp_allele(0, &cig, &REF29[..need - 1], REF29, HP_A, HP_E), None);
        let cig = cs(vec![Cigar::Match(need as u32)]);
        assert_eq!(hp_allele(0, &cig, &REF29[..need], REF29, HP_A, HP_E), Some(0));
    }

    #[test]
    fn hp_allele_is_undefined_when_the_anchor_base_is_not_aligned() {
        // 5M 5D 19M deletes REF29[5..10), which swallows the anchor: there is no aligned base to hang the canonical rewrite off.
        let mut seq = REF29[..5].to_vec();
        seq.extend_from_slice(&REF29[10..]);
        let cig = cs(vec![Cigar::Match(5), Cigar::Del(5), Cigar::Match(19)]);
        assert_eq!(hp_allele(0, &cig, &seq, REF29, HP_A, HP_E), None);
    }

    #[test]
    fn hp_allele_is_undefined_when_an_n_op_crosses_the_run() {
        let mut seq = REF29[..10].to_vec();
        seq.extend_from_slice(&REF29[12..]);
        let cig = cs(vec![Cigar::Match(10), Cigar::RefSkip(2), Cigar::Match(17)]);
        assert_eq!(hp_allele(0, &cig, &seq, REF29, HP_A, HP_E), None);
    }
}

// ===========================================================================
/// One read's reassignment at one normalized indel site.
#[derive(Debug, Clone)]
pub struct IndelEdit {
    /// The coordinate `rewrite_read` looks this edit up by.
    pub ref_pos: u32,
    pub from: Allele,
    pub to: Allele,
}

/// Which of the three edit varieties an `IndelEdit` is, discriminated by `from`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum EditKind {
    Ref,
    Ins,
    Del,
}

pub(crate) fn edit_kind(a: &Allele) -> EditKind {
    match a {
        Allele::Ref => EditKind::Ref,
        Allele::Ins(_) => EditKind::Ins,
        Allele::Del(_) => EditKind::Del,
    }
}

/// All edits for one read, keyed by REFERENCE coordinate so substitutions.
#[derive(Debug, Default, Clone)]
pub struct ReadEdits {
    pub subs: Vec<(u32, u8)>,
    pub indels: Vec<IndelEdit>,
}

impl ReadEdits {
    pub fn is_empty(&self) -> bool {
        self.subs.is_empty() && self.indels.is_empty()
    }
}

pub struct RewriteResult {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub cigar: CigarString,
    /// `true` iff at least one indel edit was actually applied during the walk.
    /// A substitution alone never sets this.
    pub structure_changed: bool,
    /// Keys of the edits this walk ACTUALLY applied -- gain/revert/replace `Ins`.
    /// A verbatim copy or an unsupported-target fallthrough never adds its key here.
    pub applied: HashSet<(u32, EditKind)>,
}

/// Total reference bases consumed by a CIGAR (M/=/X/D/N).
/// The rewrite invariant. Test-only: the rewrite path maintains this invariant structurally, so
/// nothing in production recomputes it.
#[cfg(test)]
fn ref_consumed(c: &CigarString) -> u32 {
    c.iter()
        .map(|op| match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) | Cigar::RefSkip(n) => *n,
            _ => 0,
        })
        .sum()
}

fn push_op(ops: &mut Vec<Cigar>, op: Cigar) {
    let len = match op {
        Cigar::Match(n) | Cigar::Ins(n) | Cigar::Del(n) | Cigar::RefSkip(n)
        | Cigar::SoftClip(n) | Cigar::HardClip(n) | Cigar::Pad(n) | Cigar::Equal(n)
        | Cigar::Diff(n) => n,
    };
    if len == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        match (&last, &op) {
            (Cigar::Match(a), Cigar::Match(b)) => { *last = Cigar::Match(a + b); return; }
            (Cigar::Ins(a), Cigar::Ins(b)) => { *last = Cigar::Ins(a + b); return; }
            (Cigar::Del(a), Cigar::Del(b)) => { *last = Cigar::Del(a + b); return; }
            (Cigar::SoftClip(a), Cigar::SoftClip(b)) => { *last = Cigar::SoftClip(a + b); return; }
            _ => {}
        }
    }
    ops.push(op);
}

/// Rewrite one read's SEQ/QUAL/CIGAR against `edits`, which are keyed by REFERENCE coordinate.
/// Contract: `refseq` MUST be the full contig the read is aligned to, indexed by absolute reference
/// coordinate (not a windowed slice around the read).
pub fn rewrite_read(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    qual: &[u8],
    edits: &ReadEdits,
    refseq: &[u8],
) -> RewriteResult {
    if edits.is_empty() {
        // Common fast path: no edits, so nothing about the read changes.
        return RewriteResult {
            seq: seq.to_vec(),
            qual: qual.to_vec(),
            cigar: cigar.clone(),
            structure_changed: false,
            applied: HashSet::new(),
        };
    }

    if edits.indels.is_empty() {
        // Substitutions-only fast path: apply the subs directly to a copy of SEQ.
        let subs: HashMap<u32, u8> = edits.subs.iter().copied().collect();
        let mut out_seq = seq.to_vec();
        let mut rp: u32 = pos.max(0) as u32;
        let mut qp: usize = 0;
        for op in cigar.iter() {
            match *op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for k in 0..n {
                        if let Some(&b) = subs.get(&(rp + k)) {
                            out_seq[qp + k as usize] = b;
                        }
                    }
                    qp += n as usize;
                    rp += n;
                }
                Cigar::Ins(n) | Cigar::SoftClip(n) => qp += n as usize,
                Cigar::Del(n) | Cigar::RefSkip(n) => rp += n,
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        return RewriteResult {
            seq: out_seq,
            qual: qual.to_vec(),
            cigar: cigar.clone(),
            structure_changed: false,
            applied: HashSet::new(),
        };
    }

    let subs: HashMap<u32, u8> = edits.subs.iter().copied().collect();
    // Keyed by (ref_pos, kind_of(from)); on an identical key the first edit wins, the rest are reported as unapplied.
    let mut indels: HashMap<(u32, EditKind), &IndelEdit> = HashMap::new();
    for e in &edits.indels {
        let key = (e.ref_pos, edit_kind(&e.from));
        indels.entry(key).or_insert(e);
    }

    let mut out_seq: Vec<u8> = Vec::with_capacity(seq.len() + 8);
    let mut out_qual: Vec<u8> = Vec::with_capacity(qual.len() + 8);
    let mut ops: Vec<Cigar> = Vec::new();
    let mut structure_changed = false;
    // Edits genuinely applied, keyed exactly like `indels` above.
    let mut applied: HashSet<(u32, EditKind)> = HashSet::new();

    let mut rp: u32 = pos.max(0) as u32; // reference position of the next ref-consuming base
    let mut qp: usize = 0; // query position of the next original base
    // Ref position of the most recent *aligned* base emitted so far — the pileup anchor convention.
    let mut last_ref: Option<u32> = None;

    // Quality assigned to synthesized bases: never more confident than the flanks.
    let fill_qual = |out_qual: &Vec<u8>, next: Option<u8>| -> u8 {
        match (out_qual.last().copied(), next) {
            (Some(l), Some(r)) => l.min(r),
            (Some(l), None) => l,
            (None, Some(r)) => r,
            (None, None) => 0,
        }
    };

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let mut k: u32 = 0;
                while k < n {
                    let cur = rp + k;
                    let qi = qp + k as usize;
                    out_seq.push(*subs.get(&cur).unwrap_or(&seq[qi]));
                    out_qual.push(qual[qi]);
                    push_op(&mut ops, Cigar::Match(1));
                    last_ref = Some(cur);
                    k += 1;

                    // An indel anchored on this aligned base?
                    if let Some(e) = indels.get(&(cur, EditKind::Ref)) {
                        match (&e.from, &e.to) {
                            (Allele::Ref, Allele::Ins(s)) => {
                                let next = qual.get(qp + k as usize).copied();
                                let f = fill_qual(&out_qual, next);
                                out_seq.extend_from_slice(s);
                                out_qual.extend(std::iter::repeat(f).take(s.len()));
                                push_op(&mut ops, Cigar::Ins(s.len() as u32));
                                applied.insert((cur, EditKind::Ref));
                                structure_changed = true;
                            }
                            (Allele::Ref, Allele::Del(m)) => {
                                // The following m reference bases become a deletion.
                                let m = *m;
                                if k + m <= n {
                                    push_op(&mut ops, Cigar::Del(m));
                                    k += m; // skip both the ref bases and their query partners
                                    applied.insert((cur, EditKind::Ref));
                                    structure_changed = true;
                                }
                            }
                            _ => {}
                        }
                    }
                }
                qp += n as usize;
                rp += n;
            }

            Cigar::Ins(n) => {
                // While `last_ref` is `None` there is no aligned base before this insertion to anchor an edit to.
                let key = last_ref.map(|lr| (lr, EditKind::Ins));
                let e = key.and_then(|k| indels.remove(&k));
                match e.map(|e| &e.to) {
                    // Reverted: drop the inserted bases entirely.
                    Some(Allele::Ref) => {
                        applied.insert(key.expect("key set whenever `e` matched"));
                        structure_changed = true;
                    }
                    // Replaced: emit the site's consensus inserted bases instead.
                    Some(Allele::Ins(s)) => {
                        let next = qual.get(qp + n as usize).copied();
                        let f = fill_qual(&out_qual, next);
                        out_seq.extend_from_slice(s);
                        out_qual.extend(std::iter::repeat(f).take(s.len()));
                        push_op(&mut ops, Cigar::Ins(s.len() as u32));
                        applied.insert(key.expect("key set whenever `e` matched"));
                        structure_changed = true;
                    }
                    // No edit, or an unsupported cross-kind target: copy verbatim.
                    _ => {
                        out_seq.extend_from_slice(&seq[qp..qp + n as usize]);
                        out_qual.extend_from_slice(&qual[qp..qp + n as usize]);
                        push_op(&mut ops, Cigar::Ins(n));
                    }
                }
                qp += n as usize;
            }

            Cigar::Del(n) => {
                // See the `Ins` arm above: skip the lookup entirely while there is no preceding aligned base to anchor to.
                let key = last_ref.map(|lr| (lr, EditKind::Del));
                let e = key.and_then(|k| indels.remove(&k));
                let splice = |out_seq: &mut Vec<u8>, out_qual: &mut Vec<u8>, from: u32, count: u32| {
                    let f = match (out_qual.last().copied(), qual.get(qp).copied()) {
                        (Some(l), Some(r)) => l.min(r),
                        (Some(l), None) => l,
                        (None, Some(r)) => r,
                        (None, None) => 0,
                    };
                    for j in 0..count {
                        let idx = (from + j) as usize;
                        debug_assert!(
                            idx < refseq.len(),
                            "refseq index {idx} out of range (refseq.len() = {}); \
                             refseq must be the full contig indexed by absolute \
                             reference coordinate",
                            refseq.len()
                        );
                        out_seq.push(refseq.get(idx).copied().unwrap_or(b'N'));
                        out_qual.push(f);
                    }
                };
                match e.map(|e| &e.to) {
                    // Reverted: the deleted reference bases come back as matches.
                    Some(Allele::Ref) => {
                        splice(&mut out_seq, &mut out_qual, rp, n);
                        push_op(&mut ops, Cigar::Match(n));
                        applied.insert(key.expect("key set whenever `e` matched"));
                        structure_changed = true;
                    }
                    // Shrunk: left-aligned, so the shorter deletion comes first and the freed reference bases follow as matches.
                    Some(Allele::Del(m)) if *m <= n => {
                        push_op(&mut ops, Cigar::Del(*m));
                        let back = n - *m;
                        if back > 0 {
                            splice(&mut out_seq, &mut out_qual, rp + *m, back);
                            push_op(&mut ops, Cigar::Match(back));
                        }
                        applied.insert(key.expect("key set whenever `e` matched"));
                        structure_changed = true;
                    }
                    // No edit, or an unsupported target: copy verbatim.
                    _ => push_op(&mut ops, Cigar::Del(n)),
                }
                rp += n;
            }

            Cigar::SoftClip(n) => {
                out_seq.extend_from_slice(&seq[qp..qp + n as usize]);
                out_qual.extend_from_slice(&qual[qp..qp + n as usize]);
                push_op(&mut ops, Cigar::SoftClip(n));
                qp += n as usize;
            }
            Cigar::RefSkip(n) => {
                push_op(&mut ops, Cigar::RefSkip(n));
                rp += n;
            }
            Cigar::HardClip(n) => push_op(&mut ops, Cigar::HardClip(n)),
            Cigar::Pad(n) => push_op(&mut ops, Cigar::Pad(n)),
        }
    }

    RewriteResult {
        seq: out_seq,
        qual: out_qual,
        cigar: CigarString(ops),
        structure_changed,
        applied,
    }
}

/// This is the rewrite primitive the generic HP length model needs and `rewrite_read` cannot provide.
pub fn rewrite_hp_run(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    qual: &[u8],
    a: u32,
    e: u32,
    b: u8,
    target_k: i32,
) -> Option<(CigarString, Vec<u8>, Vec<u8>)> {
    if e < a + 2 {
        return None;
    }
    let run_len = (e - a - 1) as i32;
    if target_k < -run_len {
        return None;
    }
    if match_run_len(pos, cigar, a) == 0 {
        return None;
    }

    // Everything strictly before the run and everything from `e` on is copied verbatim.
    let mut pre_ops: Vec<Cigar> = Vec::new();
    let mut pre_seq: Vec<u8> = Vec::new();
    let mut pre_qual: Vec<u8> = Vec::new();
    let mut post_ops: Vec<Cigar> = Vec::new();
    let mut post_seq: Vec<u8> = Vec::new();
    let mut post_qual: Vec<u8> = Vec::new();
    // Only the run's QUALITIES are accumulated.
    let mut run_qual: Vec<u8> = Vec::new();
    // Reference bases of the run the read actually accounts for.
    let mut run_ref: u32 = 0;

    let mut rp = pos.max(0) as u32;
    let mut qp = 0usize;
    let mut last_aligned: Option<u32> = None;
    // True once a reference position at or past `e` has been reached: the run is a contiguous reference interval.
    let mut past = false;

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for i in 0..n {
                    let cur = rp + i;
                    let qi = qp + i as usize;
                    let base = *seq.get(qi)?;
                    let q = *qual.get(qi)?;
                    if cur <= a {
                        push_op(&mut pre_ops, Cigar::Match(1));
                        pre_seq.push(base);
                        pre_qual.push(q);
                    } else if cur < e {
                        run_ref += 1;
                        run_qual.push(q);
                    } else {
                        past = true;
                        push_op(&mut post_ops, Cigar::Match(1));
                        post_seq.push(base);
                        post_qual.push(q);
                    }
                    last_aligned = Some(cur);
                }
                rp += n;
                qp += n as usize;
            }
            Cigar::Ins(n) => {
                let n = n as usize;
                let bases = seq.get(qp..qp + n)?;
                let quals = qual.get(qp..qp + n)?;
                // Anchored on the anchor base or on any base of the run: part of the run's length, so absorbed.
                if last_aligned.is_some_and(|l| l >= a && l < e) {
                    run_qual.extend_from_slice(quals);
                } else if past {
                    push_op(&mut post_ops, Cigar::Ins(n as u32));
                    post_seq.extend_from_slice(bases);
                    post_qual.extend_from_slice(quals);
                } else {
                    push_op(&mut pre_ops, Cigar::Ins(n as u32));
                    pre_seq.extend_from_slice(bases);
                    pre_qual.extend_from_slice(quals);
                }
                qp += n;
            }
            Cigar::Del(n) => {
                // A deletion can straddle the run's edges; each part goes where its reference positions do.
                let (s, t) = (rp, rp + n);
                let head = (a + 1).min(t).saturating_sub(s);
                let inside = overlap(s, t, a + 1, e);
                let tail = t.saturating_sub(e.max(s));
                push_op(&mut pre_ops, Cigar::Del(head));
                run_ref += inside;
                if tail > 0 {
                    past = true;
                    push_op(&mut post_ops, Cigar::Del(tail));
                }
                rp = t;
            }
            Cigar::RefSkip(n) => {
                if overlap(rp, rp + n, a + 1, e) > 0 {
                    return None;
                }
                if rp >= e {
                    past = true;
                    push_op(&mut post_ops, Cigar::RefSkip(n));
                } else {
                    push_op(&mut pre_ops, Cigar::RefSkip(n));
                }
                rp += n;
            }
            Cigar::SoftClip(n) => {
                let n = n as usize;
                let bases = seq.get(qp..qp + n)?;
                let quals = qual.get(qp..qp + n)?;
                if past {
                    push_op(&mut post_ops, Cigar::SoftClip(n as u32));
                    post_seq.extend_from_slice(bases);
                    post_qual.extend_from_slice(quals);
                } else {
                    push_op(&mut pre_ops, Cigar::SoftClip(n as u32));
                    pre_seq.extend_from_slice(bases);
                    pre_qual.extend_from_slice(quals);
                }
                qp += n;
            }
            Cigar::HardClip(n) => {
                let ops = if past { &mut post_ops } else { &mut pre_ops };
                push_op(ops, Cigar::HardClip(n));
            }
            Cigar::Pad(n) => {
                let ops = if past { &mut post_ops } else { &mut pre_ops };
                push_op(ops, Cigar::Pad(n));
            }
        }
    }

    if run_ref != run_len as u32 {
        return None;
    }

    // Canonical body: the run becomes exactly `run_len + target_k` copies of `b`.
    let target_len = (run_len + target_k) as usize;
    let body_seq: Vec<u8> = vec![b; target_len];
    // QUALITY is not rewritten wholesale: the read's own qualities over the run are kept.
    let mut body_qual: Vec<u8> = Vec::with_capacity(target_len);
    if let Some(add) = target_len.checked_sub(run_qual.len()).filter(|n| *n > 0) {
        let right = run_qual.first().copied().or_else(|| post_qual.first().copied());
        let f = match (pre_qual.last().copied(), right) {
            (Some(l), Some(r)) => l.min(r),
            (Some(l), None) => l,
            (None, Some(r)) => r,
            (None, None) => 0,
        };
        body_qual.extend(std::iter::repeat(f).take(add));
    }
    let drop = run_qual.len().saturating_sub(target_len);
    body_qual.extend_from_slice(&run_qual[drop..]);

    let mut ops = pre_ops;
    if target_k > 0 {
        push_op(&mut ops, Cigar::Ins(target_k as u32));
        push_op(&mut ops, Cigar::Match(run_len as u32));
    } else {
        push_op(&mut ops, Cigar::Del((-target_k) as u32));
        push_op(&mut ops, Cigar::Match((run_len + target_k) as u32));
    }
    for op in post_ops {
        push_op(&mut ops, op);
    }

    let mut out_seq = pre_seq;
    out_seq.extend_from_slice(&body_seq);
    out_seq.extend_from_slice(&post_seq);
    let mut out_qual = pre_qual;
    out_qual.extend_from_slice(&body_qual);
    out_qual.extend_from_slice(&post_qual);

    Some((CigarString(ops), out_seq, out_qual))
}

/// Remove tags invalidated by a structural change.
pub fn strip_stale_tags(rec: &mut rust_htslib::bam::Record) {
    for tag in [
        b"MM".as_slice(),
        b"ML".as_slice(),
        b"MN".as_slice(),
        b"NM".as_slice(),
        b"MD".as_slice(),
        b"cs".as_slice(),
    ] {
        // Absent tags are expected; removal failure just means "wasn't there".
        let _ = rec.remove_aux(tag);
    }
}

#[cfg(test)]
mod rewrite_tests {
    use super::*;

    // Reference used by every test below (index 0-based): 0123456789...
    const REF: &[u8] = b"ACGTAAAAACGTACGT";

    fn cig(ops: Vec<Cigar>) -> CigarString {
        CigarString(ops)
    }

    /// Assert the design's core invariant plus SEQ/QUAL length agreement.
    fn assert_invariants(before: &CigarString, r: &RewriteResult) {
        assert_eq!(
            ref_consumed(before),
            ref_consumed(&r.cigar),
            "reference consumption changed: {before:?} -> {:?}",
            r.cigar
        );
        assert_eq!(r.seq.len(), r.qual.len(), "SEQ and QUAL lengths disagree");
        let query: u32 = r
            .cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) | Cigar::Equal(n)
                | Cigar::Diff(n) => *n,
                _ => 0,
            })
            .sum();
        assert_eq!(query as usize, r.seq.len(), "CIGAR query length != SEQ length");
    }

    #[test]
    fn no_edits_is_a_faithful_copy() {
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let r = rewrite_read(0, &c, &seq, &qual, &ReadEdits::default(), REF);
        assert_eq!(r.seq, seq);
        assert_eq!(r.qual, qual);
        assert_eq!(r.cigar, c);
        assert!(!r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn substitution_only_preserves_length_and_cigar() {
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits { subs: vec![(2, b'T')], indels: vec![] };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACTTAAAA".to_vec());
        assert_eq!(r.qual, qual);
        assert!(!r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn reverting_an_insertion_removes_query_bases() {
        // Read = ACGT + "GG" inserted after ref index 3 + AAAA.
        let c = cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)]);
        let seq = b"ACGTGGAAAA".to_vec();
        let qual = vec![30u8; 10];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ins(b"GG".to_vec()), to: Allele::Ref }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(8)]));
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn reverting_a_deletion_splices_reference_bases_back() {
        // Read skips ref indices 4..6 (AA), i.e. 4M 2D 3M over ACGT|AA|AAA...
        let c = cig(vec![Cigar::Match(4), Cigar::Del(2), Cigar::Match(3)]);
        let seq = b"ACGTAAA".to_vec();
        let qual = vec![20u8, 21, 22, 23, 24, 25, 26];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Del(2), to: Allele::Ref }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(9)]));
        // Spliced bases take min(left flank, right flank) = min(23, 24) = 23.
        assert_eq!(r.qual, vec![20, 21, 22, 23, 23, 23, 24, 25, 26]);
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn gaining_an_insertion_splits_a_match_run() {
        // Plain 8M read gains a consensus "A" insertion after ref index 3.
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![20u8, 21, 22, 23, 24, 25, 26, 27];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4)]));
        // Inserted base takes min(23, 24) = 23.
        assert_eq!(r.qual, vec![20, 21, 22, 23, 23, 24, 25, 26, 27]);
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn gaining_a_deletion_converts_match_to_del() {
        // Plain 8M read gains a 2bp deletion after ref index 3.
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Del(2) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Del(2), Cigar::Match(2)]));
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn replacing_an_insertion_swaps_the_inserted_bases() {
        let c = cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)]);
        let seq = b"ACGTGGAAAA".to_vec();
        let qual = vec![30u8; 10];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit {
                ref_pos: 3,
                from: Allele::Ins(b"GG".to_vec()),
                to: Allele::Ins(b"A".to_vec()),
            }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4)]));
        assert_invariants(&c, &r);
    }

    #[test]
    fn shrinking_a_deletion_emits_del_then_match() {
        // 3bp deletion shrinks to 1bp: Del(1) first, then the two reference bases that stop being deleted come back as matches.
        let c = cig(vec![Cigar::Match(4), Cigar::Del(3), Cigar::Match(2)]);
        let seq = b"ACGTAC".to_vec();
        let qual = vec![30u8; 6];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Del(3), to: Allele::Del(1) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        // ref[5..7] = "AA" return to the read.
        assert_eq!(r.seq, b"ACGTAAAC".to_vec());
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(4)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn one_read_carrying_sub_insertion_and_deletion_edits() {
        // 12M read: substitution at ref 1, gained insertion after ref 3, gained deletion after ref 8.
        let c = cig(vec![Cigar::Match(12)]);
        let seq = b"ACGTAAAAACGT".to_vec();
        let qual = vec![30u8; 12];
        let edits = ReadEdits {
            subs: vec![(1, b'T')],
            indels: vec![
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) },
                IndelEdit { ref_pos: 8, from: Allele::Ref, to: Allele::Del(2) },
            ],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ATGTAAAAAAT".to_vec());
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(5), Cigar::Del(2), Cigar::Match(1)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn soft_clips_are_copied_verbatim() {
        let c = cig(vec![Cigar::SoftClip(3), Cigar::Match(8), Cigar::SoftClip(2)]);
        let seq = b"TTTACGTAAAAGG".to_vec();
        let qual = vec![30u8; 13];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(&r.seq[..3], b"TTT");
        assert_eq!(&r.seq[r.seq.len() - 2..], b"GG");
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::SoftClip(3), Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4), Cigar::SoftClip(2)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn read_starting_at_nonzero_pos_uses_absolute_reference_coordinates() {
        // Read aligned at POS 4 (the AAAAA run).
        let c = cig(vec![Cigar::Match(6)]);
        let seq = b"AAAAAC".to_vec();
        let qual = vec![30u8; 6];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 5, from: Allele::Ref, to: Allele::Del(1) }],
        };
        let r = rewrite_read(4, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"AAAAC".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(2), Cigar::Del(1), Cigar::Match(3)]));
        assert_invariants(&c, &r);
    }

    #[test]
    fn two_edits_sharing_a_ref_pos_but_differing_in_kind_are_both_applied() {
        // A single read can legitimately carry two IndelEdits with the SAME ref_pos: one keyed at the read's own CIGAR anchor.
        let c = cig(vec![Cigar::Match(16), Cigar::Ins(1)]);
        let mut seq = REF.to_vec();
        seq.extend_from_slice(b"A");
        let qual = vec![30u8; 17];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![
                // Gained: at ref_pos 15 (from: Ref), splice in a 2bp consensus insertion "AT" right after the last matched base.
                IndelEdit {
                    ref_pos: 15,
                    from: Allele::Ref,
                    to: Allele::Ins(b"AT".to_vec()),
                },
                // Reverted: this read's OWN carried insertion (from: Ins), also anchored at ref_pos 15, is dropped entirely.
                IndelEdit {
                    ref_pos: 15,
                    from: Allele::Ins(b"A".to_vec()),
                    to: Allele::Ref,
                },
            ],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        // The gained "AT" is spliced in (proving the Ref-kind lookup fired).
        let mut expected_seq = REF.to_vec();
        expected_seq.extend_from_slice(b"AT");
        assert_eq!(r.seq, expected_seq);
        assert_eq!(r.cigar, cig(vec![Cigar::Match(16), Cigar::Ins(2)]));
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn duplicate_ref_pos_and_kind_keeps_the_first_edit_and_never_panics() {
        // two edits with the IDENTICAL pair -- both `from: Ref` at position 3 -- are a LEGITIMATE input.
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) },
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Del(1) },
            ],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        // The FIRST edit (Ref -> Ins) is the one applied...
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4)]));
        // ...and its key appears exactly once in `applied` -- the SECOND edit (Ref -> Del) never gets a chance to run at all.
        assert_eq!(r.applied.len(), 1);
        assert!(r.applied.contains(&(3, EditKind::Ref)));
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn gained_deletion_overrunning_its_match_op_is_refused_not_truncated() {
        // a gained deletion whose span runs past the remaining length of its Match op must be REFUSED outright.
        let c = cig(vec![Cigar::Match(4)]);
        let seq = b"ACGT".to_vec();
        let qual = vec![30u8; 4];
        let edits = ReadEdits {
            subs: vec![],
            // Anchored at ref 0, the first base of the only Match op: after that anchor base is emitted.
            indels: vec![IndelEdit { ref_pos: 0, from: Allele::Ref, to: Allele::Del(4) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, seq, "the read must be unchanged when the deletion cannot fit");
        assert_eq!(r.cigar, c, "CIGAR must be unchanged");
        assert!(r.applied.is_empty(), "a refused edit must not appear in `applied`");
        assert!(!r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn substitution_only_with_equal_and_diff_ops_preserves_cigar_verbatim() {
        // the substitutions-only fast path must return the ORIGINAL CIGAR byte-for-byte, `=`/`X` ops included.
        let c = cig(vec![Cigar::Equal(2), Cigar::Diff(1), Cigar::Equal(5)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits { subs: vec![(2, b'T')], indels: vec![] };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.cigar, c, "CIGAR must be byte-identical to the input, `=`/`X` ops included");
        assert_eq!(r.seq, b"ACTTAAAA".to_vec());
        assert_eq!(r.qual, qual);
        assert!(!r.structure_changed);
        assert!(r.applied.is_empty());
        assert_invariants(&c, &r);
    }

    #[test]
    fn one_edit_shared_by_two_ins_ops_separated_only_by_a_deletion_is_applied_once() {
        // CIGAR 10M 1I 2D 1I 10M -- both insertions share the SAME lookup key.
        let c = cig(vec![
            Cigar::Match(10),
            Cigar::Ins(1),
            Cigar::Del(2),
            Cigar::Ins(1),
            Cigar::Match(10),
        ]);
        let mut seq = vec![b'A'; 10];
        seq.push(b'T'); // first inserted base
        seq.push(b'C'); // second inserted base
        seq.extend(vec![b'A'; 10]);
        let qual = vec![30u8; 22];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 9, from: Allele::Ins(b"T".to_vec()), to: Allele::Ref }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);

        // The FIRST Ins op reaches the key first and is reverted.
        let mut expected_seq = vec![b'A'; 10];
        expected_seq.push(b'C');
        expected_seq.extend(vec![b'A'; 10]);
        assert_eq!(
            r.seq, expected_seq,
            "exactly one Ins op must be dropped -- the other must survive untouched"
        );
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(10), Cigar::Del(2), Cigar::Ins(1), Cigar::Match(10)])
        );
        assert_eq!(
            r.applied.len(),
            1,
            "the shared key must be recorded as applied exactly once, not twice"
        );
        assert!(r.applied.contains(&(9, EditKind::Ins)));
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn one_edit_shared_by_two_del_ops_separated_only_by_an_insertion_is_applied_once() {
        // Mirror of the Ins-side test above: CIGAR 10M 2D 1I 2D 10M.
        let c = cig(vec![
            Cigar::Match(10),
            Cigar::Del(2),
            Cigar::Ins(1),
            Cigar::Del(2),
            Cigar::Match(10),
        ]);
        let mut seq = vec![b'A'; 10];
        seq.push(b'T'); // the lone inserted base, between the two deletions
        seq.extend(vec![b'A'; 10]);
        let qual = vec![30u8; 21];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 9, from: Allele::Del(2), to: Allele::Ref }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);

        // The FIRST deletion is reverted: REF[10..12] ("GT") splices back in as matches.
        let mut expected_seq = vec![b'A'; 10];
        expected_seq.extend_from_slice(b"GT");
        expected_seq.push(b'T');
        expected_seq.extend(vec![b'A'; 10]);
        assert_eq!(
            r.seq, expected_seq,
            "exactly one deletion must be reverted -- the other must survive as a real deletion"
        );
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(12), Cigar::Ins(1), Cigar::Del(2), Cigar::Match(10)])
        );
        assert_eq!(
            r.applied.len(),
            1,
            "the shared key must be recorded as applied exactly once, not twice"
        );
        assert!(r.applied.contains(&(9, EditKind::Del)));
        assert!(r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn offsetting_insertion_and_deletion_leave_length_unchanged_but_set_structure_changed() {
        // 12M read gains a 1bp insertion after ref index 3 and a 1bp deletion after ref index 8.
        let c = cig(vec![Cigar::Match(12)]);
        let seq = b"ACGTAAAAACGT".to_vec();
        let qual = vec![30u8; 12];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) },
                IndelEdit { ref_pos: 8, from: Allele::Ref, to: Allele::Del(1) },
            ],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(
            r.seq.len(),
            seq.len(),
            "offsetting insertion+deletion should leave SEQ length unchanged"
        );
        assert_eq!(r.seq, b"ACGTAAAAAAGT".to_vec());
        assert!(
            r.structure_changed,
            "an indel edit was applied even though output length is unchanged"
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn leading_insertion_is_unaffected_by_an_edit_anchored_before_any_aligned_base() {
        // `last_ref` starts as `None`: a read whose first ref-relevant op is an insertion has no aligned base preceding it.
        let c = cig(vec![Cigar::Ins(2), Cigar::Match(10)]);
        let seq = b"GGACGTAAAAAC".to_vec();
        let qual = vec![30u8; 12];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit {
                ref_pos: 0,
                from: Allele::Ins(b"GG".to_vec()),
                to: Allele::Ref,
            }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, seq, "leading insertion must be copied verbatim");
        assert_eq!(r.cigar, c);
        assert!(!r.structure_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn no_edits_preserves_equal_and_diff_ops_verbatim() {
        // 's fast path: with zero edits, the original CIGAR must come back byte-identical, not folded into `M`.
        let c = cig(vec![Cigar::Equal(4), Cigar::Diff(1), Cigar::Equal(3)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let r = rewrite_read(0, &c, &seq, &qual, &ReadEdits::default(), REF);
        assert_eq!(r.cigar, c, "CIGAR must be byte-identical when there are no edits");
        assert_eq!(r.seq, seq);
        assert_eq!(r.qual, qual);
        assert!(!r.structure_changed);
        assert_invariants(&c, &r);
    }

    // ---- Plan 2026-09-26 T8d: canonical homopolymer-run rewrite ----
    const HPR_A: u32 = 3;
    const HPR_E: u32 = 9;
    const HPR_B: u8 = b'A';

    /// `rewrite_hp_run` returns the same three pieces `rewrite_read` does.
    fn hp_result(t: (CigarString, Vec<u8>, Vec<u8>)) -> RewriteResult {
        RewriteResult {
            cigar: t.0,
            seq: t.1,
            qual: t.2,
            structure_changed: true,
            applied: HashSet::new(),
        }
    }

    #[test]
    fn hp_rewrite_to_the_same_length_is_a_faithful_copy() {
        let c = cig(vec![Cigar::Match(16)]);
        let seq = REF.to_vec();
        let qual = vec![30u8; 16];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 0).expect("must rewrite"),
        );
        assert_eq!(r.cigar, c);
        assert_eq!(r.seq, seq);
        assert_eq!(r.qual, qual);
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_shrinks_a_reference_length_run_to_a_deletion() {
        // Ref -> Del(1): the canonical form is left-aligned.
        let c = cig(vec![Cigar::Match(16)]);
        let seq = REF.to_vec();
        let qual = vec![30u8; 16];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, -1).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(11)]));
        assert_eq!(r.seq.len(), 15, "one base must leave the read");
        assert_eq!(r.seq, b"ACGTAAAACGTACGT".to_vec());
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_grows_a_run_with_a_left_aligned_insertion() {
        let c = cig(vec![Cigar::Match(16)]);
        let seq = REF.to_vec();
        let qual = vec![30u8; 16];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 1).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(12)]));
        assert_eq!(r.seq, b"ACGTAAAAAACGTACGT".to_vec(), "the run must grow to 6 A's");
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_absorbs_the_reads_own_insertion() {
        // 9M 1I 6M: the extra A is anchored on REF[8], the last base of the run.
        let c = cig(vec![Cigar::Match(9), Cigar::Ins(1), Cigar::Match(7)]);
        let mut seq = REF[..9].to_vec();
        seq.push(b'A');
        seq.extend_from_slice(&REF[9..]);
        let qual = vec![30u8; 17];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 0).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(16)]));
        assert_eq!(r.seq, REF.to_vec());
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_turns_an_insertion_into_a_deletion() {
        // The transition `rewrite_read` cannot make at all (cross-kind).
        let c = cig(vec![Cigar::Match(9), Cigar::Ins(1), Cigar::Match(7)]);
        let mut seq = REF[..9].to_vec();
        seq.push(b'A');
        seq.extend_from_slice(&REF[9..]);
        let qual = vec![30u8; 17];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, -1).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(11)]));
        assert_eq!(r.seq, b"ACGTAAAACGTACGT".to_vec());
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_grows_the_reads_own_deletion() {
        // The other transition `rewrite_read` refuses: Del(1) -> Del(3) would consume reference from the following op there.
        let c = cig(vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(11)]);
        let mut seq = REF[..4].to_vec();
        seq.extend_from_slice(&REF[5..]);
        let qual = vec![30u8; 15];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, -3).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Del(3), Cigar::Match(9)]));
        assert_eq!(r.seq.len(), 13);
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_can_delete_the_whole_run() {
        let c = cig(vec![Cigar::Match(16)]);
        let seq = REF.to_vec();
        let qual = vec![30u8; 16];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, -5).expect("must rewrite"),
        );
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(4), Cigar::Del(5), Cigar::Match(7)]),
            "a zero-length match after the deletion must not be emitted"
        );
        assert_eq!(r.seq, b"ACGTCGTACGT".to_vec());
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_grows_a_run_the_read_deleted_entirely() {
        // The run's base is a property of the REFERENCE.
        let c = cig(vec![Cigar::Match(4), Cigar::Del(5), Cigar::Match(7)]);
        let mut seq = REF[..4].to_vec();
        seq.extend_from_slice(&REF[9..]);
        let qual = vec![30u8; 11];

        // ...back to the reference length: the whole run is synthesized as `b`.
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 0).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(16)]));
        assert_eq!(r.seq, REF.to_vec(), "the restored run must be five A's");
        assert_invariants(&c, &r);

        // ...and past it, to a run two bases longer than the reference.
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 2).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(12)]));
        assert_eq!(r.seq, b"ACGTAAAAAAACGTACGT".to_vec());
        assert_invariants(&c, &r);

        assert!(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, -5).is_some(),
            "the read already carries this length, so no base has to be synthesized"
        );
    }

    #[test]
    fn hp_rewrite_synthesizes_the_run_base_not_whatever_the_read_carries() {
        // Every base the rewrite emits over the run must be `b`.
        let c = cig(vec![Cigar::Match(16)]);
        let mut seq = REF.to_vec();
        seq[4] = b'C';
        let qual = vec![30u8; 16];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 1).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(12)]));
        assert_eq!(r.seq[4], b'A', "the synthesized base is the run base, not the read's C");
        assert!(
            !r.seq[4..10].contains(&b'C'),
            "the rewritten run must be a pure run of the tract base"
        );
        assert_eq!(r.seq, b"ACGTAAAAAACGTACGT".to_vec());
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_insertion_is_never_a_slice_of_the_read() {
        // Building the body as "synthesized `b` x add.
        let c = cig(vec![Cigar::Match(9), Cigar::Ins(1), Cigar::Match(7)]);
        let mut seq = REF[..9].to_vec();
        seq[4] = b'C';
        seq.push(b'A');
        seq.extend_from_slice(&REF[9..]);
        let qual = vec![30u8; 17];
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 2).expect("must rewrite"),
        );
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(12)]));
        assert_eq!(
            &r.seq[4..6],
            b"AA",
            "the insertion must be 2 copies of the run base, not a slice of the read"
        );
        assert_eq!(r.seq, b"ACGTAAAAAAACGTACGT".to_vec());
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_refuses_a_target_shorter_than_nothing() {
        let c = cig(vec![Cigar::Match(16)]);
        let seq = REF.to_vec();
        let qual = vec![30u8; 16];
        assert!(rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, -6).is_none());
    }

    #[test]
    fn hp_rewrite_refuses_a_read_that_does_not_cover_the_whole_run() {
        // Ends at reference 6, inside the run: there is no whole run to replace.
        let c = cig(vec![Cigar::Match(6)]);
        let seq = REF[..6].to_vec();
        let qual = vec![30u8; 6];
        assert!(rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, -1).is_none());
        // Starts at reference 5, inside the run: same, from the other side.
        let c = cig(vec![Cigar::Match(11)]);
        let seq = REF[5..].to_vec();
        let qual = vec![30u8; 11];
        assert!(rewrite_hp_run(5, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, -1).is_none());
    }

    #[test]
    fn hp_rewrite_keeps_soft_clips_and_absolute_coordinates() {
        // 2S 12M from POS 2: the clip is copied verbatim and the run still resolves by ABSOLUTE reference coordinate.
        let c = cig(vec![Cigar::SoftClip(2), Cigar::Match(12)]);
        let mut seq = b"NN".to_vec();
        seq.extend_from_slice(&REF[2..14]);
        let qual = vec![30u8; 14];
        let r = hp_result(
            rewrite_hp_run(2, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 1).expect("must rewrite"),
        );
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::SoftClip(2), Cigar::Match(2), Cigar::Ins(1), Cigar::Match(10)])
        );
        assert_eq!(r.seq, b"NNGTAAAAAACGTAC".to_vec());
        assert_invariants(&c, &r);
    }

    #[test]
    fn hp_rewrite_gives_synthesized_bases_the_minimum_flanking_quality() {
        // Same rule as `fill_qual`: never more confident than the neighbours.
        let c = cig(vec![Cigar::Match(16)]);
        let seq = REF.to_vec();
        let mut qual = vec![40u8; 16];
        qual[3] = 20;
        qual[4] = 30;
        let r = hp_result(
            rewrite_hp_run(0, &c, &seq, &qual, HPR_A, HPR_E, HPR_B, 1).expect("must rewrite"),
        );
        assert_eq!(r.qual[4], 20, "min(20, 30), not 40 and not an invented value");
        assert_eq!(r.qual[3], 20, "the anchor's own quality must survive");
        assert_eq!(r.qual[5], 30, "the original run bases keep their qualities");
        assert_invariants(&c, &r);
    }

    #[test]
    fn stale_tags_are_stripped_and_others_survive() {
        use rust_htslib::bam::record::{Aux, Record};

        let mut rec = Record::new();
        rec.set(
            b"r1",
            Some(&cig(vec![Cigar::Match(4)])),
            b"ACGT",
            &[30u8; 4],
        );
        rec.push_aux(b"MM", Aux::String("C+m,0;")).unwrap();
        rec.push_aux(b"ML", Aux::U8(200)).unwrap();
        rec.push_aux(b"MN", Aux::I32(4)).unwrap();
        rec.push_aux(b"NM", Aux::I32(2)).unwrap();
        rec.push_aux(b"MD", Aux::String("4")).unwrap();
        rec.push_aux(b"cs", Aux::String(":4")).unwrap();
        rec.push_aux(b"XY", Aux::I32(7)).unwrap();

        strip_stale_tags(&mut rec);

        assert!(rec.aux(b"MM").is_err(), "MM must be stripped");
        assert!(rec.aux(b"ML").is_err(), "ML must be stripped");
        assert!(rec.aux(b"MN").is_err(), "MN must be stripped");
        assert!(rec.aux(b"NM").is_err(), "NM must be stripped");
        assert!(rec.aux(b"MD").is_err(), "MD must be stripped");
        assert!(rec.aux(b"cs").is_err(), "cs must be stripped");
        assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)), "unrelated tags must survive");
    }

    #[test]
    fn stripping_absent_tags_is_a_no_op() {
        use rust_htslib::bam::record::{Aux, Record};

        let mut rec = Record::new();
        rec.set(b"r1", Some(&cig(vec![Cigar::Match(4)])), b"ACGT", &[30u8; 4]);
        rec.push_aux(b"XY", Aux::I32(7)).unwrap();
        strip_stale_tags(&mut rec); // must not panic on missing tags
        assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)));
    }
}
