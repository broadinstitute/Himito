//! Small-indel site model and length-changing read rewrite for the ONT denoiser.
//! Every rewrite preserves the reference bases a read consumes, so POS, reference end and
//! coordinate-sort order are invariant and the output never needs re-sorting or re-indexing.

use std::collections::{HashMap, HashSet};

use rust_htslib::bam::record::{Cigar, CigarString};

use crate::denoise::{HOMOPLASMIC_VAF, MIN_STRAND, STRAND_BIAS_P};

/// An indel allele at a left-normalized site.
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
}

/// Minimum plain match-context flank required on each side of an indel edit.
pub const INDEL_FLANK: usize = 5;
/// Exclusive upper bound on modelled indel length: 1..=4 bp are correctable, longer ones are never touched.
pub const MAX_INDEL_LEN: u32 = 5;
// Per-read indel error rate in repeat context L: ERR0 * ERR_SCALE^(L-1), capped at ERR_CAP.
const ERR0: f64 = 0.01;
const ERR_SCALE: f64 = 1.5;
const ERR_CAP: f64 = 0.4;
/// Minimum VAF for an indel alt to be kept.
const VAF_ABS_MIN: f64 = 0.05;
/// Significance of the binomial excess-over-error test, Bonferroni-corrected over sites.
const EXCESS_ALPHA: f64 = 1e-3;
/// Per-base decay of misobservation mass with net-length distance in `assign_allele`.
const DELTA: f64 = 0.3;

/// User-settable indel denoising options.
#[derive(Debug, Clone)]
pub struct IndelOpts {
    pub enabled: bool,
    /// OFF by default (CLI `--indel-hp-length-model`) until an eval shows it helps.
    pub hp_length_model: bool,
    /// An allele carried by at least this fraction of reads at a site is never corrected away.
    pub protect_vaf: f64,
}

impl Default for IndelOpts {
    fn default() -> Self {
        IndelOpts { enabled: false, hp_length_model: false, protect_vaf: 0.2 }
    }
}

/// Shift an indel maximally left, so one event gets one site however the aligner anchored it.
/// `anchor` is the last reference base before the event; an insertion rotates as it shifts.
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

/// Reference repeat context `L` at a normalized site: the number of tandem copies of the event's
/// primitive unit around `norm_pos + 1`, at least 1.
pub fn repeat_context(refseq: &[u8], norm_pos: u32, allele: &Allele) -> u32 {
    let start = norm_pos as usize + 1;
    let unit: &[u8] = match allele {
        Allele::Ref => return 1,
        Allele::Ins(s) => s,
        Allele::Del(n) => match refseq.get(start..start + *n as usize) {
            Some(w) => w,
            None => return 1,
        },
    };
    if unit.is_empty() {
        return 1;
    }
    let unit = primitive_period(unit);
    let u = unit.len();

    // Copies running rightward from the event position, then leftward.
    let mut copies = 0u32;
    let mut p = start;
    while refseq.get(p..p + u).map_or(false, |w| w == unit) {
        copies += 1;
        p += u;
    }
    let mut p = start;
    while p >= u && refseq.get(p - u..p).map_or(false, |w| w == unit) {
        copies += 1;
        p -= u;
    }
    copies.max(1)
}

/// Exclusive reference end of the homopolymer run starting at `anchor + 1`, or `anchor + 1` itself
/// when that base does not repeat (a lone base is not a run).
pub fn hp_tract_end(refseq: &[u8], anchor: u32) -> u32 {
    let start = anchor as usize + 1;
    let Some(&base) = refseq.get(start) else {
        return start as u32;
    };
    let mut p = start + 1;
    while refseq.get(p) == Some(&base) {
        p += 1;
    }
    if p - start >= 2 {
        p as u32
    } else {
        start as u32
    }
}

/// If `anchor` is the left-normalized anchor of a homopolymer run (a run of >= 2 starting at
/// `anchor + 1`, preceded by a different base), returns the run's exclusive end and its base.
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

/// The allele that changes a homopolymer run of base `b` by `k` bases.
pub fn hp_length_allele(k: i32, b: u8) -> Allele {
    match k {
        0 => Allele::Ref,
        k if k > 0 => Allele::Ins(vec![b; k as usize]),
        k => Allele::Del((-k) as u32),
    }
}

/// Per-read indel error probability in repeat context `l`.
pub fn error_rate(l: u32) -> f64 {
    let exp = (l as i32 - 1).max(0);
    (ERR0 * ERR_SCALE.powi(exp)).min(ERR_CAP)
}

/// One-sided binomial P(X >= n) over `depth` reads at per-read error rate `eps`.
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
    /// Hashed qnames of reads that already voted at this site: one vote per read.
    voted: HashSet<u64>,
    /// At an HP site this, and NOT `depth - alt_total`, is the reference support.
    hp_ref: u32,
}

impl IndelSite {
    /// Count one ALT observation; returns `false`, counting nothing, if this read already voted here.
    pub fn add_from_read(&mut self, qname: &[u8], allele: Allele, reverse: bool) -> bool {
        if !self.mark_voted(qname) {
            return false;
        }
        let e = self.obs.entry(allele).or_insert((0, 0));
        if reverse {
            e.1 += 1;
        } else {
            e.0 += 1;
        }
        true
    }

    /// Count one HP-site reference vote, under the same one-vote-per-read cap as `add_from_read`.
    pub fn add_hp_ref(&mut self, qname: &[u8]) {
        if self.mark_voted(qname) {
            self.hp_ref += 1;
        }
    }

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

/// The decision for one site. All allele lists are sorted, for deterministic output.
#[derive(Debug, Clone)]
pub struct IndelSiteModel {
    /// Alleles passing candidacy (REF always) with frequencies renormalized over this set.
    pub kept: Vec<(Allele, f64)>,
    pub eps: f64,
    /// The repeat context length that produced `eps`.
    pub context_len: u32,
    pub strand_biased: usize,
    /// Raw observed VAF of every allele at the site, REF included. Never renormalized, so candidacy
    /// rejection does not affect the `protect_vaf` check in `assign_allele`.
    pub obs_vaf: Vec<(Allele, f64)>,
    /// Alleles rejected by a strand gate (per-strand minimum or strand bias).
    pub strand_rejected: Vec<Allele>,
    /// The site's homopolymer run `(exclusive end, base)` when the HP length model is on; `None`
    /// when the model is off or the site is not a homopolymer.
    pub hp_run: Option<(u32, u8)>,
}

/// Decide which alleles survive at a site and compute their frequencies.
pub fn fit_indel_site(
    site: &IndelSite,
    refseq: &[u8],
    norm_pos: u32,
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

    // One eps per site, from the deepest context of any observed allele.
    let context_len = site
        .obs
        .keys()
        .map(|a| repeat_context(refseq, norm_pos, a))
        .max()
        .unwrap_or(1);
    let eps = error_rate(context_len);

    // Reads that cast no vote at this site account for any slack.
    debug_assert!(
        site.hp_ref + site.alt_total() <= site.depth,
        "indel site accounting: hp_ref ({}) + alt_total ({}) exceeds depth ({})",
        site.hp_ref,
        site.alt_total(),
        site.depth
    );

    let hp_run = if o.hp_length_model { hp_site(refseq, norm_pos) } else { None };
    // At an HP site REF support is voted explicitly; elsewhere REF gets whatever depth the events did not claim.
    let ref_support = if hp_run.is_some() {
        site.hp_ref
    } else {
        depth.saturating_sub(site.alt_total())
    };
    // Raw counts over kept alleles.
    let mut raw = vec![(Allele::Ref, ref_support as f64)];
    let mut obs_vaf = vec![(Allele::Ref, (ref_support as f64 / depth as f64).min(1.0))];
    let mut strand_rejected: Vec<Allele> = Vec::new();
    let mut strand_biased = 0usize;

    for (allele, &(fwd, rev)) in site.obs.iter() {
        let n = fwd + rev;
        let vaf = (n as f64 / depth as f64).min(1.0);
        obs_vaf.push((allele.clone(), vaf));
        if n == 0 {
            continue;
        }
        if fwd < MIN_STRAND || rev < MIN_STRAND {
            strand_rejected.push(allele.clone());
            continue;
        }
        let eps_l = error_rate(repeat_context(refseq, norm_pos, allele));
        if vaf < VAF_ABS_MIN.max(1.5 * eps_l)
            || indel_excess_pvalue(n, depth, eps_l) > EXCESS_ALPHA / n_sites.max(1) as f64
        {
            continue;
        }
        // Near-homoplasmic exemption: a single-strand artifact cannot reach that VAF.
        if vaf < HOMOPLASMIC_VAF
            && crate::call::is_strand_skewed(fwd as usize, rev as usize, expected_fwd_frac, STRAND_BIAS_P)
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

    // `site.obs` is a HashMap with a randomly-seeded hasher.
    kept.sort_by(|(a, _), (b, _)| a.cmp(b));
    obs_vaf.sort_by(|(a, _), (b, _)| a.cmp(b));
    strand_rejected.sort();

    IndelSiteModel { kept, eps, context_len, strand_biased, obs_vaf, strand_rejected, hp_run }
}

/// Outcome of MAP assignment for one read at one site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Assignment {
    /// The read already carries the MAP allele, its allele is protected, or the site is degenerate.
    Keep,
    /// Reassign the read to this allele.
    Move(Allele),
    /// The MAP allele is only reachable through a transition the rewrite walk does not support.
    Refused,
}

/// Maximum-a-posteriori allele for a read observing `observed`.
pub fn assign_allele(observed: &Allele, m: &IndelSiteModel, o: &IndelOpts) -> Assignment {
    let eps = m.eps.clamp(0.0, 0.5);
    let decay = |a: &Allele, b: &Allele| DELTA.powi((a.net_len() - b.net_len()).abs());

    // Given truth k, error mass spreads over the kept alleles plus `observed`, decaying in net-length distance.
    let mut best: Option<(&Allele, f64)> = None;
    for (k, f) in &m.kept {
        let like = if k == observed {
            1.0 - eps
        } else {
            let mut z: f64 = m.kept.iter().filter(|(a, _)| a != k).map(|(a, _)| decay(a, k)).sum();
            if !m.kept.iter().any(|(a, _)| a == observed) {
                z += decay(observed, k);
            }
            if z <= 0.0 {
                continue;
            }
            eps * decay(observed, k) / z
        };
        let val = f * like;
        if best.map_or(true, |(_, b)| val > b) {
            best = Some((k, val));
        }
    }

    let Some((target, _)) = best else {
        return Assignment::Keep;
    };
    if target == observed {
        return Assignment::Keep;
    }
    // Outside an HP site the rewrite walk can neither cross Ins <-> Del nor grow a deletion.
    if m.hp_run.is_none() {
        let reachable = match (observed, target) {
            (Allele::Ref, _) | (_, Allele::Ref) | (Allele::Ins(_), Allele::Ins(_)) => true,
            (Allele::Del(n), Allele::Del(t)) => t <= n,
            _ => false,
        };
        if !reachable {
            return Assignment::Refused;
        }
    }
    // A substantial allele is never corrected away, even when it lost candidacy -- unless a strand gate rejected it.
    let observed_vaf = m
        .obs_vaf
        .iter()
        .find(|(a, _)| a == observed)
        .map_or(0.0, |(_, v)| *v);
    if observed_vaf >= o.protect_vaf && !m.strand_rejected.contains(observed) {
        return Assignment::Keep;
    }
    Assignment::Move(target.clone())
}

/// One indel carried by a read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadEvent {
    /// Left-normalized anchor: the site key.
    pub norm_pos: u32,
    /// The read's own CIGAR anchor: the last aligned base BEFORE the event.
    pub anchor: u32,
    /// `None` for an indel of `MAX_INDEL_LEN` or more (never corrected); `norm_pos` is then `anchor`.
    pub allele: Option<Allele>,
}

/// Aligned (`M`/`=`/`X`) bases left in the CIGAR op covering reference position `target`, counting
/// `target` itself; 0 when `target` is deleted, skipped, or outside the read.
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
pub fn read_events(pos: i64, cigar: &CigarString, seq: &[u8], refseq: &[u8]) -> Vec<ReadEvent> {
    let mut out = Vec::new();
    let mut rp = pos.max(0) as u32; // next reference position
    let mut qp = 0usize; // next query position
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
            Cigar::Ins(n) | Cigar::Del(n) => {
                let is_ins = matches!(op, Cigar::Ins(_));
                // An indel before any aligned base has no anchor and is skipped.
                if let Some(anchor) = last_aligned {
                    out.push(if n < MAX_INDEL_LEN {
                        let allele = if is_ins {
                            Allele::Ins(seq[qp..qp + n as usize].to_vec())
                        } else {
                            Allele::Del(n)
                        };
                        let (norm_pos, allele) = normalize_left(refseq, anchor, &allele);
                        ReadEvent { norm_pos, anchor, allele: Some(allele) }
                    } else {
                        ReadEvent { norm_pos: anchor, anchor, allele: None }
                    });
                }
                if is_ins {
                    qp += n as usize;
                } else {
                    rp += n;
                }
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

/// Net length change a read applies to the homopolymer run with anchor `a` and exclusive end `e`.
/// `None` unless the read aligns the anchor and reaches `INDEL_FLANK` bases past the run, and also
/// when it inserts a base other than the run's or an `N` op crosses the run.
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
                // An insertion anchored on the anchor base or on any base of the run belongs to the run.
                if last_aligned.is_some_and(|l| l >= a && l < e) {
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
    (rp >= e + INDEL_FLANK as u32).then_some(k)
}

/// One read's reassignment at one normalized indel site.
#[derive(Debug, Clone)]
pub struct IndelEdit {
    /// The aligned base `rewrite_read` hangs this edit off: the site's normalized anchor for a gained
    /// indel (`from == Ref`), the read's own CIGAR anchor otherwise.
    pub ref_pos: u32,
    pub from: Allele,
    pub to: Allele,
}

/// The variety of an `IndelEdit`'s `from` allele, which selects where the rewrite walk looks it up.
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

/// All edits for one read, keyed by reference coordinate; `subs` are `(ref_pos, new base)`.
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
    /// `true` iff at least one indel edit was applied (`!applied.is_empty()`).
    pub structure_changed: bool,
    /// `(ref_pos, edit_kind(from))` keys of the indel edits actually applied.
    pub applied: HashSet<(u32, EditKind)>,
}

/// Quality for synthesized bases: never more confident than the flanking bases.
fn flank_qual(left: Option<u8>, right: Option<u8>) -> u8 {
    match (left, right) {
        (Some(l), Some(r)) => l.min(r),
        (Some(q), None) | (None, Some(q)) => q,
        (None, None) => 0,
    }
}

/// Append `op`, merging it into a preceding op of the same kind (M/I/D/S); zero-length ops are dropped.
fn push_op(ops: &mut Vec<Cigar>, op: Cigar) {
    if op.len() == 0 {
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

/// Rewrite one read's SEQ/QUAL/CIGAR against `edits`. `refseq` must be the full contig the read is
/// aligned to, indexed by absolute reference coordinate (reverted deletions splice bases from it).
pub fn rewrite_read(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    qual: &[u8],
    edits: &ReadEdits,
    refseq: &[u8],
) -> RewriteResult {
    let subs: HashMap<u32, u8> = edits.subs.iter().copied().collect();

    if edits.indels.is_empty() {
        // Substitutions only: the CIGAR, `=`/`X` ops included, is returned verbatim.
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

    // On a duplicate (ref_pos, kind) key the first edit wins; the rest stay unapplied.
    let mut indels: HashMap<(u32, EditKind), &IndelEdit> = HashMap::new();
    for e in &edits.indels {
        indels.entry((e.ref_pos, edit_kind(&e.from))).or_insert(e);
    }

    let mut out_seq: Vec<u8> = Vec::with_capacity(seq.len() + 8);
    let mut out_qual: Vec<u8> = Vec::with_capacity(qual.len() + 8);
    let mut ops: Vec<Cigar> = Vec::new();
    let mut applied: HashSet<(u32, EditKind)> = HashSet::new();

    let mut rp: u32 = pos.max(0) as u32; // reference position of the next ref-consuming base
    let mut qp: usize = 0; // query position of the next original base
    // Reference position of the last aligned base emitted: the pileup anchor convention.
    let mut last_ref: Option<u32> = None;

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

                    // A gained indel anchored on this aligned base.
                    match indels.get(&(cur, EditKind::Ref)).map(|e| &e.to) {
                        Some(Allele::Ins(s)) => {
                            let f = flank_qual(out_qual.last().copied(), qual.get(qp + k as usize).copied());
                            out_seq.extend_from_slice(s);
                            out_qual.extend(std::iter::repeat(f).take(s.len()));
                            push_op(&mut ops, Cigar::Ins(s.len() as u32));
                            applied.insert((cur, EditKind::Ref));
                        }
                        // The following m reference bases become a deletion, if they lie in this op.
                        Some(&Allele::Del(m)) if k + m <= n => {
                            push_op(&mut ops, Cigar::Del(m));
                            k += m; // skip both the ref bases and their query partners
                            applied.insert((cur, EditKind::Ref));
                        }
                        _ => {}
                    }
                }
                qp += n as usize;
                rp += n;
            }

            Cigar::Ins(n) => {
                // Before the first aligned base there is no anchor to look an edit up by.
                let edit = last_ref.and_then(|lr| indels.remove(&(lr, EditKind::Ins)).map(|e| (lr, &e.to)));
                match edit {
                    // Reverted: drop the inserted bases entirely.
                    Some((lr, Allele::Ref)) => {
                        applied.insert((lr, EditKind::Ins));
                    }
                    // Replaced: emit the site's consensus inserted bases instead.
                    Some((lr, Allele::Ins(s))) => {
                        let f = flank_qual(out_qual.last().copied(), qual.get(qp + n as usize).copied());
                        out_seq.extend_from_slice(s);
                        out_qual.extend(std::iter::repeat(f).take(s.len()));
                        push_op(&mut ops, Cigar::Ins(s.len() as u32));
                        applied.insert((lr, EditKind::Ins));
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
                let edit = last_ref.and_then(|lr| indels.remove(&(lr, EditKind::Del)).map(|e| (lr, &e.to)));
                // Emit reference bases [from, from + count) at flank quality.
                let splice = |out_seq: &mut Vec<u8>, out_qual: &mut Vec<u8>, from: u32, count: u32| {
                    let f = flank_qual(out_qual.last().copied(), qual.get(qp).copied());
                    for idx in from as usize..(from + count) as usize {
                        debug_assert!(idx < refseq.len(), "refseq must be the full contig (index {idx})");
                        out_seq.push(refseq.get(idx).copied().unwrap_or(b'N'));
                        out_qual.push(f);
                    }
                };
                match edit {
                    // Reverted: the deleted reference bases come back as matches.
                    Some((lr, Allele::Ref)) => {
                        splice(&mut out_seq, &mut out_qual, rp, n);
                        push_op(&mut ops, Cigar::Match(n));
                        applied.insert((lr, EditKind::Del));
                    }
                    // Shrunk: left-aligned, so the shorter deletion comes first and the freed reference bases follow as matches.
                    Some((lr, &Allele::Del(m))) if m <= n => {
                        push_op(&mut ops, Cigar::Del(m));
                        splice(&mut out_seq, &mut out_qual, rp + m, n - m);
                        push_op(&mut ops, Cigar::Match(n - m));
                        applied.insert((lr, EditKind::Del));
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
        structure_changed: !applied.is_empty(),
        applied,
    }
}

/// CIGAR ops together with the query bases and qualities they emit.
#[derive(Default)]
struct Piece {
    ops: Vec<Cigar>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl Piece {
    fn push(&mut self, op: Cigar, seq: &[u8], qual: &[u8]) {
        push_op(&mut self.ops, op);
        self.seq.extend_from_slice(seq);
        self.qual.extend_from_slice(qual);
    }
}

/// Rewrite the homopolymer run of base `b` with anchor `a` and exclusive end `e` to `run_len +
/// target_k` bases in canonical left-aligned form, absorbing any indels the read carries in it.
/// Unlike `rewrite_read`, every target length is reachable. `None` when the read does not align the
/// anchor and account for the whole run, or an `N` op crosses the run.
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
    if target_k < -run_len || match_run_len(pos, cigar, a) == 0 {
        return None;
    }

    // Everything up to the anchor and everything from `e` on is copied verbatim.
    let (mut pre, mut post) = (Piece::default(), Piece::default());
    // Of the run itself only the read's qualities are kept.
    let mut run_qual: Vec<u8> = Vec::new();
    // Reference bases of the run the read actually accounts for.
    let mut run_ref: u32 = 0;

    let mut rp = pos.max(0) as u32;
    let mut qp = 0usize;
    let mut last_aligned: Option<u32> = None;
    // Set once the walk reaches the run's end: later ops belong after it.
    let mut past = false;

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for i in 0..n {
                    let cur = rp + i;
                    let qi = qp + i as usize;
                    let (base, q) = (*seq.get(qi)?, *qual.get(qi)?);
                    if cur <= a {
                        pre.push(Cigar::Match(1), &[base], &[q]);
                    } else if cur < e {
                        run_ref += 1;
                        run_qual.push(q);
                    } else {
                        past = true;
                        post.push(Cigar::Match(1), &[base], &[q]);
                    }
                    last_aligned = Some(cur);
                }
                rp += n;
                qp += n as usize;
            }
            Cigar::Ins(n) => {
                let bases = seq.get(qp..qp + n as usize)?;
                let quals = qual.get(qp..qp + n as usize)?;
                // Anchored on the anchor base or on any base of the run: part of the run's length, so absorbed.
                if last_aligned.is_some_and(|l| l >= a && l < e) {
                    run_qual.extend_from_slice(quals);
                } else {
                    (if past { &mut post } else { &mut pre }).push(Cigar::Ins(n), bases, quals);
                }
                qp += n as usize;
            }
            Cigar::Del(n) => {
                // A deletion can straddle the run's edges; each part goes where its reference positions do.
                let (s, t) = (rp, rp + n);
                pre.push(Cigar::Del((a + 1).min(t).saturating_sub(s)), &[], &[]);
                run_ref += overlap(s, t, a + 1, e);
                let tail = t.saturating_sub(e.max(s));
                if tail > 0 {
                    past = true;
                    post.push(Cigar::Del(tail), &[], &[]);
                }
                rp = t;
            }
            Cigar::RefSkip(n) => {
                if overlap(rp, rp + n, a + 1, e) > 0 {
                    return None;
                }
                if rp >= e {
                    past = true;
                    post.push(Cigar::RefSkip(n), &[], &[]);
                } else {
                    pre.push(Cigar::RefSkip(n), &[], &[]);
                }
                rp += n;
            }
            Cigar::SoftClip(n) => {
                let bases = seq.get(qp..qp + n as usize)?;
                let quals = qual.get(qp..qp + n as usize)?;
                (if past { &mut post } else { &mut pre }).push(Cigar::SoftClip(n), bases, quals);
                qp += n as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {
                (if past { &mut post } else { &mut pre }).push(*op, &[], &[]);
            }
        }
    }

    if run_ref != run_len as u32 {
        return None;
    }

    // Canonical body: exactly `run_len + target_k` copies of `b`. The read's own run qualities are
    // kept, dropped from the left when shrinking or padded on the left at flank quality when growing.
    let target_len = (run_len + target_k) as usize;
    let right = run_qual.first().copied().or_else(|| post.qual.first().copied());
    let pad = flank_qual(pre.qual.last().copied(), right);
    let mut body_qual = vec![pad; target_len.saturating_sub(run_qual.len())];
    body_qual.extend_from_slice(&run_qual[run_qual.len().saturating_sub(target_len)..]);

    let mut ops = pre.ops;
    if target_k > 0 {
        push_op(&mut ops, Cigar::Ins(target_k as u32));
        push_op(&mut ops, Cigar::Match(run_len as u32));
    } else {
        push_op(&mut ops, Cigar::Del((-target_k) as u32));
        push_op(&mut ops, Cigar::Match((run_len + target_k) as u32));
    }
    for op in post.ops {
        push_op(&mut ops, op);
    }
    let out_seq = [pre.seq, vec![b; target_len], post.seq].concat();
    let out_qual = [pre.qual, body_qual, post.qual].concat();
    Some((CigarString(ops), out_seq, out_qual))
}

/// Remove tags a structural change invalidates: base modifications (MM/ML/MN) and alignment
/// differences (NM/MD/cs).
pub fn strip_stale_tags(rec: &mut rust_htslib::bam::Record) {
    for tag in [b"MM", b"ML", b"MN", b"NM", b"MD", b"cs"] {
        // Absent tags are expected; removal failure just means "wasn't there".
        let _ = rec.remove_aux(tag);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Cigar::{
        Del as D, Diff as X, Equal as E, Ins as I, Match as M, RefSkip as N, SoftClip as S,
    };

    fn cig(ops: &[Cigar]) -> CigarString {
        CigarString(ops.to_vec())
    }

    fn ins(s: &[u8]) -> Allele {
        Allele::Ins(s.to_vec())
    }

    // C + 7 A + G: a 1bp A event normalizes to anchor 0, the run's HP site (exclusive end 8).
    const HP_REF: &[u8] = b"CAAAAAAAG";
    // A 12-mer A-run, deep enough to cap `error_rate`.
    const DEEP_HP_REF: &[u8] = b"CAAAAAAAAAAAAG";
    const UNIQ_REF: &[u8] = b"ACGTACGTACGTACGT";
    const C5_REF: &[u8] = b"TACCCCCTCTA";
    // The A-run is 8..=15: HP anchor 7, exclusive end 16.
    const REF29: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC";

    #[test]
    fn normalize_left_moves_to_the_leftmost_equivalent_anchor() {
        // Every placement of a 1bp event inside the A-run collapses to anchor 0.
        for anchor in 0..=6 {
            assert_eq!(normalize_left(HP_REF, anchor, &Allele::Del(1)), (0, Allele::Del(1)));
        }
        for anchor in 0..=7 {
            assert_eq!(normalize_left(HP_REF, anchor, &ins(b"A")), (0, ins(b"A")));
        }
        let cases: &[(&[u8], u32, Allele, (u32, Allele))] = &[
            // An insertion rotates as it shifts: "TA" after 1 is "AT" after 0.
            (b"CATATATG", 1, ins(b"TA"), (0, ins(b"AT"))),
            (b"ACGTACGT", 1, Allele::Del(1), (1, Allele::Del(1))), // unique sequence: no move
            (b"AAAAA", 3, Allele::Del(1), (0, Allele::Del(1))),    // stops at the contig start
            (HP_REF, 4, Allele::Ref, (4, Allele::Ref)),
        ];
        for (r, anchor, a, want) in cases {
            assert_eq!(&normalize_left(r, *anchor, a), want, "{a:?} at {anchor}");
        }
    }

    #[test]
    fn repeat_context_counts_copies_of_the_primitive_unit() {
        let cases: &[(&[u8], u32, Allele, u32)] = &[
            (HP_REF, 0, Allele::Del(1), 7),
            (HP_REF, 0, ins(b"A"), 7),
            (b"CATATATATG", 0, ins(b"AT"), 4),
            (b"CATATATATG", 0, Allele::Del(2), 4),
            // Multi-copy events reduce to their period.
            (b"CATATATATG", 0, ins(b"ATAT"), 4),
            (b"CATATATATG", 0, Allele::Del(4), 4),
            (b"CAAAAAAAAG", 0, Allele::Del(3), 8),
            (b"CAAAAAAAAG", 0, ins(b"AA"), 8),
            (C5_REF, 1, Allele::Del(1), 5),
            (b"ACGTACGT", 1, Allele::Del(1), 1),
            (b"ACGTACGT", 1, Allele::Ref, 1),
        ];
        for (r, pos, a, want) in cases {
            assert_eq!(repeat_context(r, *pos, a), *want, "{a:?} at {pos}");
        }
    }

    #[test]
    fn hp_tract_end_and_hp_site() {
        assert_eq!(hp_tract_end(REF29, 7), 16);
        assert_eq!(hp_tract_end(REF29, 11), 16); // from inside the run
        assert_eq!(hp_tract_end(REF29, 20), 21); // REF29[21] is a lone G
        assert_eq!(hp_tract_end(REF29, 28), 29); // past the contig end
        // Only the run's leftmost anchor is an HP site; 15 is followed by a lone C.
        assert_eq!(hp_site(REF29, 7), Some((16, b'A')));
        for anchor in (8..=15).chain([20, 28]) {
            assert_eq!(hp_site(REF29, anchor), None, "anchor {anchor}");
        }
        assert_eq!(hp_length_allele(0, b'A'), Allele::Ref);
        assert_eq!(hp_length_allele(2, b'A'), ins(b"AA"));
        assert_eq!(hp_length_allele(-3, b'A'), Allele::Del(3));
    }

    #[test]
    fn error_rate_grows_with_context_and_is_capped() {
        assert!((error_rate(1) - 0.01).abs() < 1e-12);
        assert!((error_rate(2) - 0.015).abs() < 1e-12);
        for l in 1..40 {
            assert!(error_rate(l + 1) >= error_rate(l), "not monotonic at {l}");
            assert!(error_rate(l) <= ERR_CAP, "cap exceeded at {l}");
        }
        assert_eq!(error_rate(30), ERR_CAP);
    }

    /// A site of `depth` reads with `events` (allele, fwd, rev); the remaining reads are REF, strand-balanced.
    fn site(depth: u32, events: &[(Allele, u32, u32)]) -> IndelSite {
        let mut s = IndelSite { depth, ..Default::default() };
        let (mut f, mut r) = (0, 0);
        for (a, fwd, rev) in events {
            s.obs.insert(a.clone(), (*fwd, *rev));
            f += fwd;
            r += rev;
        }
        let rest = depth - f - r;
        s.col_fwd = f + rest / 2;
        s.col_rev = r + rest - rest / 2;
        s
    }

    fn fit(r: &[u8], pos: u32, s: &IndelSite, n_sites: usize, hp: bool) -> IndelSiteModel {
        let o = IndelOpts { hp_length_model: hp, ..Default::default() };
        fit_indel_site(s, r, pos, n_sites, &o)
    }

    #[test]
    fn candidacy_gates() {
        let del = Allele::Del(1);
        // (case, refseq, norm_pos, depth, fwd, rev, column strand totals, n_sites, kept, strand-biased count)
        let cases: &[(&str, &[u8], u32, u32, u32, u32, Option<(u32, u32)>, usize, bool, usize)] = &[
            ("balanced 20% in unique sequence", UNIQ_REF, 1, 100, 10, 10, None, 1, true, 0),
            ("same evidence in a 7-mer run", HP_REF, 0, 100, 10, 10, None, 1, false, 0),
            ("one reverse read fails min_strand", UNIQ_REF, 1, 100, 19, 1, None, 1, false, 0),
            ("18/2 against a balanced column", UNIQ_REF, 1, 100, 18, 2, None, 1, false, 1),
            ("80% skewed: near-homoplasmic exemption", UNIQ_REF, 1, 100, 75, 5, None, 1, true, 0),
            ("deep, significant skew but only a 0.11 effect", UNIQ_REF, 1, 5240, 328, 196, Some((2693, 2547)), 1, true, 0),
            ("deep single-strand artifact, 0.50 effect", UNIQ_REF, 1, 4640, 449, 15, Some((2171, 2469)), 1, false, 1),
            ("10% at depth 431 in a C5 run", C5_REF, 1, 431, 22, 21, None, 1, true, 0),
            ("6% at depth 60 in a C5 run", C5_REF, 1, 60, 2, 2, None, 1, false, 0),
            ("10% at depth 60 in a C5 run: excess not significant", C5_REF, 1, 60, 3, 3, None, 1, false, 0),
            ("6% at depth 60 in unique sequence", UNIQ_REF, 1, 60, 2, 2, None, 1, false, 0),
            ("43/431 does not survive 16569 tests", C5_REF, 1, 431, 22, 21, None, 16569, false, 0),
            ("419/5334 survives 16569 tests", C5_REF, 1, 5334, 210, 209, None, 16569, true, 0),
        ];
        for &(case, r, pos, depth, fwd, rev, col, n_sites, kept, biased) in cases {
            let mut s = site(depth, &[(del.clone(), fwd, rev)]);
            if let Some((cf, cr)) = col {
                (s.col_fwd, s.col_rev) = (cf, cr);
            }
            let m = fit(r, pos, &s, n_sites, false);
            assert_eq!(m.kept.iter().any(|(a, _)| *a == del), kept, "{case}");
            assert!(m.kept.iter().any(|(a, _)| *a == Allele::Ref), "{case}: REF is always kept");
            assert_eq!(m.strand_biased, biased, "{case}");
        }
        // The 0.11-effect case is saved by the effect-size gate, not by an insignificant p-value.
        assert!(crate::call::strand_bias_pvalue(328, 196, 2693.0 / 5240.0) < STRAND_BIAS_P);
    }

    #[test]
    fn kept_and_obs_vaf_are_sorted() {
        let m = fit(UNIQ_REF, 1, &site(100, &[(ins(b"A"), 10, 10), (Allele::Del(1), 10, 10)]), 1, false);
        let want = [Allele::Ref, ins(b"A"), Allele::Del(1)];
        assert!(m.kept.iter().map(|(a, _)| a).eq(want.iter()), "{:?}", m.kept);
        assert!(m.obs_vaf.iter().map(|(a, _)| a).eq(want.iter()), "{:?}", m.obs_vaf);
    }

    #[test]
    fn assign_allele_cases() {
        use Allele::{Del, Ref};
        use Assignment::{Keep, Move, Refused};
        // (case, refseq, norm_pos, site events at depth 100, HP length model, observed, observed kept, expected)
        let cases: Vec<(&str, &[u8], u32, Vec<(Allele, u32, u32)>, bool, Allele, bool, Assignment)> = vec![
            ("40% INS heteroplasmy: INS read", HP_REF, 0, vec![(ins(b"A"), 20, 20)], false, ins(b"A"), true, Keep),
            ("40% INS heteroplasmy: REF read", HP_REF, 0, vec![(ins(b"A"), 20, 20)], false, Ref, true, Keep),
            ("90% INS: a REF read dropped the indel", HP_REF, 0, vec![(ins(b"A"), 45, 45)], false, Ref, true, Move(ins(b"A"))),
            ("90% INS: INS read", HP_REF, 0, vec![(ins(b"A"), 45, 45)], false, ins(b"A"), true, Keep),
            ("REF-only site: REF read", UNIQ_REF, 1, vec![(Del(1), 1, 1)], false, Ref, true, Keep),
            ("REF-only site: stray Del", UNIQ_REF, 1, vec![(Del(1), 1, 1)], false, Del(1), false, Move(Ref)),
            ("REF-only site: unobserved stray Ins", UNIQ_REF, 1, vec![(Del(1), 1, 1)], false, ins(b"T"), false, Move(Ref)),
            ("failed Del(2) snaps to kept Del(1)", UNIQ_REF, 1, vec![(Del(1), 40, 40), (Del(2), 1, 1)], false, Del(2), false, Move(Del(1))),
            ("Del -> Ins crosses kinds", UNIQ_REF, 1, vec![(ins(b"T"), 45, 45), (Del(1), 1, 1)], false, Del(1), false, Refused),
            ("Del(1) -> Del(3) grows a deletion", UNIQ_REF, 1, vec![(Del(3), 45, 45), (Del(1), 1, 1)], false, Del(1), false, Refused),
            // At an HP site every kept length is reachable.
            ("Ins -> Del without the HP model", HP_REF, 0, vec![(Del(1), 45, 45), (ins(b"A"), 1, 1)], false, ins(b"A"), false, Refused),
            ("Ins -> Del with the HP model", HP_REF, 0, vec![(Del(1), 45, 45), (ins(b"A"), 1, 1)], true, ins(b"A"), false, Move(Del(1))),
            ("Del(1) -> Del(3) with the HP model", HP_REF, 0, vec![(Del(3), 45, 45), (Del(1), 1, 1)], true, Del(1), false, Move(Del(3))),
            // protect_vaf shields floor-rejected alleles, but not rare or strand-rejected ones.
            ("5% floor-rejected Del", DEEP_HP_REF, 0, vec![(Del(1), 3, 2)], false, Del(1), false, Move(Ref)),
            ("30% floor-rejected Del", DEEP_HP_REF, 0, vec![(Del(1), 15, 15)], false, Del(1), false, Keep),
            ("30% strand-rejected Ins", UNIQ_REF, 1, vec![(ins(b"T"), 30, 0)], false, ins(b"T"), false, Move(Ref)),
            ("HP model: 30% floor-rejected Ins", DEEP_HP_REF, 0, vec![(Del(1), 30, 30), (ins(b"A"), 15, 15)], true, ins(b"A"), false, Keep),
            ("HP model: 30% strand-rejected Ins", DEEP_HP_REF, 0, vec![(Del(1), 30, 30), (ins(b"A"), 30, 0)], true, ins(b"A"), false, Move(Ref)),
        ];
        for (case, r, pos, events, hp, observed, observed_kept, want) in cases {
            let m = fit(r, pos, &site(100, &events), 1, hp);
            assert_eq!(m.hp_run.is_some(), hp, "{case}");
            assert_eq!(m.kept.iter().any(|(a, _)| *a == observed), observed_kept, "{case}: {:?}", m.kept);
            assert_eq!(assign_allele(&observed, &m, &IndelOpts::default()), want, "{case}");
        }
    }

    #[test]
    fn match_run_len_counts_aligned_bases_left_in_the_op() {
        let cases: &[(&[Cigar], i64, u32, u32)] = &[
            (&[M(10)], 0, 0, 10),
            (&[M(10)], 0, 3, 7),
            (&[M(10)], 0, 9, 1),
            (&[M(10)], 0, 10, 0), // past the end
            (&[M(10)], 5, 3, 0),  // before POS
            (&[M(5), D(3), M(5)], 0, 5, 0),
            (&[M(5), D(3), M(5)], 0, 7, 0),
            (&[M(5), D(3), M(5)], 0, 8, 5),
            (&[M(5), N(3), M(5)], 0, 6, 0),
            // An insertion splits ops without consuming reference.
            (&[M(5), I(2), M(5)], 0, 4, 1),
            (&[M(5), I(2), M(5)], 0, 5, 5),
        ];
        for (ops, pos, target, want) in cases {
            assert_eq!(match_run_len(*pos, &cig(ops), *target), *want, "{ops:?} at {pos}, target {target}");
        }
    }

    #[test]
    fn read_events_reports_left_normalized_events() {
        let ev = |norm_pos, anchor, allele| ReadEvent { norm_pos, anchor, allele };
        let (a150, a41) = (vec![b'A'; 150], vec![b'A'; 41]);
        let cases: Vec<(&[u8], i64, Vec<Cigar>, &[u8], Vec<ReadEvent>)> = vec![
            (REF29, 0, vec![M(29)], REF29, vec![]),
            // An extra / missing A: the aligner's anchor and the normalized site differ.
            (REF29, 0, vec![M(16), I(1), M(13)], b"ACGTACGTAAAAAAAAACGTACGTACGTAC", vec![ev(7, 15, Some(ins(b"A")))]),
            (REF29, 0, vec![M(15), D(1), M(13)], b"ACGTACGTAAAAAAACGTACGTACGTAC", vec![ev(7, 14, Some(Allele::Del(1)))]),
            (REF29, 0, vec![M(19), I(1), M(10)], b"ACGTACGTAAAAAAAACGTGACGTACGTAC", vec![ev(18, 18, Some(ins(b"G")))]),
            (REF29, 0, vec![M(19), I(6), M(10)], b"ACGTACGTAAAAAAAACGTGGGGGGACGTACGTAC", vec![ev(18, 18, None)]),
            (REF29, 0, vec![M(16), I(1), M(4), D(1), M(8)], b"ACGTACGTAAAAAAAAACGTCGTACGTAC",
             vec![ev(7, 15, Some(ins(b"A"))), ev(19, 19, Some(Allele::Del(1)))]),
            (REF29, 0, vec![S(3), M(19), I(1), M(10)], b"TTTACGTACGTAAAAAAAACGTGACGTACGTAC", vec![ev(18, 18, Some(ins(b"G")))]),
            (REF29, 16, vec![M(3), I(1), M(10)], b"CGTGACGTACGTAC", vec![ev(18, 18, Some(ins(b"G")))]),
            // An insertion right after a deletion anchors at the last ALIGNED base, not the last deleted one.
            (&a150, 100, vec![M(20), D(2), I(1), M(20)], &a41,
             vec![ev(0, 119, Some(Allele::Del(2))), ev(0, 119, Some(ins(b"A")))]),
        ];
        for (r, pos, ops, seq, want) in cases {
            assert_eq!(read_events(pos, &cig(&ops), seq, r), want, "{ops:?}");
        }
    }

    #[test]
    fn hp_allele_is_the_net_run_length_change() {
        let r = REF29;
        let cases: Vec<(&str, Vec<Cigar>, Vec<u8>, Option<i32>)> = vec![
            ("pure match", vec![M(29)], r.to_vec(), Some(0)),
            ("extra A anchored on the run's last base", vec![M(16), I(1), M(13)], [&r[..16], b"A", &r[16..]].concat(), Some(1)),
            ("2bp deleted inside the run", vec![M(10), D(2), M(17)], [&r[..10], &r[12..]].concat(), Some(-2)),
            ("offsetting ins and del", vec![M(10), I(1), M(2), D(1), M(16)], [&r[..10], b"A", &r[10..12], &r[13..]].concat(), Some(0)),
            ("inserted G breaks the run", vec![M(16), I(1), M(13)], [&r[..16], b"G", &r[16..]].concat(), None),
            ("read ends inside the run", vec![M(13)], r[..13].to_vec(), None),
            // e + INDEL_FLANK = 21.
            ("one base short of the right flank", vec![M(20)], r[..20].to_vec(), None),
            ("exactly the right flank", vec![M(21)], r[..21].to_vec(), Some(0)),
            ("anchor deleted", vec![M(5), D(5), M(19)], [&r[..5], &r[10..]].concat(), None),
            ("N op crosses the run", vec![M(10), N(2), M(17)], [&r[..10], &r[12..]].concat(), None),
        ];
        for (case, ops, seq, want) in cases {
            assert_eq!(hp_allele(0, &cig(&ops), &seq, r, 7, 16), want, "{case}");
        }
    }

    // Reference for the rewrite tests; its A-run 4..9 has HP anchor 3.
    const REF: &[u8] = b"ACGTAAAAACGTACGT";

    fn ref_consumed(c: &CigarString) -> u32 {
        c.iter().filter(|op| matches!(op, M(_) | E(_) | X(_) | D(_) | N(_))).map(|op| op.len()).sum()
    }

    /// Reference consumption is preserved, and SEQ, QUAL and the CIGAR's query length agree.
    fn assert_invariants(before: &CigarString, cigar: &CigarString, seq: &[u8], qual: &[u8]) {
        assert_eq!(ref_consumed(before), ref_consumed(cigar), "{before} -> {cigar}");
        assert_eq!(seq.len(), qual.len(), "SEQ and QUAL lengths disagree");
        let query: u32 = cigar.iter().filter(|op| matches!(op, M(_) | E(_) | X(_) | I(_) | S(_))).map(|op| op.len()).sum();
        assert_eq!(query as usize, seq.len(), "CIGAR query length != SEQ length");
    }

    fn ed(ref_pos: u32, from: Allele, to: Allele) -> IndelEdit {
        IndelEdit { ref_pos, from, to }
    }

    /// Run `rewrite_read` with qualities 20, 21, ... and check SEQ, CIGAR and the applied edit keys.
    fn check_rewrite(
        pos: i64,
        ops: &[Cigar],
        seq: &[u8],
        subs: &[(u32, u8)],
        indels: Vec<IndelEdit>,
        want_seq: &[u8],
        want_ops: &[Cigar],
        want_applied: &[(u32, EditKind)],
    ) -> RewriteResult {
        let c = cig(ops);
        let qual: Vec<u8> = (0..seq.len()).map(|i| 20 + i as u8).collect();
        let edits = ReadEdits { subs: subs.to_vec(), indels };
        let r = rewrite_read(pos, &c, seq, &qual, &edits, REF);
        assert_eq!(r.seq, want_seq, "{c}");
        assert_eq!(r.cigar, cig(want_ops), "{c}");
        assert_eq!(r.applied, want_applied.iter().copied().collect::<HashSet<_>>(), "{c}");
        assert_eq!(r.structure_changed, !want_applied.is_empty(), "{c}");
        if want_applied.is_empty() {
            assert_eq!(r.qual, qual, "{c}");
        }
        assert_invariants(&c, &r.cigar, &r.seq, &r.qual);
        r
    }

    #[test]
    fn rewrite_read_cases() {
        use Allele::{Del, Ref};
        use EditKind as K;
        // No indel edits: the CIGAR comes back verbatim, `=`/`X` ops included.
        check_rewrite(0, &[M(8)], b"ACGTAAAA", &[], vec![], b"ACGTAAAA", &[M(8)], &[]);
        check_rewrite(0, &[E(4), X(1), E(3)], b"ACGTAAAA", &[], vec![], b"ACGTAAAA", &[E(4), X(1), E(3)], &[]);
        check_rewrite(0, &[M(8)], b"ACGTAAAA", &[(2, b'T')], vec![], b"ACTTAAAA", &[M(8)], &[]);
        check_rewrite(0, &[E(2), X(1), E(5)], b"ACGTAAAA", &[(2, b'T')], vec![], b"ACTTAAAA", &[E(2), X(1), E(5)], &[]);

        // Revert, gain and replace an insertion; inserted bases take min(left, right) flank quality.
        check_rewrite(0, &[M(4), I(2), M(4)], b"ACGTGGAAAA", &[], vec![ed(3, ins(b"GG"), Ref)], b"ACGTAAAA", &[M(8)], &[(3, K::Ins)]);
        let r = check_rewrite(0, &[M(8)], b"ACGTAAAA", &[], vec![ed(3, Ref, ins(b"A"))], b"ACGTAAAAA", &[M(4), I(1), M(4)], &[(3, K::Ref)]);
        assert_eq!(r.qual, [20, 21, 22, 23, 23, 24, 25, 26, 27]);
        check_rewrite(0, &[M(4), I(2), M(4)], b"ACGTGGAAAA", &[], vec![ed(3, ins(b"GG"), ins(b"A"))], b"ACGTAAAAA", &[M(4), I(1), M(4)], &[(3, K::Ins)]);

        // Revert, gain and shrink a deletion; spliced reference bases take min(left, right) flank quality.
        let r = check_rewrite(0, &[M(4), D(2), M(3)], b"ACGTAAA", &[], vec![ed(3, Del(2), Ref)], b"ACGTAAAAA", &[M(9)], &[(3, K::Del)]);
        assert_eq!(r.qual, [20, 21, 22, 23, 23, 23, 24, 25, 26]);
        check_rewrite(0, &[M(8)], b"ACGTAAAA", &[], vec![ed(3, Ref, Del(2))], b"ACGTAA", &[M(4), D(2), M(2)], &[(3, K::Ref)]);
        // A shrunk deletion stays left-aligned: Del(1) first, then the freed bases as matches.
        check_rewrite(0, &[M(4), D(3), M(2)], b"ACGTAC", &[], vec![ed(3, Del(3), Del(1))], b"ACGTAAAC", &[M(4), D(1), M(4)], &[(3, K::Del)]);

        // A substitution, a gained insertion and a gained deletion in one read.
        check_rewrite(0, &[M(12)], b"ACGTAAAAACGT", &[(1, b'T')], vec![ed(3, Ref, ins(b"A")), ed(8, Ref, Del(2))],
                      b"ATGTAAAAAAT", &[M(4), I(1), M(5), D(2), M(1)], &[(3, K::Ref), (8, K::Ref)]);
        // Soft clips are copied verbatim; POS shifts edits to absolute reference coordinates.
        check_rewrite(0, &[S(3), M(8), S(2)], b"TTTACGTAAAAGG", &[], vec![ed(3, Ref, ins(b"A"))],
                      b"TTTACGTAAAAAGG", &[S(3), M(4), I(1), M(4), S(2)], &[(3, K::Ref)]);
        check_rewrite(4, &[M(6)], b"AAAAAC", &[], vec![ed(5, Ref, Del(1))], b"AAAAC", &[M(2), D(1), M(3)], &[(5, K::Ref)]);

        // Same ref_pos, different kinds: a gained insertion and the read's own reverted one both apply.
        check_rewrite(0, &[M(16), I(1)], &[REF, b"A"].concat(), &[], vec![ed(15, Ref, ins(b"AT")), ed(15, ins(b"A"), Ref)],
                      &[REF, b"AT"].concat(), &[M(16), I(2)], &[(15, K::Ref), (15, K::Ins)]);
        // Same (ref_pos, kind): the first edit wins.
        check_rewrite(0, &[M(8)], b"ACGTAAAA", &[], vec![ed(3, Ref, ins(b"A")), ed(3, Ref, Del(1))],
                      b"ACGTAAAAA", &[M(4), I(1), M(4)], &[(3, K::Ref)]);
        // A gained deletion running past its match op is refused, not truncated.
        check_rewrite(0, &[M(4)], b"ACGT", &[], vec![ed(0, Ref, Del(4))], b"ACGT", &[M(4)], &[]);
        // Two Ins (Del) ops sharing an anchor across a Del (Ins): one edit applies to the first only.
        let a10: &[u8] = b"AAAAAAAAAA";
        check_rewrite(0, &[M(10), I(1), D(2), I(1), M(10)], &[a10, b"TC", a10].concat(), &[], vec![ed(9, ins(b"T"), Ref)],
                      &[a10, b"C", a10].concat(), &[M(10), D(2), I(1), M(10)], &[(9, K::Ins)]);
        check_rewrite(0, &[M(10), D(2), I(1), D(2), M(10)], &[a10, b"T", a10].concat(), &[], vec![ed(9, Del(2), Ref)],
                      &[a10, b"GT", b"T", a10].concat(), &[M(12), I(1), D(2), M(10)], &[(9, K::Del)]);
        // A leading insertion has no anchor, so no edit can reach it.
        check_rewrite(0, &[I(2), M(10)], b"GGACGTAAAAAC", &[], vec![ed(0, ins(b"GG"), Ref)], b"GGACGTAAAAAC", &[I(2), M(10)], &[]);
    }

    /// `rewrite_hp_run` on REF's A-run (anchor 3, end 9) with uniform quality.
    fn check_hp(pos: i64, ops: &[Cigar], seq: &[u8], target_k: i32, want: Option<(&[Cigar], &[u8])>) {
        let c = cig(ops);
        let qual = vec![30u8; seq.len()];
        let got = rewrite_hp_run(pos, &c, seq, &qual, 3, 9, b'A', target_k);
        match (got, want) {
            (Some((gc, gs, gq)), Some((wc, ws))) => {
                assert_eq!((&gc, gs.as_slice()), (&cig(wc), ws), "{c} to k={target_k}");
                assert_invariants(&c, &gc, &gs, &gq);
            }
            (None, None) => {}
            (got, _) => panic!("{c} to k={target_k}: got {:?}", got.map(|g| g.0.to_string())),
        }
    }

    #[test]
    fn rewrite_hp_run_cases() {
        let with_ins = [&REF[..9], b"A", &REF[9..]].concat(); // own extra A, anchored on the run's last base
        let no_run = [&REF[..4], &REF[9..]].concat();
        let mut c_in_run = REF.to_vec();
        c_in_run[4] = b'C';
        let c_in_run_with_ins = [&REF[..4], b"C", &REF[5..9], b"A", &REF[9..]].concat();

        check_hp(0, &[M(16)], REF, 0, Some((&[M(16)], REF)));
        check_hp(0, &[M(16)], REF, -1, Some((&[M(4), D(1), M(11)], b"ACGTAAAACGTACGT")));
        check_hp(0, &[M(16)], REF, 1, Some((&[M(4), I(1), M(12)], b"ACGTAAAAAACGTACGT")));
        // Deleting the whole run emits no zero-length match.
        check_hp(0, &[M(16)], REF, -5, Some((&[M(4), D(5), M(7)], b"ACGTCGTACGT")));
        // The read's own indels are absorbed, including the cross-kind and growing-deletion
        // transitions `rewrite_read` cannot make.
        check_hp(0, &[M(9), I(1), M(7)], &with_ins, 0, Some((&[M(16)], REF)));
        check_hp(0, &[M(9), I(1), M(7)], &with_ins, -1, Some((&[M(4), D(1), M(11)], b"ACGTAAAACGTACGT")));
        check_hp(0, &[M(4), D(1), M(11)], &[&REF[..4], &REF[5..]].concat(), -3, Some((&[M(4), D(3), M(9)], b"ACGTAACGTACGT")));
        // A run the read deleted entirely is resynthesized from the reference base.
        check_hp(0, &[M(4), D(5), M(7)], &no_run, 0, Some((&[M(16)], REF)));
        check_hp(0, &[M(4), D(5), M(7)], &no_run, 2, Some((&[M(4), I(2), M(12)], b"ACGTAAAAAAACGTACGT")));
        check_hp(0, &[M(4), D(5), M(7)], &no_run, -5, Some((&[M(4), D(5), M(7)], &no_run)));
        // Every emitted run base is `b`, never the read's own (wrong) base nor a slice of the read.
        check_hp(0, &[M(16)], &c_in_run, 1, Some((&[M(4), I(1), M(12)], b"ACGTAAAAAACGTACGT")));
        check_hp(0, &[M(9), I(1), M(7)], &c_in_run_with_ins, 2, Some((&[M(4), I(2), M(12)], b"ACGTAAAAAAACGTACGT")));
        // Soft clips are kept and the run resolves by absolute coordinate.
        check_hp(2, &[S(2), M(12)], &[b"NN", &REF[2..14]].concat(), 1, Some((&[S(2), M(2), I(1), M(10)], b"NNGTAAAAAACGTAC")));
        // Refused: a negative target length, or a read that does not cover the whole run.
        check_hp(0, &[M(16)], REF, -6, None);
        check_hp(0, &[M(6)], &REF[..6], -1, None);
        check_hp(5, &[M(11)], &REF[5..], -1, None);
    }

    #[test]
    fn rewrite_hp_run_synthesized_bases_take_the_minimum_flanking_quality() {
        let mut qual = vec![40u8; 16];
        qual[3] = 20;
        qual[4] = 30;
        let (_, _, q) = rewrite_hp_run(0, &cig(&[M(16)]), REF, &qual, 3, 9, b'A', 1).unwrap();
        // Anchor keeps its own; the synthesized base gets min(20, 30); the run bases keep theirs.
        assert_eq!(q[3..6], [20, 20, 30]);
    }

    #[test]
    fn strip_stale_tags_removes_only_invalidated_tags() {
        use rust_htslib::bam::record::{Aux, Record};
        let stale: [&[u8]; 6] = [b"MM", b"ML", b"MN", b"NM", b"MD", b"cs"];
        let mut rec = Record::new();
        rec.set(b"r1", Some(&cig(&[M(4)])), b"ACGT", &[30u8; 4]);
        for tag in stale {
            rec.push_aux(tag, Aux::I32(1)).unwrap();
        }
        rec.push_aux(b"XY", Aux::I32(7)).unwrap();
        strip_stale_tags(&mut rec);
        strip_stale_tags(&mut rec); // absent tags are a no-op
        for tag in stale {
            assert!(rec.aux(tag).is_err(), "{} must be stripped", String::from_utf8_lossy(tag));
        }
        assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)));
    }
}
