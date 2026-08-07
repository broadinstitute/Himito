//! Small-indel site model for the ONT denoiser.
//!
//! Indels have no per-base quality analogue, so this module supplies the missing
//! per-read likelihood via a repeat-context-scaled error rate `eps(L)`. That one
//! function drives both the candidacy gate and the MAP assignment.

use std::collections::HashMap;

use rust_htslib::bam::record::{Cigar, CigarString};

/// An indel allele at a normalized site, indexed conceptually by its net length
/// change: `Ref` = 0, `Ins(s)` = +s.len(), `Del(n)` = -n.
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

    /// True when both alleles are the same variety (both insertions, both deletions,
    /// or both reference). Cross-kind reassignment is not supported by the rewrite
    /// walk, so the assignment step uses this to refuse those transitions.
    pub fn kind_matches(&self, other: &Allele) -> bool {
        matches!(
            (self, other),
            (Allele::Ref, Allele::Ref)
                | (Allele::Ins(_), Allele::Ins(_))
                | (Allele::Del(_), Allele::Del(_))
        )
    }
}

/// Tunable parameters for indel denoising. `Default` is DISABLED, so every existing
/// caller that passes `IndelOpts::default()` keeps byte-identical v1 behavior.
#[derive(Debug, Clone)]
pub struct IndelOpts {
    pub enabled: bool,
    /// Exclusive upper bound: with 5, lengths 1..=4 are correctable.
    pub max_len: u32,
    pub vaf: f64,
    pub err0: f64,
    pub err_scale: f64,
    pub err_cap: f64,
    pub floor_mult: f64,
    pub delta: f64,
    /// Minimum plain match-context flank required on each side of an edit.
    /// Consumers must clamp this to at least 1 wherever it is used: the FIX-4
    /// match-run bound for a gained deletion is `m + flank + 1` (anchor base +
    /// `m` deleted bases + `flank` genuine right-hand flank matches), and at
    /// `flank == 0` that bound degenerates to one short of the true
    /// no-truncation requirement, which trips a `debug_assert` in debug builds
    /// and silently truncates the deletion in release. `apply_corrections`
    /// enforces the clamp via `iopts.flank.max(1)`; this field itself is left
    /// unvalidated (a `pub` field with no constructor) so any future direct
    /// caller must not assume 0 is safe here.
    pub flank: usize,
    /// An allele carried by at least this fraction of reads at a site is never
    /// corrected away, regardless of candidacy. This guards against `vaf_floor`
    /// exceeding 1.0 in deep repeat contexts (which would otherwise make it
    /// impossible for ANY alt allele to clear candidacy, so a genuine
    /// heteroplasmy would be silently rewritten to REF read by read). Applies
    /// symmetrically: a substantial REF allele is equally protected from being
    /// converted to an indel.
    pub protect_vaf: f64,
}

impl Default for IndelOpts {
    fn default() -> Self {
        IndelOpts {
            enabled: false,
            max_len: 5,
            vaf: 0.05,
            err0: 0.01,
            err_scale: 1.5,
            err_cap: 0.4,
            floor_mult: 3.0,
            delta: 0.3,
            flank: 5,
            protect_vaf: 0.2,
        }
    }
}

/// Shift an indel maximally left so the same biological event, however the aligner
/// anchored it, collapses to one canonical `(anchor, allele)` key.
///
/// `anchor` is the reference position of the last base BEFORE the event, which is
/// exactly where htslib's pileup reports it.
pub fn normalize_left(refseq: &[u8], anchor: u32, allele: &Allele) -> (u32, Allele) {
    match allele {
        Allele::Ref => (anchor, Allele::Ref),

        // Deleted bases are refseq[a+1 ..= a+n]. Shifting the anchor left by one
        // makes them refseq[a ..= a+n-1], which is the same edit exactly when the
        // base entering the deletion equals the base leaving it.
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

        // `s` is inserted after refseq[a]. Shifting left by one puts refseq[a] after
        // the insertion instead, which is the same edit exactly when the last base of
        // `s` equals refseq[a]; the inserted string rotates right to stay equivalent.
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

/// Reference repeat context `L` at a normalized site: how many tandem copies of the
/// event's own sequence surround it. For a 1bp event this is the homopolymer run
/// length. Unique sequence gives 1. This is where ONT's indel error mass concentrates,
/// so it scales both the error rate and the candidacy floor.
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
    let u = unit.len();
    let start = norm_pos as usize + 1;

    // Copies running rightward from the event position...
    let mut copies = 0u32;
    let mut p = start;
    while refseq.get(p..p + u).map_or(false, |w| w == &unit[..]) {
        copies += 1;
        p += u;
    }
    // ...and leftward. (After left-normalization there is usually nothing to the
    // left, but counting both directions keeps the measure independent of which
    // end the aligner happened to anchor.)
    let mut p = start;
    while p >= u && refseq.get(p - u..p).map_or(false, |w| w == &unit[..]) {
        copies += 1;
        p -= u;
    }
    copies.max(1)
}

/// Per-junction indel error probability, growing with repeat context and capped.
/// This is the likelihood term that base quality supplies for substitutions and that
/// indels have no equivalent for: without it, MAP assignment degenerates to
/// "everything becomes the most frequent allele" and heteroplasmy is annihilated.
pub fn error_rate(l: u32, o: &IndelOpts) -> f64 {
    let exp = (l as i32 - 1).max(0);
    (o.err0 * o.err_scale.powi(exp)).min(o.err_cap)
}

/// Candidacy VAF floor for an allele in context `L`: an absolute minimum, raised to a
/// multiple of the local error rate wherever that is higher.
pub fn vaf_floor(l: u32, o: &IndelOpts) -> f64 {
    o.vaf.max(o.floor_mult * error_rate(l, o))
}

/// Accumulated observations at one left-normalized indel site.
#[derive(Debug, Default, Clone)]
pub struct IndelSite {
    /// Non-reference events only: allele -> (forward count, reverse count).
    pub obs: HashMap<Allele, (u32, u32)>,
    /// Total reads spanning this site (the VAF denominator).
    pub depth: u32,
    /// Column strand totals, used as the null for the strand-bias test. Mito coverage
    /// is frequently asymmetric, so the null is the column's OWN forward fraction
    /// rather than a library-wide constant.
    pub col_fwd: u32,
    pub col_rev: u32,
    /// Qnames of reads that have already cast an ALT vote at this site (see FIX 1
    /// in the task-8 review round). Left-normalization can pull two of one read's
    /// OWN events -- e.g. an insertion and a nearby deletion, exactly what
    /// coexists in a short tandem repeat -- to the SAME site; without this cap
    /// that read would contribute two ALT votes while depth counts it once,
    /// letting `alt_total` exceed `depth`.
    voted: std::collections::HashSet<Vec<u8>>,
}

impl IndelSite {
    pub fn add(&mut self, allele: Allele, reverse: bool) {
        let e = self.obs.entry(allele).or_insert((0, 0));
        if reverse {
            e.1 += 1;
        } else {
            e.0 += 1;
        }
    }

    /// Like `add`, but caps a single read to at most one ALT vote at this site.
    /// Returns `false` (and drops the vote) if `qname` has already voted here.
    pub fn add_from_read(&mut self, qname: Vec<u8>, allele: Allele, reverse: bool) -> bool {
        if !self.voted.insert(qname) {
            return false;
        }
        self.add(allele, reverse);
        true
    }

    fn alt_total(&self) -> u32 {
        self.obs.values().map(|&(f, r)| f + r).sum()
    }
}

/// The decision for one site: which alleles are kept, their renormalized frequencies,
/// and the context error rate used for MAP assignment.
#[derive(Debug, Clone)]
pub struct IndelSiteModel {
    pub kept: Vec<(Allele, f64)>,
    pub eps: f64,
    /// The repeat context length that produced `eps`; reported in
    /// `reassignments_by_hp_length` so correction inside long runs is auditable.
    pub context_len: u32,
    pub strand_biased: usize,
    /// The RAW observed VAF of every allele at the site (kept or not), including
    /// `Allele::Ref`. Unlike `kept`, this is never renormalized, so it is the
    /// correct source for the `protect_vaf` substantiality check in
    /// `assign_allele` -- candidacy rejection must not affect protection. Values
    /// are clamped to 1.0 (an allele's count can exceed `depth` on corrupted or
    /// overlapping input); the clamp bounds the raw VAF that the protection
    /// comparison above reads and makes the floor check in `fit_indel_site`
    /// slightly stricter on such input. It does not gate the near-homoplasmic
    /// exemption: an unclamped VAF above 1.0 would already exceed
    /// `homoplasmic_vaf`, so the exemption fires either way.
    pub obs_vaf: Vec<(Allele, f64)>,
    /// Alleles rejected specifically by a strand gate: the per-strand count check
    /// (`min_strand`) or the strand-bias binomial test. These are suspected
    /// single-strand artifacts -- exactly what those gates exist to catch -- so
    /// `assign_allele` must not let `protect_vaf` shield them even at a
    /// substantial VAF. An allele rejected only by the context floor (or one that
    /// is kept) is not included here and remains eligible for protection. Sorted
    /// for the same determinism reason as `kept` and `obs_vaf`.
    pub strand_rejected: Vec<Allele>,
}

/// Decide which alleles survive at a site and compute their frequencies.
///
/// REF is always kept, mirroring "the reference allele is always eligible" in the
/// substitution model. A non-reference allele must clear three gates: per-strand
/// counts, the context-scaled VAF floor, and the strand-bias binomial test (with the
/// same near-homoplasmic exemption the SNV model and `call::filter_strand_bias` use).
pub fn fit_indel_site(
    site: &IndelSite,
    refseq: &[u8],
    norm_pos: u32,
    min_strand: u32,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
    o: &IndelOpts,
) -> IndelSiteModel {
    let depth = site.depth.max(1);
    let col_total = site.col_fwd + site.col_rev;
    let expected_fwd_frac = if col_total == 0 {
        0.5
    } else {
        site.col_fwd as f64 / col_total as f64
    };

    // One eps per site. Alleles at the same site can in principle sit in different
    // repeat contexts; take the deepest. This is conservative for CANDIDACY (it
    // raises `vaf_floor`, so fewer alleles clear the gate) but anti-conservative
    // for ASSIGNMENT: the crossover threshold `eps/(1-eps)` grows with eps, so a
    // higher eps makes MAP assignment MORE willing to move reads, not less.
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

    let mut strand_biased = 0usize;
    // Raw counts over kept alleles; REF gets whatever depth the events did not claim.
    let mut raw: Vec<(Allele, f64)> = Vec::new();
    let ref_support = depth.saturating_sub(site.alt_total());
    raw.push((Allele::Ref, ref_support as f64));

    // RAW observed VAF for every allele, whether kept or not, for the protection
    // check in `assign_allele` -- rejection by candidacy must not strip protection
    // from an allele that a substantial fraction of reads actually carry.
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
        if vaf < vaf_floor(l, o) {
            continue;
        }
        // Near-homoplasmic exemption: a single-strand artifact is absent from the
        // other strand's reads, so it can never reach a high apparent frequency.
        // Testing an allele carried by most reads would correct away a real variant.
        if strand_bias_p > 0.0
            && vaf < homoplasmic_vaf
            && crate::call::strand_bias_pvalue(fwd as usize, rev as usize, expected_fwd_frac)
                < strand_bias_p
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

    // Deterministic order: `site.obs` is a HashMap with a randomly-seeded hasher,
    // so without this the two Vecs (and any downstream tie-break in
    // `assign_allele`) would vary between runs on identical input.
    kept.sort_by(|(a, _), (b, _)| a.cmp(b));
    obs_vaf.sort_by(|(a, _), (b, _)| a.cmp(b));
    strand_rejected.sort();

    IndelSiteModel { kept, eps, context_len, strand_biased, obs_vaf, strand_rejected }
}

/// Outcome of MAP assignment for one read at one site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Assignment {
    /// The read already carries the MAP allele, or the site is degenerate.
    Keep,
    /// Reassign the read to this allele.
    Move(Allele),
    /// The MAP allele is only reachable through a transition the rewrite walk does
    /// not support (cross-kind, or growing a deletion). The read is left untouched,
    /// but counted separately so the restriction's cost is auditable.
    Refused,
}

/// Maximum-a-posteriori allele for a read observing `observed`.
///
/// The `eps` likelihood is what makes this safe to run in both directions. With
/// frequency alone, `argmax_k f_k` is identical for every read and the whole site
/// collapses to its most frequent allele, annihilating heteroplasmy. Here a read only
/// leaves its own allele when `f_own / f_other < eps / (1 - eps)`.
pub fn assign_allele(observed: &Allele, m: &IndelSiteModel, o: &IndelOpts) -> Assignment {
    let eps = m.eps.clamp(0.0, 0.5);

    // Error mass spreads over the alleles that could have been observed given truth
    // k, decaying in net-length distance.
    let mut best: Option<(&Allele, f64)> = None;
    for (k, f) in &m.kept {
        let like = if k == observed {
            1.0 - eps
        } else {
            // Z normalizes over the kept set PLUS the observed allele. The observed
            // allele is very often NOT kept -- that is precisely the allele we are
            // correcting away. Summing over kept-only makes Z zero whenever REF is
            // the sole survivor, which is the single most common case in real data:
            // one read carrying a lone spurious indel. That read must still be
            // corrected, so it has to appear in its own normalization.
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
    // Restriction 1: the rewrite walk cannot turn an insertion into a deletion or
    // vice versa. Restriction 2: growing a deletion would consume reference from the
    // following CIGAR op. Both leave the read untouched, counted as Refused.
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
    // Protection: an allele carried by a substantial fraction of reads is never
    // corrected away, even when it lost candidacy (e.g. `vaf_floor` exceeding 1.0
    // in a deep repeat context, which would otherwise make every alt allele
    // uncorrectable-into and let a genuine heteroplasmy be erased read by read).
    // This is symmetric: it protects a substantial observed REF just as much as
    // a substantial observed indel. It does NOT apply to an allele rejected by a
    // strand gate (`min_strand` or the strand-bias test): that rejection means the
    // allele is a suspected single-strand artifact, which is exactly what those
    // gates exist to catch, so protection must not shield it.
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
    /// The reference position of the last aligned base BEFORE the event, exactly
    /// as the read's own CIGAR places it (before any left-normalization). This is
    /// the coordinate `rewrite_read` actually keys its `Ins`/`Del` lookups on, via
    /// `last_ref` -- it is NOT necessarily `norm_pos`. Left-normalization only
    /// ever moves an event left, so `anchor >= norm_pos`; they coincide only when
    /// the aligner already left-aligned the indel (e.g. in unique sequence).
    pub anchor: u32,
    /// `None` means the indel is at or above `max_len`: untouchable, and a marker
    /// that this read must not be assigned at nearby sites either. Omitting it
    /// entirely would make the read indistinguishable from a clean REF observation
    /// and expose it to being "repaired" into an allele it cannot represent.
    pub allele: Option<Allele>,
}

/// Length of the contiguous `M`/`=`/`X` run starting AT `target`, walking a read's
/// own CIGAR from its POS. Returns 0 when `target` is not covered by such a run at
/// all -- inside a `D`/`N`, before the read's start, or past its end.
///
/// This is the read's actual "aligned query base" test, and is stricter than
/// checking only the read's overall `[pos, end_pos())` span: that span includes any
/// `D`/`N` the read carries, but a read has no query base -- and hence no depth or
/// rewrite context -- inside those. Used both to gate an ALT vote in pass 1 (FIX 2:
/// a read that deletes the normalized column must not cast a vote there even though
/// its overall span covers it) and to size-check a gained deletion's flank in pass 2
/// (FIX 3: the deletion and its right-hand flank must fit inside ONE such run, not
/// merely inside the read's overall span).
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
///
/// The pileup reports what a read has AT A COLUMN, which cannot answer "does this
/// read carry the event belonging to the site normalization moved left to?" — at the
/// normalized column the carrier reports `Indel::None` just like a non-carrier. Going
/// through the read's own CIGAR answers it exactly, per record.
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

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                rp += n;
                qp += n as usize;
            }
            Cigar::Ins(n) => {
                // The anchor is the last aligned reference base before the event.
                // An insertion before any aligned base has no anchor and is skipped.
                if rp > 0 {
                    let anchor = rp - 1;
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
                if rp > 0 {
                    let anchor = rp - 1;
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

#[cfg(test)]
mod tests {
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
        assert_eq!(o.max_len, 5);
    }

    // ref:  0123456789
    //       C A A A A A A A G   (index 0 = C, 1..=7 = A x7, 8 = G)
    const HP_REF: &[u8] = b"CAAAAAAAG";

    #[test]
    fn deletion_in_homopolymer_normalizes_to_one_anchor() {
        // A 1bp deletion of an A can be anchored anywhere in the run (1..=6);
        // every placement must collapse to anchor 0 (the C before the run).
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
        // ref: C A T A T A T G ; inserting "TA" after index 1 (the A at 1) is the
        // same event as inserting "AT" after index 0, reached by rotating right.
        // (Verified by direct read reconstruction: ref[0..=1]+"TA"+ref[2..] ==
        // ref[0..=0]+"AT"+ref[1..] == "CATATATATG".)
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

    // Build a site: `events` is (allele, fwd_count, rev_count). REF support is
    // whatever depth remains after the events, split evenly across strands.
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
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_some(), "balanced 20% del in unique seq must be kept");
        assert!(kept_of(&m, &Allele::Ref).is_some(), "REF is always kept");
    }

    #[test]
    fn same_indel_inside_long_homopolymer_is_rejected() {
        // Identical evidence, but now in a 7-mer A-run: floor = 3 * eps(7) which is
        // far above 0.20, so the context gate rejects it.
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 10, 10)]);
        let o = IndelOpts::default();
        assert!(vaf_floor(7, &o) > 0.20, "test premise: floor must exceed the 20% VAF");
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_none(), "20% del in a 7-mer run must be rejected");
    }

    #[test]
    fn strand_skewed_indel_is_rejected_and_counted() {
        // 20 alt reads, 18/2 split, against a balanced column -> strand-biased.
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 18, 2)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_none());
        assert_eq!(m.strand_biased, 1);
    }

    #[test]
    fn near_homoplasmic_indel_is_exempt_from_strand_bias() {
        // 80% of the column, badly skewed -- but a single-strand artifact cannot
        // reach 80%, so the exemption keeps it (mirrors the SNV model).
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 75, 5)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_some(), "near-homoplasmic allele must survive");
        assert_eq!(m.strand_biased, 0);
    }

    #[test]
    fn low_count_indel_fails_the_per_strand_gate() {
        // Only one read on the reverse strand: fails min_strand=2 regardless of VAF.
        let ins = Allele::Ins(b"T".to_vec());
        let s = site(100, &[(ins.clone(), 19, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &ins).is_none());
    }

    #[test]
    fn kept_and_obs_vaf_are_returned_in_deterministic_sorted_order() {
        // Three alleles survive candidacy in unique sequence (L=1, floor 0.05): REF,
        // an insertion, and a deletion. `site.obs` is a HashMap with a randomly
        // seeded hasher, so without an explicit sort the emitted order (and any
        // downstream `>` tie-break in `assign_allele`) would vary between runs on
        // identical input, making the emitted BAM nondeterministic.
        let ins = Allele::Ins(b"A".to_vec());
        let del = Allele::Del(1);
        let s = site(100, &[(ins.clone(), 10, 10), (del.clone(), 10, 10)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert_eq!(m.kept.len(), 3, "REF, Ins, and Del must all survive candidacy here");
        // Assert the EXACT expected sequence rather than comparing the list to a
        // sorted copy of itself: `Ref` is always pushed first and is always the
        // Ord-minimum, so the only HashMap-order variability is between the two
        // remaining alleles, and one of those two permutations is already sorted --
        // a self-comparison would only catch a missing `sort_by` about half the time.
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
        // 40% INS / 60% REF in a long homopolymer (large eps). Nothing may flip.
        let ins = Allele::Ins(b"A".to_vec());
        let s = site(100, &[(ins.clone(), 20, 20)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0; // isolate assignment from the gate for this test
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, &o);
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
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, &o);
        assert!(m.eps > 0.1, "test premise: long homopolymer must have a large eps");
        assert_eq!(assign_allele(&Allele::Ref, &m, &o), Assignment::Move(ins.clone()));
        assert_eq!(assign_allele(&ins, &m, &o), Assignment::Keep);
    }

    #[test]
    fn failed_allele_snaps_to_the_kept_allele_of_the_same_kind() {
        // Site keeps a 1bp deletion at 80%; a stray 2bp deletion fails candidacy and
        // must snap onto the kept 1bp allele rather than reverting all the way to REF.
        // (The kept allele has to DOMINATE for snapping to beat reverting: at a site
        // that is mostly REF, reverting a stray deletion to REF is the correct MAP
        // answer, not a bug.)
        let d1 = Allele::Del(1);
        let d2 = Allele::Del(2);
        let s = site(100, &[(d1.clone(), 40, 40), (d2.clone(), 1, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &d1).is_some());
        assert!(kept_of(&m, &d2).is_none());
        assert_eq!(assign_allele(&d2, &m, &o), Assignment::Move(d1.clone()));
    }

    #[test]
    fn cross_kind_reassignment_is_refused() {
        // Site keeps only an insertion; a failed deletion would MAP across kinds,
        // which the rewrite walk does not support -- the read is left untouched.
        let ins = Allele::Ins(b"T".to_vec());
        let del = Allele::Del(1);
        let s = site(100, &[(ins.clone(), 45, 45), (del.clone(), 1, 1)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_none());
        // The MAP target here is Ins (it dominates the site), which is a cross-kind
        // transition the rewrite walk does not support, so it must be Refused --
        // not silently reinterpreted as Keep.
        assert_eq!(assign_allele(&del, &m, &o), Assignment::Refused);
    }

    #[test]
    fn growing_a_deletion_is_refused() {
        // Site keeps a 3bp deletion; a read carrying a 1bp deletion would have to
        // grow, consuming reference from the next CIGAR op. Refused.
        let d1 = Allele::Del(1);
        let d3 = Allele::Del(3);
        let s = site(100, &[(d3.clone(), 45, 45), (d1.clone(), 1, 1)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        // d3 dominates the site and is the MAP target, but growing Del(1) into
        // Del(3) is a transition the rewrite walk does not support -- Refused,
        // not silently reinterpreted as Keep.
        assert_eq!(assign_allele(&d1, &m, &o), Assignment::Refused);
    }

    #[test]
    fn lone_spurious_indel_at_a_ref_only_site_is_still_corrected() {
        // THE most common real case: a single read carries a random indel error, so
        // no alt allele clears candidacy and REF is the sole survivor. The stray read
        // must still be reverted. (An earlier design returned Keep whenever fewer
        // than two alleles survived, which silently disabled the denoiser's primary
        // job; the Z normalization must therefore include the observed allele.)
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 1, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
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

    // ref: C AAAAAAAAAAAA G (12 A's) -> a 1bp event here sits in L=12, where
    // vaf_floor(12) = floor_mult(3.0) * err_cap(0.4) = 1.2 > 1.0: no alt allele can
    // ever clear candidacy at this depth, regardless of true frequency.
    const DEEP_HP_REF: &[u8] = b"CAAAAAAAAAAAAG";

    #[test]
    fn substantial_allele_erased_from_candidacy_is_still_protected_at_assignment() {
        // A genuine 50/50 Del(1) heteroplasmy in a 12-mer homopolymer. Candidacy
        // genuinely rejects the allele (that is fix 1's premise, not a bug to work
        // around) -- but `protect_vaf` must stop `assign_allele` from rewriting
        // every read carrying it back to REF, which would destroy a real 50/50
        // heteroplasmy read by read (mirrors mtDNA's m.303-315 / m.16184-16193
        // poly-C heteroplasmy hotspots).
        let del = Allele::Del(1);
        let o = IndelOpts::default();
        assert!(vaf_floor(12, &o) > 1.0, "test premise: floor must exceed 1.0 at L=12");
        let s = site(100, &[(del.clone(), 25, 25)]);
        let m = fit_indel_site(&s, DEEP_HP_REF, 0, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_none(), "candidacy must genuinely reject the allele");
        assert_eq!(
            assign_allele(&del, &m, &o),
            Assignment::Keep,
            "a substantial allele must be protected from erasure, not destroyed"
        );
    }

    #[test]
    fn low_frequency_indel_below_protect_vaf_is_still_reverted() {
        // Same deep-homopolymer context (candidacy rejects the allele for the same
        // reason), but now the indel is a rare ~5% event, well below the default
        // `protect_vaf` of 0.2. Protection must NOT kick in here, proving fix 1
        // did not disable ordinary cleaning of spurious low-frequency indels.
        let del = Allele::Del(1);
        let o = IndelOpts::default();
        let s = site(100, &[(del.clone(), 3, 2)]); // 5/100 = 0.05 < protect_vaf 0.2
        let m = fit_indel_site(&s, DEEP_HP_REF, 0, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_none(), "candidacy must reject the allele");
        assert_eq!(
            assign_allele(&del, &m, &o),
            Assignment::Move(Allele::Ref),
            "a low-frequency indel below protect_vaf must still be reverted"
        );
    }

    #[test]
    fn strand_rejected_allele_is_not_shielded_by_protect_vaf() {
        // 30% VAF in unique sequence, but entirely on the forward strand (0 reverse
        // reads): a classic single-strand artifact, which `min_strand` exists to
        // catch. 30% clears `protect_vaf` (0.2), but protection must NOT apply --
        // shielding a suspected strand artifact would defeat the denoiser's main
        // purpose of removing exactly this kind of error.
        let ins = Allele::Ins(b"T".to_vec());
        let s = site(100, &[(ins.clone(), 30, 0)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
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
        // Same 30% VAF, but now strand-balanced (15f/15r) in a deep homopolymer
        // where `vaf_floor` exceeds 1.0: candidacy rejects it purely on the context
        // floor, not on a strand gate. This is not a suspected artifact, so
        // `protect_vaf` must still shield it.
        let del = Allele::Del(1);
        let o = IndelOpts::default();
        assert!(vaf_floor(12, &o) > 1.0, "test premise: floor must exceed 1.0 at L=12");
        let s = site(100, &[(del.clone(), 15, 15)]);
        let m = fit_indel_site(&s, DEEP_HP_REF, 0, 2, SB_P, HOM_VAF, &o);
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

    // Reference shared by the read_events tests. Index map:
    //   0-7   ACGTACGT
    //   8-15  AAAAAAAA   (an 8bp homopolymer run)
    //   16-28 CGTACGTACGTAC
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
        // 5M 3D 5M: targets 5..8 fall inside the deletion, and have no aligned
        // query base at all.
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
        // An insertion consumes no reference, so the two Match ops it splits are
        // still contiguous in REFERENCE space -- but `rewrite_read`'s Match arm only
        // sees ONE op at a time, so this helper (used by FIX 2/3, both of which
        // reason about a single contiguous CIGAR M/=/X op) correctly reports the
        // SECOND op's own length once the target falls past the insertion, not the
        // combined length of both flanking runs.
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
        // Read carries an extra A, anchored by the aligner at reference index 15
        // (the last A of the run). Normalization must pull it back to index 7.
        let seq = b"ACGTACGTAAAAAAAAACGTACGTACGTAC"; // 30bp: 16 + 1 inserted + 13
        let cig = cs(vec![Cigar::Match(16), Cigar::Ins(1), Cigar::Match(13)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        // The contrast this test exists for: the aligner's own CIGAR anchor (15)
        // and the left-normalized site (7) genuinely differ. `rewrite_read` keys
        // its Ins/Del lookups on `anchor`, not `norm_pos` -- see CRITICAL 1.
        assert_eq!(evs[0].anchor, 15, "the aligner's own CIGAR anchor must be reported");
        assert_eq!(evs[0].norm_pos, 7, "left-normalization must still pull the site to 7");
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"A".to_vec())));
    }

    #[test]
    fn deletion_event_is_normalized_out_of_the_homopolymer() {
        // Read is missing one A from the run; the aligner anchored the deletion at
        // reference index 14, normalization pulls it to 7.
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
        // Insert "G" after reference index 18 (a T). Nothing to shift through.
        let seq = b"ACGTACGTAAAAAAAACGTGACGTACGTAC"; // 30bp
        let cig = cs(vec![Cigar::Match(19), Cigar::Ins(1), Cigar::Match(10)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        // In unique sequence anchor and normalized site coincide -- unlike the
        // homopolymer fixtures above.
        assert_eq!(evs[0].anchor, 18);
        assert_eq!(evs[0].norm_pos, 18);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"G".to_vec())));
    }

    #[test]
    fn oversized_indel_is_reported_as_untouchable() {
        // A 6bp insertion with max_len 5: recorded with allele None so the caller
        // knows this read must not be assigned at nearby sites. Reporting nothing at
        // all would make the read look like a clean REF observation and expose it to
        // being "repaired" into some other allele.
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
        // Leading soft clip consumes query but no reference; the event must still
        // land at reference index 18.
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
}
