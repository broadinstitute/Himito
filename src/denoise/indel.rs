//! Small-indel site model for the ONT denoiser.
//!
//! Indels have no per-base quality analogue, so this module supplies the missing
//! per-read likelihood via a repeat-context-scaled error rate `eps(L)`. That one
//! function drives both the candidacy gate and the MAP assignment.

use std::collections::HashMap;

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
    /// `assign_allele` -- candidacy rejection must not affect protection.
    pub obs_vaf: Vec<(Allele, f64)>,
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

    for (allele, &(fwd, rev)) in site.obs.iter() {
        let n = fwd + rev;
        obs_vaf.push((allele.clone(), (n as f64 / depth as f64).min(1.0)));
        if n == 0 {
            continue;
        }
        if fwd < min_strand || rev < min_strand {
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

    IndelSiteModel { kept, eps, context_len, strand_biased, obs_vaf }
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
    // a substantial observed indel.
    let observed_vaf = m
        .obs_vaf
        .iter()
        .find(|(a, _)| a == observed)
        .map_or(0.0, |(_, v)| *v);
    if observed_vaf >= o.protect_vaf {
        return Assignment::Keep;
    }
    Assignment::Move(target.clone())
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
        let alleles: Vec<&Allele> = m.kept.iter().map(|(a, _)| a).collect();
        let mut sorted = alleles.clone();
        sorted.sort();
        assert_eq!(alleles, sorted, "kept must be in deterministic sorted order");
        let obs_alleles: Vec<&Allele> = m.obs_vaf.iter().map(|(a, _)| a).collect();
        let mut obs_sorted = obs_alleles.clone();
        obs_sorted.sort();
        assert_eq!(obs_alleles, obs_sorted, "obs_vaf must be in deterministic sorted order");
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
}
