//! Small-indel site model for the ONT denoiser.
//!
//! Indels have no per-base quality analogue, so this module supplies the missing
//! per-read likelihood via a repeat-context-scaled error rate `eps(L)`. That one
//! function drives both the candidacy gate and the MAP assignment.

use std::collections::HashMap;

/// An indel allele at a normalized site, indexed conceptually by its net length
/// change: `Ref` = 0, `Ins(s)` = +s.len(), `Del(n)` = -n.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
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
}
