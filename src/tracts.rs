//! Known hypervariable homopolymer / short-tandem tracts of the rCRS.
//!
//! Plan 2026-09-26 T5. Calls inside one of these runs are spelled in several
//! equivalent ways against the reference (`m.16183A>AC` and `m.16189T>TC` are
//! the same C-run length once `m.16189T>C` is in the background), so they are
//! collapsed to a single sequence-level allele per read before the tree search.

use crate::scite::VariantSpan;

/// A hypervariable tract on rCRS, 1-based inclusive.
///
/// `ref_seq` is the rCRS sequence over `[start, end]`; it is checked against
/// every variant's REF before collapsing, so a non-rCRS reference silently
/// disables this. `start` is the anchor base *before* the run, so left-aligned
/// insertions (`m.302A>AC`, `m.16183A>AC`, `m.513G>GCA`) fall inside.
pub struct Tract {
    pub name: &'static str,
    pub start: u64,
    pub end: u64,
    pub ref_seq: &'static str,
}

pub const KNOWN_TRACTS: &[Tract] = &[
    Tract { name: "D310", start: 302, end: 315, ref_seq: "ACCCCCCCTCCCCC" },
    Tract { name: "HVS1-C", start: 16183, end: 16193, ref_seq: "ACCCCCTCCCC" },
    Tract { name: "CA514", start: 513, end: 525, ref_seq: "GCACACACACACC" },
];

/// A tract as a half-open interval, matching [`VariantSpan`] conventions: the
/// REF interval of a variant at the last base of the tract is `[end, end + 1)`.
fn tract_span(tract: &Tract) -> VariantSpan {
    VariantSpan { start: tract.start, end: tract.end + 1 }
}

/// Index into `tracts` of the tract that fully contains `span`, if any.
pub fn containing_tract_in(tracts: &[Tract], span: &VariantSpan) -> Option<usize> {
    tracts.iter().position(|t| {
        let t = tract_span(t);
        span.start >= t.start && span.end <= t.end
    })
}

/// Index of the [`KNOWN_TRACTS`] entry that fully contains `span`, if any.
#[allow(dead_code)] // the KNOWN_TRACTS shorthand; T5 goes through `*_in`.
pub fn containing_tract(span: &VariantSpan) -> Option<usize> {
    containing_tract_in(KNOWN_TRACTS, span)
}

/// Whether `span` overlaps any tract in `tracts`.
pub fn overlaps_tract_in(tracts: &[Tract], span: &VariantSpan) -> bool {
    tracts.iter().any(|t| tract_span(t).overlaps(span))
}

/// Whether `span` overlaps any [`KNOWN_TRACTS`] entry (used by T6).
#[allow(dead_code)] // consumed by plan T6 (per-locus SCITE error rates).
pub fn overlaps_tract(span: &VariantSpan) -> bool {
    overlaps_tract_in(KNOWN_TRACTS, span)
}

/// The matrix row name for tract allele `alt`: `m.<start><ref_seq>><alt>`.
///
/// `scite::variant_spans` parses this as `[start, end + 1)`, i.e. the whole
/// tract, so two alleles of one tract are automatically mutually exclusive in
/// the tree search.
pub fn tract_allele_name(tract: &Tract, alt: &str) -> String {
    format!("m.{}{}>{}", tract.start, tract.ref_seq, alt)
}

/// Whether `name` is a tract-allele row name produced by [`tract_allele_name`].
pub fn is_tract_allele(name: &str, tracts: &[Tract]) -> bool {
    tracts.iter().any(|t| {
        name.strip_prefix(&format!("m.{}{}>", t.start, t.ref_seq))
            .is_some_and(|alt| !alt.is_empty())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn span(variant: &str) -> VariantSpan {
        crate::scite::variant_spans(&[variant.to_string()])[0]
    }

    #[test]
    fn known_tract_ref_seqs_match_rcrs() {
        let fasta = std::fs::read_to_string("test_data/rCRS.fasta")
            .expect("test_data/rCRS.fasta must be readable");
        let seq: String = fasta
            .lines()
            .filter(|l| !l.starts_with('>'))
            .flat_map(|l| l.trim().chars())
            .collect();
        assert_eq!(seq.len(), 16569, "rCRS is 16569 bp");
        for tract in KNOWN_TRACTS {
            let expected = &seq[(tract.start as usize - 1)..(tract.end as usize)];
            assert_eq!(
                tract.ref_seq, expected,
                "{} ({}-{}) does not match rCRS",
                tract.name, tract.start, tract.end
            );
        }
    }

    #[test]
    fn containing_tract_accepts_variants_inside_and_rejects_those_outside() {
        // Anchor insertion at the first base, SNV in the middle, SNV at the
        // last base: all inside HVS1-C (index 1).
        assert_eq!(containing_tract(&span("m.16183A>AC")), Some(1));
        assert_eq!(containing_tract(&span("m.16189T>C")), Some(1));
        assert_eq!(containing_tract(&span("m.16193C>T")), Some(1));
        // A deletion wholly inside the tract.
        assert_eq!(containing_tract(&span("m.16192CC>C")), Some(1));
        // The other two tracts.
        assert_eq!(containing_tract(&span("m.302A>AC")), Some(0));
        assert_eq!(containing_tract(&span("m.513G>GCA")), Some(2));
        // One base before the anchor, one base past the end, and far away.
        assert_eq!(containing_tract(&span("m.16182A>C")), None);
        assert_eq!(containing_tract(&span("m.16194A>G")), None);
        assert_eq!(containing_tract(&span("m.3243A>G")), None);
        // Straddling the 3' boundary: REF covers 16193 and 16194.
        assert_eq!(containing_tract(&span("m.16193CA>C")), None);
    }

    #[test]
    fn overlaps_tract_is_true_for_straddlers_and_false_outside() {
        assert!(overlaps_tract(&span("m.16189T>C")));
        assert!(overlaps_tract(&span("m.16193CA>C")));
        assert!(overlaps_tract(&span("m.302A>AC")));
        assert!(!overlaps_tract(&span("m.16182A>C")));
        assert!(!overlaps_tract(&span("m.16194A>G")));
        assert!(!overlaps_tract(&span("m.3243A>G")));
    }

    #[test]
    fn is_tract_allele_recognizes_generated_names_only() {
        let name = tract_allele_name(&KNOWN_TRACTS[1], "ACCCCCCCCCC");
        assert_eq!(name, "m.16183ACCCCCTCCCC>ACCCCCCCCCC");
        assert!(is_tract_allele(&name, KNOWN_TRACTS));
        assert!(is_tract_allele(
            &tract_allele_name(&KNOWN_TRACTS[0], "ACCCCCCCTCCCC"),
            KNOWN_TRACTS
        ));
        // Constituent calls inside a tract are not tract alleles.
        assert!(!is_tract_allele("m.16189T>C", KNOWN_TRACTS));
        assert!(!is_tract_allele("m.16183A>AC", KNOWN_TRACTS));
        // Right position, wrong REF: not a tract allele either.
        assert!(!is_tract_allele("m.16183ACCCCCTCCC>ACCCC", KNOWN_TRACTS));
    }
}
