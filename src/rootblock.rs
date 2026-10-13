//! Partition variants into an unordered root block and an informative set.

use crate::lineage::BinaryMatrix;
use crate::scite::fisher_greater;
use crate::scite::CleanedMatrix;
use adjustp::{adjust, Procedure};
use std::collections::HashMap;

/// Dropout-corrected HF a variant needs to be a candidate at all (a gate, not the decision).
const MIN_HF: f64 = 0.80;
/// Benjamini-Hochberg q above which a candidate's absences count as unstructured.
const MAX_Q: f64 = 0.05;
/// Below this many absent calls a candidate joins the block untested.
const MIN_ABSENT: usize = 10;
/// Effect-size floor: a partner's enrichment odds ratio must reach this to count as structure.
const MIN_OR: f64 = 2.0;

/// Why a candidate ended up on its side of the partition (the audit table's `reason` column).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockReason {
    /// Admitted without a partner test, reported as `high_raw_vaf`.
    HighRawVaf,
    /// Tested: every partner's q exceeded `MAX_Q`.
    UnstructuredAbsences,
    /// Tested: some partner's q was at or below `MAX_Q`.
    StructuredAbsences,
    /// Absences exceed what dropout explains.
    ExcessAbsences,
    /// Would have been blocked, but pulled back because its REF span overlaps another variant.
    SpanConflict,
}

impl BlockReason {
    pub fn as_str(self) -> &'static str {
        match self {
            BlockReason::HighRawVaf => "high_raw_vaf",
            BlockReason::UnstructuredAbsences => "unstructured_absences",
            BlockReason::StructuredAbsences => "structured_absences",
            BlockReason::ExcessAbsences => "excess_absences",
            BlockReason::SpanConflict => "span_conflict",
        }
    }
}

/// One row per candidate, i.e. per variant clearing `MIN_HF`.
/// Variants below the gate are never candidates and produce no row.
#[derive(Debug, Clone)]
pub struct RootBlockAudit {
    pub variant: usize,
    pub hf: f64,
    pub n_absent: usize,
    /// `None` when admitted without testing (too few absences).
    pub min_q: Option<f64>,
    /// The variant achieving `min_q`; `None` for `HighRawVaf` and for `ExcessAbsences`.
    pub partner: Option<usize>,
    pub in_block: bool,
    pub reason: BlockReason,
}

#[derive(Debug, Clone)]
pub struct RootBlock {
    /// Variant indices into the original matrix, ascending.
    pub block: Vec<usize>,
    /// Variant indices into the original matrix, ascending.
    pub informative: Vec<usize>,
    pub audit: Vec<RootBlockAudit>,
}

/// One-sided Fisher (odds ratio, p) for "reads absent at `v` are enriched for alt calls at `u`"; `None` if untestable.
fn absence_vs_alt(matrix: &BinaryMatrix, v: usize, u: usize) -> Option<(f64, f64)> {
    let (mut a, mut b, mut c, mut d) = (0usize, 0usize, 0usize, 0usize);
    for r in 0..matrix.reads.len() {
        let (Some(vv), Some(uu)) = (matrix.data[v][r], matrix.data[u][r]) else {
            continue;
        };
        match (vv, uu) {
            (0, 1) => a += 1,
            (0, _) => b += 1,
            (_, 1) => c += 1,
            _ => d += 1,
        }
    }
    if a + b == 0 || a + c == 0 {
        return None;
    }
    Some(fisher_greater(a, b, c, d))
}

/// Split the variants into the root block and the informative set.
/// `dropout_rate` is the expected per-read dropout; 0 disables the dropout correction and the excess-absence rule.
pub fn partition(matrix: &BinaryMatrix, dropout_rate: f64) -> RootBlock {
    let n = matrix.variants.len();
    let mut block = Vec::new();
    let mut informative = Vec::new();
    let mut audit = Vec::new();

    // Pass 1: the HF gate and untested admission; the rest need the excess-absence check and the partner test.
    let mut pending: Vec<(usize, f64, usize, usize)> = Vec::new();
    for v in 0..n {
        let (present, n_absent, hf) = crate::lineage::row_hf(&matrix.data[v]);
        let covered = present + n_absent;

        // Dropout inflates the absent count, which drags observed HF below the gate as depth grows.
        let expected_dropout = dropout_rate * covered as f64;
        let corrected_absent = (n_absent as f64 - expected_dropout).max(0.0);
        let denom = present as f64 + corrected_absent;
        let corrected_hf = if denom > 0.0 { present as f64 / denom } else { 0.0 };

        if corrected_hf < MIN_HF {
            informative.push(v);
        } else if n_absent < MIN_ABSENT {
            block.push(v);
            audit.push(RootBlockAudit {
                variant: v, hf, n_absent, min_q: None, partner: None, in_block: true,
                reason: BlockReason::HighRawVaf,
            });
        } else {
            pending.push((v, hf, n_absent, covered));
        }
    }

    // Excess absences versus dropout, before any partner test.
    if dropout_rate > 0.0 {
        pending.retain(|&(v, hf, n_absent, covered)| {
            let p_excess = crate::scite::binomial_upper_tail(covered, n_absent, dropout_rate);
            let excess = p_excess <= MAX_Q && n_absent as f64 >= MIN_OR * dropout_rate * covered as f64;
            if excess {
                informative.push(v);
                audit.push(RootBlockAudit {
                    variant: v,
                    hf,
                    n_absent,
                    min_q: Some(p_excess),
                    partner: None,
                    in_block: false,
                    reason: BlockReason::ExcessAbsences,
                });
            }
            !excess
        });
    }

    // Pass 2: every p-value and odds ratio, then ONE Benjamini-Hochberg correction across all p-values together.
    let mut pairs: Vec<(usize, usize, f64, f64)> = Vec::new();
    for &(v, ..) in &pending {
        for u in 0..n {
            if u == v {
                continue;
            }
            if let Some((or, p)) = absence_vs_alt(matrix, v, u) {
                pairs.push((v, u, or, p));
            }
        }
    }
    let qs = if pairs.is_empty() {
        Vec::new()
    } else {
        let pvals: Vec<f64> = pairs.iter().map(|&(_, _, _, p)| p).collect();
        adjust(&pvals, Procedure::BenjaminiHochberg)
    };

    // Two best-partner tables, both built in ascending `u` order with strict-< updates.
    let mut best: HashMap<usize, (f64, usize)> = HashMap::new();
    let mut best_material: HashMap<usize, (f64, usize)> = HashMap::new();
    for (idx, &(v, u, or, _)) in pairs.iter().enumerate() {
        let q = qs[idx];
        let slot = best.entry(v).or_insert((q, u));
        if q < slot.0 {
            *slot = (q, u);
        }
        if or >= MIN_OR {
            let mslot = best_material.entry(v).or_insert((q, u));
            if q < mslot.0 {
                *mslot = (q, u);
            }
        }
    }

    // Pass 3: decide each pending candidate, preserving ascending order.
    for (v, hf, n_absent, _) in pending {
        let (min_q, partner, in_block, reason) = match (best.get(&v), best_material.get(&v)) {
            // No testable partner: cannot assess structure, so admit.
            (None, _) => (None, None, true, BlockReason::HighRawVaf),
            // Report the partner that carried the structured signal, not merely the smallest-q one.
            (_, Some(&(mq, mu))) if mq <= MAX_Q => (Some(mq), Some(mu), false, BlockReason::StructuredAbsences),
            (Some(&(bq, bu)), _) => (Some(bq), Some(bu), true, BlockReason::UnstructuredAbsences),
        };
        if in_block {
            block.push(v);
        } else {
            informative.push(v);
        }
        audit.push(RootBlockAudit { variant: v, hf, n_absent, min_q, partner, in_block, reason });
    }

    block.sort_unstable();
    informative.sort_unstable();
    audit.sort_unstable_by_key(|a| a.variant);
    RootBlock { block, informative, audit }
}

/// A `BinaryMatrix` holding only `keep`, in the order given.
pub fn restrict(matrix: &BinaryMatrix, keep: &[usize]) -> BinaryMatrix {
    BinaryMatrix {
        variants: keep.iter().map(|&v| matrix.variants[v].clone()).collect(),
        reads: matrix.reads.clone(),
        data: keep.iter().map(|&v| matrix.data[v].clone()).collect(),
    }
}

/// Append `block_names` to a cleaned matrix as all-present rows; returns the row count before appending.
pub fn append_block_to_cleaned(cleaned: &mut CleanedMatrix, block_names: &[String]) -> usize {
    let n_informative = cleaned.variants.len();
    let n_reads = cleaned.reads.len();
    for name in block_names {
        cleaned.variants.push(name.clone());
        cleaned.data.push(vec![1u8; n_reads]);
    }
    n_informative
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `data` is `[variant][read]`.
    fn matrix(rows: Vec<Vec<i8>>) -> BinaryMatrix {
        let n_reads = rows[0].len();
        BinaryMatrix {
            variants: (0..rows.len()).map(|i| format!("m.{i}A>G")).collect(),
            reads: (0..n_reads).map(|r| format!("r{r}")).collect(),
            data: rows
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|&c| if c < 0 { None } else { Some(c as u8) })
                        .collect()
                })
                .collect(),
        }
    }

    /// `n` alt calls then `m` ref calls.
    fn run(n_alt: usize, n_ref: usize) -> Vec<i8> {
        let mut v = vec![1i8; n_alt];
        v.extend(vec![0i8; n_ref]);
        v
    }

    #[test]
    fn variant_below_hf_gate_is_informative_and_unaudited() {
        // hf = 50/100 = 0.50, below the 0.80 gate
        let m = matrix(vec![run(50, 50)]);
        let rb = partition(&m, 0.0);
        assert_eq!(rb.informative, vec![0]);
        assert!(rb.block.is_empty());
        assert!(rb.audit.is_empty(), "sub-gate variants produce no audit row");
    }

    /// Fewer absences than dropout predicts still get the partner test: `MIN_ABSENT` counts raw absences,
    /// only the HF gate is dropout-corrected.
    #[test]
    fn subclone_below_dropout_rate_is_tested_when_min_absent_met() {
        let n = 1000;
        let v0 = run(960, 40);
        let v1: Vec<i8> = (0..n).map(|r| if r >= 960 { 1 } else { 0 }).collect();
        let rb = partition(&matrix(vec![v0, v1]), 0.05);
        assert_eq!(rb.informative, vec![0, 1]);
        assert_eq!(rb.audit[0].reason, BlockReason::StructuredAbsences);
    }

    #[test]
    fn high_hf_variant_with_too_few_absences_joins_block_as_high_raw_vaf() {
        // hf = 95/100 = 0.95 (>= 0.80), n_absent = 5 (< 10) -> admitted untested
        let m = matrix(vec![run(95, 5)]);
        let rb = partition(&m, 0.0);
        assert_eq!(rb.block, vec![0]);
        assert!(rb.informative.is_empty());
        assert_eq!(rb.audit.len(), 1);
        assert!(rb.audit[0].in_block);
        assert_eq!(rb.audit[0].n_absent, 5);
        assert!(rb.audit[0].min_q.is_none());
        assert!(rb.audit[0].partner.is_none());
        assert_eq!(rb.audit[0].reason, BlockReason::HighRawVaf);
    }

    /// A germline variant whose absences are scattered.
    #[test]
    fn germline_with_scattered_absences_joins_block() {
        // v0: 90 alt then 10 ref -> hf 0.90, n_absent 10. v1 alternates, so it is unrelated to v0's absences.
        let v0 = run(90, 10);
        let v1: Vec<i8> = (0..100).map(|r| if r % 2 == 0 { 1 } else { 0 }).collect();
        let rb = partition(&matrix(vec![v0, v1]), 0.0);
        assert_eq!(rb.block, vec![0], "unstructured absences -> block");
        assert_eq!(rb.informative, vec![1]);
        let a = &rb.audit[0];
        assert!(a.in_block);
        assert!(a.min_q.unwrap() > 0.05, "min_q was {:?}", a.min_q);
        assert_eq!(a.reason, BlockReason::UnstructuredAbsences);
    }

    /// A real subclone: the reads lacking v0 all carry a private allele.
    #[test]
    fn true_subclone_with_private_allele_stays_informative() {
        // v0: 90 alt then 10 ref -> hf 0.90, n_absent 10. v1: ref on 0..90, alt on 90..100 -> perfectly marks v0's absences.
        let v0 = run(90, 10);
        let v1: Vec<i8> = (0..100).map(|r| if r >= 90 { 1 } else { 0 }).collect();
        let rb = partition(&matrix(vec![v0, v1]), 0.0);
        assert!(rb.block.is_empty(), "structured absences must not be collapsed");
        assert_eq!(rb.informative, vec![0, 1]);
        let a = &rb.audit[0];
        assert!(!a.in_block);
        assert!(a.min_q.unwrap() < 0.05, "min_q was {:?}", a.min_q);
        assert_eq!(a.partner, Some(1));
        assert_eq!(a.reason, BlockReason::StructuredAbsences);
    }

    /// The effect-size floor is what stops coverage from eroding the homoplasmy set.
    #[test]
    fn effect_size_floor_keeps_significant_but_weak_enrichment_in_block() {
        // v0: hf 0.90, n_absent 300 (well past MIN_ABSENT, so it is tested).
        let v0 = run(2700, 300);
        // v1: alt on 540 of v0's 2700 present reads and on 84 of its 300 absent reads.
        let v1: Vec<i8> = (0..3000).map(|r| i8::from(r < 540 || (2700..2784).contains(&r))).collect();
        let rb = partition(&matrix(vec![v0, v1]), 0.0);
        assert_eq!(rb.block, vec![0], "weak enrichment must stay homoplasmic");
        assert_eq!(rb.informative, vec![1]);
        let a = &rb.audit[0];
        assert!(a.in_block);
        assert_eq!(a.reason, BlockReason::UnstructuredAbsences);
        assert!(
            a.min_q.unwrap() <= 0.05,
            "the enrichment IS significant (min_q {:?}) — only the odds-ratio floor blocks it",
            a.min_q
        );
    }

    /// Read-level dropout makes two germline variants drop out on the SAME bad reads.
    #[test]
    fn read_level_dropout_confounder_still_joins_block() {
        // Reads 90..100 are poor quality: both variants drop out there.
        let v0 = run(90, 10);
        let v1 = run(90, 10);
        let rb = partition(&matrix(vec![v0, v1]), 0.0);
        assert_eq!(rb.block, vec![0, 1], "correlated dropout is not a subclone");
        assert!(rb.informative.is_empty());
        assert!(
            rb.audit[0].min_q.unwrap() > 0.05,
            "the confounder must be TESTED and survive, not skipped as untestable"
        );
    }

    /// A candidate with no partner offering any alt call among jointly-covered reads cannot be tested.
    #[test]
    fn candidate_with_no_testable_partner_joins_block() {
        // v1 is all-ref, so cell `a` (v0 absent & v1 alt) can never be non-zero and the pair is skipped entirely.
        let v0 = run(90, 10);
        let v1 = vec![0i8; 100];
        let rb = partition(&matrix(vec![v0, v1]), 0.0);
        assert_eq!(rb.block, vec![0]);
        assert!(rb.audit[0].min_q.is_none());
        // Untestable for lack of a partner, not for too few absences.
        assert_eq!(rb.audit[0].reason, BlockReason::HighRawVaf);
    }

    /// No partner marks the absent reads, so only the excess-absence rule can keep `v0` informative.
    #[test]
    fn excess_absences_keep_variant_informative_without_a_partner() {
        let n = 1000;
        let v0 = run(850, 150);
        let v1 = vec![0i8; n];
        let rb = partition(&matrix(vec![v0, v1]), 0.05);
        assert_eq!(rb.informative, vec![0, 1]);
        assert!(rb.block.is_empty());
        let a = &rb.audit[0];
        assert!(!a.in_block);
        assert_eq!(a.reason.as_str(), "excess_absences");
        assert!(a.partner.is_none());
        let q = a.min_q.expect("excess absences record the binomial tail");
        assert!(q <= 0.05, "min_q was {q}");
    }

    #[test]
    fn modest_absence_count_stays_blocked() {
        let n = 1000;
        let v0 = run(945, 55);
        let v1 = vec![0i8; n];
        let rb = partition(&matrix(vec![v0, v1]), 0.05);
        assert_eq!(rb.block, vec![0]);
        assert_eq!(rb.informative, vec![1]);
        assert!(rb.audit[0].in_block);
        assert_ne!(rb.audit[0].reason.as_str(), "excess_absences");
    }

    /// `dropout_rate == 0` skips the excess-absence rule.
    #[test]
    fn zero_dropout_skips_the_excess_absence_rule() {
        let n = 1000;
        let v0 = run(850, 150);
        let v1 = vec![0i8; n];
        let rb = partition(&matrix(vec![v0, v1]), 0.0);
        assert_eq!(rb.block, vec![0]);
        assert_eq!(rb.informative, vec![1]);
        assert_eq!(rb.audit[0].reason, BlockReason::HighRawVaf);
    }

    /// `v0` and `v1` each have exactly one partner producing a moderate, non-trivial p.
    #[test]
    fn pooled_bh_blocks_what_per_candidate_bh_would_keep_informative() {
        // A 300-read row: uncovered except for `(call, length)` runs starting at read `start`.
        let row = |start: usize, runs: &[(i8, usize)]| {
            let mut r = vec![-1i8; start];
            for &(call, len) in runs {
                r.extend(vec![call; len]);
            }
            r.resize(300, -1);
            r
        };
        let pad = row(200, &[(1, 10), (0, 90)]);
        let rows = vec![
            // v0, v1, v2: hf 0.90, n_absent 10, each on its own 100-read block.
            row(0, &[(1, 90), (0, 10)]),
            row(100, &[(1, 90), (0, 10)]),
            row(200, &[(1, 90), (0, 10)]),
            // u0: v0's only partner; u1: the same shape on v1's block.
            row(0, &[(1, 15), (0, 75), (1, 5), (0, 5)]),
            row(100, &[(1, 15), (0, 75), (1, 5), (0, 5)]),
            // Padding partners for v2: a=0, b=10, c=10, d=80 -> p = 1.0 exactly.
            pad.clone(),
            pad.clone(),
            pad,
        ];

        let rb = partition(&matrix(rows), 0.0);
        assert_eq!(
            rb.block,
            vec![0, 1, 2],
            "pooling all 5 pairs must push v0's and v1's q past MAX_Q"
        );
        let q0 = rb.audit[0].min_q.unwrap();
        let q1 = rb.audit[1].min_q.unwrap();
        assert!(q0 > 0.05, "q0 was {q0}");
        assert!(q1 > 0.05, "q1 was {q1}");
        // Confirms the underlying per-pair p-value really is moderate.
        assert!(q0 < 0.5, "q0 was {q0}, expected a moderate (not saturated) value");
    }

    #[test]
    fn restrict_keeps_only_named_variants_in_order_and_all_reads() {
        let m = matrix(vec![run(1, 2), run(2, 1), run(3, 0)]);
        let r = restrict(&m, &[2, 0]);
        assert_eq!(r.variants, vec!["m.2A>G".to_string(), "m.0A>G".to_string()]);
        assert_eq!(r.reads.len(), 3);
        assert_eq!(r.data[0], vec![Some(1), Some(1), Some(1)]);
        assert_eq!(r.data[1], vec![Some(1), Some(0), Some(0)]);
        // Reads survive an empty informative set.
        let r = restrict(&m, &[]);
        assert!(r.variants.is_empty() && r.data.is_empty());
        assert_eq!(r.reads.len(), 3);
    }

    #[test]
    fn append_block_adds_all_present_columns_at_the_end() {
        let mut c = CleanedMatrix {
            variants: vec!["m.1A>G".to_string()],
            reads: vec!["r0".to_string(), "r1".to_string()],
            data: vec![vec![1, 0]],
            attachment: vec![0, 1],
        };
        let n_informative = append_block_to_cleaned(&mut c, &["m.9A>G".to_string()]);

        assert_eq!(n_informative, 1, "returns the informative column count");
        assert_eq!(c.variants, vec!["m.1A>G".to_string(), "m.9A>G".to_string()]);
        assert_eq!(c.data[0], vec![1, 0], "informative column untouched");
        assert_eq!(c.data[1], vec![1, 1], "block variants are present in every read");
        assert_eq!(c.attachment, vec![0, 1], "attachments are not renumbered");
    }
}
