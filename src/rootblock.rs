//! Partition variants into an unordered root block and an informative set.

use crate::lineage::BinaryMatrix;
use crate::scite::fisher_greater;
use crate::scite::CleanedMatrix;
use adjustp::{adjust, Procedure};
use std::collections::HashMap;

/// `min_hf` is a candidate *gate*, not the decision.
#[derive(Debug, Clone)]
pub struct RootBlockConfig {
    pub min_hf: f64,
    pub max_q: f64,
    pub min_absent: usize,
    /// Effect-size floor for rule 3: the partner's enrichment odds ratio must reach this.
    pub min_or: f64,
}

impl Default for RootBlockConfig {
    fn default() -> Self {
        Self { min_hf: 0.80, max_q: 0.05, min_absent: 10, min_or: 2.0 }
    }
}

/// Why a candidate ended up on its side of the partition, for the audit table's `reason` column.
/// Purely descriptive — records an existing decision, it never participates in making one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockReason {
    /// Admitted without a partner test, reported as `high_raw_vaf`.
    HighRawVaf,
    /// Tested: every partner's q exceeded `max_q`.
    UnstructuredAbsences,
    /// Tested: some partner's q was at or below `max_q`.
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

/// One row per candidate, i.e. per variant clearing `min_hf`.
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

impl RootBlock {
    /// The partition that leaves every variant in the search — what `--no-root-block` produces.
    pub fn disabled(n_variants: usize) -> Self {
        RootBlock { block: Vec::new(), informative: (0..n_variants).collect(), audit: Vec::new() }
    }
}

/// Heteroplasmic frequency and absent-call count for one variant, over covered reads only.
pub(crate) fn hf_and_absent(matrix: &BinaryMatrix, v: usize) -> (f64, usize) {
    let present = matrix.data[v].iter().filter(|c| **c == Some(1)).count();
    let absent = matrix.data[v].iter().filter(|c| **c == Some(0)).count();
    let covered = present + absent;
    let hf = if covered > 0 { present as f64 / covered as f64 } else { 0.0 };
    (hf, absent)
}

/// One-sided Fisher p for "reads absent at `v` are enriched for alt calls at `u`".
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

/// `0.0` reproduces the absolute, coverage-sensitive behaviour.
pub fn partition(matrix: &BinaryMatrix, cfg: &RootBlockConfig, dropout_rate: f64) -> RootBlock {
    let n = matrix.variants.len();
    let mut block = Vec::new();
    let mut informative = Vec::new();
    let mut audit = Vec::new();

    // Pass 1: apply rules 1 and 2, and collect the candidates needing the excess-absence check and rule 3.
    let mut pending: Vec<(usize, f64, usize, usize)> = Vec::new();
    for v in 0..n {
        let present = matrix.data[v].iter().filter(|c| **c == Some(1)).count();
        let (hf, n_absent) = hf_and_absent(matrix, v);
        let covered = present + n_absent;

        // Dropout inflates the absent count, which drags observed HF below the gate as depth grows.
        let expected_dropout = dropout_rate * covered as f64;
        let corrected_absent = (n_absent as f64 - expected_dropout).max(0.0);
        let denom = present as f64 + corrected_absent;
        let corrected_hf = if denom > 0.0 { present as f64 / denom } else { 0.0 };

        if corrected_hf < cfg.min_hf {
            informative.push(v);
        } else if n_absent < cfg.min_absent {
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
    let pending: Vec<(usize, f64, usize)> = if dropout_rate > 0.0 {
        let mut rest = Vec::new();
        for (v, hf, n_absent, covered) in pending {
            let p_excess = crate::scite::binomial_upper_tail(covered, n_absent, dropout_rate);
            let effect_floor = cfg.min_or * dropout_rate * covered as f64;
            if p_excess <= cfg.max_q && n_absent as f64 >= effect_floor {
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
            } else {
                rest.push((v, hf, n_absent));
            }
        }
        rest
    } else {
        pending.into_iter().map(|(v, hf, n_absent, _covered)| (v, hf, n_absent)).collect()
    };

    // Pass 2: every p-value and odds ratio, then ONE Benjamini-Hochberg correction across all p-values together.
    let mut pairs: Vec<(usize, usize, f64, f64)> = Vec::new();
    for &(v, _, _) in &pending {
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
        if or >= cfg.min_or {
            let mslot = best_material.entry(v).or_insert((q, u));
            if q < mslot.0 {
                *mslot = (q, u);
            }
        }
    }

    // Pass 3: decide each pending candidate, preserving ascending order.
    for (v, hf, n_absent) in pending {
        let structured = best_material.get(&v).is_some_and(|&(q, _)| q <= cfg.max_q);
        let (min_q, partner, in_block, reason) = match best.get(&v) {
            // No testable partner: cannot assess structure, so admit.
            None => (None, None, true, BlockReason::HighRawVaf),
            Some(&(bq, bu)) => {
                if structured {
                    // Report the partner that carried the structured signal, not merely the smallest-q one.
                    let &(mq, mu) = best_material.get(&v).unwrap();
                    (Some(mq), Some(mu), false, BlockReason::StructuredAbsences)
                } else {
                    (Some(bq), Some(bu), true, BlockReason::UnstructuredAbsences)
                }
            }
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

/// Append the root-block variants to a cleaned matrix as all-present columns.
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
    use crate::lineage::BinaryMatrix;

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
    fn hf_uses_covered_reads_only_ignoring_uncovered_cells() {
        // 6 alt, 2 ref, 92 uncovered -> hf = 6/8 = 0.75, n_absent = 2
        let mut row = vec![1i8; 6];
        row.extend(vec![0i8; 2]);
        row.extend(vec![-1i8; 92]);
        let m = matrix(vec![row]);
        let (hf, n_absent) = hf_and_absent(&m, 0);
        assert!((hf - 0.75).abs() < 1e-12, "hf was {hf}");
        assert_eq!(n_absent, 2);
    }

    #[test]
    fn variant_below_hf_gate_is_informative_and_unaudited() {
        // hf = 50/100 = 0.50, below the 0.80 gate
        let m = matrix(vec![run(50, 50)]);
        let rb = partition(&m, &RootBlockConfig::default(), 0.0);
        assert_eq!(rb.informative, vec![0]);
        assert!(rb.block.is_empty());
        assert!(rb.audit.is_empty(), "sub-gate variants produce no audit row");
    }

    /// Fixed `min_absent` would skip the structure test for real subclones at high depth, so it scales.
    /// Rule 2 must use raw `min_absent`; dropout correction applies to the HF gate only.
    #[test]
    fn subclone_below_dropout_rate_is_tested_when_min_absent_met() {
        let n = 1000;
        let v0 = run(960, 40);
        let v1: Vec<i8> = (0..n).map(|r| if r >= 960 { 1 } else { 0 }).collect();
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.05);
        assert_eq!(rb.informative, vec![0, 1]);
        assert_eq!(rb.audit[0].reason, BlockReason::StructuredAbsences);
    }

    #[test]
    fn high_hf_variant_with_too_few_absences_joins_block_as_high_raw_vaf() {
        // hf = 95/100 = 0.95 (>= 0.80), n_absent = 5 (< 10) -> rule 2
        let m = matrix(vec![run(95, 5)]);
        let rb = partition(&m, &RootBlockConfig::default(), 0.0);
        assert_eq!(rb.block, vec![0]);
        assert!(rb.informative.is_empty());
        assert_eq!(rb.audit.len(), 1);
        assert!(rb.audit[0].in_block);
        assert_eq!(rb.audit[0].n_absent, 5);
        assert!(rb.audit[0].min_q.is_none(), "rule 2 admits without testing");
        assert!(rb.audit[0].partner.is_none());
        assert_eq!(rb.audit[0].reason, BlockReason::HighRawVaf);
    }

    /// A germline variant whose absences are scattered.
    #[test]
    fn germline_with_scattered_absences_joins_block() {
        // v0: 90 alt then 10 ref -> hf 0.90, n_absent 10. v1.
        let v0 = run(90, 10);
        let v1: Vec<i8> = (0..100).map(|r| if r % 2 == 0 { 1 } else { 0 }).collect();
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.0);
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
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.0);
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
        let n = 3000;
        let mut v0 = vec![0i8; n];
        for r in 0..2700 {
            v0[r] = 1;
        }
        // v0: hf 0.90, n_absent 300 (well past min_absent, so it is tested).
        let mut v1 = vec![0i8; n];
        // Among v0's present reads (0..2700): 540 alt.
        for r in 0..540 {
            v1[r] = 1;
        }
        for r in 2700..2784 {
            v1[r] = 1;
        }
        let m = matrix(vec![v0, v1]);

        let rb = partition(&m, &RootBlockConfig::default(), 0.0);
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

        // Same data, floor removed: significance alone now reclassifies it.
        let no_floor = RootBlockConfig { min_or: 1.0, ..RootBlockConfig::default() };
        let rb = partition(&m, &no_floor, 0.0);
        assert_eq!(rb.informative, vec![0, 1], "without the floor, q alone wins");
        assert!(rb.block.is_empty());
        assert_eq!(rb.audit[0].reason, BlockReason::StructuredAbsences);
    }

    /// Read-level dropout makes two germline variants drop out on the SAME bad reads.
    #[test]
    fn read_level_dropout_confounder_still_joins_block() {
        // Reads 90..100 are poor quality: both variants drop out there.
        let v0 = run(90, 10);
        let v1 = run(90, 10);
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.0);
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
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.0);
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
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.05);
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
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.05);
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
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.0);
        assert_eq!(rb.block, vec![0]);
        assert_eq!(rb.informative, vec![1]);
        assert_ne!(rb.audit[0].reason.as_str(), "excess_absences");
        assert_eq!(rb.audit[0].reason, BlockReason::HighRawVaf);
    }

    /// `v0` and `v1` each have exactly one partner producing a moderate, non-trivial p.
    #[test]
    fn pooled_bh_blocks_what_per_candidate_bh_would_keep_informative() {
        let n_reads = 300;
        let mut rows: Vec<Vec<i8>> = vec![vec![-1i8; n_reads]; 8];

        // v0 (row 0): block A, reads 0..100. hf 0.90, n_absent 10.
        for r in 0..90 {
            rows[0][r] = 1;
        }
        for r in 90..100 {
            rows[0][r] = 0;
        }
        // v1 (row 1): block B, reads 100..200.
        for r in 100..190 {
            rows[1][r] = 1;
        }
        for r in 190..200 {
            rows[1][r] = 0;
        }
        // v2 (row 2): block C, reads 200..300.
        for r in 200..290 {
            rows[2][r] = 1;
        }
        for r in 290..300 {
            rows[2][r] = 0;
        }

        // u0 (row 3): v0's only partner.
        for r in 0..15 {
            rows[3][r] = 1;
        }
        for r in 15..90 {
            rows[3][r] = 0;
        }
        for r in 90..95 {
            rows[3][r] = 1;
        }
        for r in 95..100 {
            rows[3][r] = 0;
        }
        // u1 (row 4): the same shape as u0, shifted onto v1's block (100..200).
        for r in 100..115 {
            rows[4][r] = 1;
        }
        for r in 115..190 {
            rows[4][r] = 0;
        }
        for r in 190..195 {
            rows[4][r] = 1;
        }
        for r in 195..200 {
            rows[4][r] = 0;
        }

        // u_p1..u_p3 (rows 5,6,7): padding partners for v2. a=0, b=10, c=10, d=80 -> p = 1.0 exactly.
        for k in 5..8 {
            for r in 200..210 {
                rows[k][r] = 1;
            }
            for r in 210..300 {
                rows[k][r] = 0;
            }
        }

        let rb = partition(&matrix(rows), &RootBlockConfig::default(), 0.0);
        assert_eq!(
            rb.block,
            vec![0, 1, 2],
            "pooling all 5 pairs must push v0's and v1's q past max_q"
        );
        let q0 = rb.audit[0].min_q.unwrap();
        let q1 = rb.audit[1].min_q.unwrap();
        assert!(q0 > 0.05, "q0 was {q0}");
        assert!(q1 > 0.05, "q1 was {q1}");
        // Confirms the underlying per-pair p-value really is moderate.
        assert!(q0 < 0.5, "q0 was {q0}, expected a moderate (not saturated) value");
    }

    use crate::scite::CleanedMatrix;

    #[test]
    fn restrict_keeps_only_named_variants_in_order_and_all_reads() {
        let m = matrix(vec![run(1, 2), run(2, 1), run(3, 0)]);
        let r = restrict(&m, &[2, 0]);
        assert_eq!(r.variants, vec!["m.2A>G".to_string(), "m.0A>G".to_string()]);
        assert_eq!(r.reads.len(), 3);
        assert_eq!(r.data[0], vec![Some(1), Some(1), Some(1)]);
        assert_eq!(r.data[1], vec![Some(1), Some(0), Some(0)]);
    }

    #[test]
    fn restrict_to_nothing_yields_zero_variants_but_keeps_reads() {
        let m = matrix(vec![run(1, 2)]);
        let r = restrict(&m, &[]);
        assert!(r.variants.is_empty());
        assert!(r.data.is_empty());
        assert_eq!(r.reads.len(), 3, "reads survive an empty informative set");
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

    #[test]
    fn append_block_preserves_index_alignment_with_tree_node_ids() {
        // write_mutation_tree does `variants[node_id]`, so node ids must index the appended block.
        let mut c = CleanedMatrix {
            variants: vec!["m.1A>G".to_string(), "m.2A>G".to_string()],
            reads: vec!["r0".to_string()],
            data: vec![vec![1], vec![0]],
            attachment: vec![0],
        };
        append_block_to_cleaned(&mut c, &["m.9A>G".to_string()]);
        assert_eq!(c.variants[0], "m.1A>G", "node 0 still resolves to its variant");
        assert_eq!(c.variants[1], "m.2A>G", "node 1 still resolves to its variant");
    }

    #[test]
    fn default_config_matches_the_documented_values() {
        // These four numbers are quoted in docs/root-block-collapse.md and in the CLI help.
        let c = RootBlockConfig::default();
        assert!((c.min_hf - 0.80).abs() < 1e-12);
        assert!((c.max_q - 0.05).abs() < 1e-12);
        assert_eq!(c.min_absent, 10);
    }
}
