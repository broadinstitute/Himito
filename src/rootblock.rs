//! Partition variants into an unordered root block (homoplasmic: the molecules
//! lacking them are not distinguished by any other allele) and an informative
//! set (everything the mutation-tree search should actually order).
//!
//! See `docs/root-block-collapse.md` for why HF alone cannot make this call.

use crate::lineage::BinaryMatrix;
use crate::scite::fisher_greater;
use crate::scite::CleanedMatrix;
use adjustp::{adjust, Procedure};
use std::collections::HashMap;

/// `min_hf` is a candidate *gate*, not the decision — the decision is the
/// absence-structure test in [`partition`].
#[derive(Debug, Clone)]
pub struct RootBlockConfig {
    pub min_hf: f64,
    pub max_q: f64,
    pub min_absent: usize,
    /// Effect-size floor for rule 3: a partner's enrichment odds ratio must
    /// reach this before its absences count as *structured*. A one-sided Fisher
    /// p shrinks with read count, so significance alone would move true
    /// homoplasmies into the search as coverage rises; demanding a real odds
    /// ratio (which a genuine subclone produces but dropout does not) decouples
    /// the call from depth.
    pub min_or: f64,
}

impl Default for RootBlockConfig {
    fn default() -> Self {
        Self { min_hf: 0.80, max_q: 0.05, min_absent: 10, min_or: 2.0 }
    }
}

/// Why a candidate ended up on its side of the partition, for the audit
/// table's `reason` column. Purely descriptive — records an existing
/// decision, it never participates in making one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockReason {
    /// Admitted without testing. Two distinct causes reach this: rule 2 (too
    /// few jointly-covered absences to test at all), and no candidate partner
    /// offering a testable alt call (`partition`'s pass 3,
    /// `best.get(&v) == None`). Both mean "cannot assess structure, so admit",
    /// and both print `NA` for `min_q`/`partner`, so the audit row does not
    /// distinguish them — deliberately named for the shared outcome rather
    /// than for either cause, since naming it after one would misreport the
    /// other.
    Untested,
    /// Tested: every partner's q exceeded `max_q`. No evidence the absences
    /// are structured, so blocked.
    UnstructuredAbsences,
    /// Tested: some partner's q was at or below `max_q`. The absences are
    /// structured, so kept in the search.
    StructuredAbsences,
    /// Would have been blocked, but pulled back because its REF span
    /// overlaps another variant (`resolve_span_conflicts` in `scite.rs`).
    SpanConflict,
}

impl BlockReason {
    pub fn as_str(self) -> &'static str {
        match self {
            BlockReason::Untested => "untested",
            BlockReason::UnstructuredAbsences => "unstructured_absences",
            BlockReason::StructuredAbsences => "structured_absences",
            BlockReason::SpanConflict => "span_conflict",
        }
    }
}

/// One row per candidate, i.e. per variant clearing `min_hf`. Variants below the
/// gate are never candidates and produce no row.
#[derive(Debug, Clone)]
pub struct RootBlockAudit {
    pub variant: usize,
    pub hf: f64,
    pub n_absent: usize,
    /// `None` when admitted without testing (too few absences).
    pub min_q: Option<f64>,
    /// The variant achieving `min_q`; `None` when untested.
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
    /// The partition that leaves every variant in the search — what
    /// `--no-root-block` produces.
    pub fn disabled(n_variants: usize) -> Self {
        RootBlock { block: Vec::new(), informative: (0..n_variants).collect(), audit: Vec::new() }
    }
}

/// Heteroplasmic frequency and absent-call count for one variant, over covered
/// reads only. Mirrors `lineage.rs:515` so the two stages agree.
pub(crate) fn hf_and_absent(matrix: &BinaryMatrix, v: usize) -> (f64, usize) {
    let present = matrix.data[v].iter().filter(|c| **c == Some(1)).count();
    let absent = matrix.data[v].iter().filter(|c| **c == Some(0)).count();
    let covered = present + absent;
    let hf = if covered > 0 { present as f64 / covered as f64 } else { 0.0 };
    (hf, absent)
}

/// One-sided Fisher p for "reads absent at `v` are enriched for alt calls at
/// `u`", over reads jointly covered at both.
///
/// The 2x2, with `a` as the enriched cell `fisher_greater` tests:
///
/// |             | `u` alt | `u` ref |
/// |-------------|---------|---------|
/// | `v` absent  | a       | b       |
/// | `v` present | c       | d       |
///
/// The question is whether the molecules that LACK `v` are distinguished by
/// carrying some other allele. A real subclone at `1 - HF` answers yes; dropout
/// answers no.
///
/// Deliberately tested against `u`'s ALT calls rather than `u`'s absences:
/// dropout is partly read-level, so two germline variants drop out on the same
/// poor-quality reads and an absence-vs-absence test would call that a subclone.
/// Demanding positive distinguishing evidence is what dropout cannot fake. See
/// `read_level_dropout_confounder_still_joins_block`.
///
/// `None` when the pair carries no information: no jointly-covered absence at
/// `v`, or no alt call at `u` among jointly-covered reads. Otherwise
/// `Some((odds_ratio, p))` for the enrichment of the `v`-absent / `u`-alt cell.
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

/// `dropout_rate` is the per-read dropout (false-negative) rate used to make the
/// HF gate depth-relative: expected dropout absences (`dropout_rate x covered`)
/// are discounted before the gate is judged. `0.0` reproduces the absolute,
/// coverage-sensitive behaviour. It is a parameter rather than a
/// `RootBlockConfig` field because it is not a root-block knob: it is the SCITE
/// false-negative rate, and taking it from the one place that owns it
/// (`run_scite_pipeline`) means the two can never disagree.
pub fn partition(matrix: &BinaryMatrix, cfg: &RootBlockConfig, dropout_rate: f64) -> RootBlock {
    let n = matrix.variants.len();
    let mut block = Vec::new();
    let mut informative = Vec::new();
    let mut audit = Vec::new();

    // Pass 1: apply rules 1 and 2, and collect the candidates needing rule 3.
    // `pending` holds (variant, hf, n_absent) in ascending variant order.
    let mut pending: Vec<(usize, f64, usize)> = Vec::new();
    for v in 0..n {
        let present = matrix.data[v].iter().filter(|c| **c == Some(1)).count();
        let (hf, n_absent) = hf_and_absent(matrix, v);
        let covered = present + n_absent;

        // Dropout inflates the absent count, which drags observed HF below the
        // gate as depth grows and sheds true homoplasmies with coverage.
        // Discount the absences a pure-dropout homoplasmy would produce before
        // judging the HF gate. With `dropout_rate == 0` `corrected_hf == hf`.
        // Review 2026-09-25 T4: rule 2 uses raw `min_absent`; testability is a
        // sample-size question, not depth-relative.
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
                reason: BlockReason::Untested,
            });
        } else {
            pending.push((v, hf, n_absent));
        }
    }

    // Pass 2: every (candidate, partner) p-value and odds ratio, then ONE
    // Benjamini-Hochberg correction across all p-values together.
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

    // Two best-partner tables, both built in ascending `u` order with strict-<
    // updates so ties resolve to the lowest partner index deterministically:
    //   `best`          — smallest q over ALL partners, for the audit's min_q.
    //   `best_material`  — smallest q among partners whose odds ratio clears
    //                      `min_or`, for the decision. Only a materially large
    //                      enrichment (which a real subclone produces and
    //                      dropout cannot) can pull a variant out of the block,
    //                      so the call no longer tightens with coverage.
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
            None => (None, None, true, BlockReason::Untested),
            Some(&(bq, bu)) => {
                if structured {
                    // Report the partner that carried the structured signal, not
                    // merely the smallest-q one (they coincide unless a
                    // sub-`min_or` partner happens to have an even smaller q).
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

/// A `BinaryMatrix` holding only `keep`, in the order given. All reads survive.
pub fn restrict(matrix: &BinaryMatrix, keep: &[usize]) -> BinaryMatrix {
    BinaryMatrix {
        variants: keep.iter().map(|&v| matrix.variants[v].clone()).collect(),
        reads: matrix.reads.clone(),
        data: keep.iter().map(|&v| matrix.data[v].clone()).collect(),
    }
}

/// Append the root-block variants to a cleaned matrix as all-present columns,
/// and return the number of informative columns that preceded them.
///
/// Block variants are homoplasmic by construction, so every read carries them —
/// which is the point: the block designation corrects exactly the dropout that
/// flagged the variant in the first place.
///
/// They go at the END because `write_mutation_tree` resolves node labels with
/// `variants[node_id]` (`scite.rs:1246`) and node ids index the informative
/// matrix. Prepending would silently relabel every node.
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

    /// `data` is `[variant][read]`. `1` = alt, `0` = ref, `-1` = uncovered.
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

    /// At high depth, depth-inflated `min_absent` wrongly skips the structure
    /// test for a real subclone below the dropout rate. Rule 2 must use raw
    /// `min_absent`; dropout correction applies to the HF gate only.
    /// Review 2026-09-25 T4.
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
    fn high_hf_variant_with_too_few_absences_joins_block_untested() {
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
        assert_eq!(rb.audit[0].reason, BlockReason::Untested);
    }

    /// A germline variant whose absences are scattered: the reads lacking it
    /// are not distinguished by carrying anything else.
    #[test]
    fn germline_with_scattered_absences_joins_block() {
        // v0: 90 alt then 10 ref -> hf 0.90, n_absent 10 (testable).
        // v1: alt on every even-indexed read -> among v0's absent reads
        //     (90..100) exactly 5 are even, matching the global rate.
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
        // v0: 90 alt then 10 ref -> hf 0.90, n_absent 10 (testable).
        // v1: ref on 0..90, alt on 90..100 -> perfectly marks v0's absences.
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

    /// The effect-size floor is what stops coverage from eroding the homoplasmy
    /// set. At high depth a one-sided Fisher p reaches significance on a *small*
    /// enrichment, so q alone would move a true homoplasmy into the search. Here
    /// the absences are enriched for a partner allele (q well below `max_q`) but
    /// the odds ratio is only ~1.56: with the default `min_or = 2.0` the variant
    /// stays blocked; drop the floor to 1.0 and the same data reclassifies it.
    #[test]
    fn effect_size_floor_keeps_significant_but_weak_enrichment_in_block() {
        let n = 3000;
        let mut v0 = vec![0i8; n];
        for r in 0..2700 {
            v0[r] = 1;
        }
        // v0: hf 0.90, n_absent 300 (well past min_absent, so it is tested).
        let mut v1 = vec![0i8; n];
        // Among v0's present reads (0..2700): 540 alt. Among v0's absent reads
        // (2700..3000): 84 alt. -> a=84, b=216, c=540, d=2160: OR 1.556, p ~1e-3.
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

    /// Read-level dropout makes two germline variants drop out on the SAME bad
    /// reads. An absence-vs-absence test would read that as a subclone and keep
    /// both in the search. Testing against the partner's ALT calls does not.
    ///
    /// DO NOT "simplify" absence_vs_alt into an absence-vs-absence test; this
    /// test is the reason it is written the way it is.
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

    /// A candidate with no partner offering any alt call among jointly-covered
    /// reads cannot be tested, and is admitted rather than left in the search.
    #[test]
    fn candidate_with_no_testable_partner_joins_block() {
        // v1 is all-ref, so cell `a` (v0 absent & v1 alt) can never be non-zero
        // and the pair is skipped entirely.
        let v0 = run(90, 10);
        let v1 = vec![0i8; 100];
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default(), 0.0);
        assert_eq!(rb.block, vec![0]);
        assert!(rb.audit[0].min_q.is_none());
        // Untestable for lack of a partner, not for too few absences (v0 has
        // 10, clearing `min_absent`) — both are grouped under the same
        // "untested, admitted" reason; see `BlockReason::Untested`.
        assert_eq!(rb.audit[0].reason, BlockReason::Untested);
    }

    /// Regression guard for "one Benjamini-Hochberg correction across ALL
    /// pairs together", not per-candidate. `v0` and `v1` each have exactly
    /// one partner producing a moderate, non-trivial p (~0.0255 — neither
    /// ~1.0 nor a vanishing extreme). In isolation (m=1 each) that p survives
    /// BH untouched and both would stay informative. `v2` contributes three
    /// more pairs, all non-significant (p = 1.0), that carry no signal of
    /// their own but enlarge the pool to m=5, which is enough to push q for
    /// v0 and v1 above `max_q`.
    ///
    /// If a refactor ever moved the `adjust` call inside the per-candidate
    /// loop (pooling only each candidate's own partners), this test would
    /// start failing: v0 and v1 would wrongly stay informative. Every other
    /// fixture in this file has m=1 either way, or ties at p=1.0 exactly, so
    /// none of them would catch that regression.
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
        // v1 (row 1): block B, reads 100..200. Same shape as v0.
        for r in 100..190 {
            rows[1][r] = 1;
        }
        for r in 190..200 {
            rows[1][r] = 0;
        }
        // v2 (row 2): block C, reads 200..300. Same shape; exists only to
        // supply padding pairs below.
        for r in 200..290 {
            rows[2][r] = 1;
        }
        for r in 290..300 {
            rows[2][r] = 0;
        }

        // u0 (row 3): v0's only partner. Among v0's 10 absent reads (90..100):
        // 5 alt, 5 ref. Among v0's 90 present reads (0..90): 15 alt, 75 ref.
        // -> a=5, b=5, c=15, d=75, fisher_greater p = 0.02546 (moderate).
        // u0's own hf = 20/100 = 0.20, below the gate, so it is never itself
        // a candidate and contributes no row of its own.
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

        // u_p1..u_p3 (rows 5,6,7): padding partners for v2. a=0, b=10, c=10,
        // d=80 -> p = 1.0 exactly (a=0 saturates the one-sided Fisher test).
        // Their own hf = 10/100 = 0.10, below the gate.
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
        // Confirms the underlying per-pair p-value really is moderate, not
        // one of the trivial extremes (~0 or exactly 1.0) used elsewhere in
        // this file — this test is about BH pooling, not about the Fisher
        // test itself.
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
        // write_mutation_tree does `variants[node_id]`; node ids index the
        // INFORMATIVE matrix, so appending must not move them.
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
        // These four numbers are quoted in docs/root-block-collapse.md and in
        // the CLI help. Changing one without the others is a bug.
        let c = RootBlockConfig::default();
        assert!((c.min_hf - 0.80).abs() < 1e-12);
        assert!((c.max_q - 0.05).abs() < 1e-12);
        assert_eq!(c.min_absent, 10);
    }
}
