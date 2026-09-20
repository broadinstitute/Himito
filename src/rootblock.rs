//! Partition variants into an unordered root block (homoplasmic: the molecules
//! lacking them are not distinguished by any other allele) and an informative
//! set (everything the mutation-tree search should actually order).
//!
//! See `docs/root-block-collapse.md` for why HF alone cannot make this call.

use crate::lineage::BinaryMatrix;
use crate::scite::fisher_greater;
use adjustp::{adjust, Procedure};
use std::collections::HashMap;

/// `min_hf` is a candidate *gate*, not the decision — the decision is the
/// absence-structure test in [`partition`].
#[derive(Debug, Clone)]
pub struct RootBlockConfig {
    pub min_hf: f64,
    pub max_q: f64,
    pub min_absent: usize,
}

impl Default for RootBlockConfig {
    fn default() -> Self {
        Self { min_hf: 0.80, max_q: 0.05, min_absent: 10 }
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
/// `v`, or no alt call at `u` among jointly-covered reads.
fn absence_vs_alt_p(matrix: &BinaryMatrix, v: usize, u: usize) -> Option<f64> {
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
    Some(fisher_greater(a, b, c, d).1)
}

pub fn partition(matrix: &BinaryMatrix, cfg: &RootBlockConfig) -> RootBlock {
    let n = matrix.variants.len();
    let mut block = Vec::new();
    let mut informative = Vec::new();
    let mut audit = Vec::new();

    // Pass 1: apply rules 1 and 2, and collect the candidates needing rule 3.
    // `pending` holds (variant, hf, n_absent) in ascending variant order.
    let mut pending: Vec<(usize, f64, usize)> = Vec::new();
    for v in 0..n {
        let (hf, n_absent) = hf_and_absent(matrix, v);
        if hf < cfg.min_hf {
            informative.push(v);
        } else if n_absent < cfg.min_absent {
            block.push(v);
            audit.push(RootBlockAudit {
                variant: v, hf, n_absent, min_q: None, partner: None, in_block: true,
            });
        } else {
            pending.push((v, hf, n_absent));
        }
    }

    // Pass 2: every (candidate, partner) p-value, then ONE Benjamini-Hochberg
    // correction across all of them together.
    let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
    for &(v, _, _) in &pending {
        for u in 0..n {
            if u == v {
                continue;
            }
            if let Some(p) = absence_vs_alt_p(matrix, v, u) {
                pairs.push((v, u, p));
            }
        }
    }
    let qs = if pairs.is_empty() {
        Vec::new()
    } else {
        let pvals: Vec<f64> = pairs.iter().map(|&(_, _, p)| p).collect();
        adjust(&pvals, Procedure::BenjaminiHochberg)
    };

    // Best (smallest q) partner per candidate. `pairs` is built in ascending
    // `u` order, and only a strictly smaller q displaces the incumbent, so ties
    // resolve to the lowest partner index and the result is deterministic.
    let mut best: HashMap<usize, (f64, usize)> = HashMap::new();
    for (idx, &(v, u, _)) in pairs.iter().enumerate() {
        let q = qs[idx];
        let slot = best.entry(v).or_insert((q, u));
        if q < slot.0 {
            *slot = (q, u);
        }
    }

    // Pass 3: decide each pending candidate, preserving ascending order.
    for (v, hf, n_absent) in pending {
        let (min_q, partner, in_block) = match best.get(&v) {
            // No testable partner: cannot assess structure, so admit.
            None => (None, None, true),
            Some(&(q, u)) => (Some(q), Some(u), q > cfg.max_q),
        };
        if in_block {
            block.push(v);
        } else {
            informative.push(v);
        }
        audit.push(RootBlockAudit { variant: v, hf, n_absent, min_q, partner, in_block });
    }

    block.sort_unstable();
    informative.sort_unstable();
    audit.sort_unstable_by_key(|a| a.variant);
    RootBlock { block, informative, audit }
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
        let rb = partition(&m, &RootBlockConfig::default());
        assert_eq!(rb.informative, vec![0]);
        assert!(rb.block.is_empty());
        assert!(rb.audit.is_empty(), "sub-gate variants produce no audit row");
    }

    #[test]
    fn high_hf_variant_with_too_few_absences_joins_block_untested() {
        // hf = 95/100 = 0.95 (>= 0.80), n_absent = 5 (< 10) -> rule 2
        let m = matrix(vec![run(95, 5)]);
        let rb = partition(&m, &RootBlockConfig::default());
        assert_eq!(rb.block, vec![0]);
        assert!(rb.informative.is_empty());
        assert_eq!(rb.audit.len(), 1);
        assert!(rb.audit[0].in_block);
        assert_eq!(rb.audit[0].n_absent, 5);
        assert!(rb.audit[0].min_q.is_none(), "rule 2 admits without testing");
        assert!(rb.audit[0].partner.is_none());
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
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default());
        assert_eq!(rb.block, vec![0], "unstructured absences -> block");
        assert_eq!(rb.informative, vec![1]);
        let a = &rb.audit[0];
        assert!(a.in_block);
        assert!(a.min_q.unwrap() > 0.05, "min_q was {:?}", a.min_q);
    }

    /// A real subclone: the reads lacking v0 all carry a private allele.
    #[test]
    fn true_subclone_with_private_allele_stays_informative() {
        // v0: 90 alt then 10 ref -> hf 0.90, n_absent 10 (testable).
        // v1: ref on 0..90, alt on 90..100 -> perfectly marks v0's absences.
        let v0 = run(90, 10);
        let v1: Vec<i8> = (0..100).map(|r| if r >= 90 { 1 } else { 0 }).collect();
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default());
        assert!(rb.block.is_empty(), "structured absences must not be collapsed");
        assert_eq!(rb.informative, vec![0, 1]);
        let a = &rb.audit[0];
        assert!(!a.in_block);
        assert!(a.min_q.unwrap() < 0.05, "min_q was {:?}", a.min_q);
        assert_eq!(a.partner, Some(1));
    }

    /// Read-level dropout makes two germline variants drop out on the SAME bad
    /// reads. An absence-vs-absence test would read that as a subclone and keep
    /// both in the search. Testing against the partner's ALT calls does not.
    ///
    /// DO NOT "simplify" absence_vs_alt_p into an absence-vs-absence test; this
    /// test is the reason it is written the way it is.
    #[test]
    fn read_level_dropout_confounder_still_joins_block() {
        // Reads 90..100 are poor quality: both variants drop out there.
        let v0 = run(90, 10);
        let v1 = run(90, 10);
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default());
        assert_eq!(rb.block, vec![0, 1], "correlated dropout is not a subclone");
        assert!(rb.informative.is_empty());
    }

    /// A candidate with no partner offering any alt call among jointly-covered
    /// reads cannot be tested, and is admitted rather than left in the search.
    #[test]
    fn candidate_with_no_testable_partner_joins_block() {
        // v1 is all-ref, so cell `a` (v0 absent & v1 alt) can never be non-zero
        // and the pair is skipped entirely.
        let v0 = run(90, 10);
        let v1 = vec![0i8; 100];
        let rb = partition(&matrix(vec![v0, v1]), &RootBlockConfig::default());
        assert_eq!(rb.block, vec![0]);
        assert!(rb.audit[0].min_q.is_none());
    }
}
