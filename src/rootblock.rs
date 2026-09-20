//! Partition variants into an unordered root block (homoplasmic: the molecules
//! lacking them are not distinguished by any other allele) and an informative
//! set (everything the mutation-tree search should actually order).
//!
//! See `docs/root-block-collapse.md` for why HF alone cannot make this call.

use crate::lineage::BinaryMatrix;

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

pub fn partition(matrix: &BinaryMatrix, cfg: &RootBlockConfig) -> RootBlock {
    let mut block = Vec::new();
    let mut informative = Vec::new();
    let mut audit = Vec::new();

    for v in 0..matrix.variants.len() {
        let (hf, n_absent) = hf_and_absent(matrix, v);

        // Rule 1: below the gate, never a candidate.
        if hf < cfg.min_hf {
            informative.push(v);
            continue;
        }

        // Rule 2: too few absences to assess structure, and too few to inform
        // topology either way.
        if n_absent < cfg.min_absent {
            block.push(v);
            audit.push(RootBlockAudit {
                variant: v,
                hf,
                n_absent,
                min_q: None,
                partner: None,
                in_block: true,
            });
            continue;
        }

        // Rule 3 (absence-structure test) arrives in Task 2. Until then a
        // testable candidate stays in the search.
        informative.push(v);
        audit.push(RootBlockAudit {
            variant: v,
            hf,
            n_absent,
            min_q: None,
            partner: None,
            in_block: false,
        });
    }

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
}
