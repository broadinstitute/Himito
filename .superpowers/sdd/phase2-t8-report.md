# ONT-recall Phase 2 T8 — max-partner co-occurrence (experiment, report only)

Hidden flag `--cooccurrence-stat {sum-jaccard,max-partner}` on `call` and
`quick-start`, default `sum-jaccard` (unchanged behaviour).

## Statistic

For each tested variant `i` (same tested set as today: admitted by
`permutation_test` and matrix alt-frequency <= `permutation_frequency_threshold`),
over every other matrix row `j` with `shared = |alt_i & alt_j| >= 2`:

    p_ij = P(X >= shared),  X ~ Hypergeom(N = |cov_i & cov_j|,
                                          K = |alt_j & cov_i|, n = |alt_i & cov_j|)

`p_i = min(1, min_j p_ij * #partners)` (Bonferroni), `p_i = 1` with no partner, then
BH at `p_value_threshold` as usual. No permutation is drawn, so `--permutation-mode`
is ignored in this mode. The math is the same as `perm_retest.py --mode max-partner`,
and it uses `linkage::log10_hyper_sf`.

## Oracle (`perm_retest.py`, alpha 0.01)

| matrix | sum-Jaccard adaptive: truth / non-truth | max-partner: truth / non-truth |
|---|---|---|
| `r3_final` (merge215k ont-r10, 48k reads, 554 tested, 66 truth) | not run (too slow at 48k reads) | **66 / 164** |
| `HG00097_r10` (HapMap ~420x, 574 tested, 64 truth) | 55 / 3 | **63 / 89** |
| `r3_t4_den` (merge215k ont-denoised, 46 tested, 29 truth) | 28 / 5 | **29 / 6** |

## Rust end-to-end (merge215k `ont-r10`, `--genotype-from-bam`)

`BAM=HapMap_ont_merge.bam run_hapmap_bench.sh ont-r10 p2t8_maxp --cooccurrence-stat max-partner`

| run | TP / FP / FN | recall | precision |
|---|---|---|---|
| `r3_final` (sum-Jaccard adaptive) | 114 / 6 / 7 | 0.942 | 0.950 |
| `p2t8_maxp` (max-partner) | 112 / 18 / 9 | 0.926 | 0.862 |

By truth HL: 0.02 rose from 42 to 44 TP, which is the HG002-style recovery. 0.005 fell
from 12 to 8 TP.
New FP: 545G>A, 955AC>A, 1800GA>G, 1896TA>T, 3298CT>C, 4357T>C, 6666TC>T, 6691GA>G,
13365CG>C, 15150G>A, 16162A>G, 16386TA>T. Most are low-HF homopolymer indels and SNPs
that ride on the same error-prone reads.

## Reading

Max-partner recovers almost every truth allele, but on deep ONT data it does not
separate real lineages from correlated systematic errors. Errors that co-occur on the
same bad reads look like a strong single partner, and with 10^4 to 10^5 co-covered
reads even a modest enrichment drives the hypergeometric p to ~0. The summed Jaccard
dilutes those pairs, and that dilution is exactly what makes it a useful filter here.
Downstream filters cut the oracle's 164 non-truth passes to 18 FP, but that is still
3x the baseline. It is not worth promoting as-is. A possible follow-up: require the
partner itself to be a confident or PASS call, or require per-strand linkage the way
`linkage::rescue` does. The owner decides.

Unit tests: `sum_jaccard_rejects_a_small_lineage`,
`max_partner_keeps_a_small_lineage_and_rejects_noise`,
`max_partner_p_is_partner_count_times_min_hypergeometric_p`.
`CARGO_TARGET_DIR=target cargo test`: **442 passed / 1 ignored**.
