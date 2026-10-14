# T6 PacBio precision evaluation

Commit: 22f30270d3d1b08e6f7e972f40732d7669da1a8c
Date: 2026-10-14T20:33:31Z

## Full depth (HapMap_pb_merge.bam, ~6494×)

| metric | value | gate |
|---|---|---|
| Precision | 0.918 | ≥ 0.91 |
| Sensitivity | 0.926 | ≥ 0.92 |
| Padj (7 corroborated as TP) | 0.975 | ≥ 0.97 |
| TP / FP / FN | 112 / 10 / 9 | |

FPs: 7 corroborated (5532, 9342, 10092, 12731, 15150, 15622, 16162) + {1199, 3976, 4357}.

Truth called but not PASS (indel_false_threshold 0.1): 247, 285, 302, 16182 (Potential_Artifact).

## Depth series (downsample.py → quick-start -d pacbio)

See `2026-10-14-t6-depth-series.md`. At ~1000×, S matches the plan gray-zone table (0.818). At ≤100×, gray zone is empty (T=3) so PacBio SNP admission matches the prior exemption for AC≥3.

## ONT non-regression

`scripts/pb_precision/ont_regression.sh check --full` → **ONT regression: PASS** (call-only + quick-start).
