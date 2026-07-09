#!/usr/bin/env bash
# Sweep SCITE fp/fn on a fixed matrix; re-runs only `Himito lineage` per cell.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HIMITO="${HIMITO:-/Users/suhang/Analysis/Himito/target/release/Himito}"

OUTDIR="" PROFILE="ont-r10"
FP_GRID="0.0005 0.001 0.005 0.01"
FN_GRID="0.02 0.05 0.1 0.2"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --fp-grid) FP_GRID="$2"; shift 2;;
    --fn-grid) FN_GRID="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || { echo "usage: --outdir DIR [--profile P] [--fp-grid \"...\"] [--fn-grid \"...\"]" >&2; exit 1; }

MATRIX="$OUTDIR/himito/sim.matrix.csv"
VCF="$OUTDIR/himito/sim.vcf"
SWEEP="$OUTDIR/sweep_metrics.tsv"
[[ -f "$MATRIX" && -f "$VCF" ]] || { echo "run run_eval.sh first (missing matrix/vcf)" >&2; exit 1; }
: > "$SWEEP"

for fp in $FP_GRID; do
  for fn in $FN_GRID; do
    pfx="$OUTDIR/himito/sweep_fp${fp}_fn${fn}"
    if "$HIMITO" lineage -m "$MATRIX" -v "$VCF" \
         --fp-rate "$fp" --fn-rate "$fn" --min-hf 0.01 --max-hf 0.95 \
         -o "$pfx" >/dev/null 2>&1; then
      python "$HERE/score_lineage.py" \
        --truth-tree "$OUTDIR/truth/truth_mutation_tree.tsv" \
        --recon-tree "${pfx}.mutation_tree.tsv" \
        --truth-variants "$OUTDIR/truth/truth_variants.txt" \
        --vcf "$VCF" --profile "$PROFILE" --fp "$fp" --fn "$fn" \
        --metrics-out "$SWEEP" >/dev/null
    else
      echo "lineage failed at fp=$fp fn=$fn (skipped)" >&2
    fi
  done
done

echo "=== sweep results ($SWEEP) ==="
column -t "$SWEEP"
echo "=== best cell by ad_f1 (tie-break var_f1) ==="
[[ $(wc -l < "$SWEEP") -gt 1 ]] || { echo "no successful sweep cells (all lineage runs failed?)" >&2; exit 1; }
# columns: 1=profile 2=fp 3=fn ... 9=var_f1 ... 12=ad_f1
tail -n +2 "$SWEEP" | sort -t$'\t' -k12,12gr -k9,9gr | head -1 \
  | awk -F'\t' '{printf "profile=%s fp=%s fn=%s  ad_f1=%s var_f1=%s\n",$1,$2,$3,$12,$9}'
