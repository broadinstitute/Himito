#!/usr/bin/env bash
# Sweep SCITE fp/fn on a fixed matrix; re-runs only `Himito lineage` per cell.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
HIMITO="${HIMITO:-$REPO/target/release/Himito}"

OUTDIR="" PROFILE="ont-r10"
FP_GRID="0.0005 0.001 0.005 0.01"
FN_GRID="0.02 0.05 0.1 0.2"
# HF band for `Himito lineage` only, same as run_eval.sh / run_himito.sh.
# score_lineage.py does not band its detected-variant set, so this does not move
# var_precision across sweep cells -- only which variants enter the SCITE matrix.
MIN_HF=0.01 MAX_HF=0.95
# Sweeps re-run full MCMC per cell; use a lighter budget than the CLI defaults
# (10000×4). Override for a final high-quality pass if needed.
MCMC_ITERS="${MCMC_ITERS:-2000}"
MCMC_CHAINS="${MCMC_CHAINS:-2}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --fp-grid) FP_GRID="$2"; shift 2;;
    --fn-grid) FN_GRID="$2"; shift 2;;
    --min-hf) MIN_HF="$2"; shift 2;;
    --max-hf) MAX_HF="$2"; shift 2;;
    --mcmc-iterations) MCMC_ITERS="$2"; shift 2;;
    --mcmc-chains) MCMC_CHAINS="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || {
  echo "usage: --outdir DIR [--profile P] [--fp-grid \"...\"] [--fn-grid \"...\"] [--min-hf F] [--max-hf F] [--mcmc-iterations N] [--mcmc-chains N]" >&2
  exit 1
}

MATRIX="$OUTDIR/himito/sim.matrix.csv"
VCF="$OUTDIR/himito/sim.vcf"
SWEEP="$OUTDIR/sweep_metrics.tsv"
[[ -f "$MATRIX" && -f "$VCF" ]] || { echo "run run_eval.sh first (missing matrix/vcf)" >&2; exit 1; }
[[ -x "$HIMITO" ]] || { echo "missing Himito binary: $HIMITO (cargo build --release)" >&2; exit 1; }
: > "$SWEEP"

n_fp=$(wc -w <<<"$FP_GRID"); n_fn=$(wc -w <<<"$FN_GRID")
n_cells=$((n_fp * n_fn))
echo "sweep: ${n_cells} cells (${n_fp}×${n_fn}), MCMC ${MCMC_CHAINS}×${MCMC_ITERS}, HF=[${MIN_HF},${MAX_HF}), himito=$HIMITO" >&2

i=0
for fp in $FP_GRID; do
  for fn in $FN_GRID; do
    i=$((i + 1))
    pfx="$OUTDIR/himito/sweep_fp${fp}_fn${fn}"
    echo "[$i/$n_cells] fp=$fp fn=$fn" >&2
    if "$HIMITO" lineage -m "$MATRIX" -v "$VCF" \
         --fp-rate "$fp" --fn-rate "$fn" --min-hf "$MIN_HF" --max-hf "$MAX_HF" \
         --mcmc-iterations "$MCMC_ITERS" --mcmc-chains "$MCMC_CHAINS" \
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
