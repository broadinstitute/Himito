#!/usr/bin/env bash
# One full simulate -> reconstruct -> score cycle for a single read profile.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF=/Users/suhang/Analysis/Himito/rCRS.fasta

OUTDIR="" PROFILE="ont-r10" NMUT=12 DEPTH=300 SEED=1 FP=0.001 FN=0.05
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --n-mutations) NMUT="$2"; shift 2;;
    --total-depth) DEPTH="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --fp) FP="$2"; shift 2;;
    --fn) FN="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || { echo "usage: --outdir DIR [--profile ...] [--n-mutations N] [--total-depth N] [--seed N] [--fp F] [--fn F]" >&2; exit 1; }
: "${PBSIM_MODEL_DIR:?set PBSIM_MODEL_DIR}"

mkdir -p "$OUTDIR"
python "$HERE/simulate_tree.py" --reference "$REF" --n-mutations "$NMUT" --seed "$SEED" --outdir "$OUTDIR"
"$HERE/simulate_reads.sh" --outdir "$OUTDIR" --profile "$PROFILE" --total-depth "$DEPTH" --seed "$SEED"
"$HERE/run_himito.sh" --outdir "$OUTDIR" --profile "$PROFILE" --sample SIM --fp "$FP" --fn "$FN"
python "$HERE/score_lineage.py" \
  --truth-tree "$OUTDIR/truth/truth_mutation_tree.tsv" \
  --recon-tree "$OUTDIR/himito/sim_lineage.mutation_tree.tsv" \
  --truth-variants "$OUTDIR/truth/truth_variants.txt" \
  --vcf "$OUTDIR/himito/sim.vcf" \
  --profile "$PROFILE" --fp "$FP" --fn "$FN" \
  --metrics-out "$OUTDIR/metrics.tsv"

echo "=== metrics ($OUTDIR/metrics.tsv) ==="
column -t "$OUTDIR/metrics.tsv"
