#!/usr/bin/env bash
# One full simulate -> reconstruct -> score cycle for a single read profile.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
# Export so child scripts (run_himito.sh) see the same paths. A plain shell
# assignment is NOT inherited by subprocesses — that previously made Linux
# runs fall back to a hardcoded macOS REF and die at samtools sort.
export REF="${REF:-$REPO/rCRS.fasta}"
export HIMITO="${HIMITO:-$REPO/target/release/Himito}"

OUTDIR="" PROFILE="ont-r10" NMUT=12 DEPTH=300 SEED=1 FP=0.001 FN=0.05
# HF band for `Himito lineage` only (forwarded via run_himito.sh). score_lineage.py
# no longer bands its detected-variant set -- var_precision counts every PASS/. call.
MIN_HF=0.1 MAX_HF=0.99
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --n-mutations) NMUT="$2"; shift 2;;
    --total-depth) DEPTH="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --fp) FP="$2"; shift 2;;
    --fn) FN="$2"; shift 2;;
    --min-hf) MIN_HF="$2"; shift 2;;
    --max-hf) MAX_HF="$2"; shift 2;;
    --ref) export REF="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || {
  echo "usage: --outdir DIR [--profile ...] [--n-mutations N] [--total-depth N] [--seed N] [--fp F] [--fn F] [--min-hf F] [--max-hf F] [--ref FASTA]" >&2
  exit 1
}
# Resolve in simulate_reads.sh; export a sane default here for visibility.
if [[ -z "${PBSIM_MODEL_DIR:-}" || ! -d "${PBSIM_MODEL_DIR:-}" ]]; then
  export PBSIM_MODEL_DIR="$HERE/pbsim3_models"
fi
[[ -f "$REF" ]] || { echo "missing reference FASTA: $REF" >&2; exit 1; }

mkdir -p "$OUTDIR"
echo "$OUTDIR"
python "$HERE/simulate_tree.py" --reference "$REF" --n-mutations "$NMUT" --seed "$SEED" --outdir "$OUTDIR" --min-hf "$MIN_HF" --max-hf "$MAX_HF"
"$HERE/simulate_reads.sh" --outdir "$OUTDIR" --profile "$PROFILE" --total-depth "$DEPTH" --seed "$SEED"
"$HERE/run_himito.sh" --outdir "$OUTDIR" --profile "$PROFILE" --sample SIM \
  --fp "$FP" --fn "$FN" --min-hf "$MIN_HF" --max-hf "$MAX_HF"
python "$HERE/score_lineage.py" \
  --truth-tree "$OUTDIR/truth/truth_mutation_tree.tsv" \
  --recon-tree "$OUTDIR/himito/sim_lineage.mutation_tree.tsv" \
  --truth-variants "$OUTDIR/truth/truth_variants.txt" \
  --vcf "$OUTDIR/himito/sim.vcf" \
  --profile "$PROFILE" --fp "$FP" --fn "$FN" \
  --metrics-out "$OUTDIR/metrics.tsv"

echo "=== metrics ($OUTDIR/metrics.tsv) ==="
column -t "$OUTDIR/metrics.tsv"
