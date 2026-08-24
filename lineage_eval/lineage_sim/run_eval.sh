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
# Call p-value forwarded to run_himito.sh (Himito call -p); default matches that
# script, where the reason for 1 rather than 0.1 is documented: at 0.1 the caller
# emits an empty VCF for these simulated cells and every metric collapses to zero.
PVAL=1
# HF band for `Himito lineage` only (forwarded via run_himito.sh). score_lineage.py
# no longer bands its detected-variant set -- var_precision counts every PASS/. call.
# 0.01 matches the committed benchmarks, whose truth SNVs are *observed* at
# HF 0.06-0.21 and would be mostly filtered away by a 0.1 floor.
MIN_HF=0.01 MAX_HF=0.99
# Band for simulate_tree.py's *truth* frequencies. This is a different constraint
# from the lineage filter above -- one is what we simulate, the other is how we
# filter what got called -- and the two used to share a single flag, which made
# the committed benchmarks impossible to reproduce in one invocation.
#
# Every truth mutation must land above this floor to stay callable. 0.05 keeps the
# deepest n=10 mutation near 8% while still leaving internal clones big enough for
# mutation *order* to be identifiable; the old 0.1 forced near-uniform frequency
# splits, which is what starved those internal clones (see DEFAULT_INTERNAL_KEEP
# in simulate_tree.py) and also made n=12 infeasible entirely.
SIM_MIN_HF=0.05 SIM_MAX_HF=0.99
# Internal-clone mass held back by simulate_tree.py; controls whether mutation
# order along a chain is identifiable at all. Empty = use the simulator default.
INTERNAL_KEEP=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --n-mutations) NMUT="$2"; shift 2;;
    --total-depth) DEPTH="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --fp) FP="$2"; shift 2;;
    --fn) FN="$2"; shift 2;;
    --pval) PVAL="$2"; shift 2;;
    --min-hf) MIN_HF="$2"; shift 2;;
    --max-hf) MAX_HF="$2"; shift 2;;
    --sim-min-hf) SIM_MIN_HF="$2"; shift 2;;
    --sim-max-hf) SIM_MAX_HF="$2"; shift 2;;
    --internal-keep) INTERNAL_KEEP="$2"; shift 2;;
    --ref) export REF="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || {
  echo "usage: --outdir DIR [--profile ...] [--n-mutations N] [--total-depth N] [--seed N] [--fp F] [--fn F] [--pval P] [--min-hf F] [--max-hf F] [--sim-min-hf F] [--sim-max-hf F] [--internal-keep F] [--ref FASTA]" >&2
  exit 1
}
# Resolve in simulate_reads.sh; export a sane default here for visibility.
if [[ -z "${PBSIM_MODEL_DIR:-}" || ! -d "${PBSIM_MODEL_DIR:-}" ]]; then
  export PBSIM_MODEL_DIR="$HERE/pbsim3_models"
fi
[[ -f "$REF" ]] || { echo "missing reference FASTA: $REF" >&2; exit 1; }

mkdir -p "$OUTDIR"
echo "$OUTDIR"
SIM_TREE_ARGS=()
[[ -n "$INTERNAL_KEEP" ]] && SIM_TREE_ARGS+=(--internal-keep "$INTERNAL_KEEP")
python "$HERE/simulate_tree.py" --reference "$REF" --n-mutations "$NMUT" --seed "$SEED" --outdir "$OUTDIR" --min-hf "$SIM_MIN_HF" --max-hf "$SIM_MAX_HF" "${SIM_TREE_ARGS[@]+"${SIM_TREE_ARGS[@]}"}"
"$HERE/simulate_reads.sh" --outdir "$OUTDIR" --profile "$PROFILE" --total-depth "$DEPTH" --seed "$SEED"
"$HERE/run_himito.sh" --outdir "$OUTDIR" --profile "$PROFILE" --sample SIM \
  --fp "$FP" --fn "$FN" --pval "$PVAL" --min-hf "$MIN_HF" --max-hf "$MAX_HF"
python "$HERE/score_lineage.py" \
  --truth-tree "$OUTDIR/truth/truth_mutation_tree.tsv" \
  --recon-tree "$OUTDIR/himito/sim_lineage.mutation_tree.tsv" \
  --truth-variants "$OUTDIR/truth/truth_variants.txt" \
  --vcf "$OUTDIR/himito/sim.vcf" \
  --profile "$PROFILE" --fp "$FP" --fn "$FN" \
  --metrics-out "$OUTDIR/metrics.tsv"

echo "=== metrics ($OUTDIR/metrics.tsv) ==="
