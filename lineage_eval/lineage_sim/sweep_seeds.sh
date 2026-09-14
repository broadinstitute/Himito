#!/usr/bin/env bash
# Run one full eval per seed and aggregate into a single seed_metrics.tsv.
#
# Purpose: establish the *distribution* of var_/ad_/hap_ F1 at a fixed
# configuration. A single seed cannot distinguish "the method has a systematic
# weakness" from "this instance happened to contain one coin-flip edge", and
# every conclusion about a proposed improvement needs that baseline first.
#
# No `set -e`: a seed whose frequency assignment fails in simulate_tree.py, or
# whose lineage step finds no informative variants, must not abort the sweep.
# Failures are recorded per seed and the run continues.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

OUTDIR="" PROFILE="ont-denoised" NMUT=10 DEPTH=1000
SEEDS="1 2 3 4 5 6 7 8 9 10"
# `call -f`: the 0.2 data-type default drops every truth SNV on simulated data
# and leaves an empty VCF. See README "The -f gate silently empties the VCF".
FREQ_THRESHOLD=0.05
REF_ARG="${REF:-$REPO/rCRS.fasta}"
KEEP_INTERMEDIATE=0
# Floor on simulated truth frequencies. Empty = run_eval.sh's default (0.05).
# simulate_tree.py cannot satisfy 0.05 above ~12 mutations ("Could not find valid
# frequency assignment"); n=15 needs <=0.03 and n=20 needs <=0.02. Lowering it
# pushes the deepest clones toward the caller's detection floor, so raising
# --n-mutations without lowering this is the only way to add tree difficulty
# while holding per-variant detectability fixed.
SIM_MIN_HF=""
# Forwarded to run_eval.sh. Must sit below --sim-min-hf: denoise's keep threshold
# is a frequency floor that silently removes rarer truth variants before they ever
# reach the caller, and no amount of depth compensates.
DENOISE_VAF=""
# Tree shape forwarded to run_eval.sh -> simulate_tree.py (random|chain|star).
TOPOLOGY=""
# Terminal mass each internal node retains. This, not tree shape, is what makes an
# individual ordering call hard; see DEFAULT_INTERNAL_KEEP in simulate_tree.py.
# Lowering it passes more mass downward, so it raises ordering difficulty without
# costing detection. Empty = simulator default (0.20).
INTERNAL_KEEP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --n-mutations) NMUT="$2"; shift 2;;
    --total-depth) DEPTH="$2"; shift 2;;
    --seeds) SEEDS="$2"; shift 2;;
    --frequency-threshold) FREQ_THRESHOLD="$2"; shift 2;;
    --sim-min-hf) SIM_MIN_HF="$2"; shift 2;;
    --denoise-vaf) DENOISE_VAF="$2"; shift 2;;
    --topology) TOPOLOGY="$2"; shift 2;;
    --internal-keep) INTERNAL_KEEP="$2"; shift 2;;
    --ref) REF_ARG="$2"; shift 2;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || {
  echo "usage: --outdir DIR [--profile P] [--n-mutations N] [--total-depth N] [--seeds \"1 2 3\"] [--frequency-threshold F] [--sim-min-hf F] [--denoise-vaf F] [--topology random|chain|star] [--internal-keep F] [--ref FASTA] [--keep-intermediate]" >&2
  exit 1
}
[[ -f "$REF_ARG" ]] || { echo "missing reference FASTA: $REF_ARG (pass --ref)" >&2; exit 1; }

mkdir -p "$OUTDIR"
AGG="$OUTDIR/seed_metrics.tsv"
: > "$AGG"
FAILED=""

n_seeds=$(wc -w <<<"$SEEDS" | tr -d ' ')
echo "multi-seed: ${n_seeds} seeds, profile=$PROFILE n_mutations=$NMUT depth=$DEPTH -f=$FREQ_THRESHOLD" >&2

i=0
for s in $SEEDS; do
  i=$((i + 1))
  D="$OUTDIR/seed$s"
  echo "[$i/$n_seeds] seed=$s -> $D" >&2
  rm -rf "$D"
  EVAL_ARGS=(--outdir "$D" --profile "$PROFILE"
             --n-mutations "$NMUT" --total-depth "$DEPTH" --seed "$s"
             --ref "$REF_ARG" --frequency-threshold "$FREQ_THRESHOLD")
  [[ -n "$SIM_MIN_HF" ]] && EVAL_ARGS+=(--sim-min-hf "$SIM_MIN_HF")
  [[ -n "$DENOISE_VAF" ]] && EVAL_ARGS+=(--denoise-vaf "$DENOISE_VAF")
  [[ -n "$TOPOLOGY" ]] && EVAL_ARGS+=(--topology "$TOPOLOGY")
  [[ -n "$INTERNAL_KEEP" ]] && EVAL_ARGS+=(--internal-keep "$INTERNAL_KEEP")
  if ! "$HERE/run_eval.sh" "${EVAL_ARGS[@]}" >"$OUTDIR/seed$s.log" 2>&1; then
    echo "  FAILED (see $OUTDIR/seed$s.log)" >&2
    FAILED="$FAILED $s"
    continue
  fi
  # Prepend a seed column; take the header from the first successful seed only.
  if [[ ! -s "$AGG" ]]; then
    { printf 'seed\t'; head -1 "$D/metrics.tsv"; } >> "$AGG"
  fi
  { printf '%s\t' "$s"; tail -n +2 "$D/metrics.tsv"; } >> "$AGG"
  # Reads and BAM/GFA dominate disk; drop them unless explicitly kept.
  [[ "$KEEP_INTERMEDIATE" == "1" ]] || rm -rf "$D/reads" "$D/himito/aln.sorted.bam" \
      "$D/himito/aln.denoised.bam" "$D/himito/sim.gfa"
done

[[ -s "$AGG" ]] || { echo "no seed completed successfully" >&2; exit 1; }
echo "=== per-seed metrics ($AGG) ==="
column -t "$AGG"
[[ -n "$FAILED" ]] && echo "failed seeds:$FAILED" >&2

echo "=== summary (mean / min / max over successful seeds) ==="
awk -F'\t' 'NR==1{for(j=1;j<=NF;j++) h[j]=$j; next}
  { n++; for(j=1;j<=NF;j++) if ($j ~ /^[0-9.]+$/) { s[j]+=$j; if(!(j in mn)||$j<mn[j]) mn[j]=$j; if(!(j in mx)||$j>mx[j]) mx[j]=$j } }
  END{ printf "n_seeds=%d\n", n
       for(j=1;j<=NF;j++) if (h[j] ~ /_(f1|precision|recall)$/)
         printf "  %-16s mean=%.4f  min=%.4f  max=%.4f\n", h[j], s[j]/n, mn[j], mx[j] }' "$AGG"
