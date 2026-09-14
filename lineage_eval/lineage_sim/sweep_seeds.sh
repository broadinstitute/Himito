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

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --n-mutations) NMUT="$2"; shift 2;;
    --total-depth) DEPTH="$2"; shift 2;;
    --seeds) SEEDS="$2"; shift 2;;
    --frequency-threshold) FREQ_THRESHOLD="$2"; shift 2;;
    --ref) REF_ARG="$2"; shift 2;;
    --keep-intermediate) KEEP_INTERMEDIATE=1; shift;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || {
  echo "usage: --outdir DIR [--profile P] [--n-mutations N] [--total-depth N] [--seeds \"1 2 3\"] [--frequency-threshold F] [--ref FASTA] [--keep-intermediate]" >&2
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
  if ! "$HERE/run_eval.sh" --outdir "$D" --profile "$PROFILE" \
        --n-mutations "$NMUT" --total-depth "$DEPTH" --seed "$s" \
        --ref "$REF_ARG" --frequency-threshold "$FREQ_THRESHOLD" \
        >"$OUTDIR/seed$s.log" 2>&1; then
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
