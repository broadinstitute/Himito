#!/usr/bin/env bash
# Cartesian sweep of run_eval.sh over seeds × n-mutations × depths.
# Per-cell intermediates under $OUTDIR/seed{S}_mut{M}_depth{D}/; combined
# metrics (with identifying columns) at $OUTDIR/metrics.tsv.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUTDIR="" PROFILE="ont-r10" PVAL=0.1
SEEDS="" NMUTS="" DEPTHS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --seeds) SEEDS="$2"; shift 2;;
    --n-mutations) NMUTS="$2"; shift 2;;
    --depths) DEPTHS="$2"; shift 2;;
    --pval) PVAL="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" && -n "$SEEDS" && -n "$NMUTS" && -n "$DEPTHS" ]] || {
  echo "usage: --outdir DIR --seeds \"...\" --n-mutations \"...\" --depths \"...\" [--profile P] [--pval P]" >&2
  exit 1
}

mkdir -p "$OUTDIR"
COMBINED="$OUTDIR/metrics.tsv"
: > "$COMBINED"

n_s=$(wc -w <<<"$SEEDS")
n_m=$(wc -w <<<"$NMUTS")
n_d=$(wc -w <<<"$DEPTHS")
n_cells=$((n_s * n_m * n_d))
echo "sweep: ${n_cells} cells (${n_s}×${n_m}×${n_d}), profile=$PROFILE pval=$PVAL" >&2

i=0
n_ok=0
for seed in $SEEDS; do
  for nmut in $NMUTS; do
    for depth in $DEPTHS; do
      i=$((i + 1))
      cell="$OUTDIR/seed${seed}_mut${nmut}_depth${depth}"
      echo "[$i/$n_cells] seed=$seed mut=$nmut depth=$depth" >&2
      if "$HERE/run_eval.sh" \
           --outdir "$cell" \
           --profile "$PROFILE" \
           --n-mutations "$nmut" \
           --total-depth "$depth" \
           --seed "$seed" \
           --pval "$PVAL"; then
        cell_metrics="$cell/metrics.tsv"
        if [[ ! -f "$cell_metrics" || ! -s "$cell_metrics" ]]; then
          echo "missing metrics at $cell_metrics (skipped)" >&2
          continue
        fi
        header=$(head -n 1 "$cell_metrics")
        row=$(tail -n 1 "$cell_metrics")
        if [[ ! -s "$COMBINED" ]]; then
          printf 'seed\tn_mutations\ttotal_depth\t%s\n' "$header" > "$COMBINED"
        fi
        printf '%s\t%s\t%s\t%s\n' "$seed" "$nmut" "$depth" "$row" >> "$COMBINED"
        n_ok=$((n_ok + 1))
      else
        echo "run_eval failed at seed=$seed mut=$nmut depth=$depth (skipped)" >&2
      fi
    done
  done
done

[[ "$n_ok" -gt 0 ]] || {
  echo "no successful sweep cells (all run_eval runs failed?)" >&2
  exit 1
}

echo "=== combined metrics ($COMBINED) ==="
column -t "$COMBINED"
