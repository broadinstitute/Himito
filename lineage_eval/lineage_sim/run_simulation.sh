#!/usr/bin/env bash
# Cartesian sweep of run_eval.sh over seeds × n-mutations × depths.
# Per-cell intermediates under $OUTDIR/seed{S}_mut{M}_depth{D}/; combined
# metrics (with identifying columns) at $OUTDIR/metrics.tsv.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# PVAL default matches run_eval.sh / run_himito.sh, where the choice of 1 over
# 0.1 is documented: at 0.1 the caller emits an empty VCF for these cells.
OUTDIR="" PROFILE="ont-r10" PVAL=1
# SCITE rates forwarded to run_eval.sh / Himito lineage (same defaults as run_eval.sh).
FP=0.001 FN=0.05
SEEDS="" NMUTS="" DEPTHS=""
# Forwarded to run_eval.sh when set; empty = use each script's own default.
MIN_HF="" MAX_HF="" SIM_MIN_HF="" SIM_MAX_HF="" INTERNAL_KEEP=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --seeds) SEEDS="$2"; shift 2;;
    --n-mutations) NMUTS="$2"; shift 2;;
    --depths) DEPTHS="$2"; shift 2;;
    --pval) PVAL="$2"; shift 2;;
    --fp) FP="$2"; shift 2;;
    --fn) FN="$2"; shift 2;;
    --min-hf) MIN_HF="$2"; shift 2;;
    --max-hf) MAX_HF="$2"; shift 2;;
    --sim-min-hf) SIM_MIN_HF="$2"; shift 2;;
    --sim-max-hf) SIM_MAX_HF="$2"; shift 2;;
    --internal-keep) INTERNAL_KEEP="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
EVAL_ARGS=()
[[ -n "$MIN_HF" ]] && EVAL_ARGS+=(--min-hf "$MIN_HF")
[[ -n "$MAX_HF" ]] && EVAL_ARGS+=(--max-hf "$MAX_HF")
[[ -n "$SIM_MIN_HF" ]] && EVAL_ARGS+=(--sim-min-hf "$SIM_MIN_HF")
[[ -n "$SIM_MAX_HF" ]] && EVAL_ARGS+=(--sim-max-hf "$SIM_MAX_HF")
[[ -n "$INTERNAL_KEEP" ]] && EVAL_ARGS+=(--internal-keep "$INTERNAL_KEEP")
[[ -n "$OUTDIR" && -n "$SEEDS" && -n "$NMUTS" && -n "$DEPTHS" ]] || {
  echo "usage: --outdir DIR --seeds \"...\" --n-mutations \"...\" --depths \"...\" [--profile P] [--pval P] [--fp F] [--fn F] [--min-hf F] [--max-hf F] [--sim-min-hf F] [--sim-max-hf F] [--internal-keep F]" >&2
  exit 1
}

mkdir -p "$OUTDIR"
COMBINED="$OUTDIR/metrics.tsv"
: > "$COMBINED"

n_s=$(wc -w <<<"$SEEDS")
n_m=$(wc -w <<<"$NMUTS")
n_d=$(wc -w <<<"$DEPTHS")
n_cells=$((n_s * n_m * n_d))
echo "sweep: ${n_cells} cells (${n_s}×${n_m}×${n_d}), profile=$PROFILE pval=$PVAL fp=$FP fn=$FN" >&2

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
           --pval "$PVAL" \
           --fp "$FP" \
           --fn "$FN" \
           "${EVAL_ARGS[@]+"${EVAL_ARGS[@]}"}"; then
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
