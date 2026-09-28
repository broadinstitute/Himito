#!/usr/bin/env bash
# Sweep Himito lineage --hp-fp-rate x --hp-fn-rate on a fixed matrix.
# Plan 2026-09-27 T11. Same shape as sweep_fpfn.sh: simulation, alignment,
# denoise, build, and call run once per seed; only lineage and score_lineage.py
# re-run per grid cell.
#
# Defaults are multiples of the profile's resolved global fp/fn:
#   {1,2,4,8,16} x fp and {1,2,4,8,16,32} x fn, each capped at 0.5.
# Global --fp-rate/--fn-rate stay at that profile default unless --fp/--fn
# are given. Pass --seeds to loop seed*/ cells from sweep_seeds.sh
# --keep-intermediate (or from run_eval.sh --stop-after call, one directory
# per seed). Without --seeds, --outdir itself is the one cell.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
HIMITO="${HIMITO:-$REPO/target/release/Himito}"

# Plan 2026-09-27 T6. truth/background.txt stores Himito SNV names. The scorer
# wants the <REF><pos><ALT> tokens those names came from.
background_spec_from_truth() {
  local file="$1"
  [[ -f "$file" ]] || return 0
  local tokens=() line
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" ]] && continue
    if [[ "$line" =~ ^m\.([0-9]+)([ACGT])\>([ACGT])$ ]]; then
      tokens+=("${BASH_REMATCH[2]}${BASH_REMATCH[1]}${BASH_REMATCH[3]}")
    fi
  done < "$file"
  (IFS=,; echo "${tokens[*]}")
}

# Mirrors lineage::resolve_error_rates (src/lineage.rs) and run_eval.sh.
resolve_scite_rates() {
  case "$1" in
    hifi) echo "0.005 0.05";;
    ont-denoised) echo "0.0001 0.01";;
    *) echo "0.001 0.05";;  # ont-r10, ont-r9, fallback
  esac
}

profile_dtype() {
  case "$1" in
    hifi) echo "pacbio";;
    ont-r10) echo "ont-r10";;
    ont-denoised) echo "ont-denoised";;
    *) echo "profile must be hifi or ont-r10 or ont-denoised" >&2; return 1;;
  esac
}

# $1 base rate, remaining args are integer multipliers. Cap each product at 0.5
# and drop duplicates so a multiplier that overflows does not repeat 0.5.
default_grid() {
  python - "$1" "${@:2}" <<'PY'
import sys
base = float(sys.argv[1])
seen = []
for raw in sys.argv[2:]:
    value = min(0.5, base * int(raw))
    text = f"{value:.10g}"
    if text not in seen:
        seen.append(text)
print(" ".join(seen))
PY
}

OUTDIR="" PROFILE="ont-denoised"
HP_FP_GRID="" HP_FN_GRID=""
FP="" FN=""
MIN_HF=0.01 MAX_HF=0.95
BACKGROUND=""
COLLAPSE_HP=0
SEEDS=""
OBJECTIVE="hp_intrusion"
SECONDARY="ad_f1"
MCMC_ITERS="${MCMC_ITERS:-2000}"
MCMC_CHAINS="${MCMC_CHAINS:-2}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --hp-fp-grid) HP_FP_GRID="$2"; shift 2;;
    --hp-fn-grid) HP_FN_GRID="$2"; shift 2;;
    --fp) FP="$2"; shift 2;;
    --fn) FN="$2"; shift 2;;
    --min-hf) MIN_HF="$2"; shift 2;;
    --max-hf) MAX_HF="$2"; shift 2;;
    --mcmc-iterations) MCMC_ITERS="$2"; shift 2;;
    --mcmc-chains) MCMC_CHAINS="$2"; shift 2;;
    --background) BACKGROUND="$2"; shift 2;;
    --collapse-hp) COLLAPSE_HP=1; shift;;
    --seeds) SEEDS="$2"; shift 2;;
    --objective) OBJECTIVE="$2"; shift 2;;
    --secondary) SECONDARY="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || {
  echo "usage: --outdir DIR [--profile P] [--hp-fp-grid \"...\"] [--hp-fn-grid \"...\"] [--fp F] [--fn F] [--min-hf F] [--max-hf F] [--mcmc-iterations N] [--mcmc-chains N] [--background SPEC] [--collapse-hp] [--seeds \"1 2\"] [--objective hp_intrusion] [--secondary ad_f1]" >&2
  exit 1
}

read -r DEFAULT_FP DEFAULT_FN < <(resolve_scite_rates "$PROFILE")
DTYPE="$(profile_dtype "$PROFILE")"
GRID_FP="${FP:-$DEFAULT_FP}"
GRID_FN="${FN:-$DEFAULT_FN}"
[[ -n "$HP_FP_GRID" ]] || HP_FP_GRID="$(default_grid "$GRID_FP" 1 2 4 8 16)"
[[ -n "$HP_FN_GRID" ]] || HP_FN_GRID="$(default_grid "$GRID_FN" 1 2 4 8 16 32)"
SCORE_FP="${FP:-$DEFAULT_FP}"
SCORE_FN="${FN:-$DEFAULT_FN}"
if [[ "$COLLAPSE_HP" == 1 ]]; then COLLAPSE_TXT=true; else COLLAPSE_TXT=false; fi

[[ -x "$HIMITO" ]] || { echo "missing Himito binary: $HIMITO (cargo build --release)" >&2; exit 1; }

# seed|directory. Without --seeds the outdir is the single cell (seed NA).
CELLS=()
if [[ -n "$SEEDS" ]]; then
  for s in $SEEDS; do
    CELLS+=("$s|$OUTDIR/seed$s")
  done
else
  CELLS+=("NA|$OUTDIR")
fi

n_fp=$(wc -w <<<"$HP_FP_GRID"); n_fn=$(wc -w <<<"$HP_FN_GRID")
n_rates=$((n_fp * n_fn))
n_cells=$((n_rates * ${#CELLS[@]}))
echo "hp sweep: ${n_cells} lineage runs (${#CELLS[@]} seed(s) x ${n_fp} x ${n_fn}), global fp=$SCORE_FP fn=$SCORE_FN, MCMC ${MCMC_CHAINS}x${MCMC_ITERS}, collapse_hp=$COLLAPSE_TXT" >&2
echo "  hp_fp_grid: $HP_FP_GRID" >&2
echo "  hp_fn_grid: $HP_FN_GRID" >&2

SWEEP="$OUTDIR/hp_sweep_metrics.tsv"
: > "$SWEEP"
first=1

for entry in "${CELLS[@]}"; do
  seed="${entry%%|*}"
  cell="${entry#*|}"
  MATRIX="$cell/himito/sim.matrix.csv"
  VCF="$cell/himito/sim.vcf"
  if [[ ! -f "$MATRIX" || ! -f "$VCF" ]]; then
    echo "missing matrix/vcf in $cell (run run_eval.sh --stop-after call, or sweep_seeds.sh --keep-intermediate)" >&2
    exit 1
  fi
  cell_background="$BACKGROUND"
  if [[ -z "$cell_background" || "$cell_background" == "none" ]]; then
    derived="$(background_spec_from_truth "$cell/truth/background.txt")"
    [[ -n "$derived" ]] && cell_background="$derived"
  fi
  for hp_fp in $HP_FP_GRID; do
    for hp_fn in $HP_FN_GRID; do
      pfx="$cell/himito/sweep_hpfp${hp_fp}_hpfn${hp_fn}"
      echo "seed=$seed hp_fp=$hp_fp hp_fn=$hp_fn" >&2
      LINEAGE_ARGS=(
        -m "$MATRIX" -v "$VCF"
        --hp-fp-rate "$hp_fp" --hp-fn-rate "$hp_fn"
        --min-hf "$MIN_HF" --max-hf "$MAX_HF"
        --mcmc-iterations "$MCMC_ITERS" --mcmc-chains "$MCMC_CHAINS"
        -d "$DTYPE"
        -o "$pfx"
      )
      # Global rates stay at the profile default unless the caller set them.
      [[ -n "$FP" ]] && LINEAGE_ARGS+=(--fp-rate "$FP")
      [[ -n "$FN" ]] && LINEAGE_ARGS+=(--fn-rate "$FN")
      [[ "$COLLAPSE_HP" == 1 ]] && LINEAGE_ARGS+=(--collapse-hp)
      if ! "$HIMITO" lineage "${LINEAGE_ARGS[@]}" >/dev/null 2>&1; then
        echo "lineage failed at seed=$seed hp_fp=$hp_fp hp_fn=$hp_fn (skipped)" >&2
        continue
      fi
      HAP_ARGS=()
      if [[ -s "${pfx}.cleaned_matrix_all_variants.csv" && -s "$cell/truth/clones.tsv" ]]; then
        HAP_ARGS=(--recon-matrix "${pfx}.cleaned_matrix_all_variants.csv"
                  --truth-clones "$cell/truth/clones.tsv")
      fi
      SCORE_ARGS=()
      [[ -f "$cell/truth/background.txt" ]] && SCORE_ARGS+=(--truth-background "$cell/truth/background.txt")
      [[ -f "$cell/truth/hp_tracts.tsv" ]] && SCORE_ARGS+=(--truth-hp-tracts "$cell/truth/hp_tracts.tsv")
      [[ -f "$cell/truth/hp_lineage.tsv" ]] && SCORE_ARGS+=(--truth-hp-lineage "$cell/truth/hp_lineage.tsv")
      if [[ -n "$cell_background" && "$cell_background" != "none" ]]; then
        SCORE_ARGS+=(--background "$cell_background")
      fi
      tmp="$cell/himito/.hp_score_${hp_fp}_${hp_fn}.tsv"
      rm -f "$tmp"
      python "$HERE/score_lineage.py" \
        --truth-tree "$cell/truth/truth_mutation_tree.tsv" \
        --recon-tree "${pfx}.mutation_tree.tsv" \
        --truth-variants "$cell/truth/truth_variants.txt" \
        --vcf "$VCF" --profile "$PROFILE" --fp "$SCORE_FP" --fn "$SCORE_FN" \
        "${HAP_ARGS[@]+"${HAP_ARGS[@]}"}" \
        "${SCORE_ARGS[@]+"${SCORE_ARGS[@]}"}" \
        --metrics-out "$tmp" >/dev/null
      header=$(head -n 1 "$tmp")
      row=$(tail -n 1 "$tmp")
      rm -f "$tmp"
      if [[ "$first" == 1 ]]; then
        printf 'seed\thp_fp_rate\thp_fn_rate\tcollapse_hp\t%s\n' "$header" > "$SWEEP"
        first=0
      fi
      printf '%s\t%s\t%s\t%s\t%s\n' "$seed" "$hp_fp" "$hp_fn" "$COLLAPSE_TXT" "$row" >> "$SWEEP"
    done
  done
done

[[ $(wc -l < "$SWEEP") -gt 1 ]] || { echo "no successful sweep cells (all lineage runs failed?)" >&2; exit 1; }
python "$HERE/pick_best_hp_rates.py" \
  --metrics "$SWEEP" \
  --summary "$OUTDIR/hp_sweep_summary.tsv" \
  --json "$OUTDIR/best_hp_rates.json" \
  --report "$OUTDIR/hp_sweep_report.md" \
  --objective "$OBJECTIVE" \
  --secondary "$SECONDARY"
echo "=== hp sweep ($SWEEP) ==="
echo "summary: $OUTDIR/hp_sweep_summary.tsv"
echo "best:    $OUTDIR/best_hp_rates.json"
echo "report:  $OUTDIR/hp_sweep_report.md"
