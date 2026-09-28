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

OUTDIR="" PROFILE="ont-r10" NMUT=12 DEPTH=300 SEED=1 FP="" FN=""
# Himito call -f override forwarded to run_himito.sh. Empty = omit `-f` so
# call::resolve_thresholds applies the data-type default (0.2 here). `-p` is
# never forwarded: run_himito.sh lets the same resolver pick the p-value.
FREQ_THRESHOLD=""
# denoise's keep threshold, forwarded to run_himito.sh (which defaults it to 0.03).
# That default was swept on n=10 cells whose rarest truth variant sits at HF ~0.05.
# It is a *frequency* floor, not a coverage floor, so depth cannot compensate: at
# --n-mutations 15 --sim-min-hf 0.03 the six variants at HF 0.033-0.045 never reach
# raw_matrix.csv and var_recall caps at 0.53 even at --total-depth 3000. Lowering it
# to 0.01 restores them (0.53 -> 0.87). Keep this below --sim-min-hf.
DENOISE_VAF=""
# Shared HF floor for `Himito call -v` and `Himito lineage --min-hf` (forwarded
# via run_himito.sh). score_lineage.py does not re-band: it counts every PASS/.
# call, which is evaluation at this floor once the caller already cut there.
# 0.01 matches the committed benchmarks, whose truth SNVs are *observed* at
# HF 0.06-0.21 and would be mostly filtered away by a 0.1 floor.
# max-hf 0.95 matches the Himito lineage CLI (excludes near-fixed sites).
# Empty --fp/--fn: run_himito.sh passes -d so lineage::resolve_error_rates applies.
MIN_HF=0.01 MAX_HF=0.95
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
# Tree shape, forwarded to simulate_tree.py. 'chain' is the one knob that raises
# mutation-ordering difficulty WITHOUT costing variant recall: a unary path never
# splits mass between siblings, so the minimum clone frequency is higher than a
# random tree of the same --n-mutations. Empty = simulator default (random).
TOPOLOGY=""
# Plan 2026-09-26 T9. Off unless passed; forwarded to run_himito.sh's lineage step.
COLLAPSE_HP=0
# Plan 2026-09-27 T3. Empty = omit the flag (simulate_tree.py default is 0).
HP_LINEAGE_INDELS=""
# Plan 2026-09-27 T3. Forwarded only to Himito lineage, via run_himito.sh.
HP_FP_RATE=""
HP_FN_RATE=""
# Plan 2026-09-27 T4. Empty = omit --hp-unstable. Depth is forwarded only with it.
HP_UNSTABLE=""
HP_MIN_SUBCLONE_FREQ=""
# Plan 2026-09-27 T5. Empty = omit; simulate_tree.py defaults to none.
BACKGROUND=""
# Plan 2026-09-27 T8. One truth, then hifi and ont-denoised. Linux-only: hifi
# needs ccs. Off unless passed; --profile is ignored while this is set.
PAIR=0
# Plan 2026-09-27 T11. "call" stops once sim.matrix.csv and sim.vcf exist.
# Lineage and score_lineage.py are a later sweep step. Reads and truth stay.
STOP_AFTER=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --n-mutations) NMUT="$2"; shift 2;;
    --total-depth) DEPTH="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --fp) FP="$2"; shift 2;;
    --fn) FN="$2"; shift 2;;
    --frequency-threshold) FREQ_THRESHOLD="$2"; shift 2;;
    --denoise-vaf) DENOISE_VAF="$2"; shift 2;;
    --min-hf) MIN_HF="$2"; shift 2;;
    --max-hf) MAX_HF="$2"; shift 2;;
    --sim-min-hf) SIM_MIN_HF="$2"; shift 2;;
    --sim-max-hf) SIM_MAX_HF="$2"; shift 2;;
    --internal-keep) INTERNAL_KEEP="$2"; shift 2;;
    --topology) TOPOLOGY="$2"; shift 2;;
    --ref) export REF="$2"; shift 2;;
    --collapse-hp) COLLAPSE_HP=1; shift;;
    --hp-lineage-indels) HP_LINEAGE_INDELS="$2"; shift 2;;
    --hp-fp-rate) HP_FP_RATE="$2"; shift 2;;
    --hp-fn-rate) HP_FN_RATE="$2"; shift 2;;
    --hp-unstable) HP_UNSTABLE="$2"; shift 2;;
    --hp-min-subclone-freq) HP_MIN_SUBCLONE_FREQ="$2"; shift 2;;
    --background) BACKGROUND="$2"; shift 2;;
    --pair) PAIR=1; shift;;
    --stop-after)
      STOP_AFTER="$2"; shift 2
      [[ "$STOP_AFTER" == "call" ]] || { echo "--stop-after must be call (got $STOP_AFTER)" >&2; exit 1; }
      ;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || {
  echo "usage: --outdir DIR [--profile ...] [--n-mutations N] [--total-depth N] [--seed N] [--fp F] [--fn F] [--frequency-threshold F] [--denoise-vaf F] [--min-hf F] [--max-hf F] [--sim-min-hf F] [--sim-max-hf F] [--internal-keep F] [--topology random|chain|star] [--ref FASTA] [--collapse-hp] [--hp-lineage-indels K] [--hp-fp-rate F] [--hp-fn-rate F] [--hp-unstable SPEC] [--hp-min-subclone-freq F] [--background SPEC] [--pair] [--stop-after call]" >&2
  exit 1
}
# Plan 2026-09-27 T8. Fail before any simulation. An ONT-only run is not a pair.
if [[ "$PAIR" == 1 ]]; then
  command -v ccs >/dev/null 2>&1 || {
    echo "--pair needs HiFi/ccs and is Linux-only for now" >&2
    exit 1
  }
fi
# Resolve in simulate_reads.sh; export a sane default here for visibility.
if [[ -z "${PBSIM_MODEL_DIR:-}" || ! -d "${PBSIM_MODEL_DIR:-}" ]]; then
  export PBSIM_MODEL_DIR="$HERE/pbsim3_models"
fi
[[ -f "$REF" ]] || { echo "missing reference FASTA: $REF" >&2; exit 1; }

# Mirrors lineage::resolve_error_rates (src/lineage.rs) so metrics.tsv records
# the rates Himito lineage actually uses when --fp/--fn are left unset.
resolve_scite_rates() {
  case "$1" in
    hifi) echo "0.005 0.05";;
    ont-denoised) echo "0.0001 0.01";;
    *) echo "0.001 0.05";;  # ont-r10, ont-r9, fallback
  esac
}

mkdir -p "$OUTDIR"
echo "$OUTDIR"
SIM_TREE_ARGS=()
[[ -n "$INTERNAL_KEEP" ]] && SIM_TREE_ARGS+=(--internal-keep "$INTERNAL_KEEP")
[[ -n "$TOPOLOGY" ]] && SIM_TREE_ARGS+=(--topology "$TOPOLOGY")
# Plan 2026-09-27 T3. Empty omits the flag.
[[ -n "$HP_LINEAGE_INDELS" ]] && SIM_TREE_ARGS+=(--hp-lineage-indels "$HP_LINEAGE_INDELS")
# Plan 2026-09-27 T4. Empty SPEC omits --hp-unstable. --total-depth reaches
# simulate_tree.py only together with that flag.
[[ -n "$HP_UNSTABLE" ]] && SIM_TREE_ARGS+=(--hp-unstable "$HP_UNSTABLE")
[[ -n "$HP_MIN_SUBCLONE_FREQ" ]] && SIM_TREE_ARGS+=(--hp-min-subclone-freq "$HP_MIN_SUBCLONE_FREQ")
[[ -n "$HP_UNSTABLE" ]] && SIM_TREE_ARGS+=(--total-depth "$DEPTH")
# Plan 2026-09-27 T5. Empty omits the flag.
[[ -n "$BACKGROUND" ]] && SIM_TREE_ARGS+=(--background "$BACKGROUND")
python "$HERE/simulate_tree.py" --reference "$REF" --n-mutations "$NMUT" --seed "$SEED" --outdir "$OUTDIR" --min-hf "$SIM_MIN_HF" --max-hf "$SIM_MAX_HF" "${SIM_TREE_ARGS[@]+"${SIM_TREE_ARGS[@]}"}"
# Plan 2026-09-27 T6. Pass truth files when this run wrote them, and the
# simulator background when it is not `none`, so scoring uses the same tracts.
SCORE_ARGS=()
[[ -f "$OUTDIR/truth/background.txt" ]] && SCORE_ARGS+=(--truth-background "$OUTDIR/truth/background.txt")
[[ -f "$OUTDIR/truth/hp_tracts.tsv" ]] && SCORE_ARGS+=(--truth-hp-tracts "$OUTDIR/truth/hp_tracts.tsv")
[[ -f "$OUTDIR/truth/hp_lineage.tsv" ]] && SCORE_ARGS+=(--truth-hp-lineage "$OUTDIR/truth/hp_lineage.tsv")
if [[ -n "$BACKGROUND" && "$BACKGROUND" != "none" ]]; then
  SCORE_ARGS+=(--background "$BACKGROUND")
fi

# $1 profile, $2 reconstruction directory name under OUTDIR.
# "himito" is where run_himito.sh writes. Any other name is a rename so a
# second platform (Plan 2026-09-27 T8) cannot clobber the first.
run_one_profile() {
  local profile="$1" dest="$2"
  "$HERE/simulate_reads.sh" --outdir "$OUTDIR" --profile "$profile" --total-depth "$DEPTH" --seed "$SEED"
  # Not `local -a`: macOS bash 3.2 rejects compound assignment on local.
  HIMITO_ARGS=(
    --outdir "$OUTDIR" --profile "$profile" --sample SIM
    --vaf "$MIN_HF" --min-hf "$MIN_HF" --max-hf "$MAX_HF"
  )
  [[ -n "$FP" ]] && HIMITO_ARGS+=(--fp "$FP")
  [[ -n "$FN" ]] && HIMITO_ARGS+=(--fn "$FN")
  [[ -n "$FREQ_THRESHOLD" ]] && HIMITO_ARGS+=(--frequency-threshold "$FREQ_THRESHOLD")
  [[ -n "$DENOISE_VAF" ]] && HIMITO_ARGS+=(--denoise-vaf "$DENOISE_VAF")
  [[ "$COLLAPSE_HP" == 1 ]] && HIMITO_ARGS+=(--collapse-hp)
  # Plan 2026-09-27 T3. run_himito.sh forwards these to Himito lineage only.
  [[ -n "$HP_FP_RATE" ]] && HIMITO_ARGS+=(--hp-fp-rate "$HP_FP_RATE")
  [[ -n "$HP_FN_RATE" ]] && HIMITO_ARGS+=(--hp-fn-rate "$HP_FN_RATE")
  # Plan 2026-09-27 T11. run_himito.sh exits after Himito call.
  [[ -n "$STOP_AFTER" ]] && HIMITO_ARGS+=(--stop-after "$STOP_AFTER")
  "$HERE/run_himito.sh" "${HIMITO_ARGS[@]}"
  if [[ "$dest" != "himito" ]]; then
    rm -rf "$OUTDIR/$dest"
    mv "$OUTDIR/himito" "$OUTDIR/$dest"
  fi
}

# $1 profile, $2 reconstruction directory, $3 metrics TSV path.
score_reconstruction() {
  local profile="$1" hdir="$2" metrics="$3"
  local default_fp default_fn
  read -r default_fp default_fn < <(resolve_scite_rates "$profile")
  python "$HERE/score_lineage.py" \
    --truth-tree "$OUTDIR/truth/truth_mutation_tree.tsv" \
    --recon-tree "$OUTDIR/$hdir/sim_lineage.mutation_tree.tsv" \
    --truth-variants "$OUTDIR/truth/truth_variants.txt" \
    --vcf "$OUTDIR/$hdir/sim.vcf" \
    --recon-matrix "$OUTDIR/$hdir/sim_lineage.cleaned_matrix_all_variants.csv" \
    --truth-clones "$OUTDIR/truth/clones.tsv" \
    --profile "$profile" --fp "${FP:-$default_fp}" --fn "${FN:-$default_fn}" \
    "${SCORE_ARGS[@]+"${SCORE_ARGS[@]}"}" \
    --metrics-out "$metrics"
}

# Plan 2026-09-27 T8. ONT tree scored with the PacBio tree as truth.
# pair_ad_f1 is that run's ad_f1. pair_root_fraction_delta is ONT minus PacBio.
write_pair_metrics() {
  local pb_vars="$OUTDIR/himito_pb/called_variants.txt"
  local pair_score="$OUTDIR/.pair_tree_score.tsv"
  python - "$HERE" "$OUTDIR/himito_pb/sim.vcf" "$pb_vars" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from score_lineage import detected_variants_with_hf_from_vcf
names = sorted(detected_variants_with_hf_from_vcf(sys.argv[2]))
with open(sys.argv[3], "w") as fh:
    for name in names:
        fh.write(name + "\n")
PY
  rm -f "$pair_score"
  python "$HERE/score_lineage.py" \
    --truth-tree "$OUTDIR/himito_pb/sim_lineage.mutation_tree.tsv" \
    --recon-tree "$OUTDIR/himito_ont/sim_lineage.mutation_tree.tsv" \
    --truth-variants "$pb_vars" \
    --vcf "$OUTDIR/himito_ont/sim.vcf" \
    --profile pair \
    "${SCORE_ARGS[@]+"${SCORE_ARGS[@]}"}" \
    --metrics-out "$pair_score"
  python - "$OUTDIR/metrics_ont.tsv" "$OUTDIR/metrics_pb.tsv" "$pair_score" "$OUTDIR/pair_metrics.tsv" <<'PY'
import csv, sys
ont_f, pb_f, pair_f, out_f = sys.argv[1:]

def last(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    if not rows:
        raise SystemExit(f"no data rows in {path}")
    return rows[-1]

ont, pb, pair = last(ont_f), last(pb_f), last(pair_f)
ont_rf, pb_rf = ont["root_fraction"], pb["root_fraction"]
if ont_rf in ("", "NA") or pb_rf in ("", "NA"):
    delta = "NA"
else:
    delta = f"{float(ont_rf) - float(pb_rf):.4f}"
with open(out_f, "w") as fh:
    fh.write("pair_ad_f1\tpair_root_fraction_delta\n")
    fh.write(f"{pair['ad_f1']}\t{delta}\n")
PY
  rm -f "$pair_score"
}

# Plan 2026-09-27 T11. sim.matrix.csv and sim.vcf are enough for sweep_hp_rates.sh.
# Do not score: sim_lineage.* does not exist, and truth/ is left in place.
stop_after_call() {
  echo "stop-after call: skipped lineage and scoring (truth/ and himito/sim.vcf kept)" >&2
  exit 0
}

if [[ "$PAIR" == 1 ]]; then
  # hifi first, then ont-denoised. Each reconstruction is renamed before the
  # next simulate_reads/run_himito, which both write under this same OUTDIR.
  run_one_profile hifi himito_pb
  if [[ "$STOP_AFTER" == "call" ]]; then
    run_one_profile ont-denoised himito_ont
    stop_after_call
  fi
  score_reconstruction hifi himito_pb "$OUTDIR/metrics_pb.tsv"
  run_one_profile ont-denoised himito_ont
  score_reconstruction ont-denoised himito_ont "$OUTDIR/metrics_ont.tsv"
  write_pair_metrics
  echo "=== pair metrics ($OUTDIR/pair_metrics.tsv) ==="
  echo "PacBio metrics: $OUTDIR/metrics_pb.tsv"
  echo "ONT metrics: $OUTDIR/metrics_ont.tsv"
else
  run_one_profile "$PROFILE" himito
  if [[ "$STOP_AFTER" == "call" ]]; then
    stop_after_call
  fi
  score_reconstruction "$PROFILE" himito "$OUTDIR/metrics.tsv"
  echo "=== metrics ($OUTDIR/metrics.tsv) ==="
fi
