#!/usr/bin/env bash
# THE frozen homopolymer-length benchmark. Run this to get a number that is
# comparable to baseline_hp.tsv and to other people's HP runs.
#
# Plan 2026-09-27 T8. bench_standard.sh is a different experiment (SNV-only)
# and is left unchanged. Do not score an HP run against baseline_standard.tsv.
#
# ---------------------------------------------------------------------------
# FROZEN PARAMETERS — do not edit without regenerating the baseline.
# ---------------------------------------------------------------------------
# The call/tree settings match bench_standard.sh. The three HP flags are what
# make this run a length-heteroplasmy benchmark. Changing any value below
# makes the result incomparable to baseline_hp.tsv.
#
#   --profile ont-denoised    ont-r10 skips denoise entirely and yields an empty
#                             VCF. hifi is untested on Apple Silicon.
#   --n-mutations 10          NOT a difficulty knob: total clone mass is
#                             conserved, so raising it only costs variant recall.
#                             With chain + keep 0.10 this puts min truth HF at
#                             0.85 * 0.9^10 = 0.296, clear of every gate.
#   --total-depth 1000        3000 removes essentially all headroom (ad_f1 0.996,
#                             1/10 imperfect); 300 wrecks detection (var_recall
#                             0.750) and makes the tree metrics low-power.
#   --topology chain          One unary path: 45 truth ancestral pairs instead of
#                             ~7, and *better* detection than a random tree
#                             (0.990 vs 0.920) because a chain never splits clone
#                             mass between siblings.
#   --internal-keep 0.10      The actual ordering-difficulty knob. Chain shape
#                             alone gives no headroom (ad_f1 0.987); 0.10 puts the
#                             order margin near the ONT dropout rate, where the
#                             call is genuinely hard. See DEFAULT_INTERNAL_KEEP.
#   --frequency-threshold .05 call -f defaults to 0.2, which is ABOVE every truth
#                             variant here and empties the VCF. Must stay below
#                             min truth HF (0.296 in this config).
#   --denoise-vaf 0.03        Pinned, not inherited. run_himito.sh no longer
#                             hardcodes this (it follows --vaf when unset), so the
#                             benchmark has to state it or the result would move
#                             with an unrelated default. 0.03 is the value the
#                             committed SNV baseline was produced under, and this
#                             config has ~10x headroom above it (min truth HF
#                             0.296), so it truncates nothing here.
#   --background hvs1-polyC   T16189C,A16183C. Extends the HVS1 C-run the way the
#                             real PacBio/ONT pairs do.
#   --hp-lineage-indels 2     Two of the ten edges are inherited length changes.
#   --hp-unstable common      D310 / HVS1 / C573 length states that are not tree
#                             edges (lineage-independent heteroplasmy).
#   seeds 1..30               10 seeds cannot separate a real effect from one
#                             coin-flip edge.
#
# Scored with ad_recall measured against ALL truth pairs. Results produced before
# that change are inflated and must not be compared against this baseline.
# baseline_hp.tsv is empty until --update-baseline. A compare with no data rows
# prints "no baseline yet" and exits 0.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

PROFILE=ont-denoised
NMUT=10
DEPTH=1000
TOPOLOGY=chain
INTERNAL_KEEP=0.10
FREQ_THRESHOLD=0.05
DENOISE_VAF=0.03
# Plan 2026-09-27 T8.
BACKGROUND=hvs1-polyC
HP_LINEAGE_INDELS=2
HP_UNSTABLE=common
SEEDS="$(seq -s' ' 1 30)"

BASELINE="$HERE/baseline_hp.tsv"
OUTDIR="" REF_ARG="${REF:-$REPO/rCRS.fasta}" UPDATE=0 COMPARE_ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --ref) REF_ARG="$2"; shift 2;;
    --baseline) BASELINE="$2"; shift 2;;
    --update-baseline) UPDATE=1; shift;;
    --compare) COMPARE_ONLY="$2"; shift 2;;   # score an existing seed_metrics.tsv
    --seeds) SEEDS="$2"; shift 2;;            # smoke runs only; breaks comparability
    -h|--help)
      echo "usage: $0 --outdir DIR [--ref FASTA] [--update-baseline] [--compare FILE]" >&2
      exit 0;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done

compare() {  # $1 = seed_metrics.tsv to score against $BASELINE
  python3 - "$1" "$BASELINE" <<'PY'
import csv, math, os, sys
cur_f, base_f = sys.argv[1], sys.argv[2]
# Standard nine, plus the HP columns. hp_len_l1 may be NA; those cells are
# omitted from the mean (see load).
FOCUS = ["var_precision","var_recall","var_f1",
         "ad_precision","ad_recall","ad_f1",
         "hap_precision","hap_recall","hap_f1",
         "hp_intrusion","hp_backbone_nodes","hp_len_l1",
         "root_fraction","hp_var_f1"]

def load(f):
    with open(f) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = {}
    for k in FOCUS:
        vals = []
        for r in rows:
            raw = r.get(k)
            if raw in (None, "", "NA"):
                continue
            vals.append(float(raw))
        cols[k] = vals
    return cols, len(rows)

def mean(v): return sum(v)/len(v) if v else float("nan")
def sd(v):
    if len(v) < 2: return 0.0
    m = mean(v); return math.sqrt(sum((x-m)**2 for x in v)/(len(v)-1))

def show_current(cur):
    print(f"{'metric':<20}{'mean':>8}{'sd':>8}{'min':>8}{'max':>8}")
    for k in FOCUS:
        v = cur[k]
        if not v:
            print(f"{k:<20}{'NA':>8}{'NA':>8}{'NA':>8}{'NA':>8}")
            continue
        print(f"{k:<20}{mean(v):>8.4f}{sd(v):>8.4f}{min(v):>8.4f}{max(v):>8.4f}")

cur, n_cur = load(cur_f)
print(f"seeds: {n_cur}")
# Missing file and a header-only file are both "no data rows yet".
n_base = 0
if os.path.exists(base_f) and os.path.getsize(base_f) > 0:
    _, n_base = load(base_f)
if n_base == 0:
    print("\nno baseline yet\n")
    show_current(cur)
    sys.exit(0)

base, n_base = load(base_f)
print(f"baseline: {n_base} seeds ({base_f})\n")
print(f"{'metric':<20}{'baseline':>10}{'current':>10}{'delta':>9}{'2*SE':>8}  verdict")
print("-"*70)
regressed = 0
for k in FOCUS:
    b, c = base[k], cur[k]
    if not b or not c:
        print(f"{k:<20}{'NA':>10}{'NA':>10}")
        continue
    d = mean(c) - mean(b)
    # Pooled noise band: a delta inside 2 standard errors is seed scatter.
    se = math.sqrt((sd(b)**2/max(len(b),1)) + (sd(c)**2/max(len(c),1)))
    band = 2*se
    if abs(d) <= band or band == 0 and d == 0:
        verdict = "same"
    elif d > 0:
        verdict = "IMPROVED"
    else:
        verdict = "REGRESSED"; regressed += 1
    print(f"{k:<20}{mean(b):>10.4f}{mean(c):>10.4f}{d:>+9.4f}{band:>8.4f}  {verdict}")
print()
print(f"{regressed} metric(s) regressed beyond seed noise." if regressed
      else "No metric regressed beyond seed noise.")
PY
}

if [[ -n "$COMPARE_ONLY" ]]; then
  compare "$COMPARE_ONLY"; exit $?
fi

[[ -n "$OUTDIR" ]] || { echo "usage: $0 --outdir DIR [--ref FASTA] [--update-baseline]" >&2; exit 1; }
[[ -f "$REF_ARG" ]] || { echo "missing reference FASTA: $REF_ARG (pass --ref)" >&2; exit 1; }

echo "=== HP benchmark ===" >&2
echo "profile=$PROFILE n_mutations=$NMUT depth=$DEPTH topology=$TOPOLOGY" >&2
echo "internal_keep=$INTERNAL_KEEP -f=$FREQ_THRESHOLD denoise_vaf=$DENOISE_VAF seeds=$(wc -w <<<"$SEEDS" | tr -d ' ')" >&2
echo "background=$BACKGROUND hp_lineage_indels=$HP_LINEAGE_INDELS hp_unstable=$HP_UNSTABLE" >&2

# Plan 2026-09-27 T8. sweep_seeds.sh already accepts these flags.
"$HERE/sweep_seeds.sh" --outdir "$OUTDIR" --profile "$PROFILE" \
  --n-mutations "$NMUT" --total-depth "$DEPTH" --seeds "$SEEDS" \
  --topology "$TOPOLOGY" --internal-keep "$INTERNAL_KEEP" \
  --frequency-threshold "$FREQ_THRESHOLD" --denoise-vaf "$DENOISE_VAF" \
  --background "$BACKGROUND" --hp-lineage-indels "$HP_LINEAGE_INDELS" \
  --hp-unstable "$HP_UNSTABLE" \
  --ref "$REF_ARG" || true

RESULT="$OUTDIR/seed_metrics.tsv"
[[ -s "$RESULT" ]] || { echo "no seed completed; nothing to score" >&2; exit 1; }

echo
echo "=== vs baseline ==="
compare "$RESULT"

if [[ "$UPDATE" == "1" ]]; then
  cp "$RESULT" "$BASELINE"
  echo
  echo "baseline updated: $BASELINE" >&2
fi
