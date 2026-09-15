#!/usr/bin/env bash
# THE standard Himito-lineage benchmark. Run this to get a number that is
# comparable to the committed baseline and to other people's runs.
#
# ---------------------------------------------------------------------------
# FROZEN PARAMETERS — do not edit without regenerating the baseline.
# ---------------------------------------------------------------------------
# Every value below was chosen from a measured sweep; changing any one of them
# makes the result incomparable to baseline_standard.tsv.
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
#                             committed baseline was produced under, and this
#                             config has ~10x headroom above it (min truth HF
#                             0.296), so it truncates nothing here.
#   seeds 1..30               10 seeds cannot separate a real effect from one
#                             coin-flip edge.
#
# Scored with ad_recall measured against ALL truth pairs. Results produced before
# that change are inflated and must not be compared against this baseline.
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
SEEDS="$(seq -s' ' 1 30)"

BASELINE="$HERE/baseline_standard.tsv"
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
FOCUS = ["var_precision","var_recall","var_f1",
         "ad_precision","ad_recall","ad_f1",
         "hap_precision","hap_recall","hap_f1"]

def load(f):
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    return {k: [float(r[k]) for r in rows if r.get(k) not in (None,"","NA")] for k in FOCUS}, len(rows)

def mean(v): return sum(v)/len(v) if v else float("nan")
def sd(v):
    if len(v) < 2: return 0.0
    m = mean(v); return math.sqrt(sum((x-m)**2 for x in v)/(len(v)-1))

cur, n_cur = load(cur_f)
print(f"seeds: {n_cur}")
if not os.path.exists(base_f):
    print(f"\nno baseline at {base_f} — run with --update-baseline to create it\n")
    print(f"{'metric':<16}{'mean':>8}{'sd':>8}{'min':>8}{'max':>8}")
    for k in FOCUS:
        v = cur[k]
        print(f"{k:<16}{mean(v):>8.4f}{sd(v):>8.4f}{min(v):>8.4f}{max(v):>8.4f}")
    sys.exit(0)

base, n_base = load(base_f)
print(f"baseline: {n_base} seeds ({base_f})\n")
print(f"{'metric':<16}{'baseline':>10}{'current':>10}{'delta':>9}{'2*SE':>8}  verdict")
print("-"*66)
regressed = 0
for k in FOCUS:
    b, c = base[k], cur[k]
    if not b or not c:
        print(f"{k:<16}{'NA':>10}{'NA':>10}"); continue
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
    print(f"{k:<16}{mean(b):>10.4f}{mean(c):>10.4f}{d:>+9.4f}{band:>8.4f}  {verdict}")
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

echo "=== standard benchmark ===" >&2
echo "profile=$PROFILE n_mutations=$NMUT depth=$DEPTH topology=$TOPOLOGY" >&2
echo "internal_keep=$INTERNAL_KEEP -f=$FREQ_THRESHOLD denoise_vaf=$DENOISE_VAF seeds=$(wc -w <<<"$SEEDS" | tr -d ' ')" >&2

"$HERE/sweep_seeds.sh" --outdir "$OUTDIR" --profile "$PROFILE" \
  --n-mutations "$NMUT" --total-depth "$DEPTH" --seeds "$SEEDS" \
  --topology "$TOPOLOGY" --internal-keep "$INTERNAL_KEEP" \
  --frequency-threshold "$FREQ_THRESHOLD" --denoise-vaf "$DENOISE_VAF" \
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
