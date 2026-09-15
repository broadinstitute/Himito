#!/usr/bin/env bash
# Sweep `Himito denoise --vaf` across replicate seeds and report every focus
# metric per threshold.
#
# Replaces an earlier version that hardcoded paths under ./tmp/simulation_ont_r10,
# selected cells by globbing `*mut10*`, and called `call` with ont-r10-era flags
# (-f 0.2 -p 1 --permutation-rounds 1). None of that runs against the current
# harness, and -f 0.2 empties the VCF outright. This drives `sweep_seeds.sh`
# instead, so every cell goes through exactly the same path as a normal eval.
#
# --vaf is a FREQUENCY floor applied before the graph is built. Its effect
# depends entirely on how much headroom the simulated truth has above it:
#
#   * margin >> 1 (e.g. the standard config, min truth HF 0.296, ~10x the 0.03
#     default): nothing is truncated and the threshold only trades variant
#     precision against noise.
#   * margin ~ 1 (e.g. --n-mutations 15 --sim-min-hf 0.03, min truth HF 0.033):
#     the gate deletes truth variants before they reach raw_matrix.csv and
#     var_recall is capped no matter how much depth is added.
#
# Sweep the regime you actually care about; the optimum is not shared between
# them, which is why a single hardcoded default was the wrong design.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

OUTDIR=""
VAFS="0.005 0.01 0.03 0.05"
SEEDS="1 2 3 4 5 6 7 8 9 10"
# Defaults mirror bench_standard.sh so the sweep is comparable to the baseline.
PROFILE=ont-denoised NMUT=10 DEPTH=1000 TOPOLOGY=chain INTERNAL_KEEP=0.10
FREQ_THRESHOLD=0.05 SIM_MIN_HF=""
REF_ARG="${REF:-$REPO/rCRS.fasta}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --vafs) VAFS="$2"; shift 2;;
    --seeds) SEEDS="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --n-mutations) NMUT="$2"; shift 2;;
    --total-depth) DEPTH="$2"; shift 2;;
    --topology) TOPOLOGY="$2"; shift 2;;
    --internal-keep) INTERNAL_KEEP="$2"; shift 2;;
    --frequency-threshold) FREQ_THRESHOLD="$2"; shift 2;;
    --sim-min-hf) SIM_MIN_HF="$2"; shift 2;;
    --ref) REF_ARG="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" ]] || {
  echo "usage: $0 --outdir DIR [--vafs \"0.005 0.01 0.03\"] [--seeds \"1 2 3\"] [--n-mutations N] [--topology T] [--internal-keep F] [--sim-min-hf F] [--ref FASTA]" >&2
  exit 1
}

# Must exist before the per-cell log redirect below; sweep_seeds.sh only creates
# its own $OUTDIR/v<vaf> subdirectory, not this parent.
mkdir -p "$OUTDIR"

echo "denoise --vaf sweep: [$VAFS] x $(wc -w <<<"$SEEDS" | tr -d ' ') seeds" >&2
echo "  config: $PROFILE n=$NMUT depth=$DEPTH topology=$TOPOLOGY keep=$INTERNAL_KEEP -f=$FREQ_THRESHOLD" >&2

for v in $VAFS; do
  echo "=== denoise --vaf $v ===" >&2
  ARGS=(--outdir "$OUTDIR/v$v" --profile "$PROFILE" --n-mutations "$NMUT"
        --total-depth "$DEPTH" --seeds "$SEEDS" --topology "$TOPOLOGY"
        --internal-keep "$INTERNAL_KEEP" --frequency-threshold "$FREQ_THRESHOLD"
        --denoise-vaf "$v" --ref "$REF_ARG")
  [[ -n "$SIM_MIN_HF" ]] && ARGS+=(--sim-min-hf "$SIM_MIN_HF")
  "$HERE/sweep_seeds.sh" "${ARGS[@]}" >"$OUTDIR/v$v.log" 2>&1 \
    || echo "  cell --vaf $v failed (see $OUTDIR/v$v.log)" >&2
done

echo
echo "=== summary ==="
python3 - "$OUTDIR" $VAFS <<'PY'
import csv, os, statistics as st, sys
out, vafs = sys.argv[1], sys.argv[2:]
FOCUS = ["var_precision","var_recall","var_f1","ad_f1","hap_f1"]
print(f"{'--vaf':>8}{'seeds':>6}" + "".join(f"{k:>15}" for k in FOCUS))
print("-"*(14+15*len(FOCUS)))
rows=[]
for v in vafs:
    f=os.path.join(out,f"v{v}","seed_metrics.tsv")
    if not os.path.exists(f):
        print(f"{v:>8}{'-':>6}  (no result)"); continue
    R=list(csv.DictReader(open(f),delimiter="\t"))
    if not R: continue
    m={k:st.mean([float(r[k]) for r in R if r.get(k) not in (None,"","NA")]) for k in FOCUS}
    rows.append((v,m))
    print(f"{v:>8}{len(R):>6}" + "".join(f"{m[k]:>15.4f}" for k in FOCUS))
if rows:
    best_ad=max(rows,key=lambda r:r[1]["ad_f1"])
    best_vp=max(rows,key=lambda r:r[1]["var_precision"])
    print(f"\nbest ad_f1:         --vaf {best_ad[0]} ({best_ad[1]['ad_f1']:.4f})")
    print(f"best var_precision: --vaf {best_vp[0]} ({best_vp[1]['var_precision']:.4f})")
    print("\nCheck var_recall first: if it falls as --vaf rises, the gate is")
    print("truncating truth variants and any precision gain is bought with recall.")
PY
