#!/usr/bin/env bash
# Simulate long reads from every clone genome, handling mtDNA circularity by
# pooling reads from two rotations. Requires an active himito-eval env.
# PBSIM_MODEL_DIR defaults to ./pbsim3_models next to this script.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_MODEL_DIR="$HERE/pbsim3_models"
# Prefer PBSIM_MODEL_DIR if it actually has models; otherwise use local pbsim3_models
# (avoids stale placeholders like /Path/To/Folder/... from old docs).
if [[ -n "${PBSIM_MODEL_DIR:-}" && -d "$PBSIM_MODEL_DIR" ]]; then
  :
else
  export PBSIM_MODEL_DIR="$DEFAULT_MODEL_DIR"
fi

OUTDIR="" PROFILE="" TOTAL_DEPTH=300 SEED=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --total-depth) TOTAL_DEPTH="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" && -n "$PROFILE" ]] || { echo "usage: --outdir DIR --profile {hifi,ont-r10} [--total-depth N] [--seed N]" >&2; exit 1; }

case "$PROFILE" in
  hifi)     MODEL_NAME="ERRHMM-SEQUEL.model"; PASS=10;;
  ont-r10)  MODEL_NAME="ERRHMM-ONT-HQ.model"; PASS=1;;
  *) echo "profile must be hifi or ont-r10" >&2; exit 1;;
esac
MODEL="$PBSIM_MODEL_DIR/$MODEL_NAME"
if [[ ! -f "$MODEL" && -f "$DEFAULT_MODEL_DIR/$MODEL_NAME" ]]; then
  echo "warning: $MODEL missing; falling back to $DEFAULT_MODEL_DIR" >&2
  export PBSIM_MODEL_DIR="$DEFAULT_MODEL_DIR"
  MODEL="$PBSIM_MODEL_DIR/$MODEL_NAME"
fi
[[ -f "$MODEL" ]] || {
  echo "missing model: $MODEL" >&2
  echo "hint: run ./setup_env.sh to fetch models into $DEFAULT_MODEL_DIR" >&2
  exit 1
}
command -v pbsim >/dev/null || {
  echo "pbsim not in PATH (activate the himito-eval env: conda activate himito-eval)" >&2
  exit 1
}

GENOMES="$OUTDIR/truth/clone_genomes.fa"
CLONES="$OUTDIR/truth/clones.tsv"
[[ -s "$GENOMES" && -s "$CLONES" ]] || {
  echo "missing truth outputs: $GENOMES / $CLONES (run simulate_tree.py first)" >&2
  exit 1
}
WORK="$OUTDIR/reads/work"
rm -rf "$WORK"; mkdir -p "$WORK"
OUT_FQ="$OUTDIR/reads/reads.fastq"
: > "$OUT_FQ"

ROT=8284  # floor(16569/2)

# Split multi-FASTA into per-clone files (clone id from header ">clone_<id>").
awk -v dir="$WORK" '/^>/{id=substr($0,8); f=dir"/clone_"id".fa"; print > f; next}{print >> f}' "$GENOMES"


N_CLONES=$(tail -n +2 "$CLONES" | wc -l | tr -d ' ')
N_JOBS=$((N_CLONES * 2))
echo "simulating reads: $N_CLONES clones x 2 rotations = $N_JOBS jobs (profile=$PROFILE, pass-num=$PASS, total-depth=$TOTAL_DEPTH)" >&2
if [[ "$PROFILE" == "hifi" ]]; then
  echo "HiFi is slow: each job runs pbsim (multi-pass CLR) then ccs; progress lines below — not hung." >&2
  command -v ccs >/dev/null || { echo "ccs not in PATH (HiFi requires pbccs on Linux x86)" >&2; exit 1; }
fi

job=0
# Process-substitution (not a pipe) so set -e / exit work inside the loop.
while IFS=$'\t' read -r clone_id freq; do
  fa="$WORK/clone_${clone_id}.fa"
  [[ -f "$fa" ]] || { echo "missing $fa" >&2; exit 1; }
  # depth for this clone, halved across two rotations (integer, min 1)
  depth=$(python3 -c "d=${TOTAL_DEPTH}*float('${freq}')/2; print(max(1,int(round(d))))")

  for rotidx in 0 1; do
    job=$((job + 1))
    echo "[$job/$N_JOBS] clone=$clone_id rot=$rotidx depth=$depth" >&2
    tmpl="$WORK/clone_${clone_id}_rot${rotidx}.fa"
    if [[ "$rotidx" == "0" ]]; then
      cp "$fa" "$tmpl"
    else
      # rotate sequence by $ROT bp; keep header
      python3 - "$fa" "$tmpl" "$ROT" <<'PY'
import sys
inp, out, rot = sys.argv[1], sys.argv[2], int(sys.argv[3])
hdr=None; seq=[]
for ln in open(inp):
    if ln.startswith(">"): hdr=ln.rstrip("\n")
    else: seq.append(ln.strip())
s="".join(seq); s=s[rot:]+s[:rot]
with open(out,"w") as fh:
    fh.write(hdr+"_rot1\n")
    for i in range(0,len(s),60): fh.write(s[i:i+60]+"\n")
PY
    fi

    pfx="$WORK/sim_${clone_id}_rot${rotidx}"
    # pbsim3 wants an integer seed; derive one deterministically.
    seed=$(python3 -c "import hashlib; h=int(hashlib.md5(f'${clone_id}-${rotidx}-${SEED}'.encode()).hexdigest(),16); print(h % 2000000000 + 1)")

    # Keep pbsim log on failure — success output is noisy and previously hid
    # "command not found" / bad-model / bad-genome errors behind a generic message.
    pbsim_log="${pfx}.pbsim.log"
    if ! pbsim --strategy wgs --method errhmm --errhmm "$MODEL" \
          --depth "$depth" --genome "$tmpl" --prefix "$pfx" \
          --pass-num "$PASS" --seed "$seed" >"$pbsim_log" 2>&1; then
      echo "pbsim failed for clone $clone_id rot $rotidx (genome=$tmpl depth=$depth model=$MODEL)" >&2
      echo "---- pbsim log ($pbsim_log) ----" >&2
      cat "$pbsim_log" >&2 || true
      echo "---- end pbsim log ----" >&2
      exit 1
    fi
    rm -f "$pbsim_log"

    if [[ "$PROFILE" == "hifi" ]]; then
      # Multi-pass -> HiFi via ccs. pbsim3 emits <pfx>_0001.bam (subreads).
      bam="${pfx}_0001.bam"
      [[ -f "$bam" ]] || { echo "no pbsim bam: $bam" >&2; exit 1; }
      samtools sort -n -o "${pfx}.subreads.bam" "$bam"
      ccs "${pfx}.subreads.bam" "${pfx}.hifi.fastq.gz" --min-passes 3 >/dev/null 2>&1 || {
        echo "ccs failed for clone $clone_id rot $rotidx" >&2; exit 1; }
      gzip -dc "${pfx}.hifi.fastq.gz" | awk -v c="$clone_id" -v r="$rotidx" \
        'NR%4==1{print "@clone_"c"_rot"r"_"substr($0,2); next}{print}' >> "$OUT_FQ"
    else
      # Single-pass ONT: pbsim emits <pfx>_0001.fq.gz directly.
      fq="${pfx}_0001.fq.gz"
      [[ -f "$fq" ]] || { echo "no pbsim fq.gz: $fq" >&2; exit 1; }
      gzip -dc "$fq" | awk -v c="$clone_id" -v r="$rotidx" \
        'NR%4==1{print "@clone_"c"_rot"r"_"substr($0,2); next}{print}' >> "$OUT_FQ"
    fi
  done
done < <(tail -n +2 "$CLONES" | awk -F'\t' '{print $1 "\t" $4}')

gzip -f "$OUT_FQ"
n=$(gzip -dc "$OUT_FQ.gz" | awk 'END{print NR/4}')
echo "wrote $OUT_FQ.gz ($n reads, profile=$PROFILE, depth=$TOTAL_DEPTH)"
