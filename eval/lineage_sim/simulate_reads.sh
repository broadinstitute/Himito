#!/usr/bin/env bash
# Simulate long reads from every clone genome, handling mtDNA circularity by
# pooling reads from two rotations. Requires an active himito-eval env and
# PBSIM_MODEL_DIR set.
set -euo pipefail

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
: "${PBSIM_MODEL_DIR:?set PBSIM_MODEL_DIR (see setup_env.sh output)}"

case "$PROFILE" in
  hifi)     MODEL="$PBSIM_MODEL_DIR/ERRHMM-SEQUEL.model"; PASS=10;;
  ont-r10)  MODEL="$PBSIM_MODEL_DIR/ERRHMM-ONT-HQ.model"; PASS=1;;
  *) echo "profile must be hifi or ont-r10" >&2; exit 1;;
esac
[[ -f "$MODEL" ]] || { echo "missing model: $MODEL" >&2; exit 1; }

GENOMES="$OUTDIR/truth/clone_genomes.fa"
CLONES="$OUTDIR/truth/clones.tsv"
WORK="$OUTDIR/reads/work"
rm -rf "$WORK"; mkdir -p "$WORK"
OUT_FQ="$OUTDIR/reads/reads.fastq"
: > "$OUT_FQ"

ROT=8284  # floor(16569/2)

# Split multi-FASTA into per-clone files (clone id from header ">clone_<id>").
awk -v dir="$WORK" '/^>/{id=substr($0,8); f=dir"/clone_"id".fa"; print > f; next}{print >> f}' "$GENOMES"

# Read clone frequencies (skip header).
# Use awk to emit fields so empty variant_path (ref clone) doesn't collapse under IFS.
tail -n +2 "$CLONES" | awk -F'\t' '{print $1 "\t" $4}' | while IFS=$'\t' read -r clone_id freq; do
  fa="$WORK/clone_${clone_id}.fa"
  [[ -f "$fa" ]] || { echo "missing $fa" >&2; exit 1; }
  # depth for this clone, halved across two rotations (integer, min 1)
  depth=$(python3 -c "import math,sys; d=${TOTAL_DEPTH}*float('${freq}')/2; print(max(1,int(round(d))))")

  for rotidx in 0 1; do
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
    seed=$(python3 -c "print((hash(('${clone_id}',${rotidx},${SEED})) % 2000000000) + 1)")

    pbsim --strategy wgs --method errhmm --errhmm "$MODEL" \
          --depth "$depth" --genome "$tmpl" --prefix "$pfx" \
          --pass-num "$PASS" --seed "$seed" >/dev/null 2>&1 || {
            echo "pbsim failed for clone $clone_id rot $rotidx" >&2; exit 1; }

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
done

gzip -f "$OUT_FQ"
n=$(gzip -dc "$OUT_FQ.gz" | awk 'END{print NR/4}')
echo "wrote $OUT_FQ.gz ($n reads, profile=$PROFILE, depth=$TOTAL_DEPTH)"
