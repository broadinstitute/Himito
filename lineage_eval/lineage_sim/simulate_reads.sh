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
[[ -n "$OUTDIR" && -n "$PROFILE" ]] || { echo "usage: --outdir DIR --profile {hifi,ont-r10,ont-denoised} [--total-depth N] [--seed N]" >&2; exit 1; }

# ONT uses pbsim3's qshmm method so reads carry realistic per-base quality
# scores (errhmm emits placeholder Q0 '!' quals, which break any quality-aware
# caller/denoiser). HiFi stays on errhmm (its reads are consensus-polished by
# ccs, and no QSHMM-SEQUEL model ships with pbsim3).
#
# ont-denoised has no distinct read model -- denoising happens later in
# run_himito.sh -- so its reads are simulated exactly like ont-r10.
case "$PROFILE" in
  hifi)                  METHOD="errhmm"; MODEL_NAME="ERRHMM-SEQUEL.model"; PASS=10;;
  ont-r10|ont-denoised)  METHOD="qshmm";  MODEL_NAME="QSHMM-ONT-HQ.model"; PASS=1;;
  *) echo "profile must be hifi or ont-r10 or ont-denoised" >&2; exit 1;;
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
mkdir -p "$OUTDIR/reads"
# Refuse to run two sims against the same --outdir: each does `rm -rf "$WORK"`
# below, so overlapping runs delete/recreate files under each other mid-stream,
# which surfaces much later as a bare "gzip: unexpected end of file".
LOCK="$OUTDIR/reads/.sim.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "another simulate_reads.sh appears to be running for --outdir=$OUTDIR" >&2
  echo "  (lock: $LOCK; if stale, remove with: rmdir '$LOCK')" >&2
  exit 1
fi
trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
rm -rf "$WORK"; mkdir -p "$WORK"
OUT_FQ="$OUTDIR/reads/reads.fastq"
: > "$OUT_FQ"

# Append a gzipped FASTQ into OUT_FQ (rewriting read headers), but only once the
# stream is fully written. pbsim3 returns exit 0 BEFORE its gzip output is flushed
# (confirmed: an immediate `gzip -t` fails, a few ms later it passes), so reading
# eagerly hits a truncated stream and would otherwise abort the whole sim with a
# bare "gzip: unexpected end of file". We poll `gzip -t` until the member is
# complete, and only treat a persistent failure (~10s) as a real corruption — e.g.
# a second run sharing --outdir deleting files mid-write. $4 is an optional log.
append_fq_gz() {  # path clone_id rotidx [producer_log]
  local fq="$1" cid="$2" rid="$3" diag="${4:-}" tries=0 max=100  # 100 * 0.1s = ~10s
  until gzip -t "$fq" 2>/dev/null; do
    tries=$((tries + 1))
    if (( tries >= max )); then
      echo "gzip never became intact after ${max} tries: $fq ($(wc -c <"$fq" 2>/dev/null | tr -d ' ') bytes) for clone $cid rot $rid" >&2
      echo "  the producer never finished flushing, or a second run sharing --outdir=$OUTDIR" >&2
      echo "  deleted/overwrote it mid-read (its 'rm -rf work/' races this)." >&2
      if [[ -n "$diag" && -f "$diag" ]]; then
        echo "---- producer log ($diag) ----" >&2; cat "$diag" >&2 || true; echo "---- end log ----" >&2
      fi
      exit 1
    fi
    sleep 0.1
  done
  gzip -dc "$fq" | awk -v c="$cid" -v r="$rid" \
    'NR%4==1{print "@clone_"c"_rot"r"_"substr($0,2); next}{print}' >> "$OUT_FQ"
}

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
    if ! pbsim --strategy wgs --method "$METHOD" --"$METHOD" "$MODEL" \
          --depth "$depth" --genome "$tmpl" --prefix "$pfx" \
          --pass-num "$PASS" --seed "$seed" >"$pbsim_log" 2>&1; then
      echo "pbsim failed for clone $clone_id rot $rotidx (genome=$tmpl depth=$depth model=$MODEL)" >&2
      echo "---- pbsim log ($pbsim_log) ----" >&2
      cat "$pbsim_log" >&2 || true
      echo "---- end pbsim log ----" >&2
      exit 1
    fi
    # Keep $pbsim_log until the reads are safely appended, so an integrity
    # failure below still has the producer's log for context.
    if [[ "$PROFILE" == "hifi" ]]; then
      # Multi-pass -> HiFi via ccs. pbsim3 emits <pfx>_0001.bam (subreads).
      bam="${pfx}_0001.bam"
      [[ -f "$bam" ]] || { echo "no pbsim bam: $bam" >&2; exit 1; }
      samtools sort -n -o "${pfx}.subreads.bam" "$bam"
      ccs "${pfx}.subreads.bam" "${pfx}.hifi.fastq.gz" --min-passes 3 >/dev/null 2>&1 || {
        echo "ccs failed for clone $clone_id rot $rotidx" >&2; exit 1; }
      append_fq_gz "${pfx}.hifi.fastq.gz" "$clone_id" "$rotidx" "$pbsim_log"
    else
      # Single-pass ONT: pbsim emits <pfx>_0001.fq.gz directly.
      fq="${pfx}_0001.fq.gz"
      [[ -f "$fq" ]] || { echo "no pbsim fq.gz: $fq" >&2; exit 1; }
      append_fq_gz "$fq" "$clone_id" "$rotidx" "$pbsim_log"
    fi
    rm -f "$pbsim_log"
  done
done < <(tail -n +2 "$CLONES" | awk -F'\t' '{print $1 "\t" $4}')

gzip -f "$OUT_FQ"
n=$(gzip -dc "$OUT_FQ.gz" | awk 'END{print NR/4}')
echo "wrote $OUT_FQ.gz ($n reads, profile=$PROFILE, depth=$TOTAL_DEPTH)"
