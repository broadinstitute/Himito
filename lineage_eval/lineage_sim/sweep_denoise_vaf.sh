#!/usr/bin/env bash
# Sweep denoise's keep threshold (--denoise-vaf) over the mut10 cells and score each.
#
# `call` runs with --permutation-rounds 1: at PVAL=1 the permutation test can never
# exclude anything (exclusion needs a BH q-value > 1), so the round count cannot
# change the VCF. Verified byte-identical against the 1000-round default, and it
# takes the per-cell cost from ~230s to ~2s, which is what makes this sweep feasible.
# No `set -e`: a single cell that panics or produces no callable variant must not
# abort the whole sweep. Failures are reported per cell instead.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
H="${HIMITO:-$REPO/target/release/Himito}"
REF="${REF:-$REPO/rCRS.fasta}"
SW="$HERE/tmp/simulation_ont_r10"
OUT="${OUT:-/tmp/vafsweep}"
VAFS="${VAFS:-0.01 0.02 0.03 0.04}"
# depth-100 cells yield too few callable variants for SCITE in every arm tested, so
# they only add runtime and a constant zero to every threshold's mean.
MIN_DEPTH="${MIN_DEPTH:-300}"

for v in $VAFS; do
  for src in "$SW"/*mut10*; do
    c=$(basename "$src")
    [ -s "$src/reads/reads.fastq.gz" ] || continue
    [ "${c##*depth}" -ge "$MIN_DEPTH" ] || continue
    D="$OUT/v$v/$c"; HD="$D/himito"
    rm -rf "$D"; mkdir -p "$HD"
    cp -R "$src/truth" "$D/" 2>/dev/null || continue

    # Reuse the alignment; only denoise onward depends on the swept threshold.
    BAM="$src/himito/aln.sorted.bam"
    [ -s "$BAM" ] || { echo "skip $c (no bam)"; continue; }

    log="$D/run.log"
    # --indels matches run_himito.sh: without it, uncorrected indel artifacts crowd
    # out the SNVs entirely and every threshold scores the same (and badly).
    "$H" denoise -i "$BAM" -o "$HD/dn.bam" -r "$REF" -d ont-denoised --indels \
      --vaf "$v" --min-strand 2 --stats "$HD/denoise_stats.json" >>"$log" 2>&1 \
      || { echo "DENOISEFAIL v$v $c"; continue; }
    samtools index "$HD/dn.bam" >>"$log" 2>&1
    "$H" build -i "$HD/dn.bam" -r "$REF" -k 21 -o "$HD/sim.gfa" -l 3000 \
      --min-edge-reads 2 >>"$log" 2>&1 || { echo "BUILDFAIL v$v $c"; continue; }
    "$H" call -g "$HD/sim.gfa" -r "$REF" -s SIM -d ont-denoised -o "$HD/sim.vcf" -k 21 \
      --input-bam "$HD/dn.bam" -m 2 -v 0.01 -p 1 -f 0.2 \
      --permutation-frequency-threshold 0.7 --strand-bias-threshold 0.05 \
      --indel-false-threshold 0.1 --permutation-rounds 1 >>"$log" 2>&1 \
      || { echo "CALLFAIL v$v $c"; continue; }
    if "$H" lineage -m "$HD/sim.matrix.csv" -v "$HD/sim.vcf" \
        --fp-rate 0.001 --fn-rate 0.05 --min-hf 0.01 --max-hf 0.95 \
        -o "$HD/sim_lineage" >/dev/null 2>&1; then
      python "$HERE/score_lineage.py" --truth-tree "$D/truth/truth_mutation_tree.tsv" \
        --recon-tree "$HD/sim_lineage.mutation_tree.tsv" \
        --truth-variants "$D/truth/truth_variants.txt" --vcf "$HD/sim.vcf" \
        --profile ont-r10 --fp 0.001 --fn 0.05 --metrics-out "$D/metrics.tsv" >/dev/null 2>&1 \
        || echo "SCOREFAIL v$v $c"
    else
      echo "TREEFAIL v$v $c"
    fi
  done
  echo "done vaf=$v"
done
echo "SWEEPDONE"
