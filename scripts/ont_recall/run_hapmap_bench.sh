#!/usr/bin/env bash
# HapMap ONT mito benchmark for Himito quick-start.
#
# Runs quick-start, then scores the VCF this run wrote (<prefix>.vcf).
#
# usage (from ont_benchmark/):
#   ../scripts/ont_recall/run_hapmap_bench.sh MODE PREFIX [extra quick-start flags...]
#   ../scripts/ont_recall/run_hapmap_bench.sh ont-r10 r10_default
#   ../scripts/ont_recall/run_hapmap_bench.sh ont-r10 r10_perm10k --permutation-rounds 10000
# env overrides: HIMITO RTG REF SDF TRUTH BAM OUT
set -euo pipefail
[[ $# -ge 2 ]] || { sed -n 2,10p "$0"; exit 1; }
MODE=$1; PREFIX=$2; shift 2
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HIMITO=${HIMITO:-$HERE/../../target/release/Himito}
RTG=${RTG:-/Users/suhang/Analysis/rtg-tools-3.13/rtg}
REF=${REF:-../rCRS.fasta}
SDF=${SDF:-../rtg_ref}
TRUTH=${TRUTH:-SMHTHAPMAP6_GRCh38_v1.0.0_mito_benchmark.vcf.gz}
BAM=${BAM:-HapMap_ont.bam}
OUT=${OUT:-output}

mtime() { stat -f %m "$1" 2>/dev/null || stat -c %Y "$1"; }

t0=$(date +%s)
"$HIMITO" quick-start -i "$BAM" -c chrM -r "$REF" -s HG00097_ont -d "$MODE" -o "$PREFIX" \
  "$@" > "$PREFIX.quickstart.log" 2>&1
t1=$(date +%s)
echo "quick-start: $((t1 - t0))s total; strand-bias + permutation step: $(( $(mtime "$PREFIX.vcf") - $(mtime "$PREFIX.raw_matrix.csv") ))s"

bcftools view -i 'FILTER="PASS"' -Oz -o "$PREFIX.pass.vcf.gz" "$PREFIX.vcf"
bcftools index -t -f "$PREFIX.pass.vcf.gz"
rm -rf "$OUT/$PREFIX.eval"   # vcfeval refuses to write into an existing directory
"$RTG" vcfeval -b "$TRUTH" -c "$PREFIX.pass.vcf.gz" -o "$OUT/$PREFIX.eval" -t "$SDF" \
  --squash-ploidy --sample ALT,ALT > /dev/null
cat "$OUT/$PREFIX.eval/summary.txt"
python3 "$HERE/score_vs_truth.py" --truth "$TRUTH" "$PREFIX.vcf"
