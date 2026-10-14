#!/bin/bash
# HapMap ONT mito benchmark. Scores THE VCF THIS RUN WROTE.
#
# usage (from ont_benchmark/):
#   bash benchmarking.sh ont-r10              # → PREFIX=HG00097_r10
#   bash benchmarking.sh ont-denoised         # → PREFIX=HG00097_ont
#   bash benchmarking.sh ont-r10 my_prefix
#   bash benchmarking.sh ont-r10 my_prefix --extra quick-start flags...
set -euo pipefail

HIMITO_PATH="${HIMITO_PATH:-/Users/suhang/Analysis/Himito/target/release/Himito}"
VCFEVAL_PATH="${VCFEVAL_PATH:-/Users/suhang/Analysis/rtg-tools-3.13/rtg}"
TRUTH_VCF_PATH="${TRUTH_VCF_PATH:-/Users/suhang/Analysis/Himito/ont_benchmark/SMHTHAPMAP6_GRCh38_v1.0.0_mito_benchmark.vcf.gz}"
OUTPUT_DIR="${OUTPUT_DIR:-/Users/suhang/Analysis/Himito/pb_benchmark/output}"
HERE="$(cd "$(dirname "$0")" && pwd)"

MODE=${1:-pacbio}
shift || true
case "$MODE" in
  ont-r10)      DEFAULT_PREFIX=HG00097_r10 ;;
  ont-denoised) DEFAULT_PREFIX=HG00097_ont ;;
  ont-r9|pacbio) DEFAULT_PREFIX=HG00097_${MODE} ;;
  *) echo "unknown mode: $MODE (expected ont-r10|ont-denoised|…)" >&2; exit 1 ;;
esac
PREFIX=${1:-$DEFAULT_PREFIX}
if [[ $# -gt 0 ]]; then shift; fi

"$HIMITO_PATH" quick-start \
  -i HapMap_pb_merge.bam -c chrM -r ../rCRS.fasta -s HG00097_pb \
  -d "$MODE" -o "$PREFIX" "$@"

# PASS-only call set, named so it cannot collide with the full VCF
bcftools view -i 'FILTER="PASS"' -Oz -o "${PREFIX}.pass.vcf.gz" "${PREFIX}.vcf"
bcftools index -t -f "${PREFIX}.pass.vcf.gz"

rm -rf "${OUTPUT_DIR}/${PREFIX}.eval"
"$VCFEVAL_PATH" vcfeval \
  -b "$TRUTH_VCF_PATH" -c "${PREFIX}.pass.vcf.gz" \
  -o "${OUTPUT_DIR}/${PREFIX}.eval" -t ../rtg_ref \
  --squash-ploidy --sample ALT,ALT
cat "${OUTPUT_DIR}/${PREFIX}.eval/summary.txt"

if [[ -x "$HERE/../scripts/ont_recall/score_vs_truth.py" ]]; then
  python3 "$HERE/../scripts/ont_recall/score_vs_truth.py" \
    --truth "$TRUTH_VCF_PATH" "${PREFIX}.vcf"
fi
