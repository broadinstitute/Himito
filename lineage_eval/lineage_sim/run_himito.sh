#!/usr/bin/env bash
# Align pooled simulated reads to rCRS and run Himito build -> call -> lineage.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
HIMITO="${HIMITO:-$REPO/target/release/Himito}"
REF="${REF:-$REPO/rCRS.fasta}"
OUTDIR="" PROFILE="" SAMPLE="SIM" FP=0.001 FN=0.05 KMER=21
# Call/build-tuning defaults: match `Himito quick-start`'s own defaults, NOT
# tuned-for-simulated-data overrides, so this harness measures the pipeline
# users actually run. Values below are quick-start's CLI defaults (main.rs
# QuickStart: minimal_ac=2, vaf_threshold=0.01, strand_bias_threshold=0.05,
# indel_false_threshold=0.1, build's min_edge_reads hardcoded to 2) or, where
# quick-start leaves a threshold as Option<None>, the same
# call::resolve_thresholds(data_type, None, None, None) fallback quick-start
# gets for a data_type other than ont-r9/ont-r10 (p=0.01, f=0.2, perm=0.7) --
# CALL_DATATYPE here is always "pacbio" or "ont-denoised", so that's the
# branch that applies for every profile this script supports.
MINIMAL_AC=2 VAF=0.01 PVAL=0.01 FREQ_THRESHOLD=0.2 PERM_FREQ_THRESHOLD=0.7
STRAND_BIAS_THRESHOLD=0.05 INDEL_FALSE_THRESHOLD=0.1 MIN_EDGE_READS=2
# HF band forwarded to `Himito lineage` (same band as run_eval.sh / sweep_fpfn.sh / score_lineage.py).
# Unrelated to quick-start, which has no lineage/SCITE step of its own.
MIN_HF=0.01 MAX_HF=0.95
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --fp) FP="$2"; shift 2;;
    --fn) FN="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --himito) HIMITO="$2"; shift 2;;
    --minimal-ac) MINIMAL_AC="$2"; shift 2;;
    --vaf) VAF="$2"; shift 2;;
    --pval) PVAL="$2"; shift 2;;
    --frequency-threshold) FREQ_THRESHOLD="$2"; shift 2;;
    --permutation-frequency-threshold) PERM_FREQ_THRESHOLD="$2"; shift 2;;
    --strand-bias-threshold) STRAND_BIAS_THRESHOLD="$2"; shift 2;;
    --indel-false-threshold) INDEL_FALSE_THRESHOLD="$2"; shift 2;;
    --min-edge-reads) MIN_EDGE_READS="$2"; shift 2;;
    --kmer) KMER="$2"; shift 2;;
    --min-hf) MIN_HF="$2"; shift 2;;
    --max-hf) MAX_HF="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" && -n "$PROFILE" ]] || {
  echo "usage: --outdir DIR --profile {hifi,ont-r10} [--sample S] [--fp F] [--fn F]" >&2
  echo "  [--minimal-ac N] [--vaf V] [--pval P] [--frequency-threshold F]" >&2
  echo "  [--permutation-frequency-threshold F] [--strand-bias-threshold F]" >&2
  echo "  [--indel-false-threshold F] [--min-edge-reads N] [--min-hf F] [--max-hf F]" >&2
  exit 1
}

case "$PROFILE" in
  hifi)    MMPRESET="map-hifi"; DTYPE="pacbio";;
  ont-r10) MMPRESET="lr:hq";  DTYPE="ont-denoised";;
  *) echo "profile must be hifi or ont-r10" >&2; exit 1;;
esac

HDIR="$OUTDIR/himito"; mkdir -p "$HDIR"
FQ="$OUTDIR/reads/reads.fastq.gz"
BAM="$HDIR/aln.sorted.bam"

[[ -s "$FQ" ]] || { echo "missing or empty reads file: $FQ (run simulate_reads.sh first)" >&2; exit 1; }
[[ -f "$REF" ]] || { echo "missing reference FASTA: $REF (set REF or --ref)" >&2; exit 1; }
[[ -x "$HIMITO" ]] || { echo "missing Himito binary: $HIMITO (set HIMITO or --himito; cargo build --release)" >&2; exit 1; }
command -v minimap2 >/dev/null || { echo "minimap2 not in PATH" >&2; exit 1; }
command -v samtools >/dev/null || { echo "samtools not in PATH" >&2; exit 1; }

# Do not swallow minimap2 stderr — a missing REF/broken install otherwise
# looks like "samtools sort: failed to read header from '-'".
minimap2 -ax "$MMPRESET" -t 4 "$REF" "$FQ" \
  | samtools sort -o "$BAM" -
samtools index "$BAM"

# Denoise ONT reads before graph construction (no-op for hifi/pacbio).
# quick-start ties denoise's --vaf to the SAME vaf_threshold used for calling
# (main.rs: `vaf_threshold as f64`); --min-strand 2 and the strand-bias
# defaults (0.01 p / 0.7 homoplasmic-exempt) are quick-start's own hardcoded
# values, and match the Denoise subcommand's CLI defaults, so they're left
# implicit below.
BUILD_BAM="$BAM"
CALL_DATATYPE="$DTYPE"
if [[ "$DTYPE" == ont-* ]]; then
  DENOISED="$HDIR/aln.denoised.bam"
  "$HIMITO" denoise -i "$BAM" -o "$DENOISED" -r "$REF" -d "$DTYPE" \
    --vaf "$VAF" --min-strand 2 --stats "$HDIR/denoise_stats.json"
  samtools index "$DENOISED"
  BUILD_BAM="$DENOISED"
  CALL_DATATYPE="ont-denoised"
fi

# Build anchor graph (input can be a BAM). --min-edge-reads matches
# quick-start's hardcoded edge-read gate (main.rs QuickStart: `build::start(...,
# 2)`); the CLI's own default is 1, which is why this must be passed explicitly.
"$HIMITO" build -i "$BUILD_BAM" -r "$REF" -k "$KMER" -o "$HDIR/sim.gfa" -l 3000 \
  --min-edge-reads "$MIN_EDGE_READS"

# Call variants: -o is the VCF; matrix.csv is derived as <o>.matrix.csv.
"$HIMITO" call -g "$HDIR/sim.gfa" -r "$REF" -s "$SAMPLE" -d "$CALL_DATATYPE" \
  -o "$HDIR/sim.vcf" -k "$KMER" --input-bam "$BUILD_BAM" \
  -m "$MINIMAL_AC" -v "$VAF" -p "$PVAL" -f "$FREQ_THRESHOLD" \
  --permutation-frequency-threshold "$PERM_FREQ_THRESHOLD" \
  --strand-bias-threshold "$STRAND_BIAS_THRESHOLD" \
  --indel-false-threshold "$INDEL_FALSE_THRESHOLD"

# quick-start also runs a NUMT/methylation filter before build and asm/methyl
# steps after call; both are no-ops for this harness's simulated reads (no
# MM/ML tags to classify on, well under the 50000-read downsample cap) and
# produce artifacts (fasta, methylation bed) this harness doesn't score, so
# they're intentionally not reproduced here.

# Lineage: SCITE mutation-tree reconstruction.
"$HIMITO" lineage -m "$HDIR/sim.matrix.csv" -v "$HDIR/sim.vcf" \
  --fp-rate "$FP" --fn-rate "$FN" --min-hf "$MIN_HF" --max-hf "$MAX_HF" \
  -o "$HDIR/sim_lineage"

echo "himito done: $HDIR/sim_lineage.mutation_tree.tsv"
