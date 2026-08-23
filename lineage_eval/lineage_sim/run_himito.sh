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

MINIMAL_AC=2 VAF=0.01 PVAL=1 FREQ_THRESHOLD=0.2 PERM_FREQ_THRESHOLD=0.7
STRAND_BIAS_THRESHOLD=0.05 INDEL_FALSE_THRESHOLD=0.1
# denoise's keep threshold, decoupled from the caller's -v. quick-start ties the two
# together, but they answer different questions: -v is a raw HF cut on a called
# variant, while this one is compared against the site model's EM estimate, which
# discounts alt observations explainable as quality-weighted error. That makes it a
# sharper filter on noise sites whose raw HF overlaps real low-frequency
# heteroplasmies.
#
# 0.03 is the sweep optimum over the n=10 depth>=300 cells (see sweep_denoise_vaf.sh
# and DENOISE_KEEP_VAF in src/main.rs): variant precision 0.837 -> 0.967 with recall
# held at 1.000 and ad_f1 within 0.003 of best. Empty = follow --vaf.
DENOISE_VAF=0.03
# PVAL=1 (permutation test effectively disabled) is required for this harness to
# call anything at all. At PVAL=0.1 the simulated ont-r10 n=10 cells produce an
# EMPTY VCF -- not a single one of the ten truth SNVs survives, even though the
# reads plainly carry them at 4-17% HF and `-p 1` calls all ten as PASS. That
# makes every downstream lineage metric zero, and it is why the committed
# benchmarks under lineage_eval/benchmark_results (which have var_f1 ~ 1.0)
# cannot be reproduced with a 0.1 default. Lower this only alongside a check
# that truth SNVs are still called.
# MIN_EDGE_READS deliberately diverges from quick-start's hardcoded 2: the gate is
# now inclusive (>= N reads), and 1 CIGARs every read-supported edge. At 2, ~98% of
# edges on this ONT graph go un-CIGARed and the reads reaching a bubble through them
# are recorded as missing rather than ref, which is what starved the lineage matrix.
MIN_EDGE_READS=2

# HF band forwarded to `Himito lineage` (same band as run_eval.sh / sweep_fpfn.sh).
# Unrelated to quick-start, which has no lineage/SCITE step of its own, and no longer
# related to scoring: score_lineage.py counts every PASS/. call regardless of HF.

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
    --denoise-vaf) DENOISE_VAF="$2"; shift 2;;
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
  echo "  [--minimal-ac N] [--vaf V] [--denoise-vaf V] [--pval P] [--frequency-threshold F]" >&2
  echo "  [--permutation-frequency-threshold F] [--strand-bias-threshold F]" >&2
  echo "  [--indel-false-threshold F] [--min-edge-reads N] [--min-hf F] [--max-hf F]" >&2
  exit 1
}

case "$PROFILE" in
  hifi)    MMPRESET="map-hifi"; DTYPE="pacbio";;
  ont-r10) MMPRESET="lr:hq";  DTYPE="ont-r10";;
  ont-denoised) MMPRESET="lr:hq";  DTYPE="ont-denoised";;
  *) echo "profile must be hifi or ont-r10 or ont-denoised" >&2; exit 1;;
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
#
# Indel denoising is ON by default. It used to be documented as off with a
# DENOISE_INDELS=1 opt-in, but `--indels` was ALSO appended unconditionally to the
# command below, so the switch was dead and indel denoising ran every time -- which
# is the configuration every committed benchmark here was produced under. Turning it
# genuinely off costs almost everything: on seed1_mut10_depth300 the VCF goes from 19
# PASS SNVs / 1 indel to 1 PASS SNV / 10 indels, because uncorrected indel artifacts
# crowd the graph and matrix. Set DENOISE_INDELS=0 to run the off-vs-on comparison.
DENOISE_INDELS="${DENOISE_INDELS:-1}"

BUILD_BAM="$BAM"
CALL_DATATYPE="$DTYPE"
if [[ "$DTYPE" == ont-denoised ]]; then
  DENOISED="$HDIR/aln.denoised.bam"
  DENOISE_ARGS=(--vaf "${DENOISE_VAF:-$VAF}" --min-strand 2 --stats "$HDIR/denoise_stats.json")
  if [[ "$DENOISE_INDELS" == "1" ]]; then
    DENOISE_ARGS+=(--indels)
  fi
  "$HIMITO" denoise -i "$BAM" -o "$DENOISED" -r "$REF" -d "$DTYPE" "${DENOISE_ARGS[@]}"
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
  --fp-rate "$FP" --fn-rate "$FN" --min-hf 0.01 --max-hf 0.95 \
  -o "$HDIR/sim_lineage"

echo "himito done: $HDIR/sim_lineage.mutation_tree.tsv"
