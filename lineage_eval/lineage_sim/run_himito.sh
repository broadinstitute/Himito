#!/usr/bin/env bash
# Align pooled simulated reads to rCRS and run Himito build -> call -> lineage.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
HIMITO="${HIMITO:-$REPO/target/release/Himito}"
REF="${REF:-$REPO/rCRS.fasta}"
OUTDIR="" PROFILE="" SAMPLE="SIM" FP="" FN="" KMER=21

# Call/build-tuning defaults: match `Himito quick-start`'s own defaults, NOT
# tuned-for-simulated-data overrides, so this harness measures the pipeline
# users actually run. Values below are quick-start's CLI defaults (main.rs
# QuickStart: minimal_ac=2, vaf_threshold=0.01, strand_bias_threshold=0.05,
# indel_false_threshold=0.1, build's min_edge_reads hardcoded to 2) or, where
# quick-start leaves a threshold as Option<None>, the same
# call::resolve_thresholds(data_type, None, None, None) fallback (p=0.01,
# f=0.2, perm=0.7 for every profile this script supports). `--p-value-threshold`
# is deliberately omitted so that data-type default runs; pass
# --frequency-threshold to override `-f` without touching `-p`.

# Call -v tracks lineage --min-hf unless --vaf is passed explicitly, so the
# VCF that score_lineage.py reads is already cut at the same floor lineage uses.
# Empty here; filled with MIN_HF after argv parsing.
MINIMAL_AC=2 VAF="" FREQ_THRESHOLD="" PERM_FREQ_THRESHOLD=0.7
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
# MIN_EDGE_READS deliberately diverges from quick-start's hardcoded 2: the gate is
# now inclusive (>= N reads), and 1 CIGARs every read-supported edge. At 2, ~98% of
# edges on this ONT graph go un-CIGARed and the reads reaching a bubble through them
# are recorded as missing rather than ref, which is what starved the lineage matrix.
MIN_EDGE_READS=2

# HF band for `Himito lineage` and, unless --vaf is set, for `Himito call -v`.
# score_lineage.py still counts every PASS/. call (no second HF band): with call
# already cut at this floor, that is evaluation at the same min-hf.
# Empty --fp/--fn: Himito lineage -d uses resolve_error_rates (src/lineage.rs).

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
  echo "  [--minimal-ac N] [--vaf V] [--denoise-vaf V] [--frequency-threshold F]" >&2
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

# After profile/min-hf are known: call -v follows --min-hf unless --vaf was set.
VAF="${VAF:-$MIN_HF}"

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
# (main.rs: `vaf_threshold as f64`). The strand gates are no longer flags at all:
# min-strand (2) and the strand-bias p (0.01) are the DENOISE_MIN_STRAND /
# DENOISE_STRAND_BIAS_P constants in main.rs, identical on both paths. The
# near-homoplasmic exemption does NOT match: quick-start uses 0.7, the standalone
# denoise subcommand uses 0.95 (DENOISE_HOMOPLASMIC_VAF). This script runs
# `denoise` directly, so it benchmarks at 0.95.
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
  DENOISE_ARGS=(--vaf "${DENOISE_VAF:-$VAF}" --stats "$HDIR/denoise_stats.json")
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
# `-p` / `--p-value-threshold` is omitted so call::resolve_thresholds applies
# the data-type default (0.01 for every profile this script supports). `-f` is
# only forwarded when --frequency-threshold was set; otherwise that default
# (0.2 here) applies too.
CALL_ARGS=(
  -g "$HDIR/sim.gfa" -r "$REF" -s "$SAMPLE" -d "$CALL_DATATYPE"
  -o "$HDIR/sim.vcf" -k "$KMER" --input-bam "$BUILD_BAM"
  -m "$MINIMAL_AC" -v "$VAF"
  --permutation-frequency-threshold "$PERM_FREQ_THRESHOLD"
  --strand-bias-threshold "$STRAND_BIAS_THRESHOLD"
  --indel-false-threshold "$INDEL_FALSE_THRESHOLD"
)
[[ -n "$FREQ_THRESHOLD" ]] && CALL_ARGS+=(-f "$FREQ_THRESHOLD")
"$HIMITO" call "${CALL_ARGS[@]}"

# quick-start also runs a NUMT/methylation filter before build and asm/methyl
# steps after call; both are no-ops for this harness's simulated reads (no
# MM/ML tags to classify on, well under the 50000-read downsample cap) and
# produce artifacts (fasta, methylation bed) this harness doesn't score, so
# they're intentionally not reproduced here.

# Lineage: SCITE mutation-tree reconstruction. -d selects resolve_error_rates
# presets (pacbio 0.005/0.05, ont-r10 0.001/0.05, ont-denoised 0.0001/0.01)
# unless --fp/--fn were passed.
LINEAGE_ARGS=(
  -m "$HDIR/sim.matrix.csv" -v "$HDIR/sim.vcf"
  --min-hf "$MIN_HF" --max-hf "$MAX_HF"
  -d "$DTYPE"
  -o "$HDIR/sim_lineage"
)
[[ -n "$FP" ]] && LINEAGE_ARGS+=(--fp-rate "$FP")
[[ -n "$FN" ]] && LINEAGE_ARGS+=(--fn-rate "$FN")
"$HIMITO" lineage "${LINEAGE_ARGS[@]}"

echo "himito done: $HDIR/sim_lineage.mutation_tree.tsv"
