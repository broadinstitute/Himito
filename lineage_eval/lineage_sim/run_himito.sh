#!/usr/bin/env bash
# Align pooled simulated reads to rCRS and run Himito build -> call -> lineage.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
HIMITO="${HIMITO:-$REPO/target/release/Himito}"
REF="${REF:-$REPO/rCRS.fasta}"
OUTDIR="" PROFILE="" SAMPLE="SIM" FP="" FN="" KMER=21
HIMITO_THREADS="${HIMITO_THREADS:-4}"

# Reuse the caller and denoiser presets unless the benchmark explicitly overrides
# them. The standalone caller has a different strand-bias default from quick-start.
MINIMAL_AC=2 VAF="" FREQ_THRESHOLD="" PERM_FREQ_THRESHOLD=""
STRAND_BIAS_THRESHOLD=0.05 INDEL_FALSE_THRESHOLD=""
DENOISE_VAF=""
MIN_EDGE_READS=2

# The lineage HF band is independent of the caller and denoiser floors.
MIN_HF=0.01 MAX_HF=0.95
# Himito lineage is SNV-only by default. --include-non-snvs admits indels (needed
# for the HP benchmark); --collapse-hp requires it and turns it on.
INCLUDE_NON_SNVS=0
COLLAPSE_HP=0
# Forwarded only to `Himito lineage`, not `Himito call`.
HP_FP_RATE=""
HP_FN_RATE=""
# "call" skips Himito lineage after sim.matrix.csv and sim.vcf exist.
STOP_AFTER=""
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
    --include-non-snvs) INCLUDE_NON_SNVS=1; shift;;
    --collapse-hp) COLLAPSE_HP=1; INCLUDE_NON_SNVS=1; shift;;
    --hp-fp-rate) HP_FP_RATE="$2"; shift 2;;
    --hp-fn-rate) HP_FN_RATE="$2"; shift 2;;
    --stop-after)
      STOP_AFTER="$2"; shift 2
      [[ "$STOP_AFTER" == "call" ]] || { echo "--stop-after must be call (got $STOP_AFTER)" >&2; exit 1; }
      ;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" && -n "$PROFILE" ]] || {
  echo "usage: --outdir DIR --profile {hifi,ont-r10,ont-denoised} [--sample S] [--fp F] [--fn F]" >&2
  echo "  [--minimal-ac N] [--vaf V] [--denoise-vaf V] [--frequency-threshold F]" >&2
  echo "  [--permutation-frequency-threshold F] [--strand-bias-threshold F]" >&2
  echo "  [--indel-false-threshold F] [--min-edge-reads N] [--min-hf F] [--max-hf F]" >&2
  echo "  [--include-non-snvs] [--collapse-hp] [--hp-fp-rate F] [--hp-fn-rate F] [--stop-after call]" >&2
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
minimap2 -ax "$MMPRESET" -t "$HIMITO_THREADS" "$REF" "$FQ" \
  | samtools sort -o "$BAM" -
samtools index "$BAM"

# Denoising changes the graph input. quick-start still passes the original
# BAM to call for strand evidence, linkage rescue, and genotyping.
DENOISE_INDELS="${DENOISE_INDELS:-1}"
BUILD_BAM="$BAM"
if [[ "$DTYPE" == ont-denoised ]]; then
  DENOISED="$HDIR/aln.denoised.bam"
  DENOISE_ARGS=(--stats "$HDIR/denoise_stats.json")
  [[ -n "$DENOISE_VAF" ]] && DENOISE_ARGS+=(--vaf "$DENOISE_VAF")
  [[ "$DENOISE_INDELS" == 1 ]] && DENOISE_ARGS+=(--indels)
  "$HIMITO" denoise -i "$BAM" -o "$DENOISED" -r "$REF" -d "$DTYPE" "${DENOISE_ARGS[@]}"
  samtools index "$DENOISED"
  BUILD_BAM="$DENOISED"
fi

# Build anchor graph (input can be a BAM). --min-edge-reads matches
# quick-start's hardcoded edge-read gate (main.rs QuickStart: `build::start(...,
# 2)`); the CLI's own default is 1, which is why this must be passed explicitly.
"$HIMITO" build -i "$BUILD_BAM" -r "$REF" -k "$KMER" -o "$HDIR/sim.gfa" -l 3000 \
  --min-edge-reads "$MIN_EDGE_READS"

# Leave platform-specific VAF, indel and permutation thresholds to Himito.
CALL_ARGS=(
  -g "$HDIR/sim.gfa" -r "$REF" -s "$SAMPLE" -d "$DTYPE"
  -o "$HDIR/sim.vcf" -k "$KMER" --input-bam "$BAM"
  -m "$MINIMAL_AC" --strand-bias-threshold "$STRAND_BIAS_THRESHOLD"
)
[[ -n "$VAF" ]] && CALL_ARGS+=(-v "$VAF")
[[ -n "$FREQ_THRESHOLD" ]] && CALL_ARGS+=(-f "$FREQ_THRESHOLD")
[[ -n "$PERM_FREQ_THRESHOLD" ]] && CALL_ARGS+=(--permutation-frequency-threshold "$PERM_FREQ_THRESHOLD")
[[ -n "$INDEL_FALSE_THRESHOLD" ]] && CALL_ARGS+=(--indel-false-threshold "$INDEL_FALSE_THRESHOLD")
"$HIMITO" call "${CALL_ARGS[@]}"

# The HP-rate sweep re-runs lineage itself. Leave the matrix, VCF, reads, and truth in place.
if [[ "$STOP_AFTER" == "call" ]]; then
  echo "stop-after call: skipped Himito lineage ($HDIR/sim.matrix.csv, $HDIR/sim.vcf)" >&2
  exit 0
fi

# quick-start also runs a NUMT/methylation filter before build and asm/methyl
# steps after call; both are no-ops for this harness's simulated reads (no
# MM/ML tags to classify on, well under the 50000-read downsample cap) and
# produce artifacts (fasta, methylation bed) this harness doesn't score, so
# they're intentionally not reproduced here.

# Lineage: SCITE mutation-tree reconstruction over heteroplasmic SNVs in
# [--min-hf, --max-hf). Unless --fp/--fn are passed, the rates are the -d preset;
# variants overlapping D310/HVS1-C/CA514 use twice those unless --hp-fp/fn-rate are set.
LINEAGE_ARGS=(
  -m "$HDIR/sim.matrix.csv" -v "$HDIR/sim.vcf"
  --min-hf "$MIN_HF" --max-hf "$MAX_HF"
  -d "$DTYPE"
  -o "$HDIR/sim_lineage"
)
[[ -n "$FP" ]] && LINEAGE_ARGS+=(--fp-rate "$FP")
[[ -n "$FN" ]] && LINEAGE_ARGS+=(--fn-rate "$FN")
[[ "$INCLUDE_NON_SNVS" == 1 ]] && LINEAGE_ARGS+=(--include-non-snvs)
[[ "$COLLAPSE_HP" == 1 ]] && LINEAGE_ARGS+=(--collapse-hp)
[[ -n "$HP_FP_RATE" ]] && LINEAGE_ARGS+=(--hp-fp-rate "$HP_FP_RATE")
[[ -n "$HP_FN_RATE" ]] && LINEAGE_ARGS+=(--hp-fn-rate "$HP_FN_RATE")
"$HIMITO" lineage "${LINEAGE_ARGS[@]}"

echo "himito done: $HDIR/sim_lineage.mutation_tree.tsv"
