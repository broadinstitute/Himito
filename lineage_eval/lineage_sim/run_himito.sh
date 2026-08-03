#!/usr/bin/env bash
# Align pooled simulated reads to rCRS and run Himito build -> call -> lineage.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
HIMITO="${HIMITO:-$REPO/target/release/Himito}"
REF="${REF:-$REPO/rCRS.fasta}"
OUTDIR="" PROFILE="" SAMPLE="SIM" FP=0.001 FN=0.05 KMER=21
# Call-tuning defaults: permissive for simulated data
# -m 0: keep all reads in matrix (simulated reads carry few alts each; m>=1 empties the matrix)
# -v 0.005: capture low-frequency clones (default 0.01 loses sub-1% clones)
# -p 1.0: disable permutation test (clean simulated variants have balanced strands)
MINIMAL_AC=0 VAF=0.001 PVAL=1
# HF band forwarded to `Himito lineage` (same band as run_eval.sh / sweep_fpfn.sh / score_lineage.py).
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
    --kmer) KMER="$2"; shift 2;;
    --min-hf) MIN_HF="$2"; shift 2;;
    --max-hf) MAX_HF="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" && -n "$PROFILE" ]] || { echo "usage: --outdir DIR --profile {hifi,ont-r10} [--sample S] [--fp F] [--fn F] [--minimal-ac N] [--vaf V] [--pval P] [--min-hf F] [--max-hf F]" >&2; exit 1; }

case "$PROFILE" in
  hifi)    MMPRESET="map-hifi"; DTYPE="pacbio";;
  ont-r10) MMPRESET="map-ont";  DTYPE="ont-r10";;
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

# Build anchor graph (input can be a BAM).
"$HIMITO" build -i "$BAM" -r "$REF" -k "$KMER" -o "$HDIR/sim.gfa" -l 3000

# Call variants: -o is the VCF; matrix.csv is derived as <o>.matrix.csv.
"$HIMITO" call -g "$HDIR/sim.gfa" -r "$REF" -s "$SAMPLE" -d "$DTYPE" \
  -o "$HDIR/sim.vcf" -k "$KMER" --input-bam "$BAM" \
  -m "$MINIMAL_AC" -v "$VAF" -p "$PVAL"

# Lineage: SCITE mutation-tree reconstruction.
"$HIMITO" lineage -m "$HDIR/sim.matrix.csv" -v "$HDIR/sim.vcf" \
  --fp-rate "$FP" --fn-rate "$FN" --min-hf "$MIN_HF" --max-hf "$MAX_HF" \
  -o "$HDIR/sim_lineage"

echo "himito done: $HDIR/sim_lineage.mutation_tree.tsv"
