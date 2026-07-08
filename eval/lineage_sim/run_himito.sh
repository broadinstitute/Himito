#!/usr/bin/env bash
# Align pooled simulated reads to rCRS and run Himito build -> call -> lineage.
set -euo pipefail

HIMITO=/Users/suhang/Analysis/Himito/target/release/Himito
REF=/Users/suhang/Analysis/Himito/rCRS.fasta
OUTDIR="" PROFILE="" SAMPLE="SIM" FP=0.001 FN=0.05
# Call-tuning defaults: permissive for simulated data
# -m 0: keep all reads in matrix (simulated reads carry few alts each; m>=1 empties the matrix)
# -v 0.005: capture low-frequency clones (default 0.01 loses sub-1% clones)
# -p 1.0: disable permutation test (clean simulated variants have balanced strands)
MINIMAL_AC=0 VAF=0.005 PVAL=1.0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --fp) FP="$2"; shift 2;;
    --fn) FN="$2"; shift 2;;
    --minimal-ac) MINIMAL_AC="$2"; shift 2;;
    --vaf) VAF="$2"; shift 2;;
    --pval) PVAL="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 1;;
  esac
done
[[ -n "$OUTDIR" && -n "$PROFILE" ]] || { echo "usage: --outdir DIR --profile {hifi,ont-r10} [--sample S] [--fp F] [--fn F] [--minimal-ac N] [--vaf V] [--pval P]" >&2; exit 1; }

case "$PROFILE" in
  hifi)    MMPRESET="map-hifi"; DTYPE="pacbio";;
  ont-r10) MMPRESET="map-ont";  DTYPE="ont-r10";;
  *) echo "profile must be hifi or ont-r10" >&2; exit 1;;
esac

HDIR="$OUTDIR/himito"; mkdir -p "$HDIR"
FQ="$OUTDIR/reads/reads.fastq.gz"
BAM="$HDIR/aln.sorted.bam"

minimap2 -ax "$MMPRESET" -t 4 "$REF" "$FQ" 2>/dev/null \
  | samtools sort -o "$BAM" -
samtools index "$BAM"

# Build anchor graph (input can be a BAM).
"$HIMITO" build -i "$BAM" -r "$REF" -k 21 -o "$HDIR/sim.gfa" -l 3000

# Call variants: -o is the VCF; matrix.csv is derived as <o>.matrix.csv.
"$HIMITO" call -g "$HDIR/sim.gfa" -r "$REF" -s "$SAMPLE" -d "$DTYPE" \
  -o "$HDIR/sim.vcf" -k 21 --input-bam "$BAM" \
  -m "$MINIMAL_AC" -v "$VAF" -p "$PVAL"

# Lineage: SCITE mutation-tree reconstruction.
"$HIMITO" lineage -m "$HDIR/sim.matrix.csv" -v "$HDIR/sim.vcf" \
  --fp-rate "$FP" --fn-rate "$FN" --min-hf 0.01 --max-hf 0.95 \
  -o "$HDIR/sim_lineage"

echo "himito done: $HDIR/sim_lineage.mutation_tree.tsv"
