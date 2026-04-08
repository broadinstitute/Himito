#bin/bash

# Usage:
#   ./Mitorsaw.sh <input_bam> <ref_fa> <output_prefix>
#
# Example:
#   ./Mitorsaw.sh ./sample.bam ./rCRS.fasta ./results/sample

input_file=$1
ref_fa=$2
output_prefix=$3

mitorsaw haplotype \
            --reference $ref_fa \
            --bam $input_file \
            --minimum-maf 0.01 \
            --output-vcf $output_prefix.mitorsaw.vcf.gz \
            --output-hap-stats $output_prefix.mitorsaw.stat
