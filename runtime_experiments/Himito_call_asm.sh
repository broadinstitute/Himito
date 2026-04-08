#bin/bash

# Usage:
#   ./Himito_call_asm.sh <input_bam> <ref_fa> <output_prefix>
#
# Example:
#   ./Himito_call_asm.sh ./sample.bam ./rCRS.fasta ./results/sample

input_file=$1
ref_fa=$2
output_prefix=$3

Himito filter -i $input_file -c chrM -m $output_prefix.mt.bam -n $output_prefix.numts.bam
Himito build -i $output_prefix.mt.bam -k 21 -r $ref_fa -o $output_prefix.gfa 
Himito call -g $output_prefix.gfa -r $ref_fa -k 21 -s $output_prefix -o $output_prefix.vcf
Himito asm -g $output_prefix.gfa -s $output_prefix -o $output_prefix.fasta        

