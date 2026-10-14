#!/bin/bash
# usage: ev.sh in.vcf tag   -> PASS-only vcfeval vs HapMap truth, prints summary line + FP/FN counts
set -euo pipefail
RTG=${RTG:-/Users/suhang/Analysis/rtg-tools-3.13/rtg}
W=${W:-$PWD}
TRUTH=${TRUTH:-$W/truth.vcf.gz}
SDF=${SDF:-../rtg_ref}
IN=$1; TAG=$2
awk -F'\t' '/^#/ || $7=="PASS"' "$IN" > $W/$TAG.pass.vcf
rm -f $W/$TAG.pass.vcf.gz*
$RTG bgzip $W/$TAG.pass.vcf >/dev/null 2>&1
$RTG index -f vcf $W/$TAG.pass.vcf.gz >/dev/null 2>&1
rm -rf $W/$TAG.eval
$RTG vcfeval -b $TRUTH -c $W/$TAG.pass.vcf.gz -o $W/$TAG.eval -t $SDF --squash-ploidy --sample ALT,ALT >/dev/null 2>&1
printf "%-28s " "$TAG"; tail -1 $W/$TAG.eval/summary.txt
