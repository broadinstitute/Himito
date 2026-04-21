#bin/bash

#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./Himito_quickstart.sh <input_bam> <ref_fa> <output_prefix> <sample_id> <data_type>
#
# Example:
#   ./Himito_quickstart.sh ./test_data/HG002_chrM_100x.bam ./test_data/rCRS.fasta ./results/HG002_chrM_100x HG002 pacbio

input_bam="${1:?input_file is required}"
ref_fa="${2:?ref_fa is required}"          # kept for interface compatibility
output_prefix="${3:?output_prefix is required}"
sample_id="${4:?sample_id is required}"
data_type="${5:?data_type is required}"


echo "Processing: $input_bam"


IMG="us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1.1.0"
WD="$PWD"

# /usr/bin/time writes metrics to stderr; capture into per-sample file
time_log="${output_prefix}.time.txt"
set +e
/usr/bin/time -f "elapsed=%E user=%U sys=%S maxrss_kb=%M exit=%x" \
docker run --rm -u "$(id -u):$(id -g)" -v "$WD":/work -w /work "$IMG" \
  /Himito/target/release/Himito quick-start \
  -i "${input_bam}" \
  -o "${output_prefix}" \
  -r "${ref_fa}" \
  -s "${sample_id}" \
  -d "${data_type}" \
    2> "$time_log"


