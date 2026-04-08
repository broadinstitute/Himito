#bin/bash

#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./MitoHiFi.sh <input_folder> <ref_fa> <ref_gb> <output_prefix_dir>
#
# Example:
#   ./MitoHiFi.sh ./bam_inputs ./rCRS.fasta ./rCRS.gb ./results

input_folder="${1:?input_folder is required}"
ref_fa="${2:?ref_fa is required}"          # kept for interface compatibility
ref_gb="${3:?ref_gb is required}"
output_prefix_dir="${4:?output_prefix_dir is required}"

mkdir -p "$output_prefix_dir"

# TSV log with runtime + memory from /usr/bin/time
runtime_log="${output_prefix_dir}/mitohifi_all_runtime.tsv"
echo -e "bam_file\telapsed\tuser_sec\tsys_sec\tmax_rss_kb\texit_code" > "$runtime_log"

shopt -s nullglob
fastq_files=("$input_folder"/*.fastq)

if [[ ${#fastq_files[@]} -eq 0 ]]; then
  echo "No .fastq files found in: $input_folder" >&2
  exit 1
fi

for input_file in "${fastq_files[@]}"; do
  base="$(basename "$input_file" .fastq)"
  out_prefix="${output_prefix_dir}/${base}"
  mkdir -p "./${base}"
  cd "./${base}"
  cp "$input_file" ${base}.fq
  cp "$ref_fa" "rCRS.fa"
  cp "$ref_gb" "rCRS.gb"
  echo "Processing: $input_file"

  # /usr/bin/time writes metrics to stderr; capture into per-sample file
  time_log="${output_prefix_dir}/${base}.time.txt"
  set +e
  /usr/bin/time -f "elapsed=%E user=%U sys=%S maxrss_kb=%M exit=%x" \
      docker run --rm -it \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  ghcr.io/marcelauliano/mitohifi:master \
  mitohifi.py \
    -f "rCRS.fa" \
    -r "${base}.fq" \
    -g "rCRS.gb" \
    -t 8 \
    2> "$time_log"
  exit_code=$?
  set -e

  # Parse timing fields from the time log
  elapsed="$(awk -F'[ =]+' '/^elapsed=/{print $2}' "$time_log" | tail -n1)"
  user_sec="$(awk -F'[ =]+' '/^elapsed=/{print $4}' "$time_log" | tail -n1)"
  sys_sec="$(awk -F'[ =]+' '/^elapsed=/{print $6}' "$time_log" | tail -n1)"
  max_rss_kb="$(awk -F'[ =]+' '/^elapsed=/{print $8}' "$time_log" | tail -n1)"

  echo -e "${input_file}\t${elapsed:-NA}\t${user_sec:-NA}\t${sys_sec:-NA}\t${max_rss_kb:-NA}\t${exit_code}" >> "$runtime_log"

  if [[ $exit_code -ne 0 ]]; then
    echo "FAILED: $input_file (exit $exit_code). See: $time_log" >&2
  fi
  cd ..

done

echo "Done. Runtime summary: $runtime_log"
