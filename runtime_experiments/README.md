# Runtime Experiments

This folder contains batch scripts for benchmarking runtime and memory usage of Himito and comparison tools on multiple inputs.

Each script:

- scans an input folder (`*.bam` or `*.gfa` depending on the script),
- runs one tool/workflow per file,
- records `/usr/bin/time` metrics (`elapsed`, `user`, `sys`, `maxrss_kb`, `exit`),
- writes a summary TSV in the output folder.

## Prerequisites

- Bash (`set -euo pipefail` compatible)
- `/usr/bin/time`
- Tool binaries on `PATH` as needed:
  - `Himito`
  - `mitorsaw`
  - `nextflow` (for mtDNA-Server script)

## Scripts

- `Himito_filter.sh`  
  Runs `Himito filter` for each BAM.

- `Himito_build.sh`  
  Runs `Himito build` for each BAM.

- `Himito_asm.sh`  
  Runs `Himito asm` for each GFA.

- `Himito_call.sh`  
  Runs `Himito call` for each GFA.

- `Himito_call_asm.sh`  
  Convenience wrapper to run call + asm sequence for a single sample.

- `Himito_all.sh`  
  Batch wrapper that calls `Himito_call_asm.sh` for each BAM.

- `Mitorsaw.sh` / `Mitorsaw_all.sh`  
  Runs MitoRSaw haplotype calling and logs runtime.

- `mtdnaserver.sh`  
  Runs `nextflow run genepi/mtdna-server-2` per BAM (fusion mode) with per-sample generated config.

- `MitoHiFi.sh`  
  Runtime batch helper for MitoHiFi experiments.

## Common Usage Pattern

Most scripts expect:

```bash
bash <script>.sh <input_folder> <reference_fasta> <output_folder>
```

Examples:

```bash
bash Himito_filter.sh ./bam_inputs ./NC_012920.1.fasta ./runtime_out/filter
bash Himito_call.sh ./gfa_inputs ./NC_012920.1.fasta ./runtime_out/call
bash Mitorsaw_all.sh ./bam_inputs ./NC_012920.1.fasta ./runtime_out/mitorsaw
bash mtdnaserver.sh ./bam_inputs ./NC_012920.1.fasta ./runtime_out/mtdnaserver
```

## Output Logs

Each run writes:

- per-sample time log: `<output_folder>/<sample>.time.txt`
- aggregated table: `*_runtime.tsv` in `<output_folder>`

`maxrss_kb` is peak resident memory in KB from `/usr/bin/time`.

## Notes

- Several script headers still use old comments like `run_himito_build_batch.sh`; use the actual file names in this folder.
- For Nextflow/docker workflows, ensure your runtime environment has working container execution setup.
