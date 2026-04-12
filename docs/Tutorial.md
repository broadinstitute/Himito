# Himito Docker Tutorial (Non-WDL)

This tutorial is for users who want to run Himito directly from the command line (without WDL/Cromwell).
It starts from a BAM file and runs the full pipeline in Docker:

`Filter -> Build -> Asm -> Call -> Methyl`

It is designed to be easy to reproduce and adapt into Snakemake/Nextflow later.

## 1) Prerequisites

- Docker installed and running.
- Input BAM and index (`.bai`).
- Mitochondrial reference FASTA (for example `NC_012920.1.fasta`).

Himito Docker image used in this tutorial:

```bash
docker pull us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1.1.0
```

## 2) Prepare a working directory or use the `test_data` folder as as an example

```bash
git clone https://github.com/broadinstitute/Himito.git
cd Himito/test_data
```

Put your files here (or use the test dataset):

- `HG002_chrM.bam`
- `HG002_chrM.bam.bai`
- `rCRS.fasta`

Set variables once:

```bash
IMG="us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1.1.0"
WD="$PWD"
SAMPLE="HG002_chrM"
REF="rCRS.fasta"
CHR="chrM"
K=21
```

Optional helper alias to shorten Docker commands:

```bash
alias himito_docker='docker run --rm -u "$(id -u):$(id -g)" -v "$WD":/work -w /work '"$IMG"
```

If you do not use the alias, each command below has a full Docker form.

## 3) One-command quick run 

If you only want a quick end-to-end test:

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$WD":/work -w /work "$IMG" \
  /Himito/target/release/Himito quick-start \
  -i "${SAMPLE}.bam" \
  -o "${SAMPLE}" \
  -r "$REF" \
  -s "$SAMPLE" \
  -d pacbio
```

## 4) Step-by-step pipeline

### Step A: Filter (remove NUMTs-like reads)

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$WD":/work -w /work "$IMG" \
  /Himito/target/release/Himito filter \
  -i "${SAMPLE}.bam" \
  -c "$CHR" \
  -m "${SAMPLE}.mt.bam" \
  -n "${SAMPLE}.numts.bam"
```

Expected outputs:

- `${SAMPLE}.mt.bam`
- `${SAMPLE}.numts.bam`

### Step B: Build (construct mitochondrial graph)

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$WD":/work -w /work "$IMG" \
  /Himito/target/release/Himito build \
  -i "${SAMPLE}.mt.bam" \
  -k "$K" \
  -r "$REF" \
  -o "${SAMPLE}.gfa"
```

Expected output:

- `${SAMPLE}.gfa`

### Step C: Asm (extract major haplotype)

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$WD":/work -w /work "$IMG" \
  /Himito/target/release/Himito asm \
  -g "${SAMPLE}.gfa" \
  -o "${SAMPLE}.fasta" \
  -s "$SAMPLE"
```

Expected output:

- `${SAMPLE}.fasta`

### Step D: Call (variant calling)

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$WD":/work -w /work "$IMG" \
  /Himito/target/release/Himito call \
  -g "${SAMPLE}.gfa" \
  -r "$REF" \
  -k "$K" \
  -s "$SAMPLE" \
  -d pacbio \
  -o "${SAMPLE}.vcf"
```

Expected outputs include:

- `${SAMPLE}.vcf`
- `${SAMPLE}.matrix.csv`
- `${SAMPLE}.raw_matrix.csv`

For ONT reads, change `-d pacbio` to your `-d ont-r9` or `-d ont-r10` used in your dataset setup.

### Step E: Methyl (methylation calls)

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$WD":/work -w /work "$IMG" \
  /Himito/target/release/Himito methyl \
  -g "${SAMPLE}.gfa" \
  -p 0.7 \
  -b "${SAMPLE}.mt.bam" \
  -o "${SAMPLE}.bed"
```

Expected output:

- `${SAMPLE}.bed`

## 5) Summary of inputs and outputs

```text
  input/
    ${SAMPLE}.bam
    ${SAMPLE}.bam.bai
    rCRS.fasta
  output/
    ${SAMPLE}.mt.bam
    ${SAMPLE}.numts.bam
    ${SAMPLE}.gfa
    ${SAMPLE}.methyl.gfa
    ${SAMPLE}.fasta
    ${SAMPLE}.vcf
    ${SAMPLE}.matrix.csv
    ${SAMPLE}.raw_matrix.csv
    ${SAMPLE}.bed
    ${SAMPLE}.methylation_per_read.csv

```