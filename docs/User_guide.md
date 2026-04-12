# Himito User Guide

## Overview

Himito is a specialized tool for building mitochondrial anchor-based graphical genomes from long-read sequencing data. It provides a comprehensive pipeline for mitochondrial genome analysis, including filtering nuclear mitochondrial sequences (NUMTs), assembling major haplotypes, calling variants, and analyzing methylation signals.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Commands](#commands)
- [Workflows](#workflows)
- [Input Requirements](#input-requirements)
- [Output Files](#output-files)
- [VCF output schema](#vcf-output-schema)
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)

## Installation

### Prerequisites

- Rust programming language
- Long-read sequencing data (BAM format)
- Mitochondrial reference genome (FASTA format)
- (Optional) Paired short-read sequencing data (BAM or FASTA format)

### Step 1: Install Rust

### install rust
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
### Step 2: Install Himito

### install and run Himito
```
git clone https://github.com/broadinstitute/Himito.git
cd Himito
cargo build --release
```

## Quick Start

Here's a basic workflow to get started with Himito:

```
# filter NUMTs-derived reads
./target/release/Himito filter -i <input.bam> -c <chromosome in bam, e.g. "chrM"> -m <mt_output.bam> -n <numts_output.bam>

# construct graph
./target/release/Himito build -i <mt_output.bam> -k <kmer_size> -r <NC_012920.1.fasta> -o <output.gfa> 

# (Optional) prune graph using short read sequence data
msbwt2-build -o sr_msbwt.npy <srWGS.chrM.fasta.gz>
./target/release/Himito correct -g <output.gfa> -b <bwt_file, e.g. sr_msbwt.npy> -o <corrected.gfa> -m <minimal_supporting_sr> -q <query_length, should be less than short read length>

# call variants from graph
./target/release/Himito call -g <output.gfa> -r <NC_012920.1.fasta> -k <kmer_size> -s <sampleid> -o <output.vcf>

# extract major haplotype from graph
./target/release/Himito asm -g <output.gfa>  -o <output.majorhaplotpe.fasta> -s <header string, e.g. "HG002 major haplotype">

# call methylation signals
./target/release/Himito methyl -g <output.annotated.gfa> -p <min_prob> -b <mt_test.bam> -o <methyl.bed>
```

Himito Docker can be downloaded at docker hub: [us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1.1.0](us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1.1.0).

## Commands
### filter - Remove NUMTs Reads
Separates genuine mitochondrial reads from nuclear mitochondrial sequences (NUMTs).

```
./target/release/Himito filter [OPTIONS]
```
Required Parameters:
- -i, --input <FILE>: Input BAM file containing long reads
- -c, --chromosome <STRING>: Mitochondrial chromosome name (e.g., "chrM", "MT", "M")
- -m, --mt-output <FILE>: Output BAM file for mitochondrial reads
- -n, --numts-output <FILE>: Output BAM file for NUMTs reads

Example:
```
./target/release/Himito filter -i HG002.bam -c chrM -m HG002_mt.bam -n HG002_numts.bam
```

### build - Construct Mitochondrial Graph
Creates a graph representation of the mitochondrial genome from filtered reads.

```
./target/release/Himito build [OPTIONS]
```
Required Parameters:
- -i, --input <FILE>: Input BAM file (mitochondrial reads only)
- -k, --kmer-size <INT>: K-mer size for graph construction, please choose your k carefully: too small k will give you a hair ball and mess up the downstream assembly; too large k will not compress the reads. Maybe draw your graph and make sure it is a circle.
- -r, --reference <FILE>: Mitochondrial reference genome (FASTA)
- -o, --output <FILE>: Output graph file (GFA format)

Example:
```
./target/release/Himito build -i HG002_mt.bam -k 21 -r rCRS.fasta -o HG002_mt.gfa
```
### correct - Prune graph paths using short read sequence data
Prune paths in the raw graph based on its kmer-counts in the short reads data

```
msbwt2-build -o sr_msbwt.npy <srWGS.chrM.fasta.gz>
./target/release/Himito correct [OPTIONS]
```
Required Parameters:
- -g --graphfile <File>: Input GFA file for short read correction
- -b --bwt_file <File>: msBWT file (.npy) constructed using srWGS reads aligned to chrM
- -o --outputfile <File>: Output GFA file for downstream analysis
- -m --min_support_counts <INT>: the minimal number of short reads supporting each path, paths with supports less than this threshold will be trimmed.
- -q --query_length <INT>: the Kmer length for querying the short read msBWT, path with length less than K will not be kmerized. K should be less than the short read length.

Example:
```
msbwt2-build -o sr_msbwt.npy srWGS.chrM.fasta.gz

./target/release/Himito correct -g HG002_mt.gfa -b sr_msbwt.npy -o HG002_mt_sr.gfa -m 10 -q 99

```

### call - Variant Calling
Identifies homoplasmic and heteroplasmic variants from the mitochondrial graph.

```
./target/release/Himito call [OPTIONS]
```

Required Parameters:
- -g, --graph <FILE>: Input graph file (GFA format)
- -r, --reference <FILE>: Mitochondrial reference genome (FASTA)
- -k, --kmer-size <INT>: K-mer size (should match build step)
- -s, --sample-id <STRING>: Sample identifier
- -d, --data-type <String>: Long Read Sequencing technology (pacbio or ont)
- -o, --output <FILE>: Output VCF file

Example:

```
./target/release/Himito call -g HG002_mt.gfa -r rCRS.fasta -k 21 -s HG002 -o HG002_variants.vcf
```

### asm - Primary Assembly Extraction
Extracts the major haplotype sequence from the mitochondrial graph.

```
./target/release/Himito asm [OPTIONS]
```

Required Parameters:
- -g, --graph <FILE>: Input graph file (GFA format)
- -o, --output <FILE>: Output FASTA file
- -s, --header <STRING>: FASTA header string

Example:
```
./target/release/Himito asm -g HG002_mt.gfa -o HG002_major_hap.fasta -s "HG002 mitochondrial major haplotype"
```

### methyl - Methylation Analysis
Analyzes methylation signals from the mitochondrial reads.

```
./target/release/Himito methyl [OPTIONS]
```

Required Parameters:
- -g, --graph <FILE>: Input annotated graph file (GFA format), should be the output of Himito call
- -p, --min-prob <FLOAT>: Minimum probability threshold for methylation calls
- -b, --bam <FILE>: Input BAM file with methylation tags
- -o, --output <FILE>: Output BED file

Example:
```
./target/release/Himito methyl -g HG002_mt.annotated.gfa -p 0.7 -b HG002_mt.bam -o HG002_methylation.bed
```

## Workflows
### WDL Workflows

Himito includes WDL (Workflow Description Language) workflows for running on cloud platforms. Check the wdl/ directory for available workflows.

Basic Analysis Pipeline (wdl/Himito_methyl.wdl)

- Quality Control: Ensure your BAM file is properly indexed and contains long reads
- Filtering: Remove NUMTs contamination using the filter command
- Graph Construction: Build the mitochondrial graph with appropriate k-mer size
- Variant Calling: Identify variants and estimate heteroplasmy levels
- Assembly: Extract consensus sequences for downstream analysis
- Methylation: Analyze epigenetic modifications (if data available)


## Input Requirements
### lrWGS BAM File 
- Must be coordinate-sorted and indexed (.bai file)

- Should contain long reads (PacBio, Oxford Nanopore)

- Reads should be mapped to a reference genome including mitochondrial chromosome

- For methylation analysis: BAM should contain modification tags (MM, ML tags)

### (Optional) srWGS BAM File
- BAM file should be sorted and indexed

- (Recommended) Subset to reads aligned to chrM, convert BAM to FASTA/FASTA.gz

- Construct msBWT using [rust-msBWT](https://github.com/HudsonAlpha/rust-msbwt)

- Run ```Himito correct``` with proper parameters according to the srWGS data features




### Reference Genome
- Mitochondrial reference in FASTA format
- Common references: NC_012920.1 (revised Cambridge Reference Sequence)
- Should match the reference used for initial read mapping

## Output Files

### Default file naming (`quick-start` and stepwise runs)

Many outputs share the same basename you pass with `-o` / `--output-prefix` (`quick-start`) or the explicit paths you pass to each subcommand. Typical artifacts:

| Extension / pattern | Description |
| --------------------- | ----------- |
| `.mt.bam`, `.numts.bam` | Mitochondrial and NUMT-partitioned reads (`filter` / `quick-start`) |
| `.gfa` | Sequence graph (`build` / `quick-start`; methylation step may update annotated graph) |
| `.fasta` | Major-haplotype assembly (`asm` / `quick-start`) |
| `.vcf` | Variant calls after graph-based calling and permutation filtering (`call` / `quick-start`) |
| `.bed` | Methylation summary (`methyl` / `quick-start`) |
| `.raw_matrix.csv`, `.matrix.csv` | Per-read variant support matrices written next to the VCF basename (`call` / `quick-start`) |

The **matrix CSV** has header `variant`, then one column per read name. Each row is a variant identifier; cell values are non-negative counts of support in that read (see `call` in the codebase). `.raw_matrix.csv` is produced before permutation-based filtering; `.matrix.csv` matches the variants retained in the final VCF.

### Graph Files (.gfa)
- Contains the mitochondrial genome graph structure
- Can be visualized with tools like Bandage
- Used as input for variant calling and assembly and methylation analysis

Example (GFA):

```text
H	VN:Z:1.0
S	A000022	TAACCACTCACGGGAGCTCTC	PG:J:{"pos":22}
S	A000044	ATGCATTTGGTATTTTCGTCT	PG:J:{"pos":44}
S	A000066	GGGGGTATGCACGCGATAGCA	PG:J:{"pos":66}
S	A000088	TGCGAGACGCTGGAGCCGGAG	PG:J:{"pos":88}
S	A000110	ACCCTATGTCGCAGTATCTGT	PG:J:{"pos":110}
...
S	E00001.0000	AGAAGTTATTATCTCGAACTGACACTGA	PG:J:{"dst":["A012540"],"reads":["m64012_190920_173625/161154834/ccs"],"src":["SOURCE"]}	RC:i:1
S	E00001.0001	CTCTCCCTAAGCTTCAACTAGA	PG:J:{"dst":["A012584"],"reads":["m64012_190920_173625/161154834/ccs","m64011_190901_095311/93847814/ccs"], "src":["A012540"],"variants":"36=1D7="}	RC:i:15
...
L	A010010	+	E00084.0016	+	0M
L	E00084.0016	+	A010054	+	0M
L	A010384	+	E00084.0017	+	0M

```

In this example, `S` lines are graph segments (Anchor nodes, e.g. A010010; edge nodes e.g. E00084.0016) and `L` lines connect segments with an overlap (`0M`).



### Variant Files (.vcf)
Himito writes **VCF v4.2** (uncompressed). Records are **sorted by position** (then type, REF, ALT for ties). The final list passed to `write_vcf` is produced after **allele-count / heteroplasmy filtering** and **permutation-based filtering** when you use `call` or `quick-start`.

Standard tools can read the file as a normal single-sample VCF. If your reference contig name in the FASTA header is not `chrM`, note that **data lines still use `CHROM=chrM`** in the emitted VCF body while the header **`##contig=<ID=…>`** reflects the first sequence name and length from your reference FASTA.

### VCF output schema

**Header (meta-lines)**  
Himito emits at least:

- `##fileformat=VCFv4.2`
- `##reference=<reference_sequence_name>`
- `##contig=<ID=…,length=…>` for the mitochondrial reference contig (from your FASTA’s first sequence name and length)
- `##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">`
- `##FORMAT=<ID=GT,…>`, `##FORMAT=<ID=AD,…>`, `##FORMAT=<ID=HF,…>` (genotype, allele depth, heteroplasmic fraction)

**Column header**

`#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  <sample_id>`

`<sample_id>` is the value you pass to `-s` / `--sample-id` (`call` or `quick-start`).

**Body fields**

| Field | Content |
| ----- | ------- |
| CHROM | `chrM` on each data line |
| POS | 1-based position on the reference. |
| ID | `.` |
| REF / ALT | Reference and alternate alleles (SNP or symbolic indel representation as produced by the caller) |
| QUAL | `.` |
| FILTER | `.` |
| INFO | `DP=<integer>` total read depth at the variant **start** position used for depth lookup |
| FORMAT | Fixed string `GT:AD:HF` |
| Sample | `GT:AD:HF` with **GT** `1` for called alternate, **AD** = alternate **allele read count** (supporting reads), **HF** = `AD / DP` at that position (heteroplasmic fraction, 0–1 float).|

Example (Vcf):

```text
##fileformat=VCFv4.2
##reference=chrM
##contig=<ID=chrM,length=16569>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=HF,Number=1,Type=Float,Description="Heteroplasmic Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG002
chrM	263	.	A	G	.	.	DP=1792	GT:AD:HF	1:1792:1
chrM	309	.	C	CCT	.	.	DP=1752	GT:AD:HF	1:1229:0.701484
chrM	310	.	T	C	.	.	DP=1752	GT:AD:HF	1:1636:0.93378997
chrM	456	.	C	T	.	.	DP=1753	GT:AD:HF	1:1709:0.9749002
chrM	457	.	C	T	.	.	DP=1753	GT:AD:HF	1:19:0.010838563
chrM	750	.	A	G	.	.	DP=1825	GT:AD:HF	1:1815:0.99452055
...
```

### Binary matrix for Presence of variants in each read (.csv)
- Rows: variant identifiers; columns: read IDs; values encode per-read support (see table above). Two files: pre- and post-permutation filter (`.raw_matrix.csv` vs `.matrix.csv`).

Example (read matrix):

```text
variant,m64011_190830_220126/100204788/ccs,m64011_190830_220126/100206196/ccs, ...
m.13376T>C,0,1,...
m.15175C>T,1,1,...
```

### Assembly Files (.fasta)
- Primary mitochondrial genome sequence
- Represents the major haplotype
- Can be used for phylogenetic analysis

Example (assembly fasta)
```text
>HG002 
TAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGA
TAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCT
GCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAA
GTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCC
...
```

### Methylation Files (.bed)
- BED format with methylation sites
- Includes probability scores and coverage information
- Compatible with methylation analysis tools
### read-level methylation matrix (.csv)
- Element is the methylation likelihood of each CpG site

```text
##fileformat=BED
##haplotype=majorhaplotype
#CHROM	Ref_start	Ref_end	Asm_start	Asm_end	Motif	Mod_rate	Unmod_rate	Cov	Mod_count	Unmod_count
ChrM	33	34	11	12	CG	0.23796498905908095	0.7620350109409191	1828	435	1393
ChrM	61	62	39	40	CG	0.20044419766796223	0.7995558023320377	1801	361	1440
ChrM	78	79	56	57	CG	0.047592695074709465	0.9524073049252906	1807	86	1721
...
```

## Best Practices
### Data Preparation
1. Use whole-genome long-read data (>50X mitochondrial coverage)
- Use recent mitochondrial reference sequences (rCRS)
### Parameter Optimization
- Test different k-mer sizes for your specific dataset
- Adjust methylation thresholds based on your coverage and accuracy requirements
- Consider running multiple iterations with different parameters
### Quality Control
- Check the number of reads filtered as NUMTs vs. genuine mitochondrial
- Verify graph connectivity and complexity
- Validate variant calls against known mitochondrial polymorphisms

## Troubleshooting
### Common Issues
#### Low mitochondrial read count after filtering
- Check chromosome naming convention (chrM vs. MT vs. M)
- Verify input BAM contains mitochondrial reads
- Consider adjusting filtering parameters
    - increase parameter -f in Himito filter
#### Graph construction fails
- Check reference genome format and completeness
- Check reference genome naming convention
- Ensure sufficient coverage of mitochondrial genome
- Try different k-mer sizes
#### No variants called
- Verify graph file integrity
- Check reference file is identical to the Himito Build process
- Consider lowering variant calling thresholds
#### Too many false positive (e.g. short indels)
- Adjust parameter --heteroplasmic-frequency-threshold to a larger value (e.g. 0.9 or 0.95)
- Adjust parameter --p-value-threshold to a smaller value (e.g. 0.00001)

## Getting Help
For additional support: Please submit an issue in the Himito GitHub repository

## Citation
The preprint is coming soon...
