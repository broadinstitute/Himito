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
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)

## Installation

### Prerequisites

- Rust programming language
- Long-read sequencing data (BAM format)
- Mitochondrial reference genome (FASTA format)

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

Himito Docker can be downloaded at docker hub: [hangsuunc/Himito:v1](https://hub.docker.com/r/hangsuunc/himito).

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
- -k, --kmer-size <INT>: K-mer size for graph construction
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
- -q --query_length <INT>: the Kmer length for querying the short read msBWT, path with length less than K will not be kmerized

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
- -g, --graph <FILE>: Input annotated graph file (GFA format)
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
### BAM File Requirements
- Must be coordinate-sorted and indexed (.bai file)

- Should contain long reads (PacBio, Oxford Nanopore)

- Reads should be mapped to a reference genome including mitochondrial chromosome

- For methylation analysis: BAM should contain modification tags (MM, ML tags)

### Reference Genome
- Mitochondrial reference in FASTA format
- Common references: NC_012920.1 (revised Cambridge Reference Sequence)
- Should match the reference used for initial read mapping

## Output Files
### Graph Files (.gfa)
- Contains the mitochondrial genome graph structure
- Can be visualized with tools like Bandage
- Used as input for variant calling and assembly and methylation analysis
### Variant Files (.vcf)
- Standard VCF format with mitochondrial variants
- Includes heteroplasmy frequency estimates
- Compatible with standard VCF analysis tools
### Binary matrix for Presence of variants in each read (.csv)
- presence-absence matrix representing variants in each read
### Assembly Files (.fasta)
- Primary mitochondrial genome sequence
- Represents the major haplotype
- Can be used for phylogenetic analysis
### Methylation Files (.bed)
- BED format with methylation sites
- Includes probability scores and coverage information
- Compatible with methylation analysis tools
### read-level methylation matrix (.csv)
- Element is the methylation likelihood of each CpG site

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

## Getting Help
For additional support: Please submit an issue in the Himito GitHub repository

## Citation
The preprint is coming soon...