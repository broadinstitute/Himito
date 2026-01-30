
<p align="center">
<img src="https://github.com/user-attachments/assets/3129d44f-f69a-4491-a0cf-29fab5ee88a9" alt="logo" style="width: 30%; height: auto;">
</p>

# Himito 
Building Mitochondrial anchor-based graphical genome from long reads. Filter Numts reads, assemble major haplotypes, call homoplasmic and heteroplasmic variants, analyze methylation signals

## Documentation
[User Guide](https://github.com/broadinstitute/Himito/blob/main/docs/User_guide.md)

## Usage
### install rust
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
### install Himito
```
git clone https://github.com/broadinstitute/Himito.git
cd Himito
cargo build --release
```
### Quick start (Filter-Build-Asm-Call-Methyl)
```
./target/release/Himito quick-start -i <lrWGS.bam> \
                                    -o <output_prefix> \
                                    -r <rCRS.fa> \
                                    -s <HG002> \
                                    -d pacbio

```
## Details
### build Himito graph
```
# filter NUMTs-derived reads
./target/release/Himito filter -i <input.bam> -c <chromosome in bam, e.g. "chrM"> -m <mt_output.bam> -n <numts_output.bam>

# construct graph
./target/release/Himito build -i <mt_output.bam> -k <kmer_size> -r <NC_012920.1.fasta> -o <output.gfa>

```
### (Optional) Prune the Himito graph using paired srWGS data
The graph served as the foundation for downstream assembly, variant calling and methylation analysis.

If you have paired short reads data, you can refine the graph by using ```Himito correct```.

```Himito correct``` trims graph paths according to the occurrence counts of path kmers in the paired srWGS data. We recommend to first subset the srWGS BAM file to reads aligned to chrM, then compressed these reads into a [msBWT](https://github.com/HudsonAlpha/rust-msbwt).
```
msbwt2-build -o sr_msbwt.npy <srWGS.chrM.fasta.gz>
./target/release/Himito correct -g <output.gfa> -b <bwt_file, e.g. sr_msbwt.npy> -o <corrected.gfa> -m <minimal_supporting_sr> -q <query_length, should be less than short read length>
```
### Downstream analysis
```
# call variants from graph, default for pacbio, change <-d ont> to ont data
./target/release/Himito call -g <output.gfa> -r <NC_012920.1.fasta> -k <kmer_size> -s <sampleid> -o <output.vcf>

# extract major haplotype from graph
./target/release/Himito asm -g <output.gfa>  -o <output.majorhaplotpe.fasta> -s <header string, e.g. "HG002 major haplotype">

# call methylation signals
./target/release/Himito methyl -g <output.annotated.gfa> -p <min_prob> -b <mt_test.bam> -o <methyl.bed>

# enumerate all possible haplotypes within windows
./target/release/Himito minorhap -g <output.gfa> -o <output.allhaplotype.fasta> -s <sample_id>
```

### Jan 20th 2026 Updates: NUMTs Breakpoint Calling

Himito can identify and call NUMTs (Nuclear Mitochondrial DNA segments) breakpoints from long-read BAM files. NUMTs are segments of mitochondrial DNA that have been inserted into the nuclear genome. The NUMTs calling functionality (`callnumts.rs`) works by:

1. **Identifying NUMT reads**: Parses supplementary alignment (SA) tags in BAM files to identify reads that have:
   - Primary alignment to mitochondrial DNA (chrM) and supplementary alignment to nuclear chromosomes, OR
   - Primary alignment to nuclear chromosomes and supplementary alignment to mitochondrial DNA

2. **Merging breakpoints**: Groups breakpoints by chromosome and merges consecutive breakpoints within a specified gap threshold (`max_gap_threshold`) to create intervals. This helps reduce noise and identify true NUMT insertion sites.

3. **Writing BND records**: Outputs breakend (BND) structural variant records in VCF format, representing the NUMT breakpoints. Each BND record includes:
   - Breakpoint positions on both the nuclear chromosome and mitochondrial genome
   - Strand orientation information
   - Supporting read counts (allele depth)
   - Properly formatted BND ALT fields based on strand orientations

**Usage:**
```bash
./target/release/Himito call-numts -i sample.bam \
                               -c chrM \
                               -m 10000 \
                               -r hg38.fa \
                               -o numts_breakpoints.vcf \
                               -s HG002 \
                               -a 2
```
**Output format:**
The output VCF file contains BND (breakend) records in standard VCFv4.3 format. Each record represents a NUMT breakpoint with:
- `CHROM`: Nuclear chromosome where the breakpoint occurs
- `POS`: Breakpoint position (1-based)
- `ALT`: BND format string indicating the connection to mitochondrial DNA
- `INFO`: Contains `SVTYPE=BND`
- `FORMAT`: GT (genotype) and AD (allele depth/supporting read count)

This will identify NUMT breakpoints where reads have alignments spanning both mitochondrial and nuclear genomes, merge breakpoints within 10kb, and output only those with at least 2 supporting reads. If you want to examine all possible NUMTs BND, adjust -a to <=1.

ToDo: may perform local assembly to find sequence resolved NUMTs insertions