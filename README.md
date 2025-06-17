
<p align="center">
<img src="https://github.com/user-attachments/assets/3129d44f-f69a-4491-a0cf-29fab5ee88a9" alt="logo" style="width: 30%; height: auto;">
</p>

# Himito 
Building Mitochondrial anchor-based graphical genome from long reads. Filter Numts reads, assemble major haplotypes, call homoplasmic and heteroplasmic variants, analyze methylation signals

## Documentation
[User Guide](https://github.com/broadinstitute/Himito/blob/main/docs/Userguid.md)

## Usage
### install rust
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### install and build Himito graph
```
git clone https://github.com/broadinstitute/Himito.git
cd Himito
cargo build --release

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
# call variants from graph
./target/release/Himito call -g <output.gfa> -r <NC_012920.1.fasta> -k <kmer_size> -s <sampleid> -o <output.vcf>

# extract major haplotype from graph
./target/release/Himito asm -g <output.gfa>  -o <output.majorhaplotpe.fasta> -s <header string, e.g. "HG002 major haplotype">

# call methylation signals
./target/release/Himito methyl -g <output.annotated.gfa> -p <min_prob> -b <mt_test.bam> -o <methyl.bed>
```
