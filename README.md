
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

### install and run Himito
```
git clone https://github.com/broadinstitute/Himito.git
cd Himito
cargo build --release

# filter NUMTs-derived reads
./target/release/Himito filter -i <input.bam> -c <chromosome in bam, e.g. "chrM"> -m <mt_output.bam> -n <numts_output.bam>

# construct graph
./target/release/Himito build -i <mt_output.bam> -k <kmer_size> -r <NC_012920.1.fasta> -o <output.gfa> 

# call variants from graph
./target/release/Himito call -g <output.gfa> -r <NC_012920.1.fasta> -k <kmer_size> -s <sampleid> -o <output.vcf>

# extract major haplotype from graph
./target/release/Himito asm -g <output.gfa>  -o <output.majorhaplotpe.fasta> -s <header string, e.g. "HG002 major haplotype">

# call methylation signals
./target/release/Himito methyl -g <output.annotated.gfa> -p <min_prob> -b <mt_test.bam> -o <methyl.bed>
```
