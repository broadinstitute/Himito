
<p align="center">
<img src="https://github.com/user-attachments/assets/3129d44f-f69a-4491-a0cf-29fab5ee88a9" alt="logo" style="width: 30%; height: auto;">
</p>

# Himito

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21522816-blue?style=flat&logo=zenodo)](https://doi.org/10.5281/zenodo.21522816)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2025.11.03.686348-blue)](https://doi.org/10.1101/2025.11.03.686348)

Building Mitochondrial anchor-based graphical genome from long reads. Filter Numts reads, assemble major haplotypes, call homoplasmic and heteroplasmic variants, analyze methylation signals

## Documentation
[User Guide](https://github.com/broadinstitute/Himito/blob/main/docs/User_guide.md)

## Installation

### Standard install

**Rust:**
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

**System dependencies** (needed to compile `rust-htslib`):
```bash
# Debian/Ubuntu
sudo apt-get install -y build-essential cmake pkg-config libssl-dev libclang-dev zlib1g-dev
```

**Build Himito from source:**
```bash
git clone https://github.com/broadinstitute/Himito.git
cd Himito
cargo build --release
```

**Or download a prebuilt Linux binary** from [GitHub Releases](https://github.com/broadinstitute/Himito/releases):
```bash
wget https://github.com/broadinstitute/Himito/releases/download/v1.1.2/Himito
chmod +x Himito
./Himito --help
```

### Restricted network environments (incl. China)

If GitHub, Docker Hub, or `crates.io` are slow or unreachable, use the **Zenodo source archive**, Cargo/rustup mirrors, or an offline `cargo build --offline` workflow. See **[Installation in restricted network environments](docs/install_restricted_network.md)**.

## Usage
### Quick start (Filter-Build-Asm-Call-Methyl)
```
./target/release/Himito quick-start -i <lrWGS.bam> \
                                    -o <output_prefix> \
                                    -r <rCRS.fa> \
                                    -s <HG002> \
                                    -d pacbio

```
If you are analyzing high coverage samples, remember to set ```--maximal-mt-depth``` to a reasonable number. -d accepts: pacbio (default), ont-r9, ont-r10, ont-denoised

`-h` lists only the common options of `quick-start`, `call` and `lineage`. The tuning
options (for example `-k`, `-v`, `--p-value-threshold`, `--permutation-rounds`) still work
but are listed only by `--help`. The full help text of each command is kept in
[`docs/cli/`](docs/cli/).

`quick-start` and `call` left-align indels by default; `--no-left-align` keeps the graph's
indel spellings. See [call variants from graph](#call-variants-from-graph-default-for-pacbio) below.


## Details

### Filter NUMTs-derived reads
```
# filter NUMTs-derived reads
./target/release/Himito filter -i <input.bam> -c <chromosome in bam, e.g. "chrM"> -m <mt_output.bam> -n <numts_output.bam>
```

### Update: Denoise ONT reads before graph construction
For ONT R9 and R10 data, we recommend running Himito denoise to correct ONT reads in a BAM before graph construction. Denoising is a per-base model: it removes the *random* ONT error it can see, but leaves *systematic* error, which carries high base qualities, untouched. `ont-denoised` therefore keeps the ONT thresholds and the permutation co-occurrence test on SNPs (unlike `pacbio`, which exempts SNPs), and additionally makes the strand-bias filter **drop** strand-skewed low-VAF SNPs instead of only tagging them `FILTER=Strand_bias`. A call that is still strand-skewed after denoising has survived a filter designed to catch exactly that artifact, so it is discarded.

```
./target/release/Himito denoise -i <mt_output_bam> -r <NC_012920.1.fasta> -d <ont-r10/ont-r9> -o <denoised_bam> 
```

### build Himito graph
```
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
#### call variants from graph, default for pacbio.
```
./target/release/Himito call -g <output.gfa> -r <NC_012920.1.fasta> -s <sampleid> -d <data_type> -o <output.vcf> -b <raw_mt.bam>
```

`-b/--bam` takes the raw (pre-denoise) BAM of the same reads. It drives the strand-bias
filter, ONT linkage rescue and `<output>.matrix.csv`. `--input-bam` is accepted as an alias.
Without `-b`, the matrix comes from the graph and the strand filter and linkage rescue are off.
The strand-bias action, permutation mode and genotyping gates are set by the data type.

Indels are left-aligned against the reference by default, the way `bcftools norm`
does, so that equivalent spellings of one event in a homopolymer or tandem repeat
collapse into a single record with their supporting reads pooled. Pass `--no-left-align`
to keep the graph's indel spellings instead.

This changes the reported POS/REF/ALT, AC and HF of indels, and the variant names used
as row labels in the output matrix, so it also shifts the coordinates seen by downstream
lineage analysis. The VCF records it with a `##HimitoNormalization=left-aligned` header
line. An indel that would shift past the start of the reference stops there rather than
wrapping around the circular origin. Note that graphs built by `Himito build` already
carry left-anchored indel CIGARs in practice, so on such graphs left-alignment is a
guarantee rather than a transformation and typically changes nothing.
#### extract major haplotype from graph
```
./target/release/Himito asm -g <output.gfa>  -o <output.majorhaplotpe.fasta> -s <header string, e.g. "HG002 major haplotype">
```

#### call methylation signals
```
./target/release/Himito methyl -g <output.annotated.gfa> -p <min_prob> -b <mt_test.bam> -o <methyl.bed>
```
#### enumerate all possible haplotypes within windows (deprecated)
```
./target/release/Himito minorhap -g <output.gfa> -o <output.allhaplotype.fasta> -s <sample_id>
```

### July 14th 2026 Updates: Mitochondrial lineage inference (Under-construction)

Himito includes a new function for inferring mitochondrial mutation history from heteroplasmic variants given the Infinite Sites assumption (ISA) that every variant can only mutate once and will exist on the lineage over the lifetime. 

Given a Himito read–variant matrix (`<prefix>.matrix.csv` from `call`) and the corresponding VCF, the `lineage` command filters to informative heteroplasmic sites, deduplicates reads into haplotypes, and runs an MCMC search over mutation trees and identify a ML tree (algorithm adapted from SCITE, Kuipers J et al. Single-cell sequencing data reveals widespread recurrence and loss of mutational hits in the life histories of tumorus. Genome Research 2017; 27:1885-1894.). The initial tree is a Neighbor–Joining tree constructed from reads instead of a random tree though. Sequencing errors are modeled with false-positive and false-negative rates preset by `-d` (pacbio 0.005 / 0.05, ONT 0.005 / 0.02) and set with `--fp-rate` / `--fn-rate`. 

Outputs include a mutation tree, read-level lineage Newick tree, the cleaned read–variant matrix `<prefix>.cleaned_matrix_all_variants.csv`, the variant co-occurrence table `variant_cooccurrence.tsv`, and haplotype maps. The raw co-occurrence table and `raw_haplotype_map.tsv` are written only with `--debug`.

```
./target/release/Himito lineage -m <prefix>.matrix.csv \
                                 -v <output.vcf> \
                                 -o <output_prefix> \
                                 --min-hf 0.01 \
                                 --max-hf 0.95
```

Key outputs: `<output_prefix>.mutation_tree.tsv`, `<output_prefix>.read_lineage.nwk`, `<output_prefix>.cleaned_matrix_all_variants.csv`, `<output_prefix>.variant_cooccurrence.tsv`, `<output_prefix>.cleaned_haplotype_map.tsv`. `--debug` also writes `<output_prefix>.raw_haplotype_map.tsv` and `<output_prefix>.raw_variant_cooccurrence.tsv`.

`lineage` shows `-m -v -d -o --min-hf --max-hf` in `-h`; the error-rate, MCMC and `--debug`
options are listed by `--help`.

To Do: 
1. add four gamete test
2. Add posteriors, convergence and confidence metrics