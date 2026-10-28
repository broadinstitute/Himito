
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
filter, linkage rescue and `<output>.matrix.csv`. `--input-bam` is accepted as an alias.
Without `-b`, the matrix comes from the graph and the strand filter and linkage rescue are off.
At indel sites the matrix genotyper scores each read against haplotypes that carry the sample's
majority alleles (VCF HF > 0.5) inside the site window, not bare rCRS. This matters where the sample
rewrites a homopolymer, e.g. HVS1 with 16182C/16183C/16189C, where `m.16181A>AC` lengthens a 12-C run.
SNP sites keep bare-reference haplotypes.
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

Given a Himito read–variant matrix (`<prefix>.matrix.csv` from `call`) and the corresponding VCF, the `lineage` command filters to informative heteroplasmic sites, deduplicates reads into haplotypes, and runs an MCMC search over mutation trees and identify a ML tree (algorithm adapted from SCITE, Kuipers J et al. Single-cell sequencing data reveals widespread recurrence and loss of mutational hits in the life histories of tumorus. Genome Research 2017; 27:1885-1894.). The initial tree is a Neighbor–Joining tree constructed from reads instead of a random tree though. Sequencing errors are modeled with false-positive and false-negative rates, measured per sample where the data allow:

- SNV false negatives: estimated from the sample's near-fixed SNV rows, with a floor of 0.005. If there are too few such rows, the `-d` preset is used (pacbio 0.05, ONT 0.02). The SNV false-positive rate is the `-d` preset (0.005). `--fp-rate` / `--fn-rate` override both.
- D310 and HVS1-C length alleles (alleles overlapping the tract, and indels whose run on the sample's majority alleles reaches its C run, such as `m.16181A>AC`): rates come from the sample's homopolymer length confusion, `<prefix>.hp_confusion.tsv`. `call --bam` writes this file, and `Himito hp-error` writes it for existing outputs. `lineage` uses `--hp-error`, or the file beside `<prefix>.matrix.csv`, or else an embedded HG00097 table.
  - chrM has no C run longer than 6 bp outside these tracts, so rates for longer runs are lower bounds.
  - With `--tract-alleles auto` (the default), every such allele is placed in the tree on every platform. An allele whose run is longer than the table calibrates is flagged `beyond_calibration` in `tract_annotation.tsv` rather than removed, so its placement is judged on its margin. That file lists every tract allele that reaches the tree search (`placed` = yes/no) with its best observed partner and whether it is significantly linked to it. A tract allele that went into the root block has no row there; `error_model.tsv` carries its run class, run length and `beyond_calibration` status.
  - Placed tract alleles form multi-allelic sites. The alleles at one tract position (e.g. `m.302A>ACC`/`m.302A>ACCC`) are one categorical character, scored with the sample's length-change kernel from `hp_confusion.tsv`, so a 1-bp slip costs one error term.
  - A length allele that occurs in several lineages (e.g. an HVS1 length independent of D310) is placed as recurrent leaf copies, named `<allele>#2`, `#3`, …. Each copy must raise the tree log-likelihood by more than ln(number of reads), and the reads it takes must be called that length significantly more often than the lineage's own length explains, assuming a slip rate of at least 5% (binomial, Bonferroni). A recurrent allele never has children, so it never defines a clade. Matrix outputs fold the copies into the allele's row.
  - Other `--tract-alleles` values: `site` is a synonym of `auto`, kept so the behaviour can be named explicitly, `search` gives each allele its own binary row (the earlier model), and `annotate` reports every tract allele in `tract_annotation.tsv` (`placed` = no) without placing it.
  - Any `--hp-fp-rate` / `--hp-fn-rate` given: no sites (binary rows with that rate); the alleles are still placed.
  - A site needs pure length alleles that agree on the reference run. Its root length is the commonest called length there that is not a site allele; with near-equal counts (e.g. HG00097 PacBio D310: C8 24 vs C11 23 reads) it rests on few reads.
  - `--hp-fp-rate` / `--hp-fn-rate` override the table.
- Every rate used, and where it came from, is written to `<prefix>.error_model.tsv`. 

`<prefix>.variant_disposition.tsv` accounts for every candidate variant that reached the matrix or is a PASS site inside the VCF's HF window: one row per in-window PASS site and per matrix row, with its VCF HF, matrix HF, read counts and an `outcome` of `informative`, `root_block`, `tract_annotated`, `dropped_hf_window`, `dropped_min_presence`, `dropped_min_absence` (only reachable with `--max-hf` above 1), `dropped_not_whitelisted` (the VCF does not list it as a PASS site inside the HF window; see `detail`) or `no_matrix_row` (a PASS site inside the HF window that `matrix.csv` has no row for), plus a `detail`. A non-PASS site, a PASS site with neither an HF nor an AF value, or a PASS site whose VCF HF falls outside the window produces no row at all when it has no matrix row, so the table is complete with respect to in-window PASS sites and matrix rows, not with respect to the VCF as a whole. A variant labelled `dropped_hf_window` or `dropped_not_whitelisted` whose matrix HF is at or above `--max-hf`, with at least 2 reads carrying it, may still appear in `cleaned_matrix_all_variants.csv` as a fixed homoplasmy (with `-v`, a site the VCF itself calls fixed is `dropped_not_whitelisted`, so that label is the common route and `dropped_hf_window` the rare one). Diffing two platforms' tables answers why their trees differ. `<prefix>.mutation_tree.tsv` accounts for every read: a trailing `UNASSIGNED` row counts the reads that call nothing at any informative variant. They are attached to no node and left out of the cleaned matrix, haplotype tips and Newick, rather than imputed as a wild-type clone.

Outputs include a mutation tree, read-level lineage Newick tree, the cleaned read–variant matrix `<prefix>.cleaned_matrix_all_variants.csv`, the variant co-occurrence table `variant_cooccurrence.tsv`, and haplotype maps. `variant_cooccurrence.tsv` counts observed calls only: each pair of informative variants is tabulated over the reads that call both sites, never over the tree-imputed matrix. `--debug` also writes the co-occurrence table over every filtered variant (`raw_variant_cooccurrence.tsv`, root-block variants included), the same table on the tree-imputed matrix (`cleaned_variant_cooccurrence.tsv`, which restates the tree), and `raw_haplotype_map.tsv`.

```
# only for outputs made before `call` wrote it:
./target/release/Himito hp-error -b <reads.bam> -r rCRS.fasta -o <prefix>

./target/release/Himito lineage -m <prefix>.matrix.csv \
                                 -v <output.vcf> \
                                 -o <output_prefix> \
                                 --min-hf 0.01 \
                                 --max-hf 0.95
```

Key outputs: `<output_prefix>.mutation_tree.tsv`, `<output_prefix>.read_lineage.nwk`, `<output_prefix>.cleaned_matrix_all_variants.csv`, `<output_prefix>.variant_cooccurrence.tsv`, `<output_prefix>.cleaned_haplotype_map.tsv`, `<output_prefix>.root_block.tsv`, `<output_prefix>.error_model.tsv`, `<output_prefix>.tract_annotation.tsv`, `<output_prefix>.variant_disposition.tsv`. `--debug` also writes `<output_prefix>.raw_haplotype_map.tsv`, `<output_prefix>.raw_variant_cooccurrence.tsv` and `<output_prefix>.cleaned_variant_cooccurrence.tsv`.

`lineage` shows `-m -v -d -o --min-hf --max-hf` in `-h`. The error-rate, `--hp-error`, `--tract-alleles`, MCMC and `--debug`
options are listed by `--help`.

To Do: 
1. add four gamete test
2. Add posteriors, convergence and confidence metrics
