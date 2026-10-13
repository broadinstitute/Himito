use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};
use anyhow::{Context, Result as AnyhowResult};
use log::{info};

mod agg;
mod asm;
mod build;
mod call;
mod filter;
mod methyl;
mod correct;
mod minorhap;
mod callnumts;
mod lineage;
mod rootblock;
mod tracts;
mod scite;
mod denoise;
mod denoise_indel;
mod linkage;

// Thresholds: see dev_docs/superpowers/plans/2026-10-11-ont-recall-high-coverage.md
pub const DENOISE_KEEP_VAF: f64 = 0.03;
pub const PENDING_FLOOR: f64 = 0.005;
pub const PERM_STOP_AFTER: u32 = 10;
pub const PERM_HEADROOM: f64 = 5.0;
pub const PERM_MAX_ROUNDS: usize = 100_000;
pub const PERM_SEED: u64 = 42;
pub const STRAND_BIAS_MIN_EFFECT: f64 = 0.20;
pub const STRAND_LOCAL_MIN_READS: usize = 20;
pub const GENOTYPE_MAX_ED_ONT: u32 = 4;
pub const GENOTYPE_FLANK: usize = 10;
pub const COOC_EXEMPT_HF: f64 = 0.10;
pub const COOC_EXEMPT_MIN_AC: usize = 20;
pub const LINK_MIN_HF: f64 = 0.004;
pub const LINK_MAX_HF: f64 = 0.20;
pub const LINK_MIN_AC: u32 = 10;
pub const LINK_MIN_DIST: u32 = 20;
pub const LINK_ERR_QUANTILE: f64 = 0.95;
pub const LINK_MIN_SHARED: u32 = 5;
pub const LINK_MIN_SHARED_STRAND: u32 = 3;
pub const LINK_MIN_COCOVER: u32 = 50;
pub const LINK_MIN_COCOVER_STRAND: u32 = 30;
pub const LINK_MIN_ENRICH: f64 = 10.0;
pub const LINK_LOG10P: f64 = -10.0;
pub const LINK_STRAND_LOG10P: f64 = -5.0;
pub const DENOISE_MIN_STRAND: u32 = 2;
pub const DENOISE_STRAND_BIAS_P: f64 = 0.01;
pub const DENOISE_HOMOPLASMIC_VAF: f64 = 0.95;

/// The data types Himito accepts wherever a `-d/--data-type` is taken.
pub const DATA_TYPES: [&str; 4] = ["pacbio", "ont-r9", "ont-r10", "ont-denoised"];

/// Minimal reads a variant must be present in / absent from for `lineage` to treat it as informative.
pub const LINEAGE_MIN_PRESENCE: usize = 2;
pub const LINEAGE_MIN_ABSENCE: usize = 1;

/// Minimal reads backing a haplotype for `lineage` to report it.
pub const LINEAGE_MIN_READS: usize = 3;

/// RNG seed for the SCITE MCMC search, fixed so runs are reproducible.
pub const LINEAGE_MCMC_SEED: u64 = 42;

/// Lineage root block (always on): minimum HF for a variant to be CONSIDERED homoplasmic.
pub const LINEAGE_ROOT_BLOCK_MIN_HF: f64 = 0.80;

/// Lineage root block: Benjamini-Hochberg q above which a variant's absences count as unstructured.
pub const LINEAGE_ROOT_BLOCK_MAX_Q: f64 = 0.05;

/// Lineage root block: below this many absent calls a high-HF variant is admitted without testing.
pub const LINEAGE_ROOT_BLOCK_MIN_ABSENT: usize = 10;

/// Lineage root block: minimum enrichment odds ratio for a partner's absences to be structured.
pub const LINEAGE_ROOT_BLOCK_MIN_OR: f64 = 2.0;

/// Quick-start: longest read span (bp) the graph build aligns.
pub const QUICK_START_LENGTH_MAX: usize = 3000;

/// Quick-start: RNG seed of the mt-read downsampling, fixed for reproducible subsets.
pub const QUICK_START_DOWNSAMPLE_SEED: u64 = 42;

/// Strand-skewed SNPs below this HF are dropped instead of flagged.
pub const STRAND_BIAS_DROP_MAX_HF: f64 = 0.10;

#[derive(Debug, Parser)]
#[clap(name = "Himito")]
#[clap(version = env!("CARGO_PKG_VERSION"))]
#[clap(about = "Analysis of mitochondrial genome using long reads.", long_about = None)]

struct Cli {
    #[command(flatten)]
    global: GlobalOpts,

    #[clap(subcommand)]
    command: Commands,
}

#[derive(Args, Debug)]
pub struct GlobalOpts {
    /// worker threads for parallel stages
    #[arg(long, global = true)]
    pub threads: Option<usize>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Quick start pipeline for mitochondrial assembly.
    #[clap(arg_required_else_help = true)]
    QuickStart {
        /// input BAM of reads, coordinate-sorted
        #[clap(short, long, value_parser, required = true)]
        input_bam: PathBuf,

        /// mitochondrial contig name in the BAM
        #[clap(short, long, value_parser, default_value = "chrM")]
        chromo: String,

        /// reference FASTA (rCRS)
        #[clap(short, long, value_parser, required = true)]
        reference_path: PathBuf,

        /// sample name written to the outputs
        #[clap(short, long, value_parser, required = true)]
        sample_id: String,

        /// sequencing data type; sets the filter and denoise presets
        #[clap(short, long, default_value = "pacbio", value_parser = clap::builder::PossibleValuesParser::new(DATA_TYPES))]
        data_type: String,

        /// output prefix for the BAM, GFA, VCF, matrix and BED outputs
        #[clap(short, long, value_parser, required = true)]
        output_prefix: PathBuf,

        /// maximum mtDNA reads kept after downsampling
        #[clap(long, value_parser, default_value_t = 50000)]
        maximal_mt_depth: usize,

        /// minimum probability to call a C methylated
        #[clap(short, long, value_parser, default_value_t = 0.5, hide_short_help = true)]
        threshold_methylation_prob: f64,

        /// maximum methylated fraction for a read to count as mtDNA
        #[clap(short, long, value_parser, default_value_t = 0.2, hide_short_help = true)]
        filter_max_methylation: f64,

        /// k-mer size of the anchor graph
        #[clap(short, long, value_parser, default_value_t = 21, hide_short_help = true)]
        kmer_size: usize,

        /// minimum allele count for a variant
        #[clap(short, long, value_parser, default_value_t = 2, hide_short_help = true)]
        minimal_ac: usize,

        /// minimum heteroplasmic frequency for a variant
        #[clap(short, long, value_parser, default_value_t = 0.01, hide_short_help = true)]
        vaf_threshold: f32,

        /// permutation-test p-value cutoff; 1 disables it (preset by data type)
        #[clap(long, value_parser, hide_short_help = true)]
        p_value_threshold: Option<f64>,

        /// permutation-test frequency threshold (preset by data type)
        #[clap(long, value_parser, hide_short_help = true)]
        heteroplasmic_frequency_threshold: Option<f64>,

        /// p-value cutoff of the strand-bias filter
        #[clap(long, value_parser, default_value_t = 0.05, hide_short_help = true)]
        strand_bias_threshold: f64,

        /// indels below this HF are flagged Potential_Artifact (preset by data type)
        #[clap(long, value_parser, hide_short_help = true)]
        indel_false_threshold: Option<f64>,

        /// variants above this matrix frequency skip the permutation test (preset)
        #[clap(long, value_parser, hide_short_help = true)]
        permutation_frequency_threshold: Option<f64>,

        /// minimum permutation rounds per tested variant
        #[clap(long, value_parser, default_value_t = 1000, hide_short_help = true)]
        permutation_rounds: usize,

        // `--left-align` is a no-op kept so old command lines still parse.
        /// no-op; indels are left-aligned by default
        #[clap(long = "left-align", action, hide = true)]
        left_align_compat: bool,

        /// do not left-align indels (reproduces pre-2026-09 variant names)
        #[clap(long = "no-left-align", action, conflicts_with = "left_align_compat", hide_short_help = true)]
        no_left_align: bool,

        /// turn off read-linkage rescue of low-frequency SNVs (ONT only)
        #[clap(long = "no-linkage-rescue", action, hide_short_help = true)]
        no_linkage_rescue: bool,

        /// keep threshold of the denoise site model (ont-denoised only)
        #[clap(long, value_parser, default_value_t = DENOISE_KEEP_VAF, hide_short_help = true)]
        denoise_keep_vaf: f64,
    },

    /// Filter reads derived from Numts
    #[clap(arg_required_else_help = true)]
    Filter {
        /// input path for bam file.
        #[clap(short, long, value_parser, required = true)]
        input_bam: PathBuf,

        /// contig name in the bam file
        #[clap(short, long, value_parser)]
        chromo: String,

        /// output path for mtDNA bam file
        #[clap(short, long, value_parser, required = true)]
        mt_output: PathBuf,

        /// output path for numts bam file
        #[clap(short, long, required = true, value_parser)]
        numts_output: PathBuf,

        /// min_probability to determine a C is methylated
        #[clap(short, long, value_parser, default_value_t = 0.5)]
        threshold_methylation_prob: f64,

        /// max fraction to keep a read as mtDNA derived
         #[clap(short, long, value_parser, default_value_t = 0.2)]
        fraction_max_methylation: f64,

        /// maximal reads to keep from mtDNA
        #[clap(long, value_parser, default_value_t = 50000)]
        maximal_mt_depth: usize,

        /// RNG seed for mt read downsampling (reproducible subsampling)
        #[clap(long, value_parser, default_value_t = 42)]
        downsample_seed: u64,

    },

    /// Build graph from long-read data in FASTA or Bam file.
    #[clap(arg_required_else_help = true)]
    Build {
        /// bam or fasta file with reads spanning locus of interest.
        #[clap(short, long, value_parser,required = true)]
        input_read_path: PathBuf,

        /// Reference Fasta file, rCRS
        #[clap(short, long, value_parser, required = true)]
        reference_path: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = 21)]
        kmer_size: usize,

        /// Output path for anchor graph.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// max Length to do alignment
        #[clap(short, long, value_parser, default_value_t = 3000)]
        length_max: usize,

        /// minimal reads on a graph edge before it is aligned into a CIGAR.
        #[clap(long, value_parser, default_value_t = 2)]
        min_edge_reads: usize,

    },

    /// Denoise ONT reads in a BAM (SNV error correction) before graph construction.
    #[clap(arg_required_else_help = true)]
    Denoise {
        /// input BAM (coordinate-sorted)
        #[clap(short, long, value_parser, required = true)]
        input: PathBuf,

        /// output denoised BAM
        #[clap(short, long, value_parser, required = true)]
        output: PathBuf,

        /// reference FASTA (rCRS)
        #[clap(short, long, value_parser, required = true)]
        reference: PathBuf,

        /// data type; only ont-* is denoised, others pass through unchanged
        #[clap(short, long, value_parser, default_value = "ont-r10")]
        data_type: String,

        /// minimal VAF for an allele to be kept (below this it is corrected away)
        #[clap(long, value_parser, default_value_t = DENOISE_KEEP_VAF)]
        vaf: f64,

        /// enable small-indel (<5bp) correction in addition to substitutions
        #[clap(long, value_parser, default_value_t = false)]
        indels: bool,

        /// an allele at or above this fraction of reads at a site is never corrected away.
        #[clap(long, value_parser, default_value_t = 0.2)]
        indel_protect_vaf: f64,

        /// EXPERIMENTAL: treat a homopolymer as a single LENGTH.
        #[clap(long, value_parser, default_value_t = false)]
        indel_hp_length_model: bool,

        /// optional path to write denoise statistics as JSON
        #[clap(long, value_parser)]
        stats: Option<PathBuf>,
    },

    /// Correct graph based on srWGS data
    #[clap(arg_required_else_help = true)]
    Correct {
        /// path for anchor graph.
        #[clap(short, long, value_parser)]
        graphfile: PathBuf,
        /// path for bwt file of srWGS reads, fasta or fasta.gz file
        #[clap(short, long, value_parser)]
        bwt_file: String,
        /// path for corrected graph gfa file
        #[clap(short, long, value_parser)]
        outputfile: PathBuf,
        /// path for standard linear reference FASTA file
        #[clap(short, long, value_parser)]
        reference_file: PathBuf,
        /// query length for kmerizing graph, should be less than the short read length
        #[clap(short, long, value_parser, default_value_t = 99)]
        query_length: usize,
        /// min support counts for a read to be considered
        #[clap(short, long, value_parser, default_value_t = 10)]
        min_support_counts: usize,
    },

    /// Call Variants from Sequence Graph
    #[clap(arg_required_else_help = true)]
    Call {
        /// anchor graph (GFA)
        #[clap(short, long, value_parser)]
        graphfile: PathBuf,

        /// reference FASTA (rCRS)
        #[clap(short, long, value_parser, required = true)]
        reference_fasta: PathBuf,

        /// sample name written to the VCF
        #[clap(short, long, value_parser, required = true)]
        sample_id: String,

        /// sequencing data type; sets the filter presets
        #[clap(short, long, default_value = "pacbio", value_parser = clap::builder::PossibleValuesParser::new(DATA_TYPES))]
        data_type: String,

        /// output VCF path
        #[clap(short, long, value_parser, required = true)]
        output_file: PathBuf,

        /// raw pre-denoise BAM: strand filter, linkage rescue, matrix.csv
        #[clap(short = 'b', long = "bam", alias = "input-bam", value_parser)]
        bam: Option<PathBuf>,

        /// k-mer size of the anchor graph
        #[clap(short, long, value_parser, default_value_t = 21, hide_short_help = true)]
        k: usize,

        /// minimum allele count for a variant
        #[clap(short, long, value_parser, default_value_t = 2, hide_short_help = true)]
        minimal_ac: usize,

        /// minimum heteroplasmic frequency for a variant
        #[clap(short, long, value_parser, default_value_t = 0.01, hide_short_help = true)]
        vaf_threshold: f32,

        /// permutation-test p-value cutoff; 1 disables it (preset by data type)
        #[clap(short, long, value_parser, hide_short_help = true)]
        p_value_threshold: Option<f64>,

        /// permutation-test frequency threshold (preset by data type)
        #[clap(short, long, value_parser, hide_short_help = true)]
        frequency_threshold: Option<f64>,

        /// p-value cutoff of the strand-bias filter (only used with --bam)
        #[clap(long, value_parser, default_value_t = 0.01, hide_short_help = true)]
        strand_bias_threshold: f64,

        /// indels below this HF are flagged Potential_Artifact (preset by data type)
        #[clap(long, value_parser, hide_short_help = true)]
        indel_false_threshold: Option<f64>,

        /// variants above this matrix frequency skip the permutation test (preset)
        #[clap(long, value_parser, hide_short_help = true)]
        permutation_frequency_threshold: Option<f64>,

        /// minimum permutation rounds per tested variant
        #[clap(long, value_parser, default_value_t = 1000, hide_short_help = true)]
        permutation_rounds: usize,

        // `--left-align` is a no-op kept so old command lines still parse.
        /// no-op; indels are left-aligned by default
        #[clap(long = "left-align", action, hide = true)]
        left_align_compat: bool,

        /// do not left-align indels (reproduces pre-2026-09 variant names)
        #[clap(long = "no-left-align", action, conflicts_with = "left_align_compat", hide_short_help = true)]
        no_left_align: bool,
    },

    /// Extract Major Haplotype as Fasta file from Graph
    #[clap(arg_required_else_help = true)]
    Asm {
        /// path for anchor graph.
        #[clap(short, long, value_parser)]
        graphfile: PathBuf,
        /// path for output fasta file
        #[clap(short, long, value_parser)]
        outputfile: PathBuf,

        /// header for the major haplotype, usually the sample name
        #[clap(short, long, value_parser)]
        sample: String,

    },

    /// Annotate Methylation signals to the Graph
    #[clap(arg_required_else_help = true)]
    Methyl {
        /// path for cigar annotated anchor graph.
        #[clap(short, long, value_parser)]
        graphfile: PathBuf,
        /// path for bam file with MM/ML tags.
        #[clap(short, long, value_parser)]
        bamfile: PathBuf,
        /// path for output methylation bed file
        #[clap(short, long, value_parser)]
        outputfile: PathBuf,
        /// min_probability to determine a C is methylated
        #[clap(short, long, value_parser, default_value_t = 0.5)]
        prob_min: f64,
        /// extract per-read methylation signals for the major haplotype or all reads.
        #[clap(short, long, value_parser, default_value_t = false)]
        major_haplotype: bool
    },

    /// Extract Minor Haplotype as Fasta file from Graph
    #[clap(arg_required_else_help = true)]
    Minorhap {
        /// path for anchor graph.
        #[clap(short, long, value_parser)]
        graphfile: PathBuf,
        /// length of the linear reference genome, e.g.rCRS length is 16569 (default)
        #[clap(short, long, value_parser, default_value_t = 16569)]
        ref_length: i64,
        /// bin size for the each haplotype window, default is 1000bp
        #[clap(short, long, value_parser, default_value_t = 1000)]
        bin_size: i32,
        /// pad size for the each haplotype window, default is 100bp
        #[clap(short, long, value_parser, default_value_t = 100)]
        pad_size: i64,

        /// minimal read ratio supporting each haplotype in a given window.
        #[clap(short, long, value_parser, default_value_t = 0.01)]
        min_read_ratio: f64,

        /// minimal read count supporting each haplotype in a given window, default is 1
        #[clap(short, long, value_parser, default_value_t = 1)]
        count_support: usize,

        /// path for output fasta file
        #[clap(short, long, value_parser)]
        outputfile: PathBuf,

        /// header for the major haplotype, usually the sample name
        #[clap(short, long, value_parser)]
        sample: String,
    },

    /// Call Numts from BAM file
    #[clap(arg_required_else_help = true)]
    CallNumts {
        /// input path for bam file.
        #[clap(short, long, value_parser, required = true)]
        input_bam: PathBuf,
        /// mitochondrial contig name in the bam file
        #[clap(short, long, value_parser)]
        chromo: String,
        /// output path for numts vcf file
        #[clap(short, long, value_parser, required = true)]
        output_vcf: PathBuf,
        /// path for standard linear reference FASTA file
        #[clap(short, long, value_parser, required = true)]
        reference_file: PathBuf,
        /// sample name of the bam file
        #[clap(short, long, value_parser, required = true)]
        sample_name: String,
        /// junction-signature clustering tolerance (bp)
        #[clap(long, value_parser, default_value_t = 100)]
        tolerance: i64,
        /// minimum MAPQ for both aligned segments of a junction
        #[clap(long, value_parser, default_value_t = 20)]
        min_mapq: u8,
        /// minimum aligned segment length (bp)
        #[clap(long, value_parser, default_value_t = 100)]
        min_seg_len: i64,
        /// minimal supporting-read count (AC) for a PASS variant
        #[clap(short, long, value_parser, default_value_t = 2)]
        ac_threshold: usize,
        /// base allowance for supplementary alignments per read (chimera filter)
        #[clap(long, value_parser, default_value_t = 3)]
        max_splits: i64,
        /// extra split allowance per read-kb (chimera filter)
        #[clap(long, value_parser, default_value_t = 0.1)]
        max_splits_per_kb: f64,
        /// optional BED of reference NUMTs; overlapping calls get a REFNUMT filter
        #[clap(long, value_parser)]
        ref_numt_bed: Option<PathBuf>,
        /// emit breakend (BND) records for all calls instead of sequence-resolved INS
        #[clap(long, action)]
        emit_bnd: bool,
    },

    /// Infer mitochondrial lineage tree from Himito read_var matrix CSV.
    #[clap(arg_required_else_help = true)]
    Lineage {
        /// Himito matrix CSV (<prefix>.matrix.csv)
        #[clap(short, long, value_parser, required = true)]
        matrix_file: PathBuf,

        /// Himito VCF; restricts the analysis to its heteroplasmic sites
        #[clap(short, long, value_parser)]
        vcf_file: Option<PathBuf>,

        /// sequencing data type; sets the SCITE error-rate preset
        #[clap(
            short,
            long,
            default_value = "pacbio",
            value_parser = clap::builder::PossibleValuesParser::new(DATA_TYPES)
        )]
        data_type: String,

        /// output prefix; writes <prefix>.haplotype_map.tsv and the trees
        #[clap(short, long, value_parser, required = true)]
        output_prefix: String,

        /// minimum heteroplasmic frequency (inclusive) of a variant
        #[clap(long, value_parser, default_value_t = 0.01)]
        min_hf: f64,

        /// maximum heteroplasmic frequency (exclusive); excludes fixed variants
        #[clap(long, value_parser, default_value_t = 0.95)]
        max_hf: f64,

        /// SCITE false-positive rate (preset by data type)
        #[clap(long, value_parser, hide_short_help = true)]
        fp_rate: Option<f64>,

        /// SCITE false-negative rate (preset by data type)
        #[clap(long, value_parser, hide_short_help = true)]
        fn_rate: Option<f64>,

        /// SCITE false-positive rate in HP tracts (twice --fp-rate, capped at 0.5)
        #[clap(long, value_parser, hide_short_help = true)]
        hp_fp_rate: Option<f64>,

        /// SCITE false-negative rate in HP tracts (twice --fn-rate, capped at 0.5)
        #[clap(long, value_parser, hide_short_help = true)]
        hp_fn_rate: Option<f64>,

        /// merge equivalent spellings inside known rCRS homopolymer tracts
        // Tracts: D310 302-315, HVS1-C 16183-16193, CA514 513-525.
        #[clap(long, action, hide_short_help = true)]
        collapse_hp: bool,

        /// MCMC iterations per chain
        #[clap(long, value_parser, default_value_t = 5000, hide_short_help = true)]
        mcmc_iterations: usize,

        /// independent MCMC chains; the best tree is kept
        #[clap(long, value_parser, default_value_t = 3, hide_short_help = true)]
        mcmc_chains: usize,

        /// also write raw_haplotype_map.tsv and raw_variant_cooccurrence.tsv
        #[clap(long, action, hide_short_help = true)]
        debug: bool,
    }
}

/// Validate the user-settable `--indel-*` options at the CLI boundary.
fn validate_indel_opts(o: &denoise_indel::IndelOpts) -> AnyhowResult<()> {
    anyhow::ensure!(
        o.protect_vaf > 0.0 && o.protect_vaf <= 1.0,
        "--indel-protect-vaf must be in (0, 1] (got {})", o.protect_vaf
    );
    Ok(())
}

fn validate_lineage_opts(mcmc_chains: usize) -> AnyhowResult<()> {
    if mcmc_chains == 0 {
        anyhow::bail!("--mcmc-chains must be at least 1 (got {mcmc_chains})");
    }
    Ok(())
}

pub fn init_rayon_threads(threads: Option<usize>) -> AnyhowResult<()> {
    let Some(n) = threads else {
        return Ok(());
    };
    anyhow::ensure!(n >= 1, "--threads must be at least 1 (got {n})");
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
        .context("failed to initialize rayon thread pool (already initialized?)")?;
    info!("Rayon thread pool: {n} worker(s)");
    Ok(())
}

fn main() {
    env_logger::init();
    let args = Cli::parse();
    init_rayon_threads(args.global.threads);
    match args.command {
        Commands::QuickStart {
            input_bam,
            chromo,
            sample_id,
            reference_path,
            data_type,
            output_prefix,
            threshold_methylation_prob,
            filter_max_methylation,
            kmer_size,
            minimal_ac,
            vaf_threshold,
            p_value_threshold,
            heteroplasmic_frequency_threshold,
            maximal_mt_depth,
            strand_bias_threshold,
            indel_false_threshold,
            permutation_frequency_threshold,
            permutation_rounds,
            no_left_align,
            left_align_compat: _,
            no_linkage_rescue,
            denoise_keep_vaf,
        } => {
            let left_align = !no_left_align;
            let (p_value_threshold, frequency_threshold, permutation_frequency_threshold_) =
                call::resolve_thresholds(&data_type, p_value_threshold, heteroplasmic_frequency_threshold, permutation_frequency_threshold);
            let indel_false_threshold = call::resolve_indel_false_threshold(&data_type, indel_false_threshold);
            let mt_output = output_prefix.with_extension("mt.bam");
            let numts_output = output_prefix.with_extension("numts.bam");
            let _ = filter::start(
                &input_bam,
                &chromo,
                &mt_output,
                &numts_output,
                threshold_methylation_prob,
                filter_max_methylation,
                maximal_mt_depth,
                QUICK_START_DOWNSAMPLE_SEED,
            );
            
            // For ONT, denoise the mt BAM so reads consolidate; PacBio is untouched.
            let (build_bam, data_type_) = if data_type == "ont-denoised" {
                let denoised = output_prefix.with_extension("mt.denoised.bam");
                // Denoise settings for the combined pipeline.
                if let Err(e) = denoise::start(
                    &mt_output, &denoised, &reference_path, &data_type,
                    // Denoise's keep threshold, NOT `vaf_threshold`.
                    denoise_keep_vaf, 0.7,
                    // Indel correction is on in QuickStart.
                    &denoise_indel::IndelOpts {
                        enabled: true,
                        hp_length_model: false,
                        ..Default::default()
                    },
                    None,
                ) {
                    eprintln!("Warning: denoise failed ({e:#}); falling back to raw mt BAM.");
                    (mt_output.clone(), data_type.as_str())
                } else {
                    (denoised, "ont-denoised")
                }
            } else {
                (mt_output.clone(), data_type.as_str())
            };
            info!("Data Type for variant calling: {data_type_}");

            let graph_output = output_prefix.with_extension("gfa");
            let _ = build::start(&graph_output, kmer_size, &build_bam, &reference_path, QUICK_START_LENGTH_MAX, 2);

            let assemble_output = output_prefix.with_extension("fasta");
            asm::start(&graph_output, &assemble_output, &sample_id);

            let vcf_output = output_prefix.with_extension("vcf");
            println!("{}", data_type_);
            call::start(
                &graph_output,
                &reference_path,
                kmer_size,
                minimal_ac,
                &vcf_output,
                &sample_id,
                vaf_threshold,
                &data_type_,
                p_value_threshold,
                frequency_threshold,
                // Raw mt BAM, not build_bam (denoise turns uncertain ALTs into REF).
                Some(&mt_output),
                strand_bias_threshold,
                indel_false_threshold,
                permutation_frequency_threshold_,
                permutation_rounds,
                // Adaptive permutation test; Fixed is only reachable from unit tests.
                call::PermutationMode::Adaptive,
                // None: the data-type preset (drop for ont-denoised, flag otherwise).
                None,
                STRAND_BIAS_DROP_MAX_HF,
                left_align,
                no_linkage_rescue,
            );
            let annotated_graph_output = output_prefix.with_extension("gfa");
            let methyl_output = output_prefix.with_extension("bed");
            methyl::start(&annotated_graph_output, &mt_output, &methyl_output, threshold_methylation_prob, false);
        }
        Commands::Filter {
            input_bam,
            chromo,
            mt_output,
            numts_output,
            threshold_methylation_prob,
            fraction_max_methylation,
            maximal_mt_depth,
            downsample_seed,
        } => {
            let _ = filter::start(
                &input_bam,
                &chromo,
                &mt_output,
                &numts_output,
                threshold_methylation_prob,
                fraction_max_methylation,
                maximal_mt_depth,
                downsample_seed,
            );
        }

        Commands::Build {
            output,
            kmer_size,
            input_read_path,
            reference_path,
            length_max,
            min_edge_reads,
        } => {
            build::start(&output, kmer_size, &input_read_path, &reference_path, length_max, min_edge_reads);
        }

        Commands::Denoise {
            input,
            output,
            reference,
            data_type,
            vaf,
            indels,
            indel_protect_vaf,
            indel_hp_length_model,
            stats,
        } => {
            let iopts = denoise_indel::IndelOpts {
                enabled: indels,
                protect_vaf: indel_protect_vaf,
                hp_length_model: indel_hp_length_model,
                ..Default::default()
            };
            // Only validate when --indels is actually in effect.
            if iopts.enabled {
                if let Err(e) = validate_indel_opts(&iopts) {
                    eprintln!("Error: {e:#}");
                    std::process::exit(1);
                }
            }
            if let Err(e) = denoise::start(
                &input, &output, &reference, &data_type, vaf,
                DENOISE_HOMOPLASMIC_VAF, &iopts, stats.as_ref(),
            ) {
                eprintln!("Error running denoise: {e:#}");
                std::process::exit(1);
            }
        }

        Commands::Correct {
            graphfile,
            bwt_file,
            outputfile,
            reference_file,
            query_length,
            min_support_counts
        } => {
            let _ = correct::start(&graphfile, &bwt_file, &reference_file, &outputfile, query_length, min_support_counts);
        }

        Commands::Call {
            graphfile,
            reference_fasta,
            k,
            minimal_ac,
            output_file,
            sample_id,
            vaf_threshold,
            data_type,
            p_value_threshold,
            frequency_threshold,
            bam,
            strand_bias_threshold,
            indel_false_threshold,
            permutation_frequency_threshold,
            permutation_rounds,
            no_left_align,
            left_align_compat: _,
        } => {
            let left_align = !no_left_align;
            let (p_value_threshold, frequency_threshold, permutation_frequency_threshold) =
                call::resolve_thresholds(&data_type, p_value_threshold, frequency_threshold, permutation_frequency_threshold);
            let indel_false_threshold = call::resolve_indel_false_threshold(&data_type, indel_false_threshold);
            call::start(
                &graphfile,
                &reference_fasta,
                k,
                minimal_ac,
                &output_file,
                &sample_id,
                vaf_threshold,
                &data_type,
                p_value_threshold,
                frequency_threshold,
                bam.as_deref(),
                strand_bias_threshold,
                indel_false_threshold,
                permutation_frequency_threshold,
                permutation_rounds,
                // Adaptive permutation test; Fixed is only reachable from unit tests.
                call::PermutationMode::Adaptive,
                // None: the data-type preset (drop for ont-denoised, flag otherwise).
                None,
                STRAND_BIAS_DROP_MAX_HF,
                left_align,
                false,
            );
        }

        Commands::Asm {
            graphfile,
            outputfile,
            sample,
        } => {
            asm::start(&graphfile, &outputfile, &sample);
        }

        Commands::Methyl {
            graphfile,
            bamfile,
            outputfile,
            prob_min,
            major_haplotype
        } => {
            methyl::start(&graphfile, &bamfile, &outputfile, prob_min, major_haplotype);
        }
        Commands::Minorhap {
                graphfile,
                ref_length,
                bin_size,
                pad_size,
                min_read_ratio,
                count_support,
                outputfile,
                sample,
            } => {
                minorhap::start(&graphfile, ref_length, bin_size,  pad_size, min_read_ratio, count_support,&outputfile, &sample);
        }
        Commands::CallNumts {
            input_bam,
            chromo,
            output_vcf,
            reference_file,
            sample_name,
            tolerance,
            min_mapq,
            min_seg_len,
            ac_threshold,
            max_splits,
            max_splits_per_kb,
            ref_numt_bed,
            emit_bnd,
        } => {
            let config = callnumts::CallNumtsConfig {
                tolerance,
                min_mapq,
                min_seg_len,
                ac_threshold,
                max_splits,
                max_splits_per_kb,
                emit_bnd,
                ref_numt_bed,
            };
            if let Err(e) = callnumts::start(
                &input_bam,
                &chromo,
                &output_vcf,
                &reference_file,
                &sample_name,
                &config,
            ) {
                eprintln!("Error calling NUMTs: {:#}", e);
                std::process::exit(2);
            }
        }
        Commands::Lineage {
            matrix_file,
            vcf_file,
            min_hf,
            max_hf,
            data_type,
            fp_rate,
            fn_rate,
            hp_fp_rate,
            hp_fn_rate,
            collapse_hp,
            mcmc_iterations,
            mcmc_chains,
            debug,
            output_prefix,
        } => {
            let matrix_file = matrix_file.to_str().expect("matrix-file path is not valid UTF-8");
            let vcf_file = vcf_file
                .as_deref()
                .map(|p| p.to_str().expect("vcf-file path is not valid UTF-8"));
            let (fp_rate, fn_rate) =
                lineage::resolve_error_rates(&data_type, fp_rate, fn_rate);
            let hp_rates = lineage::resolve_hp_rates(fp_rate, fn_rate, hp_fp_rate, hp_fn_rate);
            // Empty slice: collapse is the identity.
            let tracts = if collapse_hp {
                crate::tracts::KNOWN_TRACTS
            } else {
                &[]
            };
            // Root block is always on.
            let root_block = Some(rootblock::RootBlockConfig {
                min_hf: LINEAGE_ROOT_BLOCK_MIN_HF,
                max_q: LINEAGE_ROOT_BLOCK_MAX_Q,
                min_absent: LINEAGE_ROOT_BLOCK_MIN_ABSENT,
                min_or: LINEAGE_ROOT_BLOCK_MIN_OR,
            });
            // `run_mcmc_multichain` panics on an empty chain list.
            if let Err(e) = validate_lineage_opts(mcmc_chains) {
                eprintln!("Error: {e:#}");
                std::process::exit(1);
            }
            if let Err(e) = lineage::start(
                matrix_file,
                vcf_file,
                min_hf,
                max_hf,
                LINEAGE_MIN_PRESENCE,
                LINEAGE_MIN_ABSENCE,
                LINEAGE_MIN_READS,
                fp_rate,
                fn_rate,
                hp_rates,
                tracts,
                mcmc_iterations,
                mcmc_chains,
                LINEAGE_MCMC_SEED,
                root_block,
                debug,
                &output_prefix,
            ) {
                eprintln!("Error running lineage analysis: {:#}", e);
                std::process::exit(1);
            }
        }
    }


}

#[cfg(test)]
mod indel_opts_validation_tests {
    use super::*;
    use denoise_indel::IndelOpts;

    fn opts() -> IndelOpts {
        IndelOpts { enabled: true, ..Default::default() }
    }

    #[test]
    fn shipped_defaults_pass_validation() {
        assert!(validate_indel_opts(&opts()).is_ok());
    }

    #[test]
    fn protect_vaf_zero_is_rejected() {
        // --indel-protect-vaf 0 protects every non-strand-rejected read, so it is the exact inverse of "protect less.
        let mut o = opts();
        o.protect_vaf = 0.0;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn protect_vaf_one_is_accepted() {
        let mut o = opts();
        o.protect_vaf = 1.0;
        assert!(validate_indel_opts(&o).is_ok());
    }

    #[test]
    fn protect_vaf_above_one_is_rejected() {
        let mut o = opts();
        o.protect_vaf = 1.5;
        assert!(validate_indel_opts(&o).is_err());
    }

    // The former per-field range tests were deleted along with the flags that set those fields.

    #[test]
    fn hp_length_model_is_off_unless_the_flag_is_passed() {
        use clap::Parser;
        let base = [
            "Himito", "denoise", "-i", "in.bam", "-o", "out.bam", "-r", "ref.fa", "--indels",
        ];
        let hp_of = |args: &[&str]| match Cli::try_parse_from(args).unwrap().command {
            Commands::Denoise { indel_hp_length_model, .. } => indel_hp_length_model,
            other => panic!("expected the Denoise subcommand, got {other:?}"),
        };
        assert!(!hp_of(&base), "--indels alone must not enable the HP length model");
        let mut on = base.to_vec();
        on.push("--indel-hp-length-model");
        assert!(hp_of(&on), "the flag must enable it");
    }
}

#[cfg(test)]
mod lineage_opts_validation_tests {
    use super::*;

    #[test]
    fn mcmc_chains_zero_is_rejected() {
        let err = validate_lineage_opts(0).unwrap_err();
        assert_eq!(
            format!("{err:#}"),
            "--mcmc-chains must be at least 1 (got 0)"
        );
    }

    #[test]
    fn mcmc_chains_one_is_accepted() {
        assert!(validate_lineage_opts(1).is_ok());
    }

    #[test]
    fn collapse_hp_is_off_unless_the_flag_is_passed() {
        use clap::Parser;
        let base = ["Himito", "lineage", "-m", "matrix.csv", "-o", "out"];
        let collapse_of = |args: &[&str]| match Cli::try_parse_from(args).unwrap().command {
            Commands::Lineage { collapse_hp, .. } => collapse_hp,
            other => panic!("expected the Lineage subcommand, got {other:?}"),
        };
        assert!(!collapse_of(&base), "collapse must be off by default");
        let mut on = base.to_vec();
        on.push("--collapse-hp");
        assert!(collapse_of(&on), "--collapse-hp must enable it");
    }
}

#[cfg(test)]
mod linkage_cli_tests {
    use super::*;

    #[test]
    fn linkage_flags_parse() {
        let cli = Cli::try_parse_from([
            "Himito", "quick-start", "-i", "in.bam", "-r", "ref.fa", "-s", "S", "-o", "out", "--no-linkage-rescue",
        ])
        .unwrap();
        assert!(matches!(cli.command, Commands::QuickStart { no_linkage_rescue: true, .. }));
    }
}

#[cfg(test)]
mod bam_cli_tests {
    use super::*;

    const CALL: [&str; 9] = ["Himito", "call", "-g", "g.gfa", "-r", "ref.fa", "-s", "S", "-o"];

    fn call_args(extra: &[&str]) -> Vec<String> {
        let mut v: Vec<String> = CALL.iter().map(|s| s.to_string()).collect();
        v.push("out.vcf".into());
        v.extend(extra.iter().map(|s| s.to_string()));
        v
    }

    fn bam_of(extra: &[&str]) -> Option<PathBuf> {
        match Cli::try_parse_from(call_args(extra)).unwrap().command {
            Commands::Call { bam, .. } => bam,
            other => panic!("unexpected subcommand: {other:?}"),
        }
    }

    #[test]
    fn call_short_flag_sets_bam() {
        assert_eq!(bam_of(&["-b", "raw.bam"]), Some(PathBuf::from("raw.bam")));
        assert_eq!(bam_of(&["--bam", "raw.bam"]), Some(PathBuf::from("raw.bam")));
    }

    #[test]
    fn call_input_bam_is_an_alias_of_bam() {
        assert_eq!(bam_of(&["--input-bam", "raw.bam"]), Some(PathBuf::from("raw.bam")));
    }

    #[test]
    fn call_bam_is_optional() {
        assert_eq!(bam_of(&[]), None);
    }

    #[test]
    fn call_rejects_the_removed_bam_flags() {
        for flag in ["--genotype-bam", "--linkage-bam"] {
            let err = Cli::try_parse_from(call_args(&[flag, "x.bam"])).unwrap_err();
            assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument, "{flag}");
        }
    }

    #[test]
    fn quickstart_rejects_genotype_from_bam() {
        let err = Cli::try_parse_from([
            "Himito", "quick-start", "-i", "in.bam", "-r", "ref.fa", "-s", "S", "-o", "out", "--genotype-from-bam",
        ])
        .unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument);
    }
}

#[cfg(test)]
mod denoise_keep_vaf_cli_tests {
    use super::*;

    fn keep_vaf_of(extra: &[&str]) -> f64 {
        let mut args = vec!["Himito", "quick-start", "-i", "in.bam", "-r", "ref.fa", "-s", "S", "-o", "out"];
        args.extend_from_slice(extra);
        match Cli::try_parse_from(args).unwrap().command {
            Commands::QuickStart { denoise_keep_vaf, vaf_threshold, .. } => {
                assert_eq!(vaf_threshold, 0.01, "-v must stay independent of --denoise-keep-vaf");
                denoise_keep_vaf
            }
            other => panic!("expected the QuickStart subcommand, got {other:?}"),
        }
    }

    #[test]
    fn quickstart_denoise_keep_vaf_defaults_to_the_constant() {
        assert_eq!(keep_vaf_of(&[]), DENOISE_KEEP_VAF);
        assert_eq!(DENOISE_KEEP_VAF, 0.03);
    }

    #[test]
    fn quickstart_denoise_keep_vaf_explicit_value_wins() {
        assert_eq!(keep_vaf_of(&["--denoise-keep-vaf", "0.01"]), 0.01);
    }
}

// Effective left-alignment is the default; `--left-align` remains accepted and `--no-left-align` turns it off.
#[cfg(test)]
mod left_align_cli_tests {
    use super::*;

    fn parse(args: &[&str]) -> clap::error::Result<Cli> {
        Cli::try_parse_from(args)
    }

    fn effective_left_align(cli: Cli) -> bool {
        match cli.command {
            Commands::Call { no_left_align, .. } | Commands::QuickStart { no_left_align, .. } => {
                !no_left_align
            }
            other => panic!("unexpected subcommand: {other:?}"),
        }
    }

    #[test]
    fn call_left_aligns_by_default() {
        let cli = parse(&["Himito", "call", "-g", "g.gfa", "-r", "ref.fa", "-s", "S", "-o", "out.vcf"]).unwrap();
        assert!(effective_left_align(cli));
    }

    #[test]
    fn call_left_align_flag_still_left_aligns() {
        let cli = parse(&[
            "Himito", "call", "-g", "g.gfa", "-r", "ref.fa", "-s", "S", "-o", "out.vcf", "--left-align",
        ])
        .unwrap();
        assert!(effective_left_align(cli));
    }

    #[test]
    fn call_no_left_align_disables_it() {
        let cli = parse(&[
            "Himito", "call", "-g", "g.gfa", "-r", "ref.fa", "-s", "S", "-o", "out.vcf", "--no-left-align",
        ])
        .unwrap();
        assert!(!effective_left_align(cli));
    }

    #[test]
    fn call_rejects_both_left_align_flags() {
        let err = parse(&[
            "Himito", "call", "-g", "g.gfa", "-r", "ref.fa", "-s", "S", "-o", "out.vcf",
            "--left-align", "--no-left-align",
        ])
        .unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::ArgumentConflict);
    }

    #[test]
    fn quickstart_left_aligns_by_default() {
        let cli = parse(&[
            "Himito", "quick-start", "-i", "in.bam", "-r", "ref.fa", "-s", "S", "-o", "out",
        ])
        .unwrap();
        assert!(effective_left_align(cli));
    }

    #[test]
    fn quickstart_left_align_flag_still_left_aligns() {
        let cli = parse(&[
            "Himito", "quick-start", "-i", "in.bam", "-r", "ref.fa", "-s", "S", "-o", "out", "--left-align",
        ])
        .unwrap();
        assert!(effective_left_align(cli));
    }

    #[test]
    fn quickstart_no_left_align_disables_it() {
        let cli = parse(&[
            "Himito", "quick-start", "-i", "in.bam", "-r", "ref.fa", "-s", "S", "-o", "out", "--no-left-align",
        ])
        .unwrap();
        assert!(!effective_left_align(cli));
    }

    #[test]
    fn quickstart_rejects_both_left_align_flags() {
        let err = parse(&[
            "Himito", "quick-start", "-i", "in.bam", "-r", "ref.fa", "-s", "S", "-o", "out",
            "--left-align", "--no-left-align",
        ])
        .unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::ArgumentConflict);
    }
}

#[cfg(test)]
mod option_split_tests {
    use super::*;
    use clap::CommandFactory;

    struct Split {
        sub: &'static str,
        visible: &'static [&'static str],
        // (long flag, value to parse it with; "" for a switch)
        advanced: &'static [(&'static str, &'static str)],
        gone: &'static [&'static str],
    }

    const QUICK_START: Split = Split {
        sub: "quick-start",
        visible: &[
            "--input-bam", "--chromo", "--reference-path", "--sample-id", "--data-type",
            "--output-prefix", "--maximal-mt-depth", "--threads",
        ],
        advanced: &[
            ("--kmer-size", "21"), ("--minimal-ac", "2"), ("--vaf-threshold", "0.01"),
            ("--threshold-methylation-prob", "0.5"), ("--filter-max-methylation", "0.2"),
            ("--p-value-threshold", "0.01"), ("--heteroplasmic-frequency-threshold", "0.5"),
            ("--strand-bias-threshold", "0.05"), ("--indel-false-threshold", "0.1"),
            ("--permutation-frequency-threshold", "0.5"), ("--permutation-rounds", "1000"),
            ("--no-left-align", ""), ("--no-linkage-rescue", ""), ("--denoise-keep-vaf", "0.03"),
        ],
        gone: &[
            "--length-max", "--downsample-seed", "--strand-bias-drop-max-hf",
            "--strand-bias-action", "--permutation-mode", "--cooccurrence-stat",
            "--genotype-from-bam", "--left-align",
        ],
    };

    const CALL: Split = Split {
        sub: "call",
        visible: &[
            "--graphfile", "--reference-fasta", "--sample-id", "--data-type", "--output-file",
            "--bam", "--threads",
        ],
        advanced: &[
            ("--k", "21"), ("--minimal-ac", "2"), ("--vaf-threshold", "0.01"),
            ("--p-value-threshold", "0.01"), ("--frequency-threshold", "0.5"),
            ("--strand-bias-threshold", "0.01"), ("--indel-false-threshold", "0.1"),
            ("--permutation-frequency-threshold", "0.5"), ("--permutation-rounds", "1000"),
            ("--no-left-align", ""),
        ],
        gone: &[
            "--strand-bias-action", "--strand-bias-drop-max-hf", "--permutation-mode",
            "--genotype-max-ed", "--genotype-flank", "--cooccurrence-stat",
            "--genotype-bam", "--linkage-bam", "--left-align", "--input-bam",
        ],
    };

    const LINEAGE: Split = Split {
        sub: "lineage",
        visible: &[
            "--matrix-file", "--vcf-file", "--data-type", "--output-prefix", "--min-hf",
            "--max-hf", "--threads",
        ],
        advanced: &[
            ("--fp-rate", "0.005"), ("--fn-rate", "0.02"), ("--hp-fp-rate", "0.01"),
            ("--hp-fn-rate", "0.04"), ("--collapse-hp", ""), ("--mcmc-iterations", "100"),
            ("--mcmc-chains", "2"), ("--debug", ""),
        ],
        gone: &[
            "--no-root-block", "--root-block-min-hf", "--root-block-max-q",
            "--root-block-min-absent", "--root-block-min-or",
        ],
    };

    /// (short help, long help) of a subcommand, as `-h` and `--help` print them.
    fn helps(sub: &str) -> (String, String) {
        let mut root = Cli::command();
        root.build();
        let cmd = root.find_subcommand_mut(sub).expect("subcommand exists");
        (cmd.render_help().to_string(), cmd.render_long_help().to_string())
    }

    /// Long flags that head an option entry (not ones merely mentioned in prose).
    fn listed(help: &str) -> Vec<String> {
        help.lines()
            .map(str::trim_start)
            .filter(|l| l.starts_with('-'))
            .flat_map(|l| {
                l.split_whitespace()
                    .take_while(|t| t.starts_with('-'))
                    .filter(|t| t.starts_with("--"))
                    .map(|t| t.trim_end_matches(',').to_string())
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    fn cmdline(split: &Split, extra: &[&str]) -> Vec<String> {
        let base: &[&str] = match split.sub {
            "quick-start" => &["Himito", "quick-start", "-i", "in.bam", "-r", "ref.fa", "-s", "S", "-o", "out"],
            "call" => &["Himito", "call", "-g", "g.gfa", "-r", "ref.fa", "-s", "S", "-o", "out.vcf"],
            "lineage" => &["Himito", "lineage", "-m", "m.csv", "-o", "out"],
            other => panic!("no base command line for {other}"),
        };
        base.iter().chain(extra).map(|s| s.to_string()).collect()
    }

    fn check(split: &Split) {
        let (short, long) = helps(split.sub);
        let (short, long) = (listed(&short), listed(&long));
        let has = |v: &[String], f: &str| v.iter().any(|x| x == f);
        for f in split.visible {
            assert!(has(&short, f), "{} -h must list {f}", split.sub);
            assert!(has(&long, f), "{} --help must list {f}", split.sub);
        }
        for (f, _) in split.advanced {
            assert!(!has(&short, f), "{} -h must NOT list advanced {f}", split.sub);
            assert!(has(&long, f), "{} --help must list advanced {f}", split.sub);
        }
        for f in split.gone {
            assert!(!has(&short, f), "{} -h must not list {f}", split.sub);
            assert!(!has(&long, f), "{} --help must not list {f}", split.sub);
        }
    }

    fn check_parse(split: &Split) {
        for (f, v) in split.advanced {
            let extra: Vec<&str> = if v.is_empty() { vec![f] } else { vec![f, v] };
            if let Err(e) = Cli::try_parse_from(cmdline(split, &extra)) {
                panic!("{} {f}: advanced flag must still parse: {e}", split.sub);
            }
        }
        for f in split.gone {
            if *f == "--input-bam" || *f == "--left-align" {
                continue;
            }
            let err = Cli::try_parse_from(cmdline(split, &[f, "1"])).expect_err(f);
            assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument, "{} {f}", split.sub);
        }
    }

    #[test]
    fn quick_start_help_split() {
        check(&QUICK_START);
    }

    #[test]
    fn call_help_split() {
        check(&CALL);
    }

    #[test]
    fn lineage_help_split() {
        check(&LINEAGE);
    }

    #[test]
    fn quick_start_advanced_parse_and_removed_rejected() {
        check_parse(&QUICK_START);
    }

    #[test]
    fn call_advanced_parse_and_removed_rejected() {
        check_parse(&CALL);
    }

    #[test]
    fn lineage_advanced_parse_and_removed_rejected() {
        check_parse(&LINEAGE);
    }

    #[test]
    fn left_align_is_a_hidden_noop_on_quick_start_and_call() {
        for split in [&QUICK_START, &CALL] {
            Cli::try_parse_from(cmdline(split, &["--left-align"])).unwrap();
            let (short, long) = helps(split.sub);
            assert!(!listed(&short).iter().any(|f| f == "--left-align"));
            assert!(!listed(&long).iter().any(|f| f == "--left-align"));
        }
    }

    #[test]
    fn call_input_bam_alias_still_parses_but_is_not_listed() {
        Cli::try_parse_from(cmdline(&CALL, &["--input-bam", "raw.bam"])).unwrap();
        let (short, long) = helps("call");
        assert!(!listed(&short).iter().any(|f| f == "--input-bam"));
        assert!(!listed(&long).iter().any(|f| f == "--input-bam"));
    }

    #[test]
    fn data_type_is_validated_on_quick_start_and_call() {
        for split in [&QUICK_START, &CALL] {
            let err = Cli::try_parse_from(cmdline(split, &["-d", "ont-10"])).unwrap_err();
            assert_eq!(err.kind(), clap::error::ErrorKind::InvalidValue, "{}", split.sub);
            for dt in DATA_TYPES {
                Cli::try_parse_from(cmdline(split, &["-d", dt])).unwrap();
            }
        }
    }

    #[test]
    fn help_strings_are_one_short_line_without_defaults() {
        let mut root = Cli::command();
        root.build();
        for sub in ["quick-start", "call", "lineage"] {
            for arg in root.find_subcommand(sub).unwrap().get_arguments() {
                let id = arg.get_id();
                // clap's own -h/--help carries a long form by design.
                if arg.is_hide_set() || id == "help" || id == "version" {
                    continue;
                }
                assert!(arg.get_long_help().is_none(), "{sub} --{id}: long help set");
                let help = arg.get_help().map(|h| h.to_string()).unwrap_or_default();
                assert!(!help.is_empty(), "{sub} --{id}: no help");
                assert!(!help.contains('\n'), "{sub} --{id}: help spans lines");
                assert!(help.chars().count() <= 80, "{sub} --{id}: {} chars", help.chars().count());
                assert!(!help.to_lowercase().contains("default"), "{sub} --{id}: names a default");
            }
        }
    }

    #[test]
    fn lineage_root_block_is_always_on() {
        assert_eq!(LINEAGE_ROOT_BLOCK_MIN_HF, 0.80);
        assert_eq!(LINEAGE_ROOT_BLOCK_MAX_Q, 0.05);
        assert_eq!(LINEAGE_ROOT_BLOCK_MIN_ABSENT, 10);
        assert_eq!(LINEAGE_ROOT_BLOCK_MIN_OR, 2.0);
    }
}
