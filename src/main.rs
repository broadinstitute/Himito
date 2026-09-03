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
mod scite;
mod denoise;
mod denoise_indel;

/// Keep threshold for denoise's per-site substitution model: an eligible non-ref
/// allele survives only if the site model's EM frequency estimate reaches this.
///
/// Deliberately NOT the caller's `--vaf`, which the two used to share. That one is a
/// raw HF cut on an already-called variant; this one is compared against a
/// quality-weighted estimate that discounts alt observations explainable as error,
/// which makes it far sharper on noise whose raw HF overlaps real heteroplasmies.
///
/// 0.03 was chosen by sweeping 0.01-0.04 over the 20 ont-r10 n=10 depth>=300 cells in
/// lineage_eval: it lifts variant precision 0.837 -> 0.967 while holding variant
/// recall at 1.000 and ancestor-descendant F1 within 0.003 of the optimum. Above
/// 0.04 recall starts falling (0.965) as real low-frequency heteroplasmies drop out.
pub const DENOISE_KEEP_VAF: f64 = 0.03;

/// Minimal observations required on EACH strand for an allele to be a candidate.
///
/// Was `--min-strand`. Every caller in the repo — QuickStart, both eval scripts,
/// and every unit test — used 2, so the flag never did anything but restate its
/// own default.
pub const DENOISE_MIN_STRAND: u32 = 2;

/// p-value threshold for the per-allele strand-bias test (binomial, against the
/// column's own forward fraction); alleles below it lose candidacy.
///
/// Was `--strand-bias-p`. `0.0` disables the test entirely — a path `fit_site`'s
/// unit tests still exercise, which is why the plumbing below keeps this as a
/// parameter even though no CLI flag sets it any more.
pub const DENOISE_STRAND_BIAS_P: f64 = 0.01;

/// Within-column frequency at or above which an allele skips the strand-bias
/// test (a single-strand artifact cannot reach a near-homoplasmic frequency).
///
/// Was `--homoplasmic-vaf`. Note this is the *standalone denoise* value; the
/// QuickStart pipeline deliberately passes 0.7 instead (see the `QuickStart`
/// arm), matching `call.rs`'s own homoplasmic threshold.
pub const DENOISE_HOMOPLASMIC_VAF: f64 = 0.95;

/// The data types Himito accepts wherever a `-d/--data-type` is taken.
///
/// Used as the clap `PossibleValuesParser` set so an unrecognised value is
/// rejected before any work starts. This matters because the fp/fn preset
/// lookup in `lineage::resolve_error_rates` has a catch-all arm: without this
/// list a typo'd `-d ont-10` would silently run to completion with raw-ONT
/// rates and exit 0.
pub const DATA_TYPES: [&str; 4] = ["pacbio", "ont-r9", "ont-r10", "ont-denoised"];

/// Minimal reads a variant must be present in / absent from for `lineage` to
/// treat it as informative (it needs both to place a bifurcation).
///
/// Were `--min-presence` / `--min-absence`. Nothing in the repo — no WDL, no
/// eval script, no doc — ever set either, and `git log -S` finds no commit that
/// ever did.
pub const LINEAGE_MIN_PRESENCE: usize = 2;
pub const LINEAGE_MIN_ABSENCE: usize = 1;

/// Minimal reads backing a haplotype for `lineage` to report it.
///
/// Was `--min-reads`. Drops low-count haplotypes from the haplotype maps and
/// from the tips of `<prefix>.read_lineage.nwk`.
pub const LINEAGE_MIN_READS: usize = 3;

/// RNG seed for the SCITE MCMC search, fixed so runs are reproducible.
///
/// Was `--mcmc-seed`. Every invocation in the repo used the default 42;
/// `--mcmc-iterations` and `--mcmc-chains` remain flags because
/// `lineage_eval/lineage_sim/sweep_fpfn.sh` does vary those two.
pub const LINEAGE_MCMC_SEED: u64 = 42;

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
    /// Global options
    #[arg(long, global = true)]
    pub threads: Option<usize>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Quick start pipeline for mitochondrial assembly, variant calling and methylation aggregation
    #[clap(arg_required_else_help = true)]
    QuickStart {
        /// input path for bam file.
        #[clap(short, long, value_parser, required = true)]
        input_bam: PathBuf,

        /// contig name in the bam file
        #[clap(short, long, value_parser, default_value = "chrM")]
        chromo: String,

        /// Reference Fasta file, rCRS
        #[clap(short, long, value_parser, required = true)]
        reference_path: PathBuf,

        /// sample name of the bam file
        #[clap(short, long, value_parser, required = true)]
        sample_id: String,

        /// data type, pacbio, ont-r9, ont-r10
        #[clap(short, long, value_parser, default_value = "pacbio")]
        data_type: String,

        /// output prefix
        #[clap(short, long, value_parser, required = true)]
        output_prefix: PathBuf,

        /// min_probability to determine a C is methylated
        #[clap(short, long, value_parser, default_value_t = 0.5)]
        threshold_methylation_prob: f64,

        /// max fraction to keep a read as mtDNA derived
         #[clap(short, long, value_parser, default_value_t = 0.2)]
        filter_max_methylation: f64,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = 21)]
        kmer_size: usize,

        /// minimal allele count for variants
        #[clap(short, long, value_parser, default_value_t = 2)]
        minimal_ac: usize,

        /// max Length to do alignment
        #[clap(short, long, value_parser, default_value_t = 3000)]
        length_max: usize,

        /// minimal heteroplasmic frequency for variants
        #[clap(short, long, value_parser, default_value_t = 0.01)]
        vaf_threshold: f32,

        /// p-value threshold for permutation test (optional, preset by data-type,(optional, set as 1 to disable permutation test))
        #[clap(long, value_parser)]
        p_value_threshold: Option<f64>,

        /// heteroplasmic frequency threshold for permutation test (optional, preset by data-type,higher is more stringent)
        #[clap(long, value_parser)]
        heteroplasmic_frequency_threshold: Option<f64>,

        /// maximal reads to keep from mtDNA 
        #[clap(long, value_parser, default_value_t = 50000)]
        maximal_mt_depth: usize,

        /// RNG seed for mt read downsampling (reproducible subsampling)
        #[clap(long, value_parser, default_value_t = 42)]
        downsample_seed: u64,

        /// strand bias threshold for filtering variants
        #[clap(long, value_parser, default_value_t = 0.05)]
        strand_bias_threshold: f64,

        /// indel false threshold for filtering variants
        #[clap(long, value_parser, default_value_t = 0.1)]
        indel_false_threshold: f64,

        /// heteroplasmic frequency threshold for permutation test
        #[clap(long, value_parser )]
        permutation_frequency_threshold: Option<f64>,

        /// permutation rounds per variant. The smallest attainable p-value is
        /// 1/(rounds+1), so this must leave headroom below --p-value-threshold after
        /// Benjamini-Hochberg. Cost scales linearly; lower it to speed up large graphs.
        #[clap(long, value_parser, default_value_t = 1000)]
        permutation_rounds: usize,

        /// what to do with strand-skewed variants: flag them in FILTER, or additionally
        /// drop low-VAF SNPs outright (default: drop for ont-denoised, flag otherwise)
        #[clap(long, value_enum)]
        strand_bias_action: Option<call::StrandBiasAction>,

        /// heteroplasmic-frequency ceiling below which a strand-skewed SNP is dropped
        /// rather than merely flagged (only used with --strand-bias-action drop)
        #[clap(long, value_parser, default_value_t = 0.10)]
        strand_bias_drop_max_hf: f64,

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
        /// Variants are read only from CIGAR-bearing edges, and reads reaching a
        /// bubble solely through an un-CIGARed edge are not counted as covering it
        /// (their genotype stays missing rather than ref). The default of 1 aligns
        /// every read-supported edge; raise it to trade recall of low-frequency
        /// variants for fewer false positives on noisy long reads.
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

        /// an allele carried by at least this fraction of reads at a site is never
        /// corrected away, whatever candidacy says (protects real heteroplasmy in
        /// deep repeat contexts, where the candidacy floor can exceed 1.0)
        #[clap(long, value_parser, default_value_t = 0.2)]
        indel_protect_vaf: f64,

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

    ///Call Variants from Sequence Graph
    #[clap(arg_required_else_help = true)]
    Call {
        /// path for anchor graph.
        #[clap(short, long, value_parser)]
        graphfile: PathBuf,

        /// Reference Fasta path
        #[clap(short, long, value_parser, required = true)]
        reference_fasta: PathBuf,

        /// sample name of the bam file
        #[clap(short, long, value_parser, required = true)]
        sample_id: String,

        /// data type, pacbio, ont-r9, ont-r10, ont-denoised
        #[clap(short, long, value_parser, default_value = "pacbio")]
        data_type: String,

        /// output file name
        #[clap(short, long, value_parser, required = true)]
        output_file: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = 21)]
        k: usize,

        /// minimal allele count for variants
        #[clap(short, long, value_parser, default_value_t = 2)]
        minimal_ac: usize,

        /// minimal heteroplasmic frequency for variants
        #[clap(short, long, value_parser, default_value_t = 0.01)]
        vaf_threshold: f32,

        /// p-value threshold for permutation test (optional, set as 1 to disable permutation test)
        #[clap(short, long, value_parser)]
        p_value_threshold: Option<f64>,

        /// frequency threshold for permutation test (optional)
        #[clap(short, long, value_parser)]
        frequency_threshold: Option<f64>,

        /// optional BAM file; when provided, variants are additionally filtered by strand bias
        #[clap(long, value_parser)]
        input_bam: Option<PathBuf>,

        /// p-value threshold for the strand-bias filter (only used with --input-bam)
        #[clap(long, value_parser, default_value_t = 0.01)]
        strand_bias_threshold: f64,

        /// indel false threshold for filtering variants
        #[clap(long, value_parser, default_value_t = 0.1)]
        indel_false_threshold: f64,

        /// heteroplasmic frequency threshold for permutation test (optional, preset by data-type)
        #[clap(long, value_parser)]
        permutation_frequency_threshold: Option<f64>,

        /// permutation rounds per variant. The smallest attainable p-value is
        /// 1/(rounds+1), so this must leave headroom below --p-value-threshold after
        /// Benjamini-Hochberg. Cost scales linearly; lower it to speed up large graphs.
        #[clap(long, value_parser, default_value_t = 1000)]
        permutation_rounds: usize,

        /// what to do with strand-skewed variants: flag them in FILTER, or additionally
        /// drop low-VAF SNPs outright (default: drop for ont-denoised, flag otherwise)
        #[clap(long, value_enum)]
        strand_bias_action: Option<call::StrandBiasAction>,

        /// heteroplasmic-frequency ceiling below which a strand-skewed SNP is dropped
        /// rather than merely flagged (only used with --strand-bias-action drop)
        #[clap(long, value_parser, default_value_t = 0.10)]
        strand_bias_drop_max_hf: f64,
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
        /// extract the per-read level methylation signals on major haplotype or all the reads 
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

        /// minimal read ratio supporting each haplotype in a given window, default is 0.01
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

    /// Infer mitochondrial lineage tree from Himito read_var matrix CSV and heteroplasmic variants
    #[clap(arg_required_else_help = true)]
    Lineage {
        /// path for the Himito matrix CSV (<prefix>.matrix.csv)
        #[clap(short, long, value_parser, required = true)]
        matrix_file: PathBuf,

        /// path for the Himito VCF; used to read each variant's HF and
        /// restrict the analysis to heteroplasmic (non-fixed) sites
        #[clap(short, long, value_parser)]
        vcf_file: Option<PathBuf>,

        /// minimal heteroplasmic frequency (inclusive) for a variant to be considered
        #[clap(long, value_parser, default_value_t = 0.01)]
        min_hf: f64,

        /// maximal heteroplasmic frequency (exclusive) for a variant to be
        /// considered; the default excludes fixed/homoplasmic variants (HF >= 0.95)
        #[clap(long, value_parser, default_value_t = 0.95)]
        max_hf: f64,

        /// data type, pacbio, ont-r9, ont-r10, ont-denoised; selects the
        /// default SCITE fp/fn rates
        #[clap(
            short,
            long,
            default_value = "pacbio",
            value_parser = clap::builder::PossibleValuesParser::new(DATA_TYPES)
        )]
        data_type: String,

        /// SCITE false-positive rate (alpha): P(observed=1 | true=0)
        /// (optional; preset by data type)
        #[clap(long, value_parser)]
        fp_rate: Option<f64>,

        /// SCITE false-negative rate (beta): P(observed=0 | true=1)
        /// (optional; preset by data type)
        #[clap(long, value_parser)]
        fn_rate: Option<f64>,

        /// number of MCMC iterations per chain
        #[clap(long, value_parser, default_value_t = 5000)]
        mcmc_iterations: usize,

        /// number of independent MCMC chains (best tree across all is kept)
        #[clap(long, value_parser, default_value_t = 3)]
        mcmc_chains: usize,

        /// output prefix; writes <prefix>.haplotype_map.tsv
        #[clap(short, long, value_parser, required = true)]
        output_prefix: String,
    }
}

/// Validate the user-settable `--indel-*` options at the CLI boundary. Only
/// called when `--indels` is passed (see the `Denoise` arm below): a user who is
/// not using the feature must never be blocked by it.
///
/// This validation deliberately lives here, not inside `IndelOpts` or any
/// library function: the CLI is the boundary where a human-typed flag value
/// first exists, and the error message needs the flag's own name, which only
/// the CLI layer knows.
///
/// Only `protect_vaf` is checked, because it is the only numeric field a flag
/// can still set. The rest of `IndelOpts` now always holds
/// `IndelOpts::default()`, whose values are in range by construction — checking
/// a number the user cannot type would be dead code.
///
/// `--indel-protect-vaf 0` is load-bearing, not cosmetic: it protects EVERY
/// non-strand-rejected read (`observed_vaf >= 0.0` is always true), so the only
/// reads ever corrected are ones that failed a strand gate — the exact inverse
/// of the "protect less, correct more" a value of 0 looks like it should mean,
/// and it would otherwise silently disable the feature while still exiting 0
/// with a valid BAM.
fn validate_indel_opts(o: &denoise_indel::IndelOpts) -> AnyhowResult<()> {
    anyhow::ensure!(
        o.protect_vaf > 0.0 && o.protect_vaf <= 1.0,
        "--indel-protect-vaf must be in (0, 1] (got {})", o.protect_vaf
    );
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
            length_max,
            kmer_size,
            minimal_ac,
            vaf_threshold,
            p_value_threshold,
            heteroplasmic_frequency_threshold,
            maximal_mt_depth,
            downsample_seed,
            strand_bias_threshold,
            indel_false_threshold,
            permutation_frequency_threshold,
            permutation_rounds,
            strand_bias_action,
            strand_bias_drop_max_hf,
        } => {
            let (p_value_threshold, frequency_threshold, permutation_frequency_threshold_) =
                call::resolve_thresholds(&data_type, p_value_threshold, heteroplasmic_frequency_threshold, permutation_frequency_threshold);
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
                downsample_seed,
            );
            
            // For ONT, denoise the mt BAM so reads consolidate; PacBio is untouched.
            // Methylation (below) keeps reading the ORIGINAL mt_output.
            // data_type_ only flips to "ont-denoised" when denoise actually succeeds;
            // on fallback the raw reads keep their original type so the caller still
            // applies the ONT permutation-test SNP filtering those reads need.
            let (build_bam, data_type_) = if data_type == "ont-denoised" {
                let denoised = output_prefix.with_extension("mt.denoised.bam");
                // Denoise settings for the combined pipeline. The strand gates are
                // the shared constants; the near-homoplasmic exemption is
                // deliberately 0.7 here — matching `call.rs`'s own
                // permutation/homoplasmic threshold — and NOT the standalone
                // denoise value in `DENOISE_HOMOPLASMIC_VAF` (0.95). The two paths
                // have always differed; an earlier version of this comment claimed
                // they matched, which was wrong.
                if let Err(e) = denoise::start(
                    &mt_output, &denoised, &reference_path, &data_type,
                    // Denoise's keep threshold, NOT `vaf_threshold`. Passing the
                    // caller's HF cut here made the two impossible to tune apart, and
                    // at 0.01 the site model keeps marginal noise alleles whose reads
                    // then reach the VCF (precision 0.837 vs 0.967 at 0.03).
                    DENOISE_KEEP_VAF, DENOISE_MIN_STRAND, DENOISE_STRAND_BIAS_P, 0.7,
                    // Indel correction is on in QuickStart. Every other field keeps
                    // its `Default` value, which `validate_indel_opts` accepts as-is,
                    // so there is nothing to validate here.
                    &denoise_indel::IndelOpts { enabled: true, ..Default::default() },
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
            // Preserve the historical default edge-read gate (2) for this combined pipeline.
            let _ = build::start(&graph_output, kmer_size, &build_bam, &reference_path, length_max, 2);

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
                Some(&build_bam),
                strand_bias_threshold,
                indel_false_threshold,
                permutation_frequency_threshold_,
                permutation_rounds,
                strand_bias_action,
                strand_bias_drop_max_hf,
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
            stats,
        } => {
            // Every field except `enabled` and `protect_vaf` keeps its `Default`
            // value: the flags that used to set them all shipped defaults
            // identical to `IndelOpts::default()` and nothing ever overrode them.
            let iopts = denoise_indel::IndelOpts {
                enabled: indels,
                protect_vaf: indel_protect_vaf,
                ..Default::default()
            };
            // Only validate when --indels is actually in effect: a user who
            // never opted into the feature must never be blocked by it.
            if iopts.enabled {
                if let Err(e) = validate_indel_opts(&iopts) {
                    eprintln!("Error: {e:#}");
                    std::process::exit(1);
                }
            }
            if let Err(e) = denoise::start(
                &input, &output, &reference, &data_type, vaf, DENOISE_MIN_STRAND,
                DENOISE_STRAND_BIAS_P, DENOISE_HOMOPLASMIC_VAF, &iopts, stats.as_ref(),
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
            input_bam,
            strand_bias_threshold,
            indel_false_threshold,
            permutation_frequency_threshold,
            permutation_rounds,
            strand_bias_action,
            strand_bias_drop_max_hf,
        } => {
            let (p_value_threshold, frequency_threshold, permutation_frequency_threshold) =
                call::resolve_thresholds(&data_type, p_value_threshold, frequency_threshold, permutation_frequency_threshold);
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
                input_bam.as_ref(),
                strand_bias_threshold,
                indel_false_threshold,
                permutation_frequency_threshold,
                permutation_rounds,
                strand_bias_action,
                strand_bias_drop_max_hf,
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
            mcmc_iterations,
            mcmc_chains,
            output_prefix,
        } => {
            let matrix_file = matrix_file.to_str().expect("matrix-file path is not valid UTF-8");
            let vcf_file = vcf_file
                .as_deref()
                .map(|p| p.to_str().expect("vcf-file path is not valid UTF-8"));
            let (fp_rate, fn_rate) =
                lineage::resolve_error_rates(&data_type, fp_rate, fn_rate);
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
                mcmc_iterations,
                mcmc_chains,
                LINEAGE_MCMC_SEED,
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
        // FIX 4: --indel-protect-vaf 0 protects every non-strand-rejected
        // read, so it is the exact inverse of "protect less, correct more"
        // and must be rejected, not silently accepted with a zero exit code.
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

    // The former per-field range tests (err_cap, err0, err_scale, floor_mult,
    // delta, vaf, max_len, flank) were deleted along with the flags that set
    // those fields. They now always hold `IndelOpts::default()`, which
    // `shipped_defaults_pass_validation` above covers; a test that pokes a
    // value no user can supply would only be testing the poke.
}
