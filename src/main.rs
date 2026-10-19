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

pub const DENOISE_KEEP_VAF: f64 = 0.03;

/// The data types Himito accepts wherever a `-d/--data-type` is taken.
pub const DATA_TYPES: [&str; 4] = ["pacbio", "ont-r9", "ont-r10", "ont-denoised"];

/// Quick-start: longest read span (bp) the graph build aligns.
pub const QUICK_START_LENGTH_MAX: usize = 3000;

/// Quick-start: RNG seed of the mt-read downsampling, fixed for reproducible subsets.
pub const QUICK_START_DOWNSAMPLE_SEED: u64 = 42;

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
    #[arg(long, global = true, value_parser = parse_positive)]
    pub threads: Option<usize>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Quick start pipeline for mitochondrial assembly, variant calling and methylation aggregation
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

        /// sequencing data type; sets the filter and denoise presets (pacbio, ont-r9, ont-r10, ont-denoised)
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

        /// minimum heteroplasmic frequency for a variant [default: 0.01]
        #[clap(short, long, value_parser, hide_short_help = true)]
        vaf_threshold: Option<f32>,

        /// permutation-test p-value cutoff; 1 disables it (preset by data type)
        #[clap(long, value_parser, hide_short_help = true)]
        p_value_threshold: Option<f64>,

        /// permutation-test frequency threshold (preset by data type)
        #[clap(long, value_parser, hide_short_help = true)]
        heteroplasmic_frequency_threshold: Option<f64>,

        /// p-value cutoff of the strand-bias filter
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

        /// keep the graph's indel spellings instead of left-aligning them
        #[clap(long = "no-left-align", action, hide_short_help = true)]
        no_left_align: bool,

        /// turn off read-linkage rescue of low-frequency SNVs
        #[clap(long = "no-linkage-rescue", action, hide_short_help = true)]
        no_linkage_rescue: bool,

        /// keep threshold of the denoise site model (ont-denoised only)
        #[clap(long, value_parser, default_value_t = DENOISE_KEEP_VAF, hide_short_help = true)]
        denoise_keep_vaf: f64,

        /// PacBio: out-of-tract small indels with HF from this up to --indel-false-threshold PASS only if they pass the linkage test; a value at or above --indel-false-threshold turns this off [default: 0.015]
        #[clap(long, value_parser, hide_short_help = true)]
        indel_linked_min_hf: Option<f64>,
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

        /// data type; only ont-* is denoised, pacbio passes through unchanged
        #[clap(short, long, default_value = "ont-r10", value_parser = clap::builder::PossibleValuesParser::new(DATA_TYPES))]
        data_type: String,

        /// minimal VAF for an allele to be kept (below this it is corrected away)
        #[clap(long, value_parser, default_value_t = DENOISE_KEEP_VAF)]
        vaf: f64,

        /// enable small-indel (<5bp) correction in addition to substitutions
        #[clap(long, value_parser, default_value_t = false)]
        indels: bool,

        /// an allele at or above this fraction of reads at a site is never corrected away.
        #[clap(long, value_parser = parse_protect_vaf, default_value_t = 0.2)]
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

        /// minimum heteroplasmic frequency for a variant [default: 0.01]
        #[clap(short, long, value_parser, hide_short_help = true)]
        vaf_threshold: Option<f32>,

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

        /// keep the graph's indel spellings instead of left-aligning them
        #[clap(long = "no-left-align", action, hide_short_help = true)]
        no_left_align: bool,

        /// PacBio: out-of-tract small indels with HF from this up to --indel-false-threshold PASS only if they pass the linkage test; a value at or above --indel-false-threshold turns this off [default: 0.015]
        #[clap(long, value_parser, hide_short_help = true)]
        indel_linked_min_hf: Option<f64>,
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

    /// Infer mitochondrial lineage tree from Himito read_var matrix CSV and heteroplasmic variants
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
        #[clap(long, value_parser = parse_positive, default_value_t = 3, hide_short_help = true)]
        mcmc_chains: usize,

        /// also write raw_haplotype_map.tsv and raw_variant_cooccurrence.tsv
        #[clap(long, action, hide_short_help = true)]
        debug: bool,
    }
}

/// `--indel-protect-vaf` must be in (0, 1]; 0 would protect every allele.
/// Counts that must be at least 1.
fn parse_positive(s: &str) -> Result<usize, String> {
    match s.parse::<usize>() {
        Ok(n) if n >= 1 => Ok(n),
        _ => Err(format!("must be a positive integer, got {s}")),
    }
}

fn parse_protect_vaf(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|e| format!("{e}"))?;
    if v > 0.0 && v <= 1.0 { Ok(v) } else { Err(format!("must be in (0, 1], got {v}")) }
}

pub fn init_rayon_threads(threads: Option<usize>) -> AnyhowResult<()> {
    let Some(n) = threads else {
        return Ok(());
    };
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
    if let Err(e) = init_rayon_threads(args.global.threads) {
        eprintln!("Error: {e:#}");
        std::process::exit(1);
    }
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
            no_linkage_rescue,
            denoise_keep_vaf,
            indel_linked_min_hf,
        } => {
            let left_align = !no_left_align;
            let (p_value_threshold, frequency_threshold, permutation_frequency_threshold) =
                call::resolve_thresholds(&data_type, p_value_threshold, heteroplasmic_frequency_threshold, permutation_frequency_threshold);
            let indel_false_threshold = call::resolve_indel_false_threshold(&data_type, indel_false_threshold);
            let vaf_threshold = call::resolve_vaf(&data_type, vaf_threshold);
            let indel_linked_min_hf = call::resolve_indel_linked_min_hf(indel_linked_min_hf);
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
            
            // ont-denoised builds the graph from a denoised copy of the mt BAM (indel correction on).
            let mut build_bam = mt_output.clone();
            if data_type == "ont-denoised" {
                let denoised = output_prefix.with_extension("mt.denoised.bam");
                let iopts = denoise_indel::IndelOpts { enabled: true, ..Default::default() };
                // `denoise_keep_vaf` is the denoiser's own keep threshold, not `vaf_threshold`.
                match denoise::start(&mt_output, &denoised, &reference_path, &data_type, denoise_keep_vaf, &iopts, None) {
                    Ok(()) => build_bam = denoised,
                    Err(e) => eprintln!("Warning: denoise failed ({e:#}); falling back to raw mt BAM."),
                }
            }

            let graph_output = output_prefix.with_extension("gfa");
            let _ = build::start(&graph_output, kmer_size, &build_bam, &reference_path, QUICK_START_LENGTH_MAX, 2);

            let assemble_output = output_prefix.with_extension("fasta");
            asm::start(&graph_output, &assemble_output, &sample_id);

            let vcf_output = output_prefix.with_extension("vcf");
            call::start(
                &graph_output,
                &reference_path,
                kmer_size,
                minimal_ac,
                &vcf_output,
                &sample_id,
                vaf_threshold,
                &data_type,
                p_value_threshold,
                frequency_threshold,
                // Raw mt BAM, not build_bam (denoise turns uncertain ALTs into REF).
                Some(&mt_output),
                strand_bias_threshold,
                indel_false_threshold,
                permutation_frequency_threshold,
                permutation_rounds,
                left_align,
                no_linkage_rescue,
                indel_linked_min_hf,
            );
            let methyl_output = output_prefix.with_extension("bed");
            methyl::start(&graph_output, &mt_output, &methyl_output, threshold_methylation_prob, false);
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
            };
            if let Err(e) = denoise::start(&input, &output, &reference, &data_type, vaf, &iopts, stats.as_ref()) {
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
            indel_linked_min_hf,
        } => {
            let left_align = !no_left_align;
            let (p_value_threshold, frequency_threshold, permutation_frequency_threshold) =
                call::resolve_thresholds(&data_type, p_value_threshold, frequency_threshold, permutation_frequency_threshold);
            let indel_false_threshold = call::resolve_indel_false_threshold(&data_type, indel_false_threshold);
            let vaf_threshold = call::resolve_vaf(&data_type, vaf_threshold);
            let indel_linked_min_hf = call::resolve_indel_linked_min_hf(indel_linked_min_hf);
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
                left_align,
                false,
                indel_linked_min_hf,
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
            if let Err(e) = lineage::start(
                matrix_file,
                vcf_file,
                min_hf,
                max_hf,
                fp_rate,
                fn_rate,
                hp_rates,
                tracts,
                mcmc_iterations,
                mcmc_chains,
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
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn value_parsers_enforce_ranges() {
        for (v, ok) in [("0", false), ("1e-9", true), ("1", true), ("1.5", false), ("NaN", false)] {
            assert_eq!(parse_protect_vaf(v).is_ok(), ok, "{v}");
        }
        for (v, ok) in [("0", false), ("1", true), ("-1", false), ("x", false)] {
            assert_eq!(parse_positive(v).is_ok(), ok, "{v}");
        }
    }

    /// docs/cli/ holds the `-h` and `--help` text of the commands that hide advanced options from `-h`.
    #[test]
    fn docs_cli_match_help() {
        let mut root = Cli::command();
        root.build();
        for sub in ["quick-start", "call", "lineage"] {
            let cmd = root.find_subcommand_mut(sub).unwrap();
            for (kind, help) in [("short", cmd.render_help()), ("long", cmd.render_long_help())] {
                let path = format!("{}/docs/cli/{sub}.{kind}.txt", env!("CARGO_MANIFEST_DIR"));
                assert_eq!(help.to_string(), std::fs::read_to_string(&path).unwrap(), "{path} is out of date");
            }
        }
    }
}
