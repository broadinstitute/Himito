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
        #[clap(long, value_parser, default_value_t = 0.01)]
        vaf: f64,

        /// minimal observations required on EACH strand for allele candidacy
        #[clap(long, value_parser, default_value_t = 2)]
        min_strand: u32,

        /// p-value threshold for the per-allele strand-bias test (binomial, against the
        /// column's own forward fraction); alleles below it lose candidacy. 0 disables.
        #[clap(long, value_parser, default_value_t = 0.01)]
        strand_bias_p: f64,

        /// alleles at or above this within-column frequency skip the strand-bias test
        /// (a single-strand artifact cannot reach a near-homoplasmic frequency)
        #[clap(long, value_parser, default_value_t = 0.95)]
        homoplasmic_vaf: f64,

        /// enable small-indel (<5bp) correction in addition to substitutions
        #[clap(long, value_parser, default_value_t = false)]
        indels: bool,

        /// exclusive upper bound on correctable indel length (5 = lengths 1..=4)
        #[clap(long, value_parser, default_value_t = 5)]
        indel_max_len: u32,

        /// absolute minimal VAF floor for an indel allele to be kept
        #[clap(long, value_parser, default_value_t = 0.05)]
        indel_vaf: f64,

        /// per-junction indel error probability in unique sequence
        #[clap(long, value_parser, default_value_t = 0.01)]
        indel_err0: f64,

        /// multiplicative growth of the indel error rate per extra repeat copy
        #[clap(long, value_parser, default_value_t = 1.5)]
        indel_err_scale: f64,

        /// ceiling on the context-scaled indel error rate
        #[clap(long, value_parser, default_value_t = 0.4)]
        indel_err_cap: f64,

        /// candidacy floor is at least this multiple of the local indel error rate
        #[clap(long, value_parser, default_value_t = 3.0)]
        indel_floor_mult: f64,

        /// decay of error mass per unit of net-length distance between alleles
        #[clap(long, value_parser, default_value_t = 0.3)]
        indel_delta: f64,

        /// reads must extend this many bases past a site on both sides to be reassigned
        #[clap(long, value_parser, default_value_t = 5)]
        indel_flank: usize,

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

        /// maximal heteroplasmic frequency (exclusive) for a variant to be considered;
        /// default of 1.0 excludes fixed/homoplasmic variants (HF == 0.95)
        #[clap(long, value_parser, default_value_t = 0.95)]
        max_hf: f64,

        /// minimal number of reads a variant must be present in to be informative
        #[clap(long, value_parser, default_value_t = 2)]
        min_presence: usize,

        /// minimal number of reads a variant must be absent from to be informative
        #[clap(long, value_parser, default_value_t = 1)]
        min_absence: usize,

        /// minimal number of reads required to report a haplotype
        #[clap(long, value_parser, default_value_t = 3)]
        min_reads: usize,

        /// data type, pacbio, ont-r9, ont-r10, ont-denoised; selects the
        /// default SCITE fp/fn rates
        #[clap(short, long, value_parser, default_value = "pacbio")]
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

        /// RNG seed for the MCMC search (reproducible runs)
        #[clap(long, value_parser, default_value_t = 42)]
        mcmc_seed: u64,

        /// output prefix; writes <prefix>.haplotype_map.tsv
        #[clap(short, long, value_parser, required = true)]
        output_prefix: String,
    }
}

/// Validate the numeric `--indel-*` options at the CLI boundary. Only called
/// when `--indels` is passed (see the `Denoise` arm below): a user who is not
/// using the feature must never be blocked by it, however out-of-range its
/// unused defaults might theoretically be set to by some other caller.
///
/// This validation deliberately lives here, not inside `IndelOpts` or any
/// library function: the CLI is the boundary where a human-typed flag value
/// first exists, and the error message needs the flag's own name, which only
/// the CLI layer knows. Two of these are load-bearing, not cosmetic:
/// `--indel-protect-vaf 0` protects EVERY non-strand-rejected read
/// (`observed_vaf >= 0.0` is always true), so the only reads ever corrected
/// are ones that failed a strand gate -- the exact inverse of "protect less,
/// correct more" that a value of 0 looks like it should mean, and it would
/// otherwise silently disable the feature while still exiting 0 with a valid
/// BAM. `--indel-err-cap` above 0.5 is honoured by the candidacy floor
/// (`vaf_floor` via `floor_mult * error_rate`) but `assign_allele` clamps
/// `eps` to `[0.0, 0.5]` before using it for MAP assignment, so a cap above
/// 0.5 would make the gate and the assignment step disagree about the error
/// rate actually in effect.
fn validate_indel_opts(o: &denoise::indel::IndelOpts) -> AnyhowResult<()> {
    fn open_closed(flag: &str, v: f64, hi: f64) -> AnyhowResult<()> {
        anyhow::ensure!(
            v > 0.0 && v <= hi,
            "--{flag} must be in (0, {hi}] (got {v})"
        );
        Ok(())
    }
    open_closed("indel-protect-vaf", o.protect_vaf, 1.0)?;
    open_closed("indel-err0", o.err0, 0.5)?;
    open_closed("indel-err-cap", o.err_cap, 0.5)?;
    // Upper bounds matter as much as lower ones here: an unbounded multiplier can
    // push `vaf_floor` past 1.0 at every context, leaving `--indel-protect-vaf` as
    // the only thing keeping any allele -- a silent near-disable of exactly the kind
    // the `--indel-protect-vaf 0` check above exists to prevent.
    anyhow::ensure!(
        o.err_scale >= 1.0 && o.err_scale <= 10.0,
        "--indel-err-scale must be in [1.0, 10.0] (got {})", o.err_scale
    );
    anyhow::ensure!(
        o.floor_mult >= 0.0 && o.floor_mult <= 100.0,
        "--indel-floor-mult must be in [0.0, 100.0] (got {})", o.floor_mult
    );
    open_closed("indel-delta", o.delta, 1.0)?;
    // `vaf` is only ever used as `max(vaf, floor_mult * error_rate(L))`, so 0 is a
    // coherent setting meaning "no absolute floor, use the context-scaled one only".
    anyhow::ensure!(
        o.vaf >= 0.0 && o.vaf <= 1.0,
        "--indel-vaf must be in [0, 1] (got {})", o.vaf
    );
    // The feature is specified and documented as small-indel (<5bp) correction; a
    // large bound would let the rewrite walk emit correspondingly large gained
    // deletions, well outside anything the site model was designed or validated for.
    anyhow::ensure!(
        (1..=20).contains(&o.max_len),
        "--indel-max-len must be in [1, 20] (got {}); this feature targets small indels",
        o.max_len
    );
    anyhow::ensure!(
        o.flank >= 1,
        "--indel-flank must be >= 1 (got {})", o.flank
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
                // Denoise defaults for the combined pipeline: min_strand 2, strand-bias
                // p 0.01, near-homoplasmic exemption at 0.7 (matching the Denoise CLI
                // and call.rs's own permutation/homoplasmic threshold).
                if let Err(e) = denoise::start(
                    &mt_output, &denoised, &reference_path, &data_type,
                    vaf_threshold as f64, 2, 0.01, 0.7,
                    // Indel correction is on in QuickStart. Every other field keeps
                    // its `Default` value, which `validate_indel_opts` accepts as-is,
                    // so there is nothing to validate here.
                    &denoise::indel::IndelOpts { enabled: true, ..Default::default() },
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
            min_strand,
            strand_bias_p,
            homoplasmic_vaf,
            indels,
            indel_max_len,
            indel_vaf,
            indel_err0,
            indel_err_scale,
            indel_err_cap,
            indel_floor_mult,
            indel_delta,
            indel_flank,
            indel_protect_vaf,
            stats,
        } => {
            let iopts = denoise::indel::IndelOpts {
                enabled: indels,
                max_len: indel_max_len,
                vaf: indel_vaf,
                err0: indel_err0,
                err_scale: indel_err_scale,
                err_cap: indel_err_cap,
                floor_mult: indel_floor_mult,
                delta: indel_delta,
                flank: indel_flank,
                protect_vaf: indel_protect_vaf,
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
                &input, &output, &reference, &data_type, vaf, min_strand,
                strand_bias_p, homoplasmic_vaf, &iopts, stats.as_ref(),
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
            min_presence,
            min_absence,
            min_reads,
            data_type,
            fp_rate,
            fn_rate,
            mcmc_iterations,
            mcmc_chains,
            mcmc_seed,
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
                min_presence,
                min_absence,
                min_reads,
                &data_type,
                fp_rate,
                fn_rate,
                mcmc_iterations,
                mcmc_chains,
                mcmc_seed,
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
    use denoise::indel::IndelOpts;

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

    #[test]
    fn err_cap_above_half_is_rejected() {
        // FIX 4: the candidacy floor honours err_cap above 0.5, but
        // `assign_allele` clamps eps to [0.0, 0.5] -- a cap above 0.5 would
        // make the gate and the assignment step disagree.
        let mut o = opts();
        o.err_cap = 0.6;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn err0_above_half_is_rejected() {
        let mut o = opts();
        o.err0 = 0.6;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn err_scale_below_one_is_rejected() {
        let mut o = opts();
        o.err_scale = 0.5;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn floor_mult_negative_is_rejected() {
        let mut o = opts();
        o.floor_mult = -1.0;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn floor_mult_zero_is_accepted() {
        let mut o = opts();
        o.floor_mult = 0.0;
        assert!(validate_indel_opts(&o).is_ok());
    }

    #[test]
    fn delta_zero_is_rejected() {
        let mut o = opts();
        o.delta = 0.0;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn delta_above_one_is_rejected() {
        let mut o = opts();
        o.delta = 1.5;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn vaf_above_one_is_rejected() {
        let mut o = opts();
        o.vaf = 1.5;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn max_len_zero_is_rejected() {
        let mut o = opts();
        o.max_len = 0;
        assert!(validate_indel_opts(&o).is_err());
    }

    #[test]
    fn flank_zero_is_rejected() {
        let mut o = opts();
        o.flank = 0;
        assert!(validate_indel_opts(&o).is_err());
    }
}
