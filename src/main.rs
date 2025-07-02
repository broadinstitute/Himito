use std::path::PathBuf;

use clap::{Parser, Subcommand};
use crate::agg::GraphicalGenome;

mod agg;
mod asm;
mod build;
mod call;
mod filter;
mod methyl;
mod correct;
mod minorhap;

#[derive(Debug, Parser)]
#[clap(name = "Himito")]
#[clap(about = "Analysis of mitochondrial genome using long reads.", long_about = None)]

struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

const DEFAULT_KMER_SIZE: usize = 21;

#[derive(Debug, Subcommand)]

enum Commands {
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
        prob_min: f64,

        /// max fraction to keep a read as mtDNA derived
         #[clap(short, long, value_parser, default_value_t = 0.2)]
        fraction_max_methylation: f64,

    },

    /// Build graph from long-read data in FASTA or Bam file.
    #[clap(arg_required_else_help = true)]
    Build {
        /// Output path for anchor graph.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        kmer_size: usize,

        /// Reference Fasta file, rCRS
        #[clap(short, long, value_parser, required = true)]
        reference_path: PathBuf,

        /// bam or fasta file with reads spanning locus of interest.
        #[clap(short, long, value_parser,required = true)]
        input_read_path: PathBuf,
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

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        k: usize,
        /// max Length to do alignment
        #[clap(short, long, value_parser, default_value_t = 3000)]
        length_max: usize,

        /// minimal allele count for variants
        #[clap(short, long, value_parser, default_value_t = 1)]
        minimal_ac: usize,

        /// minimal heteroplasmic frequency for variants
        #[clap(short, long, value_parser, default_value_t = 0.01)]
        vaf_threshold: f32,

        /// output file name
        #[clap(short, long, value_parser, required = true)]
        output_file: PathBuf,

        /// sample name of the bam file
        #[clap(short, long, value_parser, required = true)]
        sample_id: String,
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
        major_haplotype:bool
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

        /// path for output fasta file
        #[clap(short, long, value_parser)]
        outputfile: PathBuf,

        /// header for the major haplotype, usually the sample name
        #[clap(short, long, value_parser)]
        sample: String,
    },
}

fn main() {
    let args = Cli::parse();
    match args.command {
        Commands::Filter {
            input_bam,
            chromo,
            mt_output,
            numts_output,
            prob_min,
            fraction_max_methylation
        } => {
            let _ = filter::start(&input_bam, &chromo, &mt_output, &numts_output, prob_min, fraction_max_methylation);
        }

        Commands::Build {
            output,
            kmer_size,
            input_read_path,
            reference_path,
        } => {
            build::start(&output, kmer_size, &input_read_path, &reference_path);
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
            length_max,
            minimal_ac,
            output_file,
            sample_id,
            vaf_threshold,
        } => {
            call::start(
                &graphfile,
                &reference_fasta,
                k,
                length_max,
                minimal_ac,
                &output_file,
                &sample_id,
                vaf_threshold,
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
                outputfile,
                sample,
            } => {
                minorhap::start(&graphfile, ref_length, bin_size,  pad_size, min_read_ratio,&outputfile, &sample);
        }
    }

    
}
