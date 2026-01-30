use std::{collections::{HashMap, HashSet}, path::PathBuf, io::Write};
use rust_htslib::bam::{self, Read, Writer,IndexedReader, Header, Record, record::{Aux, AuxArray}};
use bio::io::fasta::{Reader as FastaReader, Record as FastaRecord};



pub fn find_numts(bam_file: &PathBuf, chromo: &str) -> Result<Vec<(String, i32, String, String, i32, String, String)>, Box<dyn std::error::Error>> {
    let mut numts_mapping_info: Vec<(String, i32, String, String, i32, String, String)> = Vec::new();
    let mut bam = IndexedReader::from_path(bam_file)?;

    // Get the chromosome ID from the header
    let tid = bam.header().tid(chromo.as_bytes())
        .ok_or("Chromosome not found in BAM header")?;
    
    // Get the chromosome length
    let chrom_length = bam.header().target_len(tid)
        .ok_or("Could not get chromosome length")?;


    // Set the region to fetch
    bam.fetch((tid, 0, chrom_length))?;
    let header = bam.header().clone();
    
    for read in bam.records() {
        let record = read?;
        // Skip unmapped reads
        if record.is_unmapped() {
            continue;
        }
        // Get the target ID (tid) for the primary alignment
        let primary_tid = record.tid();
        // Convert tid to chromosome name using the header
        let primary_mapping_chromosome = String::from_utf8_lossy(&header.tid2name(primary_tid as u32)).to_string();
        // Get the primary alignment position (0-based)
        let primary_mapping_position = record.pos();
        // Get the primary alignment strand
        let primary_mapping_strand = record.strand().to_string();
        // Check for supplementary alignments
        if let Ok(sa_tag) = record.aux(b"SA") {
            if let Aux::String(sa_str) = sa_tag {
                let supplementary_positions: Vec<&str> = sa_str.split(';')
                    .filter(|s| !s.is_empty())
                    .collect();                
                for sa in supplementary_positions {
                    let parts: Vec<&str> = sa.split(',').collect();
                    // println!("{:?}", parts);
                    if parts.len() >= 6 {
                        let chrom = parts[0];
                        let pos = parts[1].parse::<i32>().unwrap();
                        let strand = parts[2].parse::<char>().unwrap();
                        let cigar = parts[3];
                        if chrom != chromo.to_string() {
                            numts_mapping_info.push((primary_mapping_chromosome.clone(), primary_mapping_position as i32, primary_mapping_strand.to_string(), chrom.to_string(), pos, strand.to_string(), cigar.to_string()));
                        }
                    }
                }
            }
        }

    }
    
    Ok(numts_mapping_info)
}

/// Merges consecutive values that are within max_gap of each other into intervals
/// Returns a vector of (start, end) tuples representing merged intervals
pub fn merge_by_gap(vals: &[(i32, String, String, i32, String)], max_gap: i32) -> Vec<((i32, String, String, i32, String), (i32, String, String, i32, String), i32)> {
    if vals.is_empty() {
        return Vec::new();
    }
    
    // Sort and deduplicate (equivalent to sorted(set(vals)))
    let mut sorted_vals: Vec<(i32, String, String, i32, String)> = vals.to_vec();
    sorted_vals.sort_by_key(|k| k.0);
    sorted_vals.dedup();
    
    let mut intervals = Vec::new();
    let mut start = sorted_vals[0].clone();
    let mut prev = sorted_vals[0].clone();
    
    // Iterate through sorted values starting from the second one
    let mut count = 0;
    for x in sorted_vals.clone().iter().skip(1) {
        if x.0 - prev.0 <= max_gap {
            // Extend the current interval
            prev = x.clone();
            count += 1;
        } else {
            // Save the current interval and start a new one
            intervals.push((start.clone(), prev.clone(), count+1));
            start = x.clone();
            prev = x.clone();
            count = 0;
        }
    }
    intervals.push((start, prev, count + 1));
    
    intervals
}
/// Groups breakpoints by chromosome and merges positions within max_gap_threshold
/// Returns intervals as (chromosome, start, end) tuples
pub fn get_numts_intervals(
    numts_mapping_info: &[(String, i32, String, String, i32, String, String)], 
    max_gap_threshold: i32
) -> Vec<(String, (i32, String, String, i32, String), (i32, String, String, i32, String), i32)> {
    let mut numts_intervals: Vec<(String, (i32, String, String, i32, String), (i32, String, String, i32, String), i32)> = Vec::new();
    let mut stack: HashMap<String, Vec<(i32, String, String, i32, String)>> = HashMap::new();
    
    // Group positions by chromosome
    for (_mt_chrom, _mt_start, _mt_strand, numt_chrom, numt_start, _numt_strand, _cigar) in numts_mapping_info {
        stack.entry(numt_chrom.clone())
            .or_insert_with(Vec::new)
            .push((*numt_start as i32,_numt_strand.clone(), _mt_chrom.clone(), *_mt_start, _mt_strand.clone()));
    }
    
    // Merge positions by gap for each chromosome and create intervals
    for (chrom, positionlist) in stack {
        let merged = merge_by_gap(&positionlist, max_gap_threshold);
        for (start, end, counts) in merged {
            numts_intervals.push((chrom.clone(), start, end, counts));
        }
    }
    numts_intervals
}

/// Formats BND ALT field based on strand orientations
/// Returns the BND ALT string in VCF format
fn format_bnd_alt(
    ref_base: &str,
    mt_chrom: &str,
    mt_pos: i32,
    mt_strand: &str,
    sa_chrom: &str,
    sa_pos: i32,
    sa_strand: &str,
) -> String {
    // BND format: 
    // - N[chr2:pos2[ means N connects to chr2:pos2 from left (both forward)
    // - N]chr2:pos2] means N connects to chr2:pos2 from right (primary forward, SA reverse)
    // - ]chr2:pos2]N means chr2:pos2 connects to N from left (primary reverse, SA forward)
    // - [chr2:pos2[N means chr2:pos2 connects to N from right (both reverse)
    
    let primary_forward = sa_strand == "+";
    let mt_forward = mt_strand == "+";
    
    if primary_forward && mt_forward {
        // Both forward: N[chr2:pos2[
        format!("{}[{}:{}[", ref_base, mt_chrom, mt_pos)
    } else if primary_forward && !mt_forward {
        // Primary forward, SA reverse: N]chr2:pos2]
        format!("{}]{}:{}]", ref_base, mt_chrom, mt_pos)
    } else if !primary_forward && mt_forward {
        // Primary reverse, SA forward: ]chr2:pos2]N
        format!("]{}:{}]{}", mt_chrom, mt_pos, ref_base)
    } else {
        // Both reverse: [chr2:pos2[N
        format!("[{}:{}[{}", mt_chrom, mt_pos, ref_base)
    }
}

/// Writes BND records to a VCF file from SA tag information
pub fn write_bnd_vcf(
    output_vcf: &PathBuf,
    numts_intervals: &Vec<(String, (i32, String, String, i32, String), (i32, String, String, i32, String), i32)>,
    reference_file: &PathBuf,
    sample_name: &str,
    minimal_ac: i32,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    let mut vcf_file = File::create(output_vcf)?;
    
    // reference fasta information
    let ref_reader = FastaReader::from_file(reference_file).unwrap();
    let reference_sequence: Vec<FastaRecord> = ref_reader.records().map(|r| r.unwrap()).collect();
    let ref_seq = String::from_utf8_lossy(reference_sequence[0].seq()).to_string();
    let ref_header = reference_sequence[0].id().to_string();
    
    // Write VCF header
    writeln!(vcf_file, "##fileformat=VCFv4.3")?;
    for record in reference_sequence.iter() {
        let referencename = record.id().to_string();
        let referencelength = record.seq().len() as u64;
        writeln!(vcf_file, "##contig=<ID={},length={}>", referencename, referencelength)?;
    }
    writeln!(vcf_file, "##ALT=<ID=BND,Description=\"NUMTs breakpoints\">")?;
    writeln!(vcf_file, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")?;
    writeln!(vcf_file,  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    writeln!(vcf_file,   "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Supporting Read Counts\">")?;
    writeln!(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample_name)?;
    
    // Write BND records
    for (auto_chrom, (auto_start, auto_strand, mt_chrom, mt_pos, mt_strand), (auto_end, auto_end_strand, mt_chrom_end, mt_pos_end, mt_strand_end), counts) in numts_intervals {
        // VCF uses 1-based positions
        let pos = auto_start + 1;
        let ref_seq = reference_sequence
        .iter()
        .find(|r| r.id() == auto_chrom)
        .unwrap()
        .seq();
        let ref_seq_ = String::from_utf8_lossy(ref_seq).to_string();
        let ref_base = ref_seq_[(pos - 1) as usize..pos as usize].to_string();
        // Format BND ALT field
        let alt = format_bnd_alt(
            &ref_base,
            mt_chrom,
            *mt_pos,
            mt_strand,
            auto_chrom,
            *auto_start,
            auto_strand,
        );
        
        // Create variant ID
        let variant_id = format!("BND_{}_{}_{}_{}_{}_{}", auto_chrom, auto_start, auto_strand, mt_chrom, mt_pos, mt_strand);
        
        // INFO field
        let info = "SVTYPE=BND".to_string();
        
        // Write the record
        if counts >= &minimal_ac {
            writeln!(
                vcf_file,
                "{}\t{}\t{}\t{}\t{}\t.\t.\t{}\tGT:AD\t0/1:{}",
                auto_chrom, pos, ref_base, variant_id, alt, info, counts
            )?;
        }
        
    }
    
    Ok(())
}

pub fn start(input_bam:&PathBuf, chromo: &str, max_gap_threshold: i32, output_vcf:&PathBuf, reference_file:&PathBuf, sample_name:&str, minimal_ac: i32) -> Result<(), Box<dyn std::error::Error>> {
    
    let numts_mapping =  find_numts(input_bam, chromo).expect("Error finding numts");
    let numts_intervals = get_numts_intervals(&numts_mapping, max_gap_threshold);

    write_bnd_vcf(&output_vcf, &numts_intervals, &reference_file, &sample_name, minimal_ac).expect("Error writing BND VCF");
    
    Ok(())
}