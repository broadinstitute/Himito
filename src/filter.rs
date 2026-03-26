
use std::{collections::HashSet, path::PathBuf};
use rust_htslib::bam::{self, Read,IndexedReader, Header, record::Aux};
use rand::seq::IteratorRandom;


pub fn downsample_list(mt_readnames: &HashSet<String>, max_reads:usize) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    if mt_readnames.len() <= max_reads {
        return Ok(mt_readnames.iter().cloned().collect());
    }
    let mut rng = rand::thread_rng();
    let downsampled_mt_readnames: Vec<String> = mt_readnames
        .iter()
        .choose_multiple(&mut rng, max_reads)
        .into_iter()
        .cloned()
        .collect();
    Ok(downsampled_mt_readnames)

}


pub fn find_numts(bam_file: &PathBuf, chromo: &str, mod_char:char, min_methyl_prob:f64, fraction_threshold:f64) -> Result<(HashSet<String>, HashSet<String>), Box<dyn std::error::Error>> {
    let mut numts_readnames: HashSet<String> = HashSet::new();
    let mut mt_readnames: HashSet<String> = HashSet::new();

    let mut bam = IndexedReader::from_path(bam_file).unwrap();

    // Get the chromosome ID from the header
    let tid = bam.header().tid(chromo.as_bytes())
        .ok_or("Chromosome not found in BAM header")?;
    
    // Get the chromosome length
    let chrom_length = bam.header().target_len(tid)
        .ok_or("Could not get chromosome length")?;

    println!("{},{},{}", chromo, tid, chrom_length);

    // Set the region to fetch
    bam.fetch((tid, 0, chrom_length))?;
    
    for read in bam.records() {
        let record = read?;
        let mut is_mt = true;
        
        // Skip unmapped reads
        if record.is_unmapped() {
            continue;
        }
        
        // Convert query name to String
        let query_name = String::from_utf8_lossy(record.qname()).to_string();
        
        // Check for secondary alignments
        if record.is_secondary() {
            continue;
        }
        
        // Check for supplementary alignments
        if let Ok(sa_tag) = record.aux(b"SA") {
            if let Aux::String(sa_str) = sa_tag {
                let supplementary_positions: Vec<&str> = sa_str.split(';')
                    .filter(|s| !s.is_empty())
                    .collect();
                
                for sa in supplementary_positions {
                    let parts: Vec<&str> = sa.split(',').collect();
                    if parts.len() >= 6 {
                        let chrom = parts[0];
                        if chrom != chromo.to_string() {
                            // println!("{},{}", chrom, chromo);
                            numts_readnames.insert(query_name.clone());
                            is_mt = false;
                            break;
                        }
                    }
                }
            }
        }

        // Check for methylation signals
        if record.aux(b"Mm").is_ok() || record.aux(b"Ml").is_ok() || record.aux(b"MM").is_ok() || record.aux(b"ML").is_ok() {
            let mut methylation_signal = Vec::new();
            if let Ok(mods) = record.clone().basemods_iter() {
                // Iterate over the modification types
                for res in mods {
                    if let Ok( (position, m) ) = res {
                        if m.modified_base as u8 as char != mod_char{
                            continue
                        }
                        // let strand = mod_metadata.strand;
                        let qual = m.qual as f32 / 255.0;
                        methylation_signal.push(qual);

                    }
                }                    
            }

            if methylation_signal.len() == 0{
                mt_readnames.insert(query_name.clone());
                continue
            }

            let count = methylation_signal.iter().filter(|&&q| q > min_methyl_prob as f32).count();
            let methylated_fraction: f64 = count as f64 / methylation_signal.len() as f64;
            if methylated_fraction > fraction_threshold {
                numts_readnames.insert(query_name.clone());
                is_mt = false;
            }
        }
        if is_mt {
            mt_readnames.insert(query_name.clone());
        }
    }
    
    Ok((numts_readnames, mt_readnames))
}

fn write_bams(
    bam_file: &PathBuf, 
    mt_bam: &PathBuf, 
    numts_bam: &PathBuf, 
    max_reads:usize,
    chromo: &str,
    numts_readnames: &HashSet<String>,
    mt_readnames: &HashSet<String>
) -> Result<(), Box<dyn std::error::Error>> {
    // Open input BAM
    let mut bam = IndexedReader::from_path(bam_file).unwrap();

    // Get the header from input BAM
    let header = Header::from_template(bam.header());
    
    // Create output BAM writers
    let mut mt_out = bam::Writer::from_path(mt_bam, &header, bam::Format::Bam)?;
    let mut numt_out = bam::Writer::from_path(numts_bam, &header, bam::Format::Bam)?;
    
    // Get chromosome ID
    let tid = bam.header().tid(chromo.as_bytes())
        .ok_or("Chromosome not found in BAM header")?;
    
    // Set region to fetch
    bam.fetch((tid, 0, 17000))?;

    let downsampled_mt_readnames = downsample_list(mt_readnames, max_reads)?;

    
    for read in bam.records() {
        let record = read?;
        let query_name = String::from_utf8_lossy(record.qname()).to_string();
        
        if numts_readnames.contains(&query_name) {
            numt_out.write(&record)?;
        } else if downsampled_mt_readnames.contains(&query_name) {
            mt_out.write(&record)?;
        }
    }
    Ok(())
}


// Example usage
pub fn start(input_bam:&PathBuf, chromo: &str, mt_output:&PathBuf, numts_output:&PathBuf, min_prob:f64, fraction_max_methylation:f64, max_mt_reads:usize ) -> Result<(), Box<dyn std::error::Error>> {
    
    
    match find_numts(input_bam, chromo, 'm', min_prob, fraction_max_methylation) {
        Ok((numts, mts)) => {
            println!("Found {} potential NUMTS reads", numts.len());
            // Then split into separate BAM files
            write_bams(input_bam, mt_output, numts_output, max_mt_reads, chromo, &numts, &mts)?;
        },
        Err(e) => eprintln!("Error: {}", e),
    }
    
    Ok(())
}