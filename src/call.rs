use crate::agg;
use agg::*;
use crate::build;
use indicatif::ProgressBar;
use rayon::prelude::*;
use std::collections::{HashMap};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::{collections::HashSet, path::PathBuf};
use std::error::Error;
use csv::Writer;
use ndarray::{Array1, Array2, Axis, s};
use rand::seq::SliceRandom;
use rand::thread_rng;
use statrs::distribution::{Normal, ContinuousCDF, Binomial, DiscreteCDF};
use regex::Regex;
use adjustp::{adjust, Procedure};
use bio::io::fasta::{Reader, Record};
use rust_htslib::bam::{Read as HtsRead, Reader as BamReader};

#[derive(Debug, Clone)]
pub struct Variant {
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_type: String,
    pub allele_count: usize,
}

pub fn get_variants_from_cigar(
    cigar: &str,
    ref_seq: &str,
    alt_seq: &str,
    ref_start: usize,
    allelecount: usize,
) -> (Vec<Variant>, HashMap<usize, usize>) {
    let mut poscount = HashMap::new();
    let mut variants = Vec::new();
    let mut ref_pos = 0;
    let mut alt_pos = 0;

    let mut operations: Vec<(usize, char)> = Vec::new();
    let mut num = String::new();

    for c in cigar.trim_matches('"').chars() {
        // trim_matches removes quotes at start and end
        if c.is_digit(10) {
            num.push(c);
        } else {
            let number = num.parse::<usize>().expect("number {}");
            operations.push((number, c));
            num.clear();
        }
    }

    for (length, op) in operations {
        match op {
            '=' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                }
                ref_pos += length;
                alt_pos += length;
            }
            'X' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                    let ref_allele = &ref_seq[ref_pos + i..ref_pos + i + 1];
                    let alt_allele = &alt_seq[alt_pos + i..alt_pos + i + 1];
                    variants.push(Variant {
                        pos,
                        ref_allele: ref_allele.to_string(),
                        alt_allele: alt_allele.to_string(),
                        variant_type: "SNP".to_string(),
                        allele_count: allelecount,
                    });
                }
                ref_pos += length;
                alt_pos += length;
            }
            'I' => {
                let pos = ref_start + ref_pos;
                // *poscount.entry(pos).or_insert(0) += allelecount;
                let ref_allele = if ref_pos > 0 {
                    match ref_seq.get(ref_pos - 1..ref_pos) {
                        Some(allele) => allele,
                        None => {
                            println!("{} {} {}", ref_seq, ref_seq.len(), ref_pos);
                            "-"
                        }
                    }
                } else {
                    "-"
                };

                let alt_allele = match alt_seq.get(alt_pos - 1..alt_pos + length) {
                    Some(allele) => allele,
                    None => {
                        println!("{} {} {}", alt_seq, alt_seq.len(), alt_pos);
                        "-"
                    }
                };
                variants.push(Variant {
                    pos,
                    ref_allele: ref_allele.to_string(),
                    alt_allele: alt_allele.to_string(),
                    variant_type: "INS".to_string(),
                    allele_count: allelecount,
                });
                alt_pos += length;
            }

            'D' => {
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    *poscount.entry(pos).or_insert(0) += allelecount;
                }
                let pos = ref_start + ref_pos;
                let ref_allele = match ref_seq.get(ref_pos - 1..ref_pos + length) {
                    Some(allele) => allele,
                    None => {
                        println!("{} {} {}", ref_seq, ref_seq.len(), ref_pos);
                        "-"
                    }
                };

                let alt_allele = if alt_pos > 0 {
                    match alt_seq.get(alt_pos - 1..alt_pos) {
                        Some(allele) => allele,
                        None => {
                            println!("{} {} {}", alt_seq, alt_seq.len(), alt_pos);
                            "-"
                        }
                    }
                } else {
                    "-"
                };

                if ref_pos > 0 && alt_pos > 0 {
                    if let (Some(r), Some(a)) = (
                        ref_seq.get(ref_pos - 1..ref_pos),
                        alt_seq.get(alt_pos - 1..alt_pos),
                    ) {
                        if r != a {
                            println!("{} {} {}", r, a, cigar);
                        }
                    }
                }

                variants.push(Variant {
                    pos,
                    ref_allele: ref_allele.to_string(),
                    alt_allele: alt_allele.to_string(),
                    variant_type: "DEL".to_string(),
                    allele_count: allelecount,
                });
                ref_pos += length;
            }
            _ => (), //others skip
        }
    }
    (variants, poscount)
}
fn circuliarize_variants(variants: Vec<Variant>, ref_length: usize) -> Vec<Variant> {
    let mut circular_variants = Vec::new();
    for variant in variants {
        if variant.pos < ref_length {
            circular_variants.push(variant);
        }else{
            circular_variants.push(Variant {
                pos: variant.pos - ref_length,
                ref_allele: variant.ref_allele,
                alt_allele: variant.alt_allele,
                variant_type: variant.variant_type,
                allele_count: variant.allele_count,
            });
        }
    }
    circular_variants
}

pub fn get_variant(
    graph: &mut GraphicalGenome,
    k: usize,
    ref_name: &str,
    ref_length: usize,
) -> (Vec<Variant>, HashMap<usize, usize>, HashMap<String, Vec<serde_json::Value> > ) {
    let mut coverage = HashMap::new();
    let mut read_record = HashMap::new();
    let mut var = Vec::new();
    let mut edgelist: Vec<_> = graph.edges.keys().collect();
    edgelist.sort();
    for edge in edgelist {
        let allele_count = graph
            .edges
            .get(edge)
            .and_then(|e| e.get("reads"))
            .map_or(Vec::new(), |reads| {
                reads.as_array().unwrap_or(&Vec::new()).to_vec()
            })
            .len();

        let cigar = graph
            .edges
            .get(edge)
            .and_then(|e| e.get("variants"))
            .cloned()
            .unwrap_or_else(|| serde_json::Value::String("".to_string()));

        if cigar.as_str().unwrap_or("").is_empty() {
            continue;
        }

        let src = graph
            .incoming
            .get(edge)
            .and_then(|edges| edges.first())
            .expect("Edge should have source node");
        let dst = graph
            .outgoing
            .get(edge)
            .and_then(|edge| edge.first())
            .expect("Edge should have dst");
        if src == "SOURCE" || dst == "SINK" {
            continue;
        }
        let refstart = graph.anchor[src]["pos"]
            .as_i64()
            .expect("Position should be a integer");

        let ref_seq = build::find_ref_edge(&graph, src, dst, ref_name, k);
        let src_seq = graph
            .anchor
            .get(src)
            .and_then(|anchor| anchor.get("seq").and_then(|seq| seq.as_str()))
            .unwrap_or("");
        let alt_sequence = src_seq.to_string()
            + &graph
                .edges
                .get(edge)
                .expect("edge not found")
                .get("seq")
                .and_then(|v| v.as_str())
                .unwrap_or("");
        let (variants, poscounts) = get_variants_from_cigar(
            &cigar.to_string(),
            &ref_seq,
            &alt_sequence,
            refstart as usize,
            allele_count,
        );
        let variants_circular = circuliarize_variants(variants.clone(), ref_length);
        var.extend(variants_circular.clone());

        for (pos, count) in poscounts.iter() {
            if pos < &ref_length {
                *coverage.entry(*pos).or_insert(0) += count;
            }else{
                *coverage.entry(*pos - &ref_length).or_insert(0) += count;
            }
        }
        let readlist = graph
            .edges
            .get(edge)
            .and_then(|e| e.get("reads"))
            .map_or(Vec::new(), |reads| {
                reads.as_array().unwrap_or(&Vec::new()).to_vec()
            });
        for v in &variants_circular{
            let key = generate_variant_name(&v.clone());
            read_record.entry(key).or_insert_with(Vec::new).extend(readlist.clone());
        }
            
    }
    (var, coverage, read_record)
}

fn collapse_identical_records(variants: Vec<Variant>, ref_length: usize) -> Vec<Variant> {
    if variants.is_empty() {
        return Vec::new();
    }

    let mut collapsed = HashMap::new();

    for current_var in variants {
        let pos = if current_var.pos < ref_length {current_var.pos} else {current_var.pos - ref_length};
        let ref_allele = current_var.ref_allele;
        let alt_allele = current_var.alt_allele;
        let variant_type = current_var.variant_type;
        let allele_count = current_var.allele_count;

        let key = (
            pos,
            ref_allele.clone(),
            alt_allele.clone(),
            variant_type.clone(),
        );
        *collapsed.entry(key).or_insert(0) += allele_count;
    }
    collapsed
        .into_iter()
        .map(
            |((pos, ref_allele, alt_allele, variant_type), allele_count)| Variant {
                pos,
                ref_allele,
                alt_allele,
                variant_type,
                allele_count,
            },
        )
        .collect()
}

fn format_vcf_record(variant: &Variant, coverage: HashMap<usize, usize>, indel_false_threshold: f64) -> String {
    // Add AC (allele count) to INFO field
    let read_depth = coverage.get(&variant.pos).unwrap_or(&0);
    let allele_frequency = if *read_depth == 0 {
        0.0
    } else {
        variant.allele_count as f32 / *read_depth as f32
    };

    let info = format!("DP={}", read_depth);
    let format: String = format!("GT:AD:HF");
    let genotype: String = format!("1");
    let sample: String = format!("{}:{}:{}", genotype, variant.allele_count, allele_frequency);

    let indel_len = (variant.ref_allele.len() as i32 - variant.alt_allele.len() as i32).abs();
    let is_small_indel = variant.variant_type != "SNP" && indel_len < 5;
    let filter = if is_small_indel && allele_frequency < indel_false_threshold as f32 {
        "Potential_Artifact"
    } else {
        "PASS"
    };

    match variant.variant_type.as_str() {
        "SNP" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t{}\t{}\t{}\t{}",
            variant.pos + 1,
            variant.ref_allele,
            variant.alt_allele,
            filter,
            info,
            format,
            sample
        ),
        "INS" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t{}\t{}\t{}\t{}",
            variant.pos, variant.ref_allele, variant.alt_allele, filter, info, format, sample
        ),
        "DEL" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t{}\t{}\t{}\t{}",
            variant.pos, variant.ref_allele, variant.alt_allele, filter, info, format, sample
        ),
        _ => panic!("Unknown variant type"),
    }
}

fn filter_vcf_record(
    variants: &[Variant],
    coverage: &HashMap<usize, usize>,
    minimal_ac: usize,
    hf_threshold: f32,
) -> Vec<Variant> {
    let mut filtered_var = Vec::new();
    for variant in variants{
        let allele_count = variant.allele_count;
        if allele_count < minimal_ac + 1 {
            continue;
        }
        let read_depth = coverage.get(&variant.pos).unwrap_or(&0);
        let hf = if *read_depth == 0 {
            0.0
        } else {
            variant.allele_count as f32 / *read_depth as f32
        };
        if hf < hf_threshold as f32 {
            continue;
        }

        // remove reference Ns
        let ref_allele = variant.ref_allele.as_str();
        if ref_allele.contains("N") {
            continue;
        }
        
        filtered_var.push(variant.clone());
    }
    filtered_var
}

fn write_vcf(
    variants: &[Variant],
    coverage: &HashMap<usize, usize>,
    output_file: &PathBuf,
    sample_id: &str,
    referencename:&str,
    referencelength:usize,
    indel_false_threshold: f64
) -> std::io::Result<()> {
    let mut file = File::create(Path::new(output_file))?;

    // Write VCF header
    writeln!(file, "##fileformat=VCFv4.2")?;
    writeln!(file, "##reference={}", referencename)?;
    writeln!(file, "##contig=<ID={},length={}>", referencename, referencelength)?;
    writeln!(
        file,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
    )?;
    writeln!(
        file,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        file,
        "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">"
    )?;
    writeln!(
        file,
        "##FORMAT=<ID=HF,Number=1,Type=Float,Description=\"Heteroplasmic Frequency\">"
    )?;
    writeln!(file, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    writeln!(
        file,
        "##FILTER=<ID=Potential_Artifact,Description=\"Small indel (length < 5 bp) with heteroplasmic frequency < indel_false_threshold\">"
    )?;
    writeln!(
        file,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        sample_id
    )?;

    // Sort variants by position
    let mut sorted_variants = variants.to_vec();
    sorted_variants.sort_by(|a, b| {
        // First compare positions
        a.pos
            .cmp(&b.pos)
            // Then compare variant types (to ensure consistent ordering)
            .then(a.variant_type.cmp(&b.variant_type))
            // Then compare ref alleles
            .then(a.ref_allele.cmp(&b.ref_allele))
            // Then compare alt alleles
            .then(a.alt_allele.cmp(&b.alt_allele))
    });

    // Write variant records
    for variant in sorted_variants {

        writeln!(file, "{}", format_vcf_record(&variant, coverage.clone(), indel_false_threshold))?;
    }

    Ok(())
}

pub fn generate_variant_name(variant: &Variant) -> String {
    if variant.variant_type == "SNP" {
        format!("m.{}{}>{}", variant.pos + 1, variant.ref_allele, variant.alt_allele)
    } else {
        format!("m.{}{}>{}", variant.pos, variant.ref_allele, variant.alt_allele)
    }
}

pub fn construct_matrix (read_record:&HashMap<String, Vec<serde_json::Value>>, variants:&[Variant], minimal_ac: usize) -> (Array2<f64>, Vec<Variant>, Vec<String>) {
    let mut read_set: HashSet<String> = HashSet::new();
    for (_, readlist) in read_record {
        for read in readlist {
            if let Some(read_str) = read.as_str() {
                read_set.insert(read_str.to_string());
            }
        }
    }
    let mut read_vec: Vec<String> = read_set.into_iter().collect();
    read_vec.sort();

    let read_set_dict: HashMap<String, usize> = read_vec
        .iter()
        .enumerate()
        .map(|(i, read)| (read.clone(), i))
        .collect();

    let mut var_vec: Vec<Variant> = variants.to_vec();
    var_vec.sort_by(|a, b| {
        a.pos
            .cmp(&b.pos)
            .then(a.ref_allele.cmp(&b.ref_allele))
            .then(a.alt_allele.cmp(&b.alt_allele))
    });

    let var_record_dict: HashMap<String, usize> = var_vec
        .iter()
        .enumerate()
        .map(|(i, var)| (generate_variant_name(var), i))
        .collect();


    // Create a 2D matrix filled with zeros
    let mut matrix = Array2::<f64>::zeros((var_vec.len(), read_vec.len()));

    for var in &var_vec {
        let readlist: HashSet<String> = read_record.get(&generate_variant_name(var))
            .map(|reads| {
                reads.iter()
                    .filter_map(|read| read.as_str().map(|s| s.to_string()))
                    .collect()
            })
            .unwrap_or_else(|| HashSet::new());
        
        // Get row index for this variant
        let r_index = *var_record_dict.get(&generate_variant_name(var)).unwrap();
        
        // For each read in the readlist
        for read in readlist {
            // Get column index for this read
            if let Some(c_index) = read_set_dict.get(&read) {
                // Increment the matrix cell
                matrix[[r_index, *c_index]] += 1.0;
            }
        }
    }

    // filter matrix where the sum of the column is less than 2
    let mut col_index = Vec::new();
    for i in 0..matrix.ncols() {
        if matrix.column(i).sum() > minimal_ac as f64 {
            col_index.push(i);
        }
    }
    let filtered_matrix = matrix.select(Axis(1), &col_index);
    let filtered_read_vec: Vec<String> = read_vec.iter().filter(|read| col_index.contains(&read_set_dict.get(read.as_str()).unwrap())).map(|read| read.clone()).collect();
    (filtered_matrix, var_vec, filtered_read_vec)

}



fn write_matrix_to_csv<P: AsRef<Path>>(
    matrix: &Array2<f64>,
    var_record: &[Variant],
    read_set: &[String],
    path: P
) -> Result<(), Box<dyn Error>> {
    // Create a file and CSV writer
    let file = File::create(path)?;
    let mut writer = Writer::from_writer(file);
    
    // Prepare header row (with empty cell for the corner)
    let mut header = vec!["variant".to_string()];
    header.extend(read_set.iter().cloned());
    
    // Write header
    writer.write_record(&header)?;
    
    // Write each row with its row name
    for (row_idx, variants) in var_record.iter().enumerate() {
        let var_name = generate_variant_name(variants);
        let mut row = vec![var_name.clone()];
        
        // Add the values from the matrix
        for col_idx in 0..matrix.ncols() {
            row.push(matrix[[row_idx, col_idx]].to_string());
        }
        
        writer.write_record(&row)?;
    }
    
    // Flush and finish
    writer.flush()?;
    Ok(())
}

fn jaccard_distance(vector1: &[bool], vector2: &[bool]) -> f64 {
    assert_eq!(vector1.len(), vector2.len(), "Vectors must have the same length");
    
    let mut intersection_count = 0;
    let mut union_count = 0;
    
    for (a, b) in vector1.iter().zip(vector2.iter()) {
        if *a && *b {
            intersection_count += 1;
        }
        if *a || *b {
            union_count += 1;
        }
    }
    
    if union_count == 0 {
        return 0.0; // Both vectors are all zeros
    }
    
    1.0 - (intersection_count as f64 / union_count as f64)
}

/// Generate a null distribution through permutation testing
fn get_null_distribution(
    filtered_var: &Vec<Variant>,
    matrix: &Array2<f64>, 
    permutation_round: usize,
    threshold: f64
) -> Vec<f64> {

    let bar = ProgressBar::new(permutation_round as u64);
    let summary_statistics = (0..permutation_round).into_par_iter().flat_map(|_|  {
        bar.inc(1);
        let mut local_stats = Vec::new();
        let mut rng = thread_rng();

        for (i, variant) in filtered_var.iter().enumerate() {
            let index = generate_variant_name(variant);
            let vector = matrix.slice(s![i, ..]);
            
            // Skip vectors with frequency > threshold
            let frequency = vector.sum() / vector.len() as f64;
            if frequency > threshold {
                continue;
            }
        
            // Create a shuffled copy of the vector
            let vector_data: Vec<f64> = vector.iter().copied().collect();
            let mut shuffled_data = vector_data.clone();
            shuffled_data.shuffle(&mut rng);
            let shuffled = Array1::from(shuffled_data);

            let mut all_coefficients = Vec::new();
            
            for (j, other_variant) in filtered_var.iter().enumerate() {
                let other_index = generate_variant_name(other_variant);
                if index == other_index {
                    continue;
                }

                // let other_frequency = other_variant.allele_count as f64 / *coverage.get(&other_variant.pos).unwrap() as f64;
                // if other_frequency > threshold {
                //     continue;
                // }

                let other_vector = matrix.slice(s![j, ..]);
                
                let binary_vector: Vec<bool> = shuffled.iter().map(|&x| x > 0.5).collect();
                let binary_other: Vec<bool> = other_vector.iter().map(|&x| x > 0.5).collect();

                // Calculate Jaccard distance
                let coor = (1.0 - jaccard_distance(&binary_vector, &binary_other)).abs();
                all_coefficients.push(coor);
            }
           
            local_stats.push(all_coefficients.iter().sum());
        }
        local_stats
    }).collect::<Vec<f64>>();
    bar.finish();
    summary_statistics
    
}

/// Calculate statistics for observed data
fn calculate_observation_statistics(
    filtered_var: &Vec<Variant>,
    index: usize,
    matrix: &Array2<f64>, 
) -> f64 {

    let vector = &matrix.slice(s![index, ..]);
    let mut all_coefficients = Vec::new();
    
    for (i, variant) in filtered_var.iter().enumerate() {
        if i == index {
            continue;
        }
        // let other_frequency = variant.allele_count as f64 / *coverage.get(&variant.pos).unwrap() as f64;
        // if other_frequency > threshold {
        //     continue;
        // }
        
        let other_vector =&matrix.slice(s![i, ..]);
        
        // Convert arrays to binary vectors before calculating Jaccard distance
        let binary_vector: Vec<bool> = vector.iter().map(|&x| x > 0.5).collect();
        let binary_other: Vec<bool> = other_vector.iter().map(|&x| x > 0.5).collect();

        // Calculate Jaccard distance
        let coor = (1.0 - jaccard_distance(&binary_vector, &binary_other)).abs();
        all_coefficients.push(coor);
    }
    
    all_coefficients.iter().sum()
}

/// Calculate p-value using z-score approach
fn calculate_p_value(statistics: &[f64], observation: f64) -> f64 {
    let n = statistics.len() as f64;
    
    // Calculate mean
    let mu = statistics.iter().sum::<f64>() / n;
    
    // Calculate standard deviation
    let variance = statistics.iter()
        .map(|&x| (x - mu).powi(2))
        .sum::<f64>() / n;
    let sigma = variance.sqrt();
    // println!("{:?}, {}", statistics, observation);
    
    let z_score = (observation - mu) / sigma;
    
    // Calculate p-value using normal distribution CDF
    let normal = Normal::new(0.0, 1.0).unwrap();
    if z_score.is_nan() {
        return 1.0;
    }
    return 1.0 - normal.cdf(z_score);
}

/// Build a map from read name (BAM qname) to strand (true = reverse, false = forward),
/// using primary, mapped alignments only. Also returns the library-wide forward-strand
/// fraction, which is the expected forward proportion for an unbiased allele.
fn build_strand_map(bam_file: &PathBuf) -> (HashMap<String, bool>, f64) {
    let mut bam = BamReader::from_path(bam_file)
        .unwrap_or_else(|e| panic!("Failed to open BAM {:?}: {}", bam_file, e));
    let mut strand_map = HashMap::new();
    let (mut fwd, mut rev) = (0usize, 0usize);

    for result in bam.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        if record.is_secondary() || record.is_supplementary() || record.is_unmapped() {
            continue;
        }
        let name = String::from_utf8_lossy(record.qname()).to_string();
        let is_reverse = record.is_reverse();
        if is_reverse { rev += 1; } else { fwd += 1; }
        strand_map.insert(name, is_reverse);
    }

    let total = fwd + rev;
    let fwd_frac = if total == 0 { 0.5 } else { fwd as f64 / total as f64 };
    (strand_map, fwd_frac)
}

/// Two-sided binomial-tail p-value that a variant's alt-supporting reads' forward/reverse
/// split deviates from the library's expected forward fraction. A small p-value means the
/// allele is strand-skewed, which is characteristic of systematic (e.g. homopolymer) errors
/// rather than true heteroplasmies, which should be roughly strand-balanced.
fn strand_bias_pvalue(fwd: usize, rev: usize, expected_fwd_frac: f64) -> f64 {
    let n = (fwd + rev) as u64;
    if n == 0 {
        return 1.0;
    }
    // clamp away from 0/1 so the distribution is well defined
    let p = expected_fwd_frac.clamp(1e-6, 1.0 - 1e-6);
    let binom = Binomial::new(p, n).unwrap();
    let k = fwd as u64;
    let lower = binom.cdf(k);                                  // P(X <= k)
    let upper = if k == 0 { 1.0 } else { binom.sf(k - 1) };    // P(X >= k)
    (2.0 * lower.min(upper)).min(1.0)
}

/// Remove variants whose alt-supporting reads are significantly strand-skewed relative to
/// the library's forward/reverse composition. Variants with fewer than `min_reads`
/// strand-resolved supporting reads are kept (too little evidence to judge). Returns the
/// kept variants together with the matrix subset to the corresponding rows.
fn filter_strand_bias(
    variants: &[Variant],
    matrix: &Array2<f64>,
    read_record: &HashMap<String, Vec<serde_json::Value>>,
    strand_map: &HashMap<String, bool>,
    expected_fwd_frac: f64,
    p_threshold: f64,
    min_reads: usize,
) -> (Vec<Variant>, Array2<f64>) {
    let mut keep_idx = Vec::new();
    let mut kept = Vec::new();

    for (i, variant) in variants.iter().enumerate() {
        let name = generate_variant_name(variant);
        let (mut fwd, mut rev) = (0usize, 0usize);
        let mut seen = HashSet::new();

        if let Some(reads) = read_record.get(&name) {
            for read in reads {
                if let Some(read_name) = read.as_str() {
                    // count each supporting read once
                    if !seen.insert(read_name.to_string()) {
                        continue;
                    }
                    match strand_map.get(read_name) {
                        Some(true) => rev += 1,
                        Some(false) => fwd += 1,
                        None => {} // read absent from BAM (e.g. non-primary alignment)
                    }
                }
            }
        }

        let total = fwd + rev;
        if total < min_reads {
            // not enough strand-resolved reads to make a call: keep
            keep_idx.push(i);
            kept.push(variant.clone());
            continue;
        }

        let p_value = strand_bias_pvalue(fwd, rev, expected_fwd_frac);
        if p_value < p_threshold {
            println!(
                "Strand-biased (removed), {:?}, fwd={}, rev={}, p={:.3e}",
                name, fwd, rev, p_value
            );
        } else {
            keep_idx.push(i);
            kept.push(variant.clone());
        }
    }

    let filtered_matrix = matrix.select(Axis(0), &keep_idx);
    (kept, filtered_matrix)
}

fn permutation_test(
    matrix: &Array2<f64>,
    p_value_threshold: f64,
    permutation_round: usize,
    filtered_var: &Vec<Variant>,
    coverage: &HashMap<usize, usize>,
    permutation_frequency_threshold:f64,
    data_type: &str
) -> (Vec<Variant>, Array2<f64>, Vec<String>) {

    let bar = ProgressBar::new(filtered_var.len() as u64);
    let statistics = get_null_distribution(filtered_var, &matrix, permutation_round, permutation_frequency_threshold);
    let (indices, collected_values): (Vec<_>, Vec<_>) = (0..filtered_var.len()).into_par_iter().map(|i| {
        bar.inc(1);
        let current_variant = filtered_var[i].clone();
        let index = generate_variant_name(&current_variant);
        // let frequency = filtered_var[i].allele_count as f64 / *coverage.get(&filtered_var[i].pos).unwrap() as f64;
        // let frequency = row.sum() / row.len() as f64;
        let frequency = current_variant.allele_count as f64 / *coverage.get(&current_variant.pos).unwrap() as f64;

        if frequency > permutation_frequency_threshold {
            return (Ok(i), None);
        }
        // exclude large indel
        if (current_variant.ref_allele.len() as i32 - current_variant.alt_allele.len() as i32).abs() > 50 {
            return (Ok(i), None);
        }

        let variant = filtered_var[i].clone();
        let pos = variant.pos;
        let ref_allele = variant.ref_allele;
        let alt_allele = variant.alt_allele;

        // For PacBio: exclude SNPs from permutation test (PacBio has high accuracy for SNPs)
        // For ONT: include SNPs in permutation test but with stricter filtering
        // ONT has higher error rates, so we need the permutation test to filter false positives
        if (ref_allele.len() == 1 && alt_allele.len() == 1) && (data_type == "pacbio") {
            return (Ok(i), None);
        }
        // Note: For ONT, SNPs go through permutation test which should help filter false positives
        // If precision is still low, consider adjusting p_value_threshold or frequency_threshold

        let observation = calculate_observation_statistics(&filtered_var,  i, &matrix);
        let p_value = calculate_p_value(&statistics, observation);
        (Err(index.clone()), Some((p_value, index.clone())))

 
    }).unzip(); 


    let mut raw_p_values = Vec::new();
    let mut test_index = Vec::new();

    for item in collected_values.into_iter().flatten() {
        let (p_value, index) = item;
        raw_p_values.push(p_value);
        test_index.push(index.clone());
    }

    // adjust pvalues, create excluded_index list
    let mut excluded_index = Vec::new();
    // println!("{:?}", raw_p_values);
    let qvalues = adjust(&raw_p_values, Procedure::BenjaminiHochberg);

    for (qi, q_value) in qvalues.iter().enumerate(){
        let test_index_value = &test_index[qi];
        // println!("Tested index, {:?}, {:?}", test_index_value, q_value);
        if q_value > &p_value_threshold{
            excluded_index.push(test_index_value);
            // println!("Excluded index, {:?}, {:?}", test_index_value, q_value);
        }
    }
    // println!("Excluded index, {:?}", excluded_index);

    // println!("{:?}", excluded_index);
    bar.finish();


    // filter variants
    let mut index_list = Vec::new();
    let mut f_variant: Vec<Variant> = Vec::new();
    let mut var_list: Vec<String> = Vec::new();
    // get index list and var_list
    for (r, rvariant) in filtered_var.iter().enumerate(){
        let rindex = generate_variant_name(rvariant);
        if !excluded_index.contains(&&rindex.clone()){
            index_list.push(r);
            var_list.push(rindex.clone());
        }
    }

    for v in filtered_var {
        let key = generate_variant_name(&v.clone());
        if excluded_index.contains(&&key.clone()) {
            continue;
        }
        f_variant.push(v.clone());
    }
    // filter matrix
    let filtered_matrix = matrix.select(Axis(0), &index_list);
    // println!("{}, {}, {:?}", f_variant.len(), index_list.len(), filtered_matrix.dim());
    
    (f_variant, filtered_matrix, var_list)
    
}

pub fn resolve_thresholds(
    data_type: &str,
    p_value_threshold: Option<f64>,
    frequency_threshold: Option<f64>,
) -> (f64, f64) {
    let (default_p, default_f) = if data_type == "ont-r9" {
        (0.0001, 0.5)
    } else if data_type == "ont-r10" {
        (0.01, 0.5)
    } else {
        (0.01, 0.5)
    };

    (
        p_value_threshold.unwrap_or(default_p),
        frequency_threshold.unwrap_or(default_f),
    )
}

pub fn start(
    graph_file: &PathBuf,
    fasta_reference: &PathBuf,
    k: usize,
    minimal_ac: usize,
    output_file: &PathBuf,
    sample_id: &str,
    hf_threshold: f32,
    data_type: &str,
    p_value_threshold: f64,
    heteroplasmic_frequency_threshold: f64,
    bam_file: Option<&PathBuf>,
    strand_bias_threshold: f64,
    indel_false_threshold: f64,
    permutation_frequency_threshold: f64,
) {
    if data_type != "pacbio" && data_type != "ont-r9" && data_type != "ont-r10" {
        eprintln!("Error: data type must be pacbio or ont-r9 or ont-r10");
        std::process::exit(1);
    }
    // reference fasta information
    let ref_reader = Reader::from_file(fasta_reference).unwrap();
    let reference_sequence: Vec<Record> = ref_reader.records().map(|r| r.unwrap()).collect();
    let ref_seq = String::from_utf8_lossy(reference_sequence[0].seq()).to_string();
    let ref_header = reference_sequence[0].id().to_string();
    
    // read in graphical genome
    let graph = agg::GraphicalGenome::load_graph(graph_file).unwrap();
    
    let (variants, coverage, read_record) = get_variant(&mut graph.clone(), k, &ref_header, ref_seq.len());
    let collapsed_var = collapse_identical_records(variants, ref_seq.len());
    let filtered_var = filter_vcf_record(&collapsed_var, &coverage, minimal_ac, hf_threshold);
    // modified, exclude filtered data for FPs
    
    // modified, exclude filtered data
    let (matrix, var_record, read_set) = construct_matrix(&read_record, &filtered_var, minimal_ac);
    let matrix_output_raw = output_file.with_extension("raw_matrix.csv");
    let _ = write_matrix_to_csv(&matrix, &var_record, &read_set, matrix_output_raw);

    // use matrix information to filter vcf
    // Use stricter thresholds for ONT data due to higher error rates
    // let (p_value_threshold, frequency_threshold) = if data_type == "ont" {
    //     (0.0001, 0.15)  // Stricter: lower p-value threshold, lower frequency threshold
    // } else {
    //     (0.001, 0.2)    // PacBio: original thresholds
    // };
    let (permu_filtered_var, filtered_matrix, _filtered_name) = permutation_test(&matrix, p_value_threshold, 100, &var_record, &coverage, permutation_frequency_threshold, data_type);

    // strand-bias filter: drop variants whose alt reads are strand-skewed relative to the
    // library composition (only when a BAM is provided so per-read strand is available)
    let (permu_filtered_var, filtered_matrix) = match bam_file {
        Some(bam) => {
            let (strand_map, fwd_frac) = build_strand_map(bam);
            println!("Library forward-strand fraction: {:.3}", fwd_frac);
            filter_strand_bias(
                &permu_filtered_var,
                &filtered_matrix,
                &read_record,
                &strand_map,
                fwd_frac,
                strand_bias_threshold,
                2,
            )
        }
        None => (permu_filtered_var, filtered_matrix),
    };

    // write filtered vcf
    let _ = write_vcf(
        &permu_filtered_var,
        &coverage,
        output_file,
        sample_id,
        &ref_header,
        ref_seq.len(),
        indel_false_threshold
    );
    // write matrix
    let matrix_output = output_file.with_extension("matrix.csv");
    let _ = write_matrix_to_csv(&filtered_matrix, &permu_filtered_var, &read_set, matrix_output);
}