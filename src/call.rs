use crate::agg;
use agg::*;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
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
use statrs::distribution::{Normal, ContinuousCDF};
use regex::Regex;
use adjustp::{adjust, Procedure};
use bio::io::fasta::{Reader, Record};


pub fn find_ref_edge(
    graph: &GraphicalGenome,
    src: &str,
    dst: &str,
    refstrain: &str,
    k: usize,
) -> String {
    // Create HashSet with single refstrain
    let strains = HashSet::from([refstrain.to_string()]);

    // Find all paths between anchors
    let paths = agg::FindAllPathBetweenAnchors::new(&graph, &src, &dst, strains);
    let subpaths = &paths.subpath;

    // Return empty string if no paths found
    if subpaths.is_empty() {
        return String::new();
    }

    // Get sequence from first path
    let (path, _strain) = &subpaths[0];
    let seq = agg::reconstruct_path_seq(&graph, path);

    // Return sequence with src anchor and ref_edge sequence
    seq[..seq.len() - k].to_string()
}

fn mask_ns(seq: &str) -> String {
    seq.chars()
        .map(|c| {
            let upper_c = c.to_ascii_uppercase();
            if !['A', 'G', 'C', 'T'].contains(&upper_c) {
                c.to_ascii_lowercase()
            } else {
                upper_c
            }
        })
        .collect()
}

fn alignment_to_cigar(operations: &[AlignmentOperation]) -> String {
    let mut cigar: Vec<(usize, char)> = Vec::new();

    for op in operations {
        let cigar_op = match op {
            AlignmentOperation::Match => '=',
            AlignmentOperation::Subst => 'X',
            AlignmentOperation::Del => 'D',
            AlignmentOperation::Ins => 'I',
            AlignmentOperation::Xclip(_) => 'S',
            AlignmentOperation::Yclip(_) => 'S',
        };

        if !cigar.is_empty() && cigar.last().unwrap().1 == cigar_op {
            cigar.last_mut().unwrap().0 += 1;
        } else {
            cigar.push((1, cigar_op));
        }
    }

    cigar
        .iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect()
}

fn gap_open_aligner(reference: &str, sequence: &str) -> String {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // Create an aligner with the same scoring parameters
    let mut aligner = Aligner::with_capacity(sequence.len(), reference.len(), -5, -1, &score); // match_score=0, mismatch_score=-6, gap_open=-5, gap_extend=-3

    // Perform the alignment
    let alignment = aligner.global(sequence.as_bytes(), reference.as_bytes());

    // Get the aligned sequences
    let cigar = alignment_to_cigar(&alignment.operations);
    // println!("{:?}", cigar);

    cigar
}

pub fn generate_cigar(
    mut graph: GraphicalGenome,
    ref_strain: &str,
    k: usize,
    maxlength: usize,
    minimal_read_count: usize,
) -> GraphicalGenome {
    let mut edgelist: Vec<_> = graph.edges.keys().cloned().collect();
    edgelist.sort();
    let bar = ProgressBar::new(edgelist.len() as u64);

    let updates: Vec<_> = edgelist
        .par_iter()
        .map(|edge| {
            bar.inc(1);
            let allele_count = graph
                .edges
                .get(edge)
                .and_then(|e| e.get("reads"))
                .map_or(Vec::new(), |reads| {
                    reads.as_array().unwrap_or(&Vec::new()).to_vec()
                })
                .len();

            if allele_count <= minimal_read_count {
                return (edge.clone(), None);
            }

            let src = &graph
                .edges
                .get(edge)
                .expect("Edge not found")
                .get("src")
                .expect("No src in this edge!")
                .as_array()
                .expect("not an array")
                .first()
                .expect("no vallue in src list")
                .as_str()
                .expect("src not a string");
            let src_seq = if *src == "SOURCE" {
                ""
            } else {
                graph
                    .anchor
                    .get(*src)
                    .and_then(|anchor| anchor.get("seq").and_then(|seq| seq.as_str()))
                    .unwrap_or("")
            };

            let dst = &graph
                .edges
                .get(edge)
                .expect("Edge not found")
                .get("dst")
                .expect("No dst in this edge!")
                .as_array()
                .expect("not an array")
                .first()
                .expect("no value in dst list")
                .as_str()
                .expect("dst not a string");
            // println!("src, {}, dst, {}", src, dst);
            let ref_seq = find_ref_edge(&graph, src, dst, ref_strain, k);
            let ref_length = ref_seq.len();
            let alt_sequence = src_seq.to_string()
                + &graph
                    .edges
                    .get(edge)
                    .expect("edge not found")
                    .get("seq")
                    .and_then(|v| v.as_str())
                    .unwrap_or("");
            let alt_seq = mask_ns(&alt_sequence);
            let alt_length = alt_seq.len();

            if ref_length > maxlength || ref_length < k {
                return (edge.clone(), None);
            }

            if alt_seq == ref_seq {
                return (edge.clone(), Some(format!("{}=", ref_length)));
            }

            if alt_length > maxlength {
                return (edge.clone(), None);
            }
            let cigar = gap_open_aligner(&ref_seq, &alt_seq);
            std::thread::sleep(std::time::Duration::from_millis(50));

            (edge.clone(), Some(cigar))
        })
        .collect();

    // Apply the updates sequentially after parallel computation
    for (edge, maybe_variant) in updates {
        if let Some(variant) = maybe_variant {
            if let Some(edge_data) = graph.edges.get_mut(&edge) {
                edge_data
                    .as_object_mut()
                    .expect("Edge data is not an object")
                    .insert("variants".to_string(), serde_json::Value::String(variant));
            }
        }
    }
    bar.finish();
    graph
}

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

pub fn get_variant(
    graph: &mut GraphicalGenome,
    k: usize,
    ref_name: &str,
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

        let ref_seq = find_ref_edge(&graph, src, dst, ref_name, k);
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
        var.extend(variants.clone());

        for (pos, count) in poscounts.iter() {
            *coverage.entry(*pos).or_insert(0) += count;
        }
        let readlist = graph
            .edges
            .get(edge)
            .and_then(|e| e.get("reads"))
            .map_or(Vec::new(), |reads| {
                reads.as_array().unwrap_or(&Vec::new()).to_vec()
            });
        for v in &variants{
            if v.variant_type == "SNP"{
                let key = format!("m.{}{}>{}",
                v.pos + 1,
                v.ref_allele, 
                v.alt_allele);
                read_record.entry(key).or_insert_with(Vec::new).extend(readlist.clone());
            }else {
                let key = format!("m.{}{}>{}",
                v.pos,
                v.ref_allele, 
                v.alt_allele);
                read_record.entry(key).or_insert_with(Vec::new).extend(readlist.clone());
            }

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

fn format_vcf_record(variant: &Variant, coverage: HashMap<usize, usize>) -> String {
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
    
    
    match variant.variant_type.as_str() {
        "SNP" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t.\t{}\t{}\t{}",
            variant.pos + 1,
            variant.ref_allele,
            variant.alt_allele,
            info,
            format,
            sample
        ),
        "INS" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t.\t{}\t{}\t{}",
            variant.pos, variant.ref_allele, variant.alt_allele, info, format, sample
        ),
        "DEL" => format!(
            "chrM\t{}\t.\t{}\t{}\t.\t.\t{}\t{}\t{}",
            variant.pos, variant.ref_allele, variant.alt_allele, info, format, sample
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
    referencelength:usize
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

        writeln!(file, "{}", format_vcf_record(&variant, coverage.clone()))?;
    }

    Ok(())
}


pub fn construct_matrix (read_record:&HashMap<String, Vec<serde_json::Value>>, variants:&[Variant]) -> (Array2<f64>, Vec<String>, Vec<String>) {
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

    let mut var_record: HashSet<String> = HashSet::new();
    for v in variants {
        let key = if v.variant_type == "SNP" {
            format!("m.{}{}>{}",
                v.pos + 1,
                v.ref_allele, 
                v.alt_allele)
        } else {
            format!("m.{}{}>{}",
                v.pos,
                v.ref_allele, 
                v.alt_allele)
        };
        
        var_record.insert(key);
    }
    let mut var_vec: Vec<String> = var_record.into_iter().collect();
    var_vec.sort();

    let var_record_dict: HashMap<String, usize> = var_vec
        .iter()
        .enumerate()
        .map(|(i, var)| (var.clone(), i))
        .collect();


    // Create a 2D matrix filled with zeros
    let mut matrix = Array2::<f64>::zeros((var_vec.len(), read_vec.len()));

    for var in &var_vec {
        let readlist: HashSet<String> = read_record.get(var)
            .map(|reads| {
                reads.iter()
                    .filter_map(|read| read.as_str().map(|s| s.to_string()))
                    .collect()
            })
            .unwrap_or_else(|| HashSet::new());
        
        // Get row index for this variant
        let r_index = *var_record_dict.get(var).unwrap();
        
        // For each read in the readlist
        for read in readlist {
            // Get column index for this read
            if let Some(c_index) = read_set_dict.get(&read) {
                // Increment the matrix cell
                matrix[[r_index, *c_index]] += 1.0;
            }
        }
    }
    (matrix, var_vec, read_vec)


}



fn write_matrix_to_csv<P: AsRef<Path>>(
    matrix: &Array2<f64>,
    var_record: &[String],
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
    for (row_idx, var_name) in var_record.iter().enumerate() {
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
    records: &Vec<String>,
    matrix: &Array2<f64>, 
    permutation_round: usize,
    threshold: f64
) -> Vec<f64> {

    let bar = ProgressBar::new(permutation_round as u64);
    let summary_statistics = (0..permutation_round).into_par_iter().flat_map(|_|  {
        bar.inc(1);
        let mut local_stats = Vec::new();
        let mut rng = thread_rng();

        for (i, index) in records.iter().enumerate() {
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
            
            for (j, other_index) in records.iter().enumerate() {
                if index == other_index {
                    continue;
                }

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
    recordlist: &Vec<String>,
    index: usize,
    matrix: &Array2<f64>, 
) -> f64 {

    let vector = &matrix.slice(s![index, ..]);
    let mut all_coefficients = Vec::new();
    
    for (i, other_index) in recordlist.iter().enumerate() {
        if i == index {
            continue;
        }
        
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
    1.0 - normal.cdf(z_score)
}

fn permutation_test(
    matrix: &Array2<f64>,
    records: Vec<String>,
    p_value_threshold: f64,
    permutation_round: usize,
    filtered_var: &Vec<Variant>,
    frequency_threshold:f64,
    data_type: &str
) -> (Vec<Variant>, Array2<f64>, Vec<String>) {

    let re = Regex::new(r"m\.(\d+)([A-Za-z-]+)>([A-Za-z-]+)").unwrap();
    let bar = ProgressBar::new(records.len() as u64);
    let statistics = get_null_distribution(&records, &matrix, permutation_round, frequency_threshold);
    // Replace par_iter().enumerate() with this pattern
    let (indices, collected_values): (Vec<_>, Vec<_>) = (0..filtered_var.len()).into_par_iter().map(|i| {
        bar.inc(1);
        let index = &records[i];
        let row = matrix.slice(s![i, ..]);
        let frequency = row.sum() / row.len() as f64;
        if frequency > frequency_threshold {
            return (Ok(i), None);
        }
        if let Some(caps) = re.captures(index) {
            let pos = caps.get(1).unwrap().as_str().parse::<usize>().unwrap();
            let ref_allele = caps.get(2).unwrap().as_str();
            let alt_allele = caps.get(3).unwrap().as_str();

            if (ref_allele.len() == alt_allele.len()) && (data_type == "pacbio") {
                return (Ok(i), None);
            }

        }

        let observation = calculate_observation_statistics(&records, i, &matrix);
        let p_value = calculate_p_value(&statistics, observation);
        (Err(index.clone()), Some((p_value, index.clone())))

 
    }).unzip(); 

    let mut raw_p_values = Vec::new();
    let mut test_index = Vec::new();

    for item in collected_values.into_iter().flatten() {
        let (p_value, index) = item;
        raw_p_values.push(p_value);
        test_index.push(index);
    }

    // adjust pvalues, create excluded_index list
    let mut excluded_index = Vec::new();
    let qvalues = adjust(&raw_p_values, Procedure::BenjaminiHochberg);
    for (qi, q_value) in qvalues.iter().enumerate(){
        let test_index_value = &test_index[qi];
        if q_value > &p_value_threshold{
            excluded_index.push(test_index_value);
        }
    }

    println!("{:?}", excluded_index);
    bar.finish();


    // filter variants
    let mut index_list = Vec::new();
    let mut f_variant: Vec<Variant> = Vec::new();
    let mut var_list: Vec<String> = Vec::new();
    // get index list and var_list
    for (r, rindex) in records.iter().enumerate(){
        if !excluded_index.contains(&rindex){
            index_list.push(r);
            var_list.push(rindex.clone());
        }
    }
    // for idx in &index_list {
    //     var_list.push(records[*idx].clone())
    // }
    for v in filtered_var {
        let key = if v.variant_type == "SNP" {
            format!("m.{}{}>{}",
                v.pos + 1,
                v.ref_allele, 
                v.alt_allele)
        } else {
            format!("m.{}{}>{}",
                v.pos,
                v.ref_allele, 
                v.alt_allele)
        };
        if excluded_index.contains(&&key) {
            continue;
        }
        f_variant.push(v.clone());
    }
    // filter matrix
    let filtered_matrix = matrix.select(Axis(0), &index_list);
    // println!("{}, {}, {:?}", f_variant.len(), index_list.len(), filtered_matrix.dim());
    
    (f_variant, filtered_matrix, var_list)
    
}

pub fn start(
    graph_file: &PathBuf,
    fasta_reference: &PathBuf,
    k: usize,
    maxlength: usize,
    minimal_ac: usize,
    output_file: &PathBuf,
    sample_id: &str,
    hf_threshold: f32,
    data_type: &str,
) {
    if data_type != "pacbio" && data_type != "ont" {
        eprintln!("Error: data type must be pacbio or ont");
        std::process::exit(1);
    }
    // reference fasta information
    let ref_reader = Reader::from_file(fasta_reference).unwrap();
    let reference_sequence: Vec<Record> = ref_reader.records().map(|r| r.unwrap()).collect();
    let ref_seq = String::from_utf8_lossy(reference_sequence[0].seq()).to_string();
    let ref_header = reference_sequence[0].id().to_string();
    
    // read in graphical genome
    let graph = agg::GraphicalGenome::load_graph(graph_file).unwrap();
    // generate cigar
    let mut graph_with_cigar = generate_cigar(graph, &ref_header, k, maxlength, 2);
    let (variants, coverage, read_record) = get_variant(&mut graph_with_cigar, k, &ref_header);
    let collapsed_var = collapse_identical_records(variants, ref_seq.len());
    let filtered_var = filter_vcf_record(&collapsed_var, &coverage, minimal_ac, hf_threshold);
    // modified, exclude filtered data for FPs

    let graph_output = output_file.with_extension("annotated.gfa");
    let _ = write_graph_from_graph(graph_output.to_str().unwrap(), &graph_with_cigar);
    
    // modified, exclude filtered data
    let (matrix, var_record, read_set) = construct_matrix(&read_record, &filtered_var);

    // use matrix information to filter vcf
    let (permu_filtered_var, filtered_matrix, filtered_name) = permutation_test(&matrix, var_record, 0.001, 100, &filtered_var, 0.2, data_type);

    
    // write filtered vcf
    let _ = write_vcf(
        &permu_filtered_var,
        &coverage,
        output_file,
        sample_id,
        &ref_header,
        ref_seq.len()
    );
    // write matrix
    let matrix_output = output_file.with_extension("matrix.csv");
    let _ = write_matrix_to_csv(&filtered_matrix, &filtered_name, &read_set, matrix_output);
}

