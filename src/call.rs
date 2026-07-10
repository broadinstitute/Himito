use crate::agg;
use agg::*;
use rayon::iter::WhileSome;
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
use ndarray::{Array2, Axis, s};
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
    pub filter: Option<String>,
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
                        filter:None,
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
                    filter: None,
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
                    filter: None,
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
                filter: variant.filter,
            });
        }
    }
    circular_variants
}

fn edge_reads(graph: &GraphicalGenome, edge: &str) -> HashSet<String> {
    graph
        .edges
        .get(edge)
        .and_then(|e| e.get("reads"))
        .map(|reads| {
            reads
                .as_array()
                .unwrap_or(&Vec::new())
                .iter()
                .filter_map(|r| r.as_str().map(String::from))
                .collect()
        })
        .unwrap_or_default()
}

fn get_graph_intervals(graph:&GraphicalGenome, length: i64) -> HashMap<&String, (i64, i64)>{
    let mut graph_intervals_dict = HashMap::new();
    for edge in graph.edges.keys() {
        // Edges with no computed variant record were never CIGAR-processed
        // (get_variant skips them via this same empty-cigar check) - a read
        // reaching a bubble only through such an edge was never genotyped
        // there, and must not be credited as "covering" it, or it later
        // defaults to a false ref (0) call instead of the correct NaN.
        let has_variant = graph.edges[edge]
            .get("variants")
            .and_then(|v| v.as_str())
            .map_or(false, |s| !s.is_empty());
        if !has_variant {
            continue;
        }
        let src = graph.edges[edge].get("src").unwrap().as_array().unwrap()[0].as_str().unwrap();
        let dst = graph.edges[edge].get("dst").unwrap().as_array().unwrap()[0].as_str().unwrap();
        let startpos = graph.anchor
            .get(src)
            .and_then(|v| v.get("pos"))
            .and_then(|v| v.as_i64())
            .unwrap_or(0) as i64;  
        let endpos  = graph.anchor
            .get(dst)
            .and_then(|v| v.get("pos"))
            .and_then(|v| v.as_i64())
            .unwrap_or(length) as i64;  
        if endpos > startpos {
            graph_intervals_dict.insert(edge, (startpos, endpos));
        }
        
    }
    graph_intervals_dict
}


// /// All reads traversing any parallel edge between the same anchor pair (ref + alt paths).
fn bubble_cover_reads(graph_intervals_dict: &HashMap<&String, (i64, i64)>, pos: usize, graph: &GraphicalGenome) -> HashSet<String> {
    // let empty = Vec::new();
    let mut cover = HashSet::new();
    for (edge, (start, end)) in graph_intervals_dict.iter() {
        if pos >= *start as usize && pos < *end as usize {
            cover.extend(edge_reads(graph, edge));
        }
    }
    cover
}

/// Populates `read_record`/`cover_record` for one edge's variants.
///
/// `graph_intervals_dict` holds raw (unwrapped) anchor coordinates, including
/// the duplicated wrap-around span past `ref_length` used to circularize the
/// graph. `variants` carries that same raw position, while `variants_circular`
/// carries the wrapped position used to key/name the variant. The bubble
/// lookup must use the raw position - using the wrapped one would look up
/// the wrong bubble (or none) for any variant on the wrap-around edge.
fn record_variant_reads(
    variants: &[Variant],
    variants_circular: &[Variant],
    graph_intervals_dict: &HashMap<&String, (i64, i64)>,
    graph: &GraphicalGenome,
    readlist: &[serde_json::Value],
    read_record: &mut HashMap<String, Vec<serde_json::Value>>,
    cover_record: &mut HashMap<String, HashSet<String>>,
) {
    for (raw_v, v) in variants.iter().zip(variants_circular.iter()) {
        let key = generate_variant_name(&v.clone());
        let cover_reads = bubble_cover_reads(&graph_intervals_dict, raw_v.pos, graph);
        read_record
            .entry(key.clone())
            .or_insert_with(Vec::new)
            .extend(readlist.to_vec());
        cover_record.entry(key).or_insert_with(HashSet::new).extend(cover_reads);
    }
}

pub fn get_variant(
    graph: &mut GraphicalGenome,
    k: usize,
    ref_name: &str,
    ref_length: usize,
) -> (
    Vec<Variant>,
    HashMap<usize, usize>,
    HashMap<String, Vec<serde_json::Value>>,
    HashMap<String, HashSet<String>>,
) {
    let mut coverage = HashMap::new();
    let mut read_record = HashMap::new();
    let mut cover_record: HashMap<String, HashSet<String>> = HashMap::new();
    let mut var = Vec::new();
    let mut edgelist: Vec<_> = graph.edges.keys().collect();
    let graph_intervals_dict = get_graph_intervals(graph, ref_length as i64);
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
        
        record_variant_reads(
            &variants,
            &variants_circular,
            &graph_intervals_dict,
            graph,
            &readlist,
            &mut read_record,
            &mut cover_record,
        );
    }
    (var, coverage, read_record, cover_record)
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
                filter: None,
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
    // Combine any upstream flag (e.g. Strand_bias) with the indel-artifact check here, rather
    // than letting one filter suppress the other, so a call can carry multiple FILTER reasons.
    let mut filters: Vec<String> = Vec::new();
    if let Some(existing) = variant.filter.as_deref() {
        if !existing.is_empty() {
            filters.push(existing.to_string());
        }
    }
    if is_small_indel && allele_frequency < indel_false_threshold as f32 {
        filters.push("Potential_Artifact".to_string());
    }
    let filter = if filters.is_empty() {
        "PASS".to_string()
    } else {
        filters.join(";")
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
        "##FILTER=<ID=Strand_bias,Description=\"Alt-supporting reads are significantly skewed toward one strand relative to library composition\">"
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

pub fn construct_matrix(
    read_record: &HashMap<String, Vec<serde_json::Value>>,
    cover_record: &HashMap<String, HashSet<String>>,
    variants: &[Variant],
    minimal_ac: usize,
) -> (Array2<f64>, Vec<Variant>, Vec<String>) {
    let mut read_set: HashSet<String> = HashSet::new();
    for reads in cover_record.values() {
        read_set.extend(reads.iter().cloned());
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

    // NaN = missing coverage; 0 = ref; 1 = alt
    let mut matrix = Array2::<f64>::from_elem((var_vec.len(), read_vec.len()), f64::NAN);

    for var in &var_vec {
        let name = generate_variant_name(var);
        let r_index = *var_record_dict.get(&name).unwrap();
        let alt_reads: HashSet<String> = read_record
            .get(&name)
            .map(|reads| {
                reads
                    .iter()
                    .filter_map(|read| read.as_str().map(|s| s.to_string()))
                    .collect()
            })
            .unwrap_or_default();
        let cover_reads = cover_record.get(&name).cloned().unwrap_or_default();

        for read in cover_reads {
            if let Some(c_index) = read_set_dict.get(&read) {
                matrix[[r_index, *c_index]] = if alt_reads.contains(&read) {
                    1.0
                } else {
                    0.0
                };
            }
        }
    }

    // Keep reads with more than minimal_ac alt calls across variants
    let mut col_index = Vec::new();
    for i in 0..matrix.ncols() {
        let alt_count = matrix
            .column(i)
            .iter()
            .filter(|&&x| x == 1.0)
            .count();
        if alt_count > minimal_ac {
            col_index.push(i);
        }
    }
    let filtered_matrix = matrix.select(Axis(1), &col_index);
    let filtered_read_vec: Vec<String> = col_index
        .iter()
        .map(|&i| read_vec[i].clone())
        .collect();
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
        
        // Add the values from the matrix: empty string for missing (NaN)
        for col_idx in 0..matrix.ncols() {
            let v = matrix[[row_idx, col_idx]];
            if v.is_nan() {
                row.push(String::new());
            } else {
                row.push(v.to_string());
            }
        }
        
        writer.write_record(&row)?;
    }
    
    // Flush and finish
    writer.flush()?;
    Ok(())
}

fn alt_frequency(vector: ndarray::ArrayView1<f64>) -> f64 {
    let covered: Vec<f64> = vector.iter().filter(|x| !x.is_nan()).copied().collect();
    if covered.is_empty() {
        return 0.0;
    }
    covered.iter().sum::<f64>() / covered.len() as f64
}

fn shuffle_genotypes(data: &[f64], rng: &mut impl rand::Rng) -> Vec<f64> {
    let mut result = data.to_vec();
    let covered_indices: Vec<usize> = result
        .iter()
        .enumerate()
        .filter(|(_, x)| !x.is_nan())
        .map(|(i, _)| i)
        .collect();
    let mut values: Vec<f64> = covered_indices.iter().map(|&i| result[i]).collect();
    values.shuffle(rng);
    for (idx, val) in covered_indices.iter().zip(values.iter()) {
        result[*idx] = *val;
    }
    result
}

fn genotype_jaccard_similarity(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "Vectors must have the same length");

    let mut intersection_count = 0;
    let mut union_count = 0;

    for (x, y) in a.iter().zip(b.iter()) {
        if x.is_nan() || y.is_nan() {
            continue;
        }
        let ax = *x > 0.5;
        let by = *y > 0.5;
        if ax && by {
            intersection_count += 1;
        }
        if ax || by {
            union_count += 1;
        }
    }

    if union_count == 0 {
        return 0.0;
    }

    intersection_count as f64 / union_count as f64
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
            let frequency = alt_frequency(vector);
            if frequency > threshold {
                continue;
            }

            // Shuffle genotypes only among covered reads; keep NaN positions fixed
            let vector_data: Vec<f64> = vector.iter().copied().collect();
            let shuffled = shuffle_genotypes(&vector_data, &mut rng);

            let mut all_coefficients = Vec::new();

            for (j, other_variant) in filtered_var.iter().enumerate() {
                let other_index = generate_variant_name(other_variant);
                if index == other_index {
                    continue;
                }

                let other_vector = matrix.slice(s![j, ..]);
                let other_data: Vec<f64> = other_vector.iter().copied().collect();

                let coor = genotype_jaccard_similarity(&shuffled, &other_data);
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

    let vector = matrix.slice(s![index, ..]);
    let vector_data: Vec<f64> = vector.iter().copied().collect();
    let mut all_coefficients = Vec::new();

    for (i, variant) in filtered_var.iter().enumerate() {
        if i == index {
            continue;
        }

        let other_vector = matrix.slice(s![i, ..]);
        let other_data: Vec<f64> = other_vector.iter().copied().collect();

        let coor = genotype_jaccard_similarity(&vector_data, &other_data);
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

/// Flag (but never drop) variants whose alt-supporting reads are significantly strand-skewed
/// relative to the library's forward/reverse composition; flagged variants get
/// `filter = Some("Strand_bias")` so the call survives in the VCF with FILTER set accordingly.
/// Variants with fewer than `min_reads`
/// strand-resolved supporting reads are kept (too little evidence to judge). Near-homoplasmic
/// variants (heteroplasmic frequency >= `homoplasmic_frequency_threshold`) are also kept
/// without testing: a real single-strand-only sequencing artifact can only ever produce a
/// low-to-moderate apparent frequency (it's absent from the unaffected strand's reads), so an
/// allele present in nearly every read at the locus cannot be explained by strand-specific
/// error and any measured strand skew there reflects fragmented read attribution across
/// near-duplicate graph edges, not a real signal. Returns the kept variants together with the
/// matrix subset to the corresponding rows.
fn filter_strand_bias(
    variants: &[Variant],
    matrix: &Array2<f64>,
    read_record: &HashMap<String, Vec<serde_json::Value>>,
    strand_map: &HashMap<String, bool>,
    coverage: &HashMap<usize, usize>,
    expected_fwd_frac: f64,
    p_threshold: f64,
    min_reads: usize,
    homoplasmic_frequency_threshold: f64,
) -> (Vec<Variant>, Array2<f64>) {
    let mut keep_idx = Vec::new();
    let mut kept = Vec::new();

    for (i, variant) in variants.iter().enumerate() {
        let mut variant = variant.clone();
        let name = generate_variant_name(&variant);

        let depth = coverage.get(&variant.pos).copied().unwrap_or(0);
        let hf = if depth == 0 { 0.0 } else { variant.allele_count as f64 / depth as f64 };
        if hf >= homoplasmic_frequency_threshold {
            keep_idx.push(i);
            kept.push(variant);
            continue;
        }

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
            // not enough strand-resolved reads to make a call: keep, untested
            keep_idx.push(i);
            kept.push(variant);
            continue;
        }

        let p_value = strand_bias_pvalue(fwd, rev, expected_fwd_frac);
        if p_value < p_threshold {
            // flag but keep: downstream FILTER combines this with other filters (e.g.
            // Potential_Artifact) instead of dropping the call
            variant.filter = Some("Strand_bias".to_string());
            println!(
                "Strand-biased (flagged, kept), {:?}, fwd={}, rev={}, p={:.3e}",
                name, fwd, rev, p_value
            );
        }
        keep_idx.push(i);
        kept.push(variant);
    }

    let filtered_matrix = matrix.select(Axis(0), &keep_idx);
    (kept, filtered_matrix)
}

fn permutation_test(
    matrix: &Array2<f64>,
    p_value_threshold: f64,
    heteroplasmic_error_threshold:f64,
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

        if frequency > heteroplasmic_error_threshold {
            return (Ok(i), None);
        }
        // exclude large indel
        if (current_variant.ref_allele.len() as i32 - current_variant.alt_allele.len() as i32).abs() > 50 {
            return (Ok(i), None);
        }

        let variant = filtered_var[i].clone();
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
    permutation_frequency_threshold: Option<f64>,
) -> (f64, f64, f64) {
    let (default_p, default_f, default_perm) = if data_type == "ont-r9" {
        (0.0001, 0.5, 0.7)
    } else if data_type == "ont-r10" {
        (0.01, 0.2, 0.7)
    } else {
        (0.01, 0.2, 0.8)
    };

    (
        p_value_threshold.unwrap_or(default_p),
        frequency_threshold.unwrap_or(default_f),
        permutation_frequency_threshold.unwrap_or(default_perm),
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
    
    let (variants, coverage, read_record, cover_record) =
        get_variant(&mut graph.clone(), k, &ref_header, ref_seq.len());
    let collapsed_var = collapse_identical_records(variants, ref_seq.len());
    let filtered_var = filter_vcf_record(&collapsed_var, &coverage, minimal_ac, hf_threshold);
    // modified, exclude filtered data for FPs

    // modified, exclude filtered data
    let (matrix, var_record, read_set) =
        construct_matrix(&read_record, &cover_record, &filtered_var, minimal_ac);
    let matrix_output_raw = output_file.with_extension("raw_matrix.csv");
    let _ = write_matrix_to_csv(&matrix, &var_record, &read_set, matrix_output_raw);

    // // use matrix information to filter vcf
    // Use stricter thresholds for ONT data due to higher error rates
    // let (p_value_threshold, frequency_threshold) = if data_type == "ont" {
    //     (0.0001, 0.15)  // Stricter: lower p-value threshold, lower frequency threshold
    // } else {
    //     (0.001, 0.2)    // PacBio: original thresholds
    // };
    let (permu_filtered_var, filtered_matrix, _filtered_name) = permutation_test(&matrix, p_value_threshold, heteroplasmic_frequency_threshold, 100, &var_record, &coverage, permutation_frequency_threshold, data_type);

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
                &coverage,
                fwd_frac,
                strand_bias_threshold,
                2,
                permutation_frequency_threshold,
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

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// Two bubbles: a normal one spanning raw positions [0, 5), and a
    /// wrap-around bubble spanning raw positions [8, 12) - i.e. it crosses
    /// the circular origin for a 10bp reference.
    fn wraparound_test_graph() -> GraphicalGenome {
        let mut anchor = HashMap::new();
        anchor.insert("A1".to_string(), json!({"pos": 0}));
        anchor.insert("A2".to_string(), json!({"pos": 5}));
        anchor.insert("A3".to_string(), json!({"pos": 8}));
        anchor.insert("A4".to_string(), json!({"pos": 12}));

        let mut edges = HashMap::new();
        edges.insert(
            "E1".to_string(),
            json!({"src": ["A1"], "dst": ["A2"], "reads": ["r1", "r2"], "variants": "5="}),
        );
        edges.insert(
            "E2".to_string(),
            json!({"src": ["A1"], "dst": ["A2"], "reads": ["r3"], "variants": "2=1X2="}),
        );
        edges.insert(
            "E3".to_string(),
            json!({"src": ["A3"], "dst": ["A4"], "reads": ["r4"], "variants": "4="}),
        );
        edges.insert(
            "E4".to_string(),
            json!({"src": ["A3"], "dst": ["A4"], "reads": ["r5"], "variants": "1=1X2="}),
        );

        GraphicalGenome {
            anchor,
            edges,
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        }
    }

    #[test]
    fn record_variant_reads_uses_raw_position_for_wraparound_bubble() {
        let graph = wraparound_test_graph();
        let ref_length = 10usize;
        let graph_intervals_dict = get_graph_intervals(&graph, ref_length as i64);

        // Raw position 11 lies on the wrap-around bubble [8, 12) and
        // circularizes down to wrapped position 1.
        let raw_variant = Variant {
            pos: 11,
            ref_allele: "A".to_string(),
            alt_allele: "G".to_string(),
            variant_type: "SNP".to_string(),
            allele_count: 1,
            filter: None,
        };
        let variants = vec![raw_variant];
        let variants_circular = circuliarize_variants(variants.clone(), ref_length);
        assert_eq!(variants_circular[0].pos, 1, "sanity check: position should wrap");

        let readlist = vec![json!("r4"), json!("r5")];
        let mut read_record = HashMap::new();
        let mut cover_record = HashMap::new();

        record_variant_reads(
            &variants,
            &variants_circular,
            &graph_intervals_dict,
            &graph,
            &readlist,
            &mut read_record,
            &mut cover_record,
        );

        let key = generate_variant_name(&variants_circular[0]);
        let mut cover: Vec<String> = cover_record.get(&key).cloned().unwrap_or_default().into_iter().collect();
        cover.sort();

        // The wrap-around bubble's true reads are r4/r5. Using the wrapped
        // position (1) instead of the raw one (11) would incorrectly match
        // the unrelated bubble at [0, 5), pulling in r1/r2/r3 instead.
        assert_eq!(cover, vec!["r4".to_string(), "r5".to_string()]);
    }

    /// A bubble at raw positions [0, 10) carries the real variant (edge
    /// `E_real`, cigar `5=1X4=`), with reads r1/r2. A read that couldn't be
    /// assembled through the fine-grained anchor chain instead crosses the
    /// same span via a private "skip-over" edge with no computed CIGAR
    /// (`E_private`, read r3) - mirroring what Himito's own graph does for
    /// reads it can't cleanly anchor (see e.g. E00001050 in a real run).
    fn skip_edge_test_graph() -> GraphicalGenome {
        let mut anchor = HashMap::new();
        anchor.insert("A1".to_string(), json!({"pos": 0}));
        anchor.insert("A2".to_string(), json!({"pos": 10}));

        let mut edges = HashMap::new();
        edges.insert(
            "E_real".to_string(),
            json!({"src": ["A1"], "dst": ["A2"], "reads": ["r1", "r2"], "variants": "5=1X4="}),
        );
        edges.insert(
            "E_private".to_string(),
            json!({"src": ["A1"], "dst": ["A2"], "reads": ["r3"]}),
        );

        GraphicalGenome {
            anchor,
            edges,
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        }
    }

    #[test]
    fn bubble_cover_reads_excludes_edges_with_no_computed_variant() {
        let graph = skip_edge_test_graph();
        let graph_intervals_dict = get_graph_intervals(&graph, 10);

        let cover = bubble_cover_reads(&graph_intervals_dict, 5, &graph);
        let mut cover: Vec<String> = cover.into_iter().collect();
        cover.sort();

        // r3 only reaches this span via E_private, which was never CIGAR-
        // processed (get_variant would skip it via the empty-cigar check).
        // It must not be credited as "covering" the variant - crediting it
        // as covered-but-not-alt would wrongly call it ref (0) in the
        // matrix instead of the correct NaN (unassembled/unknown).
        assert_eq!(cover, vec!["r1".to_string(), "r2".to_string()]);
    }
}