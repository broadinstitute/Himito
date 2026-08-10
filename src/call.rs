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
use statrs::distribution::{Binomial, DiscreteCDF};
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
        // SOURCE/SINK have no anchor coordinate, and get_variant skips those edges
        // outright. Defaulting the missing endpoint to 0/length instead would hand
        // them a genome-scale interval, crediting their reads as "covering" every
        // variant in the graph.
        let anchor_pos = |name: &str| {
            graph
                .anchor
                .get(name)
                .and_then(|v| v.get("pos"))
                .and_then(|v| v.as_i64())
        };
        let (Some(startpos), Some(endpos)) = (anchor_pos(src), anchor_pos(dst)) else {
            continue;
        };
        // The graph is circular: the final anchor links straight back to the first,
        // so the closing bubble has `endpos < startpos` (e.g. A016544 -> A000022 on
        // rCRS). Unwrap its end past the origin rather than dropping the edge --
        // dropping it left every variant in that bubble with no covering reads, so
        // origin-spanning variants could be emitted to the VCF yet never genotyped
        // in the matrix. Variant positions on this edge are raw/unwrapped for the
        // same reason (circuliarize_variants folds them back afterwards), so the
        // unwrapped interval is exactly what bubble_cover_reads is queried with.
        let end_unwrapped = if endpos > startpos {
            endpos
        } else if endpos < startpos {
            endpos + length
        } else {
            // Degenerate self-loop on one anchor: zero-width, nothing to cover.
            continue;
        };
        graph_intervals_dict.insert(edge, (startpos, end_unwrapped));
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
/// `graph_intervals_dict` holds raw (unwrapped) anchor coordinates: the closing
/// bubble that crosses the origin is stored as `[src.pos, dst.pos + ref_length)`
/// so it stays ordered and half-open. `variants` carries that same raw position,
/// while `variants_circular` carries the wrapped position used to key/name the
/// variant. The bubble lookup must use the raw position - using the wrapped one
/// would look up the wrong bubble (or none) for any variant past the origin.
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
        let alt_allele = variant.alt_allele.as_str();
        if alt_allele.contains("N") {
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

    // Every covered read is kept. There used to be a column filter here dropping
    // reads with `alt_count <= minimal_ac` ALT calls, but `minimal_ac` is a
    // per-variant allele-count gate ("a variant needs >= N supporting reads"),
    // not a per-read one, and reusing it that way removed ~91% of reads on a
    // real ONT run. Worse, it selected on noise: `alt_count` is summed over ALL
    // raw variants, which on ONT are mostly indel artifacts, so retention was
    // flat (~9%) regardless of how many true mutations a read carried and the
    // reference-like reads that anchor the lineage root were dropped wholesale.
    // Downstream SCITE models missing data explicitly via its FN rate, so sparse
    // reads cost it far less than a biased subsample does.
    (matrix, var_vec, read_vec)
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

/// Sum of Jaccard coefficients between one genotype vector and every OTHER variant's
/// row. Rows are compared by index rather than by `generate_variant_name`, so two
/// variants that happen to render to the same name cannot silently exclude each other.
fn sum_jaccard_against_others(vector_data: &[f64], self_index: usize, rows: &[Vec<f64>]) -> f64 {
    rows.iter()
        .enumerate()
        .filter(|(j, _)| *j != self_index)
        .map(|(_, other)| genotype_jaccard_similarity(vector_data, other))
        .sum()
}

/// Generate a PER-VARIANT null distribution through permutation testing, keyed by
/// `generate_variant_name`.
///
/// The summary statistic is a SUM of Jaccard coefficients over every other variant, so
/// it scales with how many reads carry the variant. Pooling all variants' draws into one
/// null (as this did originally) therefore turned the test into a depth filter: a
/// variant supported by 5 reads was z-scored against a null dominated by variants
/// supported by hundreds, landing far below the pooled mean and yielding p ~ 1 no matter
/// how tightly it co-segregated. Conversely a high-depth artifact cleared the pooled bar
/// on read count alone. Since `shuffle_genotypes` permutes a variant's own values among
/// its own covered reads, it already preserves that variant's alt count -- so each
/// variant's draws are the correct null FOR THAT VARIANT and simply need to be kept
/// apart rather than concatenated.
///
/// Variants whose alt frequency exceeds `threshold` are near-homoplasmic: they are the
/// co-occurrence backbone the candidates are scored against, not candidates themselves,
/// and are omitted from the returned map (as before).
fn get_null_distribution(
    filtered_var: &Vec<Variant>,
    matrix: &Array2<f64>,
    permutation_round: usize,
    threshold: f64
) -> HashMap<String, Vec<f64>> {

    let bar = ProgressBar::new(filtered_var.len() as u64);
    // Materialize the rows ONCE. The inner loop compares against every other row on
    // every round, and re-collecting each `ArrayView1` per comparison made the
    // allocation count scale with rounds x variants^2.
    let rows: Vec<Vec<f64>> = (0..filtered_var.len())
        .map(|j| matrix.slice(s![j, ..]).iter().copied().collect())
        .collect();

    let nulls: Vec<(String, Vec<f64>)> = filtered_var
        .par_iter()
        .enumerate()
        .filter_map(|(i, variant)| {
            bar.inc(1);

            // Skip vectors with frequency > threshold
            if alt_frequency(matrix.slice(s![i, ..])) > threshold {
                return None;
            }

            let mut rng = thread_rng();
            let draws: Vec<f64> = (0..permutation_round)
                .map(|_| {
                    // Shuffle genotypes only among covered reads; keep NaN positions fixed
                    let shuffled = shuffle_genotypes(&rows[i], &mut rng);
                    sum_jaccard_against_others(&shuffled, i, &rows)
                })
                .collect();

            Some((generate_variant_name(variant), draws))
        })
        .collect();
    bar.finish();
    nulls.into_iter().collect()
}

/// Calculate statistics for observed data
fn calculate_observation_statistics(
    filtered_var: &Vec<Variant>,
    index: usize,
    matrix: &Array2<f64>, 
) -> f64 {

    let rows: Vec<Vec<f64>> = (0..filtered_var.len())
        .map(|j| matrix.slice(s![j, ..]).iter().copied().collect())
        .collect();
    sum_jaccard_against_others(&rows[index], index, &rows)
}

/// One-sided empirical (rank) p-value: the fraction of `null` draws that reach or
/// exceed `observation`, with the standard add-one correction in both numerator and
/// denominator.
///
/// This replaces a z-score against a normal CDF. The statistic is a sum of Jaccard
/// coefficients bounded below by zero and strongly right-skewed -- for a scattered
/// artifact most draws are exactly 0.0 -- so a normal approximation is a poor fit and
/// its degenerate cases needed a NaN guard (sigma == 0 whenever every draw ties, which
/// is the COMMON case for a low-VAF variant, not a rare one). A rank p-value is exact
/// under the permutation null and needs no distributional assumption.
///
/// The add-one makes the smallest attainable p-value `1 / (R + 1)`, which is the honest
/// floor: R permutations cannot supply evidence beyond that. It also means R must be
/// large enough for the floor to sit below the significance threshold with room for the
/// Benjamini-Hochberg correction -- at R = 100 the floor is 0.0099, which leaves no
/// resolution at all against a 0.01 threshold.
///
/// An empty null carries no evidence, so it yields 1.0. Callers must treat that as
/// "untestable" rather than feeding it to BH, where a p-value of 1.0 reads as a
/// confirmed artifact.
pub(crate) fn empirical_p_value(null: &[f64], observation: f64) -> f64 {
    if null.is_empty() {
        return 1.0;
    }
    let at_or_above = null.iter().filter(|&&x| x >= observation).count();
    (at_or_above + 1) as f64 / (null.len() + 1) as f64
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
pub(crate) fn strand_bias_pvalue(fwd: usize, rev: usize, expected_fwd_frac: f64) -> f64 {
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

/// What `filter_strand_bias` does with a variant it finds strand-skewed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum StrandBiasAction {
    /// Set `FILTER=Strand_bias` but keep the call in the VCF (historical behaviour).
    Flag,
    /// Additionally REMOVE strand-skewed low-VAF SNPs outright.
    Drop,
}

/// Resolve the strand-bias action for a data type, honouring an explicit CLI choice.
///
/// Only `ont-denoised` defaults to dropping. Its SNPs are the ones with no other
/// discriminative filter, and denoising has already erased the strictly single-strand
/// artifacts (`fit_site`'s per-strand candidacy gate maps them to reference), so a call
/// that is STILL strand-skewed after denoising has survived a filter designed to catch
/// it. PacBio and raw ONT keep flag-only behaviour.
pub fn resolve_strand_bias_action(
    data_type: &str,
    explicit: Option<StrandBiasAction>,
) -> StrandBiasAction {
    explicit.unwrap_or(if data_type == "ont-denoised" {
        StrandBiasAction::Drop
    } else {
        StrandBiasAction::Flag
    })
}

/// Flag (and, under `StrandBiasAction::Drop`, remove) variants whose alt-supporting reads are significantly strand-skewed
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
    action: StrandBiasAction,
    drop_max_hf: f64,
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
            // The drop is deliberately narrow: SNPs only, and only below `drop_max_hf`.
            // Indels already have `Potential_Artifact` to catch them, and above the
            // ceiling a skew is likelier to be fragmented read attribution across
            // near-duplicate graph edges than a single-strand sequencing artifact.
            if action == StrandBiasAction::Drop
                && variant.variant_type == "SNP"
                && hf < drop_max_hf
            {
                println!(
                    "Strand-biased (dropped), {:?}, fwd={}, rev={}, hf={:.4}, p={:.3e}",
                    name, fwd, rev, hf, p_value
                );
                continue;
            }
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
    let nulls = get_null_distribution(filtered_var, &matrix, permutation_round, permutation_frequency_threshold);
    let collected_values: Vec<Option<(f64, String)>> = (0..filtered_var.len()).into_par_iter().map(|i| {
        bar.inc(1);
        let current_variant = &filtered_var[i];
        let index = generate_variant_name(current_variant);
        // A variant whose position carries no coverage entry cannot be scored; the
        // original `.unwrap()` here panicked on that instead of leaving it untested.
        let depth = coverage.get(&current_variant.pos).copied().unwrap_or(0);
        let frequency = if depth == 0 {
            0.0
        } else {
            current_variant.allele_count as f64 / depth as f64
        };

        if frequency > heteroplasmic_error_threshold {
            return None;
        }
        // exclude large indel
        if (current_variant.ref_allele.len() as i32 - current_variant.alt_allele.len() as i32).abs() > 50 {
            return None;
        }

        let is_snp = current_variant.ref_allele.len() == 1 && current_variant.alt_allele.len() == 1;

        // PacBio: exclude SNPs from the permutation test -- its SNP accuracy is high
        // enough that the co-occurrence test costs more true low-VAF heteroplasmies
        // than it saves false positives.
        //
        // ONT, INCLUDING `ont-denoised`: test SNPs. Denoising is a per-base quality
        // model, so it strips the RANDOM error it can see and leaves the SYSTEMATIC
        // error -- which carries high base qualities -- untouched. Exempting denoised
        // SNPs here left them with no discriminative filter anywhere in the pipeline:
        // `filter_strand_bias` only flagged, and `Potential_Artifact` is indel-only.
        // That gap is exactly where the residual low-VAF SNP false positives came from.
        if is_snp && data_type == "pacbio" {
            return None;
        }

        // No null draws for this variant (it sits above `permutation_frequency_threshold`
        // by MATRIX alt frequency even though allele_count/coverage admitted it here --
        // the two denominators can disagree). Absent evidence is not evidence of an
        // artifact: leave it untested rather than handing BH a p-value of 1.0.
        let null = match nulls.get(&index) {
            Some(n) if !n.is_empty() => n,
            _ => return None,
        };

        let observation = calculate_observation_statistics(&filtered_var, i, &matrix);
        Some((empirical_p_value(null, observation), index))
    }).collect();


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
        (0.01, 0.2, 0.7)
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
    permutation_rounds: usize,
    strand_bias_action: Option<StrandBiasAction>,
    strand_bias_drop_max_hf: f64,
) {
    if data_type != "pacbio" && data_type != "ont-r9" && data_type != "ont-r10" && data_type != "ont-denoised" {
        eprintln!("Error: data type must be pacbio or ont-r9 or ont-r10 or ont-denoised");
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
        construct_matrix(&read_record, &cover_record, &filtered_var);
    let matrix_output_raw = output_file.with_extension("raw_matrix.csv");
    let _ = write_matrix_to_csv(&matrix, &var_record, &read_set, matrix_output_raw);

    // // use matrix information to filter vcf
    // Use stricter thresholds for ONT data due to higher error rates
    // let (p_value_threshold, frequency_threshold) = if data_type == "ont" {
    //     (0.0001, 0.15)  // Stricter: lower p-value threshold, lower frequency threshold
    // } else {
    //     (0.001, 0.2)    // PacBio: original thresholds
    // };
    // Strand-bias filtering runs BEFORE the permutation test, not after.
    //
    // The permutation test scores a variant by how strongly it co-segregates with the
    // OTHER variants in the matrix. Error-prone ONT reads carry many errors at once, so
    // systematic artifacts co-occur with each other on the same bad reads and produce a
    // high Jaccard sum that is indistinguishable from a real lineage. Leaving known
    // strand-skewed artifacts in the matrix therefore fed them into both the null and
    // every observation, teaching the test that artifacts are real evidence. Removing
    // them first means the co-occurrence backbone is built only from calls that already
    // survived an independent, read-level check.
    let action = resolve_strand_bias_action(data_type, strand_bias_action);
    let (sb_filtered_var, sb_matrix) = match bam_file {
        Some(bam) => {
            let (strand_map, fwd_frac) = build_strand_map(bam);
            println!("Library forward-strand fraction: {:.3}", fwd_frac);
            println!("Strand-bias action: {:?} (SNPs below HF {:.3})", action, strand_bias_drop_max_hf);
            filter_strand_bias(
                &var_record,
                &matrix,
                &read_record,
                &strand_map,
                &coverage,
                fwd_frac,
                strand_bias_threshold,
                2,
                permutation_frequency_threshold,
                action,
                strand_bias_drop_max_hf,
            )
        }
        None => (var_record.clone(), matrix.clone()),
    };

    let (permu_filtered_var, filtered_matrix, _filtered_name) = permutation_test(&sb_matrix, p_value_threshold, heteroplasmic_frequency_threshold, permutation_rounds, &sb_filtered_var, &coverage, permutation_frequency_threshold, data_type);

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

    /// A genuinely circular graph: the last anchor links straight back to the
    /// first, so the closing bubble has `dst.pos < src.pos`. Real Himito graphs
    /// look like this (e.g. A016544 -> A000022 on rCRS); they do NOT carry a
    /// duplicated anchor span past `ref_length`.
    fn circular_test_graph() -> GraphicalGenome {
        let mut anchor = HashMap::new();
        anchor.insert("A_first".to_string(), json!({"pos": 2}));
        anchor.insert("A_mid".to_string(), json!({"pos": 5}));
        anchor.insert("A_last".to_string(), json!({"pos": 8}));

        let mut edges = HashMap::new();
        edges.insert(
            "E_mid".to_string(),
            json!({"src": ["A_first"], "dst": ["A_mid"], "reads": ["r1"], "variants": "3="}),
        );
        // Closing bubble: spans 8 -> 10 (== ref_length) -> 2, i.e. across the origin.
        edges.insert(
            "E_wrap".to_string(),
            json!({"src": ["A_last"], "dst": ["A_first"], "reads": ["rw1", "rw2"], "variants": "2=1X1="}),
        );

        GraphicalGenome {
            anchor,
            edges,
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        }
    }

    #[test]
    fn origin_spanning_bubble_is_retained_and_covers_its_reads() {
        let graph = circular_test_graph();
        let ref_length = 10usize;
        let dict = get_graph_intervals(&graph, ref_length as i64);

        assert!(
            dict.contains_key(&"E_wrap".to_string()),
            "the origin-spanning bubble must not be dropped; variants there become uncallable"
        );

        // get_variants_from_cigar walks forward from the src anchor, so a variant
        // in this bubble carries a RAW position past ref_length (8 + 2 = 10),
        // which circuliarize_variants folds back to wrapped position 0.
        let raw_pos = 10usize;
        let mut cover: Vec<String> =
            bubble_cover_reads(&dict, raw_pos, &graph).into_iter().collect();
        cover.sort();
        assert_eq!(cover, vec!["rw1".to_string(), "rw2".to_string()]);

        // The pre-origin half of the same bubble resolves too.
        let mut cover_before: Vec<String> =
            bubble_cover_reads(&dict, 8, &graph).into_iter().collect();
        cover_before.sort();
        assert_eq!(cover_before, vec!["rw1".to_string(), "rw2".to_string()]);

        // ... and it must not bleed into the unrelated bubble at [2, 5).
        let mut cover_mid: Vec<String> =
            bubble_cover_reads(&dict, 3, &graph).into_iter().collect();
        cover_mid.sort();
        assert_eq!(cover_mid, vec!["r1".to_string()]);
    }

    #[test]
    fn source_and_sink_edges_do_not_get_genome_scale_intervals() {
        let mut anchor = HashMap::new();
        anchor.insert("A1".to_string(), json!({"pos": 900}));
        anchor.insert("A2".to_string(), json!({"pos": 910}));

        let mut edges = HashMap::new();
        edges.insert(
            "E_real".to_string(),
            json!({"src": ["A1"], "dst": ["A2"], "reads": ["r1"], "variants": "5=1X4="}),
        );
        // SOURCE/SINK have no anchor coordinate. Defaulting the missing endpoint
        // to 0/ref_length would give these edges a whole-genome interval and
        // credit their reads as covering every variant in the graph.
        edges.insert(
            "E_source".to_string(),
            json!({"src": ["SOURCE"], "dst": ["A2"], "reads": ["rs"], "variants": "5="}),
        );
        edges.insert(
            "E_sink".to_string(),
            json!({"src": ["A1"], "dst": ["SINK"], "reads": ["rk"], "variants": "5="}),
        );

        let graph = GraphicalGenome {
            anchor,
            edges,
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        let dict = get_graph_intervals(&graph, 16569);

        assert!(!dict.contains_key(&"E_source".to_string()));
        assert!(!dict.contains_key(&"E_sink".to_string()));

        // A distant position must pick up nothing at all.
        assert!(bubble_cover_reads(&dict, 42, &graph).is_empty());
        // The real bubble still resolves to exactly its own reads.
        let cover: Vec<String> = bubble_cover_reads(&dict, 905, &graph).into_iter().collect();
        assert_eq!(cover, vec!["r1".to_string()]);
    }

    #[test]
    fn construct_matrix_keeps_reads_with_few_alt_calls() {
        // Reads are the cells SCITE reasons over; dropping the ones carrying few
        // ALT calls deletes the reference-like reads that anchor the root of the
        // lineage, and on noisy data selects on error load rather than lineage.
        let variants: Vec<Variant> = ["A", "C"]
            .iter()
            .enumerate()
            .map(|(i, alt)| Variant {
                pos: 100 + i,
                ref_allele: "T".to_string(),
                alt_allele: alt.to_string(),
                variant_type: "SNP".to_string(),
                allele_count: 1,
                filter: None,
            })
            .collect();
        let n0 = generate_variant_name(&variants[0]);
        let n1 = generate_variant_name(&variants[1]);

        // r_alt2 carries both ALTs, r_alt1 carries one, r_ref carries none.
        let mut read_record: HashMap<String, Vec<serde_json::Value>> = HashMap::new();
        read_record.insert(n0.clone(), vec![json!("r_alt2"), json!("r_alt1")]);
        read_record.insert(n1.clone(), vec![json!("r_alt2")]);

        let all: HashSet<String> = ["r_alt2", "r_alt1", "r_ref"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut cover_record = HashMap::new();
        cover_record.insert(n0.clone(), all.clone());
        cover_record.insert(n1.clone(), all.clone());

        let (matrix, _vars, reads) = construct_matrix(&read_record, &cover_record, &variants);

        let mut sorted = reads.clone();
        sorted.sort();
        assert_eq!(
            sorted,
            vec!["r_alt1".to_string(), "r_alt2".to_string(), "r_ref".to_string()],
            "every covered read must survive, including the all-reference one"
        );
        assert_eq!(matrix.ncols(), 3);

        let col = |name: &str| reads.iter().position(|r| r == name).unwrap();
        assert_eq!(matrix[[0, col("r_ref")]], 0.0);
        assert_eq!(matrix[[1, col("r_ref")]], 0.0);
        assert_eq!(matrix[[0, col("r_alt1")]], 1.0);
        assert_eq!(matrix[[1, col("r_alt1")]], 0.0);
        assert_eq!(matrix[[0, col("r_alt2")]], 1.0);
        assert_eq!(matrix[[1, col("r_alt2")]], 1.0);
    }

    // ---------------------------------------------------------------------
    // Per-variant permutation null + empirical p-value
    // ---------------------------------------------------------------------

    fn snp(pos: usize, alt: &str, allele_count: usize) -> Variant {
        Variant {
            pos,
            ref_allele: "T".to_string(),
            alt_allele: alt.to_string(),
            variant_type: "SNP".to_string(),
            allele_count,
            filter: None,
        }
    }

    fn indel(pos: usize, ref_allele: &str, alt: &str, allele_count: usize) -> Variant {
        Variant {
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt.to_string(),
            variant_type: if ref_allele.len() > alt.len() { "DEL" } else { "INS" }.to_string(),
            allele_count,
            filter: None,
        }
    }

    #[test]
    fn empirical_p_value_counts_null_at_or_above_observation() {
        // (1 + #{null >= obs}) / (R + 1): the +1 keeps p strictly positive, so a
        // variant can never be credited with more evidence than R rounds can supply.
        let null = vec![0.0, 1.0, 2.0, 3.0];
        // 3 of 4 null draws reach 1.0 => (1 + 3) / 5
        assert!((empirical_p_value(&null, 1.0) - 0.8).abs() < 1e-12);
        // nothing beats an observation above the whole null => floor of 1 / (R + 1)
        assert!((empirical_p_value(&null, 99.0) - 0.2).abs() < 1e-12);
        // everything beats an observation below the whole null
        assert!((empirical_p_value(&null, -1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn empirical_p_value_of_empty_null_is_one() {
        // No null draws means no evidence, which must never look significant.
        assert_eq!(empirical_p_value(&[], 42.0), 1.0);
    }

    #[test]
    fn null_distribution_is_keyed_per_variant_and_preserves_alt_count() {
        // The summary statistic is a SUM of Jaccard coefficients, so it scales with a
        // variant's alt-read count. Pooling every variant's draws into one null makes
        // the test a depth filter: a low-count variant is scored against a null built
        // mostly from high-count variants. Each variant must get its own null.
        // Three variants, so each null is a SUM over two comparisons. With only two
        // variants the two nulls would each reduce to the single symmetric Jaccard
        // J(v0, v1) and be equal by construction, testing nothing.
        let variants = vec![snp(100, "A", 8), snp(200, "C", 2), snp(300, "G", 10)];
        let n_reads = 20;
        let mut matrix = Array2::<f64>::zeros((3, n_reads));
        for j in 0..n_reads {
            matrix[[0, j]] = if j < 8 { 1.0 } else { 0.0 };  // 0.40 alt frequency
            matrix[[1, j]] = if j < 2 { 1.0 } else { 0.0 };  // 0.10
            matrix[[2, j]] = if j < 10 { 1.0 } else { 0.0 }; // 0.50, the shared anchor
        }

        let rounds = 64;
        let nulls = get_null_distribution(&variants, &matrix, rounds, 0.9);

        let n0 = generate_variant_name(&variants[0]);
        let n1 = generate_variant_name(&variants[1]);
        assert_eq!(nulls.len(), 3, "one null per variant");
        assert_eq!(nulls[&n0].len(), rounds, "one draw per round for v0");
        assert_eq!(nulls[&n1].len(), rounds, "one draw per round for v1");

        // Shuffling preserves each variant's own alt count, so the two nulls live on
        // different scales -- exactly the confound that pooling them destroyed.
        let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
        assert!(
            mean(&nulls[&n0]) > mean(&nulls[&n1]),
            "the 8-read variant's null must sit above the 2-read variant's"
        );
    }

    #[test]
    fn null_distribution_skips_variants_above_the_frequency_threshold() {
        // Near-homoplasmic variants are the co-occurrence backbone, not candidates;
        // they are excluded from the null exactly as before the per-variant rewrite.
        let variants = vec![snp(100, "A", 9), snp(200, "C", 2)];
        let mut matrix = Array2::<f64>::zeros((2, 10));
        for j in 0..10 {
            matrix[[0, j]] = if j < 9 { 1.0 } else { 0.0 };
            matrix[[1, j]] = if j < 2 { 1.0 } else { 0.0 };
        }

        let nulls = get_null_distribution(&variants, &matrix, 8, 0.7);

        assert!(
            !nulls.contains_key(&generate_variant_name(&variants[0])),
            "a variant at 0.9 alt frequency must not get a null at threshold 0.7"
        );
        assert!(nulls.contains_key(&generate_variant_name(&variants[1])));
    }

    // ---------------------------------------------------------------------
    // SNP admission into the permutation test
    // ---------------------------------------------------------------------

    /// Two co-segregating variants carried by the same 3 of 30 reads (a real lineage),
    /// plus one SNP scattered across 3 reads that share no lineage (an artifact).
    fn lineage_and_artifact() -> (Vec<Variant>, Array2<f64>, HashMap<usize, usize>) {
        let variants = vec![
            snp(100, "A", 3),                 // lineage member
            indel(200, "TT", "T", 3),         // lineage member (indel: always tested)
            snp(300, "G", 3),                 // scattered artifact
        ];
        let n_reads = 30;
        let mut matrix = Array2::<f64>::zeros((3, n_reads));
        for j in 0..n_reads {
            // reads 0,1,2 carry the lineage
            let lin = if j < 3 { 1.0 } else { 0.0 };
            matrix[[0, j]] = lin;
            matrix[[1, j]] = lin;
            // reads 10,17,24 carry the artifact -- disjoint from the lineage
            matrix[[2, j]] = if j == 10 || j == 17 || j == 24 { 1.0 } else { 0.0 };
        }
        let coverage: HashMap<usize, usize> =
            [(100, n_reads), (200, n_reads), (300, n_reads)].into_iter().collect();
        (variants, matrix, coverage)
    }

    #[test]
    fn permutation_test_evaluates_snps_for_ont_denoised() {
        // Denoising strips the RANDOM error its quality model can see and leaves the
        // SYSTEMATIC error behind, so a denoised ONT SNP still needs the co-occurrence
        // test. Skipping SNPs here left `ont-denoised` with no SNP filter at all:
        // `filter_strand_bias` only flags, and `Potential_Artifact` is indel-only.
        let (variants, matrix, coverage) = lineage_and_artifact();
        let artifact = generate_variant_name(&variants[2]);

        let (kept, _m, _names) = permutation_test(
            &matrix, 0.5, 0.9, 200, &variants, &coverage, 0.9, "ont-denoised",
        );
        let kept_names: Vec<String> = kept.iter().map(generate_variant_name).collect();

        assert!(
            !kept_names.contains(&artifact),
            "a scattered SNP must be testable, and rejected, under ont-denoised"
        );
    }

    #[test]
    fn permutation_test_still_exempts_snps_for_pacbio() {
        // PacBio SNP accuracy is high enough that the co-occurrence test costs more
        // real low-VAF heteroplasmies than it saves false positives. Unchanged.
        let (variants, matrix, coverage) = lineage_and_artifact();
        let artifact = generate_variant_name(&variants[2]);

        let (kept, _m, _names) = permutation_test(
            &matrix, 0.5, 0.9, 200, &variants, &coverage, 0.9, "pacbio",
        );
        let kept_names: Vec<String> = kept.iter().map(generate_variant_name).collect();

        assert!(
            kept_names.contains(&artifact),
            "pacbio must keep SNPs exempt from the permutation test"
        );
    }

    #[test]
    fn permutation_test_keeps_a_variant_with_no_null_draws() {
        // `get_null_distribution` filters on the MATRIX alt frequency while the test
        // loop filters on allele_count/coverage. The two can disagree, leaving a
        // tested variant with an empty null. That is missing evidence, not evidence
        // of an artifact, so the variant must be kept untested rather than dropped.
        let variants = vec![snp(100, "A", 2), snp(200, "C", 2)];
        let mut matrix = Array2::<f64>::zeros((2, 10));
        for j in 0..10 {
            // Both variants are alt in 9/10 reads => matrix frequency 0.9, above the
            // 0.7 null threshold, so neither gets a null...
            matrix[[0, j]] = if j < 9 { 1.0 } else { 0.0 };
            matrix[[1, j]] = if j < 9 { 1.0 } else { 0.0 };
        }
        // ...but allele_count/coverage is 2/100 = 0.02, well under the test's own
        // 0.5 heteroplasmic threshold, so both are admitted for testing.
        let coverage: HashMap<usize, usize> = [(100, 100), (200, 100)].into_iter().collect();

        let (kept, _m, _names) = permutation_test(
            &matrix, 0.5, 0.5, 16, &variants, &coverage, 0.7, "ont-denoised",
        );

        assert_eq!(kept.len(), 2, "variants with an empty null must be kept, not dropped");
    }

    // ---------------------------------------------------------------------
    // Strand-bias filter: drop vs flag
    // ---------------------------------------------------------------------

    /// One low-VAF SNP supported by 10 reads, all on the forward strand, against a
    /// library that is 50/50 -- the signature of a systematic ONT artifact.
    fn strand_skewed_snp() -> (
        Vec<Variant>,
        Array2<f64>,
        HashMap<String, Vec<serde_json::Value>>,
        HashMap<String, bool>,
        HashMap<usize, usize>,
    ) {
        let variants = vec![snp(100, "A", 10)];
        let name = generate_variant_name(&variants[0]);
        let matrix = Array2::<f64>::zeros((1, 10));

        let reads: Vec<serde_json::Value> =
            (0..10).map(|i| json!(format!("r{i}"))).collect();
        let read_record: HashMap<String, Vec<serde_json::Value>> =
            [(name, reads)].into_iter().collect();
        // every supporting read is forward
        let strand_map: HashMap<String, bool> =
            (0..10).map(|i| (format!("r{i}"), false)).collect();
        let coverage: HashMap<usize, usize> = [(100, 200)].into_iter().collect();
        (variants, matrix, read_record, strand_map, coverage)
    }

    #[test]
    fn strand_bias_drops_low_vaf_snps_under_the_drop_action() {
        let (variants, matrix, read_record, strand_map, coverage) = strand_skewed_snp();

        let (kept, kept_matrix) = filter_strand_bias(
            &variants, &matrix, &read_record, &strand_map, &coverage,
            0.5, 0.01, 2, 0.7, StrandBiasAction::Drop, 0.10,
        );

        assert!(kept.is_empty(), "a strand-skewed low-VAF SNP must be dropped, not flagged");
        assert_eq!(kept_matrix.nrows(), 0, "the matrix must lose the dropped row");
    }

    #[test]
    fn strand_bias_only_flags_under_the_flag_action() {
        let (variants, matrix, read_record, strand_map, coverage) = strand_skewed_snp();

        let (kept, _m) = filter_strand_bias(
            &variants, &matrix, &read_record, &strand_map, &coverage,
            0.5, 0.01, 2, 0.7, StrandBiasAction::Flag, 0.10,
        );

        assert_eq!(kept.len(), 1, "flag mode must keep the call");
        assert_eq!(kept[0].filter.as_deref(), Some("Strand_bias"));
    }

    #[test]
    fn strand_bias_drop_action_spares_indels_and_higher_vaf_snps() {
        // The drop is deliberately narrow: it targets the low-VAF SNPs that denoising
        // leaves behind. Indels already have `Potential_Artifact`, and a SNP above the
        // drop ceiling is common enough that strand skew is likelier to be fragmented
        // read attribution than a single-strand artifact.
        let (_v, _m, read_record_snp, strand_map, _c) = strand_skewed_snp();

        // Same evidence, but the call sits at 10/50 = 0.20 VAF, above the 0.10 ceiling.
        let variants = vec![snp(100, "A", 10)];
        let matrix = Array2::<f64>::zeros((1, 10));
        let coverage: HashMap<usize, usize> = [(100, 50)].into_iter().collect();
        let (kept, _m) = filter_strand_bias(
            &variants, &matrix, &read_record_snp, &strand_map, &coverage,
            0.5, 0.01, 2, 0.7, StrandBiasAction::Drop, 0.10,
        );
        assert_eq!(kept.len(), 1, "a SNP above the drop ceiling is flagged, not dropped");
        assert_eq!(kept[0].filter.as_deref(), Some("Strand_bias"));

        // An indel with identical strand evidence at the same low VAF survives too.
        let ivariants = vec![indel(100, "TT", "T", 10)];
        let iname = generate_variant_name(&ivariants[0]);
        let ireads: Vec<serde_json::Value> = (0..10).map(|i| json!(format!("r{i}"))).collect();
        let iread_record: HashMap<String, Vec<serde_json::Value>> =
            [(iname, ireads)].into_iter().collect();
        let icoverage: HashMap<usize, usize> = [(100, 200)].into_iter().collect();
        let (ikept, _m) = filter_strand_bias(
            &ivariants, &matrix, &iread_record, &strand_map, &icoverage,
            0.5, 0.01, 2, 0.7, StrandBiasAction::Drop, 0.10,
        );
        assert_eq!(ikept.len(), 1, "indels are never dropped by the strand-bias filter");
    }

    #[test]
    fn resolve_strand_bias_action_drops_only_for_denoised_ont() {
        // PacBio and raw ONT keep the historical flag-only behaviour; the drop is
        // opt-in for the one data type whose SNPs have no other filter.
        assert_eq!(resolve_strand_bias_action("ont-denoised", None), StrandBiasAction::Drop);
        assert_eq!(resolve_strand_bias_action("pacbio", None), StrandBiasAction::Flag);
        assert_eq!(resolve_strand_bias_action("ont-r10", None), StrandBiasAction::Flag);
        // An explicit CLI choice always wins.
        assert_eq!(
            resolve_strand_bias_action("ont-denoised", Some(StrandBiasAction::Flag)),
            StrandBiasAction::Flag
        );
    }
}