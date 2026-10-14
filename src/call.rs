use crate::agg;
use agg::*;
use crate::build;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::{collections::HashSet, path::PathBuf};
use std::error::Error;
use csv::Writer;
use ndarray::{Array2, Axis, s};
use statrs::distribution::{Binomial, DiscreteCDF, Poisson};
use adjustp::{adjust, Procedure};
use bio::io::fasta::{Reader, Record};
use crate::denoise_indel::{normalize_left, Allele};
use crate::linkage::{self, LinkCall};
use rust_htslib::bam::{Read as HtsRead, Reader as BamReader, Record as BamRecord};
use rust_htslib::bam::ext::BamRecordExtensions;

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
        // Edges without a variant record were never genotyped; a read crossing a bubble only via one does not cover it.
        let has_variant = graph.edges[edge]
            .get("variants")
            .and_then(|v| v.as_str())
            .map_or(false, |s| !s.is_empty());
        if !has_variant {
            continue;
        }
        let src = graph.edges[edge].get("src").unwrap().as_array().unwrap()[0].as_str().unwrap();
        let dst = graph.edges[edge].get("dst").unwrap().as_array().unwrap()[0].as_str().unwrap();
        // SOURCE/SINK have no anchor coordinate, and get_variant skips those edges outright.
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
        // The graph is circular: the final anchor links straight back to the first.
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

/// Populates `read_record`/`cover_record` for one edge's variants: the raw position in `variants` finds the
/// bubble, the matching `variants_circular` entry names the record.
/// `graph_intervals_dict` holds raw (unwrapped) anchor coordinates: the closing bubble that crosses
/// the origin is stored as `[src.pos, dst.pos + ref_length)` so it stays ordered and half-open.
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

/// `left_align`: the reference to left-align indels against; `None` keeps the spelling of the graph.
pub fn get_variant(
    graph: &mut GraphicalGenome,
    k: usize,
    ref_name: &str,
    ref_length: usize,
    left_align: Option<&[u8]>,
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
        // Align BEFORE `var` is extended and before the read records are keyed.
        let variants_circular = circuliarize_variants(variants.clone(), ref_length);
        let variants_circular = match left_align {
            Some(refseq) => variants_circular.iter().map(|v| left_align_variant(v, refseq)).collect(),
            None => variants_circular,
        };
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

/// Shift an indel maximally leftward so equivalent alignments get one canonical representation.
fn left_align_variant(variant: &Variant, refseq: &[u8]) -> Variant {
    // `normalize_left` indexes an indel by its anchor: the last reference base BEFORE the event.
    if variant.pos == 0
        || variant.ref_allele.is_empty()
        || variant.alt_allele.is_empty()
        || variant.ref_allele == "-"
        || variant.alt_allele == "-"
    {
        return variant.clone();
    }

    let allele = match variant.variant_type.as_str() {
        // ALT is the anchor base followed by the inserted sequence.
        "INS" => Allele::Ins(variant.alt_allele.as_bytes()[1..].to_vec()),
        // REF is the anchor base followed by the deleted sequence.
        "DEL" => Allele::Del((variant.ref_allele.len() - 1) as u32),
        _ => return variant.clone(),
    };

    let (anchor, normalized) = normalize_left(refseq, (variant.pos - 1) as u32, &allele);
    let anchor = anchor as usize;

    // Rebuild both alleles from the reference at the NEW anchor rather than editing the old strings.
    let anchor_base = match refseq.get(anchor..anchor + 1) {
        Some(base) => String::from_utf8_lossy(base).to_string(),
        None => return variant.clone(),
    };
    let (ref_allele, alt_allele) = match &normalized {
        Allele::Ins(inserted) => (
            anchor_base.clone(),
            format!("{}{}", anchor_base, String::from_utf8_lossy(inserted)),
        ),
        Allele::Del(n) => match refseq.get(anchor..anchor + 1 + *n as usize) {
            Some(deleted) => (String::from_utf8_lossy(deleted).to_string(), anchor_base),
            None => return variant.clone(),
        },
        Allele::Ref => return variant.clone(),
    };

    Variant { pos: anchor + 1, ref_allele, alt_allele, ..variant.clone() }
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
    // Combine any upstream flag with the indel-artifact check here, rather than letting one filter suppress the other.
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

/// Folds accepted linkage calls into the graph calls.
pub(crate) fn merge_link_calls(variants: &[Variant], links: &[LinkCall]) -> (Vec<Variant>, Vec<LinkCall>) {
    let key = |pos: usize, r: &str, a: &str| (pos, r.to_ascii_uppercase(), a.to_ascii_uppercase());
    let by_allele: HashMap<(usize, String, String), &LinkCall> = links
        .iter()
        .map(|l| (key(l.pos, &l.ref_allele, &l.alt_allele), l))
        .collect();
    let mut seen: HashSet<(usize, String, String)> = HashSet::new();
    let mut merged = variants.to_vec();
    for v in merged.iter_mut().filter(|v| v.variant_type == "SNP") {
        let k = key(v.pos, &v.ref_allele, &v.alt_allele);
        if let Some(l) = by_allele.get(&k) {
            if let Some(old) = v.filter.take() {
                println!(
                    "Linkage rescue cleared FILTER {} on {} (partner {}, log10p {:.2})",
                    old,
                    generate_variant_name(v),
                    l.partner_name,
                    l.log10p
                );
            }
            seen.insert(k);
        }
    }
    let new_calls = links
        .iter()
        .filter(|l| !seen.contains(&key(l.pos, &l.ref_allele, &l.alt_allele)))
        .cloned()
        .collect();
    (merged, new_calls)
}

/// A rescued SNV as a matrix/genotyping `Variant`.
fn link_variant(l: &LinkCall) -> Variant {
    Variant {
        pos: l.pos,
        ref_allele: l.ref_allele.clone(),
        alt_allele: l.alt_allele.clone(),
        variant_type: "SNP".to_string(),
        allele_count: l.ac as usize,
        filter: None,
    }
}

/// A rescued SNV record; DP/AD/HF are pileup counts.
fn format_link_record(l: &LinkCall) -> String {
    let hf = if l.dp == 0 { 0.0 } else { l.ac as f32 / l.dp as f32 };
    format!(
        "chrM\t{}\t.\t{}\t{}\t.\tPASS\tDP={};LINK={},{:.2}\tGT:AD:HF\t1:{}:{}",
        l.pos + 1,
        l.ref_allele,
        l.alt_allele,
        l.dp,
        l.partner_name,
        l.log10p,
        l.ac,
        hf
    )
}

pub(crate) fn write_vcf(
    variants: &[Variant],
    coverage: &HashMap<usize, usize>,
    output_file: &PathBuf,
    sample_id: &str,
    referencename:&str,
    referencelength:usize,
    indel_false_threshold: f64,
    left_align: bool,
    links: Option<&[LinkCall]>,
) -> std::io::Result<()> {
    let mut file = File::create(Path::new(output_file))?;

    // Write VCF header
    writeln!(file, "##fileformat=VCFv4.2")?;
    writeln!(file, "##reference={}", referencename)?;
    writeln!(file, "##contig=<ID={},length={}>", referencename, referencelength)?;
    if left_align {
        writeln!(file, "##HimitoNormalization=left-aligned")?;
    }
    writeln!(
        file,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
    )?;
    if links.is_some() {
        writeln!(
            file,
            "##INFO=<ID=LINK,Number=2,Type=String,Description=\"Read-linkage rescue: the most strongly linked partner allele and the log10 hypergeometric p of their co-occurrence on co-covering reads (DP/AD/HF from the raw BAM pileup)\">"
        )?;
    }
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
        "##FILTER=<ID=Strand_bias,Description=\"Alt-supporting reads are significantly skewed toward one strand relative to the reads covering the site\">"
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

    let mut sorted_links: Vec<&LinkCall> = links.unwrap_or(&[]).iter().collect();
    sorted_links.sort_by(|a, b| a.pos.cmp(&b.pos).then(a.alt_allele.cmp(&b.alt_allele)));
    let vcf_pos = |v: &Variant| if v.variant_type == "SNP" { v.pos + 1 } else { v.pos };
    let mut li = 0;
    for variant in sorted_variants {
        while li < sorted_links.len() && sorted_links[li].pos + 1 < vcf_pos(&variant) {
            writeln!(file, "{}", format_link_record(sorted_links[li]))?;
            li += 1;
        }
        writeln!(file, "{}", format_vcf_record(&variant, coverage.clone(), indel_false_threshold))?;
    }
    for l in &sorted_links[li..] {
        writeln!(file, "{}", format_link_record(l))?;
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
        let alt = alt_reads(read_record, &name);
        let cover_reads = cover_record.get(&name).cloned().unwrap_or_default();

        for read in cover_reads {
            if let Some(c_index) = read_set_dict.get(&read) {
                matrix[[r_index, *c_index]] = if alt.contains(read.as_str()) {
                    1.0
                } else {
                    0.0
                };
            }
        }
    }

    (matrix, var_vec, read_vec)
}

// ---- BAM genotyping ----

/// `(start0, ref, alt)` where `start0` is the 0-based index of REF's first base on the reference.
fn allele_span(v: &Variant) -> (usize, &str, &str) {
    let start0 = if v.variant_type == "SNP" { v.pos } else { v.pos - 1 };
    (start0, v.ref_allele.as_str(), v.alt_allele.as_str())
}

/// Bases added on each side of a genotyping site before the homopolymer push-out.
const GENOTYPE_FLANK: usize = 10;
/// Extra read bases kept on each side of the aligned window.
const GENOTYPE_SLACK: usize = 4;
/// ONT noise gate: a read segment whose best edit distance exceeds this is unusable.
const GENOTYPE_MAX_ED_ONT: u32 = 4;

fn is_ont(data_type: &str) -> bool {
    matches!(data_type, "ont-r9" | "ont-r10" | "ont-denoised")
}

/// Reference window `[left, right)` of a site span `[start, end)`: widened by `GENOTYPE_FLANK`, pushed out of
/// any homopolymer it ends in, then padded by 2. `None` when it leaves the reference (no wrap-around).
fn genotype_window(refseq: &[u8], sites_span: (usize, usize)) -> Option<(usize, usize)> {
    let len = refseq.len() as isize;
    let mut left = sites_span.0 as isize - GENOTYPE_FLANK as isize;
    let mut right = sites_span.1 as isize + GENOTYPE_FLANK as isize;
    while 0 < left && left < len && refseq[left as usize] == refseq[left as usize - 1] {
        left -= 1;
    }
    while 0 < right && right < len && refseq[right as usize] == refseq[right as usize - 1] {
        right += 1;
    }
    let (left, right) = (left - 2, right + 2);
    if left < 0 || right > len {
        None
    } else {
        Some((left as usize, right as usize))
    }
}

/// Indices into `variants` grouped into genotyping sites by overlapping REF spans; SNPs and indels never share one.
fn group_sites(variants: &[Variant]) -> Vec<Vec<usize>> {
    let mut sites: Vec<Vec<usize>> = Vec::new();
    for want_indel in [false, true] {
        let mut idx: Vec<usize> = (0..variants.len())
            .filter(|&i| {
                let (_, r, a) = allele_span(&variants[i]);
                (r.len() != a.len()) == want_indel
            })
            .collect();
        idx.sort_by_key(|&i| allele_span(&variants[i]).0); // stable
        let mut cur: Vec<usize> = Vec::new();
        let mut cur_end = 0usize;
        for i in idx {
            let (start, r, _) = allele_span(&variants[i]);
            let end = start + r.len();
            if !cur.is_empty() && start < cur_end {
                cur.push(i);
                cur_end = cur_end.max(end);
            } else {
                if !cur.is_empty() {
                    sites.push(std::mem::take(&mut cur));
                }
                cur.push(i);
                cur_end = end;
            }
        }
        if !cur.is_empty() {
            sites.push(cur);
        }
    }
    sites
}

/// Outcome of scoring one read segment at one site.
#[derive(Debug, Clone, PartialEq, Eq)]
enum SegCall {
    /// Best distance exceeds `max_ed`: the read is too noisy here to say.
    TooNoisy,
    /// Indices into the site's haplotype list tied at the minimum distance.
    Best(BTreeSet<usize>),
}

fn unit_edit_cost(a: u8, b: u8) -> i32 {
    if a == b { 0 } else { -1 }
}

type EditAligner = bio::alignment::pairwise::Aligner<fn(u8, u8) -> i32>;

thread_local! {
    // One aligner per rayon worker thread (its DP buffers are reused).
    static EDIT_ALIGNER: std::cell::RefCell<EditAligner> = std::cell::RefCell::new(
        bio::alignment::pairwise::Aligner::with_capacity(
            64, 64, 0, -1, unit_edit_cost as fn(u8, u8) -> i32,
        )
    );
}

/// Infix edit distance of each haplotype against `seg`: the haplotype is aligned end to end.
fn score_segment(haps: &[Vec<u8>], seg: &[u8], max_ed: u32) -> SegCall {
    if haps.is_empty() {
        return SegCall::TooNoisy;
    }
    let dists: Vec<u32> = EDIT_ALIGNER.with(|al| {
        let mut al = al.borrow_mut();
        haps.iter()
            .map(|h| (-al.semiglobal(h, seg).score).max(0) as u32)
            .collect()
    });
    let best = *dists.iter().min().unwrap();
    if best > max_ed {
        return SegCall::TooNoisy;
    }
    SegCall::Best(
        dists
            .iter()
            .enumerate()
            .filter(|(_, &d)| d == best)
            .map(|(i, _)| i)
            .collect(),
    )
}

/// One genotyping site: a window on the reference and the haplotypes scored inside it.
struct GenoSite {
    left: usize,
    right: usize,
    members: Vec<usize>,
    /// REF, then the haplotype of each member in order.
    haps: Vec<Vec<u8>>,
}

/// Per (site, read) evidence merged over all of the read's usable alignments.
#[derive(Default)]
struct SiteAcc {
    /// Hap indices that were the unique best of some alignment.
    winners: BTreeSet<usize>,
    /// Union of the best sets of all usable alignments.
    any_best: BTreeSet<usize>,
}

fn build_geno_sites(refseq: &[u8], variants: &[Variant]) -> Vec<GenoSite> {
    let mut sites = Vec::new();
    for members in group_sites(variants) {
        let spans: Vec<(usize, &str, &str)> =
            members.iter().map(|&i| allele_span(&variants[i])).collect();
        let lo = spans.iter().map(|s| s.0).min().unwrap();
        let hi = spans.iter().map(|s| s.0 + s.1.len()).max().unwrap();
        let (left, right) = match genotype_window(refseq, (lo, hi)) {
            Some(w) => w,
            None => continue, // crosses the circular origin: left NaN
        };
        let mut haps = vec![refseq[left..right].to_vec()];
        for &(s, r, a) in &spans {
            let mut h = refseq[left..s].to_vec();
            h.extend_from_slice(a.as_bytes());
            h.extend_from_slice(&refseq[s + r.len()..right]);
            haps.push(h);
        }
        sites.push(GenoSite { left, right, members, haps });
    }
    sites.sort_by_key(|s| s.left);
    sites
}

/// Score one alignment against every site whose window it fully spans.
fn score_record_at_sites(
    rec: &BamRecord,
    sites: &[GenoSite],
    lefts: &[usize],
    max_ed: u32,
) -> Vec<(usize, SegCall)> {
    let mut out = Vec::new();
    let rs = rec.pos().max(0) as usize;
    let re = rec.reference_end().max(0) as usize;
    let k0 = lefts.partition_point(|&l| l < rs);
    if k0 >= sites.len() || sites[k0].left >= re {
        return out;
    }
    let pairs: Vec<[i64; 2]> = rec.aligned_pairs().collect(); // matches only: [qpos, rpos]
    if pairs.is_empty() {
        return out;
    }
    let seq = rec.seq().as_bytes().to_ascii_uppercase();
    for k in k0..sites.len() {
        let site = &sites[k];
        if site.left >= re {
            break;
        }
        if site.right > re {
            continue;
        }
        let a = pairs.partition_point(|p| (p[1] as usize) < site.left);
        let b = pairs.partition_point(|p| (p[1] as usize) < site.right);
        if b <= a {
            continue;
        }
        let q_lo = pairs[a..b].iter().map(|p| p[0]).min().unwrap() as usize;
        let q_hi = pairs[a..b].iter().map(|p| p[0]).max().unwrap() as usize;
        let qa = q_lo.saturating_sub(GENOTYPE_SLACK);
        let qb = (q_hi + 1 + GENOTYPE_SLACK).min(seq.len());
        out.push((k, score_segment(&site.haps, &seq[qa..qb], max_ed)));
    }
    out
}

/// Re-genotype `variants` in every read of `bam` by local realignment against the REF.
fn genotype_from_bam(
    bam: &Path,
    refseq: &[u8],
    ref_name: &str,
    variants: &[Variant],
    max_ed: u32,
) -> (HashMap<String, Vec<serde_json::Value>>, HashMap<String, HashSet<String>>) {
    const BATCH: usize = 4096;
    let sites = build_geno_sites(refseq, variants);
    let lefts: Vec<usize> = sites.iter().map(|s| s.left).collect();

    let mut qname_ids: HashMap<String, u32> = HashMap::new();
    let mut qnames: Vec<String> = Vec::new();
    let mut acc: HashMap<(u32, u32), SiteAcc> = HashMap::new();

    let mut reader = BamReader::from_path(bam).unwrap_or_else(|e| {
        eprintln!("Error: cannot open genotyping BAM {}: {}", bam.display(), e);
        std::process::exit(1);
    });
    // Unknown name → empty genotype maps (ref pseudo-read only), no exit.
    let target_tid: Option<i32> = match reader.header().tid(ref_name.as_bytes()) {
        Some(t) => Some(t as i32),
        None => {
            eprintln!(
                "Warning: reference contig '{}' not found in genotyping BAM {}; skipping all reads",
                ref_name,
                bam.display()
            );
            None
        }
    };

    let mut flush = |batch: &mut Vec<BamRecord>| {
        let scored: Vec<(String, Vec<(usize, SegCall)>)> = batch
            .par_iter()
            .map(|rec| {
                (
                    String::from_utf8_lossy(rec.qname()).to_string(),
                    score_record_at_sites(rec, &sites, &lefts, max_ed),
                )
            })
            .collect();
        batch.clear();
        for (qname, calls) in scored {
            if calls.is_empty() {
                continue;
            }
            let qid = *qname_ids.entry(qname.clone()).or_insert_with(|| {
                qnames.push(qname.clone());
                (qnames.len() - 1) as u32
            });
            for (k, call) in calls {
                if let SegCall::Best(set) = call {
                    let e = acc.entry((k as u32, qid)).or_default();
                    if set.len() == 1 {
                        e.winners.extend(set.iter().copied());
                    }
                    e.any_best.extend(set);
                }
            }
        }
    };

    let mut batch: Vec<BamRecord> = Vec::with_capacity(BATCH);
    for rec in reader.records() {
        let rec = rec.expect("failed to read genotyping BAM record");
        let Some(tid) = target_tid else { continue };
        if rec.is_unmapped() || rec.is_secondary() || rec.seq_len() == 0 || rec.tid() != tid {
            continue;
        }
        batch.push(rec);
        if batch.len() >= BATCH {
            flush(&mut batch);
        }
    }
    flush(&mut batch);

    // Merge per (site, read): one unique winner decides; several disagreeing winners -> missing;
    // no winner -> REF unless the variant tied for best somewhere (missing).
    let mut read_record: HashMap<String, Vec<serde_json::Value>> = HashMap::new();
    let mut cover_record: HashMap<String, HashSet<String>> = HashMap::new();
    for v in variants {
        let name = generate_variant_name(v);
        read_record.entry(name.clone()).or_default();
        cover_record.entry(name).or_default().insert(ref_name.to_string());
    }
    let mut keys: Vec<(u32, u32)> = acc.keys().copied().collect();
    keys.sort_unstable(); // deterministic output order
    for key in keys {
        let a = &acc[&key];
        let site = &sites[key.0 as usize];
        let qname = &qnames[key.1 as usize];
        let call_for = |j: usize| -> Option<bool> {
            if a.winners.is_empty() {
                if a.any_best.contains(&(j + 1)) { None } else { Some(false) }
            } else if a.winners.len() > 1 {
                None
            } else {
                Some(a.winners.contains(&(j + 1)))
            }
        };
        for (j, &vi) in site.members.iter().enumerate() {
            if let Some(is_alt) = call_for(j) {
                let name = generate_variant_name(&variants[vi]);
                cover_record.get_mut(&name).unwrap().insert(qname.clone());
                if is_alt {
                    read_record.get_mut(&name).unwrap().push(serde_json::json!(qname));
                }
            }
        }
    }
    (read_record, cover_record)
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

/// Build a map from read name (BAM qname) to strand, using primary, mapped alignments only.
fn build_strand_map(bam_file: &Path) -> (HashMap<String, bool>, f64) {
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

/// Two-sided binomial p-value that the alt reads' fwd/rev split deviates from the expected forward fraction.
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

/// Strand skew: the binomial test against `null_frac` is significant AND the forward fraction is at least
/// `STRAND_BIAS_MIN_EFFECT` away from it. Used by both denoisers; `filter_strand_bias` and the indel
/// recount apply the same test via `StrandEvidence::skew`, whose values they also log.
pub(crate) fn is_strand_skewed(fwd: usize, rev: usize, null_frac: f64, p_threshold: f64) -> bool {
    strand_bias_pvalue(fwd, rev, null_frac) < p_threshold
        && (fwd as f64 / (fwd + rev) as f64 - null_frac).abs() >= STRAND_BIAS_MIN_EFFECT
}

/// Minimum |alt forward fraction - null| for a significant strand split to count as skew.
const STRAND_BIAS_MIN_EFFECT: f64 = 0.20;

/// Fewer strand-resolved alt reads than this: too little evidence, the variant is kept untested.
const STRAND_MIN_READS: usize = 2;
/// Fewest strand-resolved covering reads for the site-local null; below it the library fraction is used.
const STRAND_LOCAL_MIN_READS: usize = 20;
/// Under `StrandBiasAction::Drop`, strand-skewed SNPs below this HF are dropped instead of flagged.
const STRAND_BIAS_DROP_MAX_HF: f64 = 0.10;

/// What `filter_strand_bias` does with a variant it finds strand-skewed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum StrandBiasAction {
    Flag,
    /// Additionally REMOVE strand-skewed low-VAF SNPs outright.
    Drop,
}

/// Names of the reads carrying `name`'s alt.
fn alt_reads<'a>(read_record: &'a HashMap<String, Vec<serde_json::Value>>, name: &str) -> HashSet<&'a str> {
    read_record
        .get(name)
        .map(|reads| reads.iter().filter_map(|r| r.as_str()).collect())
        .unwrap_or_default()
}

/// `(forward, reverse)` counts of `reads`; reads absent from `strand_map` are not counted.
fn strand_counts<'a>(reads: impl Iterator<Item = &'a str>, strand_map: &HashMap<String, bool>) -> (usize, usize) {
    let (mut fwd, mut rev) = (0, 0);
    for r in reads {
        match strand_map.get(r) {
            Some(true) => rev += 1,
            Some(false) => fwd += 1,
            None => {}
        }
    }
    (fwd, rev)
}

/// Strand split of a variant's alt reads, and the forward fraction an unbiased alt would show.
struct StrandEvidence {
    fwd: usize,
    rev: usize,
    null_frac: f64,
}

impl StrandEvidence {
    /// The null is the strand mix of the covering reads WITHOUT the alt, or `library_frac` when they are too few.
    fn new(
        alt: &HashSet<&str>,
        cover: Option<&HashSet<String>>,
        strand_map: &HashMap<String, bool>,
        library_frac: f64,
    ) -> Self {
        let (fwd, rev) = strand_counts(alt.iter().copied(), strand_map);
        let others = cover.into_iter().flatten().map(String::as_str).filter(|r| !alt.contains(r));
        let (cfwd, crev) = strand_counts(others, strand_map);
        let null_frac = if cfwd + crev >= STRAND_LOCAL_MIN_READS {
            cfwd as f64 / (cfwd + crev) as f64
        } else {
            library_frac
        };
        StrandEvidence { fwd, rev, null_frac }
    }

    fn total(&self) -> usize {
        self.fwd + self.rev
    }

    /// `(alt_frac, effect, p)`: the alt forward fraction, its distance from the null, and the binomial p.
    fn skew(&self) -> (f64, f64, f64) {
        let alt_frac = self.fwd as f64 / self.total() as f64;
        let effect = (alt_frac - self.null_frac).abs();
        (alt_frac, effect, strand_bias_pvalue(self.fwd, self.rev, self.null_frac))
    }
}

/// Flags variants whose alt reads are strand-skewed (significant and at least `STRAND_BIAS_MIN_EFFECT`);
/// under `Drop`, skewed SNPs below `STRAND_BIAS_DROP_MAX_HF` are removed instead.
/// Near-homoplasmic variants (HF >= `homoplasmic_frequency_threshold`) are kept untested.
fn filter_strand_bias(
    variants: &[Variant],
    matrix: &Array2<f64>,
    read_record: &HashMap<String, Vec<serde_json::Value>>,
    strand_map: &HashMap<String, bool>,
    cover_record: &HashMap<String, HashSet<String>>,
    coverage: &HashMap<usize, usize>,
    library_frac: f64,
    p_threshold: f64,
    homoplasmic_frequency_threshold: f64,
    action: StrandBiasAction,
) -> (Vec<Variant>, Array2<f64>) {
    let mut keep_idx = Vec::new();
    let mut kept = Vec::new();

    for (i, variant) in variants.iter().enumerate() {
        let mut variant = variant.clone();
        let name = generate_variant_name(&variant);

        let depth = coverage.get(&variant.pos).copied().unwrap_or(0);
        let hf = if depth == 0 { 0.0 } else { variant.allele_count as f64 / depth as f64 };
        if hf < homoplasmic_frequency_threshold {
            let ev = StrandEvidence::new(&alt_reads(read_record, &name), cover_record.get(&name), strand_map, library_frac);
            if ev.total() >= STRAND_MIN_READS {
                let (alt_frac, effect, p_value) = ev.skew();
                if p_value < p_threshold && effect >= STRAND_BIAS_MIN_EFFECT {
                    if action == StrandBiasAction::Drop
                        && variant.variant_type == "SNP"
                        && hf < STRAND_BIAS_DROP_MAX_HF
                    {
                        println!(
                            "Strand-biased (dropped), {:?}, fwd={}, rev={}, alt_frac={:.3}, null_frac={:.3}, hf={:.4}, p={:.3e}",
                            name, ev.fwd, ev.rev, alt_frac, ev.null_frac, hf, p_value
                        );
                        continue;
                    }
                    // Flag but keep: the VCF FILTER combines this with the other filters.
                    variant.filter = Some("Strand_bias".to_string());
                    println!(
                        "Strand-biased (flagged, kept), {:?}, fwd={}, rev={}, alt_frac={:.3}, null_frac={:.3}, p={:.3e}",
                        name, ev.fwd, ev.rev, alt_frac, ev.null_frac, p_value
                    );
                }
            }
        }
        keep_idx.push(i);
        kept.push(variant);
    }

    let filtered_matrix = matrix.select(Axis(0), &keep_idx);
    (kept, filtered_matrix)
}

fn is_strand_flagged_indel(v: &Variant) -> bool {
    v.variant_type != "SNP" && v.filter.as_deref() == Some("Strand_bias")
}

/// Retest `Strand_bias`-flagged indels on BAM-genotyped alt/cover reads; a balanced split clears the flag.
fn recount_indel_strand_bias(
    variants: &[Variant],
    bam_read_record: &HashMap<String, Vec<serde_json::Value>>,
    bam_cover_record: &HashMap<String, HashSet<String>>,
    strand_map: &HashMap<String, bool>,
    library_frac: f64,
    p_threshold: f64,
) -> Vec<Variant> {
    let mut out = variants.to_vec();
    for variant in out.iter_mut().filter(|v| is_strand_flagged_indel(v)) {
        let name = generate_variant_name(variant);
        let alt = alt_reads(bam_read_record, &name);
        let ev = StrandEvidence::new(&alt, bam_cover_record.get(&name), strand_map, library_frac);
        if ev.total() < STRAND_MIN_READS {
            println!(
                "Strand-bias BAM recount, {:?}, BAM alt reads={} < {}: kept flagged",
                name, ev.total(), STRAND_MIN_READS
            );
            continue;
        }
        let (alt_frac, effect, p_value) = ev.skew();
        let cleared = p_value >= p_threshold || effect < STRAND_BIAS_MIN_EFFECT;
        println!(
            "Strand-bias BAM recount, {:?}, BAM fwd={}, rev={}, alt_frac={:.3}, null_frac={:.3}, effect={:.3}, p={:.3e}: {}",
            name, ev.fwd, ev.rev, alt_frac, ev.null_frac, effect, p_value,
            if cleared { "cleared" } else { "kept flagged" }
        );
        if cleared {
            variant.filter = None;
        }
    }
    out
}

/// Genotype only the `Strand_bias`-flagged indels in `bam` and retest their strand split there.
fn recount_flagged_indels_from_bam(
    variants: &[Variant],
    bam: &Path,
    strand_map: &HashMap<String, bool>,
    library_frac: f64,
    refseq: &[u8],
    ref_name: &str,
    p_threshold: f64,
) -> Vec<Variant> {
    let flagged: Vec<Variant> = variants.iter().filter(|v| is_strand_flagged_indel(v)).cloned().collect();
    if flagged.is_empty() {
        return variants.to_vec();
    }
    // Ungated: the matrix noise gate must not move the VCF.
    let (read_record, cover_record) = genotype_from_bam(bam, refseq, ref_name, &flagged, u32::MAX);
    println!(
        "Strand-bias BAM recount: {} flagged indel(s) from {} (forward fraction {:.3})",
        flagged.len(),
        bam.display(),
        library_frac
    );
    recount_indel_strand_bias(variants, &read_record, &cover_record, strand_map, library_frac, p_threshold)
}

// ---- adaptive permutation ----

/// Stop drawing for a variant once this many null draws reached its observation.
const PERM_STOP_AFTER: u32 = 10;
/// Floor-resolving target: `PERM_HEADROOM * m / (alpha * k0)` rounds, so the floor p clears BH.
const PERM_HEADROOM: f64 = 5.0;
/// Cap on the rounds of any one variant.
const PERM_MAX_ROUNDS: usize = 100_000;
/// Base RNG seed; each variant's stream is this XOR a hash of its name.
const PERM_SEED: u64 = 42;
/// Draws per batch; the early-stop test runs between batches.
const PERM_BATCH: u32 = 100;
/// Most floor-resolving passes (each only ever raises the target).
const PERM_MAX_STAGE2_PASSES: usize = 5;
/// Majority-lineage exemption: SNPs at or above this HF and `COOC_EXEMPT_MIN_AC` skip the test.
const COOC_EXEMPT_HF: f64 = 0.10;
const COOC_EXEMPT_MIN_AC: usize = 20;

/// PacBio hot-allele error rate and noise-tail for the depth-aware AC floor T(DP).
const PB_HOT_ERR: f64 = 0.002;
const PB_NOISE_TAIL: f64 = 0.002;

/// Smallest k with P(Poisson(PB_HOT_ERR · DP) ≥ k) < PB_NOISE_TAIL, floored at 3.
/// Only for PacBio paths (see gray-zone admission in T2).
pub fn pb_noise_threshold(dp: usize) -> usize {
    let lambda = PB_HOT_ERR * dp as f64;
    // λ = 0 → Poisson is a point mass at 0; T collapses to the floor.
    if lambda == 0.0 {
        return 3;
    }
    let pois = Poisson::new(lambda).expect("PB_HOT_ERR * dp is a valid Poisson rate");
    let mut k = 1usize;
    loop {
        // P(X >= k) = 1 - F(k-1)
        let tail = 1.0 - pois.cdf((k - 1) as u64);
        if tail < PB_NOISE_TAIL {
            return k.max(3);
        }
        k += 1;
        // Defensive: λ·DP is tiny in practice; never spin forever.
        if k > dp.max(3) + 1000 {
            return k.max(3);
        }
    }
}

/// A genotype row as bitsets over reads.
struct BitRow {
    alt: Vec<u64>,
    cov: Vec<u64>,
    /// Indices of the covered reads: the pool a shuffle draws its alt reads from.
    covered_idx: Vec<u32>,
    n_alt: u32,
}

impl BitRow {
    fn from_row(row: &[f64]) -> BitRow {
        let words = (row.len() + 63) / 64;
        let mut alt = vec![0u64; words];
        let mut cov = vec![0u64; words];
        let mut covered_idx = Vec::new();
        let mut n_alt = 0u32;
        for (j, &v) in row.iter().enumerate() {
            if v.is_nan() {
                continue;
            }
            cov[j >> 6] |= 1u64 << (j & 63);
            covered_idx.push(j as u32);
            if v > 0.5 {
                alt[j >> 6] |= 1u64 << (j & 63);
                n_alt += 1;
            }
        }
        BitRow { alt, cov, covered_idx, n_alt }
    }
}

/// Jaccard of the alt reads of a (possibly shuffled) row `x` and row `y`, over the reads both cover.
#[inline]
fn jaccard_bits(x_alt: &[u64], x_cov: &[u64], y: &BitRow) -> f64 {
    let mut inter = 0u32;
    let mut union = 0u32;
    for w in 0..x_alt.len() {
        inter += (x_alt[w] & y.alt[w]).count_ones();
        union += ((x_alt[w] & y.cov[w]) | (y.alt[w] & x_cov[w])).count_ones();
    }
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

/// Mitochondrial reference length (circular) for partner distance.
const MT_LEN: usize = 16569;

#[inline]
fn circ_dist(a: usize, b: usize) -> usize {
    let d = a.abs_diff(b);
    let wrap = MT_LEN - d.min(MT_LEN);
    d.min(wrap)
}

/// Who may serve as a co-occurrence partner for the permutation statistic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PartnerRule {
    /// Every matrix row (today's ONT behaviour).
    All,
    /// PacBio: SNPs with AC ≥ T(DP) that are not `Strand_bias`-flagged.
    PbNoiseThreshold,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PermOpts {
    pub partner_snp_only: bool,
    pub partner_min_dist: usize,
    pub partner_rule: PartnerRule,
}

impl PermOpts {
    pub fn all_rows() -> Self {
        Self {
            partner_snp_only: false,
            partner_min_dist: 0,
            partner_rule: PartnerRule::All,
        }
    }
}

pub fn resolve_perm_opts(data_type: &str) -> PermOpts {
    if data_type == "pacbio" {
        PermOpts {
            partner_snp_only: true,
            partner_min_dist: 20,
            partner_rule: PartnerRule::PbNoiseThreshold,
        }
    } else {
        PermOpts::all_rows()
    }
}

fn is_snp_variant(v: &Variant) -> bool {
    v.ref_allele.len() == 1 && v.alt_allele.len() == 1
}

/// Indices that form the partner pool under `opts` (self / min-dist filtered later per query).
fn partner_pool(
    variants: &[Variant],
    coverage: &HashMap<usize, usize>,
    opts: PermOpts,
) -> Vec<usize> {
    (0..variants.len())
        .filter(|&i| {
            let v = &variants[i];
            let snp = is_snp_variant(v);
            if opts.partner_snp_only && !snp {
                return false;
            }
            match opts.partner_rule {
                PartnerRule::All => true,
                PartnerRule::PbNoiseThreshold => {
                    if !snp {
                        return false;
                    }
                    if v.filter.as_deref() == Some("Strand_bias") {
                        return false;
                    }
                    let depth = coverage.get(&v.pos).copied().unwrap_or(0);
                    v.allele_count >= pb_noise_threshold(depth)
                }
            }
        })
        .collect()
}

/// The co-occurrence statistic: summed Jaccard of row `x` against partners,
/// excluding `self_index` and rows within `min_dist` bp (circular, `MT_LEN`).
#[inline]
fn sum_jaccard_bits(
    x_alt: &[u64],
    x_cov: &[u64],
    self_index: usize,
    self_pos: usize,
    partners: &[usize],
    positions: &[usize],
    min_dist: usize,
    rows: &[BitRow],
) -> f64 {
    let mut sum = 0.0;
    for &j in partners {
        if j == self_index {
            continue;
        }
        if min_dist > 0 && circ_dist(self_pos, positions[j]) < min_dist {
            continue;
        }
        sum += jaccard_bits(x_alt, x_cov, &rows[j]);
    }
    sum
}

/// Shuffle a row's alt calls among its covered reads (uncovered reads stay uncovered) into `out_alt`.
fn shuffle_bits(
    row: &BitRow,
    scratch: &mut Vec<u32>,
    out_alt: &mut [u64],
    rng: &mut impl rand::Rng,
) {
    scratch.clear();
    scratch.extend_from_slice(&row.covered_idx);
    out_alt.fill(0);
    let len = scratch.len();
    for t in 0..row.n_alt as usize {
        let j = rng.random_range(t..len);
        scratch.swap(t, j);
        let r = scratch[t] as usize;
        out_alt[r >> 6] |= 1u64 << (r & 63);
    }
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Null draws of one tested variant: how many reached the observation, out of how many.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct PermOutcome {
    exceed: u32,
    draws: u32,
}

impl PermOutcome {
    fn p_value(&self) -> f64 {
        (self.exceed + 1) as f64 / (self.draws + 1) as f64
    }
}

struct PermState {
    rng: rand::rngs::StdRng,
    exceed: u32,
    draws: u32,
}

/// Permutation null of each `tested` row (one outcome each, same order). Stage 1 draws up to `min_rounds`
/// per variant, stopping early after `PERM_STOP_AFTER` exceedances; stage 2 extends the variants still at
/// the p floor (no exceedance) until that floor can pass BH at `alpha`.
///
/// `partners` is the partner pool (built once); `positions[i]` is the genomic position of matrix row `i`;
/// `min_dist` excludes partners closer than that many bp (circular).
fn adaptive_permutation(
    matrix: &Array2<f64>,
    tested: &[usize],
    names: &[String],
    partners: &[usize],
    positions: &[usize],
    min_dist: usize,
    min_rounds: usize,
    alpha: f64,
) -> Vec<PermOutcome> {
    use rand::SeedableRng;
    let started = std::time::Instant::now();
    let m = tested.len();
    if m == 0 || min_rounds == 0 {
        return vec![PermOutcome { exceed: 0, draws: 0 }; m];
    }
    let bits: Vec<BitRow> = matrix.outer_iter().map(|r| BitRow::from_row(&r.to_vec())).collect();
    let words = bits.first().map_or(0, |b| b.alt.len());
    let obs: Vec<f64> = tested
        .iter()
        .map(|&i| {
            sum_jaccard_bits(
                &bits[i].alt,
                &bits[i].cov,
                i,
                positions[i],
                partners,
                positions,
                min_dist,
                &bits,
            )
        })
        .collect();
    let mut states: Vec<PermState> = tested
        .iter()
        .map(|&i| PermState {
            rng: rand::rngs::StdRng::seed_from_u64(PERM_SEED ^ fnv1a64(names[i].as_bytes())),
            exceed: 0,
            draws: 0,
        })
        .collect();

    // Extend the variants in `which` (indices into `tested`) to `target` draws.
    let extend = |states: &mut Vec<PermState>, which: &[usize], target: usize| {
        let target = target as u32;
        // Disjoint &mut access to the chosen states, one rayon task per variant.
        let mut picked: Vec<(usize, &mut PermState)> = Vec::with_capacity(which.len());
        let mut want = which.iter().copied().peekable();
        for (k, st) in states.iter_mut().enumerate() {
            if want.peek() == Some(&k) {
                picked.push((k, st));
                want.next();
            }
        }
        picked.into_par_iter().for_each_init(
            || (Vec::<u32>::new(), vec![0u64; words]),
            |(scratch, buf), (k, st)| {
                let i = tested[k];
                let row = &bits[i];
                let threshold = obs[k] - 1e-9;
                while st.draws < target && st.exceed < PERM_STOP_AFTER {
                    let n = PERM_BATCH.min(target - st.draws);
                    for _ in 0..n {
                        shuffle_bits(row, scratch, buf, &mut st.rng);
                        if sum_jaccard_bits(
                            buf,
                            &row.cov,
                            i,
                            positions[i],
                            partners,
                            positions,
                            min_dist,
                            &bits,
                        ) >= threshold
                        {
                            st.exceed += 1;
                        }
                    }
                    st.draws += n;
                }
            },
        );
    };
    let total = |states: &Vec<PermState>| states.iter().map(|s| s.draws as u64).sum::<u64>();

    let all: Vec<usize> = (0..m).collect();
    extend(&mut states, &all, min_rounds);
    let mut target = min_rounds;
    let k0 = states.iter().filter(|s| s.exceed == 0).count();
    println!(
        "Permutation stage 1: m={} k0={} R_target={} draws={} time={:.1}s",
        m, k0, target, total(&states), started.elapsed().as_secs_f64()
    );

    for pass in 0..PERM_MAX_STAGE2_PASSES {
        let floor: Vec<usize> = (0..m).filter(|&k| states[k].exceed == 0).collect();
        if floor.is_empty() {
            break;
        }
        let new_target = PERM_MAX_ROUNDS
            .min((PERM_HEADROOM * m as f64 / (alpha * floor.len() as f64)).ceil() as usize);
        if new_target <= target {
            break;
        }
        target = new_target;
        extend(&mut states, &floor, target);
        println!(
            "Permutation stage 2.{}: m={} k0={} R_target={} draws={} time={:.1}s",
            pass + 1, m, floor.len(), target, total(&states), started.elapsed().as_secs_f64()
        );
    }

    states.iter().map(|s| PermOutcome { exceed: s.exceed, draws: s.draws }).collect()
}

/// Co-occurrence permutation test: drops admitted variants whose BH q exceeds `p_value_threshold`.
fn permutation_test(
    matrix: &Array2<f64>,
    p_value_threshold: f64,
    heteroplasmic_error_threshold: f64,
    permutation_round: usize,
    filtered_var: &[Variant],
    coverage: &HashMap<usize, usize>,
    permutation_frequency_threshold: f64,
    data_type: &str,
    indel_false_threshold: f64,
) -> (Vec<Variant>, Array2<f64>) {
    // Admission: which variants are candidates for the permutation test at all.
    let admitted = |current_variant: &Variant| -> bool {
        // A variant whose position carries no coverage entry cannot be scored.
        let depth = coverage.get(&current_variant.pos).copied().unwrap_or(0);
        let frequency = if depth == 0 {
            0.0
        } else {
            current_variant.allele_count as f64 / depth as f64
        };

        if frequency > heteroplasmic_error_threshold {
            return false;
        }
        let indel_len = (current_variant.ref_allele.len() as i32 - current_variant.alt_allele.len() as i32).abs();
        // exclude large indel
        if indel_len > 50 {
            return false;
        }

        let is_snp = is_snp_variant(current_variant);

        // Jaccard co-occurrence penalizes variants on the majority lineage, so exempt high-HF, high-AC variants.
        if is_snp
            && current_variant.filter.is_none()
            && frequency >= COOC_EXEMPT_HF
            && current_variant.allele_count >= COOC_EXEMPT_MIN_AC
        {
            return false;
        }

        // PacBio: test SNPs only in the gray zone 3 ≤ AC < T(DP); AC ≥ T(DP) are exempt.
        if is_snp && data_type == "pacbio" {
            return current_variant.allele_count < pb_noise_threshold(depth);
        }
        // PacBio: small indels below indel_false_threshold can never PASS, so skip the test.
        if !is_snp && data_type == "pacbio" && indel_len < 5 && frequency < indel_false_threshold {
            return false;
        }
        true
    };

    let names: Vec<String> = filtered_var.iter().map(generate_variant_name).collect();
    let positions: Vec<usize> = filtered_var.iter().map(|v| v.pos).collect();
    let opts = resolve_perm_opts(data_type);
    let partners = partner_pool(filtered_var, coverage, opts);
    // No null for a variant above `permutation_frequency_threshold` by MATRIX alt frequency.
    let tested: Vec<usize> = (0..filtered_var.len())
        .filter(|&i| {
            admitted(&filtered_var[i])
                && alt_frequency(matrix.slice(s![i, ..])) <= permutation_frequency_threshold
        })
        .collect();
    let outcomes = adaptive_permutation(
        matrix,
        &tested,
        &names,
        &partners,
        &positions,
        opts.partner_min_dist,
        permutation_round,
        p_value_threshold,
    );
    // A variant without draws is untested: missing evidence, not evidence of an artifact, so it skips BH.
    let (test_index, raw_p_values): (Vec<&str>, Vec<f64>) = tested
        .iter()
        .zip(&outcomes)
        .filter(|(_, o)| o.draws > 0)
        .map(|(&i, o)| (names[i].as_str(), o.p_value()))
        .unzip();
    let qvalues = adjust(&raw_p_values, Procedure::BenjaminiHochberg);
    let excluded: HashSet<&str> = test_index
        .into_iter()
        .zip(qvalues)
        .filter(|(_, q)| *q > p_value_threshold)
        .map(|(name, _)| name)
        .collect();

    let keep: Vec<usize> = (0..filtered_var.len()).filter(|&i| !excluded.contains(names[i].as_str())).collect();
    let kept = keep.iter().map(|&i| filtered_var[i].clone()).collect();
    (kept, matrix.select(Axis(0), &keep))
}

pub fn resolve_thresholds(
    data_type: &str,
    p_value_threshold: Option<f64>,
    frequency_threshold: Option<f64>,
    permutation_frequency_threshold: Option<f64>,
) -> (f64, f64, f64) {
    let (default_p, default_f, default_perm) = if data_type == "ont-r9" {
        (0.0001, 0.5, 0.7)
    } else {
        (0.01, 0.2, 0.7)
    };

    (
        p_value_threshold.unwrap_or(default_p),
        frequency_threshold.unwrap_or(default_f),
        permutation_frequency_threshold.unwrap_or(default_perm),
    )
}

pub fn resolve_indel_false_threshold(data_type: &str, explicit: Option<f64>) -> f64 {
    explicit.unwrap_or(if is_ont(data_type) { 0.015 } else { 0.1 })
}

/// Floor on heteroplasmic frequency for calling (`-v`). PacBio defaults to 0.002; others to 0.01.
pub fn resolve_vaf(data_type: &str, explicit: Option<f32>) -> f32 {
    explicit.unwrap_or(if data_type == "pacbio" { 0.002 } else { 0.01 })
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
    // Raw (pre-denoise) BAM of the same reads: strand map, linkage rescue and the matrix.
    bam: Option<&Path>,
    strand_bias_threshold: f64,
    indel_false_threshold: f64,
    permutation_frequency_threshold: f64,
    permutation_rounds: usize,
    left_align: bool,
    no_linkage_rescue: bool,
) {
    // reference fasta information
    let ref_reader = Reader::from_file(fasta_reference).unwrap();
    let reference_sequence: Vec<Record> = ref_reader.records().map(|r| r.unwrap()).collect();
    let ref_seq = String::from_utf8_lossy(reference_sequence[0].seq()).to_string();
    let ref_header = reference_sequence[0].id().to_string();

    // read in graphical genome
    let graph = agg::GraphicalGenome::load_graph(graph_file).unwrap();

    let (variants, coverage, read_record, cover_record) =
        get_variant(
            &mut graph.clone(),
            k,
            &ref_header,
            ref_seq.len(),
            if left_align { Some(ref_seq.as_bytes()) } else { None },
        );
    let collapsed_var = collapse_identical_records(variants, ref_seq.len());
    let filtered_var = filter_vcf_record(&collapsed_var, &coverage, minimal_ac, hf_threshold);
    // modified, exclude filtered data for FPs

    // modified, exclude filtered data
    let (matrix, var_record, read_set) =
        construct_matrix(&read_record, &cover_record, &filtered_var);
    let matrix_output_raw = output_file.with_extension("raw_matrix.csv");
    let _ = write_matrix_to_csv(&matrix, &var_record, &read_set, matrix_output_raw);

    // Strand-bias filtering runs before the permutation test so artifacts do not feed its null.
    let action = if data_type == "ont-denoised" { StrandBiasAction::Drop } else { StrandBiasAction::Flag };
    // Linkage rescue (and with it the indel strand recount) runs on ONT only; PacBio is not evaluated.
    let linkage_bam = bam.filter(|_| is_ont(data_type) && !no_linkage_rescue);
    let (sb_filtered_var, sb_matrix) = match bam {
        Some(bam) => {
            let (strand_map, fwd_frac) = build_strand_map(bam);
            println!("Library forward-strand fraction: {:.3}", fwd_frac);
            println!("Strand-bias action: {:?} (SNPs below HF {:.3})", action, STRAND_BIAS_DROP_MAX_HF);
            let (sb_var, sb_matrix) = filter_strand_bias(
                &var_record,
                &matrix,
                &read_record,
                &strand_map,
                &cover_record,
                &coverage,
                fwd_frac,
                strand_bias_threshold,
                permutation_frequency_threshold,
                action,
            );
            let sb_var = match linkage_bam {
                Some(raw) => recount_flagged_indels_from_bam(
                    &sb_var,
                    raw,
                    &strand_map,
                    fwd_frac,
                    ref_seq.as_bytes(),
                    &ref_header,
                    strand_bias_threshold,
                ),
                None => sb_var,
            };
            (sb_var, sb_matrix)
        }
        None => (var_record.clone(), matrix.clone()),
    };

    let (permu_filtered_var, filtered_matrix) = permutation_test(
        &sb_matrix,
        p_value_threshold,
        heteroplasmic_frequency_threshold,
        permutation_rounds,
        &sb_filtered_var,
        &coverage,
        permutation_frequency_threshold,
        data_type,
        indel_false_threshold,
    );

    let link_calls: Option<Vec<LinkCall>> = linkage_bam.and_then(|bam| {
        match linkage::rescue(bam, ref_seq.as_bytes(), &ref_header) {
            Ok((calls, _summary)) => Some(calls),
            Err(e) => {
                eprintln!("Warning: linkage rescue failed ({e:#}); continuing without it.");
                None
            }
        }
    });
    let (permu_filtered_var, new_links) = match &link_calls {
        Some(calls) => merge_link_calls(&permu_filtered_var, calls),
        None => (permu_filtered_var, Vec::new()),
    };
    for l in &new_links {
        println!(
            "LINK m.{}{}>{} AC {} DP {} ({}/{}) partner {} shared {} enrich {:.1} log10p {:.2}",
            l.pos + 1, l.ref_allele, l.alt_allele, l.ac, l.dp, l.fwd, l.rev,
            l.partner_name, l.shared, l.enrich, l.log10p
        );
    }
    let rescued_var: Vec<Variant> = new_links.iter().map(link_variant).collect();

    // write filtered vcf
    let _ = write_vcf(
        &permu_filtered_var,
        &coverage,
        output_file,
        sample_id,
        &ref_header,
        ref_seq.len(),
        indel_false_threshold,
        left_align,
        link_calls.as_ref().map(|_| new_links.as_slice()),
    );
    // The VCF and every filter above run on the graph; only the lineage matrix is BAM-genotyped.
    let genotype_src = bam.map(|bam| GenotypeSource {
        bam,
        refseq: ref_seq.as_bytes(),
        ref_name: &ref_header,
        max_ed: if is_ont(data_type) { GENOTYPE_MAX_ED_ONT } else { u32::MAX },
    });
    if let Err(e) = write_final_matrix(
        output_file,
        &filtered_matrix,
        &permu_filtered_var,
        &rescued_var,
        &read_set,
        genotype_src.as_ref(),
    ) {
        eprintln!("Error writing matrix: {e}");
    }
}

/// Where to re-genotype the final variants from.
struct GenotypeSource<'a> {
    bam: &'a Path,
    refseq: &'a [u8],
    ref_name: &'a str,
    max_ed: u32,
}

/// Writes `<output>.matrix.csv`: re-genotyped from the BAM when one is given, else the graph matrix.
fn write_final_matrix(
    output_file: &Path,
    graph_matrix: &Array2<f64>,
    variants: &[Variant],
    rescued: &[Variant],
    graph_reads: &[String],
    genotype: Option<&GenotypeSource>,
) -> Result<(), Box<dyn Error>> {
    let matrix_output = output_file.with_extension("matrix.csv");
    match genotype {
        None => {
            eprintln!("Warning: no --bam given; matrix.csv comes from the graph, not from re-genotyped reads.");
            write_matrix_to_csv(graph_matrix, variants, graph_reads, matrix_output)
        }
        Some(src) => {
            let all: Vec<Variant> = variants.iter().chain(rescued).cloned().collect();
            let (read_record, cover_record) =
                genotype_from_bam(src.bam, src.refseq, src.ref_name, &all, src.max_ed);
            let (matrix, var_record, read_set) =
                construct_matrix(&read_record, &cover_record, &all);
            println!(
                "Genotyped {} variants x {} reads from {}",
                var_record.len(),
                read_set.len(),
                src.bam.display()
            );
            write_matrix_to_csv(&matrix, &var_record, &read_set, matrix_output)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use serde_json::json;

    /// A graph from named anchors and edges (each edge JSON carries `src`/`dst`); adjacency follows the edges.
    fn graph(anchors: &[(&str, serde_json::Value)], edges: &[(&str, serde_json::Value)]) -> GraphicalGenome {
        let mut g = GraphicalGenome {
            anchor: anchors.iter().map(|(n, a)| (n.to_string(), a.clone())).collect(),
            edges: HashMap::new(),
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };
        for (name, e) in edges {
            let src = e["src"][0].as_str().unwrap().to_string();
            let dst = e["dst"][0].as_str().unwrap().to_string();
            g.outgoing.entry(src.clone()).or_default().push(name.to_string());
            g.incoming.entry(dst.clone()).or_default().push(name.to_string());
            g.incoming.insert(name.to_string(), vec![src]);
            g.outgoing.insert(name.to_string(), vec![dst]);
            g.edges.insert(name.to_string(), e.clone());
        }
        g
    }

    fn sorted(reads: HashSet<String>) -> Vec<String> {
        let mut v: Vec<String> = reads.into_iter().collect();
        v.sort();
        v
    }

    /// Typed the way `get_variants_from_cigar` types it: equal lengths SNP, longer REF DEL, else INS.
    fn var(pos: usize, ref_allele: &str, alt: &str, allele_count: usize) -> Variant {
        let variant_type = match ref_allele.len().cmp(&alt.len()) {
            std::cmp::Ordering::Equal => "SNP",
            std::cmp::Ordering::Greater => "DEL",
            std::cmp::Ordering::Less => "INS",
        };
        Variant {
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt.to_string(),
            variant_type: variant_type.to_string(),
            allele_count,
            filter: None,
        }
    }

    fn snp(pos: usize, alt: &str, allele_count: usize) -> Variant {
        var(pos, "T", alt, allele_count)
    }

    #[test]
    fn record_variant_reads_uses_raw_position_for_wraparound_bubble() {
        // Bubbles [0, 5) and [8, 12): in raw coordinates the second crosses the origin of a 10 bp reference.
        let g = graph(
            &[("A1", json!({"pos": 0})), ("A2", json!({"pos": 5})), ("A3", json!({"pos": 8})), ("A4", json!({"pos": 12}))],
            &[
                ("E1", json!({"src": ["A1"], "dst": ["A2"], "reads": ["r1", "r2"], "variants": "5="})),
                ("E2", json!({"src": ["A1"], "dst": ["A2"], "reads": ["r3"], "variants": "2=1X2="})),
                ("E3", json!({"src": ["A3"], "dst": ["A4"], "reads": ["r4"], "variants": "4="})),
                ("E4", json!({"src": ["A3"], "dst": ["A4"], "reads": ["r5"], "variants": "1=1X2="})),
            ],
        );
        let dict = get_graph_intervals(&g, 10);
        // Raw position 11 wraps to 1, which lies in the other bubble.
        let raw = vec![snp(11, "G", 1)];
        let wrapped = circuliarize_variants(raw.clone(), 10);
        assert_eq!(wrapped[0].pos, 1);

        let (mut read_record, mut cover_record) = (HashMap::new(), HashMap::new());
        record_variant_reads(&raw, &wrapped, &dict, &g, &[json!("r4"), json!("r5")], &mut read_record, &mut cover_record);
        assert_eq!(sorted(cover_record[&generate_variant_name(&wrapped[0])].clone()), ["r4", "r5"]);
    }

    #[test]
    fn bubble_cover_reads_excludes_edges_with_no_computed_variant() {
        // r3 crosses the bubble only via E_private, which has no CIGAR and so was never genotyped.
        let g = graph(
            &[("A1", json!({"pos": 0})), ("A2", json!({"pos": 10}))],
            &[
                ("E_real", json!({"src": ["A1"], "dst": ["A2"], "reads": ["r1", "r2"], "variants": "5=1X4="})),
                ("E_private", json!({"src": ["A1"], "dst": ["A2"], "reads": ["r3"]})),
            ],
        );
        let dict = get_graph_intervals(&g, 10);
        assert_eq!(sorted(bubble_cover_reads(&dict, 5, &g)), ["r1", "r2"]);
    }

    #[test]
    fn origin_spanning_bubble_covers_its_reads_on_both_sides_of_the_origin() {
        // E_wrap closes the circle: 8 -> 10 (== ref_length) -> 2.
        let g = graph(
            &[("A_first", json!({"pos": 2})), ("A_mid", json!({"pos": 5})), ("A_last", json!({"pos": 8}))],
            &[
                ("E_mid", json!({"src": ["A_first"], "dst": ["A_mid"], "reads": ["r1"], "variants": "3="})),
                ("E_wrap", json!({"src": ["A_last"], "dst": ["A_first"], "reads": ["rw1", "rw2"], "variants": "2=1X1="})),
            ],
        );
        let dict = get_graph_intervals(&g, 10);
        for raw_pos in [8, 10] {
            assert_eq!(sorted(bubble_cover_reads(&dict, raw_pos, &g)), ["rw1", "rw2"], "raw pos {raw_pos}");
        }
        assert_eq!(sorted(bubble_cover_reads(&dict, 3, &g)), ["r1"]);
    }

    #[test]
    fn source_and_sink_edges_get_no_interval() {
        // SOURCE/SINK have no anchor coordinate.
        let g = graph(
            &[("A1", json!({"pos": 900})), ("A2", json!({"pos": 910}))],
            &[
                ("E_real", json!({"src": ["A1"], "dst": ["A2"], "reads": ["r1"], "variants": "5=1X4="})),
                ("E_source", json!({"src": ["SOURCE"], "dst": ["A2"], "reads": ["rs"], "variants": "5="})),
                ("E_sink", json!({"src": ["A1"], "dst": ["SINK"], "reads": ["rk"], "variants": "5="})),
            ],
        );
        let dict = get_graph_intervals(&g, 16569);
        assert_eq!(dict.len(), 1);
        assert!(bubble_cover_reads(&dict, 42, &g).is_empty());
        assert_eq!(sorted(bubble_cover_reads(&dict, 905, &g)), ["r1"]);
    }

    // ---- permutation test ----

    #[test]
    fn permutation_test_rejects_a_scattered_snp_except_on_pacbio() {
        // A SNP and an indel carried by the same 3 of 30 reads, and a SNP scattered over 3 other reads.
        let variants = vec![snp(100, "A", 3), var(200, "TT", "T", 3), snp(300, "G", 3)];
        let mut matrix = Array2::<f64>::zeros((3, 30));
        matrix.slice_mut(s![0..2, 0..3]).fill(1.0);
        for j in [10, 17, 24] {
            matrix[[2, j]] = 1.0;
        }
        let coverage: HashMap<usize, usize> = variants.iter().map(|v| (v.pos, 30)).collect();
        let artifact = generate_variant_name(&variants[2]);
        // Denoising leaves systematic error behind, so ont-denoised tests SNPs; PacBio SNPs are exempt.
        for (data_type, artifact_kept) in [("ont-denoised", false), ("pacbio", true)] {
            let (kept, _) = permutation_test(&matrix, 0.5, 0.9, 200, &variants, &coverage, 0.9, data_type, resolve_indel_false_threshold(data_type, None));
            let names: Vec<String> = kept.iter().map(generate_variant_name).collect();
            assert_eq!(names.contains(&artifact), artifact_kept, "{data_type}");
            assert_eq!(names.len(), 2 + artifact_kept as usize, "{data_type}: the lineage pair survives");
        }
    }

    #[test]
    fn permutation_test_keeps_variants_it_cannot_test() {
        // Alt in 9 of 10 matrix reads: above the 0.7 permutation frequency threshold, so no null is drawn,
        // although allele_count / depth = 0.02 admits both variants.
        let variants = vec![snp(100, "A", 2), snp(200, "C", 2)];
        let mut matrix = Array2::<f64>::zeros((2, 10));
        matrix.slice_mut(s![.., ..9]).fill(1.0);
        let coverage: HashMap<usize, usize> = [(100, 100), (200, 100)].into_iter().collect();
        let (kept, _) = permutation_test(&matrix, 0.5, 0.5, 16, &variants, &coverage, 0.7, "ont-denoised", 0.015);
        assert_eq!(kept.len(), 2, "an untested variant must be kept, not dropped");
    }

    /// X sits on 45 random reads of the majority lineage A (reads 0..240), next to minor lineages B and C.
    fn run_majority_lineage(x: Variant) -> Vec<String> {
        let n_reads = 300;
        let lineages = [(0usize, 240usize, 3usize), (240, 276, 15), (276, 300, 15)];
        let n_markers: usize = lineages.iter().map(|l| l.2).sum();
        let mut matrix = Array2::<f64>::zeros((1 + n_markers, n_reads));
        let mut rng = StdRng::seed_from_u64(12969);
        let mut placed = 0;
        while placed < 45 {
            let j = rng.random_range(0..240);
            if matrix[[0, j]] == 0.0 {
                matrix[[0, j]] = 1.0;
                placed += 1;
            }
        }
        let mut variants = vec![x];
        for (lo, hi, n) in lineages {
            for _ in 0..n {
                matrix.slice_mut(s![variants.len(), lo..hi]).fill(1.0);
                variants.push(snp(1000 + 10 * (variants.len() - 1), "G", hi - lo));
            }
        }
        let coverage: HashMap<usize, usize> = variants.iter().map(|v| (v.pos, n_reads)).collect();
        let (kept, _) = permutation_test(&matrix, 0.01, 0.2, 1000, &variants, &coverage, 0.7, "ont-r10", 0.015);
        // A untested, B/C tested: all lineage markers survive.
        assert_eq!(kept.iter().filter(|v| v.pos >= 1000).count(), n_markers, "lineage markers must survive");
        kept.iter().map(generate_variant_name).collect()
    }

    #[test]
    fn cooccurrence_test_exempts_strong_balanced_snps() {
        // The co-occurrence statistic punishes a variant that sits on the MAJORITY lineage.
        let x = snp(500, "C", 45);
        let name = generate_variant_name(&x);
        assert!(
            run_majority_lineage(x.clone()).contains(&name),
            "a strong strand-balanced SNP must be exempt from the co-occurrence veto"
        );
        let biased = Variant { filter: Some("Strand_bias".into()), ..x };
        assert!(!run_majority_lineage(biased).contains(&name), "a strand-biased SNP must still be tested and rejected");
        let indel = var(500, "TT", "T", 45);
        let name = generate_variant_name(&indel);
        assert!(!run_majority_lineage(indel).contains(&name), "an indel with the same evidence must still be tested");
    }

    /// `n` random genotype cells: ~20% NaN (uncovered), ~10% alt.
    fn random_genotype_row(rng: &mut StdRng, n: usize) -> Vec<f64> {
        (0..n)
            .map(|_| match rng.random::<f64>() {
                u if u < 0.20 => f64::NAN,
                u if u < 0.30 => 1.0,
                _ => 0.0,
            })
            .collect()
    }

    #[test]
    fn bit_jaccard_matches_f64_jaccard() {
        // Reference: alt = > 0.5, reads uncovered in either row skipped.
        let jaccard_f64 = |a: &[f64], b: &[f64]| {
            let (mut inter, mut union) = (0, 0);
            for (x, y) in a.iter().zip(b).filter(|(x, y)| !x.is_nan() && !y.is_nan()) {
                inter += (*x > 0.5 && *y > 0.5) as u32;
                union += (*x > 0.5 || *y > 0.5) as u32;
            }
            if union == 0 { 0.0 } else { inter as f64 / union as f64 }
        };
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..200 {
            let a = random_genotype_row(&mut rng, 300);
            let b = random_genotype_row(&mut rng, 300);
            let (ba, bb) = (BitRow::from_row(&a), BitRow::from_row(&b));
            assert_eq!(jaccard_bits(&ba.alt, &ba.cov, &bb), jaccard_f64(&a, &b));
        }
    }

    #[test]
    fn bit_shuffle_preserves_alt_count_and_coverage() {
        let mut rng = StdRng::seed_from_u64(11);
        let br = BitRow::from_row(&random_genotype_row(&mut rng, 300));
        let mut scratch: Vec<u32> = Vec::new();
        let mut alt = vec![0u64; br.alt.len()];
        for _ in 0..1000 {
            shuffle_bits(&br, &mut scratch, &mut alt, &mut rng);
            assert_eq!(alt.iter().map(|x| x.count_ones()).sum::<u32>(), br.n_alt, "shuffle must preserve the alt count");
            assert!(alt.iter().zip(&br.cov).all(|(a, c)| a & !c == 0), "shuffled alts must stay inside the covered reads");
        }
    }

    #[test]
    fn adaptive_permutation_resolves_the_bh_floor() {
        // 10 identical lineage rows (alt in the same 24 of 400 reads) plus 120 independent 8-read noise rows.
        let n_reads = 400;
        let mut rng = StdRng::seed_from_u64(2026);
        let mut variants = Vec::new();
        let mut matrix = Array2::<f64>::zeros((130, n_reads));
        for i in 0..130 {
            variants.push(snp(100 + i, "A", if i < 10 { 24 } else { 8 }));
            if i < 10 {
                matrix.slice_mut(s![i, ..24]).fill(1.0);
            } else {
                let mut placed = 0;
                while placed < 8 {
                    let j = rng.random_range(0..n_reads);
                    if matrix[[i, j]] == 0.0 {
                        matrix[[i, j]] = 1.0;
                        placed += 1;
                    }
                }
            }
        }
        let coverage: HashMap<usize, usize> = variants.iter().map(|v| (v.pos, n_reads)).collect();

        // 1000 rounds alone put the lineage rows at p = 1/1001, BH q = 0.013 > 0.01: stage 2 must lift that floor.
        let (kept, _) = permutation_test(&matrix, 0.01, 0.2, 1000, &variants, &coverage, 0.7, "ont-r10", 0.015);
        let lineage = kept.iter().filter(|v| v.pos < 110).count();
        assert_eq!(lineage, 10, "every lineage row must be kept");
        assert!(kept.len() - lineage <= 1, "at most one noise row may survive, got {}", kept.len() - lineage);

        let tested: Vec<usize> = (0..130).collect();
        let names: Vec<String> = variants.iter().map(generate_variant_name).collect();
        let positions: Vec<usize> = variants.iter().map(|v| v.pos).collect();
        let partners: Vec<usize> = (0..130).collect(); // PartnerRule::All
        let outcomes = adaptive_permutation(&matrix, &tested, &names, &partners, &positions, 0, 1000, 0.01);
        let needed = (PERM_HEADROOM * 130.0 / (0.01 * 10.0)).ceil() as u32;
        let max_draws = outcomes.iter().map(|o| o.draws).max().unwrap();
        assert!(max_draws >= needed, "{needed} rounds needed, drew at most {max_draws}");
        // Clear nulls stop early.
        let clear: Vec<&PermOutcome> = outcomes.iter().filter(|o| o.p_value() > 0.1).collect();
        assert!(clear.len() > 50, "expected many clear nulls, saw {}", clear.len());
        assert!(clear.iter().all(|o| o.draws < 1000), "clear nulls must stop before 1000 draws");
        // Per-variant seeds make the outcome independent of thread scheduling.
        assert_eq!(
            adaptive_permutation(&matrix, &tested, &names, &partners, &positions, 0, 1000, 0.01),
            outcomes
        );
    }

    // ---- PacBio precision plan T1: depth-aware noise threshold ----

    #[test]
    fn pb_noise_threshold_matches_poisson_table() {
        assert_eq!(pb_noise_threshold(0), 3);
        assert_eq!(pb_noise_threshold(6), 3);
        assert_eq!(pb_noise_threshold(88), 3);
        assert_eq!(pb_noise_threshold(268), 5);
        assert_eq!(pb_noise_threshold(970), 8);
        assert_eq!(pb_noise_threshold(6494), 25);
    }

    // ---- PacBio precision plan T2: gray-zone admission ----

    #[test]
    fn pacbio_snps_are_tested_only_in_gray_zone() {
        // At DP 6500, T=25: AC 20 is gray-zone (admitted), AC 30 is exempt.
        assert!(20 < pb_noise_threshold(6500));
        assert!(!(30 < pb_noise_threshold(6500)));
        // Linked AC-20 SNPs share reads with an AC-30 backbone (the only partner under PbNoiseThreshold);
        // a scattered AC-20 does not. AC-30 itself is exempt (not tested).
        let variants = vec![
            snp(1000, "A", 20),
            snp(3000, "A", 20),
            snp(5000, "C", 30),
            snp(7000, "G", 20),
        ];
        let mut matrix = Array2::<f64>::zeros((4, 200));
        for j in 0..20 {
            matrix[[0, j]] = 1.0;
            matrix[[1, j]] = 1.0;
            matrix[[2, j]] = 1.0; // backbone overlaps the linked pair
        }
        for j in 20..30 {
            matrix[[2, j]] = 1.0; // AC 30
        }
        // Scattered AC-20: disjoint from the backbone.
        for j in 50..70 {
            matrix[[3, j]] = 1.0;
        }
        let coverage: HashMap<usize, usize> = variants.iter().map(|v| (v.pos, 6500)).collect();
        let (kept, _) = permutation_test(&matrix, 0.01, 0.2, 1000, &variants, &coverage, 0.7, "pacbio", 0.1);
        let pos: Vec<usize> = kept.iter().map(|v| v.pos).collect();
        assert!(pos.contains(&1000) && pos.contains(&3000), "linked gray-zone SNPs kept: {pos:?}");
        assert!(pos.contains(&5000), "AC≥T SNP is exempt and kept: {pos:?}");
        assert!(!pos.contains(&7000), "scattered gray-zone SNP is tested and dropped: {pos:?}");

        // At DP 90, T=3: no SNP with AC≥3 is admitted (all exempt).
        assert_eq!(pb_noise_threshold(90), 3);
        let variants90 = vec![snp(100, "A", 3), snp(200, "C", 5)];
        let mut matrix90 = Array2::<f64>::zeros((2, 30));
        for j in [10, 17, 24] {
            matrix90[[0, j]] = 1.0;
        }
        for j in 0..5 {
            matrix90[[1, j]] = 1.0;
        }
        let coverage90: HashMap<usize, usize> = [(100, 90), (200, 90)].into_iter().collect();
        let (kept90, _) = permutation_test(&matrix90, 0.01, 0.2, 200, &variants90, &coverage90, 0.7, "pacbio", 0.1);
        assert_eq!(kept90.len(), 2, "at DP 90 every AC≥3 SNP is exempt: {:?}", kept90.iter().map(|v| v.pos).collect::<Vec<_>>());
    }

    // ---- frozen ONT baseline (PacBio precision plan, T0) ----
    // Later tasks APPEND rows/assertions; the ONT rows and literals below never change.

    #[test]
    fn ont_presets_are_frozen() {
        for dt in ["ont-r9", "ont-r10", "ont-denoised"] {
            let want = if dt == "ont-r9" { (0.0001, 0.5, 0.7) } else { (0.01, 0.2, 0.7) };
            assert_eq!(resolve_thresholds(dt, None, None, None), want, "{dt}");
            assert_eq!(resolve_indel_false_threshold(dt, None), 0.015, "{dt}");
            assert!(is_ont(dt), "{dt}");
            assert_eq!(resolve_perm_opts(dt), PermOpts::all_rows(), "{dt}");
            assert_eq!(resolve_vaf(dt, None), 0.01, "{dt}");
        }
    }

    #[test]
    fn resolve_vaf_presets_and_explicit_override() {
        assert_eq!(resolve_vaf("pacbio", None), 0.002);
        assert_eq!(resolve_vaf("ont-r10", None), 0.01);
        assert_eq!(resolve_vaf("pacbio", Some(0.05)), 0.05);
        assert_eq!(resolve_vaf("ont-r10", Some(0.0)), 0.0);
    }

    #[test]
    fn pacbio_small_indel_below_false_threshold_is_not_admitted() {
        // 5% indel (HF < 0.1) is not admitted → kept even with no partners.
        // 12% indel is admitted and, with no partners, dropped.
        let variants = vec![
            var(1000, "TT", "T", 50),  // HF 0.05 at DP 1000
            var(3000, "TT", "T", 120), // HF 0.12
        ];
        let mut matrix = Array2::<f64>::zeros((2, 200));
        for j in 0..50 {
            matrix[[0, j]] = 1.0;
        }
        for j in 60..180 {
            matrix[[1, j]] = 1.0;
        }
        let coverage: HashMap<usize, usize> = [(1000, 1000), (3000, 1000)].into_iter().collect();
        let (kept, _) = permutation_test(&matrix, 0.01, 0.2, 500, &variants, &coverage, 0.7, "pacbio", 0.1);
        let pos: Vec<usize> = kept.iter().map(|v| v.pos).collect();
        assert!(pos.contains(&1000), "5% indel not admitted → kept: {pos:?}");
        assert!(!pos.contains(&3000), "12% indel admitted and fails alone: {pos:?}");
    }

    /// 200 reads, depth 6500 everywhere: straddles every PacBio rule (T = 25 there, all variants gray-zone).
    fn frozen_ont_fixture() -> (Vec<Variant>, Array2<f64>, HashMap<usize, usize>) {
        let mut variants = Vec::new();
        let mut rows: Vec<Vec<usize>> = Vec::new();
        let mut add = |v: Variant, reads: Vec<usize>| {
            variants.push(v);
            rows.push(reads);
        };
        for pos in [1000, 3000, 5000] {
            add(snp(pos, "A", 20), (0..20).collect()); // lineage A
        }
        for pos in [7000, 9000] {
            add(snp(pos, "C", 6), (20..26).collect()); // lineage B
        }
        add(snp(11000, "G", 3), (30..33).collect()); // adjacent error pair
        add(snp(11001, "G", 3), (30..33).collect());
        add(snp(13000, "A", 3), vec![40, 90, 150]); // scattered
        add(var(15000, "TT", "T", 10), (0..10).collect()); // 1-bp deletion, HF 0.05 by matrix
        add(snp(17000, "C", 30), (100..130).collect()); // AC 30
        let mut matrix = Array2::<f64>::zeros((variants.len(), 200));
        for (i, reads) in rows.iter().enumerate() {
            for &j in reads {
                matrix[[i, j]] = 1.0;
            }
        }
        let coverage = variants.iter().map(|v| (v.pos, 6500)).collect();
        (variants, matrix, coverage)
    }

    #[test]
    fn ont_permutation_is_frozen() {
        let (variants, matrix, coverage) = frozen_ont_fixture();
        let names: Vec<String> = variants.iter().map(generate_variant_name).collect();
        let tested: Vec<usize> = (0..variants.len()).collect();
        // (kept variant positions, (exceed, draws) per variant, in fixture order); literals captured at the starting commit.
        // Fixture order: 1000, 3000, 5000 | 7000, 9000 | 11000, 11001 | 13000 (scattered) | 15000 del | 17000 (AC 30).
        // The scattered SNP and the AC-30 SNP are indistinguishable from the null (100 exceedances in the first
        // 100 draws, p ~ 1) and are rejected; every other row sits at the p floor and is kept.
        let kept_all: Vec<usize> = vec![1000, 3000, 5000, 7000, 9000, 11000, 11001, 15000];
        let outcomes_all: Vec<(u32, u32)> = vec![
            (0, 1000), (0, 1000), (0, 1000), (0, 1000), (0, 1000), (0, 1000), (0, 1000), (100, 100), (0, 1000), (100, 100),
        ];
        let frozen: [(&str, Vec<usize>, Vec<(u32, u32)>); 3] = [
            ("ont-r9", kept_all.clone(), outcomes_all.clone()),
            ("ont-r10", kept_all.clone(), outcomes_all.clone()),
            ("ont-denoised", kept_all, outcomes_all),
        ];
        for (dt, want_kept, want_outcomes) in frozen {
            let (p, f, perm_freq) = (0.01, 0.2, 0.7);
            let (kept, _) = permutation_test(&matrix, p, f, 1000, &variants, &coverage, perm_freq, dt, resolve_indel_false_threshold(dt, None));
            let got_kept: Vec<usize> = kept.iter().map(|v| v.pos).collect();
            let positions: Vec<usize> = variants.iter().map(|v| v.pos).collect();
            let partners: Vec<usize> = (0..variants.len()).collect(); // ONT: PartnerRule::All
            let outcomes = adaptive_permutation(&matrix, &tested, &names, &partners, &positions, 0, 1000, p);
            let got: Vec<(u32, u32)> = outcomes.iter().map(|o| (o.exceed, o.draws)).collect();
            assert_eq!(got_kept, want_kept, "{dt}: kept");
            assert_eq!(got, want_outcomes, "{dt}: (exceed, draws)");
        }
        // T2/T3 discrimination: PacBio gray-zone + partner rules change the kept set vs ont-r10.
        // AC 30 (17000) is exempt (AC ≥ T); adjacent pair 11000/11001 has no partner ≥20 bp away → dropped.
        let (pb_kept, _) = permutation_test(&matrix, 0.01, 0.2, 1000, &variants, &coverage, 0.7, "pacbio", 0.1);
        let pb_pos: Vec<usize> = pb_kept.iter().map(|v| v.pos).collect();
        assert_ne!(pb_pos, vec![1000, 3000, 5000, 7000, 9000, 11000, 11001, 15000], "pacbio kept must differ from ont-r10");
        assert!(pb_pos.contains(&17000), "AC≥T SNP stays on PacBio: {pb_pos:?}");
        assert!(!pb_pos.contains(&11000) && !pb_pos.contains(&11001), "adjacent pair dropped on PacBio: {pb_pos:?}");
    }

    #[test]
    fn pacbio_adjacent_pair_has_null_p_without_distant_partner() {
        // pos 100/101, same reads, no other partner above T(DP): obs ≈ 0 after min-dist, p ≈ 1.
        let variants = vec![snp(100, "A", 5), snp(101, "C", 5), snp(5000, "G", 30)];
        let mut matrix = Array2::<f64>::zeros((3, 100));
        for j in 0..5 {
            matrix[[0, j]] = 1.0;
            matrix[[1, j]] = 1.0;
        }
        for j in 20..50 {
            matrix[[2, j]] = 1.0;
        }
        // Depth high enough that T > 5, so only the AC-30 row is a partner.
        let coverage: HashMap<usize, usize> = variants.iter().map(|v| (v.pos, 6500)).collect();
        assert_eq!(resolve_perm_opts("pacbio").partner_min_dist, 20);
        let (kept, _) = permutation_test(&matrix, 0.01, 0.2, 500, &variants, &coverage, 0.7, "pacbio", 0.1);
        let pos: Vec<usize> = kept.iter().map(|v| v.pos).collect();
        assert!(!pos.contains(&100) && !pos.contains(&101), "adjacent pair must be rejected: {pos:?}");
        assert!(pos.contains(&5000), "backbone partner itself is exempt: {pos:?}");
    }

    #[test]
    fn pacbio_noise_overlapping_only_subthreshold_partners_is_rejected() {
        // Gray-zone noise shares reads only with other sub-threshold SNPs → empty effective partner set after
        // PbNoiseThreshold, obs = 0, p ≈ 1.
        let variants = vec![
            snp(1000, "A", 6),
            snp(3000, "C", 6),
            snp(5000, "G", 6),
        ];
        let mut matrix = Array2::<f64>::zeros((3, 100));
        for j in 0..6 {
            matrix[[0, j]] = 1.0;
            matrix[[1, j]] = 1.0;
            matrix[[2, j]] = 1.0;
        }
        let coverage: HashMap<usize, usize> = variants.iter().map(|v| (v.pos, 6500)).collect();
        // T(6500)=25: all three are gray-zone (admitted) but none qualify as partners.
        let (kept, _) = permutation_test(&matrix, 0.01, 0.2, 500, &variants, &coverage, 0.7, "pacbio", 0.1);
        assert!(kept.is_empty(), "noise with only sub-threshold partners must all be dropped: {:?}", kept.iter().map(|v| v.pos).collect::<Vec<_>>());
    }

    // ---- strand bias ----

    /// One variant's strand evidence: alt reads per strand, plus further covering reads per strand.
    struct SbCase {
        variants: Vec<Variant>,
        matrix: Array2<f64>,
        read_record: HashMap<String, Vec<serde_json::Value>>,
        strand_map: HashMap<String, bool>,
        cover_record: HashMap<String, HashSet<String>>,
        coverage: HashMap<usize, usize>,
    }

    fn sb_case(variant: Variant, alt_fwd: usize, alt_rev: usize, cover_fwd: usize, cover_rev: usize, depth: usize) -> SbCase {
        let name = generate_variant_name(&variant);
        let (mut strand_map, mut alt, mut cover) = (HashMap::new(), Vec::new(), HashSet::new());
        for (prefix, n, is_rev, is_alt) in
            [("af", alt_fwd, false, true), ("ar", alt_rev, true, true), ("cf", cover_fwd, false, false), ("cr", cover_rev, true, false)]
        {
            for i in 0..n {
                let read = format!("{prefix}{i}");
                strand_map.insert(read.clone(), is_rev);
                if is_alt {
                    alt.push(json!(read));
                }
                // As in `get_variant`, the cover set includes the alt reads.
                cover.insert(read);
            }
        }
        SbCase {
            matrix: Array2::zeros((1, 1)),
            read_record: [(name.clone(), alt)].into_iter().collect(),
            cover_record: [(name, cover)].into_iter().collect(),
            coverage: [(variant.pos, depth)].into_iter().collect(),
            variants: vec![variant],
            strand_map,
        }
    }

    fn run_sb(case: &SbCase, library_frac: f64, action: StrandBiasAction) -> (Vec<Variant>, Array2<f64>) {
        filter_strand_bias(
            &case.variants, &case.matrix, &case.read_record, &case.strand_map, &case.cover_record, &case.coverage,
            library_frac, 0.01, 0.7, action,
        )
    }

    #[test]
    fn strand_bias_drop_removes_only_low_hf_snps() {
        // 10 forward-only alt reads against a balanced library.
        let flagged = |kept: &[Variant]| kept.len() == 1 && kept[0].filter.as_deref() == Some("Strand_bias");
        let low = sb_case(snp(100, "A", 10), 10, 0, 0, 0, 200); // HF 0.05
        assert!(flagged(&run_sb(&low, 0.5, StrandBiasAction::Flag).0));
        let (kept, kept_matrix) = run_sb(&low, 0.5, StrandBiasAction::Drop);
        assert!(kept.is_empty() && kept_matrix.nrows() == 0, "dropped together with its matrix row");
        // At HF 0.20 (above STRAND_BIAS_DROP_MAX_HF), or as an indel: flagged, never dropped.
        let high = sb_case(snp(100, "A", 10), 10, 0, 0, 0, 50);
        assert!(flagged(&run_sb(&high, 0.5, StrandBiasAction::Drop).0));
        let indel = sb_case(var(100, "TT", "T", 10), 10, 0, 0, 0, 200);
        assert!(flagged(&run_sb(&indel, 0.5, StrandBiasAction::Drop).0));
    }

    #[test]
    fn strand_bias_tests_against_the_site_local_null() {
        assert!(strand_bias_pvalue(600, 400, 0.5) < 0.01);
        // (case, library forward fraction, flagged)
        let cases = [
            // 7805G>A: alt 45F/27R is skewed against the library but not against its 54%-forward site.
            (sb_case(snp(7805, "A", 72), 45, 27, 54, 46, 1000), 0.498, false),
            // 600F/400R is significant, but its 0.10 effect is below STRAND_BIAS_MIN_EFFECT.
            (sb_case(snp(500, "A", 1000), 600, 400, 50, 50, 5000), 0.5, false),
            // 33CG>C: the 34 reverse-only alt reads are left out of the null, or they would hide their own skew.
            (sb_case(var(33, "CG", "C", 34), 0, 34, 2, 2, 302), 0.498, true),
            // Only 10 covering reads without the alt: too few for a local null, so the library fraction is used.
            (sb_case(snp(100, "A", 10), 10, 0, 10, 0, 200), 0.5, true),
        ];
        for (i, (case, library_frac, flagged)) in cases.iter().enumerate() {
            let kept = run_sb(case, *library_frac, StrandBiasAction::Flag).0;
            assert_eq!(kept.len(), 1, "case {i}");
            assert_eq!(kept[0].filter.is_some(), *flagged, "case {i}");
        }
    }

    #[test]
    fn bam_recount_retests_only_strand_flagged_indels() {
        let sb = Some("Strand_bias");
        // (variant, FILTER before, BAM alt fwd/rev at a balanced site, FILTER after)
        let cases = [
            (var(302, "A", "ACC", 20), sb, 9, 11, None), // balanced in the BAM: cleared
            (var(33, "CG", "C", 20), sb, 1, 19, sb),     // still skewed
            (var(302, "A", "ACC", 20), sb, 1, 0, sb),    // too few BAM alt reads to retest
            (snp(500, "A", 20), sb, 9, 11, sb),          // SNVs are not recounted
            (var(302, "A", "ACC", 20), Some("Potential_Artifact"), 9, 11, Some("Potential_Artifact")),
        ];
        for (v, before, fwd, rev, after) in cases {
            let v = Variant { filter: before.map(String::from), ..v };
            let bam = sb_case(v.clone(), fwd, rev, 50, 50, 1000);
            let out = recount_indel_strand_bias(&[v], &bam.read_record, &bam.cover_record, &bam.strand_map, 0.5, 0.01);
            assert_eq!(out[0].filter.as_deref(), after, "{} {fwd}/{rev}", generate_variant_name(&out[0]));
        }
    }

    // ---- left alignment ----

    /// `TCAAAAGCTT` - a 4bp A-run at indices 2..6 preceded by a C.
    const HOMOPOLYMER_REF: &[u8] = b"TCAAAAGCTT";

    /// Build a DEL the way `get_variants_from_cigar` does.
    fn del_at(refseq: &[u8], pos: usize, len: usize, allele_count: usize) -> Variant {
        let s = |r: std::ops::Range<usize>| std::str::from_utf8(&refseq[r]).unwrap().to_string();
        var(pos, &s(pos - 1..pos + len), &s(pos - 1..pos), allele_count)
    }

    /// Build an INS the way `get_variants_from_cigar` does.
    fn ins_at(refseq: &[u8], pos: usize, inserted: &str, allele_count: usize) -> Variant {
        let anchor = std::str::from_utf8(&refseq[pos - 1..pos]).unwrap();
        var(pos, anchor, &format!("{anchor}{inserted}"), allele_count)
    }

    #[test]
    fn left_align_variant_shifts_indels_to_the_leftmost_anchor() {
        let hp = HOMOPOLYMER_REF;
        // rCRS 560..=580, covering the poly-C run at 568-573: slice index `i` is 1-based coordinate 560 + i.
        let rcrs: &[u8] = b"CAAAGACACCCCCCACAGTTT";
        let anchorless = Variant { variant_type: "INS".into(), ..var(3, "-", "-", 2) };
        // (reference, input, expected (pos, REF, ALT))
        let cases = [
            // The deleted / inserted A moves to the base before the A-run, as `bcftools norm` reports it.
            (hp, del_at(hp, 5, 1, 3), (2, "CA", "C")),
            (hp, ins_at(hp, 6, "A", 3), (2, "C", "CA")),
            // A multi-base insertion rotates across an AG repeat.
            (b"TTAGAGAGC" as &[u8], ins_at(b"TTAGAGAGC", 8, "AG", 2), (2, "T", "TAG")),
            // A run reaching the reference start stops at the first possible anchor.
            (b"AAAAC", del_at(b"AAAAC", 3, 1, 2), (1, "AA", "A")),
            // Unique sequence, SNPs, anchorless and origin indels stay put.
            (b"ACGTTGCA", del_at(b"ACGTTGCA", 5, 1, 4), (5, "TG", "T")),
            (hp, var(4, "A", "G", 7), (4, "A", "G")),
            (hp, anchorless, (3, "-", "-")),
            (hp, var(0, "TC", "T", 5), (0, "TC", "T")),
            // rCRS 572 CC>C and 573 C>CC both normalize to 567 (cross-checked with `bcftools norm`).
            (rcrs, del_at(rcrs, 13, 1, 4), (567 - 559, "AC", "A")),
            (rcrs, ins_at(rcrs, 14, "C", 4), (567 - 559, "A", "AC")),
        ];
        for (refseq, v, want) in cases {
            let got = left_align_variant(&v, refseq);
            assert_eq!((got.pos, got.ref_allele.as_str(), got.alt_allele.as_str()), want, "{v:?}");
            assert_eq!((got.variant_type, got.allele_count), (v.variant_type, v.allele_count), "evidence is untouched");
        }
    }

    /// One bubble with two alt edges spelling the same 1 bp A-run deletion at different anchors.
    fn homopolymer_bubble_graph() -> GraphicalGenome {
        graph(
            &[("A1", json!({"pos": 0, "seq": "TCA"})), ("A2", json!({"pos": 7, "seq": "CTT"}))],
            &[
                ("E_ref", json!({"src": ["A1"], "dst": ["A2"], "seq": "AAAG", "reads": ["chrM_ref"], "variants": "7="})),
                ("E_alt1", json!({"src": ["A1"], "dst": ["A2"], "seq": "AAG", "reads": ["r1", "r2", "r3"], "variants": "5=1D1="})),
                ("E_alt2", json!({"src": ["A1"], "dst": ["A2"], "seq": "AAG", "reads": ["r4", "r5", "r6", "r7"], "variants": "3=1D3="})),
            ],
        )
    }

    #[test]
    fn get_variant_left_align_merges_equivalent_spellings_and_their_reads() {
        let n = HOMOPOLYMER_REF.len();
        // Without alignment the two spellings stay separate records at their original anchors.
        let (variants, _, read_record, _) = get_variant(&mut homopolymer_bubble_graph(), 3, "chrM_ref", n, None);
        let mut names: Vec<String> = variants.iter().map(generate_variant_name).collect();
        names.sort();
        assert_eq!(names, ["m.3AA>A", "m.5AA>A"]);
        assert!(read_record.keys().all(|k| names.contains(k)));

        // Aligned, both records are keyed by the one name `construct_matrix` looks rows up by.
        let (variants, _, read_record, cover_record) =
            get_variant(&mut homopolymer_bubble_graph(), 3, "chrM_ref", n, Some(HOMOPOLYMER_REF));
        assert!(variants.iter().all(|v| generate_variant_name(v) == "m.2CA>C"));
        assert!(read_record.keys().chain(cover_record.keys()).all(|k| k == "m.2CA>C"));
        let mut alt: Vec<&str> = read_record["m.2CA>C"].iter().map(|r| r.as_str().unwrap()).collect();
        alt.sort();
        assert_eq!(alt, ["r1", "r2", "r3", "r4", "r5", "r6", "r7"]);
        let collapsed = collapse_identical_records(variants, n);
        assert_eq!((collapsed.len(), collapsed[0].allele_count), (1, 7), "evidence from both spellings sums");
    }

    #[test]
    fn vcf_header_declares_left_alignment_only_when_on() {
        for left_align in [true, false] {
            let path = std::env::temp_dir().join(format!("himito_left_align_{left_align}_{}.vcf", std::process::id()));
            let coverage: HashMap<usize, usize> = [(5, 30)].into_iter().collect();
            write_vcf(&[del_at(HOMOPOLYMER_REF, 5, 1, 3)], &coverage, &path, "s1", "chrM", 10, 0.1, left_align, None).unwrap();
            let vcf = std::fs::read_to_string(&path).unwrap();
            assert_eq!(vcf.contains("##HimitoNormalization=left-aligned"), left_align);
        }
    }

    #[test]
    fn small_indel_at_ont_threshold_passes_above_and_is_artifact_below() {
        let cov: HashMap<usize, usize> = [(100, 1000)].into_iter().collect();
        let rec = format_vcf_record(&var(100, "GA", "G", 18), cov.clone(), 0.015);
        assert!(rec.contains("\tPASS\t"), "HF 0.018 should pass at 0.015: {rec}");
        let rec = format_vcf_record(&var(100, "GA", "G", 11), cov, 0.015);
        assert!(rec.contains("Potential_Artifact"), "HF 0.011 should be flagged: {rec}");
    }

    #[test]
    fn indel_false_threshold_preset_by_data_type() {
        for (dt, preset) in [("ont-r9", 0.015), ("ont-r10", 0.015), ("ont-denoised", 0.015), ("pacbio", 0.1)] {
            assert_eq!(resolve_indel_false_threshold(dt, None), preset, "{dt}");
            assert_eq!(resolve_indel_false_threshold(dt, Some(0.05)), 0.05, "{dt}: an explicit value wins");
        }
    }

    // ---- BAM genotyping ----

    #[test]
    fn allele_span_starts_at_the_first_ref_base() {
        for (v, want) in [
            (snp(100, "A", 5), (100, "T", "A")),
            (del_at(HOMOPOLYMER_REF, 5, 1, 3), (4, "AA", "A")),
            (ins_at(HOMOPOLYMER_REF, 6, "A", 3), (5, "A", "AA")),
        ] {
            assert_eq!(allele_span(&v), want);
            // The variant name carries the same position, 1-based.
            assert!(generate_variant_name(&v).starts_with(&format!("m.{}", want.0 + 1)));
        }
    }

    #[test]
    fn genotype_window_flanks_pushes_out_of_homopolymers_and_pads() {
        let unique: Vec<u8> = b"ACGT".iter().cycle().take(40).copied().collect();
        // SNP at 20: [20 - 10, 21 + 10), then 2 more on each side.
        assert_eq!(genotype_window(&unique, (20, 21)), Some((8, 33)));
        // Leaving the reference on either side: no window.
        assert_eq!(genotype_window(&unique, (3, 4)), None);
        assert_eq!(genotype_window(&unique, (36, 37)), None);
        // Window ends inside A-runs are pushed out of them before padding.
        let mut hp = unique.clone();
        hp[8..12].fill(b'A');
        hp[30..34].fill(b'A');
        assert_eq!(genotype_window(&hp, (20, 21)), Some((6, 36)));
    }

    #[test]
    fn group_sites_merges_overlapping_ref_spans_but_never_a_snp_with_an_indel() {
        for (vs, want) in [
            // SNPs at the same base; m.302A>AC and m.302A>ACC, both anchored on base 302.
            (vec![snp(100, "A", 5), snp(100, "C", 4)], vec![vec![0, 1]]),
            (vec![var(302, "A", "AC", 3), var(302, "A", "ACC", 3)], vec![vec![0, 1]]),
            // Adjacent SNPs, and DELs whose REF spans [9, 11) and [11, 13) merely touch.
            (vec![snp(100, "A", 5), snp(101, "A", 5)], vec![vec![0], vec![1]]),
            (vec![var(10, "CA", "C", 3), var(12, "GT", "G", 3)], vec![vec![0], vec![1]]),
            // ONT: m.8170AT>A is the deletion error of the very T in m.8171T>A.
            (vec![var(8170, "AT", "A", 9), snp(8170, "A", 7)], vec![vec![1], vec![0]]),
        ] {
            assert_eq!(group_sites(&vs), want, "{vs:?}");
        }
    }

    #[test]
    fn score_segment_calls_the_best_haplotypes_and_gates_noise() {
        const REF: &[u8] = b"GATCCTAGGCATTGCAAGTCCGATAC";
        let mut alt = REF.to_vec();
        alt[12] = b'A'; // m.13T>A
        let haps = vec![REF.to_vec(), alt.clone()];
        let mut del = REF.to_vec();
        del.remove(12); // the SNP base is deleted from the read
        let mut noisy = REF.to_vec(); // 6 substitutions away from REF
        for i in [0, 3, 6, 9, 15, 18] {
            noisy[i] = if noisy[i] == b'A' { b'C' } else { b'A' };
        }
        // Junk read bases around the window must not count against any haplotype.
        let seg = |core: &[u8]| [&b"GGAC"[..], core, &b"CATT"[..]].concat();
        let best = |s: &[usize]| SegCall::Best(s.iter().copied().collect());
        assert_eq!(score_segment(&haps, &seg(REF), u32::MAX), best(&[0]));
        assert_eq!(score_segment(&haps, &seg(&alt), u32::MAX), best(&[1]));
        // A tie comes back as both indices: missing downstream, not REF.
        assert_eq!(score_segment(&haps, &seg(&del), u32::MAX), best(&[0, 1]));
        assert_eq!(score_segment(&haps, &seg(&noisy), 4), SegCall::TooNoisy);
        assert!(matches!(score_segment(&haps, &seg(&noisy), u32::MAX), SegCall::Best(_)));
    }

    const GENO_REF: &str = "GCTAGACATACATACATACACGTCAGCACGACTGTGCAGTGTGATCGCTAGTAGTAGTGT";

    /// GENO_REF[30] == 'A': m.31A>T.
    fn geno_snp() -> Variant {
        var(30, "A", "T", 5)
    }

    fn geno_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_data/genotype").join(name)
    }

    fn run_geno(
        sam: &str,
        variants: &[Variant],
    ) -> (HashMap<String, Vec<serde_json::Value>>, HashMap<String, HashSet<String>>) {
        genotype_from_bam(&geno_path(sam), GENO_REF.as_bytes(), "chrM", variants, u32::MAX)
    }

    #[test]
    fn genotype_from_bam_calls_ref_alt_and_leaves_ambiguous_reads_missing() {
        let name = "m.31A>T";
        let (rr, cr) = run_geno("reads.sam", &[geno_snp()]);
        let alt = alt_reads(&rr, name);
        let cover = &cr[name];

        assert!(cover.contains("r_ref") && !alt.contains("r_ref"), "r_ref is a decisive 0");
        assert!(cover.contains("r_alt") && alt.contains("r_alt"), "r_alt is a decisive 1");
        assert!(!cover.contains("r_del"), "SNP base deleted: REF/ALT tie -> NaN, not 0");
        assert!(!cover.contains("r_short"), "read ends inside the window -> not covered");
        // Primary + supplementary agree (both ALT): one call.
        assert!(cover.contains("r_split") && alt.contains("r_split"));
        // Reference pseudo-read: covered, never alt (an all-0 column).
        assert!(cover.contains("chrM") && !alt.contains("chrM"));
    }

    #[test]
    fn genotype_from_bam_disagreeing_alignments_of_one_read_are_missing() {
        let name = "m.31A>T";
        let (rr, cr) = run_geno("split_disagree.sam", &[geno_snp()]);
        assert!(!cr[name].contains("r_split"));
        assert!(!alt_reads(&rr, name).contains("r_split"));
        assert!(cr[name].contains("chrM"));
    }

    #[test]
    fn genotype_from_bam_scores_an_indel_site_and_feeds_construct_matrix() {
        let del = del_at(GENO_REF.as_bytes(), 30, 1, 4); // r_del carries it
        let (sname, dname) = ("m.31A>T", generate_variant_name(&del));
        let variants = vec![geno_snp(), del];
        let (rr, cr) = run_geno("reads.sam", &variants);
        assert!(alt_reads(&rr, &dname).contains("r_del"));
        assert!(cr[&dname].contains("r_ref") && !alt_reads(&rr, &dname).contains("r_ref"));

        let (m, vars, reads) = construct_matrix(&rr, &cr, &variants);
        let col = |n: &str| reads.iter().position(|r| r == n).unwrap();
        let row = |n: &str| vars.iter().position(|v| generate_variant_name(v) == n).unwrap();
        assert_eq!(m[[row(sname), col("r_alt")]], 1.0);
        assert_eq!(m[[row(sname), col("r_ref")]], 0.0);
        assert!(m[[row(sname), col("r_del")]].is_nan());
        assert_eq!(m[[row(&dname), col("r_del")]], 1.0);
        assert_eq!(m[[row(sname), col("chrM")]], 0.0);
        assert_eq!(m[[row(&dname), col("chrM")]], 0.0);
    }

    #[test]
    fn genotype_from_bam_variant_without_a_window_is_only_the_ref_column() {
        // The window would cross the origin: no read can call it.
        let v = var(3, "A", "T", 5);
        let name = generate_variant_name(&v);
        let (rr, cr) = run_geno("reads.sam", &[v]);
        assert!(alt_reads(&rr, &name).is_empty());
        assert_eq!(sorted(cr[&name].clone()), ["chrM"]);
    }

    #[test]
    fn genotype_from_bam_matches_ref_name_when_chrm_is_not_tid_zero() {
        let name = "m.31A>T";
        let (rr, cr) = run_geno("reads_second_contig.sam", &[geno_snp()]);
        let (rr0, cr0) = run_geno("reads.sam", &[geno_snp()]);
        assert!(alt_reads(&rr, name).contains("r_alt"));
        assert_eq!(alt_reads(&rr, name), alt_reads(&rr0, name));
        assert_eq!(cr[name], cr0[name]);
    }

    #[test]
    fn genotype_from_bam_unknown_ref_name_returns_ref_column_only() {
        let name = "m.31A>T";
        let (rr, cr) = genotype_from_bam(
            &geno_path("reads_second_contig.sam"), GENO_REF.as_bytes(), "not_a_contig", &[geno_snp()], u32::MAX,
        );
        assert!(alt_reads(&rr, name).is_empty());
        assert_eq!(sorted(cr[name].clone()), ["not_a_contig"]);
    }

    #[test]
    fn write_final_matrix_regenotypes_from_the_bam_when_given() {
        let dir = std::env::temp_dir().join(format!("himito_final_matrix_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("out.vcf");
        let graph_matrix = Array2::from_shape_vec((1, 2), vec![1.0, f64::NAN]).unwrap();
        let graph_reads = ["g1".to_string(), "g2".to_string()];
        let matrix_csv = || std::fs::read_to_string(out.with_extension("matrix.csv")).unwrap();

        write_final_matrix(&out, &graph_matrix, &[geno_snp()], &[], &graph_reads, None).unwrap();
        assert_eq!(matrix_csv(), "variant,g1,g2\nm.31A>T,1,\n");

        let bam = geno_path("reads.sam");
        let src = GenotypeSource { bam: &bam, refseq: GENO_REF.as_bytes(), ref_name: "chrM", max_ed: u32::MAX };
        write_final_matrix(&out, &graph_matrix, &[geno_snp()], &[], &graph_reads, Some(&src)).unwrap();
        // Reads sorted by name, the ref pseudo-read all 0, r_del (tie) and r_short (partial) missing.
        assert_eq!(matrix_csv(), "variant,chrM,r_alt,r_ref,r_split\nm.31A>T,0,1,0,1\n");
    }
}
