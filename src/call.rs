use crate::agg;
use agg::*;
use rayon::iter::WhileSome;
use crate::build;
use indicatif::ProgressBar;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};
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

/// Populates `read_record`/`cover_record` for one edge's variants.
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

/// `left_align` carries the reference sequence when `--left-align` is on, and is `None` otherwise.
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
            Some(refseq) => left_align_variants(&variants_circular, refseq),
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

/// Shift every indel maximally leftward so equivalent alignments get one canonical representation.
fn left_align_variants(variants: &[Variant], refseq: &[u8]) -> Vec<Variant> {
    variants
        .iter()
        .map(|v| left_align_variant(v, refseq))
        .collect()
}

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
        "INS" => match variant.alt_allele.as_bytes().get(1..) {
            Some(inserted) => Allele::Ins(inserted.to_vec()),
            None => return variant.clone(),
        },
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

    Variant {
        pos: anchor + 1,
        ref_allele,
        alt_allele,
        variant_type: variant.variant_type.clone(),
        allele_count: variant.allele_count,
        filter: variant.filter.clone(),
    }
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

    // Every covered read is kept.
    (matrix, var_vec, read_vec)
}

// ---- BAM genotyping (Plan 2026-10-11) ----

/// `(start0, ref, alt)` where `start0` is the 0-based index of REF's first base on the reference.
fn allele_span(v: &Variant) -> (usize, &str, &str) {
    let start0 = if v.variant_type == "SNP" { v.pos } else { v.pos - 1 };
    (start0, v.ref_allele.as_str(), v.alt_allele.as_str())
}

/// Mirrors `site_window` in the oracle.
fn genotype_window(
    refseq: &[u8],
    sites_span: (usize, usize),
    flank: usize,
) -> Option<(usize, usize)> {
    let len = refseq.len() as isize;
    let mut left = sites_span.0 as isize - flank as isize;
    let mut right = sites_span.1 as isize + flank as isize;
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

/// Mirrors `group_sites` in the oracle: indices into `variants`, grouped into genotyping sites.
fn group_sites(variants: &[Variant]) -> Vec<Vec<usize>> {
    let mut sites: Vec<Vec<usize>> = Vec::new();
    for want_indel in [false, true] {
        let mut idx: Vec<usize> = (0..variants.len())
            .filter(|&i| {
                let (_, r, a) = allele_span(&variants[i]);
                (r.len() != a.len()) == want_indel
            })
            .collect();
        idx.sort_by_key(|&i| allele_span(&variants[i]).0); // stable, like Python
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
    /// Indices into the site's haplotype list within `margin` of the minimum.
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
fn score_segment(haps: &[(String, Vec<u8>)], seg: &[u8], margin: u32, max_ed: u32) -> SegCall {
    if haps.is_empty() {
        return SegCall::TooNoisy;
    }
    let dists: Vec<u32> = EDIT_ALIGNER.with(|al| {
        let mut al = al.borrow_mut();
        haps.iter()
            .map(|(_, h)| (-al.semiglobal(h, seg).score).max(0) as u32)
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
            .filter(|(_, &d)| d - best < margin)
            .map(|(i, _)| i)
            .collect(),
    )
}

/// Knobs of BAM genotyping (oracle: `--flank/--slack/--margin/--max-ed`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GenotypeOpts {
    /// Bases added on each side of a site before the homopolymer push-out.
    pub flank: usize,
    /// Extra read bases kept on each side of the aligned window.
    pub slack: usize,
    /// A haplotype is "best" when within `margin` edits of the minimum.
    pub margin: u32,
    /// Noise gate: a segment whose best distance exceeds this is unusable.
    pub max_ed: u32,
}

impl Default for GenotypeOpts {
    fn default() -> Self {
        GenotypeOpts { flank: 10, slack: 4, margin: 1, max_ed: u32::MAX }
    }
}

/// Matrix genotyping knobs by data type: ONT gets the noise gate, PacBio none.
pub fn resolve_genotype_opts(data_type: &str) -> GenotypeOpts {
    let max_ed = match data_type {
        "ont-r9" | "ont-r10" | "ont-denoised" => crate::GENOTYPE_MAX_ED_ONT,
        _ => u32::MAX,
    };
    GenotypeOpts { flank: crate::GENOTYPE_FLANK, max_ed, ..GenotypeOpts::default() }
}

/// Linkage rescue runs on ONT data unless switched off; PacBio is not evaluated.
pub fn linkage_rescue_enabled(data_type: &str, no_linkage_rescue: bool) -> bool {
    matches!(data_type, "ont-r9" | "ont-r10" | "ont-denoised") && !no_linkage_rescue
}

/// One genotyping site: a window on the reference and the haplotypes scored inside it.
struct GenoSite {
    left: usize,
    right: usize,
    members: Vec<usize>,
    haps: Vec<(String, Vec<u8>)>,
}

/// Per (site, read) evidence merged over all of the read's alignments.
#[derive(Default)]
struct SiteAcc {
    any_usable: bool,
    /// Hap indices that were the unique best of some alignment.
    winners: BTreeSet<usize>,
    /// Union of the best sets of all usable alignments.
    any_best: BTreeSet<usize>,
}

fn build_geno_sites(refseq: &[u8], variants: &[Variant], opts: &GenotypeOpts) -> Vec<GenoSite> {
    let mut sites = Vec::new();
    for members in group_sites(variants) {
        let spans: Vec<(usize, &str, &str)> =
            members.iter().map(|&i| allele_span(&variants[i])).collect();
        let lo = spans.iter().map(|s| s.0).min().unwrap();
        let hi = spans.iter().map(|s| s.0 + s.1.len()).max().unwrap();
        let (left, right) = match genotype_window(refseq, (lo, hi), opts.flank) {
            Some(w) => w,
            None => continue, // crosses the circular origin: left NaN (plan decision 5)
        };
        let mut haps: Vec<(String, Vec<u8>)> = vec![("REF".to_string(), refseq[left..right].to_vec())];
        for (&i, &(s, r, a)) in members.iter().zip(&spans) {
            let mut h = refseq[left..s].to_vec();
            h.extend_from_slice(a.as_bytes());
            h.extend_from_slice(&refseq[s + r.len()..right]);
            haps.push((generate_variant_name(&variants[i]), h));
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
    opts: &GenotypeOpts,
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
        let qa = q_lo.saturating_sub(opts.slack);
        let qb = (q_hi + 1 + opts.slack).min(seq.len());
        out.push((k, score_segment(&site.haps, &seq[qa..qb], opts.margin, opts.max_ed)));
    }
    out
}

/// Re-genotype `variants` in every read of `bam` by local realignment against the REF.
pub fn genotype_from_bam(
    bam: &Path,
    refseq: &[u8],
    ref_name: &str,
    variants: &[Variant],
    opts: &GenotypeOpts,
) -> (HashMap<String, Vec<serde_json::Value>>, HashMap<String, HashSet<String>>) {
    const BATCH: usize = 4096;
    let sites = build_geno_sites(refseq, variants, opts);
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
                    score_record_at_sites(rec, &sites, &lefts, opts),
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
                    e.any_usable = true;
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

    // Merge per (site, qname): decisive calls that disagree -> NaN; with none decisive.
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
        if !a.any_usable {
            continue;
        }
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

/// Sum of Jaccard coefficients between one genotype vector and every OTHER variant's row.
fn sum_jaccard_against_others(vector_data: &[f64], self_index: usize, rows: &[Vec<f64>]) -> f64 {
    rows.iter()
        .enumerate()
        .filter(|(j, _)| *j != self_index)
        .map(|(_, other)| genotype_jaccard_similarity(vector_data, other))
        .sum()
}

/// Generate a PER-VARIANT null distribution through permutation testing.
fn get_null_distribution(
    filtered_var: &Vec<Variant>,
    matrix: &Array2<f64>,
    permutation_round: usize,
    threshold: f64
) -> HashMap<String, Vec<f64>> {

    let bar = ProgressBar::new(filtered_var.len() as u64);
    // Materialize the rows ONCE.
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

/// This replaces a z-score against a normal CDF.
/// Callers must treat that as "untestable" rather than feeding it to BH, where a p-value of 1.0 reads
/// as a confirmed artifact.
pub(crate) fn empirical_p_value(null: &[f64], observation: f64) -> f64 {
    if null.is_empty() {
        return 1.0;
    }
    let at_or_above = null.iter().filter(|&&x| x >= observation).count();
    (at_or_above + 1) as f64 / (null.len() + 1) as f64
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

/// What `filter_strand_bias` does with a variant it finds strand-skewed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StrandBiasAction {
    Flag,
    /// Additionally REMOVE strand-skewed low-VAF SNPs outright.
    Drop,
}

/// Resolve the strand-bias action for a data type, honouring an explicit CLI choice.
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

/// Variants with fewer than `min_reads` strand-resolved supporting reads are kept.
fn filter_strand_bias(
    variants: &[Variant],
    matrix: &Array2<f64>,
    read_record: &HashMap<String, Vec<serde_json::Value>>,
    strand_map: &HashMap<String, bool>,
    cover_record: &HashMap<String, HashSet<String>>,
    coverage: &HashMap<usize, usize>,
    expected_fwd_frac: f64,
    p_threshold: f64,
    min_reads: usize,
    homoplasmic_frequency_threshold: f64,
    action: StrandBiasAction,
    drop_max_hf: f64,
    min_effect: f64,
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

        let (mut cfwd, mut crev) = (0usize, 0usize);
        if let Some(cover) = cover_record.get(&name) {
            for read_name in cover {
                // The null describes the reads that did NOT carry the alt.
                if seen.contains(read_name) {
                    continue;
                }
                match strand_map.get(read_name) {
                    Some(true) => crev += 1,
                    Some(false) => cfwd += 1,
                    None => {}
                }
            }
        }
        let null_frac = if cfwd + crev >= crate::STRAND_LOCAL_MIN_READS {
            cfwd as f64 / (cfwd + crev) as f64
        } else {
            expected_fwd_frac
        };
        let alt_frac = fwd as f64 / total as f64;
        let effect = (alt_frac - null_frac).abs();
        let p_value = strand_bias_pvalue(fwd, rev, null_frac);
        if p_value < p_threshold && effect >= min_effect {
            // The drop is deliberately narrow: SNPs only, and only below `drop_max_hf`.
            if action == StrandBiasAction::Drop
                && variant.variant_type == "SNP"
                && hf < drop_max_hf
            {
                println!(
                    "Strand-biased (dropped), {:?}, fwd={}, rev={}, alt_frac={:.3}, null_frac={:.3}, hf={:.4}, p={:.3e}",
                    name, fwd, rev, alt_frac, null_frac, hf, p_value
                );
                continue;
            }
            // flag but keep: downstream FILTER combines this with other filters instead of dropping the call
            variant.filter = Some("Strand_bias".to_string());
            println!(
                "Strand-biased (flagged, kept), {:?}, fwd={}, rev={}, alt_frac={:.3}, null_frac={:.3}, p={:.3e}",
                name, fwd, rev, alt_frac, null_frac, p_value
            );
        }
        keep_idx.push(i);
        kept.push(variant);
    }

    let filtered_matrix = matrix.select(Axis(0), &keep_idx);
    (kept, filtered_matrix)
}

/// Retest `Strand_bias`-flagged indels with alt/cover reads from BAM genotyping.
#[allow(clippy::too_many_arguments)]
fn recount_indel_strand_bias(
    variants: &[Variant],
    graph_counts: &HashMap<String, (usize, usize)>,
    bam_read_record: &HashMap<String, Vec<serde_json::Value>>,
    bam_cover_record: &HashMap<String, HashSet<String>>,
    strand_map: &HashMap<String, bool>,
    expected_fwd_frac: f64,
    p_threshold: f64,
    min_reads: usize,
    min_effect: f64,
) -> Vec<Variant> {
    let mut out = variants.to_vec();
    for variant in out.iter_mut() {
        if variant.variant_type == "SNP" || variant.filter.as_deref() != Some("Strand_bias") {
            continue;
        }
        let name = generate_variant_name(variant);
        let (gfwd, grev) = graph_counts.get(&name).copied().unwrap_or((0, 0));
        let alt: HashSet<&str> = bam_read_record
            .get(&name)
            .map(|reads| reads.iter().filter_map(|r| r.as_str()).collect())
            .unwrap_or_default();
        let (mut fwd, mut rev) = (0usize, 0usize);
        for r in &alt {
            match strand_map.get(*r) {
                Some(true) => rev += 1,
                Some(false) => fwd += 1,
                None => {}
            }
        }
        let total = fwd + rev;
        if total < min_reads {
            println!(
                "Strand-bias BAM recount, {:?}, graph fwd={}, rev={}; BAM alt reads={} < {}: kept flagged",
                name, gfwd, grev, total, min_reads
            );
            continue;
        }
        let (mut cfwd, mut crev) = (0usize, 0usize);
        if let Some(cover) = bam_cover_record.get(&name) {
            for r in cover.iter().filter(|r| !alt.contains(r.as_str())) {
                match strand_map.get(r) {
                    Some(true) => crev += 1,
                    Some(false) => cfwd += 1,
                    None => {}
                }
            }
        }
        let null_frac = if cfwd + crev >= crate::STRAND_LOCAL_MIN_READS {
            cfwd as f64 / (cfwd + crev) as f64
        } else {
            expected_fwd_frac
        };
        let alt_frac = fwd as f64 / total as f64;
        let effect = (alt_frac - null_frac).abs();
        let p_value = strand_bias_pvalue(fwd, rev, null_frac);
        let cleared = p_value >= p_threshold || effect < min_effect;
        println!(
            "Strand-bias BAM recount, {:?}, graph fwd={}, rev={}; BAM fwd={}, rev={}, alt_frac={:.3}, null_frac={:.3}, effect={:.3}, p={:.3e}: {}",
            name, gfwd, grev, fwd, rev, alt_frac, null_frac, effect, p_value,
            if cleared { "cleared" } else { "kept flagged" }
        );
        if cleared {
            variant.filter = None;
        }
    }
    out
}

/// Genotype only the `Strand_bias`-flagged indels in `bam`.
#[allow(clippy::too_many_arguments)]
fn recount_flagged_indels_from_bam(
    variants: &[Variant],
    graph_read_record: &HashMap<String, Vec<serde_json::Value>>,
    graph_strand_map: &HashMap<String, bool>,
    bam: &Path,
    refseq: &[u8],
    ref_name: &str,
    opts: &GenotypeOpts,
    p_threshold: f64,
) -> Vec<Variant> {
    let flagged: Vec<Variant> = variants
        .iter()
        .filter(|v| v.variant_type != "SNP" && v.filter.as_deref() == Some("Strand_bias"))
        .cloned()
        .collect();
    if flagged.is_empty() {
        return variants.to_vec();
    }
    let graph_counts: HashMap<String, (usize, usize)> = flagged
        .iter()
        .map(|v| {
            let name = generate_variant_name(v);
            let reads: HashSet<&str> = graph_read_record
                .get(&name)
                .map(|rs| rs.iter().filter_map(|r| r.as_str()).collect())
                .unwrap_or_default();
            let rev = reads.iter().filter(|r| graph_strand_map.get(**r) == Some(&true)).count();
            let fwd = reads.iter().filter(|r| graph_strand_map.get(**r) == Some(&false)).count();
            (name, (fwd, rev))
        })
        .collect();
    let (read_record, cover_record) = genotype_from_bam(bam, refseq, ref_name, &flagged, opts);
    let (strand_map, fwd_frac) = build_strand_map(bam);
    println!(
        "Strand-bias BAM recount: {} flagged indel(s) from {} (forward fraction {:.3})",
        flagged.len(),
        bam.display(),
        fwd_frac
    );
    recount_indel_strand_bias(
        variants,
        &graph_counts,
        &read_record,
        &cover_record,
        &strand_map,
        fwd_frac,
        p_threshold,
        2,
        crate::STRAND_BIAS_MIN_EFFECT,
    )
}

// ---- adaptive permutation (ONT-recall T1) ----

/// How `permutation_test` draws its null.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PermutationMode {
    /// Draw `--permutation-rounds` per tested variant.
    Adaptive,
    Fixed,
}

/// A genotype row as bitsets over reads.
pub(crate) struct BitRow {
    pub alt: Vec<u64>,
    pub cov: Vec<u64>,
    /// Indices of the covered reads: the pool a shuffle draws its alt reads from.
    pub covered_idx: Vec<u32>,
    pub n_alt: u32,
}

impl BitRow {
    pub(crate) fn from_row(row: &[f64]) -> BitRow {
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

/// Jaccard of a (possibly shuffled) row `x` against row `y`.
#[inline]
pub(crate) fn jaccard_bits(x_alt: &[u64], x_cov: &[u64], y: &BitRow) -> f64 {
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

/// Bitset analogue of `sum_jaccard_against_others`.
#[inline]
fn sum_jaccard_bits(x_alt: &[u64], x_cov: &[u64], self_index: usize, rows: &[BitRow]) -> f64 {
    let mut sum = 0.0;
    for (j, y) in rows.iter().enumerate() {
        if j != self_index {
            sum += jaccard_bits(x_alt, x_cov, y);
        }
    }
    sum
}

/// Bitset analogue of `shuffle_genotypes`.
pub(crate) fn shuffle_bits(
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

/// Materialize the matrix rows once, as the permutation code compares each row against every other.
fn matrix_rows(matrix: &Array2<f64>) -> Vec<Vec<f64>> {
    (0..matrix.nrows())
        .map(|j| matrix.slice(s![j, ..]).iter().copied().collect())
        .collect()
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
pub(crate) struct PermOutcome {
    pub exceed: u32,
    pub draws: u32,
}

impl PermOutcome {
    pub(crate) fn p_value(&self) -> f64 {
        (self.exceed + 1) as f64 / (self.draws + 1) as f64
    }
}

pub(crate) struct AdaptiveResult {
    /// One outcome per entry of `tested`, in the same order.
    pub outcomes: Vec<PermOutcome>,
    /// Final per-variant round target of the floor-resolving stage.
    pub r_target: usize,
    pub total_draws: u64,
}

/// Draws per batch; the early-stop test runs between batches.
const PERM_BATCH: u32 = 100;
/// Most floor-resolving passes (each only ever raises the target).
const PERM_MAX_STAGE2_PASSES: usize = 5;

struct PermState {
    rng: rand::rngs::StdRng,
    exceed: u32,
    draws: u32,
}

/// Stage 1 draws up to `min_rounds` per tested variant.
pub(crate) fn adaptive_permutation(
    rows: &[Vec<f64>],
    tested: &[usize],
    names: &[String],
    min_rounds: usize,
    alpha: f64,
) -> AdaptiveResult {
    use rand::SeedableRng;
    let started = std::time::Instant::now();
    let m = tested.len();
    if m == 0 || min_rounds == 0 {
        return AdaptiveResult {
            outcomes: vec![PermOutcome { exceed: 0, draws: 0 }; m],
            r_target: min_rounds,
            total_draws: 0,
        };
    }
    let stop_after = crate::PERM_STOP_AFTER;
    let bits: Vec<BitRow> = rows.iter().map(|r| BitRow::from_row(r)).collect();
    let words = bits.first().map_or(0, |b| b.alt.len());
    let obs: Vec<f64> = tested
        .iter()
        .map(|&i| sum_jaccard_bits(&bits[i].alt, &bits[i].cov, i, &bits))
        .collect();
    let mut states: Vec<PermState> = tested
        .iter()
        .map(|&i| PermState {
            rng: rand::rngs::StdRng::seed_from_u64(crate::PERM_SEED ^ fnv1a64(names[i].as_bytes())),
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
                while st.draws < target && st.exceed < stop_after {
                    let n = PERM_BATCH.min(target - st.draws);
                    for _ in 0..n {
                        shuffle_bits(row, scratch, buf, &mut st.rng);
                        if sum_jaccard_bits(buf, &row.cov, i, &bits) >= threshold {
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
        let new_target = crate::PERM_MAX_ROUNDS
            .min((crate::PERM_HEADROOM * m as f64 / (alpha * floor.len() as f64)).ceil() as usize);
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

    AdaptiveResult {
        outcomes: states
            .iter()
            .map(|s| PermOutcome { exceed: s.exceed, draws: s.draws })
            .collect(),
        r_target: target,
        total_draws: total(&states),
    }
}

fn permutation_test(
    matrix: &Array2<f64>,
    p_value_threshold: f64,
    heteroplasmic_error_threshold:f64,
    permutation_round: usize,
    filtered_var: &Vec<Variant>,
    coverage: &HashMap<usize, usize>,
    permutation_frequency_threshold:f64,
    data_type: &str,
    mode: PermutationMode,
) -> (Vec<Variant>, Array2<f64>, Vec<String>) {

    let bar = ProgressBar::new(filtered_var.len() as u64);

    // Admission: which variants are candidates for the permutation test at all.
    let admitted = |i: usize| -> bool {
        let current_variant = &filtered_var[i];
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
        // exclude large indel
        if (current_variant.ref_allele.len() as i32 - current_variant.alt_allele.len() as i32).abs() > 50 {
            return false;
        }

        let is_snp = current_variant.ref_allele.len() == 1 && current_variant.alt_allele.len() == 1;

        // Jaccard co-occurrence penalizes variants on the majority lineage, so exempt high-HF, high-AC variants.
        if is_snp
            && current_variant.filter.is_none()
            && frequency >= crate::COOC_EXEMPT_HF
            && current_variant.allele_count >= crate::COOC_EXEMPT_MIN_AC
        {
            return false;
        }

        // PacBio: exclude SNPs from the permutation test (SNP accuracy is high enough).
        !(is_snp && data_type == "pacbio")
    };

    // Fixed: a per-variant null of `permutation_round` draws for every variant.
    let precomputed_p = mode == PermutationMode::Adaptive;
    let (nulls, adaptive_p): (HashMap<String, Vec<f64>>, HashMap<usize, f64>) = match mode {
        PermutationMode::Fixed => (
            get_null_distribution(filtered_var, &matrix, permutation_round, permutation_frequency_threshold),
            HashMap::new(),
        ),
        PermutationMode::Adaptive => {
            let rows = matrix_rows(matrix);
            let names: Vec<String> = filtered_var.iter().map(generate_variant_name).collect();
            // No null for a variant above `permutation_frequency_threshold` by MATRIX alt frequency.
            let tested: Vec<usize> = (0..filtered_var.len())
                .filter(|&i| {
                    admitted(i)
                        && alt_frequency(matrix.slice(s![i, ..])) <= permutation_frequency_threshold
                })
                .collect();
            let res = adaptive_permutation(&rows, &tested, &names, permutation_round, p_value_threshold);
            let p: HashMap<usize, f64> = tested
                .iter()
                .zip(res.outcomes.iter())
                .filter(|(_, o)| o.draws > 0)
                .map(|(&i, o)| (i, o.p_value()))
                .collect();
            (HashMap::new(), p)
        }
    };
    let collected_values: Vec<Option<(f64, String)>> = (0..filtered_var.len()).into_par_iter().map(|i| {
        bar.inc(1);
        let current_variant = &filtered_var[i];
        let index = generate_variant_name(current_variant);
        if !admitted(i) {
            return None;
        }

        if precomputed_p {
            // Untested (no draws) is missing evidence, not evidence of an artifact.
            return adaptive_p.get(&i).map(|&p| (p, index));
        }

        // No null draws for this variant.
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

pub fn resolve_indel_false_threshold(data_type: &str, explicit: Option<f64>) -> f64 {
    explicit.unwrap_or(match data_type {
        "ont-r9" | "ont-r10" | "ont-denoised" => 0.015,
        _ => 0.1,
    })
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
    permutation_mode: PermutationMode,
    strand_bias_action: Option<StrandBiasAction>,
    strand_bias_drop_max_hf: f64,
    left_align: bool,
    no_linkage_rescue: bool,
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
    let action = resolve_strand_bias_action(data_type, strand_bias_action);
    // The matrix gate must not move the VCF, so the indel recount below is ungated.
    let recount_opts = GenotypeOpts { flank: crate::GENOTYPE_FLANK, ..GenotypeOpts::default() };
    let linkage_bam = bam.filter(|_| linkage_rescue_enabled(data_type, no_linkage_rescue));
    let (sb_filtered_var, sb_matrix) = match bam {
        Some(bam) => {
            let (strand_map, fwd_frac) = build_strand_map(bam);
            println!("Library forward-strand fraction: {:.3}", fwd_frac);
            println!("Strand-bias action: {:?} (SNPs below HF {:.3})", action, strand_bias_drop_max_hf);
            let (sb_var, sb_matrix) = filter_strand_bias(
                &var_record,
                &matrix,
                &read_record,
                &strand_map,
                &cover_record,
                &coverage,
                fwd_frac,
                strand_bias_threshold,
                2,
                permutation_frequency_threshold,
                action,
                strand_bias_drop_max_hf,
                crate::STRAND_BIAS_MIN_EFFECT,
            );
            let sb_var = match linkage_bam {
                Some(raw) => recount_flagged_indels_from_bam(
                    &sb_var,
                    &read_record,
                    &strand_map,
                    raw,
                    ref_seq.as_bytes(),
                    &ref_header,
                    &recount_opts,
                    strand_bias_threshold,
                ),
                None => sb_var,
            };
            (sb_var, sb_matrix)
        }
        None => (var_record.clone(), matrix.clone()),
    };

    let (permu_filtered_var, filtered_matrix, _filtered_name) = permutation_test(&sb_matrix, p_value_threshold, heteroplasmic_frequency_threshold, permutation_rounds, &sb_filtered_var, &coverage, permutation_frequency_threshold, data_type, permutation_mode);

    let link_calls: Option<Vec<LinkCall>> = linkage_bam.and_then(|bam| {
        match linkage::rescue(bam, ref_seq.as_bytes(), &ref_header, &linkage::LinkOpts::default()) {
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
        opts: resolve_genotype_opts(data_type),
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
pub struct GenotypeSource<'a> {
    pub bam: &'a Path,
    pub refseq: &'a [u8],
    pub ref_name: &'a str,
    pub opts: GenotypeOpts,
}

/// Where `<output>.matrix.csv` came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MatrixSource {
    Bam,
    Graph,
}

const GRAPH_MATRIX_WARNING: &str =
    "Warning: no --bam given; matrix.csv comes from the graph, not from re-genotyped reads.";

/// Writes `<output>.matrix.csv`: re-genotyped from the BAM.
fn write_final_matrix(
    output_file: &Path,
    graph_matrix: &Array2<f64>,
    variants: &[Variant],
    rescued: &[Variant],
    graph_reads: &[String],
    genotype: Option<&GenotypeSource>,
) -> Result<MatrixSource, Box<dyn Error>> {
    let matrix_output = output_file.with_extension("matrix.csv");
    match genotype {
        None => {
            eprintln!("{GRAPH_MATRIX_WARNING}");
            write_matrix_to_csv(graph_matrix, variants, graph_reads, matrix_output)?;
            Ok(MatrixSource::Graph)
        }
        Some(src) => {
            let all: Vec<Variant> = variants.iter().chain(rescued).cloned().collect();
            let (read_record, cover_record) =
                genotype_from_bam(src.bam, src.refseq, src.ref_name, &all, &src.opts);
            let (matrix, var_record, read_set) =
                construct_matrix(&read_record, &cover_record, &all);
            println!(
                "Genotyped {} variants x {} reads from {}",
                var_record.len(),
                read_set.len(),
                src.bam.display()
            );
            write_matrix_to_csv(&matrix, &var_record, &read_set, matrix_output)?;
            Ok(MatrixSource::Bam)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// Two bubbles: a normal one and a wrap-around one crossing the circular origin.
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

        // Raw position 11 lies on the wrap-around bubble [8, 12) and circularizes down to wrapped position 1.
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

        // The wrap-around bubble's true reads are r4/r5.
        assert_eq!(cover, vec!["r4".to_string(), "r5".to_string()]);
    }

    /// A bubble carrying the real variant, plus a read that crosses it via a skip-over edge with no CIGAR.
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

        // r3 only reaches this span via E_private, which was never CIGAR-processed.
        assert_eq!(cover, vec!["r1".to_string(), "r2".to_string()]);
    }

    /// A genuinely circular graph: the last anchor links straight back to the first.
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

        // get_variants_from_cigar walks forward from the src anchor.
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
        // SOURCE/SINK have no anchor coordinate.
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
        // Reads are the cells SCITE reasons over.
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
        // p = (1 + #{null >= obs}) / (rounds + 1): the +1 keeps p strictly positive.
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
        // The summary statistic is a SUM of Jaccard coefficients, so it scales with a variant's alt-read count.
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

        // Shuffling preserves each variant's own alt count.
        let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
        assert!(
            mean(&nulls[&n0]) > mean(&nulls[&n1]),
            "the 8-read variant's null must sit above the 2-read variant's"
        );
    }

    #[test]
    fn null_distribution_skips_variants_above_the_frequency_threshold() {
        // Near-homoplasmic variants are the co-occurrence backbone, not candidates.
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

    /// Two co-segregating variants carried by the same 3 of 30 reads.
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
        // Denoising strips the RANDOM error its quality model can see and leaves the SYSTEMATIC error behind.
        let (variants, matrix, coverage) = lineage_and_artifact();
        let artifact = generate_variant_name(&variants[2]);

        let (kept, _m, _names) = permutation_test(
            &matrix, 0.5, 0.9, 200, &variants, &coverage, 0.9, "ont-denoised",
            PermutationMode::Fixed,
        );
        let kept_names: Vec<String> = kept.iter().map(generate_variant_name).collect();

        assert!(
            !kept_names.contains(&artifact),
            "a scattered SNP must be testable, and rejected, under ont-denoised"
        );
    }

    #[test]
    fn permutation_test_still_exempts_snps_for_pacbio() {
        // PacBio SNPs are excluded from the permutation test.
        let (variants, matrix, coverage) = lineage_and_artifact();
        let artifact = generate_variant_name(&variants[2]);

        let (kept, _m, _names) = permutation_test(
            &matrix, 0.5, 0.9, 200, &variants, &coverage, 0.9, "pacbio",
            PermutationMode::Fixed,
        );
        let kept_names: Vec<String> = kept.iter().map(generate_variant_name).collect();

        assert!(
            kept_names.contains(&artifact),
            "pacbio must keep SNPs exempt from the permutation test"
        );
    }

    /// X is a real heteroplasmy on the majority lineage (three exclusive lineages A, B, C).
    fn majority_lineage_heteroplasmy(x: Variant) -> (Vec<Variant>, Array2<f64>, HashMap<usize, usize>) {
        let n_reads = 300;
        let mut rng = StdRng::seed_from_u64(12969);
        let mut variants = vec![x];
        let mut rows: Vec<Vec<f64>> = vec![vec![0.0; n_reads]];
        // 45 random reads among A (reads 0..240)
        let mut placed = 0;
        while placed < 45 {
            let j = rng.random_range(0..240);
            if rows[0][j] == 0.0 {
                rows[0][j] = 1.0;
                placed += 1;
            }
        }
        let lineages = [(0usize, 240usize, 3usize), (240, 276, 15), (276, 300, 15)];
        let mut pos = 1000;
        for (lo, hi, n_markers) in lineages {
            for _ in 0..n_markers {
                variants.push(snp(pos, "G", hi - lo));
                let mut row = vec![0.0; n_reads];
                for cell in row.iter_mut().take(hi).skip(lo) {
                    *cell = 1.0;
                }
                rows.push(row);
                pos += 10;
            }
        }
        let mut matrix = Array2::<f64>::zeros((rows.len(), n_reads));
        for (i, row) in rows.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                matrix[[i, j]] = *v;
            }
        }
        let coverage: HashMap<usize, usize> =
            variants.iter().map(|v| (v.pos, n_reads)).collect();
        (variants, matrix, coverage)
    }

    fn run_majority_lineage(x: Variant) -> Vec<String> {
        let (variants, matrix, coverage) = majority_lineage_heteroplasmy(x);
        let (kept, _m, _n) = permutation_test(
            &matrix, 0.01, 0.2, 1000, &variants, &coverage, 0.7, "ont-r10",
            PermutationMode::Adaptive,
        );
        let kept_names: Vec<String> = kept.iter().map(generate_variant_name).collect();
        // sanity: all 33 lineage markers (A untested, B/C tested) must survive
        assert_eq!(
            kept.iter().filter(|v| v.pos >= 1000).count(),
            33,
            "lineage markers must survive the test"
        );
        kept_names
    }

    #[test]
    fn cooccurrence_test_exempts_strong_balanced_snps() {
        // The co-occurrence statistic punishes a variant that sits on the MAJORITY lineage.
        let x = |filter: Option<&str>| {
            let mut v = snp(500, "C", 45);
            v.filter = filter.map(String::from);
            v
        };

        let kept = run_majority_lineage(x(None));
        assert!(
            kept.contains(&generate_variant_name(&x(None))),
            "a strong strand-balanced SNP must be exempt from the co-occurrence veto"
        );

        let kept = run_majority_lineage(x(Some("Strand_bias")));
        assert!(
            !kept.contains(&generate_variant_name(&x(None))),
            "a strand-biased SNP must still be tested and rejected"
        );

        let mut indel_x = indel(500, "TT", "T", 45);
        indel_x.filter = None;
        let name = generate_variant_name(&indel_x);
        let kept = run_majority_lineage(indel_x);
        assert!(!kept.contains(&name), "an indel with the same evidence must still be tested");
    }

    #[test]
    fn permutation_test_keeps_a_variant_with_no_null_draws() {
        // `get_null_distribution` filters on the MATRIX alt frequency while the test loop filters on allele_count/coverage.
        let variants = vec![snp(100, "A", 2), snp(200, "C", 2)];
        let mut matrix = Array2::<f64>::zeros((2, 10));
        for j in 0..10 {
            // Both variants are alt in 9/10 reads => matrix frequency 0.9, above the 0.7 null threshold, so neither gets a null...
            matrix[[0, j]] = if j < 9 { 1.0 } else { 0.0 };
            matrix[[1, j]] = if j < 9 { 1.0 } else { 0.0 };
        }
        // ...but allele_count/coverage is 2/100 = 0.02, well under the test's own 0.5 heteroplasmic threshold.
        let coverage: HashMap<usize, usize> = [(100, 100), (200, 100)].into_iter().collect();

        let (kept, _m, _names) = permutation_test(
            &matrix, 0.5, 0.5, 16, &variants, &coverage, 0.7, "ont-denoised",
            PermutationMode::Fixed,
        );

        assert_eq!(kept.len(), 2, "variants with an empty null must be kept, not dropped");
    }

    // ---------------------------------------------------------------------

    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    /// `n` reads x `rows` random genotype rows: ~20% NaN (uncovered), ~10% alt.
    fn random_genotype_row(rng: &mut StdRng, n: usize) -> Vec<f64> {
        (0..n)
            .map(|_| {
                let u: f64 = rng.random();
                if u < 0.20 {
                    f64::NAN
                } else if u < 0.30 {
                    1.0
                } else {
                    0.0
                }
            })
            .collect()
    }

    #[test]
    fn bit_jaccard_matches_f64_jaccard() {
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..200 {
            let a = random_genotype_row(&mut rng, 300);
            let b = random_genotype_row(&mut rng, 300);
            let (ba, bb) = (BitRow::from_row(&a), BitRow::from_row(&b));
            let want = genotype_jaccard_similarity(&a, &b);
            let got = jaccard_bits(&ba.alt, &ba.cov, &bb);
            assert_eq!(got, want, "bitset Jaccard must equal the f64 Jaccard exactly");
        }
    }

    #[test]
    fn bit_shuffle_preserves_alt_count_and_coverage() {
        let mut rng = StdRng::seed_from_u64(11);
        let row = random_genotype_row(&mut rng, 300);
        let br = BitRow::from_row(&row);
        let mut scratch: Vec<u32> = Vec::new();
        let mut alt = vec![0u64; br.alt.len()];
        let popcnt = |w: &[u64]| w.iter().map(|x| x.count_ones()).sum::<u32>();
        for _ in 0..1000 {
            shuffle_bits(&br, &mut scratch, &mut alt, &mut rng);
            assert_eq!(popcnt(&alt), br.n_alt, "shuffle must preserve the alt count");
            assert!(
                alt.iter().zip(&br.cov).all(|(a, c)| a & !c == 0),
                "shuffled alts must stay inside the covered reads"
            );
        }
    }

    /// 10 identical lineage rows (alt in the same 24 of 400 reads) plus 120 independent noise rows.
    fn floor_matrix() -> (Vec<Variant>, Array2<f64>, HashMap<usize, usize>) {
        let n_reads = 400;
        let mut rng = StdRng::seed_from_u64(2026);
        let mut variants = Vec::new();
        let mut matrix = Array2::<f64>::zeros((130, n_reads));
        for i in 0..130 {
            variants.push(snp(100 + i, "A", if i < 10 { 24 } else { 8 }));
            if i < 10 {
                for j in 0..24 {
                    matrix[[i, j]] = 1.0;
                }
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
        let coverage: HashMap<usize, usize> = (0..130).map(|i| (100 + i, n_reads)).collect();
        (variants, matrix, coverage)
    }

    fn lineage_kept(kept: &[Variant]) -> (usize, usize) {
        let lineage = kept.iter().filter(|v| v.pos < 110).count();
        (lineage, kept.len() - lineage)
    }

    #[test]
    fn fixed_permutation_cannot_resolve_the_bh_floor() {
        let (variants, matrix, coverage) = floor_matrix();
        let (kept, _m, _n) = permutation_test(
            &matrix, 0.01, 0.2, 1000, &variants, &coverage, 0.7, "ont-r10",
            PermutationMode::Fixed,
        );
        let (lineage, _noise) = lineage_kept(&kept);
        assert_eq!(
            lineage, 0,
            "R = 1000 puts 10 lineage rows at p = 1/1001, so BH q = 0.013 > 0.01: all 10 excluded"
        );
    }

    #[test]
    fn adaptive_permutation_resolves_the_bh_floor() {
        let (variants, matrix, coverage) = floor_matrix();
        let (kept, _m, _n) = permutation_test(
            &matrix, 0.01, 0.2, 1000, &variants, &coverage, 0.7, "ont-r10",
            PermutationMode::Adaptive,
        );
        let (lineage, noise) = lineage_kept(&kept);
        assert_eq!(lineage, 10, "adaptive must keep every lineage row");
        assert!(noise <= 1, "at most one noise row may survive, got {noise}");

        let rows = matrix_rows(&matrix);
        let tested: Vec<usize> = (0..variants.len()).collect();
        let names: Vec<String> = variants.iter().map(generate_variant_name).collect();
        let res = adaptive_permutation(&rows, &tested, &names, 1000, 0.01);
        let expected =
            (crate::PERM_HEADROOM * 130.0 / (0.01 * 10.0)).ceil() as usize;
        assert!(
            res.r_target >= expected,
            "PERM_HEADROOM * 130 / (0.01 * 10) = {expected} rounds needed, reported {}",
            res.r_target
        );
    }

    #[test]
    fn adaptive_permutation_stops_clear_nulls_early() {
        let (variants, matrix, _cov) = floor_matrix();
        let rows = matrix_rows(&matrix);
        let tested: Vec<usize> = (0..variants.len()).collect();
        let names: Vec<String> = variants.iter().map(generate_variant_name).collect();
        let min_rounds = 1000;
        let res = adaptive_permutation(&rows, &tested, &names, min_rounds, 0.01);
        let mut checked = 0;
        for (k, o) in res.outcomes.iter().enumerate() {
            let p = (o.exceed + 1) as f64 / (o.draws + 1) as f64;
            if p > 0.1 {
                checked += 1;
                assert!(
                    (o.draws as usize) < min_rounds,
                    "variant {k} with p = {p:.3} should have stopped early, drew {}",
                    o.draws
                );
            }
        }
        assert!(checked > 50, "expected many clear nulls, saw {checked}");
    }

    #[test]
    fn adaptive_permutation_is_deterministic() {
        let (variants, matrix, _cov) = floor_matrix();
        let rows = matrix_rows(&matrix);
        let tested: Vec<usize> = (0..variants.len()).collect();
        let names: Vec<String> = variants.iter().map(generate_variant_name).collect();
        let a = adaptive_permutation(&rows, &tested, &names, 1000, 0.01);
        let b = adaptive_permutation(&rows, &tested, &names, 1000, 0.01);
        assert_eq!(a.r_target, b.r_target);
        let pv = |r: &AdaptiveResult| -> Vec<u64> {
            r.outcomes
                .iter()
                .map(|o| ((o.exceed + 1) as f64 / (o.draws + 1) as f64).to_bits())
                .collect()
        };
        assert_eq!(pv(&a), pv(&b), "p-values must be byte-identical across runs");
    }

    // ---------------------------------------------------------------------

    /// The noise drowns the pair's summed Jaccard.
    fn small_lineage_matrix() -> (Vec<Variant>, Array2<f64>, HashMap<usize, usize>) {
        let (n_reads, n_noise) = (200, 120);
        let mut rng = StdRng::seed_from_u64(2002);
        let mut variants = vec![snp(100, "A", 6), snp(110, "A", 6)];
        let mut matrix = Array2::<f64>::zeros((2 + n_noise, n_reads));
        for j in [0, 1, 2, 3, 4, 5] {
            matrix[[0, j]] = 1.0;
        }
        for j in [0, 1, 2, 3, 6, 7] {
            matrix[[1, j]] = 1.0;
        }
        for i in 2..2 + n_noise {
            variants.push(snp(200 + i, "A", 8));
            let mut placed = 0;
            while placed < 8 {
                let j = rng.random_range(0..n_reads);
                if matrix[[i, j]] == 0.0 {
                    matrix[[i, j]] = 1.0;
                    placed += 1;
                }
            }
        }
        let coverage: HashMap<usize, usize> = variants.iter().map(|v| (v.pos, n_reads)).collect();
        (variants, matrix, coverage)
    }

    fn run_small_lineage(mode: PermutationMode) -> (usize, usize) {
        let (variants, matrix, coverage) = small_lineage_matrix();
        let (kept, _m, _n) = permutation_test(
            &matrix, 0.01, 0.2, 200, &variants, &coverage, 0.7, "ont-r10", mode,
        );
        let lineage = kept.iter().filter(|v| v.pos < 200).count();
        (lineage, kept.len() - lineage)
    }

    #[test]
    fn sum_jaccard_rejects_a_small_lineage() {
        let (lineage, _) = run_small_lineage(PermutationMode::Fixed);
        assert_eq!(lineage, 0, "the summed Jaccard of a 2-allele lineage is drowned by noise");
    }

    // ---------------------------------------------------------------------

    /// One low-VAF SNP supported by 10 reads, all on the forward strand.
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
            &variants, &matrix, &read_record, &strand_map, &HashMap::new(), &coverage,
            0.5, 0.01, 2, 0.7, StrandBiasAction::Drop, 0.10, crate::STRAND_BIAS_MIN_EFFECT,
        );

        assert!(kept.is_empty(), "a strand-skewed low-VAF SNP must be dropped, not flagged");
        assert_eq!(kept_matrix.nrows(), 0, "the matrix must lose the dropped row");
    }

    #[test]
    fn strand_bias_only_flags_under_the_flag_action() {
        let (variants, matrix, read_record, strand_map, coverage) = strand_skewed_snp();

        let (kept, _m) = filter_strand_bias(
            &variants, &matrix, &read_record, &strand_map, &HashMap::new(), &coverage,
            0.5, 0.01, 2, 0.7, StrandBiasAction::Flag, 0.10, crate::STRAND_BIAS_MIN_EFFECT,
        );

        assert_eq!(kept.len(), 1, "flag mode must keep the call");
        assert_eq!(kept[0].filter.as_deref(), Some("Strand_bias"));
    }

    #[test]
    fn strand_bias_drop_action_spares_indels_and_higher_vaf_snps() {
        // The drop is deliberately narrow: it targets the low-VAF SNPs that denoising leaves behind.
        let (_v, _m, read_record_snp, strand_map, _c) = strand_skewed_snp();

        // Same evidence, but the call sits at 10/50 = 0.20 VAF, above the 0.10 ceiling.
        let variants = vec![snp(100, "A", 10)];
        let matrix = Array2::<f64>::zeros((1, 10));
        let coverage: HashMap<usize, usize> = [(100, 50)].into_iter().collect();
        let (kept, _m) = filter_strand_bias(
            &variants, &matrix, &read_record_snp, &strand_map, &HashMap::new(), &coverage,
            0.5, 0.01, 2, 0.7, StrandBiasAction::Drop, 0.10, crate::STRAND_BIAS_MIN_EFFECT,
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
            &ivariants, &matrix, &iread_record, &strand_map, &HashMap::new(), &icoverage,
            0.5, 0.01, 2, 0.7, StrandBiasAction::Drop, 0.10, crate::STRAND_BIAS_MIN_EFFECT,
        );
        assert_eq!(ikept.len(), 1, "indels are never dropped by the strand-bias filter");
    }

    // ---- column-local strand-bias null (ONT-recall T3) ----

    /// Build a one-SNP scenario: `fwd`/`rev` alt reads, and covering reads in addition to the alt reads.
    #[allow(clippy::type_complexity)]
    fn local_null_case(
        variant: Variant,
        alt_fwd: usize,
        alt_rev: usize,
        cover_fwd: usize,
        cover_rev: usize,
        depth: usize,
    ) -> (
        Vec<Variant>,
        Array2<f64>,
        HashMap<String, Vec<serde_json::Value>>,
        HashMap<String, bool>,
        HashMap<String, HashSet<String>>,
        HashMap<usize, usize>,
    ) {
        let name = generate_variant_name(&variant);
        let mut strand_map = HashMap::new();
        let mut alt_reads = Vec::new();
        for i in 0..alt_fwd {
            strand_map.insert(format!("a{i}"), false);
            alt_reads.push(json!(format!("a{i}")));
        }
        for i in 0..alt_rev {
            strand_map.insert(format!("ar{i}"), true);
            alt_reads.push(json!(format!("ar{i}")));
        }
        // the real `cover_record` set contains the alt reads too
        let mut cover: HashSet<String> = strand_map.keys().cloned().collect();
        for i in 0..cover_fwd {
            strand_map.insert(format!("c{i}"), false);
            cover.insert(format!("c{i}"));
        }
        for i in 0..cover_rev {
            strand_map.insert(format!("cr{i}"), true);
            cover.insert(format!("cr{i}"));
        }
        let read_record: HashMap<String, Vec<serde_json::Value>> =
            [(name.clone(), alt_reads)].into_iter().collect();
        let cover_record: HashMap<String, HashSet<String>> =
            [(name, cover)].into_iter().collect();
        let coverage: HashMap<usize, usize> = [(variant.pos, depth)].into_iter().collect();
        let n = alt_fwd + alt_rev;
        (vec![variant], Array2::<f64>::zeros((1, n.max(1))), read_record, strand_map, cover_record, coverage)
    }

    fn run_sb(
        case: &(
            Vec<Variant>,
            Array2<f64>,
            HashMap<String, Vec<serde_json::Value>>,
            HashMap<String, bool>,
            HashMap<String, HashSet<String>>,
            HashMap<usize, usize>,
        ),
        library_frac: f64,
        action: StrandBiasAction,
    ) -> Vec<Variant> {
        filter_strand_bias(
            &case.0, &case.1, &case.2, &case.3, &case.4, &case.5,
            library_frac, 0.01, 2, 0.7, action, 0.10, crate::STRAND_BIAS_MIN_EFFECT,
        )
        .0
    }

    #[test]
    fn strand_bias_uses_column_local_null() {
        // 7805G>A: alt 45F/27R, column 54% fwd.
        let case = local_null_case(snp(7805, "A", 72), 45, 27, 54, 46, 1000);
        let local_p = strand_bias_pvalue(45, 27, 0.54);
        assert!(local_p >= 0.01, "local null p={local_p}");
        for action in [StrandBiasAction::Flag, StrandBiasAction::Drop] {
            let kept = run_sb(&case, 0.498, action);
            assert_eq!(kept.len(), 1);
            assert_eq!(kept[0].filter, None, "{action:?}: must not be flagged against the local null");
        }
    }

    #[test]
    fn strand_bias_needs_an_effect_not_just_a_p_value() {
        // 1000 alt reads at 600/400 vs a 50% local null: p~1e-10, effect 0.10.
        let case = local_null_case(snp(500, "A", 1000), 600, 400, 50, 50, 5000);
        assert!(strand_bias_pvalue(600, 400, 0.5) < 1e-8);
        let kept = run_sb(&case, 0.5, StrandBiasAction::Flag);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].filter, None, "a 0.10 effect is below STRAND_BIAS_MIN_EFFECT");
    }

    #[test]
    fn strand_bias_still_flags_single_strand_artifacts() {
        // 33CG>C: alt 1F/54R against a local 0.536 fwd null.
        let case = local_null_case(indel(33, "CG", "C", 55), 1, 54, 536, 464, 2000);
        let kept = run_sb(&case, 0.5, StrandBiasAction::Flag);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].filter.as_deref(), Some("Strand_bias"));
        // an SNP with the same evidence is flagged under Flag and removed under Drop
        let case = local_null_case(snp(33, "C", 55), 1, 54, 536, 464, 2000);
        assert_eq!(run_sb(&case, 0.5, StrandBiasAction::Flag)[0].filter.as_deref(), Some("Strand_bias"));
        assert!(run_sb(&case, 0.5, StrandBiasAction::Drop).is_empty());
    }

    // ---- BAM recount of flagged indels (ONT-recall R3-T3) ----

    /// Flag `v` from graph counts 10/0 against a balanced column.
    fn graph_flagged_then_bam_recount(v: Variant, bam_fwd: usize, bam_rev: usize) -> Vec<Variant> {
        let graph = local_null_case(v.clone(), 10, 0, 50, 50, 1000);
        let flagged = run_sb(&graph, 0.5, StrandBiasAction::Flag);
        assert_eq!(flagged[0].filter.as_deref(), Some("Strand_bias"), "graph 10/0 must flag");
        let bam = local_null_case(v, bam_fwd, bam_rev, 50, 50, 1000);
        recount_indel_strand_bias(
            &flagged, &HashMap::new(), &bam.2, &bam.4, &bam.3, 0.5, 0.01, 2,
            crate::STRAND_BIAS_MIN_EFFECT,
        )
    }

    #[test]
    fn bam_recount_clears_strand_bias_on_a_balanced_indel() {
        // 302A>ACC: the graph edges split the alt reads by strand, the BAM reads don't.
        let out = graph_flagged_then_bam_recount(indel(302, "A", "ACC", 20), 9, 11);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].filter, None, "9/11 vs a balanced column must clear the flag");
    }

    #[test]
    fn bam_recount_keeps_strand_bias_on_a_skewed_indel() {
        let out = graph_flagged_then_bam_recount(indel(33, "CG", "C", 20), 1, 19);
        assert_eq!(out[0].filter.as_deref(), Some("Strand_bias"), "1/19 stays flagged");
    }

    #[test]
    fn bam_recount_leaves_snvs_and_other_filters_alone() {
        let out = graph_flagged_then_bam_recount(snp(500, "A", 20), 9, 11);
        assert_eq!(out[0].filter.as_deref(), Some("Strand_bias"), "SNVs are not recounted");

        let mut other = indel(302, "A", "ACC", 20);
        other.filter = Some("Potential_Artifact".into());
        let bam = local_null_case(other.clone(), 9, 11, 50, 50, 1000);
        let out = recount_indel_strand_bias(
            &[other], &HashMap::new(), &bam.2, &bam.4, &bam.3, 0.5, 0.01, 2,
            crate::STRAND_BIAS_MIN_EFFECT,
        );
        assert_eq!(out[0].filter.as_deref(), Some("Potential_Artifact"));
    }

    #[test]
    fn bam_recount_keeps_the_flag_when_the_bam_has_too_few_alt_reads() {
        let out = graph_flagged_then_bam_recount(indel(302, "A", "ACC", 20), 1, 0);
        assert_eq!(out[0].filter.as_deref(), Some("Strand_bias"));
    }

    #[test]
    fn strand_bias_local_null_excludes_the_alt_reads_themselves() {
        // The edge's cover set can be dominated by reverse-only alt reads.
        let case = local_null_case(indel(33, "CG", "C", 34), 0, 34, 2, 2, 302);
        assert_eq!(case.4.values().next().unwrap().len(), 38);
        let kept = run_sb(&case, 0.498, StrandBiasAction::Flag);
        assert_eq!(kept[0].filter.as_deref(), Some("Strand_bias"));
    }

    #[test]
    fn strand_bias_falls_back_to_library_fraction_with_few_covering_reads() {
        // Only 10 strand-resolved covering reads: the local fraction is ignored and the library fraction is the null.
        let case = local_null_case(snp(100, "A", 10), 10, 0, 10, 0, 200);
        let kept = run_sb(&case, 0.5, StrandBiasAction::Flag);
        assert_eq!(kept[0].filter.as_deref(), Some("Strand_bias"));
    }

    // ===================================================================

    /// `TCAAAAGCTT` - a 4bp A-run at indices 2..6 preceded by a C.
    const HOMOPOLYMER_REF: &[u8] = b"TCAAAAGCTT";

    /// Build a DEL the way `get_variants_from_cigar` does.
    fn del_at(refseq: &[u8], pos: usize, len: usize, allele_count: usize) -> Variant {
        Variant {
            pos,
            ref_allele: String::from_utf8(refseq[pos - 1..pos + len].to_vec()).unwrap(),
            alt_allele: String::from_utf8(refseq[pos - 1..pos].to_vec()).unwrap(),
            variant_type: "DEL".to_string(),
            allele_count,
            filter: None,
        }
    }

    /// Build an INS the way `get_variants_from_cigar` does.
    fn ins_at(refseq: &[u8], pos: usize, inserted: &str, allele_count: usize) -> Variant {
        let anchor = String::from_utf8(refseq[pos - 1..pos].to_vec()).unwrap();
        Variant {
            pos,
            ref_allele: anchor.clone(),
            alt_allele: format!("{}{}", anchor, inserted),
            variant_type: "INS".to_string(),
            allele_count,
            filter: None,
        }
    }

    #[test]
    fn left_align_shifts_homopolymer_deletion_to_start_of_run() {
        // The aligner put the deleted A at the END of the run. bcftools norm reports this at the base before the run: POS 2.
        let v = del_at(HOMOPOLYMER_REF, 5, 1, 3);
        let aligned = left_align_variants(&[v], HOMOPOLYMER_REF);

        assert_eq!(aligned[0].pos, 2);
        assert_eq!(aligned[0].ref_allele, "CA");
        assert_eq!(aligned[0].alt_allele, "C");
        assert_eq!(aligned[0].variant_type, "DEL");
        assert_eq!(aligned[0].allele_count, 3, "alignment must not touch evidence");
    }

    #[test]
    fn left_align_shifts_homopolymer_insertion_to_start_of_run() {
        // An extra A appended to the run: REF A, ALT AA at pos 6.
        let v = ins_at(HOMOPOLYMER_REF, 6, "A", 3);
        let aligned = left_align_variants(&[v], HOMOPOLYMER_REF);

        assert_eq!(aligned[0].pos, 2);
        assert_eq!(aligned[0].ref_allele, "C");
        assert_eq!(aligned[0].alt_allele, "CA");
        assert_eq!(aligned[0].variant_type, "INS");
    }

    #[test]
    fn left_align_rotates_multi_base_insertion_across_a_tandem_repeat() {
        // `TTAGAGAGC`: an AG dinucleotide repeat at 2..8.
        let refseq: &[u8] = b"TTAGAGAGC";
        let v = ins_at(refseq, 8, "AG", 2);
        let aligned = left_align_variants(&[v], refseq);

        assert_eq!(aligned[0].pos, 2);
        assert_eq!(aligned[0].ref_allele, "T");
        assert_eq!(aligned[0].alt_allele, "TAG");
    }

    #[test]
    fn left_align_leaves_indels_in_unique_sequence_untouched() {
        // Deleting the G of...TTGCA... has nowhere to go: the base entering the deletion differs from the base leaving it (G).
        let refseq: &[u8] = b"ACGTTGCA";
        let v = del_at(refseq, 5, 1, 4);
        let aligned = left_align_variants(&[v.clone()], refseq);

        assert_eq!(aligned[0].pos, v.pos);
        assert_eq!(aligned[0].ref_allele, v.ref_allele);
        assert_eq!(aligned[0].alt_allele, v.alt_allele);
    }

    #[test]
    fn left_align_leaves_snps_untouched() {
        // A SNP inside the homopolymer would shift if it were treated as an indel.
        let v = Variant {
            pos: 4,
            ref_allele: "A".to_string(),
            alt_allele: "G".to_string(),
            variant_type: "SNP".to_string(),
            allele_count: 7,
            filter: None,
        };
        let aligned = left_align_variants(&[v.clone()], HOMOPOLYMER_REF);

        assert_eq!(aligned[0].pos, 4);
        assert_eq!(aligned[0].ref_allele, "A");
        assert_eq!(aligned[0].alt_allele, "G");
    }

    #[test]
    fn left_align_preserves_length_and_order_for_anchorless_variants() {
        // `record_variant_reads` zips the raw and circularized variant vectors BY INDEX.
        let anchorless_ins = Variant {
            pos: 3,
            ref_allele: "-".to_string(),
            alt_allele: "-".to_string(),
            variant_type: "INS".to_string(),
            allele_count: 2,
            filter: None,
        };
        let at_origin = Variant {
            pos: 0,
            ref_allele: "TC".to_string(),
            alt_allele: "T".to_string(),
            variant_type: "DEL".to_string(),
            allele_count: 5,
            filter: None,
        };
        let shiftable = del_at(HOMOPOLYMER_REF, 5, 1, 3);

        let input = vec![anchorless_ins.clone(), shiftable, at_origin.clone()];
        let aligned = left_align_variants(&input, HOMOPOLYMER_REF);

        assert_eq!(aligned.len(), 3, "1:1 map: no filtering, no deduplication");
        // Index 0 and 2 are untouched, and stay at their original indices.
        assert_eq!(aligned[0].ref_allele, anchorless_ins.ref_allele);
        assert_eq!(aligned[0].allele_count, anchorless_ins.allele_count);
        assert_eq!(aligned[2].pos, at_origin.pos);
        assert_eq!(aligned[2].ref_allele, at_origin.ref_allele);
        // Index 1 is the one that moved.
        assert_eq!(aligned[1].pos, 2);
    }

    #[test]
    fn left_align_does_not_shift_past_the_reference_start() {
        // `AAAAC`: the run reaches the start of the reference.
        let refseq: &[u8] = b"AAAAC";
        let v = del_at(refseq, 3, 1, 2);
        let aligned = left_align_variants(&[v], refseq);

        assert_eq!(aligned[0].pos, 1, "shifts as far as an anchor base allows");
        assert_eq!(aligned[0].ref_allele, "AA");
        assert_eq!(aligned[0].alt_allele, "A");
    }

    #[test]
    fn left_align_lets_collapse_merge_equivalent_spellings_and_sum_evidence() {
        // This is the point of doing alignment before collapsing.
        let spelling_a = del_at(HOMOPOLYMER_REF, 5, 1, 3);
        let spelling_b = del_at(HOMOPOLYMER_REF, 3, 1, 4);
        let ref_length = HOMOPOLYMER_REF.len();

        let raw = collapse_identical_records(
            vec![spelling_a.clone(), spelling_b.clone()],
            ref_length,
        );
        assert_eq!(raw.len(), 2, "without alignment they stay separate records");

        let aligned = left_align_variants(&[spelling_a, spelling_b], HOMOPOLYMER_REF);
        let collapsed = collapse_identical_records(aligned, ref_length);

        assert_eq!(collapsed.len(), 1);
        assert_eq!(collapsed[0].pos, 2);
        assert_eq!(collapsed[0].ref_allele, "CA");
        assert_eq!(collapsed[0].alt_allele, "C");
        assert_eq!(collapsed[0].allele_count, 7, "evidence from both spellings sums");
    }

    #[test]
    fn matrix_keys_follow_the_aligned_variant_not_the_raw_one() {
        // The read/cover records are keyed by `generate_variant_name` of the circularized variant.
        let graph = wraparound_test_graph();
        let dict = get_graph_intervals(&graph, 10);

        // Raw position 11 sits in the wrap-around bubble [8, 12).
        let raw = vec![Variant { pos: 11, ..del_at(HOMOPOLYMER_REF, 5, 1, 3) }];
        let aligned = left_align_variants(&[del_at(HOMOPOLYMER_REF, 5, 1, 3)], HOMOPOLYMER_REF);

        let mut read_record = HashMap::new();
        let mut cover_record = HashMap::new();
        record_variant_reads(
            &raw,
            &aligned,
            &dict,
            &graph,
            &[json!("r4")],
            &mut read_record,
            &mut cover_record,
        );

        let aligned_key = generate_variant_name(&aligned[0]);
        assert_eq!(aligned_key, "m.2CA>C");
        assert!(read_record.contains_key(&aligned_key), "keyed by aligned name");
        assert!(
            !read_record.contains_key(&generate_variant_name(&raw[0])),
            "the raw spelling must not appear as a key"
        );
        // ...and the raw position still resolved the covering bubble.
        let mut cover: Vec<String> =
            cover_record.get(&aligned_key).cloned().unwrap_or_default().into_iter().collect();
        cover.sort();
        assert_eq!(cover, vec!["r4".to_string(), "r5".to_string()]);
    }

    /// One bubble with two alt edges spelling the same 1bp deletion at different anchors.
    fn homopolymer_bubble_graph() -> GraphicalGenome {
        let mut anchor = HashMap::new();
        anchor.insert("A1".to_string(), json!({"pos": 0, "seq": "TCA"}));
        anchor.insert("A2".to_string(), json!({"pos": 7, "seq": "CTT"}));

        let mut edges = HashMap::new();
        edges.insert(
            "E_ref".to_string(),
            json!({"src": ["A1"], "dst": ["A2"], "seq": "AAAG", "reads": ["chrM_ref"], "variants": "7="}),
        );
        edges.insert(
            "E_alt1".to_string(),
            json!({"src": ["A1"], "dst": ["A2"], "seq": "AAG", "reads": ["r1", "r2", "r3"], "variants": "5=1D1="}),
        );
        edges.insert(
            "E_alt2".to_string(),
            json!({"src": ["A1"], "dst": ["A2"], "seq": "AAG", "reads": ["r4", "r5", "r6", "r7"], "variants": "3=1D3="}),
        );

        let mut outgoing = HashMap::new();
        outgoing.insert(
            "A1".to_string(),
            vec!["E_ref".to_string(), "E_alt1".to_string(), "E_alt2".to_string()],
        );
        for e in ["E_ref", "E_alt1", "E_alt2"] {
            outgoing.insert(e.to_string(), vec!["A2".to_string()]);
        }

        let mut incoming = HashMap::new();
        for e in ["E_ref", "E_alt1", "E_alt2"] {
            incoming.insert(e.to_string(), vec!["A1".to_string()]);
        }
        incoming.insert(
            "A2".to_string(),
            vec!["E_ref".to_string(), "E_alt1".to_string(), "E_alt2".to_string()],
        );

        GraphicalGenome { anchor, edges, outgoing, incoming }
    }

    #[test]
    fn get_variant_keys_read_records_by_the_left_aligned_name() {
        // `construct_matrix` looks rows up by `generate_variant_name` of the returned variant.
        let mut graph = homopolymer_bubble_graph();
        let (variants, _coverage, read_record, cover_record) =
            get_variant(&mut graph, 3, "chrM_ref", HOMOPOLYMER_REF.len(), Some(HOMOPOLYMER_REF));

        let names: HashSet<String> = variants.iter().map(generate_variant_name).collect();
        assert_eq!(
            names,
            HashSet::from(["m.2CA>C".to_string()]),
            "both spellings must canonicalize to the same aligned name"
        );
        for key in read_record.keys() {
            assert!(names.contains(key), "read_record key {key} has no matching variant");
        }
        for key in cover_record.keys() {
            assert!(names.contains(key), "cover_record key {key} has no matching variant");
        }

        // Merging the spellings also merges their evidence: all seven alt reads land under the one canonical key.
        let mut alt_reads: Vec<String> = read_record["m.2CA>C"]
            .iter()
            .filter_map(|r| r.as_str().map(str::to_string))
            .collect();
        alt_reads.sort();
        assert_eq!(alt_reads, vec!["r1", "r2", "r3", "r4", "r5", "r6", "r7"]);
    }

    #[test]
    fn get_variant_without_left_align_keeps_both_raw_spellings() {
        // The flag is opt-in: with `None` the two spellings stay separate records at their original anchors.
        let mut graph = homopolymer_bubble_graph();
        let (variants, _coverage, read_record, _cover) =
            get_variant(&mut graph, 3, "chrM_ref", HOMOPOLYMER_REF.len(), None);

        let mut names: Vec<String> = variants.iter().map(generate_variant_name).collect();
        names.sort();
        assert_eq!(names, vec!["m.3AA>A".to_string(), "m.5AA>A".to_string()]);
        for key in read_record.keys() {
            assert!(names.contains(key), "read_record key {key} has no matching variant");
        }
    }

    /// rCRS 560..=580 (1-based), covering the poly-C run at 568-573.
    /// Slice index `i` is 1-based reference coordinate `RCRS_OFFSET + i`.
    const RCRS_560_580: &[u8] = b"CAAAGACACCCCCCACAGTTT";
    const RCRS_OFFSET: usize = 560;

    #[test]
    fn left_align_matches_bcftools_norm_on_the_rcrs_poly_c_run() {
        // Cross-checked against `bcftools norm -f rCRS.fasta`.
        let to_vcf_pos = |slice_pos: usize| slice_pos + RCRS_OFFSET - 1;
        // Slice index of the last C of the poly-C run (1-based 573).
        let last_c = 573 - RCRS_OFFSET;

        // Deletion of that last C: emitted as `572 CC>C` before alignment.
        let del = del_at(RCRS_560_580, last_c, 1, 4);
        assert_eq!(to_vcf_pos(del.pos), 572, "sanity check on the input spelling");
        let aligned = left_align_variants(&[del], RCRS_560_580);
        assert_eq!(to_vcf_pos(aligned[0].pos), 567);
        assert_eq!(aligned[0].ref_allele, "AC");
        assert_eq!(aligned[0].alt_allele, "A");

        // Insertion of a further C after that last C: emitted as `573 C>CC`.
        let ins = ins_at(RCRS_560_580, last_c + 1, "C", 4);
        assert_eq!(to_vcf_pos(ins.pos), 573, "sanity check on the input spelling");
        let aligned = left_align_variants(&[ins], RCRS_560_580);
        assert_eq!(to_vcf_pos(aligned[0].pos), 567);
        assert_eq!(aligned[0].ref_allele, "A");
        assert_eq!(aligned[0].alt_allele, "AC");
    }

    fn write_vcf_to_temp(name: &str, left_align: bool) -> String {
        let path = std::env::temp_dir().join(name);
        let variants = vec![del_at(HOMOPOLYMER_REF, 5, 1, 3)];
        let coverage: HashMap<usize, usize> = [(5, 30)].into_iter().collect();
        write_vcf(&variants, &coverage, &path, "s1", "chrM", 10, 0.1, left_align, None).unwrap();
        std::fs::read_to_string(&path).unwrap()
    }

    #[test]
    fn vcf_header_records_that_indels_were_left_aligned() {
        // Whether a VCF is normalized changes how it may be compared or merged.
        let vcf = write_vcf_to_temp("himito_left_align_on.vcf", true);
        assert!(
            vcf.contains("##HimitoNormalization=left-aligned"),
            "header must declare the normalization; got:\n{vcf}"
        );
    }

    #[test]
    fn vcf_header_omits_the_normalization_line_when_not_left_aligning() {
        let vcf = write_vcf_to_temp("himito_left_align_off.vcf", false);
        assert!(!vcf.contains("HimitoNormalization"));
    }

    #[test]
    fn resolve_strand_bias_action_drops_only_for_denoised_ont() {
        assert_eq!(resolve_strand_bias_action("ont-denoised", None), StrandBiasAction::Drop);
        assert_eq!(resolve_strand_bias_action("pacbio", None), StrandBiasAction::Flag);
        assert_eq!(resolve_strand_bias_action("ont-r10", None), StrandBiasAction::Flag);
        // An explicit CLI choice always wins.
        assert_eq!(
            resolve_strand_bias_action("ont-denoised", Some(StrandBiasAction::Flag)),
            StrandBiasAction::Flag
        );
    }

    // ===================================================================

    #[test]
    fn allele_span_snp_uses_pos_as_start() {
        let v = snp(100, "A", 5);
        assert_eq!(allele_span(&v), (100, "T", "A"));
        // Invariant the Python oracle relies on: `m.<start0 + 1>`.
        assert!(generate_variant_name(&v).starts_with(&format!("m.{}", 100 + 1)));
    }

    #[test]
    fn allele_span_del_starts_at_anchor() {
        let v = del_at(HOMOPOLYMER_REF, 5, 1, 3);
        let (start0, r, a) = allele_span(&v);
        assert_eq!((start0, r, a), (4, "AA", "A"));
        assert!(generate_variant_name(&v).starts_with(&format!("m.{}", start0 + 1)));
    }

    #[test]
    fn allele_span_ins_starts_at_anchor() {
        let v = ins_at(HOMOPOLYMER_REF, 6, "A", 3);
        let (start0, r, a) = allele_span(&v);
        assert_eq!((start0, r, a), (5, "A", "AA"));
        assert!(generate_variant_name(&v).starts_with(&format!("m.{}", start0 + 1)));
    }

    #[test]
    fn genotype_window_in_unique_sequence_is_flank_plus_two() {
        // 40 distinct-neighbour bases: no homopolymer anywhere.
        let refseq: Vec<u8> = b"ACGT".iter().cycle().take(40).copied().collect();
        // SNP at 20: span (20, 21), flank 3 -> [17, 24) then +-2 -> (15, 26).
        assert_eq!(genotype_window(&refseq, (20, 21), 3), Some((15, 26)));
    }

    #[test]
    fn genotype_window_pushes_out_of_homopolymer_then_pads() {
        // DEL span, flank 1: left starts at 3 inside the A-run, is pushed to 2.
        let v = del_at(HOMOPOLYMER_REF, 5, 1, 3);
        let (start0, r, _) = allele_span(&v);
        assert_eq!(
            genotype_window(HOMOPOLYMER_REF, (start0, start0 + r.len()), 1),
            Some((0, 9))
        );
    }

    #[test]
    fn genotype_window_outside_reference_is_none() {
        let refseq: Vec<u8> = b"ACGT".iter().cycle().take(40).copied().collect();
        assert_eq!(genotype_window(&refseq, (3, 4), 10), None);
        assert_eq!(genotype_window(&refseq, (36, 37), 10), None);
    }

    fn sites_as_names(vs: &[Variant], sites: &[Vec<usize>]) -> Vec<Vec<String>> {
        sites
            .iter()
            .map(|s| s.iter().map(|&i| generate_variant_name(&vs[i])).collect())
            .collect()
    }

    #[test]
    fn group_sites_merges_snps_at_the_same_base() {
        let vs = vec![snp(100, "A", 5), snp(100, "C", 4)];
        let sites = group_sites(&vs);
        assert_eq!(sites, vec![vec![0, 1]]);
    }

    #[test]
    fn group_sites_merges_overlapping_indels() {
        // m.302A>AC and m.302A>ACC: both anchored on base 302.
        let a = Variant {
            pos: 302, ref_allele: "A".into(), alt_allele: "AC".into(),
            variant_type: "INS".into(), allele_count: 3, filter: None,
        };
        let b = Variant { alt_allele: "ACC".into(), ..a.clone() };
        let sites = group_sites(&[a, b]);
        assert_eq!(sites, vec![vec![0, 1]]);
    }

    #[test]
    fn group_sites_never_groups_a_snp_with_an_overlapping_indel() {
        // The real ONT case: m.8170AT>A is the deletion ERROR of the very T in m.8171T>A.
        let del = Variant {
            pos: 8170, ref_allele: "AT".into(), alt_allele: "A".into(),
            variant_type: "DEL".into(), allele_count: 9, filter: None,
        };
        let sub = Variant {
            pos: 8170, ref_allele: "T".into(), alt_allele: "A".into(),
            variant_type: "SNP".into(), allele_count: 7, filter: None,
        };
        let vs = vec![del, sub];
        assert_eq!(sites_as_names(&vs, &[vec![0], vec![1]]), vec![vec!["m.8170AT>A"], vec!["m.8171T>A"]]);
        let sites = group_sites(&vs);
        assert_eq!(sites.len(), 2, "SNP and indel must be separate sites: {sites:?}");
        assert!(sites.contains(&vec![0]) && sites.contains(&vec![1]));
    }

    #[test]
    fn group_sites_keeps_adjacent_non_overlapping_variants_apart() {
        let vs = vec![snp(100, "A", 5), snp(101, "A", 5)];
        assert_eq!(group_sites(&vs), vec![vec![0], vec![1]]);
        // Two indels whose REF spans merely touch are also separate.
        let d1 = Variant {
            pos: 10, ref_allele: "CA".into(), alt_allele: "C".into(),
            variant_type: "DEL".into(), allele_count: 3, filter: None,
        }; // span [9, 11)
        let d2 = Variant {
            pos: 12, ref_allele: "GT".into(), alt_allele: "G".into(),
            variant_type: "DEL".into(), allele_count: 3, filter: None,
        }; // span [11, 13)
        assert_eq!(group_sites(&[d1, d2]), vec![vec![0], vec![1]]);
    }

    const SCORE_REF: &[u8] = b"GATCCTAGGCATTGCAAGTCCGATAC";

    /// (REF, ALT) haplotypes: ALT carries a SNP at index 12 (T -> A).
    fn score_haps() -> Vec<(String, Vec<u8>)> {
        let mut alt = SCORE_REF.to_vec();
        assert_eq!(alt[12], b'T');
        alt[12] = b'A';
        vec![
            ("REF".to_string(), SCORE_REF.to_vec()),
            ("m.13T>A".to_string(), alt),
        ]
    }

    fn with_junk(core: &[u8]) -> Vec<u8> {
        let mut v = b"GGAC".to_vec();
        v.extend_from_slice(core);
        v.extend_from_slice(b"CATT");
        v
    }

    #[test]
    fn score_segment_exact_ref_with_junk_ends_is_ref_only() {
        let haps = score_haps();
        let seg = with_junk(SCORE_REF);
        assert_eq!(
            score_segment(&haps, &seg, 1, u32::MAX),
            SegCall::Best([0usize].into_iter().collect())
        );
    }

    #[test]
    fn score_segment_alt_segment_is_alt_only() {
        let haps = score_haps();
        let seg = with_junk(&haps[1].1);
        assert_eq!(
            score_segment(&haps, &seg, 1, u32::MAX),
            SegCall::Best([1usize].into_iter().collect())
        );
    }

    #[test]
    fn score_segment_deleted_snp_base_ties_ref_and_alt() {
        let haps = score_haps();
        let mut core = SCORE_REF.to_vec();
        core.remove(12); // the SNP base is deleted from the read
        let seg = with_junk(&core);
        // A tie must come back as both indices: NaN downstream, not 0.
        assert_eq!(
            score_segment(&haps, &seg, 1, u32::MAX),
            SegCall::Best([0usize, 1].into_iter().collect())
        );
    }

    #[test]
    fn score_segment_gates_on_max_ed() {
        let haps = score_haps();
        // 6 substitutions away from REF (and 6 or 7 from ALT).
        let mut core = SCORE_REF.to_vec();
        for i in [0usize, 3, 6, 9, 15, 18] {
            core[i] = if core[i] == b'A' { b'C' } else { b'A' };
        }
        let seg = with_junk(&core);
        assert_eq!(score_segment(&haps, &seg, 1, 4), SegCall::TooNoisy);
        // The same read is decisive-or-tied, not noisy, when the gate is off.
        assert!(matches!(score_segment(&haps, &seg, 1, u32::MAX), SegCall::Best(_)));
    }

    #[test]
    fn score_segment_margin_widens_the_best_set() {
        let haps = score_haps();
        let seg = with_junk(SCORE_REF);
        // ALT is 1 edit worse than REF; margin 2 keeps it in the best set.
        assert_eq!(
            score_segment(&haps, &seg, 2, u32::MAX),
            SegCall::Best([0usize, 1].into_iter().collect())
        );
    }

    const GENO_REF: &str = "GCTAGACATACATACATACACGTCAGCACGACTGTGCAGTGTGATCGCTAGTAGTAGTGT";

    fn geno_snp() -> Variant {
        // GENO_REF[30] == 'A'; m.31A>T
        Variant {
            pos: 30, ref_allele: "A".into(), alt_allele: "T".into(),
            variant_type: "SNP".into(), allele_count: 5, filter: None,
        }
    }

    fn geno_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_data/genotype").join(name)
    }

    fn run_geno(
        sam: &str,
        variants: &[Variant],
    ) -> (HashMap<String, Vec<serde_json::Value>>, HashMap<String, HashSet<String>>) {
        genotype_from_bam(
            &geno_path(sam), GENO_REF.as_bytes(), "chrM", variants, &GenotypeOpts::default(),
        )
    }

    fn alt_reads(rr: &HashMap<String, Vec<serde_json::Value>>, name: &str) -> HashSet<String> {
        rr.get(name).map(|v| v.iter().map(|x| x.as_str().unwrap().to_string()).collect()).unwrap_or_default()
    }

    #[test]
    fn genotype_from_bam_calls_ref_alt_and_leaves_ambiguous_reads_missing() {
        let v = geno_snp();
        let name = generate_variant_name(&v);
        assert_eq!(name, "m.31A>T");
        let (rr, cr) = run_geno("reads.sam", &[v]);
        let alt = alt_reads(&rr, &name);
        let cover = &cr[&name];

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
        let v = geno_snp();
        let name = generate_variant_name(&v);
        let (rr, cr) = run_geno("split_disagree.sam", &[v]);
        assert!(!cr[&name].contains("r_split"));
        assert!(!alt_reads(&rr, &name).contains("r_split"));
        assert!(cr[&name].contains("chrM"));
    }

    #[test]
    fn genotype_from_bam_scores_an_indel_site_and_feeds_construct_matrix() {
        let refb = GENO_REF.as_bytes();
        let snp_v = geno_snp();
        let del = del_at(refb, 30, 1, 4); // m.30CA>C ; r_del carries it
        let dname = generate_variant_name(&del);
        let sname = generate_variant_name(&snp_v);
        let variants = vec![snp_v, del];
        let (rr, cr) = run_geno("reads.sam", &variants);

        assert!(alt_reads(&rr, &dname).contains("r_del"));
        assert!(cr[&dname].contains("r_ref") && !alt_reads(&rr, &dname).contains("r_ref"));

        // Same shapes construct_matrix already consumes.
        let (m, vars, reads) = construct_matrix(&rr, &cr, &variants);
        let col = |n: &str| reads.iter().position(|r| r == n).unwrap();
        let row = |n: &str| vars.iter().position(|v| generate_variant_name(v) == n).unwrap();
        assert_eq!(m[[row(&sname), col("r_alt")]], 1.0);
        assert_eq!(m[[row(&sname), col("r_ref")]], 0.0);
        assert!(m[[row(&sname), col("r_del")]].is_nan());
        assert_eq!(m[[row(&dname), col("r_del")]], 1.0);
        assert_eq!(m[[row(&sname), col("chrM")]], 0.0);
        assert_eq!(m[[row(&dname), col("chrM")]], 0.0);
    }

    #[test]
    fn genotype_from_bam_variant_without_a_window_is_only_the_ref_column() {
        // Window would cross the origin: no read can call it.
        let v = Variant {
            pos: 3, ref_allele: "A".into(), alt_allele: "T".into(),
            variant_type: "SNP".into(), allele_count: 5, filter: None,
        };
        let name = generate_variant_name(&v);
        let (rr, cr) = run_geno("reads.sam", &[v]);
        assert!(alt_reads(&rr, &name).is_empty());
        assert_eq!(cr[&name], ["chrM".to_string()].into_iter().collect());
    }

    #[test]
    fn genotype_from_bam_matches_ref_name_when_chrm_is_not_tid_zero() {
        let v = geno_snp();
        let name = generate_variant_name(&v);
        let (rr, cr) = run_geno("reads_second_contig.sam", &[v.clone()]);
        let alt = alt_reads(&rr, &name);
        let cover = &cr[&name];
        assert!(cover.contains("r_ref") && !alt.contains("r_ref"));
        assert!(cover.contains("r_alt") && alt.contains("r_alt"));
        assert!(cover.contains("r_split") && alt.contains("r_split"));
        assert!(cover.contains("chrM") && !alt.contains("chrM"));
        // Same decisive calls as the tid-0 fixture.
        let (rr0, cr0) = run_geno("reads.sam", &[v]);
        assert_eq!(alt_reads(&rr, &name), alt_reads(&rr0, &name));
        assert_eq!(cr[&name], cr0[&name]);
    }

    #[test]
    fn genotype_from_bam_unknown_ref_name_returns_ref_column_only() {
        let v = geno_snp();
        let name = generate_variant_name(&v);
        let (rr, cr) = genotype_from_bam(
            &geno_path("reads_second_contig.sam"),
            GENO_REF.as_bytes(),
            "not_a_contig",
            &[v],
            &GenotypeOpts::default(),
        );
        assert!(alt_reads(&rr, &name).is_empty());
        assert_eq!(cr[&name], ["not_a_contig".to_string()].into_iter().collect());
    }

    fn tmp_out(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("himito_geno_{}_{}", tag, std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d.join("out.vcf")
    }

    /// A stand-in for the graph matrix `start` has in hand after the permutation test.
    fn graph_matrix_fixture(v: &Variant) -> (Array2<f64>, Vec<Variant>, Vec<String>) {
        let m = Array2::from_shape_vec((1, 2), vec![1.0, f64::NAN]).unwrap();
        (m, vec![v.clone()], vec!["g1".to_string(), "g2".to_string()])
    }

    #[test]
    fn write_final_matrix_without_genotype_bam_is_unchanged() {
        let out = tmp_out("none");
        let v = geno_snp();
        let (m, vars, reads) = graph_matrix_fixture(&v);
        let source = write_final_matrix(&out, &m, &vars, &[], &reads, None).unwrap();

        assert_eq!(source, MatrixSource::Graph);
        assert!(GRAPH_MATRIX_WARNING.contains("--bam"), "{GRAPH_MATRIX_WARNING}");
        let text = std::fs::read_to_string(out.with_extension("matrix.csv")).unwrap();
        assert_eq!(text, "variant,g1,g2\nm.31A>T,1,\n");
        assert!(!out.with_extension("graph_matrix.csv").exists());
    }

    #[test]
    fn write_final_matrix_with_bam_writes_no_graph_copy() {
        let out = tmp_out("swap");
        let v = geno_snp();
        let (m, vars, reads) = graph_matrix_fixture(&v);
        let bam = geno_path("reads.sam");
        let src = GenotypeSource {
            bam: &bam, refseq: GENO_REF.as_bytes(), ref_name: "chrM", opts: GenotypeOpts::default(),
        };
        let source = write_final_matrix(&out, &m, &vars, &[], &reads, Some(&src)).unwrap();

        // matrix.csv is BAM-genotyped: reads sorted by name, ref pseudo-read all 0, r_del / r_short missing.
        assert_eq!(source, MatrixSource::Bam);
        let text = std::fs::read_to_string(out.with_extension("matrix.csv")).unwrap();
        assert_eq!(text, "variant,chrM,r_alt,r_ref,r_split\nm.31A>T,0,1,0,1\n");
        assert!(!out.with_extension("graph_matrix.csv").exists());
    }

    #[test]
    fn resolve_genotype_opts_gates_ont_noise_only() {
        for dt in ["ont-r9", "ont-r10", "ont-denoised"] {
            let o = resolve_genotype_opts(dt);
            assert_eq!(o.max_ed, 4, "{dt}");
            assert_eq!((o.flank, o.margin), (10, 1), "{dt}");
        }
        let o = resolve_genotype_opts("pacbio");
        assert_eq!(o.max_ed, u32::MAX);
        assert_eq!((o.flank, o.margin), (10, 1));
    }

    #[test]
    fn linkage_rescue_is_on_for_ont_only() {
        for dt in ["ont-r9", "ont-r10", "ont-denoised"] {
            assert!(linkage_rescue_enabled(dt, false), "{dt}");
            assert!(!linkage_rescue_enabled(dt, true), "{dt}");
        }
        assert!(!linkage_rescue_enabled("pacbio", false));
    }

    #[test]
    fn indel_false_threshold_preset_by_data_type() {
        for dt in ["ont-r9", "ont-r10", "ont-denoised"] {
            assert_eq!(resolve_indel_false_threshold(dt, None), 0.015, "{dt}");
        }
        assert_eq!(resolve_indel_false_threshold("pacbio", None), 0.1);
        // Explicit CLI value always wins, for every data type.
        for dt in ["ont-r9", "ont-r10", "ont-denoised", "pacbio"] {
            assert_eq!(resolve_indel_false_threshold(dt, Some(0.05)), 0.05, "{dt}");
        }
    }

    #[test]
    fn small_indel_at_ont_threshold_passes_above_and_is_artifact_below() {
        let del = indel(100, "GA", "G", 18);
        let cov: HashMap<usize, usize> = [(100usize, 1000usize)].into_iter().collect();
        let rec = format_vcf_record(&del, cov.clone(), 0.015);
        assert!(rec.contains("\tPASS\t"), "HF 0.018 should pass at 0.015: {rec}");

        let low = indel(100, "GA", "G", 11);
        let rec = format_vcf_record(&low, cov, 0.015);
        assert!(rec.contains("Potential_Artifact"), "HF 0.011 should be flagged: {rec}");
    }
}
