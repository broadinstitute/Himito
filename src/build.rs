// Import necessary standard library modules
use std::collections::HashSet;
use std::path::PathBuf;
use bio::io::fasta::{Reader, Record};
use rust_htslib::bam::{self, Read};
use indicatif::ProgressBar;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, Write};
use agg::*;
use crate::agg;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use rayon::prelude::*;
use rustc_hash::FxHashMap;


pub fn reverse_complement(kmer: &str) -> String {
    // Byte-level: ACGTN is pure ASCII, so this avoids UTF-8 char decoding/boundary
    // checks that `.chars()` pays on every element.
    let bytes = kmer.as_bytes();
    let mut out = vec![0u8; bytes.len()];
    for (dst, &b) in out.iter_mut().rev().zip(bytes.iter()) {
        *dst = match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'N' => b'N',
            _ => panic!("Unexpected character: {}", b as char),
        };
    }
    // Safety: every byte written above is a valid single-byte ASCII code point.
    String::from_utf8(out).unwrap()
}

/// 2-bit encoding for a single base; `None` for anything ambiguous (N, lowercase soft-mask, etc).
#[inline]
fn base_to_bits(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Encode a short (<=32bp) sequence into a 2-bit-packed u64. Returns `None` if it
/// contains any non-ACGT base or exceeds 32bp (u64 capacity).
fn encode_kmer(bytes: &[u8]) -> Option<u64> {
    if bytes.len() > 32 {
        return None;
    }
    let mut code: u64 = 0;
    for &b in bytes {
        code = (code << 2) | base_to_bits(b)?;
    }
    Some(code)
}

/// Hashable identity for a kmer window. `Bits` is the fast, allocation-free path for
/// pure-ACGT windows (the overwhelming majority in real sequence data). `Raw` is an
/// exact byte-slice fallback for anything that can't be 2-bit packed (an ambiguous
/// base like the single historical `N` at rCRS position 3107, or k > 32) so uniqueness
/// is decided by the same exact-match semantics as the original `String` comparison.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum KmerKey {
    Bits(u64),
    Raw(Box<[u8]>),
}

fn encode_kmer_key(bytes: &[u8]) -> KmerKey {
    match encode_kmer(bytes) {
        Some(code) => KmerKey::Bits(code),
        None => KmerKey::Raw(bytes.into()),
    }
}

/// Rolling k-mer scan over `contig`, calling `visit(key, pos)` for every window,
/// including ones spanning an ambiguous base. Runs in O(n) with zero per-position
/// allocation for the common all-ACGT case, only falling back to an owned byte slice
/// for the rare window that can't be 2-bit packed.
fn scan_kmers<F: FnMut(KmerKey, usize)>(contig: &[u8], k: usize, mut visit: F) {
    if k == 0 || contig.len() < k {
        return;
    }
    if k > 32 {
        // 2-bit packing needs <=32 bases per u64; fall back to exact byte comparison.
        for pos in 0..=contig.len() - k {
            visit(KmerKey::Raw(contig[pos..pos + k].into()), pos);
        }
        return;
    }

    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let mut code: u64 = 0;
    let mut run_len: usize = 0;

    for (i, &b) in contig.iter().enumerate() {
        match base_to_bits(b) {
            Some(bits) => {
                code = ((code << 2) | bits) & mask;
                run_len += 1;
            }
            None => {
                run_len = 0;
                code = 0;
            }
        }
        if i + 1 >= k {
            let pos = i + 1 - k;
            if run_len >= k {
                visit(KmerKey::Bits(code), pos);
            } else {
                visit(KmerKey::Raw(contig[pos..pos + k].into()), pos);
            }
        }
    }
}

pub fn map_to_genome(contig: &str, k: usize) -> FxHashMap<KmerKey, Vec<usize>> {
    let mut position_dict: FxHashMap<KmerKey, Vec<usize>> = FxHashMap::default();

    scan_kmers(contig.as_bytes(), k, |key, pos| {
        position_dict.entry(key).or_insert_with(Vec::new).push(pos);
    });

    position_dict
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnchorInfo {
    seq: String,
    pos: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeInfo {
    seq: String,
    reads: Vec<String>,
    samples: HashSet<String>,
    src: String,
    dst: String,
}

impl AnchorInfo {
    // Getter method for `seq`
    pub fn get_seq(&self) -> &String {
        &self.seq
    }
}

pub fn create_anchors(
    position_dict: &FxHashMap<KmerKey, Vec<usize>>,
    contig: &str,
    k: usize,
) -> HashMap<String, AnchorInfo> {
    let mut anchor_info = HashMap::new();

    // Kmers that appear exactly once become anchors; the sequence is recovered by
    // slicing the original contig at that position rather than storing it as a key.
    for positions in position_dict.values() {
        if positions.len() == 1 {
            let pos = positions[0];
            let anchor_name = format!("A{:06}", pos);
            anchor_info.insert(
                anchor_name,
                AnchorInfo {
                    seq: contig[pos..pos + k].to_string(),
                    pos,
                },
            );
        }
    }

    anchor_info
}

pub fn get_final_anchor(
    anchor_info: &HashMap<String, AnchorInfo>,
    k: usize,
) -> HashMap<String, AnchorInfo> {
    let mut final_anchor = HashMap::new();
    let mut anchornames: Vec<&String> = anchor_info.keys().collect();
    anchornames.sort();
    // println!("{:?}", anchornames);

    let mut anchor_unadjacent_list = Vec::new();
    let mut last_position = 0;

    for anchor in anchornames {
        let position = anchor_info[anchor].pos;
        // assert!(position > last_position);
        if position > last_position + k {
            anchor_unadjacent_list.push(anchor.clone());
            last_position = position;
        }
    }

    for anchor_name in &anchor_unadjacent_list {
        if let Some(anchor) = anchor_info.get(anchor_name) {
            final_anchor.insert(anchor_name.clone(), anchor.clone());
        }
    }

    final_anchor
}

pub fn mapping_info(
    anchor_info: &HashMap<String, &AnchorInfo>,
    contig: &str,
    k: usize,
) -> (HashMap<String, usize>, HashMap<String, Vec<usize>>) {
    // Seed the position dict with the encoded forward and reverse-complement kmer for
    // every anchor (anchor kmers are few, so encoding them isn't hot).
    let mut position_dict: FxHashMap<KmerKey, Vec<usize>> = FxHashMap::default();
    for info in anchor_info.values() {
        let anchor_seq = &info.seq;
        let anchor_rev = reverse_complement(anchor_seq);
        position_dict.entry(encode_kmer_key(anchor_seq.as_bytes())).or_insert_with(Vec::new);
        position_dict.entry(encode_kmer_key(anchor_rev.as_bytes())).or_insert_with(Vec::new);
    }

    // Single O(n) rolling scan of the contig instead of allocating a String per kmer.
    scan_kmers(contig.as_bytes(), k, |key, pos| {
        if let Some(v) = position_dict.get_mut(&key) {
            v.push(pos);
        }
    });

    let mut a = HashMap::new();
    let mut svs = HashMap::new();

    // Process each anchor
    for (anchor, info) in anchor_info {
        let anchor_seq = &info.seq;
        let anchor_rev = reverse_complement(anchor_seq);

        let empty = Vec::new();
        let fwd_positions = position_dict.get(&encode_kmer_key(anchor_seq.as_bytes())).unwrap_or(&empty);
        let rev_positions = position_dict.get(&encode_kmer_key(anchor_rev.as_bytes())).unwrap_or(&empty);

        // Combine positions from forward and reverse sequences
        let position_list = [fwd_positions.clone(), rev_positions.clone()].concat();
        // ! should change into position_dict.get(anchor_seq).len() == 1 and position_dict.get(&anchor_rev) == 0
        if position_list.len() == 1 {
            a.insert(anchor.clone(), position_list[0]);
        } else {
            svs.insert(anchor.clone(), position_list);
        }
    }

    (a, svs)
}


pub fn construct_edges(
    src_pos: usize,
    dst_pos: usize,
    k: usize,
    contig: &str,
    contigname: String,
    sample: String,
    anchorseq: &FxHashMap<String, String>,
) -> EdgeInfo {
    let src_seq: String;
    let mut pr = false;
    let mut src = "".to_string();

    if src_pos == 0 {
        src = "SOURCE".to_string();
        src_seq = "".to_string();
        pr = false;
    } else {
        src_seq = contig[src_pos..src_pos+k].to_string();
        src = match anchorseq.get(&src_seq) {
            Some(value) => value.clone(),
            None => {
                let reversed_src_seq = reverse_complement(&src_seq);
                pr = true;
                anchorseq
                    .get(&reversed_src_seq)
                    .unwrap_or(&"".to_string())
                    .clone()
            }
        };
    }

    let mut sr = false;
    let mut dst = "".to_string();
    let dst_seq;
    if dst_pos == contig.len() {
        dst = "SINK".to_string();
        dst_seq = "".to_string();
        sr = true;
    } else {
        dst_seq = contig[dst_pos..dst_pos+k].to_string();
        dst = match anchorseq.get(&dst_seq) {
            Some(value) => value.clone(),
            None => {
                let reversed_dst_seq = reverse_complement(&dst_seq);
                sr = true;
                anchorseq
                    .get(&reversed_dst_seq)
                    .unwrap_or(&"".to_string())
                    .clone()
            }
        };
    }

    let mut edge_seq = if src_pos == 0 {
        contig.get(0..dst_pos).unwrap_or_default().to_string()
    } else {
        contig
            .get(src_pos + k..dst_pos)
            .unwrap_or_default()
            .to_string()
    };

    if pr && sr {
        edge_seq = reverse_complement(&edge_seq);
        let node = src.clone();
        src = dst.clone();
        dst = node.clone();
    }

    EdgeInfo {
        seq: edge_seq,
        src,
        dst,
        reads: vec![contigname],
        samples: vec![sample].into_iter().collect(),
    }
}



/// Build the edges contributed by a single contig, with no shared/cross-contig state.
/// Kept separate so `create_edge_file` can run this per-contig work in parallel.
fn edges_for_contig(
    contig_name: &str,
    contig: &str,
    final_anchor: &HashMap<String, &AnchorInfo>,
    anchorseq: &FxHashMap<String, String>,
    k: usize,
    threshold: usize,
) -> Vec<EdgeInfo> {
    let mut local_edges = Vec::new();
    if contig.len() < k + 1 {
        return local_edges;
    }

    let sample_name = if contig_name.contains('|') {
        contig_name.split('|').last().unwrap_or("").to_string()
    } else {
        "".to_string()
    };

    let (a, _svs) = mapping_info(final_anchor, contig, k);
    if a.len() < threshold {
        return local_edges;
    }

    let mut splitposlist: Vec<_> = a.values().copied().collect();
    splitposlist.sort();
    let mut src_pos = 0;

    // Process all positions except the last
    for &dst_pos in &splitposlist {
        if dst_pos - src_pos < k + 1 {
            continue;
        }
        local_edges.push(construct_edges(
            src_pos,
            dst_pos,
            k,
            contig,
            contig_name.to_string(),
            sample_name.clone(),
            anchorseq,
        ));
        src_pos = dst_pos;
    }

    // Process final edge to end of contig
    let dst_pos = contig.len();
    if dst_pos - src_pos > k + 1 {
        local_edges.push(construct_edges(
            src_pos,
            dst_pos,
            k,
            contig,
            contig_name.to_string(),
            sample_name.clone(),
            anchorseq,
        ));
    }

    local_edges
}

pub fn create_edge_file(
    all_seq: &HashMap<String, String>,
    final_anchor: &HashMap<String, &AnchorInfo>,
    k: usize,
    threshold: usize
) -> (HashMap<String, EdgeInfo>, HashMap<String, Vec<String>>) {
    let anchorseq: FxHashMap<_, _> = final_anchor
        .iter()
        .map(|(anchor, info)| (info.seq.clone(), anchor.clone()))
        .collect();

    // Sort for deterministic iteration order (a plain HashMap iteration order isn't
    // stable, and we want reproducible output regardless of thread scheduling).
    let mut contigs: Vec<(&String, &String)> = all_seq.iter().collect();
    contigs.sort_by(|a, b| a.0.cmp(b.0));

    let bar = ProgressBar::new(contigs.len() as u64);

    // Each contig's edges are computed independently in parallel; no shared state is
    // touched until the sequential merge below.
    let per_contig_edges: Vec<Vec<EdgeInfo>> = contigs
        .par_iter()
        .map(|(contig_name, contig)| {
            let edges = edges_for_contig(contig_name, contig, final_anchor, &anchorseq, k, threshold);
            bar.inc(1);
            edges
        })
        .collect();
    bar.finish();

    // Sequential merge: dedup edges globally by (src, dst, seq) in O(1) amortized per
    // edge via a lookup map, instead of an O(m) linear scan of each src's edge list.
    let mut edge_info: HashMap<String, EdgeInfo> = HashMap::new();
    let mut outgoing: HashMap<String, Vec<String>> = HashMap::new();
    let mut edge_key_to_name: FxHashMap<(String, String, String), String> = FxHashMap::default();
    let mut edge_counter: usize = 0;

    for edges in per_contig_edges {
        for e in edges {
            let key = (e.src.clone(), e.dst.clone(), e.seq.clone());
            if let Some(existing_name) = edge_key_to_name.get(&key) {
                let existing_edge = edge_info.get_mut(existing_name).unwrap();
                existing_edge.reads.extend(e.reads);
                existing_edge.samples.extend(e.samples);
            } else {
                let edgename = format!("E{:08}", edge_counter);
                edge_counter += 1;
                outgoing.entry(e.src.clone()).or_default().push(edgename.clone());
                edge_key_to_name.insert(key, edgename.clone());
                edge_info.insert(edgename, e);
            }
        }
    }

    (edge_info, outgoing)
}

pub fn write_gfa(
    final_anchor: &HashMap<String, AnchorInfo>,
    edge_info: &HashMap<String, EdgeInfo>,
    output_filename: &str,
) -> std::result::Result<(), Box<dyn Error>> {
    let mut file = File::create(output_filename)?;

    writeln!(file, "H\tVN:Z:1.0")?;

    let mut anchor_output = Vec::new();
    let mut keys: Vec<_> = final_anchor.keys().collect();
    keys.sort();
    for anchor in keys.iter() {
        let info = &final_anchor[*anchor];
        let seq = &info.seq;
        let mut anchor_info_clone = HashMap::new();
        anchor_info_clone.insert("pos".to_string(), info.pos);
        // anchor_info_clone.seq = String::new();
        let json_string =
            serde_json::to_string(&anchor_info_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = format!("S\t{}\t{}\tPG:J:{}", anchor, seq, json_string);
        anchor_output.push(formatted_string);
        // writeln!(file,"S\t{}\t{}\tPG:J:{}\n", anchor, seq, json_string)?;
    }
    let mut edge_output = Vec::new();
    let mut link_output = Vec::new();

    let mut edge_keys: Vec<_> = edge_info.keys().collect();
    edge_keys.sort();
    for edge in edge_keys {
        let edge_data = &edge_info[edge];
        let seq = &edge_data.seq;
        let src = &edge_data.src;
        let dst = &edge_data.dst;
        let mut edge_data_clone = HashMap::new();
        edge_data_clone.insert("reads".to_string(), &edge_data.reads);
        let sample_vec: Vec<String> = edge_data.samples.iter().cloned().collect();
        let src_vec = vec![edge_data.src.clone()];
        let dst_vec = vec![edge_data.dst.clone()];
        edge_data_clone.insert("samples".to_string(), &sample_vec);
        edge_data_clone.insert("src".to_string(), &src_vec);
        edge_data_clone.insert("dst".to_string(), &dst_vec);

        let json_string =
            serde_json::to_string(&edge_data_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = if !edge_data.reads.is_empty() {
            format!(
                "S\t{}\t{}\tPG:J:{}\tRC:i:{}",
                edge,
                seq,
                json_string,
                edge_data.reads.len()
            )
        } else {
            format!(
                "S\t{}\t{}\tPG:J:{}\tRC:i:{}",
                edge,
                seq,
                json_string,
                edge_data.reads.len()
            )
        };
        edge_output.push(formatted_string);
        
        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", src, edge));
        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", edge, dst));

    }

    for s in anchor_output {
        writeln!(file, "{}", s)?;
    }
    for s in edge_output {
        writeln!(file, "{}", s)?;
    }
    for l in link_output {
        writeln!(file, "{}", l)?;
    }
    Ok(())
}

pub fn write_anchor_json(final_anchor: &HashMap<String, AnchorInfo>, output_path: &str) -> io::Result<()> {
    let json = serde_json::to_string_pretty(final_anchor)?;
    let mut file = File::create(output_path)?;
    file.write_all(json.as_bytes())?;
    Ok(())
}



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
    // Byte-level: avoids `.chars()` UTF-8 decoding since ACGTN is pure ASCII.
    let out: Vec<u8> = seq
        .as_bytes()
        .iter()
        .map(|&b| {
            let upper = b.to_ascii_uppercase();
            if matches!(upper, b'A' | b'G' | b'C' | b'T') {
                upper
            } else {
                b.to_ascii_lowercase()
            }
        })
        .collect();
    // Safety: every output byte is a valid single-byte ASCII code point.
    String::from_utf8(out).unwrap()
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
    // match=0, mismatch=-6, gap_open=-5, gap_extend=-3.
    //
    // The key relation is mismatch (-6) being *dearer* than a gap open (-5).
    // Inside a homopolymer/tandem-repeat run that is what lets the correct
    // two-gap edit outscore the one-gap-plus-mismatch shortcut: at chrM:310
    // (a lone T flanked by poly-C) representing the extra C's as 3I..1I costs
    // -22, while folding one of them into a T>C substitution costs -23. Under
    // the previous match=+1/mismatch=-1/gap_extend=-1 settings the shortcut
    // won, and chrM:310 was emitted as a false near-homoplasmic T>C "SNP"
    // alongside a same-run C insertion.
    let score = |a: u8, b: u8| if a == b { 0i32 } else { -6i32 };
    let mut aligner = Aligner::with_capacity(sequence.len(), reference.len(), -5, -3, &score);

    // Perform the alignment
    let alignment = aligner.global(sequence.as_bytes(), reference.as_bytes());

    // Get the aligned sequences
    let cigar = alignment_to_cigar(&alignment.operations);
    // println!("{:?}", cigar);

    // Safety net. The scoring above already prefers the clean two-gap edit in
    // a homopolymer run, so the chrM:310-style artifact no longer reaches here.
    // It still fires where a spurious substitution survives scoring - a length-1
    // X against an I/D run whose orphaned base occurs inside that run can always
    // be re-anchored losslessly, and is not a real substitution.
    normalize_homopolymer_indels(&cigar, reference, sequence)
}

/// Rewrites a CIGAR to eliminate a length-1 `X` sitting directly against an
/// `I` or `D` run whenever the X's "orphaned" base (the one with no partner
/// on the gapped side) also occurs inside that run. In that case the run can
/// be losslessly re-anchored so the X becomes a plain match and the indel is
/// redistributed around it: same total ref/alt content, no substitution.
/// Left untouched when the orphaned base doesn't occur in the adjacent run -
/// that's very likely a genuine substitution, not a homopolymer-registration
/// artifact.
fn normalize_homopolymer_indels(cigar: &str, ref_seq: &str, alt_seq: &str) -> String {
    let ops = parse_cigar(cigar);
    let ref_bytes = ref_seq.as_bytes();
    let alt_bytes = alt_seq.as_bytes();

    let mut result: Vec<(usize, char)> = Vec::new();
    let mut ref_pos = 0usize;
    let mut alt_pos = 0usize;
    let mut i = 0usize;

    let push_op = |result: &mut Vec<(usize, char)>, len: usize, op: char| {
        if len == 0 {
            return;
        }
        if let Some(last) = result.last_mut() {
            if last.1 == op {
                last.0 += len;
                return;
            }
        }
        result.push((len, op));
    };

    while i < ops.len() {
        let (len, op) = ops[i];

        // I(n) immediately followed by a single mismatch: the mismatch's ref
        // base may really belong inside the insertion run.
        if op == 'I' && i + 1 < ops.len() && ops[i + 1] == (1, 'X') {
            let ins = &alt_bytes[alt_pos..alt_pos + len];
            let ref_x = ref_bytes[ref_pos];
            if let Some(p) = ins.iter().rposition(|&b| b == ref_x) {
                push_op(&mut result, p, 'I');
                push_op(&mut result, 1, '=');
                push_op(&mut result, len - p, 'I'); // trailing insertion piece + the X's alt base
                ref_pos += 1;
                alt_pos += len + 1;
                i += 2;
                continue;
            }
        }

        // A single mismatch immediately followed by I(n): mirror of the above.
        if op == 'X' && len == 1 && i + 1 < ops.len() && ops[i + 1].1 == 'I' {
            let n = ops[i + 1].0;
            let ins = &alt_bytes[alt_pos + 1..alt_pos + 1 + n];
            let ref_x = ref_bytes[ref_pos];
            if let Some(p) = ins.iter().position(|&b| b == ref_x) {
                push_op(&mut result, 1 + p, 'I'); // the X's alt base + leading insertion piece
                push_op(&mut result, 1, '=');
                push_op(&mut result, n - p - 1, 'I');
                ref_pos += 1;
                alt_pos += 1 + n;
                i += 2;
                continue;
            }
        }

        // D(n) immediately followed by a single mismatch: mirror of I+X with
        // ref/alt roles swapped.
        if op == 'D' && i + 1 < ops.len() && ops[i + 1] == (1, 'X') {
            let del = &ref_bytes[ref_pos..ref_pos + len];
            let alt_x = alt_bytes[alt_pos];
            if let Some(q) = del.iter().rposition(|&b| b == alt_x) {
                push_op(&mut result, q, 'D');
                push_op(&mut result, 1, '=');
                push_op(&mut result, len - q, 'D');
                ref_pos += len + 1;
                alt_pos += 1;
                i += 2;
                continue;
            }
        }

        // A single mismatch immediately followed by D(n): mirror of X+I.
        if op == 'X' && len == 1 && i + 1 < ops.len() && ops[i + 1].1 == 'D' {
            let n = ops[i + 1].0;
            let del = &ref_bytes[ref_pos + 1..ref_pos + 1 + n];
            let alt_x = alt_bytes[alt_pos];
            if let Some(q) = del.iter().position(|&b| b == alt_x) {
                push_op(&mut result, 1 + q, 'D');
                push_op(&mut result, 1, '=');
                push_op(&mut result, n - q - 1, 'D');
                ref_pos += 1 + n;
                alt_pos += 1;
                i += 2;
                continue;
            }
        }

        match op {
            '=' | 'X' => {
                ref_pos += len;
                alt_pos += len;
            }
            'I' => alt_pos += len,
            'D' => ref_pos += len,
            _ => {}
        }
        push_op(&mut result, len, op);
        i += 1;
    }

    result
        .iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect()
}

/// Parses a run-length CIGAR string ("23=4I1X20=") into `(length, op)` pairs.
fn parse_cigar(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut num = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            ops.push((num.parse().expect("cigar length"), c));
            num.clear();
        }
    }
    ops
}

pub fn generate_cigar(
    graph: &mut GraphicalGenome,
    ref_strain: &str,
    k: usize,
    maxlength: usize,
    minimal_read_count: usize,
) {
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

            // Inclusive gate: an edge needs at least `minimal_read_count` reads.
            // This was `<=` (i.e. strictly more than N), which with the default of
            // 1 meant >=2 reads and left 108k of 109k edges on a real ONT graph
            // un-CIGARed. Reads reaching a bubble only through those edges are not
            // credited as covering it, so their genotypes became NaN rather than a
            // ref call -- ~92% of depth at a given site.
            if allele_count < minimal_read_count {
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
            // std::thread::sleep(std::time::Duration::from_millis(50));

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
}


/// Read support for the edge that circularizes the reference path.
///
/// The mitochondrial genome is circular, but the graph is built linearly, so the
/// reference enters at SOURCE and leaves at SINK. Circularizing splices the
/// SINK-bound tail edge and the SOURCE-bound head edge into one edge whose
/// sequence is `tail.seq + head.seq`. A read only supports that spliced edge if
/// it traverses *both* halves - i.e. it actually wraps the origin - so the
/// support set is the intersection, not the union: a read that stops at the
/// linear end never enters the head portion.
///
/// Duplicates within `head` collapse, matching the `HashSet` semantics the
/// caller already relies on for the tail half.
fn circularized_edge_reads(tail: &HashSet<String>, head: &[String]) -> HashSet<String> {
    let head_set: HashSet<&String> = head.iter().collect();
    tail.iter()
        .filter(|r| head_set.contains(r))
        .cloned()
        .collect()
}

pub fn start(output: &PathBuf, k: usize, read_path: &PathBuf, reference_path: &PathBuf, maxlength: usize, min_edge_reads: usize) {
    // Read reference records into a vector
    let ref_reader = Reader::from_file(reference_path).unwrap();
    let reference_sequence: Vec<Record> = ref_reader.records().map(|r| r.unwrap()).collect();
    let ref_seq = String::from_utf8_lossy(reference_sequence[0].seq()).to_string();
    let ref_header = reference_sequence[0].id().to_string();
    
    let position_dict = map_to_genome(&ref_seq, k);
    let anchors = create_anchors(&position_dict, &ref_seq, k);
    println!("anchors, {:?}", anchors.len());

    let unadjacent_anchor = get_final_anchor(&anchors, k);
    println!("unadjacent_anchor, {:?}", unadjacent_anchor.len());

    // Read the reads records (name and sequence) into a vector.
    let mut read_dictionary = HashMap::new();
    match read_path.extension().and_then(|ext| ext.to_str()) {
        Some("fasta") | Some("fa") | Some("fna") => {
            // Handle FASTA files
            println!("Process FASTA file");
            let reader = Reader::from_file(read_path).unwrap();
            let all_reads: Vec<Record> = reader.records().map(|r| r.unwrap()).collect();
            for record in all_reads {
                let contig_name = record.id().to_string();
                let contig = String::from_utf8_lossy(record.seq()).to_string();
                read_dictionary.insert(contig_name, contig);
            }
        }
        Some("bam") => {
            // Handle BAM files
            println!("Processing BAM file");
            let mut bam = bam::Reader::from_path(read_path).unwrap();
            for record in bam.records(){
                let r = record.expect("Failed to read BAM record");
                let contig = String::from_utf8_lossy(&r.seq().as_bytes()).to_string();
                let contig_name = String::from_utf8_lossy(&r.qname()).to_string();
                read_dictionary.insert(contig_name, contig);
            }
        },
        _ => {
            // Handle unknown file types
            println!("Unsupported file format. Please provide a FASTA or BAM file.");
        }
    }
    // insert reference sequence to the read dictionary
    let reference_id = String::from_utf8_lossy(reference_sequence[0].id().as_bytes()).to_string();
    read_dictionary.insert(reference_id.clone(), ref_seq);

    // construct edge information
    let dereferenced_anchor: HashMap<String, &AnchorInfo> = unadjacent_anchor
        .iter()
        .map(|(k, v)| (k.clone(), v)) 
        .collect();
    let (n_edge_info, _outgoing) = create_edge_file(&read_dictionary, &dereferenced_anchor, k, 1);
    let mut edge_info = n_edge_info;
    // find final anchor set in the src and dst of edges
    let mut final_anchor_list = HashSet::new();
    
    let mut f_ref_edge = String::new();
    let mut f_dst_anchor = String::new();
    let mut l_ref_edge = String::new();
    let mut f_src_anchor = String::new();
    for (edgename, edge_info) in edge_info.iter(){
        let src = edge_info.src.clone();
        let dst = edge_info.dst.clone();
        let reads = edge_info.reads.clone();
        if reads.contains(&reference_id) {
            if src == "SOURCE".to_string() {
                f_ref_edge = edgename.clone();
                f_dst_anchor = dst.clone();
            }
            if dst == "SINK".to_string() {
                l_ref_edge = edgename.clone();
                f_src_anchor = src.clone();
            }
        }
        final_anchor_list.insert(src);
        final_anchor_list.insert(dst);
    }
    println!("final_anchor number: {}", final_anchor_list.len());
    println!("final_edge number: {}", edge_info.len());
    
    let mut final_anchor = HashMap::new();
    for anchor in final_anchor_list.iter() {
        if let Some(info) = unadjacent_anchor.get(anchor) {
            final_anchor.insert(anchor.clone(), info.clone());
        }
    }

    // circularize the reference path
    let mut l_ref_seq = String::new();
    let mut l_ref_reads = HashSet::new();
    if let Some(edge) = edge_info.get(&l_ref_edge) {
        l_ref_seq.push_str(&edge.seq.clone());
        l_ref_reads.extend(edge.reads.clone().iter().cloned());
    }
    let mut head_read_count = 0usize;
    let tail_read_count = l_ref_reads.len();
    if let Some(edge) = edge_info.get(&f_ref_edge) {
        l_ref_seq.push_str(&edge.seq.clone());
        head_read_count = edge.reads.len();
        l_ref_reads = circularized_edge_reads(&l_ref_reads, &edge.reads);
    }
    println!(
        "circularized reference edge: {} reads (tail {}, head {})",
        l_ref_reads.len(),
        tail_read_count,
        head_read_count
    );
    edge_info.remove(&l_ref_edge);
    edge_info.remove(&f_ref_edge);
    edge_info.insert(l_ref_edge.clone(), EdgeInfo {
        seq: l_ref_seq,
        src: f_src_anchor.clone(),
        dst: f_dst_anchor.clone(),
        // Sorted so the emitted GFA is byte-reproducible: HashSet iteration
        // order is not stable across runs.
        reads: {
            let mut r: Vec<String> = l_ref_reads.iter().cloned().collect();
            r.sort();
            r
        },
        samples: HashSet::new(),
    });
    // Write final graph to disk.
    let _ = write_gfa(
        &final_anchor,
        &edge_info,
        "tmp.gfa",
    );
    let tmp_gfa_path = PathBuf::from("tmp.gfa");
    let mut graph = agg::GraphicalGenome::load_graph(&tmp_gfa_path).unwrap();

    // generate cigar. `min_edge_reads` is the `minimal_read_count` gate: an edge
    // is aligned (and can yield variants) only if it carries at least this many
    // reads. Lower it to recover low-frequency variants whose reads are
    // fragmented across many low-count edges on noisy long reads.
    generate_cigar(&mut graph, &ref_header, k, maxlength, min_edge_reads);
    let graph_output = output.with_extension("gfa");
    let _ = write_graph_from_graph(graph_output.to_str().unwrap(), &graph);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reconstructs the alt sequence a CIGAR implies, given the same ref_seq
    /// and alt_seq the CIGAR was computed from, so tests can prove a
    /// rewritten CIGAR still encodes exactly the same edit as the original
    /// rather than merely "having no X".
    fn reconstruct_alt(cigar: &str, _ref_seq: &str, alt_seq: &str) -> String {
        let alt_bytes = alt_seq.as_bytes();
        let mut ref_pos = 0usize;
        let mut alt_pos = 0usize;
        let mut out = Vec::new();
        let mut num = String::new();
        for c in cigar.chars() {
            if c.is_ascii_digit() {
                num.push(c);
                continue;
            }
            let len: usize = num.parse().unwrap();
            num.clear();
            match c {
                '=' | 'X' => {
                    out.extend_from_slice(&alt_bytes[alt_pos..alt_pos + len]);
                    ref_pos += len;
                    alt_pos += len;
                }
                'I' => {
                    out.extend_from_slice(&alt_bytes[alt_pos..alt_pos + len]);
                    alt_pos += len;
                }
                'D' => {
                    ref_pos += len;
                }
                _ => {}
            }
        }
        String::from_utf8(out).unwrap()
    }

    fn ref_consumed(cigar: &str) -> usize {
        let mut total = 0usize;
        let mut num = String::new();
        for c in cigar.chars() {
            if c.is_ascii_digit() {
                num.push(c);
                continue;
            }
            let len: usize = num.parse().unwrap();
            num.clear();
            if matches!(c, '=' | 'X' | 'D') {
                total += len;
            }
        }
        total
    }

    fn rs(names: &[&str]) -> HashSet<String> {
        names.iter().map(|s| s.to_string()).collect()
    }
    fn sorted(s: &HashSet<String>) -> Vec<String> {
        let mut v: Vec<String> = s.iter().cloned().collect();
        v.sort();
        v
    }

    /// The spliced edge is supported only by reads that wrap the origin, i.e.
    /// that appear on BOTH the SINK-bound tail and the SOURCE-bound head. The
    /// original code called `.intersection()` and discarded the iterator, so the
    /// head half was ignored entirely and the tail's reads passed through whole.
    #[test]
    fn circularized_edge_reads_keeps_only_reads_on_both_halves() {
        let tail = rs(&["r_wrap1", "r_tail_only", "r_wrap2"]);
        let head = ["r_wrap1".to_string(), "r_head_only".to_string(), "r_wrap2".to_string()];

        let got = circularized_edge_reads(&tail, &head);

        assert_eq!(sorted(&got), vec!["r_wrap1", "r_wrap2"]);
        // The pre-fix no-op behaviour was to return the tail unchanged.
        assert_ne!(sorted(&got), sorted(&tail));
    }

    #[test]
    fn circularized_edge_reads_handles_disjoint_and_empty_halves() {
        // No read wraps the origin: the spliced edge has no support.
        assert!(circularized_edge_reads(&rs(&["a", "b"]), &["c".to_string()]).is_empty());
        // Either half missing collapses to empty, never to the other half.
        assert!(circularized_edge_reads(&rs(&[]), &["a".to_string()]).is_empty());
        assert!(circularized_edge_reads(&rs(&["a"]), &[]).is_empty());
    }

    #[test]
    fn circularized_edge_reads_collapses_duplicate_head_entries() {
        // EdgeInfo.reads is a Vec, so the head half can repeat a read name.
        let got = circularized_edge_reads(&rs(&["a", "b"]), &["a".to_string(), "a".to_string()]);
        assert_eq!(sorted(&got), vec!["a"]);
    }

    /// The scoring parameters must make the clean two-gap edit outscore the
    /// one-gap-plus-mismatch shortcut inside a homopolymer run. This is the
    /// property the aligner's match/mismatch/gap settings exist to provide, and
    /// it is what the normalizer tests below can't see (they start from a
    /// hardcoded CIGAR). Real chrM 287-330 poly-C window from NA12877.
    #[test]
    fn aligner_emits_no_spurious_snp_in_polyc_run() {
        let ref_seq = "AAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGC";
        let alt_seq = "AAAAATTTCCACCAAACCCCCCCCCCTCCCCCCGCTTCTGGCCACAGC";

        // Insertion direction, and its exact ref<->alt mirror.
        for (r, a) in [(ref_seq, alt_seq), (alt_seq, ref_seq)] {
            let cigar = gap_open_aligner(r, a);
            assert!(
                !cigar.contains('X'),
                "poly-C indel must not be reported as a substitution, got {}",
                cigar
            );
            assert_eq!(ref_consumed(&cigar), r.len());
            assert_eq!(reconstruct_alt(&cigar, r, a), a);
        }
    }

    /// The gap-friendly scoring must not start swallowing real substitutions:
    /// an isolated SNP stays a single X, and an isolated deletion stays a D.
    #[test]
    fn aligner_preserves_isolated_snp_and_deletion() {
        let snp = gap_open_aligner("ACGTACGTACGTACGT", "ACGTACGAACGTACGT");
        assert_eq!(snp, "7=1X8=");

        let del = gap_open_aligner("ACGTACGTTTTACGTACGT", "ACGTACGTACGTACGT");
        assert_eq!(ref_consumed(&del), 19);
        assert_eq!(reconstruct_alt(&del, "ACGTACGTTTTACGTACGT", "ACGTACGTACGTACGT"), "ACGTACGTACGTACGT");
        assert!(del.contains('D') && !del.contains('X'), "got {}", del);
    }

    /// Identical sequences, and the degenerate empty-alt edge.
    #[test]
    fn aligner_handles_trivial_inputs() {
        let same = gap_open_aligner("ACGTACGTAC", "ACGTACGTAC");
        assert_eq!(same, "10=");

        let empty = gap_open_aligner("ACGTACGTAC", "");
        assert_eq!(empty, "10D");
    }

    #[test]
    fn homopolymer_insertion_absorbs_adjacent_mismatch() {
        // Real chrM 287-330 window from an NA12877 poly-C tract (7 C's, a
        // lone T at what would be chrM:310, then 5 more C's). The sample's
        // edge inserts 4 extra C's in that stretch; the raw aligner's optimal
        // scoring path represents 3 of them as a clean "I" but folds the 4th
        // together with the flanking T into a spurious "1X" (T>C) - a false
        // homoplasmic SNP that doesn't exist in any real read.
        let ref_seq = "AAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGC";
        let alt_seq = "AAAAATTTCCACCAAACCCCCCCCCCTCCCCCCGCTTCTGGCCACAGC";
        let raw_cigar = "23=4I1X20=";
        assert_eq!(reconstruct_alt(raw_cigar, ref_seq, alt_seq), alt_seq);

        let normalized = normalize_homopolymer_indels(raw_cigar, ref_seq, alt_seq);

        assert!(
            !normalized.contains('X'),
            "expected the homopolymer-adjacent mismatch to be absorbed, got {}",
            normalized
        );
        assert_eq!(ref_consumed(&normalized), ref_seq.len());
        assert_eq!(reconstruct_alt(&normalized, ref_seq, alt_seq), alt_seq);
    }

    #[test]
    fn genuine_snp_next_to_unrelated_indel_is_left_alone() {
        // A real substitution (ref 'G', alt 'A') sitting right next to an
        // insertion whose inserted base never equals the mismatched ref base
        // - there is no way to re-anchor this without changing content, so
        // it must be left as a real X.
        let ref_seq = "ACGTACGTG";
        let alt_seq = "ACGTACGTTTA"; // insert "TT" then mismatch G->A
        let raw_cigar = "8=2I1X";
        assert_eq!(reconstruct_alt(raw_cigar, ref_seq, alt_seq), alt_seq);

        let normalized = normalize_homopolymer_indels(raw_cigar, ref_seq, alt_seq);

        assert_eq!(normalized, raw_cigar);
    }

    #[test]
    fn homopolymer_deletion_absorbs_adjacent_mismatch() {
        // Exact role-swap (ref<->alt, I<->D) of homopolymer_insertion_absorbs_adjacent_mismatch:
        // a deletion is just an insertion viewed from the other sequence, so
        // this is the same real chrM poly-C edit, mirrored.
        let ref_seq = "AAAAATTTCCACCAAACCCCCCCCCCTCCCCCCGCTTCTGGCCACAGC";
        let alt_seq = "AAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGC";
        let raw_cigar = "23=4D1X20=";
        assert_eq!(reconstruct_alt(raw_cigar, ref_seq, alt_seq), alt_seq);

        let normalized = normalize_homopolymer_indels(raw_cigar, ref_seq, alt_seq);

        assert!(
            !normalized.contains('X'),
            "expected the homopolymer-adjacent mismatch to be absorbed, got {}",
            normalized
        );
        assert_eq!(ref_consumed(&normalized), ref_seq.len());
        assert_eq!(reconstruct_alt(&normalized, ref_seq, alt_seq), alt_seq);
    }
}
