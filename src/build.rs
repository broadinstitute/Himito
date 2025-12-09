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
use std::io::{self, BufRead, BufReader, Write};


pub fn reverse_complement(kmer: &str) -> String {
    kmer.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => panic!("Unexpected character: {}", c),
        })
        .collect()
}

pub fn map_to_genome(contig: &str, k: usize) -> HashMap<String, Vec<usize>> {
    let mut position_dict: HashMap<String, Vec<usize>> = HashMap::new();
    
    for i in 0..contig.len() - k + 1 {
        let kmer = contig[i..i+k].to_string();
        position_dict.entry(kmer)
            .or_insert_with(Vec::new)
            .push(i);
    }
    
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
    position_dict: &HashMap<String, Vec<usize>>,
    k: usize
) -> HashMap<String, AnchorInfo> {
    let mut anchor_updated_list = Vec::new();
    
    // Collect kmers that appear exactly once
    for (kmer, positions) in position_dict {
        if positions.len() == 1 {
            anchor_updated_list.push(kmer);
        }
    }

    let mut anchor_info = HashMap::new();
    for kmer in anchor_updated_list {
        if let Some(pos) = position_dict.get(kmer).and_then(|v| v.first()) {
            let anchor_name = format!("A{:06}", pos);
            anchor_info.insert(
                anchor_name,
                AnchorInfo {
                    seq: kmer.clone(),
                    pos: *pos,
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
    contig: String,
    k: usize,
) -> (HashMap<String, usize>, HashMap<String, Vec<usize>>) {
    // Create anchor_seq to anchor_name mapping
    let anchor_seq_map: HashMap<String, String> = anchor_info
        .iter()
        .map(|(anchor, info)| (info.seq.clone(), anchor.clone()))
        .collect();

    // Initialize position dictionary
    let mut position_dict: HashMap<String, Vec<usize>> = HashMap::new();
    for anchor_seq in anchor_seq_map.keys() {
        let anchor_rev = reverse_complement(anchor_seq);
        position_dict.insert(anchor_seq.clone(), Vec::new());
        position_dict.insert(anchor_rev, Vec::new());
    }

    // Find positions of kmers in contig
    for i in 0..contig.len() - k + 1 {
        let kmer = contig[i..i+k].to_string();
        if position_dict.contains_key(&kmer) {
            position_dict.get_mut(&kmer).unwrap().push(i);
        }
    }

    let mut a = HashMap::new();
    let mut svs = HashMap::new();

    // Process each anchor
    for (anchor, info) in anchor_info {
        let anchor_seq = &info.seq;
        let anchor_rev = reverse_complement(anchor_seq);
        
        // Combine positions from forward and reverse sequences
        let position_list = [
            position_dict.get(anchor_seq).unwrap_or(&Vec::new()).clone(),
            position_dict.get(&anchor_rev).unwrap_or(&Vec::new()).clone()
        ].concat();
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
    contig: String,
    contigname: String,
    sample: String,
    anchorseq: &HashMap<String, String>,
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



pub fn create_edge_file(
    all_seq: &HashMap<String, String>,
    final_anchor: &HashMap<String, &AnchorInfo>,
    k: usize,
    threshold: usize
) -> (HashMap<String, EdgeInfo>, HashMap<String, Vec<String>>) {
    let mut edge_info: HashMap<String, EdgeInfo> = HashMap::new();
    let mut outgoing: HashMap<String, Vec<String>> = HashMap::new();
    
    let anchorseq: HashMap<_, _> = final_anchor
        .iter()
        .map(|(anchor, info)| (info.seq.clone(), anchor.clone()))
        .collect();

    let mut contig_index = 0;
    let bar = ProgressBar::new(all_seq.len() as u64);

    for (contig_name, contig) in all_seq.iter() {
        if contig.len() < k + 1{
            continue
        }
        contig_index += 1;
        bar.inc(1);
        // let contig_name = record.id().to_string();
        // let contig = String::from_utf8_lossy(record.seq()).to_string();
        // let sample_name = contig_name.split('|').last().unwrap_or("").to_string();
        let sample_name = if contig_name.contains('|') {
                contig_name.split('|').last().unwrap_or("").to_string()
            } else {
                "".to_string()
            };
        
        let (a, _svs) = mapping_info(final_anchor, contig.to_string(), k);

        if a.len() < threshold {
            continue;
        }

        let mut splitposlist: Vec<_> = a.values().copied().collect();
        splitposlist.sort();
        let mut edgeindex = 0;
        let mut src_pos = 0;

        // Process all positions except the last
        for &dst_pos in &splitposlist {
            if dst_pos - src_pos < k+1 {
                continue
            }
            let e = construct_edges(
                src_pos,
                dst_pos,
                k,
                contig.to_string(),
                contig_name.to_string(),
                sample_name.clone(),
                &anchorseq,
            );
            
            let src = &e.src;
            let edgelist = outgoing.entry(src.clone()).or_default();
            
            // Try to find matching existing edge
            let mut found_match = false;
            for edge_name in edgelist.iter() {
                let existing_edge = edge_info.get_mut(edge_name).unwrap();
                if existing_edge.dst == e.dst && existing_edge.seq == e.seq {
                    existing_edge.reads.extend(e.reads.clone());
                    existing_edge.samples = existing_edge.samples.union(&e.samples).cloned().collect();
                    found_match = true;
                    break;
                }
            }
            
            // If no match found, create new edge
            if !found_match {
                let edgename = format!("E{:05}.{:04}", contig_index, edgeindex);
                edge_info.insert(edgename.clone(), e.clone());
                edgelist.push(edgename);
                edgeindex += 1;
            }
            // debuging why edge sequence are empty
            let mut edge_seq = if src_pos == 0 {
            let seq = contig.get(0..dst_pos).unwrap_or_default().to_string();
                if seq.is_empty() {
                    println!("Warning: Empty edge sequence created for SOURCE. dst_pos: {}, contig_len: {}", dst_pos, contig.len());
                }
                seq
            } else {
                let seq = contig.get(src_pos + k..dst_pos).unwrap_or_default().to_string();
                if seq.is_empty() {
                    println!("Warning: Empty edge sequence created. src_pos: {}, dst_pos: {}, k: {}, contig_len: {}", 
                            src_pos, dst_pos, k, contig.len());
                }
                seq
            };
            
            src_pos = dst_pos;
        }

        // Process final edge to end of contig
        let dst_pos = contig.len();
        if dst_pos - src_pos > k+ 1{
            let e = construct_edges(
                src_pos,
                dst_pos,
                k,
                contig.to_string(),
                contig_name.to_string(),
                sample_name.clone(),
                &anchorseq,
            );
            
            let src = &e.src;
            let edgelist = outgoing.entry(src.clone()).or_default();
        
            // Try to find matching existing edge for final segment
            let mut found_match = false;
            for edge_name in edgelist.iter() {
                let existing_edge = edge_info.get_mut(edge_name).unwrap();
                if existing_edge.dst == e.dst && existing_edge.seq == e.seq {
                    existing_edge.reads.extend(e.reads.clone());
                    existing_edge.samples = existing_edge.samples.union(&e.samples).cloned().collect();
                    found_match = true;
                    break;
                }
            }
            
            // If no match found for final segment, create new edge
            if !found_match {
                let edgename = format!("E{:05}.{:04}", contig_index, edgeindex);
                edge_info.insert(edgename.clone(), e.clone());
                edgelist.push(edgename);
            }
        }
        std::thread::sleep(std::time::Duration::from_millis(50));

    }
    bar.finish();

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

pub fn start(output: &PathBuf, k: usize, read_path: &PathBuf, reference_path: &PathBuf) {
    // Read reference records into a vector
    let ref_reader = Reader::from_file(reference_path).unwrap();
    let reference_sequence: Vec<Record> = ref_reader.records().map(|r| r.unwrap()).collect();
    let ref_seq = String::from_utf8_lossy(reference_sequence[0].seq()).to_string();
    
    let position_dict = map_to_genome(&ref_seq, k);
    let anchors = create_anchors(&position_dict, k);
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
    let mut edge_info = n_edge_info.clone();
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
    if let Some(edge) = edge_info.get(&f_ref_edge) {
        l_ref_seq.push_str(&edge.seq.clone());
        l_ref_reads.intersection(&edge.reads.clone().iter().cloned().collect::<HashSet<_>>());
    }
    edge_info.remove(&l_ref_edge);
    edge_info.remove(&f_ref_edge);
    edge_info.insert(l_ref_edge.clone(), EdgeInfo {
        seq: l_ref_seq,
        src: f_src_anchor.clone(),
        dst: f_dst_anchor.clone(),
        reads: l_ref_reads.iter().cloned().collect(),
        samples: HashSet::new(),
    });

    // Write final graph to disk.
    let _ = write_gfa(
        &final_anchor,
        &edge_info,
        output.to_str().unwrap(),
    );

}
