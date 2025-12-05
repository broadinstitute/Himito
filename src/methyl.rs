use ndarray::Array2;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Read, Reader, Record};
use std::path::Path;
use std::{path::PathBuf, fs::File, io::{self, Write}};
use crate::{agg::*, methyl};
use std::collections::{HashMap, HashSet};
use serde_json::{json, Value};
use std::error::Error;
use csv::Writer;

pub fn find_path_on_graph(graph: &GraphicalGenome, read_name: &str) -> Vec<String> {
    let mut node = "SOURCE".to_string();
    let mut anchor_list: Vec<String> = graph.anchor.keys().cloned().collect();
    anchor_list.sort(); // Sort the anchor list
    let mut path_items: Vec<String> = Vec::new();

    while !path_items.contains(&node) && node != "SINK" {
        let mut found_edge = false;
        
        if let Some(outgoing_edges) = graph.outgoing.get(&node) {
            for edge in outgoing_edges {
                if let Some(edge_data) = graph.edges.get(edge) {
                    if let Some(reads) = edge_data.get("reads") {
                        if let Some(reads_array) = reads.as_array() {
                            if reads_array.iter().any(|v| v.as_str().map_or(false, |s| s == read_name)) {
                                path_items.push(node.clone());
                                path_items.push(edge.clone());
                                node = edge_data.get("dst")
                                    .and_then(|dst| dst.as_array())
                                    .and_then(|arr| arr.get(0))
                                    .and_then(|v| v.as_str())
                                    .unwrap_or("SINK")
                                    .to_string();
                                found_edge = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
        
        if !found_edge {
            if node == "SOURCE" {
                node = anchor_list.first().cloned().unwrap_or_else(|| "SINK".to_string());
            } else {
                let mut next_node_found = false;
                for (i, anchor) in anchor_list.iter().enumerate() {
                    if *anchor == node && i + 1 < anchor_list.len() {
                        node = anchor_list[i + 1].clone();
                        next_node_found = true;
                        break;
                    }
                }
                if !next_node_found {
                    node = "SINK".to_string();
                }
            }
        }
    }
    
    path_items
}
pub fn find_mapping_position(graph:&GraphicalGenome, r: &Record) -> HashMap<usize, (String, usize)> {
    let read_name = String::from_utf8_lossy(&r.qname()).to_string();
    // Get forward sequence
    // let forward_contig = if r.is_reverse() {
    //     reverse_complement(&String::from_utf8_lossy(&r.seq().as_bytes()).to_string())
    // } else {
    //     String::from_utf8_lossy(&r.seq().as_bytes()).to_string()
    // };
    let forward_contig = String::from_utf8_lossy(&r.seq().as_bytes()).to_string();
    let path_items = find_path_on_graph(&graph, &read_name);
    let mut anchor_list = Vec::new();
    for item in path_items.iter(){
        if item.starts_with("A"){
            anchor_list.push(item.clone());
        }
    }
    if anchor_list.is_empty(){
        return HashMap::new();
    }

    // Get the length of anchor sequence (k)
    let k = graph.anchor.get(&anchor_list[0])
        .and_then(|anchor| anchor.get("seq"))
        .and_then(|seq| seq.as_str())
        .map_or(0, |s| s.len());

    // Create anchor_seq to anchor mapping
    let mut anchor_seq_map: HashMap<String, String> = HashMap::new();
    for anchor in &anchor_list {
        if let Some(anchor_data) = graph.anchor.get(anchor) {
            if let Some(seq) = anchor_data.get("seq").and_then(|s| s.as_str()) {
                anchor_seq_map.insert(seq.to_string(), anchor.clone());
            }
        }
    }

    // Mapping anchor to reads
    let mut position_dict: HashMap<String, Vec<usize>> = HashMap::new();
    for (anchor_seq, _) in &anchor_seq_map {
        let anchor_rev = reverse_complement(anchor_seq);
        position_dict.entry(anchor_seq.clone()).or_insert_with(Vec::new);
        position_dict.entry(anchor_rev).or_insert_with(Vec::new);
    }

    // Find positions of anchor sequences in the read
    if forward_contig.len() > k {
        for i in 1..(forward_contig.len() - k + 1) {
            if let Some(kmer) = forward_contig.get(i..i+k) {
                if position_dict.contains_key(kmer) {
                    position_dict.entry(kmer.to_string())
                        .or_insert_with(Vec::new)
                        .push(i);
                }
            }
        }
    }

    // Exclude multiple mapping
    let mut final_anchor: HashMap<String, usize> = HashMap::new();
    for anchor in &anchor_list {
        if let Some(anchor_data) = graph.anchor.get(anchor) {
            if let Some(seq) = anchor_data.get("seq").and_then(|s| s.as_str()) {
                let anchor_rev = reverse_complement(seq);
                let pos_seq = position_dict.get(seq).map_or(Vec::new(), |v| v.clone());
                let pos_rev = position_dict.get(&anchor_rev).map_or(Vec::new(), |v| v.clone());
                let position_list = [pos_seq.clone(), pos_rev.clone()].concat();
                // if r.is_reverse(){
                //     println!("{}, {}", pos_seq.len(), pos_rev.len());
                // }
                if position_list.len() == 1 && pos_rev.is_empty() {
                    final_anchor.insert(anchor.clone(), pos_seq[0]);
                }
            }
        }
    }

    // Find base pair level mapping between reads and graph entities
    let mut read_position_mapping: HashMap<usize, (String, usize)> = HashMap::new();
    
    // Sort final anchors by position
    let mut sorted_final_anchor: Vec<(usize, String)> = final_anchor
        .iter()
        .map(|(k, v)| (*v, k.clone()))
        .collect();
    sorted_final_anchor.sort_by_key(|item| item.0);
    
    if sorted_final_anchor.len() < 1 {
        // println!("{}, {}", read_name, final_anchor.len());
        return read_position_mapping;
    }

    // Extract first anchor information
    let (first_anchor_position, first_anchor) = &sorted_final_anchor[0];
    // Find the first edge
    if let Some(incoming_edges) = graph.incoming.get(first_anchor) {
        for edge in incoming_edges {
            if let Some(edge_data) = graph.edges.get(edge) {
                if let Some(reads) = edge_data.get("reads").and_then(|r| r.as_array()) {
                    if reads.iter().any(|v| v.as_str().map_or(false, |s| s == read_name)) {
                        // Check sequence
                        if let Some(edge_seq) = edge_data.get("seq").and_then(|s| s.as_str()) {
                            if let Some(contig_prefix) = forward_contig.get(0..*first_anchor_position) {
                                if edge_seq == contig_prefix {
                                    for first_pos in 0..*first_anchor_position {
                                        read_position_mapping.insert(first_pos, (edge.clone(), first_pos));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }   

    // Middle part
    for (read_position, anchor) in &sorted_final_anchor {
        // Map read positions to anchor offsets
        for i in 0..k {
            read_position_mapping.insert(read_position + i, (anchor.clone(), i));
        }

        // Process connected edges
        if let Some(outgoing_edges) = graph.outgoing.get(anchor) {
            for edge in outgoing_edges {
                if let Some(edge_data) = graph.edges.get(edge) {
                    if let Some(reads) = edge_data.get("reads").and_then(|r| r.as_array()) {
                        if reads.iter().any(|v| v.as_str().map_or(false, |s| s == read_name)) {
                            if let Some(dst_node) = graph.outgoing.get(edge).and_then(|nodes| nodes.get(0)) {
                                if dst_node == "SINK" {
                                    if let Some(edge_seq) = edge_data.get("seq").and_then(|s| s.as_str()) {
                                        if let Some(contig_suffix) = forward_contig.get(read_position + k..) {
                                            if edge_seq == contig_suffix {
                                                let edge_length = edge_seq.len();
                                                for j in 0..edge_length {
                                                    read_position_mapping.insert(read_position + k + j, (edge.clone(), j));
                                                }
                                            }
                                        }
                                    }
                                } else if let Some(&dst_pos) = final_anchor.get(dst_node) {
                                    if let Some(edge_seq) = edge_data.get("seq").and_then(|s| s.as_str()) {
                                        if let Some(contig_middle) = forward_contig.get(read_position + k..dst_pos) {
                                            if edge_seq == contig_middle {
                                                let edge_length = edge_seq.len();
                                                for j in 0..edge_length {
                                                    read_position_mapping.insert(read_position + k + j, (edge.clone(), j));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
    }

    read_position_mapping

}

pub fn get_methylation_read(r: &Record, mod_char: char) -> HashMap<usize, f32> {
    let mut methyl_pos_dict: HashMap<usize, f32> = HashMap::new();
    
    // Get forward sequence
    // let forward_sequence = if r.is_reverse() {
    //     reverse_complement(&String::from_utf8_lossy(&r.seq().as_bytes()).to_string())
    // } else {
    //     String::from_utf8_lossy(&r.seq().as_bytes()).to_string()
    // };
    let forward_sequence = String::from_utf8_lossy(&r.seq().as_bytes()).to_string();
    
    // Check for modification data
    if let Ok(mods) = r.basemods_iter() {
        // Iterate over the modification types
        for res in mods {
            if let Ok( (position, m) ) = res {
                if m.modified_base as u8 as char != mod_char{
                    continue
                }
                // let strand = mod_metadata.strand;
                let pos_usize = position as usize;
                let qual = m.qual as f32 / 255.0;
                let motif: String = if r.is_reverse(){
                    forward_sequence[pos_usize-1..pos_usize +1 ].to_string()
                    
                }else{
                    forward_sequence[pos_usize..pos_usize +2 ].to_string()
                };
                    
                if motif == "CG" {
                    if r.is_reverse(){
                        methyl_pos_dict.insert(pos_usize-1, qual);

                    }else{
                        methyl_pos_dict.insert(pos_usize, qual);

                    }
                    
                }else{
                    println!("{},{},{}", pos_usize, qual, motif);
                    continue
                }

            }
        }                    

    }
    
    methyl_pos_dict
}

pub fn validation(graph: &GraphicalGenome, read_position_mapping:&HashMap<usize, (String, usize)>, r:&Record){
    let read_sequence = String::from_utf8_lossy(&r.seq().as_bytes()).to_string();
    for (position, (item, offset)) in read_position_mapping {
        if item.starts_with("A") {
            let anchor_data = graph.anchor.get(item).expect(&format!("Anchor {} not found", item));
            let seq = anchor_data.get("seq").and_then(|s| s.as_str()).expect(&format!("No seq for anchor {}", item));
            
            let query_char = read_sequence.chars().nth(*position).expect(&format!("Position {} out of bounds in query", position));
            let anchor_char = seq.chars().nth(*offset).expect(&format!("Offset {} out of bounds in anchor {}", offset, item));
            if query_char != anchor_char{
                println!("Mismatch at position {} for item {} offset {}", position, item, offset);
            }
            // assert_eq!(query_char, anchor_char, "Mismatch at position {} for item {} offset {}", position, item, offset);
        } else if item.starts_with("E") {
            let edge_data = graph.edges.get(item).expect(&format!("Edge {} not found", item));
            let seq = edge_data.get("seq").and_then(|s| s.as_str()).expect(&format!("No seq for edge {}", item));
            
            let query_char = read_sequence.chars().nth(*position).expect(&format!("Position {} out of bounds in query", position));
            let edge_char = seq.chars().nth(*offset).expect(&format!("Offset {} out of bounds in edge {}", offset, item));
            if query_char != edge_char{
                println!("Mismatch at position {} for item {} offset {}", position, item, offset);
            }
            // assert_eq!(query_char, edge_char, "Mismatch at position {} for item {} offset {}", position, item, offset);
        }
    }
}

/// Maps positions from alternate sequence to reference sequence based on CIGAR string
pub fn get_reference_coordinates(cigar: &str, ref_start: usize) -> HashMap<usize, usize> {
    let mut ref_pos = 0;
    let mut alt_pos = 0;
    
    // Parse CIGAR string into operations
    let mut operations = Vec::new();
    let mut num = String::new();
    
    for c in cigar.chars() {
        if c.is_digit(10) {
            num.push(c);
        } else {
            if !num.is_empty() {
                let length = num.parse::<usize>().unwrap();
                operations.push((length, c));
                num.clear();
            }
        }
    }
    
    // Process each operation
    let mut position_mapping = HashMap::new();
    
    for (length, op) in operations {
        match op {
            '=' | 'M' => {  // Match
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    position_mapping.insert(alt_pos + i, pos);
                }
                ref_pos += length;
                alt_pos += length;
            },
            'X' => {  // Mismatch
                for i in 0..length {
                    let pos = ref_start + ref_pos + i;
                    position_mapping.insert(alt_pos + i, pos);
                }
                ref_pos += length;
                alt_pos += length;
            },
            'I' => {  // Insertion
                let pos = ref_start + ref_pos - 1;  // Mapping to one base pair before
                for i in 0..length {
                    position_mapping.insert(alt_pos + i, pos);
                }
                alt_pos += length;
            },
            'D' => {  // Deletion, nothing will map back to reference coordinates
                ref_pos += length;
            },
            _ => {
                // Handle other CIGAR operations or error cases
                // You might want to add logging or proper error handling here
            }
        }
    }
    
    position_mapping
}

pub fn add_methylation_to_graph(
    mut graph: &mut GraphicalGenome, 
    methyl_pos_dict: &HashMap<usize, f32>, 
    mapping_position: &HashMap<usize, (String, usize)>, 
    read_name: &str
) {
    for (pos, likelihood) in methyl_pos_dict {
        if let Some((item, offset)) = mapping_position.get(pos) {
            if item.starts_with("A") {
                if let Some(anchor_data) = graph.anchor.get_mut(item) {
                    // Assert for debugging
                    let seq = anchor_data.get("seq")
                        .and_then(|s| s.as_str())
                        .expect(&format!("No seq for anchor {}", item));
                    let c = seq.chars().nth(*offset)
                        .expect(&format!("Offset {} out of bounds in anchor {}", offset, item));
                    assert_eq!(c, 'C', "Expected C at position {} in anchor {}, found {}", offset, item, c);
                    
                    // Update methylation data
                    if !anchor_data.get("methyl").is_some() {
                        anchor_data["methyl"] = json!({});
                    }
                    
                    if let Some(methyl) = anchor_data.get_mut("methyl").and_then(|m| m.as_object_mut()) {
                        if !methyl.contains_key(read_name) {
                            methyl.insert(read_name.to_string(), json!([]));
                        }
                        
                        if let Some(read_methyl) = methyl.get_mut(read_name).and_then(|rm| rm.as_array_mut()) {
                            if !read_methyl.contains(&json!([offset, likelihood])){
                                read_methyl.push(json!([offset, likelihood]));
                            }
                            
                        }
                    }
                }
            } else if item.starts_with("E") {
                if let Some(edge_data) = graph.edges.get_mut(item) {
                    // Assert for debugging
                    let seq = edge_data.get("seq")
                        .and_then(|s| s.as_str())
                        .expect(&format!("No seq for edge {}", item));
                    let c = seq.chars().nth(*offset)
                        .expect(&format!("Offset {} out of bounds in edge {}", offset, item));
                    assert_eq!(c, 'C', "Expected C at position {} in edge {}, found {}", offset, item, c);
                    
                    // Update methylation data
                    if !edge_data.get("methyl").is_some() {
                        edge_data["methyl"] = json!({});
                    }
                    
                    if let Some(methyl) = edge_data.get_mut("methyl").and_then(|m| m.as_object_mut()) {
                        if !methyl.contains_key(read_name) {
                            methyl.insert(read_name.to_string(), json!([]));
                        }
                        
                        if let Some(read_methyl) = methyl.get_mut(read_name).and_then(|rm| rm.as_array_mut()) {
                            if !read_methyl.contains(&json!([offset, likelihood])){
                                read_methyl.push(json!([offset, likelihood]));
                            }
                        }
                    }
                }
            } else {
                println!("{}, {}", item, pos);
            }
        }
    }
}

fn find_methylation_signal_on_major_haplotype(graph: &GraphicalGenome) -> HashMap<(usize, usize, String), HashMap<String, f64>> {
    let mut methyl: HashMap<(usize, usize, String), HashMap<String, f64>> = HashMap::new();
    let major_haplotype = construct_major_haplotype_entitylist(graph);
    let major_haplotype_sequence = construct_major_haplotype(graph);
    let mut startpos = 0;
    
    let mut anchor_keys: Vec<_> = graph.anchor.keys().cloned().collect();
    anchor_keys.sort();
    
    let firstanchor = &anchor_keys[0];
    let anchorseq = graph.anchor[firstanchor]["seq"].as_str().unwrap();
    let k = anchorseq.len();
    
    for item in major_haplotype {
        let (methyl_info_list, seq, reference_position_mapping) = if item.starts_with("A") {
            // Anchor case
            let methyl_info = match graph.anchor[&item].get("methyl") {
                Some(m) => m.as_object().unwrap_or(&serde_json::Map::new()).clone(),
                None => serde_json::Map::new()
            };
            
            let seq = graph.anchor[&item]["seq"].as_str().unwrap().to_string();
            let cigar = format!("{}=", k);
            let ref_start = graph.anchor[&item]["pos"].as_u64().unwrap() as usize;
            let reference_mapping = get_reference_coordinates(&cigar, ref_start);
            
            (methyl_info, seq, reference_mapping)
        } else if item.starts_with("E") {
            // Edge case
            let methyl_info = match graph.edges[&item].get("methyl") {
                Some(m) => m.as_object().unwrap_or(&serde_json::Map::new()).clone(),
                None => serde_json::Map::new()
            };
            
            let seq = graph.edges[&item]["seq"].as_str().unwrap().to_string();
            let cigar = match graph.edges[&item].get("variants") {
                Some(v) => v.as_str().unwrap_or("").to_string(),
                None => "".to_string()
            };
            
            let src = graph.incoming[&item][0].as_str();
            let ref_start = if src == "SOURCE" {
                0
            } else {
                graph.anchor[src]["pos"].as_u64().unwrap() as usize + k
            };
            
            let reference_mapping = if cigar.is_empty() {
                HashMap::new()
            } else {
                get_reference_coordinates(&cigar, ref_start)
            };
            
            (methyl_info, seq, reference_mapping)
        } else {
            // Other case
            (serde_json::Map::new(), "".to_string(), HashMap::new())
        };
        
        if methyl_info_list.is_empty() {
            startpos += seq.len();
            continue;
        }
        
        for (read, info_list) in methyl_info_list {
            let info_array = info_list.as_array().unwrap();
            for info in info_array {
                let pos = info[0].as_u64().unwrap() as usize;
                let likelihood = info[1].as_f64().unwrap();
                
                let currentpos = startpos + pos;
                let referencepos = reference_position_mapping.get(&pos).cloned();
                if currentpos >= major_haplotype_sequence.len() || 
                   major_haplotype_sequence.chars().nth(currentpos).unwrap() != 'C' {
                    continue;
                }
                if currentpos + 1 >= major_haplotype_sequence.len() || 
                   major_haplotype_sequence.chars().nth(currentpos + 1).unwrap() != 'G' {
                    continue;
                }

                let motif = major_haplotype_sequence[currentpos..currentpos + 2].to_string();
                
                if let Some(referencepos_value) = referencepos {
                    methyl.entry((currentpos, referencepos_value, motif)).or_insert_with(HashMap::new)
                     .insert(read.clone(), likelihood);
                }

            }
        }
        
        startpos += seq.len();
    }
    
    methyl
}

fn write_bed(
    methyl: HashMap<(usize, usize, String), HashMap<String, f64>>,
    output_file: &PathBuf,
    min_prob: f64
) -> std::io::Result<()> {
    let mut file = File::create(Path::new(output_file))?;

    // Write BED header
    writeln!(file, "##fileformat=BED")?;
    writeln!(file, "##haplotype=majorhaplotype")?;
    writeln!(
        file,
        "#CHROM\tRef_start\tRef_end\tAsm_start\tAsm_end\tMotif\tMod_rate\tUnmod_rate\tCov\tMod_count\tUnmod_count"
    )?;
    // let min_prob = 0.5;
    // let mut methylation_signal: HashMap<(usize, Option<usize>), (f64, f64)> = HashMap::new();
    let mut methyl_info = Vec::new();
    for ((pos, refpos, motif), d) in methyl.iter() {
        let mut methyl_count = 0;
        let mut unmethyl_count = 0;
        let total_count = d.len();
        
        for (_, &likelihood) in d.iter() {
            if likelihood > min_prob {
                methyl_count += 1;
            } else if likelihood < 1.0 - min_prob {
                unmethyl_count += 1;
            }
        }
    
        let methyl_rate = methyl_count as f64 / total_count as f64;
        let unmethyl_rate = unmethyl_count as f64 / total_count as f64;
        // 1-based coordinates
        let ref_pos_start = refpos + 1;
        let ref_pos_end = refpos + 2;
        let asm_pos_start = pos + 1;
        let asm_pos_end = pos + 2;
        methyl_info.push((ref_pos_start, ref_pos_end, asm_pos_start, asm_pos_end, motif, methyl_rate,unmethyl_rate,total_count, methyl_count, unmethyl_count));

    }
    methyl_info.sort_by_key(|item| item.0);
    for methyl_list in methyl_info.iter(){
        let (ref_pos_start, ref_pos_end, asm_pos_start, asm_pos_end, motif, methyl_rate,unmethyl_rate,total_count, methyl_count, unmethyl_count) = methyl_list;
        writeln!(file, 
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            "ChrM", ref_pos_start, ref_pos_end, asm_pos_start, asm_pos_end, motif, methyl_rate,unmethyl_rate,total_count,methyl_count, unmethyl_count,
            )?;

    }


    Ok(())
}

fn write_methylation_to_csv<P: AsRef<Path>>(
    methyl: HashMap<(usize, usize, String), HashMap<String, f64>>,
    min_prob:f64,
    path: P
) -> Result<(), Box<dyn Error>> {
    // Create a file and CSV writer
    let file = File::create(path)?;
    let mut writer = Writer::from_writer(file);

    // find all the readsets and positions
    let mut refpos_list = Vec::new();
    let mut read_set = HashSet::new();
    for ((pos, refpos, motif), d) in methyl.iter() {
        refpos_list.push(refpos);
        for (read_name, &likelihood) in d.iter() {
            read_set.insert(read_name.clone());
        }
    }
    refpos_list.sort();
    // println!("{}", refpos_list.len());

    let ref_pos_dict: HashMap<usize, usize> = refpos_list
        .iter()
        .enumerate()
        .map(|(i, rpos)| (**rpos, i))
        .collect();

    let mut read_vec: Vec<String> = read_set.into_iter().collect();
    read_vec.sort();
    let read_set_dict: HashMap<String, usize> = read_vec
        .iter()
        .enumerate()
        .map(|(i, read)| (read.clone(), i))
        .collect();

    // construct matrix
    let mut matrix = Array2::<f64>::zeros((refpos_list.len(), read_vec.len()));
    for ((pos, refpos, motif), d) in methyl.iter() {
        let row_index = ref_pos_dict.get(&refpos).unwrap();
        for (read_name, &likelihood) in d.iter() {
            let col_index = read_set_dict.get(read_name).unwrap();
            // matrix[[*row_index, *col_index]] = likelihood;
            if likelihood > min_prob {
                matrix[[*row_index, *col_index]] = 1.0;
            }else if likelihood < 1.0-min_prob {
                 matrix[[*row_index, *col_index]] = -1.0;
            }
        }
    }
    // Prepare header row (with empty cell for the corner)
    let mut header = vec!["methylation".to_string()];
    header.extend(read_vec.iter().cloned());
    
    // Write header
    writer.write_record(&header)?;
    
    // Write each row with its row name
    for (row_idx, refpos_name) in refpos_list.iter().enumerate() {
        let mut row = vec![refpos_name.to_string().clone()];
        
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

fn write_methylation_alignment_to_csv<P: AsRef<Path>>(
    methyl: HashMap<String, HashMap<usize, Vec<f32>>>,
    min_prob:f32,
    strand_dict: HashMap<String, bool>,
    path: P,
    fraction_threshold:f64
) -> Result<(), Box<dyn Error>> {
    // Create a file and CSV writer
    let file = File::create(path)?;
    let mut writer = Writer::from_writer(file);

    // find all the readsets and positions
    let mut refposlist = HashSet::new();
    let mut read_set = HashSet::new();
    let mut coverage_dict:HashMap<usize, usize> = HashMap::new();
    for (read_name, d) in methyl.iter() {
        read_set.insert(read_name.clone());
        let strand = strand_dict.get(read_name).unwrap_or(&false);
        for (refpos, likelihood) in d.iter() {
            refposlist.insert(refpos.clone());
            let count = coverage_dict.get(&refpos).unwrap_or(&0) + 1;
            coverage_dict.insert(*refpos, count);
        }
    }
    let mut refpos_list : Vec<usize> = Vec::new();
    let total_read_number = read_set.len();
    for (selected_pos, c) in coverage_dict.iter(){
        if (fraction_threshold * total_read_number as f64) < *c as f64 {
            refpos_list.push(selected_pos.clone());
        }
    }
    refpos_list.sort();
    // println!("{}", refpos_list.len());

    let ref_pos_dict: HashMap<usize, usize> = refpos_list
        .iter()
        .enumerate()
        .map(|(i, rpos)| (rpos.clone(), i))
        .collect();

    let mut read_vec: Vec<String> = read_set.into_iter().collect();
    read_vec.sort();
    let read_set_dict: HashMap<String, usize> = read_vec
        .iter()
        .enumerate()
        .map(|(i, read)| (read.clone(), i))
        .collect();
    println!("read number {}, position number, {}", read_vec.len(),refpos_list.len());
    // construct matrix
    let mut matrix = Array2::<f64>::zeros((refpos_list.len(), read_vec.len()));
    println!("{:?}", matrix.shape());
    for (read_name, d) in methyl.iter() {
        let col_index = read_set_dict.get(read_name).unwrap();
        // let strand = strand_dict.get(read_name).unwrap_or(&false);
        for (refpos, likelihood) in d.iter() {
            if !ref_pos_dict.contains_key(&refpos){
                continue
            }
            let row_index = ref_pos_dict.get(&refpos).unwrap();
            // matrix[[*row_index, *col_index]] = likelihood;
            if likelihood.len() > 1{
                println!("{}, {}, {:?}", refpos, read_name, likelihood);
                continue
            }
            let likelihood_value = likelihood[0];
    
            if likelihood_value > min_prob {
                matrix[[*row_index, *col_index]] = 1.0;
            }else if likelihood_value < 1.0-min_prob {
                matrix[[*row_index, *col_index]] = -1.0;
            }
        }
        
    }
    // Prepare header row (with empty cell for the corner)
    let mut header = vec!["methylation".to_string()];
    header.extend(read_vec.iter().cloned());
    
    // Write header
    writer.write_record(&header)?;
    
    // Write each row with its row name
    for (row_idx, refpos_name) in refpos_list.iter().enumerate() {
        let mut row = vec![refpos_name.to_string().clone()];
        
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

pub fn start (graph_file: &PathBuf, bam_file: &PathBuf, output_file: &PathBuf, min_prob: f64, major_haploype: bool) {
    println!("Add Methylation Signals!");
    // annotate graph
    let mut graph = GraphicalGenome::load_graph(graph_file).unwrap();
    println!("Processing BAM file");
    let mut bam = Reader::from_path(bam_file).unwrap();
    let mut read_name_set = HashSet::new();
    let mut methyl_single_read: HashMap<String,HashMap<usize, Vec<f32>> >= HashMap::new(); // extract directly from reference map
    let mut strand_dict:HashMap<String, bool> = HashMap::new();
    for record in bam.records(){
        let r = record.expect("Failed to read BAM record");
        // let mut read_sequence = String::from_utf8_lossy(&r.seq().as_bytes()).to_string();
        // if r.is_reverse(){
        //     read_sequence = reverse_complement(&read_sequence);
        // }
        if r.is_unmapped(){
            continue
        }
        if r.is_supplementary(){
            continue
        }
        let read_name = String::from_utf8_lossy(&r.qname()).to_string();
        // if read_name_set.contains(&read_name) {
        //     continue
        // }
        strand_dict.insert(read_name.clone(), r.is_reverse());
        read_name_set.insert(read_name.clone());
        let read_position_mapping = find_mapping_position(&graph, &r);
        // Validate read_position_mapping
        validation(&graph, &read_position_mapping, &r);
        let methyl_pos_dict = get_methylation_read(&r, 'm');
        
        // map to reference coordinates
        let aligned_pairs:HashMap<usize, usize> = r.aligned_pairs_full()
                                .into_iter()
                                .filter_map(|pair| {
                                // Only include pairs where both positions are defined
                                match pair {
                                    [Some(q), Some(r)] => Some((q as usize, r as usize)),
                                    _ => None
                                }
                            })
                            .collect();

        for (read_pos, likelihood) in methyl_pos_dict.iter(){
            if aligned_pairs.contains_key(read_pos){
                let ref_pos = aligned_pairs.get(read_pos).copied().unwrap_or(0);
                methyl_single_read
                                    .entry(read_name.clone())
                                    .or_insert_with(HashMap::new)
                                    .entry(ref_pos)
                                    .or_insert_with(Vec::new)
                                    .push(*likelihood);
            }

        }

        if read_position_mapping.len() == 0 {
            continue
        }
        if methyl_pos_dict.len() == 0 {
            continue
        }
        // println!("{}", read_position_mapping.len());
        add_methylation_to_graph(&mut graph, &methyl_pos_dict, &read_position_mapping, &read_name);  
    }
    let graph_output = output_file.with_extension("methyl.gfa");
    let _ = write_graph_from_graph(graph_output.to_str().unwrap(), &graph);
    

    // extract major haplotype path and then construct a matrix
    let methyl_dict =  find_methylation_signal_on_major_haplotype(&graph);
    let _ = write_bed(methyl_dict.clone(), output_file, min_prob);

    // extranct read level methylation patterns
    // println!("{}", methyl_single_read.len());
    let matrix_output = output_file.with_extension("methylation_per_read.csv");
    if major_haploype {
        let _ = write_methylation_to_csv(methyl_dict, min_prob, matrix_output);
    }else{
        let _ = write_methylation_alignment_to_csv(methyl_single_read, min_prob as f32, strand_dict, matrix_output, 0.1);

    }
    
    

}