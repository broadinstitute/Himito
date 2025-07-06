use crate::agg::*;
use std::{path::PathBuf, fs::File, io::{self, Write}};
use std::collections::{HashMap, HashSet};
use intervals::{Interval, Open};

fn find_all_reads(graph:&GraphicalGenome) -> HashSet<String> {
    let mut readset = HashSet::new();
    for edge in graph.edges.keys() {
        if let Some(reads) = graph.edges[edge].get("reads").and_then(|v| v.as_array()) {
            for read in reads.iter() {
                if let Some(read_str) = read.as_str() {
                    readset.insert(read_str.to_string());
                }
            }
        }
    }
    readset
}

fn search_bounding_anchors(
    start_pos: i32,
    end_pos: i32,
    sp_anchors: &Vec<&String>,
    graph: &GraphicalGenome,
) -> Option<(String, String)> {

    // Collect positions for each anchor
    let poslist: Vec<i32> = sp_anchors
        .iter()
        .map(|anchor| {
            graph.anchor.get(anchor.as_str())
                .and_then(|v| v.get("pos"))
                .and_then(|v| v.as_i64())
                .expect("Anchor or pos not found") as i32
        })
        .collect();

    // Ensure poslist is sorted
    let mut zipped: Vec<(&String, i32)> = sp_anchors.iter().zip(poslist.iter()).map(|(a, p)| (*a, *p)).collect();    zipped.sort_by_key(|&(_, pos)| pos);

    let (sorted_anchors, sorted_poslist): (Vec<_>, Vec<_>) = zipped.into_iter().unzip();

    // Find start_index: last index where pos <= start_pos
    let start_index = match sorted_poslist.binary_search(&start_pos) {
        Ok(idx) => idx,
        Err(idx) => if idx == 0 { 0 } else { idx - 1 },
    };

    // Find end_index: first index where pos >= end_pos
    let end_index = match sorted_poslist.binary_search(&end_pos) {
        Ok(idx) => idx,
        Err(idx) => idx,
    };

    // Return the corresponding anchors
    Some((
        sorted_anchors[start_index].clone(),
        sorted_anchors[end_index].clone(),
    ))
}

fn get_graph_intervals(graph:&GraphicalGenome, length: i64) -> HashMap<&String, (i64, i64)>{
    let mut graph_intervals_dict = HashMap::new();
    for edge in graph.edges.keys() {
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


fn get_path(
    target_interval: (i64, i64),
    graph_intervals_dict: &HashMap<&String, (i64, i64)>,
    graph: &GraphicalGenome,
) -> HashMap<String, Vec<String>> {
    // Find overlapping edges
    let mut overlap_edgelist = Vec::new();
    let t_I = Interval::open(target_interval.0, target_interval.1).unwrap();
    for (edge, edge_intervals) in graph_intervals_dict {
        let e_I = Interval::open(edge_intervals.0, edge_intervals.1).unwrap();
        if e_I.intersect(t_I).is_some() {
            overlap_edgelist.push(*edge);
        }
    }

    // Build path dictionary
    let mut path_dict: HashMap<String, Vec<String>> = HashMap::new();
    for edge in &overlap_edgelist {
        if let Some(reads) = graph.edges[*edge].get("reads").and_then(|v| v.as_array()) {
            for read in reads {
                if let Some(read_str) = read.as_str() {
                    let edge_src = graph.edges[*edge].get("src").unwrap().as_array().unwrap()[0].as_str().unwrap();
                    path_dict
                        .entry(read_str.to_string())
                        .or_insert_with(Vec::new)
                        .push(edge_src.to_string());
                    path_dict
                        .entry(read_str.to_string())
                        .or_insert_with(Vec::new)
                        .push(edge.to_string());
                }
            }
        }
    }

    // Group reads by sorted, unique edge paths
    let mut path_dict_final: HashMap<String, Vec<String>> = HashMap::new();
    for (read, edge_list) in path_dict {
        if !edge_list.is_empty() {
            // Create sorted, unique edge path
            let mut unique_edges: Vec<String> = edge_list.into_iter().collect::<HashSet<_>>().into_iter().collect();
            unique_edges.sort();
            let path_key = unique_edges.join(">");
            
            path_dict_final
                .entry(path_key)
                .or_insert_with(Vec::new)
                .push(read);
        }
    }

    path_dict_final
}


fn dfs(
    haplotype_dict: &HashMap<i64, HashMap<String, Vec<String>>>,
    min_reads: usize,
) -> Vec<(Vec<String>, HashSet<String>)> {
    let mut results = Vec::new();
    let mut positions: Vec<i64> = haplotype_dict.keys().cloned().collect();
    positions.sort(); // Ensure consistent order

    fn dfs_helper(
        pos_idx: usize,
        positions: &Vec<i64>,
        haplotype_dict: &HashMap<i64, HashMap<String, Vec<String>>>,
        current_paths: Vec<String>,
        current_reads: HashSet<String>,
        min_reads: usize,
        results: &mut Vec<(Vec<String>, HashSet<String>)>,
    ) {
        if pos_idx == positions.len() {
            if current_reads.len() >= min_reads {
                results.push((current_paths, current_reads));
            }
            return;
        }
        let pos = positions[pos_idx];
        if let Some(path_map) = haplotype_dict.get(&pos) {
            for (path_name, reads_vec) in path_map {
                let reads_set: HashSet<String> = reads_vec.iter().cloned().collect();
                if pos_idx == 0 {
                    dfs_helper(
                        pos_idx + 1,
                        positions,
                        haplotype_dict,
                        vec![path_name.clone()],
                        reads_set,
                        min_reads,
                        results,
                    );
                } else {
                    let shared_reads: HashSet<String> =
                        current_reads.intersection(&reads_set).cloned().collect();
                    if shared_reads.len() >= min_reads {
                        let mut new_paths = current_paths.clone();
                        new_paths.push(path_name.clone());
                        dfs_helper(
                            pos_idx + 1,
                            positions,
                            haplotype_dict,
                            new_paths,
                            reads_set,
                            min_reads,
                            results,
                        );
                    }
                }
            }
        }
    }

    dfs_helper(
        0,
        &positions,
        haplotype_dict,
        Vec::new(),
        HashSet::new(),
        min_reads,
        &mut results,
    );
    results
}

fn find_shared_anchors(path0: &str, path1: &str) -> HashSet<String> {
    // Split paths and filter for anchors starting with "A"
    let anchors_items: HashSet<String> = path0
        .split(">")
        .filter(|p| p.starts_with("A"))
        .map(|s| s.to_string())
        .collect();
    
    let next_anchors: HashSet<String> = path1
        .split(">")
        .filter(|p| p.starts_with("A"))
        .map(|s| s.to_string())
        .collect();
    
    // Return intersection of the two sets
    anchors_items.intersection(&next_anchors).cloned().collect()
}

fn get_seq(
    itemlist: &[String],
    graph: &GraphicalGenome,
    first_anchor: &str,
    last_anchor: &str,
) -> Option<String> {
    let mut seq = String::new();
    
    // Get edges from itemlist (items starting with 'E')
    let edgelist: HashSet<&str> = itemlist
        .iter()
        .filter(|item| item.starts_with('E'))
        .map(|s| s.as_str())
        .collect();
    
    let mut anchor = first_anchor;
    
    while anchor != last_anchor {
        // Add anchor sequence
        if let Some(anchor_data) = graph.anchor.get(anchor) {
            if let Some(seq_val) = anchor_data.get("seq") {
                if let Some(seq_str) = seq_val.as_str() {
                    seq.push_str(seq_str);
                } else {
                    println!("Anchor sequence is not a string");
                    return None;
                }
            } else {
                println!("No 'seq' field found in anchor");
                return None;
            }
        } else {
            println!("Anchor not found: {}", anchor);
            return None;
        }
        
        // Find outgoing edge that's in edgelist
        let mut found_edge = false;
        if let Some(outgoing_edges) = graph.outgoing.get(anchor) {
            for edge in outgoing_edges {
                if edgelist.contains(edge.as_str()) {
                    // Add edge sequence
                    if let Some(edge_data) = graph.edges.get(edge) {
                        if let Some(seq_val) = edge_data.get("seq") {
                            if let Some(seq_str) = seq_val.as_str() {
                                seq.push_str(seq_str);
                            }
                        }
                    }
                    
                    // Get destination anchor
                    if let Some(edge_data) = graph.edges.get(edge) {
                        if let Some(dst_array) = edge_data.get("dst").and_then(|v| v.as_array()) {
                            if let Some(dst_val) = dst_array.get(0) {
                                if let Some(dst_str) = dst_val.as_str() {
                                    anchor = dst_str;
                                    found_edge = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        if !found_edge {
            println!("No outgoing edge found from anchor {} in edgelist", anchor);
            return None;
        }
    }
    Some(seq)
}


fn reconstruct_haplotypes(
    final_dict: &Vec<(Vec<String>, HashSet<String>)>, 
    graph: &GraphicalGenome,
) -> Vec<String> {
    let mut reconstructed_haplotypes = Vec::new();

    for (path, _reads) in final_dict {
        let mut haplotype = String::new();
        let mut last_anchor = String::new();

        for i in 0..path.len() - 1 {
            // Split and sort items
            let mut itemlist: Vec<String> = path[i].split('>').map(|s| s.to_string()).collect();
            itemlist.sort();
            
            let first_anchor = if i == 0 {
                itemlist[0].clone()
            } else {
                last_anchor.clone()
            };

            let shared_anchors = find_shared_anchors(&path[i], &path[i + 1]);
            if !shared_anchors.is_empty() {
                let mut shared_sorted: Vec<String> = shared_anchors.into_iter().collect();
                shared_sorted.sort();
                last_anchor = shared_sorted.last().unwrap().clone();
            } else {
                println!("no shared anchors");
                break;
            }

            if let Some(seq) = get_seq(&itemlist, graph, &first_anchor, &last_anchor) {
                haplotype.push_str(&seq);
            } else {
                println!("no seq");
                break;
            }
        }

        // Handle the last segment
        let mut itemlist: Vec<String> = path.last().unwrap().split('>').map(|s| s.to_string()).collect();
        itemlist.sort();
        let anchor_list: Vec<String> = itemlist.iter().filter(|item| item.starts_with('A')).cloned().collect();
        let first_anchor = last_anchor.clone();
        let last_anchor = anchor_list.last().unwrap().clone();

        if let Some(seq) = get_seq(&itemlist, graph, &first_anchor, &last_anchor) {
            haplotype.push_str(&seq);
        } else {
            println!("no seq");
            continue;
        }

        println!("{}", haplotype.len());
        reconstructed_haplotypes.push(haplotype);
    }

    reconstructed_haplotypes
}

fn write_multipleseq_fasta(outputfile: &PathBuf, sequence: Vec<String>, header: &str) -> io::Result<()>{
    let mut file = File::create(outputfile)?;
    for (i, haplotype) in sequence.iter().enumerate(){
        let header = format!(">{} Haplotype {} \n", header, i.to_string());
        file.write_all(header.as_bytes())?;
        let chars_per_line = 60;
        let sequence_len = haplotype.len();
        let full_lines = sequence_len / chars_per_line;
        for i in 0..full_lines {
            let start = i * chars_per_line;
            let end = start + chars_per_line;
            writeln!(file, "{}", &haplotype[start..end])?;
        }

          // Write any remaining characters that didn't make up a full line
        if sequence_len % chars_per_line != 0 {
            writeln!(file, "{}", &haplotype[full_lines * chars_per_line..])?;
        }
        
    }
  
    Ok(())
}


pub fn start (graph_file: &PathBuf, ref_length:i64, bin_size:i32, pad_size:i64, min_read_ratio:f64, count_support:usize, output_file: &PathBuf, header:&str) {
    let graph = GraphicalGenome::load_graph(graph_file).unwrap();
    let readset = find_all_reads(&graph);

    let mut minimal_read_number = readset.len() as f64 * min_read_ratio;
    if minimal_read_number < count_support as f64{
        minimal_read_number = count_support as f64;
    }
    println!("Total read number: {}, Minimal read number: {}", readset.len(), minimal_read_number);
    // let mut sp_anchor = graph.anchor.keys().collect::<Vec<&String>>();
    // sp_anchor.sort();
    let graph_intervals_dict = get_graph_intervals(&graph, ref_length);
    let interval_list: Vec<i64> = (1..ref_length).step_by(bin_size as usize).map(|x| x as i64).collect();
    let mut haplotype_dict: HashMap<i64, HashMap<String, Vec<String>>> = HashMap::new();
    let mut failed_interval = 0;
    for i in 1..interval_list.len()-1 {
        let startpos = interval_list[i] - pad_size;
        let endpos = interval_list[i+1] + pad_size;
        let target_interval = (startpos, endpos);
        let path_dict = get_path(target_interval, &graph_intervals_dict, &graph);
        let mut survival_path_dict = HashMap::new();
        for (key, value) in path_dict {
            if value.len() > minimal_read_number as usize {
                survival_path_dict.insert(key, value);
            }
        }
        println!("{}, Number of paths:{}", startpos, survival_path_dict.len());
        if survival_path_dict.len() == 0 {
            failed_interval += 1;
        }
        haplotype_dict.insert(startpos, survival_path_dict);
    }
    println!("Number of failed intervals: {}", failed_interval);
    if failed_interval > 0 {
        println!("Failed to assemble haplotypes for {} intervals", failed_interval);
        return;
    }


    let results = dfs(&haplotype_dict, minimal_read_number as usize);
    let haplotypes = reconstruct_haplotypes(&results, &graph);
    let _ =write_multipleseq_fasta(output_file, haplotypes, header).unwrap();
    // write fasta
    // write_fasta(output_file, haplotype, header).unwrap();

}