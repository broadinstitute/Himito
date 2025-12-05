use core::panic;
use flate2::read::GzDecoder;
use ndarray::Array2;
use serde_json::Value;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::{path::PathBuf};
use std::io::Write;


#[derive(Debug)]
pub struct GraphicalGenome {
    pub anchor: HashMap<String, Value>,
    pub edges: HashMap<String, Value>,
    pub outgoing: HashMap<String, Vec<String>>,
    pub incoming: HashMap<String, Vec<String>>,
}

pub fn add_unique(vec: &mut Vec<String>, item: String) {
    if !vec.contains(&item) {
        vec.push(item);
    }
}

/// Reverse complement of a k-mer
///
/// # Arguments
///
/// * `kmer` - A string representing the k-mer.
///
/// # Returns
///
/// A string containing the reverse complement of the k-mer.
///
/// # Panics
///
/// If the k-mer contains an unexpected character.
#[must_use]
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

impl GraphicalGenome {
    /// Load a graphical genome from a file
    ///
    /// # Arguments
    ///
    /// * `filename` - A string slice that holds the name of the file to be loaded
    ///
    /// # Returns
    ///
    /// A Result containing the graphical genome if the file was successfully loaded
    ///
    /// # Errors
    ///
    /// * The file cannot be opened (e.g., permission issues, file not found).
    /// * The file format is invalid (e.g., not a valid GZIP file or has unexpected content).
    ///
    /// # Panics
    ///
    /// - The file cannot be opened.
    /// - The file is not a valid GZIP-compressed file.
    /// - The annotation string is not valid JSON.
    /// - There is a memory allocation error or other unexpected issue with the `HashMap` operations.
    pub fn load_graph(filename: &PathBuf) -> io::Result<GraphicalGenome> {
        let file = File::open(filename)?;
        let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        let mut anchor_dict = HashMap::new();
        let mut edge_dict = HashMap::new();
        let mut outgoing: HashMap<String, Vec<String>> = HashMap::new();
        let mut incoming: HashMap<String, Vec<String>> = HashMap::new();
        for line in reader.lines() {
            let line = line?.trim_end().to_string();
            if line.starts_with('S') {
                let itemlist: Vec<&str> = line.split('\t').collect();
                if itemlist.len() < 4 {
                    println!("Invalid line: {}", line);
                    continue;
                }
                let name = itemlist[1];
                let seq = itemlist[2];
                let annotation = &itemlist[3][5..];

                let json_value: Value = serde_json::from_str(annotation).unwrap();

                let mut value = json_value;

                value["seq"] = Value::String(seq.to_string());

                if name.starts_with('A') {
                    anchor_dict.insert(name.to_string(), value);
                } else if name.starts_with('E') {
                    edge_dict.insert(name.to_string(), value);
                }
            } else if line.starts_with('L') {
                let itemlist: Vec<&str> = line.split('\t').collect();
                let src = itemlist[1];
                let dst = itemlist[3];
                outgoing
                    .entry(src.to_string())
                    .or_default()
                    .push(dst.to_string());
                incoming
                    .entry(dst.to_string())
                    .or_default()
                    .push(src.to_string());
            }
        }
        Ok(GraphicalGenome {
            anchor: anchor_dict,
            edges: edge_dict,
            outgoing: outgoing,
            incoming: incoming,
        })
    }

    /// Method to extract a single sample graph
    ///
    /// # Arguments
    ///
    /// * `df_single_sample` - A 2D array containing the imputed data matrix
    /// * `anchorlist` - A vector of strings containing the anchor names
    /// * `readset` - A vector of strings containing the read names
    /// * `sample` - A string slice containing the sample name
    ///
    /// # Returns
    ///
    /// A Result containing the graphical genome of the single sample
    ///
    /// # Errors
    ///
    /// * The `df_single_sample` array has an invalid shape.
    /// * The anchor or read indices are invalid.
    /// * There is a memory allocation error.
    /// * The JSON data is invalid.
    /// * There is an unexpected data structure or other internal error.
    ///
    /// # Panics
    ///
    /// If the edge do not have outgoing anchor
    pub fn extract_single_sample_graph(
        &self,
        df_single_sample: &Array2<f64>, 
        anchorlist: &Vec<String>,
        readset:&Vec<String>,
        sample: &str,
    ) -> io::Result<GraphicalGenome> {

        let mut new_edges: HashMap<String, serde_json::Value> = HashMap::new();
        let mut new_incoming = HashMap::new();
        let mut new_outgoing = HashMap::new();

        // be careful about the anchor and read index, 
        // make sure it matches with the imputed data matrix
        for (anchorindex, anchor) in anchorlist.iter().enumerate(){
            let mut d: HashMap<usize, String> = HashMap::new();
            if let Some(outgoinglist) = self.outgoing.get(anchor) {
                // Update the existing `d` variable instead of declaring a new one
                d = outgoinglist
                    .iter()
                    .enumerate()
                    .map(|(index, value)| (index, value.to_string()))
                    .collect();
            }

            for (read_index, read) in readset.iter().enumerate(){
                let edgeindex = df_single_sample[[read_index, anchorindex]].round() - 1.0;
                #[allow(clippy::cast_possible_truncation)]#[allow(clippy::cast_sign_loss)]
                let usize_index = edgeindex as usize;
                let edgename = d.get(&usize_index);
                // println!("usize_index: {}, d: {:?}, edgename: {:?}", usize_index, &d, edgename); // Assuming new_edges is a HashMap<String, SomeType>
                new_edges.entry(edgename.unwrap().to_string()).or_insert_with(|| serde_json::json!({}));
                if let Some(edge_value) = self.edges.get(&(*edgename.as_ref().unwrap()).to_string()) {
                    if let Some(seq_value) = edge_value.get("seq") {
                        let new_edge_value = new_edges.entry((*edgename.as_ref().unwrap()).to_string()).or_insert_with(|| serde_json::json!({}));
                        new_edge_value["seq"] = seq_value.clone();
                    }
                }

                let edgename_str = (*edgename.as_ref().unwrap()).to_string();
                let new_edge_value = new_edges.entry(edgename_str.clone()).or_insert_with(|| serde_json::json!({}));
                if !new_edge_value.get("reads").is_some() {
                    new_edge_value["reads"] = serde_json::Value::Array(vec![]);
                }
                let reads_vec = new_edge_value.get_mut("reads").unwrap().as_array_mut().unwrap();
                reads_vec.push(serde_json::Value::String(read.to_string()));


                // let mut new_edge_value = new_edges.entry(edgename_str.clone()).or_insert_with(|| serde_json::json!({}));
                if !new_edge_value.get("strain").is_some() {
                    new_edge_value["strain"] = serde_json::Value::Array(vec![]);
                }
                let strain_vec = new_edge_value.get_mut("strain").unwrap().as_array_mut().unwrap();
                                if !strain_vec.iter().any(|x| x == &serde_json::Value::String(sample.to_string())) {
                    strain_vec.push(serde_json::Value::String(sample.to_string()));
                }


                if let Some(outgoing_list) = self.outgoing.get(&(*edgename.as_ref().unwrap()).to_string()) {
                    if let Some(dst) = outgoing_list.get(0){
                        let incoming_list = new_incoming.entry((*edgename.as_ref().unwrap()).to_string()).or_default();
                        add_unique(incoming_list, anchor.to_string());

                        let incoming_dst_list = new_incoming.entry(dst.to_string()).or_default();
                        add_unique(incoming_dst_list, (*edgename.as_ref().unwrap()).to_string());
                        let outgoing_list = new_outgoing.entry(anchor.to_string()).or_default();
                        add_unique(outgoing_list, (*edgename.as_ref().unwrap()).to_string());
                        let outgoing_edgename_list = new_outgoing.entry((*edgename.as_ref().unwrap()).to_string()).or_default();
                        add_unique(outgoing_edgename_list, dst.to_string());
                    }
                    else{
                        panic!("edge do not have outgoing anchor")
                            
                        }

                }
            }

        }
    let new_graph = GraphicalGenome {
            anchor: self.anchor.clone(), 
            edges: new_edges,
            outgoing: new_outgoing,
            incoming: new_incoming,
            };
    Ok(new_graph) 
    }
}


pub struct FindAllPathBetweenAnchors {
    pub subpath: Vec<(Vec<String>, HashSet<String>)>,
}

impl FindAllPathBetweenAnchors {
    #[must_use]
    pub fn new(graph: &GraphicalGenome, start: &str, end: &str, read_sets: HashSet<String>) -> Self {
        let mut finder = FindAllPathBetweenAnchors {
            subpath: Vec::new(),
        };
        finder.find_path(graph, start, end, &Vec::new(), 0, read_sets);
        finder
    }

    pub fn find_path(&mut self, g: &GraphicalGenome, start: &str, end: &str, sofar: &Vec<String>, depth: usize, readset: HashSet<String>) {
        if start == end {
            let mut sofar1 = sofar.clone();
            sofar1.push(end.to_string());
            if !readset.is_empty() {
                self.subpath.push((sofar1, readset));
            }
            return;
        }

        if readset.is_empty() || start == "SINK" {
            return;
        }

        if !g.outgoing.contains_key(start) {
            return;
        }

        let depth1 = depth + 1;

        if let Some(outgoing) = g.outgoing.get(start) {
            
            for dst in outgoing {
                let mut readset1 = readset.clone();
                if dst.starts_with("E") {
                    if let Some(edge_reads) = g.edges.get(dst).and_then(|e| e.get("reads").and_then(|r| r.as_array())) {
                        readset1.retain(|read| edge_reads.iter().any(|r| r.as_str() == Some(read)));
                    }
                    
                }
                let mut sofar1 = sofar.clone();
                sofar1.push(start.to_string());
                self.find_path(g, dst, end, &sofar1, depth1, readset1);
            }
        }
    }
}

// Series parallele graph
#[must_use]
pub fn reconstruct_path_seq(graph: &GraphicalGenome, path: &[String]) -> String {
    let mut seq = String::new();
    for item in path {
        if item.starts_with('A') {
            if let Some(anchor) = graph.anchor.get(item) {
                seq += &anchor["seq"].as_str().unwrap_or_default(); // Assuming `anchor` is a HashMap and "seq" is a key
                // println!("{:?}", anchor["seq"].as_str().unwrap_or_default());
            }
        } else if item.starts_with("E") {
            if let Some(edge) = graph.edges.get(item) {
                seq += &edge["seq"].as_str().unwrap_or_default(); // Assuming `edges` is a HashMap and "seq" is a key
            }
        }
    }
    seq
}

pub fn write_graph_from_graph(filename: &str, graph: &GraphicalGenome) -> std::io::Result<()> {
    let mut file = File::create(filename)?;

    writeln!(file, "H\tVN:Z:1.0")?;

    let mut keys: Vec<_> = graph.anchor.keys().collect();
    keys.sort();
    for anchor in keys.iter() {
        let data = &graph.anchor[*anchor];
        let seq = data["seq"].as_str().unwrap_or_default();
        let mut data_clone = data.clone();
        data_clone.as_object_mut().unwrap().remove("seq");
        let json_string = serde_json::to_string(&data_clone).unwrap_or_else(|_| "{}".to_string());
        writeln!(file, "S\t{}\t{}\tPG:J:{}", anchor, seq, json_string)?;
    }
    let mut edge_output = Vec::new();
    let mut link_output = Vec::new();

    let mut edge_keys: Vec<_> = graph.edges.keys().collect();
    edge_keys.sort();
    for edge in edge_keys.iter() {
        let edge_data = &graph.edges[*edge];
        let seq = edge_data["seq"].as_str().unwrap_or_default();
        let src = graph.incoming[*edge][0].clone();
        let dst = graph.outgoing[*edge][0].clone();
        let mut edge_data_clone = edge_data.clone();
        edge_data_clone.as_object_mut().unwrap().remove("seq");
        let json_string =
            serde_json::to_string(&edge_data_clone).unwrap_or_else(|_| "{}".to_string());
        let formatted_string = if edge_data
            .get("reads")
            .and_then(|r| r.as_array())
            .map_or(false, |arr| !arr.is_empty())
        {
            format!(
                "S\t{}\t{}\tPG:J:{}\tRC:i:{}",
                edge,
                seq,
                json_string,
                edge_data
                    .get("reads")
                    .and_then(|r| r.as_array())
                    .map_or(0, |arr| arr.len())
            )
        } else {
            format!("S\t{}\t{}", edge, seq)
        };
        edge_output.push(formatted_string);

        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", src, edge));
        link_output.push(format!("L\t{}\t+\t{}\t+\t0M", edge, dst));
    }
    for s in edge_output {
        writeln!(file, "{}", s)?;
    }
    for l in link_output {
        writeln!(file, "{}", l)?;
    }
    Ok(())
}

pub fn find_most_supported_edge(graph: &GraphicalGenome, src: String) -> String {
    let empty_vec = Vec::new();
    let outgoinglist = &graph.outgoing.get(&src).unwrap_or(&empty_vec);
    let mut m = 0;
    let mut most_supported_edge = "".to_string();
    for edge in outgoinglist.iter(){
        let edge_dst = graph.edges[edge]["dst"].as_array().unwrap().first().unwrap().as_str().unwrap().to_string();
        if edge_dst == "SINK".to_string(){
            continue
        }
        let read_count = graph.edges[edge]["reads"].as_array().unwrap_or(&Vec::new()).len();
        if read_count > m {
            m = read_count;
            most_supported_edge = edge.clone();
        }
    }
    most_supported_edge

}

pub fn construct_major_haplotype(graph:&GraphicalGenome) -> String {
    let mut anchorlist: Vec<_> = graph.anchor.keys().collect();
    anchorlist.sort();
    let mut src = anchorlist.first().unwrap().to_string();
    let mut next_edge = "".to_string();
    let mut dst = anchorlist.last().unwrap().to_string();
    let mut haplotype = String::new();
    // println!("{}", src);
    let mut entity_set = HashSet::new();

    while dst > src {
        next_edge = find_most_supported_edge(&graph, src.clone().to_string());
        if next_edge == "".to_string(){
            break
        }
        if entity_set.contains(&src) {
            break
        }
        entity_set.insert(src.clone());
        let anchor_seq = graph.anchor.get(&src)
            .and_then(|v| v.get("seq"))
            .and_then(|v| v.as_str())
            .unwrap_or(""); 
        let edge_seq = graph.edges.get(&next_edge)
            .and_then(|v| v.get("seq"))
            .and_then(|v| v.as_str())
            .unwrap_or("");       
        haplotype.push_str(anchor_seq);
        haplotype.push_str(edge_seq);
        src = graph.edges[&next_edge]["dst"].as_array()
            .unwrap()
            .first()
            .unwrap()
            .as_str()
            .unwrap()
            .to_string();        
        // println!("{}, {}, {}", src, dst, next_edge);
    }

    next_edge = find_most_supported_edge(&graph, src.clone().to_string());
    if next_edge == "".to_string(){
        return haplotype
    }
    let anchor_seq = graph.anchor.get(&src)
        .and_then(|v| v.get("seq"))
        .and_then(|v| v.as_str())
        .unwrap_or(""); 
    let edge_seq = graph.edges.get(&next_edge)
        .and_then(|v| v.get("seq"))
        .and_then(|v| v.as_str())
        .unwrap_or("");
    haplotype.push_str(anchor_seq);
    haplotype.push_str(edge_seq);
    haplotype

}

pub fn construct_major_haplotype_entitylist(graph:&GraphicalGenome) -> Vec<String> {
    let mut anchorlist: Vec<_> = graph.anchor.keys().collect();
    anchorlist.sort();
    let mut src = anchorlist.first().unwrap().to_string();
    let mut next_edge = "".to_string();
    let mut dst = anchorlist.last().unwrap().to_string();
    let mut haplotype = Vec::new();
    // println!("{}", src);
    let mut entity_set = HashSet::new();

    while dst > src {
        next_edge = find_most_supported_edge(&graph, src.clone().to_string());
        if next_edge == "".to_string(){
            break
        }
        if entity_set.contains(&src) {
            break
        }
        entity_set.insert(src.clone());
    
        haplotype.push(src.clone());
        haplotype.push(next_edge.clone());
        src = graph.edges[&next_edge]["dst"].as_array()
            .unwrap()
            .first()
            .unwrap()
            .as_str()
            .unwrap()
            .to_string();        
        // println!("{}, {}, {}", src, dst, next_edge);
    }

    next_edge = find_most_supported_edge(&graph, src.clone().to_string());
    if next_edge == "".to_string(){
        return haplotype
    }
    haplotype.push(src.clone());
    haplotype.push(next_edge.clone());
    haplotype

}