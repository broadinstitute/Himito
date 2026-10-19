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


#[derive(Debug, Clone)]
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

/// Edge endpoints stored in PG:J as `"src":["A…"]` / `"dst":["A…"]` (see `build::write_gfa`).
fn edge_endpoints(edge_data: &Value) -> Option<(String, String)> {
    let src = edge_data
        .get("src")?
        .as_array()?
        .first()?
        .as_str()?
        .to_string();
    let dst = edge_data
        .get("dst")?
        .as_array()?
        .first()?
        .as_str()?
        .to_string();
    Some((src, dst))
}

/// Ensure every edge has incoming/outgoing adjacency, recovering from PG:J when L lines were missing.
fn ensure_edge_topology(
    edges: &HashMap<String, Value>,
    outgoing: &mut HashMap<String, Vec<String>>,
    incoming: &mut HashMap<String, Vec<String>>,
) {
    for (edge, data) in edges {
        if incoming.contains_key(edge) && outgoing.contains_key(edge) {
            continue;
        }
        let Some((src, dst)) = edge_endpoints(data) else {
            continue;
        };
        if !incoming.contains_key(edge) {
            incoming.insert(edge.clone(), vec![src.clone()]);
            add_unique(outgoing.entry(src).or_default(), edge.clone());
        }
        if !outgoing.contains_key(edge) {
            outgoing.insert(edge.clone(), vec![dst.clone()]);
            add_unique(incoming.entry(dst).or_default(), edge.clone());
        }
    }
}


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
        // Truncated GFAs (e.g. write failed after S lines, before L) still carry src/dst in PG:J.
        ensure_edge_topology(&edge_dict, &mut outgoing, &mut incoming);
        Ok(GraphicalGenome {
            anchor: anchor_dict,
            edges: edge_dict,
            outgoing: outgoing,
            incoming: incoming,
        })
    }

 
    /// The graphical genome of a single sample.
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

        // be careful about the anchor and read index, make sure it matches with the imputed data matrix
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
        let mut visited = HashSet::new();
        finder.find_path(graph, start, end, &Vec::new(), 0, read_sets, &mut visited, 1000);
        finder
    }

    pub fn find_path(
        &mut self, 
        g: &GraphicalGenome, 
        start: &str, 
        end: &str, 
        sofar: &Vec<String>, 
        depth: usize, 
        readset: HashSet<String>,
        visited: &mut HashSet<String>,
        max_depth: usize,
    ) {
        // Check depth limit to prevent stack overflow
        if depth > max_depth {
            return;
        }

        // Cycle detection: if we've visited this node in the current path, skip it
        if visited.contains(start) {
            return;
        }

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

        // Mark current node as visited
        visited.insert(start.to_string());
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
                self.find_path(g, dst, end, &sofar1, depth1, readset1, visited, max_depth);
            }
        }

        // Unmark current node when backtracking
        visited.remove(start);
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
        let (src, dst) = edge_endpoints(edge_data)
            .or_else(|| {
                let src = graph.incoming.get(*edge)?.first()?.clone();
                let dst = graph.outgoing.get(*edge)?.first()?.clone();
                Some((src, dst))
            })
            .unwrap_or_else(|| {
                panic!(
                    "edge {} missing src/dst (no PG:J endpoints and no L-derived adjacency)",
                    edge
                )
            });
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
    // Links before edge segments: edge S lines can be huge (embedded read lists).
    for l in link_output {
        writeln!(file, "{}", l)?;
    }
    for s in edge_output {
        writeln!(file, "{}", s)?;
    }
    Ok(())
}

pub fn edit_distance(cigar: &String) -> usize {
    let mut edit_distance = 0;
    let mut num: usize = 0;
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + ch.to_digit(10).unwrap() as usize;
            continue;
        }
        // ch is an op
        match ch {
            // count edits
            'I' | 'D' | 'X' | 'N' => {
                edit_distance += num;
            }
            'M' | '=' | 'S' | 'H' | 'P' => {}
            _ => {}
        }
        num = 0;
    }
    edit_distance
}

pub fn find_most_supported_edge(graph: &GraphicalGenome, src: String) -> String {
    let empty_vec = Vec::new();
    let outgoinglist = &graph.outgoing.get(&src).unwrap_or(&empty_vec);
    let mut m = 0;
    let mut most_supported_edge = "".to_string();
    let mut nm_value = usize::MAX;
    for edge in outgoinglist.iter(){
        let edge_dst = graph.edges[edge]["dst"].as_array().unwrap().first().unwrap().as_str().unwrap().to_string();
        if edge_dst == "SINK".to_string(){
            continue
        }
        let read_count = graph.edges[edge]["reads"].as_array().unwrap_or(&Vec::new()).len();
        let cigar = graph.edges[edge]["variants"].as_str().unwrap_or("");
        if read_count > m {
            m = read_count;
            most_supported_edge = edge.clone();
            nm_value = edit_distance(&cigar.to_string());
        } else if read_count == m {
            /// define lexicographical order
            if edge.to_string() < most_supported_edge.to_string() {
                most_supported_edge = edge.clone();
            }else{
                continue
            }
        }else{
            continue
        }
    }
    most_supported_edge

}

pub fn construct_major_haplotype(graph:&GraphicalGenome) -> String {
    let mut anchorlist: Vec<_> = graph.anchor.keys().collect();
    anchorlist.sort();
    let mut src = anchorlist.first().unwrap().to_string();
    let mut next_edge = "".to_string();
    let dst = anchorlist.last().unwrap().to_string();
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
    let dst = anchorlist.last().unwrap().to_string();
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// GFA with edge src/dst in PG:J but no L lines (truncated write / ignored IO error).
    fn gfa_edges_without_links(path: &std::path::Path) {
        let mut f = File::create(path).unwrap();
        writeln!(f, "H\tVN:Z:1.0").unwrap();
        writeln!(f, "S\tA000001\tACGT\tPG:J:{{\"pos\":1}}").unwrap();
        writeln!(f, "S\tA000002\tTGCA\tPG:J:{{\"pos\":10}}").unwrap();
        writeln!(
            f,
            "S\tE00000001\tNNNN\tPG:J:{{\"src\":[\"A000001\"],\"dst\":[\"A000002\"],\"reads\":[\"r1\"]}}\tRC:i:1"
        )
        .unwrap();
    }

    #[test]
    fn load_graph_recovers_topology_from_edge_json_without_l_lines() {
        let dir = std::env::temp_dir().join(format!("himito_agg_nolinks_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("nolinks.gfa");
        gfa_edges_without_links(&path);

        let graph = GraphicalGenome::load_graph(&path).unwrap();
        assert!(graph.edges.contains_key("E00000001"));
        assert_eq!(graph.incoming["E00000001"], vec!["A000001".to_string()]);
        assert_eq!(graph.outgoing["E00000001"], vec!["A000002".to_string()]);
        assert_eq!(graph.outgoing["A000001"], vec!["E00000001".to_string()]);
        assert_eq!(graph.incoming["A000002"], vec!["E00000001".to_string()]);
    }

    #[test]
    fn write_graph_from_graph_uses_edge_json_when_incoming_missing() {
        let dir = std::env::temp_dir().join(format!("himito_agg_write_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("out.gfa");

        let mut edges = HashMap::new();
        edges.insert(
            "E00000001".to_string(),
            serde_json::json!({
                "seq": "NNNN",
                "src": ["A000001"],
                "dst": ["A000002"],
                "reads": ["r1"]
            }),
        );
        let mut anchor = HashMap::new();
        anchor.insert("A000001".to_string(), serde_json::json!({"seq": "ACGT", "pos": 1}));
        anchor.insert("A000002".to_string(), serde_json::json!({"seq": "TGCA", "pos": 10}));

        let graph = GraphicalGenome {
            anchor,
            edges,
            // Simulate load after a truncated GFA with no L lines.
            outgoing: HashMap::new(),
            incoming: HashMap::new(),
        };

        write_graph_from_graph(out.to_str().unwrap(), &graph).unwrap();
        let text = std::fs::read_to_string(&out).unwrap();
        assert!(
            text.contains("L\tA000001\t+\tE00000001\t+\t0M"),
            "missing src->edge link:\n{text}"
        );
        assert!(
            text.contains("L\tE00000001\t+\tA000002\t+\t0M"),
            "missing edge->dst link:\n{text}"
        );
    }
}