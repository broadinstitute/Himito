// correct graph based on srWGS data
use std::process::Command;
use msbwt2::dynamic_bwt::{create_from_fastx,DynamicBWT};
use msbwt2::msbwt_core::BWT;
use msbwt2::string_util;
use msbwt2::rle_bwt::RleBWT;

use crate::{agg::*};
use std::{path::PathBuf, fs::File, io::{self, Write}};
use std::collections::{HashMap, HashSet};



pub fn start (graph_file: &PathBuf, bwt_file: &String, output_file: &PathBuf, query_length:usize, min_support_counts: usize) {
    println!("Correcting graph based on srWGS data");
    
    // construct msbwt
    println!("Loading msbwt");
    let mut bwt = RleBWT::new();
    let _ = bwt.load_numpy_file(&bwt_file);
    
    // load graph
    let mut graph = GraphicalGenome::load_graph(graph_file).unwrap();
    let mut discarded_edges: Vec<String> = Vec::new();

    // kmerize graph and check sr support counts
    println!("Kmerizing graph and checking sr support counts");
    let edgelist = graph.edges.keys().collect::<Vec<&String>>();
    let mut kmer_graph: HashMap<String, String> = HashMap::new(); // kmer_graph[kmer] = edge_id
    for edge in edgelist {
        let mut status: bool = false; // not discarded
        let edge_data = graph.edges.get(edge).unwrap();
        let src_list = edge_data.get("src").unwrap().as_array().expect("src is not an array");
        let src = src_list[0].as_str().expect("src[0] is not a string").to_string();
        let mut src_seq = "".to_string();

        if src == "SOURCE" {
            src_seq = "".to_string();
        }else if src == "SINK" {
            src_seq = "".to_string();
        }else{
            let src_anchor_data = graph.anchor.get(&src).expect(&format!("Anchor not found for src: {}", src));
            let src_seq_val = src_anchor_data.get("seq").expect(&format!("'seq' not found for anchor: {}", src));
            src_seq = src_seq_val.as_str().expect("src_seq_val is not a string").to_string();
        }
        
        let edge_seq = edge_data.get("seq").unwrap().as_str().expect("edge seq is not a string").to_string();
        let query_seq = src_seq + &edge_seq;   

        if query_seq.len() >= query_length{
            for i in 0..query_seq.len() - query_length + 1{
                let kmer = query_seq[i..i+query_length].to_string();
                let rev_kmer = reverse_complement(&kmer);
                let count = bwt.count_kmer(&string_util::convert_stoi(&kmer)) + bwt.count_kmer(&string_util::convert_stoi(&rev_kmer));
                if count < min_support_counts as u64{
                    status = true; //discarded
                    break;
                }
            }
        }else {
            let reverse_seq = reverse_complement(&query_seq);
            let count = bwt.count_kmer(&string_util::convert_stoi(&query_seq)) + bwt.count_kmer(&string_util::convert_stoi(&reverse_seq));
            if count < min_support_counts as u64{
                status = true; //discarded
            }
        }

        if status == true{
            discarded_edges.push(edge.clone());
        }
    }

    // remove discarded edges from graph
    println!("Removing discarded edges from graph");
    println!("Number of discarded edges: {}", discarded_edges.len());
    //compute total number of bases in discarded edges
    let mut total_bp = 0;
    for edge in discarded_edges.clone(){
        let edge_data = graph.edges.get(&edge).unwrap();
        let edge_seq = edge_data.get("seq").unwrap().as_str().expect("edge seq is not a string").to_string();
        total_bp += edge_seq.len();
    }
    println!("Total number of bases in discarded edges: {}", total_bp);
    // remove from graph
    for dedge in discarded_edges {
        graph.edges.remove(&dedge);
    }

    // write graph to gfa file
    let _ = write_graph_from_graph(output_file.to_str().unwrap(), &graph);

}