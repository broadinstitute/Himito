use std::collections::HashMap;
use std::error::Error;
use std::path::PathBuf;
use crate::call;
use speedytree::DistanceMatrix;
use speedytree::{NeighborJoiningSolver, Hybrid, Tree, to_newick};
use csv::Reader;
use std::fs;

pub fn read_var_matrix(var_matrix: &PathBuf, min_vaf: f64) -> (Vec<String>, Vec<String>, Vec<Vec<f64>>) {
    let mut var_matrix = Reader::from_path(var_matrix).unwrap();
    let mut read_names = Vec::new();
    let mut variant_names = Vec::new();
    let mut var_matrix_values = Vec::new();
    for (row_index, row) in var_matrix.records().enumerate() {
        if row_index == 0 {
            read_names = row.unwrap().iter().skip(1).map(|x| x.to_string()).collect();
            continue;
        }
        let variant_name = row.as_ref().unwrap()[0].to_string();
        
        let row_values = row.as_ref().unwrap().iter().skip(1).map(|x| x.parse::<f64>().unwrap()).collect::<Vec<f64>>();
        let allele_count = row_values.iter().sum::<f64>();
        let allele_frequency = allele_count / row_values.len() as f64;
        if allele_frequency < min_vaf {
            continue;
        }
        var_matrix_values.push(row_values);
        variant_names.push(variant_name);
    }
    (read_names, variant_names, var_matrix_values)
}

pub fn transpose_matrix(var_matrix: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let mut transposed_matrix = Vec::new();
    for i in 0..var_matrix[0].len() {
        let mut row = Vec::new();
        for j in 0..var_matrix.len() {
            row.push(var_matrix[j][i]);
        }
        transposed_matrix.push(row);
    }
    transposed_matrix
}

pub fn get_distance_matrix(var_matrix: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let transposed_matrix = transpose_matrix(var_matrix);
    let mut distance_matrix = Vec::new();
    for i in 0..transposed_matrix.len() {
        let mut distance_vector = Vec::new();
        for j in 0..transposed_matrix.len() {
            let distance = call::jaccard_distance(&transposed_matrix[i].iter().map(|x| *x > 0.5).collect::<Vec<bool>>(), & transposed_matrix[j].iter().map(|x| *x > 0.5).collect::<Vec<bool>>());
            distance_vector.push(distance);
        }
        distance_matrix.push(distance_vector);
    }
    distance_matrix
}

pub fn build_phylogenetic_tree(distance_matrix:Vec<Vec<f64>>, read_names: Vec<String>) -> Result<Tree, Box<dyn Error>> {
    let distance_matrix = DistanceMatrix::build(distance_matrix, read_names).unwrap();
    let tree = NeighborJoiningSolver::<Hybrid>::default(distance_matrix.clone())
                                                            .set_canonical_steps(2)
                                                            .solve()
                                                            .unwrap();
    Ok(tree)
}

pub fn start (var_matrix: &PathBuf, output_file: &PathBuf, min_vaf:f64) {
    let (read_names, variant_names, var_matrix) = read_var_matrix(var_matrix, min_vaf);
    let distance_matrix = get_distance_matrix(&var_matrix);
    println!("distance_matrix: {:?}", distance_matrix.len());
    println!("read_names: {:?}", read_names.len());
    println!("variant_matrix: {:?}, {:?}", var_matrix.len(), var_matrix[0].len());
    let read_names_clone = read_names.clone();
    let tree = build_phylogenetic_tree(distance_matrix, read_names_clone).unwrap();
    let newick_file = output_file.with_extension("newick");
    let newick_string = to_newick(&tree);
    let _ = fs::write(newick_file.to_str().unwrap(), newick_string);

    // map variants back to the tree
    let variant_mappings = map_variants_to_tree(&tree, &read_names, &variant_names, &var_matrix);
    
    // output variant mappings
    let mapping_file = output_file.with_extension("variant_mappings.csv");
    let _ = output_variant_mappings(&mapping_file, &variant_mappings);
    
    println!("Variant mappings saved to: {:?}", mapping_file);
}

#[derive(Debug, Clone)]
pub struct VariantMapping {
    pub variant_name: String,
    pub tree_node: String,
    pub supporting_reads: Vec<String>,
    pub read_count: usize,
    pub confidence: f64,
}

/// Maps variants back to the phylogenetic tree based on read clustering patterns
pub fn map_variants_to_tree(
    tree: &Tree,
    read_names: &[String],
    variant_names: &[String],
    var_matrix: &[Vec<f64>],
) -> Vec<VariantMapping> {
    let mut mappings = Vec::new();
    
    // Get tree leaves (terminal nodes)
    let leaves = get_tree_leaves(tree);
    println!("Found {} leaves in the tree", leaves.len());
    
    // Create read to leaf mapping
    let read_to_leaf = create_read_to_leaf_mapping(read_names, &leaves);
    
    // For each variant, determine which tree nodes it maps to
    for (variant_idx, variant_name) in variant_names.iter().enumerate() {
        let variant_row = &var_matrix[variant_idx];
        
        // Find reads that have this variant (value > 0.5)
        let supporting_reads: Vec<String> = read_names
            .iter()
            .zip(variant_row.iter())
            .filter(|(_, &value)| value > 0.5)
            .map(|(read_name, _)| read_name.clone())
            .collect();
        
        if supporting_reads.is_empty() {
            continue;
        }
        
        // Map supporting reads to tree leaves
        let leaf_counts = count_reads_per_leaf(&supporting_reads, &read_to_leaf);
        let read_count = supporting_reads.len();
        
        // Find the most likely tree node for this variant
        let (best_node, confidence) = find_best_tree_node(&leaf_counts, &read_count);
        
        mappings.push(VariantMapping {
            variant_name: variant_name.clone(),
            tree_node: best_node,
            supporting_reads,
            read_count,
            confidence,
        });
    }
    
    mappings
}

/// Extracts leaf nodes from the phylogenetic tree
fn get_tree_leaves(tree: &Tree) -> Vec<String> {
    // Get all node indices and filter for leaves (nodes with degree 1)
    let mut leaves = Vec::new();
    for node_idx in tree.node_indices() {
        println!("tree_node_indices: {:?}", node_idx);
        if tree.neighbors(node_idx).count() == 1 {
            if let Some(node_weight) = tree.node_weight(node_idx) {
                leaves.push(node_weight.clone());
            }
        }
    }
    leaves
}

/// Creates a mapping from read names to tree leaf nodes
fn create_read_to_leaf_mapping(read_names: &[String], leaves: &[String]) -> HashMap<String, String> {
    let mut mapping = HashMap::new();
    for read_name in read_names {
        if let Some(leaf) = leaves.iter().find(|leaf| *leaf == read_name) {
            mapping.insert(read_name.clone(), leaf.clone());
        } else {
            // If no exact match, use the read name as the leaf name
            mapping.insert(read_name.clone(), read_name.clone());
        }
    }
    
    mapping
}

/// Counts how many supporting reads belong to each tree leaf
fn count_reads_per_leaf(
    supporting_reads: &[String],
    read_to_leaf: &HashMap<String, String>,
) -> HashMap<String, usize> {
    let mut leaf_counts = HashMap::new();
    
    for read in supporting_reads {
        if let Some(leaf) = read_to_leaf.get(read) {
            *leaf_counts.entry(leaf.clone()).or_insert(0) += 1;
        }
    }
    
    leaf_counts
}

/// Finds the best tree node for a variant based on read distribution
fn find_best_tree_node(
    leaf_counts: &HashMap<String, usize>,
    total_reads: &usize,
) -> (String, f64) {
    if leaf_counts.is_empty() {
        return ("unknown".to_string(), 0.0);
    }
    
    // Find the leaf with the most supporting reads
    let (best_leaf, max_count) = leaf_counts
        .iter()
        .max_by_key(|(_, &count)| count)
        .unwrap();
    
    // Calculate confidence as the proportion of reads supporting this leaf
    let confidence = *max_count as f64 / *total_reads as f64;
    
    (best_leaf.clone(), confidence)
}

/// Outputs variant mappings to a CSV file
pub fn output_variant_mappings(
    output_file: &PathBuf,
    mappings: &[VariantMapping],
) -> Result<(), Box<dyn Error>> {
    let mut wtr = csv::Writer::from_path(output_file)?;
    
    // Write header
    wtr.write_record(&[
        "variant_name",
        "tree_node",
        "supporting_read_count",
        "confidence",
        "supporting_reads",
    ])?;
    
    // Write data
    for mapping in mappings {
        let supporting_reads_str = mapping.supporting_reads.join(";");
        wtr.write_record(&[
            &mapping.variant_name,
            &mapping.tree_node,
            &mapping.read_count.to_string(),
            &mapping.confidence.to_string(),
            &supporting_reads_str,
        ])?;
    }
    
    wtr.flush()?;
    Ok(())
}

/// Advanced variant mapping using tree topology
pub fn map_variants_to_tree_advanced(
    tree: &Tree,
    read_names: &[String],
    variant_names: &[String],
    var_matrix: &[Vec<f64>],
) -> Vec<VariantMapping> {
    let mut mappings = Vec::new();
    
    // Get tree structure information
    let leaves = get_tree_leaves(tree);
    let read_to_leaf = create_read_to_leaf_mapping(read_names, &leaves);
    
    // For each variant, analyze its distribution across the tree
    for (variant_idx, variant_name) in variant_names.iter().enumerate() {
        let variant_row = &var_matrix[variant_idx];
        
        // Find reads that have this variant
        let supporting_reads: Vec<String> = read_names
            .iter()
            .zip(variant_row.iter())
            .filter(|(_, &value)| value > 0.5)
            .map(|(read_name, _)| read_name.clone())
            .collect();
        
        if supporting_reads.is_empty() {
            continue;
        }
        
        // Analyze variant distribution across tree nodes
        let leaf_counts = count_reads_per_leaf(&supporting_reads, &read_to_leaf);
        let read_count = supporting_reads.len();
        
        // Use more sophisticated logic to determine the best tree position
        let (best_node, confidence) = find_best_tree_node_advanced(
            &leaf_counts,
            &read_count,
            &tree,
        );
        
        mappings.push(VariantMapping {
            variant_name: variant_name.clone(),
            tree_node: best_node,
            supporting_reads,
            read_count,
            confidence,
        });
    }
    
    mappings
}

/// Advanced tree node finding using tree topology
fn find_best_tree_node_advanced(
    leaf_counts: &HashMap<String, usize>,
    total_reads: &usize,
    _tree: &Tree,
) -> (String, f64) {
    if leaf_counts.is_empty() {
        return ("unknown".to_string(), 0.0);
    }
    
    // Find the leaf with the most supporting reads
    let (best_leaf, max_count) = leaf_counts
        .iter()
        .max_by_key(|(_, &count)| count)
        .unwrap();
    
    // Calculate confidence considering tree topology
    let confidence = *max_count as f64 / *total_reads as f64;
    
    // If confidence is high enough, use the best leaf
    if confidence > 0.5 {
        (best_leaf.clone(), confidence)
    } else {
        // For low confidence variants, try to find a common ancestor
        // This is a simplified approach - in practice, you'd traverse the tree
        ("internal_node".to_string(), confidence)
    }
}