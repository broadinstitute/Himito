use std::collections::HashMap;
use std::error::Error;
use std::path::PathBuf;
use ndarray::{Array2};
use crate::call;
use speedytree::DistanceMatrix;
use speedytree::robinson_foulds;
use speedytree::{NeighborJoiningSolver, Canonical, RapidBtrees, Hybrid, Tree, to_newick};
use rayon;
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
    let mut distance_matrix = DistanceMatrix::build(distance_matrix, read_names).unwrap();
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
    let tree = build_phylogenetic_tree(distance_matrix, read_names).unwrap();
    let newick_file = output_file.with_extension("newick");
    let newick_string = to_newick(&tree);
    let _ = fs::write(newick_file.to_str().unwrap(), newick_string);
    
}