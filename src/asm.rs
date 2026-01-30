use crate::agg::*;
use std::{path::PathBuf, fs::File, io::{self, Write}};



fn write_fasta(outputfile: &PathBuf, sequence: String, header: &str) -> io::Result<()>{
    let index = 0;
    let mut file = File::create(outputfile)?;

    let header = format!(">{} \n", header);
    file.write_all(header.as_bytes())?;

    let chars_per_line = 60;
    let sequence_len = sequence.len();
    let full_lines = sequence_len / chars_per_line;

    for i in 0..full_lines {
        let start = i * chars_per_line;
        let end = start + chars_per_line;
        writeln!(file, "{}", &sequence[start..end])?;
    }

    // Write any remaining characters that didn't make up a full line
    if sequence_len % chars_per_line != 0 {
        writeln!(file, "{}", &sequence[full_lines * chars_per_line..])?;
    }

    Ok(())
}

pub fn start (graph_file: &PathBuf, output_file: &PathBuf, header:&str) {
    let graph = GraphicalGenome::load_graph(graph_file).unwrap();
    let haplotype = construct_major_haplotype(&graph);
    println!("{}", haplotype.len());

    // write fasta
    write_fasta(output_file, haplotype, header).unwrap();

}