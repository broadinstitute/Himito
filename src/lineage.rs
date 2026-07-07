use std::collections::HashMap;
use anyhow::{bail,Context, Result};
use rust_htslib::bcf::{Read, Reader}; // Read trait must be in scope for .records()
use std::io::{BufWriter, Write};
use std::fs::File;

// ─── Tree data structure ──────────────────────────────────────────────────────

pub struct Node {
    pub id: usize,
    pub label: String,
    pub is_leaf: bool,
    /// IDs of child nodes (empty for leaves)
    pub children: Vec<usize>,
    /// ID of the parent node (None for root)
    pub parent: Option<usize>,
    /// Branch length to parent (0.0 for root)
    pub branch_length: f64,
    /// Number of reads that map to this leaf (0 for internal nodes)
    pub read_count: usize,
}

pub struct Tree {
    pub nodes: Vec<Node>,
    pub root: usize,
}

impl Tree {
    pub fn leaves(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .filter(|n| n.is_leaf)
            .map(|n| n.id)
            .collect()
    }

    /// Iterative post-order traversal (children before parent).
    pub fn post_order(&self) -> Vec<usize> {
        let mut result = Vec::with_capacity(self.nodes.len());
        let mut stack: Vec<(usize, bool)> = vec![(self.root, false)];
        while let Some((id, done)) = stack.pop() {
            if done {
                result.push(id);
            } else {
                stack.push((id, true));
                for &child in self.nodes[id].children.iter().rev() {
                    stack.push((child, false));
                }
            }
        }
        result
    }

    /// Iterative pre-order traversal (parent before children).
    pub fn pre_order(&self) -> Vec<usize> {
        let mut result = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![self.root];
        while let Some(id) = stack.pop() {
            result.push(id);
            for &child in self.nodes[id].children.iter().rev() {
                stack.push(child);
            }
        }
        result
    }

    /// Sum of `read_count` over all leaves in the subtree rooted at `node_id`.
    pub fn subtree_read_count(&self, node_id: usize) -> usize {
        let mut total = 0usize;
        let mut stack = vec![node_id];
        while let Some(id) = stack.pop() {
            if self.nodes[id].is_leaf {
                total += self.nodes[id].read_count;
            }
            for &child in &self.nodes[id].children {
                stack.push(child);
            }
        }
        total
    }
}

pub struct DistMatrix {
    /// Haplotype labels in row/column order (same order as `HaplotypeMatrix::haplotypes`)
    pub labels: Vec<String>,
    /// `data[i][j]` = normalised Hamming distance between haplotype i and j.
    /// Values are in [0.0, 1.0]. Diagonal is 0.0.
    pub data: Vec<Vec<f64>>,
}

/// Build a pairwise Hamming distance matrix from haplotype binary profiles.
///
/// Distance = (number of differing positions) / (total number of variant positions).
/// If there are no variants (shouldn't happen after filtering), all distances are 0.
pub fn hamming_distance_matrix(matrix: &HaplotypeMatrix) -> DistMatrix {
    let n = matrix.haplotypes.len();
    let labels: Vec<String> = matrix.haplotypes.iter().map(|h| h.id.clone()).collect();
    let mut data = vec![vec![0.0f64; n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let mut diff = 0usize;
            let mut compared = 0usize;
            for (a, b) in matrix.haplotypes[i].profile.iter().zip(&matrix.haplotypes[j].profile) {
                if let (Some(a), Some(b)) = (a, b) {
                    compared += 1;
                    if a != b {
                        diff += 1;
                    }
                }
            }
            let d = if compared > 0 { diff as f64 / compared as f64 } else { 0.0 };
            data[i][j] = d;
            data[j][i] = d;
        }
    }

    DistMatrix { labels, data }
}


// ─── Neighbor-Joining ─────────────────────────────────────────────────────────

/// Build a rooted Neighbor-Joining tree from `dist`.
///
/// Rooting strategy:
/// * If a haplotype with the all-zero profile exists (no mutations → reference
///   ancestor), the tree is re-rooted at that haplotype's parent so the
///   outgroup hangs directly from the root.
/// * Otherwise the NJ root (midpoint of the final two-taxon join) is kept.
pub fn neighbor_joining(dist: &DistMatrix, hap_matrix: &HaplotypeMatrix) -> Result<Tree> {
    let n = dist.labels.len();
    if n < 2 {
        bail!("Need at least 2 haplotypes to build a tree, found {n}");
    }

    // ── Initialise leaf nodes ────────────────────────────────────────────────
    let mut nodes: Vec<Node> = dist
        .labels
        .iter()
        .enumerate()
        .map(|(i, label)| Node {
            id: i,
            label: label.clone(),
            is_leaf: true,
            children: vec![],
            parent: None,
            branch_length: 0.0,
            read_count: hap_matrix.haplotypes[i].count,
        })
        .collect();

    // `active[pos]` = node ID; `d[pos_i][pos_j]` = distance between them.
    let mut active: Vec<usize> = (0..n).collect();
    let mut d: Vec<Vec<f64>> = dist.data.clone();

    // ── Main NJ loop ─────────────────────────────────────────────────────────
    while active.len() > 2 {
        let m = active.len();

        // Row sums (r_i = sum_k d[i][k], d[i][i]=0 so the diagonal is free)
        let row_sums: Vec<f64> = (0..m)
            .map(|i| (0..m).map(|j| d[i][j]).sum::<f64>())
            .collect();

        // Find (best_i, best_j) minimising Q[i][j] = (m-2)*d[i][j] - r_i - r_j
        let (mut best_i, mut best_j) = (0, 1);
        let mut best_q = f64::MAX;
        for i in 0..m {
            for j in (i + 1)..m {
                let q = (m as f64 - 2.0) * d[i][j] - row_sums[i] - row_sums[j];
                if q < best_q {
                    best_q = q;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        // Branch lengths from new internal node u to best_i and best_j
        let d_ij    = d[best_i][best_j];
        let delta_i = (d_ij / 2.0
            + (row_sums[best_i] - row_sums[best_j]) / (2.0 * (m as f64 - 2.0)))
            .max(0.0);
        let delta_j = (d_ij - delta_i).max(0.0);

        // Create internal node u
        let u_id = nodes.len();
        nodes[active[best_i]].branch_length = delta_i;
        nodes[active[best_i]].parent        = Some(u_id);
        nodes[active[best_j]].branch_length = delta_j;
        nodes[active[best_j]].parent        = Some(u_id);
        nodes.push(Node {
            id:            u_id,
            label:         format!("INT{u_id}"),
            is_leaf:       false,
            children:      vec![active[best_i], active[best_j]],
            parent:        None,
            branch_length: 0.0,
            read_count:    0,
        });

        // New distances: d[u][k] = (d[best_i][k] + d[best_j][k] - d_ij) / 2
        let remaining: Vec<usize> = (0..m)
            .filter(|&k| k != best_i && k != best_j)
            .collect();
        let new_row: Vec<f64> = remaining
            .iter()
            .map(|&k| (0.5 * (d[best_i][k] + d[best_j][k] - d_ij)).max(0.0))
            .collect();

        let new_m = remaining.len() + 1;
        let mut new_d = vec![vec![0.0f64; new_m]; new_m];
        for (ni, &ri) in remaining.iter().enumerate() {
            for (nj, &rj) in remaining.iter().enumerate() {
                new_d[ni][nj] = d[ri][rj];
            }
            new_d[ni][new_m - 1] = new_row[ni];
            new_d[new_m - 1][ni] = new_row[ni];
        }

        let new_active: Vec<usize> = remaining
            .iter()
            .map(|&ri| active[ri])
            .chain(std::iter::once(u_id))
            .collect();

        active = new_active;
        d = new_d;
    }

    // ── Final join: connect the last two active nodes under a root ───────────
    let root_id  = nodes.len();
    let d_last   = d[0][1];
    let half     = (d_last / 2.0).max(0.0);

    nodes[active[0]].branch_length = half;
    nodes[active[0]].parent        = Some(root_id);
    nodes[active[1]].branch_length = half;
    nodes[active[1]].parent        = Some(root_id);

    nodes.push(Node {
        id:            root_id,
        label:         "ROOT".to_string(),
        is_leaf:       false,
        children:      vec![active[0], active[1]],
        parent:        None,
        branch_length: 0.0,
        read_count:    0,
    });

    let mut tree = Tree { nodes, root: root_id };

    // ── Re-root at all-zero haplotype's parent (outgroup rooting) ────────────
    let zero_profile = vec![Some(0u8); hap_matrix.variants.len()];
    if let Some(og_hap) = hap_matrix
        .haplotypes
        .iter()
        .find(|h| h.profile == zero_profile)
    {
        reroot_at_outgroup_parent(&mut tree, &og_hap.id);
        eprintln!(
            "      Rooted at parent of all-zero haplotype ({}) — outgroup rooting",
            og_hap.id
        );
    }

    Ok(tree)
}

/// Re-root `tree` so that the parent of `outgroup_label` becomes the new root.
///
/// Algorithm: collect the path from the outgroup's parent up to the current
/// root, then reverse all parent–child edges along that path.
fn reroot_at_outgroup_parent(tree: &mut Tree, outgroup_label: &str) {
    // Find the outgroup leaf
    let og_id = match tree.nodes.iter().find(|n| n.label == outgroup_label) {
        Some(n) => n.id,
        None    => return,
    };
    // The new root will be the outgroup's direct parent
    let new_root = match tree.nodes[og_id].parent {
        Some(p) => p,
        None    => return, // outgroup is already the root
    };
    if new_root == tree.root {
        // Already correctly rooted (outgroup's parent IS the old root)
        return;
    }

    // Collect path from new_root up to (and including) old root,
    // capturing branch lengths BEFORE any modification.
    let mut path: Vec<usize>  = vec![new_root];
    let mut bls:  Vec<f64>    = vec![tree.nodes[new_root].branch_length];
    let mut cur = new_root;
    loop {
        match tree.nodes[cur].parent {
            Some(p) => {
                path.push(p);
                bls.push(tree.nodes[p].branch_length);
                if p == tree.root {
                    break;
                }
                cur = p;
            }
            None => break,
        }
    }

    // Reverse edges: for each consecutive pair (child, parent) in path,
    // make child the new parent.
    for i in 0..path.len().saturating_sub(1) {
        let child  = path[i];
        let parent = path[i + 1];

        tree.nodes[parent].children.retain(|&c| c != child);
        tree.nodes[child].children.push(parent);
        tree.nodes[parent].parent        = Some(child);
        tree.nodes[parent].branch_length = bls[i]; // edge length stays the same
    }

    tree.nodes[new_root].parent        = None;
    tree.nodes[new_root].branch_length = 0.0;
    tree.root = new_root;

    // If the old root has become a unary node (only one remaining child),
    // collapse it: attach its single child directly to new_root, summing
    // branch lengths so no information is lost.
    let old_root = *path.last().unwrap();
    if tree.nodes[old_root].children.len() == 1 {
        let sole_child = tree.nodes[old_root].children[0];
        let combined_bl =
            tree.nodes[old_root].branch_length + tree.nodes[sole_child].branch_length;

        // Replace old_root with sole_child in its parent's children
        let old_root_parent = tree.nodes[old_root].parent.unwrap();
        let pos = tree.nodes[old_root_parent]
            .children
            .iter()
            .position(|&c| c == old_root)
            .unwrap();
        tree.nodes[old_root_parent].children[pos] = sole_child;
        tree.nodes[sole_child].parent        = Some(old_root_parent);
        tree.nodes[sole_child].branch_length = combined_bl;
    }
}

/// One unique binary profile and the set of reads that share it.
pub struct Haplotype {
    pub id: String,
    /// State over filtered variants (`profile[v]` = `Some(0)`, `Some(1)`, or `None`)
    pub profile: Vec<Option<u8>>,
    pub reads: Vec<String>,
    pub count: usize,
}

/// Haplotypes together with the (shared) variant list.
pub struct HaplotypeMatrix {
    pub variants: Vec<String>,
    pub haplotypes: Vec<Haplotype>,
}

/// HF (heteroplasmic fraction) keyed by Himito variant ID (e.g. "m.13376T>C").
pub type HfMap = HashMap<String, f64>;

/// Binary matrix after filtering.
/// `data[variant_idx][read_idx]` is `Some(0)` (ref), `Some(1)` (alt), or
/// `None` (uncovered/missing in the source matrix.csv).
pub struct BinaryMatrix {
    pub variants: Vec<String>,
    pub reads: Vec<String>,
    pub data: Vec<Vec<Option<u8>>>,
}

/// Parse a Himito VCF file and return a map from variant ID to HF value.
///
/// Variant IDs are reconstructed as `m.<POS><REF>><ALT>` (the format Himito
/// uses in the matrix CSV), e.g. `m.13376T>C`.
pub fn parse_vcf(vcf_path: &str, min_hf: f64, max_hf: f64) -> Result<HfMap> {
    let mut reader = Reader::from_path(vcf_path)
        .with_context(|| format!("Cannot read VCF: {vcf_path}"))?;
    let mut map = HfMap::new();

    for result in reader.records() {
        #[allow(unused_mut)] // rec.format(b"HF") requires &mut self in rust_htslib
        let mut rec = result?;

        // Only lineage-informative if the call itself is trustworthy.
        if !rec.has_filter("PASS".as_bytes()) {
            continue;
        }

        // Extract position (0-based in htslib → 1-based in VCF/Himito IDs)
        let pos = rec.pos() + 1;

        // rec.alleles() borrows rec immutably; collect into owned Strings so
        // we can release the borrow before calling rec.format() (mut borrow).
        let (ref_allele, alt_allele) = {
            let alleles = rec.alleles();
            (
                std::str::from_utf8(alleles[0]).unwrap_or("").to_owned(),
                std::str::from_utf8(alleles[1]).unwrap_or("").to_owned(),
            )
        };

        let vid = format!("m.{pos}{ref_allele}>{alt_allele}");

        // Read the HF FORMAT field (float, first sample, first value)
        if let Ok(hf_data) = rec.format(b"HF").float() {
            if let Some(&hf_val) = hf_data.get(0).and_then(|s| s.first()) {
                let hf = hf_val as f64;
                if hf >= min_hf && hf < max_hf {
                    map.insert(vid, hf);
                }
            }
        }
    }
    Ok(map)
}

/// Read `matrix_path` (Himito `.matrix.csv`) and apply filters:
///
/// * **HF bounds** — keep rows where `min_hf ≤ HF < max_hf`.
/// * **Prevalence** — keep rows where the variant is present in ≥ `min_presence`
///   reads AND absent from ≥ `min_absence` reads (guarantees a bifurcation).
///
/// Counts are binarised: any count ≥ 1 becomes 1.
pub fn load_and_filter_matrix(
    matrix_path: &str,
    hf_map: &HfMap,   // contains only variants that already passed HF filtering
    min_presence: usize,
    min_absence: usize,
) -> Result<BinaryMatrix> {
    let file = std::fs::File::open(matrix_path)
        .with_context(|| format!("Cannot read matrix CSV: {matrix_path}"))?;
    parse_binary_matrix(file, hf_map, min_presence, min_absence)
}

/// `Read`-generic core of [`load_and_filter_matrix`] so it can be exercised
/// without touching disk. `Some(0)` = explicit ref call, `Some(1)` = alt,
/// `None` = missing/uncovered (empty cell in the CSV).
fn parse_binary_matrix<R: std::io::Read>(
    reader: R,
    hf_map: &HfMap,
    min_presence: usize,
    min_absence: usize,
) -> Result<BinaryMatrix> {
    let mut rdr = csv::Reader::from_reader(reader);
    let headers = rdr.headers()?.clone();
    let reads: Vec<String> = headers.iter().skip(1).map(String::from).collect();

    let mut variants: Vec<String> = Vec::new();
    let mut data: Vec<Vec<Option<u8>>> = Vec::new();

    for result in rdr.records() {
        let rec = result?;
        let vid = rec[0].to_string();

        if !hf_map.is_empty() && !hf_map.contains_key(&vid) {
            continue;
        }

        let row: Vec<Option<u8>> = rec
            .iter()
            .skip(1)
            .map(|v| {
                if v.is_empty() {
                    None
                } else {
                    v.parse::<u32>().ok().map(|n| u8::from(n >= 1))
                }
            })
            .collect();

        let present = row.iter().filter(|c| **c == Some(1)).count();
        let absent = row.iter().filter(|c| **c == Some(0)).count();
        if present < min_presence || absent < min_absence {
            continue;
        }

        variants.push(vid);
        data.push(row);
    }

    Ok(BinaryMatrix { variants, reads, data })
}



/// Collapse reads with identical binary profiles into unique haplotypes.
///
/// Haplotypes are sorted by ascending mutation count (number of 1-bits),
/// then lexicographically by profile, so that the all-zero haplotype (if
/// present) is always `H0000`.
pub fn deduplicate(matrix: &BinaryMatrix, min_reads: usize) -> HaplotypeMatrix {
    let n_reads    = matrix.reads.len();
    let n_variants = matrix.variants.len();

    // profile (per-read, across variants) → list of read indices
    let mut profile_map: HashMap<Vec<Option<u8>>, Vec<usize>> = HashMap::new();

    for read_idx in 0..n_reads {
        let profile: Vec<Option<u8>> = (0..n_variants)
            .map(|v| matrix.data[v][read_idx])
            .collect();
        profile_map.entry(profile).or_default().push(read_idx);
    }

    let mut haplotypes: Vec<Haplotype> = profile_map
        .into_iter()
        .map(|(profile, indices)| {
            let reads = indices.iter().map(|&r| matrix.reads[r].clone()).collect();
            let count = indices.len();
            Haplotype {
                id: String::new(), // assigned after sort
                profile,
                reads,
                count,
            }
        })
        .collect();
    // filter out haplotypes with less than min_reads
    haplotypes = haplotypes.into_iter().filter(|h| h.count >= min_reads).collect();
    // Deterministic ordering: fewest mutations first, then lexicographic profile
    haplotypes.sort_by(|a, b| {
        let ma: usize = a.profile.iter().filter_map(|&x| x).map(usize::from).sum();
        let mb: usize = b.profile.iter().filter_map(|&x| x).map(usize::from).sum();
        ma.cmp(&mb).then_with(|| a.profile.cmp(&b.profile))
    });

    for (i, h) in haplotypes.iter_mut().enumerate() {
        h.id = format!("H{i:04}");
    }

    HaplotypeMatrix {
        variants: matrix.variants.clone(),
        haplotypes,
    }
}


/// Run the four-gamete test on every pair of variants.
///
/// Returns a list of `(variant_a, variant_b)` pairs where all four
/// combinations {00, 01, 10, 11} are observed across haplotypes — a
/// violation of the Infinite Sites Assumption (ISA).
pub fn four_gamete_test(matrix: &HaplotypeMatrix) -> Vec<(String, String)> {
    let n_var = matrix.variants.len();
    let profiles: Vec<Vec<u8>> = matrix
        .haplotypes
        .iter()
        .map(|h| h.profile.iter().map(|c| c.unwrap_or(0)).collect())
        .collect();
    let mut violations = Vec::new();

    for i in 0..n_var {
        for j in (i + 1)..n_var {
            let mut seen = [false; 4]; // index = a*2 + b for states a,b ∈ {0,1}
            for profile in &profiles {
                let idx = (profile[i] * 2 + profile[j]) as usize;
                seen[idx] = true;
                if seen.iter().all(|&s| s) {
                    break;
                }
            }
            if seen.iter().all(|&s| s) {
                violations.push((
                    matrix.variants[i].clone(),
                    matrix.variants[j].clone(),
                ));
            }
        }
    }

    violations
}

pub fn write_haplotype_map(hap_matrix: &HaplotypeMatrix, path: &str) -> Result<()> {
    let mut w = BufWriter::new(
        File::create(path).with_context(|| format!("Cannot create {path}"))?,
    );
    writeln!(w, "haplotype_id\tn_mutations\tn_reads\tread_name")?;
    for hap in &hap_matrix.haplotypes {
        let n_mut: usize = hap.profile.iter().filter_map(|&b| b).map(usize::from).sum();
        let readlist = hap.reads.join(",");
        assert_eq!(hap.count, hap.reads.len());
        writeln!(w, "{}\t{n_mut}\t{}\t{readlist}", hap.id, hap.count)?;
        
    }
    Ok(())
}


pub fn start(
    matrix_file: &str,
    vcf_file: Option<&str>,
    min_hf: f64,
    max_hf: f64,
    min_presence: usize,
    min_absence: usize,
    min_reads: usize,
    fp_rate: f64,
    fn_rate: f64,
    mcmc_iterations: usize,
    mcmc_chains: usize,
    mcmc_seed: u64,
    output_prefix: &str,
) -> Result<()> {
    println!("Starting lineage analysis...");

    // ── Step 1: load VCF HF values, then filter the binary matrix ─────────────
    println!("[1/6] Parsing VCF: {}", vcf_file.as_ref().unwrap_or(&""));
    let hf_map = parse_vcf(vcf_file.as_ref().unwrap_or(&""), min_hf, max_hf)?;

    eprintln!("[1/6] Loading and filtering matrix: {}", matrix_file);
    let binary = load_and_filter_matrix(
        &matrix_file,
        &hf_map,
        min_presence,
        min_absence,
    )?;

    if binary.variants.is_empty() {
        anyhow::bail!(
            "No informative variants remain after filtering. \
             Try adjusting --min-hf / --max-hf thresholds."
        );
    }

    // ── Step 2: deduplicate reads into haplotypes ───────────────────────────────
    eprintln!("[2/6] Deduplicating reads into haplotypes...");
    let hap_matrix = deduplicate(&binary, min_reads);
    // write how many haplotypes and how many heteroplasmic variants on each haplotype
    eprintln!("[2/6] Found {} haplotypes across {} variants.", hap_matrix.haplotypes.len(), hap_matrix.variants.len());
    if hap_matrix.haplotypes.is_empty() {
        anyhow::bail!("No haplotypes remain after --min-reads filtering.");
    }
    // Global haplotype map (all variants)
    let hmap_path = format!("{}.raw_haplotype_map.tsv", output_prefix);
    write_haplotype_map(&hap_matrix, &hmap_path)?;

    // ── Step 3: SCITE mutation-tree search (NJ-informed MCMC) ─────────────────
    eprintln!("[3/3] Running SCITE mutation-tree search...");
    crate::scite::run_scite_pipeline(
        &binary,
        &hap_matrix,
        fp_rate,
        fn_rate,
        mcmc_iterations,
        mcmc_chains,
        mcmc_seed,
        output_prefix,
    )?;

    Ok(())

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn hf_map(names: &[&str]) -> HfMap {
        names.iter().map(|n| (n.to_string(), 0.5)).collect()
    }

    #[test]
    fn load_and_filter_matrix_preserves_missing_values() {
        let csv = "variant,r1,r2,r3\nm.100A>G,1,0,\nm.200C>T,0,1,1\n";
        let hf = hf_map(&["m.100A>G", "m.200C>T"]);
        let matrix = parse_binary_matrix(Cursor::new(csv), &hf, 1, 1).unwrap();

        assert_eq!(matrix.variants, vec!["m.100A>G", "m.200C>T"]);
        assert_eq!(matrix.reads, vec!["r1", "r2", "r3"]);
        assert_eq!(
            matrix.data,
            vec![
                vec![Some(1), Some(0), None],
                vec![Some(0), Some(1), Some(1)],
            ]
        );
    }

    #[test]
    fn load_and_filter_matrix_prevalence_filter_does_not_count_missing_as_absent() {
        // present=1 (r1), absent=0 (no explicit ref call) -> fails min_absence=1
        let csv = "variant,r1,r2\nm.100A>G,1,\n";
        let hf = hf_map(&["m.100A>G"]);
        let matrix = parse_binary_matrix(Cursor::new(csv), &hf, 1, 1).unwrap();
        assert!(matrix.variants.is_empty());
    }

    #[test]
    fn deduplicate_treats_different_missing_patterns_as_different_haplotypes() {
        let matrix = BinaryMatrix {
            variants: vec!["m.100A>G".to_string()],
            reads: vec!["r1".to_string(), "r2".to_string()],
            // Both reads are "ref-or-missing" if you ignore the distinction,
            // but r1 is an explicit ref call and r2 is genuinely uncovered.
            data: vec![vec![Some(0), None]],
        };
        let hap_matrix = deduplicate(&matrix, 1);
        assert_eq!(hap_matrix.haplotypes.len(), 2);
    }

    #[test]
    fn hamming_distance_matrix_normalizes_by_jointly_covered_sites_only() {
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["A".to_string(), "B".to_string(), "C".to_string()],
            haplotypes: vec![
                Haplotype {
                    id: "H0".to_string(),
                    profile: vec![Some(1), Some(0), None],
                    reads: vec![],
                    count: 1,
                },
                Haplotype {
                    id: "H1".to_string(),
                    profile: vec![Some(1), Some(1), Some(1)],
                    reads: vec![],
                    count: 1,
                },
            ],
        };
        let dist = hamming_distance_matrix(&hap_matrix);
        // Jointly covered sites are A and B only (H0 is missing at C).
        // They differ at B only -> distance = 1/2.
        assert!((dist.data[0][1] - 0.5).abs() < 1e-12);
        assert!((dist.data[1][0] - 0.5).abs() < 1e-12);
    }
}