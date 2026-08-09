use std::collections::HashMap;
use anyhow::{bail,Context, Result};
use rust_htslib::bcf::{Read, Reader}; // Read trait must be in scope for .records()
use std::io::{BufWriter, Write};
use std::fs::File;
use log::info;

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

    // ── Re-root at ancestral haplotype's parent (outgroup rooting) ───────────
    // The ancestral (reference) haplotype carries no alt calls. Match on the
    // absence of any `Some(1)` rather than an exact all-`Some(0)` profile, so an
    // otherwise-reference haplotype with uncovered sites (`None`) still roots the
    // tree. Haplotypes are sorted fewest-mutations-first, so the first match is
    // the best-covered ancestral candidate.
    if let Some(og_hap) = hap_matrix
        .haplotypes
        .iter()
        .find(|h| !h.profile.iter().any(|&c| c == Some(1)))
    {
        reroot_at_outgroup_parent(&mut tree, &og_hap.id);
        eprintln!(
            "      Rooted at parent of ancestral haplotype ({}) — outgroup rooting",
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

        // rec.alleles() borrows rec immutably; collect into owned Strings so we
        // can release the borrow before calling rec.format() (mut borrow). Keep
        // every ALT so multiallelic sites contribute one variant id per ALT
        // (skipping alleles[1..] gracefully handles ALT-less records).
        let (ref_allele, alt_alleles): (String, Vec<String>) = {
            let alleles = rec.alleles();
            let ref_allele = std::str::from_utf8(alleles[0]).unwrap_or("").to_owned();
            let alts = alleles
                .iter()
                .skip(1)
                .map(|a| std::str::from_utf8(a).unwrap_or("").to_owned())
                .collect();
            (ref_allele, alts)
        };

        // Read the HF FORMAT field (float, first sample). Prefer Himito's HF;
        // fall back to AF when HF is absent. HF/AF are Number=A, so the i-th
        // value lines up with the i-th ALT; fall back to the first value for
        // single-valued encodings.
        let hf_data = rec
            .format(b"HF")
            .float()
            .or_else(|_| rec.format(b"AF").float())
            .ok();
        
        for (alt_idx, alt_allele) in alt_alleles.iter().enumerate() {
            let hf_val = hf_data
                .as_ref()
                .and_then(|d| d.get(0).and_then(|s| s.get(alt_idx).or_else(|| s.first())));
            let Some(&hf_val) = hf_val else { continue };
            let hf = hf_val as f64;
            // println!("hf_data: {:?} {} {}", hf, min_hf, max_hf);
            if hf >= min_hf && hf < max_hf {
                let vid = format!("m.{pos}{ref_allele}>{alt_allele}");
                // Keep the first occurrence of a variant id (deterministic across
                // duplicate records) rather than silently overwriting.
                map.entry(vid).or_insert(hf);
            }
        }
    }
    Ok(map)
}

/// Read `matrix_path` (Himito `.matrix.csv`) and apply filters:
///
/// * **Prevalence** — keep rows where the variant is present in ≥ `min_presence`
///   reads AND absent from ≥ `min_absence` reads (guarantees a bifurcation).
/// * **HF bounds** — keep rows with `min_hf ≤ HF < max_hf`. When a VCF was
///   provided (`hf_map` non-empty) the HF comes from the VCF and this function
///   only checks membership. When no VCF was provided (`hf_map` empty) the HF is
///   recomputed from the matrix as `present / (present + absent)` — alt calls
///   over *covered* reads — matching Himito's `HF = allele_count / read_depth`,
///   so the with-/without-VCF variant sets agree.
///
/// Counts are binarised: any count ≥ 1 becomes 1.
pub fn load_and_filter_matrix(
    matrix_path: &str,
    hf_map: &HfMap,   // contains only variants that already passed HF filtering
    min_presence: usize,
    min_absence: usize,
    min_hf: f64,
    max_hf: f64,
) -> Result<BinaryMatrix> {
    let file = std::fs::File::open(matrix_path)
        .with_context(|| format!("Cannot read matrix CSV: {matrix_path}"))?;
    parse_binary_matrix(file, hf_map, min_presence, min_absence, min_hf, max_hf)
}

/// `Read`-generic core of [`load_and_filter_matrix`] so it can be exercised
/// without touching disk. `Some(0)` = explicit ref call, `Some(1)` = alt,
/// `None` = missing/uncovered (empty cell in the CSV).
fn parse_binary_matrix<R: std::io::Read>(
    reader: R,
    hf_map: &HfMap,
    min_presence: usize,
    min_absence: usize,
    min_hf: f64,
    max_hf: f64,
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

        // No-VCF fallback HF: alt calls over *covered* reads (present + absent),
        // NOT row.len(). This mirrors Himito's VCF `HF = allele_count / read_depth`
        // so the homoplasmy filter keeps the same variants with and without a VCF.
        // Uncovered (`None`) cells must be excluded from the denominator.
        let covered = present + absent;
        let freq = if covered > 0 { present as f64 / covered as f64 } else { 0.0 };

        if hf_map.is_empty() && (freq < min_hf || freq >= max_hf) {
            continue;
        }

        // println!("vid: {}, freq: {}", vid, freq);

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
    let profiles: Vec<&Vec<Option<u8>>> =
        matrix.haplotypes.iter().map(|h| &h.profile).collect();
    let mut violations = Vec::new();

    for i in 0..n_var {
        for j in (i + 1)..n_var {
            let mut seen = [false; 4]; // index = a*2 + b for states a,b ∈ {0,1}
            for profile in &profiles {
                // A haplotype only informs the pair if it covers BOTH sites;
                // missing calls (`None`) are excluded rather than read as ref.
                if let (Some(a), Some(b)) = (profile[i], profile[j]) {
                    let idx = (a * 2 + b) as usize;
                    seen[idx] = true;
                    if seen.iter().all(|&s| s) {
                        break;
                    }
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
    // Always represent the ancestral haplotype (the tree root: no heteroplasmic
    // variants). If the deduplicated haplotypes already contain an all-reference
    // profile it is kept as-is; otherwise a synthetic zero-read root row is
    // appended so the ancestral state is consistent with the lineage tree.
    let has_ancestral = hap_matrix
        .haplotypes
        .iter()
        .any(|h| !h.profile.iter().any(|&b| b == Some(1)));
    if !has_ancestral {
        writeln!(w, "Hroot\t0\t0\t")?;
    }
    Ok(())
}


/// Per-data-type SCITE error rates, with explicit values always winning.
///
/// `None` means "not supplied on the command line", so the data-type default
/// applies; `Some(v)` is used verbatim. Mirrors `call::resolve_thresholds`.
///
/// The defaults come from the fp/fn sweeps in `lineage_eval`
/// (`sweep_fpfn.sh`, scored by `score_lineage.py`):
///
/// * `pacbio` (HiFi) - every cell of the 5x6 grid (fp 0.0005-0.02,
///   fn 0.01-0.2) reached `ad_f1 == 1.0`, so HiFi is insensitive here and any
///   in-grid value works. `fp = 0.005` / `fn = 0.05` sits comfortably inside.
/// * `ont-denoised` - denoising removes most per-read noise, so the rates drop
///   accordingly.
/// * `ont-r9` / `ont-r10` - the raw-ONT sweep is driven almost entirely by
///   `fn`: `ad_f1` holds at its 0.375 plateau for `fn <= 0.05` and collapses to
///   0.118 at `fn >= 0.1`, while `fp` is near-irrelevant across 0-0.02. The
///   default therefore sits at the top of the safe `fn` range.
pub fn resolve_error_rates(
    data_type: &str,
    fp_rate: Option<f64>,
    fn_rate: Option<f64>,
) -> (f64, f64) {
    let (default_fp, default_fn) = match data_type {
        "pacbio" => (0.005, 0.05),
        "ont-denoised" => (0.0001, 0.001),
        // ont-r9 / ont-r10 and anything else: raw long-read noise levels.
        _ => (0.001, 0.05),
    };

    (
        fp_rate.unwrap_or(default_fp),
        fn_rate.unwrap_or(default_fn),
    )
}

pub fn start(
    matrix_file: &str,
    vcf_file: Option<&str>,
    min_hf: f64,
    max_hf: f64,
    min_presence: usize,
    min_absence: usize,
    min_reads: usize,
    data_type: &str,
    fp_rate: f64,
    fn_rate: f64,
    mcmc_iterations: usize,
    mcmc_chains: usize,
    mcmc_seed: u64,
    output_prefix: &str,
) -> Result<()> {
    env_logger::init();
    if !matches!(data_type, "pacbio" | "ont-r9" | "ont-r10" | "ont-denoised") {
        anyhow::bail!(
            "data type must be pacbio, ont-r9, ont-r10 or ont-denoised (got {data_type})"
        );
    }
    info!("Starting lineage analysis...");
    info!("Data type {data_type}: SCITE fp={fp_rate}, fn={fn_rate}");
    //ont-denoise fp=0.0001 fn=0.001

    // ── Step 1: load VCF HF values, then filter the binary matrix ─────────────
    // The VCF is optional: without it, no HF filter is applied and every variant
    // in the matrix is retained (an empty HfMap disables HF filtering downstream).
    let hf_map = match vcf_file {
        Some(path) => {
            info!("[1/6] Parsing VCF: {}", path);
            parse_vcf(path, min_hf, max_hf)?
        }
        None => {
            info!("[1/6] No VCF provided; skipping HF filter (all matrix variants retained).");
            HfMap::new()
        }
    };
    // println!("hf_map: {:?}", hf_map);

    info!("[1/6] Loading and filtering matrix: {}", matrix_file);
    let binary = load_and_filter_matrix(
        &matrix_file,
        &hf_map,
        min_presence,
        min_absence,
        min_hf,
        max_hf,
    )?;
    if binary.variants.is_empty() {
        anyhow::bail!(
            "No informative variants remain after filtering. \
             Try adjusting --min-hf / --max-hf thresholds."
        );
    }

    // ── Step 2: deduplicate reads into haplotypes ───────────────────────────────
    info!("[2/6] Deduplicating reads into haplotypes...");
    let hap_matrix = deduplicate(&binary, min_reads);
    // write how many haplotypes and how many heteroplasmic variants on each haplotype
    info!("[2/6] Found {} haplotypes across {} variants.", hap_matrix.haplotypes.len(), hap_matrix.variants.len());
    if hap_matrix.haplotypes.is_empty() {
        anyhow::bail!("No haplotypes remain after --min-reads filtering.");
    }
    // Diagnostic: variant pairs violating the infinite-sites assumption (all four
    // gametes observed) hint at recurrent/back mutation or recombination.
    let violations = four_gamete_test(&hap_matrix);
    if !violations.is_empty() {
        info!(
            "[2/6] Four-gamete test: {} variant pair(s) violate the infinite-sites assumption.",
            violations.len()
        );
    }
    // Global haplotype map (all variants)
    let hmap_path = format!("{}.raw_haplotype_map.tsv", output_prefix);
    write_haplotype_map(&hap_matrix, &hmap_path)?;

    // ── Step 3: SCITE mutation-tree search (NJ-informed MCMC) ─────────────────
    info!("[3/6] Running SCITE mutation-tree search...");
    crate::scite::run_scite_pipeline(
        &binary,
        &hap_matrix,
        fp_rate,
        fn_rate,
        mcmc_iterations,
        mcmc_chains,
        mcmc_seed,
        min_reads,
        output_prefix,
    )?;

    // ── Step 4: write final haplotype map (filtered by SCITE) ─────────────────
    info!("[4/6] Deduplicating cleaned reads into haplotypes...");
    // The cleaned matrix already contains exactly the variants the tree was built
    // on; re-applying prevalence/HF filters here would re-drop variants after
    // imputation (which removed the `None` cells and can push a variant to
    // absent==0 or freq==1.0), making the no-VCF output diverge from the VCF one.
    // Pass fully-permissive thresholds so this reload only deduplicates.
    let cleaned_binary = load_and_filter_matrix(
        &format!("{}.cleaned_matrix.csv", output_prefix),
        &hf_map,
        0,
        0,
        0.0,
        f64::INFINITY,
    )?;
    let cleaned_hap_matrix = deduplicate(&cleaned_binary, min_reads);
    // write how many haplotypes and how many heteroplasmic variants on each haplotype
    info!("[4/6] Found {} haplotypes across {} variants.", cleaned_hap_matrix.haplotypes.len(), cleaned_hap_matrix.variants.len());
    if cleaned_hap_matrix.haplotypes.is_empty() {
        anyhow::bail!("No haplotypes remain after --min-reads filtering.");
    }
    // Global haplotype map (all variants)
    let hmap_path = format!("{}.cleaned_haplotype_map.tsv", output_prefix);
    write_haplotype_map(&cleaned_hap_matrix, &hmap_path)?;

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
    fn resolve_error_rates_presets_differ_per_data_type() {
        assert_eq!(resolve_error_rates("pacbio", None, None), (0.005, 0.05));
        assert_eq!(resolve_error_rates("ont-denoised", None, None), (0.0001, 0.001));
        assert_eq!(resolve_error_rates("ont-r10", None, None), (0.001, 0.05));
        assert_eq!(resolve_error_rates("ont-r9", None, None), (0.001, 0.05));
        // pacbio and ont-denoised must not collapse onto the raw-ONT fallback.
        assert_ne!(
            resolve_error_rates("pacbio", None, None),
            resolve_error_rates("ont-r10", None, None)
        );
    }

    #[test]
    fn resolve_error_rates_explicit_values_override_the_preset() {
        // Either side can be overridden independently; the other keeps its preset.
        assert_eq!(resolve_error_rates("pacbio", Some(0.02), None), (0.02, 0.05));
        assert_eq!(resolve_error_rates("pacbio", None, Some(0.2)), (0.005, 0.2));
        assert_eq!(
            resolve_error_rates("ont-denoised", Some(0.01), Some(0.1)),
            (0.01, 0.1)
        );
        // 0.0 is a real user choice, not "unset".
        assert_eq!(resolve_error_rates("ont-r10", Some(0.0), Some(0.0)), (0.0, 0.0));
    }

    #[test]
    fn load_and_filter_matrix_preserves_missing_values() {
        let csv = "variant,r1,r2,r3\nm.100A>G,1,0,\nm.200C>T,0,1,1\n";
        let hf = hf_map(&["m.100A>G", "m.200C>T"]);
        let min_presence = 1;
        let min_absence = 1;
        let min_hf = 0.0;
        let max_hf = 1.0;
        let matrix = parse_binary_matrix(Cursor::new(csv), &hf, min_presence, min_absence, min_hf, max_hf).unwrap();

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
        let matrix = parse_binary_matrix(Cursor::new(csv), &hf, 1, 1, 0.0, 1.0).unwrap();
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
    fn write_haplotype_map_appends_synthetic_root_when_no_ancestral_haplotype() {
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["m.100A>G".to_string()],
            haplotypes: vec![Haplotype {
                id: "H0000".to_string(),
                profile: vec![Some(1)],
                reads: vec!["r1".to_string(), "r2".to_string()],
                count: 2,
            }],
        };
        let path = std::env::temp_dir().join("himito_test_hapmap_synthetic_root.tsv");
        write_haplotype_map(&hap_matrix, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(
            content,
            "haplotype_id\tn_mutations\tn_reads\tread_name\n\
             H0000\t1\t2\tr1,r2\n\
             Hroot\t0\t0\t\n"
        );
    }

    #[test]
    fn write_haplotype_map_keeps_existing_ancestral_haplotype_without_duplicating() {
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["m.100A>G".to_string()],
            haplotypes: vec![Haplotype {
                id: "H0000".to_string(),
                profile: vec![Some(0)],
                reads: vec!["r1".to_string()],
                count: 1,
            }],
        };
        let path = std::env::temp_dir().join("himito_test_hapmap_existing_root.tsv");
        write_haplotype_map(&hap_matrix, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(
            content,
            "haplotype_id\tn_mutations\tn_reads\tread_name\n\
             H0000\t0\t1\tr1\n"
        );
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

    fn hap2(id: &str, a: Option<u8>, b: Option<u8>) -> Haplotype {
        Haplotype { id: id.into(), profile: vec![a, b], reads: vec![], count: 1 }
    }

    #[test]
    fn four_gamete_test_flags_all_four_gametes() {
        let matrix = HaplotypeMatrix {
            variants: vec!["A".into(), "B".into()],
            haplotypes: vec![
                hap2("H0", Some(0), Some(0)),
                hap2("H1", Some(0), Some(1)),
                hap2("H2", Some(1), Some(0)),
                hap2("H3", Some(1), Some(1)),
            ],
        };
        assert_eq!(four_gamete_test(&matrix), vec![("A".to_string(), "B".to_string())]);
    }

    #[test]
    fn four_gamete_test_missing_calls_do_not_fabricate_the_ref_gamete() {
        // Only 01, 10, 11 are truly observed. The "00" gamete would only appear
        // if a None were read as ref (the old behavior). With None excluded there
        // is no violation.
        let matrix = HaplotypeMatrix {
            variants: vec!["A".into(), "B".into()],
            haplotypes: vec![
                hap2("H0", None, None),
                hap2("H1", Some(0), Some(1)),
                hap2("H2", Some(1), Some(0)),
                hap2("H3", Some(1), Some(1)),
            ],
        };
        assert!(four_gamete_test(&matrix).is_empty());
    }

    #[test]
    fn reroot_at_outgroup_parent_collapses_old_root_and_preserves_path_lengths() {
        // ROOT(4) -> {INT3(3) -> {H0000(0), H1(1)}, H2(2)}
        let mut tree = Tree {
            nodes: vec![
                Node { id: 0, label: "H0000".into(), is_leaf: true, children: vec![], parent: Some(3), branch_length: 1.0, read_count: 1 },
                Node { id: 1, label: "H1".into(), is_leaf: true, children: vec![], parent: Some(3), branch_length: 2.0, read_count: 1 },
                Node { id: 2, label: "H2".into(), is_leaf: true, children: vec![], parent: Some(4), branch_length: 5.0, read_count: 1 },
                Node { id: 3, label: "INT3".into(), is_leaf: false, children: vec![0, 1], parent: Some(4), branch_length: 3.0, read_count: 0 },
                Node { id: 4, label: "ROOT".into(), is_leaf: false, children: vec![3, 2], parent: None, branch_length: 0.0, read_count: 0 },
            ],
            root: 4,
        };
        reroot_at_outgroup_parent(&mut tree, "H0000");

        assert_eq!(tree.root, 3);
        assert_eq!(tree.nodes[3].parent, None);
        // Old root (4) collapses; H2 hangs directly under the new root with its
        // branch length summed through the old root (5.0 + 3.0).
        assert_eq!(tree.nodes[2].parent, Some(3));
        assert!((tree.nodes[2].branch_length - 8.0).abs() < 1e-12);
        // Outgroup and sibling keep their original lengths to INT3.
        assert!((tree.nodes[0].branch_length - 1.0).abs() < 1e-12);
        assert!((tree.nodes[1].branch_length - 2.0).abs() < 1e-12);
        // New root's children are exactly {H0000, H1, H2}.
        let mut kids = tree.nodes[3].children.clone();
        kids.sort();
        assert_eq!(kids, vec![0, 1, 2]);
    }

    #[test]
    fn parse_vcf_splits_multiallelic_sites_and_keeps_first_duplicate() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##FORMAT=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmic fraction\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chrM\t100\t.\tA\tG,T\t.\tPASS\t.\tHF\t0.3,0.6
chrM\t200\t.\tC\tT\t.\tPASS\t.\tHF\t0.9
chrM\t100\t.\tA\tG\t.\tPASS\t.\tHF\t0.99
";
        let path = std::env::temp_dir().join("himito_test_parse_vcf_multiallelic.vcf");
        std::fs::write(&path, vcf).unwrap();
        let map = parse_vcf(path.to_str().unwrap(), 0.0, 1.0).unwrap();
        std::fs::remove_file(&path).ok();

        // Both ALTs of the multiallelic site become distinct variant ids.
        assert_eq!(map.get("m.100A>G").copied(), Some(0.3_f32 as f64));
        assert_eq!(map.get("m.100A>T").copied(), Some(0.6_f32 as f64));
        assert_eq!(map.get("m.200C>T").copied(), Some(0.9_f32 as f64));
        // The later duplicate record (HF 0.99) does not overwrite the first.
        assert_eq!(map.len(), 3);
    }
}

/// Parity checks: the with-VCF and without-VCF variant filtering must select the
/// same variants, so a lineage run produces identical output either way (assuming
/// a VCF whose HF matches the matrix and which applies no extra filters).
#[cfg(test)]
mod vcf_parity_tests {
    use super::*;
    use std::io::Cursor;

    /// Build a `variant,r1..r10` CSV. `cells` are raw tokens ("" = uncovered).
    fn build_csv(rows: &[(&str, [&str; 10])]) -> String {
        let mut csv = String::from("variant,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10\n");
        for (vid, cells) in rows {
            csv.push_str(vid);
            for c in cells {
                csv.push(',');
                csv.push_str(c);
            }
            csv.push('\n');
        }
        csv
    }

    #[test]
    fn no_vcf_matrix_frequency_filter_matches_vcf_hf_filter() {
        // present/absent/uncovered per variant, thresholds min_hf=0.2 max_hf=0.85:
        //   m.100A>G  present=3 absent=3 uncovered=4 -> HF 3/6  = 0.50  KEEP
        //   m.200C>T  present=1 absent=3 uncovered=6 -> HF 1/4  = 0.25  KEEP
        //             (row.len() denom would give 1/10 = 0.10 -> wrongly DROPPED,
        //              so this row exercises the covered-read denominator)
        //   m.300G>A  present=9 absent=1 uncovered=0 -> HF 9/10 = 0.90  DROP (homoplasmic)
        //   m.400T>C  present=1 absent=9 uncovered=0 -> HF 1/10 = 0.10  DROP (too low)
        let csv = build_csv(&[
            ("m.100A>G", ["1", "1", "1", "0", "0", "0", "", "", "", ""]),
            ("m.200C>T", ["1", "0", "0", "0", "", "", "", "", "", ""]),
            ("m.300G>A", ["1", "1", "1", "1", "1", "1", "1", "1", "1", "0"]),
            ("m.400T>C", ["1", "0", "0", "0", "0", "0", "0", "0", "0", "0"]),
        ]);

        let (min_presence, min_absence, min_hf, max_hf) = (1usize, 1usize, 0.2f64, 0.85f64);

        // "With VCF": hf_map holds exactly the variants whose HF is in range, as
        // parse_vcf would produce (HF computed as alt/coverage).
        let mut with_vcf_hf = HfMap::new();
        with_vcf_hf.insert("m.100A>G".to_string(), 0.50);
        with_vcf_hf.insert("m.200C>T".to_string(), 0.25);

        let with_vcf = parse_binary_matrix(
            Cursor::new(csv.clone()),
            &with_vcf_hf,
            min_presence,
            min_absence,
            min_hf,
            max_hf,
        )
        .unwrap();

        // "Without VCF": empty hf_map -> frequency recomputed from the matrix.
        let without_vcf = parse_binary_matrix(
            Cursor::new(csv),
            &HfMap::new(),
            min_presence,
            min_absence,
            min_hf,
            max_hf,
        )
        .unwrap();

        // Both keep exactly the two heteroplasmic variants, same order, same data.
        assert_eq!(without_vcf.variants, vec!["m.100A>G", "m.200C>T"]);
        assert_eq!(with_vcf.variants, without_vcf.variants);
        assert_eq!(with_vcf.data, without_vcf.data);
        assert_eq!(with_vcf.reads, without_vcf.reads);
    }

    #[test]
    fn no_vcf_denominator_excludes_uncovered_reads() {
        // Heteroplasmic among covered reads but "rare" if uncovered cells are
        // counted: present=2 absent=2 uncovered=6 -> covered HF 2/4 = 0.50,
        // row.len() HF 2/10 = 0.20. A threshold between the two isolates the bug.
        let csv = build_csv(&[("m.500A>C", ["1", "1", "0", "0", "", "", "", "", "", ""])]);
        let kept = parse_binary_matrix(
            Cursor::new(csv),
            &HfMap::new(),
            1,
            1,
            0.30f64,
            0.85f64,
        )
        .unwrap();
        // Covered-read denominator (0.50) retains it; the old row.len() denominator
        // (0.20 < 0.30) would have dropped it.
        assert_eq!(kept.variants, vec!["m.500A>C"]);
    }
}