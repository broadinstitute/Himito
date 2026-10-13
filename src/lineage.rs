use std::collections::HashMap;
use anyhow::{bail,Context, Result};
use rust_htslib::bcf::{Read, Reader}; // Read trait must be in scope for .records()
use std::io::{BufWriter, Write};
use std::fs::File;
use log::{info, warn};

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
    pub data: Vec<Vec<f64>>,
}

/// Build a pairwise Hamming distance matrix from haplotype binary profiles.
/// The matrix never contains NaN.
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
            let d = if compared > 0 {
                diff as f64 / compared as f64
            } else {
                f64::NAN
            };
            data[i][j] = d;
            data[j][i] = d;
        }
    }

    let mut sum = 0.0;
    let mut count = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            if data[i][j].is_finite() {
                sum += data[i][j];
                count += 1;
            }
        }
    }
    let mean = if count > 0 { sum / count as f64 } else { 0.0 };
    for i in 0..n {
        for j in (i + 1)..n {
            if !data[i][j].is_finite() {
                data[i][j] = mean;
                data[j][i] = mean;
            }
        }
    }

    DistMatrix { labels, data }
}


// ─── Neighbor-Joining ─────────────────────────────────────────────────────────

/// Build a rooted Neighbor-Joining tree from `dist`.
pub fn neighbor_joining(hap_matrix: &HaplotypeMatrix) -> Result<Tree> {
    // The distance matrix is a pure function of `hap_matrix`.
    let dist = hamming_distance_matrix(hap_matrix);
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
    if let Some(og_hap) = pick_outgroup(hap_matrix) {
        reroot_at_outgroup_parent(&mut tree, &og_hap.id);
        eprintln!(
            "      Rooted at parent of ancestral haplotype ({}) — outgroup rooting",
            og_hap.id
        );
    }

    Ok(tree)
}

/// Choose the ancestral (no-alt) haplotype to use as the NJ outgroup.
fn pick_outgroup(hm: &HaplotypeMatrix) -> Option<&Haplotype> {
    hm.haplotypes
        .iter()
        .filter(|h| !h.profile.iter().any(|&c| c == Some(1)))
        .filter(|h| h.profile.iter().any(|c| c.is_some()))
        .max_by(|a, b| {
            let cov = |h: &Haplotype| h.profile.iter().filter(|c| c.is_some()).count();
            cov(a)
                .cmp(&cov(b))
                .then(a.count.cmp(&b.count))
                .then(std::cmp::Ordering::Less)
        })
}

/// Re-root `tree` so that the parent of `outgroup_label` becomes the new root.
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

    // Collect path from new_root up to (and including) old root, capturing branch lengths BEFORE any modification.
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

    // Reverse edges: for each consecutive pair (child, parent) in path, make child the new parent.
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

    // If the old root has become a unary node, collapse it: attach its single child directly to new_root.
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
pub struct BinaryMatrix {
    pub variants: Vec<String>,
    pub reads: Vec<String>,
    pub data: Vec<Vec<Option<u8>>>,
}

/// Every row of a Himito `.matrix.csv`, binarised but **unfiltered**.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RawMatrix {
    pub variants: Vec<String>,
    pub reads: Vec<String>,
    pub data: Vec<Vec<Option<u8>>>,
}

/// Parse a Himito VCF file and return a map from variant ID to HF value.
pub fn parse_vcf(vcf_path: &str, min_hf: f64, max_hf: f64) -> Result<HfMap> {
    let mut reader = Reader::from_path(vcf_path)
        .with_context(|| format!("Cannot read VCF: {vcf_path}"))?;
    // `call` writes this only when indels were left-aligned.
    let left_aligned = reader.header().header_records().iter().any(|rec| {
        matches!(
            rec,
            rust_htslib::bcf::header::HeaderRecord::Generic { key, .. }
                if key == "HimitoNormalization"
        )
    });
    if !left_aligned {
        warn!(
            "VCF was not left-aligned; homopolymer indels may be split across several variant names; re-run `Himito call` without `--no-left-align`"
        );
    }
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

        // rec.alleles() borrows rec immutably; collect into owned Strings.
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

        // Read the HF FORMAT field (float, first sample).
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
                // Keep the first occurrence of a variant id (deterministic across duplicate records) rather than silently overwriting.
                map.entry(vid).or_insert(hf);
            }
        }
    }
    Ok(map)
}

/// How a variant's heteroplasmic frequency is judged during matrix filtering.
pub enum HfFilter {
    /// A VCF was supplied.
    FromVcf { map: HfMap, min: f64, max: f64 },
    /// No VCF was supplied.
    FromMatrix { min: f64, max: f64 },
}

impl HfFilter {
    /// `tract_rows` holds the names of the collapsed tract-allele rows.
    fn keeps(&self, vid: &str, freq: f64, tract_rows: &std::collections::HashSet<String>) -> bool {
        match self {
            HfFilter::FromVcf { map, min, max } => {
                if tract_rows.contains(vid) {
                    freq >= *min && freq < *max
                } else {
                    map.contains_key(vid)
                }
            }
            HfFilter::FromMatrix { min, max } => freq >= *min && freq < *max,
        }
    }

    /// Whether the decision for `vid` can be made from its name alone, before the row is counted.
    fn rejects_by_name(&self, vid: &str, tract_rows: &std::collections::HashSet<String>) -> bool {
        matches!(self, HfFilter::FromVcf { .. })
            && !tract_rows.contains(vid)
            && !self.keeps(vid, 0.0, tract_rows)
    }
}

/// Read `matrix_path` (Himito `.matrix.csv`) and apply filters:
#[allow(dead_code)]
pub fn load_and_filter_matrix(
    matrix_path: &str,
    hf: &HfFilter,
    min_presence: usize,
    min_absence: usize,
) -> Result<BinaryMatrix> {
    let file = std::fs::File::open(matrix_path)
        .with_context(|| format!("Cannot read matrix CSV: {matrix_path}"))?;
    let raw = read_raw_matrix(file)?;
    Ok(filter_matrix(&raw, hf, min_presence, min_absence, &Default::default()))
}

/// Parse every row of a Himito `.matrix.csv` into a [`RawMatrix`], applying no filter at all.
pub fn read_raw_matrix<R: std::io::Read>(reader: R) -> Result<RawMatrix> {
    let mut rdr = csv::Reader::from_reader(reader);
    let headers = rdr.headers()?.clone();
    let reads: Vec<String> = headers.iter().skip(1).map(String::from).collect();

    let mut variants: Vec<String> = Vec::new();
    let mut data: Vec<Vec<Option<u8>>> = Vec::new();

    for result in rdr.records() {
        let rec = result?;
        variants.push(rec[0].to_string());
        data.push(
            rec.iter()
                .skip(1)
                .map(|v| {
                    if v.is_empty() {
                        None
                    } else {
                        v.parse::<u32>().ok().map(|n| u8::from(n >= 1))
                    }
                })
                .collect(),
        );
    }

    Ok(RawMatrix { variants, reads, data })
}

/// Apply the prevalence and heteroplasmy filters to a [`RawMatrix`].
pub fn filter_matrix(
    raw: &RawMatrix,
    hf: &HfFilter,
    min_presence: usize,
    min_absence: usize,
    tract_rows: &std::collections::HashSet<String>,
) -> BinaryMatrix {
    let mut variants: Vec<String> = Vec::new();
    let mut data: Vec<Vec<Option<u8>>> = Vec::new();

    for (vid, row) in raw.variants.iter().zip(raw.data.iter()) {
        // On the VCF path the verdict needs only the name, so rows that cannot survive are rejected before anything is counted.
        if hf.rejects_by_name(vid, tract_rows) {
            continue;
        }

        let present = row.iter().filter(|c| **c == Some(1)).count();
        let absent = row.iter().filter(|c| **c == Some(0)).count();
        if present < min_presence || absent < min_absence {
            continue;
        }

        // Matrix-derived HF: alt calls over *covered* reads (present + absent), NOT row.len().
        let covered = present + absent;
        let freq = if covered > 0 { present as f64 / covered as f64 } else { 0.0 };

        if !hf.keeps(vid, freq, tract_rows) {
            continue;
        }

        variants.push(vid.clone());
        data.push(row.clone());
    }

    BinaryMatrix { variants, reads: raw.reads.clone(), data }
}

/// Collapse every known homopolymer tract in `raw` into sequence-level alleles.
pub fn collapse_tracts(raw: RawMatrix, tracts: &[crate::tracts::Tract]) -> RawMatrix {
    let spans = crate::scite::variant_spans(&raw.variants);
    let n_reads = raw.reads.len();

    // Rows inside each tract, with the edit they apply, in row order.
    let mut members: Vec<Vec<(usize, TractOp)>> = (0..tracts.len()).map(|_| Vec::new()).collect();
    let mut disabled = vec![false; tracts.len()];
    for (i, vid) in raw.variants.iter().enumerate() {
        let Some(t) = crate::tracts::containing_tract_in(tracts, &spans[i]) else {
            continue;
        };
        match tract_op(vid, &tracts[t]) {
            // A name that is not a plain SNV/insertion/deletion is treated like a straddling variant: it keeps its own row.
            ParsedOp::NotCollapsible => {}
            ParsedOp::RefMismatch => {
                if !disabled[t] {
                    warn!(
                        "[lineage] {} REF does not match the {} tract reference ({}); \
                         this is not rCRS, so the tract is left uncollapsed",
                        vid, tracts[t].name, tracts[t].ref_seq
                    );
                }
                disabled[t] = true;
            }
            ParsedOp::Op(op) => members[t].push((i, op)),
        }
    }

    // Per tract, the sequence each read carries over the tract.
    let mut collapsed: Vec<Option<String>> = Vec::new();
    let mut tract_rows: Vec<Vec<(String, Vec<Option<u8>>)>> =
        (0..tracts.len()).map(|_| Vec::new()).collect();
    for (t, tract) in tracts.iter().enumerate() {
        if disabled[t] || members[t].is_empty() {
            members[t].clear();
            continue;
        }
        let (mut n_partial, mut n_conflict) = (0usize, 0usize);
        collapsed.clear();
        for read in 0..n_reads {
            let mut missing = 0usize;
            let mut ops: Vec<&TractOp> = Vec::new();
            for (row, op) in &members[t] {
                match raw.data[*row][read] {
                    None => missing += 1,
                    Some(1) => ops.push(op),
                    _ => {}
                }
            }
            if missing == members[t].len() {
                // Uncovered everywhere in the tract; nothing new is lost.
                collapsed.push(None);
            } else if missing > 0 {
                n_partial += 1;
                collapsed.push(None);
            } else {
                match apply_tract_ops(tract.ref_seq, &ops) {
                    Some(seq) => collapsed.push(Some(seq)),
                    None => {
                        n_conflict += 1;
                        collapsed.push(None);
                    }
                }
            }
        }

        // One row per distinct non-reference sequence, sorted by that sequence.
        let alts: std::collections::BTreeSet<&String> = collapsed
            .iter()
            .flatten()
            .filter(|s| s.as_str() != tract.ref_seq)
            .collect();
        tract_rows[t] = alts
            .iter()
            .map(|alt| {
                let row = collapsed
                    .iter()
                    .map(|s| s.as_ref().map(|s| u8::from(s == *alt)))
                    .collect();
                (crate::tracts::tract_allele_name(tract, alt), row)
            })
            .collect();

        info!(
            "[lineage] Tract {} ({}-{}): {} row(s) collapsed into {} allele(s); \
             {} read(s) set to missing for partial coverage, {} for conflicting calls.",
            tract.name,
            tract.start,
            tract.end,
            members[t].len(),
            tract_rows[t].len(),
            n_partial,
            n_conflict
        );
        if n_conflict > 0 {
            warn!(
                "[lineage] Tract {}: {} read(s) carry calls that cannot be applied together \
                 (overlapping deletions, two insertions at one base, or contradictory \
                 substitutions) and were set to missing.",
                tract.name, n_conflict
            );
        }
    }

    // Reassemble in input order, putting each tract's alleles where its first row was.
    let owner: HashMap<usize, usize> = members
        .iter()
        .enumerate()
        .flat_map(|(t, rows)| rows.iter().map(move |(row, _)| (*row, t)))
        .collect();
    let mut variants = Vec::with_capacity(raw.variants.len());
    let mut data = Vec::with_capacity(raw.data.len());
    for (i, (vid, row)) in raw.variants.iter().zip(raw.data.iter()).enumerate() {
        match owner.get(&i) {
            Some(&t) => {
                if members[t][0].0 == i {
                    for (name, cells) in &tract_rows[t] {
                        variants.push(name.clone());
                        data.push(cells.clone());
                    }
                }
            }
            None => {
                variants.push(vid.clone());
                data.push(row.clone());
            }
        }
    }

    RawMatrix { variants, reads: raw.reads, data }
}

/// One edit a carried variant applies to a tract's reference sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
enum TractOp {
    /// Replace the base at `offset`.
    Sub { offset: usize, base: u8 },
    /// Insert `seq` after `offset`.
    Ins { offset: usize, seq: Vec<u8> },
    /// Delete `len` bases starting at `offset` (the base after the anchor).
    Del { offset: usize, len: usize },
}

/// What a variant name means for a tract.
enum ParsedOp {
    Op(TractOp),
    /// Not a plain SNV, insertion or deletion — left as its own row.
    NotCollapsible,
    /// The REF spelled in the name is not what the tract reference has there.
    RefMismatch,
}

/// Split `m.<pos><ref>><alt>` into its three fields, requiring both alleles to be plain base strings.
fn parse_variant_name(variant: &str) -> Option<(u64, &str, &str)> {
    let rest = variant.strip_prefix("m.")?;
    let digits_len = rest.chars().take_while(|c| c.is_ascii_digit()).count();
    let pos = rest[..digits_len].parse::<u64>().ok()?;
    let (ref_allele, alt_allele) = rest[digits_len..].split_once('>')?;
    let bases = |a: &str| !a.is_empty() && a.bytes().all(|b| b.is_ascii_alphabetic());
    (bases(ref_allele) && bases(alt_allele)).then_some((pos, ref_allele, alt_allele))
}

/// Turn one in-tract variant name into the edit it applies to `tract.ref_seq`.
fn tract_op(variant: &str, tract: &crate::tracts::Tract) -> ParsedOp {
    let Some((pos, ref_allele, alt_allele)) = parse_variant_name(variant) else {
        return ParsedOp::NotCollapsible;
    };
    // `containing_tract_in` already placed the whole REF span inside the tract.
    let offset = (pos - tract.start) as usize;
    let (r, a) = (ref_allele.as_bytes(), alt_allele.as_bytes());
    if !tract.ref_seq.as_bytes()[offset..offset + r.len()].eq_ignore_ascii_case(r) {
        return ParsedOp::RefMismatch;
    }

    if r.len() == 1 && a.len() == 1 {
        ParsedOp::Op(TractOp::Sub { offset, base: a[0].to_ascii_uppercase() })
    } else if r.len() == 1 && a.len() > 1 && a[0].eq_ignore_ascii_case(&r[0]) {
        ParsedOp::Op(TractOp::Ins {
            offset,
            seq: a[1..].to_ascii_uppercase(),
        })
    } else if a.len() == 1 && r.len() > 1 && r[0].eq_ignore_ascii_case(&a[0]) {
        // The anchor is kept; the bases after it go away.
        ParsedOp::Op(TractOp::Del { offset: offset + 1, len: r.len() - 1 })
    } else {
        ParsedOp::NotCollapsible
    }
}

/// Apply every operation a read carries to `ref_seq` in one left-to-right walk.
fn apply_tract_ops(ref_seq: &str, ops: &[&TractOp]) -> Option<String> {
    let bases = ref_seq.as_bytes();
    let n = bases.len();
    let mut sub: Vec<Option<u8>> = vec![None; n];
    let mut ins: Vec<Option<&[u8]>> = vec![None; n];
    let mut del = vec![false; n];

    for op in ops {
        match op {
            TractOp::Sub { offset, base } => {
                // Two *different* substitutions at one offset contradict; the same one twice is just a duplicated row.
                if sub[*offset].is_some_and(|b| b != *base) {
                    return None;
                }
                sub[*offset] = Some(*base);
            }
            TractOp::Ins { offset, seq } => {
                if ins[*offset].is_some() {
                    return None;
                }
                ins[*offset] = Some(seq);
            }
            TractOp::Del { offset, len } => {
                for o in *offset..(*offset + *len) {
                    if o >= n || del[o] {
                        return None;
                    }
                    del[o] = true;
                }
            }
        }
    }
    // A deletion that swallows a substituted base is a contradiction too.
    if (0..n).any(|o| del[o] && sub[o].is_some()) {
        return None;
    }

    let mut out: Vec<u8> = Vec::with_capacity(n + 4);
    for o in 0..n {
        if !del[o] {
            out.push(sub[o].unwrap_or(bases[o]));
        }
        if let Some(seq) = ins[o] {
            out.extend_from_slice(seq);
        }
    }
    String::from_utf8(out).ok()
}

/// Near-fixed variants dropped by the HF ceiling, not for lack of support.
pub fn collect_fixed_homoplasmies(
    raw: &RawMatrix,
    max_hf: f64,
    min_presence: usize,
    exclude: &std::collections::HashSet<String>,
) -> Vec<String> {
    let mut fixed = Vec::new();
    for (vid, row) in raw.variants.iter().zip(raw.data.iter()) {
        if exclude.contains(vid) {
            continue;
        }
        let present = row.iter().filter(|c| **c == Some(1)).count();
        let absent = row.iter().filter(|c| **c == Some(0)).count();
        let covered = present + absent;
        let freq = if covered > 0 { present as f64 / covered as f64 } else { 0.0 };
        if present >= min_presence && freq >= max_hf {
            fixed.push(vid.clone());
        }
    }
    fixed
}

/// Write `cleaned` to `out_path`, then append one all-present row per name in `fixed`.
fn write_cleaned_matrix_all_variants(
    cleaned: &crate::scite::CleanedMatrix,
    fixed: &[String],
    out_path: &str,
) -> Result<()> {
    let mut matrix = cleaned.clone();
    let n_reads = matrix.reads.len();
    for name in fixed {
        matrix.variants.push(name.clone());
        matrix.data.push(vec![1u8; n_reads]);
    }
    crate::scite::write_cleaned_matrix(&matrix, out_path)
}


/// Collapse reads with identical binary profiles into unique haplotypes.
/// Haplotypes are sorted by ascending mutation count (number of 1-bits), then lexicographically by
/// profile, so that the all-zero haplotype (if present) is always `H0000`.
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
pub fn four_gamete_test(matrix: &HaplotypeMatrix) -> Vec<(String, String)> {
    let n_var = matrix.variants.len();
    let profiles: Vec<&Vec<Option<u8>>> =
        matrix.haplotypes.iter().map(|h| &h.profile).collect();
    let mut violations = Vec::new();

    for i in 0..n_var {
        for j in (i + 1)..n_var {
            let mut seen = [false; 4]; // index = a*2 + b for states a,b ∈ {0,1}
            for profile in &profiles {
                // A haplotype only informs the pair if it covers BOTH sites; missing calls are excluded rather than read as ref.
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
    // Always represent the ancestral haplotype (the tree root: no heteroplasmic variants).
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
pub fn resolve_error_rates(
    data_type: &str,
    fp_rate: Option<f64>,
    fn_rate: Option<f64>,
) -> (f64, f64) {
    let (default_fp, default_fn) = match data_type {
        "pacbio" => (0.005, 0.05),
        "ont-r9" | "ont-r10" | "ont-denoised" => (0.005, 0.02),
        // Unknown data types fall back to the ONT BAM-matrix rates.
        _ => (0.005, 0.02),
    };

    (
        fp_rate.unwrap_or(default_fp),
        fn_rate.unwrap_or(default_fn),
    )
}

/// An explicit flag is used verbatim, including 0.
pub fn resolve_hp_rates(
    fp_rate: f64,
    fn_rate: f64,
    hp_fp_rate: Option<f64>,
    hp_fn_rate: Option<f64>,
) -> (f64, f64) {
    (
        hp_fp_rate.unwrap_or((2.0 * fp_rate).min(0.5)),
        hp_fn_rate.unwrap_or((2.0 * fn_rate).min(0.5)),
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
    fp_rate: f64,
    fn_rate: f64,
    hp_rates: (f64, f64),
    tracts: &[crate::tracts::Tract],
    mcmc_iterations: usize,
    mcmc_chains: usize,
    mcmc_seed: u64,
    root_block: Option<crate::rootblock::RootBlockConfig>,
    debug: bool,
    output_prefix: &str,
) -> Result<()> {
    // No `data_type` parameter: the rates arrive already resolved by `resolve_error_rates`.
    info!("Starting lineage analysis...");
    info!("SCITE fp={fp_rate}, fn={fn_rate}, hp_fp={}, hp_fn={}", hp_rates.0, hp_rates.1);

    // ── Step 1: decide how HF is judged, then filter the binary matrix ────────
    let hf = match vcf_file {
        Some(path) => {
            info!("[1/6] Parsing VCF: {}", path);
            HfFilter::FromVcf { map: parse_vcf(path, min_hf, max_hf)?, min: min_hf, max: max_hf }
        }
        None => {
            info!("[1/6] No VCF provided; HF recomputed from the matrix itself.");
            HfFilter::FromMatrix { min: min_hf, max: max_hf }
        }
    };

    info!("[1/6] Loading and filtering matrix: {}", matrix_file);
    // The filter below and the fixed-homoplasmy scan further down both read these same rows.
    let raw = read_raw_matrix(
        File::open(matrix_file)
            .with_context(|| format!("Cannot read matrix CSV: {matrix_file}"))?,
    )?;
    // An empty slice is the identity and leaves `tract_rows` empty.
    if tracts.is_empty() {
        info!("[lineage] HP tract collapse: off (pass --collapse-hp to enable)");
    } else {
        info!("[lineage] HP tract collapse: on ({} tracts)", tracts.len());
    }
    let raw = collapse_tracts(raw, tracts);
    let tract_rows: std::collections::HashSet<String> = raw
        .variants
        .iter()
        .filter(|v| crate::tracts::is_tract_allele(v, tracts))
        .cloned()
        .collect();
    let binary = filter_matrix(&raw, &hf, min_presence, min_absence, &tract_rows);
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
        anyhow::bail!(
            "No haplotypes remain: every distinct read profile is backed by fewer \
             than {min_reads} reads. This usually means the input is too shallow \
             or too noisy; widening --min-hf / --max-hf admits more variants and \
             so more distinguishable haplotypes."
        );
    }
    // Diagnostic: variant pairs violating the infinite-sites assumption hint at recurrent/back mutation or recombination.
    let violations = four_gamete_test(&hap_matrix);
    if !violations.is_empty() {
        info!(
            "[2/6] Four-gamete test: {} variant pair(s) violate the infinite-sites assumption.",
            violations.len()
        );
    }
    // The pre-SCITE haplotype map is a debug artifact.
    if debug {
        let hmap_path = format!("{}.raw_haplotype_map.tsv", output_prefix);
        write_haplotype_map(&hap_matrix, &hmap_path)?;
    }

    // ── Step 3: SCITE mutation-tree search (NJ-informed MCMC) ─────────────────
    info!("[3/6] Running SCITE mutation-tree search...");
    // Returns the cleaned reads already deduplicated into haplotypes.
    let out = crate::scite::run_scite_pipeline(
        &binary,
        fp_rate,
        fn_rate,
        hp_rates,
        mcmc_iterations,
        mcmc_chains,
        mcmc_seed,
        min_reads,
        root_block.as_ref(),
        debug,
        output_prefix,
    )?;

    // Near-fixed homoplasmies were removed by the HF ceiling before the tree was built.
    let kept: std::collections::HashSet<String> = binary.variants.iter().cloned().collect();
    let fixed_homoplasmies = collect_fixed_homoplasmies(&raw, max_hf, min_presence, &kept);
    if !fixed_homoplasmies.is_empty() {
        info!(
            "[3/6] Recovered {} near-fixed homoplasmic variant(s) excluded by the HF ceiling.",
            fixed_homoplasmies.len()
        );
    }
    write_cleaned_matrix_all_variants(
        &out.cleaned,
        &fixed_homoplasmies,
        &format!("{output_prefix}.cleaned_matrix_all_variants.csv"),
    )?;
    let cleaned_hap_matrix = out.cleaned_haplotypes;

    // ── Step 4: write final haplotype map (filtered by SCITE) ─────────────────
    info!(
        "[4/6] Cleaned reads deduplicated into {} haplotypes across {} variants.",
        cleaned_hap_matrix.haplotypes.len(),
        cleaned_hap_matrix.variants.len()
    );
    if cleaned_hap_matrix.haplotypes.is_empty() {
        anyhow::bail!(
            "No haplotypes remain: every distinct read profile is backed by fewer \
             than {min_reads} reads. This usually means the input is too shallow \
             or too noisy; widening --min-hf / --max-hf admits more variants and \
             so more distinguishable haplotypes."
        );
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

    /// An `HfFilter::FromVcf` holding exactly `names`.
    fn from_vcf(names: &[&str]) -> HfFilter {
        HfFilter::FromVcf {
            map: names.iter().map(|n| (n.to_string(), 0.5)).collect(),
            min: 0.0,
            max: 1.0,
        }
    }

    /// Parse `csv` and filter it, i.e. what the old `parse_binary_matrix` did in one call.
    fn filter_from_csv(
        csv: &str,
        hf: &HfFilter,
        min_presence: usize,
        min_absence: usize,
    ) -> BinaryMatrix {
        let raw = read_raw_matrix(Cursor::new(csv.to_string())).unwrap();
        filter_matrix(&raw, hf, min_presence, min_absence, &Default::default())
    }

    fn hp(id: &str, profile: Vec<Option<u8>>, count: usize) -> Haplotype {
        Haplotype {
            id: id.to_string(),
            profile,
            reads: vec![],
            count,
        }
    }

    #[test]
    fn pick_outgroup_prefers_best_covered_ancestral_haplotype() {
        let hm = HaplotypeMatrix {
            variants: vec!["v0".into(), "v1".into()],
            haplotypes: vec![
                hp("H0000", vec![None, None], 1),
                hp("H0001", vec![Some(0), Some(0)], 5),
                hp("H0002", vec![Some(1), Some(1)], 5),
            ],
        };
        assert_eq!(pick_outgroup(&hm).map(|h| h.id.as_str()), Some("H0001"));
    }

    #[test]
    fn pick_outgroup_rejects_zero_coverage_ancestor() {
        let hm = HaplotypeMatrix {
            variants: vec!["v0".into(), "v1".into()],
            haplotypes: vec![
                hp("H0000", vec![None, None], 1),
                hp("H0001", vec![Some(1), Some(1)], 5),
            ],
        };
        assert_eq!(pick_outgroup(&hm).map(|h| h.id.as_str()), None);
    }

    #[test]
    fn pick_outgroup_tie_keeps_earlier_haplotype() {
        let hm = HaplotypeMatrix {
            variants: vec!["v0".into(), "v1".into()],
            haplotypes: vec![
                hp("H0000", vec![Some(0), Some(0)], 5),
                hp("H0001", vec![Some(0), Some(0)], 5),
                hp("H0002", vec![Some(1), Some(1)], 5),
            ],
        };
        assert_eq!(pick_outgroup(&hm).map(|h| h.id.as_str()), Some("H0000"));
    }

    #[test]
    fn resolve_error_rates_presets_differ_per_data_type() {
        assert_eq!(resolve_error_rates("pacbio", None, None), (0.005, 0.05));
        let ont = (0.005, 0.02);
        assert_eq!(resolve_error_rates("ont-denoised", None, None), ont);
        assert_eq!(resolve_error_rates("ont-r10", None, None), ont);
        assert_eq!(resolve_error_rates("ont-r9", None, None), ont);
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
    fn resolve_hp_rates_defaults_to_twice_the_global_rate_capped_at_half() {
        assert_eq!(resolve_hp_rates(0.005, 0.05, None, None), (0.01, 0.1));
        // 2 × 0.4 = 0.8, which the cap pulls back to 0.5.
        assert_eq!(resolve_hp_rates(0.4, 0.01, None, None), (0.5, 0.02));
        // An explicit flag is used verbatim, including a value above the cap and 0.
        assert_eq!(resolve_hp_rates(0.005, 0.05, Some(0.9), None), (0.9, 0.1));
        assert_eq!(resolve_hp_rates(0.005, 0.05, None, Some(0.0)), (0.01, 0.0));
    }

    #[test]
    fn load_and_filter_matrix_preserves_missing_values() {
        let csv = "variant,r1,r2,r3\nm.100A>G,1,0,\nm.200C>T,0,1,1\n";
        let hf = from_vcf(&["m.100A>G", "m.200C>T"]);
        let matrix = filter_from_csv(csv, &hf, 1, 1);

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
        let hf = from_vcf(&["m.100A>G"]);
        let matrix = filter_from_csv(csv, &hf, 1, 1);
        assert!(matrix.variants.is_empty());
    }

    #[test]
    fn collect_fixed_homoplasmies_finds_near_fixed_sites_not_already_kept() {
        // m.100: freq 4/4 = 1.0, supported -> recovered. m.200: freq 3/4 = 0.75 -> heteroplasmic, not recovered. m.300.
        let csv = "variant,r1,r2,r3,r4\n\
                   m.100A>G,1,1,1,1\n\
                   m.200C>T,1,1,1,0\n\
                   m.300G>A,1,,,\n\
                   m.400T>C,1,1,1,1\n";
        let raw = read_raw_matrix(Cursor::new(csv.to_string())).unwrap();
        let mut kept = std::collections::HashSet::new();
        kept.insert("m.400T>C".to_string());
        let fixed = collect_fixed_homoplasmies(&raw, 0.95, 2, &kept);
        assert_eq!(fixed, vec!["m.100A>G".to_string()]);
    }

    #[test]
    fn write_cleaned_matrix_all_variants_appends_all_present_rows() {
        // Built in memory; the expected text is unchanged.
        let cleaned = crate::scite::CleanedMatrix {
            variants: vec!["m.200C>T".to_string()],
            reads: vec!["r1".to_string(), "r2".to_string(), "r3".to_string()],
            data: vec![vec![1, 0, 1]],
            attachment: vec![0, 0, 0],
        };
        let out = std::env::temp_dir().join("himito_test_cleaned_all.csv");
        write_cleaned_matrix_all_variants(
            &cleaned,
            &["m.100A>G".to_string()],
            out.to_str().unwrap(),
        )
        .unwrap();
        let content = std::fs::read_to_string(&out).unwrap();
        std::fs::remove_file(&out).ok();
        assert_eq!(
            content,
            "variant,r1,r2,r3\nm.200C>T,1,0,1\nm.100A>G,1,1,1\n"
        );
    }

    #[test]
    fn deduplicate_treats_different_missing_patterns_as_different_haplotypes() {
        let matrix = BinaryMatrix {
            variants: vec!["m.100A>G".to_string()],
            reads: vec!["r1".to_string(), "r2".to_string()],
            // Both reads are "ref-or-missing" if you ignore the distinction.
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
        assert!((dist.data[0][1] - 0.5).abs() < 1e-12);
        assert!((dist.data[1][0] - 0.5).abs() < 1e-12);
    }

    fn hap2(id: &str, a: Option<u8>, b: Option<u8>) -> Haplotype {
        Haplotype { id: id.into(), profile: vec![a, b], reads: vec![], count: 1 }
    }

    #[test]
    fn hamming_distance_matrix_fills_uncovered_pairs_with_mean_distance() {
        // H0 and H1 share no covered site.
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["A".into(), "B".into()],
            haplotypes: vec![
                hap2("H0", Some(1), None),
                hap2("H1", None, Some(1)),
                hap2("H2", Some(0), Some(0)),
            ],
        };
        let dist = hamming_distance_matrix(&hap_matrix);
        assert_eq!(dist.data[0][2], 1.0);
        assert_eq!(dist.data[2][0], 1.0);
        assert_eq!(dist.data[1][2], 1.0);
        assert_eq!(dist.data[2][1], 1.0);
        assert_eq!(dist.data[0][1], 1.0);
        assert_eq!(dist.data[1][0], 1.0);
        for row in &dist.data {
            for &v in row {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn hamming_distance_matrix_zero_fills_when_no_pair_has_shared_coverage() {
        // With no finite off-diagonal distance, the fill is 0.0.
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["A".into(), "B".into()],
            haplotypes: vec![
                hap2("H0", Some(1), None),
                hap2("H1", None, Some(1)),
            ],
        };
        let dist = hamming_distance_matrix(&hap_matrix);
        assert_eq!(dist.data[0][1], 0.0);
        assert_eq!(dist.data[1][0], 0.0);
        assert!(dist.data.iter().flatten().all(|v| v.is_finite()));
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
        // Only 01, 10, 11 are truly observed.
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
        // Old root (4) collapses; H2 hangs directly under the new root with its branch length summed through the old root.
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

    // Captures `warn!` so the missing-header case is observable.
    struct NormalizationLog;

    static NORMALIZATION_LOG: std::sync::Mutex<Vec<(std::thread::ThreadId, String)>> =
        std::sync::Mutex::new(Vec::new());

    impl log::Log for NormalizationLog {
        fn enabled(&self, metadata: &log::Metadata) -> bool {
            metadata.level() <= log::Level::Warn
        }

        fn log(&self, record: &log::Record) {
            if self.enabled(record.metadata()) {
                NORMALIZATION_LOG.lock().unwrap().push((
                    std::thread::current().id(),
                    record.args().to_string(),
                ));
            }
        }

        fn flush(&self) {}
    }

    fn install_normalization_log() {
        static LOGGER: NormalizationLog = NormalizationLog;
        let _ = log::set_logger(&LOGGER);
        log::set_max_level(log::LevelFilter::Warn);
    }

    fn logs_during(f: impl FnOnce()) -> Vec<String> {
        let me = std::thread::current().id();
        let start = NORMALIZATION_LOG.lock().unwrap().len();
        f();
        NORMALIZATION_LOG.lock().unwrap()[start..]
            .iter()
            .filter(|(id, _)| *id == me)
            .map(|(_, line)| line.clone())
            .collect()
    }

    fn minimal_pass_vcf(extra_header: &str) -> String {
        format!(
            "\
##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
{extra_header}##FORMAT=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmic fraction\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chrM\t100\t.\tA\tG\t.\tPASS\t.\tHF\t0.3
"
        )
    }

    /// A VCF without `##HimitoNormalization` must still parse.
    #[test]
    fn parse_vcf_warns_once_when_header_lacks_himito_normalization() {
        install_normalization_log();
        let missing = std::env::temp_dir().join("himito_test_parse_vcf_no_norm.vcf");
        std::fs::write(&missing, minimal_pass_vcf("")).unwrap();
        let logs = logs_during(|| {
            parse_vcf(missing.to_str().unwrap(), 0.0, 1.0).unwrap();
        });
        std::fs::remove_file(&missing).ok();
        let hits: Vec<_> = logs
            .iter()
            .filter(|l| l.contains("VCF was not left-aligned"))
            .collect();
        assert_eq!(
            hits.len(),
            1,
            "expected one left-align warning, got {logs:?}"
        );
        assert!(
            hits[0].contains("without `--no-left-align`"),
            "warning should name --no-left-align, got {}",
            hits[0]
        );

        let present = std::env::temp_dir().join("himito_test_parse_vcf_with_norm.vcf");
        std::fs::write(
            &present,
            minimal_pass_vcf("##HimitoNormalization=left-aligned\n"),
        )
        .unwrap();
        let logs = logs_during(|| {
            parse_vcf(present.to_str().unwrap(), 0.0, 1.0).unwrap();
        });
        std::fs::remove_file(&present).ok();
        assert!(
            logs.iter().all(|l| !l.contains("VCF was not left-aligned")),
            "left-aligned header must not warn, got {logs:?}"
        );
    }

    /// Temp prefixes are stable across runs.
    fn prepare_debug_output_prefix(dir_name: &str) -> (String, String) {
        let dir = std::env::temp_dir().join(dir_name);
        std::fs::create_dir_all(&dir).unwrap();
        let matrix_path = dir.join("matrix.csv");
        std::fs::write(
            &matrix_path,
            "variant,r1,r2,r3,r4\n\
             m.100A>G,1,1,0,0\n\
             m.200C>T,1,0,1,0\n",
        )
        .unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        for suffix in [
            ".cleaned_matrix.csv",
            ".raw_haplotype_map.tsv",
            ".raw_variant_cooccurrence.tsv",
            ".cleaned_matrix_all_variants.csv",
            ".cleaned_haplotype_map.tsv",
            ".mutation_tree.tsv",
            ".variant_cooccurrence.tsv",
            ".molecule_summary.tsv",
            ".read_lineage.nwk",
            ".root_block.tsv",
        ] {
            std::fs::remove_file(format!("{prefix}{suffix}")).ok();
        }
        (matrix_path.to_str().unwrap().to_string(), prefix)
    }

    fn run_small_lineage(
        matrix: &str,
        prefix: &str,
        debug: bool,
        tracts: &[crate::tracts::Tract],
    ) {
        start(
            matrix,
            None,
            0.01,
            0.95,
            1,
            1,
            1,
            0.01,
            0.1,
            (0.01, 0.1),
            tracts,
            20,
            1,
            7,
            None,
            debug,
            prefix,
        )
        .unwrap();
    }

    #[test]
    fn debug_false_omits_raw_files_and_cleaned_matrix() {
        let (matrix, prefix) = prepare_debug_output_prefix("himito_t3_debug_off");
        run_small_lineage(&matrix, &prefix, false, &[]);
        for suffix in [
            ".cleaned_matrix.csv",
            ".raw_haplotype_map.tsv",
            ".raw_variant_cooccurrence.tsv",
        ] {
            let path = format!("{prefix}{suffix}");
            assert!(
                !std::path::Path::new(&path).exists(),
                "{path} must not be written when debug is false"
            );
        }
    }

    #[test]
    fn debug_true_writes_both_raw_files() {
        let (matrix, prefix) = prepare_debug_output_prefix("himito_t3_debug_on");
        run_small_lineage(&matrix, &prefix, true, &[]);
        for suffix in [".raw_haplotype_map.tsv", ".raw_variant_cooccurrence.tsv"] {
            let path = format!("{prefix}{suffix}");
            assert!(
                std::path::Path::new(&path).exists(),
                "{path} must be written when debug is true"
            );
        }
    }

    /// Two HVS1-C rows plus one variant outside every tract.
    fn write_hvs1_opt_in_matrix(dir_name: &str) -> (String, String) {
        let dir = std::env::temp_dir().join(dir_name);
        std::fs::create_dir_all(&dir).unwrap();
        let matrix_path = dir.join("matrix.csv");
        std::fs::write(
            &matrix_path,
            "variant,r1,r2,r3,r4\n\
             m.16183A>AC,1,0,0,0\n\
             m.16189T>TC,0,1,0,0\n\
             m.100A>G,1,1,0,0\n",
        )
        .unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        std::fs::remove_file(format!("{prefix}.cleaned_matrix_all_variants.csv")).ok();
        (matrix_path.to_str().unwrap().to_string(), prefix)
    }

    fn cleaned_names(prefix: &str) -> String {
        std::fs::read_to_string(format!("{prefix}.cleaned_matrix_all_variants.csv")).unwrap()
    }

    #[test]
    fn empty_tracts_keep_hvs1_names_and_known_tracts_collapse_them() {
        let (matrix, off) = write_hvs1_opt_in_matrix("himito_t9_collapse_off");
        run_small_lineage(&matrix, &off, false, &[]);
        let off_csv = cleaned_names(&off);
        assert!(
            off_csv.contains("m.16183A>AC"),
            "empty tracts must keep m.16183A>AC, got {off_csv}"
        );
        assert!(
            off_csv.contains("m.16189T>TC"),
            "empty tracts must keep m.16189T>TC, got {off_csv}"
        );
        assert!(
            !off_csv.contains("m.16183ACCCCCTCCCC>"),
            "empty tracts must not emit a tract allele, got {off_csv}"
        );

        let (matrix, on) = write_hvs1_opt_in_matrix("himito_t9_collapse_on");
        run_small_lineage(&matrix, &on, false, crate::tracts::KNOWN_TRACTS);
        let on_csv = cleaned_names(&on);
        assert!(
            !on_csv.contains("m.16183A>AC"),
            "KNOWN_TRACTS must drop the constituent name m.16183A>AC, got {on_csv}"
        );
        assert!(
            !on_csv.contains("m.16189T>TC"),
            "KNOWN_TRACTS must drop the constituent name m.16189T>TC, got {on_csv}"
        );
        assert!(
            on_csv.contains("m.16183ACCCCCTCCCC>"),
            "KNOWN_TRACTS must emit an HVS1-C tract allele, got {on_csv}"
        );
    }

    #[test]
    fn empty_tracts_drop_a_tract_looking_row_absent_from_the_vcf() {
        // With no tracts, `tract_rows` is empty, so `HfFilter::FromVcf` drops a name that is not in the VCF map.
        let dir = std::env::temp_dir().join("himito_t9_vcf_empty_tracts");
        std::fs::create_dir_all(&dir).unwrap();
        let matrix_path = dir.join("matrix.csv");
        std::fs::write(
            &matrix_path,
            "variant,r1,r2,r3,r4\n\
             m.16183A>AC,1,0,0,0\n\
             m.100A>G,1,1,0,0\n",
        )
        .unwrap();
        let vcf_path = dir.join("sample.vcf");
        std::fs::write(
            &vcf_path,
            "\
##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##HimitoNormalization=left-aligned
##FORMAT=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmic fraction\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chrM\t100\t.\tA\tG\t.\tPASS\t.\tHF\t0.5
",
        )
        .unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        std::fs::remove_file(format!("{prefix}.cleaned_matrix_all_variants.csv")).ok();
        start(
            matrix_path.to_str().unwrap(),
            Some(vcf_path.to_str().unwrap()),
            0.01,
            0.95,
            1,
            1,
            1,
            0.01,
            0.1,
            (0.01, 0.1),
            &[],
            20,
            1,
            7,
            None,
            false,
            &prefix,
        )
        .unwrap();
        let csv = cleaned_names(&prefix);
        assert!(
            csv.contains("m.100A>G"),
            "the VCF-listed variant must be kept, got {csv}"
        );
        assert!(
            !csv.contains("m.16183A>AC"),
            "a tract-looking row absent from the VCF must be dropped, got {csv}"
        );
        assert!(
            !csv.contains("m.16183ACCCCCTCCCC>"),
            "empty tracts must not rename the dropped row, got {csv}"
        );
    }
}

/// Parity checks: the with-VCF and without-VCF variant filtering must select the same variants.
#[cfg(test)]
mod tract_collapse_tests {
    use super::*;
    use crate::tracts::KNOWN_TRACTS;
    use std::collections::HashSet;

    /// rCRS 16183-16193: `A` + `C5` + `T` + `C4`.
    const HVS1: &str = "ACCCCCTCCCC";
    /// The tract with T16189C only: a single run of 10 C's.
    const C10: &str = "ACCCCCCCCCC";
    /// T16189C plus one extra C, spelled either side of the interrupting T.
    const C11: &str = "ACCCCCCCCCCC";

    /// `1` = alt, `0` = ref, `.` = uncovered.
    fn cells(spec: &str) -> Vec<Option<u8>> {
        spec.chars()
            .map(|c| match c {
                '1' => Some(1),
                '0' => Some(0),
                '.' => None,
                other => panic!("bad cell {other:?}"),
            })
            .collect()
    }

    fn raw(variants: &[&str], rows: &[&str]) -> RawMatrix {
        let n_reads = rows[0].chars().count();
        RawMatrix {
            variants: variants.iter().map(|v| v.to_string()).collect(),
            reads: (1..=n_reads).map(|i| format!("r{i}")).collect(),
            data: rows.iter().map(|r| cells(r)).collect(),
        }
    }

    fn hvs1(alt: &str) -> String {
        format!("m.16183{HVS1}>{alt}")
    }

    #[test]
    fn collapse_tracts_merges_equivalent_spellings_of_one_c_run() {
        // Case 1. r1 carries the extra C before the interrupting T, r2 after it.
        let out = collapse_tracts(
            raw(
                &["m.16183A>AC", "m.16189T>TC", "m.16189T>C"],
                &["1000", "0100", "1111"],
            ),
            KNOWN_TRACTS,
        );
        assert_eq!(out.variants, vec![hvs1(C10), hvs1(C11)]);
        assert_eq!(out.data, vec![cells("0011"), cells("1100")]);
        assert_eq!(out.reads, vec!["r1", "r2", "r3", "r4"]);
    }

    #[test]
    fn collapse_tracts_empty_slice_returns_the_input_unchanged() {
        // The case-1 fixture, with no tracts to collapse.
        let input = raw(
            &["m.16183A>AC", "m.16189T>TC", "m.16189T>C"],
            &["1000", "0100", "1111"],
        );
        let expected = input.clone();
        let out = collapse_tracts(input, &[]);
        assert_eq!(out.variants, expected.variants);
        assert_eq!(out.reads, expected.reads);
        assert_eq!(out.data, expected.data);
    }

    #[test]
    fn collapse_tracts_gives_reference_reads_zero_and_emits_no_reference_row() {
        // Case 2. r3 carries nothing inside the tract.
        let out = collapse_tracts(
            raw(&["m.16183A>AC", "m.16189T>C"], &["100", "010"]),
            KNOWN_TRACTS,
        );
        // r1 keeps the T and gains a C; r2 is C10.
        assert_eq!(out.variants, vec![hvs1(C10), hvs1("ACCCCCCTCCCC")]);
        assert_eq!(out.data, vec![cells("010"), cells("100")]);
    }

    #[test]
    fn collapse_tracts_marks_partially_covered_reads_missing() {
        // Case 3. r2 is uncovered on one constituent row only, r3 on both.
        let out = collapse_tracts(
            raw(&["m.16183A>AC", "m.16189T>C"], &["1..", "01."]),
            KNOWN_TRACTS,
        );
        assert_eq!(out.variants, vec![hvs1("ACCCCCCTCCCC")]);
        assert_eq!(out.data, vec![cells("1..")]);
    }

    #[test]
    fn collapse_tracts_sets_conflicting_reads_missing() {
        // Case 4, all three conflict kinds.

        // Two insertions after the same offset.
        let out = collapse_tracts(
            raw(&["m.16183A>ACC", "m.16183A>AC"], &["10", "11"]),
            KNOWN_TRACTS,
        );
        assert_eq!(out.variants, vec![hvs1("ACCCCCCTCCCC")]);
        assert_eq!(out.data, vec![vec![None, Some(1)]]);

        // A deletion whose deleted base is also substituted.
        let out = collapse_tracts(
            raw(&["m.16188CT>C", "m.16189T>C"], &["10", "11"]),
            KNOWN_TRACTS,
        );
        assert_eq!(out.variants, vec![hvs1(C10)]);
        assert_eq!(out.data, vec![vec![None, Some(1)]]);

        // Two different substitutions at one offset.
        let out = collapse_tracts(
            raw(&["m.16189T>A", "m.16189T>C"], &["10", "11"]),
            KNOWN_TRACTS,
        );
        assert_eq!(out.variants, vec![hvs1(C10)]);
        assert_eq!(out.data, vec![vec![None, Some(1)]]);
    }

    #[test]
    fn collapse_tracts_leaves_variants_straddling_the_boundary_alone() {
        // Case 5.
        let out = collapse_tracts(
            raw(&["m.16189T>C", "m.16193CA>C"], &["11", "01"]),
            KNOWN_TRACTS,
        );
        assert_eq!(out.variants, vec![hvs1(C10), "m.16193CA>C".to_string()]);
        assert_eq!(out.data, vec![cells("11"), cells("01")]);
    }

    #[test]
    fn collapse_tracts_leaves_non_indel_non_snv_names_alone() {
        // Case 5, second half of the rule: an MNV is not collapsible.
        let out = collapse_tracts(
            raw(&["m.16188CT>GA", "m.16189T>C"], &["10", "11"]),
            KNOWN_TRACTS,
        );
        assert_eq!(out.variants, vec!["m.16188CT>GA".to_string(), hvs1(C10)]);
        assert_eq!(out.data, vec![cells("10"), cells("11")]);
    }

    #[test]
    fn collapse_tracts_skips_a_tract_whose_reference_does_not_match() {
        // Case 6. rCRS has C at 16185, so `m.16185G>A` says this is not rCRS; the whole tract passes through untouched.
        let input = raw(&["m.16185G>A", "m.16189T>C"], &["10", "11"]);
        let expected = input.clone();
        assert_eq!(collapse_tracts(input, KNOWN_TRACTS), expected);
    }

    #[test]
    fn collapse_tracts_keeps_other_rows_in_place() {
        // Case 7.
        let out = collapse_tracts(
            raw(
                &["m.100A>G", "m.16189T>C", "m.200C>T", "m.16183A>AC", "m.300G>A"],
                &["10", "11", "01", "10", "11"],
            ),
            KNOWN_TRACTS,
        );
        assert_eq!(
            out.variants,
            vec![
                "m.100A>G".to_string(),
                hvs1(C10),
                hvs1(C11),
                "m.200C>T".to_string(),
                "m.300G>A".to_string(),
            ]
        );
        assert_eq!(
            out.data,
            vec![cells("10"), cells("01"), cells("10"), cells("01"), cells("11")]
        );
    }

    /// The rows `lineage::start` feeds the filter: three tract calls plus one ordinary variant.
    fn collapsed_fixture() -> (RawMatrix, HashSet<String>) {
        let csv = "variant,r1,r2,r3,r4\n\
                   m.16183A>AC,1,0,0,0\n\
                   m.16189T>TC,0,1,0,0\n\
                   m.16189T>C,1,1,1,1\n\
                   m.3243A>G,1,1,0,0\n";
        let raw = collapse_tracts(
            read_raw_matrix(std::io::Cursor::new(csv)).unwrap(),
            KNOWN_TRACTS,
        );
        let tract_rows = raw
            .variants
            .iter()
            .filter(|v| crate::tracts::is_tract_allele(v, KNOWN_TRACTS))
            .cloned()
            .collect();
        (raw, tract_rows)
    }

    #[test]
    fn tract_rows_are_judged_on_matrix_frequency_on_the_vcf_path() {
        // Step 5.
        let (raw, tract_rows) = collapsed_fixture();
        assert_eq!(tract_rows.len(), 2);
        let mut map = HfMap::new();
        map.insert("m.3243A>G".to_string(), 0.5);

        let kept = filter_matrix(
            &raw,
            &HfFilter::FromVcf { map: map.clone(), min: 0.2, max: 0.95 },
            1,
            1,
            &tract_rows,
        );
        assert_eq!(
            kept.variants,
            vec![hvs1(C10), hvs1(C11), "m.3243A>G".to_string()]
        );

        // Both tract alleles sit at 2/4 = 0.50.
        let kept = filter_matrix(
            &raw,
            &HfFilter::FromVcf { map, min: 0.2, max: 0.4 },
            1,
            1,
            &tract_rows,
        );
        assert_eq!(kept.variants, vec!["m.3243A>G".to_string()]);
    }

    #[test]
    fn tract_rows_are_judged_on_matrix_frequency_on_the_no_vcf_path_too() {
        let (raw, tract_rows) = collapsed_fixture();
        let kept = filter_matrix(
            &raw,
            &HfFilter::FromMatrix { min: 0.2, max: 0.95 },
            1,
            1,
            &tract_rows,
        );
        assert_eq!(
            kept.variants,
            vec![hvs1(C10), hvs1(C11), "m.3243A>G".to_string()]
        );
        let kept = filter_matrix(
            &raw,
            &HfFilter::FromMatrix { min: 0.2, max: 0.4 },
            1,
            1,
            &tract_rows,
        );
        assert!(kept.variants.is_empty());
    }
}

#[cfg(test)]
mod vcf_parity_tests {
    use super::*;
    use std::io::Cursor;

    /// Build a `variant,r1..r10` CSV.
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

    /// Parse `csv` and filter it, i.e. what the old `parse_binary_matrix` did in one call.
    fn filter_from_csv(
        csv: &str,
        hf: &HfFilter,
        min_presence: usize,
        min_absence: usize,
    ) -> BinaryMatrix {
        let raw = read_raw_matrix(Cursor::new(csv.to_string())).unwrap();
        filter_matrix(&raw, hf, min_presence, min_absence, &Default::default())
    }

    #[test]
    fn no_vcf_matrix_frequency_filter_matches_vcf_hf_filter() {
        // present/absent/uncovered per variant, thresholds min_hf=0.2 max_hf=0.85.
        let csv = build_csv(&[
            ("m.100A>G", ["1", "1", "1", "0", "0", "0", "", "", "", ""]),
            ("m.200C>T", ["1", "0", "0", "0", "", "", "", "", "", ""]),
            ("m.300G>A", ["1", "1", "1", "1", "1", "1", "1", "1", "1", "0"]),
            ("m.400T>C", ["1", "0", "0", "0", "0", "0", "0", "0", "0", "0"]),
        ]);

        let (min_presence, min_absence, min_hf, max_hf) = (1usize, 1usize, 0.2f64, 0.85f64);

        // "With VCF": the map holds exactly the variants whose HF is in range, as parse_vcf would produce.
        let mut with_vcf_hf = HfMap::new();
        with_vcf_hf.insert("m.100A>G".to_string(), 0.50);
        with_vcf_hf.insert("m.200C>T".to_string(), 0.25);

        let with_vcf = filter_from_csv(
            &csv,
            &HfFilter::FromVcf { map: with_vcf_hf, min: min_hf, max: max_hf },
            min_presence,
            min_absence,
        );

        // "Without VCF": frequency recomputed from the matrix against the band.
        let without_vcf = filter_from_csv(
            &csv,
            &HfFilter::FromMatrix { min: min_hf, max: max_hf },
            min_presence,
            min_absence,
        );

        // Both keep exactly the two heteroplasmic variants, same order, same data.
        assert_eq!(without_vcf.variants, vec!["m.100A>G", "m.200C>T"]);
        assert_eq!(with_vcf.variants, without_vcf.variants);
        assert_eq!(with_vcf.data, without_vcf.data);
        assert_eq!(with_vcf.reads, without_vcf.reads);
    }

    #[test]
    fn no_vcf_denominator_excludes_uncovered_reads() {
        // Heteroplasmic among covered reads but "rare" if uncovered cells are counted.
        let csv = build_csv(&[("m.500A>C", ["1", "1", "0", "0", "", "", "", "", "", ""])]);
        let kept = filter_from_csv(
            &csv,
            &HfFilter::FromMatrix { min: 0.30, max: 0.85 },
            1,
            1,
        );
        // Covered-read denominator (0.50) retains it; the old row.len() denominator (0.20 < 0.30) would have dropped it.
        assert_eq!(kept.variants, vec!["m.500A>C"]);
    }

    #[test]
    fn a_vcf_matching_no_variants_keeps_nothing_rather_than_falling_back() {
        // BEHAVIOR CHANGE, deliberate.
        let csv = build_csv(&[("m.100A>G", ["1", "1", "1", "0", "0", "0", "", "", "", ""])]);
        let kept = filter_from_csv(
            &csv,
            &HfFilter::FromVcf { map: HfMap::new(), min: 0.2, max: 0.85 },
            1,
            1,
        );
        assert!(
            kept.variants.is_empty(),
            "an empty VCF map must keep nothing, got {:?}",
            kept.variants
        );
    }
}