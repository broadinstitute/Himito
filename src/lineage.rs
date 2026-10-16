use std::collections::HashMap;
use anyhow::{bail,Context, Result};
use rust_htslib::bcf::{Read, Reader}; // Read trait must be in scope for .records()
use std::io::{BufWriter, Write};
use std::fs::File;
use log::{info, warn};

/// Minimal reads a variant must be present in / absent from to be informative.
const MIN_PRESENCE: usize = 2;
const MIN_ABSENCE: usize = 1;
/// Minimal reads backing a reported haplotype.
const MIN_READS: usize = 3;

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

/// Pairwise Hamming distance over jointly covered sites; a pair sharing no covered site gets the mean of the others (0 if none).
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

/// Neighbor-Joining tree over the haplotypes, rooted at the ancestral haplotype's parent when there is one.
pub fn neighbor_joining(hap_matrix: &HaplotypeMatrix) -> Result<Tree> {
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

/// The ancestral (covered, no-alt) haplotype used as the NJ outgroup: best covered, then most reads, then earliest.
fn pick_outgroup(hm: &HaplotypeMatrix) -> Option<&Haplotype> {
    let covered = |h: &Haplotype| h.profile.iter().filter(|c| c.is_some()).count();
    hm.haplotypes
        .iter()
        .filter(|h| !h.profile.contains(&Some(1)) && covered(h) > 0)
        .min_by_key(|h| std::cmp::Reverse((covered(h), h.count)))
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

        // HF (else AF) of the first sample: one value per ALT, or a single shared value.
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
            if hf >= min_hf && hf < max_hf {
                let vid = format!("m.{pos}{ref_allele}>{alt_allele}");
                // The first record wins on a duplicated variant id.
                map.entry(vid).or_insert(hf);
            }
        }
    }
    Ok(map)
}

/// How a variant's heteroplasmic frequency is judged during matrix filtering.
pub enum HfFilter {
    /// Keep the variants in the (already HF-filtered) VCF map; tract-allele rows, which no VCF lists,
    /// are judged on their matrix frequency against `[min, max)`.
    FromVcf { map: HfMap, min: f64, max: f64 },
    /// Keep variants whose matrix frequency lies in `[min, max)`.
    FromMatrix { min: f64, max: f64 },
}

impl HfFilter {
    fn keeps(&self, vid: &str, freq: f64, tracts: &[crate::tracts::Tract]) -> bool {
        match self {
            HfFilter::FromVcf { map, min, max } => {
                if crate::tracts::is_tract_allele(vid, tracts) {
                    freq >= *min && freq < *max
                } else {
                    map.contains_key(vid)
                }
            }
            HfFilter::FromMatrix { min, max } => freq >= *min && freq < *max,
        }
    }
}

/// Alt and ref call counts of one matrix row, and the alt fraction over *covered* reads (0 if none).
pub(crate) fn row_hf(row: &[Option<u8>]) -> (usize, usize, f64) {
    let present = row.iter().filter(|c| **c == Some(1)).count();
    let absent = row.iter().filter(|c| **c == Some(0)).count();
    let covered = present + absent;
    let hf = if covered > 0 { present as f64 / covered as f64 } else { 0.0 };
    (present, absent, hf)
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

/// Apply the prevalence and heteroplasmy filters to a [`RawMatrix`]; `tracts` names the collapsed tract-allele rows.
pub fn filter_matrix(raw: &RawMatrix, hf: &HfFilter, tracts: &[crate::tracts::Tract]) -> BinaryMatrix {
    let mut variants: Vec<String> = Vec::new();
    let mut data: Vec<Vec<Option<u8>>> = Vec::new();

    for (vid, row) in raw.variants.iter().zip(raw.data.iter()) {
        let (present, absent, freq) = row_hf(row);
        if present < MIN_PRESENCE || absent < MIN_ABSENCE || !hf.keeps(vid, freq, tracts) {
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
            let seq = if missing > 0 {
                // A read uncovered everywhere in the tract loses nothing new; only partial coverage is counted.
                n_partial += usize::from(missing < members[t].len());
                None
            } else {
                let seq = apply_tract_ops(tract.ref_seq, &ops);
                n_conflict += usize::from(seq.is_none());
                seq
            };
            collapsed.push(seq);
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
    for (i, (vid, row)) in raw.variants.into_iter().zip(raw.data).enumerate() {
        match owner.get(&i) {
            Some(&t) if members[t][0].0 == i => {
                for (name, cells) in std::mem::take(&mut tract_rows[t]) {
                    variants.push(name);
                    data.push(cells);
                }
            }
            Some(_) => {}
            None => {
                variants.push(vid);
                data.push(row);
            }
        }
    }

    RawMatrix { variants, reads: raw.reads, data }
}

/// One edit a carried variant applies to a tract's reference sequence.
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

/// Near-fixed variants (not in `kept`) dropped by the HF ceiling, not for lack of support.
pub fn collect_fixed_homoplasmies(raw: &RawMatrix, max_hf: f64, kept: &[String]) -> Vec<String> {
    let kept: std::collections::HashSet<&String> = kept.iter().collect();
    raw.variants
        .iter()
        .zip(raw.data.iter())
        .filter(|(vid, row)| {
            let (present, _, freq) = row_hf(row);
            !kept.contains(vid) && present >= MIN_PRESENCE && freq >= max_hf
        })
        .map(|(vid, _)| vid.clone())
        .collect()
}


/// Collapse reads with identical binary profiles into unique haplotypes backed by at least `MIN_READS` reads.
/// Haplotypes are sorted by ascending mutation count (number of 1-bits), then lexicographically by
/// profile, so that the all-zero haplotype (if present) is always `H0000`.
pub fn deduplicate(matrix: &BinaryMatrix) -> HaplotypeMatrix {
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
    haplotypes.retain(|h| h.count >= MIN_READS);
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


/// SCITE (fp, fn) rates: explicit values verbatim (0 included), else the data type's preset.
pub fn resolve_error_rates(
    data_type: &str,
    fp_rate: Option<f64>,
    fn_rate: Option<f64>,
) -> (f64, f64) {
    // Every non-PacBio type uses the ONT preset.
    let default_fn = if data_type == "pacbio" { 0.05 } else { 0.02 };
    (fp_rate.unwrap_or(0.005), fn_rate.unwrap_or(default_fn))
}

/// HP-tract (fp, fn) rates: explicit values verbatim (0 included), else twice the global rates capped at 0.5.
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

/// Errors when no read profile is backed by at least `MIN_READS` reads (all were dropped).
fn ensure_haplotypes(hap_matrix: &HaplotypeMatrix) -> Result<()> {
    if hap_matrix.haplotypes.is_empty() {
        bail!(
            "No haplotypes remain: every distinct read profile is backed by fewer \
             than {MIN_READS} reads. This usually means the input is too shallow \
             or too noisy; widening --min-hf / --max-hf admits more variants and \
             so more distinguishable haplotypes."
        );
    }
    Ok(())
}

pub fn start(
    matrix_file: &str,
    vcf_file: Option<&str>,
    min_hf: f64,
    max_hf: f64,
    fp_rate: f64,
    fn_rate: f64,
    hp_rates: (f64, f64),
    tracts: &[crate::tracts::Tract],
    mcmc_iterations: usize,
    mcmc_chains: usize,
    debug: bool,
    output_prefix: &str,
) -> Result<()> {
    info!("Starting lineage analysis...");
    info!("SCITE fp={fp_rate}, fn={fn_rate}, hp_fp={}, hp_fn={}", hp_rates.0, hp_rates.1);

    // ── Step 1: decide how HF is judged, then filter the binary matrix ────────
    let hf = match vcf_file {
        Some(path) => {
            info!("[1/4] Parsing VCF: {}", path);
            HfFilter::FromVcf { map: parse_vcf(path, min_hf, max_hf)?, min: min_hf, max: max_hf }
        }
        None => {
            info!("[1/4] No VCF provided; HF recomputed from the matrix itself.");
            HfFilter::FromMatrix { min: min_hf, max: max_hf }
        }
    };

    info!("[1/4] Loading and filtering matrix: {}", matrix_file);
    let raw = read_raw_matrix(
        File::open(matrix_file)
            .with_context(|| format!("Cannot read matrix CSV: {matrix_file}"))?,
    )?;
    // An empty slice is the identity.
    if tracts.is_empty() {
        info!("[lineage] HP tract collapse: off (pass --collapse-hp to enable)");
    } else {
        info!("[lineage] HP tract collapse: on ({} tracts)", tracts.len());
    }
    let raw = collapse_tracts(raw, tracts);
    let binary = filter_matrix(&raw, &hf, tracts);
    if binary.variants.is_empty() {
        bail!(
            "No informative variants remain after filtering. \
             Try adjusting --min-hf / --max-hf thresholds."
        );
    }

    // ── Step 2: deduplicate reads into haplotypes ───────────────────────────────
    info!("[2/4] Deduplicating reads into haplotypes...");
    let hap_matrix = deduplicate(&binary);
    info!("[2/4] Found {} haplotypes across {} variants.", hap_matrix.haplotypes.len(), hap_matrix.variants.len());
    ensure_haplotypes(&hap_matrix)?;
    // Variant pairs violating the infinite-sites assumption hint at recurrent/back mutation or recombination.
    let violations = four_gamete_test(&hap_matrix);
    if !violations.is_empty() {
        info!(
            "[2/4] Four-gamete test: {} variant pair(s) violate the infinite-sites assumption.",
            violations.len()
        );
    }
    if debug {
        write_haplotype_map(&hap_matrix, &format!("{output_prefix}.raw_haplotype_map.tsv"))?;
    }

    // ── Step 3: SCITE mutation-tree search (NJ-informed MCMC) ─────────────────
    info!("[3/4] Running SCITE mutation-tree search...");
    let out = crate::scite::run_scite_pipeline(
        &binary,
        fp_rate,
        fn_rate,
        hp_rates,
        mcmc_iterations,
        mcmc_chains,
        debug,
        output_prefix,
    )?;

    // Near-fixed homoplasmies were removed by the HF ceiling before the tree was built; report them as all-present rows.
    let fixed_homoplasmies = collect_fixed_homoplasmies(&raw, max_hf, &binary.variants);
    if !fixed_homoplasmies.is_empty() {
        info!(
            "[3/4] Recovered {} near-fixed homoplasmic variant(s) excluded by the HF ceiling.",
            fixed_homoplasmies.len()
        );
    }
    let mut cleaned = out.cleaned;
    crate::rootblock::append_block_to_cleaned(&mut cleaned, &fixed_homoplasmies);
    crate::scite::write_cleaned_matrix(&cleaned, &format!("{output_prefix}.cleaned_matrix_all_variants.csv"))?;

    // ── Step 4: write final haplotype map (filtered by SCITE) ─────────────────
    let cleaned_hap_matrix = out.cleaned_haplotypes;
    info!(
        "[4/4] Cleaned reads deduplicated into {} haplotypes across {} variants.",
        cleaned_hap_matrix.haplotypes.len(),
        cleaned_hap_matrix.variants.len()
    );
    ensure_haplotypes(&cleaned_hap_matrix)?;
    write_haplotype_map(&cleaned_hap_matrix, &format!("{output_prefix}.cleaned_haplotype_map.tsv"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tracts::KNOWN_TRACTS;
    use std::io::Cursor;
    use std::path::Path;

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

    fn hap(id: &str, profile: &str, count: usize) -> Haplotype {
        Haplotype { id: id.to_string(), profile: cells(profile), reads: vec![], count }
    }

    /// Haplotypes over variants `v0, v1, …`.
    fn haps(haplotypes: Vec<Haplotype>) -> HaplotypeMatrix {
        let n = haplotypes[0].profile.len();
        HaplotypeMatrix { variants: (0..n).map(|i| format!("v{i}")).collect(), haplotypes }
    }

    /// An `HfFilter::FromVcf` holding exactly `names`.
    fn from_vcf(names: &[&str]) -> HfFilter {
        HfFilter::FromVcf { map: names.iter().map(|n| (n.to_string(), 0.5)).collect(), min: 0.0, max: 1.0 }
    }

    /// Parse `csv` and filter it without tract collapse.
    fn filter_from_csv(csv: &str, hf: &HfFilter) -> BinaryMatrix {
        filter_matrix(&read_raw_matrix(Cursor::new(csv)).unwrap(), hf, &[])
    }

    /// Write a one-sample chrM VCF with HF in FORMAT to the temp dir.
    fn write_vcf(name: &str, extra_header: &str, records: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::write(
            &path,
            format!(
                "##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n{extra_header}\
                 ##FORMAT=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmic fraction\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n{records}"
            ),
        )
        .unwrap();
        path
    }

    #[test]
    fn pick_outgroup_prefers_a_covered_ancestor_and_keeps_the_earlier_on_ties() {
        let pick = |hs: Vec<Haplotype>| pick_outgroup(&haps(hs)).map(|h| h.id.clone());
        assert_eq!(pick(vec![hap("H0", "..", 1), hap("H1", "00", 5), hap("H2", "11", 5)]), Some("H1".into()));
        // An ancestor without coverage is no outgroup.
        assert_eq!(pick(vec![hap("H0", "..", 1), hap("H1", "11", 5)]), None);
        assert_eq!(pick(vec![hap("H0", "00", 5), hap("H1", "00", 5), hap("H2", "11", 5)]), Some("H0".into()));
    }

    #[test]
    fn resolve_error_rates_uses_presets_unless_overridden() {
        assert_eq!(resolve_error_rates("pacbio", None, None), (0.005, 0.05));
        for ont in ["ont-r9", "ont-r10", "ont-denoised"] {
            assert_eq!(resolve_error_rates(ont, None, None), (0.005, 0.02));
        }
        // Each side overrides independently; 0 is a real choice, not "unset".
        assert_eq!(resolve_error_rates("pacbio", Some(0.02), None), (0.02, 0.05));
        assert_eq!(resolve_error_rates("pacbio", None, Some(0.2)), (0.005, 0.2));
        assert_eq!(resolve_error_rates("ont-r10", Some(0.0), Some(0.0)), (0.0, 0.0));
    }

    #[test]
    fn resolve_hp_rates_defaults_to_twice_the_global_rate_capped_at_half() {
        assert_eq!(resolve_hp_rates(0.005, 0.05, None, None), (0.01, 0.1));
        assert_eq!(resolve_hp_rates(0.4, 0.01, None, None), (0.5, 0.02));
        // An explicit flag is used verbatim, including a value above the cap and 0.
        assert_eq!(resolve_hp_rates(0.005, 0.05, Some(0.9), None), (0.9, 0.1));
        assert_eq!(resolve_hp_rates(0.005, 0.05, None, Some(0.0)), (0.01, 0.0));
    }

    #[test]
    fn filter_matrix_preserves_missing_values() {
        let csv = "variant,r1,r2,r3,r4\nm.100A>G,1,1,0,\nm.200C>T,0,1,1,1\n";
        let matrix = filter_from_csv(csv, &from_vcf(&["m.100A>G", "m.200C>T"]));
        assert_eq!(matrix.variants, vec!["m.100A>G", "m.200C>T"]);
        assert_eq!(matrix.reads, vec!["r1", "r2", "r3", "r4"]);
        assert_eq!(matrix.data, vec![cells("110."), cells("0111")]);
    }

    #[test]
    fn prevalence_filter_does_not_count_missing_as_absent() {
        // Two alt calls and no explicit ref call: fails MIN_ABSENCE.
        let matrix = filter_from_csv("variant,r1,r2,r3\nm.100A>G,1,1,\n", &from_vcf(&["m.100A>G"]));
        assert!(matrix.variants.is_empty());
    }

    #[test]
    fn matrix_frequency_filter_matches_the_vcf_filter() {
        // HF over covered reads: 0.50, 0.25, 0.90, 0.20 against [0.25, 0.85).
        // m.200C>T passes only because uncovered reads stay out of the denominator (2/8, not 2/10).
        let csv = "variant,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10\n\
                   m.100A>G,1,1,1,0,0,0,,,,\n\
                   m.200C>T,1,1,0,0,0,0,0,0,,\n\
                   m.300G>A,1,1,1,1,1,1,1,1,1,0\n\
                   m.400T>C,1,1,0,0,0,0,0,0,0,0\n";
        let with_vcf = filter_from_csv(csv, &from_vcf(&["m.100A>G", "m.200C>T"]));
        let without_vcf = filter_from_csv(csv, &HfFilter::FromMatrix { min: 0.25, max: 0.85 });
        assert_eq!(without_vcf.variants, vec!["m.100A>G", "m.200C>T"]);
        assert_eq!(with_vcf.variants, without_vcf.variants);
        assert_eq!(with_vcf.data, without_vcf.data);
    }

    #[test]
    fn a_vcf_matching_no_variants_keeps_nothing() {
        let kept = filter_from_csv("variant,r1,r2,r3\nm.100A>G,1,1,0\n", &from_vcf(&[]));
        assert!(kept.variants.is_empty(), "{:?}", kept.variants);
    }

    #[test]
    fn collect_fixed_homoplasmies_finds_supported_near_fixed_sites_not_already_kept() {
        // m.100: 4/4, recovered. m.200: 3/4, heteroplasmic. m.300: one alt call only. m.400: already kept.
        let csv = "variant,r1,r2,r3,r4\n\
                   m.100A>G,1,1,1,1\n\
                   m.200C>T,1,1,1,0\n\
                   m.300G>A,1,,,\n\
                   m.400T>C,1,1,1,1\n";
        let raw = read_raw_matrix(Cursor::new(csv)).unwrap();
        assert_eq!(collect_fixed_homoplasmies(&raw, 0.95, &["m.400T>C".to_string()]), vec!["m.100A>G"]);
    }

    #[test]
    fn deduplicate_treats_different_missing_patterns_as_different_haplotypes() {
        let matrix = BinaryMatrix {
            variants: vec!["m.100A>G".to_string()],
            reads: (1..=6).map(|r| format!("r{r}")).collect(),
            // Every read is "ref-or-missing" if you ignore the distinction; three per group clear `MIN_READS`.
            data: vec![cells("000...")],
        };
        assert_eq!(deduplicate(&matrix).haplotypes.len(), 2);
    }

    #[test]
    fn write_haplotype_map_adds_a_synthetic_root_only_without_an_ancestral_haplotype() {
        let path = std::env::temp_dir().join("himito_test_hapmap.tsv");
        for (profile, rows) in [("1", "H0000\t1\t2\tr1,r2\nHroot\t0\t0\t\n"), ("0", "H0000\t0\t2\tr1,r2\n")] {
            let hap_matrix = HaplotypeMatrix {
                variants: vec!["m.100A>G".to_string()],
                haplotypes: vec![Haplotype {
                    id: "H0000".to_string(),
                    profile: cells(profile),
                    reads: vec!["r1".to_string(), "r2".to_string()],
                    count: 2,
                }],
            };
            write_haplotype_map(&hap_matrix, path.to_str().unwrap()).unwrap();
            assert_eq!(
                std::fs::read_to_string(&path).unwrap(),
                format!("haplotype_id\tn_mutations\tn_reads\tread_name\n{rows}")
            );
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn hamming_distance_matrix_normalizes_by_jointly_covered_sites_only() {
        // Jointly covered sites are v0 and v1; they differ at v1.
        let dist = hamming_distance_matrix(&haps(vec![hap("H0", "10.", 1), hap("H1", "111", 1)]));
        assert_eq!(dist.data, vec![vec![0.0, 0.5], vec![0.5, 0.0]]);
    }

    #[test]
    fn hamming_distance_matrix_fills_uncovered_pairs_with_the_mean_distance() {
        // H0/H1 share no covered site; the other pairs are at 1.0 and 0.0.
        let dist = hamming_distance_matrix(&haps(vec![hap("H0", "1.", 1), hap("H1", ".1", 1), hap("H2", "01", 1)]));
        assert_eq!(dist.data, vec![vec![0.0, 0.5, 1.0], vec![0.5, 0.0, 0.0], vec![1.0, 0.0, 0.0]]);
        // With no finite distance at all, the fill is 0.
        let dist = hamming_distance_matrix(&haps(vec![hap("H0", "1.", 1), hap("H1", ".1", 1)]));
        assert_eq!(dist.data, vec![vec![0.0, 0.0], vec![0.0, 0.0]]);
    }

    #[test]
    fn four_gamete_test_needs_all_four_gametes_observed() {
        let gametes = |first: &str| haps(vec![hap("H0", first, 1), hap("H1", "01", 1), hap("H2", "10", 1), hap("H3", "11", 1)]);
        assert_eq!(four_gamete_test(&gametes("00")), vec![("v0".to_string(), "v1".to_string())]);
        // Missing calls do not fabricate the ref gamete.
        assert!(four_gamete_test(&gametes("..")).is_empty());
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
        let mut kids = tree.nodes[3].children.clone();
        kids.sort();
        assert_eq!(kids, vec![0, 1, 2]);
    }

    #[test]
    fn parse_vcf_splits_multiallelic_sites_and_keeps_first_duplicate() {
        let path = write_vcf(
            "himito_test_parse_vcf_multiallelic.vcf",
            "",
            "chrM\t100\t.\tA\tG,T\t.\tPASS\t.\tHF\t0.3,0.6\n\
             chrM\t200\t.\tC\tT\t.\tPASS\t.\tHF\t0.9\n\
             chrM\t100\t.\tA\tG\t.\tPASS\t.\tHF\t0.99\n",
        );
        let map = parse_vcf(path.to_str().unwrap(), 0.0, 1.0).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(map.get("m.100A>G").copied(), Some(0.3_f32 as f64));
        assert_eq!(map.get("m.100A>T").copied(), Some(0.6_f32 as f64));
        assert_eq!(map.get("m.200C>T").copied(), Some(0.9_f32 as f64));
        // The later duplicate record (HF 0.99) does not overwrite the first.
        assert_eq!(map.len(), 3);
    }

    /// Captures `warn!` lines with the thread that logged them, since tests run in parallel.
    struct WarnLog;

    static WARNINGS: std::sync::Mutex<Vec<(std::thread::ThreadId, String)>> = std::sync::Mutex::new(Vec::new());

    impl log::Log for WarnLog {
        fn enabled(&self, metadata: &log::Metadata) -> bool {
            metadata.level() <= log::Level::Warn
        }

        fn log(&self, record: &log::Record) {
            if self.enabled(record.metadata()) {
                WARNINGS.lock().unwrap().push((std::thread::current().id(), record.args().to_string()));
            }
        }

        fn flush(&self) {}
    }

    /// The left-alignment warnings `parse_vcf` logs on this thread for a VCF carrying `extra_header`.
    fn left_align_warnings(extra_header: &str) -> Vec<String> {
        static LOGGER: WarnLog = WarnLog;
        let _ = log::set_logger(&LOGGER);
        log::set_max_level(log::LevelFilter::Warn);
        let path = write_vcf("himito_test_parse_vcf_norm.vcf", extra_header, "chrM\t100\t.\tA\tG\t.\tPASS\t.\tHF\t0.3\n");
        let (me, start) = (std::thread::current().id(), WARNINGS.lock().unwrap().len());
        parse_vcf(path.to_str().unwrap(), 0.0, 1.0).unwrap();
        std::fs::remove_file(&path).ok();
        WARNINGS.lock().unwrap()[start..]
            .iter()
            .filter(|(id, line)| *id == me && line.contains("VCF was not left-aligned"))
            .map(|(_, line)| line.clone())
            .collect()
    }

    #[test]
    fn parse_vcf_warns_once_when_header_lacks_himito_normalization() {
        let hits = left_align_warnings("");
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert!(hits[0].contains("without `--no-left-align`"), "{}", hits[0]);
        assert!(left_align_warnings("##HimitoNormalization=left-aligned\n").is_empty());
    }

    /// A matrix CSV in which every cell of `rows` stands for three identical reads, so each haplotype clears `MIN_READS`.
    fn tripled_csv(rows: &[(&str, &str)]) -> String {
        let mut csv = String::from("variant");
        for r in 1..=3 * rows[0].1.len() {
            csv += &format!(",r{r}");
        }
        for (vid, row) in rows {
            csv += &format!("\n{vid}");
            for c in row.chars() {
                csv += &format!(",{c},{c},{c}");
            }
        }
        csv + "\n"
    }

    /// Run `start` without a VCF on `csv` in a fresh temp dir `dir`; returns the output prefix.
    fn run_small_lineage(dir: &str, csv: &str, debug: bool, tracts: &[crate::tracts::Tract]) -> String {
        let dir = std::env::temp_dir().join(dir);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let matrix = dir.join("matrix.csv");
        std::fs::write(&matrix, csv).unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        start(matrix.to_str().unwrap(), None, 0.01, 0.95, 0.01, 0.1, (0.01, 0.1), tracts, 20, 1, debug, &prefix).unwrap();
        prefix
    }

    #[test]
    fn debug_gates_the_raw_outputs() {
        let csv = tripled_csv(&[("m.100A>G", "1100"), ("m.200C>T", "1010")]);
        for debug in [false, true] {
            let prefix = run_small_lineage(&format!("himito_lineage_debug_{debug}"), &csv, debug, &[]);
            for suffix in [".raw_haplotype_map.tsv", ".raw_variant_cooccurrence.tsv"] {
                assert_eq!(Path::new(&format!("{prefix}{suffix}")).exists(), debug, "{suffix} with debug={debug}");
            }
            assert!(!Path::new(&format!("{prefix}.cleaned_matrix.csv")).exists());
        }
    }

    #[test]
    fn collapse_hp_replaces_hvs1_rows_with_a_tract_allele() {
        let csv = tripled_csv(&[("m.16183A>AC", "1000"), ("m.16189T>TC", "0100"), ("m.100A>G", "1100")]);
        let off: &[crate::tracts::Tract] = &[];
        for (tracts, on) in [(off, false), (KNOWN_TRACTS, true)] {
            let prefix = run_small_lineage(&format!("himito_lineage_collapse_{on}"), &csv, false, tracts);
            let out = std::fs::read_to_string(format!("{prefix}.cleaned_matrix_all_variants.csv")).unwrap();
            assert_eq!(out.contains("m.16183A>AC"), !on, "{out}");
            assert_eq!(out.contains("m.16189T>TC"), !on, "{out}");
            assert_eq!(out.contains("m.16183ACCCCCTCCCC>"), on, "{out}");
        }
    }

    /// rCRS 16183-16193: `A` + `C5` + `T` + `C4`.
    const HVS1: &str = "ACCCCCTCCCC";
    /// The tract with T16189C only: a single run of 10 C's.
    const C10: &str = "ACCCCCCCCCC";
    /// T16189C plus one extra C, spelled either side of the interrupting T.
    const C11: &str = "ACCCCCCCCCCC";

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
        // r1 carries the extra C before the interrupting T, r2 after it.
        let out = collapse_tracts(
            raw(&["m.16183A>AC", "m.16189T>TC", "m.16189T>C"], &["1000", "0100", "1111"]),
            KNOWN_TRACTS,
        );
        assert_eq!(out.variants, vec![hvs1(C10), hvs1(C11)]);
        assert_eq!(out.data, vec![cells("0011"), cells("1100")]);
        assert_eq!(out.reads, vec!["r1", "r2", "r3", "r4"]);
    }

    #[test]
    fn collapse_tracts_empty_slice_returns_the_input_unchanged() {
        let input = raw(&["m.16183A>AC", "m.16189T>TC", "m.16189T>C"], &["1000", "0100", "1111"]);
        assert_eq!(collapse_tracts(input.clone(), &[]), input);
    }

    #[test]
    fn collapse_tracts_gives_reference_reads_zero_and_emits_no_reference_row() {
        // r3 carries nothing inside the tract; r1 keeps the T and gains a C; r2 is C10.
        let out = collapse_tracts(raw(&["m.16183A>AC", "m.16189T>C"], &["100", "010"]), KNOWN_TRACTS);
        assert_eq!(out.variants, vec![hvs1(C10), hvs1("ACCCCCCTCCCC")]);
        assert_eq!(out.data, vec![cells("010"), cells("100")]);
    }

    #[test]
    fn collapse_tracts_marks_partially_covered_reads_missing() {
        // r2 is uncovered on one constituent row only, r3 on both.
        let out = collapse_tracts(raw(&["m.16183A>AC", "m.16189T>C"], &["1..", "01."]), KNOWN_TRACTS);
        assert_eq!(out.variants, vec![hvs1("ACCCCCCTCCCC")]);
        assert_eq!(out.data, vec![cells("1..")]);
    }

    #[test]
    fn collapse_tracts_sets_conflicting_reads_missing() {
        for (variants, allele) in [
            // Two insertions after the same offset.
            (["m.16183A>ACC", "m.16183A>AC"], "ACCCCCCTCCCC"),
            // A deletion whose deleted base is also substituted.
            (["m.16188CT>C", "m.16189T>C"], C10),
            // Two different substitutions at one offset.
            (["m.16189T>A", "m.16189T>C"], C10),
        ] {
            let out = collapse_tracts(raw(&variants, &["10", "11"]), KNOWN_TRACTS);
            assert_eq!(out.variants, vec![hvs1(allele)], "{variants:?}");
            assert_eq!(out.data, vec![cells(".1")], "{variants:?}");
        }
    }

    #[test]
    fn collapse_tracts_leaves_variants_straddling_the_boundary_alone() {
        let out = collapse_tracts(raw(&["m.16189T>C", "m.16193CA>C"], &["11", "01"]), KNOWN_TRACTS);
        assert_eq!(out.variants, vec![hvs1(C10), "m.16193CA>C".to_string()]);
        assert_eq!(out.data, vec![cells("11"), cells("01")]);
    }

    #[test]
    fn collapse_tracts_leaves_non_indel_non_snv_names_alone() {
        // An MNV is not collapsible.
        let out = collapse_tracts(raw(&["m.16188CT>GA", "m.16189T>C"], &["10", "11"]), KNOWN_TRACTS);
        assert_eq!(out.variants, vec!["m.16188CT>GA".to_string(), hvs1(C10)]);
        assert_eq!(out.data, vec![cells("10"), cells("11")]);
    }

    #[test]
    fn collapse_tracts_skips_a_tract_whose_reference_does_not_match() {
        // rCRS has C at 16185, so `m.16185G>A` says this is not rCRS; the whole tract passes through untouched.
        let input = raw(&["m.16185G>A", "m.16189T>C"], &["10", "11"]);
        assert_eq!(collapse_tracts(input.clone(), KNOWN_TRACTS), input);
    }

    #[test]
    fn collapse_tracts_keeps_other_rows_in_place() {
        let out = collapse_tracts(
            raw(
                &["m.100A>G", "m.16189T>C", "m.200C>T", "m.16183A>AC", "m.300G>A"],
                &["10", "11", "01", "10", "11"],
            ),
            KNOWN_TRACTS,
        );
        assert_eq!(
            out.variants,
            vec!["m.100A>G".to_string(), hvs1(C10), hvs1(C11), "m.200C>T".to_string(), "m.300G>A".to_string()]
        );
        assert_eq!(out.data, vec![cells("10"), cells("01"), cells("10"), cells("01"), cells("11")]);
    }

    #[test]
    fn tract_alleles_are_judged_on_matrix_frequency_with_or_without_a_vcf() {
        // Two tract alleles at 2/4 = 0.50, plus m.3243A>G, which the VCF lists.
        let csv = "variant,r1,r2,r3,r4\n\
                   m.16183A>AC,1,0,0,0\n\
                   m.16189T>TC,0,1,0,0\n\
                   m.16189T>C,1,1,1,1\n\
                   m.3243A>G,1,1,0,0\n";
        let raw = collapse_tracts(read_raw_matrix(Cursor::new(csv)).unwrap(), KNOWN_TRACTS);
        let kept = |hf: HfFilter| filter_matrix(&raw, &hf, KNOWN_TRACTS).variants;
        let vcf = |max| HfFilter::FromVcf { map: [("m.3243A>G".to_string(), 0.5)].into(), min: 0.2, max };
        let all = vec![hvs1(C10), hvs1(C11), "m.3243A>G".to_string()];
        assert_eq!(kept(vcf(0.95)), all);
        assert_eq!(kept(HfFilter::FromMatrix { min: 0.2, max: 0.95 }), all);
        // A ceiling below 0.50 drops both tract alleles even on the VCF path.
        assert_eq!(kept(vcf(0.4)), vec!["m.3243A>G"]);
        assert!(kept(HfFilter::FromMatrix { min: 0.2, max: 0.4 }).is_empty());
    }
}
