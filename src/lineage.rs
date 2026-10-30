use std::collections::HashMap;
use anyhow::{bail,Context, Result};
use rust_htslib::bcf::{Read, Reader}; // Read trait must be in scope for .records()
use std::io::{BufWriter, Write};
use std::fs::File;
use std::path::Path;
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
#[derive(Clone)]
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
///
/// The map is a name whitelist (see [`HfFilter::FromVcf`]). `[min_hf, max_hf)` is applied here to the
/// VCF HF and then again to the matrix HF in [`filter_matrix`], on purpose: a site the VCF itself calls
/// fixed is not a lineage candidate, and dropping this window would widen the whitelist.
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
    /// Keep the variants the VCF names, each still judged on its own matrix frequency against
    /// `[min, max)`. The VCF is a whitelist only: its HF reflects pileup depth, while everything
    /// downstream (root block, SCITE, co-occurrence) counts matrix rows, and the two disagree
    /// badly at homopolymer tracts.
    FromVcf { map: HfMap, min: f64, max: f64 },
    /// Keep variants whose matrix frequency lies in `[min, max)`.
    FromMatrix { min: f64, max: f64 },
}

impl HfFilter {
    fn keeps(&self, vid: &str, freq: f64, tracts: &[crate::tracts::Tract]) -> bool {
        match self {
            HfFilter::FromVcf { map, min, max } => {
                // Tract-allele rows are collapsed spellings no VCF lists; they are judged on matrix HF alone.
                (map.contains_key(vid) || crate::tracts::is_tract_allele(vid, tracts))
                    && freq >= *min
                    && freq < *max
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

/// Why the matrix row `vid` is not a candidate for the tree: `(outcome, detail)`, or `None` if it passes every
/// gate. The one place the gates live, so [`filter_matrix`] and [`variant_dispositions`] cannot disagree.
fn screen_variant(
    vid: &str,
    (n_present, n_absent, matrix_hf): (usize, usize, f64),
    hf: &HfFilter,
    tracts: &[crate::tracts::Tract],
) -> Option<(&'static str, String)> {
    if n_present < MIN_PRESENCE {
        return Some(("dropped_min_presence", format!("{n_present} present < {MIN_PRESENCE}")));
    }
    // The HF window is judged before absence: a variant with no ref call at all has matrix HF 1.0, and "above the
    // window" is the useful reason for it. `collect_fixed_homoplasmies` recovers such variants by matrix HF alone,
    // without consulting the whitelist, so with `-v` a fixed site is usually labelled `dropped_not_whitelisted`
    // here (its VCF HF is already at or above `--max-hf`) and only a whitelisted site labelled `dropped_hf_window`.
    if !hf.keeps(vid, matrix_hf, tracts) {
        return Some(match hf {
            HfFilter::FromVcf { map, min, max }
                if !map.contains_key(vid) && !crate::tracts::is_tract_allele(vid, tracts) =>
            {
                (
                    "dropped_not_whitelisted",
                    format!("not in the VCF whitelist (not a PASS site, or its VCF HF is outside [{min}, {max}))"),
                )
            }
            HfFilter::FromVcf { min, max, .. } | HfFilter::FromMatrix { min, max } => {
                ("dropped_hf_window", format!("matrix HF {matrix_hf:.4} outside [{min}, {max})"))
            }
        });
    }
    if n_absent < MIN_ABSENCE {
        return Some(("dropped_min_absence", format!("{n_absent} absent < {MIN_ABSENCE}")));
    }
    None
}

/// Apply the prevalence and heteroplasmy filters to a [`RawMatrix`]; `tracts` names the collapsed tract-allele rows.
pub fn filter_matrix(raw: &RawMatrix, hf: &HfFilter, tracts: &[crate::tracts::Tract]) -> BinaryMatrix {
    let mut variants: Vec<String> = Vec::new();
    let mut data: Vec<Vec<Option<u8>>> = Vec::new();

    for (vid, row) in raw.variants.iter().zip(raw.data.iter()) {
        if screen_variant(vid, row_hf(row), hf, tracts).is_some() {
            continue;
        }
        variants.push(vid.clone());
        data.push(row.clone());
    }

    BinaryMatrix { variants, reads: raw.reads.clone(), data }
}

/// Where one candidate variant ended up. One row per in-window PASS site (the whitelist
/// `parse_vcf` builds) and per matrix row, so an ONT run and a PacBio run can be compared with a single
/// `diff`. A non-PASS site, a PASS site with neither an HF nor an AF value, or a PASS site whose VCF HF is outside the
/// window gets no row at all when it has no matrix row.
///
/// `outcome` is one of:
/// - `informative`: placed in the mutation tree;
/// - `pending`: written only in the provisional table a failed run leaves behind; the variant survived filtering
///   but the tree search had not run, so its final outcome is unknown. A table with any `pending` row is from a
///   failed run and must not be compared against another platform's;
/// - `root_block`: near-fixed, removed from the search and reported unordered at the root;
/// - `tract_annotated`: a homopolymer-tract allele reported in `tract_annotation.tsv` but kept out of the search;
/// - `dropped_hf_window`: matrix HF outside `[--min-hf, --max-hf)`;
/// - `dropped_min_presence`: fewer than [`MIN_PRESENCE`] reads carry it;
/// - `dropped_min_absence`: fewer than [`MIN_ABSENCE`] reads call it absent;
/// - `dropped_not_whitelisted`: the VCF does not list it as a PASS site inside the HF window (see `detail`);
/// - `no_matrix_row`: a PASS site inside the HF window that the matrix has no row for.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantDisposition {
    pub variant: String,
    /// `None` when the site is not in the in-window whitelist (a collapsed tract allele, or a no-VCF run).
    pub vcf_hf: Option<f64>,
    pub matrix_hf: f64,
    pub n_present: usize,
    pub n_absent: usize,
    pub outcome: &'static str,
    pub detail: String,
}

/// Classify every candidate by the same gates as [`filter_matrix`]. `blocked` names the variants the root
/// block took; the pre-search pass has no root block yet and uses [`provisional_dispositions`] instead.
pub fn variant_dispositions(
    raw: &RawMatrix,
    hf: &HfFilter,
    tracts: &[crate::tracts::Tract],
    blocked: &[String],
) -> Vec<VariantDisposition> {
    let vcf_hf = |vid: &str| match hf {
        HfFilter::FromVcf { map, .. } => map.get(vid).copied(),
        HfFilter::FromMatrix { .. } => None,
    };
    let mut rows: Vec<VariantDisposition> = raw
        .variants
        .iter()
        .zip(raw.data.iter())
        .map(|(vid, row)| {
            let stats @ (n_present, n_absent, matrix_hf) = row_hf(row);
            let (outcome, detail) = screen_variant(vid, stats, hf, tracts).unwrap_or_else(|| {
                if blocked.contains(vid) {
                    ("root_block", "near-fixed: reported unordered at the root".to_string())
                } else {
                    ("informative", String::new())
                }
            });
            VariantDisposition { variant: vid.clone(), vcf_hf: vcf_hf(vid), matrix_hf, n_present, n_absent, outcome, detail }
        })
        .collect();
    // Whitelisted names the matrix never genotyped: silent otherwise, and the reason an ONT/PacBio
    // comparison can look like a lineage difference when it is a calling difference.
    if let HfFilter::FromVcf { map, .. } = hf {
        for (vid, hfv) in map {
            if !raw.variants.iter().any(|v| v == vid) {
                rows.push(VariantDisposition {
                    variant: vid.clone(),
                    vcf_hf: Some(*hfv),
                    matrix_hf: 0.0,
                    n_present: 0,
                    n_absent: 0,
                    outcome: "no_matrix_row",
                    detail: "PASS in the VCF but absent from matrix.csv".to_string(),
                });
            }
        }
    }
    rows.sort_by(|a, b| a.variant.cmp(&b.variant));
    rows
}

/// The table written before the tree search: every variant that survived filtering reads `pending`, so a run that
/// fails before the final rewrite leaves a file no reader can mistake for a successful run's. Rows a filter already
/// dropped are final and are left as [`variant_dispositions`] reports them.
pub fn provisional_dispositions(
    raw: &RawMatrix,
    hf: &HfFilter,
    tracts: &[crate::tracts::Tract],
) -> Vec<VariantDisposition> {
    let mut rows = variant_dispositions(raw, hf, tracts, &[]);
    for r in rows.iter_mut().filter(|r| r.outcome == "informative") {
        r.outcome = "pending";
        r.detail = "the tree search had not run".to_string();
    }
    rows
}

/// Relabel the candidates the search withheld as tract annotations (`annotated` comes from the search itself).
pub fn mark_tract_annotated(rows: &mut [VariantDisposition], annotated: &[String]) {
    for r in rows.iter_mut().filter(|r| r.outcome == "informative" && annotated.contains(&r.variant)) {
        r.outcome = "tract_annotated";
        r.detail = "reported in tract_annotation.tsv, kept out of the tree search".to_string();
    }
}

pub fn write_variant_disposition(rows: &[VariantDisposition], path: &str) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(w, "variant\tvcf_hf\tmatrix_hf\tn_present\tn_absent\toutcome\tdetail")?;
    for r in rows {
        let vcf = r.vcf_hf.map_or("NA".to_string(), |h| format!("{h:.6}"));
        let detail = if r.detail.is_empty() { "NA" } else { r.detail.as_str() };
        writeln!(
            w,
            "{}\t{vcf}\t{:.6}\t{}\t{}\t{}\t{detail}",
            r.variant, r.matrix_hf, r.n_present, r.n_absent, r.outcome
        )?;
    }
    Ok(())
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


/// Floor on the per-sample SNV false-negative rate.
pub const SNV_FN_FLOOR: f64 = 0.005;
/// Fewer near-fixed SNV sites (or calls) than this and the data-type preset is used instead.
const DROPOUT_MIN_SITES: usize = 5;
const DROPOUT_MIN_CALLS: usize = 500;
/// SNVs this close to a known tract are genotyped inside its window and are left out of the estimate.
const DROPOUT_TRACT_MARGIN: u64 = 5;
/// A site whose absences are this unlikely under the pooled rate is real near-fixed heteroplasmy, not dropout.
const DROPOUT_OUTLIER_P: f64 = 0.001;
/// Near-fixed gate of the dropout estimate: `min(--max-hf, DROPOUT_MAX_GATE)`, so `--max-hf 1.0` still measures dropout.
const DROPOUT_MAX_GATE: f64 = 0.95;

/// The reference pseudo-read `call` adds to `matrix.csv` (named after the contig, all-ref); never a real molecule.
pub fn is_reference_pseudo_read(name: &str) -> bool {
    matches!(name, "chrM" | "MT" | "chrMT" | "NC_012920.1")
}

/// Per-sample SNV dropout measured on near-fixed SNV rows.
#[derive(Debug, Clone, PartialEq)]
pub struct DropoutEstimate {
    /// `max(raw, SNV_FN_FLOOR)`.
    pub rate: f64,
    /// Pooled absent / called over the kept sites.
    pub raw: f64,
    pub n_sites: usize,
    pub n_calls: usize,
    pub n_absent: usize,
    /// Sites dropped as real near-fixed heteroplasmy.
    pub outliers: Vec<String>,
}

/// The 1-based position of a single-base SNV name, else `None`.
fn snv_position(variant: &str) -> Option<u64> {
    let (pos, r, a) = parse_variant_name(variant)?;
    let base = |x: &str| x.len() == 1 && matches!(x.as_bytes()[0], b'A' | b'C' | b'G' | b'T' | b'N');
    (base(r) && base(a)).then_some(pos)
}

/// SNV false-negative rate from near-fixed SNV rows: pooled ref calls / calls over SNVs with matrix HF >=
/// `min(max_hf, DROPOUT_MAX_GATE)`,
/// ignoring the reference pseudo-read and SNVs within `DROPOUT_TRACT_MARGIN` bp of a known tract. Sites whose
/// absences are significant against the pooled rate are dropped once as real heteroplasmy. `None` when too few remain.
pub fn estimate_snv_dropout(raw: &RawMatrix, max_hf: f64) -> Option<DropoutEstimate> {
    let gate = max_hf.min(DROPOUT_MAX_GATE);
    let keep: Vec<usize> = raw
        .reads
        .iter()
        .enumerate()
        .filter(|(_, name)| !is_reference_pseudo_read(name))
        .map(|(i, _)| i)
        .collect();
    let near_tract = |pos: u64| {
        crate::tracts::KNOWN_TRACTS
            .iter()
            .any(|t| pos + DROPOUT_TRACT_MARGIN >= t.start && pos <= t.end + DROPOUT_TRACT_MARGIN)
    };
    let mut sites: Vec<(&String, usize, usize)> = Vec::new();
    for (vid, row) in raw.variants.iter().zip(raw.data.iter()) {
        let Some(pos) = snv_position(vid) else { continue };
        if near_tract(pos) {
            continue;
        }
        let present = keep.iter().filter(|&&r| row[r] == Some(1)).count();
        let absent = keep.iter().filter(|&&r| row[r] == Some(0)).count();
        if present >= MIN_PRESENCE && present as f64 / (present + absent) as f64 >= gate {
            sites.push((vid, present, absent));
        }
    }
    let pooled = |s: &[(&String, usize, usize)]| {
        let calls: usize = s.iter().map(|&(_, p, a)| p + a).sum();
        let absent: usize = s.iter().map(|&(_, _, a)| a).sum();
        (calls, absent, if calls > 0 { absent as f64 / calls as f64 } else { 0.0 })
    };
    let (_, _, rate0) = pooled(&sites);
    let (outlier_sites, sites): (Vec<_>, Vec<_>) = sites
        .into_iter()
        .partition(|&(_, p, a)| crate::scite::binomial_upper_tail(p + a, a, rate0) <= DROPOUT_OUTLIER_P);
    let (n_calls, n_absent, raw_rate) = pooled(&sites);
    if sites.len() < DROPOUT_MIN_SITES || n_calls < DROPOUT_MIN_CALLS {
        return None;
    }
    Some(DropoutEstimate {
        rate: raw_rate.max(SNV_FN_FLOOR),
        raw: raw_rate,
        n_sites: sites.len(),
        n_calls,
        n_absent,
        outliers: outlier_sites.iter().map(|&(v, _, _)| v.clone()).collect(),
    })
}

/// What to do with homopolymer-tract alleles (D310, HVS1-C) in the tree search.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum TractAlleleMode {
    /// Tract alleles are placed as multi-allelic sites (`Site`); an allele whose run is longer than the
    /// confusion table calibrates is placed too, flagged `beyond_calibration` in `tract_annotation.tsv`,
    /// so its placement is judged on its margin rather than deleted.
    #[default]
    Auto,
    /// Every tract allele is placed, with the alleles at one tract position scored as one categorical site whose
    /// states may recur on several branches.
    Site,
    /// Every tract allele is placed as its own binary row (the model before multi-allelic sites).
    Search,
    /// Every tract allele is reported (tract_annotation.tsv) but never placed.
    Annotate,
}

impl std::str::FromStr for TractAlleleMode {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, String> {
        match s {
            "auto" => Ok(Self::Auto),
            "site" => Ok(Self::Site),
            "search" => Ok(Self::Search),
            "annotate" => Ok(Self::Annotate),
            other => Err(format!("unknown --tract-alleles {other:?} (auto, site, search, annotate)")),
        }
    }
}

/// How `lineage` should pick its error rates; every `None` falls back to a default.
#[derive(Debug, Clone, Default)]
pub struct RateOptions {
    pub data_type: String,
    pub fp_rate: Option<f64>,
    pub fn_rate: Option<f64>,
    pub hp_fp_rate: Option<f64>,
    pub hp_fn_rate: Option<f64>,
    /// `hp_confusion.tsv` to take tract-allele rates from; `None` looks beside the matrix, then the embedded default.
    pub hp_error: Option<String>,
    pub tract_alleles: TractAlleleMode,
}

/// One row of `<prefix>.error_model.tsv`.
#[derive(Debug, Clone, PartialEq)]
pub struct ErrorModelRow {
    pub scope: String,
    pub variant: String,
    pub param: &'static str,
    pub value: f64,
    pub source: String,
    pub detail: String,
}

pub fn write_error_model(rows: &[ErrorModelRow], path: &str) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(w, "scope\tvariant\tparam\tvalue\tsource\tdetail")?;
    for r in rows {
        writeln!(w, "{}\t{}\t{}\t{:.6}\t{}\t{}", r.scope, r.variant, r.param, r.value, r.source, r.detail)?;
    }
    Ok(())
}

/// What the `hp *` rows of `error_model.tsv` apply to; D310/HVS1-C alleles get their own rows.
const HP_DEFAULT_DETAIL: &str = "default for CA514 alleles and tract alleles without a homopolymer run";

/// Resolved SNV and default HP rates, with the rows that document them.
#[derive(Debug, Clone)]
pub struct ResolvedRates {
    pub fp_rate: f64,
    pub fn_rate: f64,
    pub hp_rates: (f64, f64),
    pub rows: Vec<ErrorModelRow>,
}

/// SNV fp: `--fp-rate` or the preset. SNV fn: `--fn-rate`, else the per-sample estimate, else the preset.
/// Default HP rates are twice `--fp-rate`/`--fn-rate` when given, else twice the preset, never twice the estimate.
pub fn resolve_rates(opts: &RateOptions, raw: &RawMatrix, max_hf: f64) -> ResolvedRates {
    let (preset_fp, preset_fn) = resolve_error_rates(&opts.data_type, None, None);
    let preset = format!("preset:{}", opts.data_type);
    let estimate = estimate_snv_dropout(raw, max_hf);
    let detail = match &estimate {
        Some(e) => format!(
            "sites={} calls={} absent={} raw={:.6} outliers={}",
            e.n_sites,
            e.n_calls,
            e.n_absent,
            e.raw,
            if e.outliers.is_empty() { "-".to_string() } else { e.outliers.join(",") }
        ),
        None => format!("fewer than {DROPOUT_MIN_SITES} near-fixed SNV sites or {DROPOUT_MIN_CALLS} calls"),
    };
    let (fp_rate, fp_source) = match opts.fp_rate {
        Some(v) => (v, "cli".to_string()),
        None => (preset_fp, preset.clone()),
    };
    let (fn_rate, fn_source) = match (opts.fn_rate, &estimate) {
        (Some(v), _) => (v, "cli".to_string()),
        (None, Some(e)) => (e.rate, if e.raw < SNV_FN_FLOOR { "estimated(floor)" } else { "estimated" }.to_string()),
        (None, None) => (preset_fn, preset.clone()),
    };
    let base_fp = opts.fp_rate.unwrap_or(preset_fp);
    let base_fn = opts.fn_rate.unwrap_or(preset_fn);
    let hp_rates = resolve_hp_rates(base_fp, base_fn, opts.hp_fp_rate, opts.hp_fn_rate);
    let hp_source = |given: Option<f64>, base_given: Option<f64>| match (given, base_given) {
        (Some(_), _) => "cli".to_string(),
        (None, Some(_)) => "2x-cli".to_string(),
        (None, None) => "2x-preset".to_string(),
    };
    let row = |scope: &str, param: &'static str, value: f64, source: String, detail: String| ErrorModelRow {
        scope: scope.to_string(),
        variant: "*".to_string(),
        param,
        value,
        source,
        detail,
    };
    let rows = vec![
        row("snv", "fp", fp_rate, fp_source, String::new()),
        row("snv", "fn", fn_rate, fn_source, detail),
        row("hp", "fp", hp_rates.0, hp_source(opts.hp_fp_rate, opts.fp_rate), HP_DEFAULT_DETAIL.to_string()),
        row("hp", "fn", hp_rates.1, hp_source(opts.hp_fn_rate, opts.fn_rate), HP_DEFAULT_DETAIL.to_string()),
    ];
    ResolvedRates { fp_rate, fn_rate, hp_rates, rows }
}

// ─── Tract-allele rates from the homopolymer confusion table ─────────────────────

/// Rows with fewer observations than this are not used.
const HP_TABLE_MIN_N: u64 = 200;
/// Tract-allele rates are clamped into [floor, cap].
const HP_RATE_FLOOR: f64 = 0.005;
const HP_RATE_CAP: f64 = 0.5;
/// Background variants within this many bp of an allele shape the run it sits in.
const RUN_BACKGROUND_BP: i64 = 10;

// Embedded data lives under src/data/ because resource/ and test_data/ are git-ignored.
const RCRS_FASTA: &str = include_str!("data/NC_012920.1.fasta");
const DEFAULT_HP_TABLE_ONT: &str = include_str!("data/hp_error/HG00097_ont-r10.hp_confusion.tsv");
const DEFAULT_HP_TABLE_PACBIO: &str = include_str!("data/hp_error/HG00097_pacbio.hp_confusion.tsv");

/// The embedded rCRS (NC_012920.1), uppercased.
pub fn rcrs() -> &'static [u8] {
    static RCRS: std::sync::OnceLock<Vec<u8>> = std::sync::OnceLock::new();
    RCRS.get_or_init(|| {
        RCRS_FASTA
            .lines()
            .filter(|l| !l.starts_with('>'))
            .flat_map(|l| l.trim().bytes())
            .map(|b| b.to_ascii_uppercase())
            .collect()
    })
}

/// Base class of a homopolymer run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RunClass {
    AT,
    CG,
}

impl RunClass {
    fn of(base: u8) -> Self {
        if matches!(base, b'C' | b'G') { RunClass::CG } else { RunClass::AT }
    }
}

impl std::fmt::Display for RunClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            RunClass::AT => "A/T",
            RunClass::CG => "C/G",
        })
    }
}

/// `m.<pos><REF>><ALT>` with plain ACGTN alleles as (0-based pos, REF, ALT).
fn plain_edit(name: &str) -> Option<(i64, &str, &str)> {
    let (pos, r, a) = parse_variant_name(name)?;
    let acgtn = |x: &str| x.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'));
    (pos >= 1 && acgtn(r) && acgtn(a)).then_some((pos as i64 - 1, r, a))
}

/// The homopolymer run an allele sits in (see [`implied_run_at`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImpliedRun {
    pub class: RunClass,
    pub len: usize,
    /// The run's extent on the reference, 1-based inclusive. Inserted bases count at the last REF base of their edit,
    /// so an insertion anchored just left of a run is placed on that run.
    pub ref_start: u64,
    pub ref_end: u64,
}

impl ImpliedRun {
    /// The run's extent as a [`crate::scite::VariantSpan`] (1-based, half-open).
    pub fn span(&self) -> crate::scite::VariantSpan {
        crate::scite::VariantSpan { start: self.ref_start, end: self.ref_end + 1 }
    }
}

/// [`implied_run_at`] without the run's extent.
#[cfg(test)]
pub fn implied_run(name: &str, refseq: &[u8], background: &[String]) -> Option<(RunClass, usize)> {
    implied_run_at(name, refseq, background).map(|r| (r.class, r.len))
}

/// The homopolymer run allele `name` sits in once `background` (other variant names) is applied to `refseq`.
/// Background edits within `RUN_BACKGROUND_BP` that do not overlap the allele are applied left to right with it
/// (an edit overlapping an applied one is skipped). The run is taken at the edit itself: REF and ALT are trimmed
/// of their shared prefix and suffix, and the run is the longest single-base run intersecting what the ALT
/// changed (leftmost on ties). For a pure deletion it is the run beside the junction whose base is the deleted
/// base (else the longer of the two runs beside the junction). `None` if the allele's REF does not match `refseq`.
/// Mirrors `review_2026-10-20/reference/implied_run.py`.
pub fn implied_run_at(name: &str, refseq: &[u8], background: &[String]) -> Option<ImpliedRun> {
    let matches_ref = |pos: i64, r: &str| {
        let (s, e) = (pos as usize, pos as usize + r.len());
        e <= refseq.len() && &refseq[s..e] == r.as_bytes()
    };
    // Edits drop the bases REF and ALT share at the front: a VCF indel's anchor base is not part of what it changes,
    // so a background allele may change that base (m.16189T>TC on 16189T>C). An insertion has an empty REF.
    fn core(pos: i64, r: &str, a: &str) -> (i64, String, String) {
        let pre = r.bytes().zip(a.bytes()).take_while(|(x, y)| x == y).count();
        (pos + pre as i64, r[pre..].to_string(), a[pre..].to_string())
    }
    // Half-open REF spans; an empty span is an insertion point, which clashes with an insertion at the same point
    // or an edit strictly containing it.
    fn clash(a: (i64, i64), b: (i64, i64)) -> bool {
        match (a.0 == a.1, b.0 == b.1) {
            (false, false) => a.0 < b.1 && b.0 < a.1,
            (true, true) => a.0 == b.0,
            (true, false) => b.0 < a.0 && a.0 < b.1,
            (false, true) => a.0 < b.0 && b.0 < a.1,
        }
    }
    let (vpos, vref, valt) = plain_edit(name)?;
    if !matches_ref(vpos, vref) {
        return None;
    }
    let (cpos, cref, calt) = core(vpos, vref, valt);
    let vspan = (cpos, cpos + cref.len() as i64);
    // (pos, ref, alt, is_v)
    let mut edits: Vec<(i64, String, String, bool)> = vec![(cpos, cref.clone(), calt.clone(), true)];
    for b in background {
        let Some((pos, r, a)) = plain_edit(b) else { continue };
        if !matches_ref(pos, r) {
            continue;
        }
        let (bpos, bref, balt) = core(pos, r, a);
        let bspan = (bpos, bpos + bref.len() as i64);
        if clash(vspan, bspan) {
            continue; // changes what the allele changes
        }
        if bspan.1 <= vspan.0 - RUN_BACKGROUND_BP || bspan.0 >= vspan.1 + RUN_BACKGROUND_BP {
            continue; // too far to matter
        }
        edits.push((bpos, bref, balt, false));
    }
    // By position; an insertion goes in before an edit starting at the same base; the allele before background.
    edits.sort_by_key(|e| (e.0, !e.1.is_empty(), !e.3));
    let mut out: Vec<u8> = Vec::with_capacity(refseq.len() + 16);
    // 0-based reference position of each base of `out`; inserted bases count at the base before them.
    let mut src: Vec<i64> = Vec::with_capacity(refseq.len() + 16);
    let mut applied: Vec<(i64, i64)> = Vec::new();
    let mut cur: i64 = 0;
    let mut alt_iv = None;
    for (pos, r, a, is_v) in edits {
        let span = (pos, pos + r.len() as i64);
        if pos < cur || applied.iter().any(|&t| clash(t, span)) {
            if is_v {
                return None;
            }
            continue;
        }
        applied.push(span);
        out.extend_from_slice(&refseq[cur as usize..pos as usize]);
        src.extend(cur..pos);
        let start = out.len();
        out.extend_from_slice(a.as_bytes());
        let last = if r.is_empty() { (pos - 1).max(0) } else { pos + r.len() as i64 - 1 };
        src.extend((0..a.len() as i64).map(|j| (pos + j).min(last)));
        if is_v {
            alt_iv = Some((start, out.len()));
        }
        cur = pos + r.len() as i64;
    }
    out.extend_from_slice(&refseq[cur as usize..]);
    src.extend(cur..refseq.len() as i64);
    let (alt_start, _) = alt_iv?;
    let (vref, valt) = (cref.as_str(), calt.as_str());
    // What the allele changed: REF/ALT minus their shared prefix and suffix.
    let (r, a) = (vref.as_bytes(), valt.as_bytes());
    let max_trim = r.len().min(a.len());
    let pre = r.iter().zip(a).take_while(|(x, y)| x == y).count().min(max_trim);
    let suf = r.iter().rev().zip(a.iter().rev()).take_while(|(x, y)| x == y).count().min(max_trim - pre);
    let (lo, hi) = (alt_start + pre, alt_start + a.len() - suf);
    // Maximal runs of `out` as (start, end, base).
    let mut runs: Vec<(usize, usize, u8)> = Vec::new();
    let mut i = 0;
    while i < out.len() {
        let mut j = i + 1;
        while j < out.len() && out[j] == out[i] {
            j += 1;
        }
        runs.push((i, j, out[i]));
        i = j;
    }
    let longest = |cands: Vec<(usize, usize, u8)>| {
        cands
            .into_iter()
            .fold(None::<(usize, usize, u8)>, |best, c| match best {
                Some(b) if b.1 - b.0 >= c.1 - c.0 => Some(b),
                _ => Some(c),
            })
            .map(|(s, e, b)| ImpliedRun {
                class: RunClass::of(b),
                len: e - s,
                ref_start: src[s] as u64 + 1,
                ref_end: src[e - 1] as u64 + 1,
            })
    };
    if lo < hi {
        return longest(runs.iter().copied().filter(|&(s, e, _)| s < hi && lo < e).collect());
    }
    // Pure deletion: the junction sits between out[lo - 1] and out[lo].
    let deleted = &r[pre..r.len() - suf];
    let beside: Vec<(usize, usize, u8)> = runs
        .iter()
        .copied()
        .filter(|&(s, e, _)| (lo > 0 && s <= lo - 1 && lo - 1 < e) || (s <= lo && lo < e))
        .collect();
    let same_base: Vec<(usize, usize, u8)> = match deleted.first() {
        Some(&d) if deleted.iter().all(|&x| x == d) => beside.iter().copied().filter(|&(_, _, b)| b == d).collect(),
        _ => Vec::new(),
    };
    longest(if same_base.is_empty() { beside } else { same_base })
}

/// Whether `name` is a homopolymer-tract allele: its REF overlaps D310 or HVS1-C, or it is an indel whose implied
/// run on the sample background reaches into a tract's body (past its anchor A). The second case catches length
/// alleles that VCF left-normalisation anchors outside the tract, e.g. `m.16181A>AC` once 16182C and 16183C join
/// the HVS1 C run, without taking in the A-run indels next to the tracts.
pub fn is_homopolymer_tract_allele(name: &str, span: &crate::scite::VariantSpan, run: Option<&ImpliedRun>) -> bool {
    if crate::tracts::overlaps_homopolymer_tract(span) {
        return true;
    }
    let indel = plain_edit(name).is_some_and(|(_, r, a)| r.len() != a.len());
    indel && run.is_some_and(|r| crate::tracts::overlaps_homopolymer_tract_body(&r.span()))
}

/// Rates for one tract allele from the confusion table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TractRates {
    pub fp: f64,
    pub fn_: f64,
    /// Longest run length of the allele's class with at least `HP_TABLE_MIN_N` observations.
    pub calibrated_max: u32,
    /// The allele's run is longer than `calibrated_max`: the rates are lower bounds.
    pub beyond_calibration: bool,
}

/// Per (class, reference run length): `[n, n(obs = L), n(obs = L - 1), n(obs = L + 1)]`, and the count of every
/// observed length change; strands pooled.
#[derive(Debug, Clone, PartialEq)]
pub struct HpRateTable {
    rows: std::collections::BTreeMap<(RunClass, u32), [u64; 4]>,
    hist: std::collections::BTreeMap<(RunClass, u32), std::collections::BTreeMap<i64, u64>>,
}

impl HpRateTable {
    /// `None` when no (class, length) has `HP_TABLE_MIN_N` observations.
    pub fn from_confusion(rows: &[crate::hperror::ConfusionRow]) -> Option<Self> {
        let mut agg: std::collections::BTreeMap<(RunClass, u32), [u64; 4]> = std::collections::BTreeMap::new();
        let mut hist: std::collections::BTreeMap<(RunClass, u32), std::collections::BTreeMap<i64, u64>> =
            std::collections::BTreeMap::new();
        for &(base, ref_len, obs_len, _strand, n) in rows {
            let key = (RunClass::of(base as u8), ref_len);
            let e = agg.entry(key).or_default();
            e[0] += n;
            let delta = obs_len as i64 - ref_len as i64;
            match delta {
                0 => e[1] += n,
                -1 => e[2] += n,
                1 => e[3] += n,
                _ => {}
            }
            *hist.entry(key).or_default().entry(delta).or_default() += n;
        }
        agg.retain(|_, v| v[0] >= HP_TABLE_MIN_N);
        hist.retain(|k, _| agg.contains_key(k));
        (!agg.is_empty()).then_some(HpRateTable { rows: agg, hist })
    }

    /// P(observed − true length = Δ) for a `class` run of length `len`, from the length-change counts at the nearest
    /// calibrated length at or below it (the shortest if none; a class with no rows pools both classes).
    pub fn kernel(&self, class: RunClass, len: u32) -> Vec<(i64, f64)> {
        let mut by_len: std::collections::BTreeMap<u32, std::collections::BTreeMap<i64, u64>> =
            self.hist.iter().filter(|((c, _), _)| *c == class).map(|((_, l), h)| (*l, h.clone())).collect();
        if by_len.is_empty() {
            for ((_, l), h) in &self.hist {
                let e = by_len.entry(*l).or_default();
                for (d, n) in h {
                    *e.entry(*d).or_default() += n;
                }
            }
        }
        let h = by_len.range(..=len).next_back().or_else(|| by_len.iter().next()).map(|(_, h)| h).unwrap();
        let total: u64 = h.values().sum();
        h.iter().map(|(&d, &n)| (d, n as f64 / total as f64)).collect()
    }

    pub fn from_text(text: &str) -> Result<Option<Self>> {
        Ok(Self::from_confusion(&crate::hperror::read_confusion(text)?))
    }

    /// Rows usable for `class`; a class with none pools both classes per length.
    fn class_rows(&self, class: RunClass) -> std::collections::BTreeMap<u32, [u64; 4]> {
        let own: std::collections::BTreeMap<u32, [u64; 4]> =
            self.rows.iter().filter(|((c, _), _)| *c == class).map(|((_, l), v)| (*l, *v)).collect();
        if !own.is_empty() {
            return own;
        }
        let mut pooled: std::collections::BTreeMap<u32, [u64; 4]> = std::collections::BTreeMap::new();
        for ((_, l), v) in &self.rows {
            let e = pooled.entry(*l).or_default();
            for k in 0..4 {
                e[k] += v[k];
            }
        }
        pooled
    }

    /// Calibrated length range `[min, max]` of `class`.
    pub fn calibrated_range(&self, class: RunClass) -> (u32, u32) {
        let rows = self.class_rows(class);
        (*rows.keys().next().unwrap(), *rows.keys().last().unwrap())
    }

    /// fn = 1 - P0(L); fp = max(P(-1 | L+1), P(+1 | L-1)). Each length is looked up at the nearest calibrated
    /// length at or below it (the shortest calibrated length if there is none), so gaps in the table are bridged.
    pub fn allele_rates(&self, class: RunClass, len: u32) -> TractRates {
        let rows = self.class_rows(class);
        let lmax = *rows.keys().last().unwrap();
        let frac = |l: i64, k: usize| {
            let l = l.max(0) as u32;
            let v = rows.range(..=l).next_back().or_else(|| rows.iter().next()).map(|(_, v)| *v).unwrap();
            v[k] as f64 / v[0] as f64
        };
        let l = len as i64;
        let clamp = |x: f64| x.clamp(HP_RATE_FLOOR, HP_RATE_CAP);
        TractRates {
            fp: clamp(frac(l + 1, 2).max(frac(l - 1, 3))),
            fn_: clamp(1.0 - frac(l, 1)),
            calibrated_max: lmax,
            beyond_calibration: len > lmax,
        }
    }
}

/// The confusion table and where it came from: `--hp-error`, `<X>.hp_confusion.tsv` beside `<X>.matrix.csv`,
/// else the embedded HG00097 table of the data type.
pub fn load_hp_table(opts: &RateOptions, matrix_file: &str) -> Result<(HpRateTable, String)> {
    if let Some(path) = &opts.hp_error {
        let text = std::fs::read_to_string(path).with_context(|| format!("Cannot read --hp-error {path}"))?;
        let table = HpRateTable::from_text(&text)
            .with_context(|| format!("Cannot parse --hp-error {path}"))?
            .with_context(|| format!("--hp-error {path}: no run length with {HP_TABLE_MIN_N} observations"))?;
        return Ok((table, format!("table:{path}")));
    }
    if let Some(stem) = matrix_file.strip_suffix(".matrix.csv") {
        let sibling = format!("{stem}.hp_confusion.tsv");
        if Path::new(&sibling).exists() {
            match std::fs::read_to_string(&sibling).map_err(anyhow::Error::from).and_then(|t| HpRateTable::from_text(&t)) {
                Ok(Some(table)) => {
                    info!("[rates] Tract-allele rates from {sibling}");
                    return Ok((table, format!("table:{sibling}")));
                }
                Ok(None) => warn!("{sibling}: no run length with {HP_TABLE_MIN_N} observations; using the default table"),
                Err(e) => warn!("{sibling}: {e:#}; using the default table"),
            }
        }
    }
    let (text, name) = if opts.data_type == "pacbio" {
        (DEFAULT_HP_TABLE_PACBIO, "HG00097_pacbio")
    } else {
        (DEFAULT_HP_TABLE_ONT, "HG00097_ont-r10")
    };
    warn!(
        "No hp_confusion.tsv for this sample: tract-allele rates come from the embedded {name} table (one sample). \
         `call --bam` writes <X>.hp_confusion.tsv beside the matrix, or pass that file with --hp-error."
    );
    let table = HpRateTable::from_text(text)?.context("embedded HP table has no usable rows")?;
    Ok((table, format!("default:{name}")))
}

/// Implied runs are read on rows with matrix HF above this: the sample's majority alleles, the same rule `call`
/// uses (on VCF HF) for the haplotypes it genotypes length alleles against.
const RUN_BACKGROUND_MIN_HF: f64 = 0.5;

/// Majority rows (matrix HF > [`RUN_BACKGROUND_MIN_HF`] over real reads, present >= `MIN_PRESENCE`): the background
/// an allele's implied run is read on.
fn majority_names(raw: &RawMatrix) -> Vec<String> {
    let keep: Vec<usize> =
        raw.reads.iter().enumerate().filter(|(_, n)| !is_reference_pseudo_read(n)).map(|(i, _)| i).collect();
    raw.variants
        .iter()
        .zip(raw.data.iter())
        .filter(|(_, row)| {
            let present = keep.iter().filter(|&&r| row[r] == Some(1)).count();
            let absent = keep.iter().filter(|&&r| row[r] == Some(0)).count();
            present >= MIN_PRESENCE && present as f64 / (present + absent) as f64 > RUN_BACKGROUND_MIN_HF
        })
        .map(|(v, _)| v.clone())
        .collect()
}

/// The run-length change of a pure homopolymer indel: REF and ALT minus their shared prefix and suffix leave one
/// side empty and the other one base repeated. `None` for substitutions and compound edits.
fn homopolymer_length_change(name: &str) -> Option<i64> {
    let (_, r, a) = plain_edit(name)?;
    let (r, a) = (r.as_bytes(), a.as_bytes());
    let pre = r.iter().zip(a).take_while(|(x, y)| x == y).count();
    let (r, a) = (&r[pre..], &a[pre..]);
    let suf = r.iter().rev().zip(a.iter().rev()).take_while(|(x, y)| x == y).count();
    let (r, a) = (&r[..r.len() - suf], &a[..a.len() - suf]);
    let run = |x: &[u8]| !x.is_empty() && x.iter().all(|&b| b == x[0]);
    match (r.is_empty(), a.is_empty()) {
        (true, false) if run(a) => Some(a.len() as i64),
        (false, true) if run(r) => Some(-(r.len() as i64)),
        _ => None,
    }
}

/// The site's root length: the most common length called at `anchor` over the real reads of `raw`, among the
/// states that are not `members` (REF at `ref_len`, or a length allele the filters dropped). Ties go to REF, then the
/// shorter length; with no such call it is `ref_len`.
fn site_root_length(raw: &RawMatrix, anchor: u64, ref_len: i64, members: &[&str]) -> i64 {
    let spans = crate::scite::variant_spans(&raw.variants);
    let rows: Vec<usize> = (0..raw.variants.len()).filter(|&i| spans[i].start == anchor).collect();
    let mut counts: std::collections::BTreeMap<i64, usize> = std::collections::BTreeMap::new();
    for (r, read) in raw.reads.iter().enumerate() {
        if is_reference_pseudo_read(read) {
            continue;
        }
        let ones: Vec<usize> = rows.iter().copied().filter(|&i| raw.data[i][r] == Some(1)).collect();
        let len = match ones[..] {
            [] if rows.iter().all(|&i| raw.data[i][r] == Some(0)) => Some(ref_len),
            [i] if !members.contains(&raw.variants[i].as_str()) => {
                homopolymer_length_change(&raw.variants[i]).map(|d| ref_len + d)
            }
            _ => None,
        };
        if let Some(len) = len {
            *counts.entry(len).or_default() += 1;
        }
    }
    let best = counts.values().copied().max().unwrap_or(0);
    if best == 0 || counts.get(&ref_len) == Some(&best) {
        return ref_len;
    }
    counts.iter().find(|(_, &n)| n == best).map_or(ref_len, |(&len, _)| len)
}

/// Multi-allelic sites from the placed homopolymer-tract alleles `(variant, class, run length)`. Pure length
/// alleles anchored at the same base (their REF spans start together, so they are mutually exclusive) form one site
/// when they change a homopolymer and agree on the reference run length (run length less length change; the
/// members' majority, ties to the longer); substitutions, compound edits, disagreeing alleles and a second allele
/// of the same length stay binary rows. The root length comes from `raw` ([`site_root_length`]); the kernel from
/// `table`, for every length among the states.
fn tract_site_specs(
    variants: &[String],
    spans: &[crate::scite::VariantSpan],
    candidates: &[(usize, RunClass, usize)],
    table: &HpRateTable,
    raw: &RawMatrix,
) -> Vec<crate::scite::SiteSpec> {
    let mut by_anchor: std::collections::BTreeMap<u64, Vec<(usize, RunClass, i64, i64)>> =
        std::collections::BTreeMap::new();
    for &(v, class, len) in candidates {
        if let Some(d) = homopolymer_length_change(&variants[v]) {
            by_anchor.entry(spans[v].start).or_default().push((v, class, len as i64, d));
        }
    }
    by_anchor
        .into_iter()
        .filter_map(|(anchor, group)| {
            // The reference run length the members imply (run length less change): the majority, ties to the longer.
            let mut votes: std::collections::BTreeMap<i64, usize> = std::collections::BTreeMap::new();
            for m in &group {
                *votes.entry(m.2 - m.3).or_default() += 1;
            }
            let top = votes.values().copied().max()?;
            let ref_len = votes.iter().rev().find(|(_, &n)| n == top).map(|(&l, _)| l)?;
            let mut seen_len: Vec<i64> = Vec::new();
            let mut agree = Vec::new();
            for m in group {
                // A member changes a homopolymer (run >= 2; the implied run is the run its own bases join) on the
                // shared reference run, at a length no other member has.
                if m.2 >= 2 && m.2 - m.3 == ref_len && ref_len >= 1 && !seen_len.contains(&m.2) {
                    seen_len.push(m.2);
                    agree.push(m);
                } else {
                    warn!(
                        "[rates] {} (run {}, change {:+}) does not fit the length site at this position \
                         (reference run {ref_len}): kept as its own row",
                        variants[m.0],
                        m.2,
                        m.3
                    );
                }
            }
            if agree.is_empty() {
                return None;
            }
            let class = agree[0].1;
            let names: Vec<&str> = agree.iter().map(|m| variants[m.0].as_str()).collect();
            let root_len = site_root_length(raw, anchor, ref_len, &names);
            let member_len: Vec<i64> = agree.iter().map(|m| m.2).collect();
            let mut lens: Vec<i64> = member_len.clone();
            lens.push(root_len);
            lens.sort_unstable();
            lens.dedup();
            let kernel = lens.iter().map(|&l| (l, table.kernel(class, l.max(1) as u32))).collect();
            Some(crate::scite::SiteSpec { members: agree.iter().map(|m| m.0).collect(), member_len, root_len, kernel })
        })
        .collect()
}

/// Per-variant rates (and the report on each tract allele) for `binary`, with their `error_model.tsv` rows.
/// D310/HVS1-C alleles ([`is_homopolymer_tract_allele`]) take table rates for their implied run (CLI `--hp-*-rate`
/// wins), other tract alleles the default HP rates, everything else the SNV rates.
pub fn build_variant_model(
    binary: &BinaryMatrix,
    raw: &RawMatrix,
    resolved: &ResolvedRates,
    opts: &RateOptions,
    table: Option<(&HpRateTable, &str)>,
) -> (crate::scite::VariantModel, Vec<ErrorModelRow>) {
    let background = majority_names(raw);
    let spans = crate::scite::variant_spans(&binary.variants);
    let mut rates = Vec::with_capacity(binary.variants.len());
    let mut tract = Vec::with_capacity(binary.variants.len());
    let mut rows = Vec::new();
    // (variant, run class, run length) of the placed homopolymer-tract alleles.
    let mut site_candidates: Vec<(usize, RunClass, usize)> = Vec::new();
    for (name, span) in binary.variants.iter().zip(spans.iter()) {
        let at = implied_run_at(name, rcrs(), &background);
        let run = is_homopolymer_tract_allele(name, span, at.as_ref()).then_some(at).flatten().map(|r| (r.class, r.len));
        let Some((class, len)) = run else {
            let in_tract = crate::tracts::overlaps_tract(span);
            rates.push(if in_tract { resolved.hp_rates } else { (resolved.fp_rate, resolved.fn_rate) });
            tract.push(None);
            continue;
        };
        // Without a table both HP rates come from the command line (the table is not loaded then).
        let tr = table.map(|(t, _)| t.allele_rates(class, len as u32));
        let fp = opts.hp_fp_rate.or(tr.map(|t| t.fp)).unwrap_or(resolved.hp_rates.0);
        let fn_ = opts.hp_fn_rate.or(tr.map(|t| t.fn_)).unwrap_or(resolved.hp_rates.1);
        let detail = match tr {
            Some(t) => format!(
                "run={class}:{len} calibrated_max={} beyond_calibration={}",
                t.calibrated_max,
                if t.beyond_calibration { "yes" } else { "no" }
            ),
            None => format!("run={class}:{len}"),
        };
        let table_source = table.map_or("cli", |(_, src)| src);
        let source = |cli: Option<f64>| if cli.is_some() { "cli".to_string() } else { table_source.to_string() };
        for (param, value, src) in [("fp", fp, source(opts.hp_fp_rate)), ("fn", fn_, source(opts.hp_fn_rate))] {
            rows.push(ErrorModelRow {
                scope: "hp".to_string(),
                variant: name.clone(),
                param,
                value,
                source: src,
                detail: detail.clone(),
            });
        }
        rates.push((fp, fn_));
        let (placed, reason) = match opts.tract_alleles {
            TractAlleleMode::Annotate => (false, "--tract-alleles annotate".to_string()),
            TractAlleleMode::Auto | TractAlleleMode::Site | TractAlleleMode::Search => match tr {
                Some(t) if t.beyond_calibration => {
                    (true, format!("beyond_calibration (run {len} > calibrated {})", t.calibrated_max))
                }
                Some(t) => (true, format!("calibrated (run {len} <= calibrated {})", t.calibrated_max)),
                None => (true, "no confusion table".to_string()),
            },
        };
        if placed {
            site_candidates.push((rates.len() - 1, class, len));
        }
        tract.push(Some(crate::scite::TractReport { run_class: class.to_string(), run_length: len, placed, reason }));
    }
    // Site kernels come from the table; any explicit `--hp-*-rate` keeps the alleles as binary rows (so it is used).
    let any_cli = opts.hp_fp_rate.is_some() || opts.hp_fn_rate.is_some();
    let sites = match (opts.tract_alleles, table) {
        (TractAlleleMode::Auto | TractAlleleMode::Site, Some((t, src))) if !any_cli => {
            let sites = tract_site_specs(&binary.variants, &spans, &site_candidates, t, raw);
            for site in &sites {
                let names: Vec<&str> = site.members.iter().map(|&v| binary.variants[v].as_str()).collect();
                let lens: Vec<String> = site.member_len.iter().map(|l| l.to_string()).collect();
                rows.push(ErrorModelRow {
                    scope: "site".to_string(),
                    variant: names.join(","),
                    param: "root_len",
                    value: site.root_len as f64,
                    source: src.to_string(),
                    detail: format!(
                        "member_len={} (one categorical character; root = the commonest other length; states may recur)",
                        lens.join(",")
                    ),
                });
            }
            sites
        }
        _ => Vec::new(),
    };
    for (class, (t, src)) in [RunClass::AT, RunClass::CG].into_iter().flat_map(|c| table.map(|t| (c, t))) {
        let (lo, hi) = t.calibrated_range(class);
        rows.push(ErrorModelRow {
            scope: "hp_table".to_string(),
            variant: "*".to_string(),
            param: if class == RunClass::AT { "calibrated_max_AT" } else { "calibrated_max_CG" },
            value: hi as f64,
            source: src.to_string(),
            detail: format!("calibrated run lengths {lo}..{hi} (n >= {HP_TABLE_MIN_N})"),
        });
    }
    (crate::scite::VariantModel { rates, tract, sites }, rows)
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
    rate_options: &RateOptions,
    tracts: &[crate::tracts::Tract],
    mcmc_iterations: usize,
    mcmc_chains: usize,
    debug: bool,
    output_prefix: &str,
) -> Result<()> {
    info!("Starting lineage analysis...");

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
    // Rates come from the raw matrix: the SNV fn is measured on near-fixed SNV rows the filters below drop.
    let resolved = resolve_rates(rate_options, &raw, max_hf);
    let (fp_rate, fn_rate, hp_rates) = (resolved.fp_rate, resolved.fn_rate, resolved.hp_rates);
    info!("SCITE fp={fp_rate}, fn={fn_rate}, default hp_fp={}, hp_fn={}", hp_rates.0, hp_rates.1);
    // The table is only needed when an HP rate is left to it.
    let hp_table = if rate_options.hp_fp_rate.is_some() && rate_options.hp_fn_rate.is_some() {
        None
    } else {
        Some(load_hp_table(rate_options, matrix_file)?)
    };
    // An empty slice is the identity.
    if tracts.is_empty() {
        info!("[lineage] HP tract collapse: off (pass --collapse-hp to enable)");
    } else {
        info!("[lineage] HP tract collapse: on ({} tracts)", tracts.len());
    }
    let raw = collapse_tracts(raw, tracts);
    let binary = filter_matrix(&raw, &hf, tracts);
    // With --debug, written before anything can fail, then rewritten once the search has named its root block
    // and annotations.
    let disposition_path = format!("{output_prefix}.variant_disposition.tsv");
    if debug {
        write_variant_disposition(&provisional_dispositions(&raw, &hf, tracts), &disposition_path)?;
    }
    if binary.variants.is_empty() {
        let where_ = if debug {
            format!("{disposition_path} lists every candidate and why it was dropped; ")
        } else {
            "Re-run with --debug to write variant_disposition.tsv; ".to_string()
        };
        bail!("No variant survives filtering. {where_}try relaxing --min-hf / --max-hf.");
    }
    let (variant_model, allele_rows) =
        build_variant_model(&binary, &raw, &resolved, rate_options, hp_table.as_ref().map(|(t, s)| (t, s.as_str())));
    let error_model: Vec<ErrorModelRow> = resolved.rows.iter().cloned().chain(allele_rows).collect();
    for r in &error_model {
        info!("[rates] {} {} {} = {:.6} ({}) {}", r.scope, r.variant, r.param, r.value, r.source, r.detail);
    }
    write_error_model(&error_model, &format!("{output_prefix}.error_model.tsv"))?;

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
        Some(&variant_model),
        mcmc_iterations,
        mcmc_chains,
        debug,
        output_prefix,
    )?;

    if debug {
        let mut dispositions = variant_dispositions(&raw, &hf, tracts, &out.block_names);
        mark_tract_annotated(&mut dispositions, &out.annotated_names);
        write_variant_disposition(&dispositions, &disposition_path)?;
    }

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

    /// A `matrix.csv` with `n_reads` reads plus a trailing `chrM` pseudo-read column (always ref);
    /// each row is `(name, n_alt, n_ref)`, alt calls first, the rest uncovered.
    fn csv_with_pseudo_read(n_reads: usize, rows: &[(&str, usize, usize)]) -> String {
        let mut out = String::from("variant");
        for r in 0..n_reads {
            out.push_str(&format!(",r{r}"));
        }
        out.push_str(",chrM\n");
        for &(name, alt, rf) in rows {
            out.push_str(name);
            for r in 0..n_reads {
                out.push(',');
                if r < alt {
                    out.push('1');
                } else if r < alt + rf {
                    out.push('0');
                }
            }
            out.push_str(",0\n");
        }
        out
    }

    const FIXED_SNVS: [&str; 6] = ["m.1000A>G", "m.2000A>G", "m.3000A>G", "m.4000A>G", "m.5000A>G", "m.6000A>G"];

    fn fixed_rows(alt: usize, rf: usize) -> Vec<(&'static str, usize, usize)> {
        FIXED_SNVS.iter().map(|&n| (n, alt, rf)).collect()
    }

    fn raw_from(csv: &str) -> RawMatrix {
        read_raw_matrix(Cursor::new(csv)).unwrap()
    }

    #[test]
    fn estimate_snv_dropout_pools_snv_rows_and_ignores_the_reference_pseudo_read() {
        let est = estimate_snv_dropout(&raw_from(&csv_with_pseudo_read(200, &fixed_rows(198, 2))), 0.95).unwrap();
        assert_eq!((est.n_sites, est.n_calls, est.n_absent), (6, 1200, 12));
        assert!((est.raw - 0.01).abs() < 1e-12 && (est.rate - 0.01).abs() < 1e-12, "{est:?}");
        assert!(est.outliers.is_empty());
    }

    #[test]
    fn estimate_snv_dropout_drops_a_near_fixed_heteroplasmy_as_an_outlier() {
        // Six sites at 1/200 absent, one at 10/200 (HF 0.95, a real ~5% heteroplasmy).
        let mut rows = fixed_rows(199, 1);
        rows.push(("m.7000C>T", 190, 10));
        let est = estimate_snv_dropout(&raw_from(&csv_with_pseudo_read(200, &rows)), 0.95).unwrap();
        assert_eq!(est.outliers, vec!["m.7000C>T"]);
        assert_eq!((est.n_sites, est.n_absent), (6, 6));
        assert!((est.raw - 0.005).abs() < 1e-12, "{est:?}");
    }

    #[test]
    fn estimate_snv_dropout_skips_indels_and_tract_adjacent_snvs() {
        let mut rows = fixed_rows(198, 2);
        // m.16182A>C sits next to HVS1-C; m.8270C>CA is an indel. Counted, either would raise the rate.
        rows.push(("m.16182A>C", 190, 9));
        rows.push(("m.8270C>CA", 190, 9));
        let est = estimate_snv_dropout(&raw_from(&csv_with_pseudo_read(200, &rows)), 0.95).unwrap();
        assert_eq!(est.n_sites, 6);
        assert!((est.raw - 0.01).abs() < 1e-12, "{est:?}");
    }

    #[test]
    fn estimate_snv_dropout_floors_at_0_005_and_needs_five_sites() {
        let est = estimate_snv_dropout(&raw_from(&csv_with_pseudo_read(200, &fixed_rows(200, 0))), 0.95).unwrap();
        assert_eq!((est.raw, est.rate), (0.0, SNV_FN_FLOOR));
        let four: Vec<_> = fixed_rows(198, 2).into_iter().take(4).collect();
        assert!(estimate_snv_dropout(&raw_from(&csv_with_pseudo_read(200, &four)), 0.95).is_none());
        // Five sites but under 500 calls.
        let thin: Vec<_> = fixed_rows(80, 1).into_iter().take(5).collect();
        assert!(estimate_snv_dropout(&raw_from(&csv_with_pseudo_read(100, &thin)), 0.95).is_none());
    }

    #[test]
    fn resolve_rates_prefers_cli_then_estimate_then_preset_and_never_doubles_the_estimate() {
        let raw = raw_from(&csv_with_pseudo_read(200, &fixed_rows(198, 2)));
        let opts = |fp: Option<f64>, fnr: Option<f64>| RateOptions {
            data_type: "ont-r10".to_string(),
            fp_rate: fp,
            fn_rate: fnr,
            ..Default::default()
        };
        let r = resolve_rates(&opts(None, None), &raw, 0.95);
        assert_eq!((r.fp_rate, r.fn_rate), (0.005, 0.01));
        assert_eq!(r.hp_rates, (0.01, 0.04), "HP defaults stay at twice the preset, not twice the estimate");
        let src = |r: &ResolvedRates, scope: &str, param: &str| {
            r.rows.iter().find(|x| x.scope == scope && x.param == param).unwrap().source.clone()
        };
        assert_eq!(src(&r, "snv", "fn"), "estimated");
        assert_eq!(src(&r, "snv", "fp"), "preset:ont-r10");
        assert_eq!(src(&r, "hp", "fn"), "2x-preset");

        let r = resolve_rates(&opts(None, Some(0.03)), &raw, 0.95);
        assert_eq!((r.fn_rate, r.hp_rates.1), (0.03, 0.06));
        assert_eq!((src(&r, "snv", "fn"), src(&r, "hp", "fn")), ("cli".to_string(), "2x-cli".to_string()));

        let few = raw_from(&csv_with_pseudo_read(200, &fixed_rows(198, 2)[..3]));
        let r = resolve_rates(&opts(None, None), &few, 0.95);
        assert_eq!((r.fn_rate, src(&r, "snv", "fn")), (0.02, "preset:ont-r10".to_string()));
    }

    #[test]
    fn start_writes_the_error_model_with_its_sources() {
        let mut rows = fixed_rows(198, 2);
        rows.push(("m.9000A>G", 100, 100));
        rows.push(("m.9100C>T", 50, 150));
        let dir = std::env::temp_dir();
        let matrix = dir.join("himito_test_error_model.matrix.csv");
        std::fs::write(&matrix, csv_with_pseudo_read(200, &rows)).unwrap();
        for (fn_rate, want) in [(None, "estimated"), (Some(0.02), "cli")] {
            let prefix = dir.join("himito_test_error_model").to_str().unwrap().to_string();
            let opts = RateOptions { data_type: "pacbio".to_string(), fn_rate, ..Default::default() };
            start(matrix.to_str().unwrap(), None, 0.01, 0.95, &opts, &[], 200, 1, false, &prefix).unwrap();
            let table = std::fs::read_to_string(format!("{prefix}.error_model.tsv")).unwrap();
            let fn_row = table.lines().find(|l| l.starts_with("snv\t*\tfn\t")).unwrap();
            assert_eq!(fn_row.split('\t').nth(4), Some(want), "{table}");
        }
    }

    #[test]
    fn implied_run_matches_the_reference_cases() {
        let bg = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        let cases: [(&str, Vec<String>, (RunClass, usize)); 10] = [
            ("m.302A>ACC", bg(&[]), (RunClass::CG, 9)),
            ("m.302A>ACCC", bg(&[]), (RunClass::CG, 10)),
            ("m.302A>AC", bg(&[]), (RunClass::CG, 8)),
            ("m.310T>TC", bg(&[]), (RunClass::CG, 6)),
            ("m.16189T>C", bg(&[]), (RunClass::CG, 10)),
            ("m.16183A>C", bg(&[]), (RunClass::CG, 6)),
            ("m.16189T>C", bg(&["m.16182A>C", "m.16183A>C"]), (RunClass::CG, 12)),
            ("m.16183A>C", bg(&["m.16182A>C"]), (RunClass::CG, 7)),
            ("m.302A>ACC", bg(&["m.310T>TC"]), (RunClass::CG, 9)),
            ("m.3243A>G", bg(&[]), (RunClass::CG, 3)),
        ];
        for (name, background, want) in cases {
            assert_eq!(implied_run(name, rcrs(), &background), Some(want), "{name} {background:?}");
        }
        // Deletions take the run beside the junction that lost the base, not the anchor's run.
        assert_eq!(implied_run("m.302AC>A", rcrs(), &[]), Some((RunClass::CG, 6)));
        assert_eq!(implied_run("m.310TC>T", rcrs(), &[]), Some((RunClass::CG, 4)));
        assert_eq!(implied_run("m.16183AC>A", rcrs(), &[]), Some((RunClass::CG, 4)));
        assert_eq!(implied_run("m.16183AC>A", rcrs(), &bg(&["m.16189T>C"])), Some((RunClass::CG, 9)));
        // A collapsed (--collapse-hp) tract allele is scored at what it changed: this one is m.310T>TC.
        assert_eq!(implied_run("m.302ACCCCCCCTCCCCC>ACCCCCCCTCCCCCC", rcrs(), &[]), Some((RunClass::CG, 6)));
        assert_eq!(implied_run("m.302ACCCCCCCTCCCCC>ACCCCCCCCCTCCCCC", rcrs(), &[]), Some((RunClass::CG, 9)));
        // A REF that is not the reference base is no allele.
        assert_eq!(implied_run("m.302C>CC", rcrs(), &[]), None);
        // A background SNV on the insertion's anchor base is applied: the anchor is not part of the edit, so 310C
        // joins C7 and C5 and the extra C makes 14.
        assert_eq!(implied_run("m.310T>TC", rcrs(), &bg(&["m.310T>C"])), Some((RunClass::CG, 14)));
        // Another insertion at the same point does clash with the allele, and is skipped.
        assert_eq!(implied_run("m.310T>TC", rcrs(), &bg(&["m.310T>TCC"])), Some((RunClass::CG, 6)));
    }

    #[test]
    fn implied_run_at_reports_where_the_run_lies_on_the_reference() {
        let bg = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        let extent = |name: &str, background: &[String]| {
            let r = implied_run_at(name, rcrs(), background).unwrap();
            (r.len, r.ref_start, r.ref_end)
        };
        // Inserted bases sit at their anchor, so the D310 run of m.302A>ACC starts at 302.
        assert_eq!(extent("m.302A>ACC", &[]), (9, 302, 309));
        assert_eq!(extent("m.310T>TC", &[]), (6, 310, 315));
        assert_eq!(extent("m.310TC>T", &[]), (4, 312, 315));
        // HVS1 on a 16182C/16183C/16189C background: one more C, anchored at 16181, runs to 16193.
        let hvs1 = bg(&["m.16182A>C", "m.16183A>C", "m.16189T>C"]);
        assert_eq!(extent("m.16181A>AC", &hvs1), (13, 16181, 16193));
        // On bare rCRS the same insertion is a lone C between As.
        assert_eq!(extent("m.16181A>AC", &[]), (1, 16181, 16181));
        // Insertions anchored on a base the background changes keep that change (16189C, 16183C).
        assert_eq!(extent("m.16189T>TC", &hvs1), (13, 16182, 16193));
        assert_eq!(extent("m.16183A>AC", &hvs1), (13, 16182, 16193));
    }

    #[test]
    fn tract_membership_follows_the_implied_run_for_indels_only() {
        let bg = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        let hvs1 = bg(&["m.16182A>C", "m.16183A>C", "m.16189T>C"]);
        let is_tract = |name: &str, background: &[String]| {
            let span = crate::scite::variant_spans(&[name.to_string()])[0];
            is_homopolymer_tract_allele(name, &span, implied_run_at(name, rcrs(), background).as_ref())
        };
        assert!(is_tract("m.16181A>AC", &hvs1), "the insertion lengthens the HVS1 C run");
        assert!(is_tract("m.16181A>ACCC", &hvs1));
        assert!(!is_tract("m.16181A>AC", &[]), "without the background it does not touch the tract");
        // An SNV that only borders the run keeps its own rates, even when the run reaches the tract.
        assert!(!is_tract("m.16182A>C", &bg(&["m.16183A>C"])));
        // REF overlap is enough on its own, as before.
        assert!(is_tract("m.302A>ACC", &[]) && is_tract("m.16183A>C", &[]) && is_tract("m.16189T>C", &[]));
        assert!(!is_tract("m.3243A>G", &[]) && !is_tract("m.3243A>AG", &[]));
        // Indels in the A runs beside the tracts reach the anchor A at 302 or 16183, not the C run: not tract alleles.
        for name in ["m.16179CA>C", "m.16179C>CA", "m.299CA>C", "m.299C>CA"] {
            assert!(!is_tract(name, &[]), "{name}");
        }
    }

    #[test]
    fn embedded_rcrs_matches_test_data() {
        let path = format!("{}/test_data/rCRS.fasta", env!("CARGO_MANIFEST_DIR"));
        let text = std::fs::read_to_string(path).unwrap();
        let seq: Vec<u8> = text
            .lines()
            .filter(|l| !l.starts_with('>'))
            .flat_map(|l| l.trim().bytes())
            .map(|b| b.to_ascii_uppercase())
            .collect();
        assert_eq!(rcrs().len(), 16569);
        assert_eq!(rcrs(), seq.as_slice());
    }

    #[test]
    fn confusion_rates_use_empirical_rows_and_flat_extension() {
        // C runs only: L4 n=1000, L5 n=500, L6 n=100 (too few, ignored).
        let rows = vec![
            ('C', 4, 4, '+', 900), ('C', 4, 3, '+', 60), ('C', 4, 5, '-', 40),
            ('C', 5, 5, '+', 400), ('C', 5, 4, '+', 75), ('C', 5, 6, '+', 25),
            ('C', 6, 6, '+', 100),
        ];
        let t = HpRateTable::from_confusion(&rows).unwrap();
        assert_eq!(t.calibrated_range(RunClass::CG), (4, 5));
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12;
        let r5 = t.allele_rates(RunClass::CG, 5);
        // fn = 1 - 400/500; fp = max(P(-1 | 6 -> 5) = 75/500, P(+1 | 4) = 40/1000)
        assert!(close(r5.fn_, 0.2) && close(r5.fp, 0.15) && !r5.beyond_calibration, "{r5:?}");
        let r9 = t.allele_rates(RunClass::CG, 9);
        assert!(close(r9.fn_, 0.2) && close(r9.fp, 0.15) && r9.beyond_calibration && r9.calibrated_max == 5, "{r9:?}");
        // No A/T rows: that class pools both classes, which here are the C rows.
        assert_eq!(t.allele_rates(RunClass::AT, 5), r5);
        // Clamped into [0.005, 0.5].
        let perfect = HpRateTable::from_confusion(&[('A', 4, 4, '+', 1000)]).unwrap();
        let r = perfect.allele_rates(RunClass::AT, 4);
        assert_eq!((r.fp, r.fn_), (HP_RATE_FLOOR, HP_RATE_FLOOR));
        assert!(HpRateTable::from_confusion(&[('A', 4, 4, '+', 199)]).is_none());
    }

    #[test]
    fn confusion_rates_bridge_a_gap_with_the_nearest_shorter_length() {
        // C/G L4 and L6 are calibrated, L5 is not: L5 reads L4, L6 reads L6, longer runs read L6.
        let rows = vec![
            ('C', 4, 4, '+', 900), ('C', 4, 3, '+', 100),
            ('C', 5, 5, '+', 50),
            ('C', 6, 6, '+', 300), ('C', 6, 5, '+', 100),
        ];
        let t = HpRateTable::from_confusion(&rows).unwrap();
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12;
        assert!(close(t.allele_rates(RunClass::CG, 5).fn_, 0.1));
        assert!(close(t.allele_rates(RunClass::CG, 6).fn_, 0.25));
        let r9 = t.allele_rates(RunClass::CG, 9);
        assert!(close(r9.fn_, 0.25) && r9.beyond_calibration && r9.calibrated_max == 6, "{r9:?}");
        assert!(close(t.allele_rates(RunClass::CG, 2).fn_, 0.1), "below the table: the shortest calibrated length");
    }

    #[test]
    fn estimate_snv_dropout_still_measures_dropout_at_max_hf_one() {
        let est = estimate_snv_dropout(&raw_from(&csv_with_pseudo_read(200, &fixed_rows(198, 2))), 1.0).unwrap();
        assert!((est.raw - 0.01).abs() < 1e-12, "gate is min(max_hf, 0.95): {est:?}");
    }

    #[test]
    fn explicit_hp_rates_skip_the_table() {
        let dir = std::env::temp_dir().join("himito_test_no_table");
        std::fs::create_dir_all(&dir).unwrap();
        let mut rows = fixed_rows(198, 2);
        rows.push(("m.302A>ACC", 100, 100));
        rows.push(("m.9100C>T", 50, 150));
        let matrix = dir.join("m.matrix.csv");
        std::fs::write(&matrix, csv_with_pseudo_read(200, &rows)).unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        let opts = RateOptions {
            data_type: "ont-r10".into(),
            hp_fp_rate: Some(0.2),
            hp_fn_rate: Some(0.3),
            hp_error: Some("/nonexistent/table.tsv".into()),
            ..Default::default()
        };
        // An unusable --hp-error does not matter when both HP rates are given.
        start(matrix.to_str().unwrap(), None, 0.01, 0.95, &opts, &[], 200, 1, false, &prefix).unwrap();
        let table = std::fs::read_to_string(format!("{prefix}.error_model.tsv")).unwrap();
        assert!(table.contains("hp\tm.302A>ACC\tfn\t0.300000\tcli\trun=C/G:9\n"), "{table}");
        assert!(!table.contains("hp_table"), "{table}");
    }

    #[test]
    fn embedded_default_tables_parse() {
        for text in [DEFAULT_HP_TABLE_ONT, DEFAULT_HP_TABLE_PACBIO] {
            let t = HpRateTable::from_text(text).unwrap().unwrap();
            assert_eq!(t.calibrated_range(RunClass::CG), (4, 6));
        }
    }

    #[test]
    fn load_hp_table_prefers_the_flag_then_the_sibling_then_the_default() {
        let dir = std::env::temp_dir().join("himito_test_hp_table");
        std::fs::create_dir_all(&dir).unwrap();
        let matrix = dir.join("s.matrix.csv").to_str().unwrap().to_string();
        let sibling = dir.join("s.hp_confusion.tsv");
        let _ = std::fs::remove_file(&sibling);
        let opts = |hp_error: Option<String>| RateOptions { data_type: "pacbio".into(), hp_error, ..Default::default() };
        assert_eq!(load_hp_table(&opts(None), &matrix).unwrap().1, "default:HG00097_pacbio");
        std::fs::write(&sibling, DEFAULT_HP_TABLE_ONT).unwrap();
        assert_eq!(load_hp_table(&opts(None), &matrix).unwrap().1, format!("table:{}", sibling.display()));
        let flag = dir.join("flag.tsv");
        std::fs::write(&flag, DEFAULT_HP_TABLE_PACBIO).unwrap();
        let flag = flag.to_str().unwrap().to_string();
        assert_eq!(load_hp_table(&opts(Some(flag.clone())), &matrix).unwrap().1, format!("table:{flag}"));
        // A flag pointing at an unusable file is an error, not a silent fallback.
        let bad = dir.join("bad.tsv");
        std::fs::write(&bad, "base\tref_len\tobs_len\tstrand\tn\nA\t4\t4\t+\t3\n").unwrap();
        assert!(load_hp_table(&opts(Some(bad.to_str().unwrap().into())), &matrix).is_err());
    }

    /// A filtered matrix over `names` (4 reads, two alt) and a raw matrix with no near-fixed rows.
    fn model_inputs(names: &[&str]) -> (BinaryMatrix, RawMatrix) {
        let data: Vec<Vec<Option<u8>>> = names.iter().map(|_| cells("1100")).collect();
        let reads: Vec<String> = (1..=4).map(|r| format!("r{r}")).collect();
        let variants: Vec<String> = names.iter().map(|s| s.to_string()).collect();
        (
            BinaryMatrix { variants: variants.clone(), reads: reads.clone(), data: data.clone() },
            RawMatrix { variants, reads, data },
        )
    }

    fn model_for(names: &[&str], opts: &RateOptions) -> (crate::scite::VariantModel, Vec<ErrorModelRow>) {
        let (binary, raw) = model_inputs(names);
        let resolved = resolve_rates(opts, &raw, 0.95);
        let table = HpRateTable::from_text(DEFAULT_HP_TABLE_ONT).unwrap().unwrap();
        build_variant_model(&binary, &raw, &resolved, opts, Some((&table, "table:test")))
    }

    #[test]
    fn tract_alleles_take_table_rates_and_others_keep_theirs() {
        let opts = RateOptions { data_type: "ont-r10".into(), tract_alleles: TractAlleleMode::Search, ..Default::default() };
        let (model, rows) = model_for(&["m.302A>ACC", "m.3243A>G", "m.514CA>C"], &opts);
        let close = |a: f64, b: f64| (a - b).abs() < 1e-6;
        assert!(close(model.rates[0].0, 0.081778) && close(model.rates[0].1, 0.118002), "{:?}", model.rates[0]);
        assert_eq!(model.rates[1], (0.005, 0.02), "an SNV keeps the SNV rates (preset fn: no near-fixed SNVs here)");
        assert_eq!(model.rates[2], (0.01, 0.04), "CA514 keeps twice the preset");
        assert!(model.tract.iter().flatten().all(|t| t.placed), "Search places every tract allele");
        assert!(model.tract[0].is_some() && model.tract[1].is_none(), "D310 is recognised, the SNV is not");
        let detail = &rows.iter().find(|r| r.variant == "m.302A>ACC").unwrap().detail;
        assert_eq!(detail, "run=C/G:9 calibrated_max=6 beyond_calibration=yes");
    }

    #[test]
    fn hvs1_insertion_anchored_before_the_tract_takes_tract_rates() {
        // With 16182C, 16183C and 16189C near-fixed, the HVS1 C run spans 16182-16193 and VCF left-normalisation
        // anchors a one-C insertion at 16181, outside HVS1-C (16183-16193).
        let reads: Vec<String> = (1..=4).map(|r| format!("r{r}")).collect();
        let ins = "m.16181A>AC";
        let binary = BinaryMatrix { variants: vec![ins.into()], reads: reads.clone(), data: vec![cells("1100")] };
        let raw_with = |background: &[&str]| {
            let mut variants: Vec<String> = background.iter().map(|s| s.to_string()).collect();
            let mut data: Vec<Vec<Option<u8>>> = background.iter().map(|_| cells("1111")).collect();
            variants.push(ins.into());
            data.push(cells("1100"));
            RawMatrix { variants, reads: reads.clone(), data }
        };
        let opts = RateOptions { data_type: "ont-r10".into(), ..Default::default() };
        let table = HpRateTable::from_text(DEFAULT_HP_TABLE_ONT).unwrap().unwrap();
        let model = |raw: &RawMatrix| {
            let resolved = resolve_rates(&opts, raw, 0.95);
            let (m, rows) = build_variant_model(&binary, raw, &resolved, &opts, Some((&table, "table:test")));
            (m, rows, resolved)
        };

        let (m, rows, _) = model(&raw_with(&["m.16182A>C", "m.16183A>C", "m.16189T>C"]));
        let close = |a: f64, b: f64| (a - b).abs() < 1e-6;
        assert!(close(m.rates[0].0, 0.081778) && close(m.rates[0].1, 0.118002), "{:?}", m.rates[0]);
        let detail = &rows.iter().find(|r| r.variant == ins).unwrap().detail;
        assert_eq!(detail, "run=C/G:13 calibrated_max=6 beyond_calibration=yes");
        let a = m.tract[0].clone().expect("a C13 insertion at HVS1-C is a tract allele");
        assert_eq!((a.run_class.as_str(), a.run_length, a.placed), ("C/G", 13, true));

        // Without that background the insertion is an ordinary indel with the SNV rates.
        let (m, rows, resolved) = model(&raw_with(&[]));
        assert_eq!(m.rates[0], (resolved.fp_rate, resolved.fn_rate));
        assert!(m.tract[0].is_none() && rows.iter().all(|r| r.variant != ins));
    }

    #[test]
    fn hp_table_kernel_is_the_length_change_histogram() {
        let table = HpRateTable::from_text(DEFAULT_HP_TABLE_PACBIO).unwrap().unwrap();
        let k6 = table.kernel(RunClass::CG, 6);
        assert!((k6.iter().map(|(_, p)| p).sum::<f64>() - 1.0).abs() < 1e-12);
        let p0 = k6.iter().find(|(d, _)| *d == 0).unwrap().1;
        assert!((p0 - (1.0 - table.allele_rates(RunClass::CG, 6).fn_)).abs() < 1e-12, "P(Δ = 0) is 1 - fn");
        assert_eq!(table.kernel(RunClass::CG, 13), k6, "lengths past calibration use the longest calibrated one");
        assert_eq!(table.kernel(RunClass::CG, 2), table.kernel(RunClass::CG, 4), "and shorter ones the shortest");
    }

    #[test]
    fn placed_tract_alleles_at_one_anchor_form_one_site() {
        let names = ["m.302A>ACC", "m.302A>ACCC", "m.310T>TC", "m.3243A>G"];
        let opts = |dt: &str, mode| RateOptions { data_type: dt.into(), tract_alleles: mode, ..Default::default() };
        let (model, rows) = model_for(&names, &opts("pacbio", TractAlleleMode::Auto));
        let shape: Vec<(Vec<usize>, Vec<i64>, i64)> =
            model.sites.iter().map(|s| (s.members.clone(), s.member_len.clone(), s.root_len)).collect();
        assert_eq!(shape, [(vec![0, 1], vec![9, 10], 7), (vec![2], vec![6], 5)]);
        let lens: Vec<i64> = model.sites[0].kernel.iter().map(|k| k.0).collect();
        assert_eq!(lens, [7, 9, 10], "a kernel for every state's length");
        assert!(rows.iter().any(|r| r.scope == "site" && r.variant == "m.302A>ACC,m.302A>ACCC" && r.value == 7.0));
        // ONT auto places the C9/C10 alleles too: the data-type string no longer decides it.
        let ont = model_for(&names, &opts("ont-r10", TractAlleleMode::Auto)).0;
        assert_eq!(ont.sites.iter().map(|s| s.members.clone()).collect::<Vec<_>>(), [vec![0, 1], vec![2]]);
        let ont_site = model_for(&names, &opts("ont-r10", TractAlleleMode::Site)).0;
        assert_eq!(ont_site.sites.len(), 2);
        assert!(ont_site.tract.iter().flatten().all(|t| t.placed));
        // search keeps one binary row per allele.
        assert!(model_for(&names, &opts("pacbio", TractAlleleMode::Search)).0.sites.is_empty());
        assert_eq!("site".parse::<TractAlleleMode>(), Ok(TractAlleleMode::Site));
    }

    #[test]
    fn only_pure_length_alleles_are_site_members() {
        for (name, want) in [
            ("m.302A>ACC", Some(2)),
            ("m.310TC>T", Some(-1)),
            ("m.16183A>AC", Some(1)),
            ("m.16189T>C", None),
            ("m.310T>C", None),
            // Collapsed tract alleles: one run changed is a length allele, both runs is not.
            ("m.302ACCCCCCCTCCCCC>ACCCCCCCCTCCCCC", Some(1)),
            ("m.302ACCCCCCCTCCCCC>ACCCCCCCCTCCCCCC", None),
        ] {
            assert_eq!(homopolymer_length_change(name), want, "{name}");
        }
        let opts = RateOptions { data_type: "pacbio".into(), tract_alleles: TractAlleleMode::Site, ..Default::default() };
        let (model, _) = model_for(&["m.16183A>AC", "m.16183A>C", "m.16189T>C"], &opts);
        assert_eq!(model.sites.iter().map(|s| s.members.clone()).collect::<Vec<_>>(), [vec![0]]);
    }

    #[test]
    fn length_sites_do_not_depend_on_allele_order_and_take_one_allele_per_length() {
        let opts = RateOptions { data_type: "pacbio".into(), tract_alleles: TractAlleleMode::Site, ..Default::default() };
        // m.302A>AT inserts a lone T (no homopolymer): listed first, it must not set the reference run.
        for names in [["m.302A>AT", "m.302A>ACC", "m.302A>ACCC"], ["m.302A>ACC", "m.302A>ACCC", "m.302A>AT"]] {
            let (model, rows) = model_for(&names, &opts);
            let members: Vec<Vec<&str>> =
                model.sites.iter().map(|s| s.members.iter().map(|&v| names[v]).collect()).collect();
            assert_eq!(members, [vec!["m.302A>ACC", "m.302A>ACCC"]], "{names:?}");
            assert_eq!(model.sites[0].root_len, 7);
            assert!(rows.iter().filter(|r| r.scope == "site").all(|r| !r.variant.is_empty()));
        }
        // Two spellings of one length: the first is the site's, the other stays its own row.
        let (model, _) = model_for(&["m.302A>ACC", "m.302AC>ACCC"], &opts);
        assert_eq!(model.sites.iter().map(|s| s.members.clone()).collect::<Vec<_>>(), [vec![0]]);
    }

    #[test]
    fn one_cli_hp_rate_turns_sites_off() {
        for (fp, fnr) in [(Some(0.2), None), (None, Some(0.3)), (Some(0.2), Some(0.3))] {
            let opts = RateOptions { data_type: "pacbio".into(), hp_fp_rate: fp, hp_fn_rate: fnr, ..Default::default() };
            assert!(model_for(&["m.302A>ACC", "m.302A>ACCC"], &opts).0.sites.is_empty(), "{fp:?} {fnr:?}");
        }
    }

    #[test]
    fn the_site_root_length_reaches_the_model() {
        // C8 reads (m.302A>AC, a raw row the filters dropped) outnumber REF calls: the D310 site roots at C8.
        let reads: Vec<String> = (1..=8).map(|r| format!("r{r}")).collect();
        let binary = BinaryMatrix {
            variants: vec!["m.302A>ACC".into(), "m.302A>ACCC".into()],
            reads: reads.clone(),
            data: vec![cells("11000000"), cells("00110000")],
        };
        let raw = RawMatrix {
            variants: vec!["m.302A>AC".into(), "m.302A>ACC".into(), "m.302A>ACCC".into()],
            reads,
            data: vec![cells("00001110"), cells("11000000"), cells("00110000")],
        };
        let opts = RateOptions { data_type: "pacbio".into(), ..Default::default() };
        let resolved = resolve_rates(&opts, &raw, 0.95);
        let table = HpRateTable::from_text(DEFAULT_HP_TABLE_PACBIO).unwrap().unwrap();
        let (model, _) = build_variant_model(&binary, &raw, &resolved, &opts, Some((&table, "table:test")));
        assert_eq!(model.sites[0].root_len, 8);
    }

    #[test]
    fn a_site_roots_at_the_commonest_length_that_is_not_one_of_its_alleles() {
        // Reads at D310: C8 (m.302A>AC, filtered out of the search) x5, C9 (the site's allele) x4, REF C7 x2.
        let csv = "variant,chrM,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11\n\
                   m.302A>AC,0,1,1,1,1,1,0,0,0,0,0,0\n\
                   m.302A>ACC,0,0,0,0,0,0,1,1,1,1,0,0\n";
        let raw = read_raw_matrix(Cursor::new(csv)).unwrap();
        assert_eq!(site_root_length(&raw, 302, 7, &["m.302A>ACC"]), 8);
        // With C8 at two reads REF wins, and the pseudo-read never counts.
        let csv = "variant,chrM,r1,r2,r3,r4,r5\nm.302A>AC,1,1,1,0,0,0\nm.302A>ACC,1,0,0,0,0,0\n";
        let raw = read_raw_matrix(Cursor::new(csv)).unwrap();
        assert_eq!(site_root_length(&raw, 302, 7, &["m.302A>ACC"]), 7);
    }

    #[test]
    fn ont_site_mode_runs_end_to_end_on_the_embedded_table() {
        let dir = std::env::temp_dir().join("himito_lineage_ont_site");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let matrix = dir.join("matrix.csv");
        std::fs::write(
            &matrix,
            tripled_csv(&[("m.302A>ACC", "1100"), ("m.302A>ACCC", "0011"), ("m.16181A>AC", "1010"), ("m.100A>G", "1000")]),
        )
        .unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        let opts = RateOptions { data_type: "ont-r10".into(), tract_alleles: TractAlleleMode::Site, ..Default::default() };
        start(matrix.to_str().unwrap(), None, 0.01, 0.95, &opts, &[], 50, 1, false, &prefix).unwrap();
        let model = std::fs::read_to_string(format!("{prefix}.error_model.tsv")).unwrap();
        assert!(model.contains("site\tm.302A>ACC,m.302A>ACCC\troot_len\t"), "{model}");
        let tree = std::fs::read_to_string(format!("{prefix}.mutation_tree.tsv")).unwrap();
        assert!(tree.contains("m.302A>ACC") && !tree.contains("TRACT_ANNOTATED"), "{tree}");
    }

    #[test]
    fn cli_hp_rates_override_the_table() {
        let opts = RateOptions {
            data_type: "ont-r10".into(),
            hp_fp_rate: Some(0.2),
            hp_fn_rate: Some(0.3),
            ..Default::default()
        };
        let (model, rows) = model_for(&["m.302A>ACC"], &opts);
        assert_eq!(model.rates[0], (0.2, 0.3));
        assert!(rows.iter().filter(|r| r.variant == "m.302A>ACC").all(|r| r.source == "cli"));
        assert!(model.tract[0].as_ref().is_some_and(|t| t.placed), "auto places alleles with explicit HP rates");
    }

    #[test]
    fn auto_places_beyond_calibration_alleles_on_every_platform() {
        let opts = |dt: &str| RateOptions { data_type: dt.into(), tract_alleles: TractAlleleMode::Auto, ..Default::default() };
        for dt in ["ont-r10", "pacbio"] {
            let (model, _rows) = model_for(&["m.302A>ACC"], &opts(dt));
            let t = model.tract[0].as_ref().expect("a D310 length allele is a tract allele");
            assert!(t.placed, "{dt}: auto must place the allele");
            assert_eq!(t.reason, "beyond_calibration (run 9 > calibrated 6)", "{dt}");
        }
    }

    #[test]
    fn annotate_mode_still_withholds_tract_alleles_on_every_platform() {
        let opts = |dt: &str| RateOptions {
            data_type: dt.into(),
            tract_alleles: TractAlleleMode::Annotate,
            ..Default::default()
        };
        for dt in ["ont-r10", "pacbio"] {
            let (model, _rows) = model_for(&["m.302A>ACC"], &opts(dt));
            let t = model.tract[0].as_ref().expect("a D310 length allele is a tract allele");
            assert!(!t.placed, "{dt}: annotate must withhold the allele");
            assert_eq!(t.reason, "--tract-alleles annotate", "{dt}");
        }
    }

    #[test]
    fn a_calibrated_tract_allele_is_placed_and_says_so() {
        let opts = RateOptions { data_type: "ont-r10".into(), tract_alleles: TractAlleleMode::Auto, ..Default::default() };
        // m.310T>TC implies a C/G run of 6, which the table calibrates.
        let (model, _rows) = model_for(&["m.310T>TC"], &opts);
        let t = model.tract[0].as_ref().expect("D310 allele");
        assert!(t.placed);
        assert_eq!(t.reason, "calibrated (run 6 <= calibrated 6)");
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
        let opts = RateOptions {
            data_type: "pacbio".to_string(),
            fp_rate: Some(0.01),
            fn_rate: Some(0.1),
            hp_fp_rate: Some(0.01),
            hp_fn_rate: Some(0.1),
            ..Default::default()
        };
        start(matrix.to_str().unwrap(), None, 0.01, 0.95, &opts, tracts, 20, 1, debug, &prefix).unwrap();
        prefix
    }

    #[test]
    fn debug_gates_the_raw_outputs() {
        let csv = tripled_csv(&[("m.100A>G", "1100"), ("m.200C>T", "1010")]);
        for debug in [false, true] {
            let prefix = run_small_lineage(&format!("himito_lineage_debug_{debug}"), &csv, debug, &[]);
            for suffix in [
                ".raw_haplotype_map.tsv",
                ".raw_variant_cooccurrence.tsv",
                ".variant_disposition.tsv",
                ".tract_annotation.tsv",
            ] {
                assert_eq!(Path::new(&format!("{prefix}{suffix}")).exists(), debug, "{suffix} with debug={debug}");
            }
            assert!(Path::new(&format!("{prefix}.cleaned_variant_cooccurrence.tsv")).exists());
            assert!(!Path::new(&format!("{prefix}.cleaned_matrix.csv")).exists());
        }
    }

    #[test]
    fn start_reads_tract_membership_on_the_majority_background() {
        let hp_row = |prefix: &str| {
            std::fs::read_to_string(format!("{prefix}.error_model.tsv"))
                .unwrap()
                .lines()
                .find(|l| l.starts_with("hp\tm.16181A>AC\tfn"))
                .map(str::to_string)
        };
        // 16182C/16183C/16189C in 3 of 4 reads (HF 0.75: majority, not near-fixed) make 16181A>AC an HVS1 allele.
        let csv = tripled_csv(&[
            ("m.100A>G", "1010"),
            ("m.16181A>AC", "1100"),
            ("m.16182A>C", "1110"),
            ("m.16183A>C", "1110"),
            ("m.16189T>C", "1110"),
        ]);
        let prefix = run_small_lineage("himito_lineage_hvs1_majority", &csv, false, &[]);
        let row = hp_row(&prefix).expect("16181A>AC takes tract rates");
        assert!(row.ends_with("\tcli\trun=C/G:13"), "{row}");
        // Without that background it is an ordinary indel.
        let csv = tripled_csv(&[("m.100A>G", "1010"), ("m.16181A>AC", "1100")]);
        let prefix = run_small_lineage("himito_lineage_hvs1_bare", &csv, false, &[]);
        assert_eq!(hp_row(&prefix), None);
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
        // A ceiling below 0.50 drops both tract alleles and the whitelisted SNV (matrix HF 0.50) on the VCF path.
        assert!(kept(vcf(0.4)).is_empty());
        assert!(kept(HfFilter::FromMatrix { min: 0.2, max: 0.4 }).is_empty());
    }

    #[test]
    fn a_whitelisted_snv_still_has_to_clear_the_matrix_hf_window() {
        // One row, 99 alt of 100 covered reads: matrix HF 0.99, above max_hf.
        let mut row: Vec<Option<u8>> = vec![Some(1); 99];
        row.push(Some(0));
        let raw = RawMatrix {
            variants: vec!["m.73A>G".to_string()],
            reads: (0..100).map(|i| format!("r{i}")).collect(),
            data: vec![row],
        };
        // The VCF lists it (its own VCF HF was 0.948, inside the window) but the matrix says near-fixed.
        let hf = HfFilter::FromVcf { map: [("m.73A>G".to_string(), 0.948)].into(), min: 0.01, max: 0.95 };
        assert!(filter_matrix(&raw, &hf, KNOWN_TRACTS).variants.is_empty());

        // The same row inside the window is kept.
        let hf = HfFilter::FromVcf { map: [("m.73A>G".to_string(), 0.948)].into(), min: 0.01, max: 0.999 };
        assert_eq!(filter_matrix(&raw, &hf, KNOWN_TRACTS).variants, vec!["m.73A>G".to_string()]);
    }

    #[test]
    fn a_variant_outside_the_whitelist_is_dropped_whatever_its_matrix_hf() {
        let raw = RawMatrix {
            variants: vec!["m.3243A>G".to_string()],
            reads: (0..10).map(|i| format!("r{i}")).collect(),
            data: vec![(0..10).map(|i| Some(u8::from(i < 5))).collect()],
        };
        let hf = HfFilter::FromVcf { map: [("m.73A>G".to_string(), 0.5)].into(), min: 0.01, max: 0.95 };
        assert!(filter_matrix(&raw, &hf, KNOWN_TRACTS).variants.is_empty());
    }

    #[test]
    fn disposition_names_why_each_candidate_left_the_analysis() {
        let raw = RawMatrix {
            variants: vec!["m.100A>G".to_string(), "m.200A>G".to_string(), "m.300A>G".to_string()],
            reads: (0..10).map(|i| format!("r{i}")).collect(),
            data: vec![
                (0..10).map(|i| Some(u8::from(i < 5))).collect(), // HF 0.5, kept
                vec![Some(1); 10],                                // HF 1.0, above the window
                (0..10).map(|i| Some(u8::from(i < 1))).collect(), // 1 present, below MIN_PRESENCE
            ],
        };
        let hf = HfFilter::FromVcf {
            map: [("m.100A>G".to_string(), 0.5), ("m.200A>G".to_string(), 0.9), ("m.300A>G".to_string(), 0.1)].into(),
            min: 0.01,
            max: 0.95,
        };
        let rows = variant_dispositions(&raw, &hf, KNOWN_TRACTS, &[]);
        let by = |n: &str| rows.iter().find(|r| r.variant == n).unwrap().outcome;
        assert_eq!(by("m.100A>G"), "informative");
        assert_eq!(by("m.200A>G"), "dropped_hf_window");
        assert_eq!(by("m.300A>G"), "dropped_min_presence");
    }

    #[test]
    fn disposition_reports_a_whitelisted_site_with_no_matrix_row() {
        let raw = RawMatrix { variants: Vec::new(), reads: vec!["r0".to_string()], data: Vec::new() };
        let hf = HfFilter::FromVcf { map: [("m.16368T>C".to_string(), 0.09)].into(), min: 0.01, max: 0.95 };
        let rows = variant_dispositions(&raw, &hf, KNOWN_TRACTS, &[]);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].outcome, "no_matrix_row");
        assert_eq!(rows[0].vcf_hf, Some(0.09));
    }

    /// Ten reads, one row per variant: `n_alt` alt calls, then `n_ref` ref calls, the rest uncovered.
    fn counted_raw(rows: &[(&str, usize, usize)]) -> RawMatrix {
        RawMatrix {
            variants: rows.iter().map(|r| r.0.to_string()).collect(),
            reads: (0..10).map(|i| format!("r{i}")).collect(),
            data: rows
                .iter()
                .map(|&(_, alt, rf)| {
                    (0..10).map(|i| if i < alt { Some(1) } else if i < alt + rf { Some(0) } else { None }).collect()
                })
                .collect(),
        }
    }

    #[test]
    fn disposition_covers_every_remaining_outcome() {
        let raw = counted_raw(&[
            ("m.100A>G", 5, 5),  // whitelisted, HF 0.5
            ("m.200A>G", 5, 5),  // same counts, but the root block takes it
            ("m.300A>G", 5, 5),  // same counts, not in the VCF
            ("m.400A>G", 10, 0), // whitelisted, no ref call at all (HF 1.0 is inside the window below)
            ("m.500A>G", 0, 0),  // whitelisted, never covered
        ]);
        let map: HfMap = ["m.100A>G", "m.200A>G", "m.400A>G", "m.500A>G"].map(|n| (n.to_string(), 0.5)).into();
        let hf = HfFilter::FromVcf { map, min: 0.01, max: 1.01 };
        let rows = variant_dispositions(&raw, &hf, KNOWN_TRACTS, &["m.200A>G".to_string()]);
        let row = |n: &str| rows.iter().find(|r| r.variant == n).unwrap();

        assert_eq!(row("m.100A>G").outcome, "informative");
        assert_eq!((row("m.100A>G").n_present, row("m.100A>G").n_absent, row("m.100A>G").matrix_hf), (5, 5, 0.5));
        assert_eq!(row("m.100A>G").vcf_hf, Some(0.5));
        assert_eq!(row("m.200A>G").outcome, "root_block");
        assert_eq!(row("m.300A>G").outcome, "dropped_not_whitelisted");
        assert_eq!(row("m.300A>G").vcf_hf, None);
        assert_eq!(row("m.400A>G").outcome, "dropped_min_absence");
        assert_eq!(row("m.500A>G").outcome, "dropped_min_presence");
        assert_eq!(rows.len(), 5);
        assert!(rows.windows(2).all(|w| w[0].variant <= w[1].variant), "sorted by variant");
    }

    #[test]
    fn disposition_marks_tract_alleles_the_search_withheld() {
        let raw = counted_raw(&[("m.100A>G", 5, 5), ("m.16189T>TC", 5, 5)]);
        let hf = HfFilter::FromMatrix { min: 0.01, max: 0.95 };
        let mut rows = variant_dispositions(&raw, &hf, &[], &[]);
        mark_tract_annotated(&mut rows, &["m.16189T>TC".to_string()]);
        let outcomes: Vec<(&str, &str)> = rows.iter().map(|r| (r.variant.as_str(), r.outcome)).collect();
        assert_eq!(outcomes, [("m.100A>G", "informative"), ("m.16189T>TC", "tract_annotated")]);
    }

    #[test]
    fn disposition_agrees_with_filter_matrix_on_every_boundary() {
        // Every gate at its edge: presence 1 and 2, absence 0 and 1, HF exactly at min and at max, a tract
        // allele and an unlisted name, under both HF modes.
        let raw = counted_raw(&[
            ("m.100A>G", 1, 9),
            ("m.101A>G", 2, 8),
            ("m.102A>G", 9, 0),
            ("m.103A>G", 9, 1),
            ("m.104A>G", 5, 5),
            ("m.105A>G", 2, 0),
            ("m.106A>G", 0, 10),
            ("m.16189T>TC", 5, 5),
        ]);
        let listed: HfMap = raw.variants.iter().filter(|v| !v.starts_with("m.104")).map(|v| (v.clone(), 0.5)).collect();
        let filters = [
            HfFilter::FromVcf { map: listed, min: 0.2, max: 0.9 },
            HfFilter::FromMatrix { min: 0.5, max: 0.9 },
            HfFilter::FromMatrix { min: 0.1, max: 0.5 },
        ];
        for hf in &filters {
            for tracts in [KNOWN_TRACTS, &[][..]] {
                let kept = filter_matrix(&raw, hf, tracts).variants;
                let rows = variant_dispositions(&raw, hf, tracts, &[]);
                let informative: Vec<String> =
                    rows.iter().filter(|r| r.outcome == "informative").map(|r| r.variant.clone()).collect();
                let mut kept_sorted = kept.clone();
                kept_sorted.sort();
                assert_eq!(informative, kept_sorted, "{:?}", rows);
                assert!(rows.iter().all(|r| r.outcome != "no_matrix_row"));
            }
        }
    }

    #[test]
    fn write_variant_disposition_writes_one_tab_separated_row_each() {
        let rows = vec![
            VariantDisposition {
                variant: "m.100A>G".to_string(),
                vcf_hf: Some(0.25),
                matrix_hf: 0.5,
                n_present: 5,
                n_absent: 5,
                outcome: "informative",
                detail: String::new(),
            },
            VariantDisposition {
                variant: "m.200A>G".to_string(),
                vcf_hf: None,
                matrix_hf: 1.0,
                n_present: 10,
                n_absent: 0,
                outcome: "dropped_min_absence",
                detail: "0 absent < 1".to_string(),
            },
        ];
        let path = std::env::temp_dir().join("himito_test_variant_disposition.tsv");
        write_variant_disposition(&rows, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(
            content,
            "variant\tvcf_hf\tmatrix_hf\tn_present\tn_absent\toutcome\tdetail\n\
             m.100A>G\t0.250000\t0.500000\t5\t5\tinformative\tNA\n\
             m.200A>G\tNA\t1.000000\t10\t0\tdropped_min_absence\t0 absent < 1\n"
        );
    }

    #[test]
    fn start_with_no_surviving_variant_still_writes_the_disposition_table() {
        let dir = std::env::temp_dir().join("himito_test_disposition_empty");
        std::fs::create_dir_all(&dir).unwrap();
        let matrix = dir.join("m.matrix.csv");
        // One variant at HF 1.0: above the window, so nothing survives.
        std::fs::write(&matrix, csv_with_pseudo_read(20, &[("m.9000A>G", 20, 0)])).unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        let opts = RateOptions { data_type: "pacbio".to_string(), ..Default::default() };
        let err = start(matrix.to_str().unwrap(), None, 0.01, 0.95, &opts, &[], 20, 1, true, &prefix).unwrap_err();
        assert!(err.to_string().contains("variant_disposition.tsv"), "{err}");
        let table = std::fs::read_to_string(format!("{prefix}.variant_disposition.tsv")).unwrap();
        assert!(table.contains("m.9000A>G\tNA\t"), "{table}");
        assert!(table.contains("dropped_min_absence") || table.contains("dropped_hf_window"), "{table}");
    }

    #[test]
    fn a_run_that_fails_after_the_provisional_write_leaves_survivors_marked_pending() {
        let dir = std::env::temp_dir().join("himito_test_disposition_pending");
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();
        let mut rows = fixed_rows(198, 2);
        rows.push(("m.9000A>G", 100, 100));
        let matrix = dir.join("m.matrix.csv");
        std::fs::write(&matrix, csv_with_pseudo_read(200, &rows)).unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        // A directory where error_model.tsv belongs: the run gets past filtering and the provisional write,
        // then fails before the tree search.
        std::fs::create_dir(format!("{prefix}.error_model.tsv")).unwrap();
        let opts = RateOptions { data_type: "pacbio".to_string(), ..Default::default() };
        start(matrix.to_str().unwrap(), None, 0.01, 0.95, &opts, &[], 100, 1, true, &prefix).unwrap_err();
        let table = std::fs::read_to_string(format!("{prefix}.variant_disposition.tsv")).unwrap();
        let row = |n: &str| -> Vec<String> {
            table.lines().find(|l| l.starts_with(&format!("{n}\t"))).unwrap().split('\t').map(String::from).collect()
        };
        assert_eq!(row("m.9000A>G")[5..], ["pending", "the tree search had not run"], "{table}");
        assert!(!table.contains("informative"), "{table}");
        // A row a filter dropped is already final and reads as it does after the search.
        assert_eq!(row("m.1000A>G")[5], "dropped_hf_window", "{table}");
        assert_ne!(row("m.1000A>G")[6], "the tree search had not run", "{table}");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn start_writes_a_disposition_table_that_names_the_root_block() {
        let dir = std::env::temp_dir().join("himito_test_disposition_run");
        std::fs::create_dir_all(&dir).unwrap();
        let mut rows = fixed_rows(198, 2);
        rows.push(("m.9000A>G", 100, 100));
        rows.push(("m.9100C>T", 50, 150));
        let matrix = dir.join("m.matrix.csv");
        std::fs::write(&matrix, csv_with_pseudo_read(200, &rows)).unwrap();
        let prefix = dir.join("out").to_str().unwrap().to_string();
        let opts = RateOptions { data_type: "pacbio".to_string(), ..Default::default() };
        start(matrix.to_str().unwrap(), None, 0.01, 0.95, &opts, &[], 100, 1, true, &prefix).unwrap();
        let table = std::fs::read_to_string(format!("{prefix}.variant_disposition.tsv")).unwrap();
        let outcome = |n: &str| {
            table.lines().find(|l| l.starts_with(&format!("{n}\t"))).unwrap().split('\t').nth(5).unwrap().to_string()
        };
        assert_eq!(outcome("m.9000A>G"), "informative", "{table}");
        assert_eq!(outcome("m.9100C>T"), "informative", "{table}");
        assert!(!table.contains("pending"), "{table}");
        // The six fixed SNVs sit at HF 0.99: above the 0.95 ceiling.
        assert_eq!(outcome("m.1000A>G"), "dropped_hf_window", "{table}");
    }
}
