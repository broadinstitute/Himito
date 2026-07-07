use std::io::{BufWriter, Write};
use std::fs::File;
use anyhow::{Context, Result};
use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::lineage::{self, BinaryMatrix, HaplotypeMatrix, load_and_filter_matrix};

/// A SCITE-style mutation tree: `n_mutations` mutation nodes (ids
/// `0..n_mutations`, in the same order as `BinaryMatrix::variants`) plus one
/// root node (id `n_mutations`, the unmutated germline state).
/// `parent[i]` is node `i`'s parent; the root is its own parent.
#[derive(Debug, Clone)]
pub struct MutationTree {
    pub n_mutations: usize,
    pub parent: Vec<usize>,
}

impl MutationTree {
    pub fn root(&self) -> usize {
        self.n_mutations
    }

    /// Bitmask of mutation ids on the path from `node` to the root,
    /// inclusive of `node` itself if it is a mutation node. Root's mask is 0.
    pub fn ancestor_mask(&self, node: usize) -> u64 {
        let mut mask = 0u64;
        let mut cur = node;
        while cur != self.root() {
            mask |= 1u64 << cur;
            cur = self.parent[cur];
        }
        mask
    }

    /// Bitmask of all proper descendants of `node` (not including itself).
    ///
    /// Note: `swap_subtrees` (Task 8) builds its own "descendants" list using
    /// `ancestor_mask` with a filter that *does* include `node` itself — those
    /// two conventions are intentionally different. This method is only used in
    /// `propose_prune_reattach` to exclude a node's own subtree from reattachment targets.
    fn descendant_mask(&self, node: usize) -> u64 {
        let children = self.children_of();
        let mut mask = 0u64;
        let mut stack = vec![node];
        while let Some(cur) = stack.pop() {
            for &c in &children[cur] {
                if c < self.n_mutations {
                    mask |= 1u64 << c;
                }
                stack.push(c);
            }
        }
        mask
    }

    fn children_of(&self) -> Vec<Vec<usize>> {
        let mut children = vec![Vec::new(); self.parent.len()];
        for (node, &p) in self.parent.iter().enumerate() {
            if node != p {
                children[p].push(node);
            }
        }
        children
    }

    /// Returns a new tree identical to `self` except `node`'s parent is
    /// `new_parent`. Does not validate acyclicity — callers must only pass a
    /// `new_parent` outside `node`'s own subtree (see `descendant_mask`).
    fn with_parent(&self, node: usize, new_parent: usize) -> MutationTree {
        let mut parent = self.parent.clone();
        parent[node] = new_parent;
        MutationTree { n_mutations: self.n_mutations, parent }
    }

    /// A random valid tree: process mutations in random order, attaching
    /// each one to a uniformly random *already-placed* node (root counts as
    /// placed from the start). This is a "random recursive tree" — every
    /// node it produces is acyclic by construction. It is not uniformly
    /// distributed over all labeled trees the way reference SCITE's
    /// Prüfer-sequence generator is, but that distinction doesn't matter for
    /// MCMC correctness (see the plan header) — only for how fast chains
    /// that DO start randomly mix, which the multi-chain restarts (Task 10)
    /// already compensate for.
    pub fn random(n_mutations: usize, rng: &mut impl Rng) -> Self {
        let root = n_mutations;
        let mut parent = vec![root; n_mutations + 1];
        parent[root] = root;

        let mut order: Vec<usize> = (0..n_mutations).collect();
        for i in (1..order.len()).rev() {
            let j = rng.random_range(0..=i);
            order.swap(i, j);
        }

        let mut placed: Vec<usize> = vec![root];
        for &node in &order {
            let p = placed[rng.random_range(0..placed.len())];
            parent[node] = p;
            placed.push(node);
        }

        MutationTree { n_mutations, parent }
    }
}

/// Fixed sequencing/genotyping error rates used in the SCITE likelihood model.
pub struct ErrorRates {
    /// alpha: P(observed = 1 | true = 0) — false positive rate.
    pub fp_rate: f64,
    /// beta: P(observed = 0 | true = 1) — false negative rate.
    pub fn_rate: f64,
}

/// Log-likelihood of `profile` (one read's observed calls, one per variant,
/// `None` = missing/uncovered and contributes nothing) given that
/// `ancestor_mask` bit `i` set means variant `i` is present under the
/// candidate tree attachment being scored.
pub fn attachment_log_likelihood(
    profile: &[Option<u8>],
    ancestor_mask: u64,
    rates: &ErrorRates,
) -> f64 {
    let mut ll = 0.0;
    for (i, call) in profile.iter().enumerate() {
        let Some(observed) = call else { continue };
        let expected_mutated = (ancestor_mask & (1u64 << i)) != 0;
        let p = match (*observed, expected_mutated) {
            (1, false) => rates.fp_rate,
            (0, false) => 1.0 - rates.fp_rate,
            (0, true) => rates.fn_rate,
            (1, true) => 1.0 - rates.fn_rate,
            _ => unreachable!("observed genotype must be 0 or 1"),
        };
        ll += p.ln();
    }
    ll
}

/// The tree node (mutation node or root) whose implied genotype best
/// explains `profile` under `rates`, together with that best log-likelihood.
pub fn best_attachment(
    profile: &[Option<u8>],
    tree: &MutationTree,
    rates: &ErrorRates,
) -> (usize, f64) {
    let mut best_node = tree.root();
    let mut best_ll = f64::NEG_INFINITY;
    for node in 0..=tree.n_mutations {
        let mask = tree.ancestor_mask(node);
        let ll = attachment_log_likelihood(profile, mask, rates);
        if ll > best_ll {
            best_ll = ll;
            best_node = node;
        }
    }
    (best_node, best_ll)
}

/// Total log-likelihood of `matrix` under `tree`: sum, over every read, of
/// that read's best-attachment log-likelihood.
pub fn tree_log_likelihood(matrix: &BinaryMatrix, tree: &MutationTree, rates: &ErrorRates) -> f64 {
    let n_reads = matrix.reads.len();
    (0..n_reads)
        .map(|r| {
            let profile: Vec<Option<u8>> = matrix.data.iter().map(|row| row[r]).collect();
            best_attachment(&profile, tree, rates).1
        })
        .sum()
}

fn subtree_all_carry(tree: &lineage::Tree, hap_matrix: &HaplotypeMatrix, node_id: usize, mutation: usize) -> bool {
    let mut stack = vec![node_id];
    while let Some(id) = stack.pop() {
        let node = &tree.nodes[id];
        if node.is_leaf {
            if hap_matrix.haplotypes[node.id].profile[mutation] != Some(1) {
                return false;
            }
        } else {
            stack.extend(node.children.iter().copied());
        }
    }
    true
}

/// Build an initial mutation tree from the Neighbor-Joining haplotype tree:
/// each mutation is ranked by the size (in attached reads) of the largest NJ
/// subtree entirely carrying it, then chained in descending order.
pub fn from_nj_tree(hap_matrix: &HaplotypeMatrix, nj_tree: &lineage::Tree) -> MutationTree {
    let n = hap_matrix.variants.len();
    let mut best_size = vec![0usize; n];

    for node in &nj_tree.nodes {
        let size = nj_tree.subtree_read_count(node.id);
        for m in 0..n {
            if size > best_size[m] && subtree_all_carry(nj_tree, hap_matrix, node.id, m) {
                best_size[m] = size;
            }
        }
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| best_size[b].cmp(&best_size[a]));

    let root = n;
    let mut parent = vec![root; n + 1];
    parent[root] = root;
    let mut prev = root;
    for &m in &order {
        parent[m] = prev;
        prev = m;
    }

    MutationTree { n_mutations: n, parent }
}

fn sample_two_distinct(n: usize, rng: &mut impl Rng) -> (usize, usize) {
    let first = rng.random_range(0..n);
    let mut second = rng.random_range(0..n);
    while second == first {
        second = rng.random_range(0..n);
    }
    (first, second)
}

/// Pick a uniformly random non-root node and reattach it to a uniformly
/// random node outside its own subtree (excluding itself), so the result is
/// always acyclic.
///
/// This proposal is symmetric: the number of valid reattachment targets for
/// the chosen node depends only on the size of its own subtree, which this
/// move does not change (moving a subtree doesn't change what's inside it,
/// only where it's attached) — so plain Metropolis-Hastings acceptance
/// (Task 9) needs no extra correction term for this move.
pub fn propose_prune_reattach(tree: &MutationTree, rng: &mut impl Rng) -> MutationTree {
    let n = tree.n_mutations;
    let node = rng.random_range(0..n);
    let forbidden = tree.descendant_mask(node) | (1u64 << node);
    let valid_targets: Vec<usize> = (0..=n).filter(|&t| forbidden & (1u64 << t) == 0).collect();
    let new_parent = valid_targets[rng.random_range(0..valid_targets.len())];
    tree.with_parent(node, new_parent)
}

/// Relabel two mutation nodes: whichever tree positions `first` and `second`
/// occupied, they now swap. Symmetric proposal (no Hastings correction),
/// ported from reference SCITE's `get_new_parent_vec_swap_fast`.
fn swap_labels(tree: &MutationTree, first: usize, second: usize) -> MutationTree {
    let n = tree.n_mutations;
    let mut parent = tree.parent.clone();

    for i in 0..n {
        if parent[i] == first && i != second {
            parent[i] = second;
        } else if parent[i] == second && i != first {
            parent[i] = first;
        }
    }
    parent.swap(first, second);
    if parent[first] == first {
        parent[first] = second;
    }
    if parent[second] == second {
        parent[second] = first;
    }

    MutationTree { n_mutations: n, parent }
}

pub fn propose_swap_labels(tree: &MutationTree, rng: &mut impl Rng) -> MutationTree {
    let (first, second) = sample_two_distinct(tree.n_mutations, rng);
    swap_labels(tree, first, second)
}

/// Swap the tree positions of nodes `a` and `b`. When they're in different
/// lineages this is a direct mutual reattachment. When one is an ancestor of
/// the other, the ancestor is instead reinserted under a member of the
/// descendant's (updated) subtree — `choice_index` selects which one, and
/// must be `< ` the number of such candidates (see the two unit tests above
/// for how to construct a case where that count is exactly 1). Returns the
/// Metropolis-Hastings neighborhood-correction factor to multiply into the
/// acceptance ratio (`1.0` for the different-lineages case).
///
/// Ported from reference SCITE's `propose_new_tree` movetype 3.
fn swap_subtrees(tree: &MutationTree, a: usize, b: usize, choice_index: usize) -> (MutationTree, f64) {
    let n = tree.n_mutations;
    let (node, next_node) = if (tree.ancestor_mask(b) & (1u64 << a)) != 0 { (b, a) } else { (a, b) };

    let mut parent = tree.parent.clone();

    if (tree.ancestor_mask(node) & (1u64 << next_node)) == 0 {
        parent[node] = tree.parent[next_node];
        parent[next_node] = tree.parent[node];
        return (MutationTree { n_mutations: n, parent }, 1.0);
    }

    parent[node] = tree.parent[next_node];
    let descendants: Vec<usize> = (0..n)
        .filter(|&i| (tree.ancestor_mask(i) & (1u64 << node)) != 0)
        .collect();
    parent[next_node] = descendants[choice_index];

    let proposal = MutationTree { n_mutations: n, parent };
    let next_descendants = (0..n)
        .filter(|&i| (proposal.ancestor_mask(i) & (1u64 << next_node)) != 0)
        .count();

    let nbh_correction = descendants.len() as f64 / next_descendants as f64;
    (proposal, nbh_correction)
}

pub fn propose_swap_subtrees(tree: &MutationTree, rng: &mut impl Rng) -> (MutationTree, f64) {
    let n = tree.n_mutations;
    let (a, b) = sample_two_distinct(n, rng);
    let (node, next_node) = if (tree.ancestor_mask(b) & (1u64 << a)) != 0 { (b, a) } else { (a, b) };
    if (tree.ancestor_mask(node) & (1u64 << next_node)) == 0 {
        return swap_subtrees(tree, node, next_node, 0); // choice_index unused in this branch
    }
    let descendants_count = (0..n).filter(|&i| (tree.ancestor_mask(i) & (1u64 << node)) != 0).count();
    let choice_index = rng.random_range(0..descendants_count);
    swap_subtrees(tree, node, next_node, choice_index)
}

/// Weighted move dispatcher, matching reference SCITE's default move-type
/// probabilities for the non-error-rate moves (`0.55` prune-and-reattach,
/// `0.4` swap-labels, `0.05` swap-subtrees), renormalized to sum to 1 since
/// this plan has no error-rate move at all.
pub fn propose_move(tree: &MutationTree, rng: &mut impl Rng) -> (MutationTree, f64) {
    let r: f64 = rng.random();
    if r < 0.55 {
        (propose_prune_reattach(tree, rng), 1.0)
    } else if r < 0.95 {
        (propose_swap_labels(tree, rng), 1.0)
    } else {
        propose_swap_subtrees(tree, rng)
    }
}

/// Run one MCMC chain of `n_iterations` moves, starting from `initial_tree`
/// if given, otherwise a random tree, and return the highest-likelihood tree
/// seen at any point in the chain (MCMC explores, it doesn't monotonically
/// improve, so the running best must be tracked separately from the current
/// state).
pub fn run_mcmc(
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
    n_iterations: usize,
    initial_tree: Option<&MutationTree>,
    rng: &mut impl Rng,
) -> (MutationTree, f64) {
    let n = matrix.variants.len();
    let mut current = initial_tree.cloned().unwrap_or_else(|| MutationTree::random(n, rng));
    let mut current_ll = tree_log_likelihood(matrix, &current, rates);
    let mut best = current.clone();
    let mut best_ll = current_ll;

    for _ in 0..n_iterations {
        let (proposal, nbh_correction) = propose_move(&current, rng);
        let proposal_ll = tree_log_likelihood(matrix, &proposal, rates);

        let acceptance = nbh_correction * (proposal_ll - current_ll).exp();
        if rng.random::<f64>() < acceptance {
            current = proposal;
            current_ll = proposal_ll;
            if current_ll > best_ll {
                best = current.clone();
                best_ll = current_ll;
            }
        }
    }

    (best, best_ll)
}

/// Run `n_chains` independent MCMC chains (each with its own derived seed)
/// and return the single best tree found across all of them. `initial_tree`
/// is used only for the first chain; the rest always start random.
pub fn run_mcmc_multichain(
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
    n_iterations: usize,
    n_chains: usize,
    initial_tree: Option<&MutationTree>,
    seed: u64,
) -> (MutationTree, f64) {
    let mut best: Option<(MutationTree, f64)> = None;
    for chain in 0..n_chains {
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(chain as u64));
        let chain_initial = if chain == 0 { initial_tree } else { None };
        let (tree, ll) = run_mcmc(matrix, rates, n_iterations, chain_initial, &mut rng);
        best = match best {
            Some((_, best_ll)) if best_ll >= ll => best,
            _ => Some((tree, ll)),
        };
    }
    best.expect("n_chains must be >= 1")
}

/// The fully-imputed, denoised genotype matrix implied by attaching every
/// read to its single best-fitting node in `tree`.
#[derive(Debug)]
pub struct CleanedMatrix {
    pub variants: Vec<String>,
    pub reads: Vec<String>,
    /// `data[variant_idx][read_idx]`, always `0` or `1` (no missing values).
    pub data: Vec<Vec<u8>>,
    /// `attachment[read_idx]` = the tree node id that read attached to.
    pub attachment: Vec<usize>,
}

pub fn attach_all_reads(matrix: &BinaryMatrix, tree: &MutationTree, rates: &ErrorRates) -> CleanedMatrix {
    let n_variants = matrix.variants.len();
    let n_reads = matrix.reads.len();
    let mut data = vec![vec![0u8; n_reads]; n_variants];
    let mut attachment = vec![0usize; n_reads];

    for r in 0..n_reads {
        let profile: Vec<Option<u8>> = matrix.data.iter().map(|row| row[r]).collect();
        let (node, _ll) = best_attachment(&profile, tree, rates);
        let mask = tree.ancestor_mask(node);
        attachment[r] = node;
        for v in 0..n_variants {
            data[v][r] = u8::from((mask & (1u64 << v)) != 0);
        }
    }

    CleanedMatrix {
        variants: matrix.variants.clone(),
        reads: matrix.reads.clone(),
        data,
        attachment,
    }
}

pub fn write_cleaned_matrix(matrix: &CleanedMatrix, path: &str) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    write!(w, "variant")?;
    for read in &matrix.reads {
        write!(w, ",{read}")?;
    }
    writeln!(w)?;
    for (v_idx, variant) in matrix.variants.iter().enumerate() {
        write!(w, "{variant}")?;
        for r_idx in 0..matrix.reads.len() {
            write!(w, ",{}", matrix.data[v_idx][r_idx])?;
        }
        writeln!(w)?;
    }
    Ok(())
}

pub fn write_variant_cooccurrence(matrix: &CleanedMatrix, path: &str) -> Result<()> {
    let n = matrix.variants.len();
    let mut counts = vec![vec![0usize; n]; n];
    for i in 0..n {
        for j in 0..n {
            counts[i][j] = (0..matrix.reads.len())
                .filter(|&r| matrix.data[i][r] == 1 && matrix.data[j][r] == 1)
                .count();
        }
    }

    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    write!(w, "variant")?;
    for v in &matrix.variants {
        write!(w, "\t{v}")?;
    }
    writeln!(w)?;
    for i in 0..n {
        write!(w, "{}", matrix.variants[i])?;
        for j in 0..n {
            write!(w, "\t{}", counts[i][j])?;
        }
        writeln!(w)?;
    }
    Ok(())
}

pub fn write_molecule_summary(matrix: &CleanedMatrix, path: &str) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(w, "read_name\tn_variants\tvariants")?;
    for (r_idx, read) in matrix.reads.iter().enumerate() {
        let present: Vec<&str> = matrix
            .variants
            .iter()
            .enumerate()
            .filter(|(v_idx, _)| matrix.data[*v_idx][r_idx] == 1)
            .map(|(_, name)| name.as_str())
            .collect();
        writeln!(w, "{read}\t{}\t{}", present.len(), present.join(","))?;
    }
    Ok(())
}

pub fn write_mutation_tree(tree: &MutationTree, matrix: &CleanedMatrix, path: &str) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(w, "node_id\tvariant\tparent_id\tparent_variant\tn_reads_attached")?;
    let node_name = |id: usize| -> String {
        if id == tree.root() { "ROOT".to_string() } else { matrix.variants[id].clone() }
    };
    for node in 0..=tree.n_mutations {
        let parent = tree.parent[node];
        let n_reads = matrix.attachment.iter().filter(|&&a| a == node).count();
        writeln!(w, "{node}\t{}\t{parent}\t{}\t{n_reads}", node_name(node), node_name(parent))?;
    }
    Ok(())
}

/// Recursively emit the Newick token for the subtree at `node`. `edge_muts` are
/// the mutations accumulated (ancestral→derived) on the branch leading into the
/// node that will actually be emitted here — pass-through mutation nodes fold
/// their mutation into this list rather than becoming their own node.
fn emit_lineage_node(
    node: usize,
    edge_muts: &[String],
    tree: &MutationTree,
    children: &[Vec<usize>],
    haps_by_node: &[Vec<(String, usize)>],
    variants: &[String],
) -> String {
    let is_root = node == tree.root();
    let child_nodes = &children[node];
    let haps = &haps_by_node[node];
    let n_effective = child_nodes.len() + haps.len();

    // Collapse: a non-root node with a single effective child is dissolved.
    if !is_root && n_effective == 1 {
        if child_nodes.len() == 1 {
            // Single child-lineage: recurse, appending this node's edge mutations.
            let c = child_nodes[0];
            let mut muts = edge_muts.to_vec();
            muts.push(variants[c].clone());
            return emit_lineage_node(c, &muts, tree, children, haps_by_node, variants);
        } else {
            // Single haplotype tip: the mutation(s) ride on that tip's branch.
            let (hid, count) = &haps[0];
            return format!(
                "{hid}_n{count}:{}[&&NHX:mutation={}:reads={count}]",
                edge_muts.len(),
                edge_muts.join(",")
            );
        }
    }

    // Otherwise this node is emitted (branch point, root, or dead-end).
    let mut tokens: Vec<String> = Vec::new();
    for &c in child_nodes {
        let muts = vec![variants[c].clone()];
        tokens.push(emit_lineage_node(c, &muts, tree, children, haps_by_node, variants));
    }
    for (hid, count) in haps {
        // Haplotype attached directly to a surviving node: zero-length tip, reads only.
        tokens.push(format!("{hid}_n{count}:0[&&NHX:reads={count}]"));
    }

    if is_root {
        if tokens.is_empty() {
            "ROOT".to_string()
        } else {
            format!("({})ROOT", tokens.join(","))
        }
    } else {
        let blen = edge_muts.len();
        let muts = edge_muts.join(",");
        if tokens.is_empty() {
            // Dead-end ancestral node (no descendants, no reads) — keep the mutation.
            format!("anc{node}:{blen}[&&NHX:mutation={muts}]")
        } else {
            format!("({})anc{node}:{blen}[&&NHX:mutation={muts}]", tokens.join(","))
        }
    }
}

/// Write the read lineage tree (Newick + NHX branch mutations, tsinfer-style)
/// to `path`. Tips are haplotypes (labeled `H<id>_n<reads>`, with a machine-
/// readable `reads` NHX field); branches carry the mutations acquired, and
/// pass-through mutation nodes are collapsed so mutations ride on a single
/// branch (comma-joined, ancestral→derived).
pub fn write_read_lineage_newick(
    tree: &MutationTree,
    hap_matrix: &HaplotypeMatrix,
    rates: &ErrorRates,
    path: &str,
) -> Result<()> {
    let children = tree.children_of();
    let mut haps_by_node: Vec<Vec<(String, usize)>> = vec![Vec::new(); tree.n_mutations + 1];
    for hap in &hap_matrix.haplotypes {
        let (node, _ll) = best_attachment(&hap.profile, tree, rates);
        haps_by_node[node].push((hap.id.clone(), hap.count));
    }

    let body = emit_lineage_node(
        tree.root(),
        &[],
        tree,
        &children,
        &haps_by_node,
        &hap_matrix.variants,
    );

    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(w, "{body};")?;
    Ok(())
}

/// Run the MCMC mutation-tree search (starting from an NJ-tree-derived guess
/// when one can be built), attach every read to its best-fitting node, and
/// write all four SCITE output files with the given `output_prefix`.
pub fn run_scite_pipeline(
    binary: &BinaryMatrix,
    hap_matrix: &HaplotypeMatrix,
    fp_rate: f64,
    fn_rate: f64,
    n_iterations: usize,
    n_chains: usize,
    seed: u64,
    min_reads: usize,
    output_prefix: &str,
) -> Result<()> {
    if binary.variants.len() < 2 {
        anyhow::bail!(
            "SCITE requires at least 2 variants (got {}); \
             relax --min-hf / --min-presence / --min-absence thresholds.",
            binary.variants.len()
        );
    }
    if binary.variants.len() > 63 {
        anyhow::bail!(
            "SCITE mutation-tree search supports at most 63 variants (got {}); \
             tighten the lineage HF/prevalence thresholds.",
            binary.variants.len()
        );
    }

    let dist = lineage::hamming_distance_matrix(hap_matrix);
    let initial_tree = lineage::neighbor_joining(&dist, hap_matrix)
        .ok()
        .map(|nj_tree| from_nj_tree(hap_matrix, &nj_tree));
    eprintln!(
        "[SCITE] Initial tree: {}",
        if initial_tree.is_some() { "derived from Neighbor-Joining haplotype tree" } else { "random (NJ tree unavailable)" }
    );

    eprintln!(
        "[SCITE] {} variants, {} reads. Running MCMC: {n_chains} chains x {n_iterations} iterations",
        binary.variants.len(),
        binary.reads.len()
    );

    let rates = ErrorRates { fp_rate, fn_rate };
    let (tree, ll) = run_mcmc_multichain(binary, &rates, n_iterations, n_chains, initial_tree.as_ref(), seed);
    eprintln!("[SCITE] Best tree log-likelihood: {ll:.3}");

    let cleaned = attach_all_reads(binary, &tree, &rates);

    write_cleaned_matrix(&cleaned, &format!("{output_prefix}.cleaned_matrix.csv"))?;
    write_variant_cooccurrence(&cleaned, &format!("{output_prefix}.variant_cooccurrence.tsv"))?;
    write_molecule_summary(&cleaned, &format!("{output_prefix}.molecule_summary.tsv"))?;
    write_mutation_tree(&tree, &cleaned, &format!("{output_prefix}.mutation_tree.tsv"))?;

    // Deduplicate the SCITE-cleaned reads into haplotypes for the lineage tree.
    // Keep the tree's variant set/order (no prevalence re-filter) so the mutation
    // indices still line up with the tree nodes.
    let cleaned_binary = BinaryMatrix {
        variants: cleaned.variants.clone(),
        reads: cleaned.reads.clone(),
        data: cleaned
            .data
            .iter()
            .map(|row| row.iter().map(|&b| Some(b)).collect())
            .collect(),
    };
    let cleaned_hap_matrix = lineage::deduplicate(&cleaned_binary, min_reads);
    write_read_lineage_newick(&tree, &cleaned_hap_matrix, &rates, &format!("{output_prefix}.read_lineage.nwk"))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chain_tree() -> MutationTree {
        // 0 -> 1 -> 2 -> root(3): node0's parent is 1, node1's parent is 2,
        // node2's parent is root, root's parent is itself.
        MutationTree { n_mutations: 3, parent: vec![1, 2, 3, 3] }
    }

    fn assert_valid_tree(tree: &MutationTree) {
        let n = tree.n_mutations;
        assert_eq!(tree.parent.len(), n + 1);
        assert_eq!(tree.parent[tree.root()], tree.root());
        for start in 0..=n {
            let mut cur = start;
            let mut steps = 0;
            while cur != tree.root() {
                cur = tree.parent[cur];
                steps += 1;
                assert!(steps <= n, "cycle detected starting from node {start}");
            }
        }
    }

    #[test]
    fn ancestor_mask_includes_self_and_all_ancestors_up_to_root() {
        let tree = chain_tree();
        assert_eq!(tree.ancestor_mask(0), 0b111);
        assert_eq!(tree.ancestor_mask(1), 0b110);
        assert_eq!(tree.ancestor_mask(2), 0b100);
        assert_eq!(tree.ancestor_mask(3), 0);
    }

    #[test]
    fn descendant_mask_excludes_self_and_includes_transitive_children() {
        let tree = chain_tree();
        assert_eq!(tree.descendant_mask(2), 0b011);
        assert_eq!(tree.descendant_mask(1), 0b001);
        assert_eq!(tree.descendant_mask(0), 0);
        assert_eq!(tree.descendant_mask(3), 0b111);
    }

    #[test]
    fn with_parent_changes_only_the_targeted_node() {
        let tree = chain_tree();
        let moved = tree.with_parent(0, 3);
        assert_eq!(moved.parent, vec![3, 2, 3, 3]);
        assert_eq!(tree.parent, vec![1, 2, 3, 3]);
    }

    #[test]
    fn attachment_log_likelihood_matches_hand_computed_value() {
        let profile = vec![Some(1), Some(0), None];
        let ancestor_mask = 0b011;
        let rates = ErrorRates { fp_rate: 0.1, fn_rate: 0.2 };

        let expected = (1.0 - rates.fn_rate).ln() + rates.fn_rate.ln();
        let actual = attachment_log_likelihood(&profile, ancestor_mask, &rates);
        assert!((actual - expected).abs() < 1e-12);
    }

    #[test]
    fn attachment_log_likelihood_scores_true_negative_and_false_positive() {
        let profile = vec![Some(0), Some(1)];
        let ancestor_mask = 0;
        let rates = ErrorRates { fp_rate: 0.1, fn_rate: 0.2 };

        let expected = (1.0 - rates.fp_rate).ln() + rates.fp_rate.ln();
        let actual = attachment_log_likelihood(&profile, ancestor_mask, &rates);
        assert!((actual - expected).abs() < 1e-12);
    }

    #[test]
    fn best_attachment_picks_the_node_matching_observed_mutations() {
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 2, 2] };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.3 };
        let profile = vec![Some(1), Some(0)];

        let expected_ll = (1.0 - rates.fn_rate).ln() + (1.0 - rates.fp_rate).ln();
        let (node, ll) = best_attachment(&profile, &tree, &rates);
        assert_eq!(node, 0);
        assert!((ll - expected_ll).abs() < 1e-12);
    }

    #[test]
    fn tree_log_likelihood_sums_best_attachment_over_all_reads() {
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 2, 2] };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.3 };
        let matrix = BinaryMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            reads: vec!["r1".to_string(), "r2".to_string()],
            data: vec![vec![Some(1), Some(0)], vec![Some(0), Some(1)]],
        };

        let per_read_ll = (1.0 - rates.fn_rate).ln() + (1.0 - rates.fp_rate).ln();
        let expected = 2.0 * per_read_ll;
        let actual = tree_log_likelihood(&matrix, &tree, &rates);
        assert!((actual - expected).abs() < 1e-12);
    }

    #[test]
    fn random_tree_is_always_valid() {
        let mut rng = StdRng::seed_from_u64(1);
        for n in [1, 2, 5, 10] {
            for _ in 0..20 {
                let tree = MutationTree::random(n, &mut rng);
                assert_valid_tree(&tree);
            }
        }
    }

    #[test]
    fn from_nj_tree_orders_mutations_by_largest_fully_carrying_subtree() {
        use lineage::{Haplotype, HaplotypeMatrix, Node, Tree};

        // 3 leaf haplotypes, all children of one root:
        //   H0 (5 reads): neither A nor B
        //   H1 (3 reads): A only
        //   H2 (2 reads): A and B
        // A's largest fully-carrying subtree is {H1} (3 reads) since H0 lacks A.
        // B's largest fully-carrying subtree is {H2} (2 reads) since H0, H1 lack B.
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            haplotypes: vec![
                Haplotype { id: "H0".to_string(), profile: vec![Some(0), Some(0)], reads: vec![], count: 5 },
                Haplotype { id: "H1".to_string(), profile: vec![Some(1), Some(0)], reads: vec![], count: 3 },
                Haplotype { id: "H2".to_string(), profile: vec![Some(1), Some(1)], reads: vec![], count: 2 },
            ],
        };
        let nj_tree = Tree {
            nodes: vec![
                Node { id: 0, label: "H0".to_string(), is_leaf: true, children: vec![], parent: Some(3), branch_length: 0.0, read_count: 5 },
                Node { id: 1, label: "H1".to_string(), is_leaf: true, children: vec![], parent: Some(3), branch_length: 0.0, read_count: 3 },
                Node { id: 2, label: "H2".to_string(), is_leaf: true, children: vec![], parent: Some(3), branch_length: 0.0, read_count: 2 },
                Node { id: 3, label: "ROOT".to_string(), is_leaf: false, children: vec![0, 1, 2], parent: None, branch_length: 0.0, read_count: 0 },
            ],
            root: 3,
        };

        let tree = from_nj_tree(&hap_matrix, &nj_tree);
        // A ranked ahead of B -> A's parent is root, B's parent is A.
        assert_eq!(tree.parent, vec![2, 0, 2]);
    }

    #[test]
    fn propose_prune_reattach_always_produces_a_valid_tree() {
        let mut rng = StdRng::seed_from_u64(2);
        let tree = MutationTree { n_mutations: 4, parent: vec![2, 2, 4, 4, 4] };
        for _ in 0..200 {
            let proposal = propose_prune_reattach(&tree, &mut rng);
            assert_valid_tree(&proposal);
            assert_eq!(proposal.n_mutations, tree.n_mutations);
        }
    }

    #[test]
    fn propose_prune_reattach_never_attaches_a_node_under_its_own_descendant() {
        let mut rng = StdRng::seed_from_u64(3);
        let tree = MutationTree { n_mutations: 4, parent: vec![2, 2, 4, 4, 4] };
        for _ in 0..200 {
            let proposal = propose_prune_reattach(&tree, &mut rng);
            assert_ne!(proposal.parent[2], 0);
            assert_ne!(proposal.parent[2], 1);
        }
    }

    #[test]
    fn swap_labels_exchanges_two_nodes_positions_in_the_chain() {
        let tree = chain_tree(); // parent = [1, 2, 3, 3]: chain 0->1->2->root(3)
        let swapped = swap_labels(&tree, 0, 2);
        // Relabeling 0<->2 turns "0->1->2->root" into "2->1->0->root".
        assert_eq!(swapped.parent, vec![3, 0, 1, 3]);
    }

    #[test]
    fn propose_swap_labels_always_produces_a_valid_tree() {
        let mut rng = StdRng::seed_from_u64(4);
        let tree = MutationTree { n_mutations: 4, parent: vec![2, 2, 4, 4, 4] };
        for _ in 0..200 {
            let proposal = propose_swap_labels(&tree, &mut rng);
            assert_valid_tree(&proposal);
        }
    }

    #[test]
    fn swap_subtrees_different_lineages_swaps_parents_directly() {
        let tree = MutationTree { n_mutations: 4, parent: vec![2, 2, 4, 4, 4] };
        let (proposal, correction) = swap_subtrees(&tree, 0, 3, 0);
        assert_eq!(proposal.parent, vec![4, 2, 4, 2, 4]);
        assert_eq!(correction, 1.0);
    }

    #[test]
    fn swap_subtrees_same_lineage_reinserts_ancestor_with_hastings_correction() {
        let tree = MutationTree { n_mutations: 4, parent: vec![2, 2, 4, 4, 4] };
        // node=1, next_node=2: node1's only current "descendant" (self-inclusive)
        // is itself, so this is fully deterministic regardless of choice_index.
        let (proposal, correction) = swap_subtrees(&tree, 1, 2, 0);
        assert_eq!(proposal.parent, vec![2, 4, 1, 4, 4]);
        assert!((correction - 0.5).abs() < 1e-12);
    }

    #[test]
    fn swap_labels_star_topology_does_not_double_fire_self_loop_fixup() {
        // Star tree: all leaves are direct children of root(4).
        // Swapping any two leaves must not break the star structure or
        // double-apply the self-loop fixup. This is the topology the
        // chain_tree test misses.
        let tree = MutationTree { n_mutations: 4, parent: vec![4, 4, 4, 4, 4] };
        let swapped = swap_labels(&tree, 0, 3);
        assert_valid_tree(&swapped);
        // Labels 0 and 3 are leaves with no children — the parent array is
        // unchanged for all other nodes; only positions 0 and 3 swap.
        assert_eq!(swapped.parent, vec![4, 4, 4, 4, 4]);
    }

    #[test]
    fn propose_move_always_produces_a_valid_tree_with_a_positive_finite_correction() {
        let mut rng = StdRng::seed_from_u64(5);
        let tree = MutationTree { n_mutations: 4, parent: vec![2, 2, 4, 4, 4] };
        for _ in 0..500 {
            let (proposal, correction) = propose_move(&tree, &mut rng);
            assert_valid_tree(&proposal);
            assert!(correction.is_finite() && correction > 0.0);
        }
    }

    #[test]
    fn run_mcmc_recovers_a_true_ancestor_relationship() {
        let mut data_a = Vec::new();
        let mut data_b = Vec::new();
        for _ in 0..6 {
            data_a.push(Some(1));
            data_b.push(Some(1));
        }
        for _ in 0..3 {
            data_a.push(Some(1));
            data_b.push(Some(0));
        }
        data_a.push(Some(0)); // simulated false negative
        data_b.push(Some(1));

        let n_reads = data_a.len();
        let matrix = BinaryMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            reads: (0..n_reads).map(|i| format!("r{i}")).collect(),
            data: vec![data_a, data_b],
        };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.1 };

        let mut rng = StdRng::seed_from_u64(42);
        let (tree, _ll) = run_mcmc(&matrix, &rates, 3000, None, &mut rng);

        assert_eq!(tree.parent[1], 0);
    }

    #[test]
    fn run_mcmc_multichain_finds_the_true_tree_with_fewer_iterations_per_chain() {
        let mut data_a = Vec::new();
        let mut data_b = Vec::new();
        for _ in 0..6 {
            data_a.push(Some(1));
            data_b.push(Some(1));
        }
        for _ in 0..3 {
            data_a.push(Some(1));
            data_b.push(Some(0));
        }
        data_a.push(Some(0));
        data_b.push(Some(1));

        let n_reads = data_a.len();
        let matrix = BinaryMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            reads: (0..n_reads).map(|i| format!("r{i}")).collect(),
            data: vec![data_a, data_b],
        };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.1 };

        let (tree, _ll) = run_mcmc_multichain(&matrix, &rates, 1000, 4, None, 7);
        assert_eq!(tree.parent[1], 0);
    }

    #[test]
    fn attach_all_reads_imputes_missing_and_corrects_false_negatives() {
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] }; // A ancestor of B
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.1 };
        let matrix = BinaryMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            reads: vec!["r_full".to_string(), "r_missingB".to_string(), "r_errorA".to_string()],
            data: vec![
                vec![Some(1), Some(1), Some(0)],
                vec![Some(1), None, Some(1)],
            ],
        };

        let cleaned = attach_all_reads(&matrix, &tree, &rates);

        assert_eq!(cleaned.variants, vec!["A", "B"]);
        assert_eq!(cleaned.reads, vec!["r_full", "r_missingB", "r_errorA"]);
        assert_eq!(cleaned.data[0], vec![1, 1, 1]); // A: unchanged, unchanged, false negative flipped
        assert_eq!(cleaned.data[1], vec![1, 0, 1]); // B: unchanged, imputed absent, unchanged
        assert_eq!(cleaned.attachment, vec![1, 0, 1]);
    }

    fn small_cleaned_matrix() -> CleanedMatrix {
        CleanedMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            reads: vec!["r1".to_string(), "r2".to_string()],
            data: vec![vec![1, 0], vec![1, 1]],
            attachment: vec![1, 0],
        }
    }

    #[test]
    fn write_cleaned_matrix_writes_matrix_csv_shape() {
        let matrix = small_cleaned_matrix();
        let path = std::env::temp_dir().join("himito_test_write_cleaned_matrix.csv");
        write_cleaned_matrix(&matrix, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(content, "variant,r1,r2\nA,1,0\nB,1,1\n");
    }

    #[test]
    fn write_variant_cooccurrence_counts_reads_with_both_variants() {
        let matrix = small_cleaned_matrix();
        let path = std::env::temp_dir().join("himito_test_write_cooccurrence.tsv");
        write_variant_cooccurrence(&matrix, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(content, "variant\tA\tB\nA\t1\t1\nB\t1\t2\n");
    }

    #[test]
    fn write_molecule_summary_lists_variant_count_and_names_per_read() {
        let matrix = small_cleaned_matrix();
        let path = std::env::temp_dir().join("himito_test_write_molecule_summary.tsv");
        write_molecule_summary(&matrix, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(content, "read_name\tn_variants\tvariants\nr1\t2\tA,B\nr2\t1\tB\n");
    }

    #[test]
    fn write_mutation_tree_lists_every_node_with_parent_and_attached_read_count() {
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        let matrix = small_cleaned_matrix();
        let path = std::env::temp_dir().join("himito_test_write_mutation_tree.tsv");
        write_mutation_tree(&tree, &matrix, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(
            content,
            "node_id\tvariant\tparent_id\tparent_variant\tn_reads_attached\n\
             0\tA\t2\tROOT\t1\n\
             1\tB\t0\tA\t1\n\
             2\tROOT\t2\tROOT\t0\n"
        );
    }

    #[test]
    fn run_scite_pipeline_writes_all_four_output_files() {
        let matrix = BinaryMatrix {
            variants: vec!["m.A".to_string(), "m.B".to_string()],
            reads: (0..10).map(|i| format!("r{i}")).collect(),
            data: vec![
                vec![Some(1); 9].into_iter().chain([Some(0)]).collect(),
                vec![Some(1); 6].into_iter().chain(vec![Some(0); 3]).chain([Some(1)]).collect(),
            ],
        };
        let hap_matrix = lineage::deduplicate(&matrix, 1);

        let prefix = std::env::temp_dir()
            .join("himito_test_scite_pipeline_out")
            .to_str()
            .unwrap()
            .to_string();

        run_scite_pipeline(&matrix, &hap_matrix, 0.01, 0.1, 500, 2, 123, 1, &prefix).unwrap();

        for suffix in [
            ".cleaned_matrix.csv",
            ".variant_cooccurrence.tsv",
            ".molecule_summary.tsv",
            ".mutation_tree.tsv",
            ".read_lineage.nwk",
        ] {
            let path = format!("{prefix}{suffix}");
            assert!(std::path::Path::new(&path).exists(), "expected output file {path} to exist");
            std::fs::remove_file(&path).ok();
        }
    }

    #[test]
    fn write_read_lineage_newick_collapses_and_annotates_reads_and_mutations() {
        use lineage::{Haplotype, HaplotypeMatrix};

        // root(4) -> A(0) -> B(1) -> {C(2), D(3)}
        //   A: no haplotypes, single child B      -> collapses onto B's branch
        //   B: branch point (C, D, and hap H3)    -> ancestral node anc1
        //   C: single haplotype H0, no children   -> collapses onto H0's tip branch
        //   D: two haplotypes H1, H2              -> ancestral node anc3
        //   H3 attaches to internal node B        -> length-0 tip, reads only
        let tree = MutationTree { n_mutations: 4, parent: vec![4, 0, 1, 1, 4] };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.05 };

        let hap_matrix = HaplotypeMatrix {
            variants: vec!["mA".into(), "mB".into(), "mC".into(), "mD".into()],
            haplotypes: vec![
                Haplotype { id: "H0".into(), profile: vec![Some(1), Some(1), Some(1), Some(0)], reads: vec![], count: 9 },
                Haplotype { id: "H1".into(), profile: vec![Some(1), Some(1), Some(0), Some(1)], reads: vec![], count: 7 },
                Haplotype { id: "H2".into(), profile: vec![Some(1), Some(1), None, Some(1)], reads: vec![], count: 4 },
                Haplotype { id: "H3".into(), profile: vec![Some(1), Some(1), Some(0), Some(0)], reads: vec![], count: 2 },
            ],
        };

        let path = std::env::temp_dir().join("himito_test_read_lineage.nwk");
        write_read_lineage_newick(&tree, &hap_matrix, &rates, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(
            content,
            "((H0_n9:1[&&NHX:mutation=mC:reads=9],\
             (H1_n7:0[&&NHX:reads=7],H2_n4:0[&&NHX:reads=4])anc3:1[&&NHX:mutation=mD],\
             H3_n2:0[&&NHX:reads=2])anc1:2[&&NHX:mutation=mA,mB])ROOT;\n"
        );
    }
}
