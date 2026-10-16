use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::fs::File;
use anyhow::{Context, Result};
use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;
use log::{info, warn};
use statrs::function::gamma::ln_gamma;
use adjustp::{adjust, Procedure};

use crate::lineage::{self, BinaryMatrix, HaplotypeMatrix};
use crate::rootblock;

/// A growable bitset over mutation/variant ids, backed by `Vec<u64>` words.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Mask {
    words: Vec<u64>,
}

impl Mask {
    /// An empty mask with no bits set.
    pub fn new() -> Self {
        Mask { words: Vec::new() }
    }

    /// Set bit `i`, growing the backing storage as needed.
    #[inline]
    pub fn set(&mut self, i: usize) {
        let word = i / 64;
        if word >= self.words.len() {
            self.words.resize(word + 1, 0);
        }
        self.words[word] |= 1u64 << (i % 64);
    }

    /// Test whether bit `i` is set.
    #[inline]
    pub fn get(&self, i: usize) -> bool {
        let word = i / 64;
        word < self.words.len() && (self.words[word] & (1u64 << (i % 64))) != 0
    }

    /// Number of set bits (mutations carried).
    #[inline]
    pub fn count_ones(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }
}

/// Log-likelihood differences below this are treated as exact ties.
const LL_TIE_EPS: f64 = 1e-9;

/// A SCITE-style mutation tree: mutation nodes `0..n_mutations` (indexed like the matrix variants)
/// plus the root `n_mutations`, which is its own parent.
#[derive(Debug, Clone)]
pub struct MutationTree {
    pub n_mutations: usize,
    pub parent: Vec<usize>,
}

impl MutationTree {
    pub fn root(&self) -> usize {
        self.n_mutations
    }

    /// Bitmask of mutation ids on the path from `node` (inclusive) to the root.
    pub fn ancestor_mask(&self, node: usize) -> Mask {
        let mut mask = Mask::new();
        let mut cur = node;
        while cur != self.root() {
            mask.set(cur);
            cur = self.parent[cur];
        }
        mask
    }

    /// Bitmask of all proper descendants of `node` (not including itself).
    fn descendant_mask(&self, node: usize) -> Mask {
        let children = self.children_of();
        let mut mask = Mask::new();
        let mut stack = vec![node];
        while let Some(cur) = stack.pop() {
            for &c in &children[cur] {
                if c < self.n_mutations {
                    mask.set(c);
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

    /// Returns a new tree identical to `self` except `node`'s parent is `new_parent`.
    fn with_parent(&self, node: usize, new_parent: usize) -> MutationTree {
        let mut parent = self.parent.clone();
        parent[node] = new_parent;
        MutationTree { n_mutations: self.n_mutations, parent }
    }

    /// A random tree: mutations, in random order, each attach to a uniformly random already-placed node.
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
    /// Per-variant `(fp, fn)`; columns past its end (root-block genotypes) use the global rates.
    pub per_variant: Vec<(f64, f64)>,
}

/// Error rates are clamped into `[RATE_FLOOR, 1 - RATE_FLOOR]` before any logarithm is taken.
const RATE_FLOOR: f64 = 1e-12;

impl ErrorRates {
    /// The clamped `(fp, fn)` for variant `v`.
    fn clamped_for(&self, v: usize) -> (f64, f64) {
        let (fp, fnr) = self.per_variant.get(v).copied().unwrap_or((self.fp_rate, self.fn_rate));
        (fp.clamp(RATE_FLOOR, 1.0 - RATE_FLOOR), fnr.clamp(RATE_FLOOR, 1.0 - RATE_FLOOR))
    }
}

/// Log-likelihood of `profile`; bit `i` of `ancestor_mask` means variant `i` is present.
pub fn attachment_log_likelihood(
    profile: &[Option<u8>],
    ancestor_mask: &Mask,
    rates: &ErrorRates,
) -> f64 {
    let mut ll = 0.0;
    for (i, call) in profile.iter().enumerate() {
        let Some(observed) = call else { continue };
        let (fp, fn_rate) = rates.clamped_for(i);
        let expected_mutated = ancestor_mask.get(i);
        let p = match (*observed, expected_mutated) {
            (1, false) => fp,
            (0, false) => 1.0 - fp,
            (0, true) => fn_rate,
            (1, true) => 1.0 - fn_rate,
            _ => unreachable!("observed genotype must be 0 or 1"),
        };
        ll += p.ln();
    }
    ll
}

/// The best `(node, ll)` among `(node, ll, n_mutations)` candidates; near-ties go to fewer mutations.
#[inline]
fn pick_best_node(candidates: impl Iterator<Item = (usize, f64, usize)>, root: usize) -> (usize, f64) {
    let (mut best_node, mut best_ll, mut best_muts) = (root, f64::NEG_INFINITY, usize::MAX);
    for (node, ll, n_muts) in candidates {
        let strictly_better = ll > best_ll + LL_TIE_EPS;
        let tied_and_simpler = (ll - best_ll).abs() <= LL_TIE_EPS && n_muts < best_muts;
        if strictly_better || tied_and_simpler {
            best_ll = best_ll.max(ll);
            best_node = node;
            best_muts = n_muts;
        }
    }
    (best_node, best_ll)
}

/// The tree node whose implied genotype best explains `profile` under `rates`.
/// Also the reference the tests check [`AttachmentScorer::best_attachment`] against.
pub fn best_attachment(
    profile: &[Option<u8>],
    tree: &MutationTree,
    rates: &ErrorRates,
) -> (usize, f64) {
    let candidates = (0..=tree.n_mutations).map(|node| {
        let mask = tree.ancestor_mask(node);
        (node, attachment_log_likelihood(profile, &mask, rates), mask.count_ones())
    });
    pick_best_node(candidates, tree.root())
}

/// Node visit order and mutation depth for one tree, computed once and reused across every read.
pub struct TreeLayout {
    preorder: Vec<usize>,
    depth: Vec<usize>,
}

impl TreeLayout {
    pub fn new(tree: &MutationTree) -> Self {
        let children = tree.children_of();
        let root = tree.root();
        let mut preorder = Vec::with_capacity(tree.n_mutations + 1);
        let mut depth = vec![0usize; tree.n_mutations + 1];
        let mut stack = vec![root];
        while let Some(node) = stack.pop() {
            preorder.push(node);
            for &c in &children[node] {
                depth[c] = depth[node] + 1;
                stack.push(c);
            }
        }
        TreeLayout { preorder, depth }
    }
}

/// Per-cell call codes: no call, an explicit reference call, an explicit alternate call.
const CODE_UNCOVERED: u8 = 0;
const CODE_REF: u8 = 1;
const CODE_ALT: u8 = 2;

/// Per-read log-likelihood terms, precomputed once per matrix so scoring a tree is one pass per read.
pub struct AttachmentScorer {
    /// `base[r]` = log P(read r's calls | genotype is all-reference).
    base: Vec<f64>,
    /// `code[r * n_variants + m]` — one of the three `CODE_*` values.
    code: Vec<u8>,
    /// Row stride of `code`.
    n_variants: usize,
    /// `delta[variant][code]` is log P(call | variant mutated) − log P(call | variant not mutated).
    delta: Vec<[f64; 3]>,
}

impl AttachmentScorer {
    pub fn new(matrix: &BinaryMatrix, rates: &ErrorRates) -> Self {
        let n_variants = matrix.variants.len();
        let n_reads = matrix.reads.len();
        let mut base = vec![0.0f64; n_reads];
        let mut code = vec![CODE_UNCOVERED; n_reads * n_variants];
        // delta[_][CODE_UNCOVERED] stays 0: an uncovered cell costs nothing at any node.
        let mut delta = vec![[0.0f64; 3]; n_variants];

        for m in 0..n_variants {
            let (fp, fnr) = rates.clamped_for(m);
            let ln_fp = fp.ln();
            let ln_1mfp = (1.0 - fp).ln();
            let ln_fn = fnr.ln();
            let ln_1mfn = (1.0 - fnr).ln();
            delta[m][CODE_REF as usize] = ln_fn - ln_1mfp;
            delta[m][CODE_ALT as usize] = ln_1mfn - ln_fp;

            for r in 0..n_reads {
                match matrix.data[m][r] {
                    Some(1) => {
                        base[r] += ln_fp;
                        code[r * n_variants + m] = CODE_ALT;
                    }
                    Some(0) => {
                        base[r] += ln_1mfp;
                        code[r * n_variants + m] = CODE_REF;
                    }
                    _ => {}
                }
            }
        }

        AttachmentScorer { base, code, n_variants, delta }
    }

    /// How many reads this scorer was built over.
    fn n_reads(&self) -> usize {
        self.base.len()
    }

    /// The best-fitting node for read `read` and its log-likelihood.
    pub fn best_attachment(
        &self,
        read: usize,
        tree: &MutationTree,
        layout: &TreeLayout,
        scratch: &mut Vec<f64>,
    ) -> (usize, f64) {
        let root = tree.root();
        scratch.clear();
        scratch.resize(tree.n_mutations + 1, 0.0);
        scratch[root] = self.base[read];

        let codes = &self.code[read * self.n_variants..(read + 1) * self.n_variants];
        for &node in &layout.preorder {
            if node != root {
                // `node` is the variant acquired on the edge into it.
                let delta = self.delta[node][codes[node] as usize];
                scratch[node] = scratch[tree.parent[node]] + delta;
            }
        }

        let scratch = &*scratch;
        pick_best_node((0..=tree.n_mutations).map(|node| (node, scratch[node], layout.depth[node])), root)
    }
}

/// Total log-likelihood under `tree`: each read's best-attachment log-likelihood, summed.
pub fn tree_log_likelihood_with(tree: &MutationTree, scorer: &AttachmentScorer) -> f64 {
    let layout = TreeLayout::new(tree);
    let mut scratch = Vec::new();
    (0..scorer.n_reads())
        .map(|r| scorer.best_attachment(r, tree, &layout, &mut scratch).1)
        .sum()
}

/// [`tree_log_likelihood_with`] for a one-off tree.
pub fn tree_log_likelihood(matrix: &BinaryMatrix, tree: &MutationTree, rates: &ErrorRates) -> f64 {
    let scorer = AttachmentScorer::new(matrix, rates);
    tree_log_likelihood_with(tree, &scorer)
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

/// The half-open reference interval `[start, end)` a variant's REF allele occupies.
/// Two variants whose intervals overlap are mutually-exclusive alleles — a single molecule cannot
/// carry both — so they must never lie on one root-to-leaf path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantSpan {
    pub start: u64,
    pub end: u64,
}

impl VariantSpan {
    pub fn overlaps(&self, other: &VariantSpan) -> bool {
        self.start < other.end && other.start < self.end
    }
}

/// The REF interval of an `m.<pos><ref>><alt>` name: `m.310T>C` → `[310, 311)`; an empty or `-` REF
/// counts as one base.
fn parse_variant_span(variant: &str) -> Option<VariantSpan> {
    let rest = variant.trim_start_matches("m.");
    let digits_len = rest.chars().take_while(|c| c.is_ascii_digit()).count();
    let start = rest[..digits_len].parse::<u64>().ok()?;
    let ref_len = match rest[digits_len..].split('>').next().unwrap_or("") {
        "" | "-" => 1,
        a => a.len() as u64,
    };
    Some(VariantSpan { start, end: start.saturating_add(ref_len) })
}

/// Per-variant REF intervals, indexed like `variants`. A name with no parseable position gets a
/// unique one-base sentinel interval far past any real coordinate, so it never conflicts.
pub fn variant_spans(variants: &[String]) -> Vec<VariantSpan> {
    let mut unparsed: Vec<&str> = Vec::new();
    let spans = variants
        .iter()
        .enumerate()
        .map(|(i, v)| {
            parse_variant_span(v).unwrap_or_else(|| {
                unparsed.push(v);
                let sentinel = u64::MAX - 1 - i as u64;
                VariantSpan { start: sentinel, end: sentinel + 1 }
            })
        })
        .collect();
    if !unparsed.is_empty() {
        let examples = &unparsed[..unparsed.len().min(5)];
        warn!(
            "[SCITE] {}/{} variant name(s) carry no parseable reference position \
             (expected `m.<pos><ref>><alt>`, e.g. m.310T>C); they are exempt from \
             same-locus exclusivity{}. Examples: {}",
            unparsed.len(),
            variants.len(),
            if unparsed.len() == variants.len() {
                " — the constraint is inactive for this run"
            } else {
                ""
            },
            examples.join(", ")
        );
    }
    spans
}

/// One `(fp, fn)` per search-matrix variant.
fn per_variant_rates(
    spans: &[VariantSpan],
    fp_rate: f64,
    fn_rate: f64,
    hp_rates: (f64, f64),
) -> Vec<(f64, f64)> {
    spans
        .iter()
        .map(|span| {
            if crate::tracts::overlaps_tract(span) {
                hp_rates
            } else {
                (fp_rate, fn_rate)
            }
        })
        .collect()
}

/// Whether any mutation in `tree` has a strict ancestor whose REF interval overlaps its own.
pub fn violates_position_exclusivity(tree: &MutationTree, spans: &[VariantSpan]) -> bool {
    (0..tree.n_mutations).any(|node| {
        let span = spans[node];
        let mut cur = tree.parent[node];
        while cur != tree.root() {
            if spans[cur].overlaps(&span) {
                return true;
            }
            cur = tree.parent[cur];
        }
        false
    })
}

/// Lift each mutation whose REF interval overlaps an ancestor's to be a sibling of that ancestor.
/// The reattachment target is always an ancestor of the moved node (never inside its own subtree), so
/// the result stays acyclic; each move strictly decreases that node's depth, so the loop terminates.
pub fn enforce_position_exclusivity(tree: &MutationTree, spans: &[VariantSpan]) -> MutationTree {
    let mut result = tree.clone();
    loop {
        let mut changed = false;
        for node in 0..result.n_mutations {
            let span = spans[node];
            let mut cur = result.parent[node];
            let mut shallowest_conflict = None;
            while cur != result.root() {
                if spans[cur].overlaps(&span) {
                    shallowest_conflict = Some(cur);
                }
                cur = result.parent[cur];
            }
            if let Some(ancestor) = shallowest_conflict {
                result.parent[node] = result.parent[ancestor];
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    result
}

/// Build an initial mutation tree from the Neighbor-Joining haplotype tree: mutations chained in
/// descending order of the largest NJ subtree (in reads) that entirely carries them.
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

    // The clonality-ranked chain can stack mutually-exclusive overlapping alleles on one path.
    let tree = MutationTree { n_mutations: n, parent };
    enforce_position_exclusivity(&tree, &variant_spans(&hap_matrix.variants))
}

/// Two distinct indices in `0..n` by rejection sampling (panics rather than spinning forever if `n < 2`).
fn sample_two_distinct(n: usize, rng: &mut impl Rng) -> (usize, usize) {
    assert!(n >= 2, "sample_two_distinct needs at least 2 nodes (got {n})");
    let first = rng.random_range(0..n);
    let mut second = rng.random_range(0..n);
    while second == first {
        second = rng.random_range(0..n);
    }
    (first, second)
}

/// Reattach a random non-root node to a random node outside its own subtree (stays acyclic).
/// Symmetric, so it needs no Hastings correction.
pub fn propose_prune_reattach(tree: &MutationTree, rng: &mut impl Rng) -> MutationTree {
    let n = tree.n_mutations;
    let node = rng.random_range(0..n);
    let mut forbidden = tree.descendant_mask(node);
    forbidden.set(node);
    let valid_targets: Vec<usize> = (0..=n).filter(|&t| !forbidden.get(t)).collect();
    let new_parent = valid_targets[rng.random_range(0..valid_targets.len())];
    tree.with_parent(node, new_parent)
}

/// Relabel two mutation nodes: whichever tree positions `first` and `second` occupied, they now swap.
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

/// Swap the tree positions of nodes `a` and `b`. When one descends from the other, the descendant
/// takes the ancestor's place and the ancestor moves under the `choice_index`-th node of the
/// descendant's subtree. Returns the Hastings neighbourhood correction (1.0 for unrelated nodes).
fn swap_subtrees(tree: &MutationTree, a: usize, b: usize, choice_index: usize) -> (MutationTree, f64) {
    let n = tree.n_mutations;
    let (node, next_node) = if tree.ancestor_mask(b).get(a) { (b, a) } else { (a, b) };

    let mut parent = tree.parent.clone();

    if !tree.ancestor_mask(node).get(next_node) {
        parent[node] = tree.parent[next_node];
        parent[next_node] = tree.parent[node];
        return (MutationTree { n_mutations: n, parent }, 1.0);
    }

    parent[node] = tree.parent[next_node];
    let descendants: Vec<usize> = (0..n)
        .filter(|&i| tree.ancestor_mask(i).get(node))
        .collect();
    parent[next_node] = descendants[choice_index];

    let proposal = MutationTree { n_mutations: n, parent };
    let next_descendants = (0..n)
        .filter(|&i| proposal.ancestor_mask(i).get(next_node))
        .count();

    let nbh_correction = descendants.len() as f64 / next_descendants as f64;
    (proposal, nbh_correction)
}

pub fn propose_swap_subtrees(tree: &MutationTree, rng: &mut impl Rng) -> (MutationTree, f64) {
    let n = tree.n_mutations;
    let (a, b) = sample_two_distinct(n, rng);
    let (node, next_node) = if tree.ancestor_mask(b).get(a) { (b, a) } else { (a, b) };
    if !tree.ancestor_mask(node).get(next_node) {
        return swap_subtrees(tree, node, next_node, 0); // choice_index unused in this branch
    }
    let descendants_count = (0..n).filter(|&i| tree.ancestor_mask(i).get(node)).count();
    let choice_index = rng.random_range(0..descendants_count);
    swap_subtrees(tree, node, next_node, choice_index)
}

/// Weighted move dispatcher using reference SCITE's default move probabilities (no error-rate move).
pub fn propose_move(tree: &MutationTree, rng: &mut impl Rng) -> (MutationTree, f64) {
    // Only one tree exists, and both swap moves need two distinct mutations.
    if tree.n_mutations < 2 {
        return (tree.clone(), 1.0);
    }
    let r: f64 = rng.random();
    if r < 0.55 {
        (propose_prune_reattach(tree, rng), 1.0)
    } else if r < 0.95 {
        (propose_swap_labels(tree, rng), 1.0)
    } else {
        propose_swap_subtrees(tree, rng)
    }
}

/// Run one MCMC chain of `n_iterations` moves from `initial_tree` (else a random tree); returns the
/// best tree seen and its log-likelihood.
pub fn run_mcmc(
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
    n_iterations: usize,
    initial_tree: Option<&MutationTree>,
    rng: &mut impl Rng,
) -> (MutationTree, f64) {
    let n = matrix.variants.len();
    let spans = variant_spans(&matrix.variants);
    let mut current = initial_tree.cloned().unwrap_or_else(|| MutationTree::random(n, rng));
    current = enforce_position_exclusivity(&current, &spans);
    let scorer = AttachmentScorer::new(matrix, rates);
    let mut current_ll = tree_log_likelihood_with(&current, &scorer);
    let mut best = current.clone();
    let mut best_ll = current_ll;

    // Zero or one mutation admits exactly one tree.
    if n < 2 {
        return (current, current_ll);
    }

    for _ in 0..n_iterations {
        let (proposal, nbh_correction) = propose_move(&current, rng);
        // Hard constraint: overlapping (mutually exclusive) alleles never share a root-to-leaf path.
        if violates_position_exclusivity(&proposal, &spans) {
            continue;
        }
        let proposal_ll = tree_log_likelihood_with(&proposal, &scorer);

        // Track the best proposal even if the Metropolis-Hastings test rejects it.
        if proposal_ll > best_ll {
            best = proposal.clone();
            best_ll = proposal_ll;
        }

        let acceptance = nbh_correction * (proposal_ll - current_ll).exp();
        if rng.random::<f64>() < acceptance {
            current = proposal;
            current_ll = proposal_ll;
        }
    }

    (best, best_ll)
}

/// Base MCMC seed, fixed so runs are reproducible.
const MCMC_SEED: u64 = 42;

/// Run `n_chains` independent MCMC chains (seeds `MCMC_SEED + chain`; `initial_tree` seeds only the first)
/// and return the best tree found across all of them.
pub fn run_mcmc_multichain(
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
    n_iterations: usize,
    n_chains: usize,
    initial_tree: Option<&MutationTree>,
) -> (MutationTree, f64) {
    let mut best: Option<(MutationTree, f64)> = None;
    for chain in 0..n_chains {
        let mut rng = StdRng::seed_from_u64(MCMC_SEED.wrapping_add(chain as u64));
        let chain_initial = if chain == 0 { initial_tree } else { None };
        let (tree, ll) = run_mcmc(matrix, rates, n_iterations, chain_initial, &mut rng);
        best = match best {
            Some((_, best_ll)) if best_ll >= ll => best,
            _ => Some((tree, ll)),
        };
    }
    best.expect("n_chains must be >= 1")
}

/// The imputed, denoised genotype matrix from attaching every read to its best node in `tree`.
#[derive(Debug, Clone)]
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

    let scorer = AttachmentScorer::new(matrix, rates);
    let layout = TreeLayout::new(tree);
    let mut scratch = Vec::new();

    // Per-node genotypes depend only on the tree: build them once, not per read.
    let masks: Vec<Mask> = (0..=tree.n_mutations).map(|n| tree.ancestor_mask(n)).collect();

    for r in 0..n_reads {
        let (node, _ll) = scorer.best_attachment(r, tree, &layout, &mut scratch);
        attachment[r] = node;
        for v in 0..n_variants {
            data[v][r] = u8::from(masks[node].get(v));
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

/// Lewontin LD (D', r²) and one-sided Fisher enrichment for one variant pair's 2x2 read table.
struct PairCooccurrence {
    n_co: usize,
    n11: usize,
    n10: usize,
    n01: usize,
    n00: usize,
    d_prime: f64,
    r2: f64,
    odds_ratio: f64,
    fisher_p: f64,
}

pub(crate) fn lchoose(n: f64, k: f64) -> f64 {
    if k < 0.0 || k > n {
        return f64::NEG_INFINITY;
    }
    ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0)
}

/// `(odds_ratio, p)`, with `p` the one-sided Fisher's exact p for enrichment of the `n11` cell.
pub(crate) fn fisher_greater(n11: usize, n10: usize, n01: usize, n00: usize) -> (f64, f64) {
    let (a, b, c, d) = (n11 as f64, n10 as f64, n01 as f64, n00 as f64);
    let r1 = a + b;
    let r2 = c + d;
    let c1 = a + c;
    let n = a + b + c + d;
    if r1 == 0.0 || r2 == 0.0 || c1 == 0.0 || (b + d) == 0.0 {
        return (f64::NAN, 1.0);
    }
    let denom = lchoose(n, c1);
    let hi = r1.min(c1) as i64;
    let mut p = 0.0;
    let mut x = n11 as i64;
    while x <= hi {
        p += (lchoose(r1, x as f64) + lchoose(r2, c1 - x as f64) - denom).exp();
        x += 1;
    }
    let odds_ratio = if b * c > 0.0 { (a * d) / (b * c) } else { f64::INFINITY };
    (odds_ratio, p.min(1.0))
}

/// D' and r^2 for the 2x2 table `[[n11, n10], [n01, n00]]`.
fn dprime_r2(n11: usize, n10: usize, n01: usize, n00: usize) -> (f64, f64) {
    let (a, b, c, d) = (n11 as f64, n10 as f64, n01 as f64, n00 as f64);
    let n = a + b + c + d;
    if n == 0.0 {
        return (0.0, 0.0);
    }
    let p_a = (a + b) / n;
    let p_b = (a + c) / n;
    let d_stat = a / n - p_a * p_b;
    let d_max = if d_stat >= 0.0 {
        (p_a * (1.0 - p_b)).min((1.0 - p_a) * p_b)
    } else {
        (p_a * p_b).min((1.0 - p_a) * (1.0 - p_b))
    };
    let d_prime = if d_max > 1e-12 { d_stat / d_max } else { 0.0 };
    let denom = p_a * (1.0 - p_a) * p_b * (1.0 - p_b);
    let r2 = if denom > 1e-12 { (d_stat * d_stat) / denom } else { 0.0 };
    (d_prime, r2)
}

fn pair_cooccurrence(n11: usize, n10: usize, n01: usize, n00: usize) -> PairCooccurrence {
    let (d_prime, r2) = dprime_r2(n11, n10, n01, n00);
    let (odds_ratio, fisher_p) = fisher_greater(n11, n10, n01, n00);
    PairCooccurrence { n_co: n11 + n10 + n01 + n00, n11, n10, n01, n00, d_prime, r2, odds_ratio, fisher_p }
}

/// Statistics for every variant pair over the reads calling both sites;
/// `cell(v, r)` is read `r`'s call at variant `v` (`None` = uncovered).
fn write_cooccurrence(
    variants: &[String],
    n_reads: usize,
    cell: impl Fn(usize, usize) -> Option<u8>,
    path: &str,
) -> Result<()> {
    let n = variants.len();
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let (mut n11, mut n10, mut n01, mut n00) = (0usize, 0usize, 0usize, 0usize);
            for r in 0..n_reads {
                if let (Some(vi), Some(vj)) = (cell(i, r), cell(j, r)) {
                    match (vi, vj) {
                        (1, 1) => n11 += 1,
                        (1, 0) => n10 += 1,
                        (0, 1) => n01 += 1,
                        _ => n00 += 1,
                    }
                }
            }
            pairs.push((i, j, pair_cooccurrence(n11, n10, n01, n00)));
        }
    }

    let qvalues = if pairs.is_empty() {
        Vec::new()
    } else {
        let pvals: Vec<f64> = pairs.iter().map(|(_, _, s)| s.fisher_p).collect();
        adjust(&pvals, Procedure::BenjaminiHochberg)
    };

    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(w, "variant_a\tvariant_b\tn_co\tn11\tn10\tn01\tn00\td_prime\tr2\todds_ratio\tfisher_p\tfisher_q")?;
    for (idx, (i, j, s)) in pairs.iter().enumerate() {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{:.6e}\t{:.6e}",
            variants[*i], variants[*j], s.n_co, s.n11, s.n10, s.n01, s.n00,
            s.d_prime, s.r2, s.odds_ratio, s.fisher_p, qvalues[idx]
        )?;
    }
    Ok(())
}

/// Co-occurrence on the SCITE-cleaned matrix over its first `n` (non-block) variants.
pub fn write_variant_cooccurrence(matrix: &CleanedMatrix, n: usize, path: &str) -> Result<()> {
    write_cooccurrence(&matrix.variants[..n], matrix.reads.len(), |v, r| Some(matrix.data[v][r]), path)
}

/// Same statistics on the raw (pre-SCITE) matrix.
pub fn write_raw_variant_cooccurrence(matrix: &BinaryMatrix, path: &str) -> Result<()> {
    write_cooccurrence(&matrix.variants, matrix.reads.len(), |v, r| matrix.data[v][r], path)
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

/// How much evidence pins the order of one parent→child edge.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeEvidence {
    /// Log-likelihood (nats) lost by swapping this mutation with its parent.
    pub support_ll: Option<f64>,
    /// `(n_child_without_parent, n_parent_without_child)` over reads calling both sites.
    pub counts: Option<(usize, usize)>,
}

/// Significance level for the count margin behind an edge's order.
const ORDER_SIGN_TEST_ALPHA: f64 = 0.05;
/// Minimum likelihood margin (nats) for an edge's order to count as resolved.
const ORDER_MIN_SUPPORT_LL: f64 = 1.0;

impl EdgeEvidence {
    /// Whether the data resolves this edge's direction; `None` when the parent is the root.
    pub fn order_resolved(&self) -> Option<bool> {
        let (n_child, n_parent) = self.counts?;
        if self.support_ll? < ORDER_MIN_SUPPORT_LL {
            return Some(false);
        }
        // Sign test on the larger of the two counts.
        let sign_p = binomial_upper_tail(n_child + n_parent, n_child.max(n_parent), 0.5);
        Some(sign_p < ORDER_SIGN_TEST_ALPHA)
    }
}

/// Upper tail `P(X >= k)` for `X ~ Binomial(n, p)`.
pub(crate) fn binomial_upper_tail(n: usize, k: usize, p: f64) -> f64 {
    if k == 0 || n == 0 {
        return 1.0;
    }
    if k > n || !(p > 0.0 && p < 1.0) {
        return if p >= 1.0 { 1.0 } else { 0.0 };
    }
    let ln_p = p.ln();
    let ln_q = (1.0 - p).ln();
    let nf = n as f64;
    (k..=n)
        .map(|i| {
            let ii = i as f64;
            (lchoose(nf, ii) + ii * ln_p + (nf - ii) * ln_q).exp()
        })
        .sum::<f64>()
        .min(1.0)
}

/// [`EdgeEvidence`] for every mutation node, indexed by mutation id; edges from the root carry none.
pub fn edge_evidence(
    tree: &MutationTree,
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
) -> Vec<EdgeEvidence> {
    let scorer = AttachmentScorer::new(matrix, rates);
    let tree_ll = tree_log_likelihood_with(tree, &scorer);

    (0..tree.n_mutations)
        .map(|node| {
            let parent = tree.parent[node];
            if parent == tree.root() {
                return EdgeEvidence { support_ll: None, counts: None };
            }
            let swapped = swap_labels(tree, node, parent);
            EdgeEvidence {
                support_ll: Some(tree_ll - tree_log_likelihood_with(&swapped, &scorer)),
                counts: Some(orientation_counts(matrix, node, parent)),
            }
        })
        .collect()
}

fn format_support(support: Option<f64>) -> String {
    match support {
        Some(s) => format!("{s:.3}"),
        None => "NA".to_string(),
    }
}

/// The root-block decision for every candidate (variant clearing the HF gate), one row each.
pub fn write_root_block_table(
    audit: &[rootblock::RootBlockAudit],
    variants: &[String],
    path: &str,
) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(w, "variant\thf\tn_absent\tmin_q\tpartner_variant\tin_block\treason")?;
    for a in audit {
        let min_q = match a.min_q {
            Some(q) => format!("{q:.6e}"),
            None => "NA".to_string(),
        };
        // At a saturated min_q (1.0) no partner discriminated better than any other.
        let partner = match (a.partner, a.min_q) {
            (Some(u), Some(q)) if q < 1.0 => variants[u].clone(),
            _ => "NA".to_string(),
        };
        writeln!(
            w,
            "{}\t{:.6}\t{}\t{min_q}\t{partner}\t{}\t{}",
            variants[a.variant], a.hf, a.n_absent, a.in_block, a.reason.as_str()
        )?;
    }
    Ok(())
}

pub fn write_mutation_tree(
    tree: &MutationTree,
    matrix: &CleanedMatrix,
    evidence: &[EdgeEvidence],
    block_names: &[String],
    path: &str,
) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(
        w,
        "node_id\tvariant\tparent_id\tparent_variant\tn_reads_attached\t\
         order_support_ll\tn_child_without_parent\tn_parent_without_child"
    )?;
    let node_name = |id: usize| -> String {
        if id == tree.root() { "ROOT".to_string() } else { matrix.variants[id].clone() }
    };
    for node in 0..=tree.n_mutations {
        let parent = tree.parent[node];
        let n_reads = matrix.attachment.iter().filter(|&&a| a == node).count();
        let ev = evidence.get(node);
        let support = format_support(ev.and_then(|e| e.support_ll));
        let (n_child, n_parent) = match ev.and_then(|e| e.counts) {
            Some((c, p)) => (c.to_string(), p.to_string()),
            None => ("NA".to_string(), "NA".to_string()),
        };
        writeln!(
            w,
            "{node}\t{}\t{parent}\t{}\t{n_reads}\t{support}\t{n_child}\t{n_parent}",
            node_name(node),
            node_name(parent)
        )?;
    }
    // The whole root block is one row.
    if !block_names.is_empty() {
        let root_reads = matrix.attachment.iter().filter(|&&a| a == tree.root()).count();
        writeln!(
            w,
            "ROOT_BLOCK\t{}\t{}\tROOT\t{root_reads}\tNA\tNA\tNA",
            block_names.join(","),
            tree.root()
        )?;
    }
    Ok(())
}

/// Recursively emit the Newick token for the subtree at `node`; `edge_muts` are the mutations
/// (ancestral → derived) on the branch into it, into which pass-through nodes are folded.
fn emit_lineage_node(
    node: usize,
    edge_muts: &[usize],
    tree: &MutationTree,
    children: &[Vec<usize>],
    haps_by_node: &[Vec<(String, usize)>],
    variants: &[String],
    evidence: &[EdgeEvidence],
) -> String {
    let is_root = node == tree.root();
    let child_nodes = &children[node];
    let haps = &haps_by_node[node];
    let n_effective = child_nodes.len() + haps.len();
    let mut_names = |muts: &[usize]| -> String {
        muts.iter().map(|&m| variants[m].as_str()).collect::<Vec<&str>>().join(",")
    };
    let mut_support = |muts: &[usize]| -> String {
        muts.iter()
            .map(|&m| format_support(evidence.get(m).and_then(|e| e.support_ll)))
            .collect::<Vec<String>>()
            .join(",")
    };

    // Collapse: a non-root node with a single effective child is dissolved.
    if !is_root && n_effective == 1 {
        if child_nodes.len() == 1 {
            // Single child-lineage: recurse, appending this node's edge mutations.
            let c = child_nodes[0];
            let mut muts = edge_muts.to_vec();
            muts.push(c);
            return emit_lineage_node(c, &muts, tree, children, haps_by_node, variants, evidence);
        } else {
            // Single haplotype tip: the mutation(s) ride on that tip's branch.
            let (hid, count) = &haps[0];
            return format!(
                "{hid}_n{count}:{}[&&NHX:mutation={}:order_support={}:reads={count}]",
                edge_muts.len(),
                mut_names(edge_muts),
                mut_support(edge_muts)
            );
        }
    }

    // Otherwise this node is emitted (branch point, root, or dead-end).
    let mut tokens: Vec<String> = Vec::new();
    for &c in child_nodes {
        tokens.push(emit_lineage_node(c, &[c], tree, children, haps_by_node, variants, evidence));
    }

    if is_root {
        // The ancestral (no-mutation) state is always an explicit tip, even with zero reads.
        let root_reads: usize = haps.iter().map(|(_, count)| *count).sum();
        let root_hid = if haps.len() == 1 { haps[0].0.clone() } else { "Hroot".to_string() };
        tokens.push(format!("{root_hid}_n{root_reads}:0[&&NHX:reads={root_reads}]"));
        return format!("({})ROOT", tokens.join(","));
    }

    for (hid, count) in haps {
        // Haplotype attached directly to a surviving node: zero-length tip, reads only.
        tokens.push(format!("{hid}_n{count}:0[&&NHX:reads={count}]"));
    }

    let blen = edge_muts.len();
    let muts = mut_names(edge_muts);
    let support = mut_support(edge_muts);
    if tokens.is_empty() {
        // Dead-end ancestral node (no descendants, no reads) — keep the mutation.
        format!("anc{node}:{blen}[&&NHX:mutation={muts}:order_support={support}]")
    } else {
        format!(
            "({})anc{node}:{blen}[&&NHX:mutation={muts}:order_support={support}]",
            tokens.join(",")
        )
    }
}

/// Write the read lineage tree (Newick + NHX branch mutations, tsinfer-style) to `path`.
pub fn write_read_lineage_newick(
    tree: &MutationTree,
    hap_matrix: &HaplotypeMatrix,
    rates: &ErrorRates,
    evidence: &[EdgeEvidence],
    block_names: &[String],
    path: &str,
) -> Result<()> {
    let children = tree.children_of();
    let mut haps_by_node: Vec<Vec<(String, usize)>> = vec![Vec::new(); tree.n_mutations + 1];
    for hap in &hap_matrix.haplotypes {
        // Profiles also cover the root-block columns past the tree's mutations; those score alike at every node.
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
        evidence,
    );

    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    write!(w, "{body}")?;
    // The body ends with the `)ROOT` node, which carries the root block as an NHX field.
    if !block_names.is_empty() {
        write!(w, "[&&NHX:root_block={}]", block_names.join(","))?;
    }
    writeln!(w, ";")?;
    Ok(())
}

/// Reorder the mutations on each maximal unary path towards the maximum-likelihood order: seed with
/// the read-orientation order, then hill-climb on adjacent swaps (likelihood ties go to orientation).
pub fn polish_unary_path_order(
    tree: &MutationTree,
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
) -> MutationTree {
    let mut current = tree.clone();
    let scorer = AttachmentScorer::new(matrix, rates);
    let paths = maximal_unary_paths(&current);
    info!(
        "[SCITE] Unary-path polish: {} maximal unary path(s) (len≥2) over {} variants",
        paths.len(),
        tree.n_mutations
    );
    for path in paths {
        let names = |order: &[usize]| -> String {
            order
                .iter()
                .map(|&i| matrix.variants[i].as_str())
                .collect::<Vec<&str>>()
                .join(" → ")
        };
        let start_names = names(&path);
        let ll_start = tree_log_likelihood_with(&current, &scorer);

        let preferred = order_path_by_orientation(&path, matrix);
        let rank: HashMap<usize, usize> =
            preferred.iter().enumerate().map(|(i, &m)| (m, i)).collect();

        let mut order = path.clone();
        let mut ll = ll_start;
        let mut swaps = 0usize;
        let mut seeded = false;

        // Seed with the full orientation ordering when it does not cost likelihood.
        if preferred != order {
            let candidate = apply_path_reorder(&current, &order, &preferred);
            let ll_candidate = tree_log_likelihood_with(&candidate, &scorer);
            if ll_candidate >= ll - 1e-6 {
                current = candidate;
                order = preferred.clone();
                ll = ll_candidate;
                seeded = true;
            }
        }

        // Hill climb on adjacent transpositions.
        let max_passes = order.len() * order.len() + 8;
        for _ in 0..max_passes {
            let mut improved = false;
            for i in 0..order.len() - 1 {
                let mut trial_order = order.clone();
                trial_order.swap(i, i + 1);
                let trial = apply_path_reorder(&current, &order, &trial_order);
                let trial_ll = tree_log_likelihood_with(&trial, &scorer);
                let strictly_better = trial_ll > ll + LL_TIE_EPS;
                let tied_and_preferred = (trial_ll - ll).abs() <= LL_TIE_EPS
                    && rank[&trial_order[i]] < rank[&order[i]];
                if strictly_better || tied_and_preferred {
                    current = trial;
                    order = trial_order;
                    ll = ll.max(trial_ll);
                    swaps += 1;
                    improved = true;
                }
            }
            if !improved {
                break;
            }
        }

        if order == path {
            info!(
                "[SCITE] Unary-path polish: {} mutations unchanged [{}]",
                path.len(),
                start_names
            );
        } else {
            let how = match (seeded, swaps) {
                (true, 0) => "orientation order adopted".to_string(),
                (true, n) => format!("orientation order adopted, then {n} swap(s)"),
                (false, n) => format!("{n} swap(s)"),
            };
            info!(
                "[SCITE] Unary-path polish: {} mutations reordered ({how}, ΔLL={:.3}): \
                 [{}] => [{}]",
                path.len(),
                ll - ll_start,
                start_names,
                names(&order)
            );
        }
    }
    current
}

/// Maximal chains (≥ 2 mutations, top-down) in which every node but the last has exactly one child.
fn maximal_unary_paths(tree: &MutationTree) -> Vec<Vec<usize>> {
    let children = tree.children_of();
    let root = tree.root();
    let mut paths = Vec::new();
    for start in 0..tree.n_mutations {
        let p = tree.parent[start];
        let parent_is_branch = p == root || children[p].len() != 1;
        if !parent_is_branch {
            continue;
        }
        let mut path = vec![start];
        let mut cur = start;
        while let [only_child] = children[cur][..] {
            path.push(only_child);
            cur = only_child;
        }
        if path.len() >= 2 {
            paths.push(path);
        }
    }
    paths
}

fn allele_frequency(matrix: &BinaryMatrix, m: usize) -> f64 {
    let mut n0 = 0usize;
    let mut n1 = 0usize;
    for &cell in &matrix.data[m] {
        match cell {
            Some(0) => n0 += 1,
            Some(1) => n1 += 1,
            _ => {}
        }
    }
    let denom = n0 + n1;
    if denom == 0 {
        0.0
    } else {
        n1 as f64 / denom as f64
    }
}

/// P(m | o) among jointly called reads; 0.5 if never observed together with o=1.
fn nesting_p_m_given_o(matrix: &BinaryMatrix, m: usize, o: usize) -> f64 {
    let mut n11 = 0usize;
    let mut n01 = 0usize;
    for (&cm, &co) in matrix.data[m].iter().zip(matrix.data[o].iter()) {
        let (Some(vm), Some(vo)) = (cm, co) else { continue };
        if vo == 1 && vm == 1 {
            n11 += 1;
        } else if vo == 1 && vm == 0 {
            n01 += 1;
        }
    }
    let denom = n11 + n01;
    if denom == 0 {
        0.5
    } else {
        n11 as f64 / denom as f64
    }
}

/// Coverage-matched orientation counts for the pair `(m, o)`, over reads that call *both* sites.
fn orientation_counts(matrix: &BinaryMatrix, m: usize, o: usize) -> (usize, usize) {
    let mut n10 = 0usize;
    let mut n01 = 0usize;
    for (&cm, &co) in matrix.data[m].iter().zip(matrix.data[o].iter()) {
        let (Some(vm), Some(vo)) = (cm, co) else { continue };
        match (vm, vo) {
            (1, 0) => n10 += 1,
            (0, 1) => n01 += 1,
            _ => {}
        }
    }
    (n10, n01)
}

/// `Σ_o (n10 − n01)` over the other mutations `o` on the path: positive when `m` is seen without
/// `o` more often than the reverse, i.e. `m` looks ancestral.
fn orientation_score(matrix: &BinaryMatrix, m: usize, path: &[usize]) -> i64 {
    path.iter()
        .filter(|&&o| o != m)
        .map(|&o| {
            let (n10, n01) = orientation_counts(matrix, m, o);
            n10 as i64 - n01 as i64
        })
        .sum()
}

/// Ancestral → derived ordering of the mutations on one unary path (≥ 2 mutations): by orientation
/// score, then allele frequency, then mean nesting `P(m | o)`, then id.
fn order_path_by_orientation(path: &[usize], matrix: &BinaryMatrix) -> Vec<usize> {
    let mut scored: Vec<(usize, i64, f64, f64)> = path
        .iter()
        .map(|&m| {
            let nest = path
                .iter()
                .filter(|&&o| o != m)
                .map(|&o| nesting_p_m_given_o(matrix, m, o))
                .sum::<f64>()
                / (path.len() - 1) as f64;
            (m, orientation_score(matrix, m, path), allele_frequency(matrix, m), nest)
        })
        .collect();
    scored.sort_by(|a, b| {
        b.1.cmp(&a.1)
            .then_with(|| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal))
            .then_with(|| b.3.partial_cmp(&a.3).unwrap_or(std::cmp::Ordering::Equal))
            .then_with(|| a.0.cmp(&b.0))
    });
    scored.into_iter().map(|(m, _, _, _)| m).collect()
}

/// `tree` with the unary path `old_path` rechained as `new_order` (same nodes, top-down).
fn apply_path_reorder(tree: &MutationTree, old_path: &[usize], new_order: &[usize]) -> MutationTree {
    debug_assert_eq!(old_path.len(), new_order.len());
    let children = tree.children_of();
    let old_tip = *old_path.last().unwrap();
    let new_tip = *new_order.last().unwrap();

    let mut parent = tree.parent.clone();
    parent[new_order[0]] = tree.parent[old_path[0]];
    for i in 1..new_order.len() {
        parent[new_order[i]] = new_order[i - 1];
    }
    // The old tip's children (off the path) follow the new tip.
    for &c in &children[old_tip] {
        if !old_path.contains(&c) {
            parent[c] = new_tip;
        }
    }
    MutationTree { n_mutations: tree.n_mutations, parent }
}

/// Pull every block member whose REF span overlaps any other variant back into the search: imputing
/// it on every read would put two mutually-exclusive alleles on one molecule.
/// `variants` is the full variant list `rb` indexes.
fn resolve_span_conflicts(mut rb: rootblock::RootBlock, variants: &[String]) -> rootblock::RootBlock {
    let spans = variant_spans(variants);

    let (pulled_back, kept): (Vec<usize>, Vec<usize>) = rb.block.iter().copied().partition(|&v| {
        (0..variants.len()).any(|u| u != v && spans[v].overlaps(&spans[u]))
    });

    if !pulled_back.is_empty() {
        let mut names: Vec<String> = pulled_back.iter().map(|&v| variants[v].clone()).collect();
        names.sort();
        info!(
            "[SCITE] Root block: pulled {} variant(s) back into the search — their span \
             overlaps another variant, so imputing them present on every read would assert \
             two mutually-exclusive alleles on one molecule: {}",
            pulled_back.len(),
            names.join(", ")
        );
    }

    for a in rb.audit.iter_mut() {
        if pulled_back.contains(&a.variant) {
            a.in_block = false;
            a.reason = rootblock::BlockReason::SpanConflict;
        }
    }

    rb.informative.extend(pulled_back);
    rb.block = kept;
    rb.block.sort_unstable();
    rb.informative.sort_unstable();
    rb.audit.sort_unstable_by_key(|a| a.variant);
    rb
}

/// The cleaned matrix (root-block columns appended) and its haplotypes, the tips of `<prefix>.read_lineage.nwk`.
pub struct SciteOutput {
    pub cleaned: CleanedMatrix,
    pub cleaned_haplotypes: HaplotypeMatrix,
}

/// Run the MCMC mutation-tree search, attach every read to its best-fitting node, and write the
/// SCITE outputs under `output_prefix`.
pub fn run_scite_pipeline(
    binary: &BinaryMatrix,
    fp_rate: f64,
    fn_rate: f64,
    hp_rates: (f64, f64),
    n_iterations: usize,
    n_chains: usize,
    debug: bool,
    output_prefix: &str,
) -> Result<SciteOutput> {
    if binary.variants.is_empty() {
        anyhow::bail!(
            "SCITE requires at least 1 variant (got 0); \
             relax the --min-hf / --max-hf thresholds."
        );
    }
    // TODO: dropout_rate is still the global fn, not the HP-tract fn.
    let rb = resolve_span_conflicts(rootblock::partition(binary, fn_rate), &binary.variants);
    // The raw co-occurrence and root-block tables report the full matrix; everything else the search matrix.
    let full_binary = binary;
    let block_names: Vec<String> =
        rb.block.iter().map(|&v| binary.variants[v].clone()).collect();
    if !block_names.is_empty() {
        info!(
            "[SCITE] Root block: {} homoplasmic variant(s) removed from the search \
             and reported unordered at the root: {}",
            block_names.len(),
            block_names.join(", ")
        );
    }
    let search_matrix = rootblock::restrict(binary, &rb.informative);
    let binary = &search_matrix;
    // Haplotypes over the search matrix only seed the NJ-derived initial tree.
    let hap_matrix = lineage::deduplicate(binary);
    let no_informative = binary.variants.is_empty();
    if no_informative {
        info!(
            "[SCITE] Every variant is homoplasmic: the tree is ROOT plus a \
             {}-variant block, so the MCMC search is skipped.",
            block_names.len()
        );
    }
    let single_variant = binary.variants.len() == 1;
    if single_variant {
        info!(
            "[SCITE] Only 1 variant ({}) passed filtering: the mutation tree is \
             fixed (root → variant), so the MCMC search and unary-path polish are \
             skipped and the pairwise co-occurrence tables will be empty.",
            binary.variants[0]
        );
    }
    let initial_tree = lineage::neighbor_joining(&hap_matrix)
        .ok()
        .map(|nj_tree| from_nj_tree(&hap_matrix, &nj_tree));
    info!(
        "[SCITE] Initial tree: {}",
        if initial_tree.is_some() { "derived from Neighbor-Joining haplotype tree" } else { "random (NJ tree unavailable)" }
    );

    let spans = variant_spans(&binary.variants);
    let per_variant = per_variant_rates(&spans, fp_rate, fn_rate, hp_rates);
    let n_hp = spans.iter().filter(|span| crate::tracts::overlaps_tract(span)).count();
    info!(
        "[SCITE] {n_hp} of {} variant(s) overlap a known HP tract and use fp={}, fn={}",
        spans.len(),
        hp_rates.0,
        hp_rates.1
    );
    // Overlap is not transitive, so report the constrained pairs rather than grouping into sites.
    let conflict_pairs: Vec<String> = (0..spans.len())
        .flat_map(|i| (i + 1..spans.len()).map(move |j| (i, j)))
        .filter(|&(i, j)| spans[i].overlaps(&spans[j]))
        .map(|(i, j)| format!("{}/{}", binary.variants[i], binary.variants[j]))
        .collect();
    if !conflict_pairs.is_empty() {
        info!(
            "[SCITE] Same-locus exclusivity active for {} overlapping allele pair(s): {}",
            conflict_pairs.len(),
            conflict_pairs.join(", ")
        );
    }

    info!(
        "[SCITE] {} variants, {} reads. {}",
        binary.variants.len(),
        binary.reads.len(),
        if single_variant {
            "Only one tree is possible; skipping MCMC.".to_string()
        } else {
            format!("Running MCMC: {n_chains} chains x {n_iterations} iterations")
        }
    );

    let rates = ErrorRates { fp_rate, fn_rate, per_variant };
    let (tree, ll) = if no_informative {
        (MutationTree { n_mutations: 0, parent: vec![0] }, 0.0)
    } else {
        run_mcmc_multichain(binary, &rates, n_iterations, n_chains, initial_tree.as_ref())
    };
    info!("[SCITE] Best tree log-likelihood: {ll:.3}");

    let tree = if no_informative { tree } else { polish_unary_path_order(&tree, binary, &rates) };
    let tree = enforce_position_exclusivity(&tree, &spans);
    let ll_polish = tree_log_likelihood(binary, &tree, &rates);
    info!("[SCITE] After unary-path polish log-likelihood: {ll_polish:.3}");

    let mut cleaned = attach_all_reads(binary, &tree, &rates);
    let n_informative = rootblock::append_block_to_cleaned(&mut cleaned, &block_names);

    let evidence = edge_evidence(&tree, binary, &rates);
    let unresolved: Vec<String> = (0..tree.n_mutations)
        .filter(|&m| evidence[m].order_resolved() == Some(false))
        .map(|m| {
            let (n_child, n_parent) = evidence[m].counts.unwrap_or((0, 0));
            format!(
                "{} under {} ({:.2} nats, {n_child} vs {n_parent} reads)",
                binary.variants[m],
                binary.variants[tree.parent[m]],
                evidence[m].support_ll.unwrap_or(0.0)
            )
        })
        .collect();
    if !unresolved.is_empty() {
        info!(
            "[SCITE] {} edge(s) whose order this data does not resolve — they could \
             be exchanged with their parent: {}",
            unresolved.len(),
            unresolved.join("; ")
        );
    }

    write_variant_cooccurrence(&cleaned, n_informative, &format!("{output_prefix}.variant_cooccurrence.tsv"))?;
    if debug {
        write_raw_variant_cooccurrence(full_binary, &format!("{output_prefix}.raw_variant_cooccurrence.tsv"))?;
    }
    write_molecule_summary(&cleaned, &format!("{output_prefix}.molecule_summary.tsv"))?;
    write_mutation_tree(&tree, &cleaned, &evidence, &block_names, &format!("{output_prefix}.mutation_tree.tsv"))?;
    write_root_block_table(&rb.audit, &full_binary.variants, &format!("{output_prefix}.root_block.tsv"))?;

    // Deduplicate the SCITE-cleaned reads into haplotypes for the lineage tree.
    let cleaned_binary = BinaryMatrix {
        variants: cleaned.variants.clone(),
        reads: cleaned.reads.clone(),
        data: cleaned
            .data
            .iter()
            .map(|row| row.iter().map(|&b| Some(b)).collect())
            .collect(),
    };
    let cleaned_hap_matrix = lineage::deduplicate(&cleaned_binary);
    write_read_lineage_newick(
        &tree,
        &cleaned_hap_matrix,
        &rates,
        &evidence,
        &block_names,
        &format!("{output_prefix}.read_lineage.nwk"),
    )?;

    Ok(SciteOutput {
        cleaned,
        cleaned_haplotypes: cleaned_hap_matrix,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    impl ErrorRates {
        /// The same rates for every variant.
        fn uniform(fp: f64, fnr: f64) -> Self {
            ErrorRates { fp_rate: fp, fn_rate: fnr, per_variant: Vec::new() }
        }
    }

    fn mask(bits: &[usize]) -> Mask {
        let mut m = Mask::new();
        bits.iter().for_each(|&b| m.set(b));
        m
    }

    fn chain_tree() -> MutationTree {
        // 0 -> 1 -> 2 -> root(3).
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

    /// Build a matrix from per-variant cell tokens: '1' alt, '0' ref, '.' uncovered.
    fn matrix_from_rows(names: &[&str], rows: &[&str]) -> BinaryMatrix {
        let n_reads = rows[0].len();
        BinaryMatrix {
            variants: names.iter().map(|s| s.to_string()).collect(),
            reads: (0..n_reads).map(|i| format!("r{i}")).collect(),
            data: rows
                .iter()
                .map(|row| {
                    row.chars()
                        .map(|c| match c {
                            '1' => Some(1),
                            '0' => Some(0),
                            '.' => None,
                            other => panic!("bad cell {other}"),
                        })
                        .collect()
                })
                .collect(),
        }
    }

    /// A random matrix whose cells are uncovered, ref or alt with equal probability.
    fn random_matrix(n_variants: usize, n_reads: usize, rng: &mut impl Rng) -> BinaryMatrix {
        BinaryMatrix {
            variants: (0..n_variants).map(|i| format!("m.{}A>G", 100 + i)).collect(),
            reads: (0..n_reads).map(|i| format!("r{i}")).collect(),
            data: (0..n_variants)
                .map(|_| {
                    (0..n_reads)
                        .map(|_| match rng.random_range(0..3u8) {
                            0 => None,
                            1 => Some(0),
                            _ => Some(1),
                        })
                        .collect()
                })
                .collect(),
        }
    }

    /// Read `r`'s calls at every variant.
    fn profile(matrix: &BinaryMatrix, r: usize) -> Vec<Option<u8>> {
        matrix.data.iter().map(|row| row[r]).collect()
    }

    /// Run `write` on a temp file called `name` and return what it wrote.
    fn written(name: &str, write: impl FnOnce(&str) -> Result<()>) -> String {
        let path = std::env::temp_dir().join(name);
        write(path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        content
    }

    /// Run the pipeline into the temp prefix `himito_scite_<name>`; returns its output and that prefix.
    fn run_pipeline(name: &str, binary: &BinaryMatrix, debug: bool) -> (SciteOutput, String) {
        let prefix = std::env::temp_dir().join(format!("himito_scite_{name}")).to_str().unwrap().to_string();
        let out = run_scite_pipeline(binary, 0.01, 0.05, (0.01, 0.05), 200, 1, debug, &prefix).unwrap();
        (out, prefix)
    }

    fn read_output(prefix: &str, suffix: &str) -> String {
        std::fs::read_to_string(format!("{prefix}{suffix}")).unwrap()
    }

    /// `n_present` alt calls followed by `n_absent` ref calls.
    fn present_then_absent(n_present: usize, n_absent: usize) -> String {
        format!("{}{}", "1".repeat(n_present), "0".repeat(n_absent))
    }

    #[test]
    fn ancestor_mask_includes_self_and_all_ancestors_up_to_root() {
        let tree = chain_tree();
        assert_eq!(tree.ancestor_mask(0), mask(&[0, 1, 2]));
        assert_eq!(tree.ancestor_mask(1), mask(&[1, 2]));
        assert_eq!(tree.ancestor_mask(2), mask(&[2]));
        assert_eq!(tree.ancestor_mask(3), Mask::new());
    }

    #[test]
    fn descendant_mask_excludes_self_and_includes_transitive_children() {
        let tree = chain_tree();
        assert_eq!(tree.descendant_mask(2), mask(&[0, 1]));
        assert_eq!(tree.descendant_mask(1), mask(&[0]));
        assert_eq!(tree.descendant_mask(0), Mask::new());
        assert_eq!(tree.descendant_mask(3), mask(&[0, 1, 2]));
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
        let ancestor_mask = mask(&[0, 1]);
        let rates = ErrorRates::uniform(0.1, 0.2);

        let expected = (1.0 - rates.fn_rate).ln() + rates.fn_rate.ln();
        let actual = attachment_log_likelihood(&profile, &ancestor_mask, &rates);
        assert!((actual - expected).abs() < 1e-12);
    }

    #[test]
    fn attachment_log_likelihood_scores_true_negative_and_false_positive() {
        let profile = vec![Some(0), Some(1)];
        let ancestor_mask = Mask::new();
        let rates = ErrorRates::uniform(0.1, 0.2);

        let expected = (1.0 - rates.fp_rate).ln() + rates.fp_rate.ln();
        let actual = attachment_log_likelihood(&profile, &ancestor_mask, &rates);
        assert!((actual - expected).abs() < 1e-12);
    }

    /// Columns past `per_variant` — the root-block genotypes appended beyond `n_mutations`.
    #[test]
    fn columns_past_per_variant_use_the_global_rate() {
        let profile = vec![Some(1), Some(1)];
        let rates = ErrorRates { fp_rate: 0.1, fn_rate: 0.2, per_variant: vec![(0.3, 0.4)] };
        // Variant 0 uses its own fp (0.3); the extra column falls back to 0.1.
        let expected = 0.3_f64.ln() + 0.1_f64.ln();
        let actual = attachment_log_likelihood(&profile, &Mask::new(), &rates);
        assert!((actual - expected).abs() < 1e-12, "{actual} vs {expected}");
    }

    #[test]
    fn best_attachment_picks_the_node_matching_observed_mutations() {
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 2, 2] };
        let rates = ErrorRates::uniform(0.01, 0.3);
        let profile = vec![Some(1), Some(0)];

        let expected_ll = (1.0 - rates.fn_rate).ln() + (1.0 - rates.fp_rate).ln();
        let (node, ll) = best_attachment(&profile, &tree, &rates);
        assert_eq!(node, 0);
        assert!((ll - expected_ll).abs() < 1e-12);
    }

    #[test]
    fn near_ties_attach_to_the_node_with_fewest_mutations() {
        // Chain root(2) -> A(0) -> B(1): uncovered at A and ref at B fits root and A equally well.
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        let rates = ErrorRates::uniform(0.01, 0.1);
        let (node, ll) = best_attachment(&[None, Some(0)], &tree, &rates);
        assert_eq!(node, tree.root());
        assert!((ll - (1.0 - rates.fp_rate).ln()).abs() < 1e-12);
        // A fully uncovered read scores 0 everywhere.
        assert_eq!(best_attachment(&[None, None], &tree, &rates), (tree.root(), 0.0));

        // The scorer path must not fabricate an alt call at the uncovered site either.
        let matrix = matrix_from_rows(&["A", "B"], &["1.", "10"]);
        let cleaned = attach_all_reads(&matrix, &tree, &rates);
        assert_eq!(cleaned.attachment[1], tree.root());
        assert_eq!(cleaned.data, vec![vec![1, 0], vec![1, 0]]);
    }

    #[test]
    fn absence_at_a_high_fn_variant_attaches_below_that_variant() {
        // Chain root -> A -> B; the read carries A but not B.
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        let profile = [Some(1u8), Some(0)];
        let global = ErrorRates::uniform(0.01, 0.05);
        assert_eq!(best_attachment(&profile, &tree, &global).0, 0, "global fn keeps the absence above B");

        let high = ErrorRates { fp_rate: 0.01, fn_rate: 0.05, per_variant: vec![(0.01, 0.05), (0.01, 0.999)] };
        assert_eq!(best_attachment(&profile, &tree, &high).0, 1, "high fn at B explains the absence as dropout");
        let matrix = matrix_from_rows(&["A", "B"], &["1", "0"]);
        assert_eq!(attach_all_reads(&matrix, &tree, &high).attachment, vec![1]);
    }

    /// Known-tract spans, including a whole-tract allele row and a boundary straddler, take the HP rates.
    #[test]
    fn per_variant_rates_are_higher_inside_known_hp_tracts() {
        let variants = vec![
            "m.3243A>G".to_string(),
            "m.16189T>C".to_string(),
            "m.16183ACCCCCTCCCC>ACCCCCCCCCC".to_string(),
            "m.16193CA>C".to_string(),
            "m.16182A>C".to_string(),
        ];
        let spans = variant_spans(&variants);
        let rates = per_variant_rates(&spans, 0.01, 0.05, (0.02, 0.10));
        assert_eq!(rates[0], (0.01, 0.05), "outside every tract");
        assert_eq!(rates[1], (0.02, 0.10), "HVS1-C SNV");
        assert_eq!(rates[2], (0.02, 0.10), "whole-tract allele");
        assert_eq!(rates[3], (0.02, 0.10), "straddles the tract boundary");
        assert_eq!(rates[4], (0.01, 0.05), "one base before the tract");
    }

    #[test]
    fn tree_log_likelihood_sums_best_attachment_over_all_reads() {
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 2, 2] };
        let rates = ErrorRates::uniform(0.01, 0.3);
        let matrix = matrix_from_rows(&["A", "B"], &["10", "01"]);

        let per_read_ll = (1.0 - rates.fn_rate).ln() + (1.0 - rates.fp_rate).ln();
        let actual = tree_log_likelihood(&matrix, &tree, &rates);
        assert!((actual - 2.0 * per_read_ll).abs() < 1e-12);
    }

    #[test]
    fn incremental_scorer_matches_reference_best_attachment_on_random_trees() {
        let mut rng = StdRng::seed_from_u64(0xC0FFEE);
        // The degenerate settings are the ones that matter.
        let mut settings: Vec<ErrorRates> = [
            (0.005, 0.05),  // pacbio default
            (0.0001, 0.01), // ont-denoised default
            (0.0, 0.0),     // both at the floor
            (0.0, 0.05),    // fp only at the floor
            (0.005, 0.0),   // fn only at the floor
            (1.0, 1.0),     // both at the ceiling
            (0.5, 0.5),     // uninformative
        ]
        .into_iter()
        .map(|(fp, fnr)| ErrorRates::uniform(fp, fnr))
        .collect();
        settings.push(ErrorRates {
            per_variant: (0..9).map(|i| if i % 2 == 0 { (0.2, 0.4) } else { (0.001, 0.02) }).collect(),
            ..ErrorRates::uniform(0.005, 0.05)
        });

        for rates in &settings {
            for trial in 0..300usize {
                let matrix = random_matrix(1 + trial % 9, 1 + trial % 7, &mut rng);
                let tree = MutationTree::random(matrix.variants.len(), &mut rng);
                let scorer = AttachmentScorer::new(&matrix, rates);
                let layout = TreeLayout::new(&tree);
                let mut scratch = Vec::new();
                let mut want_total = 0.0;
                for r in 0..matrix.reads.len() {
                    let (want_node, want_ll) = best_attachment(&profile(&matrix, r), &tree, rates);
                    let (got_node, got_ll) = scorer.best_attachment(r, &tree, &layout, &mut scratch);
                    assert_eq!(got_node, want_node, "node mismatch (trial {trial}, read {r})");
                    assert!(
                        (got_ll - want_ll).abs() < 1e-9,
                        "ll mismatch (trial {trial}, read {r}): {got_ll} vs {want_ll}"
                    );
                    want_total += want_ll;
                }
                let total = tree_log_likelihood(&matrix, &tree, rates);
                assert!((total - want_total).abs() < 1e-9, "trial {trial}: {total} vs {want_total}");
            }
        }
    }

    #[test]
    fn zero_error_rates_are_clamped_to_finite_scores() {
        // `--fp-rate 0` is a valid input. Star tree: r0 fits node 0 exactly; r1 (alt at both) fits no node.
        let matrix = matrix_from_rows(&["m.100A>G", "m.200C>T"], &["11", "01"]);
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 2, 2] };
        let rates = ErrorRates::uniform(0.0, 0.0);
        let scorer = AttachmentScorer::new(&matrix, &rates);
        let layout = TreeLayout::new(&tree);
        for r in 0..2 {
            let (want_node, want_ll) = best_attachment(&profile(&matrix, r), &tree, &rates);
            let (node, ll) = scorer.best_attachment(r, &tree, &layout, &mut Vec::new());
            assert!(ll.is_finite() && want_ll.is_finite(), "read {r}: scorer {ll}, reference {want_ll}");
            assert_eq!(node, want_node, "read {r}");
            assert!((ll - want_ll).abs() < 1e-9, "read {r}: {ll} vs {want_ll}");
        }
        assert_eq!(scorer.best_attachment(0, &tree, &layout, &mut Vec::new()).0, 0);
    }

    #[test]
    fn tree_layout_is_a_preorder_with_mutation_depths() {
        let mut rng = StdRng::seed_from_u64(7);
        for n in 1..12usize {
            let tree = MutationTree::random(n, &mut rng);
            let layout = TreeLayout::new(&tree);
            assert_eq!(layout.preorder.len(), n + 1);
            let mut seen = vec![false; n + 1];
            for &node in &layout.preorder {
                // The forward scoring pass needs each parent scored before its children.
                assert!(!seen[node], "node {node} visited twice");
                assert!(node == tree.root() || seen[tree.parent[node]], "node {node} visited before its parent");
                seen[node] = true;
                assert_eq!(layout.depth[node], tree.ancestor_mask(node).count_ones(), "depth of node {node}");
            }
        }
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

        // Three leaves under one root: H0 (5 reads) carries neither variant, H1 (3) A only, H2 (2) A and B.
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
    fn from_nj_tree_breaks_an_overlapping_chain_off_the_clonality_ranking() {
        // One haplotype carries m.310T>C, the other m.310T>TC.
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["m.310T>C".to_string(), "m.310T>TC".to_string()],
            haplotypes: vec![
                lineage::Haplotype { id: "H0".to_string(), profile: vec![Some(1), Some(0)], reads: vec!["r0".to_string()], count: 1 },
                lineage::Haplotype { id: "H1".to_string(), profile: vec![Some(0), Some(1)], reads: vec!["r1".to_string()], count: 1 },
            ],
        };
        let nj_tree = lineage::neighbor_joining(&hap_matrix).unwrap();
        let tree = from_nj_tree(&hap_matrix, &nj_tree);
        assert!(!violates_position_exclusivity(&tree, &variant_spans(&hap_matrix.variants)));
    }

    fn span(start: u64, end: u64) -> VariantSpan {
        VariantSpan { start, end }
    }

    #[test]
    fn variant_spans_parses_ref_interval_and_falls_back_uniquely() {
        let variants = vec![
            "m.310T>C".to_string(),   // SNP: one base
            "m.310T>TC".to_string(),  // insertion anchored at 310
            "m.310TAA>T".to_string(), // deletion of 311-312, anchored at 310
            "m.8274C>T".to_string(),
            "m.5->AC".to_string(), // Himito's "no anchor base" REF placeholder
            "m.7A>-".to_string(),  // ...and the ALT-side placeholder
            "A".to_string(),
            "B".to_string(),
        ];
        let s = variant_spans(&variants);
        assert_eq!(s[0], span(310, 311));
        assert_eq!(s[1], span(310, 311));
        assert_eq!(s[2], span(310, 313));
        assert_eq!(s[3], span(8274, 8275));
        assert_eq!(s[4], span(5, 6));
        assert_eq!(s[5], span(7, 8));
        // Non-parseable names never collide with each other or with real loci.
        assert!(!s[6].overlaps(&s[7]));
        assert!(s[6].start > 8274 && s[7].start > 8274);
    }

    #[test]
    fn variant_spans_flags_a_deletion_against_a_substitution_inside_its_span() {
        // The case a start-position-only check misses.
        let s = variant_spans(&["m.310TAA>T".to_string(), "m.311A>G".to_string()]);
        assert!(s[0].overlaps(&s[1]));
        // The base just past the deletion is untouched.
        let outside = variant_spans(&["m.310TAA>T".to_string(), "m.313G>A".to_string()]);
        assert!(!outside[0].overlaps(&outside[1]));
    }

    #[test]
    fn violates_position_exclusivity_flags_overlapping_alleles_on_one_path() {
        // spans: 16126(0) -> 310(1) -> 310(2), a chain.
        let spans = vec![span(16126, 16127), span(310, 311), span(310, 311)];
        let chain = MutationTree { n_mutations: 3, parent: vec![3, 0, 1, 3] };
        assert!(violates_position_exclusivity(&chain, &spans));

        // Same two 310 alleles as siblings under 16126: no path carries both.
        let branched = MutationTree { n_mutations: 3, parent: vec![3, 0, 0, 3] };
        assert!(!violates_position_exclusivity(&branched, &spans));

        // m.311A>G [311,312) beneath the deletion m.310TAA>T [310,313).
        let deletion = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        assert!(violates_position_exclusivity(&deletion, &[span(310, 313), span(311, 312)]));
    }

    #[test]
    fn enforce_position_exclusivity_lifts_overlapping_alleles_to_siblings() {
        // (spans, chain, expected parents): 16126 -> 310T>TC -> 310T>C; and three 310 alleles under
        // 500, which needs a second pass since lifting node 2 leaves node 3 still under it.
        let cases = [
            (vec![span(16126, 16127), span(310, 311), span(310, 311)], vec![3, 0, 1, 3], vec![3, 0, 0, 3]),
            (
                vec![span(500, 501), span(310, 311), span(310, 311), span(310, 311)],
                vec![4, 0, 1, 2, 4],
                vec![4, 0, 0, 0, 4],
            ),
        ];
        for (spans, parent, expected) in cases {
            let chain = MutationTree { n_mutations: spans.len(), parent };
            let fixed = enforce_position_exclusivity(&chain, &spans);
            assert_eq!(fixed.parent, expected);
            assert!(!violates_position_exclusivity(&fixed, &spans));
        }
    }

    #[test]
    fn enforce_position_exclusivity_always_returns_an_acyclic_tree() {
        // 12 variants over 4 loci; every fourth is a 5-base deletion straddling its neighbours.
        let spans: Vec<VariantSpan> = (0..12u64)
            .map(|i| {
                let locus = (i % 4) * 10;
                span(locus, locus + if i % 4 == 0 { 5 } else { 1 })
            })
            .collect();
        for seed in 0..500u64 {
            let mut rng = StdRng::seed_from_u64(seed);
            let tree = MutationTree::random(12, &mut rng);
            let fixed = enforce_position_exclusivity(&tree, &spans);
            assert_valid_tree(&fixed);
            assert!(!violates_position_exclusivity(&fixed, &spans), "seed {seed} left an overlapping pair on one path");
        }
    }

    #[test]
    fn run_mcmc_never_places_same_position_alleles_on_one_path() {
        // Reads: 6 with 16126 + insertion, 6 with 16126 + substitution, 3 with 16126 only, 3 germline.
        let variants = ["m.16126T>C", "m.310T>TC", "m.310T>C"];
        let matrix = matrix_from_rows(
            &variants,
            &["111111111111111000", "111111000000000000", "000000111111000000"],
        );
        let rates = ErrorRates::uniform(0.01, 0.1);
        let mut rng = StdRng::seed_from_u64(11);
        let (tree, _ll) = run_mcmc(&matrix, &rates, 3000, None, &mut rng);

        assert!(
            !violates_position_exclusivity(&tree, &variant_spans(&matrix.variants)),
            "the two position-310 alleles must never lie on one root-to-leaf path"
        );
    }

    /// Test-only RNG that replays a fixed sequence of `u64` words, cycling once exhausted.
    struct ScriptedRng { words: Vec<u64>, idx: usize }
    impl rand::RngCore for ScriptedRng {
        fn next_u32(&mut self) -> u32 { (self.next_u64() >> 32) as u32 }
        fn next_u64(&mut self) -> u64 {
            let v = self.words[self.idx % self.words.len()];
            self.idx += 1;
            v
        }
        fn fill_bytes(&mut self, dest: &mut [u8]) {
            for chunk in dest.chunks_mut(8) {
                let v = self.next_u64().to_le_bytes();
                chunk.copy_from_slice(&v[..chunk.len()]);
            }
        }
    }

    #[test]
    fn run_mcmc_records_a_rejected_but_superior_swap_subtrees_proposal() {
        let matrix = matrix_from_rows(&["m.100A>G", "m.200C>T", "m.300G>A", "m.400T>C"], &["0", "1", "0", "0"]);
        let rates = ErrorRates::uniform(0.01, 0.6);
        let start = MutationTree { n_mutations: 4, parent: vec![2, 2, 4, 4, 4] };
        let start_ll = tree_log_likelihood(&matrix, &start, &rates);

        let words = vec![
            u64::MAX,   // propose_move: random::<f64>() ~= 0.9999999999999999 >= 0.95 -> swap-subtrees.
            1u64 << 62, // sample_two_distinct: random_range(0..4) == 1 -> first = 1.
            1u64 << 63, // sample_two_distinct: random_range(0..4) == 2 -> second = 2 (!= first, no retry).
            0u64,       // swap_subtrees's only same-lineage candidate: random_range(0..1) == 0 regardless.
            u64::MAX,   // MH accept test: random::<f64>() ~= 0.9999999999999999 >= acceptance -> reject.
        ];
        let mut rng = ScriptedRng { words, idx: 0 };

        let (best, best_ll) = run_mcmc(&matrix, &rates, 1, Some(&start), &mut rng);

        let proposal = MutationTree { n_mutations: 4, parent: vec![2, 4, 1, 4, 4] };
        let proposal_ll = tree_log_likelihood(&matrix, &proposal, &rates);
        assert!(proposal_ll > start_ll, "setup: the proposal must beat the start ({proposal_ll} vs {start_ll})");
        // Mirrors `run_mcmc`: acceptance = nbh_correction * exp(proposal - current).
        let mh_draw: f64 = ScriptedRng { words: vec![u64::MAX], idx: 0 }.random();
        let acceptance = 0.5 * (proposal_ll - start_ll).exp();
        assert!(acceptance <= mh_draw, "setup: the scripted draw must reject the proposal ({acceptance} vs {mh_draw})");

        assert_eq!(best.parent, proposal.parent, "a rejected-but-superior proposal must still be returned as best");
        assert!((best_ll - proposal_ll).abs() < 1e-9, "reported {best_ll}, actual {proposal_ll}");
    }

    #[test]
    fn run_mcmc_reaches_the_optimum_on_an_exhaustively_searchable_space() {
        let matrix = matrix_from_rows(
            &["m.100A>G", "m.200C>T", "m.300G>A"],
            &["111111000", "111000000", "000111000"],
        );
        let rates = ErrorRates::uniform(0.01, 0.05);

        // Brute-force every acyclic parent vector on 3 mutations.
        let mut optimum = f64::NEG_INFINITY;
        for p0 in 0..=3usize {
            for p1 in 0..=3usize {
                for p2 in 0..=3usize {
                    let t = MutationTree { n_mutations: 3, parent: vec![p0, p1, p2, 3] };
                    let acyclic = (0..3).all(|s| {
                        let (mut cur, mut steps) = (s, 0);
                        while cur != 3 && steps <= 3 {
                            cur = t.parent[cur];
                            steps += 1;
                        }
                        cur == 3
                    });
                    if acyclic {
                        optimum = optimum.max(tree_log_likelihood(&matrix, &t, &rates));
                    }
                }
            }
        }

        for seed in [3u64, 17, 101, 2024] {
            let mut rng = StdRng::seed_from_u64(seed);
            let (_best, best_ll) = run_mcmc(&matrix, &rates, 3000, None, &mut rng);
            assert!(best_ll >= optimum - 1e-9, "seed {seed}: got {best_ll}, optimum {optimum}");
        }
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
        // Relabeling 0<->2 turns "0->1->2->root" into "2->1->0->root".
        let swapped = swap_labels(&chain_tree(), 0, 2);
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
        // node=1, next_node=2: node1's only current "descendant" is itself.
        let (proposal, correction) = swap_subtrees(&tree, 1, 2, 0);
        assert_eq!(proposal.parent, vec![2, 4, 1, 4, 4]);
        assert!((correction - 0.5).abs() < 1e-12);
    }

    #[test]
    fn swap_labels_star_topology_does_not_double_fire_self_loop_fixup() {
        // Star tree: leaves 0 and 3 have no children, so relabeling them changes nothing.
        let tree = MutationTree { n_mutations: 4, parent: vec![4, 4, 4, 4, 4] };
        let swapped = swap_labels(&tree, 0, 3);
        assert_valid_tree(&swapped);
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
    fn propose_move_leaves_a_single_mutation_tree_unchanged() {
        // `sample_two_distinct` would otherwise loop forever.
        let tree = MutationTree { n_mutations: 1, parent: vec![1, 1] };
        let mut rng = StdRng::seed_from_u64(4);
        for _ in 0..100 {
            let (proposal, correction) = propose_move(&tree, &mut rng);
            assert_eq!(proposal.parent, tree.parent);
            assert_eq!(correction, 1.0);
        }
    }

    /// A ancestral to B: 6 reads carry both, 3 only A, and one carries B with A dropped out.
    fn ancestor_pair_matrix() -> BinaryMatrix {
        matrix_from_rows(&["A", "B"], &["1111111110", "1111110001"])
    }

    #[test]
    fn run_mcmc_recovers_a_true_ancestor_relationship() {
        let mut rng = StdRng::seed_from_u64(42);
        let (tree, _ll) = run_mcmc(&ancestor_pair_matrix(), &ErrorRates::uniform(0.01, 0.1), 3000, None, &mut rng);
        assert_eq!(tree.parent[1], 0);
    }

    #[test]
    fn run_mcmc_multichain_finds_the_true_tree_with_fewer_iterations_per_chain() {
        let rates = ErrorRates::uniform(0.01, 0.1);
        let (tree, _ll) = run_mcmc_multichain(&ancestor_pair_matrix(), &rates, 1000, 4, None);
        assert_eq!(tree.parent[1], 0);
    }

    #[test]
    fn attach_all_reads_imputes_missing_and_corrects_false_negatives() {
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] }; // A ancestor of B
        let rates = ErrorRates::uniform(0.01, 0.1);
        let matrix = BinaryMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            reads: vec!["r_full".to_string(), "r_missingB".to_string(), "r_errorA".to_string()],
            data: vec![vec![Some(1), Some(1), Some(0)], vec![Some(1), None, Some(1)]],
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
        let content = written("himito_test_write_cleaned_matrix.csv", |p| write_cleaned_matrix(&matrix, p));
        assert_eq!(content, "variant,r1,r2\nA,1,0\nB,1,1\n");
    }

    #[test]
    fn write_variant_cooccurrence_reports_pairwise_ld_and_fisher_stats() {
        let matrix = small_cleaned_matrix();
        let content = written("himito_test_write_cooccurrence.tsv", |p| write_variant_cooccurrence(&matrix, 2, p));
        let mut lines = content.lines();
        assert_eq!(
            lines.next().unwrap(),
            "variant_a\tvariant_b\tn_co\tn11\tn10\tn01\tn00\td_prime\tr2\todds_ratio\tfisher_p\tfisher_q"
        );
        // r1=[A,B]=[1,1], r2=[A,B]=[0,1] -> n11=1 (r1), n01=1 (r2)
        let fields: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(&fields[0..7], &["A", "B", "2", "1", "0", "1", "0"]);
        assert_eq!(lines.next(), None);
    }

    #[test]
    fn fisher_greater_matches_known_2x2_table() {
        // Tea-tasting-style table with a clear enrichment in n11.
        let (odds_ratio, p) = fisher_greater(8, 2, 1, 9);
        assert!((odds_ratio - 36.0).abs() < 1e-9);
        assert!(p < 0.01, "expected a small one-sided p-value, got {p}");
    }

    #[test]
    fn dprime_r2_is_one_for_perfect_co_occurrence() {
        let (d_prime, r2) = dprime_r2(10, 0, 0, 10);
        assert!((d_prime - 1.0).abs() < 1e-9);
        assert!((r2 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn dprime_r2_is_zero_for_independent_variants() {
        // p_a = p_b = 0.5 and the joint frequency exactly matches p_a*p_b
        let (d_prime, r2) = dprime_r2(5, 5, 5, 5);
        assert!(d_prime.abs() < 1e-9);
        assert!(r2.abs() < 1e-9);
    }

    #[test]
    fn write_raw_variant_cooccurrence_restricts_to_jointly_covered_reads() {
        // r1: A alt, B missing (excluded); r2: A ref, B alt; r3: A alt, B alt
        let matrix = matrix_from_rows(&["A", "B"], &["101", ".11"]);
        let content = written("himito_test_write_raw_cooccurrence.tsv", |p| write_raw_variant_cooccurrence(&matrix, p));
        let fields: Vec<&str> = content.lines().nth(1).unwrap().split('\t').collect();
        // r1 dropped (B uncovered): n_co=2, n11=1 (r3), n01=1 (r2)
        assert_eq!(&fields[0..7], &["A", "B", "2", "1", "0", "1", "0"]);
    }

    #[test]
    fn write_molecule_summary_lists_variant_count_and_names_per_read() {
        let matrix = small_cleaned_matrix();
        let content = written("himito_test_write_molecule_summary.tsv", |p| write_molecule_summary(&matrix, p));
        assert_eq!(content, "read_name\tn_variants\tvariants\nr1\t2\tA,B\nr2\t1\tB\n");
    }

    #[test]
    fn write_mutation_tree_lists_every_node_with_parent_and_attached_read_count() {
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        let matrix = small_cleaned_matrix();
        let evidence = vec![
            EdgeEvidence { support_ll: None, counts: None },
            EdgeEvidence { support_ll: Some(4.25), counts: Some((3, 7)) },
        ];
        let nodes = "node_id\tvariant\tparent_id\tparent_variant\tn_reads_attached\t\
                     order_support_ll\tn_child_without_parent\tn_parent_without_child\n\
                     0\tA\t2\tROOT\t1\tNA\tNA\tNA\n\
                     1\tB\t0\tA\t1\t4.250\t3\t7\n\
                     2\tROOT\t2\tROOT\t0\tNA\tNA\tNA\n";
        let content = written("himito_test_write_mutation_tree.tsv", |p| write_mutation_tree(&tree, &matrix, &evidence, &[], p));
        assert_eq!(content, nodes);

        // The whole root block is one extra row.
        let block = ["m.750A>G".to_string(), "m.16183A>C".to_string()];
        let content =
            written("himito_test_write_mutation_tree_block.tsv", |p| write_mutation_tree(&tree, &matrix, &evidence, &block, p));
        assert_eq!(content, format!("{nodes}ROOT_BLOCK\tm.750A>G,m.16183A>C\t2\tROOT\t0\tNA\tNA\tNA\n"));
    }

    #[test]
    fn write_root_block_table_reports_each_candidate() {
        use rootblock::{BlockReason, RootBlockAudit};
        let variants = vec!["m.750A>G".to_string(), "m.310T>TC".to_string(), "m.16189T>C".to_string()];
        let audit = vec![
            RootBlockAudit {
                variant: 0, hf: 0.90, n_absent: 12,
                min_q: Some(0.42), partner: Some(1), in_block: true,
                reason: BlockReason::UnstructuredAbsences,
            },
            // No partner test: NA, not 0.
            RootBlockAudit {
                variant: 1, hf: 0.99, n_absent: 3,
                min_q: None, partner: None, in_block: true,
                reason: BlockReason::HighRawVaf,
            },
            // A saturated min_q blanks the (arbitrary) partner.
            RootBlockAudit {
                variant: 2, hf: 0.87, n_absent: 20,
                min_q: Some(1.0), partner: Some(1), in_block: true,
                reason: BlockReason::UnstructuredAbsences,
            },
        ];
        let content = written("himito_test_root_block.tsv", |p| write_root_block_table(&audit, &variants, p));
        assert_eq!(
            content,
            "variant\thf\tn_absent\tmin_q\tpartner_variant\tin_block\treason\n\
             m.750A>G\t0.900000\t12\t4.200000e-1\tm.310T>TC\ttrue\tunstructured_absences\n\
             m.310T>TC\t0.990000\t3\tNA\tNA\ttrue\thigh_raw_vaf\n\
             m.16189T>C\t0.870000\t20\t1.000000e0\tNA\ttrue\tunstructured_absences\n"
        );
    }

    #[test]
    fn write_read_lineage_newick_collapses_and_annotates_reads_and_mutations() {
        use lineage::Haplotype;

        // root(4) -> A(0) -> B(1) -> {C(2), D(3)}.
        let tree = MutationTree { n_mutations: 4, parent: vec![4, 0, 1, 1, 4] };
        let rates = ErrorRates::uniform(0.01, 0.05);
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["mA".into(), "mB".into(), "mC".into(), "mD".into()],
            haplotypes: vec![
                Haplotype { id: "H0".into(), profile: vec![Some(1), Some(1), Some(1), Some(0)], reads: vec![], count: 9 },
                Haplotype { id: "H1".into(), profile: vec![Some(1), Some(1), Some(0), Some(1)], reads: vec![], count: 7 },
                Haplotype { id: "H2".into(), profile: vec![Some(1), Some(1), None, Some(1)], reads: vec![], count: 4 },
                Haplotype { id: "H3".into(), profile: vec![Some(1), Some(1), Some(0), Some(0)], reads: vec![], count: 2 },
            ],
        };
        // mA hangs from the root (no order to support); the rest carry evidence.
        let evidence = vec![
            EdgeEvidence { support_ll: None, counts: None },
            EdgeEvidence { support_ll: Some(12.5), counts: Some((0, 9)) },
            EdgeEvidence { support_ll: Some(0.0), counts: Some((2, 2)) },
            EdgeEvidence { support_ll: Some(3.0), counts: Some((1, 4)) },
        ];
        let content = written("himito_test_read_lineage.nwk", |p| {
            write_read_lineage_newick(&tree, &hap_matrix, &rates, &evidence, &[], p)
        });

        // No haplotype has zero mutations, so the root state becomes a synthetic zero-read tip.
        assert_eq!(
            content,
            "((H0_n9:1[&&NHX:mutation=mC:order_support=0.000:reads=9],\
             (H1_n7:0[&&NHX:reads=7],H2_n4:0[&&NHX:reads=4])anc3:1[&&NHX:mutation=mD:order_support=3.000],\
             H3_n2:0[&&NHX:reads=2])anc1:2[&&NHX:mutation=mA,mB:order_support=NA,12.500],\
             Hroot_n0:0[&&NHX:reads=0])ROOT;\n"
        );
    }

    #[test]
    fn write_read_lineage_newick_uses_all_zero_haplotype_as_root_tip() {
        use lineage::Haplotype;

        // root(1) -> A(0).
        let tree = MutationTree { n_mutations: 1, parent: vec![1, 1] };
        let rates = ErrorRates::uniform(0.01, 0.05);
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["mA".into()],
            haplotypes: vec![
                Haplotype { id: "H0000".into(), profile: vec![Some(0)], reads: vec![], count: 5 },
                Haplotype { id: "H0001".into(), profile: vec![Some(1)], reads: vec![], count: 3 },
            ],
        };
        let evidence = vec![EdgeEvidence { support_ll: None, counts: None }];
        let write = |block: &[String], name: &str| {
            written(name, |p| write_read_lineage_newick(&tree, &hap_matrix, &rates, &evidence, block, p))
        };

        let body = "(H0001_n3:1[&&NHX:mutation=mA:order_support=NA:reads=3],H0000_n5:0[&&NHX:reads=5])ROOT";
        assert_eq!(write(&[], "himito_test_read_lineage_root_tip.nwk"), format!("{body};\n"));
        // The root block is an NHX field on the ROOT node.
        let block = ["m.750A>G".to_string(), "m.16183A>C".to_string()];
        assert_eq!(
            write(&block, "himito_test_read_lineage_root_block.nwk"),
            format!("{body}[&&NHX:root_block=m.750A>G,m.16183A>C];\n")
        );
    }

    #[test]
    fn polish_unary_path_order_is_noop_when_order_already_matches_af() {
        let rates = ErrorRates::uniform(0.001, 0.001);
        let matrix = matrix_from_rows(&["A", "B"], &["1111111111", "0000001111"]);
        let correct = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        assert_eq!(polish_unary_path_order(&correct, &matrix, &rates).parent, correct.parent);
    }

    #[test]
    fn path_order_uses_jointly_covered_counts_not_coverage_biased_allele_frequency() {
        // A is ancestral: among the 4 reads covering both sites it appears once without B and never the other way round.
        let matrix = matrix_from_rows(&["A", "B"], &["1100................", "10001111111111100000"]);
        assert!(allele_frequency(&matrix, 0) < allele_frequency(&matrix, 1));
        assert_eq!(order_path_by_orientation(&[0, 1], &matrix), vec![0, 1]);
        assert_eq!(order_path_by_orientation(&[1, 0], &matrix), vec![0, 1]);

        // A tree in the wrong order (B ancestral) is reordered.
        let wrong = MutationTree { n_mutations: 2, parent: vec![1, 2, 2] };
        let polished = polish_unary_path_order(&wrong, &matrix, &ErrorRates::uniform(0.001, 0.05));
        assert_eq!(polished.parent[0], 2, "A should hang from the root");
        assert_eq!(polished.parent[1], 0, "B should hang from A");
    }

    #[test]
    fn path_order_falls_back_to_allele_frequency_when_joint_counts_are_tied() {
        // A and B never appear apart among jointly covered reads (n10 = n01 = 0).
        let matrix = matrix_from_rows(&["A", "B"], &["11001111", "1100...."]);
        assert_eq!(order_path_by_orientation(&[1, 0], &matrix), vec![0, 1]);
    }

    #[test]
    fn polish_keeps_climbing_after_the_first_candidate_ordering() {
        // Adopting one candidate ordering stops short; adjacent swaps reach the optimum.
        let rates = ErrorRates::uniform(0.001, 0.05);
        let matrix = matrix_from_rows(&["m0", "m1", "m2"], &["10010001111", "0.101.0010.", "101010011.."]);
        let tree = MutationTree { n_mutations: 3, parent: vec![1, 2, 3, 3] };

        let polished = polish_unary_path_order(&tree, &matrix, &rates);

        assert_eq!(polished.parent[0], 3, "m0 should hang from the root");
        assert_eq!(polished.parent[2], 0, "m2 should hang from m0");
        assert_eq!(polished.parent[1], 2, "m1 should hang from m2");
        let ll = tree_log_likelihood(&matrix, &polished, &rates);
        assert!(ll > -6.73, "expected the exhaustive-ML likelihood, got {ll}");
    }

    #[test]
    fn polish_is_deterministic_and_terminates_on_a_likelihood_flat_path() {
        // All three variants are carried as a block.
        let rates = ErrorRates::uniform(0.001, 0.05);
        let matrix = matrix_from_rows(&["A", "B", "C"], &["111111110000"; 3]);
        let flipped = MutationTree { n_mutations: 3, parent: vec![1, 2, 3, 3] };

        let polished = polish_unary_path_order(&flipped, &matrix, &rates);

        assert_eq!(polished.parent, vec![3, 0, 1, 3]);
        let before = tree_log_likelihood(&matrix, &flipped, &rates);
        let after = tree_log_likelihood(&matrix, &polished, &rates);
        assert!((after - before).abs() < 1e-9, "flat path: {before} vs {after}");
    }

    #[test]
    fn polish_unary_path_order_moves_distal_children_to_new_tip() {
        // Path ROOT→0→1 with distal children 2, 3 under tip 1; m1 (in every read) is ancestral to m0.
        let rates = ErrorRates::uniform(0.01, 0.05);
        let matrix = matrix_from_rows(&["m0", "m1", "m2", "m3"], &["00001111", "11111111", "00000011", "00001100"]);
        let tree = MutationTree { n_mutations: 4, parent: vec![4, 0, 1, 1, 4] };
        let polished = polish_unary_path_order(&tree, &matrix, &rates);
        // 1 under ROOT, 0 under 1, distal 2 and 3 under the new tip 0.
        assert_eq!(polished.parent, vec![1, 4, 0, 0, 4]);
    }

    #[test]
    fn edge_evidence_measures_how_firmly_reads_pin_each_edge() {
        let rates = ErrorRates::uniform(0.001, 0.05);
        // Block-carried: nothing separates A, B and C, so no edge's order has any support.
        let flat = matrix_from_rows(&["A", "B", "C"], &["111111110000"; 3]);
        let chain = MutationTree { n_mutations: 3, parent: vec![3, 0, 1, 3] };
        let evidence = edge_evidence(&chain, &flat, &rates);
        assert_eq!(evidence[0], EdgeEvidence { support_ll: None, counts: None }, "A hangs from the root");
        for ev in &evidence[1..] {
            assert!(ev.support_ll.unwrap().abs() < 1e-9, "expected no evidence, got {ev:?}");
            assert_eq!(ev.counts, Some((0, 0)));
        }

        // Six reads carry A without B, none B without A.
        let pinned = matrix_from_rows(&["A", "B"], &["1111111111", "1111000000"]);
        let chain = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        let evidence = edge_evidence(&chain, &pinned, &rates);
        assert!(evidence[1].support_ll.unwrap() > 5.0, "expected decisive evidence, got {:?}", evidence[1]);
        assert_eq!(evidence[1].counts, Some((0, 6)));
    }

    #[test]
    fn order_resolved_flags_a_near_tie_even_when_the_likelihood_looks_confident() {
        // Confident by likelihood, but the count margin is one read.
        let inverted = EdgeEvidence { support_ll: Some(5.989), counts: Some((8, 9)) };
        assert_eq!(inverted.order_resolved(), Some(false));

        // A flat likelihood is unresolved whatever the counts say.
        let flat = EdgeEvidence { support_ll: Some(0.0), counts: Some((4, 5)) };
        assert_eq!(flat.order_resolved(), Some(false));

        for counts in [(3usize, 28usize), (2, 29), (7, 16), (11, 23)] {
            let ev = EdgeEvidence { support_ll: Some(26.9), counts: Some(counts) };
            assert_eq!(ev.order_resolved(), Some(true), "{counts:?} should be resolved");
        }

        let root_edge = EdgeEvidence { support_ll: None, counts: None };
        assert_eq!(root_edge.order_resolved(), None);
    }

    #[test]
    fn resolve_span_conflicts_pulls_every_overlapping_block_member_back() {
        // A blocked indel overlapping an informative SNV at 310; a block chain where B overlaps only
        // A; two block members overlapping only each other.
        let cases: [(&[&str], &[usize]); 3] = [
            (&["m.310T>TC", "m.310T>C"], &[0]),
            (&["m.300A>G", "m.300AT>A", "m.301TG>T"], &[1, 2]),
            (&["m.500A>G", "m.500A>T"], &[0, 1]),
        ];
        for (names, block) in cases {
            let variants: Vec<String> = names.iter().map(|s| s.to_string()).collect();
            let n = variants.len();
            let rb = rootblock::RootBlock {
                block: block.to_vec(),
                informative: (0..n).filter(|v| !block.contains(v)).collect(),
                audit: block
                    .iter()
                    .map(|&variant| rootblock::RootBlockAudit {
                        variant, hf: 0.9, n_absent: 12, min_q: None, partner: None, in_block: true,
                        reason: rootblock::BlockReason::HighRawVaf,
                    })
                    .collect(),
            };
            let rb = resolve_span_conflicts(rb, &variants);
            assert!(rb.block.is_empty(), "{names:?}: every overlapping block member must be pulled back");
            assert_eq!(rb.informative, (0..n).collect::<Vec<_>>());
            assert!(
                rb.audit.iter().all(|a| !a.in_block && a.reason == rootblock::BlockReason::SpanConflict),
                "{names:?}: the audit must record the span conflict"
            );
        }
    }

    #[test]
    fn run_scite_pipeline_handles_more_than_63_variants() {
        // 65 variants need a second mask word; all-ref reads past `n_var` keep every HF below the root-block gate.
        let (n_var, n_reads) = (65, 131);
        let matrix = BinaryMatrix {
            variants: (0..n_var).map(|i| format!("m.{i}")).collect(),
            reads: (0..n_reads).map(|j| format!("r{j}")).collect(),
            data: (0..n_var).map(|i| (0..n_reads).map(|j| Some(u8::from(j > i && j <= n_var))).collect()).collect(),
        };
        let (_, prefix) = run_pipeline("65var", &matrix, false);
        for suffix in [".variant_cooccurrence.tsv", ".molecule_summary.tsv", ".mutation_tree.tsv", ".read_lineage.nwk"] {
            let path = format!("{prefix}{suffix}");
            assert!(std::path::Path::new(&path).exists(), "expected output file {path} to exist");
            std::fs::remove_file(&path).ok();
        }
    }

    #[test]
    fn run_scite_pipeline_handles_a_single_variant() {
        // The only possible tree is root → m.100A>G, and there is no pair to report.
        let matrix = matrix_from_rows(&["m.100A>G"], &["111000."]);
        let (_, prefix) = run_pipeline("1var", &matrix, false);

        let tree_tsv = read_output(&prefix, ".mutation_tree.tsv");
        assert!(tree_tsv.contains("0\tm.100A>G\t1\tROOT\t"), "expected m.100A>G under ROOT, got:\n{tree_tsv}");
        let cooc = read_output(&prefix, ".variant_cooccurrence.tsv");
        assert_eq!(cooc.lines().count(), 1, "expected header only in .variant_cooccurrence.tsv");
        let nwk = read_output(&prefix, ".read_lineage.nwk");
        assert!(nwk.contains("mutation=m.100A>G"), "got newick: {nwk}");
    }

    #[test]
    fn run_scite_pipeline_with_all_variants_homoplasmic_does_not_panic() {
        let germline = present_then_absent(95, 5);
        let binary = matrix_from_rows(&["m.750A>G", "m.16183A>C"], &[&germline, &germline]);
        let (out, _) = run_pipeline("rb_allblock", &binary, false);
        assert_eq!(out.cleaned.variants.len(), 2);
    }

    #[test]
    fn one_informative_variant_after_partition_still_takes_the_single_variant_path() {
        let binary = matrix_from_rows(&["m.750A>G", "m.302A>ACC"], &[&present_then_absent(95, 5), &"10".repeat(50)]);
        let (out, prefix) = run_pipeline("rb_single", &binary, true);

        assert_eq!(out.cleaned.variants, vec!["m.302A>ACC", "m.750A>G"], "block variants stay in the cleaned matrix");
        let tsv = read_output(&prefix, ".mutation_tree.tsv");
        assert!(tsv.contains("ROOT_BLOCK\tm.750A>G\t"), "the germline variant should be blocked:\n{tsv}");
        assert!(tsv.contains("0\tm.302A>ACC\t1\tROOT\t"), "the informative variant stays a tree node:\n{tsv}");
        // The raw data still carries the blocked variant's dropouts.
        let raw = read_output(&prefix, ".raw_variant_cooccurrence.tsv");
        assert!(raw.contains("m.750A>G"), "the blocked variant must stay in the raw co-occurrence table:\n{raw}");
    }

    #[test]
    fn variant_cooccurrence_omits_block_variants_entirely() {
        let binary = matrix_from_rows(
            &["m.750A>G", "m.16183A>C", "m.9477G>A"],
            &[&present_then_absent(95, 5), &"10".repeat(50), &present_then_absent(30, 70)],
        );
        let rb = rootblock::partition(&binary, 0.05);
        assert_eq!((rb.block, rb.informative), (vec![0], vec![1, 2]), "test premise: v0 alone blocks");

        let (_, prefix) = run_pipeline("rb_cooc_omit", &binary, false);
        let cooc = read_output(&prefix, ".variant_cooccurrence.tsv");
        assert_eq!(cooc.lines().count(), 2, "expected one informative pair:\n{cooc}");
        assert!(!cooc.contains("m.750A>G"), "the block variant must be omitted entirely:\n{cooc}");
    }
}
