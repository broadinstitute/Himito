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

/// A growable bitset over mutation/variant ids, backed by `Vec<u64>` words.
///
/// Bit `i` set means mutation/variant `i` is present. Backing the mask with a
/// word vector (rather than a single `u64`) lifts the old 63-variant ceiling on
/// the SCITE search: any number of mutations is representable. An empty mask
/// (no words) represents the all-zero germline state.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Mask {
    words: Vec<u64>,
}

impl Mask {
    /// An empty mask with no bits set.
    pub fn new() -> Self {
        Mask { words: Vec::new() }
    }

    /// Build a mask from an explicit list of set bit indices (test/helper use).
    pub fn from_bits(bits: &[usize]) -> Self {
        let mut m = Mask::new();
        for &b in bits {
            m.set(b);
        }
        m
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

/// Log-likelihood differences below this are treated as exact ties. Attachment
/// and unary-path ordering both have genuinely flat regions (see
/// `polish_unary_path_order`), and resolving those by float noise makes the
/// output depend on iteration order rather than on the data.
const LL_TIE_EPS: f64 = 1e-9;

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
    ///
    /// Note: `swap_subtrees` (Task 8) builds its own "descendants" list using
    /// `ancestor_mask` with a filter that *does* include `node` itself — those
    /// two conventions are intentionally different. This method is only used in
    /// `propose_prune_reattach` to exclude a node's own subtree from reattachment targets.
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
    ancestor_mask: &Mask,
    rates: &ErrorRates,
) -> f64 {
    let mut ll = 0.0;
    for (i, call) in profile.iter().enumerate() {
        let Some(observed) = call else { continue };
        let expected_mutated = ancestor_mask.get(i);
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
///
/// Ties are broken parsimoniously: among nodes that explain the read equally
/// well, the one carrying the fewest mutations wins. A read that is uncovered
/// at a variant pays nothing for having it imputed as alt, so without this
/// tie-break every partially-covered read drifts to the deepest tied node and
/// `attach_all_reads` writes alt calls that no observation supports.
pub fn best_attachment(
    profile: &[Option<u8>],
    tree: &MutationTree,
    rates: &ErrorRates,
) -> (usize, f64) {
    let mut best_node = tree.root();
    let mut best_ll = f64::NEG_INFINITY;
    let mut best_muts = usize::MAX;
    for node in 0..=tree.n_mutations {
        let mask = tree.ancestor_mask(node);
        let ll = attachment_log_likelihood(profile, &mask, rates);
        let n_muts = mask.count_ones();
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

/// The half-open reference interval `[start, end)` a variant's REF allele
/// occupies, in the 1-based coordinates Himito variant names use.
///
/// Two variants whose intervals overlap are mutually-exclusive alleles — a
/// single molecule cannot carry both — so they must never lie on one
/// root-to-leaf path. Comparing intervals rather than start positions is what
/// catches a deletion against a substitution *inside* the deleted span:
/// `m.310TAA>T` covers `[310, 313)` and `m.311A>G` covers `[311, 312)`, which
/// share no start position but cannot coexist on one molecule.
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

/// Parse the REF interval out of an `m.<pos><ref>><alt>` name (as built by
/// `call::generate_variant_name` and `lineage.rs`): `m.310T>C` → `[310, 311)`,
/// `m.310TAA>T` → `[310, 313)`, `m.310T>TC` → `[310, 311)`.
///
/// Indel names anchor on the base *before* the event, and that anchor base is
/// included in the interval. That is deliberately conservative: it keeps an
/// indel exclusive with a substitution at its anchor (`m.310T>TC` vs
/// `m.310T>C`), which is correct under VCF allele semantics — both names assert
/// a different content for base 310 — at the cost of also separating a deletion
/// from a substitution at the retained anchor base. Splitting those onto sibling
/// branches is the safe direction to err.
///
/// A name with no parseable position gets a unique one-base sentinel interval
/// far past any real coordinate, so it can never spuriously conflict with
/// another variant. `ok` is false in that case so callers can report it.
fn parse_variant_span(variant: &str, index: usize) -> (VariantSpan, bool) {
    let rest = variant.trim_start_matches("m.");
    let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();

    let Ok(start) = digits.parse::<u64>() else {
        let sentinel = u64::MAX - 1 - index as u64;
        return (VariantSpan { start: sentinel, end: sentinel + 1 }, false);
    };

    // REF allele: everything between the position and the `>`. `-` is Himito's
    // placeholder for "no anchor base available" (an indel at the very start of
    // a reference block); treat it as one base so it still conflicts locally
    // rather than with nothing at all.
    let ref_allele = rest[digits.len()..].split('>').next().unwrap_or("");
    let ref_len = match ref_allele {
        "" | "-" => 1,
        a => a.len() as u64,
    };

    (VariantSpan { start, end: start.saturating_add(ref_len) }, true)
}

/// Per-variant REF intervals, indexed like `variants`.
///
/// Warns if any name is unparseable: the sentinel fallback makes such a variant
/// unconstrained, so a wholesale naming-convention change would silently turn
/// the exclusivity guarantee into a no-op. This makes that failure loud.
pub fn variant_spans(variants: &[String]) -> Vec<VariantSpan> {
    let mut spans = Vec::with_capacity(variants.len());
    let mut unparsed: Vec<&str> = Vec::new();
    for (i, v) in variants.iter().enumerate() {
        let (span, ok) = parse_variant_span(v, i);
        if !ok {
            unparsed.push(v.as_str());
        }
        spans.push(span);
    }
    if !unparsed.is_empty() {
        let examples: Vec<&str> = unparsed.iter().copied().take(5).collect();
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

/// Whether any mutation in `tree` has a strict ancestor whose REF interval
/// overlaps its own — i.e. two mutually-exclusive alleles placed on one lineage.
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

/// Return a tree equal to `tree` except that every mutation whose REF interval
/// overlaps one of its ancestors' is lifted to become a sibling of the
/// shallowest such ancestor, so no two overlapping alleles remain on one
/// root-to-leaf path. The reattachment target is always an ancestor of the moved
/// node (never inside its own subtree), so the result stays acyclic; each move
/// strictly decreases that node's depth, so the loop terminates.
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

    // The clonality-ranked chain can stack mutually-exclusive overlapping
    // alleles on one path; split them onto sibling branches before search.
    let tree = MutationTree { n_mutations: n, parent };
    enforce_position_exclusivity(&tree, &variant_spans(&hap_matrix.variants))
}

/// Two distinct indices in `0..n` by rejection sampling.
///
/// Panics when `n < 2`: with fewer than two mutations there is no such pair, and
/// the rejection loop would otherwise spin forever. Callers must screen for that
/// (see [`propose_move`], which skips both swap moves in that case).
fn sample_two_distinct(n: usize, rng: &mut impl Rng) -> (usize, usize) {
    assert!(n >= 2, "sample_two_distinct needs at least 2 nodes (got {n})");
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
    let mut forbidden = tree.descendant_mask(node);
    forbidden.set(node);
    let valid_targets: Vec<usize> = (0..=n).filter(|&t| !forbidden.get(t)).collect();
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

/// Weighted move dispatcher, matching reference SCITE's default move-type
/// probabilities for the non-error-rate moves (`0.55` prune-and-reattach,
/// `0.4` swap-labels, `0.05` swap-subtrees), renormalized to sum to 1 since
/// this plan has no error-rate move at all.
pub fn propose_move(tree: &MutationTree, rng: &mut impl Rng) -> (MutationTree, f64) {
    // With 0 or 1 mutations only one tree exists, and both swap moves need two
    // distinct mutation nodes to pick from. Propose the current tree unchanged
    // (always accepted, no state change) rather than sampling an impossible pair.
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
    let spans = variant_spans(&matrix.variants);
    let mut current = initial_tree.cloned().unwrap_or_else(|| MutationTree::random(n, rng));
    current = enforce_position_exclusivity(&current, &spans);
    let mut current_ll = tree_log_likelihood(matrix, &current, rates);
    let mut best = current.clone();
    let mut best_ll = current_ll;

    // A single mutation (or none) admits exactly one tree — it hangs off the
    // root — so there is nothing to search. Return it directly instead of
    // burning `n_iterations` no-op proposals.
    if n < 2 {
        return (current, current_ll);
    }

    for _ in 0..n_iterations {
        let (proposal, nbh_correction) = propose_move(&current, rng);
        // Hard constraint: never let two mutually-exclusive overlapping
        // alleles share a root-to-leaf path.
        if violates_position_exclusivity(&proposal, &spans) {
            continue;
        }
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
            data[v][r] = u8::from(mask.get(v));
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

/// Pairwise Lewontin linkage-disequilibrium (D', r^2) and one-sided Fisher's
/// exact enrichment for a variant pair's 2x2 read table
/// `[[n11, n10], [n01, n00]]`, mirroring `h1_cooccurrence` in the
/// `heteroplasmy_lineage_test.py` companion script.
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

fn lchoose(n: f64, k: f64) -> f64 {
    if k < 0.0 || k > n {
        return f64::NEG_INFINITY;
    }
    ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0)
}

/// One-sided Fisher's exact p for enrichment of the `n11` cell. Returns
/// `(odds_ratio, p_greater)`.
fn fisher_greater(n11: usize, n10: usize, n01: usize, n00: usize) -> (f64, f64) {
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

fn write_pair_cooccurrence_table(
    variants: &[String],
    pairs: &[(usize, usize, PairCooccurrence)],
    path: &str,
) -> Result<()> {
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

/// Pairwise co-occurrence stats on the SCITE-cleaned (tree-imputed, no
/// missing values) matrix: every read counts toward every pair's table.
pub fn write_variant_cooccurrence(matrix: &CleanedMatrix, path: &str) -> Result<()> {
    let n = matrix.variants.len();
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let (mut n11, mut n10, mut n01, mut n00) = (0usize, 0usize, 0usize, 0usize);
            for r in 0..matrix.reads.len() {
                match (matrix.data[i][r], matrix.data[j][r]) {
                    (1, 1) => n11 += 1,
                    (1, 0) => n10 += 1,
                    (0, 1) => n01 += 1,
                    _ => n00 += 1,
                }
            }
            pairs.push((i, j, pair_cooccurrence(n11, n10, n01, n00)));
        }
    }
    write_pair_cooccurrence_table(&matrix.variants, &pairs, path)
}

/// Same statistics on the raw (pre-SCITE) matrix. Reads with a missing call
/// (`None`) at either site are excluded from that pair's table, so `n_co`
/// counts only reads jointly covering both variants.
pub fn write_raw_variant_cooccurrence(matrix: &BinaryMatrix, path: &str) -> Result<()> {
    let n = matrix.variants.len();
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let (mut n11, mut n10, mut n01, mut n00) = (0usize, 0usize, 0usize, 0usize);
            for r in 0..matrix.reads.len() {
                if let (Some(vi), Some(vj)) = (matrix.data[i][r], matrix.data[j][r]) {
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
    write_pair_cooccurrence_table(&matrix.variants, &pairs, path)
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
///
/// Both fields are `None` when the parent is the root: there is no mutation
/// above to swap with, so the edge has no order to support.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeEvidence {
    /// Log-likelihood (nats) lost by swapping this mutation with its parent.
    /// Near zero means the two could be exchanged without the data noticing.
    pub support_ll: Option<f64>,
    /// `(n_child_without_parent, n_parent_without_child)` over reads calling
    /// both sites. A near-tie here means the order is a coin flip even when
    /// `support_ll` looks large, because each read of margin is worth
    /// `ln((1-β)/β)` nats — at β = 0.05 that is 2.94 nats for a single read.
    pub counts: Option<(usize, usize)>,
}

/// Significance level for the count margin behind an edge's order.
const ORDER_SIGN_TEST_ALPHA: f64 = 0.05;
/// Minimum likelihood margin (nats) for an edge's order to count as resolved.
const ORDER_MIN_SUPPORT_LL: f64 = 1.0;

impl EdgeEvidence {
    /// Whether the data actually resolves this edge's parent→child direction:
    /// `Some(true)` resolved, `Some(false)` arbitrary, `None` no order to resolve
    /// (the parent is the root).
    ///
    /// Both tests must pass. The likelihood margin alone is not enough, because
    /// each read of count margin is worth `ln((1-β)/β)` nats — 2.94 at β = 0.05 —
    /// so a single read can look like decisive evidence. The sign test on the
    /// counts asks the blunter question: would a margin this large arise from a
    /// coin flip? On the ont-r10 eval set this pair of tests separates every
    /// correctly-ordered edge (p ≤ 0.047) from both arbitrary ones (p = 0.5).
    pub fn order_resolved(&self) -> Option<bool> {
        let (n_child, n_parent) = self.counts?;
        let support = self.support_ll?;
        if support < ORDER_MIN_SUPPORT_LL {
            return Some(false);
        }
        Some(sign_test_p(n_child, n_parent) < ORDER_SIGN_TEST_ALPHA)
    }
}

/// Two-sided-in-spirit sign test on a count margin: `P(X ≥ max(a, b))` for
/// `X ~ Binomial(a + b, 0.5)`. Returns 1.0 for an empty table.
fn sign_test_p(a: usize, b: usize) -> f64 {
    let n = a + b;
    if n == 0 {
        return 1.0;
    }
    let k = a.max(b);
    let ln_half_pow = (n as f64) * 0.5f64.ln();
    (k..=n)
        .map(|i| (lchoose(n as f64, i as f64) + ln_half_pow).exp())
        .sum::<f64>()
        .min(1.0)
}

/// Log-likelihood lost by swapping `node` with its parent, i.e. the evidence
/// pinning this edge's direction. `None` when the parent is the root.
///
/// This is the same quantity the unary-path hill climb maximises, reported so a
/// reader can tell an order the data insists on from one that was picked
/// arbitrarily off a flat likelihood.
pub fn edge_order_support(
    tree: &MutationTree,
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
    node: usize,
) -> Option<f64> {
    let parent = tree.parent[node];
    if parent == tree.root() || node == tree.root() {
        return None;
    }
    let swapped = swap_labels(tree, node, parent);
    Some(tree_log_likelihood(matrix, tree, rates) - tree_log_likelihood(matrix, &swapped, rates))
}

/// [`EdgeEvidence`] for every mutation node, indexed by mutation id.
pub fn edge_evidence(
    tree: &MutationTree,
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
) -> Vec<EdgeEvidence> {
    (0..tree.n_mutations)
        .map(|node| {
            let parent = tree.parent[node];
            if parent == tree.root() {
                EdgeEvidence { support_ll: None, counts: None }
            } else {
                EdgeEvidence {
                    support_ll: edge_order_support(tree, matrix, rates, node),
                    counts: Some(orientation_counts(matrix, node, parent)),
                }
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

pub fn write_mutation_tree(
    tree: &MutationTree,
    matrix: &CleanedMatrix,
    evidence: &[EdgeEvidence],
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
    Ok(())
}

/// Recursively emit the Newick token for the subtree at `node`. `edge_muts` are
/// the mutation ids accumulated (ancestral→derived) on the branch leading into
/// the node that will actually be emitted here — pass-through mutation nodes fold
/// their mutation into this list rather than becoming their own node. Each branch
/// reports an `order_support` list parallel to its `mutation` list, so a collapsed
/// run whose internal order the data does not resolve is visible as such.
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
        // The root carries no heteroplasmic variants, so it is always emitted as
        // an explicit ancestral haplotype tip. When reads best-attach to the root
        // that haplotype supplies the id and read count; otherwise a synthetic
        // zero-read tip keeps the ancestral state visible in every lineage tree.
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

/// Write the read lineage tree (Newick + NHX branch mutations, tsinfer-style)
/// to `path`. Tips are haplotypes (labeled `H<id>_n<reads>`, with a machine-
/// readable `reads` NHX field); branches carry the mutations acquired, and
/// pass-through mutation nodes are collapsed so mutations ride on a single
/// branch (comma-joined, ancestral→derived). Each branch also carries an
/// `order_support` list parallel to `mutation`: the nats of likelihood pinning
/// each mutation above the next, or `NA` where there is no mutation above.
pub fn write_read_lineage_newick(
    tree: &MutationTree,
    hap_matrix: &HaplotypeMatrix,
    rates: &ErrorRates,
    evidence: &[EdgeEvidence],
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
        evidence,
    );

    let mut w = BufWriter::new(File::create(path).with_context(|| format!("Cannot create {path}"))?);
    writeln!(w, "{body};")?;
    Ok(())
}

/// Re-order the mutations on each maximal unary path towards the maximum-
/// likelihood ordering, using the orientation ranking (see
/// [`order_path_by_orientation`]) as both the starting candidate and the
/// tie-break on likelihood-flat stretches.
///
/// Why a path needs polishing at all: on a unary path only the reads that
/// attach to an *interior* node distinguish one ordering from another, and a
/// false negative at the ancestral mutation turns a "carries both" read into a
/// "carries only the child" read — precisely the observation that argues for the
/// reversed order. So the MCMC's choice of ordering along a chain is decided by
/// a handful of reads, and where no read attaches in the interior the likelihood
/// is exactly flat and the choice is pure noise.
///
/// The search is a hill climb over adjacent transpositions. Swapping two
/// neighbours on a unary path changes exactly one implied genotype (the prefix
/// between them), which makes it the natural elementary move; a full-path
/// candidate, in contrast, is all-or-nothing, so one bad pair anywhere on a long
/// path throws away every other correction. A swap is taken when it strictly
/// improves the likelihood, or when the likelihood is tied and the swap moves
/// towards the orientation-preferred order — the latter is what makes flat
/// stretches deterministic instead of iteration-order dependent.
pub fn polish_unary_path_order(
    tree: &MutationTree,
    matrix: &BinaryMatrix,
    rates: &ErrorRates,
) -> MutationTree {
    let mut current = tree.clone();
    let paths = maximal_unary_paths(&current);
    info!(
        "[SCITE] Unary-path polish: {} maximal unary path(s) (len≥2) over {} variants",
        paths.len(),
        tree.n_mutations
    );
    for path in paths {
        if path.len() < 2 {
            continue;
        }
        let names = |order: &[usize]| -> String {
            order
                .iter()
                .map(|&i| matrix.variants[i].as_str())
                .collect::<Vec<&str>>()
                .join(" → ")
        };
        let start_names = names(&path);
        let ll_start = tree_log_likelihood(matrix, &current, rates);

        // Rank preferred by the data-driven orientation statistic; used to seed
        // the climb and to break exact likelihood ties in a fixed direction.
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
            let ll_candidate = tree_log_likelihood(matrix, &candidate, rates);
            if ll_candidate >= ll - 1e-6 {
                current = candidate;
                order = preferred.clone();
                ll = ll_candidate;
                seeded = true;
            }
        }

        // Hill climb on adjacent transpositions. Each accepted move either
        // strictly raises the likelihood or strictly reduces the number of
        // inversions against `rank` at equal likelihood, so this terminates; the
        // pass cap is belt-and-braces against float noise.
        let max_passes = order.len() * order.len() + 8;
        for _ in 0..max_passes {
            let mut improved = false;
            for i in 0..order.len() - 1 {
                let mut trial_order = order.clone();
                trial_order.swap(i, i + 1);
                let trial = apply_path_reorder(&current, &order, &trial_order);
                let trial_ll = tree_log_likelihood(matrix, &trial, rates);
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

/// Maximal unary paths: each is a sequence of ≥1 mutation nodes where every
/// non-final node has exactly one child. Starts only at nodes whose parent is
/// the root or a branching node (≠1 child).
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
        while children[cur].len() == 1 {
            let nxt = children[cur][0];
            if nxt >= tree.n_mutations {
                break;
            }
            path.push(nxt);
            cur = nxt;
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

/// Coverage-matched orientation counts for the pair `(m, o)`, over reads that
/// call *both* sites: `(n10, n01)` = (m alt & o ref, m ref & o alt).
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

/// Signed orientation score of `m` against the other mutations on its path:
/// `Σ_o (n10 - n01)`, positive when `m` tends to be the one seen *without* its
/// partner and therefore the more ancestral of the pair.
///
/// Under a perfect phylogeny with false-negative rate β, if `m` is ancestral to
/// `o` then a `(m=0, o=1)` read requires a dropout at `m` (probability β) while
/// `(m=1, o=0)` reads are genuine, so the pairwise log-likelihood ratio between
/// the two orientations is `(n10 - n01) · ln((1-β)/β)`. That weight is a
/// positive constant shared by every pair on the path, so it scales the score
/// without changing the ordering and is left out. Summing the signed pairwise
/// margins and sorting (a Borda ranking of the orientation tournament) is the
/// standard cheap approximation to minimum-violation ranking; the likelihood
/// hill climb in `polish_unary_path_order` then refines it.
///
/// This is deliberately restricted to jointly covered reads: comparing raw
/// allele frequencies across sites with different depth is not meaningful, and
/// mitochondrial coverage is uneven enough (read ends, NUMT filtering) for that
/// bias to invert the order on its own.
fn orientation_score(matrix: &BinaryMatrix, m: usize, path: &[usize]) -> i64 {
    path.iter()
        .filter(|&&o| o != m)
        .map(|&o| {
            let (n10, n01) = orientation_counts(matrix, m, o);
            n10 as i64 - n01 as i64
        })
        .sum()
}

/// Ancestral → derived ordering of the mutations on one unary path.
///
/// Primary key is the coverage-matched pairwise orientation score; allele
/// frequency (a marginal statistic, which also uses reads covering only one of
/// the two sites) breaks its ties, then mean nesting, then the variant index so
/// the result is deterministic on genuinely tied data.
fn order_path_by_orientation(path: &[usize], matrix: &BinaryMatrix) -> Vec<usize> {
    let mut scored: Vec<(usize, i64, f64, f64)> = path
        .iter()
        .map(|&m| {
            let others: Vec<usize> = path.iter().copied().filter(|&o| o != m).collect();
            let nest = if others.is_empty() {
                0.5
            } else {
                others
                    .iter()
                    .map(|&o| nesting_p_m_given_o(matrix, m, o))
                    .sum::<f64>()
                    / others.len() as f64
            };
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

fn apply_path_reorder(tree: &MutationTree, old_path: &[usize], new_order: &[usize]) -> MutationTree {
    debug_assert_eq!(old_path.len(), new_order.len());
    let children = tree.children_of();
    let external_parent = tree.parent[old_path[0]];
    let old_tip = *old_path.last().unwrap();
    let new_tip = *new_order.last().unwrap();

    let mut parent = tree.parent.clone();
    parent[new_order[0]] = external_parent;
    for i in 1..new_order.len() {
        parent[new_order[i]] = new_order[i - 1];
    }
    // Distal children of the old tip (not on the unary path) follow the new tip.
    let on_path: std::collections::HashSet<usize> = old_path.iter().copied().collect();
    for &c in &children[old_tip] {
        if !on_path.contains(&c) {
            parent[c] = new_tip;
        }
    }
    MutationTree {
        n_mutations: tree.n_mutations,
        parent,
    }
}

/// Run the MCMC mutation-tree search (starting from an NJ-tree-derived guess
/// when one can be built), attach every read to its best-fitting node, and
/// write all SCITE output files with the given `output_prefix`.
///
/// Returns the SCITE-cleaned reads deduplicated into haplotypes — the same
/// matrix that backs the tips of `<output_prefix>.read_lineage.nwk`. The caller
/// needs exactly this to write `<prefix>.cleaned_haplotype_map.tsv`; handing it
/// back is what lets `lineage::start` skip re-reading `cleaned_matrix.csv` off
/// disk and deduplicating it a second time.
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
) -> Result<HaplotypeMatrix> {
    if binary.variants.is_empty() {
        anyhow::bail!(
            "SCITE requires at least 1 variant (got 0); \
             relax the --min-hf / --max-hf thresholds."
        );
    }
    // One variant is a valid, if degenerate, run: the tree can only be
    // root → variant, so there is no topology to search and no variant pair to
    // score. Every output file is still produced (the two co-occurrence tables
    // are header-only), which keeps single-variant samples in the same reporting
    // format as the rest.
    let single_variant = binary.variants.len() == 1;
    if single_variant {
        info!(
            "[SCITE] Only 1 variant ({}) passed filtering: the mutation tree is \
             fixed (root → variant), so the MCMC search and unary-path polish are \
             skipped and the pairwise co-occurrence tables will be empty.",
            binary.variants[0]
        );
    }
    let dist = lineage::hamming_distance_matrix(hap_matrix);
    let initial_tree = lineage::neighbor_joining(&dist, hap_matrix)
        .ok()
        .map(|nj_tree| from_nj_tree(hap_matrix, &nj_tree));
    info!(
        "[SCITE] Initial tree: {}",
        if initial_tree.is_some() { "derived from Neighbor-Joining haplotype tree" } else { "random (NJ tree unavailable)" }
    );

    // Overlap is not transitive (a deletion can overlap two substitutions that
    // do not overlap each other), so report the constrained pairs rather than
    // grouping into sites. Index order keeps the log deterministic.
    let spans = variant_spans(&binary.variants);
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

    let rates = ErrorRates { fp_rate, fn_rate };
    let (tree, ll) = run_mcmc_multichain(binary, &rates, n_iterations, n_chains, initial_tree.as_ref(), seed);
    info!("[SCITE] Best tree log-likelihood: {ll:.3}");

    let tree = polish_unary_path_order(&tree, binary, &rates);
    let tree = enforce_position_exclusivity(&tree, &spans);
    let ll_polish = tree_log_likelihood(binary, &tree, &rates);
    info!("[SCITE] After unary-path polish log-likelihood: {ll_polish:.3}");

    let cleaned = attach_all_reads(binary, &tree, &rates);

    // Per-edge order evidence: how much likelihood pins each mutation above its
    // parent, and the read counts behind it. Reported rather than acted on — a
    // near-tie means the pipeline had to pick one of several equally good orders.
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

    write_cleaned_matrix(&cleaned, &format!("{output_prefix}.cleaned_matrix.csv"))?;
    write_variant_cooccurrence(&cleaned, &format!("{output_prefix}.variant_cooccurrence.tsv"))?;
    write_raw_variant_cooccurrence(binary, &format!("{output_prefix}.raw_variant_cooccurrence.tsv"))?;
    write_molecule_summary(&cleaned, &format!("{output_prefix}.molecule_summary.tsv"))?;
    write_mutation_tree(&tree, &cleaned, &evidence, &format!("{output_prefix}.mutation_tree.tsv"))?;

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
    write_read_lineage_newick(
        &tree,
        &cleaned_hap_matrix,
        &rates,
        &evidence,
        &format!("{output_prefix}.read_lineage.nwk"),
    )?;

    Ok(cleaned_hap_matrix)
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
        assert_eq!(tree.ancestor_mask(0), Mask::from_bits(&[0, 1, 2]));
        assert_eq!(tree.ancestor_mask(1), Mask::from_bits(&[1, 2]));
        assert_eq!(tree.ancestor_mask(2), Mask::from_bits(&[2]));
        assert_eq!(tree.ancestor_mask(3), Mask::new());
    }

    #[test]
    fn descendant_mask_excludes_self_and_includes_transitive_children() {
        let tree = chain_tree();
        assert_eq!(tree.descendant_mask(2), Mask::from_bits(&[0, 1]));
        assert_eq!(tree.descendant_mask(1), Mask::from_bits(&[0]));
        assert_eq!(tree.descendant_mask(0), Mask::new());
        assert_eq!(tree.descendant_mask(3), Mask::from_bits(&[0, 1, 2]));
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
        let ancestor_mask = Mask::from_bits(&[0, 1]);
        let rates = ErrorRates { fp_rate: 0.1, fn_rate: 0.2 };

        let expected = (1.0 - rates.fn_rate).ln() + rates.fn_rate.ln();
        let actual = attachment_log_likelihood(&profile, &ancestor_mask, &rates);
        assert!((actual - expected).abs() < 1e-12);
    }

    #[test]
    fn attachment_log_likelihood_scores_true_negative_and_false_positive() {
        let profile = vec![Some(0), Some(1)];
        let ancestor_mask = Mask::new();
        let rates = ErrorRates { fp_rate: 0.1, fn_rate: 0.2 };

        let expected = (1.0 - rates.fp_rate).ln() + rates.fp_rate.ln();
        let actual = attachment_log_likelihood(&profile, &ancestor_mask, &rates);
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
    fn best_attachment_prefers_the_fewest_mutations_among_equally_likely_nodes() {
        // Chain root(2) -> A(0) -> B(1). The read is uncovered at A and an
        // explicit ref at B, so ROOT and A explain it identically well. ROOT is
        // the parsimonious choice; attaching to A would invent an alt call at A.
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.1 };
        let profile = vec![None, Some(0)];

        let (node, ll) = best_attachment(&profile, &tree, &rates);
        assert_eq!(node, tree.root());
        // The likelihood is unchanged by the tie-break: one true-negative call.
        assert!((ll - (1.0 - rates.fp_rate).ln()).abs() < 1e-12);
    }

    #[test]
    fn best_attachment_of_a_fully_uncovered_read_is_the_root() {
        // No observed call at all -> every node has log-likelihood 0. The read
        // carries no evidence, so it must not be assigned any mutation.
        let tree = MutationTree { n_mutations: 3, parent: vec![3, 0, 1, 3] };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.1 };
        let profile = vec![None, None, None];

        let (node, ll) = best_attachment(&profile, &tree, &rates);
        assert_eq!(node, tree.root());
        assert_eq!(ll, 0.0);
    }

    #[test]
    fn attach_all_reads_does_not_fabricate_alt_calls_for_uncovered_sites() {
        // r_partial is uncovered at A and an explicit ref at B: attaching it
        // below A would write A=1 into the cleaned matrix out of nothing.
        let tree = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.05 };
        let matrix = BinaryMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            reads: vec!["r_both".to_string(), "r_partial".to_string()],
            data: vec![vec![Some(1), None], vec![Some(1), Some(0)]],
        };

        let cleaned = attach_all_reads(&matrix, &tree, &rates);

        assert_eq!(cleaned.attachment[1], tree.root());
        assert_eq!(cleaned.data[0], vec![1, 0]);
        assert_eq!(cleaned.data[1], vec![1, 0]);
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

    fn span(start: u64, end: u64) -> VariantSpan {
        VariantSpan { start, end }
    }

    #[test]
    fn variant_spans_parses_ref_interval_and_falls_back_uniquely() {
        let variants = vec![
            "m.310T>C".to_string(),      // SNP: one base
            "m.310T>TC".to_string(),     // insertion anchored at 310
            "m.310TAA>T".to_string(),    // deletion of 311-312, anchored at 310
            "m.8274C>T".to_string(),
            "m.5->AC".to_string(),       // Himito's "no anchor base" REF placeholder
            "m.7A>-".to_string(),         // ...and the ALT-side placeholder
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
        // The gap a start-position-only check misses: 311 lies inside the
        // deletion's REF span but shares no start coordinate with it.
        let variants = vec!["m.310TAA>T".to_string(), "m.311A>G".to_string()];
        let s = variant_spans(&variants);
        assert!(s[0].overlaps(&s[1]));

        // ...while the base just past the deletion is untouched and free to
        // sit anywhere on the same lineage.
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
    }

    #[test]
    fn violates_position_exclusivity_flags_a_deletion_over_its_deleted_bases() {
        // m.310TAA>T [310,313) with m.311A>G [311,312) beneath it.
        let spans = vec![span(310, 313), span(311, 312)];
        let chain = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        assert!(violates_position_exclusivity(&chain, &spans));
    }

    #[test]
    fn enforce_position_exclusivity_splits_overlapping_alleles_to_siblings() {
        // 16126(0) -> 310T>TC(1) -> 310T>C(2): the impossible chain.
        let spans = vec![span(16126, 16127), span(310, 311), span(310, 311)];
        let chain = MutationTree { n_mutations: 3, parent: vec![3, 0, 1, 3] };
        let fixed = enforce_position_exclusivity(&chain, &spans);
        // 310T>C is lifted to become a sibling of 310T>TC, both under 16126.
        assert_eq!(fixed.parent, vec![3, 0, 0, 3]);
        assert!(!violates_position_exclusivity(&fixed, &spans));
    }

    #[test]
    fn enforce_position_exclusivity_resolves_three_overlapping_alleles() {
        // 500(0) -> 310(1) -> 310(2) -> 310(3): needs more than one pass, since
        // lifting node 2 leaves node 3 still under it.
        let spans = vec![
            span(500, 501),
            span(310, 311),
            span(310, 311),
            span(310, 311),
        ];
        let chain = MutationTree { n_mutations: 4, parent: vec![4, 0, 1, 2, 4] };
        let fixed = enforce_position_exclusivity(&chain, &spans);
        assert_eq!(fixed.parent, vec![4, 0, 0, 0, 4]);
        assert!(!violates_position_exclusivity(&fixed, &spans));
    }

    #[test]
    fn enforce_position_exclusivity_always_returns_an_acyclic_tree() {
        // Heavy overlap (12 variants over 4 loci, one of them a wide deletion)
        // is the stress case for the lift-to-sibling loop.
        let mut spans = Vec::new();
        for i in 0..12u64 {
            let locus = (i % 4) * 10;
            // Every fourth variant is a 5-base deletion straddling its neighbours.
            let width = if i % 4 == 0 { 5 } else { 1 };
            spans.push(span(locus, locus + width));
        }
        for seed in 0..500u64 {
            let mut rng = StdRng::seed_from_u64(seed);
            let tree = MutationTree::random(12, &mut rng);
            let fixed = enforce_position_exclusivity(&tree, &spans);
            assert_valid_tree(&fixed);
            assert!(
                !violates_position_exclusivity(&fixed, &spans),
                "seed {seed} left an overlapping pair on one path"
            );
        }
    }

    #[test]
    fn from_nj_tree_breaks_an_overlapping_chain_off_the_clonality_ranking() {
        // Two haplotypes: one carrying m.310T>C, one carrying m.310T>TC. The
        // clonality chain would stack them; the two alleles must end up on
        // separate branches instead.
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["m.310T>C".to_string(), "m.310T>TC".to_string()],
            haplotypes: vec![
                lineage::Haplotype {
                    id: "H0".to_string(),
                    profile: vec![Some(1), Some(0)],
                    reads: vec!["r0".to_string()],
                    count: 1,
                },
                lineage::Haplotype {
                    id: "H1".to_string(),
                    profile: vec![Some(0), Some(1)],
                    reads: vec!["r1".to_string()],
                    count: 1,
                },
            ],
        };
        let dist = lineage::hamming_distance_matrix(&hap_matrix);
        let nj_tree = lineage::neighbor_joining(&dist, &hap_matrix).unwrap();
        let tree = from_nj_tree(&hap_matrix, &nj_tree);
        assert!(!violates_position_exclusivity(
            &tree,
            &variant_spans(&hap_matrix.variants)
        ));
    }

    #[test]
    fn run_mcmc_never_places_same_position_alleles_on_one_path() {
        // m.16126T>C is the shared ancestor; m.310T>TC and m.310T>C are its two
        // mutually-exclusive descendants (never co-observed on a read).
        let variants = vec![
            "m.16126T>C".to_string(),
            "m.310T>TC".to_string(),
            "m.310T>C".to_string(),
        ];
        let mut anc = Vec::new(); // 16126
        let mut ins = Vec::new(); // 310T>TC
        let mut sub = Vec::new(); // 310T>C
        let mut push = |a: u8, i: u8, s: u8, n: usize, av: &mut Vec<Option<u8>>, iv: &mut Vec<Option<u8>>, sv: &mut Vec<Option<u8>>| {
            for _ in 0..n {
                av.push(Some(a));
                iv.push(Some(i));
                sv.push(Some(s));
            }
        };
        push(1, 1, 0, 6, &mut anc, &mut ins, &mut sub); // 16126 + insertion
        push(1, 0, 1, 6, &mut anc, &mut ins, &mut sub); // 16126 + substitution
        push(1, 0, 0, 3, &mut anc, &mut ins, &mut sub); // 16126 only
        push(0, 0, 0, 3, &mut anc, &mut ins, &mut sub); // germline

        let n_reads = anc.len();
        let matrix = BinaryMatrix {
            variants: variants.clone(),
            reads: (0..n_reads).map(|i| format!("r{i}")).collect(),
            data: vec![anc, ins, sub],
        };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.1 };
        let mut rng = StdRng::seed_from_u64(11);
        let (tree, _ll) = run_mcmc(&matrix, &rates, 3000, None, &mut rng);

        let spans = variant_spans(&variants);
        assert!(
            !violates_position_exclusivity(&tree, &spans),
            "the two position-310 alleles must never lie on one root-to-leaf path"
        );
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
    fn write_variant_cooccurrence_reports_pairwise_ld_and_fisher_stats() {
        let matrix = small_cleaned_matrix();
        let path = std::env::temp_dir().join("himito_test_write_cooccurrence.tsv");
        write_variant_cooccurrence(&matrix, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        let mut lines = content.lines();
        assert_eq!(
            lines.next().unwrap(),
            "variant_a\tvariant_b\tn_co\tn11\tn10\tn01\tn00\td_prime\tr2\todds_ratio\tfisher_p\tfisher_q"
        );
        // r1=[A,B]=[1,0], r2=[A,B]=[0,1] -> n11=1 (r1), n01=1 (r2)
        let row = lines.next().unwrap();
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(&fields[0..7], &["A", "B", "2", "1", "0", "1", "0"]);
        assert_eq!(lines.next(), None);
    }

    #[test]
    fn fisher_greater_matches_known_2x2_table() {
        // classic tea-tasting-style table with a clear enrichment in n11
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
        let matrix = BinaryMatrix {
            variants: vec!["A".to_string(), "B".to_string()],
            reads: vec!["r1".to_string(), "r2".to_string(), "r3".to_string()],
            data: vec![
                vec![Some(1), Some(0), Some(1)],
                vec![None, Some(1), Some(1)],
            ],
        };
        let path = std::env::temp_dir().join("himito_test_write_raw_cooccurrence.tsv");
        write_raw_variant_cooccurrence(&matrix, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        let row = content.lines().nth(1).unwrap();
        let fields: Vec<&str> = row.split('\t').collect();
        // r1 dropped (B uncovered): n_co=2, n11=1 (r3), n01=1 (r2)
        assert_eq!(&fields[0..7], &["A", "B", "2", "1", "0", "1", "0"]);
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
        let evidence = vec![
            EdgeEvidence { support_ll: None, counts: None },
            EdgeEvidence { support_ll: Some(4.25), counts: Some((3, 7)) },
        ];
        let path = std::env::temp_dir().join("himito_test_write_mutation_tree.tsv");
        write_mutation_tree(&tree, &matrix, &evidence, path.to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(
            content,
            "node_id\tvariant\tparent_id\tparent_variant\tn_reads_attached\t\
             order_support_ll\tn_child_without_parent\tn_parent_without_child\n\
             0\tA\t2\tROOT\t1\tNA\tNA\tNA\n\
             1\tB\t0\tA\t1\t4.250\t3\t7\n\
             2\tROOT\t2\tROOT\t0\tNA\tNA\tNA\n"
        );
    }

    #[test]
    fn run_scite_pipeline_handles_more_than_63_variants() {
        // 65 variants exceeds the old u64-bitmask limit. Reads form a nested
        // chain: read j carries variants 0..j, so every variant has both
        // present and absent reads and the deepest mutations live on bits ≥64.
        let n_var = 65;
        let n_reads = n_var + 1;
        let variants: Vec<String> = (0..n_var).map(|i| format!("m.{i}")).collect();
        let reads: Vec<String> = (0..n_reads).map(|j| format!("r{j}")).collect();
        let data: Vec<Vec<Option<u8>>> = (0..n_var)
            .map(|i| (0..n_reads).map(|j| Some(u8::from(j > i))).collect())
            .collect();
        let matrix = BinaryMatrix { variants, reads, data };
        let hap_matrix = lineage::deduplicate(&matrix, 1);

        let prefix = std::env::temp_dir()
            .join("himito_test_scite_pipeline_65var")
            .to_str()
            .unwrap()
            .to_string();

        run_scite_pipeline(&matrix, &hap_matrix, 0.01, 0.1, 50, 1, 7, 1, &prefix).unwrap();

        let nwk = format!("{prefix}.read_lineage.nwk");
        assert!(std::path::Path::new(&nwk).exists(), "expected {nwk} to exist");
        for suffix in [
            ".cleaned_matrix.csv",
            ".variant_cooccurrence.tsv",
            ".raw_variant_cooccurrence.tsv",
            ".molecule_summary.tsv",
            ".mutation_tree.tsv",
            ".read_lineage.nwk",
        ] {
            std::fs::remove_file(format!("{prefix}{suffix}")).ok();
        }
    }

    #[test]
    fn run_scite_pipeline_handles_a_single_variant() {
        // One variant: the only possible tree is root → m.A, the co-occurrence
        // tables have no pair to report, and the swap proposals (which need two
        // distinct mutation nodes) must never be reached.
        let matrix = BinaryMatrix {
            variants: vec!["m.100A>G".to_string()],
            reads: (0..10).map(|i| format!("r{i}")).collect(),
            data: vec![vec![Some(1); 7]
                .into_iter()
                .chain(vec![Some(0); 2])
                .chain([None])
                .collect()],
        };
        let hap_matrix = lineage::deduplicate(&matrix, 1);

        let prefix = std::env::temp_dir()
            .join("himito_test_scite_pipeline_1var")
            .to_str()
            .unwrap()
            .to_string();

        run_scite_pipeline(&matrix, &hap_matrix, 0.01, 0.1, 500, 2, 123, 1, &prefix).unwrap();

        let tree_tsv = std::fs::read_to_string(format!("{prefix}.mutation_tree.tsv")).unwrap();
        // The mutation hangs directly off the root, with no order to resolve.
        assert!(
            tree_tsv.contains("0\tm.100A>G\t1\tROOT\t"),
            "expected m.100A>G under ROOT, got:\n{tree_tsv}"
        );
        // Header-only co-occurrence tables: a single variant forms no pair.
        for suffix in [".variant_cooccurrence.tsv", ".raw_variant_cooccurrence.tsv"] {
            let content = std::fs::read_to_string(format!("{prefix}{suffix}")).unwrap();
            assert_eq!(content.lines().count(), 1, "expected header only in {suffix}");
        }
        let nwk = std::fs::read_to_string(format!("{prefix}.read_lineage.nwk")).unwrap();
        assert!(nwk.contains("mutation=m.100A>G"), "got newick: {nwk}");

        for suffix in [
            ".cleaned_matrix.csv",
            ".variant_cooccurrence.tsv",
            ".raw_variant_cooccurrence.tsv",
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
    fn propose_move_leaves_a_single_mutation_tree_unchanged() {
        // Guards the infinite rejection loop in `sample_two_distinct`: with one
        // mutation the swap moves have no second node to draw, so `propose_move`
        // must return the current tree instead of proposing one.
        let tree = MutationTree { n_mutations: 1, parent: vec![1, 1] };
        let mut rng = StdRng::seed_from_u64(4);
        for _ in 0..100 {
            let (proposal, correction) = propose_move(&tree, &mut rng);
            assert_eq!(proposal.parent, tree.parent);
            assert_eq!(correction, 1.0);
        }
    }

    #[test]
    fn run_scite_pipeline_writes_all_output_files() {
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
            ".raw_variant_cooccurrence.tsv",
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
        // mA hangs from the root (no order to support); the rest carry evidence.
        let evidence = vec![
            EdgeEvidence { support_ll: None, counts: None },
            EdgeEvidence { support_ll: Some(12.5), counts: Some((0, 9)) },
            EdgeEvidence { support_ll: Some(0.0), counts: Some((2, 2)) },
            EdgeEvidence { support_ll: Some(3.0), counts: Some((1, 4)) },
        ];
        write_read_lineage_newick(&tree, &hap_matrix, &rates, &evidence, path.to_str().unwrap())
            .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();

        // No haplotype carries zero mutations, so the ancestral root state is
        // materialised as a synthetic zero-read tip (Hroot) under ROOT. Each
        // branch reports the per-mutation order support alongside the mutations,
        // so a collapsed run of mutations whose order is arbitrary is visible.
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
        use lineage::{Haplotype, HaplotypeMatrix};

        // root(1) -> A(0). The all-zero haplotype (H0000, 5 reads) best-attaches
        // to the root and supplies the ancestral tip's id/count directly — no
        // synthetic Hroot, and its reads are not double-counted elsewhere.
        let tree = MutationTree { n_mutations: 1, parent: vec![1, 1] };
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.05 };
        let hap_matrix = HaplotypeMatrix {
            variants: vec!["mA".into()],
            haplotypes: vec![
                Haplotype { id: "H0000".into(), profile: vec![Some(0)], reads: vec![], count: 5 },
                Haplotype { id: "H0001".into(), profile: vec![Some(1)], reads: vec![], count: 3 },
            ],
        };

        let path = std::env::temp_dir().join("himito_test_read_lineage_root_tip.nwk");
        let evidence = vec![EdgeEvidence { support_ll: None, counts: None }];
        write_read_lineage_newick(&tree, &hap_matrix, &rates, &evidence, path.to_str().unwrap())
            .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(
            content,
            "(H0001_n3:1[&&NHX:mutation=mA:order_support=NA:reads=3],\
             H0000_n5:0[&&NHX:reads=5])ROOT;\n"
        );
    }

    #[test]
    fn polish_unary_path_order_restores_af_preferred_order() {
        // Higher AF on A, lower on B; reversed tree has B under ROOT, A under B.
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.05 };
        let matrix = BinaryMatrix {
            variants: vec!["A".into(), "B".into()],
            reads: (0..10).map(|i| format!("r{i}")).collect(),
            // 6 A-only, 4 both, 0 B-only → AF(A)=1.0, AF(B)=0.4
            data: vec![
                vec![Some(1); 10],
                vec![Some(0); 6].into_iter().chain(vec![Some(1); 4]).collect(),
            ],
        };
        let reversed = MutationTree { n_mutations: 2, parent: vec![1, 2, 2] };
        let polished = polish_unary_path_order(&reversed, &matrix, &rates);
        // A under ROOT, B under A
        assert_eq!(polished.parent[0], 2);
        assert_eq!(polished.parent[1], 0);
    }

    #[test]
    fn polish_unary_path_order_is_noop_when_order_already_matches_af() {
        // Correct A→B already matches AF order; polish must leave parents unchanged.
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.001 };
        let matrix = BinaryMatrix {
            variants: vec!["A".into(), "B".into()],
            reads: (0..10).map(|i| format!("r{i}")).collect(),
            data: vec![
                vec![Some(1); 10],
                vec![Some(0); 6].into_iter().chain(vec![Some(1); 4]).collect(),
            ],
        };
        let correct = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };
        let polished = polish_unary_path_order(&correct, &matrix, &rates);
        assert_eq!(polished.parent, correct.parent);
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

    #[test]
    fn path_order_uses_jointly_covered_counts_not_coverage_biased_allele_frequency() {
        // A is ancestral: among the 4 reads covering both sites it appears once
        // without B and never the other way round. But A is only covered in
        // those 4 reads (AF 2/4 = 0.50) while B is covered in all 20
        // (AF 12/20 = 0.60), so a global-AF sort puts B first — the per-variant
        // denominators are not comparable. The coverage-matched pairwise counts
        // (n10 = 1, n01 = 0) put A first.
        let matrix = matrix_from_rows(
            &["A", "B"],
            &[
                "1100................",
                "10001111111111100000",
            ],
        );
        assert!(allele_frequency(&matrix, 0) < allele_frequency(&matrix, 1));
        assert_eq!(order_path_by_orientation(&[0, 1], &matrix), vec![0, 1]);
        assert_eq!(order_path_by_orientation(&[1, 0], &matrix), vec![0, 1]);
    }

    #[test]
    fn path_order_falls_back_to_allele_frequency_when_joint_counts_are_tied() {
        // A and B never appear apart among jointly covered reads (n10 = n01 = 0),
        // so the pairwise statistic is uninformative. The marginal counts still
        // rank A ancestral, using reads that cover only one of the two sites.
        let matrix = matrix_from_rows(
            &["A", "B"],
            &[
                "11001111",
                "1100....",
            ],
        );
        assert_eq!(order_path_by_orientation(&[1, 0], &matrix), vec![0, 1]);
    }

    #[test]
    fn polish_reorders_a_path_that_a_global_af_sort_leaves_untouched() {
        // Same coverage-biased pair as above, with the tree already in the wrong
        // order (B ancestral). A global-AF sort agrees with the wrong order and
        // therefore proposes nothing at all.
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.05 };
        let matrix = matrix_from_rows(
            &["A", "B"],
            &[
                "1100................",
                "10001111111111100000",
            ],
        );
        let wrong = MutationTree { n_mutations: 2, parent: vec![1, 2, 2] };
        let polished = polish_unary_path_order(&wrong, &matrix, &rates);
        assert_eq!(polished.parent[0], 2, "A should hang from the root");
        assert_eq!(polished.parent[1], 0, "B should hang from A");
    }

    #[test]
    fn polish_keeps_climbing_after_the_first_candidate_ordering() {
        // Starting chain ROOT -> 2 -> 1 -> 0 (LL -15.706). The orientation sort
        // proposes 2 -> 0 -> 1, which is an improvement (LL -9.716) but not the
        // best ordering: 0 -> 2 -> 1 scores -6.722 and is the exhaustive ML over
        // all six permutations. A single accept/reject on one candidate stops at
        // the first improvement; adjacent swaps reach the optimum.
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.05 };
        let matrix = matrix_from_rows(
            &["m0", "m1", "m2"],
            &["10010001111", "0.101.0010.", "101010011.."],
        );
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
        // All three variants are carried as a block, so no read attaches to an
        // interior node and every one of the six orderings has the same
        // likelihood. The climb must settle on the orientation-preferred order
        // (index order, since every statistic ties here) instead of cycling.
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.05 };
        let matrix = matrix_from_rows(
            &["A", "B", "C"],
            &["111111110000", "111111110000", "111111110000"],
        );
        let flipped = MutationTree { n_mutations: 3, parent: vec![1, 2, 3, 3] };

        let polished = polish_unary_path_order(&flipped, &matrix, &rates);

        assert_eq!(polished.parent[0], 3);
        assert_eq!(polished.parent[1], 0);
        assert_eq!(polished.parent[2], 1);
        let before = tree_log_likelihood(&matrix, &flipped, &rates);
        let after = tree_log_likelihood(&matrix, &polished, &rates);
        assert!((after - before).abs() < 1e-9, "flat path: {before} vs {after}");
    }

    #[test]
    fn edge_order_support_is_flat_when_no_read_separates_two_mutations() {
        // Same block-carried data: nothing distinguishes A above B from B above
        // A, so the edge carries no evidence at all.
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.05 };
        let matrix = matrix_from_rows(
            &["A", "B", "C"],
            &["111111110000", "111111110000", "111111110000"],
        );
        let chain = MutationTree { n_mutations: 3, parent: vec![3, 0, 1, 3] };

        assert_eq!(edge_order_support(&chain, &matrix, &rates, 0), None, "parent is the root");
        let s1 = edge_order_support(&chain, &matrix, &rates, 1).unwrap();
        let s2 = edge_order_support(&chain, &matrix, &rates, 2).unwrap();
        assert!(s1.abs() < 1e-9, "expected no evidence, got {s1}");
        assert!(s2.abs() < 1e-9, "expected no evidence, got {s2}");
    }

    #[test]
    fn edge_order_support_is_large_when_reads_pin_the_order() {
        // Six reads carry A without B, none carry B without A, so putting B
        // above A would have to explain all six as dropouts at B.
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.05 };
        let matrix = matrix_from_rows(&["A", "B"], &["1111111111", "1111000000"]);
        let chain = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };

        let support = edge_order_support(&chain, &matrix, &rates, 1).unwrap();
        assert!(support > 5.0, "expected decisive evidence, got {support}");
    }

    #[test]
    fn order_resolved_flags_a_near_tie_even_when_the_likelihood_looks_confident() {
        // Numbers taken from the ont-r10 eval set. The inverted edge carries ~6
        // nats of likelihood — which reads as confident — but its count margin is
        // one read (sign-test p = 0.5), so the order is not actually resolved.
        let inverted = EdgeEvidence { support_ll: Some(5.989), counts: Some((8, 9)) };
        assert_eq!(inverted.order_resolved(), Some(false));

        // A flat likelihood is unresolved whatever the counts say.
        let flat = EdgeEvidence { support_ll: Some(0.0), counts: Some((4, 5)) };
        assert_eq!(flat.order_resolved(), Some(false));

        // Genuinely supported edges from the same tree: lopsided counts.
        for counts in [(3usize, 28usize), (2, 29), (7, 16), (11, 23)] {
            let ev = EdgeEvidence { support_ll: Some(26.9), counts: Some(counts) };
            assert_eq!(ev.order_resolved(), Some(true), "{counts:?} should be resolved");
        }

        // No mutation above: there is no order to resolve.
        let root_edge = EdgeEvidence { support_ll: None, counts: None };
        assert_eq!(root_edge.order_resolved(), None);
    }

    #[test]
    fn edge_evidence_reports_the_counts_behind_each_order() {
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.05 };
        let matrix = matrix_from_rows(&["A", "B"], &["110", "101"]);
        let chain = MutationTree { n_mutations: 2, parent: vec![2, 0, 2] };

        let evidence = edge_evidence(&chain, &matrix, &rates);
        assert_eq!(evidence[0].counts, None, "A hangs from the root");
        assert_eq!(evidence[0].support_ll, None);
        // B without A: read 3. A without B: read 2.
        assert_eq!(evidence[1].counts, Some((1, 1)));
        assert!(evidence[1].support_ll.is_some());
    }

    #[test]
    fn polish_unary_path_order_moves_distal_children_to_new_tip() {
        // Path ROOT→0→1 with distal children 2,3 under tip 1. AF prefers 1 ancestral to 0.
        let rates = ErrorRates { fp_rate: 0.01, fn_rate: 0.05 };
        let matrix = BinaryMatrix {
            variants: vec!["m0".into(), "m1".into(), "m2".into(), "m3".into()],
            reads: (0..8).map(|i| format!("r{i}")).collect(),
            data: vec![
                // m0: present in reads 4..7 only → lower AF
                vec![Some(0); 4].into_iter().chain(vec![Some(1); 4]).collect(),
                // m1: present in all reads 0..7 → higher AF (ancestral)
                vec![Some(1); 8],
                vec![Some(0); 6].into_iter().chain(vec![Some(1); 2]).collect(),
                vec![Some(0); 4].into_iter().chain(vec![Some(1); 2]).chain(vec![Some(0); 2]).collect(),
            ],
        };
        // parent: 0→root(4), 1→0, 2→1, 3→1
        let tree = MutationTree { n_mutations: 4, parent: vec![4, 0, 1, 1, 4] };
        let polished = polish_unary_path_order(&tree, &matrix, &rates);
        // After reorder: 1 under ROOT, 0 under 1, distal 2 and 3 under new tip 0
        assert_eq!(polished.parent[1], 4);
        assert_eq!(polished.parent[0], 1);
        assert_eq!(polished.parent[2], 0);
        assert_eq!(polished.parent[3], 0);
    }

    #[test]
    fn polish_fixes_fully_flipped_ont_eval_tree() {
        let mat = std::env::var("HIMITO_TEST_MATRIX").unwrap_or_else(|_| {
            format!("{}/lineage_eval/lineage_sim/tmp/eval_ont_r10/himito/sim.matrix.csv", env!("CARGO_MANIFEST_DIR"))
        });
        let vcf = std::env::var("HIMITO_TEST_VCF").unwrap_or_else(|_| {
            format!("{}/lineage_eval/lineage_sim/tmp/eval_ont_r10/himito/sim.vcf", env!("CARGO_MANIFEST_DIR"))
        });
        let hf = lineage::HfFilter::FromVcf(lineage::parse_vcf(&vcf, 0.01, 0.95).unwrap());
        let binary = lineage::load_and_filter_matrix(&mat, &hf, 2, 1).unwrap();
        let n = binary.variants.len();
        let root = n;
        let ix = |name: &str| {
            binary.variants.iter().position(|v| v == name).unwrap_or_else(|| panic!("missing {name} in {:?}", binary.variants))
        };
        let mut parent = vec![root; n + 1];
        parent[root] = root;
        let flipped = [
            ("m.12477T>G", "ROOT"),
            ("m.11749A>T", "m.12477T>G"),
            ("m.7102T>C", "m.11749A>T"),
            ("m.5967T>C", "ROOT"),
            ("m.12183G>C", "m.5967T>C"),
            ("m.15258A>C", "ROOT"),
            ("m.1700T>A", "m.15258A>C"),
            ("m.14577T>A", "m.15258A>C"),
            ("m.2358A>C", "m.14577T>A"),
            ("m.8171T>A", "ROOT"),
        ];
        for (child, par) in flipped {
            let c = ix(child);
            parent[c] = if par == "ROOT" { root } else { ix(par) };
        }
        let tree = MutationTree { n_mutations: n, parent };
        let rates = ErrorRates { fp_rate: 0.001, fn_rate: 0.05 };
        let ll0 = tree_log_likelihood(&binary, &tree, &rates);
        let polished = polish_unary_path_order(&tree, &binary, &rates);
        let ll1 = tree_log_likelihood(&binary, &polished, &rates);
        eprintln!("n_vars={n} ll0={ll0:.3} ll1={ll1:.3} delta={:.3}", ll1 - ll0);
        assert_eq!(polished.parent[ix("m.5967T>C")], ix("m.12183G>C"));
        assert_eq!(polished.parent[ix("m.12183G>C")], root);
        assert_eq!(polished.parent[ix("m.2358A>C")], ix("m.15258A>C"));
        assert_eq!(polished.parent[ix("m.14577T>A")], ix("m.2358A>C"));
        assert_eq!(polished.parent[ix("m.11749A>T")], root);

        // The long chain is recovered as 11749 → {12477, 7102}: 11749 ancestral to
        // both, as in the truth tree. The order *within* the 12477/7102 pair is not
        // recoverable from this data set and is deliberately not asserted:
        //   * the pairwise counts are n10 = 8 vs n01 = 9 — a one-read margin
        //     (sign-test p ≈ 0.7), so the coverage-matched evidence is a coin flip;
        //   * the truth order (12477 ancestral) needs 9 dropouts at 12477 while the
        //     reverse needs 8 at 7102, so *any* likelihood with a symmetric
        //     false-negative rate prefers the reverse, by ≈ 6 nats;
        //   * the full pipeline on this same data already returns the reverse
        //     (`sim_lineage.mutation_tree.tsv`), and did so before this polish
        //     existed — the old whole-path rule only reached the truth order from
        //     this artificially flipped starting tree, not from the MCMC's tree.
        // `edge_order_support` is what surfaces the ambiguity to the user.
        let chain = [ix("m.12477T>G"), ix("m.7102T>C")];
        let deep = if polished.parent[chain[0]] == chain[1] { chain[0] } else { chain[1] };
        let shallow = if deep == chain[0] { chain[1] } else { chain[0] };
        assert_eq!(polished.parent[shallow], ix("m.11749A>T"));
        assert_eq!(polished.parent[deep], shallow);
        let (n10, n01) = orientation_counts(&binary, chain[0], chain[1]);
        assert!(
            (n10 as i64 - n01 as i64).abs() <= 2,
            "this pair is expected to be a near-tie (n10={n10}, n01={n01})"
        );

        assert!(ll1 > ll0 + 1.0, "expected large LL gain from fixing flipped lineages");
    }
}
