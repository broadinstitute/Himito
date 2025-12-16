/// MCMC algorithm for finding optimal phylogenetic trees from binary mutation matrices.
use std::f64;
use std::sync::Arc;
use rand::Rng;
use rayon::prelude::*;
use indicatif::ProgressBar;

/// Tree and error rate combination for tracking optimal trees
#[derive(Clone, Debug)]
pub struct TreeBeta {
    pub tree: Vec<i32>,
    pub beta: f64,
}

/// Result of MCMC search
#[derive(Debug)]
pub struct MCMCResult {
    pub best_tree: Vec<i32>,
    pub best_score: f64,
    pub best_beta: f64,
    pub all_optimal_trees: Vec<TreeBeta>,
    pub n_mutations: usize,
    pub n_samples: usize,
}

/// Convert parent vector to ancestor matrix
/// Returns a matrix where [i][j] = true if i is an ancestor of j
pub fn parent_vector_to_anc_matrix(parent: &[i32], n: usize) -> Vec<Vec<bool>> {
    let mut anc_matrix = vec![vec![false; n]; n];
    let root = n as i32;
    
    for i in 0..n {
        let mut anc = i as i32;
        while anc < root {
            if parent[anc as usize] < root {
                anc_matrix[parent[anc as usize] as usize][i] = true;
            }
            anc = parent[anc as usize];
        }
    }
    
    // Each node is an ancestor of itself
    for i in 0..n {
        anc_matrix[i][i] = true;
    }
    
    anc_matrix
}

/// Get all descendants of a node
pub fn get_descendants(anc_matrix: &[Vec<bool>], node: usize, n: usize) -> Vec<usize> {
    let mut descendants = Vec::new();
    for i in 0..n {
        if anc_matrix[node][i] {
            descendants.push(i);
        }
    }
    descendants
}

/// Get all non-descendants of a node (ancestors and nodes in different branches)
pub fn get_non_descendants(anc_matrix: &[Vec<bool>], node: usize, n: usize) -> Vec<usize> {
    let mut non_descendants = Vec::new();
    for i in 0..n {
        if !anc_matrix[node][i] {
            non_descendants.push(i);
        }
    }
    non_descendants
}

/// Get child list from parent vector
pub fn get_child_list_from_parent_vector(parent: &[i32], n: usize) -> Vec<Vec<usize>> {
    let mut child_list = vec![Vec::new(); n + 1];
    for i in 0..n {
        child_list[parent[i] as usize].push(i);
    }
    child_list
}

/// Compute breadth-first traversal of tree
pub fn get_breadth_first_traversal(parent: &[i32], n: usize) -> Vec<usize> {
    let child_list = get_child_list_from_parent_vector(parent, n);
    let mut bft = vec![n]; // Start with root
    let mut i = 0;
    
    while i < bft.len() && bft.len() < n + 1 {
        let node = bft[i];
        for &child in &child_list[node] {
            if bft.len() < n + 1 {
                bft.push(child);
            }
        }
        i += 1;
    }
    
    bft
}

/// Generate random parent vector (random tree)
pub fn get_rand_parent_vec(n: usize, rng: &mut impl Rng) -> Vec<i32> {
    // Generate random Pruefer code
    let code_length = n - 1;
    let mut code = Vec::with_capacity(code_length);
    for _ in 0..code_length {
        code.push(rng.gen_range(0..=n) as i32);
    }
    
    pruefer_code_to_parent_vector(&code, n)
}

/// Convert Pruefer code to parent vector
fn pruefer_code_to_parent_vector(code: &[i32], n: usize) -> Vec<i32> {
    let node_count = n + 1;
    let mut parent = vec![0; node_count];
    let root = n as i32;
    
    // Get last occurrence of each node in code
    let mut last_occ = vec![-1; node_count + 1];
    for (i, &node) in code.iter().enumerate() {
        if node != root {
            last_occ[node as usize] = i as i32;
        }
    }
    
    // Initialize queue: nodes that are not in code
    let mut queue = vec![true; node_count + 1];
    for &node in code {
        queue[node as usize] = false;
    }
    
    let mut queue_cutter = -1;
    let mut next_idx = get_next_in_queue(&queue, 0, node_count + 1);
    
    for (i, &code_node) in code.iter().enumerate() {
        if queue_cutter >= 0 {
            parent[queue_cutter as usize] = code_node;
            queue_cutter = -1;
        } else {
            parent[next_idx] = code_node;
            next_idx = get_next_in_queue(&queue, next_idx + 1, node_count + 1);
        }
        
        if last_occ[code_node as usize] == i as i32 {
            update_queue(code_node, &mut queue, next_idx);
            queue_cutter = update_queue_cutter(code_node, &queue, next_idx);
        }
    }
    
    if queue_cutter >= 0 {
        parent[queue_cutter as usize] = root;
    } else {
        parent[next_idx] = root;
    }
    
    parent
}

fn get_next_in_queue(queue: &[bool], pos: usize, length: usize) -> usize {
    for i in pos..length {
        if queue[i] {
            return i;
        }
    }
    length
}

fn update_queue(node: i32, queue: &mut [bool], next_idx: usize) {
    if node as usize >= next_idx {
        queue[node as usize] = true;
    }
}

fn update_queue_cutter(node: i32, _queue: &[bool], next_idx: usize) -> i32 {
    if (node as usize) < next_idx {
        node
    } else {
        -1
    }
}

// ============================================================================
// Scoring Functions (Simplified: only 0->1 and 1->0 transitions)
// ============================================================================

/// Compute log scores for binary transitions only
/// Returns a 2x2 matrix: [observed][true_state]
/// observed: 0 or 1
/// true_state: 0 (absent) or 1 (present)
pub fn get_binary_log_scores(fd: f64, ad1: f64) -> Vec<Vec<f64>> {
    let mut log_scores = vec![vec![0.0; 2]; 2];
    
    // Observed 0, true 0: P(obs=0 | true=0) = 1 - fd
    log_scores[0][0] = if 1.0 - fd > 0.0 {
        (1.0 - fd).ln()
    } else {
        f64::NEG_INFINITY
    };
    
    // Observed 1, true 0: P(obs=1 | true=0) = fd
    log_scores[1][0] = if fd > 0.0 {
        fd.ln()
    } else {
        f64::NEG_INFINITY
    };
    
    // Observed 0, true 1: P(obs=0 | true=1) = ad1
    log_scores[0][1] = if ad1 > 0.0 {
        ad1.ln()
    } else {
        f64::NEG_INFINITY
    };
    
    // Observed 1, true 1: P(obs=1 | true=1) = 1 - ad1
    log_scores[1][1] = if 1.0 - ad1 > 0.0 {
        (1.0 - ad1).ln()
    } else {
        f64::NEG_INFINITY
    };
    
    log_scores
}

/// Compute root attachment score (attaching sample to root)
fn root_attachment_score(n: usize, log_scores: &[Vec<f64>], data_vector: &[i32]) -> f64 {
    let mut score = 0.0;
    for gene in 0..n {
        let obs = data_vector[gene] as usize;
        score += log_scores[obs][0]; // True state is 0 at root
    }
    score
}

/// Compute attachment scores for a sample to all nodes in the tree
fn get_attachment_scores(
    parent: &[i32],
    n: usize,
    log_scores: &[Vec<f64>],
    data_vector: &[i32],
    bft: &[usize],
) -> Vec<f64> {
    let mut attachment_score = vec![f64::NEG_INFINITY; n + 1];
    attachment_score[n] = root_attachment_score(n, log_scores, data_vector);
    
    for i in 1..bft.len() {
        let node = bft[i];
        if node < n {
            attachment_score[node] = attachment_score[parent[node] as usize];
            let obs = data_vector[node] as usize;
            attachment_score[node] -= log_scores[obs][0]; // Remove score for true=0
            attachment_score[node] += log_scores[obs][1]; // Add score for true=1
        }
    }
    
    attachment_score
}

/// Fast approximate tree scoring (max attachment per sample)
fn score_tree_fast_max(
    n: usize,
    m: usize,
    log_scores: &[Vec<f64>],
    data_matrix: &[Vec<i32>],
    parent: &[i32],
) -> f64 {
    let bft = get_breadth_first_traversal(parent, n);
    let mut tree_score = 0.0;
    
    for sample in 0..m {
        let scores = get_attachment_scores(parent, n, log_scores, &data_matrix[sample], &bft);
        let max_score = scores.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        tree_score += max_score;
    }
    
    tree_score
}

/// Accurate tree scoring (max attachment per sample)
fn score_tree_accurate_max(
    n: usize,
    m: usize,
    log_scores: &[Vec<f64>],
    data_matrix: &[Vec<i32>],
    parent: &[i32],
) -> f64 {
    let bft = get_breadth_first_traversal(parent, n);
    let mut tree_score = 0.0;
    
    for sample in 0..m {
        let scores = get_attachment_scores(parent, n, log_scores, &data_matrix[sample], &bft);
        let max_score = scores.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        tree_score += max_score;
    }
    
    tree_score
}

/// Score a tree (uses fast approximation, then accurate if close to best)
pub fn score_tree(
    n: usize,
    m: usize,
    log_scores: &[Vec<f64>],
    data_matrix: &[Vec<i32>],
    parent: &[i32],
    best_tree_log_score: f64,
) -> f64 {
    const EPSILON: f64 = 1e-12;
    let approx = score_tree_fast_max(n, m, log_scores, data_matrix, parent);
    
    if approx > best_tree_log_score - EPSILON {
        score_tree_accurate_max(n, m, log_scores, data_matrix, parent)
    } else {
        approx
    }
}

/// Accurate tree scoring (always use accurate method)
pub fn score_tree_accurate(
    n: usize,
    m: usize,
    log_scores: &[Vec<f64>],
    data_matrix: &[Vec<i32>],
    parent: &[i32],
) -> f64 {
    score_tree_accurate_max(n, m, log_scores, data_matrix, parent)
}

// ============================================================================
// Tree Move Proposals
// ============================================================================

/// Propose a new tree by applying a random move
/// Returns (proposed_tree, neighborhood_correction)
pub fn propose_new_tree(
    move_probs: &[f64],
    n: usize,
    curr_tree_anc_matrix: &[Vec<bool>],
    curr_tree_parent_vec: &[i32],
    rng: &mut impl Rng,
) -> (Vec<i32>, f64) {
    let movetype = sample_random_move(move_probs, rng);
    let mut nbh_correction = 1.0;
    
    match movetype {
        1 => {
            // Prune and re-attach
            let node_to_move = rng.gen_range(0..n);
            let possible_parents = get_non_descendants(curr_tree_anc_matrix, node_to_move, n);
            let new_parent = choose_parent(&possible_parents, n, rng);
            let prop_tree_par_vec = get_new_parent_vec_fast(
                curr_tree_parent_vec,
                node_to_move,
                new_parent,
                n,
            );
            (prop_tree_par_vec, nbh_correction)
        }
        2 => {
            // Swap two node labels
            let nodes_to_swap = sample_two_elements_without_replacement(n, rng);
            let prop_tree_par_vec = get_new_parent_vec_swap_fast(
                curr_tree_parent_vec,
                nodes_to_swap.0,
                nodes_to_swap.1,
                n,
            );
            (prop_tree_par_vec, nbh_correction)
        }
        3 => {
            // Swap two subtrees
            let nodes_to_swap = sample_two_elements_without_replacement(n, rng);
            let (node_to_move, next_node_to_move) = reorder_to_start_with_descendant(
                nodes_to_swap,
                curr_tree_anc_matrix,
            );
            
            if !curr_tree_anc_matrix[next_node_to_move][node_to_move] {
                // Nodes are in different lineages -- simple case
                let mut prop_tree_par_vec = curr_tree_parent_vec.to_vec();
                prop_tree_par_vec[node_to_move] = curr_tree_parent_vec[next_node_to_move];
                prop_tree_par_vec[next_node_to_move] = curr_tree_parent_vec[node_to_move];
                (prop_tree_par_vec, nbh_correction)
            } else {
                // Nodes are in the same lineage -- need to avoid cycles
                let mut prop_tree_par_vec = curr_tree_parent_vec.to_vec();
                prop_tree_par_vec[node_to_move] = curr_tree_parent_vec[next_node_to_move];
                
                let descendants = get_descendants(curr_tree_anc_matrix, node_to_move, n);
                let prop_tree_anc_matrix = parent_vector_to_anc_matrix(&prop_tree_par_vec, n);
                let next_descendants = get_descendants(&prop_tree_anc_matrix, next_node_to_move, n);
                
                if !descendants.is_empty() {
                    let chosen = descendants[rng.gen_range(0..descendants.len())];
                    prop_tree_par_vec[next_node_to_move] = chosen as i32;
                    nbh_correction = (descendants.len() as f64) / (next_descendants.len() as f64);
                }
                
                (prop_tree_par_vec, nbh_correction)
            }
        }
        _ => {
            // Default: return current tree
            (curr_tree_parent_vec.to_vec(), nbh_correction)
        }
    }
}

fn sample_random_move(move_probs: &[f64], rng: &mut impl Rng) -> usize {
    if move_probs.len() < 2 {
        return 1;
    }
    
    let r: f64 = rng.gen();
    let mut prob_sum = if move_probs.len() > 1 { move_probs[1] } else { 0.0 };
    
    // Start at index 1; index 0 is for changing error rate
    for i in 1..move_probs.len() - 1 {
        if r <= prob_sum {
            return i;
        }
        prob_sum += move_probs[i + 1];
    }
    
    move_probs.len() - 1
}

fn choose_parent(possible_parents: &[usize], root: usize, rng: &mut impl Rng) -> usize {
    let mut possible_parents_with_root = possible_parents.to_vec();
    possible_parents_with_root.push(root);
    let chosen_pos = rng.gen_range(0..possible_parents_with_root.len());
    possible_parents_with_root[chosen_pos]
}

fn get_new_parent_vec_fast(
    curr_tree_parent_vec: &[i32],
    node_to_move: usize,
    new_parent: usize,
    _n: usize,
) -> Vec<i32> {
    let mut prop_tree_par_vec = curr_tree_parent_vec.to_vec();
    prop_tree_par_vec[node_to_move] = new_parent as i32;
    prop_tree_par_vec
}

fn get_new_parent_vec_swap_fast(
    curr_tree_parent_vec: &[i32],
    first: usize,
    second: usize,
    n: usize,
) -> Vec<i32> {
    let mut prop_tree_par_vec = curr_tree_parent_vec.to_vec();
    
    for i in 0..n {
        if prop_tree_par_vec[i] == first as i32 && i != second {
            prop_tree_par_vec[i] = second as i32;
        } else if prop_tree_par_vec[i] == second as i32 && i != first {
            prop_tree_par_vec[i] = first as i32;
        }
    }
    
    let temp = prop_tree_par_vec[first];
    prop_tree_par_vec[first] = prop_tree_par_vec[second];
    prop_tree_par_vec[second] = temp;
    
    // Ensure tree is connected
    if prop_tree_par_vec[first] == first as i32 {
        prop_tree_par_vec[first] = second as i32;
    }
    if prop_tree_par_vec[second] == second as i32 {
        prop_tree_par_vec[second] = first as i32;
    }
    
    prop_tree_par_vec
}

fn reorder_to_start_with_descendant(
    nodes_to_swap: (usize, usize),
    curr_tree_anc_matrix: &[Vec<bool>],
) -> (usize, usize) {
    if curr_tree_anc_matrix[nodes_to_swap.0][nodes_to_swap.1] {
        // First is ancestor of second, so swap
        (nodes_to_swap.1, nodes_to_swap.0)
    } else {
        nodes_to_swap
    }
}

fn sample_two_elements_without_replacement(n: usize, rng: &mut impl Rng) -> (usize, usize) {
    let first = rng.gen_range(0..n);
    let mut second = rng.gen_range(0..n);
    while first == second {
        second = rng.gen_range(0..n);
    }
    (first, second)
}

// ============================================================================
// Tree List Management
// ============================================================================

/// Update list of optimal trees
pub fn update_tree_list(
    best_trees: &mut Vec<TreeBeta>,
    curr_tree_parent_vec: &[i32],
    n: usize,
    curr_score: f64,
    best_score: f64,
    beta: f64,
) {
    if curr_score > best_score {
        // New best score, reset list
        best_trees.clear();
        best_trees.push(TreeBeta {
            tree: curr_tree_parent_vec.to_vec(),
            beta,
        });
    } else if (curr_score - best_score).abs() < 1e-10 {
        // Same score, check if tree is duplicate
        if !is_duplicate_tree(best_trees, curr_tree_parent_vec, n) {
            best_trees.push(TreeBeta {
                tree: curr_tree_parent_vec.to_vec(),
                beta,
            });
        }
    }
}

fn is_duplicate_tree(best_trees: &[TreeBeta], new_tree: &[i32], n: usize) -> bool {
    for tree_beta in best_trees {
        let mut is_same = true;
        for i in 0..n {
            if tree_beta.tree[i] != new_tree[i] {
                is_same = false;
                break;
            }
        }
        if is_same {
            return true;
        }
    }
    false
}

// ============================================================================
// MCMC Algorithm
// ============================================================================
/// Result from a single MCMC repetition
#[derive(Clone, Debug)]
struct RepetitionResult {
    best_tree: Vec<i32>,
    best_tree_log_score: f64,
    best_score: f64,
    best_beta: f64,
    best_trees: Vec<TreeBeta>,
}

/// Run a single MCMC repetition
fn run_single_repetition(
    r: usize,
    fd: f64,
    ad1: f64,
    no_of_reps: usize,
    no_of_loops: usize,
    gamma: f64,
    move_probs: &[f64],
    n: usize,
    m: usize,
    data_matrix: &[Vec<i32>],
    log_scores: &[Vec<f64>],
    use_tree_list: bool,
    initial_tree: Option<&[i32]>,
    use_initial_tree: bool,
) -> RepetitionResult {
    let mut rng = rand::thread_rng();
    
    // Start MCMC with initial tree if provided (only for first repetition), otherwise random tree
    let mut curr_tree_parent_vec = if use_initial_tree {
        // First repetition: use initial tree if provided, otherwise random
        if let Some(init_tree) = initial_tree {
            if init_tree.len() != n {
                panic!("Initial tree length {} does not match number of mutations {}", init_tree.len(), n);
            }
            init_tree.to_vec()
        } else {
            get_rand_parent_vec(n, &mut rng)
        }
    } else {
        // Subsequent repetitions: always use random tree
        get_rand_parent_vec(n, &mut rng)
    };
    let mut curr_tree_anc_matrix = parent_vector_to_anc_matrix(&curr_tree_parent_vec, n);
    let curr_beta = ad1;
    
    let mut curr_tree_log_score = score_tree_accurate(
        n,
        m,
        log_scores,
        data_matrix,
        &curr_tree_parent_vec,
    );
    
    let mut curr_score = curr_tree_log_score;
    let mut best_tree_log_score = curr_tree_log_score;
    let mut best_score = curr_score;
    let mut best_beta = curr_beta;
    let mut best_trees: Vec<TreeBeta> = Vec::new();
    
    for it in 0..no_of_loops {
        if it % 100000 == 0 && it > 0 {
            println!(
                "MCMC repetition {}/{}, step {}: best tree score {:.4}, best overall score {:.4}",
                r + 1,
                no_of_reps,
                it,
                best_tree_log_score,
                best_score
            );
        }
        
        
        // Move changed tree
        let (prop_tree_par_vec, nbh_correction) = propose_new_tree(
            move_probs,
            n,
            &curr_tree_anc_matrix,
            &curr_tree_parent_vec,
            &mut rng,
        );
            
        let prop_tree_log_score = score_tree(
            n,
            m,
            log_scores,
            data_matrix,
            &prop_tree_par_vec,
            best_tree_log_score,
        );
            
        // Metropolis-Hastings acceptance
        let acceptance_prob = nbh_correction
            * ((prop_tree_log_score - curr_tree_log_score) * gamma).exp();
        
        if rng.gen::<f64>() < acceptance_prob {
            curr_tree_anc_matrix = parent_vector_to_anc_matrix(&prop_tree_par_vec, n);
            curr_tree_parent_vec = prop_tree_par_vec;
            curr_tree_log_score = prop_tree_log_score;
            curr_score = curr_tree_log_score;
        }
        
        
        // Update list of optimal trees
        if use_tree_list {
            update_tree_list(
                &mut best_trees,
                &curr_tree_parent_vec,
                n,
                curr_score,
                best_score,
                curr_beta,
            );
        }
        
        // Update best tree
        if curr_score > best_score {
            best_tree_log_score = curr_tree_log_score;
            best_score = curr_score;
            best_beta = curr_beta;
        }
    }
    
    // Get the best tree from this repetition
    let best_tree = if !best_trees.is_empty() {
        best_trees[0].tree.clone()
    } else {
        curr_tree_parent_vec.clone()
    };
    
    RepetitionResult {
        best_tree,
        best_tree_log_score,
        best_score,
        best_beta,
        best_trees,
    }
}

pub fn run_mcmc(
    fd: f64,
    ad1: f64,
    no_of_reps: usize,
    no_of_loops: usize,
    gamma: f64,
    move_probs: &[f64],
    n: usize,
    m: usize,
    data_matrix: &[Vec<i32>],
    _score_type: char,
    use_tree_list: bool,
    initial_tree: Option<&[i32]>,
    rng: &mut impl Rng,
) -> MCMCResult {
    const BURN_IN_PHASE: f64 = 0.25;
    let _burn_in = (no_of_loops as f64 * BURN_IN_PHASE) as usize;
    
    // Compute log scores (only binary: 0->1 and 1->0)
    let log_scores = get_binary_log_scores(fd, ad1);
    
    // Prepare data for parallel execution
    let data_matrix_ref = data_matrix;
    let move_probs_ref = move_probs;
    let initial_tree_ref = initial_tree;
    
    // Run repetitions in parallel
    let progress_bar = Arc::new(ProgressBar::new(no_of_reps as u64));
    let results: Vec<RepetitionResult> = (0..no_of_reps)
        .into_par_iter()
        .map(|r| {
            let result = run_single_repetition(
                r,
                fd,
                ad1,
                no_of_reps,
                no_of_loops,
                gamma,
                move_probs_ref,
                n,
                m,
                data_matrix_ref,
                &log_scores,
                use_tree_list,
                initial_tree_ref,
                r == 0, // Only first repetition uses initial tree
            );
            progress_bar.inc(1);
            result
        })
        .collect();
    progress_bar.finish();
    
    // Merge results from all repetitions
    let mut best_tree_log_score = f64::NEG_INFINITY;
    let mut best_score = f64::NEG_INFINITY;
    let mut best_beta = ad1;
    let mut best_trees: Vec<TreeBeta> = Vec::new();
    let mut best_tree: Vec<i32> = Vec::new();
    
    for result in results {
        // Update overall best
        if result.best_score > best_score {
            best_tree_log_score = result.best_tree_log_score;
            best_score = result.best_score;
            best_beta = result.best_beta;
            best_tree = result.best_tree.clone();
        }
        
        // Merge tree lists
        if use_tree_list {
            for tree_beta in result.best_trees {
                if tree_beta.beta == best_beta {
                    if tree_beta.tree.len() == n {
                        // Check if this tree is better or equal to current best
                        let tree_score = score_tree_accurate(
                            n,
                            m,
                            &log_scores,
                            data_matrix_ref,
                            &tree_beta.tree,
                        );
                        if tree_score > best_score {
                            // New best, reset list
                            best_trees.clear();
                            best_trees.push(tree_beta);
                            best_score = tree_score;
                            best_tree_log_score = tree_score;
                        } else if (tree_score - best_score).abs() < 1e-10 {
                            // Same score, check if duplicate
                            if !is_duplicate_tree(&best_trees, &tree_beta.tree, n) {
                                best_trees.push(tree_beta);
                            }
                        }
                    }
                }
            }
        }
    }
    
    println!("Best log score for tree: {:.4}", best_tree_log_score);
    
    if best_trees.is_empty() {
        // If no trees in list, add the best one
        best_trees.push(TreeBeta {
            tree: best_tree.clone(),
            beta: best_beta,
        });
    }
    
    // Get the best tree
    let final_best_tree = if !best_trees.is_empty() {
        best_trees[0].tree.clone()
    } else if !best_tree.is_empty() {
        best_tree
    } else {
        Vec::new()
    };
    
    MCMCResult {
        best_tree: final_best_tree,
        best_score: best_tree_log_score,
        best_beta,
        all_optimal_trees: best_trees,
        n_mutations: n,
        n_samples: m,
    }
}

/// Find optimal tree from binary matrix
/// 
/// # Arguments
/// * `data_matrix` - Binary mutation matrix (m samples x n mutations), values should be 0 or 1
/// * `fd` - False discovery rate (0→1 errors), default 0.0001
/// * `ad1` - Allelic dropout rate (1→0 errors), default 0.2
/// * `no_of_reps` - Number of MCMC repetitions, default 1
/// * `no_of_loops` - Chain length per repetition, default 100000
/// * `gamma` - Scaling factor for scores, default 1.0
/// * `move_probs` - Probabilities for moves [beta_move, prune_reattach, swap_labels, swap_subtrees]
///   Default: [0.0, 0.55, 0.4, 0.05]
/// * `score_type` - 'm' for max, 's' for sum (currently only 'm' is implemented)
/// * `use_tree_list` - Whether to maintain list of optimal trees, default true
/// * `initial_tree` - Optional initial tree to start MCMC from. If None, starts from random tree.
pub fn find_optimal_tree(
    data_matrix: &[Vec<i32>],
    fd: Option<f64>,
    ad1: Option<f64>,
    no_of_reps: Option<usize>,
    no_of_loops: Option<usize>,
    gamma: Option<f64>,
    move_probs: Option<&[f64]>,
    score_type: Option<char>,
    use_tree_list: Option<bool>,
    initial_tree: Option<&[i32]>,
    rng: &mut impl Rng,
) -> MCMCResult {
    let m = data_matrix.len();
    if m == 0 {
        panic!("Data matrix is empty");
    }
    let n = data_matrix[0].len();
    
    // Validate matrix is binary
    for row in data_matrix {
        if row.len() != n {
            panic!("Inconsistent row lengths in data matrix");
        }
        for &val in row {
            if val != 0 && val != 1 {
                panic!("Data matrix contains non-binary value: {}", val);
            }
        }
    }
    
    let fd = fd.unwrap_or(0.0001);
    let ad1 = ad1.unwrap_or(0.2);
    let no_of_reps = no_of_reps.unwrap_or(1);
    let no_of_loops = no_of_loops.unwrap_or(100000);
    let gamma = gamma.unwrap_or(1.0);
    let default_move_probs = vec![0.0, 0.55, 0.4, 0.05];
    let move_probs = move_probs.unwrap_or(&default_move_probs);
    let score_type = score_type.unwrap_or('m');
    let use_tree_list = use_tree_list.unwrap_or(true);
    
    // Validate initial tree if provided
    if let Some(init_tree) = initial_tree {
        if init_tree.len() != n {
            panic!("Initial tree length {} does not match number of taxa {}", init_tree.len(), n);
        }
        // Validate that all parent indices are valid (0 to n, where n is root)
        for (i, &parent) in init_tree.iter().enumerate() {
            if parent < 0 || parent > n as i32 {
                panic!("Invalid parent {} for node {} in initial tree (must be 0-{})", parent, i, n);
            }
            if parent == i as i32 {
                panic!("Node {} cannot be its own parent in initial tree", i);
            }
        }
    }
    
    println!("Running MCMC to find optimal tree...");
    println!("  Mutations: {}", n);
    println!("  Samples: {}", m);
    println!("  MCMC iterations: {}", no_of_loops);
    println!("  Repetitions: {}", no_of_reps);
    println!("  Error rates: fd={}, ad1={}", fd, ad1);
    println!("  Score type: {}", score_type);
    if initial_tree.is_some() {
        println!("  Starting from provided initial tree");
    } else {
        println!("  Starting from random tree");
    }
    println!();
    
    run_mcmc(
        fd,
        ad1,
        no_of_reps,
        no_of_loops,
        gamma,
        move_probs,
        n,
        m,
        data_matrix,
        score_type,
        use_tree_list,
        initial_tree,
        rng,
    )
}

/// Compute a table of log scores of observing one genotype given another.
/// Simplified version that only considers binary transitions (0→1 and 1→0).
/// 
/// # Arguments
/// * `fd` - False discovery rate (0→1 errors)
/// * `ad1` - Allelic dropout rate (1→0 errors)
/// 
/// # Returns
/// A 2x2 matrix where:
/// - First dimension (2): observed state (0 or 1)
/// - Second dimension (2): true state (0 = absent, 1 = present)
pub fn get_log_scores(fd: f64, ad1: f64) -> Vec<Vec<f64>> {
    let mut log_scores = vec![vec![0.0; 2]; 2];
    
    // Observed 0, true 0: P(obs=0 | true=0) = 1 - fd
    log_scores[0][0] = if 1.0 - fd > 0.0 {
        (1.0 - fd).ln()
    } else {
        f64::NEG_INFINITY
    };
    
    // Observed 1, true 0: P(obs=1 | true=0) = fd
    log_scores[1][0] = if fd > 0.0 {
        fd.ln()
    } else {
        f64::NEG_INFINITY
    };
    
    // Observed 0, true 1: P(obs=0 | true=1) = ad1
    log_scores[0][1] = if ad1 > 0.0 {
        ad1.ln()
    } else {
        f64::NEG_INFINITY
    };
    
    // Observed 1, true 1: P(obs=1 | true=1) = 1 - ad1
    log_scores[1][1] = if 1.0 - ad1 > 0.0 {
        (1.0 - ad1).ln()
    } else {
        f64::NEG_INFINITY
    };
    
    log_scores
}
