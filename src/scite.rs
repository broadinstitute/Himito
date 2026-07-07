use std::io::{BufWriter, Write};
use std::fs::File;
use anyhow::{Context, Result};
use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::lineage::{self, BinaryMatrix, HaplotypeMatrix};

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
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chain_tree() -> MutationTree {
        // 0 -> 1 -> 2 -> root(3): node0's parent is 1, node1's parent is 2,
        // node2's parent is root, root's parent is itself.
        MutationTree { n_mutations: 3, parent: vec![1, 2, 3, 3] }
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
}
