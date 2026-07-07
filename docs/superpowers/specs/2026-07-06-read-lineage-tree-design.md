# Read Lineage Tree Output — Design

**Date:** 2026-07-06
**Component:** `Himito lineage` / `src/scite.rs`
**Status:** Approved (pending spec review)

## Goal

Add a fifth output to the SCITE lineage pipeline: a **read lineage tree** in
Newick format, styled after `tsinfer`/`tskit` tree sequences — haplotypes are
nodes and **mutations are annotated on the branches** (not as node labels).
File: `<prefix>.read_lineage.nwk`.

## Motivation

The pipeline already emits the mutation tree (`.mutation_tree.tsv`), where each
node *is* a mutation. Users want the dual view: a tree whose **tips are observed
haplotypes** and whose **branches are labeled with the mutations acquired**, so
it loads in standard tree viewers (ete3, Dendroscope, IcyTree) and reads like a
tsinfer tree — ancestral haplotypes at internal nodes, sampled haplotypes at
tips, mutations marked on the edges where they occur.

## Semantics (tsinfer-style mapping)

The ML `MutationTree` already contains everything needed. Each internal SCITE
node is an ancestral genotype (its `ancestor_mask`), and the single mutation
acquired on the edge *into* that node is exactly the mutation the node
represents. Re-expressing that tree as a haplotype tree:

- **Tips** = observed haplotypes from `HaplotypeMatrix` (the same `H####` ids as
  `.raw_haplotype_map.tsv`). Every read in a haplotype shares one profile, so a
  haplotype attaches as a unit to a single node via
  `best_attachment(haplotype.profile, tree, rates)`.
- **Internal nodes** = ancestral haplotypes at genuine branch points, labeled
  `anc<node_id>`; the germline node is labeled `ROOT`.
- **Branches carry the mutations**, via an NHX comment:
  `[&&NHX:mutation=m.310T>TC]`.

## Unary collapse (confirmed behavior)

To match tsinfer (no pass-through internal nodes), mutation nodes that are not
real branch points are dissolved. A node's **effective children** are its
surviving child-lineages plus its attached haplotype tips (each attached
haplotype counts as one child).

- A mutation node is **kept** as an ancestral internal node (`anc<id>`) iff it
  has **≥2 effective children** — a genuine branch point.
- A mutation node with **exactly one effective child** is **collapsed**: its
  mutation is folded onto that single child's branch.
  - If the single child is a **child-lineage**, recurse into it, appending this
    node's mutation to the branch (mutations accumulate ancestral→derived).
  - If the single child is a **single haplotype tip**, the mutation rides on
    that tip's own branch: `H0002_n56:1[&&NHX:mutation=m.309C>CCC:reads=56]`.
- A mutation node with **zero effective children** (a dead-end mutation nobody's
  best-fit carries and that has no descendants) is emitted as an ancestral leaf
  `anc<id>:<len>[&&NHX:mutation=...]` so the mutation is not lost.
- `ROOT` (germline) is **never collapsed** — it is always the tree anchor, even
  if it has a single child. When `ROOT` has one child, the tree is simply
  `(<child subtree>)ROOT;`.
- Consequently a branch may carry **several** mutations, comma-joined in order
  from ancestral to derived:
  `[&&NHX:mutation=m.3144A>G,m.310T>TC]`.

## Output format details

- **Tip label:** `H<id>_n<read_count>` (e.g. `H0002_n56`). Alphanumeric +
  underscore, no quoting needed — the read count is human-readable in the label.
- **Internal label:** `anc<node_id>` for ancestral haplotypes, `ROOT` for
  germline.
- **Read support (machine-readable):** every haplotype tip carries its read
  count as an NHX field on its branch, `[&&NHX:reads=<count>]`, so viewers can
  size/color tips by support without parsing the label. When the tip's branch
  also carries mutations, both fields appear:
  `[&&NHX:mutation=m.309C>CCC:reads=56]`. Internal (ancestral) nodes and `ROOT`
  have no `reads` field (they are unobserved).
- **Branch mutations:** NHX `[&&NHX:mutation=<comma-list>]`. Mutation names
  contain `.` and `>` but never `:` or `=` (the NHX field delimiters), so they
  are safe unquoted inside NHX.
- **Branch length:** the number of mutations on that branch (`1` normally, more
  when collapsed). A haplotype tip that acquires no new mutation relative to its
  parent node gets length `0` and no `mutation` field (but still carries its
  `reads` field).
- **Root:** has no incoming branch and no mutation annotation.
- Standard Newick termination with `;` and a trailing newline.

### Edge cases

- A surviving mutation node with attached haplotypes but no child lineages emits
  those haplotypes as its tip children.
- A mutation node that is a genuine tip in the SCITE tree (no child mutations)
  *and* has no attached haplotypes is an unobserved ancestral haplotype; after
  collapse it can only persist if it is a branch point, so in practice such
  dead-end nodes are collapsed away. If one is the sole content of the tree it is
  emitted as a labeled leaf `anc<id>` so the mutation is not lost.
- `run_scite_pipeline` already guarantees ≥2 variants, so the tree always has a
  usable structure.

## Interface

```rust
/// Write the read lineage tree (Newick + NHX branch mutations, tsinfer-style)
/// to `path`. Tips are haplotypes; branches carry the mutations acquired,
/// with pass-through mutation nodes collapsed onto a single branch.
pub fn write_read_lineage_newick(
    tree: &MutationTree,
    hap_matrix: &HaplotypeMatrix,
    rates: &ErrorRates,
    path: &str,
) -> Result<()>;
```

Called in `run_scite_pipeline` after the existing four writers, writing
`<output_prefix>.read_lineage.nwk`. `run_scite_pipeline` already receives
`hap_matrix` and constructs `rates`, so no signature change is required there or
in `lineage::start` / the CLI.

## Testing (TDD)

Unit test on a small hand-built `MutationTree` + `HaplotypeMatrix` asserting the
exact Newick string, covering:
1. A branch point with ≥2 haplotype tips, each carrying `[&&NHX:...:reads=<n>]`.
2. A collapsed pass-through mutation (branch carrying a comma-joined mutation
   list), with `mutation` and `reads` fields co-present on a tip.
3. A haplotype attached to an internal (read-bearing) node → zero-length tip
   edge with no `mutation` field but still a `reads` field.

Plus the existing `run_scite_pipeline` file-existence test extended to assert
`.read_lineage.nwk` is written.

## Out of scope

- Real `tskit` `.trees` output (would add the `tskit` crate + node-time/site-
  position modeling). Explicitly deferred.
- Any change to error-rate handling, the MCMC, or the other four outputs.
