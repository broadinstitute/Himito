# Next steps: relaxing the infinite-sites assumption

Status: **not implemented.** This records the design so the decision does not
have to be re-derived.

## Where we are

Himito's lineage inference assumes infinite sites (ISA): every variant arises
exactly once and never reverts. This is not enforced by a check — it is an
invariant of the data structure. `scite::MutationTree` has one node per variant
and defines a genotype as the set of mutations on the root-to-node path, so no
search move can produce a second copy of a mutation or a back-mutation.

(The same property is what makes `scite::AttachmentScorer` possible: because
`genotype(child) = genotype(parent) + one mutation`, a read's score at a child
is its score at the parent plus a single precomputed term.)

Two consequences follow, and both are currently unaddressed:

1. **Every output tree is a perfect phylogeny regardless of the input.** Genuine
   recurrent mutation, back-mutation, or NUMT-derived recombination is absorbed
   as false-positive/false-negative noise and imputed away. Nothing reports that
   this happened.

   The size of the effect *is* recoverable today with no new code: diff
   `<prefix>.cleaned_matrix.csv` against the input `matrix.csv`. Every cell that
   differs is an observation the model had to call wrong.

2. **`<prefix>.variant_cooccurrence.tsv` cannot show a violation.** It is
   computed on the SCITE-cleaned matrix, where any two mutations are either
   nested or in disjoint subtrees. Nested implies `n01 = 0`; disjoint implies
   `n11 = 0`. One cell of the 2x2 table is always structurally zero, so `D'` is
   always exactly +/-1 and Fisher's p is a function of the topology, not the
   data. **Use `<prefix>.raw_variant_cooccurrence.tsv` for anything
   inferential.**

## Option A: quantitative ISA screen on the raw matrix

The cheapest real improvement, and a prerequisite for B and C.

`lineage::four_gamete_test` asks a binary question — were all four gametes ever
seen — over the deduplicated haplotype matrix. At `fp = 0.005` and `fn = 0.05`
over hundreds of reads, a site pair with zero real violations will routinely show
all four gametes by chance, so the count it logs is close to uninformative. Its
result is also discarded: only `.len()` reaches the log line at
`src/lineage.rs:748`; the variant-pair names are dropped.

Replace the binary test with a counted one. The 2x2 tables are already computed
in `scite::write_raw_variant_cooccurrence`. Under a fitted tree where `i` is
ancestral to `j`, the minority gamete `(i=0, j=1)` requires a dropout at `i`, so
its count is approximately `Binomial(n11 + n01, fn)`. Test the observed minority
count against that tail.

```rust
/// Pairs whose minority gamete count exceeds what fp/fn can explain.
pub fn isa_violations(
    matrix: &BinaryMatrix,
    tree: &MutationTree,
    rates: &ErrorRates,
    alpha: f64,
) -> Vec<(usize, usize, usize, f64)>; // (i, j, minority_count, p)
```

Run it on the **raw** `BinaryMatrix`, never on `CleanedMatrix` — the cleaned one
cannot violate, for the reason in point 2 above. Write
`<prefix>.isa_violations.tsv`. Apply Benjamini-Hochberg across pairs; `adjustp`
is already a dependency and `scite::write_pair_cooccurrence_table` already uses
it.

This is reporting only. It changes no tree and breaks no existing output.

## Option B: drop the offending variant

Given Option A's list, drop the lower-support variant of each failing pair before
the MCMC. One flag, a handful of lines, and it makes the ISA claim honest: the
tree is then built only on variants for which ISA is defensible.

The cost is information. A recurrent variant at a real mtDNA hotspot is
biologically interesting, and this throws it away.

## Option C: finite sites — allow a variant two nodes

The principled fix. Let a flagged variant occupy two nodes in the tree, so it can
arise independently on two lineages. This is the Dollo / persistent-phylogeny
relaxation, and it is what the SCITE successors do (SiFit, SCARLET).

Sketch:

- `MutationTree.parent` becomes indexed by *node*, with a separate
  `node_variant: Vec<usize>` mapping node to variant. Most nodes map one-to-one;
  a flagged variant gets two nodes sharing a `node_variant` value.
- `ancestor_mask` sets a bit per *variant*, so a genotype carrying either copy
  reads as carrying the variant. `attachment_log_likelihood` needs no change.
- `AttachmentScorer` still applies unchanged: the incremental identity
  `genotype(child) = genotype(parent) + one variant` holds as long as the two
  copies are never on one root-to-leaf path.
- That last condition is the new hard constraint, and it is structurally the same
  one `scite::violates_position_exclusivity` already enforces for overlapping REF
  spans. Reuse that machinery rather than writing a second copy of it.
- Model selection matters more than the search: doubling a variant always fits at
  least as well, so it needs a penalty. BIC over the number of duplicated
  variants is the obvious choice, with Option A's p-values deciding which
  variants are even eligible.

## What was considered and dropped

A Gusfield O(H*M) perfect-phylogeny construction with minimum-flip repair
(Chen, Eulenstein, Fernandez-Baca & Sanderson, *Minimum-flip supertrees:
complexity and algorithms*, IEEE/ACM TCBB 3(2):165-173, 2006) was specced as an
error-model-free alternative output and then dropped.

The reason: minimising flips is exactly maximising the existing SCITE likelihood
with `fp == fn`. With equal rates, a read's log-likelihood at a node is
`(#agreements)*ln(1-p) + (#disagreements)*ln(p)`, which is
`const + (#disagreements)*ln(p/(1-p))` — and `ln(p/(1-p)) < 0` for `p < 0.5`, so
the maximiser is the minimiser of disagreements. Running `Himito lineage` with
`--fp-rate` and `--fn-rate` set equal already performs a minimum-flip search.

What that substitution does *not* give you is Gusfield's exact answer in O(H*M)
when the data happens to be conflict-free, as opposed to MCMC's heuristic one.
If that exactness ever matters, this is the note to reopen.
