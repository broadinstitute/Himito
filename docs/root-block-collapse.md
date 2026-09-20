# Collapsing homoplasmic variants into a root block

Status: **implemented and validated on one real sample; NOT ready for
cross-depth or cross-technology comparison.** On by default;
`--no-root-block` disables. Shipped in commits `ad22eb6..e174098`.

It recovers the rCRS haplogroup backbone blindly on real PacBio data (see
"Real-data validation" below), but its decision rule is depth-dependent, which
is disqualifying for the comparison this feature exists to support. Read both
caveats before using it for anything published.

Caveats before trusting it:

* **The benchmark re-run promised below has NOT been done.** `sweep_fpfn.sh`
  needs `sim.matrix.csv` from a prior `run_eval.sh` run, which is not present,
  and regenerating it needs the simulation environment. `sweep_metrics.tsv` in
  the tree is still the pre-change baseline.
* **Block membership is depth-dependent — a known limitation, not yet
  fixed.** The decision rule is a q-value threshold (see "The decision rule"
  below), and a q-value tests *statistical power*, not effect size: with more
  reads, the same dropout pattern reaches significance and the variant is
  pulled out of the block. Measured on `Himito_benchmark/HG002.matrix.csv` by
  subsampling the same sample to different read counts: 125/250/500 reads all
  block 10 variants, 1000 reads blocks 8, 2000 reads blocks only 3. The
  planned fix is to replace the criterion with an **equivalence test** —
  blocking on positive evidence that the association is negligible, rather
  than on a failure to detect one, which would make membership converge as
  depth increases instead of shrinking. That work is not yet done. Until it
  lands, **do not use this feature to compare block composition across
  different depths or technologies** — a shrinking block at higher depth is
  an artifact of this test, not necessarily a real difference in the
  underlying variants.

Verified: the escape hatch is exact. All eight output files are byte-identical
between a binary built at the pre-feature commit `4d59bd6` and this one run
with `--no-root-block`, on real ONT data.

**Real-data validation.** `SMaHT/SMHT005-3AK-pacbio.matrix.csv`, run at stock
default thresholds, blocks exactly the nine variants that make up the rCRS
haplogroup backbone: `m.263A>G, m.750A>G, m.1438A>G, m.3010G>A, m.4769A>G,
m.8860A>G, m.9743A>G, m.15326A>G, m.16519T>C` (HF range 0.836–0.899). The one
candidate the test rejects and keeps in the search, `m.16189T>C`, is refused
at q = 1e-203 because its absences are marked by the adjacent poly-C indel
`m.16188CT>C` — a read-level dropout confounder the test correctly separates
from the real backbone. This is the feature's best validation to date: a
haplogroup-defining variant set recovered blindly, with the one genuinely
structured absence pattern in the matrix correctly excluded.

## Where we are

`Himito lineage` builds a SCITE mutation tree over every variant that clears the
`--min-hf` / `--max-hf` window. On real data the top of that tree is a long unary
run of germline, haplogroup-defining variants, and **its order is not
reproducible between sequencing technologies**.

HG00290, same individual, same six variants, both trees pure chains:

| depth | ONT (`readnum`) | PacBio (`readnum`) |
|-------|-----------------|--------------------|
| 1 | `750A>G` (240) | `13827A>G` (1) |
| 2 | `16526G>A` (213) | `16526G>A` (5) |
| 3 | `13827A>G` (774) | `750A>G` (16) |
| 4 | `9477G>A` (2148) | `9477G>A` (61) |
| 5 | `310T>TC` (1071) | `310T>TC` (138) |
| 6 | `302A>AC` (3377) | `302A>AC` (582) |

The top three are exactly reversed. The bottom three agree.

### Why the order is noise

On a unary path `root → A → B`, a read observing `A=1, B=0` costs nothing under
`A → B` (it attaches at `A`, exact match) but costs `min(ln fp, ln fn)` under
`B → A`. So the variant observed *alone* more often is placed higher, and the
variant *missing* more often is pushed to the bottom. The inferred trunk order
is a sort of the sites by explicit-REF-call rate.

Missingness here is homopolymer dropout (302, 310, 16181–16183 are poly-C
tracts), coverage, and reads not spanning a site on a 16.5 kb circle. All three
are technology-specific. Hence the flip.

### Why `--max-hf` makes it worse, not better

Matrix-derived HF is `present / (present + absent)` over covered reads only
(`lineage.rs:510-515`). A near-homoplasmic variant therefore clears `--max-hf`
**iff enough covering reads carry an erroneous REF call.** Among germline
variants the filter admits the dirty ones and excludes the clean ones, then
hands exactly those erroneous REF calls to the MCMC as the ordering signal.

The CLI default is `--max-hf 0.95` (`main.rs:566`); the runs that produced the
trees above used `0.99`, which widens the admission window further.

### Why no HF threshold can fix it

Germline `9477G>A`, one individual:

* ONT: HF ≈ 0.84
* PacBio: HF ≈ 0.97

No single cutoff separates that from a genuine 96% heteroplasmy. **HF magnitude
is the wrong discriminator.**

## Design

Partition variants before the search. Homoplasmic ones become an unordered block
annotated on the root and are removed from the matrix the MCMC sees; the tree is
searched over the informative variants only.

```
parse_binary_matrix                -> BinaryMatrix                  (unchanged)
  rootblock::partition             -> block + informative           [NEW]
  restrict to informative          -> a real sub-BinaryMatrix       [NEW]
  run_mcmc_multichain / polish /
  exclusivity / attach_all_reads   -> over the sub-matrix           (unchanged)
  re-expand block into outputs                                      [NEW]
```

### New module: `src/rootblock.rs`

`src/scite.rs` is already ~150 KB and carries matrix scoring, MCMC, polish,
imputation and six output writers. This is a self-contained statistical test
with a clean interface, so it gets its own unit.

```rust
pub struct RootBlockConfig {
    pub min_hf: f64,      // default 0.80 - candidate gate only, NOT the decision
    pub max_q: f64,       // default 0.05
    pub min_absent: usize // default 10
}

pub struct RootBlock {
    pub block: Vec<usize>,        // variant indices, ascending
    pub informative: Vec<usize>,  // variant indices, ascending
    pub audit: Vec<RootBlockAudit>,
}

/// One row per candidate, i.e. per variant clearing `min_hf`. Variants below
/// the gate are never candidates and produce no row.
pub struct RootBlockAudit {
    pub variant: usize,
    pub hf: f64,
    pub n_absent: usize,
    /// `None` when admitted by rule 2 (too few absences to test).
    pub min_q: Option<f64>,
    /// The `u` achieving `min_q`; `None` under rule 2.
    pub partner: Option<usize>,
    pub in_block: bool,
    /// Why: rule 2/no-testable-partner, tested-and-blocked, tested-and-kept,
    /// or (set later by `resolve_span_conflicts`) span-conflict pull-back.
    /// See "Outputs" below for the four values. Overwritten, not merged, by
    /// the span-conflict pass.
    pub reason: BlockReason,
}

pub fn partition(matrix: &BinaryMatrix, cfg: &RootBlockConfig) -> RootBlock;
```

It depends only on `BinaryMatrix` and the existing pair statistics. Testable in
isolation: hand it a synthetic matrix, assert the partition.

`scite::fisher_greater` (`scite.rs:955`), `scite::pair_cooccurrence`
(`scite.rs:998`) and the Benjamini-Hochberg helper used at `scite.rs:1011` are
currently private. They become `pub(crate)` so `rootblock` can reuse them rather
than reimplementing a second copy of Fisher's exact test.

**It runs on the raw `BinaryMatrix`, never the cleaned one.** Per
`docs/finite-sites-relaxation.md`, co-occurrence on the SCITE-cleaned matrix is
structurally degenerate (one cell of every 2x2 is forced to zero), so any
association test there measures topology, not data.

### The decision rule

HF is recomputed inside `partition` using the same convention as
`lineage.rs:515` (`present / (present + absent)`), so the two stages agree.

For a variant `v`:

1. `hf(v) < cfg.min_hf` -> **informative.** Not a candidate; never tested.
2. `hf(v) >= cfg.min_hf` and `n_absent(v) < cfg.min_absent` -> **block.**
   Too few absences to test structure, and too few to inform topology either
   way. This is the PacBio case (`readnum` 0-5 on the trunk).
3. Otherwise, run the absence-structure test below.
4. **Span-conflict pull-back** (`resolve_span_conflicts`, `scite.rs`), applied
   after 1-3 have produced a candidate block: any block member whose REF span
   overlaps *any* other variant in the full matrix — blocked or informative —
   is pulled back into the search, unconditionally. `append_block_to_cleaned`
   imputes every block member present on every read; wherever the overlapping
   variant is also called present, that would assert two mutually-exclusive
   alleles on one molecule, which was structurally impossible before the
   partition ran. A real heteroplasmy at that locus means the high-HF member
   is a major allele, not fixed germline, and belongs in the tree next to its
   partner. This rule overrides whatever rules 1-3 decided; it never adds a
   variant to the block, only removes one.

#### Absence-structure test

For candidate `v` and each other variant `u`, over reads jointly covered at
both:

|             | `u` alt | `u` ref |
|-------------|---------|---------|
| `v` absent  | a       | b       |
| `v` present | c       | d       |

One-sided Fisher for enrichment of cell **a** (`scite::fisher_greater`,
`scite.rs:955`), Benjamini-Hochberg across all tested pairs (the same correction
already applied in `scite.rs:1011`). Then:

> `v` joins the root block iff `min_u q_vu > cfg.max_q`.

The question it asks is: **are the molecules that lack `v` distinguished by
carrying some other allele?** A real subclone sitting at `1 - HF` answers yes for
at least one `u`. Pure dropout answers no for every `u`. This is what separates
ONT's `9477G>A` at HF 0.84 from a genuine 84% heteroplasmy, which a threshold
cannot.

#### Why the test is against `u`'s ALT calls, not `u`'s absences

**This is the load-bearing choice in the design. Do not "simplify" it.**

Dropout is partly read-level: a truncated or poor-quality read drops out
everywhere. So the absences of two germline variants are positively correlated,
and an absence-vs-absence test would read read quality as a subclone and keep
both variants in the tree. Testing against `u`'s alt calls demands positive
distinguishing evidence, which dropout does not produce.

#### Known false negative

A real subclone defined *only* by the loss of `v`, carrying no private allele
anywhere, is invisible to this test and will be folded into the block. In
mitochondrial data this is unlikely but not impossible. Mitigations: the
`--no-root-block` escape, and the per-variant audit table below, so the decision
is auditable rather than silent.

### Index remapping

The restriction produces a genuine `BinaryMatrix` carrying its own `variants`
name vector — **not** an index mapping threaded through the pipeline. Everything
downstream (`AttachmentScorer`, `run_mcmc_multichain`, `polish_unary_path_order`,
`enforce_position_exclusivity`, `attach_all_reads`, `edge_evidence`, all writers)
addresses variants by position within whatever matrix it was handed, so it works
unchanged. Only the re-expansion step ever needs to know a block exists.

### Imputation

Block variants are homoplasmic by construction, so in `CleanedMatrix` every read
gets `1` for every block variant. This is the point, not a shortcut: the block
designation corrects exactly the dropout that flagged the variant.
`cleaned.attachment` is untouched — attachments are over the informative tree.

## Outputs

**New `<prefix>.root_block.tsv`** — one row per *candidate*, i.e. per variant
clearing `--root-block-min-hf`. Variants below the gate produce no row:

```
variant  hf  n_absent  min_q  partner_variant  in_block  reason
```

`partner_variant` is the `u` achieving `min_q`. Candidates admitted by rule 2
(untestable) report `min_q = NA`, `partner_variant = NA`, `in_block = true`.
`partner_variant` is also reported as `NA` — even though `partner` is
`Some(u)` internally — whenever `min_q >= 1.0`: a saturated tie means no
partner was ever strictly better than the incumbent, so whichever `u` happens
to sit lowest in variant order wins by construction, not by evidence, and
printing it would read as "the best-discriminating partner" when it is
nothing of the kind.

`reason` disambiguates rows that would otherwise look identical or
contradictory. It takes one of four values:

* `untested` — admitted without testing: either rule 2 (too few
  jointly-covered absences), or no candidate partner offered a testable alt
  call at all. Both report `min_q = NA`, `partner_variant = NA`.
* `unstructured_absences` — tested; every partner's q exceeded `max_q`. No
  evidence of structure, so blocked.
* `structured_absences` — tested; some partner's q was at or below `max_q`.
  The absences are structured, so kept in the search (`in_block = false`).
* `span_conflict` — rule 3 (or rule 2) would have blocked this variant, but
  rule 4 pulled it back because its REF span overlaps another variant. Without
  this column, a span-conflict pull-back prints `min_q = NA`,
  `partner_variant = NA`, `in_block = false` — the same `NA`/`NA` shape rule 2
  produces, but with `in_block` flipped, and the table gives no reason why.
  `reason` is what tells the two apart.

**`<prefix>.mutation_tree.tsv`** — one extra row: reserved node id, `variant` =
comma-joined block members, `parent_variant` = `ROOT`, `n_reads_attached` =
reads attaching at root, support columns `NA`. Existing columns are unchanged so
downstream parsers and the plotting script keep working.

**`<prefix>.variant_cooccurrence.tsv`** — pairs involving a block variant are
**omitted entirely**, not printed with `NA` fields. `write_variant_cooccurrence_upto`
restricts the table to the first `n_informative` columns before building any
row, so a block variant's row simply does not exist in this file (it is still
present, with real statistics, in `<prefix>.raw_variant_cooccurrence.tsv`).
The reason is the same one that motivates the restriction in the first place:
the variant is invariant in the cleaned matrix, so `dprime_r2`'s zero-variance
guard (`scite.rs:991`) would return `0.0` for any pair involving it, which
does not crash but is a false statement about the data.

**`<prefix>.molecule_summary.tsv`** — every block variant is listed in every
read's `variants` column and counted in `n_variants`, **even for reads that
did not span the locus.** This follows directly from imputation: block
variants are homoplasmic by construction, so `CleanedMatrix` carries them as
an all-`1` column for every read regardless of what that read actually
covered, and `write_molecule_summary` reads straight off `CleanedMatrix`. This
is consistent with `<prefix>.cleaned_matrix.csv`, which shows the same all-`1`
columns, and is intended, not a bug — but it does mean per-read variant counts
in this file are not raw observations for block members.

**Newick** — the root already emits an explicit tip with no edge mutations
(`scite.rs:1329`). Annotate that token with `[&&NHX:root_block=v1,v2,...]`,
reusing the existing NHX convention. No structural change to
`emit_lineage_node`.

## CLI

| flag | default | meaning |
|------|---------|---------|
| `--no-root-block` | off (collapsing is **on**) | disable entirely; reproduces pre-change behaviour |
| `--root-block-min-hf` | `0.80` | candidate gate |
| `--root-block-max-q` | `0.05` | BH q above which absences count as unstructured |
| `--root-block-min-absent` | `10` | below this, admit without testing |

## Edge cases that must not panic

* **Empty block** -> output byte-identical to `--no-root-block`. This is a test,
  not a hope.
* **All variants homoplasmic** -> informative set empty. The guard at
  `scite.rs:1701` currently bails on zero variants telling the user to relax
  `--min-hf`/`--max-hf`. It needs a separate branch: the correct result here is a
  valid `ROOT + block` tree, not an error.
* **One informative variant** -> the existing `single_variant` path
  (`scite.rs:1712`) must still fire.

## Tests

1. Synthetic germline (absences scattered) plus a true 90% heteroplasmy
   (absences clustered with a private alt) -> germline in block, heteroplasmy
   not.
2. **Read-level dropout confounder**: two germline variants whose absences
   co-occur on the same poor-quality reads -> both still land in the block. This
   is the test that proves an absence-vs-absence test would have failed. Name it
   so it survives future refactoring.
3. Empty block -> output identical to disabled.
4. All-homoplasmic -> valid tree, no panic.
5. Round-trip: block variants are all-1 columns in the cleaned matrix.
6. Existing pipeline tests re-run under `--no-root-block`, proving the escape
   hatch is faithful.

## Consequences of shipping this on by default

The next run changes:

* `lineage_eval/lineage_sim/sweep_metrics.tsv`
* the `HG00097_30x_chrM_hg38*` example outputs
* the Mitorsaw / MitoHiFi / mtdnaserver benchmark comparisons under
  `Himito_benchmark/`

`--no-root-block` reproduces the old numbers. Re-running the benchmark suite
was intended to be part of this work, not a follow-up — but per the header,
**it has not actually been done**: `sweep_fpfn.sh` needs a `sim.matrix.csv`
this checkout does not have, and `sweep_metrics.tsv` in the tree is still the
pre-change baseline. Treat the benchmark comparison as outstanding.

## Rejected alternatives

* **Tighten `--max-hf` only.** Fewest lines, but as the `9477G>A` case shows the
  threshold has to be set per technology, and it still both misses germline
  variants and swallows real high-HF heteroplasmies.
* **Drop homoplasmic variants entirely.** Loses the haplogroup calls from every
  output and leaves nothing to show a reviewer about what was excluded.
* **Collapse at report time only.** Zero risk to the likelihood, but the MCMC
  still spends its effort ordering the trunk and reads still pass through it, so
  it fixes the figures without fixing the instability. (Note that
  `emit_lineage_node`, `scite.rs:1273`, already does a purely *topological*
  collapse of unary runs in the Newick output. That is report-only and
  support-blind, and is not a substitute for this.)
* **Post-hoc from the fitted tree** (run MCMC, collapse the low-support
  root-adjacent unary run, re-run). Needs no new variant-level statistic, but
  fails on the actual data: ONT's trunk has *high* order support — hundreds of
  dropout-driven reads — so a support-based rule flags PacBio's trunk and misses
  ONT's, which is the exact failure being fixed.
* **PhyloTree / MITOMAP annotation.** Reliable for mitochondria and defensible
  to reviewers, but adds an external database dependency and a new required
  input, and cannot catch private or novel germline variants. Worth revisiting
  as a cross-check, not as the primary rule.

## Related, out of scope

`--data-type` (`main.rs:571`) offers `pacbio`, `ont-r9`, `ont-r10`,
`ont-denoised`. There is no `pacbio-denoised` preset, so an ONT-denoised /
PacBio-raw comparison is confounded: denoising changes the explicit-REF-call
profile, which is what both this partition and the trunk order depend on.
