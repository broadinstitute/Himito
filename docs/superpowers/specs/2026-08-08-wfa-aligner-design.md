# Replace `gap_open_aligner` with a WFA-based aligner

Date: 2026-08-08
Status: approved (pending implementation)

## Problem

`gap_open_aligner` in `src/build.rs` aligns every graph edge's alt sequence
against its reference sequence using rust-bio's `Aligner::global`, a
Needleman-Wunsch/Gotoh gap-affine aligner that runs in O(n*m) time regardless of
how similar the two sequences are.

`generate_cigar` calls it once per edge, under rayon, on sequences up to
`--length-max` (default 3000 bp). Mitochondrial graph edges are near-identical to
the reference, so the quadratic aligner spends nearly all its time on a band it
never needs. Measured on 3000 bp pairs at 0.1% divergence, a batch of 100
alignments takes 22.2 s.

## Goal

Make edge alignment dramatically faster without changing which alignments are
considered optimal, and without regressing the homopolymer-artifact handling that
`normalize_homopolymer_indels` provides.

## Non-goal

Improving accuracy by retuning the scoring model. See "Why scoring cannot fix the
homopolymer artifact" below - this was investigated and rejected on evidence.

## Approach

Replace the rust-bio aligner with the wavefront alignment algorithm (WFA) via the
`rust_wfa` crate, keeping rust-bio as a bounded-memory fallback.

WFA is exact: it computes the same optimal gap-affine alignment as
Needleman-Wunsch, in O(n*s) time where `s` is the alignment score. For similar
sequences `s` is small and the win is large.

### Scoring equivalence

The current rust-bio scoring is match `M=+1`, mismatch `S=-1`, gap open `O=-5`,
gap extend `E=-1`, with a gap of length `k` costing `O + k*E`.

WFA requires the match penalty to be 0. The score-maximization problem converts to
penalty-minimization by:

```
x (mismatch penalty)     = 2*(M - S) = 4
o (gap open penalty)     = -2*O      = 10
e (gap extend penalty)   = M - 2*E   = 3
```

with the relation `penalty = M*(n + m) - 2*score` for a global alignment of
sequences of length `n` and `m`.

This was validated empirically: across 800 randomized pairs spanning lengths
50-3000 and divergences 0.1%-15%, the WFA penalty equalled the transformed
rust-bio score in 800/800 cases. **These penalties are therefore exactly
equivalent to today's scoring**, and are what the implementation will use.

### Measured speedup

100 alignment pairs per cell, release build, Apple Silicon:

| length | divergence | rust-bio | rust_wfa | speedup |
|-------:|-----------:|---------:|---------:|--------:|
| 500  | 0.1% |   480 ms |   1.9 ms |  259x |
| 500  | 1%   |   472 ms |   6.8 ms |   70x |
| 500  | 5%   |   510 ms |  93.6 ms |  5.4x |
| 500  | 15%  |   537 ms |   467 ms |  1.1x |
| 1500 | 0.1% |  4604 ms |   4.6 ms |  990x |
| 1500 | 1%   |  4755 ms |  50.9 ms |   93x |
| 1500 | 5%   |  4996 ms |   669 ms |  7.5x |
| 1500 | 15%  |  5035 ms |  4070 ms |  1.2x |
| 3000 | 0.1% | 22207 ms |   9.0 ms | 2471x |
| 3000 | 1%   | 22414 ms |   166 ms |  135x |
| 3000 | 5%   | 22752 ms |  2708 ms |  8.4x |
| 3000 | 15%  | 23516 ms | 18312 ms |  1.3x |

The benefit is concentrated below ~5% divergence, which is where mitochondrial
graph edges live. Above that, WFA stops paying for itself - which motivates the
fallback guard.

## The accuracy trap this design must avoid

WFA is exact but its traceback resolves ties differently from rust-bio. On the
chrM:310 poly-C case that `normalize_homopolymer_indels` exists to handle:

```
rust-bio: 23=4I1X20=      normalizer fires,  X removed
rust_wfa: 23=1X2=4I18=    normalizer misses, X survives
```

Both CIGARs have the same optimal penalty. But `normalize_homopolymer_indels`
only matches a length-1 `X` **directly adjacent** to an `I`/`D` run. WFA puts a
`2=` between them, so a drop-in swap would silently reintroduce the false
near-homoplasmic T>C SNP at chrM:310 that the normalizer was written to suppress.

### Fix: an indel left-alignment pre-pass

Insert a `left_align_indels` step between alignment and normalization. It
repeatedly slides each `I`/`D` run leftward through a preceding `=` run, one base
at a time, while the base leaving the right end of the gap equals the base
entering at the left end. This is the standard `bcftools norm` left-shift,
restricted to sliding across `=` only - never across an `X` or another indel - so
it is content-preserving by construction.

Verified on all three cases the normalizer's existing tests cover:

```
chrM:310 ins: 23=1X2=4I18= -> 23=1X4I20= -> 23=3I1=1I20=   (X removed)
chrM:310 del: 23=1X2=4D18= -> 23=1X4D20= -> 23=3D1=1D20=   (X removed)
genuine SNP:  8=2I1X       -> 7=2I1=1X   -> 7=2I1=1X       (X preserved)
```

All three round-trip: applying the CIGAR to the reference reproduces the alt
sequence exactly.

## Why scoring cannot fix the homopolymer artifact

Retuning penalties was considered and rejected. For the chrM:310 case, the
biologically correct edit is two separate insertions (3 C's before the lone T, 1 C
after), which pays two gap-open penalties:

```
correct 2-gap edit: (10 + 3*3) + (10 + 1*3) = 32
artifact 1-gap+X:   (10 + 4*3) + 4          = 26
```

The artifact is cheaper. Under any affine gap model with a meaningful gap-open
penalty, one gap plus one mismatch beats two gaps, so no penalty retuning makes
the correct edit score-optimal. `normalize_homopolymer_indels` performs a
correction the scoring model structurally cannot, and is therefore load-bearing
and must be preserved.

## Design

### New module: `src/wfa_align.rs`

`build.rs` is already 1136 lines and mixes graph construction with alignment.
All alignment and CIGAR post-processing moves into a focused module with a small
surface:

```rust
/// Align `alt` against `reference`, returning a CIGAR in Himito convention:
/// `I` consumes alt, `D` consumes ref, `=`/`X` consume both.
pub fn align_cigar(reference: &str, alt: &str) -> String
```

Internals, all private to the module:

- `WFA_PENALTIES` - the `x=4, o=10, e=3` constant.
- `wfa_cigar(reference, alt) -> Option<String>` - runs `wavefront_align` and
  converts its gapped `query_aligned`/`text_aligned` rows into a run-length CIGAR.
  Returns `None` on any `AlignmentError`.
- `bio_cigar(reference, alt) -> String` - the current rust-bio path, retained
  verbatim as the fallback.
- `left_align_indels(cigar, reference, alt) -> String` - new, described above.
- `normalize_homopolymer_indels`, `parse_cigar`, `alignment_to_cigar` - moved from
  `build.rs` unchanged, along with their three existing tests.

`align_cigar` composes them:

```
prefilter -> wfa_cigar (or bio_cigar fallback) -> left_align_indels
          -> normalize_homopolymer_indels
```

### Handling `rust_wfa`'s constraints

`wavefront_align(query, text, pens)` rejects two inputs:

- `ZeroLength` if either string is empty.
- `QueryTooLong` if `query.len() > text.len()`.

`wfa_cigar` handles the length constraint by passing the shorter sequence as
`query` and, when that means swapping ref and alt, mirroring `I` <-> `D` when
building the CIGAR. Empty inputs return `None` and route to the fallback.

### Fallback guard

WFA uses O(s^2) **memory**. At 3000 bp and 15% divergence a single alignment
holds roughly 80 MB of wavefront grid; `generate_cigar` runs these under rayon
across every core, so a cluster of divergent edges could spike memory into the
gigabytes.

Before aligning, `align_cigar` prefilters with
`bio::alignment::distance::simd::bounded_levenshtein(ref, alt, k)`, a SIMD-backed
bounded edit distance that returns `None` once the distance exceeds `k`.

- `k = max(1, 5% of max(ref.len(), alt.len()))`
- distance within `k` -> WFA path
- distance exceeds `k`, or WFA returns an error -> rust-bio fallback

At 5% divergence WFA is still 8.4x faster than rust-bio and needs only a few MB,
so the threshold sits right where the tradeoff turns. Above it, behaviour is
identical to today's, since rust-bio is exactly what runs now. The prefilter is
cheap relative to either aligner.

Note the two thresholds measure different things: `k` bounds edit distance
(indels and mismatches), while the WFA penalty `s` is a weighted score. The 5%
edit-distance cap is a conservative proxy that bounds `s` well below the levels
where memory becomes a concern.

### Call site

`gap_open_aligner` is removed from `build.rs`. Its single caller in
`generate_cigar` (line ~870) calls `wfa_align::align_cigar(&ref_seq, &alt_seq)`
instead. No signature change, no change to `generate_cigar`'s structure, no CLI
change.

### Dependency

Add to `Cargo.toml`, renamed because the crate's lib target is confusingly named
`lib`:

```toml
wfa = { package = "rust_wfa", version = "1.0.0" }
```

Cost, flagged for the record: `rust_wfa` was last published in 2022 and pulls in
clap 3.2, strum 0.24, num_cpus and rand 0.8, duplicating Himito's clap 4.5 and
rand 0.9. This adds roughly a minute to a cold build and nothing at runtime.
Cargo resolves the duplicate majors without conflict. It compiles clean on the
current toolchain.

## Testing

Existing, preserved unchanged (moved into `wfa_align.rs`):

- `homopolymer_insertion_absorbs_adjacent_mismatch`
- `genuine_snp_next_to_unrelated_indel_is_left_alone`
- `homopolymer_deletion_absorbs_adjacent_mismatch`

New:

1. **Left-align unit tests** - the three verified cases above, asserting both the
   resulting CIGAR and that a genuine SNP keeps its `X`.
2. **Differential optimality test** - randomized sequences across lengths
   50-3000 and divergences 0.1%-15%; assert the WFA penalty equals
   `M*(n+m) - 2*bio_score`. This is the 800/800 check from the design probe,
   promoted into the suite with a fixed RNG seed for reproducibility.
3. **Round-trip property test** - applying the final CIGAR to the reference
   reproduces the alt sequence exactly, and ref-consuming ops sum to `ref.len()`
   while alt-consuming ops sum to `alt.len()`.
4. **Swap path** - a case where `alt.len() > ref.len()`, confirming `I`/`D` are
   mirrored correctly rather than inverted.
5. **Fallback path** - a pair beyond the 5% threshold produces a valid CIGAR
   identical to the rust-bio result.
6. **Degenerate inputs** - empty alt, and `alt == ref`.

## Risks

| Risk | Mitigation |
|------|-----------|
| WFA tie-break reintroduces the chrM:310 artifact | `left_align_indels` pre-pass, verified on both mirrors; covered by tests 1 and the existing three |
| O(s^2) memory blowup under rayon | 5% bounded-levenshtein prefilter routes divergent pairs to rust-bio |
| `rust_wfa` unmaintained since 2022 | Exactness pinned by the differential test; rust-bio retained in-tree as a working fallback, so reverting is a one-line change |
| CIGARs shift for edges that were already correct | Left-alignment is a canonicalization, so downstream variant positions become more consistent, not less; round-trip test guarantees content is preserved |
