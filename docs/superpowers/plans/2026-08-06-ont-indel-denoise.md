# ONT Small-Indel Denoiser Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `Himito denoise` to correct small (< 5 bp) indel errors in ONT reads by two-directional MAP reassignment at left-normalized indel sites, behind an opt-in `--indels` flag.

**Architecture:** `src/denoise.rs` becomes a module directory. Pure model code (allele types, left-normalization, repeat context, ε/floor, site fit, MAP assignment) lives in `src/denoise/indel.rs`; the length-changing SEQ/QUAL/CIGAR rewrite lives in `src/denoise/rewrite.rs`; `src/denoise.rs` keeps orchestration (pileup collection, deferred site decision, stats, `start`). The correction payload migrates from query coordinates to reference coordinates so substitutions and indels can be applied in one left-to-right walk.

**Tech Stack:** Rust 2021, `rust_htslib` 0.49 (`bam::pileup::Indel`, `bam::record::{Cigar, CigarString}`), `clap` derive, `serde`. No new dependencies.

**Spec:** [`docs/superpowers/specs/2026-08-06-ont-indel-denoise-design.md`](../specs/2026-08-06-ont-indel-denoise-design.md)

## Global Constraints

- **Opt-in.** `IndelOpts::enabled` defaults to `false`. `QuickStart` passes the default. With indels disabled, output must be **byte-identical** to current behavior.
- **ONT only.** `data_type` not starting with `ont` remains a byte-identical passthrough copy.
- **Reference-consumption invariant.** Every rewrite transition preserves the read's total reference consumption (sum of `M`/`=`/`X`/`D`/`N` lengths). `POS`, reference end, and coordinate-sort order are invariant. Assert this in every rewrite test.
- **Exclusive length bound.** `indel_max_len` default `5` means lengths 1–4 are correctable; ≥ 5 is untouchable.
- **Defaults (all tunable, all unvalidated until Task 10):** `max_len 5`, `vaf 0.05`, `err0 0.01`, `err_scale 1.5`, `err_cap 0.4`, `floor_mult 3.0`, `delta 0.3`, `flank 5`, `protect_vaf 0.2`.
- **Substantial alleles are never corrected away.** An allele carried by `>= protect_vaf` of the reads at a site is left alone regardless of candidacy. Added after Task 3's review: `vaf_floor(L)` exceeds 1.0 once L >= 10, so in deep repeat contexts nothing clears candidacy and `assign_allele` would revert every read carrying a real indel to reference. Human mtDNA's poly-C tracts (m.303-315, m.16184-16193) are exactly that shape and are heteroplasmy hotspots. Protection applies symmetrically to REF.
- **Reused SNV parameters:** `min_strand`, `strand_bias_p`, `homoplasmic_vaf` are shared with the existing substitution model — do not duplicate them into `IndelOpts`.
- **Never invent sequence.** Inserted bases are always the site's consensus inserted bases or reference bases, never a per-read guess.
- **Depth is per-column, never a tally of observations.** A read carrying an event reports `Indel::None` at the normalized site's own column *and* its event at the anchor column; counting both double-counts that read and halves every real indel's apparent VAF. The denominator is the pileup depth at `norm_pos`.
- **The flank guard is asymmetric.** A read must extend `flank` reference bases to the left of a site and `flank + max_len` to the right — only a gained deletion consumes extra reference. A symmetric `flank + max_len` on both sides is unsatisfiable near a contig start.
- **Two documented transition restrictions** (deviations from the spec's general MAP statement, added deliberately to bound rewrite complexity):
  1. **No cross-kind reassignment.** A read observing `Ins` may only be assigned `Ref` or another `Ins`; a read observing `Del` only `Ref` or another `Del`. Cross-kind MAP results leave the read untouched, counted in `reads_skipped_unsupported_transition`.
  2. **`Del(n) → Del(m)` only when `m ≤ n`.** Growing a deletion would consume reference from the following CIGAR op. Growing results leave the read untouched, counted in the same stat.
- **Commit after every task.** Run `cargo build --release` before each commit; it must succeed.

---

### Task 1: Indel allele types, options, and left-normalization

Creates the model module and the canonicalization that makes site accumulation possible. Left-normalization is what collapses one biological event anchored at several columns into a single site; without it, support splits and the gate over-corrects.

**Files:**
- Create: `src/denoise/indel.rs`
- Modify: `src/denoise.rs` (add `mod indel;` near the top, after the `use` block ending at line 6)

**Interfaces:**
- Consumes: nothing.
- Produces: `pub enum Allele { Ref, Ins(Vec<u8>), Del(u32) }` with `pub fn net_len(&self) -> i32` and `pub fn kind_matches(&self, other: &Allele) -> bool`; `pub struct IndelOpts` with `Default`; `pub fn normalize_left(refseq: &[u8], anchor: u32, allele: &Allele) -> (u32, Allele)`.

- [ ] **Step 1: Create the module file with types and a stub, plus the failing tests**

Create `src/denoise/indel.rs`:

```rust
//! Small-indel site model for the ONT denoiser.
//!
//! Indels have no per-base quality analogue, so this module supplies the missing
//! per-read likelihood via a repeat-context-scaled error rate `eps(L)`. That one
//! function drives both the candidacy gate and the MAP assignment.

use std::collections::HashMap;

/// An indel allele at a normalized site, indexed conceptually by its net length
/// change: `Ref` = 0, `Ins(s)` = +s.len(), `Del(n)` = -n.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Allele {
    Ref,
    Ins(Vec<u8>),
    Del(u32),
}

impl Allele {
    /// Net length change this allele applies to the read relative to the reference.
    pub fn net_len(&self) -> i32 {
        match self {
            Allele::Ref => 0,
            Allele::Ins(s) => s.len() as i32,
            Allele::Del(n) => -(*n as i32),
        }
    }

    /// True when both alleles are the same variety (both insertions, both deletions,
    /// or both reference). Cross-kind reassignment is not supported by the rewrite
    /// walk, so the assignment step uses this to refuse those transitions.
    pub fn kind_matches(&self, other: &Allele) -> bool {
        matches!(
            (self, other),
            (Allele::Ref, Allele::Ref)
                | (Allele::Ins(_), Allele::Ins(_))
                | (Allele::Del(_), Allele::Del(_))
        )
    }
}

/// Tunable parameters for indel denoising. `Default` is DISABLED, so every existing
/// caller that passes `IndelOpts::default()` keeps byte-identical v1 behavior.
#[derive(Debug, Clone)]
pub struct IndelOpts {
    pub enabled: bool,
    /// Exclusive upper bound: with 5, lengths 1..=4 are correctable.
    pub max_len: u32,
    pub vaf: f64,
    pub err0: f64,
    pub err_scale: f64,
    pub err_cap: f64,
    pub floor_mult: f64,
    pub delta: f64,
    pub flank: usize,
}

impl Default for IndelOpts {
    fn default() -> Self {
        IndelOpts {
            enabled: false,
            max_len: 5,
            vaf: 0.05,
            err0: 0.01,
            err_scale: 1.5,
            err_cap: 0.4,
            floor_mult: 3.0,
            delta: 0.3,
            flank: 5,
        }
    }
}

/// Shift an indel maximally left so the same biological event, however the aligner
/// anchored it, collapses to one canonical `(anchor, allele)` key.
///
/// `anchor` is the reference position of the last base BEFORE the event, which is
/// exactly where htslib's pileup reports it.
pub fn normalize_left(refseq: &[u8], anchor: u32, allele: &Allele) -> (u32, Allele) {
    let _ = (refseq, anchor);
    (anchor, allele.clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn net_len_and_kind() {
        assert_eq!(Allele::Ref.net_len(), 0);
        assert_eq!(Allele::Ins(b"AC".to_vec()).net_len(), 2);
        assert_eq!(Allele::Del(3).net_len(), -3);
        assert!(Allele::Ins(b"A".to_vec()).kind_matches(&Allele::Ins(b"GG".to_vec())));
        assert!(!Allele::Ins(b"A".to_vec()).kind_matches(&Allele::Del(1)));
        assert!(!Allele::Ref.kind_matches(&Allele::Del(1)));
    }

    #[test]
    fn indel_opts_default_is_disabled() {
        let o = IndelOpts::default();
        assert!(!o.enabled, "indel correction must be opt-in");
        assert_eq!(o.max_len, 5);
    }

    // ref:  0123456789
    //       C A A A A A A A G   (index 0 = C, 1..=7 = A x7, 8 = G)
    const HP_REF: &[u8] = b"CAAAAAAAG";

    #[test]
    fn deletion_in_homopolymer_normalizes_to_one_anchor() {
        // A 1bp deletion of an A can be anchored anywhere in the run (1..=6);
        // every placement must collapse to anchor 0 (the C before the run).
        for anchor in 0..=6u32 {
            let (pos, allele) = normalize_left(HP_REF, anchor, &Allele::Del(1));
            assert_eq!(pos, 0, "anchor {anchor} did not normalize to 0");
            assert_eq!(allele, Allele::Del(1));
        }
    }

    #[test]
    fn insertion_in_homopolymer_normalizes_to_one_anchor() {
        // Inserting an extra A is legal after any A in the run; all collapse to 0.
        for anchor in 0..=7u32 {
            let (pos, allele) = normalize_left(HP_REF, anchor, &Allele::Ins(b"A".to_vec()));
            assert_eq!(pos, 0, "anchor {anchor} did not normalize to 0");
            assert_eq!(allele, Allele::Ins(b"A".to_vec()));
        }
    }

    #[test]
    fn insertion_rotates_when_shifting_left() {
        // ref: C A T A T A T G ; inserting "TA" after index 2 (the T at 2) is the
        // same event as inserting "AT" after index 0, reached by rotating right.
        let r = b"CATATATG";
        let (pos, allele) = normalize_left(r, 2, &Allele::Ins(b"TA".to_vec()));
        assert_eq!(pos, 0);
        assert_eq!(allele, Allele::Ins(b"AT".to_vec()));
    }

    #[test]
    fn unique_sequence_event_does_not_move() {
        // ref: A C G T A C G T ; deleting the G at index 2 is already left-aligned.
        let r = b"ACGTACGT";
        let (pos, allele) = normalize_left(r, 1, &Allele::Del(1));
        assert_eq!(pos, 1);
        assert_eq!(allele, Allele::Del(1));
    }

    #[test]
    fn normalization_stops_at_contig_start() {
        // Whole reference is one homopolymer: shifting can never pass index 0.
        let r = b"AAAAA";
        let (pos, _) = normalize_left(r, 3, &Allele::Del(1));
        assert_eq!(pos, 0);
    }

    #[test]
    fn ref_allele_never_moves() {
        let (pos, allele) = normalize_left(HP_REF, 4, &Allele::Ref);
        assert_eq!(pos, 4);
        assert_eq!(allele, Allele::Ref);
    }
}
```

Add the module declaration to `src/denoise.rs`, immediately after the existing `use rust_htslib::...` line (currently line 6):

```rust
pub(crate) mod indel;
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test denoise::indel:: 2>&1 | tail -20`
Expected: `net_len_and_kind` and `indel_opts_default_is_disabled` PASS; the five normalization tests FAIL with assertion errors like `anchor 1 did not normalize to 0` (the stub returns the anchor unchanged).

- [ ] **Step 3: Implement `normalize_left`**

Replace the stub body in `src/denoise/indel.rs`:

```rust
pub fn normalize_left(refseq: &[u8], anchor: u32, allele: &Allele) -> (u32, Allele) {
    match allele {
        Allele::Ref => (anchor, Allele::Ref),

        // Deleted bases are refseq[a+1 ..= a+n]. Shifting the anchor left by one
        // makes them refseq[a ..= a+n-1], which is the same edit exactly when the
        // base entering the deletion equals the base leaving it.
        Allele::Del(n) => {
            let n = *n as usize;
            let mut a = anchor as usize;
            while a > 0 {
                match (refseq.get(a), refseq.get(a + n)) {
                    (Some(entering), Some(leaving)) if entering == leaving => a -= 1,
                    _ => break,
                }
            }
            (a as u32, Allele::Del(n as u32))
        }

        // `s` is inserted after refseq[a]. Shifting left by one puts refseq[a] after
        // the insertion instead, which is the same edit exactly when the last base of
        // `s` equals refseq[a]; the inserted string rotates right to stay equivalent.
        Allele::Ins(s) => {
            let mut a = anchor as usize;
            let mut s = s.clone();
            while a > 0 && !s.is_empty() {
                match refseq.get(a) {
                    Some(&b) if b == *s.last().unwrap() => {
                        s.rotate_right(1);
                        a -= 1;
                    }
                    _ => break,
                }
            }
            (a as u32, Allele::Ins(s))
        }
    }
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test denoise::indel:: 2>&1 | tail -20`
Expected: all 7 tests PASS.

- [ ] **Step 5: Build and commit**

```bash
cargo build --release
git add src/denoise/indel.rs src/denoise.rs
git commit -m "feat(denoise): indel allele types, IndelOpts, left-normalization"
```

---

### Task 2: Repeat context and the context-scaled error rate

`L` is the reference repeat context at a normalized site — the homopolymer run length for 1 bp events, the tandem-copy count for longer ones. It drives `eps(L)` (the per-read likelihood in Task 3) and `vaf_floor(L)` (the candidacy gate).

**Files:**
- Modify: `src/denoise/indel.rs`

**Interfaces:**
- Consumes: `Allele`, `IndelOpts` from Task 1.
- Produces: `pub fn repeat_context(refseq: &[u8], norm_pos: u32, allele: &Allele) -> u32`; `pub fn error_rate(l: u32, o: &IndelOpts) -> f64`; `pub fn vaf_floor(l: u32, o: &IndelOpts) -> f64`.

- [ ] **Step 1: Write the failing tests**

Append inside the existing `mod tests` block in `src/denoise/indel.rs`:

```rust
    #[test]
    fn repeat_context_counts_homopolymer_run() {
        // HP_REF = C AAAAAAA G -> a 1bp A event normalized to anchor 0 sits in a 7-run.
        assert_eq!(repeat_context(HP_REF, 0, &Allele::Del(1)), 7);
        assert_eq!(repeat_context(HP_REF, 0, &Allele::Ins(b"A".to_vec())), 7);
    }

    #[test]
    fn repeat_context_counts_tandem_copies() {
        // ref: C ATATATAT G -> a 2bp "AT" event at anchor 0 sits in 4 copies.
        let r = b"CATATATATG";
        assert_eq!(repeat_context(r, 0, &Allele::Ins(b"AT".to_vec())), 4);
        assert_eq!(repeat_context(r, 0, &Allele::Del(2)), 4);
    }

    #[test]
    fn repeat_context_is_one_in_unique_sequence() {
        let r = b"ACGTACGT";
        // Deleting the single G at index 2 (anchor 1): "G" occurs once here.
        assert_eq!(repeat_context(r, 1, &Allele::Del(1)), 1);
        assert_eq!(repeat_context(r, 1, &Allele::Ref), 1);
    }

    #[test]
    fn error_rate_grows_with_context_and_is_capped() {
        let o = IndelOpts::default(); // err0 0.01, scale 1.5, cap 0.4
        assert!((error_rate(1, &o) - 0.01).abs() < 1e-12);
        assert!((error_rate(2, &o) - 0.015).abs() < 1e-12);
        assert!(error_rate(3, &o) > error_rate(2, &o));
        // Monotonic non-decreasing, never above the cap.
        for l in 1..40u32 {
            assert!(error_rate(l + 1, &o) >= error_rate(l, &o), "not monotonic at {l}");
            assert!(error_rate(l, &o) <= o.err_cap + 1e-12, "cap exceeded at {l}");
        }
        // Deep enough runs saturate exactly at the cap.
        assert!((error_rate(30, &o) - o.err_cap).abs() < 1e-12);
    }

    #[test]
    fn vaf_floor_is_the_larger_of_absolute_and_error_scaled() {
        let o = IndelOpts::default(); // vaf 0.05, floor_mult 3.0
        // L=1: 3 * 0.01 = 0.03 < 0.05, so the absolute floor wins.
        assert!((vaf_floor(1, &o) - 0.05).abs() < 1e-12);
        // Deep homopolymer: 3 * 0.4 = 1.2 dominates.
        assert!(vaf_floor(30, &o) > 1.0);
        for l in 1..30u32 {
            assert!(vaf_floor(l + 1, &o) >= vaf_floor(l, &o), "not monotonic at {l}");
        }
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test denoise::indel:: 2>&1 | tail -20`
Expected: FAIL to compile — `cannot find function 'repeat_context' in this scope` (and likewise `error_rate`, `vaf_floor`).

- [ ] **Step 3: Implement the three functions**

Append to `src/denoise/indel.rs`, after `normalize_left`:

```rust
/// Reference repeat context `L` at a normalized site: how many tandem copies of the
/// event's own sequence surround it. For a 1bp event this is the homopolymer run
/// length. Unique sequence gives 1. This is where ONT's indel error mass concentrates,
/// so it scales both the error rate and the candidacy floor.
pub fn repeat_context(refseq: &[u8], norm_pos: u32, allele: &Allele) -> u32 {
    let unit: Vec<u8> = match allele {
        Allele::Ref => return 1,
        Allele::Ins(s) => s.clone(),
        Allele::Del(n) => {
            let a = norm_pos as usize + 1;
            match refseq.get(a..a + *n as usize) {
                Some(w) => w.to_vec(),
                None => return 1,
            }
        }
    };
    if unit.is_empty() {
        return 1;
    }
    let u = unit.len();
    let start = norm_pos as usize + 1;

    // Copies running rightward from the event position...
    let mut copies = 0u32;
    let mut p = start;
    while refseq.get(p..p + u).map_or(false, |w| w == &unit[..]) {
        copies += 1;
        p += u;
    }
    // ...and leftward. (After left-normalization there is usually nothing to the
    // left, but counting both directions keeps the measure independent of which
    // end the aligner happened to anchor.)
    let mut p = start;
    while p >= u && refseq.get(p - u..p).map_or(false, |w| w == &unit[..]) {
        copies += 1;
        p -= u;
    }
    copies.max(1)
}

/// Per-junction indel error probability, growing with repeat context and capped.
/// This is the likelihood term that base quality supplies for substitutions and that
/// indels have no equivalent for: without it, MAP assignment degenerates to
/// "everything becomes the most frequent allele" and heteroplasmy is annihilated.
pub fn error_rate(l: u32, o: &IndelOpts) -> f64 {
    let exp = (l as i32 - 1).max(0);
    (o.err0 * o.err_scale.powi(exp)).min(o.err_cap)
}

/// Candidacy VAF floor for an allele in context `L`: an absolute minimum, raised to a
/// multiple of the local error rate wherever that is higher.
pub fn vaf_floor(l: u32, o: &IndelOpts) -> f64 {
    o.vaf.max(o.floor_mult * error_rate(l, o))
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test denoise::indel:: 2>&1 | tail -20`
Expected: all 12 tests PASS.

- [ ] **Step 5: Build and commit**

```bash
cargo build --release
git add src/denoise/indel.rs
git commit -m "feat(denoise): repeat-context detection and context-scaled indel error rate"
```

---

### Task 3: Site candidacy and two-directional MAP assignment

The decision core. Candidacy reuses the existing SNV gates (per-strand counts, strand-bias binomial with the near-homoplasmic exemption); assignment adds the `eps(L)` likelihood that makes symmetric correction safe.

**Files:**
- Modify: `src/denoise/indel.rs`

**Interfaces:**
- Consumes: `Allele`, `IndelOpts`, `repeat_context`, `error_rate`, `vaf_floor`; `crate::call::strand_bias_pvalue(fwd: usize, rev: usize, expected_fwd_frac: f64) -> f64`.
- Produces: `pub struct IndelSite { pub obs: HashMap<Allele, (u32, u32)>, pub depth: u32, pub col_fwd: u32, pub col_rev: u32 }` with `pub fn add(&mut self, allele: Allele, reverse: bool)`; `pub struct IndelSiteModel { pub kept: Vec<(Allele, f64)>, pub eps: f64, pub context_len: u32, pub strand_biased: usize }`; `pub enum Assignment { Keep, Move(Allele), Refused }`; `pub fn fit_indel_site(site: &IndelSite, refseq: &[u8], norm_pos: u32, min_strand: u32, strand_bias_p: f64, homoplasmic_vaf: f64, o: &IndelOpts) -> IndelSiteModel`; `pub fn assign_allele(observed: &Allele, m: &IndelSiteModel, o: &IndelOpts) -> Assignment`.

`IndelSite` carries only ALT counts from `add()`. Its `depth`, `col_fwd`, and `col_rev` are filled in by the caller from the pileup's per-column totals (Task 8) — a read must be counted **once** at a site, and a read carrying an event reports `Indel::None` at the site's own column as well as its event at the anchor column, so accumulating depth from observations would double-count it.

`Assignment` distinguishes "already MAP" from "MAP target refused by a rewrite restriction" so Task 8 can count `reads_skipped_unsupported_transition` instead of silently conflating the two.

- [ ] **Step 1: Write the failing tests**

Append inside `mod tests` in `src/denoise/indel.rs`:

```rust
    // Build a site: `events` is (allele, fwd_count, rev_count). REF support is
    // whatever depth remains after the events, split evenly across strands.
    fn site(depth: u32, events: &[(Allele, u32, u32)]) -> IndelSite {
        let mut s = IndelSite { depth, ..Default::default() };
        let mut used_f = 0;
        let mut used_r = 0;
        for (a, f, r) in events {
            for _ in 0..*f { s.add(a.clone(), false); }
            for _ in 0..*r { s.add(a.clone(), true); }
            used_f += f;
            used_r += r;
        }
        // Remaining depth votes REF, strand-balanced.
        let rest = depth - (used_f + used_r);
        s.col_fwd = used_f + rest / 2;
        s.col_rev = used_r + (rest - rest / 2);
        s
    }

    const UNIQ_REF: &[u8] = b"ACGTACGTACGTACGT";
    const SB_P: f64 = 0.01;
    const HOM_VAF: f64 = 0.7;

    fn kept_of<'a>(m: &'a IndelSiteModel, a: &Allele) -> Option<&'a (Allele, f64)> {
        m.kept.iter().find(|(k, _)| k == a)
    }

    #[test]
    fn balanced_indel_in_unique_sequence_is_kept() {
        // 20% 1bp deletion, strand-balanced, in unique sequence (L=1, floor 0.05).
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 10, 10)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_some(), "balanced 20% del in unique seq must be kept");
        assert!(kept_of(&m, &Allele::Ref).is_some(), "REF is always kept");
    }

    #[test]
    fn same_indel_inside_long_homopolymer_is_rejected() {
        // Identical evidence, but now in a 7-mer A-run: floor = 3 * eps(7) which is
        // far above 0.20, so the context gate rejects it.
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 10, 10)]);
        let o = IndelOpts::default();
        assert!(vaf_floor(7, &o) > 0.20, "test premise: floor must exceed the 20% VAF");
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_none(), "20% del in a 7-mer run must be rejected");
    }

    #[test]
    fn strand_skewed_indel_is_rejected_and_counted() {
        // 20 alt reads, 18/2 split, against a balanced column -> strand-biased.
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 18, 2)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_none());
        assert_eq!(m.strand_biased, 1);
    }

    #[test]
    fn near_homoplasmic_indel_is_exempt_from_strand_bias() {
        // 80% of the column, badly skewed -- but a single-strand artifact cannot
        // reach 80%, so the exemption keeps it (mirrors the SNV model).
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 75, 5)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_some(), "near-homoplasmic allele must survive");
        assert_eq!(m.strand_biased, 0);
    }

    #[test]
    fn low_count_indel_fails_the_per_strand_gate() {
        // Only one read on the reverse strand: fails min_strand=2 regardless of VAF.
        let ins = Allele::Ins(b"T".to_vec());
        let s = site(100, &[(ins.clone(), 19, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &ins).is_none());
    }

    // ---- assignment (the three rows of the design's table) ----

    #[test]
    fn heteroplasmic_site_preserves_both_alleles() {
        // 40% INS / 60% REF in a long homopolymer (large eps). Nothing may flip.
        let ins = Allele::Ins(b"A".to_vec());
        let s = site(100, &[(ins.clone(), 20, 20)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0; // isolate assignment from the gate for this test
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &ins).is_some() && kept_of(&m, &Allele::Ref).is_some());
        assert_eq!(assign_allele(&ins, &m, &o), Assignment::Keep, "40% allele must not flip");
        assert_eq!(assign_allele(&Allele::Ref, &m, &o), Assignment::Keep, "60% allele must not flip");
    }

    #[test]
    fn near_homoplasmic_site_repairs_reads_that_dropped_the_indel() {
        // 90% INS / 10% REF in a long homopolymer: REF reads are dropped-indel errors.
        let ins = Allele::Ins(b"A".to_vec());
        let s = site(100, &[(ins.clone(), 45, 45)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;
        let m = fit_indel_site(&s, HP_REF, 0, 2, SB_P, HOM_VAF, &o);
        assert!(m.eps > 0.1, "test premise: long homopolymer must have a large eps");
        assert_eq!(assign_allele(&Allele::Ref, &m, &o), Assignment::Move(ins.clone()));
        assert_eq!(assign_allele(&ins, &m, &o), Assignment::Keep);
    }

    #[test]
    fn failed_allele_snaps_to_the_kept_allele_of_the_same_kind() {
        // Site keeps a 1bp deletion at 80%; a stray 2bp deletion fails candidacy and
        // must snap onto the kept 1bp allele rather than reverting all the way to REF.
        // (The kept allele has to DOMINATE for snapping to beat reverting: at a site
        // that is mostly REF, reverting a stray deletion to REF is the correct MAP
        // answer, not a bug.)
        let d1 = Allele::Del(1);
        let d2 = Allele::Del(2);
        let s = site(100, &[(d1.clone(), 40, 40), (d2.clone(), 1, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &d1).is_some());
        assert!(kept_of(&m, &d2).is_none());
        assert_eq!(assign_allele(&d2, &m, &o), Assignment::Move(d1.clone()));
    }

    #[test]
    fn cross_kind_reassignment_is_refused() {
        // Site keeps only an insertion; a failed deletion would MAP across kinds,
        // which the rewrite walk does not support -- the read is left untouched.
        let ins = Allele::Ins(b"T".to_vec());
        let del = Allele::Del(1);
        let s = site(100, &[(ins.clone(), 45, 45), (del.clone(), 1, 1)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert!(kept_of(&m, &del).is_none());
        // MAP over {REF, Ins} for an observed Del: whichever wins, it is either REF
        // (allowed) or Ins (refused). Assert we never emit the cross-kind Ins.
        assert_ne!(assign_allele(&del, &m, &o), Assignment::Move(ins.clone()));
    }

    #[test]
    fn growing_a_deletion_is_refused() {
        // Site keeps a 3bp deletion; a read carrying a 1bp deletion would have to
        // grow, consuming reference from the next CIGAR op. Refused.
        let d1 = Allele::Del(1);
        let d3 = Allele::Del(3);
        let s = site(100, &[(d3.clone(), 45, 45), (d1.clone(), 1, 1)]);
        let mut o = IndelOpts::default();
        o.vaf = 0.01;
        o.floor_mult = 0.0;
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert_ne!(assign_allele(&d1, &m, &o), Assignment::Move(d3.clone()));
    }

    #[test]
    fn lone_spurious_indel_at_a_ref_only_site_is_still_corrected() {
        // THE most common real case: a single read carries a random indel error, so
        // no alt allele clears candidacy and REF is the sole survivor. The stray read
        // must still be reverted. (An earlier design returned Keep whenever fewer
        // than two alleles survived, which silently disabled the denoiser's primary
        // job; the Z normalization must therefore include the observed allele.)
        let del = Allele::Del(1);
        let s = site(100, &[(del.clone(), 1, 1)]);
        let o = IndelOpts::default();
        let m = fit_indel_site(&s, UNIQ_REF, 1, 2, SB_P, HOM_VAF, &o);
        assert_eq!(m.kept.len(), 1, "only REF should survive");
        // A read already carrying the sole kept allele is untouched...
        assert_eq!(assign_allele(&Allele::Ref, &m, &o), Assignment::Keep);
        // ...but the stray deletion is reverted.
        assert_eq!(assign_allele(&del, &m, &o), Assignment::Move(Allele::Ref));
        // Same for a stray insertion nobody else carries.
        assert_eq!(
            assign_allele(&Allele::Ins(b"T".to_vec()), &m, &o),
            Assignment::Move(Allele::Ref)
        );
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test denoise::indel:: 2>&1 | tail -20`
Expected: FAIL to compile — `cannot find type 'IndelSite' in this scope` (and `IndelSiteModel`, `fit_indel_site`, `assign_allele`).

- [ ] **Step 3: Implement the site model**

Append to `src/denoise/indel.rs`:

```rust
/// Accumulated observations at one left-normalized indel site.
#[derive(Debug, Default, Clone)]
pub struct IndelSite {
    /// Non-reference events only: allele -> (forward count, reverse count).
    pub obs: HashMap<Allele, (u32, u32)>,
    /// Total reads spanning this site (the VAF denominator).
    pub depth: u32,
    /// Column strand totals, used as the null for the strand-bias test. Mito coverage
    /// is frequently asymmetric, so the null is the column's OWN forward fraction
    /// rather than a library-wide constant.
    pub col_fwd: u32,
    pub col_rev: u32,
}

impl IndelSite {
    pub fn add(&mut self, allele: Allele, reverse: bool) {
        let e = self.obs.entry(allele).or_insert((0, 0));
        if reverse {
            e.1 += 1;
        } else {
            e.0 += 1;
        }
    }

    fn support(&self, a: &Allele) -> u32 {
        self.obs.get(a).map_or(0, |&(f, r)| f + r)
    }

    fn alt_total(&self) -> u32 {
        self.obs.values().map(|&(f, r)| f + r).sum()
    }
}

/// The decision for one site: which alleles are kept, their renormalized frequencies,
/// and the context error rate used for MAP assignment.
#[derive(Debug, Clone)]
pub struct IndelSiteModel {
    pub kept: Vec<(Allele, f64)>,
    pub eps: f64,
    /// The repeat context length that produced `eps`; reported in
    /// `reassignments_by_hp_length` so correction inside long runs is auditable.
    pub context_len: u32,
    pub strand_biased: usize,
}

impl IndelSiteModel {
    fn freq(&self, a: &Allele) -> f64 {
        self.kept.iter().find(|(k, _)| k == a).map_or(0.0, |(_, f)| *f)
    }
}

/// Decide which alleles survive at a site and compute their frequencies.
///
/// REF is always kept, mirroring "the reference allele is always eligible" in the
/// substitution model. A non-reference allele must clear three gates: per-strand
/// counts, the context-scaled VAF floor, and the strand-bias binomial test (with the
/// same near-homoplasmic exemption the SNV model and `call::filter_strand_bias` use).
pub fn fit_indel_site(
    site: &IndelSite,
    refseq: &[u8],
    norm_pos: u32,
    min_strand: u32,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
    o: &IndelOpts,
) -> IndelSiteModel {
    let depth = site.depth.max(1);
    let col_total = site.col_fwd + site.col_rev;
    let expected_fwd_frac = if col_total == 0 {
        0.5
    } else {
        site.col_fwd as f64 / col_total as f64
    };

    // One eps per site. Alleles at the same site can in principle sit in different
    // repeat contexts; take the deepest, which is the conservative choice (it raises
    // the error rate, making assignment less willing to move reads).
    let context_len = site
        .obs
        .keys()
        .map(|a| repeat_context(refseq, norm_pos, a))
        .max()
        .unwrap_or(1);
    let eps = error_rate(context_len, o);

    let mut strand_biased = 0usize;
    // Raw counts over kept alleles; REF gets whatever depth the events did not claim.
    let mut raw: Vec<(Allele, f64)> = Vec::new();
    let ref_support = depth.saturating_sub(site.alt_total());
    raw.push((Allele::Ref, ref_support as f64));

    for (allele, &(fwd, rev)) in site.obs.iter() {
        let n = fwd + rev;
        if n == 0 {
            continue;
        }
        if fwd < min_strand || rev < min_strand {
            continue;
        }
        let l = repeat_context(refseq, norm_pos, allele);
        let vaf = n as f64 / depth as f64;
        if vaf < vaf_floor(l, o) {
            continue;
        }
        // Near-homoplasmic exemption: a single-strand artifact is absent from the
        // other strand's reads, so it can never reach a high apparent frequency.
        // Testing an allele carried by most reads would correct away a real variant.
        if strand_bias_p > 0.0
            && vaf < homoplasmic_vaf
            && crate::call::strand_bias_pvalue(fwd as usize, rev as usize, expected_fwd_frac)
                < strand_bias_p
        {
            strand_biased += 1;
            continue;
        }
        raw.push((allele.clone(), n as f64));
    }

    // Renormalize over the kept set, so mass on rejected alleles does not leak.
    let sum: f64 = raw.iter().map(|(_, c)| *c).sum();
    let kept: Vec<(Allele, f64)> = if sum > 0.0 {
        raw.into_iter().map(|(a, c)| (a, c / sum)).collect()
    } else {
        vec![(Allele::Ref, 1.0)]
    };

    IndelSiteModel { kept, eps, context_len, strand_biased }
}

/// Outcome of MAP assignment for one read at one site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Assignment {
    /// The read already carries the MAP allele, or the site is degenerate.
    Keep,
    /// Reassign the read to this allele.
    Move(Allele),
    /// The MAP allele is only reachable through a transition the rewrite walk does
    /// not support (cross-kind, or growing a deletion). The read is left untouched,
    /// but counted separately so the restriction's cost is auditable.
    Refused,
}

/// Maximum-a-posteriori allele for a read observing `observed`.
///
/// The `eps` likelihood is what makes this safe to run in both directions. With
/// frequency alone, `argmax_k f_k` is identical for every read and the whole site
/// collapses to its most frequent allele, annihilating heteroplasmy. Here a read only
/// leaves its own allele when `f_own / f_other < eps / (1 - eps)`.
pub fn assign_allele(observed: &Allele, m: &IndelSiteModel, o: &IndelOpts) -> Assignment {
    let eps = m.eps.clamp(0.0, 0.5);

    // Error mass spreads over the alleles that could have been observed given truth
    // k, decaying in net-length distance.
    let mut best: Option<(&Allele, f64)> = None;
    for (k, f) in &m.kept {
        let like = if k == observed {
            1.0 - eps
        } else {
            // Z normalizes over the kept set PLUS the observed allele. The observed
            // allele is very often NOT kept -- that is precisely the allele we are
            // correcting away. Summing over kept-only makes Z zero whenever REF is
            // the sole survivor, which is the single most common case in real data:
            // one read carrying a lone spurious indel. That read must still be
            // corrected, so it has to appear in its own normalization.
            let mut z: f64 = m
                .kept
                .iter()
                .filter(|(a, _)| a != k)
                .map(|(a, _)| o.delta.powi((a.net_len() - k.net_len()).abs()))
                .sum();
            if !m.kept.iter().any(|(a, _)| a == observed) {
                z += o.delta.powi((observed.net_len() - k.net_len()).abs());
            }
            if z <= 0.0 {
                continue;
            }
            eps * o.delta.powi((observed.net_len() - k.net_len()).abs()) / z
        };
        let val = f * like;
        if best.map_or(true, |(_, b)| val > b) {
            best = Some((k, val));
        }
    }

    let target = match best {
        Some((t, _)) => t,
        None => return Assignment::Keep,
    };
    if target == observed {
        return Assignment::Keep;
    }
    // Restriction 1: the rewrite walk cannot turn an insertion into a deletion or
    // vice versa. Restriction 2: growing a deletion would consume reference from the
    // following CIGAR op. Both leave the read untouched, counted as Refused.
    match (observed, target) {
        (Allele::Ref, _) | (_, Allele::Ref) => {}
        (a, b) if !a.kind_matches(b) => return Assignment::Refused,
        _ => {}
    }
    if let (Allele::Del(n), Allele::Del(mm)) = (observed, target) {
        if mm > n {
            return Assignment::Refused;
        }
    }
    Assignment::Move(target.clone())
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test denoise::indel:: 2>&1 | tail -25`
Expected: all 23 tests PASS.

- [ ] **Step 5: Build and commit**

```bash
cargo build --release
git add src/denoise/indel.rs
git commit -m "feat(denoise): indel site candidacy and two-directional MAP assignment"
```

---

### Task 4: The reference-coordinate rewrite walk

The only length-changing code in the denoiser. It rebuilds SEQ, QUAL, and CIGAR in one left-to-right pass, and its correctness rests entirely on the reference-consumption invariant.

**Files:**
- Create: `src/denoise/rewrite.rs`
- Modify: `src/denoise.rs` (add `mod rewrite;` beside `mod indel;`)

**Interfaces:**
- Consumes: `Allele` from `super::indel`.
- Produces: `pub struct IndelEdit { pub ref_pos: u32, pub from: Allele, pub to: Allele }`; `pub struct ReadEdits { pub subs: Vec<(u32, u8)>, pub indels: Vec<IndelEdit> }` (`Default`); `pub struct RewriteResult { pub seq: Vec<u8>, pub qual: Vec<u8>, pub cigar: CigarString, pub len_changed: bool }`; `pub fn rewrite_read(pos: i64, cigar: &CigarString, seq: &[u8], qual: &[u8], edits: &ReadEdits, refseq: &[u8]) -> RewriteResult`; `pub fn ref_consumed(c: &CigarString) -> u32`.

- [ ] **Step 1: Create the module with a stub and the failing tests**

Create `src/denoise/rewrite.rs`:

```rust
//! Length-changing read rewrite for the ONT denoiser.
//!
//! Every transition here preserves the read's total reference consumption:
//! `Ref->Ins` adds query only; `Ref->Del` converts M(n) to D(n) (both consume n
//! reference bases); reverting `Ins` removes query only; reverting `Del` converts
//! D(n) back to M(n). POS, reference end, and coordinate-sort order are therefore
//! invariant, so the output never needs re-sorting or re-indexing.

use std::collections::HashMap;

use rust_htslib::bam::record::{Cigar, CigarString};

use super::indel::Allele;

/// One read's reassignment at one normalized indel site.
#[derive(Debug, Clone)]
pub struct IndelEdit {
    /// Reference position of the last base BEFORE the event (the pileup anchor).
    pub ref_pos: u32,
    pub from: Allele,
    pub to: Allele,
}

/// All edits for one read, keyed by REFERENCE coordinate so substitutions and
/// length-changing indel edits can be applied in a single ordered walk.
#[derive(Debug, Default, Clone)]
pub struct ReadEdits {
    pub subs: Vec<(u32, u8)>,
    pub indels: Vec<IndelEdit>,
}

impl ReadEdits {
    pub fn is_empty(&self) -> bool {
        self.subs.is_empty() && self.indels.is_empty()
    }
}

pub struct RewriteResult {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub cigar: CigarString,
    pub len_changed: bool,
}

/// Total reference bases consumed by a CIGAR (M/=/X/D/N). The rewrite invariant.
pub fn ref_consumed(c: &CigarString) -> u32 {
    c.iter()
        .map(|op| match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) | Cigar::RefSkip(n) => *n,
            _ => 0,
        })
        .sum()
}

fn push_op(ops: &mut Vec<Cigar>, op: Cigar) {
    let len = match op {
        Cigar::Match(n) | Cigar::Ins(n) | Cigar::Del(n) | Cigar::RefSkip(n)
        | Cigar::SoftClip(n) | Cigar::HardClip(n) | Cigar::Pad(n) | Cigar::Equal(n)
        | Cigar::Diff(n) => n,
    };
    if len == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        match (&last, &op) {
            (Cigar::Match(a), Cigar::Match(b)) => { *last = Cigar::Match(a + b); return; }
            (Cigar::Ins(a), Cigar::Ins(b)) => { *last = Cigar::Ins(a + b); return; }
            (Cigar::Del(a), Cigar::Del(b)) => { *last = Cigar::Del(a + b); return; }
            (Cigar::SoftClip(a), Cigar::SoftClip(b)) => { *last = Cigar::SoftClip(a + b); return; }
            _ => {}
        }
    }
    ops.push(op);
}

pub fn rewrite_read(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    qual: &[u8],
    edits: &ReadEdits,
    refseq: &[u8],
) -> RewriteResult {
    let _ = (pos, edits, refseq);
    RewriteResult {
        seq: seq.to_vec(),
        qual: qual.to_vec(),
        cigar: cigar.clone(),
        len_changed: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference used by every test below (index 0-based):
    //   0123456789...
    //   A C G T A A A A A C G T A C G T
    const REF: &[u8] = b"ACGTAAAAACGTACGT";

    fn cig(ops: Vec<Cigar>) -> CigarString {
        CigarString(ops)
    }

    /// Assert the design's core invariant plus SEQ/QUAL length agreement.
    fn assert_invariants(before: &CigarString, r: &RewriteResult) {
        assert_eq!(
            ref_consumed(before),
            ref_consumed(&r.cigar),
            "reference consumption changed: {before:?} -> {:?}",
            r.cigar
        );
        assert_eq!(r.seq.len(), r.qual.len(), "SEQ and QUAL lengths disagree");
        let query: u32 = r
            .cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) | Cigar::Equal(n)
                | Cigar::Diff(n) => *n,
                _ => 0,
            })
            .sum();
        assert_eq!(query as usize, r.seq.len(), "CIGAR query length != SEQ length");
    }

    #[test]
    fn no_edits_is_a_faithful_copy() {
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let r = rewrite_read(0, &c, &seq, &qual, &ReadEdits::default(), REF);
        assert_eq!(r.seq, seq);
        assert_eq!(r.qual, qual);
        assert_eq!(r.cigar, c);
        assert!(!r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn substitution_only_preserves_length_and_cigar() {
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits { subs: vec![(2, b'T')], indels: vec![] };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACTTAAAA".to_vec());
        assert_eq!(r.qual, qual);
        assert!(!r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn reverting_an_insertion_removes_query_bases() {
        // Read = ACGT + "GG" inserted after ref index 3 + AAAA.
        let c = cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)]);
        let seq = b"ACGTGGAAAA".to_vec();
        let qual = vec![30u8; 10];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ins(b"GG".to_vec()), to: Allele::Ref }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(8)]));
        assert!(r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn reverting_a_deletion_splices_reference_bases_back() {
        // Read skips ref indices 4..6 (AA), i.e. 4M 2D 3M over ACGT|AA|AAA...
        let c = cig(vec![Cigar::Match(4), Cigar::Del(2), Cigar::Match(3)]);
        let seq = b"ACGTAAA".to_vec();
        let qual = vec![20u8, 21, 22, 23, 24, 25, 26];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Del(2), to: Allele::Ref }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(9)]));
        // Spliced bases take min(left flank, right flank) = min(23, 24) = 23.
        assert_eq!(r.qual, vec![20, 21, 22, 23, 23, 23, 24, 25, 26]);
        assert!(r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn gaining_an_insertion_splits_a_match_run() {
        // Plain 8M read gains a consensus "A" insertion after ref index 3.
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![20u8, 21, 22, 23, 24, 25, 26, 27];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4)]));
        // Inserted base takes min(23, 24) = 23.
        assert_eq!(r.qual, vec![20, 21, 22, 23, 23, 24, 25, 26, 27]);
        assert!(r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn gaining_a_deletion_converts_match_to_del() {
        // Plain 8M read gains a 2bp deletion after ref index 3.
        let c = cig(vec![Cigar::Match(8)]);
        let seq = b"ACGTAAAA".to_vec();
        let qual = vec![30u8; 8];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Del(2) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Del(2), Cigar::Match(2)]));
        assert!(r.len_changed);
        assert_invariants(&c, &r);
    }

    #[test]
    fn replacing_an_insertion_swaps_the_inserted_bases() {
        let c = cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)]);
        let seq = b"ACGTGGAAAA".to_vec();
        let qual = vec![30u8; 10];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit {
                ref_pos: 3,
                from: Allele::Ins(b"GG".to_vec()),
                to: Allele::Ins(b"A".to_vec()),
            }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ACGTAAAAA".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4)]));
        assert_invariants(&c, &r);
    }

    #[test]
    fn shrinking_a_deletion_emits_del_then_match() {
        // 3bp deletion shrinks to 1bp: Del(1) first (left-aligned), then the two
        // reference bases that stop being deleted come back as matches.
        let c = cig(vec![Cigar::Match(4), Cigar::Del(3), Cigar::Match(2)]);
        let seq = b"ACGTAC".to_vec();
        let qual = vec![30u8; 6];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Del(3), to: Allele::Del(1) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        // ref[5..7] = "AA" return to the read.
        assert_eq!(r.seq, b"ACGTAAAC".to_vec());
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(4)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn one_read_carrying_sub_insertion_and_deletion_edits() {
        // 12M read: substitution at ref 1, gained insertion after ref 3, gained
        // deletion after ref 8. All applied in a single walk, left to right.
        let c = cig(vec![Cigar::Match(12)]);
        let seq = b"ACGTAAAAACGT".to_vec();
        let qual = vec![30u8; 12];
        let edits = ReadEdits {
            subs: vec![(1, b'T')],
            indels: vec![
                IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) },
                IndelEdit { ref_pos: 8, from: Allele::Ref, to: Allele::Del(2) },
            ],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"ATGTAAAAAAT".to_vec());
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(5), Cigar::Del(2), Cigar::Match(2)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn soft_clips_are_copied_verbatim() {
        let c = cig(vec![Cigar::SoftClip(3), Cigar::Match(8), Cigar::SoftClip(2)]);
        let seq = b"TTTACGTAAAAGG".to_vec();
        let qual = vec![30u8; 13];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 3, from: Allele::Ref, to: Allele::Ins(b"A".to_vec()) }],
        };
        let r = rewrite_read(0, &c, &seq, &qual, &edits, REF);
        assert_eq!(&r.seq[..3], b"TTT");
        assert_eq!(&r.seq[r.seq.len() - 2..], b"GG");
        assert_eq!(
            r.cigar,
            cig(vec![Cigar::SoftClip(3), Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4), Cigar::SoftClip(2)])
        );
        assert_invariants(&c, &r);
    }

    #[test]
    fn read_starting_at_nonzero_pos_uses_absolute_reference_coordinates() {
        // Read aligned at POS 4 (the AAAAA run). Edit anchored at absolute ref 5.
        let c = cig(vec![Cigar::Match(6)]);
        let seq = b"AAAAAC".to_vec();
        let qual = vec![30u8; 6];
        let edits = ReadEdits {
            subs: vec![],
            indels: vec![IndelEdit { ref_pos: 5, from: Allele::Ref, to: Allele::Del(1) }],
        };
        let r = rewrite_read(4, &c, &seq, &qual, &edits, REF);
        assert_eq!(r.seq, b"AAAAC".to_vec());
        assert_eq!(r.cigar, cig(vec![Cigar::Match(2), Cigar::Del(1), Cigar::Match(3)]));
        assert_invariants(&c, &r);
    }
}
```

Add to `src/denoise.rs`, beside the `mod indel;` line from Task 1:

```rust
pub(crate) mod rewrite;
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test denoise::rewrite:: 2>&1 | tail -25`
Expected: `no_edits_is_a_faithful_copy` PASSES (the stub copies); the other 10 FAIL with SEQ/CIGAR mismatches (the stub ignores all edits).

- [ ] **Step 3: Implement `rewrite_read`**

Replace the stub body in `src/denoise/rewrite.rs`:

```rust
pub fn rewrite_read(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    qual: &[u8],
    edits: &ReadEdits,
    refseq: &[u8],
) -> RewriteResult {
    let subs: HashMap<u32, u8> = edits.subs.iter().copied().collect();
    let indels: HashMap<u32, &IndelEdit> =
        edits.indels.iter().map(|e| (e.ref_pos, e)).collect();

    let mut out_seq: Vec<u8> = Vec::with_capacity(seq.len() + 8);
    let mut out_qual: Vec<u8> = Vec::with_capacity(qual.len() + 8);
    let mut ops: Vec<Cigar> = Vec::new();

    let mut rp: u32 = pos.max(0) as u32; // reference position of the next ref-consuming base
    let mut qp: usize = 0; // query position of the next original base
    let mut last_ref: u32 = rp; // ref position of the most recent aligned base

    // Quality assigned to synthesized bases: never more confident than the flanks.
    let fill_qual = |out_qual: &Vec<u8>, next: Option<u8>| -> u8 {
        match (out_qual.last().copied(), next) {
            (Some(l), Some(r)) => l.min(r),
            (Some(l), None) => l,
            (None, Some(r)) => r,
            (None, None) => 0,
        }
    };

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let mut k: u32 = 0;
                while k < n {
                    let cur = rp + k;
                    let qi = qp + k as usize;
                    out_seq.push(*subs.get(&cur).unwrap_or(&seq[qi]));
                    out_qual.push(qual[qi]);
                    push_op(&mut ops, Cigar::Match(1));
                    last_ref = cur;
                    k += 1;

                    // An indel anchored on this aligned base?
                    if let Some(e) = indels.get(&cur) {
                        match (&e.from, &e.to) {
                            (Allele::Ref, Allele::Ins(s)) => {
                                let next = qual.get(qp + k as usize).copied();
                                let f = fill_qual(&out_qual, next);
                                out_seq.extend_from_slice(s);
                                out_qual.extend(std::iter::repeat(f).take(s.len()));
                                push_op(&mut ops, Cigar::Ins(s.len() as u32));
                            }
                            (Allele::Ref, Allele::Del(m)) => {
                                // The following m reference bases become a deletion.
                                // The flank guard in denoise.rs guarantees they lie
                                // inside this same match run.
                                debug_assert!(k + m <= n, "Ref->Del ran past its match op");
                                let m = (*m).min(n - k);
                                push_op(&mut ops, Cigar::Del(m));
                                k += m; // skip both the ref bases and their query partners
                            }
                            _ => {}
                        }
                    }
                }
                qp += n as usize;
                rp += n;
            }

            Cigar::Ins(n) => {
                let e = indels
                    .get(&last_ref)
                    .filter(|e| matches!(e.from, Allele::Ins(_)));
                match e.map(|e| &e.to) {
                    // Reverted: drop the inserted bases entirely.
                    Some(Allele::Ref) => {}
                    // Replaced: emit the site's consensus inserted bases instead.
                    Some(Allele::Ins(s)) => {
                        let next = qual.get(qp + n as usize).copied();
                        let f = fill_qual(&out_qual, next);
                        out_seq.extend_from_slice(s);
                        out_qual.extend(std::iter::repeat(f).take(s.len()));
                        push_op(&mut ops, Cigar::Ins(s.len() as u32));
                    }
                    // No edit, or an unsupported cross-kind target: copy verbatim.
                    _ => {
                        out_seq.extend_from_slice(&seq[qp..qp + n as usize]);
                        out_qual.extend_from_slice(&qual[qp..qp + n as usize]);
                        push_op(&mut ops, Cigar::Ins(n));
                    }
                }
                qp += n as usize;
            }

            Cigar::Del(n) => {
                let e = indels
                    .get(&last_ref)
                    .filter(|e| matches!(e.from, Allele::Del(_)));
                let splice = |out_seq: &mut Vec<u8>, out_qual: &mut Vec<u8>, from: u32, count: u32| {
                    let f = match (out_qual.last().copied(), qual.get(qp).copied()) {
                        (Some(l), Some(r)) => l.min(r),
                        (Some(l), None) => l,
                        (None, Some(r)) => r,
                        (None, None) => 0,
                    };
                    for j in 0..count {
                        let idx = (from + j) as usize;
                        out_seq.push(refseq.get(idx).copied().unwrap_or(b'N'));
                        out_qual.push(f);
                    }
                };
                match e.map(|e| &e.to) {
                    // Reverted: the deleted reference bases come back as matches.
                    Some(Allele::Ref) => {
                        splice(&mut out_seq, &mut out_qual, rp, n);
                        push_op(&mut ops, Cigar::Match(n));
                    }
                    // Shrunk: left-aligned, so the shorter deletion comes first and
                    // the freed reference bases follow as matches.
                    Some(Allele::Del(m)) if *m <= n => {
                        push_op(&mut ops, Cigar::Del(*m));
                        let back = n - *m;
                        if back > 0 {
                            splice(&mut out_seq, &mut out_qual, rp + *m, back);
                            push_op(&mut ops, Cigar::Match(back));
                        }
                    }
                    // No edit, or an unsupported target: copy verbatim.
                    _ => push_op(&mut ops, Cigar::Del(n)),
                }
                rp += n;
            }

            Cigar::SoftClip(n) => {
                out_seq.extend_from_slice(&seq[qp..qp + n as usize]);
                out_qual.extend_from_slice(&qual[qp..qp + n as usize]);
                push_op(&mut ops, Cigar::SoftClip(n));
                qp += n as usize;
            }
            Cigar::RefSkip(n) => {
                push_op(&mut ops, Cigar::RefSkip(n));
                rp += n;
            }
            Cigar::HardClip(n) => push_op(&mut ops, Cigar::HardClip(n)),
            Cigar::Pad(n) => push_op(&mut ops, Cigar::Pad(n)),
        }
    }

    let len_changed = out_seq.len() != seq.len();
    RewriteResult {
        seq: out_seq,
        qual: out_qual,
        cigar: CigarString(ops),
        len_changed,
    }
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test denoise::rewrite:: 2>&1 | tail -25`
Expected: all 11 tests PASS.

- [ ] **Step 5: Build and commit**

```bash
cargo build --release
git add src/denoise/rewrite.rs src/denoise.rs
git commit -m "feat(denoise): reference-coordinate SEQ/QUAL/CIGAR rewrite walk"
```

---

### Task 5: Stale-tag stripping for length-changed reads

`MM`/`ML`/`MN` index a read by its own base positions, so a length change silently corrupts them. `NM`/`MD`/`cs` describe the old alignment. Both sets must go when a read's length changes.

**Files:**
- Modify: `src/denoise/rewrite.rs`

**Interfaces:**
- Consumes: nothing new.
- Produces: `pub fn strip_stale_tags(rec: &mut rust_htslib::bam::Record)`.

- [ ] **Step 1: Write the failing test**

Append inside `mod tests` in `src/denoise/rewrite.rs`:

```rust
    #[test]
    fn stale_tags_are_stripped_and_others_survive() {
        use rust_htslib::bam::record::{Aux, Record};

        let mut rec = Record::new();
        rec.set(
            b"r1",
            Some(&cig(vec![Cigar::Match(4)])),
            b"ACGT",
            &[30u8; 4],
        );
        rec.push_aux(b"MM", Aux::String("C+m,0;")).unwrap();
        rec.push_aux(b"ML", Aux::U8(200)).unwrap();
        rec.push_aux(b"NM", Aux::I32(2)).unwrap();
        rec.push_aux(b"XY", Aux::I32(7)).unwrap();

        strip_stale_tags(&mut rec);

        assert!(rec.aux(b"MM").is_err(), "MM must be stripped");
        assert!(rec.aux(b"ML").is_err(), "ML must be stripped");
        assert!(rec.aux(b"NM").is_err(), "NM must be stripped");
        assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)), "unrelated tags must survive");
    }

    #[test]
    fn stripping_absent_tags_is_a_no_op() {
        use rust_htslib::bam::record::{Aux, Record};

        let mut rec = Record::new();
        rec.set(b"r1", Some(&cig(vec![Cigar::Match(4)])), b"ACGT", &[30u8; 4]);
        rec.push_aux(b"XY", Aux::I32(7)).unwrap();
        strip_stale_tags(&mut rec); // must not panic on missing tags
        assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)));
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test denoise::rewrite::tests::stale 2>&1 | tail -15`
Expected: FAIL to compile — `cannot find function 'strip_stale_tags' in this scope`.

- [ ] **Step 3: Implement `strip_stale_tags`**

Append to `src/denoise/rewrite.rs`, after `rewrite_read`:

```rust
/// Remove tags invalidated by a length change.
///
/// `MM`/`ML`/`MN` are base-modification tags indexed by the read's own base
/// positions; a length change silently desynchronizes them. This is safe here
/// because methylation aggregation runs on the ORIGINAL mt BAM (see main.rs
/// QuickStart), not the denoised one. `NM`/`MD`/`cs` describe the pre-edit
/// alignment and are simply stale.
pub fn strip_stale_tags(rec: &mut rust_htslib::bam::Record) {
    for tag in [
        b"MM".as_slice(),
        b"ML".as_slice(),
        b"MN".as_slice(),
        b"NM".as_slice(),
        b"MD".as_slice(),
        b"cs".as_slice(),
    ] {
        // Absent tags are expected; removal failure just means "wasn't there".
        let _ = rec.remove_aux(tag);
    }
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test denoise::rewrite:: 2>&1 | tail -20`
Expected: all 13 tests PASS.

- [ ] **Step 5: Build and commit**

```bash
cargo build --release
git add src/denoise/rewrite.rs
git commit -m "feat(denoise): strip base-modification and alignment tags on length change"
```

---

### Task 6: Migrate `Corrections` to reference coordinates

Pure refactor with no behavior change. Substitutions move from query to reference coordinates and `apply_corrections` starts routing through `rewrite_read`, so Task 8 can add indel edits to the same payload.

**Files:**
- Modify: `src/denoise.rs:259` (the `Corrections` type alias), `src/denoise.rs:275-359` (`compute_corrections`), `src/denoise.rs:362-397` (`apply_corrections`), `src/denoise.rs:398+` (`start`, to pass `refs` through), and the two existing tests that construct or assert on `Corrections`.

**Interfaces:**
- Consumes: `rewrite::{ReadEdits, RewriteResult, rewrite_read, strip_stale_tags}`.
- Produces: `pub type Corrections = HashMap<Vec<u8>, rewrite::ReadEdits>;` and `fn apply_corrections(in_bam: &Path, out_bam: &Path, corrections: &Corrections, refs: &HashMap<Vec<u8>, Vec<u8>>) -> Result<(u64, u64)>`.

- [ ] **Step 1: Write the failing test**

Append inside `mod tests` in `src/denoise.rs`:

```rust
    #[test]
    fn corrections_are_keyed_by_reference_position() {
        // The single G error sits at reference position 2 (reads start at POS 0).
        // Recording it in reference coordinates is what lets Task 7 mix
        // length-changing indel edits into the same per-read payload (Task 8).
        let dir = std::env::temp_dir().join("himito_denoise_refcoord");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], bool)> = vec![];
        let clean = b"AAAAAA";
        for i in 0..9 {
            reads.push((["r0","r1","r2","r3","r4","r5","r6","r7","r8"][i], 0, clean, i % 2 == 0));
        }
        reads.push(("rerr", 0, b"AAGAAA", false));
        write_test_bam(&bam, "chrM", 6, &reads);

        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), b"AAAAAA".to_vec());
        let (corr, _) = compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF).unwrap();

        let edits = corr.get(&b"rerr".to_vec()).expect("rerr must be corrected");
        assert_eq!(edits.subs, vec![(2u32, b'A')]);
        assert!(edits.indels.is_empty());
        std::fs::remove_dir_all(&dir).ok();
    }
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test denoise::tests::corrections_are_keyed 2>&1 | tail -15`
Expected: FAIL to compile — `no field 'subs' on type '&Vec<(usize, u8)>'`.

- [ ] **Step 3: Migrate the type and both passes**

In `src/denoise.rs`, replace the type alias at line 259:

```rust
pub type Corrections = HashMap<Vec<u8>, rewrite::ReadEdits>;
```

Add to the `use` block at the top of the file:

```rust
use crate::denoise::rewrite::{rewrite_read, strip_stale_tags, ReadEdits};
```

In `compute_corrections`, replace the correction-recording block (currently lines 341-353) with reference-coordinate recording:

```rust
        let mut changed = false;
        for (qname, _qpos, allele, qual) in obs {
            let corr = map_allele(allele, qual, &model);
            if corr != allele {
                corrections
                    .entry(qname)
                    .or_default()
                    .subs
                    .push((refpos as u32, idx_to_base(corr)));
                stats.bases_modified += 1;
                stats.substitution_matrix[allele][corr] += 1;
                changed = true;
            }
        }
```

Replace `apply_corrections` in full:

```rust
fn apply_corrections(
    in_bam: &Path,
    out_bam: &Path,
    corrections: &Corrections,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
) -> Result<(u64, u64)> {
    let mut reader = Reader::from_path(in_bam)
        .with_context(|| format!("cannot open BAM {in_bam:?}"))?;
    let header_view = reader.header().to_owned();
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
        .with_context(|| format!("cannot create BAM {out_bam:?}"))?;

    let mut reads_processed = 0u64;
    let mut reads_modified = 0u64;
    let empty: Vec<u8> = Vec::new();

    for result in reader.records() {
        let mut rec = result.context("read BAM record")?;
        // Only primary mapped reads are candidates for correction.
        let correctable =
            !(rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped());
        if correctable {
            reads_processed += 1;
            if let Some(edits) = corrections.get(rec.qname()) {
                if !edits.is_empty() {
                    let refseq = refs
                        .get(header_view.tid2name(rec.tid() as u32))
                        .unwrap_or(&empty);
                    // Extract owned parts, rewrite, rebuild (aux preserved by set()).
                    let qname = rec.qname().to_vec();
                    let cigar = rec.cigar().take();
                    let seq = rec.seq().as_bytes();
                    let qual = rec.qual().to_vec();
                    let out =
                        rewrite_read(rec.pos(), &cigar, &seq, &qual, edits, refseq);
                    rec.set(&qname, Some(&out.cigar), &out.seq, &out.qual);
                    if out.len_changed {
                        strip_stale_tags(&mut rec);
                    }
                    reads_modified += 1;
                }
            }
        }
        writer.write(&rec).context("write BAM record")?;
    }
    Ok((reads_processed, reads_modified))
}
```

In `start`, pass `refs` to the apply pass:

```rust
    let (reads_processed, reads_modified) = apply_corrections(input, output, &corrections, &refs)?;
```

Finally update the two existing tests. In `compute_corrections_flags_a_single_base_error`, replace the correction assertion:

```rust
        let edits = corr.get(&b"rerr".to_vec()).expect("rerr must be corrected");
        assert_eq!(edits.subs, vec![(2u32, b'A')]);
```

In `apply_corrections_rewrites_seq_and_preserves_tags`, replace the corrections construction and call:

```rust
        let mut corr: Corrections = HashMap::new();
        corr.insert(
            b"rerr".to_vec(),
            ReadEdits { subs: vec![(2u32, b'A')], indels: vec![] },
        );
        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), b"AAAAAA".to_vec());
        let (proc, modified) = apply_corrections(&inb, &outb, &corr, &refs).unwrap();
```

- [ ] **Step 4: Run the full denoise suite to verify no behavior change**

Run: `cargo test denoise:: 2>&1 | tail -25`
Expected: all tests PASS, including the pre-existing `apply_corrections_rewrites_seq_and_preserves_tags` (the `XY` tag still survives, because a substitution-only rewrite reports `len_changed == false` and skips tag stripping) and `start_ont_denoises_and_writes_stats_json`.

- [ ] **Step 5: Build and commit**

```bash
cargo build --release
git add src/denoise.rs
git commit -m "refactor(denoise): key corrections by reference position, route writes through rewrite_read"
```

---

### Task 7: Enumerate a read's own normalized indel events from its CIGAR

The pileup tells you what a read has *at a column*; it cannot tell you whether a read carries the event belonging to a site that normalization moved left. That question — "does THIS read carry the event at site s?" — is answerable only per-record, from the read's own CIGAR. This function is what lets assignment happen in the write-back pass instead of needing a third BAM pass.

**Files:**
- Modify: `src/denoise/indel.rs`

**Interfaces:**
- Consumes: `Allele`, `normalize_left` from Task 1.
- Produces: `pub struct ReadEvent { pub norm_pos: u32, pub allele: Option<Allele> }` (`allele: None` marks an indel at or above `max_len` — untouchable, and a signal that nearby sites must not be assigned for this read); `pub fn read_events(pos: i64, cigar: &CigarString, seq: &[u8], refseq: &[u8], max_len: u32) -> Vec<ReadEvent>`.

- [ ] **Step 1: Write the failing tests**

Add this import to the top of `src/denoise/indel.rs`:

```rust
use rust_htslib::bam::record::{Cigar, CigarString};
```

Append inside `mod tests` in `src/denoise/indel.rs`:

```rust
    // Reference shared by the read_events tests. Index map:
    //   0-7   ACGTACGT
    //   8-15  AAAAAAAA   (an 8bp homopolymer run)
    //   16-28 CGTACGTACGTAC
    const REF29: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC";

    fn cs(ops: Vec<Cigar>) -> CigarString {
        CigarString(ops)
    }

    #[test]
    fn read_with_no_indels_yields_no_events() {
        let evs = read_events(0, &cs(vec![Cigar::Match(29)]), REF29, REF29, 5);
        assert!(evs.is_empty());
    }

    #[test]
    fn insertion_event_is_normalized_out_of_the_homopolymer() {
        // Read carries an extra A, anchored by the aligner at reference index 15
        // (the last A of the run). Normalization must pull it back to index 7.
        let seq = b"ACGTACGTAAAAAAAAACGTACGTACGTAC"; // 30bp: 16 + 1 inserted + 13
        let cig = cs(vec![Cigar::Match(16), Cigar::Ins(1), Cigar::Match(13)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].norm_pos, 7);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"A".to_vec())));
    }

    #[test]
    fn deletion_event_is_normalized_out_of_the_homopolymer() {
        // Read is missing one A from the run; the aligner anchored the deletion at
        // reference index 14, normalization pulls it to 7.
        let seq = b"ACGTACGTAAAAAAACGTACGTACGTAC"; // 28bp: 29 ref bases minus one A
        let cig = cs(vec![Cigar::Match(15), Cigar::Del(1), Cigar::Match(13)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].norm_pos, 7);
        assert_eq!(evs[0].allele, Some(Allele::Del(1)));
    }

    #[test]
    fn unique_sequence_insertion_keeps_its_anchor() {
        // Insert "G" after reference index 18 (a T). Nothing to shift through.
        let seq = b"ACGTACGTAAAAAAAACGTGACGTACGTAC"; // 30bp
        let cig = cs(vec![Cigar::Match(19), Cigar::Ins(1), Cigar::Match(10)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].norm_pos, 18);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"G".to_vec())));
    }

    #[test]
    fn oversized_indel_is_reported_as_untouchable() {
        // A 6bp insertion with max_len 5: recorded with allele None so the caller
        // knows this read must not be assigned at nearby sites. Reporting nothing at
        // all would make the read look like a clean REF observation and expose it to
        // being "repaired" into some other allele.
        let seq = b"ACGTACGTAAAAAAAACGTGGGGGGACGTACGTAC"; // 35bp
        let cig = cs(vec![Cigar::Match(19), Cigar::Ins(6), Cigar::Match(10)]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].allele, None);
        assert_eq!(evs[0].norm_pos, 18);
    }

    #[test]
    fn multiple_events_in_one_read_are_all_reported() {
        // An extra A in the run (normalizes to 7) plus a 1bp deletion at index 19.
        let seq = b"ACGTACGTAAAAAAAAACGTCGTACGTAC"; // 29bp
        let cig = cs(vec![
            Cigar::Match(16),
            Cigar::Ins(1),
            Cigar::Match(4),
            Cigar::Del(1),
            Cigar::Match(9),
        ]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 2);
        assert_eq!(evs[0].norm_pos, 7);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"A".to_vec())));
        assert_eq!(evs[1].norm_pos, 19);
        assert_eq!(evs[1].allele, Some(Allele::Del(1)));
    }

    #[test]
    fn soft_clips_do_not_shift_reference_coordinates() {
        // Leading soft clip consumes query but no reference; the event must still
        // land at reference index 18.
        let seq = b"TTTACGTACGTAAAAAAAACGTGACGTACGTAC"; // 3 clipped + 30
        let cig = cs(vec![
            Cigar::SoftClip(3),
            Cigar::Match(19),
            Cigar::Ins(1),
            Cigar::Match(10),
        ]);
        let evs = read_events(0, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].norm_pos, 18);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"G".to_vec())));
    }

    #[test]
    fn read_aligned_at_nonzero_pos_reports_absolute_coordinates() {
        // Read starts at reference index 16; insert "G" after index 18.
        let seq = b"CGTGACGTACGTAC"; // ref[16..19] + "G" + ref[19..29]
        let cig = cs(vec![Cigar::Match(3), Cigar::Ins(1), Cigar::Match(10)]);
        let evs = read_events(16, &cig, seq, REF29, 5);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].norm_pos, 18);
        assert_eq!(evs[0].allele, Some(Allele::Ins(b"G".to_vec())));
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test denoise::indel:: 2>&1 | tail -20`
Expected: FAIL to compile — `cannot find function 'read_events' in this scope`.

- [ ] **Step 3: Implement `read_events`**

Append to `src/denoise/indel.rs`:

```rust
/// One indel carried by a read, in normalized site coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadEvent {
    pub norm_pos: u32,
    /// `None` means the indel is at or above `max_len`: untouchable, and a marker
    /// that this read must not be assigned at nearby sites either. Omitting it
    /// entirely would make the read indistinguishable from a clean REF observation
    /// and expose it to being "repaired" into an allele it cannot represent.
    pub allele: Option<Allele>,
}

/// Enumerate a read's own indel events straight from its CIGAR, left-normalized.
///
/// The pileup reports what a read has AT A COLUMN, which cannot answer "does this
/// read carry the event belonging to the site normalization moved left to?" — at the
/// normalized column the carrier reports `Indel::None` just like a non-carrier. Going
/// through the read's own CIGAR answers it exactly, per record.
pub fn read_events(
    pos: i64,
    cigar: &CigarString,
    seq: &[u8],
    refseq: &[u8],
    max_len: u32,
) -> Vec<ReadEvent> {
    let mut out = Vec::new();
    let mut rp = pos.max(0) as u32; // next reference position
    let mut qp = 0usize; // next query position

    for op in cigar.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                rp += n;
                qp += n as usize;
            }
            Cigar::Ins(n) => {
                // The anchor is the last aligned reference base before the event.
                // An insertion before any aligned base has no anchor and is skipped.
                if rp > 0 {
                    let anchor = rp - 1;
                    if n < max_len {
                        let s = seq[qp..qp + n as usize].to_vec();
                        let (np, a) = normalize_left(refseq, anchor, &Allele::Ins(s));
                        out.push(ReadEvent { norm_pos: np, allele: Some(a) });
                    } else {
                        out.push(ReadEvent { norm_pos: anchor, allele: None });
                    }
                }
                qp += n as usize;
            }
            Cigar::Del(n) => {
                if rp > 0 {
                    let anchor = rp - 1;
                    if n < max_len {
                        let (np, a) = normalize_left(refseq, anchor, &Allele::Del(n));
                        out.push(ReadEvent { norm_pos: np, allele: Some(a) });
                    } else {
                        out.push(ReadEvent { norm_pos: anchor, allele: None });
                    }
                }
                rp += n;
            }
            Cigar::SoftClip(n) => qp += n as usize,
            Cigar::RefSkip(n) => rp += n,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    out
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test denoise::indel:: 2>&1 | tail -25`
Expected: all tests PASS (23 from Tasks 1–3, plus 8 new).

- [ ] **Step 5: Build and commit**

```bash
cargo build --release
git add src/denoise/indel.rs
git commit -m "feat(denoise): enumerate per-read normalized indel events from CIGAR"
```

---

### Task 8: Site accumulation in the pileup pass, per-read assignment in the write-back pass

The orchestration task. Two passes, as today.

**Pass 1** (the existing pileup loop) records per-column depth and strand totals, and accumulates normalized ALT event counts into sites. Indel collection must happen **before** the `non_ref() == 0` early-`continue` — that gate is about non-reference *substitutions*, and a column with a perfectly clean base composition can still carry indel events. Sites are decided after the loop closes, because normalization can move an event into a column already passed.

**Pass 2** (the existing record walk in `apply_corrections`) enumerates each read's own events via `read_events`, determines its allele at every decided site it spans, assigns, and rewrites.

**Files:**
- Modify: `src/denoise.rs` — `DenoiseStats` (lines 261-273), `compute_corrections` (lines 275-359), `apply_corrections`, and `start`

**Interfaces:**
- Consumes: `indel::{Allele, Assignment, IndelOpts, IndelSite, IndelSiteModel, ReadEvent, assign_allele, fit_indel_site, normalize_left, read_events}`; `rewrite::IndelEdit`.
- Produces: `pub type IndelModels = HashMap<Vec<u8>, Vec<(u32, indel::IndelSiteModel)>>` (per contig, **sorted ascending by position** so the record walk can binary-search a read's span); `compute_corrections(..., iopts: &IndelOpts) -> Result<(Corrections, IndelModels, DenoiseStats)>`; `apply_corrections(in_bam, out_bam, corrections, refs, models, iopts, stats: &mut DenoiseStats) -> Result<(u64, u64)>`.

**Why depth comes from the column map, not from counting observations.** A read carrying an insertion normalized to site `s` reports `Indel::None` at column `s` *and* its insertion at the anchor column. Counting both as observations of site `s` double-counts that read — depth reads high and VAF reads roughly half its true value, systematically under-calling real indels. The per-column depth at `s` counts every spanning read exactly once, which is what the VAF denominator must be. ALT support comes from normalized events; `REF` support is `depth − alt_total`.

- [ ] **Step 1: Write the failing tests**

Append inside `mod tests` in `src/denoise.rs`:

```rust
    // Reference shared by the indel wiring tests. Index map:
    //   0-7   ACGTACGT
    //   8-15  AAAAAAAA   (an 8bp homopolymer run)
    //   16-28 CGTACGTACGTAC
    // Long enough that the shipped default --indel-flank 5 is satisfiable at the
    // sites used below (normalized positions 7 and 18).
    const REF29: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC";

    // Build a BAM whose reads carry explicit CIGARs, so indel-bearing reads can be
    // constructed directly. Each read: (qname, pos, seq, cigar ops, reverse).
    fn write_cigar_bam(
        path: &Path,
        contig: &str,
        contig_len: usize,
        reads: &[(&str, i64, &[u8], Vec<Cigar>, bool)],
    ) {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", contig);
        sq.push_tag(b"LN", &contig_len.to_string());
        header.push_record(&sq);
        let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
        for (name, pos, seq, ops, reverse) in reads {
            let mut rec = Record::new();
            let quals = vec![30u8; seq.len()];
            rec.set(name.as_bytes(), Some(&CigarString(ops.clone())), seq, &quals);
            rec.set_tid(0);
            rec.set_pos(*pos);
            rec.set_mapq(60);
            rec.set_flags(if *reverse { 16 } else { 0 });
            w.write(&rec).unwrap();
        }
    }

    fn indel_refs() -> HashMap<Vec<u8>, Vec<u8>> {
        let mut refs = HashMap::new();
        refs.insert(b"chrM".to_vec(), REF29.to_vec());
        refs
    }

    fn indels_on() -> indel::IndelOpts {
        indel::IndelOpts { enabled: true, ..Default::default() }
    }

    const CLEAN_SEQ: &[u8] = b"ACGTACGTAAAAAAAACGTACGTACGTAC"; // == REF29, 29bp
    // An extra A inside the run: 16M 1I 13M. Aligner anchors at index 15;
    // normalization pulls the site to index 7.
    const HP_INS_SEQ: &[u8] = b"ACGTACGTAAAAAAAAACGTACGTACGTAC"; // 30bp
    // A spurious G after index 18, in unique sequence: 19M 1I 10M.
    const UNIQ_INS_SEQ: &[u8] = b"ACGTACGTAAAAAAAACGTGACGTACGTAC"; // 30bp

    fn clean_cigar() -> Vec<Cigar> { vec![Cigar::Match(29)] }
    fn hp_ins_cigar() -> Vec<Cigar> { vec![Cigar::Match(16), Cigar::Ins(1), Cigar::Match(13)] }
    fn uniq_ins_cigar() -> Vec<Cigar> { vec![Cigar::Match(19), Cigar::Ins(1), Cigar::Match(10)] }

    const NAMES: [&str; 9] = ["r0","r1","r2","r3","r4","r5","r6","r7","r8"];

    #[test]
    fn indels_disabled_by_default_leaves_indel_reads_untouched() {
        let dir = std::env::temp_dir().join("himito_denoise_indel_off");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, UNIQ_INS_SEQ, uniq_ins_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let off = indel::IndelOpts::default();
        assert!(!off.enabled);
        let (_corr, models, stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &off).unwrap();
        assert!(models.is_empty(), "no sites may be built while --indels is off");
        assert_eq!(stats.indel_sites_examined, 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn site_depth_counts_each_spanning_read_exactly_once() {
        // The double-count trap: a read carrying the insertion reports Indel::None at
        // the normalized column (7) AND its insertion at the anchor column (15).
        // Depth at site 7 must be 10 -- the number of spanning reads -- not 19.
        let dir = std::env::temp_dir().join("himito_denoise_indel_depth");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        reads.push(("rdrop", 0, CLEAN_SEQ, clean_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, _stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();

        let list = models.get(&b"chrM".to_vec()).expect("chrM must have sites");
        let (_, model) = list.iter().find(|(p, _)| *p == 7).expect("site at 7");
        let ins = indel::Allele::Ins(b"A".to_vec());
        let f_ins = model.kept.iter().find(|(a, _)| *a == ins).expect("insertion kept").1;
        // 9 of 10 reads carry it. A double-counted denominator would give ~0.47.
        assert!(
            (f_ins - 0.9).abs() < 1e-9,
            "insertion frequency should be 0.9, got {f_ins} (depth double-counted?)"
        );
        assert_eq!(model.context_len, 8, "site 7 sits in an 8bp homopolymer run");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn models_are_sorted_by_position_within_each_contig() {
        // The record walk binary-searches this list; unsorted input silently drops
        // sites from every read's span.
        let dir = std::env::temp_dir().join("himito_denoise_indel_sorted");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let both: Vec<Cigar> = vec![
            Cigar::Match(16), Cigar::Ins(1), Cigar::Match(4), Cigar::Ins(1), Cigar::Match(9),
        ];
        // Extra A in the run (site 7) AND a spurious G after index 19 (site 19).
        let both_seq: &[u8] = b"ACGTACGTAAAAAAAAACGTGCGTACGTAC"; // 31bp? see cigar: 16+1+4+1+9 = 31
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, both_seq, both.clone(), i % 2 == 0));
        }
        reads.push(("rclean", 0, CLEAN_SEQ, clean_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let (_corr, models, _stats) =
            compute_corrections(&bam, &indel_refs(), 2, 0.01, SB_P, HOM_VAF, &indels_on()).unwrap();
        let list = models.get(&b"chrM".to_vec()).expect("chrM must have sites");
        assert!(list.len() >= 2, "expected at least two sites, got {}", list.len());
        assert!(
            list.windows(2).all(|w| w[0].0 <= w[1].0),
            "site list must be sorted ascending by position"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn lone_spurious_insertion_is_removed_from_a_substitution_clean_column() {
        // Every aligned base matches the reference, so `counts.non_ref() == 0` at
        // every column and the SNV path skips them all. The lone spurious insertion
        // must STILL be collected and corrected -- this is the early-continue trap.
        // It is also the REF-only-site case: no alt allele clears candidacy, so the
        // site keeps REF alone and the stray read must still be reverted.
        let dir = std::env::temp_dir().join("himito_denoise_indel_cleancol");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, UNIQ_INS_SEQ, uniq_ins_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        assert_eq!(stats.columns_skipped, stats.columns_processed,
            "test premise: every column is substitution-clean");
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rins" {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "spurious insertion must be removed");
                assert_eq!(rec.cigar().to_string(), "29M");
                seen = true;
            } else {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "clean reads must be untouched");
            }
        }
        assert!(seen);
        assert_eq!(stats.reads_reassigned_indel_to_ref, 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn homoplasmic_indel_survives_and_repairs_the_read_that_dropped_it() {
        // 9 of 10 reads carry a 1bp insertion inside the 8bp A-run, so eps is large
        // enough for the repair direction to fire (in unique sequence eps is 0.01 and
        // a 10% REF allele correctly would NOT flip).
        let dir = std::env::temp_dir().join("himito_denoise_indel_repair");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        reads.push(("rdrop", 0, CLEAN_SEQ, clean_cigar(), false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        assert_eq!(stats.indel_sites_kept, 1, "the 90% insertion must survive candidacy");
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut repaired = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            assert_eq!(rec.seq().as_bytes(), HP_INS_SEQ, "every read should carry the insertion");
            if rec.qname() == b"rdrop" {
                repaired = true;
            }
        }
        assert!(repaired);
        assert_eq!(stats.reads_reassigned_ref_to_ins, 1);
        assert_eq!(stats.indel_bases_inserted, 1);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn read_ending_too_close_to_a_site_is_skipped_by_the_flank_guard() {
        // Same 90%-insertion site (normalized position 7), but the read that dropped
        // the insertion is aligned 0..14 -- it ends before 7 + flank(5) + max_len(5)
        // = 17, so the rewrite walk would lack the match context it needs.
        let dir = std::env::temp_dir().join("himito_denoise_indel_guard");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, HP_INS_SEQ, hp_ins_cigar(), i % 2 == 0));
        }
        let short_seq: &[u8] = b"ACGTACGTAAAAAA"; // ref[0..14]
        reads.push(("rshort", 0, short_seq, vec![Cigar::Match(14)], false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rshort" {
                assert_eq!(rec.seq().as_bytes(), short_seq, "guarded read must be untouched");
            }
        }
        assert!(stats.reads_skipped_span_guard >= 1, "the skip must be counted");
        assert_eq!(stats.reads_reassigned_ref_to_ins, 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn read_carrying_an_oversized_indel_is_not_assigned_nearby() {
        // A 6bp insertion exceeds max_len 5. The read must be left alone -- and in
        // particular must not be treated as a clean REF observation and "repaired"
        // into the site's kept allele.
        let dir = std::env::temp_dir().join("himito_denoise_indel_large");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let big_seq: &[u8] = b"ACGTACGTAAAAAAAACGTGGGGGGACGTACGTAC"; // 35bp
        let big_cigar = vec![Cigar::Match(19), Cigar::Ins(6), Cigar::Match(10)];
        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, UNIQ_INS_SEQ, uniq_ins_cigar(), i % 2 == 0));
        }
        reads.push(("rbig", 0, big_seq, big_cigar, false));
        write_cigar_bam(&bam, "chrM", 29, &reads);

        let refs = indel_refs();
        let iopts = indels_on();
        let (corr, models, mut stats) =
            compute_corrections(&bam, &refs, 2, 0.01, SB_P, HOM_VAF, &iopts).unwrap();
        apply_corrections(&bam, &outb, &corr, &refs, &models, &iopts, &mut stats).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rbig" {
                assert_eq!(rec.seq().as_bytes(), big_seq, "oversized indel must be untouched");
            }
        }
        std::fs::remove_dir_all(&dir).ok();
    }
```

Verify the test module already imports `rust_htslib::bam::record::{Cigar, CigarString}` (it imports `Cigar` and `CigarString` for `write_test_bam`); add `Cigar` if only `CigarString` is present. Every existing call to `compute_corrections` in the test module gains a trailing `&indel::IndelOpts::default()` and destructures three values; every existing call to `apply_corrections` gains `&models`, `&iopts`, and `&mut stats`.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test denoise::tests:: 2>&1 | tail -25`
Expected: FAIL to compile — `this function takes 6 arguments but 7 arguments were supplied`, and `no field 'indel_sites_examined' on type 'DenoiseStats'`.

- [ ] **Step 3: Extend the stats struct**

In `src/denoise.rs`, extend `DenoiseStats` (after `substitution_matrix`):

```rust
    // --- indel denoising (all zero unless --indels is enabled) ---
    pub indel_sites_examined: u64,
    pub indel_sites_kept: u64,
    pub indel_sites_strand_biased: u64,
    pub indel_events_examined: u64,
    pub reads_reassigned_ref_to_ins: u64,
    pub reads_reassigned_ref_to_del: u64,
    pub reads_reassigned_indel_to_ref: u64,
    pub reads_reassigned_indel_to_indel: u64,
    pub indel_bases_inserted: u64,
    pub indel_bases_removed: u64,
    /// Index = repeat-context length L; audits how much correction fired in long runs.
    pub reassignments_by_hp_length: Vec<u64>,
    pub reads_skipped_span_guard: u64,
    pub reads_skipped_unsupported_transition: u64,
```

- [ ] **Step 4: Implement pass 1 (site accumulation and decision)**

In `src/denoise.rs`, add the model type alias beside `Corrections`:

```rust
/// Decided indel sites, per contig, SORTED ASCENDING BY POSITION so the write-back
/// pass can binary-search the sites a read spans.
pub type IndelModels = HashMap<Vec<u8>, Vec<(u32, indel::IndelSiteModel)>>;
```

Replace `compute_corrections` in full:

```rust
fn compute_corrections(
    bam_path: &Path,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    min_strand: u32,
    vaf: f64,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
    iopts: &indel::IndelOpts,
) -> Result<(Corrections, IndelModels, DenoiseStats)> {
    let mut reader = Reader::from_path(bam_path)
        .with_context(|| format!("cannot open BAM {bam_path:?}"))?;
    // Own the tid->contig-name map so we can borrow the reader mutably for pileup.
    let header = reader.header().to_owned();

    let mut corrections: Corrections = HashMap::new();
    let mut stats = DenoiseStats::default();

    // ALT event counts per left-normalized site. Normalization can move an event into
    // a column the loop has already passed, so sites are decided after it closes.
    let mut sites: HashMap<(Vec<u8>, u32), indel::IndelSite> = HashMap::new();
    // Per-column (depth, fwd, rev). This -- NOT a tally of observations -- is the VAF
    // denominator: a read carrying an event reports Indel::None at the normalized
    // column as well as its event at the anchor column, so counting observations
    // would count that read twice and halve every real indel's apparent frequency.
    let mut depths: HashMap<(Vec<u8>, u32), (u32, u32, u32)> = HashMap::new();

    let mut plp = reader.pileup();
    plp.set_max_depth(1_000_000);
    for p in plp {
        let col = p.context("pileup error")?;
        let name = header.tid2name(col.tid());
        let refseq = match refs.get(name) {
            Some(s) => s,
            None => continue,
        };
        let refpos = col.pos() as usize;
        let ref_idx = refseq.get(refpos).and_then(|&b| base_to_idx(b));

        let mut counts = ColumnCounts::new(ref_idx);
        // Remember each usable observation to correct after the fit.
        let mut obs: Vec<(Vec<u8>, usize, usize, u8)> = Vec::new(); // (qname, qpos, allele, qual)

        // Indel collection MUST happen here, before the `non_ref() == 0` skip below:
        // that gate is about non-reference SUBSTITUTIONS, and a column with a
        // perfectly clean base composition can still carry indel events. Collecting
        // after it would silently lose every indel at a substitution-clean column.
        let mut col_depth = 0u32;
        let mut col_fwd = 0u32;
        let mut col_rev = 0u32;
        let mut col_events: Vec<(indel::Allele, bool)> = Vec::new();

        for a in col.alignments() {
            let rec = a.record();
            if rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped() {
                continue;
            }

            if iopts.enabled && !a.is_refskip() {
                if let Some(q) = a.qpos() {
                    col_depth += 1;
                    if rec.is_reverse() { col_rev += 1 } else { col_fwd += 1 }
                    match a.indel() {
                        bam::pileup::Indel::Ins(n) if n < iopts.max_len => {
                            let s: Vec<u8> =
                                (0..n as usize).map(|i| rec.seq()[q + 1 + i]).collect();
                            col_events.push((indel::Allele::Ins(s), rec.is_reverse()));
                            stats.indel_events_examined += 1;
                        }
                        bam::pileup::Indel::Del(n) if n < iopts.max_len => {
                            col_events.push((indel::Allele::Del(n), rec.is_reverse()));
                            stats.indel_events_examined += 1;
                        }
                        // Indel::None contributes depth only; an indel at or above
                        // max_len contributes depth and is never corrected.
                        _ => {}
                    }
                }
            }

            if a.is_del() || a.is_refskip() {
                continue;
            }
            let qpos = match a.qpos() {
                Some(x) => x,
                None => continue,
            };
            let base = rec.seq()[qpos];
            let allele = match base_to_idx(base) {
                Some(i) => i,
                None => continue,
            };
            let qual = rec.qual()[qpos];
            counts.add(allele, qual, rec.is_reverse());
            obs.push((rec.qname().to_vec(), qpos, allele, qual));
        }

        if iopts.enabled {
            depths.insert((name.to_vec(), refpos as u32), (col_depth, col_fwd, col_rev));
            for (allele, reverse) in col_events {
                let (norm_pos, norm_allele) =
                    indel::normalize_left(refseq, refpos as u32, &allele);
                sites
                    .entry((name.to_vec(), norm_pos))
                    .or_default()
                    .add(norm_allele, reverse);
            }
        }

        stats.columns_processed += 1;
        stats.bases_examined += obs.len() as u64;

        if counts.non_ref() == 0 {
            stats.columns_skipped += 1;
            continue;
        }

        let model = fit_site(&counts, min_strand, vaf, strand_bias_p, homoplasmic_vaf);
        stats.alt_sites_preserved += model.kept_alt_count(ref_idx) as u64;
        stats.alt_alleles_strand_biased += model.strand_biased_alts as u64;

        let mut changed = false;
        for (qname, _qpos, allele, qual) in obs {
            let corr = map_allele(allele, qual, &model);
            if corr != allele {
                corrections
                    .entry(qname)
                    .or_default()
                    .subs
                    .push((refpos as u32, idx_to_base(corr)));
                stats.bases_modified += 1;
                stats.substitution_matrix[allele][corr] += 1;
                changed = true;
            }
        }
        if changed {
            stats.columns_corrected += 1;
        }
    }

    // ---- deferred indel site decisions ----
    let mut models: IndelModels = HashMap::new();
    if iopts.enabled {
        for ((contig, norm_pos), mut site) in sites.into_iter() {
            let refseq = match refs.get(&contig) {
                Some(s) => s,
                None => continue,
            };
            // Fill the denominator from the column totals: each spanning read once.
            let (depth, fwd, rev) = depths
                .get(&(contig.clone(), norm_pos))
                .copied()
                .unwrap_or((0, 0, 0));
            if depth == 0 {
                continue; // nothing spans this site; nothing to decide
            }
            site.depth = depth;
            site.col_fwd = fwd;
            site.col_rev = rev;

            stats.indel_sites_examined += 1;
            let model = indel::fit_indel_site(
                &site, refseq, norm_pos, min_strand, strand_bias_p, homoplasmic_vaf, iopts,
            );
            stats.indel_sites_strand_biased += model.strand_biased as u64;
            // `kept` includes REF; count only surviving non-reference alleles.
            stats.indel_sites_kept +=
                model.kept.iter().filter(|(a, _)| *a != indel::Allele::Ref).count() as u64;
            models.entry(contig).or_default().push((norm_pos, model));
        }
        // The write-back pass binary-searches each contig's list.
        for list in models.values_mut() {
            list.sort_by_key(|(p, _)| *p);
        }
        for edits in corrections.values_mut() {
            edits.subs.sort_by_key(|&(p, _)| p);
        }
    }

    Ok((corrections, models, stats))
}
```

- [ ] **Step 5: Implement pass 2 (per-read assignment inside the write-back walk)**

Replace `apply_corrections` in `src/denoise.rs`:

```rust
fn apply_corrections(
    in_bam: &Path,
    out_bam: &Path,
    corrections: &Corrections,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    models: &IndelModels,
    iopts: &indel::IndelOpts,
    stats: &mut DenoiseStats,
) -> Result<(u64, u64)> {
    let mut reader = Reader::from_path(in_bam)
        .with_context(|| format!("cannot open BAM {in_bam:?}"))?;
    let header_view = reader.header().to_owned();
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
        .with_context(|| format!("cannot create BAM {out_bam:?}"))?;

    let mut reads_processed = 0u64;
    let mut reads_modified = 0u64;
    let empty: Vec<u8> = Vec::new();

    for result in reader.records() {
        let mut rec = result.context("read BAM record")?;
        // Only primary mapped reads are candidates for correction.
        let correctable =
            !(rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped());
        if correctable {
            reads_processed += 1;
            let contig = header_view.tid2name(rec.tid() as u32).to_vec();
            let refseq = refs.get(&contig).unwrap_or(&empty);
            let mut edits = corrections.get(rec.qname()).cloned().unwrap_or_default();

            if iopts.enabled {
                let cigar = rec.cigar().take();
                let seq = rec.seq().as_bytes();
                let rstart = rec.pos();
                let rend = cigar.end_pos();
                let evs = indel::read_events(rstart, &cigar, &seq, refseq, iopts.max_len);

                // What this read carries at each normalized site, and where it has an
                // indel too large to represent.
                let mut carried: HashMap<u32, indel::Allele> = HashMap::new();
                let mut blocked: Vec<i64> = Vec::new();
                for e in &evs {
                    match &e.allele {
                        Some(a) => { carried.insert(e.norm_pos, a.clone()); }
                        None => blocked.push(e.norm_pos as i64),
                    }
                }

                let flank = iopts.flank as i64;
                // A gained deletion consumes up to max_len reference bases to the
                // right, so the right-hand requirement is larger than the left.
                let reach = flank + iopts.max_len as i64;

                if let Some(list) = models.get(&contig) {
                    let lo = list.partition_point(|(p, _)| (*p as i64) < rstart);
                    for (pos, model) in &list[lo..] {
                        let s = *pos as i64;
                        if s > rend {
                            break;
                        }
                        let observed =
                            carried.get(pos).cloned().unwrap_or(indel::Allele::Ref);
                        let target = match indel::assign_allele(&observed, model, iopts) {
                            indel::Assignment::Keep => continue,
                            indel::Assignment::Refused => {
                                stats.reads_skipped_unsupported_transition += 1;
                                continue;
                            }
                            indel::Assignment::Move(t) => t,
                        };
                        // The rewrite walk needs plain match context on both sides.
                        if s - flank < rstart || s + reach > rend {
                            stats.reads_skipped_span_guard += 1;
                            continue;
                        }
                        // An oversized indel nearby means this read's alignment in
                        // this window cannot be represented by the site's alleles.
                        if blocked.iter().any(|b| (b - s).abs() <= reach) {
                            stats.reads_skipped_span_guard += 1;
                            continue;
                        }
                        match (&observed, &target) {
                            (indel::Allele::Ref, indel::Allele::Ins(s2)) => {
                                stats.reads_reassigned_ref_to_ins += 1;
                                stats.indel_bases_inserted += s2.len() as u64;
                            }
                            (indel::Allele::Ref, indel::Allele::Del(n)) => {
                                stats.reads_reassigned_ref_to_del += 1;
                                stats.indel_bases_removed += *n as u64;
                            }
                            (_, indel::Allele::Ref) => {
                                stats.reads_reassigned_indel_to_ref += 1;
                            }
                            _ => stats.reads_reassigned_indel_to_indel += 1,
                        }
                        let l = model.context_len as usize;
                        if stats.reassignments_by_hp_length.len() <= l {
                            stats.reassignments_by_hp_length.resize(l + 1, 0);
                        }
                        stats.reassignments_by_hp_length[l] += 1;

                        edits.indels.push(rewrite::IndelEdit {
                            ref_pos: *pos,
                            from: observed,
                            to: target,
                        });
                    }
                }
                edits.indels.sort_by_key(|e| e.ref_pos);
            }

            if !edits.is_empty() {
                // Extract owned parts, rewrite, rebuild (aux preserved by set()).
                let qname = rec.qname().to_vec();
                let cigar = rec.cigar().take();
                let seq = rec.seq().as_bytes();
                let qual = rec.qual().to_vec();
                let out = rewrite_read(rec.pos(), &cigar, &seq, &qual, &edits, refseq);
                rec.set(&qname, Some(&out.cigar), &out.seq, &out.qual);
                if out.len_changed {
                    strip_stale_tags(&mut rec);
                }
                reads_modified += 1;
            }
        }
        writer.write(&rec).context("write BAM record")?;
    }
    Ok((reads_processed, reads_modified))
}
```

Update `start` to thread the new values through (the `indels` parameter arrives in Task 9; until then pass `&indel::IndelOpts::default()`):

```rust
    let (corrections, models, mut stats) = compute_corrections(
        input, &refs, min_strand, vaf, strand_bias_p, homoplasmic_vaf, indels,
    )?;
    let (reads_processed, reads_modified) =
        apply_corrections(input, output, &corrections, &refs, &models, indels, &mut stats)?;
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `cargo test denoise:: 2>&1 | tail -30`
Expected: all tests PASS, including the seven new indel tests and every pre-existing test.

- [ ] **Step 7: Build and commit**

```bash
cargo build --release
git add src/denoise.rs
git commit -m "feat(denoise): accumulate indel sites in pileup, assign per read in write-back"
```

---

### Task 9: CLI surface and `IndelOpts` plumbing

Exposes the feature. `start` takes one struct rather than nine more scalars, and `QuickStart` passes the disabled default so no existing result shifts.

**Files:**
- Modify: `src/denoise.rs` (`start` signature at line ~398, and the log line), `src/main.rs:195-233` (the `Denoise` subcommand), `src/main.rs:629-647` (the `Denoise` handler), `src/main.rs:553-556` (the `QuickStart` denoise call)

**Interfaces:**
- Consumes: `indel::IndelOpts`.
- Produces: `pub fn start(input: &PathBuf, output: &PathBuf, reference: &PathBuf, data_type: &str, vaf: f64, min_strand: u32, strand_bias_p: f64, homoplasmic_vaf: f64, indels: &indel::IndelOpts, stats_out: Option<&PathBuf>) -> Result<()>`.

- [ ] **Step 1: Write the failing test**

Append inside `mod tests` in `src/denoise.rs`:

```rust
    #[test]
    fn start_with_indels_disabled_is_byte_identical_to_v1() {
        // The opt-in guarantee: with the default IndelOpts, denoised output must not
        // differ by a single byte from the substitution-only behavior.
        let dir = std::env::temp_dir().join("himito_denoise_optin");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let out_a = dir.join("a.bam");
        let out_b = dir.join("b.bam");
        let reff = dir.join("ref.fa");
        std::fs::write(&reff, ">chrM\nACGTACGTAAAAAAAACGTACGTACGTAC\n").unwrap();

        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, UNIQ_INS_SEQ, uniq_ins_cigar(), false));
        write_cigar_bam(&inb, "chrM", 29, &reads);

        let off = indel::IndelOpts::default();
        start(&inb, &out_a, &reff, "ont-r10", 0.01, 2, SB_P, HOM_VAF, &off, None).unwrap();
        start(&inb, &out_b, &reff, "ont-r10", 0.01, 2, SB_P, HOM_VAF, &off, None).unwrap();
        assert_eq!(std::fs::read(&out_a).unwrap(), std::fs::read(&out_b).unwrap());

        // And the indel-bearing read keeps its insertion.
        let mut r = Reader::from_path(&out_a).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rins" {
                assert_eq!(rec.seq().as_bytes(), UNIQ_INS_SEQ);
                seen = true;
            }
        }
        assert!(seen);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn start_with_indels_enabled_removes_the_spurious_insertion() {
        let dir = std::env::temp_dir().join("himito_denoise_optin_on");
        std::fs::create_dir_all(&dir).unwrap();
        let inb = dir.join("in.bam");
        let outb = dir.join("out.bam");
        let reff = dir.join("ref.fa");
        let statsf = dir.join("stats.json");
        std::fs::write(&reff, ">chrM\nACGTACGTAAAAAAAACGTACGTACGTAC\n").unwrap();

        let mut reads: Vec<(&str, i64, &[u8], Vec<Cigar>, bool)> = vec![];
        for (i, n) in NAMES.iter().enumerate() {
            reads.push((n, 0, CLEAN_SEQ, clean_cigar(), i % 2 == 0));
        }
        reads.push(("rins", 0, UNIQ_INS_SEQ, uniq_ins_cigar(), false));
        write_cigar_bam(&inb, "chrM", 29, &reads);

        let mut on = indel::IndelOpts::default();
        on.enabled = true;
        start(&inb, &outb, &reff, "ont-r10", 0.01, 2, SB_P, HOM_VAF, &on, Some(&statsf)).unwrap();

        let mut r = Reader::from_path(&outb).unwrap();
        let mut seen = false;
        for rec in r.records() {
            let rec = rec.unwrap();
            if rec.qname() == b"rins" {
                assert_eq!(rec.seq().as_bytes(), CLEAN_SEQ, "insertion must be removed");
                assert_eq!(rec.cigar().to_string(), "29M", "CIGAR must be rebuilt");
                seen = true;
            }
        }
        assert!(seen);
        let json = std::fs::read_to_string(&statsf).unwrap();
        assert!(json.contains("\"reads_reassigned_indel_to_ref\": 1"));
        std::fs::remove_dir_all(&dir).ok();
    }
```

Update the two existing `start(...)` calls in the test module (`start_pacbio_is_byte_identical_passthrough` and `start_ont_denoises_and_writes_stats_json`) to pass `&indel::IndelOpts::default()` before the final `stats_out` argument.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test denoise::tests::start_with 2>&1 | tail -15`
Expected: FAIL to compile — `this function takes 9 arguments but 10 arguments were supplied`.

- [ ] **Step 3: Add the parameter and the CLI flags**

In `src/denoise.rs`, add the `indels` parameter to `start` (before `stats_out`):

```rust
pub fn start(
    input: &PathBuf,
    output: &PathBuf,
    reference: &PathBuf,
    data_type: &str,
    vaf: f64,
    min_strand: u32,
    strand_bias_p: f64,
    homoplasmic_vaf: f64,
    indels: &indel::IndelOpts,
    stats_out: Option<&PathBuf>,
) -> Result<()> {
```

and replace the temporary `&indel::IndelOpts::default()` from Task 7 with `indels`. Extend the summary log line so the indel work is visible without opening the JSON:

```rust
    if indels.enabled {
        info!(
            "[denoise] indel sites {} examined / {} alt alleles kept / {} strand-biased; \
             reads reassigned ref->ins {} ref->del {} indel->ref {} indel->indel {}; \
             skipped span-guard {} unsupported-transition {}",
            stats.indel_sites_examined, stats.indel_sites_kept, stats.indel_sites_strand_biased,
            stats.reads_reassigned_ref_to_ins, stats.reads_reassigned_ref_to_del,
            stats.reads_reassigned_indel_to_ref, stats.reads_reassigned_indel_to_indel,
            stats.reads_skipped_span_guard, stats.reads_skipped_unsupported_transition,
        );
    }
```

In `src/main.rs`, add the flags to the `Denoise` subcommand, immediately before the existing `stats` field (line 230):

```rust
        /// enable small-indel (<5bp) correction in addition to substitutions
        #[clap(long, value_parser, default_value_t = false)]
        indels: bool,

        /// exclusive upper bound on correctable indel length (5 = lengths 1..=4)
        #[clap(long, value_parser, default_value_t = 5)]
        indel_max_len: u32,

        /// absolute minimal VAF floor for an indel allele to be kept
        #[clap(long, value_parser, default_value_t = 0.05)]
        indel_vaf: f64,

        /// per-junction indel error probability in unique sequence
        #[clap(long, value_parser, default_value_t = 0.01)]
        indel_err0: f64,

        /// multiplicative growth of the indel error rate per extra repeat copy
        #[clap(long, value_parser, default_value_t = 1.5)]
        indel_err_scale: f64,

        /// ceiling on the context-scaled indel error rate
        #[clap(long, value_parser, default_value_t = 0.4)]
        indel_err_cap: f64,

        /// candidacy floor is at least this multiple of the local indel error rate
        #[clap(long, value_parser, default_value_t = 3.0)]
        indel_floor_mult: f64,

        /// decay of error mass per unit of net-length distance between alleles
        #[clap(long, value_parser, default_value_t = 0.3)]
        indel_delta: f64,

        /// reads must extend this many bases past a site on both sides to be reassigned
        #[clap(long, value_parser, default_value_t = 5)]
        indel_flank: usize,

        /// an allele carried by at least this fraction of reads at a site is never
        /// corrected away, whatever candidacy says (protects real heteroplasmy in
        /// deep repeat contexts, where the candidacy floor can exceed 1.0)
        #[clap(long, value_parser, default_value_t = 0.2)]
        indel_protect_vaf: f64,
```

Replace the `Denoise` handler (lines 629-647):

```rust
        Commands::Denoise {
            input,
            output,
            reference,
            data_type,
            vaf,
            min_strand,
            strand_bias_p,
            homoplasmic_vaf,
            indels,
            indel_max_len,
            indel_vaf,
            indel_err0,
            indel_err_scale,
            indel_err_cap,
            indel_floor_mult,
            indel_delta,
            indel_flank,
            indel_protect_vaf,
            stats,
        } => {
            let iopts = denoise::indel::IndelOpts {
                enabled: indels,
                max_len: indel_max_len,
                vaf: indel_vaf,
                err0: indel_err0,
                err_scale: indel_err_scale,
                err_cap: indel_err_cap,
                floor_mult: indel_floor_mult,
                delta: indel_delta,
                flank: indel_flank,
                protect_vaf: indel_protect_vaf,
            };
            if let Err(e) = denoise::start(
                &input, &output, &reference, &data_type, vaf, min_strand,
                strand_bias_p, homoplasmic_vaf, &iopts, stats.as_ref(),
            ) {
                eprintln!("Error running denoise: {e:#}");
                std::process::exit(1);
            }
        }
```

In the `QuickStart` denoise call (lines 553-556), pass the disabled default so QuickStart results do not shift:

```rust
                if let Err(e) = denoise::start(
                    &mt_output, &denoised, &reference_path, &data_type,
                    vaf_threshold as f64, 2, 0.01, 0.7,
                    // Indel correction stays off in QuickStart until the lineage_sim
                    // validation gate in the design doc passes.
                    &denoise::indel::IndelOpts::default(),
                    None,
                ) {
```

- [ ] **Step 4: Run the full suite and check the CLI**

Run: `cargo test denoise:: 2>&1 | tail -30`
Expected: all tests PASS.

Run: `cargo build --release && ./target/release/Himito denoise --help 2>&1 | grep -c indel`
Expected: at least `10` (one line per new flag).

- [ ] **Step 5: Commit**

```bash
git add src/denoise.rs src/main.rs
git commit -m "feat(denoise): --indels CLI surface; QuickStart keeps indel correction off"
```

---

### Task 10: Validation harness wiring and user documentation

Makes the design's ship gate runnable. The harness gains an opt-in passthrough so `--indels` off vs on can be compared on identical simulated data.

**Files:**
- Modify: `lineage_eval/lineage_sim/run_himito.sh:95-104`
- Modify: `docs/User_guide.md`

**Interfaces:**
- Consumes: the `--indels` flag from Task 9.
- Produces: `DENOISE_INDELS` environment variable (default `0`) honored by `run_himito.sh`.

- [ ] **Step 1: Add the opt-in passthrough to the harness**

In `lineage_eval/lineage_sim/run_himito.sh`, replace the denoise block (lines 95-104):

```bash
# Indel denoising is OFF by default, matching quick-start. Set DENOISE_INDELS=1 to
# run the off-vs-on comparison the design's ship gate requires:
#   GFA edge/anchor count and matrix column count must drop, with no regression in
#   the lineage tree score or small-indel precision/recall.
DENOISE_INDELS="${DENOISE_INDELS:-0}"

BUILD_BAM="$BAM"
CALL_DATATYPE="$DTYPE"
if [[ "$DTYPE" == ont-* ]]; then
  DENOISED="$HDIR/aln.denoised.bam"
  DENOISE_ARGS=(--vaf "$VAF" --min-strand 2 --stats "$HDIR/denoise_stats.json")
  if [[ "$DENOISE_INDELS" == "1" ]]; then
    DENOISE_ARGS+=(--indels)
  fi
  "$HIMITO" denoise -i "$BAM" -o "$DENOISED" -r "$REF" -d "$DTYPE" "${DENOISE_ARGS[@]}"
  samtools index "$DENOISED"
  BUILD_BAM="$DENOISED"
  CALL_DATATYPE="ont-denoised"
fi
```

- [ ] **Step 2: Verify the harness still parses and defaults to off**

Run: `bash -n lineage_eval/lineage_sim/run_himito.sh && echo SYNTAX_OK`
Expected: `SYNTAX_OK`

Run: `grep -n 'DENOISE_INDELS' lineage_eval/lineage_sim/run_himito.sh`
Expected: three lines — the comment-adjacent default assignment, the `if` test, and nothing that enables it unconditionally.

- [ ] **Step 3: Document the feature and its one sharp edge**

Append a section to `docs/User_guide.md`:

```markdown
### Small-indel denoising (ONT, opt-in)

`Himito denoise --indels` extends substitution correction to indels shorter than
`--indel-max-len` (default 5, exclusive — so lengths 1–4). Reads are reassigned at
left-normalized indel sites in both directions: a read can lose a spurious indel,
gain a consensus indel it dropped, or snap from one indel allele to another of the
same kind.

Correction aggressiveness is governed by a repeat-context-scaled error rate: an indel
inside a long homopolymer must clear a much higher frequency bar than the same indel
in unique sequence, because that is where ONT's indel errors concentrate. A read only
leaves its own allele when that allele's frequency falls below `eps / (1 - eps)` of a
competing one.

**Not enabled by QuickStart.** Indel correction is off unless you pass `--indels`
explicitly.

**Reads whose length changed lose their `MM`/`ML`/`MN` and `NM`/`MD`/`cs` tags.**
Base-modification tags are indexed by the read's own base positions and cannot survive
a length change. This is safe for the standard pipeline because methylation
aggregation reads the original mt BAM, not the denoised one — but do not use an
indel-denoised BAM as a methylation input.

**Do not re-run the denoiser on its own output.** Allele frequencies are computed from
the input BAM and decisions are applied once. Feeding denoised reads back in would let
each pass shift the frequencies the next pass sees, compounding correction.

Statistics for the indel pass (sites examined/kept/strand-biased, per-transition read
counts, and `reassignments_by_hp_length`, which shows how much correction fired inside
long homopolymers) are written to the `--stats` JSON and summarized in the log.
```

- [ ] **Step 4: Verify the docs render and the full suite is green**

Run: `cargo test 2>&1 | tail -15`
Expected: `test result: ok.` with 0 failures.

Run: `grep -c 'indel' docs/User_guide.md`
Expected: a non-zero count (at least 8).

- [ ] **Step 5: Commit**

```bash
git add lineage_eval/lineage_sim/run_himito.sh docs/User_guide.md
git commit -m "feat(denoise): DENOISE_INDELS harness switch; document indel denoising"
```

---

## Post-Implementation: Run the Validation Gate

The plan delivers the feature **disabled**. Flipping the default is a separate decision that requires evidence. Run:

```bash
cd lineage_eval/lineage_sim
DENOISE_INDELS=0 ./run_himito.sh --profile ont-r10 --outdir /tmp/indel_off
DENOISE_INDELS=1 ./run_himito.sh --profile ont-r10 --outdir /tmp/indel_on
```

Compare, per the design's §9:

1. GFA edge and anchor counts (`grep -c '^L' sim.gfa`, `grep -c '^S' sim.gfa`) — should **drop**.
2. Matrix column count (`head -1 sim.vcf.matrix.csv | tr ',' '\n' | wc -l`) vs. the true variant count — should **drop toward truth**.
3. `score_lineage.py` tree score — must **not regress**.
4. Small-indel precision/recall vs. truth VCF — must **not regress**.

Only if (1) and (2) improve with no regression in (3) or (4) should `IndelOpts::default().enabled` be flipped to `true` and the QuickStart call updated. If the numbers disappoint, tune `indel_err0` / `indel_err_scale` / `indel_floor_mult` first — those set how aggressive correction is — before concluding the approach is wrong.
