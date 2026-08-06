# ONT BAM Denoiser (`Himito denoise`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an ONT-only, SNV-only BAM denoiser that corrects per-read substitution errors with a baldur-style per-site ML error model, so reads consolidate and the existing `build`/`call`/`lineage` pipeline recovers low-frequency variants.

**Architecture:** New `src/denoise.rs` module. Two passes over the BAM: (1) `pileup` → per-column ML fit → collect per-read base corrections keyed by qname; (2) rewrite each read's `SEQ` (aux tags preserved by `Record::set`). A `denoise` subcommand exposes it; `QuickStart` auto-applies it for `ont-*` (routing `build`/`call` to the denoised BAM while methylation keeps the original). PacBio is a passthrough.

**Tech Stack:** Rust, `rust-htslib` 0.49.0 (`bam::Reader`/`Writer`/pileup), `bio::io::fasta`, `serde`/`serde_json` (stats), `rayon` (already used).

## Global Constraints

- ONT-only behavior change; **PacBio path unchanged** (`-d pacbio` → byte-identical passthrough).
- **SNV-only**: never touch indel/soft-clip/ref-skip positions; base identity only, never `CIGAR`/length/positions.
- Downstream `build`/`call`/`lineage` and methylation are **unchanged**; methylation always reads the **original** (un-denoised) BAM.
- Defaults (verbatim): `--vaf 0.01`, `--min-strand 2`. EM: `max_iter = 100`, `tol = 1e-6`. Pileup `set_max_depth(1_000_000)` (avoid htslib's 8000 truncation).
- `rust-htslib` `Record::set(qname, Some(&cigar), &seq, &qual)` preserves aux data (verified: record.rs:340 "aux data is left unchanged"); use it for the SEQ rewrite. No manual aux re-attach needed.
- `Seq` indexing `rec.seq()[qpos]` is O(1) decoded base access — use it; never call `.as_bytes()` per pileup observation.
- Base index convention: `A=0, C=1, G=2, T=3` (uppercase); any other base (`N`, etc.) is "not a base" → skipped.

---

### Task 1: Numeric core — error model, per-column counts, EM fit, MAP correction

Pure functions, no I/O. This is the heart of the denoiser and is fully unit-tested.

**Files:**
- Create: `src/denoise.rs`
- Modify: `src/main.rs` (add `mod denoise;` near the other `mod` lines, ~line 16)

**Interfaces:**
- Produces (used by Tasks 2 & 4):
  - `pub(crate) fn phred_to_e(q: u8) -> f64`
  - `fn base_to_idx(b: u8) -> Option<usize>`, `fn idx_to_base(i: usize) -> u8`
  - `const NQ: usize = 94;`
  - `struct ColumnCounts { n: [[u32; NQ]; 4], strand: [[u32; 2]; 4], ref_idx: Option<usize> }` with `fn new(ref_idx: Option<usize>) -> Self`, `fn add(&mut self, allele: usize, q: u8, reverse: bool)`, `fn non_ref(&self) -> u32`
  - `struct SiteModel { freq: [f64; 4], kept: [bool; 4] }` with `fn kept_alt_count(&self, ref_idx: Option<usize>) -> usize`
  - `fn fit_site(c: &ColumnCounts, min_strand: u32, vaf: f64) -> SiteModel`
  - `fn map_allele(observed: usize, q: u8, m: &SiteModel) -> usize`

- [ ] **Step 1: Add the module declaration**

In `src/main.rs`, next to the existing `mod scite;` (~line 16), add:

```rust
mod denoise;
```

- [ ] **Step 2: Write the failing tests**

Create `src/denoise.rs` with the test module first:

```rust
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use anyhow::{Context, Result};

// ... (implementation added in later steps) ...

#[cfg(test)]
mod tests {
    use super::*;

    fn q(v: u8) -> u8 { v } // readability helper for quality values

    #[test]
    fn phred_to_e_matches_definition() {
        assert!((phred_to_e(10) - 0.1).abs() < 1e-12);
        assert!((phred_to_e(20) - 0.01).abs() < 1e-12);
        assert!((phred_to_e(0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn base_index_round_trips_and_rejects_n() {
        for (b, i) in [(b'A', 0usize), (b'C', 1), (b'G', 2), (b'T', 3)] {
            assert_eq!(base_to_idx(b), Some(i));
            assert_eq!(idx_to_base(i), b);
        }
        assert_eq!(base_to_idx(b'N'), None);
        assert_eq!(base_to_idx(b'a'), Some(0)); // lowercase accepted
    }

    // Build a column: `obs` is (allele_idx, qual, reverse) repeated.
    fn column(ref_idx: usize, obs: &[(usize, u8, bool)]) -> ColumnCounts {
        let mut c = ColumnCounts::new(Some(ref_idx));
        for &(a, ql, rev) in obs { c.add(a, ql, rev); }
        c
    }

    #[test]
    fn pure_ref_column_with_sparse_errors_keeps_only_ref() {
        // ref = A(0); 30 clean A plus 2 stray G/C at low quality, not strand-balanced.
        let mut obs = vec![];
        for i in 0..30 { obs.push((0usize, q(30), i % 2 == 0)); }
        obs.push((2, q(8), false)); // a single G error (fwd only)
        obs.push((1, q(8), true));  // a single C error (rev only)
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01);
        assert!(m.kept[0]); // ref kept
        assert_eq!(m.kept_alt_count(Some(0)), 0); // no alt survives
        // the two errors correct to ref A
        assert_eq!(map_allele(2, q(8), &m), 0);
        assert_eq!(map_allele(1, q(8), &m), 0);
    }

    #[test]
    fn strand_balanced_heteroplasmy_is_preserved() {
        // ref A(0) ~60%, alt G(2) ~40%, both strand-balanced, decent quality.
        let mut obs = vec![];
        for i in 0..60 { obs.push((0usize, q(30), i % 2 == 0)); }
        for i in 0..40 { obs.push((2usize, q(30), i % 2 == 0)); }
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01);
        assert!(m.kept[0] && m.kept[2]);
        assert_eq!(m.kept_alt_count(Some(0)), 1);
        // confident observations of each allele are preserved
        assert_eq!(map_allele(0, q(30), &m), 0);
        assert_eq!(map_allele(2, q(30), &m), 2);
        // a third-base error (T) collapses to the higher-frequency allele (ref A)
        assert_eq!(map_allele(3, q(8), &m), 0);
    }

    #[test]
    fn strand_skewed_candidate_is_corrected_away() {
        // ref A(0); an alt G(2) at 10% but ALL on the forward strand -> fails candidacy.
        let mut obs = vec![];
        for i in 0..90 { obs.push((0usize, q(30), i % 2 == 0)); }
        for _ in 0..10 { obs.push((2usize, q(30), false)); } // fwd only
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.01);
        assert!(!m.kept[2]);
        assert_eq!(map_allele(2, q(30), &m), 0); // corrected to ref
    }

    #[test]
    fn low_vaf_but_balanced_allele_below_threshold_is_dropped() {
        // alt G(2) present twice per strand (passes candidacy) but VAF ~0.02 < 0.05 threshold.
        let mut obs = vec![];
        for i in 0..196 { obs.push((0usize, q(30), i % 2 == 0)); }
        for _ in 0..2 { obs.push((2usize, q(30), false)); }
        for _ in 0..2 { obs.push((2usize, q(30), true)); }
        let c = column(0, &obs);
        let m = fit_site(&c, 2, 0.05); // higher vaf threshold
        assert!(!m.kept[2]);
    }
}
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `cargo test denoise:: 2>&1 | tail -20`
Expected: compile error (functions/types not defined).

- [ ] **Step 4: Write the implementation**

Add above the `#[cfg(test)]` block in `src/denoise.rs`:

```rust
pub(crate) const NQ: usize = 94; // Phred qualities 0..=93

#[inline]
pub(crate) fn phred_to_e(q: u8) -> f64 {
    let q = q.min((NQ - 1) as u8);
    10f64.powf(-(q as f64) / 10.0)
}

#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[inline]
fn idx_to_base(i: usize) -> u8 {
    [b'A', b'C', b'G', b'T'][i]
}

pub(crate) struct ColumnCounts {
    n: [[u32; NQ]; 4],     // n[allele][quality]
    strand: [[u32; 2]; 4], // strand[allele][0=fwd, 1=rev]
    ref_idx: Option<usize>,
}

impl ColumnCounts {
    fn new(ref_idx: Option<usize>) -> Self {
        ColumnCounts { n: [[0; NQ]; 4], strand: [[0; 2]; 4], ref_idx }
    }

    fn add(&mut self, allele: usize, q: u8, reverse: bool) {
        let q = q.min((NQ - 1) as u8) as usize;
        self.n[allele][q] += 1;
        self.strand[allele][if reverse { 1 } else { 0 }] += 1;
    }

    fn non_ref(&self) -> u32 {
        (0..4)
            .filter(|&a| Some(a) != self.ref_idx)
            .map(|a| self.strand[a][0] + self.strand[a][1])
            .sum()
    }
}

pub(crate) struct SiteModel {
    freq: [f64; 4],
    kept: [bool; 4],
}

impl SiteModel {
    fn kept_alt_count(&self, ref_idx: Option<usize>) -> usize {
        (0..4)
            .filter(|&a| self.kept[a] && Some(a) != ref_idx)
            .count()
    }
}

fn fit_site(c: &ColumnCounts, min_strand: u32, vaf: f64) -> SiteModel {
    // Eligibility: reference always eligible; a non-ref allele needs >= min_strand
    // observations on EACH strand (baldur's strand-balance candidacy).
    let mut eligible = [false; 4];
    for a in 0..4 {
        eligible[a] = Some(a) == c.ref_idx
            || (c.strand[a][0] >= min_strand && c.strand[a][1] >= min_strand);
    }

    // Precompute per-quality error and total observations.
    let e: [f64; NQ] = std::array::from_fn(|q| phred_to_e(q as u8));
    let total_obs: f64 = c.n.iter().flatten().map(|&x| x as f64).sum();

    // Initialize f from raw proportions restricted to eligible alleles.
    let mut f = [0.0f64; 4];
    let mut init_sum = 0.0;
    for a in 0..4 {
        if eligible[a] {
            let t: u32 = c.n[a].iter().sum();
            f[a] = t as f64;
            init_sum += t as f64;
        }
    }
    if init_sum <= 0.0 {
        // No eligible observations: fall back to reference-only.
        if let Some(r) = c.ref_idx {
            f[r] = 1.0;
        }
    } else {
        for a in 0..4 {
            f[a] /= init_sum;
        }
    }

    // EM. Observations of non-eligible alleles (errors) contribute their
    // responsibility to the eligible alleles; frequencies normalize over all obs.
    if total_obs > 0.0 {
        for _iter in 0..100 {
            let mut newf = [0.0f64; 4];
            for a in 0..4 {
                for qi in 0..NQ {
                    let cnt = c.n[a][qi];
                    if cnt == 0 {
                        continue;
                    }
                    let eq = e[qi];
                    // p_k ∝ f_k * P(obs=a | true=k)
                    let mut p = [0.0f64; 4];
                    let mut psum = 0.0;
                    for k in 0..4 {
                        if !eligible[k] || f[k] <= 0.0 {
                            continue;
                        }
                        let like = if k == a { 1.0 - eq } else { eq / 3.0 };
                        p[k] = f[k] * like;
                        psum += p[k];
                    }
                    if psum <= 0.0 {
                        continue;
                    }
                    let cntf = cnt as f64;
                    for k in 0..4 {
                        if p[k] > 0.0 {
                            newf[k] += cntf * p[k] / psum;
                        }
                    }
                }
            }
            for k in 0..4 {
                newf[k] /= total_obs;
            }
            let delta: f64 = (0..4).map(|k| (newf[k] - f[k]).abs()).fold(0.0, f64::max);
            f = newf;
            if delta < 1e-6 {
                break;
            }
        }
    }

    // Keep rule: reference always kept; a non-ref eligible allele kept iff f >= vaf.
    let mut kept = [false; 4];
    for a in 0..4 {
        if Some(a) == c.ref_idx {
            kept[a] = true;
        } else if eligible[a] && f[a] >= vaf {
            kept[a] = true;
        }
    }
    // Renormalize freq over kept alleles (guards map_allele against dropped mass).
    let ksum: f64 = (0..4).filter(|&a| kept[a]).map(|a| f[a]).sum();
    let mut freq = [0.0f64; 4];
    if ksum > 0.0 {
        for a in 0..4 {
            if kept[a] {
                freq[a] = f[a] / ksum;
            }
        }
    } else if let Some(r) = c.ref_idx {
        freq[r] = 1.0;
    }

    SiteModel { freq, kept }
}

fn map_allele(observed: usize, q: u8, m: &SiteModel) -> usize {
    let eq = phred_to_e(q);
    let mut best = observed;
    let mut best_val = -1.0f64;
    for k in 0..4 {
        if !m.kept[k] {
            continue;
        }
        let like = if k == observed { 1.0 - eq } else { eq / 3.0 };
        let val = m.freq[k] * like;
        if val > best_val {
            best_val = val;
            best = k;
        }
    }
    best
}
```

Note: the leading `use` lines (`HashMap`, `Path`/`PathBuf`, `anyhow`) are consumed by later tasks; if the compiler warns about unused imports at this task, keep them — Task 2/4 use them.

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test denoise:: 2>&1 | tail -20`
Expected: all Task-1 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/denoise.rs src/main.rs
git commit -m "feat(denoise): per-site ML error model core (EM fit + MAP correction)"
```

---

### Task 2: Pileup pass — build corrections + statistics from a BAM

**Files:**
- Modify: `src/denoise.rs`
- Test: `src/denoise.rs` (`#[cfg(test)]` module, add a small in-memory BAM helper)

**Interfaces:**
- Consumes: `ColumnCounts`, `fit_site`, `map_allele`, `base_to_idx`, `idx_to_base` (Task 1).
- Produces (used by Tasks 3 & 4):
  - `pub type Corrections = HashMap<Vec<u8>, Vec<(usize, u8)>>;` (qname → list of `(qpos, new_base)`)
  - `#[derive(Default, Debug, serde::Serialize)] pub struct DenoiseStats { ... }` (fields below)
  - `fn compute_corrections(bam_path: &Path, refs: &HashMap<Vec<u8>, Vec<u8>>, min_strand: u32, vaf: f64) -> Result<(Corrections, DenoiseStats)>`

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)] mod tests` in `src/denoise.rs`:

```rust
use rust_htslib::bam::{self, Format, Header, Read, Reader, Record, Writer};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};

// Build a tiny coordinate-sorted BAM. Reads are full-length matches at `pos`.
// Each read: (qname, pos, seq bytes, reverse). Quality fixed at Q30.
fn write_test_bam(path: &Path, contig: &str, contig_len: usize, reads: &[(&str, i64, &[u8], bool)]) {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", contig);
    sq.push_tag(b"LN", &contig_len.to_string());
    header.push_record(&sq);
    let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
    for (name, pos, seq, reverse) in reads {
        let mut rec = Record::new();
        let quals = vec![30u8; seq.len()];
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        rec.set(name.as_bytes(), Some(&cigar), seq, &quals);
        rec.set_tid(0);
        rec.set_pos(*pos);
        rec.set_mapq(60);
        rec.set_flags(if *reverse { 16 } else { 0 });
        w.write(&rec).unwrap();
    }
}

#[test]
fn compute_corrections_flags_a_single_base_error() {
    let dir = std::env::temp_dir().join("himito_denoise_t2");
    std::fs::create_dir_all(&dir).unwrap();
    let bam = dir.join("in.bam");
    // Reference AAAAAA at chrM. Ten reads all "AAAAAA" except one with a G at pos 2.
    let refseq = b"AAAAAA".to_vec();
    let mut reads: Vec<(&str, i64, &[u8], bool)> = vec![];
    let clean = b"AAAAAA";
    let errd = b"AAGAAA"; // error at index 2
    for i in 0..9 { reads.push((["r0","r1","r2","r3","r4","r5","r6","r7","r8"][i], 0, clean, i % 2 == 0)); }
    reads.push(("rerr", 0, errd, false));
    write_test_bam(&bam, "chrM", 6, &reads);

    let mut refs = HashMap::new();
    refs.insert(b"chrM".to_vec(), refseq);
    let (corr, stats) = compute_corrections(&bam, &refs, 2, 0.01).unwrap();

    // The single G (not strand-balanced, VAF 10% but fwd-only) is corrected to A.
    assert_eq!(corr.get(&b"rerr".to_vec()), Some(&vec![(2usize, b'A')]));
    assert_eq!(stats.bases_modified, 1);
    assert_eq!(stats.substitution_matrix[2][0], 1); // G(2) -> A(0)
    assert!(stats.columns_processed >= 1);
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test denoise::tests::compute_corrections_flags_a_single_base_error 2>&1 | tail -20`
Expected: compile error (`compute_corrections`, `DenoiseStats`, `Corrections` undefined).

- [ ] **Step 3: Write the implementation**

Add to `src/denoise.rs` (before the test module). Add imports at the top of the file:

```rust
use rust_htslib::bam::{self, Read, Reader};
```

Then:

```rust
pub type Corrections = HashMap<Vec<u8>, Vec<(usize, u8)>>;

#[derive(Default, Debug, serde::Serialize)]
pub struct DenoiseStats {
    pub reads_processed: u64,
    pub reads_modified: u64,
    pub bases_examined: u64,
    pub bases_modified: u64,
    pub columns_processed: u64,
    pub columns_corrected: u64,
    pub columns_skipped: u64,
    pub alt_sites_preserved: u64,
    pub substitution_matrix: [[u64; 4]; 4], // [from][to]
}

fn compute_corrections(
    bam_path: &Path,
    refs: &HashMap<Vec<u8>, Vec<u8>>,
    min_strand: u32,
    vaf: f64,
) -> Result<(Corrections, DenoiseStats)> {
    let mut reader = Reader::from_path(bam_path)
        .with_context(|| format!("cannot open BAM {bam_path:?}"))?;
    // Own the tid->contig-name map so we can borrow the reader mutably for pileup.
    let header = reader.header().to_owned();

    let mut corrections: Corrections = HashMap::new();
    let mut stats = DenoiseStats::default();

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

        for a in col.alignments() {
            let rec = a.record();
            if rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped() {
                continue;
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

        stats.columns_processed += 1;
        stats.bases_examined += obs.len() as u64;

        if counts.non_ref() == 0 {
            stats.columns_skipped += 1;
            continue;
        }

        let model = fit_site(&counts, min_strand, vaf);
        stats.alt_sites_preserved += model.kept_alt_count(ref_idx) as u64;

        let mut changed = false;
        for (qname, qpos, allele, qual) in obs {
            let corr = map_allele(allele, qual, &model);
            if corr != allele {
                corrections
                    .entry(qname)
                    .or_default()
                    .push((qpos, idx_to_base(corr)));
                stats.bases_modified += 1;
                stats.substitution_matrix[allele][corr] += 1;
                changed = true;
            }
        }
        if changed {
            stats.columns_corrected += 1;
        }
    }

    Ok((corrections, stats))
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test denoise::tests::compute_corrections_flags_a_single_base_error 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/denoise.rs
git commit -m "feat(denoise): pileup pass computing per-read corrections + stats"
```

---

### Task 3: Write-back pass — apply corrections to a new BAM (aux preserved)

**Files:**
- Modify: `src/denoise.rs`
- Test: `src/denoise.rs` (test module)

**Interfaces:**
- Consumes: `Corrections` (Task 2).
- Produces (used by Task 4): `fn apply_corrections(in_bam: &Path, out_bam: &Path, corrections: &Corrections) -> Result<(u64, u64)>` returning `(reads_processed, reads_modified)`.

- [ ] **Step 1: Write the failing test**

Add to the test module (reuses `write_test_bam` from Task 2):

```rust
#[test]
fn apply_corrections_rewrites_seq_and_preserves_tags() {
    use rust_htslib::bam::record::Aux;
    let dir = std::env::temp_dir().join("himito_denoise_t3");
    std::fs::create_dir_all(&dir).unwrap();
    let inb = dir.join("in.bam");
    let outb = dir.join("out.bam");

    // One read "AAGAAA" with an aux tag we must preserve.
    {
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chrM");
        sq.push_tag(b"LN", &6usize.to_string());
        header.push_record(&sq);
        let mut w = Writer::from_path(&inb, &header, Format::Bam).unwrap();
        let mut rec = Record::new();
        let seq = b"AAGAAA";
        let quals = vec![30u8; 6];
        rec.set(b"rerr", Some(&CigarString(vec![Cigar::Match(6)])), seq, &quals);
        rec.set_tid(0);
        rec.set_pos(0);
        rec.push_aux(b"XY", Aux::I32(7)).unwrap();
        w.write(&rec).unwrap();
    }

    let mut corr: Corrections = HashMap::new();
    corr.insert(b"rerr".to_vec(), vec![(2usize, b'A')]);
    let (proc, modified) = apply_corrections(&inb, &outb, &corr).unwrap();
    assert_eq!(proc, 1);
    assert_eq!(modified, 1);

    // Verify SEQ was corrected and the aux tag survived.
    let mut r = Reader::from_path(&outb).unwrap();
    let rec = r.records().next().unwrap().unwrap();
    assert_eq!(rec.seq().as_bytes(), b"AAAAAA");
    assert!(matches!(rec.aux(b"XY").unwrap(), Aux::I32(7)));
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test denoise::tests::apply_corrections_rewrites_seq_and_preserves_tags 2>&1 | tail -20`
Expected: compile error (`apply_corrections` undefined).

- [ ] **Step 3: Write the implementation**

Add to `src/denoise.rs`:

```rust
fn apply_corrections(in_bam: &Path, out_bam: &Path, corrections: &Corrections) -> Result<(u64, u64)> {
    let mut reader = Reader::from_path(in_bam)
        .with_context(|| format!("cannot open BAM {in_bam:?}"))?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
        .with_context(|| format!("cannot create BAM {out_bam:?}"))?;

    let mut reads_processed = 0u64;
    let mut reads_modified = 0u64;

    for result in reader.records() {
        let mut rec = result.context("read BAM record")?;
        // Only primary mapped reads are candidates for correction.
        let correctable =
            !(rec.is_secondary() || rec.is_supplementary() || rec.is_unmapped());
        if correctable {
            reads_processed += 1;
            if let Some(edits) = corrections.get(rec.qname()) {
                // Extract owned parts, edit SEQ, rebuild (aux preserved by set()).
                let qname = rec.qname().to_vec();
                let cigar = rec.cigar().take();
                let mut seq = rec.seq().as_bytes();
                let qual = rec.qual().to_vec();
                for &(qpos, base) in edits {
                    if qpos < seq.len() {
                        seq[qpos] = base;
                    }
                }
                rec.set(&qname, Some(&cigar), &seq, &qual);
                reads_modified += 1;
            }
        }
        writer.write(&rec).context("write BAM record")?;
    }
    Ok((reads_processed, reads_modified))
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test denoise::tests::apply_corrections_rewrites_seq_and_preserves_tags 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/denoise.rs
git commit -m "feat(denoise): write-back pass applying corrections, aux preserved"
```

---

### Task 4: Orchestrator `start()` — reference loading, passthrough, stats output

**Files:**
- Modify: `src/denoise.rs`
- Test: `src/denoise.rs` (test module)

**Interfaces:**
- Consumes: `compute_corrections`, `apply_corrections`, `DenoiseStats` (Tasks 2–3).
- Produces (used by Task 5/6): `pub fn start(input: &PathBuf, output: &PathBuf, reference: &PathBuf, data_type: &str, vaf: f64, min_strand: u32, stats_out: Option<&PathBuf>) -> Result<()>`

- [ ] **Step 1: Write the failing tests**

Add to the test module:

```rust
#[test]
fn start_pacbio_is_byte_identical_passthrough() {
    let dir = std::env::temp_dir().join("himito_denoise_t4pb");
    std::fs::create_dir_all(&dir).unwrap();
    let inb = dir.join("in.bam");
    let outb = dir.join("out.bam");
    let reff = dir.join("ref.fa");
    std::fs::write(&reff, ">chrM\nAAAAAA\n").unwrap();
    write_test_bam(&inb, "chrM", 6, &[("r0", 0, b"AAGAAA", false)]);

    start(&inb, &outb, &reff, "pacbio", 0.01, 2, None).unwrap();
    assert_eq!(std::fs::read(&inb).unwrap(), std::fs::read(&outb).unwrap());
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn start_ont_denoises_and_writes_stats_json() {
    let dir = std::env::temp_dir().join("himito_denoise_t4ont");
    std::fs::create_dir_all(&dir).unwrap();
    let inb = dir.join("in.bam");
    let outb = dir.join("out.bam");
    let reff = dir.join("ref.fa");
    let statsf = dir.join("stats.json");
    std::fs::write(&reff, ">chrM\nAAAAAA\n").unwrap();
    let clean = b"AAAAAA";
    let mut reads: Vec<(&str, i64, &[u8], bool)> = vec![];
    for i in 0..9 { reads.push((["r0","r1","r2","r3","r4","r5","r6","r7","r8"][i], 0, clean, i % 2 == 0)); }
    reads.push(("rerr", 0, b"AAGAAA", false));
    write_test_bam(&inb, "chrM", 6, &reads);

    start(&inb, &outb, &reff, "ont-r10", 0.01, 2, Some(&statsf)).unwrap();

    let mut r = Reader::from_path(&outb).unwrap();
    let mut modified_ok = false;
    for rec in r.records() {
        let rec = rec.unwrap();
        if rec.qname() == b"rerr" {
            assert_eq!(rec.seq().as_bytes(), b"AAAAAA");
            modified_ok = true;
        }
    }
    assert!(modified_ok);
    let json = std::fs::read_to_string(&statsf).unwrap();
    assert!(json.contains("\"bases_modified\":1"));
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test denoise::tests::start_ 2>&1 | tail -20`
Expected: compile error (`start` undefined).

- [ ] **Step 3: Write the implementation**

Add imports at the top of `src/denoise.rs`:

```rust
use bio::io::fasta;
use log::info;
```

Then add:

```rust
pub fn start(
    input: &PathBuf,
    output: &PathBuf,
    reference: &PathBuf,
    data_type: &str,
    vaf: f64,
    min_strand: u32,
    stats_out: Option<&PathBuf>,
) -> Result<()> {
    // PacBio (and anything non-ONT) is a byte-identical passthrough.
    if !data_type.starts_with("ont") {
        std::fs::copy(input, output)
            .with_context(|| format!("passthrough copy {input:?} -> {output:?}"))?;
        info!("[denoise] data-type '{data_type}' is not ONT; passthrough copy.");
        return Ok(());
    }

    // Load reference contigs into memory (uppercased).
    let mut refs: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    let rdr = fasta::Reader::from_file(reference)
        .with_context(|| format!("cannot open reference {reference:?}"))?;
    for rec in rdr.records() {
        let rec = rec.context("read reference record")?;
        refs.insert(rec.id().as_bytes().to_vec(), rec.seq().to_ascii_uppercase());
    }

    let (corrections, mut stats) = compute_corrections(input, &refs, min_strand, vaf)?;
    let (reads_processed, reads_modified) = apply_corrections(input, output, &corrections)?;
    stats.reads_processed = reads_processed;
    stats.reads_modified = reads_modified;

    let rate = if stats.bases_examined > 0 {
        stats.bases_modified as f64 / stats.bases_examined as f64
    } else {
        0.0
    };
    info!(
        "[denoise] reads {}/{} modified; bases {}/{} modified ({:.2}%); \
         columns {} processed / {} corrected / {} skipped; alt sites preserved {}",
        stats.reads_modified, stats.reads_processed,
        stats.bases_modified, stats.bases_examined, rate * 100.0,
        stats.columns_processed, stats.columns_corrected, stats.columns_skipped,
        stats.alt_sites_preserved,
    );

    if let Some(path) = stats_out {
        let json = serde_json::to_string_pretty(&stats).context("serialize stats")?;
        std::fs::write(path, json).with_context(|| format!("write stats {path:?}"))?;
    }
    Ok(())
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test denoise:: 2>&1 | tail -20`
Expected: all `denoise::` tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/denoise.rs
git commit -m "feat(denoise): orchestrator with reference load, passthrough, stats JSON"
```

---

### Task 5: CLI subcommand `denoise`

**Files:**
- Modify: `src/main.rs` (add `Denoise` variant to `Commands` enum ~after the `Build` variant; add a match arm in `main`)

**Interfaces:**
- Consumes: `denoise::start` (Task 4).

- [ ] **Step 1: Add the subcommand definition**

In `src/main.rs`, in the `Commands` enum, immediately after the closing `},` of the `Build { ... }` variant, add:

```rust
    /// Denoise ONT reads in a BAM (SNV error correction) before graph construction.
    #[clap(arg_required_else_help = true)]
    Denoise {
        /// input BAM (coordinate-sorted)
        #[clap(short, long, value_parser, required = true)]
        input: PathBuf,

        /// output denoised BAM
        #[clap(short, long, value_parser, required = true)]
        output: PathBuf,

        /// reference FASTA (rCRS)
        #[clap(short, long, value_parser, required = true)]
        reference: PathBuf,

        /// data type; only ont-* is denoised, others pass through unchanged
        #[clap(short, long, value_parser, default_value = "ont-r10")]
        data_type: String,

        /// minimal VAF for an allele to be kept (below this it is corrected away)
        #[clap(long, value_parser, default_value_t = 0.01)]
        vaf: f64,

        /// minimal observations required on EACH strand for allele candidacy
        #[clap(long, value_parser, default_value_t = 2)]
        min_strand: u32,

        /// optional path to write denoise statistics as JSON
        #[clap(long, value_parser)]
        stats: Option<PathBuf>,
    },
```

- [ ] **Step 2: Add the match arm**

In `src/main.rs` `main()`, add a new arm (e.g. after the `Commands::Build { .. } => { .. }` arm):

```rust
        Commands::Denoise {
            input,
            output,
            reference,
            data_type,
            vaf,
            min_strand,
            stats,
        } => {
            if let Err(e) = denoise::start(&input, &output, &reference, &data_type, vaf, min_strand, stats.as_ref()) {
                eprintln!("Error running denoise: {e:#}");
                std::process::exit(1);
            }
        }
```

- [ ] **Step 3: Build and smoke-test the CLI**

Run:
```bash
cargo build --release 2>&1 | tail -3
./target/release/Himito denoise --help 2>&1 | grep -E "min-strand|--stats|--data-type"
```
Expected: builds; help shows `--min-strand`, `--stats`, `--data-type`.

- [ ] **Step 4: Commit**

```bash
git add src/main.rs
git commit -m "feat(denoise): add 'Himito denoise' subcommand"
```

---

### Task 6: `QuickStart` wiring — auto-denoise for ONT

**Files:**
- Modify: `src/main.rs` (`Commands::QuickStart` arm, ~lines 483–525)

**Interfaces:**
- Consumes: `denoise::start` (Task 4).

- [ ] **Step 1: Insert the denoise step and route build/call to the denoised BAM**

In the `QuickStart` arm, replace the block from the `graph_output` line through the `build::start(...)` call (currently ~lines 499–501) with:

```rust
            // For ONT, denoise the mt BAM so reads consolidate; PacBio is untouched.
            // Methylation (below) keeps reading the ORIGINAL mt_output.
            let build_bam = if data_type.starts_with("ont") {
                let denoised = output_prefix.with_extension("mt.denoised.bam");
                if let Err(e) = denoise::start(
                    &mt_output, &denoised, &reference_path, &data_type,
                    vaf_threshold as f64, 2, None,
                ) {
                    eprintln!("Warning: denoise failed ({e:#}); falling back to raw mt BAM.");
                    mt_output.clone()
                } else {
                    denoised
                }
            } else {
                mt_output.clone()
            };

            let graph_output = output_prefix.with_extension("gfa");
            // Preserve the historical default edge-read gate (2) for this combined pipeline.
            let _ = build::start(&graph_output, kmer_size, &build_bam, &reference_path, length_max, 2);
```

Then, in the `call::start(...)` call, change the strand-map BAM argument from `Some(&mt_output)` to `Some(&build_bam)` (consistent with the graph). **Leave the `methyl::start(&annotated_graph_output, &mt_output, ...)` call unchanged** — methylation must keep the original BAM.

- [ ] **Step 2: Build**

Run: `cargo build --release 2>&1 | tail -3`
Expected: builds cleanly.

- [ ] **Step 3: Run the existing suite (no regressions)**

Run: `cargo test 2>&1 | tail -5`
Expected: all tests pass (`test result: ok`).

- [ ] **Step 4: Commit**

```bash
git add src/main.rs
git commit -m "feat(denoise): QuickStart auto-denoises ONT; methylation keeps original BAM"
```

---

### Task 7: Eval-harness integration + end-to-end validation

**Files:**
- Modify: `lineage_eval/lineage_sim/run_himito.sh` (insert denoise before `build` for ONT)

**Interfaces:**
- Consumes: the `denoise` subcommand (Task 5).

- [ ] **Step 1: Insert the denoise step in `run_himito.sh`**

In `lineage_eval/lineage_sim/run_himito.sh`, immediately after `samtools index "$BAM"` and before the `"$HIMITO" build ...` line, add:

```bash
# Denoise ONT reads before graph construction (no-op for hifi/pacbio).
BUILD_BAM="$BAM"
if [[ "$DTYPE" == ont-* ]]; then
  DENOISED="$HDIR/aln.denoised.bam"
  "$HIMITO" denoise -i "$BAM" -o "$DENOISED" -r "$REF" -d "$DTYPE" \
    --vaf "$MIN_HF" --min-strand 2 --stats "$HDIR/denoise_stats.json"
  samtools index "$DENOISED"
  BUILD_BAM="$DENOISED"
fi
```

Then change the `build` and `call` invocations to use `"$BUILD_BAM"` instead of `"$BAM"` (the `build -i` argument and the `call --input-bam` argument). Leave any methylation invocation, if present, on `"$BAM"`.

- [ ] **Step 2: Rebuild and run one ONT eval cycle**

Run (reusing the existing simulated data dir):
```bash
cargo build --release 2>&1 | tail -2
cd lineage_eval/lineage_sim
./run_eval.sh --outdir /tmp/eval_denoise --profile ont-r10 --n-mutations 12 --total-depth 600 --seed 1 2>&1 | tail -20
```
Expected: the pipeline completes through `score_lineage` (no "No haplotypes remain" crash); a `denoise_stats.json` is written.

- [ ] **Step 3: Compare denoise vs. no-denoise metrics**

Inspect `/tmp/eval_denoise/metrics.tsv` and `denoise_stats.json`. Expected direction (per spec §8):
- `bases_modified` / `bases_examined` ≈ ONT error rate (a few %), NOT tens of %.
- variant precision rises sharply (FPs drop) vs. the no-denoise baseline (0.06 precision), recall stays ≥ baseline (~0.92).
- `lineage` completes.

Record the before/after numbers in the run log (no code change; this is the validation gate).

- [ ] **Step 4: Commit**

```bash
git add lineage_eval/lineage_sim/run_himito.sh
git commit -m "feat(denoise): wire ONT denoise into lineage_eval harness"
```

---

## Notes for the implementer

- The whole feature is ONT-only and additive; if any task regresses `cargo test`, stop and fix before proceeding.
- Do **not** commit the design spec or eval scratch data; only the files listed per task.
- `serde` derive and `serde_json` are already dependencies (see `src/build.rs`, `src/call.rs`); no `Cargo.toml` change is expected. If the compiler reports `serde::Serialize` is unavailable, add `serde = { version = "1", features = ["derive"] }` to `Cargo.toml` — but check first, it is almost certainly already present.
- Base-quality is preserved unchanged on corrected reads (downstream `build` uses sequence, not quality).
