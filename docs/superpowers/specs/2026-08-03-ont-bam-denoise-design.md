# ONT BAM Denoiser (`Himito denoise`) — Design

**Date:** 2026-08-03
**Status:** Approved (design), pending implementation plan
**Scope:** ONT-only, SNV-only read error correction as a BAM preprocessing step. PacBio unchanged.

## 1. Motivation

Himito's variant calling for ONT has low sensitivity **and** brittle precision. Root cause (traced in `src/build.rs` / `src/call.rs`):

- The anchor graph consolidates reads into edges by **exact inter-anchor sequence identity** (`create_edge_file`, dedup by `(src, dst, seq)`).
- An edge is aligned into a CIGAR — and can therefore yield a variant — only if it has **more than `min_edge_reads` reads** (`generate_cigar`, gate `allele_count <= minimal_read_count`); variant calling skips empty-CIGAR edges (`get_variant`).

ONT per-base error (~2–5%) gives each read a unique error pattern, so reads spanning a variant shatter into many **single-read edges** that never earn a CIGAR. Confirmed empirically: variants at 4.7–11.5% VAF with 500–660× depth (e.g. `m.14595A>G`, 72 ALT reads at 11.5%) produce **zero** CIGAR-bearing edges and are dropped. Lowering `min_edge_reads` does not help — it only multiplies false positives (179 → 391 FPs in the eval, recall unchanged) and floods the read×variant matrix so `lineage` fails with *"No haplotypes remain after --min-reads filtering."*

The fragmentation is structural to exact-sequence consolidation. **Denoising reads before graph construction removes it at the source:** once per-read errors are corrected, reads from the same haplotype become (near-)identical, consolidate into shared edges, and the existing pipeline works unchanged.

## 2. Goal & scope

A new `Himito denoise` step that uses a baldur-style per-site maximum-likelihood error model to **correct substitution (SNV) errors in individual ONT reads**, writing a denoised BAM. The existing `build` → `call` → `lineage` pipeline runs on that BAM **unchanged**, and thereby produces both the confident VCF and the clean per-read matrix `lineage` needs (satisfies the "both outputs" requirement for free).

**In scope (v1):**
- ONT data types (`ont-r9`, `ont-r10`); PacBio is a no-op passthrough.
- Substitution errors only.
- Standalone subcommand; auto-applied by `QuickStart` for ONT.

**Non-goals (v1):**
- Indels / homopolymer correction (require realignment).
- Formal likelihood-ratio test (LRT) for allele significance.
- Read phasing (that remains the graph/`lineage`'s job).
- Any change to methylation, assembly, `call`, `lineage`, or PacBio behavior.

## 3. Method (baldur-adapted, "lighter ML")

Reference: baldur v1.1.8 per-site ML caller. We use its **error model + candidacy + ML frequency estimate**, but replace the full Fisher-scoring/LRT machinery with EM + a min-VAF keep rule (sufficient for denoising; the downstream `call` still emits the VCF with its own statistics). Escalation to the full LRT is future work if validation shows mis-correction.

Processed **per reference column**, independently and in parallel.

### 3.1 Column construction
For a reference position, gather observations from **primary, mapped reads** aligned there via a CIGAR match op (`M`/`=`/`X`). Each observation is `(qname, qpos, base ∈ {A,C,G,T}, e = 10^(−Q/10), strand)`. Reads with an insertion/deletion/soft-clip/ref-skip at this column are **left untouched** there (SNV-only). Columns with **zero non-reference observations are skipped** (already consensus).

### 3.2 Candidacy
The reference allele is always eligible. A non-reference allele is eligible only if it has **≥ `min_strand` (default 2) observations on *each* strand** — baldur's strand-balance filter, which rejects strand-skewed systematic ONT errors.

### 3.3 ML frequency estimate (EM)
Estimate the frequency vector `f` over the eligible allele set under the per-observation model

```
P(obs = a | true = k) = (1 − e)      if k == a
                       = e / 3        otherwise
```

Bin observations by `(allele a, quality level q)` into counts `n[a][q]` with error `e_q` (baldur's `n_ij`), then EM initialized from the raw allele proportions (`f_k^(0) = observed_count(k) / N_obs`, restricted to eligible alleles):

- **E-step:** responsibility of bin `(a,q)` to true allele `k`:
  `r_k ∝ f_k · P(a | k, e_q)` (normalized over eligible `k`).
- **M-step:** `f_k = (Σ_{a,q} n[a][q] · r_k) / N_obs`.
- Iterate to `Δf < tol` or `max_iter`.

Drop alleles with `f < --vaf` (default 0.01); renormalize over the survivors. The reference allele is never dropped. This EM is over a ≤4-letter alphabet × ≤~94 quality bins → constant cost per column, independent of depth.

### 3.4 Per-read MAP correction
For each observation `(a, e)` at the column, compute the posterior over surviving alleles and set the read's base to the maximum-a-posteriori allele:

```
correct(a) = argmax_k  f_k · P(obs = a | true = k)
```

Behavior (verified analytically):
- A random error at a site with no surviving alt (its `f ≈ 0`) snaps to the local consensus — removed.
- A true heteroplasmy is preserved: for `f_a = 0.05, e = 0.05`, `f_a(1−e) = 0.0475` beats the ref pull `f_ref·e/3 = 0.016`, so the base is kept. Only genuinely low-quality (`e` large) or non-surviving-allele bases are corrected.
- At a heteroplasmic site with two surviving alleles, confident observations of each are preserved; third-base errors go to the higher-`f` allele.

Corrections are collected as `{ qname → [(qpos, new_base)] }` (only changed bases stored).

### 3.5 Write-back
Second pass over the BAM: for each **primary** read present in the corrections map, rebuild `SEQ` with the substitutions applied (base identity only; `CIGAR`, positions, length unchanged); all other reads (and secondary/supplementary/unmapped) pass through unchanged. Output is coordinate-sorted and indexed.

**Implementation note:** modifying `SEQ` via rust_htslib rebuilds the record core (`qname`/`cigar`/`seq`/`qual`); **all auxiliary tags must be re-attached** so nothing is dropped. `QUAL` is preserved as-is (downstream `build` uses sequence, not quality).

## 4. CLI & pipeline integration

New subcommand:

```
Himito denoise -i <in.bam> -o <out.bam> -r <ref.fa> -d <data-type> \
               [--vaf 0.01] [--min-strand 2] [--stats <file.json>] [--threads N]
```

- `-d pacbio` → **no-op passthrough** (copy input to output), so the subcommand is always safe to call.
- `QuickStart` **auto-inserts** denoise for `ont-r9` / `ont-r10`, between the mt-filtered BAM and `build`; `pacbio` skips it. This is the only place `QuickStart` behavior changes, and only for ONT.
- `lineage_eval/lineage_sim/run_himito.sh` gains one conditional line (denoise before `build`) for ONT profiles.

## 5. Methylation safety

Himito's methylation reads `MM`/`ML` tags, indexed by each read's own C positions; rewriting a base near a C could desync those tags. **Decision:** the denoised BAM feeds **only** `build` → `call` → `lineage`. Methylation aggregation continues to read the **original** mt BAM, so methylation output is bit-identical to today. In `QuickStart`, the denoised BAM is a separate file used for the variant/lineage path; the methylation path keeps its current input. (The denoiser still preserves all aux tags for general correctness, per §3.5.)

## 6. Statistics & reporting

The correction pass (§3.4) visits every base and decides each one, so summary statistics are free. An **aggregate summary is always printed to the log** (stderr) at the end of a run; `--stats <file>` additionally writes the same numbers as machine-readable **JSON** (consumed by the `lineage_eval` harness to track denoise behavior across runs). Reported:

- `reads_processed` (primary mapped) and `reads_modified` (≥1 base changed).
- `bases_examined` and `bases_modified` (count **and** rate). The rate should sit near the ONT error rate (~2–5%) — a built-in sanity check; a wildly high rate signals over-correction.
- `columns_processed`, `columns_corrected` (≥1 change), `columns_skipped` (pure-reference).
- `alt_sites_preserved` — surviving non-reference alleles after the min-VAF keep rule (i.e. real heteroplasmies kept, not corrected away).
- `substitution_matrix` — the 12 directed counts (A>C, A>G, …, T>G) of applied corrections, exposing systematic ONT error directions.

Per-column and per-read breakdowns are out of scope for v1 (see §11).

## 7. Performance

No pairwise read comparison: every read contributes to a per-site `(allele × quality)` count summary and is corrected against it, so there is **no O(N²) term**.

Let `B` = total aligned bases `= N × mean_read_len`, `C` = reference columns (~16,569 for mito).

- Pileup build: `O(B)`.
- EM: `O(C · A · Q · iter)` with `A ≤ 4`, `Q ≤ ~94`, `iter ~20–50` → constant per column, independent of depth.
- Correction + write-back: `O(B)`.

**Total `O(B)` — linear in read count** for fixed read length. For the `maximal_mt_depth = 50,000`-read cap, `B ≈ 0.1–0.8 Gbp` depending on ONT read length; estimated wall time **seconds to ~30 s**, I/O-bound (BGZF decompress/recompress usually dominates) and embarrassingly parallel over columns. Peak memory a few hundred MB (transient per-column buffers + a corrections map ≈ error-rate × `B`).

## 8. Validation

Reuse `lineage_eval/lineage_sim`: `simulate → [denoise] → build/call/lineage → score_lineage`, comparing **no-denoise vs denoise** on the same simulated BAM:

- Variant detection precision / recall / F1 (from the VCF vs truth).
- Whether `lineage` completes (haplotypes consolidate rather than fragmenting to singletons).
- Tree accuracy: RF and quartet distance (already in `score_lineage.py`).
- Read-level check: fraction of reads that become identical to their clone's true sequence after denoising.

**Success criteria:** precision recovers sharply (false positives drop), recall ≥ current, and `lineage` no longer fails with "No haplotypes remain."

## 9. Files touched

- **New:** `src/denoise.rs` (column model, EM, correction, BAM I/O), plus unit tests.
- `src/main.rs`: `Denoise` subcommand + handler; `QuickStart` wiring for `ont-*`.
- `lineage_eval/lineage_sim/run_himito.sh`: conditional denoise step for ONT.
- Dependencies: `rust_htslib` (already present); no new crates expected.

## 10. Testing

Unit tests on synthetic pileup columns:
- Pure-ref column with sparse errors → all errors corrected to ref.
- Clean heteroplasmy (e.g. 60/40, both strand-balanced) → both alleles preserved; third-base errors removed.
- Low-VAF real allele passing candidacy → preserved; sub-threshold / strand-skewed candidate → corrected away (documents the over-correction boundary).
- Reads with indel/clip at the column → untouched at that column.
- Tag preservation: `MM`/`ML` and other aux tags survive a SEQ rewrite.
- `-d pacbio` → byte-identical passthrough.
- Statistics: on a known synthetic input, `bases_modified`, `alt_sites_preserved`, and the substitution matrix match hand-computed values; `--stats` JSON is well-formed and matches the logged summary.

## 11. Future work

- Indel / homopolymer correction (requires local realignment).
- Full Fisher-scoring ML + per-allele LRT (calibrated significance) if the lighter keep-rule over- or under-corrects.
- Optional single-pass streaming write-back if the two-pass I/O proves to dominate.
- Finer statistics granularity: per-column TSV (`--stats-per-site`) and/or per-read correction distributions (§6 Q6 options 2–3).
