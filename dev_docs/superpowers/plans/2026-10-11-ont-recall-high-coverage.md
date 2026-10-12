# ONT Variant Recall at High Coverage — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Do the tasks **in order**, one commit per task. Each task starts from a green `cargo test`. Do **not** start anything under "Phase 2" or "Needs a decision" without owner sign-off.

**Goal:** Raise low-HF recall for ONT `quick-start` at high depth (≈420× on the HapMap mixture) without losing precision. Today the low-frequency calls are thrown away by fixed thresholds that don't account for depth. They aren't missing for lack of evidence.

**Revision 2 (2026-10-11, after the owner's validation run).** The T1 mechanism is confirmed on the real binary: `ont-r10 --permutation-rounds 10000 --indel-false-threshold 0.05` gives **88 TP / 1 FP (recall 0.727)**, against 37 TP before. This matches the offline projection exactly. At first the run looked like "no improvement", because `ont_benchmark/benchmarking.sh` scores `HG00097_ont.vcf` (the old denoised output), whatever `-o` the run used. T0 fixes the harness. The cost is runtime: the permutation step went from 74 s to 705 s. So T1 is now "get the 10k-round result by default, at roughly the 1k-round cost", not "raise the rounds".

**Revision 3 (2026-10-12, full-depth HapMap).** T1–T5 are done (see **Results**). On the merged HapMap BAM (`merge215k_r10`, 50k-read cap, ~3,600–5,000×), `ont-r10` scores **97 TP / 2 FP / 24 FN (recall 0.802, precision 0.980)**. Read-linkage rescue of sub-2% alleles should bring it to about **113–114 TP with no new FP**. New tasks R3-T0…R3-T5 and the implementation order are in **Revision 3** at the end of this file. Do them in order; they don't depend on Phase 2.

**Starting point:** branch `hs_lineage_posterio` at `4f77488` or later. md5s at review time: `src/call.rs` `6d992919dc98ca5d576f448b0f7a6810`, `src/denoise.rs` `d3b544da05831a5522b05242c0cc0fa3`, `src/denoise_indel.rs` `3443ac99f973a48cb42854d7feee3e9a`, `src/main.rs` `411495e937779b45296abf080580f3e2`. If a md5 differs, re-locate every line number below by function name before editing.

**Tech stack:** Rust (`rayon`, `rand 0.9`, `statrs 0.16`, `adjustp`), with Python helpers (`numpy`, `pysam`) under `scripts/ont_recall/`.

**Data:** `ont_benchmark/` (untracked; do not commit data). `HapMap_ont.bam` is the SMaHT HapMap 6-line mixture, ONT R10, about 420× on chrM. Truth: `SMHTHAPMAP6_GRCh38_v1.0.0_mito_benchmark.vcf.gz`, with 121 records and the nominal mixture HL in `FORMAT/HL`.

**Oracles / tools (written during the review, in `scripts/ont_recall/`, uncommitted; T0 commits them):**

| script | what it is |
|---|---|
| `score_vs_truth.py` | Exact (POS, REF, ALT) scorer on PASS records. It reproduces `rtg vcfeval --squash-ploidy` exactly on both baseline runs (37/1/84 and 79/1/42). |
| `perm_retest.py` | Offline re-implementation of `call`'s permutation test on `<prefix>.raw_matrix.csv`. `--mode adaptive` is **the oracle for T1**. |
| `fn_pileup_audit.py` | Raw-vs-denoised read evidence at low-HL truth SNVs. |
| `run_hapmap_bench.sh` | quick-start → PASS filter → vcfeval → scorer, run on **the VCF this run wrote**, plus timing of the permutation step. |

---

## Why: the evidence

### Measured on HapMap ONT (exact scorer, which equals vcfeval here)

| run | TP | FP | FN | recall | precision | strand-bias + permutation step |
|---|---|---|---|---|---|---|
| `ont-r10` defaults | 37 | 1 | 84 | 0.306 | 0.974 | 74 s (580 tested) |
| `ont-denoised` defaults | 79 | 1 | 42 | 0.653 | 0.988 | 2 s (96 tested) |
| **`ont-r10 --permutation-rounds 10000 --indel-false-threshold 0.05`** | **88** | **1** | **33** | **0.727** | **0.989** | **705 s** |

The single FP in every run is `15622T>C`. It has 101 strand-balanced reads (54/47, about 24%) and is very likely a real heteroplasmy that the assembly-derived truth set lacks.

### Root causes (all verified in code and reads)

1. **Permutation-test p-value floor × BH (the main cause, T1).** The empirical p-value can't go below `1/(R+1)`. On raw r10 at 420×, `m = 580` candidates are tested; most are noise admitted by `minimal_ac 2` / `-v 0.01`. Exactly 54 hit the floor, all truth. BH's smallest q is then `(1/1001)·580/54 = 0.0107 > 0.01`, so **every variant under 20% HF is rejected**, including 10% variants with 50 reads. More depth means more noise candidates and a larger `m`, so recall *falls* as coverage rises. Denoising "helped" mostly by shrinking `m` to 96, not by making real variants detectable. Today's only way out is more rounds for every variant, and the cost scales linearly (705 s).
2. **The strand-bias null is library-wide (T3).** `filter_strand_bias` uses the library forward fraction (0.498). `7805G>A` (truth 10%, alt 45/27) gets p = 0.034 and is flagged. Against the reads covering the site (0.537), p = 0.156 and the effect is 0.09. At high depth, small skews become "significant". `denoise.rs::fit_site` already uses the column-local null.
3. **The co-occurrence test vetoes strong, isolated variants (T4).** The summed-Jaccard statistic is pulled *below* its null by anti-correlation with other lineages. `12969C>T` (about 15%, 54 graph reads, 39/21 strands) sits only on HG005 molecules: 59/59 of its carriers carry the HG005 marker `12618A`, and 0 carry any HG02622 marker. It gets q = 0.99. `9342G>A` (about 18%, 38/41) is removed the same way.
4. **`Potential_Artifact` re-filters indels the co-occurrence test already accepted (T5).** After T1, the only small indels below 0.05 HF that survive the test are `285CAA>C` (0.018) and `302A>ACC` (0.031). Both are truth.
5. **Fixed-fraction gates in denoise ignore depth (Phase 2: T6, T7).** `DENOISE_KEEP_VAF = 0.03` (`main.rs`, used in `denoise::fit_site`) erases every HG002 (~2%) and HG02486 (~1.5%) allele. It also erases 7 of 34 HG02257 alleles, which really sit at 3.5–4.5%, right on the line. Sites are decided independently: 13/13 reads carrying `9540C` also carry `9545G`, yet denoise kept `9540C` and rewrote `9545G` on 12 of them. Indel denoising erases `8271AC>A` (10%): `vaf_floor(L=5) = max(0.05, 3·0.0506) = 0.152`. `lineage_sim/README.md` documents the same thing on simulations: "Raising `--total-depth` from 1000 to 3000 did not move `var_recall`".

### The 33 FNs left after the validated quick fix

| class | n | example(s) | addressed by |
|---|---|---|---|
| HG00438, nominal 0.5% (2–8 raw reads, inside the noise; 99th-percentile noise count is 10 reads) | 12 | 3777, 5585 | out of scope |
| HG002 / HG02486 private, ~1.5–2% (4–12 raw reads; sum-Jaccard q 0.09–0.34) | 13 | 309, 456, 4336, 491 | Phase 2 T8 (experiment) |
| no alt reads in the data (≤1 read, or absent despite HL 0.835) | 4 | 13376, 11087, 10944, 14861 | none (truth/sample discordance) |
| REF contains `N` (`3106CN>C`; `filter_vcf_record` drops REF with N) | 1 | 3106 | report with `--exclude-pos 3106` |
| strand-bias false flag | 1 | 7805G>A | **T3** |
| co-occurrence veto | 1 | 12969C>T | **T4** |
| `Potential_Artifact` on tested indels | 2 | 285CAA>C, 302A>ACC | **T5** |
| polyC insertion, q 0.021 | 1 | 16182A>AC | none |

**Expected after T1–T5 on HapMap `ont-r10`, default flags:** TP 92, recall 0.760. Against truth, FP becomes 2: `15622T>C` plus `9342G>A` (T4), both likely real. The ceiling given the reads is about 104.

---

## Decisions (defaults the agent uses unless the owner edits this table)

| # | Item | Default in this plan |
|---|---|---|
| 1 | Permutation algorithm | **Adaptive by default** (T1). Keep today's behaviour behind `--permutation-mode fixed`, for A/B runs. |
| 2 | Adaptive constants | Stop a variant after `PERM_STOP_AFTER = 20` exceedances. Floor headroom `PERM_HEADROOM = 10`. Cap `PERM_MAX_ROUNDS = 100_000`. `--permutation-rounds` becomes the *minimum* rounds (default stays 1000). |
| 3 | RNG | Seeded per variant: `PERM_SEED = 42` XOR FNV-1a-64(variant name), using `rand::rngs::StdRng`, so results are deterministic under rayon. |
| 4 | Strand-bias effect gate | `STRAND_BIAS_MIN_EFFECT = 0.20` (absolute difference in forward fraction, alt vs covering reads). Use the column-local null when ≥ `STRAND_LOCAL_MIN_READS = 20` strand-resolved covering reads, otherwise the library fraction. |
| 5 | Co-occurrence exemption | SNPs only, `HF ≥ COOC_EXEMPT_HF = 0.10`, `AC ≥ COOC_EXEMPT_MIN_AC = 20`, and **not** strand-flagged after T3. |
| 6 | `--indel-false-threshold` | Becomes `Option<f64>` with a data-type preset: `ont-r9`/`ont-r10`/`ont-denoised` → **0.015**; `pacbio` → 0.1 (unchanged). |
| 7 | Phase 2 (T6–T8) | Not started until the owner signs off after T2. |

---

## Global constraints

- Write the test first and watch it fail, then implement, then run `cargo test` for the whole crate and get green.
- Don't weaken or delete an existing test. If an existing test calls a function whose signature you change, update **only the call site** (e.g. pass an empty map) and leave its assertions alone. If an assertion has to change because of the intended fix, say why in the commit message.
- Cite the task in comments: `// ONT-recall T3`. This avoids clashing with the 2026-10-11 genotyped-matrix plan's T0–T9.
- Every new constant goes in `main.rs` next to `DENOISE_KEEP_VAF`, with a doc comment that gives the evidence number it came from (copy it from this plan).
- `cargo build --release` must not add warnings in the files you touch.
- If something is ambiguous, stop and write `> QUESTION:` into the task. Don't guess.
- After each of T1, T3, T4 and T5, run the HapMap harness for **both** `ont-r10` and `ont-denoised` and add a row to **Results**.

---

### Task 0: Harness, oracles, baseline (no Rust changes)

- [ ] `git add scripts/ont_recall/*.py scripts/ont_recall/run_hapmap_bench.sh` and commit: `scripts: ONT recall harness and offline oracles (ONT-recall T0)`.
- [ ] In `ont_benchmark/benchmarking.sh` (untracked, owner's file), make the `bcftools view` input `${PREFIX}.vcf` and the vcfeval output dir `${OUTPUT_DIR}/${PREFIX}.eval`, or just switch to `run_hapmap_bench.sh`. Don't commit `ont_benchmark/`.
- [ ] `cargo build --release && cargo test`, and record the pass count here: `____ passed`.
- [ ] From `ont_benchmark/`, reproduce the three baseline rows of the evidence table:
  ```bash
  S=../scripts/ont_recall/run_hapmap_bench.sh
  $S ont-r10      base_r10
  $S ont-denoised base_den
  $S ont-r10      base_r10_perm10k --permutation-rounds 10000 --indel-false-threshold 0.05
  ```
  Expected: 37/1/84, 79/1/42 and 88/1/33. Permutation step about 74 s, 2 s and 705 s. If any row differs by more than ±1 TP, stop and report.
- [ ] Oracle sanity check: `python3 ../scripts/ont_recall/perm_retest.py --matrix base_r10.raw_matrix.csv --truth $TRUTH --mode fixed --rounds 1000` should print `truth 0 non-truth 0, min_q 0.0107`. With `--mode adaptive` it should print `truth 54 non-truth 0`, with about 7.2e5 total draws and final R_target about 10.7k.

### Task 1: Adaptive, floor-aware permutation test (highest impact)

**Files:** `src/call.rs`: `shuffle_genotypes` (~1186), `genotype_jaccard_similarity` (~1202), `sum_jaccard_against_others` (~1232), `get_null_distribution` (~1257), `empirical_p_value` (~1332), `permutation_test` (~1512); `src/main.rs`: the `permutation_rounds` help text in `QuickStart` (~191) and `Call` (~432), plus a new `--permutation-mode` on both.

**Algorithm.** It must match `perm_retest.py --mode adaptive`.
1. Admission is unchanged. A variant is tested when `allele_count/coverage ≤ heteroplasmic_error_threshold`, it is not an indel longer than 50 bp, it is not a PacBio SNP, and its matrix alt frequency is `≤ permutation_frequency_threshold`; the last condition is today's "has a null". `m = |tested|`.
2. **Stage 1:** for each tested variant (rayon), draw in batches of 100 up to `min_rounds = --permutation-rounds`, and **stop early** once exceedances reach `PERM_STOP_AFTER`. An exceedance is a null draw `≥ obs`, the same comparison as `empirical_p_value`.
3. **Stage 2 (resolve the floor):** loop at most 5 times. Set `k0 = #tested with 0 exceedances`. If `k0 == 0`, break. Set `R_target = min(PERM_MAX_ROUNDS, ceil(PERM_HEADROOM · m / (alpha · k0)))`, where `alpha = p_value_threshold`. If `R_target ≤` the current target, break. Otherwise extend only those `k0` variants to `R_target`, still with the early stop.
4. `p = (exceed + 1) / (draws + 1)`, then BH and exclusion exactly as today.
5. Log one line per stage, to stdout like the rest of `call`: `m`, `k0`, `R_target`, total draws and wall time.

**Speed.** Represent rows as bitsets, `struct BitRow { alt: Vec<u64>, cov: Vec<u64>, covered_idx: Vec<u32>, n_alt: u32 }`. Value `> 0.5` means alt, NaN means uncovered, 0 means covered ref; `alt ⊆ cov` always holds. For a shuffled row `x` against row `y`:
`inter = popcnt(x_alt & y.alt)`, `union = popcnt((x_alt & y.cov) | (y.alt & x_cov))`, `J = if union == 0 { 0 } else { inter / union }`. This equals `genotype_jaccard_similarity` by construction, and the test below proves it. Build a shuffle by choosing `n_alt` of `covered_idx` with a partial Fisher–Yates on a scratch buffer, keeping one per rayon task.

- [ ] **Step 1 (failing tests)** in `call.rs` `mod tests`:
  - `bit_jaccard_matches_f64_jaccard`: 200 random row pairs over 300 reads, about 20% NaN and about 10% alt, seeded. For every pair, `jaccard_bits` equals `genotype_jaccard_similarity` (exact `==` on the f64 result).
  - `bit_shuffle_preserves_alt_count_and_coverage`: `popcnt(alt) == n_alt`, `alt & !cov == 0`, over 1000 shuffles.
  - `fixed_permutation_cannot_resolve_the_bh_floor` (characterization of today's behaviour). The matrix has 400 reads, 10 identical "lineage" rows (alt in the same 24 reads) and 120 independent noise rows (8 random alts each), seeded. `permutation_test(.., mode Fixed, rounds 1000, alpha 0.01, ..)` keeps **0** of the 10 lineage rows, because the q floor is `(1/1001)·130/10 = 0.013`. (Checked offline: fixed keeps 0 of 10 with min q 0.0130; adaptive keeps 10 of 10 and 0 noise, with R_target 13000, and noise rows with p > 0.1 stop by 200 draws.)
  - `adaptive_permutation_resolves_the_bh_floor`: same matrix, adaptive mode. All 10 lineage rows are kept, ≤ 1 noise row is kept, and the reported `R_target` is ≥ 13000.
  - `adaptive_permutation_stops_clear_nulls_early`: every noise row whose final p > 0.1 used `< min_rounds` draws. Expose `draws` in a small `PermOutcome { exceed: u32, draws: u32 }`.
  - `adaptive_permutation_is_deterministic`: two runs give byte-identical p-values.
  - All existing `permutation_test_*` and `null_distribution_*` tests must stay green unchanged. They exercise the fixed path, so keep `get_null_distribution`/`empirical_p_value`, which are used by `--permutation-mode fixed`.
- [ ] **Step 2 (implement).** Add `pub enum PermutationMode { Adaptive, Fixed }` (clap `ValueEnum`, default `Adaptive`). Thread it from both CLI arms into `call::start` → `permutation_test`. Put the new code in a section `// ---- adaptive permutation (ONT-recall T1) ----` above `permutation_test`. The help text for `--permutation-rounds` should read: "minimum rounds per tested variant; with the default adaptive mode, variants at the p-value floor are extended automatically until BH can resolve them".
- [ ] **Step 3 (acceptance).** `cargo test`. Then run `$S ont-r10 t1_r10` with **no** `--permutation-rounds`:
  - With no extra flags, expect TP 86, because `247GA>G` and `8271AC>A` stay `Potential_Artifact` until T5. Adding `--indel-false-threshold 0.05` should give the validation run's 88. FP ≤ 1 either way.
  - Strand-bias + permutation step ≤ 150 s (validation run: 705 s; the 1000-round baseline: 74 s).
  - Oracle parity: the set of truth variants passing equals `perm_retest.py --mode adaptive` on `t1_r10.raw_matrix.csv`, give or take borderline variants with q in [0.005, 0.02]. List any that differ.
  - `$S ont-denoised t1_den`: TP ≥ 79, FP ≤ 1.
  - Commit: `call: adaptive floor-aware permutation test (ONT-recall T1)`.

### Task 2: Re-check `ont-r10` on the lineage simulations (report only, no code)

`lineage_sim/README.md` ("Use `ont-denoised`, not `ont-r10`") says ont-r10's VCF "comes back empty regardless of `-f`" on 94,694-edge graphs. That is the same `m`-inflation mechanism T1 fixes.

- [ ] Run `bench_standard.sh` (30 seeds, frozen config) with `--profile ont-r10` and with `--profile ont-denoised`, both on the T1 build. Also run one low-HF cell: `--n-mutations 15 --sim-min-hf 0.03 --frequency-threshold 0.02 --denoise-vaf 0.01`. Report `var_recall`, `var_precision`, `ad_f1` and the permutation-step time.
- [ ] If `ont-r10` now beats or matches `ont-denoised` on variant recall at equal precision, add the measured table to that README section. Don't change benchmark defaults. Feed the result into "Needs a decision" #1.

### Task 3: Strand bias, column-local null plus an effect-size gate

**Files:** `src/call.rs`: `filter_strand_bias` (~1425), its caller in `start` (~1725–1745), the `##FILTER=<ID=Strand_bias…>` header line (~675); `main.rs` for the constants.

- [ ] **Step 1 (failing tests):**
  - `strand_bias_uses_column_local_null`. Library fraction 0.498. The `cover_record` reads for the variant are 54% forward. The alt is 45 forward / 27 reverse, which is the `7805G>A` case. It must **not** be flagged. Today it gets p = 0.034 against the library and is flagged.
  - `strand_bias_needs_an_effect_not_just_a_p_value`. 1000 alt reads at 600/400 against a 50% local null: p ≈ 1e-10 but effect 0.10, so not flagged. This is the high-coverage case.
  - `strand_bias_still_flags_single_strand_artifacts`. Alt 1/55 against a local 0.536 (the `33CG>C` case): flagged under both `Flag` and `Drop`.
  - `strand_bias_falls_back_to_library_fraction_with_few_covering_reads`. Fewer than 20 strand-resolved covering reads means the library fraction is used.
  - Existing `strand_bias_*` tests: update call sites only, passing an empty `cover_record`. Their alt is 10/0, so the effect is 0.5 and they still flag.
- [ ] **Step 2 (implement).** Add `cover_record: &HashMap<String, HashSet<String>>` and `min_effect: f64` to `filter_strand_bias`. Local fraction = forward/(forward+reverse) over the covering reads found in `strand_map`; the reference pseudo-read isn't in `strand_map` and is skipped naturally. Flag iff `p < p_threshold && |alt_fwd_frac − null_frac| ≥ min_effect`. Log both fractions in the existing `println!`. Change the header description to "…relative to the reads covering the site".
- [ ] **Step 3 (acceptance):** `$S ont-r10 t3_r10`: `7805G>A` is PASS (+1 TP), `33CG>C` is still `Strand_bias`, and `1598G>A` (43/0, if it reaches the VCF) is flagged. `$S ont-denoised t3_den`: `7805` is PASS and nothing else changes. Commit: `call: column-local strand-bias null with effect-size gate (ONT-recall T3)`.

### Task 4: Don't let the co-occurrence test veto strong, strand-balanced SNPs

**Depends on T3** (the exemption trusts the strand flag). **File:** `src/call.rs` `permutation_test`, the per-variant closure (~1527–1570).

- [ ] **Step 1 (failing test)** `cooccurrence_test_exempts_strong_balanced_snps`. The fixture mirrors 12969C>T: X is a heteroplasmy on the **majority** lineage. Use 300 reads in three exclusive lineages: A, 240 reads with 3 markers (HF 0.8, so above the null threshold and untested); B, 36 reads with 15 markers; C, 24 reads with 15 markers. X is a SNP on 45 random A reads (HF 0.15, AC 45, `filter: None`). This was checked offline: today X gets p = 1.000, because its zero overlap with the 30 B/C markers outweighs its small excess over the 3 A markers. With the exemption it is kept. The same X with `filter = Some("Strand_bias")` is still tested and rejected, and an indel with the same evidence is still tested. (A fixture where X sits on a *minority* lineage does not fail today, so don't use one.)
- [ ] **Step 2 (implement).** In the admission closure, before the PacBio-SNP check, add `if is_snp && current_variant.filter.is_none() && frequency >= COOC_EXEMPT_HF && current_variant.allele_count >= COOC_EXEMPT_MIN_AC { return None; }`. Add a comment that the statistic punishes variants on the majority lineage, citing 12969C>T (q 0.99).
- [ ] **Step 3 (acceptance).** In the review's offline check of `ont-r10`, only two non-PASS SNPs fall in this band, `12969C>T` (truth) and `9342G>A` (not in truth), while `1598G>A` is excluded by the strand flag. Expect `12969` PASS (+1 TP) and `9342` PASS (+1 FP against truth; list it under "likely real" in Results), and **no other new PASS**. `ont-denoised` should give the same two. Commit: `call: exempt strong strand-balanced SNPs from the co-occurrence veto (ONT-recall T4)`.

### Task 5: Data-type preset for `--indel-false-threshold`

**Files:** `src/main.rs` (`QuickStart` ~184, `Call` ~426: make it `Option<f64>`), `src/call.rs` (new `resolve_indel_false_threshold(data_type, explicit) -> f64`, next to `resolve_thresholds`).

- [ ] **Step 1 (failing tests):** the resolver returns 0.015 for `ont-r9`, `ont-r10` and `ont-denoised`, and 0.1 for `pacbio`; an explicit value always wins. `format_vcf_record` with threshold 0.015 marks a 1-bp deletion at HF 0.018 as PASS, and one at HF 0.011 as `Potential_Artifact`.
- [ ] **Step 2 (implement)** and update the help text: "indels below this HF are flagged Potential_Artifact; ONT default 0.015 because small indels under 20% HF must already pass the co-occurrence test".
- [ ] **Step 3 (acceptance):** `ont-r10`: `285CAA>C`, `302A>ACC`, `247GA>G` and `8271AC>A` all PASS, so TP should reach **92**. `ont-denoised`: `247GA>G` is PASS and `540A>AC` (non-truth, HF 0.011) is still flagged. Commit: `call: data-type preset for indel_false_threshold (ONT-recall T5)`.

---

## Phase 2 (owner sign-off after T2; mainly matters if `ont-denoised` stays the discovery path)

### Task 6: Denoise keeps minor alleles that are linked to kept minor alleles

**Files:** `src/denoise.rs`: `fit_site` (~98), `map_allele` (~247), `compute_corrections` (~391; the SNV loop is ~683–712, the deferred-decision block ~722–770), `DenoiseStats`; `src/main.rs`: a new `quick-start --denoise-keep-vaf` (default `DENOISE_KEEP_VAF`), so it can be swept without rebuilding.

- [ ] **6a.** Expose `--denoise-keep-vaf`, with a resolver test and help text. This is a separate commit.
- [ ] **6b. Defer, then rescue.**
  - In `fit_site`, mark an allele **pending** when it is eligible (passes the per-strand and strand-bias gates) and `PENDING_FLOOR (0.005) ≤ f < vaf`.
  - In the SNV loop, observations of a pending allele are *not* corrected yet. Record `(qname, refpos, would-be base)` in `deferred[(contig, refpos, allele)]`.
  - For every read, record its minor alleles: kept non-ref alleles with `f < homoplasmic_vaf`, plus pending ones.
  - Record each read's aligned span, `(min, max refpos)` over the columns where it appears.
  - After the pileup, for each pending allele P, look at partners M: kept minor alleles, more than 20 bp from P, on reads spanning both. Over those reads, `a = |P ∧ M|` and the hypergeometric upper tail is `P(X ≥ a)`. Rescue P if `a ≥ 3` and `min_p · n_partners ≤ 1e-4`.
  - Repeat up to 3 rounds, adding rescued alleles to the partner set.
  - Rescued alleles keep their observations. For everything else, apply the deferred substitutions, then sort every touched read's `subs` by position. **Note:** today the sort runs only inside `if iopts.enabled`, ~767.
  - Add `stats.alt_alleles_rescued_by_linkage`.
  - Tests: a 9540/9545-style fixture (an allele at 2.8% that rides on the same reads as a kept 3.5% allele is preserved); an unlinked 2.8% allele is still corrected; a linked partner within 20 bp doesn't count; output is deterministic.
- [ ] **Acceptance:**
  - HapMap `ont-denoised`: the 7 erased HG02257 alleles (3439, 3552, 4715, 9545, 10873, 13263, 15043) survive in `mt.denoised.bam`; check with `fn_pileup_audit.py`, and the raw and denoised counts should match. TP ≥ 86, with FP unchanged.
  - `bench_standard.sh` (30 seeds): `var_precision` and `ad_f1` must not drop by more than 1 seed's worth against `baseline_standard.tsv`.
  - Low-HF cell: `var_recall` should go up.

### Task 7: Depth-aware indel candidacy in denoise

**File:** `src/denoise_indel.rs` `fit_indel_site` (~426, the gate `vaf < vaf_floor(l, o)` ~513).

- [ ] Candidacy becomes `vaf ≥ max(o.vaf_abs_min (new, 0.02), 1.5·eps)` and binomial upper tail `P(X ≥ n | depth, eps) ≤ o.excess_alpha (new, 1e-3) / n_sites`. Keep `vaf_floor` for the `protect_vaf` logic. Write the tests first: a 10% one-base deletion in a 5-C run at depth 431 survives (`8271AC>A`; today the floor is 0.152); a 6% one at depth 60 does not.
- [ ] **Acceptance:** `8271AC>A` survives denoise. `indel_alt_alleles_kept` grows by less than 2× on HapMap. HapMap `ont-denoised` FP is unchanged. The `bench_hp.sh` metrics don't regress.

### Task 8 (experiment, report only): partner-specific co-occurrence for small lineages

HG002's private alleles (~2%, 4–12 reads) get sum-Jaccard q 0.09–0.34, because the sum over about 600 rows drowns a lineage that has only a few partners.

- [ ] Add a hidden `--cooccurrence-stat max-partner`. For each tested variant, take the minimum over partners (≥ 2 shared reads) of the hypergeometric upper-tail p on co-covered reads, times the number of partners. Then BH as usual.
- [ ] Prototype it in `perm_retest.py` first. Report HapMap truth/non-truth passes and the sim metrics. The owner decides.

---

## Needs a decision (don't implement without sign-off)

1. **Which path is the default for ONT variant discovery?** After T1–T5, HapMap is expected at `ont-r10` 92 TP against `ont-denoised` about 81. Option: discover on the raw graph (`ont-r10`) and keep `--genotype-from-bam` for the lineage matrix. Denoising would then only serve graph size and lineage, and T6/T7 lose priority. Decide using T2.
2. **ONT `indel_false_threshold = 0.015`** rests on one sample. Confirm on the sims (`bench_hp.sh`) and on one HPRC ONT sample before it ships.
3. **Pending–pending linkage in T6** (needed to rescue HG002/HG02486 in denoise). It raises the risk of rescuing correlated errors on "bad reads".
4. **A per-site error model (panel of normals)** for 1–2% heteroplasmy at ≥ 400×. Count-only rules can't separate these from systematic error.
5. **Truth-set gaps to report, not tune against:** `15622T>C` (~24%, 101 balanced reads), `9342G>A` (~18%), `12969C>T` (present at ~15% on HG005 molecules, not 83.5%), `14861G>A` (absent), and `3106CN>C` (unscoreable, so use `--exclude-pos 3106`).

---

## Results

### cargo test
- baseline: **392** passed / after T5: **406** passed (1 ignored)

### HapMap ONT (`run_hapmap_bench.sh`; exact scorer; TP/FP/FN, recall, permutation-step seconds)

| step | ont-r10 | ont-denoised | notes |
|---|---|---|---|
| baseline (defaults) | 37/1/84, 0.306, 65 s | 79/1/42, 0.653, 2 s | T0 reproduce |
| validation: r10 `--permutation-rounds 10000 --indel-false-threshold 0.05` | 88/1/33, 0.727, 756 s | — | owner/T0 |
| T1 | 86/1/35, 0.711, ~1 s | 79/1/42, 0.653, ~1 s | adaptive; +ift05 → 88/1; PERM kept (20,10) |
| T3 | 87/1/34, 0.719 | 80/1/41, 0.661 | 7805 PASS; local null excludes alt reads |
| T4 | 88/2/33, 0.727 | 81/2/40, 0.669 | +12969 TP; +9342 likely-real FP |
| T5 | **92/2/29, 0.760** | **82/2/39, 0.678** | indel preset 0.015; matches plan |

Likely-real non-truth PASS after T5: `15622T>C`, `9342G>A`.

### Permutation parameter sweep (oracle on `base_r10` / `t5_r10.raw_matrix.csv`)
Focused candidates all keep truth 54 / non-truth 0. Post-T5 HapMap A/B: `(10,5)` matches end-to-end TP/FP (r10 92/2, den 82/2) with ~46% fewer floor draws (666k → 359k). **Committed `(10,5)`** in `adaa96a` (was plan default 20/10).

### T2 simulation re-check
Standard (30 seeds) and low-HF (10 seeds) after T1–T5 + perm tune (`adaa96a`). Outdir: `lineage_eval/lineage_sim/results/2026-10-12-t2-ont-r10-vs-denoised/`.

| cell | profile | var_precision | var_recall | ad_f1 | hap_f1 |
|---|---|---|---|---|---|
| standard | ont-r10 | 0.965 | 0.937 | 0.756 | 0.292 |
| standard | ont-denoised | 0.997 | 1.000 | 0.994 | 0.972 |
| low-HF | ont-r10 | 0.970 | 0.947 | 0.795 | 0.169 |
| low-HF | ont-denoised | 0.930 | 1.000 | 0.992 | 0.912 |

**Decision #1 input:** T1 fixed empty VCF on ont-r10, but denoise still wins lineage on default sims. Prefer `ont-r10` for real R10 discovery (HapMap 92 vs 82 TP); keep `ont-denoised` as sim default. Phase 2 lower priority if discovery moves off denoise.

---

## Revision 3 (2026-10-12): read-linkage rescue at full depth

**Starting point:** `hs_lineage_posterio` at `e209313` or later. md5s at revision time: `src/call.rs` `1619cc0833a4ae303c58cf6eced8d947`, `src/main.rs` `4a0b18c3b7f5d0011131568e3386673d`, `src/denoise.rs` `d3b544da05831a5522b05242c0cc0fa3`. If a md5 differs, re-locate line numbers by function name.

### Implementation recommendation (order, effort, expected gain on `merge215k_r10`)

| order | task | effort | expected effect |
|---|---|---|---|
| 1 | **R3-T1** `genotype_from_bam` contig lookup | small (~20 lines + 1 fixture) | The BAM-genotyped lineage matrix stops being empty on GRCh38 BAMs; needed by R3-T3 |
| 2 | **R3-T2** read-linkage rescue (`src/linkage.rs`) | medium (~500 lines + tests) | **+16 TP** (HG00438 ×12, HG02486 491/3537/4793/16239), 0 new FP expected |
| 3 | **R3-T3** BAM strand recount for flagged indels | small | +1 TP (302A>ACC) |
| 4 | **R3-T4** effect-size gate in denoise's strand test | small | `ont-denoised` only: stops erasing 10% alleles at depth |
| 5 | **R3-T5** acceptance and report | runs only | |

Target after R3-T3: **114/121, recall 0.942, FP 2.** Both FPs (9342G>A, 15622T>C) are well-supported and likely missing from the truth set. The other 7 truth variants can't be reached: 3106 (REF N), 14861, 10944, 11087 and 13376 (≤ 23 reads, below the candidate gates or absent), and 309/310 (HG002 alleles in the poly-C run with no linked partner).

Do R3-T2 before considering Phase 2. With linkage rescue on, ont-r10 doesn't need denoise for discovery.

### Where we are (exact scorer; vcfeval agrees)

| run | TP/FP/FN | recall | precision | notes |
|---|---|---|---|---|
| `HapMap_ont.bam` (~420×), ont-r10, after T5 | 92/2/29 | 0.760 | 0.979 | Results table above |
| `merge215k` (55,517 chrM records after the 50k-read cap), ont-denoised | 71/2/50 | 0.587 | 0.973 | All 50 FNs are absent from `raw_matrix`. Denoise erased them (see R3-T4) |
| `merge215k`, **ont-r10** | **97/2/24** | **0.802** | **0.980** | Permutation step 19 s; 597 raw-matrix rows × 48,205 reads |

**The 24 FNs of `merge215k_r10`:**

| truth | n | where lost | R3 fix |
|---|---|---|---|
| HG00438 private (nominal 0.5%; 0.6–1.7% in reads) | 12 | Never become candidates: graph HF below `-v 0.01` | R3-T2 |
| HG02486 491, 3537, 4793 (~1.7–1.8% in graph) | 3 | In `raw_matrix`; rejected by the summed-Jaccard permutation test | R3-T2 |
| HG02486 16239 (~1.1%) | 1 | Not a candidate | R3-T2 |
| 302A>ACC | 1 | `Strand_bias` from graph-assigned reads (77 alt reads) | R3-T3 |
| HG002 309, 310 | 2 | In `raw_matrix`; no linked partner (poly-C) | none |
| 3106 (REF N), 14861, 10944, 11087, 13376 | 5 | No usable reads | none |

### Evidence for the linkage rule

Measured on `ont_benchmark/HG00097_ont.mt.bam`, the same filter output as `merge215k_r10.mt.bam` (55,559 vs 55,517 records, both indexed). Reproduce with `scripts/ont_recall/linkage_oracle.py`.

- **Controls.** 400 random non-truth positions give **760 noise SNV alleles at 0.2–5% HF**, 526 of them strand-skewed. That makes 43,484 noise-involving pairs to test. Under the full rule (overall + per-strand), **0 of 760 noise alleles are rescued**. No truth pair *across* different lineages passes (0 of 389 earlier, 0 here).
- **Truth.** HG00438 **12/12**, HG02486 **4/5** (11087 has 6 reads), HG002 5/8. Those 5 are already called; 309, 310 and 13376 have no partner.
- **Why the per-strand requirement exists.** Without it, 2 noise alleles are rescued: 2542A (129/2) and 2567A (117/3). Both are forward-only, 25 bp apart, at log10 p −10.0 and 14.7× enrichment. The strongest noise pair in which **both** alleles are strand-balanced reaches only log10 p −4.5. Same-lineage truth pairs have a median log10 p of −14.5 (weakest −6.3) and a median enrichment of 52×.
- **Why there is no strand gate in front.** At ~5,000×, 6 true sub-2% alleles are strand-skewed: 3777C 16/69, 6962A 59/14, 10609C 8/22, 12882T 8/31, 16239T 18/38, 10944T 2/21. A strand gate would lose them. The per-strand linkage requirement takes its place.
- **Error-prone reads.** Read error (NM / aligned length) has median 0.013 and 95th percentile 0.129. Dropping the top 5% changed no outcome here. It stays as cheap insurance against errors co-occurring on bad reads, which `call.rs` already warns about.
- **Built-in limit.** Enrichment can't exceed 1/HF of the partner, so the rule only rescues lineages below about 8–10%. Higher ones are already called by the permutation test. For example, 7805A (10%) is not rescued by linkage.
- **Pair-level numbers** reported earlier on the same BAM: 6800–8557 shared 53 (expected 1.7, p 6e-85); 3537–4793 57 (vs 1.8, p 6e-90); 10223–10609 17 (vs 0.19, p 2e-33). Truth–noise pairs: 0 of 349 reach p < 1e-6.

**Strand-gate correction (recorded here, since it was only said in chat).** Keep `STRAND_BIAS_MIN_EFFECT = 0.20`; **don't raise it to 0.30**. At ~5,000×, real alleles reach a skew of 0.19 (8584A 230/95, 15004T 80/161). The flagged non-truth calls, counted from BAM CIGAR windows, are at 0.23 (6653CT>C), 0.28 (6650AG>A), 0.36 (14808TC>T, 14813AC>A), 0.38 (2487AC>A), 0.44 (33CG>C) and 0.53 (6460G>T). In BAM reads 302A>ACC is nearly balanced (57/81 against column 0.525; effect 0.11). The graph's alt-read assignment makes it look skewed (R3-T3).

### R3 decisions (defaults unless the owner edits)

| # | Item | Default |
|---|---|---|
| R3-1 | Source of linkage reads | **Raw, pre-denoise BAM pileup** (`mt_output` in quick-start; `call --linkage-bam <BAM>`). SNVs only. |
| R3-2 | Candidates | SNV alleles with `0.004 ≤ HF < 0.20` and `AC ≥ 10`, on the reads left after R3-3. No strand gate. The partner pool is the same set, whether or not an allele is already PASS. |
| R3-3 | Read filter | Drop reads whose error rate (max over records of NM / aligned length) is above the **95th percentile**, from every count. |
| R3-4 | Pair rule | Distance > 20 bp. Overall: co-covered ≥ 50, shared ≥ 5, enrichment ≥ 10, log10 p ≤ −10. **And each strand separately:** co-covered ≥ 30, shared ≥ 3, enrichment ≥ 10, log10 p ≤ −5. p is the hypergeometric upper tail on co-covered reads. −10 is stricter than Bonferroni 0.01 over 10⁸ pairs. |
| R3-5 | Accept | ≥ 1 passing partner. The partner doesn't itself need to be accepted. |
| R3-6 | Default on? | **On** for `ont-r9`/`ont-r10`/`ont-denoised` in quick-start, with `--no-linkage-rescue` to turn it off. **Off** for `pacbio` until evaluated (the PacBio HapMap run is 111/12/10). |
| R3-7 | VCF | New records are PASS SNVs with `INFO DP=<pileup depth>;LINK=<partner m-name>,<log10p>`, plus a new `##INFO=<ID=LINK,...>` header. `AD`/`HF` come from pileup counts. If the same allele already exists as a graph record: skip it if PASS; if filtered, clear its FILTER (SNVs only). |
| R3-8 | Lineage matrix | Rescued variants are appended to the list passed to `write_final_matrix`. They get rows only via `--genotype-from-bam`; the graph matrix is unchanged. |

All R3 constants go in `main.rs` as `LINK_*`, with the evidence numbers above in their doc comments.

### Task R3-T0: Preconditions

- [ ] Commit the agent's uncommitted `max-partner` prototype in `scripts/ont_recall/perm_retest.py` (T8 prototype) and the new `scripts/ont_recall/linkage_oracle.py`: `scripts: linkage oracle and max-partner prototype (ONT-recall R3-T0)`.
- [ ] `cargo test`: record the count (406 expected after T5).
- [ ] Baseline: `bash benchmarking.sh ont-r10 r3_base` on `HapMap_ont_merge.bam` should give 97/2/24.
- [ ] Oracle on the run's own BAM: `samtools index r3_base.mt.bam && python3 ../scripts/ont_recall/linkage_oracle.py --bam r3_base.mt.bam --truth $TRUTH --random 400` (about 6 min). Expect `noise alleles rescued: 0 of ~760`, HG00438 12/12, HG02486 4/4, with 10944, 11087 and 13376 below the candidate gates. Keep the output; R3-T2 parity is checked against it.

### Task R3-T1: `genotype_from_bam` must find the contig by name

**Bug:** `src/call.rs` ~1086, `if rec.is_unmapped() || rec.is_secondary() || rec.seq_len() == 0 || rec.tid() != 0 { continue; }`. In a GRCh38 BAM chrM is tid 24, so every read is skipped. Every HapMap run so far logs "Genotyped N variants x 1 reads", and `matrix.csv` holds only the `chrM` pseudo-read column. The sims never caught it because their BAMs have chrM as the only contig.

- [ ] **Step 1 (failing test).** Add a fixture `test_data/genotype/reads_second_contig.sam`: the same reads as `reads.sam`, with a header of `@SQ SN:chr1 LN:100` followed by `@SQ SN:chrM LN:60`, and the reads placed on `chrM` (tid 1). Assert that `genotype_from_bam(.., ref_name = "chrM", ..)` returns the same `read_record`/`cover_record` as for `reads.sam`. Also assert that an unknown `ref_name` returns maps holding only the reference pseudo-read, with one `eprintln!` and no exit.
- [ ] **Step 2.** Resolve `let tid = reader.header().tid(ref_name.as_bytes())` once, then skip records with `rec.tid() != tid as i32`. Grep for any other `tid() != 0` / `tid() == 0` assumption in `call.rs` and `build.rs`, and list what you find in the commit message.
- [ ] **Step 3 (acceptance).** On `merge215k`, the log should show "Genotyped ~97 variants x >1,000 reads", and `matrix.csv` should have more than 1,000 read columns. Commit: `call: genotype_from_bam resolves the contig by name (ONT-recall R3-T1)`.

### Task R3-T2: Read-linkage rescue of low-frequency SNVs

**New file:** `src/linkage.rs` (add `mod linkage;` in `main.rs`). **Oracle:** `scripts/ont_recall/linkage_oracle.py` (same rule, BAM pileup).

**Algorithm**
1. **Reads.** One pass over mapped, non-secondary records on the contig, with tid looked up by name as in R3-T1. Per qname keep: an id (`u32`); `reverse` from the primary record; `err = max(NM / query_alignment_length)` over the qname's records (if NM is absent, use 0.0 and `eprintln!` once); and the span list from `aligned_blocks()`, merged per qname. Compute the 95th-percentile `err` and mark reads above it as `bad`.
2. **Candidates.** One pileup pass, following the `reader.pileup()` loop in `denoise.rs::compute_corrections` (~422), with `max_depth` set to `u32::MAX`. Per column, over non-`bad` reads without deletion or refskip, count each base. Every non-reference base with `AC ≥ LINK_MIN_AC` and `LINK_MIN_HF ≤ AC/depth < LINK_MAX_HF` becomes a `Candidate { pos, ref, alt, ac, dp, fwd, rev, reads: Vec<u32> }`, with reads sorted.
3. **Pairs.** Build an inverted index from read to the candidates it carries. For every read, count every candidate pair more than `LINK_MIN_DIST` apart into `HashMap<(u32,u32), u32>`. Keep the pairs with shared ≥ `LINK_MIN_SHARED`.
4. **Co-covered counts.** "Co-covered" means the read's merged span covers both positions; the oracle uses pileup coverage, so allow small differences. `N` (overall, forward, reverse) comes from a Fenwick tree: sort the queries by `min(p1,p2)`, insert reads by span start, query span end ≥ `max(p1,p2)`. Keep separate trees for the two strands, so `N_fwd + N_rev = N`. `K` and `n` are the alt reads of each allele whose span covers the other position. `k` is the shared count.
5. **Test.** `log10_hyper_sf(k, N, K, n)` is the log-space sum of `ln_gamma` terms (`statrs::function::gamma::ln_gamma`, already used in `scite.rs`). Enrichment = `k·N/(K·n)`. Apply the overall rule and then the per-strand rule (R3-4).
6. **Accept** a candidate with ≥ 1 passing partner that is not already a PASS call (R3-5). Return `LinkCall { pos, ref, alt, ac, dp, partner_name, shared, enrich, log10p }`.

**Wiring:** `call::start` takes `linkage_bam: Option<&PathBuf>` and `link_opts`. After `permutation_test`, run the rescue and merge it into `write_vcf` per R3-7. This needs a small `format_link_record`, and the VCF sort becomes "by POS across both lists". Then append the variants for `write_final_matrix` (R3-8). In quick-start, add `--no-linkage-rescue`; pass `Some(&mt_output)` when the data type is ONT and the flag is absent (R3-6). In `call`, add `--linkage-bam`. Print one summary line: candidates, pairs tested, pairs passing, accepted, wall time.

- [ ] **Step 1 (failing tests).** Build fixtures programmatically: a `bam::Writer` into a temp dir, with header `chr1`, `chrM`, so the tid lookup is exercised too.
  - `log10_hyper_sf` test vectors, to ±0.01: `(53, 2635, 68, 65) → −84.22`, `(17, 4040, 30, 25) → −32.64`, `(11, 4000, 129, 117) → −2.92`, `(2, 1000, 20, 20) → −1.24`, `(1, 1000, 20, 20) → −0.48`.
  - The Fenwick co-cover count equals a brute-force count on 2,000 random spans × 500 random position pairs.
  - **Linked on both strands → accepted:** two alleles 500 bp apart, carried together by 20 forward and 20 reverse reads out of 2,000 spanning reads.
  - **Linked on one strand only → rejected:** same as above, but all 40 shared reads are forward (the 2542A/2567A case).
  - **Independent alleles → rejected:** 40 random carriers each, out of 2,000 reads.
  - **Linked only on bad reads → rejected:** the 40 shared carriers are the highest-error reads (NM/len 0.2), while the rest are at 0.01.
  - **Distance ≤ 20 bp → never tested.**
  - **Already-PASS allele → not emitted again**, and a filtered graph record of the same SNV gets its FILTER cleared.
- [ ] **Step 2 (implement).**
- [ ] **Step 3 (acceptance).**
  - `merge215k` ont-r10 (`bash benchmarking.sh ont-r10 r3_t2`): **TP ≥ 112** (expected 113), FP ≤ 3. List every non-truth `LINK` record with its partner and log10 p, so the owner can review it.
  - **Oracle parity:** the accepted truth alleles match R3-T0's oracle output, ±1 allele.
  - Runtime increase ≤ 60 s.
  - `HapMap_ont.bam` (~420×): TP ≥ 92, with no new non-truth PASS. Little gain is expected there, because most HG00438 alleles have fewer than 10 reads.
  - `ont-denoised` on `merge215k`: report it; linkage reads come from the raw BAM, so expect gains there too.
  - Commit: `call: read-linkage rescue of low-frequency SNVs (ONT-recall R3-T2)`.

### Task R3-T3: Recount strand bias from BAM reads for flagged small indels

**Why:** graph edges assign 302A>ACC's alt reads unevenly by strand, so it is flagged `Strand_bias`. In BAM reads it is 57/81 against column 0.525 (effect 0.11). **Depends on** R3-T1 (contig lookup) and R3-T2 (the raw BAM is threaded into `call::start`).

- [ ] **Step 1 (failing test).** An indel flagged from graph counts (10/0), whose BAM-genotyped alt reads are 9/11 against a balanced column, gets its flag cleared. One that is 1/19 in the BAM stays flagged.
- [ ] **Step 2.** After `filter_strand_bias`, collect the **indels** with `filter == Some("Strand_bias")`. Run `genotype_from_bam` on just those, against the raw BAM. Recount fwd/rev with `strand_map`, using a local null over the BAM cover set. Recompute p and effect, and clear the flag if `p ≥ strand_bias_threshold || effect < STRAND_BIAS_MIN_EFFECT`. Log the graph and BAM counts side by side. SNVs are unchanged.
- [ ] **Step 3 (acceptance)** on `merge215k`:
  - 302A>ACC is PASS (+1 TP).
  - These stay `Strand_bias`: 33CG>C, 2487AC>A, 6650AG>A, 6653CT>C, 14808TC>T, 14813AC>A. Their BAM CIGAR-window effects are 0.44, 0.38, 0.28, 0.23, 0.36 and 0.36. Realignment counts may differ slightly; report them.
  - FP unchanged.
  - Commit: `call: recount strand bias for flagged indels from BAM reads (ONT-recall R3-T3)`.

### Task R3-T4: Effect-size gate in denoise's strand-bias test

**Why:** at ~5,000×, `denoise.rs::fit_site` (~141–146) and `denoise_indel.rs::fit_indel_site` (~518–525) treat `p < DENOISE_STRAND_BIAS_P` alone as an artifact. In `merge215k` ont-denoised this erased real alleles:

| allele | strands | p | effect | reads raw → after denoise |
|---|---|---|---|---|
| 7805A (10%) | 328/196 | 3e-7 | 0.11 | 524 → 12 |
| 12372A | — | 4e-3 | 0.07 | 402 → 42 |
| 14179G | — | 7e-3 | 0.07 | 383 → 64 |
| 8584A | 230/95 | — | 0.19 | 325 → 8 |
| 15004T | 80/161 | — | 0.19 | 241 → 35 |

- [ ] **Step 1 (failing tests).** In `fit_site`, an allele at 328/196 in a column with forward fraction 0.514 stays eligible. An allele at 449/15 (the 1598G>A artifact, effect 0.50) is still rejected. Add the same pair of tests for `fit_indel_site`.
- [ ] **Step 2.** Add `&& (fwd as f64 / (fwd + rev) as f64 - expected_fwd_frac).abs() >= crate::STRAND_BIAS_MIN_EFFECT` to both strand-bias conditions. Use the same constant as `call`, so they can't drift apart.
- [ ] **Step 3 (acceptance).**
  - `merge215k` ont-denoised: `fn_pileup_audit.py --raw r3_t4.mt.bam --denoised r3_t4.mt.denoised.bam` shows the denoised/raw alt-count ratio ≥ 0.9 for the five alleles above. TP should rise from 71, with FP unchanged. Sub-3% alleles are still erased by `DENOISE_KEEP_VAF` (Phase 2 T6).
  - `lineage_sim bench_standard.sh --profile ont-denoised` (30 seeds): no drop in `var_precision`/`ad_f1`. The sims' errors carry no strand bias, so expect no change.
  - Commit: `denoise: effect-size gate on the strand-bias test (ONT-recall R3-T4)`.

### Task R3-T5: Acceptance report

- [ ] Fill in the R3 results table below, and the cargo test count.
- [ ] Add one line to the `lineage_sim/README.md` section on choosing `ont-r10` vs `ont-denoised`: with linkage rescue on, ont-r10 is the recommended discovery profile for real ONT data, and ont-denoised remains the sim default for lineage metrics (T2 numbers).

### R3: Needs a decision

1. **Default linkage rescue on for ONT** (R3-6). The evidence is one deep sample. Confirm on one HPRC ONT sample, where heteroplasmy is rare, that it adds ≤ 1 call per sample.
2. **`--maximal-mt-depth` (50,000).** The merged BAM has 215k chrM records. At 50k reads the linkage evidence is already overwhelming (truth pairs median p 3e-15). Raising the cap costs graph size (179 MB GFA) and runtime for little gain. Keep it.
3. **Truth-set gaps.** Admitting candidates down to 0.4% can surface real cell-line sub-lineage variants the assembly truth lacks, as 9342G>A and 15622T>C already are. Report linked non-truth calls separately; don't tune against them.
4. **Linkage for low-HF homopolymer indels (follow-up).** 285CAA>C is linked to HG02257 (p 4e-174) and 302A>ACC to HG02257 (p 7e-69). The earlier 10916TC>T FP (1,100×) was unlinked (best p 0.015). Indel linkage could replace the HF floor (`indel_false_threshold`). It isn't needed at full depth today, so it isn't scheduled.

### R3 results

cargo test after R3-T4: **427** passed / 1 ignored

| step | merge215k ont-r10 (TP/FP/FN) | merge215k ont-denoised | HapMap_ont ~420× ont-r10 | notes |
|---|---|---|---|---|
| R3 baseline | 97/2/24 | 71/2/50 | 92/2/29 | plan table (pre-R3) |
| R3-T1 | genotype fix | | | `r3_final`: Genotyped 127 × **45,545** reads (was 1) |
| R3-T2+T3 (`r3_final`) | **114/6/7, recall 0.942** | — | ≥92 expected | 302A>ACC PASS; HG00438 12/12 |
| R3-T4 (`r3_t4_den`) | — | **109/4/12, recall 0.901** | — | was 71/2/50; linkage on raw BAM + effect gate |

Non-truth PASS on merge215k ont-r10: `9342G>A`, `15622T>C` (graph), `10092C>T`↔`12731T>C` (LINK), `302A>AC`, `567AC>A`. Truth FNs: 309, 310, 3106, 10944, 11087, 13376, 14861 (plan unreachable set).

`fn_pileup_audit` after T4: 7805/12372/14179 preserved (ratio ≈1); 8584 and 15004 (effect ≈0.19, just under 0.20) still largely erased by denoise (keep-VAF / borderline).
