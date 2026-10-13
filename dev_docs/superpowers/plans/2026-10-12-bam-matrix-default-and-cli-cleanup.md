# BAM-Genotyped Matrix by Default, CLI Cleanup, Short Comments — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Do the tasks **in order**, one commit per task (C4 is one commit per file). Each task starts from a green `cargo test`.

**Goal (owner request, 2026-10-12):**
1. The variant × read matrix is always re-genotyped from the raw BAM. Drop `--genotype-from-bam` and the `graph_matrix.csv` output. Build the noise gate and the higher SCITE error rates into the program.
2. Cut the `quick-start`, `call` and `lineage` options down. Each option is either **visible** (in `-h` and `--help`), **advanced** (only in `--help`), or **hard-coded** (flag removed, value becomes a constant).
3. Make every comment in `src/*.rs` short.

**Starting point:** branch `hs_lineage_posterio` at `e852b1b` or later (ONT-recall R3 + Phase 2 done). md5s: `src/main.rs` `12a7cf3097dee5870a7fa6d6af2686c3`, `src/call.rs` `d8191a38cc3770033b68f2827a88cc40`, `src/lineage.rs` `d40bae082828501a6d4d784761841447`. If one differs, re-locate code by function name.

**Tools (in the repo, uncommitted; C0 commits them):**
- `scripts/ont_recall/matrix_cell_eval.py` scores matrices cell by cell on the HapMap mixture.
- `scripts/dev/check_comment_only_diff.py` fails if a diff changes anything but comments. It was self-tested on all 17 `src/*.rs` files.

---

## Why: the matrix evidence (HapMap `merge215k`, ont-r10, 2026-10-12)

`matrix_cell_eval.py`: 275,438 cells over 96 truth variants. Reads are assigned to a cell line from private alleles at other sites. The BAM matrix is the Python reference `bam_genotype_matrix.py`, which matched the Rust code 100% on the sim.

| class | matrix | correct 1s | FN | FP | missing |
|---|---|---|---|---|---|
| major SNV | graph | 69,549 | 0.37% | 0.21% | **18.6%** |
| major SNV | BAM, max-ed 4 | 83,024 | 0.26% | 0.62% | **2.9%** |
| minor SNV | graph | 6,257 | 1.32% | 0.07% | **15.9%** |
| minor SNV | BAM, no gate | 8,291 | 1.10% | 0.57% | 1.6% |
| minor SNV | BAM, **max-ed 4** | 8,176 | 0.74% | **0.34%** | **3.1%** |
| minor indel | graph | 481 | **31.4%** | 0.39% | 15.7% |
| minor indel | BAM, max-ed 4 | 658 | 11.0% | 1.20% | 4.3% |

Max-ed 6 sits between those (minor SNV FP 0.45%, missing 2.1%). Max-ed 2 cuts FP to 0.19% but leaves 6% missing. Decisive cells agree 99.9% between the two matrices. On the sim (plan 2026-10-11), the denoised graph matrix had 22% false 0s; the raw-BAM matrix had 2.2%.

The matrix now comes from **raw reads in every ONT mode**. So the `ont-denoised` SCITE preset (fp 0.0001, fn 0.01), which was tuned for the denoised graph matrix, no longer matches its input. Its measured FP is 0.3–0.6% on SNVs.

---

## Decisions (defaults unless the owner edits)

| # | Item | Default |
|---|---|---|
| D1 | Matrix source | Always `genotype_from_bam` on the **raw, pre-denoise** BAM. `quick-start` uses `<prefix>.mt.bam`. **No `graph_matrix.csv`** is written. |
| D2 | `call` BAM input | `--input-bam`, `--genotype-bam` and `--linkage-bam` become **one** `-b/--bam <raw pre-denoise BAM>`, used for the strand map, linkage rescue (ONT) and the matrix. `alias = "input-bam"` keeps old scripts parsing. Without `--bam`: no strand filter, no linkage, and `matrix.csv` comes from the graph with a one-line warning, so old WDLs keep running. |
| D3 | Noise gate | `GENOTYPE_MAX_ED_ONT = 4` for `ont-r9`/`ont-r10`/`ont-denoised`; none for `pacbio`. Resolved by a new `call::resolve_genotype_opts(data_type)`. `GENOTYPE_FLANK = 10` is unchanged. |
| D4 | SCITE presets | `ont-r9`/`ont-r10`/`ont-denoised`: **fp 0.005, fn 0.02** (one ONT preset, since the matrix source is the same). `pacbio` stays 0.005/0.05. HP-tract rates stay at 2× the resolved rates. C2's sweep can override these. |
| D5 | Option classes | Exactly the tables in C3. Advanced = `#[arg(hide_short_help = true)]` (clap 4.5). Hard-coded = flag removed, constant in `main.rs`. |
| D6 | Removed experiment | `--cooccurrence-stat` and the `max-partner` code go, along with their tests. Phase 2 T8 rejected it because it lets noise through on deep ONT data. The Python prototype stays in `scripts/ont_recall/perm_retest.py`. |
| D7 | `--permutation-mode` | The flag is removed. `PermutationMode::Fixed` stays as an internal code path, because the T1 characterization tests use it. |
| D8 | Comment scope | `src/*.rs` only, plus clap help text. Not the Python/shell scripts or the plan docs. |

---

## Global constraints

- Write a failing test first, implement, then get `cargo test` green for the whole crate.
- Don't delete or weaken a test unless this plan says so (C1 and D6 name every intended change). When an assertion changes because of an intended behaviour change, say why in the commit message.
- No new flags, except `-b/--bam` on `call`.
- New and edited comments follow the C4 rules from the start.
- After C1, C2 and C3, re-run the HapMap harness (`ont_benchmark/benchmarking.sh ont-r10 <prefix>` on `HapMap_ont_merge.bam`) and confirm the **VCF is byte-identical** to the C0 baseline. These tasks must not change calling.

---

### Task C0: Preconditions and baseline

- [ ] `.gitignore` has `scripts/*` with only `scripts/ont_recall/` re-included, so add `!scripts/dev/` and `!scripts/dev/**`. Delete `scripts/dev/__pycache__/`, which the reviewer's self-test left behind. Then commit `scripts/ont_recall/matrix_cell_eval.py`, `scripts/dev/check_comment_only_diff.py` and this plan (`git add -f`, since `superpowers` is ignored, the same as the other plans): `scripts: matrix cell eval and comment-only diff check (cleanup C0)`.
- [ ] `cargo test`: record the count.
- [ ] Baseline on `HapMap_ont_merge.bam` with `bash benchmarking.sh ont-r10 c0_r10`. Expect R3 `r3_final`: 114/6/7. Keep `c0_r10.vcf` for the byte-identity checks.
- [ ] Save the help text from before the change: `for c in quick-start call lineage; do Himito $c -h > docs/cli/$c.before.short.txt; Himito $c --help > docs/cli/$c.before.long.txt; done`. These files don't get committed until C3.
- [ ] Sim baseline: `lineage_eval/lineage_sim/bench_standard.sh` (30 seeds) for `--profile ont-r10` and `--profile ont-denoised`. Record `var_*`, `ad_f1` and `hap_f1`.

### Task C1: BAM-genotyped matrix is the only matrix

**Files:** `src/call.rs`: `write_final_matrix` (~2519), `start` (the `bam_file` / `genotype_bam` / `linkage_bam` parameters, ~2339–2356), `GenotypeOpts` (~1002). `src/main.rs`: the `QuickStart` arm (~1040–1056), the `Call` arm (~1160–1206), `quickstart_linkage_bam` (~914).

- [ ] **Step 1 (failing tests).**
  - `resolve_genotype_opts`: `ont-r10`, `ont-r9` and `ont-denoised` give `max_ed == 4`; `pacbio` gives `u32::MAX`; `flank == 10`, `margin == 1`.
  - Rewrite `write_final_matrix_with_genotype_bam_swaps_matrix_and_keeps_graph_copy` as `write_final_matrix_with_bam_writes_no_graph_copy`: `matrix.csv` holds the BAM matrix and `graph_matrix.csv` does **not** exist. This is an intended change; say so in the commit.
  - Keep `write_final_matrix_without_genotype_bam_is_unchanged`. It is now the D2 fallback, and must also check that the warning is printed. Capture it through a returned flag or a small `fn` that formats it.
  - CLI parse tests: `quick-start` no longer accepts `--genotype-from-bam`. `call -b raw.bam` and `call --input-bam raw.bam` both set `bam`. `call --genotype-bam` and `call --linkage-bam` are rejected.
- [ ] **Step 2 (implement).**
  - `call::start` takes `bam: Option<&Path>` in place of `bam_file`, `genotype_bam` and `linkage_bam`. Use it for `build_strand_map`, for `linkage::rescue` when the data type is ONT and linkage isn't disabled, for the R3-T3 indel strand recount, and for `write_final_matrix`.
  - `quick-start` passes `Some(&mt_output)`, never `build_bam`. Read names are the same in the raw and denoised BAMs, so the strand map doesn't change; check this with the byte-identity test.
  - Delete `quickstart_linkage_bam`; its logic moves into `start` with the data-type check. Keep its test assertions as tests of the new check.
  - Write `GENOTYPE_MAX_ED_ONT` and `GENOTYPE_FLANK` into `main.rs`, each with a one-line comment.
- [ ] **Step 3 (scripts).**
  - `lineage_eval/lineage_sim/run_himito.sh`: replace `--input-bam "$BUILD_BAM"` and the `GENOTYPE_BAM_MODE` case with `--bam "$BAM"` (the raw BAM). Remove the `--genotype-bam raw|denoised|none` option from `run_himito.sh`, `run_eval.sh` and the README.
  - `scripts/ont_recall/run_hapmap_bench.sh`: drop `--genotype-from-bam`.
  - `lineage_eval/lineage_sim/bam_genotype_matrix.py`: fetch the contig by name (`--contig`, default `chrM`) in place of `bam.references[0]`. It has the same bug R3-T1 fixed in Rust.
  - Tell the owner to drop `--genotype-from-bam` from `ont_benchmark/benchmarking.sh` (untracked).
- [ ] **Step 4 (acceptance).**
  - The HapMap VCF is byte-identical to `c0_r10.vcf`, and `graph_matrix.csv` is absent.
  - `matrix_cell_eval.py --bam c1_r10.mt.bam(.bai) --truth $T new=c1_r10.matrix.csv`: minor SNV FP ≤ 0.40%, missing ≤ 4%; major SNV missing ≤ 4%.
  - Sim `bench_standard.sh` (30 seeds, both profiles), with C0's presets still in place: report `ad_f1`/`hap_f1` against the C0 baseline. **Watch-out:** on the sim, max-ed 4 raised missing from 12% to 28% (plan 2026-10-11 table), because pbsim reads are noisier than real R10 (HapMap median read error 1.3%). If `hap_f1` or `ad_f1` drops by more than one seed-SD, try `GENOTYPE_MAX_ED_ONT = 6`, report both, and stop for the owner.
  - Commit: `call: BAM-genotyped matrix always, one --bam input, ONT noise gate (cleanup C1)`.

### Task C2: SCITE error-rate presets for the BAM matrix

**File:** `src/lineage.rs` `resolve_error_rates` (~1095) and its test `resolve_error_rates_presets_differ_per_data_type`.

- [ ] **Step 1 (sweep, report).** Run `lineage_sim/sweep_fpfn.sh` on the `bench_standard` config (C1 build, ≥ 10 seeds) for `ont-r10` and `ont-denoised`, over fp ∈ {0.001, 0.003, 0.005, 0.01} × fn ∈ {0.01, 0.02, 0.05}. Rank by `ad_f1`, then `hap_f1`. If (0.005, 0.02) is within one seed-SD of the best, use it (D4). Otherwise use the best pair, and record why.
- [ ] **Step 2 (failing test, then implement).** Update the preset test to the chosen values. All three ONT types must now return the **same** tuple; this is intended (D4). Keep a separate assertion for `pacbio`.
- [ ] **Step 3 (acceptance).** Re-run `bench_standard.sh` (30 seeds, both profiles). `hap_f1`/`ad_f1` must be ≥ the C1 numbers. HapMap VCF byte-identical. Commit: `lineage: ONT SCITE presets for the raw-BAM matrix (cleanup C2)`.

### Task C3: Option cleanup for `quick-start`, `call`, `lineage`

**Rules.** **Visible** options show in `-h` and `--help`. **Advanced** options get `hide_short_help = true`, so they show only in `--help`, and keep their names and short flags because scripts and WDLs use them. **Hard-coded** options lose their flag; the value becomes a `main.rs` constant with a one-line comment. Every help string is one line, ≤ 80 characters, with no defaults in the text (clap prints those). Add `value_parser = PossibleValuesParser::new(DATA_TYPES)` to `-d` on `quick-start` and `call` (`lineage` already has it).

**`quick-start`**

| option | class | why / where used |
|---|---|---|
| `-i` `-c` `-r` `-s` `-d` `-o`, `--maximal-mt-depth`, global `--threads` | visible | core inputs; `--maximal-mt-depth` used by README and `Himito_lineage.wdl` |
| `-k`, `-m`, `-v`, `-t`, `-f`, `--p-value-threshold`, `--heteroplasmic-frequency-threshold`, `--strand-bias-threshold`, `--indel-false-threshold`, `--permutation-frequency-threshold`, `--permutation-rounds`, `--no-left-align`, `--no-linkage-rescue`, `--denoise-keep-vaf` | advanced | used by WDLs, README, User_guide or the HapMap harness, or useful as A/B switches |
| `-l` (3000), `--downsample-seed` (42), `--strand-bias-drop-max-hf` (0.10) | hard-coded | never set anywhere |
| `--strand-bias-action` | hard-coded | data-type preset (`resolve_strand_bias_action(dt, None)`) |
| `--permutation-mode` | hard-coded | adaptive (D7) |
| `--cooccurrence-stat` | removed | D6 |
| `--genotype-from-bam` | removed | always on (C1) |
| `--left-align` | keep `hide = true` | no-op kept for compatibility; `Himito_lineage.wdl` passes it |

**`call`**

| option | class | why / where used |
|---|---|---|
| `-g` `-r` `-s` `-d` `-o`, `-b/--bam` | visible | core inputs (C1) |
| `-k`, `-m`, `-v`, `-p`, `-f`, `--strand-bias-threshold`, `--indel-false-threshold`, `--permutation-frequency-threshold`, `--permutation-rounds`, `--no-left-align` | advanced | `run_himito.sh`, WDLs |
| `--strand-bias-action`, `--strand-bias-drop-max-hf`, `--permutation-mode`, `--genotype-max-ed`, `--genotype-flank` | hard-coded | presets or constants (C1, D3, D7) |
| `--cooccurrence-stat`, `--genotype-bam`, `--linkage-bam` | removed | D6, D2 (`--input-bam` survives as an alias of `--bam`) |
| `--left-align` | keep `hide = true` | compatibility |

**`lineage`**

| option | class | why / where used |
|---|---|---|
| `-m` `-v` `-d` `-o`, `--min-hf`, `--max-hf` | visible | core; used across the sims and README |
| `--fp-rate`, `--fn-rate`, `--hp-fp-rate`, `--hp-fn-rate`, `--collapse-hp`, `--mcmc-iterations`, `--mcmc-chains`, `--debug` | advanced | sweeps, the hp-sweep WDL, README |
| `--no-root-block`, `--root-block-min-hf` (0.80), `--root-block-max-q` (0.05), `--root-block-min-absent` (10), `--root-block-min-or` (2.0) | hard-coded | never set anywhere; root block stays on |

- [ ] **Step 1 (failing tests).**
  - Render `Cli::command()` help for each subcommand with clap's `render_help` / `render_long_help`. Assert that the short help contains every **visible** long flag and none of the **advanced** ones, that the long help contains both, and that no **hard-coded or removed** flag appears in either.
  - Parse tests: every advanced flag still parses, removed flags are rejected, and the `--input-bam` alias works.
- [ ] **Step 2 (implement).** Remove the hard-coded and removed fields from the clap structs and their match arms, and pass the constants. Remove the `max-partner` code and its tests (D6). Fix `validate_lineage_opts` and any other validator that only checked a removed flag.
- [ ] **Step 3 (docs).** Update `README.md` and `docs/User_guide.md` for removed and renamed flags (`--genotype-from-bam`, `--genotype-bam`, `--linkage-bam`, `--input-bam` → `--bam`), and note that advanced options are in `--help`. Write `docs/cli/<cmd>.short.txt` and `docs/cli/<cmd>.long.txt` from the new binary and commit them, so CLI changes show up in review. Delete the C0 `*.before.*` files.
- [ ] **Step 4 (acceptance).** HapMap VCF byte-identical. `run_himito.sh` single-seed smoke for both profiles gives metrics identical to C2's run for the same seed. `grep` across `wdl/`, `lineage_eval/`, `scripts/`, `README.md` and `docs/` finds no hard-coded or removed flag. Commit: `cli: visible/advanced/hard-coded option split for quick-start, call, lineage (cleanup C3)`.

### Task C4: Short comments in `src/*.rs`

**Rules**
1. **Doc comment on an item:** one summary line (≤ 100 characters). Add at most 2 more lines only for a non-obvious invariant, unit, or reason.
2. **No history:** delete "Was `--x`", "Plan 2026-…", "ONT-recall …", "Review 2026-…", "previously", "an earlier version". History lives in git and `dev_docs/`.
3. **No evidence numbers** (benchmarks, sweeps, read counts). The constants block in `main.rs` gets **one** pointer line: `// Thresholds: see dev_docs/superpowers/plans/2026-10-11-ont-recall-high-coverage.md`.
4. **Inline comments** only where the code isn't self-explanatory, one line each. Delete comments that restate the code.
5. **Tests:** one line saying what the test pins, if the name doesn't already say it.
6. **clap help:** one line, ≤ 80 characters (already done for the three subcommands in C3; apply to the rest).
7. Keep `SAFETY:` notes and license headers, shortened.

- [ ] **Per file, largest first:** `scite.rs`, `call.rs`, `denoise.rs`, `denoise_indel.rs`, `lineage.rs`, `main.rs`, `build.rs`, `callnumts.rs`, then the rest.
  - Edit comments only.
  - Run `python3 scripts/dev/check_comment_only_diff.py HEAD`. It must print `comment-only` for the file.
  - Run `cargo build --release` (no new warnings) and `cargo test` (same count).
  - Commit `docs(src): shorten comments in <file> (cleanup C4)`.
- [ ] **Acceptance:**
  - `grep -c "^\s*//" src/*.rs` total goes from **5,124** (at `e852b1b`) to **≤ 2,500**.
  - `grep -rn "Plan 2026\|ONT-recall\|Review 2026\|Was \`--" src/` returns nothing.
  - `check_comment_only_diff.py <C3 commit>` over the whole C4 range prints `comment-only` for every file.

### Task C5: Final verification and report

- [ ] `cargo test` count (C3 changes it by removed and added tests; list the deltas).
- [ ] HapMap `merge215k` ont-r10 and ont-denoised: TP/FP/FN (should equal C0), plus the `matrix_cell_eval.py` table for the final `matrix.csv`.
- [ ] `bench_standard.sh` (30 seeds, both profiles) against C0.
- [ ] Fill in **Results**.

---

## Needs a decision (don't implement without sign-off)

1. **`call` without `--bam`** (D2): the graph-matrix fallback keeps old WDLs working. The alternative is to make `--bam` required and update `wdl/Himito_call.wdl`, `Downsample_experiment.wdl`, `Himito_methyl.wdl` and `AssemblyBenchmark_sr_refine.wdl`.
2. **Noise gate value** if the sims disagree with HapMap: 4 (HapMap-optimal) or 6 (safer on noisier reads). Decide on C1's report.
3. **PacBio gate and presets:** not measured here. The PacBio HapMap VCF is 111/12/10; a cell-level eval on a PacBio BAM would decide whether `max-ed` should apply there too.

---

## Results

### cargo test
- C0: ____ / C5: ____

### HapMap `merge215k`

| step | ont-r10 TP/FP/FN | VCF identical to C0? | matrix: minor SNV FP / missing | notes |
|---|---|---|---|---|
| C0 | | — | graph: 0.07% / 15.9% | |
| C1 | | | | |
| C2 | | | | presets: |
| C3 | | | | |
| C5 | | | | |

### `bench_standard.sh` (30 seeds)

| step | profile | var_f1 | ad_f1 | hap_f1 |
|---|---|---|---|---|
| C0 | ont-r10 | | | |
| C0 | ont-denoised | | | |
| C1 | | | | |
| C2 | | | | |
