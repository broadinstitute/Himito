# Simulator: Homopolymer Length Heteroplasmy in rCRS — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Do the tasks **in order**, one commit per task (or per numbered step where a task says so). Each task starts from green `pytest lineage_eval/lineage_sim/tests` and green `cargo test`.

**Goal:** Make `lineage_eval/lineage_sim/` able to simulate the homopolymer (HP) length variation that dominates real PacBio/ONT disagreement, so that HP-related changes in `Himito lineage` (T5/T6/T8/T9 of the 2026-09-26 plan, and the leaf-only / categorical-tract work proposed on 2026-09-27) can be measured instead of guessed.

**Starting point:** branch `hs_lineage_posterio` at `1fbffdd` or later. `simulate_tree.py` is SNV-only; `score_lineage.py` keys everything on exact `m.<pos><ref>><alt>` strings.

**Tech stack:** Python 3.11 (`simulate_tree.py`, `score_lineage.py`, pytest in `tests/`), Bash harness (`run_eval.sh`, `simulate_reads.sh`, `run_himito.sh`, `sweep_seeds.sh`, `bench_standard.sh`), pbsim3, minimap2, samtools, WDL 1.0 (`wdl/Himito-lineage/Himito_lineage_sim.wdl` and its Docker image from `lineage_eval/lineage_sim/Dockerfile`; validate with `womtool` or `miniwdl`). Rust only where a task says so.

**Task map:** T1–T6 simulator and scorer; T7 error-model measurement; T8 frozen local benchmark and pair mode; T9 first measurements; T10 revise the existing simulation WDL; T11 new WDL that sweeps `--hp-fp-rate` × `--hp-fn-rate` and reports the best pair.

## Review of the current simulator (2026-09-27)

What it does well: the truth/clone/read pipeline is deterministic per seed, circularity is handled by two rotations, the three frequency gates are documented, `bench_standard.sh` freezes a comparable configuration, and `score_lineage.py` has already fixed the two scoring traps (ROOT in ancestor sets, `ad_recall` over shared pairs only).

What blocks HP work:

1. **The benchmark is blind to the thing that diverges.** `AVOID_RANGES` in `simulate_tree.py` excludes 295–320 and 16180–16200 by design, and the truth contains no indels at all. Every HP tract that dominates the real `HG00xxx` PacBio/ONT trees (D310, HVS1 poly-C, 573, 455) is either excluded or never mutated. The T7 sweep of `--hp-fn-rate` over 1×/2×/4× produced identical `ad_f1` because no truth variant was ever affected by the parameter. Any HP conclusion drawn from `bench_standard.sh` is vacuous.

2. **Two different biological things are both called "HP indel", and the simulator must produce both.** (a) A *lineage-borne* length change: one clone gained a C in a run, and the allele is inherited — an ordinary ISA edge that happens to be hard to observe. (b) *Lineage-independent* length heteroplasmy: D310 / 16184 / 573 change length within an individual at a rate that makes them recur across clones, so the tract length is a per-molecule state, not a tree edge. The real trees show (b): `m.302A>AC` vs `m.310T>TC` in HG00290 PacBio has 55 child-without-parent and 138 parent-without-child reads out of ~800, far beyond fn = 0.05 dropout. A simulator with only (a) would still reward putting HP alleles on the backbone.

3. **`clone_genome` and `write_truth` cannot represent indels.** `clone_genome` asserts `g[pos-1] == ref` and replaces one character; `path_to_root` returns `(pos, ref, alt)` triples with single-base semantics; `clones.tsv` `variant_path` is a comma list of SNV names. Indels shift downstream coordinates, so edits must be applied right-to-left, and names must be produced in the form `Himito call` writes them (`call::generate_variant_name`: SNV `m.<pos+1>`, indel `m.<anchor>` with anchor base in REF and ALT, left-aligned by default since T4).

4. **Truth names are ambiguous once HP tracts exist.** In a T16189C background, a +1 C in HVS1 can be spelled `m.16183A>AC` or `m.16189T>TC`, and `--collapse-hp` renames both to `m.16183ACCCCCTCCCC>ACCCCCCCCCCC`. `score_lineage.py` compares strings, so the same reconstruction scores differently with and without `--collapse-hp` (the T9 `hap_f1` 0.8421 → 0.8000 drop was exactly a renaming, on a single CA514 row). Truth and reconstruction need a common canonical form before comparison.

5. **`score_lineage.py` has no metric for the failure mode we care about.** On the real pairs the damage is HP alleles being ordered *above* true SNV markers (HG00232: `3643G>A` and `12520A>G` hang under `302A>ACC` in PacBio and under `16182A>AC` in ONT). `ad_f1` restricted to shared truth variants can't see it when the HP alleles are not truth variants, and `var_precision` counts every non-truth call as a false positive, so a *correct* HP-length call would be penalised. Nothing reports the ROOT-attached fraction, which is the single clearest PacBio/ONT discrepancy (3–18% vs 34–64%).

6. **The read error model is unmeasured.** pbsim3's `QSHMM-ONT-HQ` is trained on real data and does emit HP-length errors, but nobody has checked whether P(observed length | true length) in the simulated reads resembles the real ONT R10 curve (roughly 40–60% under-call of a C8 in the HG00133 pair). If pbsim is too clean, HP parameters tuned on the simulation will be wrong on real data; if it is too dirty, the benchmark will punish the caller for something the real data doesn't do.

7. **Smaller harness issues** (fix in passing, own commits): `clones.tsv` `n_reads` is always 0 (dead column); `run_simulation.sh` forwards `--internal-keep` but not `--topology`, `--denoise-vaf` or `--collapse-hp`, so it drifts from `sweep_seeds.sh`; `sweep_fpfn.sh` defaults `PROFILE=ont-r10` while every other script defaults `ont-denoised`; `simulate_reads.sh` rounds clone depth to `max(1, round(depth/2))` per rotation, so a clone at frequency 0.0005 still gets 2 reads at depth 1000 — fine, but `--hp-unstable` (T4) multiplies the clone count, so the floor must not silently inflate total depth (T4 Step 5 checks this).

## Design summary

| Mode | Flag (simulate_tree.py) | What it produces | Truth representation | Scored as |
|---|---|---|---|---|
| Background homoplasmies | `--background T16189C,A16183C` | Fixed SNVs in every clone incl. ref | `truth/background.txt` | Excluded from `var_*` denominators; changes tract sequences |
| Lineage-borne HP indel | `--hp-lineage-indels K` | K of the N edges are HP length changes (ISA) | Ordinary tree node, canonical name | `var_*`, `ad_*`, `hap_*` after canonicalisation; also `hp_var_*` |
| Lineage-independent length heteroplasmy | `--hp-unstable PRESET\|SPEC` | Each clone split into sub-clones by tract length state | `truth/hp_tracts.tsv`; **not** in the tree | `hp_intrusion`, `hp_len_l1`; excluded from `var_*` |

All three compose. `bench_standard.sh` is untouched; a new `bench_hp.sh` freezes the HP configuration.

## rCRS tract catalogue

Checked against `rCRS.fasta` on 2026-09-27 (1-based, inclusive). `anchor` is the base before the run, i.e. the position `Himito call` uses for a left-aligned insertion when no background SNV extends the run.

| id | run | base | len | anchor | rCRS context | common alleles (literature) | note |
|---|---|---|---|---|---|---|---|
| D310 | 303–309 | C | 7 | 302 | `A CCCCCCC T CCCCC` | +1, +2 (C8, C9); with T310C the run is C13 | most common somatic/germline length heteroplasmy |
| HVS1 | 16184–16188 | C | 5 | 16183 | `AA CCCCC T CCCC A` | with T16189C: C10, alleles C9–C12; with A16183C the anchor is 16182; with A16182C too, 16181 | the "16184 poly-C"; `m.16182A>AC` in HG00232 is this |
| C573 | 568–573 | C | 6 | 567 | `A CCCCCC A` | +1 … +7 (573.XC) | HG00253 carries +6 and +7 |
| T455 | 452–455 | T | 4 | 451 | `A TTTT C` | +1 (455.1T) | HG00253 `m.451A>AT` |
| C960 | 956–960 | C | 5 | 955 | `A CCCCC T CCCC` | +1 … +4 (960.XC) | haplogroup-associated |
| C5899 | 5895–5899 | C | 5 | 5894 | `A CCCCC A` | +1 (5899.1C) | |
| A12418 | 12418–12425 | A | 8 | 12417 | `A AAAAAAAA` | +1, −1 (12418.1A) | longest A-run in rCRS |
| CA514 | 514–523 | CA | 5 units | 513 | `G CACACACACA CC` | ±1 unit | period 2; already in `src/tracts.rs` |
| R8281 | 8272–8289 | CCCCCTCTA | 2 units | 8271 | `CACCCCCTCTA CCCCCTCTA GAG` | −9 bp (8281–8289del) | period 9; optional, off by default |

`src/tracts.rs::KNOWN_TRACTS` covers D310 (302–315), HVS1-C (16183–16193) and CA514 (513–525) only. T1 adds a test that every catalogue entry Himito is expected to handle is present there or is explicitly listed as "not yet in KNOWN_TRACTS", so the two lists cannot drift silently.

## Global constraints

- Tests first (`pytest lineage_eval/lineage_sim/tests`), watch them fail, implement, get green. Do not weaken existing tests; if an assertion changes because of the intended change, update it and say why in the commit message.
- **Existing seeds must reproduce byte-identically** when no new flag is passed. `build_tree` only draws from `rng` for SNV placement and parent choice; every new draw must happen *after* those, from a separately seeded `random.Random(seed + 1_000_003)`, so `truth/*` for `--n-mutations 10 --seed 1` with no HP flags is `cmp`-identical to today. T2 Step 5 checks this.
- Variant names in truth files must be exactly what `Himito call` (left-aligned, default since 2026-09-26 T4) writes for the same event on the same clone background. Verified against a real VCF in T2 Step 4.
- No change to `bench_standard.sh`, `baseline_standard.tsv`, or the `FIELDS` order of existing columns in `score_lineage.py`. New columns are appended.
- Cite the task in comments, e.g. `# Plan 2026-09-27 T3`.
- If something is ambiguous, stop and write `> QUESTION:` into the task instead of guessing.

---

### Task 1: Tract catalogue module (`hp_tracts.py`)

**Files:** new `lineage_eval/lineage_sim/hp_tracts.py`, `tests/test_hp_tracts.py`.

- [ ] **Step 1 (failing tests).**
  (a) `find_homopolymer_runs(seq, min_len=4)` returns `(start, end, base, length)` 1-based inclusive; on rCRS it returns 35 runs of length ≥ 6 and includes `(303, 309, 'C', 7)`, `(568, 573, 'C', 6)`, `(12418, 12425, 'A', 8)`, `(9478, 9484, 'T', 7)`.
  (b) `CATALOGUE` (the table above, as a list of `Tract` dataclasses with `id, start, end, unit, copies, anchor`) matches `seq[start-1:end]` for every entry.
  (c) `apply_background(seq, ["T16189C", "A16183C"])` returns a sequence where `runs_after` contains `(16183, 16193, 'C', 11)` and the HVS1 tract's effective anchor is 16182. With only `T16189C` the run is `(16184, 16193, 'C', 10)`, anchor 16183.
  (d) `effective_tract(tract, background_seq)` returns the run boundaries and anchor on the background genome, extending left/right over any base equal to the run base.
  (e) Parity with Rust: parse `src/tracts.rs` with a regex for `Tract { name: "...", start: N, end: M, ref_seq: "..." }` and assert each of D310, HVS1-C, CA514 appears in `CATALOGUE` with the same `ref_seq` over `[start, end]`. Entries of `CATALOGUE` marked `in_himito=False` are asserted **absent** from `tracts.rs`; when someone later adds them to Rust this test fails and the flag gets flipped. That is the intended coupling.
- [ ] **Step 2 (implement).** Pure functions, no I/O except `load_reference` reuse from `simulate_tree.py` (move `load_reference` into `hp_tracts.py` and re-export it from `simulate_tree.py` so nothing else changes).
- [ ] **Step 3.** Commit: `sim: rCRS homopolymer tract catalogue with background-aware anchors (plan T1)`.

---

### Task 2: Indel-capable clone genomes and Himito-compatible names

**Files:** `simulate_tree.py`, `tests/test_simulate_tree.py`.

- [ ] **Step 1 (failing tests).**
  (a) `Edit(pos, ref, alt)` where `ref`/`alt` follow VCF anchored semantics. `apply_edits(seq, edits)` applies SNVs and indels; when two edits touch the same anchor (the T16189C + `m.16189T>TC` case) the result equals applying them sequentially in ascending position order on a fresh copy. Test with `[Edit(16189,'T','C'), Edit(16183,'A','AC')]` on rCRS: the result has `C11` over 16184–16194 shifted by one, and length 16570.
  (b) `himito_name(edit)` gives `m.3144A>G` for an SNV and `m.302A>AC` / `m.302ACC>A` for an insertion / deletion at anchor 302. This mirrors `call::generate_variant_name` (SNV uses `pos`, indel uses anchor). Add a fixture comparing five names with lines from `test_data/*.vcf` (or `test.vcf` at repo root, which has `m.302A>AC`-type records).
  (c) `left_align(edit, seq)`: an insertion of `C` after 16188 on rCRS becomes anchor 16183 (`m.16183A>AC`); on the `T16189C` background it stays anchor 16183; on `T16189C,A16183C` it becomes anchor 16182 (`m.16182A>AC`). Deletion of one C at 309 becomes `m.302AC>A`. This is the Python twin of `call::left_align_variant` and must produce the same names; T2 Step 4 verifies on a real VCF.
  (d) `clone_genome` still passes `test_clone_genome_applies_only_path_snvs` unchanged.
- [ ] **Step 2 (implement).** Replace the `(pos, ref, alt)` tuples with `Edit`. `path_to_root` returns `list[Edit]`. `clone_genome` applies right-to-left by anchor so upstream coordinates never shift. `Node` gains `edit: Edit | None`; keep `pos/ref/alt` as properties for the existing tests. `variant_id` becomes a thin wrapper over `himito_name`.
- [ ] **Step 3: `write_truth`.** `clones.tsv` gains no new columns yet; `variant_path` uses canonical names (T6 defines canonicalisation; until then, left-aligned Himito names). Drop the dead `n_reads` column: the header becomes `clone_id\tnode_variant\tvariant_path\tfrequency`. Grep the harness for `$5` / `n_reads` readers (`simulate_reads.sh` reads `$1` and `$4` only) and update `test_write_truth_emits_expected_files_and_headers`.
- [ ] **Step 4: name parity with `Himito call`.** Build one clone with `Edit(16183,'A','AC')` on a `T16189C` background and one with `Edit(302,'A','AC')`, simulate 200× ONT reads of that single clone with `simulate_reads.sh`, run `Himito quick-start`/`call` with defaults (left-align on), and assert the VCF contains exactly `m.16183A>AC` and `m.302A>AC` (plus the background SNV at HF ≈ 1). Record the command and VCF lines in this plan under "T2 results". If the names differ, the Python left-aligner is wrong, not Himito; fix Python.
- [ ] **Step 5: byte-identity.** `python simulate_tree.py --n-mutations 10 --seed 1 --outdir /tmp/t2_after` and `git stash`-based run of the pre-T2 script into `/tmp/t2_before`; `cmp` `truth_mutation_tree.tsv`, `truth_variants.txt`, `clone_genomes.fa`. `clones.tsv` may differ only by the dropped `n_reads` column (show with `cut -f1-4 | cmp`).
- [ ] **Step 6.** Commit: `sim: indel-capable clone genomes with left-aligned Himito names (plan T2)`.

---

### Task 3: Lineage-borne HP indels (`--hp-lineage-indels K`)

**Files:** `simulate_tree.py`, `tests/test_simulate_tree.py`, `run_eval.sh`, `sweep_seeds.sh`.

- [ ] **Step 1 (failing tests).**
  (a) `build_tree(..., hp_lineage_indels=3, hp_rng=random.Random(...))` yields a tree with `n_mutations` nodes of which exactly 3 have `edit.is_indel()`, each inside a `CATALOGUE` tract, drawn without replacement across tracts (two indels never share a tract, so the ISA truth stays a strict tree of mutually compatible alleles). Length change is drawn from `{+1: 0.6, +2: 0.2, -1: 0.2}` for HP tracts and `{±1 unit}` for CA514.
  (b) The other `n_mutations − K` SNVs still respect `AVOID_RANGES` and never fall inside a tract used by (a).
  (c) With `hp_lineage_indels=0` (default) the tree is byte-identical to today for every seed in `1..5` (compare `truth_mutation_tree.tsv` text).
  (d) A clone whose path carries both an HP indel and a background SNV in the same tract (T5) gets the name from `left_align` on the *background* sequence, e.g. `m.16182A>AC`, not `m.16183A>AC`.
- [ ] **Step 2 (implement).** After the SNV placement loop (so existing draws are untouched), pick `K` distinct node ids from `1..n_mutations` using `hp_rng`, and replace their SNV edit with an HP edit at a tract drawn from `CATALOGUE` entries with `in_himito` or `simulate=True`. Store `tract_id` on the node. `write_truth` adds `truth/hp_lineage.tsv`: `node_id\tvariant\ttract\tdelta`.
- [ ] **Step 3: harness.** `run_eval.sh` and `sweep_seeds.sh` accept `--hp-lineage-indels K` and forward it to `simulate_tree.py`. Empty = omit the flag. In the same commit, `run_himito.sh` and `run_eval.sh` accept `--hp-fp-rate F` / `--hp-fn-rate F` and forward them to `Himito lineage` only (the flags already exist in `main.rs`); T10 and T11 depend on this.
- [ ] **Step 4.** `pytest`, commit: `sim: lineage-borne homopolymer length changes as tree edges (plan T3)`.

---

### Task 4: Lineage-independent length heteroplasmy (`--hp-unstable`)

**Files:** `simulate_tree.py`, `simulate_reads.sh`, `tests/test_simulate_tree.py`.

**Model.** Each unstable tract `t` has a length distribution `π_t` over deltas relative to the (background-adjusted) reference run, e.g. D310 `{0: 0.55, +1: 0.35, +2: 0.10}`. Every clone `c` (including `ref`) is split into sub-clones, one per combination of tract states, with frequency `f_c · Π_t π_t(δ_t)`; sub-clones below `--hp-min-subclone-freq` (default 0.002) are dropped and their mass renormalised within the clone so total mass is unchanged. Tract length is therefore independent of lineage by construction — the same length appears in every clone at the same proportion. That is the homoplasy the real trees show, in its strongest form; `--hp-unstable-drift D` (optional, Step 6) lets `π_t` differ per clone by a Dirichlet perturbation so intermediate cases can be simulated too.

- [ ] **Step 1 (failing tests).**
  (a) `parse_hp_unstable("common")` returns the preset `{D310: {0:0.55,1:0.35,2:0.10}, HVS1: {0:0.6,1:0.3,2:0.1}, C573: {0:0.7,1:0.2,2:0.1}}` and `parse_hp_unstable("D310:0=0.5,1=0.5;C573:0=0.8,-1=0.2")` parses the explicit form; distributions must sum to 1 ± 1e-9 and reference every id in `CATALOGUE`, else `ValueError`.
  (b) `expand_unstable(clones, spec, min_freq)` on two clones with masses 0.6/0.4 and one tract `{0:0.5,1:0.5}` returns four sub-clones `1.h0`, `1.h1`, `2.h0`, `2.h1` with frequencies 0.3/0.3/0.2/0.2, each with `hp_state = {D310: δ}`; total frequency 1.0 ± 1e-9.
  (c) Sub-clone genomes: `1.h1` equals clone 1's genome with one extra C at the D310 anchor (apply after the clone's path edits, on the background). When a clone's lineage path already carries a lineage-borne indel in the same tract (T3), the unstable δ is applied on top and the tract is reported as `mixed` in `hp_tracts.tsv`.
  (d) With no `--hp-unstable`, `clones.tsv` and `clone_genomes.fa` are unchanged (byte-identical to T2 output).
- [ ] **Step 2 (implement).** `clones.tsv` gains a column `hp_state` (`D310:+1;C573:0`, empty for non-HP runs) appended after `frequency` so `simulate_reads.sh`'s `$1`/`$4` reads still work. `clone_genomes.fa` headers become `>clone_<id>` where `<id>` is the sub-clone id (`3.h2`); `simulate_reads.sh`'s awk splits on `>clone_` so dots in ids are fine — add a test that `awk -v dir=... '/^>/{id=substr($0,8)...}'` produces `clone_3.h2.fa` (run the awk in the test via `subprocess`).
- [ ] **Step 3: truth files.** New `truth/hp_tracts.tsv`: `tract\tanchor\tref_len\tbackground_len\tdelta\ttruth_freq\tcanonical_name`, one row per (tract, δ) with δ ≠ 0, `truth_freq` = total mass carrying that δ. `canonical_name` is defined in T6 (write the left-aligned Himito name here; T6 rewrites it). `truth_variants.txt` does **not** list these; `truth_mutation_tree.tsv` does not contain them.
- [ ] **Step 4: read simulation.** `simulate_reads.sh` needs no logic change beyond ids. Verify a `--hp-unstable common` run at depth 1000 yields total reads within 5% of the no-HP run for the same seed (the `max(1, …)` floor must not inflate depth; if it does, drop sub-clones whose per-rotation depth rounds to 0 *before* renormalising instead of flooring at 1, and log how much mass was dropped).
- [ ] **Step 5: harness.** `run_eval.sh` / `sweep_seeds.sh` accept `--hp-unstable SPEC` and `--hp-min-subclone-freq F`.
- [ ] **Step 6 (optional, own commit).** `--hp-unstable-drift D`: per clone, `π_t,c ~ Dirichlet(π_t / D)`; `D = 0` disables. Recorded per clone in `hp_state` with its own frequencies. Skip if the owner has not asked for it by the time T4 is otherwise done.
- [ ] **Step 7.** `pytest`, commit: `sim: lineage-independent homopolymer length heteroplasmy as per-clone tract states (plan T4)`.

---

### Task 5: Background homoplasmies (`--background`)

**Files:** `simulate_tree.py`, `tests/test_simulate_tree.py`, `run_eval.sh`, `sweep_seeds.sh`, `score_lineage.py` (exclusion only).

- [ ] **Step 1 (failing tests).**
  (a) `--background T16189C,A16183C` applies both SNVs to the root genome before any tree edit; `clone_genomes.fa` `clone_ref` differs from rCRS at exactly those positions.
  (b) SNV placement never picks a background position; `AVOID_RANGES` still applies.
  (c) `truth/background.txt` lists `m.16189T>C` and `m.16183A>C` (Himito names on rCRS).
  (d) Tract anchors used by T3/T4 come from `effective_tract` on the background genome (the HVS1 +1 C is named `m.16182A>AC`).
  (e) A preset `--background hvs1-polyC` expands to `T16189C,A16183C`; `--background none` is the default.
- [ ] **Step 2: scorer exclusion.** `score_lineage.py --truth-background truth/background.txt`: names in that file are removed from `detected` before `var_precision` is computed and counted in a new column `n_background_calls`. They are fixed homoplasmies, and `Himito lineage` already drops them as `fixed`; the VCF will still report them at HF ≈ 1.
- [ ] **Step 3.** Commit: `sim: fixed background homoplasmies that reshape rCRS homopolymer tracts (plan T5)`.

---

### Task 6: Canonical names and HP metrics in `score_lineage.py`

**Files:** `score_lineage.py`, `hp_tracts.py` (shared canonicaliser), `tests/test_score_lineage.py`, `run_eval.sh`, `sweep_fpfn.sh`.

**Canonical form.** For any name whose REF span lies inside an effective tract on the background genome, apply the edit to the tract's background sequence and rename it `m.<anchor><tract_ref>><tract_alt>` — exactly what `lineage::collapse_tracts` produces with `--collapse-hp`, but computed from the Python catalogue and the *background* genome rather than rCRS. Names outside tracts are unchanged. Both truth and reconstruction go through this, so a run with `--collapse-hp` and one without score the same, and `m.16183A>AC` / `m.16189T>TC` / `m.16182A>AC` spellings of one allele collide as they should.

- [ ] **Step 1 (failing tests).**
  (a) `canonical_name("m.16189T>TC", background=["T16189C"])` == `canonical_name("m.16183A>AC", background=["T16189C"])` == `m.16183ACCCCCTCCCC>ACCCCCCCCCCC`. The canonical REF is the **rCRS** tract sequence (not the background one) and the ALT is the full tract sequence on the read with the background SNV applied, matching `tracts.rs::tract_allele_name` and what `--collapse-hp` emits. A run without `--collapse-hp` therefore canonicalises to the same string as one with it. On the `T16189C,A16183C` background the tract is extended to anchor 16182 and the canonical REF becomes `m.16182AACCCCCTCCCC>…`; add that as a second assertion so the extended-anchor case (the HG00232 `m.16182A>AC` miss) is covered.
  (b) `canonical_name("m.3144A>G", ...)` is unchanged.
  (c) A tract-allele name already in collapsed form is returned unchanged.
  (d) `score()` with a truth of `{m.16183A>AC}` and a reconstruction of `{m.16189T>TC}` on the `T16189C` background gives `var_recall = 1.0` after canonicalisation, 0.0 without.
- [ ] **Step 2: apply canonicalisation** to `truth_vars`, `truth_parent` keys/values, `recon_parent` keys/values (including comma-joined `ROOT_BLOCK` members), VCF names, and matrix row names in `parse_matrix_haplotypes`. `--background` is passed to the scorer so the effective tracts are right.
- [ ] **Step 3: HP-aware detection metrics.** Read `truth/hp_tracts.tsv` (T4) via `--truth-hp-tracts`. Its `canonical_name`s form `unstable`. Then: `var_precision` denominator excludes `unstable ∪ background` (they are real alleles, not false positives); new columns `n_hp_unstable_truth`, `n_hp_unstable_called` (how many of the truth δ alleles appear in the VCF), `hp_var_precision` / `hp_var_recall` restricted to `truth/hp_lineage.tsv` names (T3).
- [ ] **Step 4: `hp_intrusion`.** From the recon tree after canonicalisation: the number of truth SNV nodes whose recon ancestor set contains any `unstable` name, divided by the number of truth SNVs in the recon tree. 0.0 is perfect. Also `hp_backbone_nodes`: count of `unstable` names that have ≥ 1 child in the recon tree. These are the two numbers that would have flagged the HG00232 PacBio/ONT disagreement.
- [ ] **Step 5: `hp_len_l1`.** For each unstable tract, the reconstruction's length distribution is the fraction of reads in `cleaned_matrix_all_variants.csv` carrying each canonical δ allele (reads carrying none are δ = 0). `hp_len_l1` is the mean over tracts of `Σ_δ |truth_freq(δ) − recon_freq(δ)|`. Report `NA` when no unstable tract exists.
- [ ] **Step 6: `root_fraction`.** Fraction of reads attached to ROOT in `mutation_tree.tsv` (`n_reads_attached` of the ROOT row over the column total). Always reported; this is the PacBio/ONT discrepancy metric and costs nothing.
- [ ] **Step 7: wiring.** `run_eval.sh` passes `--truth-background`, `--truth-hp-tracts` and `--background` when the files exist; `sweep_fpfn.sh` likewise. `FIELDS` gains the new columns at the end. Update `test_all_ones_row_outside_keep_does_not_change_score_haplotypes` only if the signature changed.
- [ ] **Step 8.** `pytest`, commit in two parts: `eval: canonicalise homopolymer tract alleles before scoring (plan T6a)` and `eval: hp_intrusion, hp_len_l1, root_fraction and HP-aware var_precision (plan T6b)`.

---

### Task 7: Measure the HP error model (report first, augment only if needed)

**Files:** new `lineage_eval/lineage_sim/measure_hp_errors.py`, optional `inject_hp_errors.py`.

- [ ] **Step 1: `measure_hp_errors.py --bam X --ref rCRS.fasta [--min-len 4] [--tracts hp_tracts]`.** For every reference HP run of length ≥ `min_len` (from `find_homopolymer_runs`), walk each primary alignment that spans `[anchor, end+2]` and compute the read's run length from the CIGAR (same rule as `denoise_indel::hp_allele`: inserted bases inside the run must equal the run base, else the read is uninformative). Output `hp_confusion.tsv`: `base\tref_len\tobs_len\tstrand\tn`. Skip runs listed in `--tracts` and any run where the modal observed length ≠ `ref_len` (the sample is not reference there). Unit test on a 3-read synthetic BAM built with `pysam` (add `pysam` to `environment.yml`).
- [ ] **Step 2: run it** on (i) a `--profile ont-denoised` cell's `aln.sorted.bam` (raw pbsim ONT-HQ), (ii) the same cell's `aln.denoised.bam`, (iii) a `--profile hifi` cell, and (iv) the real pair `lineage_eval/benchmark_results/HPRC_debug/HG00232.ont.bam` plus its PacBio counterpart if present (else `test_data/HG002_chrM_100x.bam`). Write P(obs = ref_len − 1 | ref_len) and P(obs = ref_len + 1 | ref_len) for C/G and A/T runs of length 5–9 into this plan under "T7 results", one table per BAM.
- [ ] **Step 3: decision.** If pbsim ONT-HQ under-calls a C8 at less than half the real ONT rate, implement `inject_hp_errors.py`: reads the pbsim `*.maf.gz` truth alignments (pbsim3 writes them next to the `.fq.gz`; `simulate_reads.sh` must stop deleting `$WORK` before they are used), and for each HP run in the read's true alignment resamples the observed length from the real confusion table (`--hp-error-table hp_confusion.tsv`), rewriting sequence and quality. Opt-in via `simulate_reads.sh --hp-error-table`. Otherwise record "pbsim is adequate" with the numbers and skip.
- [ ] **Step 4.** Commit: `eval: measure homopolymer length confusion from a BAM (plan T7)` (+ `sim: inject calibrated HP length errors from pbsim MAF (plan T7b)` if built).

---

### Task 8: Frozen HP benchmark and the in-silico PacBio/ONT pair

**Files:** new `bench_hp.sh`, `run_eval.sh` (one flag), `README.md`.

- [ ] **Step 1: `bench_hp.sh`.** Same shape as `bench_standard.sh` with its own `baseline_hp.tsv`. Frozen config: `--profile ont-denoised --n-mutations 10 --total-depth 1000 --topology chain --internal-keep 0.10 --frequency-threshold 0.05 --denoise-vaf 0.03 --background hvs1-polyC --hp-lineage-indels 2 --hp-unstable common`, seeds 1–30. Focus metrics add `hp_intrusion`, `hp_backbone_nodes`, `hp_len_l1`, `root_fraction`, `hp_var_f1`.
- [ ] **Step 2: pair mode.** `run_eval.sh --pair` simulates reads twice from the same `truth/` (profiles `hifi` and `ont-denoised`, same seed), runs `run_himito.sh` for each into `himito_pb/` and `himito_ont/`, scores each, and adds `pair_ad_f1` = `ad_f1` of the ONT tree scored *against the PacBio tree as if it were truth*, plus `pair_root_fraction_delta`. This is the in-silico twin of the eight real pairs; on the current code it should reproduce the ONT ROOT inflation and the HP backbone permutation. Note that `hifi` needs `ccs` on Linux x86; document that `--pair` is Linux-only for now.
- [ ] **Step 3: README.** New section "Homopolymer length heteroplasmy" describing the three modes, the catalogue, canonical names, the new metrics, and that `bench_standard.sh` is unchanged and still SNV-only on purpose. Update the "SNV-only truth" bullet to point here.
- [ ] **Step 4.** Commit: `eval: bench_hp.sh frozen HP benchmark and --pair PacBio/ONT mode (plan T8)`.

---

### Task 9: First measurements (report only, no code change)

- [ ] Run `bench_hp.sh` at the current `Himito` head, three ways: defaults; `--collapse-hp`; `--collapse-hp --hp-fn-rate 0.4` (ONT). Paste `hp_intrusion`, `hp_backbone_nodes`, `hp_len_l1`, `root_fraction`, `ad_f1`, `hap_f1` means ± 2 SE into "T9 results".
- [ ] Run `run_eval.sh --pair` on seeds 1–5 and report `pair_ad_f1` and both `root_fraction`s. If the ONT `root_fraction` is not markedly higher than PacBio's, T7's error model is too clean and T7 Step 3 must be revisited before this benchmark is trusted.
- [ ] Do **not** tune defaults in this plan; the numbers feed the 2026-09-27 lineage proposals (leaf-only HP alleles, categorical tract observations, BAM-calibrated rates).

### Task 10: Revise `wdl/Himito-lineage/Himito_lineage_sim.wdl` for the HP modes

**Files:** `wdl/Himito-lineage/Himito_lineage_sim.wdl`, `lineage_eval/lineage_sim/Dockerfile` (rebuild note only), `lineage_eval/lineage_sim/Himito_lineage_sim.inputs.json`.

The WDL is a port of `run_simulation.sh`: `ValidateConfig` → scatter `RunEvalCell` (one `run_eval.sh` per seed × n_mutations × depth) → `CombineMetrics`. It must learn the new `run_eval.sh` flags, ship the new truth files, and keep the existing inputs JSON valid. Defaults keep every new mode **off**, so a run with today's inputs produces the same cells as today (modulo the appended metric columns).

- [ ] **Step 1: workflow inputs.** Add, with the same comment style as the existing gates block:
  ```wdl
  # Plan 2026-09-27. Homopolymer length heteroplasmy; all off by default so an
  # unchanged inputs JSON reproduces the SNV-only benchmark.
  String background = "none"          # none | hvs1-polyC | comma list e.g. "T16189C,A16183C"
  Int hp_lineage_indels = 0           # K of the n_mutations edges are HP length changes
  String hp_unstable = ""             # "" | common | "D310:0=0.5,1=0.5;C573:..."
  Float hp_min_subclone_freq = 0.002
  Boolean collapse_hp = false         # Himito lineage --collapse-hp
  Float? hp_fp_rate                   # Himito lineage --hp-fp-rate (unset = 2x fp)
  Float? hp_fn_rate                   # Himito lineage --hp-fn-rate (unset = 2x fn)
  ```
  Build `background_arg`, `hp_lineage_arg`, `hp_unstable_arg`, `collapse_hp_arg`, `hp_fp_arg`, `hp_fn_arg` strings the way `fp_arg`/`fn_arg` are built today (empty string when off), and pass them into `RunEvalCell`. `run_eval.sh` must accept `--hp-fp-rate` / `--hp-fn-rate` and forward them to `run_himito.sh` → `Himito lineage` (add to T3 Step 3's flag list if not already there; `run_himito.sh` currently forwards only `--fp`/`--fn`/`--collapse-hp`).
- [ ] **Step 2: `ValidateConfig`.** Extend the Python block: `hp_unstable` must parse with `hp_tracts.parse_hp_unstable` (import from `/opt/lineage_sim`), `background` must be `none`, a preset, or a list of `<REF><pos><ALT>` tokens whose REF matches rCRS, `hp_lineage_indels` must satisfy `0 <= K <= n_mutations` for every `n_mutations` in the grid and `K <= number of simulate=True tracts` (T3 draws without replacement). Warn (not error) when `hp_unstable` is set and `collapse_hp` is false, since the unstable alleles then enter the tree under several spellings. Warn when `hp_fp_rate`/`hp_fn_rate` are set and no HP mode is on (they have nothing to act on).
- [ ] **Step 3: `RunEvalCell`.** Forward the six new args. Add the mode to the cell name so bundles from different configurations do not collide: `cell = "seed~{seed}_mut~{n_mutations}_depth~{total_depth}" + (if hp_lineage_indels > 0 || hp_unstable != "" then "_hp" else "")`. Prepend the identifying columns `background`, `hp_lineage_indels`, `hp_unstable`, `collapse_hp`, `hp_fp_rate`, `hp_fn_rate` after `internal_keep` in `cell_metrics.tsv` (write `NA` for unset optionals). Add new outputs alongside the bundle so a sweep can reuse a cell without untarring 50 MB: `File truth_dir_tar = "~{cell}_truth.tar.gz"` (all of `truth/`, which now includes `hp_tracts.tsv`, `hp_lineage.tsv`, `background.txt`), `File matrix = "~{cell}/himito/sim.matrix.csv"`, `File vcf = "~{cell}/himito/sim.vcf"`, `File cleaned_matrix = "~{cell}/himito/sim_lineage.cleaned_matrix_all_variants.csv"`. `parameter_meta` entries for every new input.
- [ ] **Step 4: `CombineMetrics`.** No change in logic, but it concatenates by position, so add a guard: every file's header must equal the first file's header, else exit 1 naming the file. Mixed-schema cells (e.g. a cell that failed before scoring) used to silently misalign columns.
- [ ] **Step 5: image and inputs.** Update the Dockerfile header comment to say the image must be rebuilt after this plan (it bakes `/opt/lineage_sim/*.py`; `hp_tracts.py`, `measure_hp_errors.py` and `pysam` from T7 are new). Add `pysam` to the image's conda/pip install. Add `Himito_lineage_sim.hp.inputs.json` next to the existing inputs JSON with the `bench_hp.sh` configuration (`background = "hvs1-polyC"`, `hp_lineage_indels = 2`, `hp_unstable = "common"`, `topology = "chain"`, `internal_keep = 0.10`, `frequency_threshold = 0.05`, `denoise_vaf = 0.03`, seeds 1–30, one `n_mutations` = 10, one depth = 1000). Leave the existing inputs JSON untouched.
- [ ] **Step 6: check.** `womtool validate` (or `miniwdl check`) on the WDL; run `miniwdl run --dir /tmp/t10 ... Himito_lineage_sim.wdl` with `seeds=[1]`, `n_mutations=[6]`, `depths=[200]` twice, once with the old inputs JSON and once with the HP inputs, using a locally built image. Paste the two `metrics.tsv` headers into this plan under "T10 results" and confirm the SNV-only run's non-HP columns match a local `run_eval.sh` at the same seed.
- [ ] **Step 7.** Commit: `wdl: forward HP simulation modes and Himito HP rates through Himito_lineage_sim.wdl (plan T10)`.

---

### Task 11: `Himito_lineage_hp_sweep.wdl` — sweep `--hp-fp-rate` × `--hp-fn-rate` and report the best pair

**Files:** new `wdl/Himito-lineage/Himito_lineage_hp_sweep.wdl`, new `lineage_eval/lineage_sim/sweep_hp_rates.sh` (the local twin, same shape as `sweep_fpfn.sh`), new `lineage_eval/lineage_sim/pick_best_hp_rates.py`, `tests/test_pick_best_hp_rates.py`, `Himito_lineage_hp_sweep.inputs.json`.

**Design.** Simulation, alignment, denoise, build and call are independent of the HP rates, so run them **once per seed** and re-run only `Himito lineage` + `score_lineage.py` per grid cell, exactly as `sweep_fpfn.sh` does for the global rates. The workflow is:

```
ValidateConfig (T10's task, imported)
  └─ scatter seed → SimulateCell        (run_eval.sh minus the lineage/score steps;
                                          emits truth.tar.gz, sim.matrix.csv, sim.vcf)
       └─ scatter (hp_fp, hp_fn) → LineageScoreCell
                                         (Himito lineage --hp-fp-rate --hp-fn-rate [--collapse-hp];
                                          score_lineage.py with the T6 truth files)
  └─ CombineMetrics (T10's, with header guard)
  └─ PickBestHpRates → best_hp_rates.json, hp_sweep_summary.tsv, hp_sweep_report.md
```

`SimulateCell` is `RunEvalCell` with a `--stop-after call` flag added to `run_eval.sh` (Step 1), so the two WDLs share the Docker image and the scripts.

- [ ] **Step 1: `run_eval.sh --stop-after call`.** After `run_himito.sh` has produced `sim.matrix.csv` and `sim.vcf`, exit 0 without running `lineage` or the scorer. Implement by forwarding `--stop-after call` to `run_himito.sh`, which skips the lineage step. Test by hand: the outdir has `himito/sim.matrix.csv`, `himito/sim.vcf`, `truth/`, and no `sim_lineage.*`.
- [ ] **Step 2: `sweep_hp_rates.sh` (local twin).** Copy `sweep_fpfn.sh`; grids `HP_FP_GRID` and `HP_FN_GRID` default to multiples of the profile's resolved global rates, `{1,2,4,8,16}×fp` and `{1,2,4,8,16,32}×fn` (for `ont-denoised` fn = 0.01 that reaches 0.32, which is where the real ONT dropout at a C8 sits; cap at 0.5). Pass `--hp-fp-rate`/`--hp-fn-rate` to `Himito lineage`; the global `--fp-rate`/`--fn-rate` stay at the profile default unless `--fp`/`--fn` are given. Score with `--truth-background`, `--truth-hp-tracts`, `--background` (T6 Step 7), and append the `hp_fp_rate`/`hp_fn_rate` columns. Accept `--collapse-hp` and `--seeds` (loop over `seed*/` cells produced by `sweep_seeds.sh --keep-intermediate`, or by Step 1) so the local script can also sweep across seeds.
- [ ] **Step 3: `pick_best_hp_rates.py` (tests first).**
  Input: the combined `metrics.tsv` (columns include `seed`, `hp_fp_rate`, `hp_fn_rate` and the metric columns). Output: `hp_sweep_summary.tsv` (one row per grid cell: mean and SE over seeds for each focus metric, `n_seeds`), `best_hp_rates.json`, and `hp_sweep_report.md`.
  Selection rule, in this order, and written into the report so the choice is auditable:
  1. Primary objective `--objective` (default `hp_intrusion`, minimised; `ad_f1`, `hap_f1`, `hp_len_l1` selectable). Cells whose mean is within `2·SE_pooled` of the best mean are **tied** with it (the same noise band `bench_standard.sh` uses).
  2. Among tied cells, secondary objective `--secondary ad_f1` (maximised), same tie rule.
  3. Among cells still tied, prefer the one closest to the global rates (smallest `log2(hp_fn/fn) + log2(hp_fp/fp)`): don't move a parameter the data can't distinguish.
  4. Refuse to declare a best cell (`"best": null`, exit 0, report says why) when the primary objective's best and worst cells are themselves within `2·SE` — the grid is flat and no parameter should change. This is the guard the T7 sweep lacked.
  Also report `root_fraction` at the chosen cell versus at the global-rate cell, since that is the PacBio/ONT symptom the rates are meant to fix.
  Tests: (a) a synthetic table where one cell is clearly best on `hp_intrusion` → chosen; (b) two cells tied on primary, one better on `ad_f1` → the latter; (c) all cells within noise → `best: null` and the report contains "flat"; (d) the summary's SE for a 3-seed cell equals `sd/√3`; (e) a table missing the `hp_fn_rate` column raises a clear error.
- [ ] **Step 4: the WDL.**
  - Inputs: everything `Himito_lineage_sim.wdl` takes (import its `ValidateConfig` and `CombineMetrics` via `import "Himito_lineage_sim.wdl" as sim`), plus `Array[Float] hp_fp_grid`, `Array[Float] hp_fn_grid`, `Boolean collapse_hp = false`, `String objective = "hp_intrusion"`, `String secondary = "ad_f1"`. One `n_mutations` and one `depth` (scalars, not arrays): a sweep is over rates and seeds, not over the simulation grid. `ValidateConfig` additionally requires an HP mode to be on (`hp_lineage_indels > 0` or `hp_unstable != ""`); otherwise the rates have nothing to act on and the sweep is guaranteed flat — fail early with that sentence.
  - `SimulateCell`: as `RunEvalCell` with `--stop-after call`; outputs `truth_tar`, `matrix`, `vcf`, `background_txt` (may be empty). 1 CPU, 4 GiB, 50 GB SSD, preemptible 2, like today.
  - `LineageScoreCell`: inputs `matrix`, `vcf`, `truth_tar`, `hp_fp_rate`, `hp_fn_rate`, `collapse_hp`, profile, `min_hf`/`max_hf`, `fp_arg`/`fn_arg`, `mcmc_iterations = 2000`, `mcmc_chains = 2` (the `sweep_fpfn.sh` budget; expose both). Runs `Himito lineage` and `score_lineage.py`, and prepends `seed`, `hp_fp_rate`, `hp_fn_rate`, `collapse_hp` columns. Small: 1 CPU, 2 GiB, 10 GB HDD, preemptible 3. Nested scatter: `scatter (seed) { call SimulateCell } scatter (seed) { scatter (fp) { scatter (fn) { call LineageScoreCell } } }`, cross-referencing `SimulateCell` outputs by index (WDL 1.0 has no `zip` over three arrays; index into `SimulateCell.matrix[i]` with `range(length(seeds))`).
  - `PickBestHpRates`: runs `pick_best_hp_rates.py` on `CombineMetrics.combined`.
  - Outputs: `File metrics`, `File summary`, `File best_hp_rates` (JSON: `{"hp_fp_rate": …, "hp_fn_rate": …, "objective": …, "value": …, "se": …, "n_seeds": …, "flat": false, "root_fraction_best": …, "root_fraction_global": …}`), `File report` (Markdown table of the grid with the chosen cell marked, plus the selection trace), `Array[File] per_cell_metrics`.
- [ ] **Step 5: inputs JSON.** `Himito_lineage_hp_sweep.inputs.json`: the `bench_hp.sh` configuration, seeds 1–10 (30 × 30 cells is too expensive for a first pass; note this in the file), `hp_fp_grid = [0.0001, 0.0002, 0.0004, 0.0008, 0.0016]`, `hp_fn_grid = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32]`, `collapse_hp = true`.
- [ ] **Step 6: check.** `womtool validate` both WDLs. `miniwdl run` with `seeds=[1,2]`, a 2×2 grid, `n_mutations=6`, `depth=200`, locally built image. Confirm `SimulateCell` ran twice and `LineageScoreCell` eight times, the combined table has 8 rows, and `best_hp_rates.json` is either a cell or `flat`. Run `sweep_hp_rates.sh` locally on the same two seeds and `cmp` the `hp_sweep_summary.tsv` from both paths (they must agree, since the same scripts run underneath). Paste into "T11 results".
- [ ] **Step 7: README.** Add `sweep_hp_rates.sh` and the two WDLs to the README section from T8 Step 3, with the selection rule spelled out and the warning that a flat grid means "leave the default alone", not "pick the top row".
- [ ] **Step 8.** Commit in two parts: `eval: sweep_hp_rates.sh and pick_best_hp_rates.py with a noise-aware selection rule (plan T11a)` and `wdl: Himito_lineage_hp_sweep.wdl sweeps --hp-fp-rate x --hp-fn-rate per seed and reports the best pair (plan T11b)`.

---

## On hold / out of scope

- Somatic *de novo* HP mutation along branches (a per-edge rate rather than a global `π_t`): T4 Step 6's drift covers the middle ground; a true birth-death process is not needed to test the current failure mode.
- NUMT-derived recombination and back-mutation at SNVs: separate plan (`finite-sites-relaxation.md`).
- 8281–8289 9-bp deletion: in the catalogue with `simulate=False`; enable once period-9 repeats are handled by `collapse_tracts`.

## Final verification

- [ ] `pytest lineage_eval/lineage_sim/tests` green; report the pass count.
- [ ] `cargo test` green (T1 Step 1e reads `src/tracts.rs` but changes nothing in Rust).
- [ ] `python simulate_tree.py --n-mutations 10 --seed 1` with no HP flags reproduces the pre-plan `truth/` byte-for-byte except the dropped `n_reads` column (T2 Step 5 output attached).
- [ ] `bench_standard.sh --seeds "1 2 3"` still runs and its `metrics.tsv` columns are a prefix of the new `FIELDS`.
- [ ] T2 Step 4 name-parity VCF lines, T7 confusion tables and T9 tables are attached to this plan.
- [ ] `womtool validate` passes on `Himito_lineage_sim.wdl` and `Himito_lineage_hp_sweep.wdl`; the T10 and T11 `miniwdl` smoke runs are attached; the existing `Himito_lineage_sim.inputs.json` still runs unchanged.
- [ ] The Docker image has been rebuilt from the post-plan `lineage_sim/` and the tag recorded here.

## T2 results

Name parity with `Himito call`, left-align on (VCF header `##HimitoNormalization=left-aligned`). Binary: `target/release/Himito`. Reads: `simulate_reads.sh` from `lineage_eval/lineage_sim/` with `PATH="/opt/miniconda3/envs/himito-eval/bin:$PATH"`. Two single-clone truths (not the random tree), each at frequency 1, 200× `--profile ont-denoised --seed 1`.

```bash
./simulate_reads.sh --outdir /tmp/t2_parity_hvs1 --profile ont-denoised --total-depth 200 --seed 1
./simulate_reads.sh --outdir /tmp/t2_parity_d310 --profile ont-denoised --total-depth 200 --seed 1
```

`/tmp/t2_parity_hvs1` genome is rCRS plus `T16189C` plus `Edit(16183,'A','AC')` after `left_align` on that background (name stays `m.16183A>AC`). `/tmp/t2_parity_d310` is rCRS plus `Edit(302,'A','AC')` (name `m.302A>AC`).

Call path matches `run_himito.sh` through `Himito call` (lineage not run). No `--no-left-align`.

```bash
minimap2 -ax lr:hq -t 4 rCRS.fasta reads.fastq.gz | samtools sort -o himito/aln.sorted.bam
Himito denoise -i himito/aln.sorted.bam -o himito/aln.denoised.bam -r rCRS.fasta -d ont-denoised --vaf 0.01 --indels
Himito build -i himito/aln.denoised.bam -r rCRS.fasta -k 21 -o himito/sim.gfa -l 3000 --min-edge-reads 2
Himito call -g himito/sim.gfa -r rCRS.fasta -s SIM -d ont-denoised -o himito/sim.vcf -k 21 \
  --input-bam himito/aln.denoised.bam -m 2 -v 0.01 \
  --permutation-frequency-threshold 0.7 --strand-bias-threshold 0.05 --indel-false-threshold 0.1
```

HVS1 VCF records (only these two):

```
chrM	16183	.	A	AC	.	PASS	DP=96	GT:AD:HF	1:47:0.48958334
chrM	16189	.	T	C	.	PASS	DP=96	GT:AD:HF	1:88:0.9166667
```

Names `m.16183A>AC` and `m.16189T>C` (HF 0.917).

D310 VCF record (only this one):

```
chrM	302	.	A	AC	.	PASS	DP=47	GT:AD:HF	1:28:0.59574467
```

Name `m.302A>AC`. Python left-align was not changed after this run.

Step 5 (`--n-mutations 10 --seed 1`): `cmp` identical for `truth_mutation_tree.tsv`, `truth_variants.txt`, and `clone_genomes.fa` between `/tmp/t2_before` (T1 script `b7842ab`) and `/tmp/t2_after`. `cut -f1-4` of both `clones.tsv` files is identical. The pre-T2 header still ends in `n_reads`; the post-T2 header is `clone_id	node_variant	variant_path	frequency`.
