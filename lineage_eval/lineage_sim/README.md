# Himito-lineage simulation evaluation

Simulate a ground-truth mitochondrial clonal mutation tree, generate long reads
(PacBio HiFi or ONT-R10), reconstruct the tree with `Himito call` + `Himito
lineage`, and score reconstruction accuracy across replicate seeds.

**Start here:** [The three frequency gates](#the-three-frequency-gates-read-this-first)
— the defaults silently empty the VCF on simulated data — and
[Why `ad_recall` is scored against *all* truth pairs](#why-ad_recall-is-scored-against-all-truth-pairs),
which governs how the tree metrics may be compared across runs.

## Truth model

The truth is a **clonal mutation tree** (SCITE-style), matching what `Himito
lineage` reconstructs: root = rCRS, each edge adds one heteroplasmic SNV, each
node has a cumulative heteroplasmy frequency that sets its simulated read depth.
mtDNA circularity is handled by simulating reads from two rotations of each
clone genome (offset 0 and offset 8284) and pooling them, so the rCRS origin is
covered by full-length reads.

## Platform requirements

**ONT (`--profile ont-denoised`, or `ont-r10`)** runs on any platform (Linux
x86, macOS ARM64). pbsim3's `QSHMM-ONT-HQ` model needs no extra tools beyond the
conda env. Prefer `ont-denoised` — `ont-r10` skips the denoise step and yields an
empty VCF; see [Use `ont-denoised`, not `ont-r10`](#use-ont-denoised-not-ont-r10).

**HiFi (`--profile hifi`)** requires `ccs` (PacBio `pbccs`). Historically `ccs`
had no macOS-ARM64 build and this profile was Linux-x86-only; as of the current
bioconda index `pbccs 4.0.0` does resolve and install on osx-arm64, so
`conda env create` succeeds there. That only means the env builds — the HiFi
path has not been run end-to-end on Apple Silicon here. On macOS ARM64, prefer
`--profile ont-denoised`.

## One-time setup

```bash
cd /Users/suhang/Analysis/Himito/lineage_eval/lineage_sim
./setup_env.sh                # creates conda env himito-eval, fetches pbsim3 models
conda activate himito-eval
# PBSIM_MODEL_DIR defaults to ./pbsim3_models; only export if models live elsewhere
# export PBSIM_MODEL_DIR="$PWD/pbsim3_models"
```

`setup_env.sh` clones pbsim3 just for its `data/*.model` files and verifies
that `pbsim`, `ccs`, `minimap2`, `samtools`, and Python (`numpy`, `dendropy`)
all resolve inside the env.

## Detection ceiling and tuning guidance

**Read this before interpreting metrics.**

### The three frequency gates (read this first)

The pipeline applies **three independent frequency floors**, and a truth variant
is destroyed by whichever is highest. A truth variant's HF is the summed
frequency of every clone carrying it — at `--n-mutations 10 --sim-min-hf 0.05`
that is roughly 0.06–0.23, and it falls as `--n-mutations` rises.

| gate | default | set with | stage |
|------|---------|----------|-------|
| `denoise --vaf` | follows `--vaf` (0.01) unless set | `--denoise-vaf` | before the graph |
| `call -v` | tracks `--min-hf` (0.01) | `--min-hf` / `--vaf` | caller |
| `call -f` | **0.2** (`call::resolve_thresholds`) | `--frequency-threshold` | permutation test |

> **INVARIANT: every one of these must sit strictly below `--sim-min-hf`.**
> Any gate above it silently truncates the truth set and caps `var_recall`
> **regardless of depth** — these are frequency floors, not coverage floors, so
> adding reads cannot compensate.

Both defaults were tuned for `--n-mutations 10`, where min truth HF is ~0.05, and
both break silently outside that.

**`call -f` = 0.2.** Measured on `--n-mutations 10 --total-depth 1000 --seed 1`,
`ont-denoised`:

| `call -f` | VCF records | PASS SNVs | truth SNVs recovered |
|-----------|-------------|-----------|----------------------|
| **0.2 (default)** | 0 | 0 | **0 / 10** |
| 0.1  | 3   | 3  | 3 / 10 |
| **0.05** | 306 | 10 | **10 / 10** |
| 0.02 | 593 | 10 | 10 / 10 |

Here all 10 truth SNVs *are* in `<prefix>.raw_matrix.csv` at every setting — the
graph finds them and the permutation gate throws them away. The symptom is an
empty VCF, then `Himito lineage` aborting with:

```
Error running lineage analysis: No informative variants remain after filtering.
```

**`denoise --vaf`.** This used to be hardcoded to 0.03 in `run_himito.sh` — the
sweep optimum over the n=10 depth>=300 cells, valid only where the truth has
headroom above it. It is now unset by default and follows `--vaf`, and the
resolved value is echoed to stderr at run time so it is never silent. The
measurement that motivated the change, with the old 0.03 in force: Measured at `--n-mutations 15 --sim-min-hf 0.03` (min truth HF 0.0330),
one seed, `-f 0.02` so the permutation gate is not the constraint:

| `--denoise-vaf` | truth recovered | var_recall |
|-----------------|-----------------|------------|
| **0.03 (default)** | 8 / 15 | **0.533** |
| 0.01 | 13 / 15 | 0.867 |
| 0.005 | 14 / 15 | 0.933 |

Unlike the `-f` case these variants are **absent from `raw_matrix.csv` entirely**
— denoise removes them before the graph is built, so nothing downstream can
recover them. Raising `--total-depth` from 1000 to 3000 did not move
`var_recall` at all.

**Working settings.** `--n-mutations 10` at `--sim-min-hf 0.05` clears both
defaults except `-f`, so it needs only `--frequency-threshold 0.05`. Anything
rarer needs both:

```bash
# n=15, min truth HF 0.033 -> every gate must be below it
--sim-min-hf 0.03 --denoise-vaf 0.01 --frequency-threshold 0.02
```

### Use `ont-denoised`, not `ont-r10`

`run_himito.sh` runs `Himito denoise` only when the profile is `ont-denoised`;
`ont-r10` skips it. Without indel denoising the graph is swamped (~108k edges
vs ~26k on the same reads). Before the 2026-10-11 adaptive permutation fix
(ONT-recall T1), that inflated candidate set emptied the VCF via the BH
p-value floor. After T1 the VCF is no longer empty, but lineage metrics still
favour denoise on the default pbsim accuracy (~89% identity). Measured on the
frozen standard config (30 seeds, Himito at `adaa96a`):

| profile | var_precision | var_recall | ad_f1 | hap_f1 | mean edges |
|---|---|---|---|---|---|
| ont-r10 | 0.965 | 0.937 | 0.756 | 0.292 | 107,637 |
| ont-denoised | 0.997 | 1.000 | 0.994 | 0.972 | 26,001 |

`simulate_reads.sh` simulates both profiles identically, so `ont-denoised`
changes only the Himito path, not the data. Keep `ont-denoised` as the
simulation default for lineage metrics; see also Needs-a-decision #1 in
`dev_docs/superpowers/plans/2026-10-11-ont-recall-high-coverage.md`. The
committed benchmarks were produced under `ont-denoised`.

With read-linkage rescue on (ONT-recall R3-T2, default for ONT quick-start),
`ont-r10` is the recommended discovery profile for real ONT data: on the
full-depth HapMap merge it recovers low-HF private alleles that denoise
erases, without needing Phase 2 denoise linkage. Keep `ont-denoised` for
simulation lineage work (T2 table above).

### Low `var_recall` is usually a gate, not the graph

Earlier revisions of this file attributed missing low-frequency variants to GFA
edges failing to form. That is not what the measurements show on the
`ont-denoised` path: the losses tracked `denoise --vaf` exactly, and lowering it
recovered them (8/15 → 14/15 above) with depth held constant. Suspect the gates
first.

Still true, and independent of the gates: clone frequency falls as
`--n-mutations` rises, because total clone mass is conserved. At
`--n-mutations 12 --total-depth 300` a smoke run produced `var_recall ≈ 0.33`,
`n_shared = 3`, and a degenerate tree.

- Check `var_recall` and `n_shared` **first**. If `n_shared` is well under
  `n_truth_vars`, fix that first: `ad_recall` is scored against all truth pairs,
  so undetected variants cap it directly.
- Verify each gate sits below `--sim-min-hf` before blaming depth.
- Raise `--total-depth` for genuine coverage limits — but note it cannot fix a
  frequency gate.
- To make the *tree* harder without hurting detection, use `--topology` and
  `--internal-keep` rather than `--n-mutations`. See
  [Difficulty knobs](#difficulty-knobs).

### Difficulty knobs

**Do not use `--n-mutations` as a difficulty knob.** Total clone mass is
conserved, so more mutations means smaller clones, lower HF and worse detection —
it trades detection difficulty for tree difficulty rather than adding the latter.
Measured at depth 1000, 10 seeds each:

| config | var_recall | truth pairs | ad_f1 | imperfect |
|--------|-----------|-------------|-------|-----------|
| n=10 random | 0.920 | 7.3 | 0.867 | 7/10 |
| n=12 random | 0.858 | 7.0 | 0.766 | 9/10 |
| n=15 random, gates cleared | 0.933 | 12.0 | 0.903 | 9/10 |
| n=15 random, depth 3000 | 0.987 | 12.0 | **0.996** | **1/10** |

The last row is the trap: give n=15 enough depth and the tree is essentially
always reconstructed correctly. **Bigger trees are not harder trees.**

Two knobs vary tree difficulty at fixed detection:

`--topology {random,chain,star}` (`simulate_tree.py`, forwarded by `run_eval.sh`
and `sweep_seeds.sh`) sets tree *shape*. A `chain` is one unary path, which is
where ordering is genuinely hard — ordering evidence is only ever a handful of
"parent without child" reads along a path; a `star` has no ancestral pairs at
all and is a control, not a difficulty setting. A chain is also the **gentlest
shape on frequency**: with one child per node no mass is split between siblings,
so `cum_freq` decays only by `internal_keep` per level.

`--internal-keep F` (default 0.20) sets how much terminal mass each internal node
retains. This, not shape, is what makes an individual ordering call hard — see
`DEFAULT_INTERNAL_KEEP` in `simulate_tree.py`. Lowering it passes *more* mass
downward, so it raises difficulty while improving detectability.

Measured, n=10 / depth 1000 / 10 seeds:

| config | var_prec | var_recall | truth pairs | ad_f1 | ad range | imperfect |
|--------|----------|-----------|-------------|-------|----------|-----------|
| `random`, keep .20 | 1.000 | 0.920 | 7.3 | 0.867 | 0.500–1.000 | 7/10 |
| `chain`, keep .20 | 0.982 | **0.990** | **45.0** | 0.987 | 0.889–1.000 | 2/10 |
| `chain`, keep .10 | 0.982 | **0.990** | **45.0** | 0.980 | 0.889–1.000 | **4/10** |

`chain` multiplies ancestral pairs 6.9× *and* improves `var_recall` — the
decoupling `--n-mutations` cannot deliver. But shape alone made `ad_f1` go **up**
(0.867 → 0.987) — partly because a chain has no starved intermediate clones, and
partly because its near-complete detection is now rewarded rather than punished. Combine it
with `--internal-keep 0.10` for headroom at clean detection.

Two caveats. With ~44 pairs a single order flip costs ~2% of `ad_f1` instead of
~12%, so `chain` gives finer resolution over a compressed range; pair it with a
raw discordant-pair count if you want visible contrast. And `chain` frequencies
are fully deterministic — seeds vary SNV positions only, not tree shape or clone
masses — so a chain that does not fit the HF band fails immediately with the
feasible `--n-mutations` rather than retrying 100 seeds.

## Run one evaluation

```bash
# ONT (works everywhere). --frequency-threshold is NOT optional here: without
# it `call -f` defaults to 0.2 and the VCF comes back empty.
./run_eval.sh --outdir /tmp/eval_ont --profile ont-denoised \
  --n-mutations 10 --total-depth 1000 --seed 1 \
  --frequency-threshold 0.05

# HiFi (Linux/x86 only)
./run_eval.sh --outdir /tmp/eval_hifi --profile hifi \
  --n-mutations 10 --total-depth 1000 --seed 1 \
  --frequency-threshold 0.05
```

`run_eval.sh` defaults `REF` to `$REPO/rCRS.fasta`; pass `--ref` if your
reference lives elsewhere (e.g. `--ref ../../test_data/rCRS.fasta`).

### `run_eval.sh` flags

| flag | default | what it does |
|------|---------|--------------|
| `--profile` | `ont-r10` | use `ont-denoised`; `ont-r10` skips denoise and empties the VCF |
| `--n-mutations` / `--total-depth` / `--seed` | 12 / 300 / 1 | simulation size. Not a difficulty knob — see [Difficulty knobs](#difficulty-knobs) |
| `--frequency-threshold F` | unset (Himito 0.2) | `call -f`. **Must be below `--sim-min-hf`** |
| `--denoise-vaf F` | unset (follows `--vaf`) | denoise keep threshold. **Must be below `--sim-min-hf`**; sweep it with `sweep_denoise_vaf.sh` |
| `--sim-min-hf F` / `--sim-max-hf F` | 0.05 / 0.99 | band for *simulated truth* frequencies (distinct from `--min-hf`) |
| `--min-hf F` / `--max-hf F` | 0.01 / 0.95 | band for `call -v` and `lineage`; **not** the truth band |
| `--topology random\|chain\|star` | `random` | tree shape; `chain` for ordering difficulty |
| `--internal-keep F` | 0.20 | terminal mass per internal node; lower = harder ordering |
| `--fp F` / `--fn F` | unset | SCITE error rates. Unset uses `lineage::resolve_error_rates` (`ont-r10` 0.001/0.05, `ont-denoised` 0.0001/0.01, `hifi` 0.005/0.05); `metrics.tsv` records the resolved values |
| `--ref FASTA` | `$REPO/rCRS.fasta` | reference |

`--sim-min-hf` and `--min-hf` are easy to confuse: the first constrains what is
*simulated*, the second filters what was *called*.

### Call-tuning defaults in `run_himito.sh`

The wrapper passes sim-appropriate `Himito call` flags by default:

| Flag | Default | Why |
|------|---------|-----|
| `--minimal-ac` | 2 | Same as Himito `call` / quick-start. |
| `--vaf` | lineage `--min-hf` (0.01) | Call floor tracks the lineage HF band unless `--vaf` is set. |
| `--frequency-threshold` | unset (Himito 0.2) | Permutation-test heteroplasmy gate (`call -f`). **Leaving this unset empties the VCF on simulated data** — see [The `-f` gate](#the-three-frequency-gates-read-this-first). Pass `0.05`. `--p-value-threshold` is never set, so `call::resolve_thresholds` supplies the data-type p-value (0.01). |

These can be overridden when calling `run_himito.sh` directly:
```bash
./run_himito.sh --outdir /tmp/eval_ont --profile ont-denoised \
  --frequency-threshold 0.05
```

## The standard benchmark (`bench_standard.sh`)

**Use this to produce a number anyone can compare against.** It runs a frozen
configuration and scores the result against a committed baseline.

```bash
./bench_standard.sh --outdir /tmp/std --ref ../../test_data/rCRS.fasta
```

It prints each focus metric next to `baseline_standard.tsv`, with a `±2 SE` noise
band, and labels every metric `same` / `IMPROVED` / `REGRESSED`. A delta inside
the band is seed scatter, not a result. To re-score an existing run without
re-running the pipeline: `--compare path/to/seed_metrics.tsv`. To move the
baseline after an intentional change: `--update-baseline`.

### The frozen parameters, and why each one

| parameter | value | why not something else |
|-----------|-------|------------------------|
| `--profile` | `ont-denoised` | `ont-r10` skips denoise and yields an empty VCF; `hifi` is untested here |
| `--n-mutations` | 10 | **not** a difficulty knob — raising it only costs variant recall. 10 keeps min truth HF at 0.091, clear of every gate |
| `--total-depth` | 1000 | 3000 removes nearly all headroom (`ad_f1` 0.996, 1/10 imperfect); 300 wrecks detection (`var_recall` 0.750) |
| `--topology` | `chain` | 45 truth ancestral pairs instead of ~7, **and** better detection (0.990 vs 0.920) — a chain never splits clone mass between siblings |
| `--internal-keep` | 0.10 | the real ordering-difficulty knob. Chain shape alone gives no headroom (`ad_f1` 0.987) |
| `--frequency-threshold` | 0.05 | `call -f` defaults to 0.2, above every truth variant here. Must stay under min truth HF (0.296) |
| `--denoise-vaf` | 0.03 | pinned, not inherited: `run_himito.sh` no longer hardcodes it, so the benchmark states it or the result drifts with an unrelated default. ~10× headroom here, so it truncates nothing |
| seeds | 1–30 | 10 seeds cannot separate a real effect from one coin-flip edge |

Changing any of these makes the result incomparable to the baseline — that is
the point of freezing them. For exploration, call `sweep_seeds.sh` directly.

### Committed baseline

`baseline_standard.tsv`, 30/30 seeds, no failures. Mean `n_truth_pairs` = 45.0,
mean `n_shared` = 9.8 of 10.

| metric | mean | sd | min | max |
|--------|------|----|-----|-----|
| `var_precision` | 0.9909 | 0.0277 | 0.9091 | 1.000 |
| `var_recall` | 0.9833 | 0.0461 | 0.8000 | 1.000 |
| `var_f1` | 0.9863 | 0.0274 | 0.8889 | 1.000 |
| `ad_precision` | 0.9932 | 0.0122 | 0.9556 | 1.000 |
| `ad_recall` | 0.9607 | 0.0892 | 0.6222 | 1.000 |
| `ad_f1` | **0.9744** | 0.0534 | 0.7671 | 1.000 |
| `hap_precision` | 0.9215 | 0.1399 | 0.3750 | 1.000 |
| `hap_recall` | 0.9200 | 0.1584 | 0.3000 | 1.000 |
| `hap_f1` | **0.9202** | 0.1489 | 0.3333 | 1.000 |

Headroom: 11/30 seeds imperfect on `ad_f1`, 13/30 on `hap_f1`, 7/30 on `var_f1`.
Detection is clean (`var_precision` 0.991) without being saturated, so all three
families can still move in either direction.

**Resolution.** The `2*SE` band is what a run has to beat to count as a change:
roughly **±0.028 on `ad_f1`** and **±0.077 on `hap_f1`** at 30 seeds. `hap_*` is
much noisier because a single mis-grouped clone moves it by ~0.1; treat small
`hap_*` deltas with suspicion and raise the seed count if you need to resolve
them.

## Homopolymer length heteroplasmy

`bench_standard.sh` stays SNV-only on purpose. Homopolymer length variation,
which dominates real PacBio/ONT disagreement, is a separate frozen config:
`bench_hp.sh`. The three simulator modes below compose. Leaving all three unset
reproduces the SNV-only truth the standard benchmark scores.

| mode | flag | what it produces | scored as |
|------|------|------------------|-----------|
| Background homoplasmies | `--background` | Fixed SNVs in every clone, including the reference. `none` (default) applies nothing. `hvs1-polyC` is `T16189C,A16183C`. Any other value is a comma-separated `<REF><pos><ALT>` list. | `truth/background.txt`. Dropped from the `var_*` denominators. They move tract anchors. |
| Lineage-borne length change | `--hp-lineage-indels K` | K of the N tree edges are inherited length changes. Unset omits the flag (K = 0). | Ordinary tree nodes, under their canonical names, in `var_*` / `ad_*` / `hap_*` and in `hp_var_*`. |
| Lineage-independent length heteroplasmy | `--hp-unstable` | Each clone is split into sub-clones by tract-length state. `common`, or `TRACT:delta=prob,...;TRACT:...`. | `truth/hp_tracts.tsv`. Not a tree edge. `hp_intrusion` and `hp_len_l1`. Dropped from `var_*`. |

`common` is D310 `{0: 0.55, +1: 0.35, +2: 0.10}`, HVS1 `{0: 0.6, +1: 0.3, +2: 0.1}`, C573 `{0: 0.7, +1: 0.2, +2: 0.1}`. Deltas are relative to the catalogue run length.

### Catalogue (`hp_tracts.py`)

`CATALOGUE` is the rCRS tract list, 1-based inclusive. `anchor` is the base before the run: the position a left-aligned insertion uses. `in_himito` is set only for the three tracts already in `src/tracts.rs` (`D310`, `HVS1`, `CA514`). `R8281` is in the catalogue with simulation off.

| id | run | unit × copies | anchor |
|----|-----|---------------|--------|
| D310 | 303–309 | C × 7 | 302 |
| HVS1 | 16184–16188 | C × 5 | 16183 |
| C573 | 568–573 | C × 6 | 567 |
| T455 | 452–455 | T × 4 | 451 |
| C960 | 956–960 | C × 5 | 955 |
| C5899 | 5895–5899 | C × 5 | 5894 |
| A12418 | 12418–12425 | A × 8 | 12417 |
| CA514 | 514–523 | CA × 5 | 513 |
| R8281 | 8272–8289 | CCCCCTCTA × 2 | 8271 |

### Canonical names

`canonical_name` in `hp_tracts.py` rewrites an allele whose REF span sits in a tract window to `m.<anchor><rCRS window>><background + edit>`. The REF side of that name is rCRS over the window, not the background bases. On a `T16189C` background, `m.16189T>TC` and `m.16183A>AC` are one name. `score_lineage.py` canonicalises both trees, the truth list, the VCF, and matrix row names before comparing, so two spellings of one allele are not a mismatch. A name outside every tract, including a typical SNV, is unchanged.

### HP metrics

Appended after `n_background_calls`. Read them by header name. `bench_hp.sh` averages `hp_intrusion`, `hp_backbone_nodes`, `hp_len_l1`, `root_fraction`, and `hp_var_f1` with the standard nine; `NA` cells are left out of the mean.

| column | meaning |
|--------|---------|
| `hp_var_precision` / `hp_var_recall` / `hp_var_f1` | Detection on the `variant` column of `truth/hp_lineage.tsv` only. `NA` when that file is absent. |
| `hp_intrusion` | Among truth SNVs that also appear in the reconstructed tree, the fraction whose reconstructed ancestors include an unstable length allele. |
| `hp_backbone_nodes` | Count of unstable names that are the parent of at least one reconstructed node. |
| `hp_len_l1` | Mean over tracts of Σ \|truth frequency − reconstructed frequency\| across length alleles. `NA` when there is no non-zero δ row or no matrix. |
| `root_fraction` | ROOT `n_reads_attached` divided by the sum of that column. A `ROOT_BLOCK` row is an annotation and is left out of the sum. `NA` when the column is absent. |

### Frozen HP benchmark (`bench_hp.sh`)

Same driver shape as `bench_standard.sh`: `--seeds` for a smoke run (breaks comparability), `--compare FILE` to re-score an existing `seed_metrics.tsv`, and `--update-baseline` to copy the result onto `baseline_hp.tsv`. A delta inside the 2·SE band is `same`. A higher `var_*`, `ad_*`, `hap_*`, or `hp_var_f1` is `IMPROVED`. A lower `hp_intrusion`, `hp_backbone_nodes`, `hp_len_l1`, or `root_fraction` is `IMPROVED`. `baseline_hp.tsv` starts with no data rows. Until it has some, compare prints `no baseline yet` and exits 0. Score HP runs against `baseline_hp.tsv`.

```bash
./bench_hp.sh --outdir /tmp/hp --ref ../../test_data/rCRS.fasta
```

Frozen on the standard call/tree settings (`ont-denoised`, 10 mutations, depth 1000, `chain`, internal-keep 0.10, `-f` 0.05, denoise-vaf 0.03, seeds 1–30):

| parameter | value |
|-----------|-------|
| `--background` | `hvs1-polyC` |
| `--hp-lineage-indels` | 2 |
| `--hp-unstable` | `common` |

### In-silico PacBio/ONT pair (`run_eval.sh --pair`)

One `simulate_tree.py` call, then reads from that same `truth/` twice: profile `hifi`, then `ont-denoised`, same seed. Reconstructions are moved to `himito_pb/` and `himito_ont/` before the next run writes `himito/`, and scored to `metrics_pb.tsv` and `metrics_ont.tsv`. `pair_metrics.tsv` has two columns:

| column | meaning |
|--------|---------|
| `pair_ad_f1` | `ad_f1` of the ONT mutation tree scored with the PacBio mutation tree as the truth tree. PacBio VCF calls are the truth variant set. |
| `pair_root_fraction_delta` | ONT `root_fraction` minus PacBio `root_fraction`. |

This is the in-silico twin of the real sample pairs. `--pair` needs `ccs` (HiFi) and is Linux-only for now. If `ccs` is not on `PATH`, the script exits and does not report an ONT-only run as a pair. `--profile` is ignored while `--pair` is set.

### Sweep HP rates (`sweep_hp_rates.sh`)

Call and the read simulation do not depend on `--hp-fp-rate` or `--hp-fn-rate`. `run_eval.sh --stop-after call` (or `sweep_seeds.sh --keep-intermediate`) writes `himito/sim.matrix.csv`, `himito/sim.vcf`, and `truth/` once per seed. `sweep_hp_rates.sh` re-runs only `Himito lineage` and `score_lineage.py` on that matrix. The global `--fp-rate` / `--fn-rate` stay at the profile default unless `--fp` / `--fn` are set. The default grid is `{1,2,4,8,16}×fp` and `{1,2,4,8,16,32}×fn`, each product capped at 0.5. For `ont-denoised` that is fp `0.0001…0.0016` and fn `0.01…0.32`.

```bash
./sweep_hp_rates.sh --outdir /tmp/hp_sweep --profile ont-denoised \
  --seeds "1 2 3 4 5 6 7 8 9 10" --collapse-hp
```

`pick_best_hp_rates.py` writes `hp_sweep_summary.tsv` (mean and SE over seeds), `best_hp_rates.json`, and `hp_sweep_report.md`. The selection rule, in order:

1. Primary objective (default `hp_intrusion`, minimised; `ad_f1`, `hap_f1`, and `hp_len_l1` are the other choices). A cell whose mean is within 2·SE of the best mean is tied with it. That is the same noise band `bench_standard.sh` uses.
2. Among tied cells, secondary objective (default `ad_f1`, maximised), same tie rule.
3. Among cells still tied, the one closest to the global rates: smallest `|log2(hp_fn/fn)| + |log2(hp_fp/fp)|`. Do not move a parameter the data cannot distinguish.
4. If the primary objective's best and worst cells are themselves within 2·SE, the grid is flat. `best` is null and the report says why.

A flat grid means leave the default alone. It does not mean pick the top row.

The report also gives `root_fraction` at the chosen cell and at the global-rate cell. That is the PacBio/ONT symptom these rates are meant to move.

### Workflows

`wdl/Himito-lineage/Himito_lineage_sim.wdl` is the full simulate → call → lineage → score grid (seeds × n_mutations × depth), including the HP modes. `wdl/Himito-lineage/Himito_lineage_hp_sweep.wdl` imports that file's `ValidateConfig` and `CombineMetrics`. It runs `SimulateCell` once per seed (`--stop-after call`), then `LineageScoreCell` for every `hp_fp_grid` × `hp_fn_grid` pair, then the picker above. `n_mutations` and `depth` are scalars. `ValidateConfig` fails the sweep unless `hp_lineage_indels > 0` or `hp_unstable` is set: otherwise the rates have nothing to act on and the sweep is guaranteed flat.

`lineage_eval/lineage_sim/Himito_lineage_hp_sweep.inputs.json` is the `bench_hp.sh` configuration with `collapse_hp` true and the ont-denoised default grid. It uses seeds 1–10. Thirty seeds times that 30-cell grid is too expensive for a first pass. MCMC in the sweep is 2 chains × 2000 iterations, the `sweep_fpfn.sh` budget. The image must be rebuilt after this plan so `/opt/lineage_sim` contains `pick_best_hp_rates.py` and `run_eval.sh --stop-after call`.

## Multi-seed baselines (`sweep_seeds.sh`)

A single seed cannot separate "the method has a systematic weakness" from "this
instance contained one coin-flip edge". `sweep_seeds.sh` runs one full
`run_eval.sh` per seed and aggregates into `seed_metrics.tsv` with a `seed`
column, tolerating per-seed failures.

```bash
./sweep_seeds.sh --outdir /tmp/ms --profile ont-denoised \
  --n-mutations 10 --total-depth 1000 --seeds "$(seq -s' ' 1 30)" \
  --frequency-threshold 0.05 --ref ../../test_data/rCRS.fasta
```

It forwards `--sim-min-hf`, `--denoise-vaf`, `--topology` and `--internal-keep`,
and prints a mean/min/max summary. Reads, BAMs and GFAs are deleted per seed
unless `--keep-intermediate` is passed (~30 MB/seed otherwise).

Measured baselines, `ont-denoised`, `-f 0.05`:

| config | seeds | var_prec | var_recall | ad_f1 | ad range | hap_f1 | imperfect `ad` |
|--------|-------|----------|-----------|-------|----------|--------|----------------|
| n=10, depth 1000 | 10 | 1.000 | 0.920 | 0.867 | 0.500–1.000 | 0.936 | 7/10 |
| n=10, depth 300 | 30 | 0.749 | 0.750 | 0.615 | 0.286–0.857 | 0.621 | 30/30 |
| n=12, depth 1000 | 10 | 0.991 | 0.858 | 0.766 | 0.500–1.000 | 0.854 | 9/10 |
| n=10, depth 1000, `chain` keep .10 | 10 | 0.982 | 0.990 | **0.980** | 0.889–1.000 | 0.935 | 4/10 |

Scored with the corrected `ad_recall` denominator; see
[Why `ad_recall` is scored against *all* truth pairs](#why-ad_recall-is-scored-against-all-truth-pairs).

Outputs land under `<outdir>/`: `truth/`, `reads/`, `himito/`, and
`metrics.tsv`.


## Sweep denoise --vaf

```bash
./sweep_denoise_vaf.sh --outdir /tmp/dv --ref ../../test_data/rCRS.fasta \
  --vafs "0.005 0.01 0.03 0.05" --seeds "1 2 3 4 5 6 7 8 9 10"
```

Defaults to the standard config and forwards `--topology` / `--internal-keep` /
`--sim-min-hf`, so a cell differs from the frozen benchmark only in `--vaf`.

Measured on the standard config (10 seeds/cell, min truth HF **0.296**, so a
~10x margin over the gate):

| `--vaf` | var_precision | var_recall | ad_f1 | hap_f1 |
|---------|---------------|------------|-------|--------|
| 0.005 | 0.8957 | **0.9900** | 0.9800 | 0.8912 |
| 0.01 | 0.9496 | **0.9900** | 0.9800 | 0.9174 |
| **0.03** | **0.9818** | **0.9900** | 0.9800 | **0.9352** |
| 0.05 | 0.9818 | **0.9900** | 0.9800 | 0.9352 |

`var_recall` is flat: with that much headroom the gate truncates nothing, and the
precision gain up to 0.03 is free. `ad_f1` does not move at all — the threshold
changes which noise reaches the matrix, not the tree. 0.03 is the optimum here
and is what `bench_standard.sh` pins.

**This is also why 0.03 cannot be a default.** The same value in a low-margin
config destroys the truth set: at `--n-mutations 15 --sim-min-hf 0.03` (min truth
HF 0.033, ~1.1x margin) it cut `var_recall` to 0.533, and 0.01 restored it to
0.867 on identical reads. The optimum is a property of the config, not of the
tool — sweep the regime you are actually in, and check `var_recall` before
reading any precision gain.

## Sweep fp/fn

After `run_eval.sh` completes (matrix and VCF are fixed), sweep SCITE error
rates by re-running only `Himito lineage` per grid cell:

```bash
./sweep_fpfn.sh --outdir /tmp/eval_ont --profile ont-denoised \
  --fp-grid "0.0005 0.001 0.005 0.01" --fn-grid "0.02 0.05 0.1 0.2"
```

Writes `sweep_metrics.tsv` and prints the best cell by `ad_f1` (tie-break:
`var_f1`).

**Check the spread before trusting that line.** When detection is clean the
matrix is unambiguous and the SCITE error priors have nothing to do: on
`n=10 / depth 1000 / seed 1 / ont-denoised` all 16 cells returned identical
scores (`var_f1` 1.000, `ad_f1` 0.875, `hap_f1` 0.900), so the reported "best
cell" was just the first row of a 16-way tie, not an optimum. To rank fp/fn you
need instances hard enough to separate the cells — lower `--total-depth`, more
mutations, and several seeds averaged.

## Metrics (in `metrics.tsv` / `sweep_metrics.tsv`)

Three metric families, each with its own denominators. Precision = "of what we
called, how much is right"; recall/sensitivity = "of the truth, how much we
recovered".

| Column | Meaning |
|--------|---------|
| `var_precision` / `recall` / `f1` | Variant *detection*: PASS SNVs in `sim.vcf` vs the truth SNV set. Independent of the tree. |
| `ad_precision` / `recall` / `f1` | **Ancestor–descendant** accuracy: are truth ancestral pairs preserved in the reconstruction? `ad_recall` is scored against *every* truth pair, `ad_precision` only over shared variants — see [below](#why-ad_recall-is-scored-against-all-truth-pairs). The primary tree-topology metric. |
| `hap_precision` / `recall` / `f1` | **Clone recovery**: a reconstructed haplotype matches a truth clone when their full variant sets are equal. Needs `--recon-matrix` + `--truth-clones`; reads `NA` otherwise. |
| `n_truth_vars` / `n_detected_vars` | Denominators for `var_recall` / `var_precision`. |
| `n_shared` | Variants present in both the truth tree and the reconstructed tree (the `ad_precision` taxon set). |
| `n_truth_pairs` | Ancestral pairs in the truth tree — the `ad_recall` denominator. Read `ad_*` next to this. |
| `n_truth_clones` / `n_recon_haps` | Denominators for `hap_recall` / `hap_precision`. Excludes the mutation-free bin. |

`ad_*` counts ordered (ancestor, descendant) pairs with ROOT excluded — ROOT is an ancestor of everything in both trees, so
counting it would hand every run |n_shared| free true positives. When both trees
are flat over the shared set they agree, and `ad_*` reads 1.0; when only one is
flat it reads 0.0.

Column order in the TSV: `profile`, `fp`, `fn`, `n_truth_vars`,
`n_detected_vars`, `n_shared`, `n_truth_pairs`, `var_precision`, `var_recall`,
`var_f1`,
`ad_precision`, `ad_recall`, `ad_f1`, `n_truth_clones`, `n_recon_haps`,
`hap_precision`, `hap_recall`, `hap_f1`, then `n_background_calls` and the
homopolymer columns. Read columns by header name, not position. Definitions
for those columns are in
[Homopolymer length heteroplasmy](#homopolymer-length-heteroplasmy).

### Why `ad_recall` is scored against *all* truth pairs

`ad_recall`'s denominator is every ancestral pair in the truth tree
(`n_truth_pairs`), **not** only pairs among `shared`. `ad_precision` stays
restricted to `shared`.

This asymmetry is deliberate and load-bearing. When both were restricted to
`shared`, a variant the caller missed left the denominator entirely, so a run
that detected *less* scored *higher* — the reconstruction was never charged for
ancestry it had no chance to recover. Correlation of `ad_f1` with `var_recall`,
before and after, over 100 runs:

| config | before | after |
|--------|--------|-------|
| n=10 depth 1000 | −0.459 | **+0.679** |
| n=10 depth 300 | −0.403 | **+0.271** |
| n=15 gates clear | −0.625 | **+0.723** |
| n=10 chain keep .20 | −0.111 | **+0.980** |
| **pooled (n=100)** | +0.049 | **+0.863** |

The clearest case: `--n-mutations 15` with `denoise --vaf` destroying 40% of the
truth set scored a **perfect `ad_f1` = 1.000 on all 10 seeds** under the old
definition — the best score in the whole benchmark, from the worst-detecting
config. It now scores 0.374, the worst.

A perfectly correct tree built from fewer detected variants now degrades as it
should (4-mutation chain truth, reconstruction correct on everything it saw):

| variants detected | `ad_precision` | `ad_recall` | `ad_f1` |
|-------------------|----------------|-------------|---------|
| 4 / 4 | 1.000 | 1.000 | 1.000 |
| 3 / 4 | 1.000 | 0.500 | 0.667 |
| 2 / 4 | 1.000 | 0.167 | 0.286 |

All four rows scored 1.000 before. Precision stays 1.000 throughout because the
tree really is correct — a *false-positive* variant is a detection error that
`var_precision` already reports, so charging it here too would double-count one
mistake.

**`ad_*` from before this change is not comparable to `ad_*` after it.** Earlier
numbers were inflated by however much detection was incomplete: near-perfect
detection barely moves (n=15 depth 3000, `var_recall` 0.987: 1.000 → 0.996),
poor detection moves a lot (n=10 depth 300, `var_recall` 0.750: 0.848 → 0.615).

## Suggested experiments

- **Run `bench_standard.sh` first.** It is the frozen config and it reports
  against the committed baseline with a noise band. Only reach for
  `sweep_seeds.sh` when you are deliberately exploring off-config.
- **Ordering difficulty:** vary `--internal-keep` (0.20/0.15/0.10) under
  `--topology chain`. This is the axis that isolates tree inference.
- **Depth titration:** vary `--total-depth` (300/1000/3000) — but confirm every
  frequency gate is below `--sim-min-hf` first, or you are titrating a gate.
- **Tree size:** vary `--n-mutations`, understanding it measures *detection*
  scaling, not tree difficulty.

Two things that did **not** work, so you can skip them:

- **`sweep_fpfn.sh` does not move `ad_*`.** Across a 16-cell grid plus six
  hand-picked rates — including the empirically measured fp/fn — every cell
  returned identical scores. Error-rate priors are not the lever on
  clean-detection instances. (The empirical rates, measured against
  truth-labelled reads, are fp ≈ 0.009 and fn ≈ 0.21, versus the `ont-denoised`
  priors of 0.0001/0.01 — mis-specified by 90× and 21×, but correcting them
  changed nothing.)
- **Raising `--n-mutations` does not make the tree harder.** See
  [Difficulty knobs](#difficulty-knobs).

## Notes / limitations

- **HiFi is untested on Apple Silicon.** `pbccs 4.0.0` now installs from
  bioconda on osx-arm64, but the HiFi path has not been run end-to-end there;
  use `--profile ont-denoised` for local development on Apple Silicon.
  `run_eval.sh --pair` uses that HiFi path, so `--pair` is Linux-only for now
  and exits when `ccs` is missing.
- **Three frequency gates must all sit below `--sim-min-hf`.** `call -f` (0.2),
  `denoise --vaf` (0.03) and `call -v`. The first two defaults were tuned for
  `--n-mutations 10` and break silently otherwise; no amount of depth
  compensates. See [The three frequency gates](#the-three-frequency-gates-read-this-first).
- **`ad_*` changed definition.** `ad_recall` is now scored against all truth
  pairs, so numbers from before that change are inflated and not comparable; see
  [Why `ad_recall` is scored against *all* truth pairs](#why-ad_recall-is-scored-against-all-truth-pairs).
- **`bench_standard.sh` is SNV-only on purpose.** It still uses `AVOID_RANGES`
  in `simulate_tree.py`, so the committed baseline stays on those sites.
  Length heteroplasmy is the separate frozen config in
  [Homopolymer length heteroplasmy](#homopolymer-length-heteroplasmy).
- **Very low depth can yield zero HiFi reads.** Raise `--total-depth` if `ccs`
  produces empty output.
- **Read-level truth labels** are embedded in FASTQ headers
  (`@clone_<id>_rot<n>_...`) for debugging; Himito ignores them.
- **Frequency assignment can fail.** If `simulate_tree.py` raises
  "Could not find valid frequency assignment", reduce `--n-mutations`, lower
  `--sim-min-hf`, or try a different `--seed`. At `--sim-min-hf 0.05` the
  practical ceiling is `--n-mutations 12`; n=15 needs ≤0.03 and n=20 needs ≤0.02.
  Under `--topology chain` frequencies are deterministic, so reseeding cannot
  help — that path fails immediately and reports the feasible `--n-mutations`.
