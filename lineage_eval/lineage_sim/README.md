# Himito-lineage simulation evaluation

Simulate a ground-truth mitochondrial clonal mutation tree, generate long reads
(PacBio HiFi or ONT-R10), reconstruct the tree with `Himito call` + `Himito
lineage`, and score reconstruction accuracy across replicate seeds.

**Start here:** [The three frequency gates](#the-three-frequency-gates-read-this-first)
— the defaults silently empty the VCF on simulated data — and
[`ad_*` is anti-correlated with detection](#ad_-is-anti-correlated-with-detection),
which governs how the tree metrics may be compared.

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
| `denoise --vaf` | **0.03** (hardcoded in `run_himito.sh`) | `--denoise-vaf` | before the graph |
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

**`denoise --vaf` = 0.03.** `run_himito.sh:37` documents this as "the sweep
optimum over the n=10 depth>=300 cells". It bites whenever clones are rarer than
that. Measured at `--n-mutations 15 --sim-min-hf 0.03` (min truth HF 0.0330),
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
`ont-r10` skips it. Without indel denoising the graph is swamped (94,694 edges
vs 25,503 on the same reads) and the VCF comes back empty regardless of `-f`.
`simulate_reads.sh` simulates both profiles identically, so `ont-denoised`
changes only the Himito path, not the data. The committed benchmarks were
produced under `ont-denoised`.

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
  `n_truth_vars`, fix that before reading any tree metric — see
  [`ad_*` is anti-correlated with detection](#ad_-is-anti-correlated-with-detection).
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
| n=10 random | 0.920 | 6.4 | 0.975 | 2/10 |
| n=12 random | 0.858 | 10.3 | 0.915 | 3/10 |
| n=15 random, gates cleared | 0.933 | — | 0.978 | 3/10 |
| n=15 random, depth 3000 | 0.987 | — | **1.000** | **0/10** |

The last row is the trap: give n=15 enough depth and the tree is reconstructed
perfectly every time. **Bigger trees are not harder trees.**

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
| `random`, keep .20 | 1.000 | 0.920 | 6.4 | 0.975 | 0.875–1.000 | 2/10 |
| `chain`, keep .20 | 0.982 | **0.990** | **44.1** | 0.998 | 0.978–1.000 | 1/10 |
| `chain`, keep .10 | 0.982 | **0.990** | **44.1** | 0.991 | 0.956–1.000 | **3/10** |

`chain` multiplies ancestral pairs 6.9× *and* improves `var_recall` — the
decoupling `--n-mutations` cannot deliver. But shape alone made `ad_f1` go **up**
(0.975 → 0.998), because a chain has no starved intermediate clones. Combine it
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
| `--denoise-vaf F` | unset (`run_himito.sh` 0.03) | denoise keep threshold. **Must be below `--sim-min-hf`** |
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

| config | seeds | var_prec | var_recall | ad_f1 | hap_f1 | imperfect `ad` |
|--------|-------|----------|-----------|-------|--------|----------------|
| n=10, depth 1000 | 10 | 1.000 | 0.920 | 0.975 | 0.936 | 2/10 |
| n=10, depth 300 | 30 | 0.749 | 0.750 | 0.848 | 0.621 | 16/30 |
| n=12, depth 1000 | 10 | 0.991 | 0.858 | 0.915 | 0.854 | 3/10 |
| n=10, depth 1000, `chain` keep .10 | 10 | 0.982 | 0.990 | 0.991 | 0.935 | 3/10 |

Outputs land under `<outdir>/`: `truth/`, `reads/`, `himito/`, and
`metrics.tsv`.


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
| `ad_precision` / `recall` / `f1` | **Ancestor–descendant** accuracy: over shared variants, are truth ancestral pairs preserved in the reconstruction? The primary tree-topology metric. |
| `hap_precision` / `recall` / `f1` | **Clone recovery**: a reconstructed haplotype matches a truth clone when their full variant sets are equal. Needs `--recon-matrix` + `--truth-clones`; reads `NA` otherwise. |
| `n_truth_vars` / `n_detected_vars` | Denominators for `var_recall` / `var_precision`. |
| `n_shared` | Variants present in both the truth tree and the reconstructed tree (the `ad_*` taxon set). |
| `n_truth_clones` / `n_recon_haps` | Denominators for `hap_recall` / `hap_precision`. Excludes the mutation-free bin. |

`ad_*` counts ordered (ancestor, descendant) pairs over `n_shared` variants,
with ROOT excluded — ROOT is an ancestor of everything in both trees, so
counting it would hand every run |n_shared| free true positives. When both trees
are flat over the shared set they agree, and `ad_*` reads 1.0; when only one is
flat it reads 0.0.

Column order in the TSV: `profile`, `fp`, `fn`, `n_truth_vars`,
`n_detected_vars`, `n_shared`, `var_precision`, `var_recall`, `var_f1`,
`ad_precision`, `ad_recall`, `ad_f1`, `n_truth_clones`, `n_recon_haps`,
`hap_precision`, `hap_recall`, `hap_f1`. Read columns by header name, not
position.

### `ad_*` is anti-correlated with detection

`ad_*` is restricted to `shared` (variants in *both* trees), so a variant the
caller misses leaves the denominator entirely and the reconstruction is never
charged for it. Missing variants therefore produce a smaller, easier tree and a
**higher** `ad_f1`. Correlation with `var_recall` over 60 runs:

| metric | depth 300 (n=30) | depth 1000 (n=10) | pooled (n=60) |
|--------|------------------|-------------------|---------------|
| `ad_f1` | **−0.403** | **−0.459** | **−0.228** |
| `hap_f1` | +0.315 | +0.641 | +0.579 |

At depth 300 the 14 seeds scoring a perfect `ad_f1` had *worse* detection than
the 16 imperfect ones (mean `var_recall` 0.686 vs 0.806). The extreme case:
`--n-mutations 15` with `denoise --vaf` blocking 40% of truth variants scored
`ad_f1 = 1.000` on **all 10 seeds**.

Consequences:

- **A genuine improvement to variant detection will lower `ad_f1`.** Do not
  compare `ad_*` across runs with different `n_shared`.
- Always read `ad_*` next to `n_shared`. Five of the 30 depth-300 seeds scored
  `ad` on ≤2 truth pairs.
- `hap_*` is correctly signed at every depth and is the better headline metric
  when detection is incomplete.

## Suggested experiments

- **Replicates first.** Always use `sweep_seeds.sh` (10+ seeds) before drawing a
  conclusion; single-seed differences of 0.1 in `ad_f1` are one coin-flip edge.
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
- **Three frequency gates must all sit below `--sim-min-hf`.** `call -f` (0.2),
  `denoise --vaf` (0.03) and `call -v`. The first two defaults were tuned for
  `--n-mutations 10` and break silently otherwise; no amount of depth
  compensates. See [The three frequency gates](#the-three-frequency-gates-read-this-first).
- **`ad_*` rises when detection falls.** Never compare it across runs with
  different `n_shared`; see [`ad_*` is anti-correlated with detection](#ad_-is-anti-correlated-with-detection).
- **SNV-only truth (no indels).** rCRS homopolymer/control-region positions are
  excluded (see `AVOID_RANGES` in `simulate_tree.py`) to keep the benchmark on
  cleanly callable sites.
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
