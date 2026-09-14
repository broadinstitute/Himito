# Himito-lineage simulation evaluation

Simulate a ground-truth mitochondrial clonal mutation tree, generate long reads
(PacBio HiFi or ONT-R10), reconstruct the tree with `Himito call` + `Himito
lineage`, and score reconstruction accuracy. Then sweep the SCITE `fp`/`fn`
rates to find the best settings per read type.

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

### The `-f` gate silently empties the VCF

`Himito call -f` (`--frequency-threshold`, the permutation-test heteroplasmy
gate) defaults to **0.2** via `call::resolve_thresholds`, and `run_himito.sh`
does not override it. Simulated truth variants land well below that: a truth
variant's HF is the summed frequency of every clone carrying it, which at
`--n-mutations 10` is roughly **0.06–0.23**. The default gate therefore discards
nearly all of them.

Measured on `--n-mutations 10 --total-depth 1000 --seed 1`, profile
`ont-denoised`, holding everything else at the wrapper defaults:

| `call -f` | VCF records | PASS SNVs | truth SNVs recovered |
|-----------|-------------|-----------|----------------------|
| **0.2 (default)** | 0 | 0 | **0 / 10** |
| 0.1  | 3   | 3  | 3 / 10 |
| **0.05** | 306 | 10 | **10 / 10** |
| 0.02 | 593 | 10 | 10 / 10 |

This is *not* a detection failure: all 10 truth SNVs are present in
`<prefix>.raw_matrix.csv` at every setting. The graph finds them and the gate
throws them away. The symptom is an empty VCF, and then `Himito lineage`
aborting with:

```
Error running lineage analysis: No informative variants remain after filtering.
```

**Pass `--frequency-threshold 0.05`** for simulated data. At 0.05 the extra
records are indel artifacts that the caller marks `Potential_Artifact` (296 of
the 306 above), so they never reach `var_precision`, which scored 1.000.

### Use `ont-denoised`, not `ont-r10`

`run_himito.sh` runs `Himito denoise` only when the profile is `ont-denoised`;
`ont-r10` skips it. Without indel denoising the graph is swamped (94,694 edges
vs 25,503 on the same reads) and the VCF comes back empty regardless of `-f`.
`simulate_reads.sh` simulates both profiles identically, so `ont-denoised`
changes only the Himito path, not the data. The committed benchmarks were
produced under `ont-denoised`.

### Clone frequency vs. the GFA threshold

Independently of the gate above, low-frequency clones (~3% of reads) may not
form GFA edges at all and are then invisible to the caller. At
`--n-mutations 12 --total-depth 300` a smoke run produced `var_recall ≈ 0.33`,
`n_shared = 3`, and a degenerate tree — a real property of the graph step.

- Raise `--total-depth` (e.g. 600–1000) so low-frequency clones accumulate
  enough reads to form GFA edges.
- Reduce `--n-mutations` (e.g. 6–8) so clone frequencies are higher on average.
- Check `var_recall` and `n_shared` first; if `n_shared < 3`, the `ad_*` metrics
  are low-power regardless of fp/fn.

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

Optional flags: `--fp F` (SCITE fp rate) and `--fn F` (SCITE fn rate). Left
unset, `lineage::resolve_error_rates` supplies the per-profile defaults
(`ont-r10` 0.001/0.05, `ont-denoised` 0.0001/0.01, `hifi` 0.005/0.05) and
`metrics.tsv` records the resolved values.

Outputs land under `<outdir>/`: `truth/`, `reads/`, `himito/`, and
`metrics.tsv`.

### Call-tuning defaults in `run_himito.sh`

The wrapper passes sim-appropriate `Himito call` flags by default:

| Flag | Default | Why |
|------|---------|-----|
| `--minimal-ac` | 2 | Same as Himito `call` / quick-start. |
| `--vaf` | lineage `--min-hf` (0.01) | Call floor tracks the lineage HF band unless `--vaf` is set. |
| `--frequency-threshold` | unset (Himito 0.2) | Permutation-test heteroplasmy gate (`call -f`). **Leaving this unset empties the VCF on simulated data** — see [The `-f` gate](#the--f-gate-silently-empties-the-vcf). Pass `0.05`. `--p-value-threshold` is never set, so `call::resolve_thresholds` supplies the data-type p-value (0.01). |

These can be overridden when calling `run_himito.sh` directly:
```bash
./run_himito.sh --outdir /tmp/eval_ont --profile ont-denoised \
  --frequency-threshold 0.05
```

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

## Suggested experiments

- **Depth titration:** vary `--total-depth` (e.g. 300/600/1000) to find where
  heteroplasmy resolution breaks down.
- **Tree size:** vary `--n-mutations` (6/12/24) to test scaling.
- **Replicates:** vary `--seed` and average metrics per profile.
- **fp/fn refinement:** run `sweep_fpfn.sh` per profile; the recommended
  `--fp`/`--fn` for HiFi vs ONT-R10 are the best-cell values, averaged across
  seeds.

## Notes / limitations

- **HiFi is untested on Apple Silicon.** `pbccs 4.0.0` now installs from
  bioconda on osx-arm64, but the HiFi path has not been run end-to-end there;
  use `--profile ont-denoised` for local development on Apple Silicon.
- **`call -f` must be lowered for simulated data.** The 0.2 default drops every
  truth SNV and leaves an empty VCF; pass `--frequency-threshold 0.05`.
- **SNV-only truth (no indels).** rCRS homopolymer/control-region positions are
  excluded (see `AVOID_RANGES` in `simulate_tree.py`) to keep the benchmark on
  cleanly callable sites.
- **Very low depth can yield zero HiFi reads.** Raise `--total-depth` if `ccs`
  produces empty output.
- **Read-level truth labels** are embedded in FASTQ headers
  (`@clone_<id>_rot<n>_...`) for debugging; Himito ignores them.
- **Frequency assignment can fail.** If `simulate_tree.py` raises
  "Could not find valid frequency assignment", reduce `--n-mutations` or try a
  different `--seed`.
