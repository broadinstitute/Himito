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

**ONT-R10 (`--profile ont-r10`)** runs on any platform (Linux x86, macOS
ARM64). pbsim3's `ERRHMM-ONT-HQ` model needs no extra tools beyond the conda
env.

**HiFi (`--profile hifi`)** requires `ccs` (PacBio `pbccs`). `ccs` has no
macOS-ARM64 build; the HiFi profile must be run on a **Linux x86** machine (or
inside a `linux/amd64` container). Attempting it on macOS ARM64 will fail when
`ccs` is not available. On macOS ARM64, use `--profile ont-r10`.

## One-time setup

```bash
cd /Users/suhang/Analysis/Himito/eval/lineage_sim
./setup_env.sh                # creates conda env himito-eval, fetches pbsim3 models
conda activate himito-eval
export PBSIM_MODEL_DIR=/Users/suhang/Analysis/Himito/eval/lineage_sim/pbsim3_models
```

`setup_env.sh` clones pbsim3 just for its `data/*.model` files and verifies
that `pbsim`, `ccs`, `minimap2`, `samtools`, and Python (`numpy`, `dendropy`)
all resolve inside the env.

## Detection ceiling and tuning guidance

**Read this before interpreting metrics.** At the default simulation settings
(`--n-mutations 12`, `--total-depth 300`), Himito's graph step surfaces only a
subset of truth variants: low-frequency clones (~3% of reads) do not form GFA
edges and are invisible to the caller. In a smoke run this produced
`var_recall ≈ 0.33`, `n_shared = 3`, and a degenerate tree (`ad_f1 = 0`) —
not a harness bug, but a real property of the graph/detection step.

For a scientifically useful benchmark, **tune so each truth variant sits well
above the GFA detection threshold and several detected variants share ancestry**:

- Raise `--total-depth` (e.g. 600–1000) so low-frequency clones accumulate
  enough reads to form GFA edges.
- Reduce `--n-mutations` (e.g. 6–8) so clone frequencies are higher on average.
- Check `var_recall` and `n_shared` first; if `n_shared < 3`, ad/pc metrics
  are low-power regardless of fp/fn.

## Run one evaluation

```bash
# ONT (works everywhere)
./run_eval.sh --outdir /tmp/eval_ont  --profile ont-r10 \
  --n-mutations 12 --total-depth 600 --seed 1

# HiFi (Linux/x86 only)
./run_eval.sh --outdir /tmp/eval_hifi --profile hifi \
  --n-mutations 12 --total-depth 600 --seed 1
```

Optional flags: `--fp F` (SCITE fp rate, default 0.001) and `--fn F` (SCITE fn
rate, default 0.05).

Outputs land under `<outdir>/`: `truth/`, `reads/`, `himito/`, and
`metrics.tsv`.

### Call-tuning defaults in `run_himito.sh`

The wrapper passes sim-appropriate `Himito call` flags by default:

| Flag | Default | Why |
|------|---------|-----|
| `--minimal-ac 0` | 0 | Matrix construction keeps only reads whose alt-count exceeds `minimal-ac`. Simulated reads carry few alts each; the clinical default (1) empties the matrix so `Himito lineage` has nothing to work with. |
| `--vaf 0.005` | 0.005 | Captures low-frequency clones (~0.5%). The clinical default (0.01) drops sub-1% clones. |
| `--pval 1.0` | 1.0 | Disables the permutation test. Clean simulated variants have balanced strand coverage; the test rejects them at clinical thresholds. |

These can be overridden when calling `run_himito.sh` directly:
```bash
./run_himito.sh --outdir /tmp/eval_ont --profile ont-r10 \
  --minimal-ac 0 --vaf 0.005 --pval 1.0
```

## Sweep fp/fn

After `run_eval.sh` completes (matrix and VCF are fixed), sweep SCITE error
rates by re-running only `Himito lineage` per grid cell:

```bash
./sweep_fpfn.sh --outdir /tmp/eval_ont --profile ont-r10 \
  --fp-grid "0.0005 0.001 0.005 0.01" --fn-grid "0.02 0.05 0.1 0.2"
```

Writes `sweep_metrics.tsv` and prints the best cell by `ad_f1` (tie-break:
`var_f1`).

## Metrics (in `metrics.tsv` / `sweep_metrics.tsv`)

| Column | Meaning |
|--------|---------|
| `var_precision` / `recall` / `f1` | Variant *detection*: PASS SNVs in `sim.vcf` vs the truth SNV set. Independent of the tree. |
| `ad_precision` / `recall` / `f1` | **Ancestor–descendant** accuracy: over shared variants, are truth ancestral pairs preserved in the reconstruction? The primary tree-topology metric. |
| `pc_recall` | **Parent–child** edge recall: stricter — fraction of exact truth parent→child mutation edges recovered. |
| `n_shared` | Variants present in both the truth tree and the reconstructed tree (denominator for tree metrics). |

Precision = "of what we called, how much is right" (specificity of calls);
recall/sensitivity = "of the truth, how much we recovered".

Column order in the TSV: `profile`, `fp`, `fn`, `n_truth_vars`,
`n_detected_vars`, `n_shared`, `var_precision`, `var_recall`, `var_f1`,
`ad_precision`, `ad_recall`, `ad_f1`, `pc_recall`.

## Suggested experiments

- **Depth titration:** vary `--total-depth` (e.g. 300/600/1000) to find where
  heteroplasmy resolution breaks down.
- **Tree size:** vary `--n-mutations` (6/12/24) to test scaling.
- **Replicates:** vary `--seed` and average metrics per profile.
- **fp/fn refinement:** run `sweep_fpfn.sh` per profile; the recommended
  `--fp`/`--fn` for HiFi vs ONT-R10 are the best-cell values, averaged across
  seeds.

## Notes / limitations

- **HiFi is Linux/x86 only.** `ccs` (pbccs) has no macOS-ARM64 binary; use
  `--profile ont-r10` for local development on Apple Silicon.
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
