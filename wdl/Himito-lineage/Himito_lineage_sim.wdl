version 1.0

# WDL port of run_simulation.sh: a Cartesian sweep of the simulate -> reconstruct
# -> score cycle (run_eval.sh) over seeds x n-mutations x depths. Each grid cell
# runs one full RunEvalCell task; CombineMetrics concatenates the per-cell metrics
# (with identifying seed/n_mutations/total_depth columns) into a single metrics.tsv.
#
# ---------------------------------------------------------------------------
# THE THREE FREQUENCY GATES
# ---------------------------------------------------------------------------
# The pipeline applies three independent frequency floors, and a truth variant is
# destroyed by whichever is highest:
#
#   denoise --vaf   (before the graph; follows --vaf when unset)
#   call -v         (= min_hf)
#   call -f         (= frequency_threshold; Himito's own default is 0.2)
#
# EVERY one of them must sit below the LOWEST TRUTH HF, or truth variants are
# deleted before they can be called and var_recall is capped no matter how much
# depth is used -- these are frequency floors, not coverage floors.
#
# Note that is the lowest *truth* HF, not sim_min_hf. sim_min_hf is only the floor
# the simulator may descend to; the achieved minimum is usually well above it (the
# standard chain config floors at 0.296 under sim_min_hf=0.05), so comparing gates
# against sim_min_hf rejects perfectly good configs. ValidateConfig computes the
# exact minimum for topology=chain and falls back to sim_min_hf as a lower bound
# otherwise, then fails the workflow before a single expensive cell runs.
#
# Defaults changed from the first version of this WDL, all of which produced an
# EMPTY VCF and a failed `Himito lineage` on every cell:
#   profile              ont-r10 -> ont-denoised  (ont-r10 skips denoise entirely)
#   frequency_threshold  unset (0.2) -> 0.05      (0.2 is above every truth variant)
#   sim_min_hf           0.01 -> 0.05             (0.01 sits under the gates)
#   sim_max_hf           0.95 -> 0.99             (matches run_eval.sh)
#
# NOTE: topology/denoise_vaf require an image built from the current
# lineage_sim/ scripts. Rebuild the docker image before using them.

workflow HimitoLineageSim {
    input {
        Array[Int] seeds
        Array[Int] n_mutations
        Array[Int] depths
        File reference_fa

        # ont-denoised | ont-r10 | hifi. ont-r10 SKIPS the denoise step, which
        # leaves indel artifacts crowding the graph and yields an empty VCF;
        # hifi needs ccs -> linux/amd64.
        String profile = "ont-denoised"

        # Himito call -f. Himito's own default is 0.2, which is above essentially
        # every simulated truth variant and empties the VCF. Must be < sim_min_hf.
        Float frequency_threshold = 0.05

        # Himito denoise --vaf. Unset follows --vaf (= min_hf). Previously
        # hardcoded to 0.03 inside run_himito.sh, which silently deleted truth
        # variants in any config whose clones were rarer than that. Must be
        # < sim_min_hf when set.
        Float? denoise_vaf

        Float? fp
        Float? fn

        Float min_hf = 0.01
        Float max_hf = 0.95

        # Band for the SIMULATED truth frequencies -- a different thing from
        # min_hf/max_hf, which filter what was called.
        Float sim_min_hf = 0.05
        Float sim_max_hf = 0.99

        # Tree shape. "chain" is the knob that raises mutation-ORDERING difficulty
        # without costing detection (a unary path never splits clone mass between
        # siblings, so its minimum HF is higher than a random tree of the same
        # size). "star" has no ancestral pairs and is a control. Raising
        # n_mutations is NOT a difficulty knob: it only shrinks clones and costs
        # variant recall.
        String topology = "random"

        # Terminal mass each internal node keeps. THIS is what makes an individual
        # ordering call hard; lower is harder, and it costs no detection.
        Float internal_keep = 0.20

        String docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/himito-lineage-sim:dev"
        RuntimeAttr? runtime_attr_override
    }

    String fp_arg = if defined(fp) then "--fp " + select_first([fp]) else ""
    String fn_arg = if defined(fn) then "--fn " + select_first([fn]) else ""
    String denoise_vaf_arg = if defined(denoise_vaf) then "--denoise-vaf " + select_first([denoise_vaf]) else ""
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, object {}])

    # Fail the whole workflow on a misconfigured gate before any cell burns a VM.
    call ValidateConfig {
        input:
            profile = profile,
            frequency_threshold = frequency_threshold,
            denoise_vaf = denoise_vaf,
            min_hf = min_hf,
            sim_min_hf = sim_min_hf,
            topology = topology,
            internal_keep = internal_keep,
            n_mutations = n_mutations,
            docker = docker
    }

    scatter (seed in seeds) {
        scatter (nmut in n_mutations) {
            scatter (depth in depths) {
                call RunEvalCell {
                    input:
                        seed = seed,
                        n_mutations = nmut,
                        total_depth = depth,
                        profile = profile,
                        fp_arg = fp_arg,
                        fn_arg = fn_arg,
                        denoise_vaf_arg = denoise_vaf_arg,
                        frequency_threshold = frequency_threshold,
                        min_hf = min_hf,
                        max_hf = max_hf,
                        sim_min_hf = sim_min_hf,
                        sim_max_hf = sim_max_hf,
                        topology = topology,
                        internal_keep = internal_keep,
                        reference_fa = reference_fa,
                        validated = ValidateConfig.ok,
                        docker = docker,
                        runtime_attr_override = runtime_attr
                }
            }
        }
    }

    Array[File] cell_metrics = flatten(flatten(RunEvalCell.cell_metrics))
    Array[File] cell_bundles = flatten(flatten(RunEvalCell.cell_bundle))

    call CombineMetrics {
        input:
            cell_metrics = cell_metrics,
            docker = docker
    }

    output {
        File metrics = CombineMetrics.combined
        Array[File] per_cell_metrics = cell_metrics
        Array[File] per_cell_bundles = cell_bundles
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task ValidateConfig {

    meta {
        description: "Check the frequency-gate invariant and topology feasibility before the grid runs. Cheap, and it turns a whole sweep of empty VCFs into one clear error."
    }

    input {
        String profile
        Float frequency_threshold
        Float? denoise_vaf
        Float min_hf
        Float sim_min_hf
        String topology
        Float internal_keep
        Array[Int] n_mutations
        String docker
    }

    # Unset denoise_vaf follows --vaf, which run_eval.sh sets from min_hf.
    Float effective_denoise_vaf = select_first([denoise_vaf, min_hf])

    command <<<
        set -euo pipefail
        python3 <<'PY'
import sys

profile   = "~{profile}"
freq      = float("~{frequency_threshold}")
denoise   = float("~{effective_denoise_vaf}")
min_hf    = float("~{min_hf}")
sim_min   = float("~{sim_min_hf}")
topology  = "~{topology}"
keep      = float("~{internal_keep}")
nmuts     = [int(x) for x in "~{sep=' ' n_mutations}".split()]
REF_FRACTION = 0.15   # simulate_tree.py default; not exposed by run_eval.sh

errors, warnings = [], []

if topology not in ("random", "chain", "star"):
    errors.append(f"topology must be random|chain|star, got {topology!r}")

# --- the gate invariant -------------------------------------------------------
# Compare gates against the LOWEST TRUTH HF, not against sim_min_hf. sim_min_hf is
# only the floor the simulator is permitted to descend to; the frequency actually
# achieved is usually well above it (the standard chain config floors at 0.296
# under sim_min_hf=0.05), so checking against sim_min_hf rejects working configs.
#
# For a chain the minimum is exact and deterministic. For random/star it is
# seed-dependent and unknowable here, so sim_min_hf is the only guaranteed bound
# and the check is necessarily weaker.
GATES = (("call -f (frequency_threshold)", freq),
         ("denoise --vaf", denoise),
         ("call -v (min_hf)", min_hf))

def check_gates(min_truth_hf, exact, ctx):
    for label, value in GATES:
        if value >= min_truth_hf:
            errors.append(
                f"{label} = {value} is at or above the {'' if exact else 'lowest possible '}"
                f"truth HF {min_truth_hf:.4f}{ctx}. This gate deletes truth variants "
                "before they can be called, and depth cannot compensate."
            )
        elif value > min_truth_hf / 2.0:
            warnings.append(
                f"{label} = {value} has under 2x margin on the {'' if exact else 'lowest possible '}"
                f"truth HF {min_truth_hf:.4f}{ctx}. Measured: denoise --vaf 0.03 against a "
                "truth HF of 0.033 cut var_recall to 0.533, and 0.01 restored it to 0.867 "
                "on the same reads. Check var_recall in the output."
            )

if topology == "chain":
    # A chain never splits mass between siblings, so the deepest clone is exactly
    # (1 - ref_fraction) * (1 - internal_keep)^n -- deterministic, so an
    # infeasible chain cannot be rescued by reseeding.
    for n in sorted(set(nmuts)):
        deepest = (1.0 - REF_FRACTION) * (1.0 - keep) ** n
        if deepest <= sim_min:
            errors.append(
                f"topology=chain with n_mutations={n} puts the deepest clone at "
                f"{deepest:.4f}, at or below sim_min_hf={sim_min}. Chain "
                "frequencies are deterministic, so no seed fixes this: lower "
                "sim_min_hf or internal_keep, or reduce n_mutations."
            )
        else:
            check_gates(deepest, True, f" for n_mutations={n}")
else:
    check_gates(sim_min, False, "")

# --- non-fatal ----------------------------------------------------------------
if profile == "ont-r10":
    warnings.append(
        "profile=ont-r10 SKIPS Himito denoise, so indel artifacts crowd the graph "
        "and the VCF usually comes back empty. Use ont-denoised unless you are "
        "deliberately measuring the un-denoised path."
    )
if topology == "star":
    warnings.append(
        "topology=star produces no ancestral pairs, so ad_recall has an empty "
        "denominator and ad_* is uninformative. It is a control, not a config."
    )

for w in warnings:
    print(f"WARNING: {w}")
if errors:
    print()
    for e in errors:
        print(f"ERROR: {e}")
    sys.exit(1)

print(f"config OK: profile={profile} topology={topology} internal_keep={keep}")
if topology == "chain":
    basis = ", ".join(
        f"n={n}:{(1.0 - REF_FRACTION) * (1.0 - keep) ** n:.4f}" for n in sorted(set(nmuts)))
    print(f"  gates checked against the exact chain minimum ({basis})")
else:
    print(f"  gates checked against sim_min_hf={sim_min} (lowest possible truth HF; "
          "the achieved minimum is seed-dependent and usually higher)")
print(f"  call -f={freq}  denoise --vaf={denoise}  call -v={min_hf}")
PY
        echo ok > ok.txt
    >>>

    output {
        String ok = read_string("ok.txt")
    }

    runtime {
        cpu:    1
        memory: "1 GiB"
        disks:  "local-disk 10 HDD"
        docker: docker
    }
}

task RunEvalCell {

    meta {
        description: "One simulate -> reconstruct -> score cycle (run_eval.sh) for a single (seed, n_mutations, total_depth) grid cell."
    }

    parameter_meta {
        seed:          "RNG seed for simulate_tree.py / simulate_reads.sh"
        n_mutations:   "number of heteroplasmic SNVs in the truth tree. NOT a difficulty knob: clone mass is conserved, so raising it only costs variant recall"
        total_depth:   "total simulated read depth across all clones"
        profile:       "read profile: ont-denoised (recommended), ont-r10 (skips denoise; usually empty VCF), or hifi"
        fp_arg:        "optional '--fp <rate>' override; empty uses the profile's Himito default"
        fn_arg:        "optional '--fn <rate>' override; empty uses the profile's Himito default"
        denoise_vaf_arg: "optional '--denoise-vaf <F>' override; empty follows --vaf. Must stay below sim_min_hf"
        frequency_threshold: "Himito call -f. Must stay below sim_min_hf; Himito's own default (0.2) empties the VCF"
        min_hf:        "shared HF floor for Himito call -v and lineage --min-hf"
        sim_min_hf:    "floor for the SIMULATED truth frequencies -- distinct from min_hf, which filters what was called"
        topology:      "truth tree shape: random | chain | star. chain raises ordering difficulty without costing detection"
        internal_keep: "terminal mass kept by each internal node; lower is harder to order and costs no detection"
        reference_fa:  "mitochondrial reference FASTA (required)"
        validated:     "sentinel from ValidateConfig; forces gate validation to run before any cell"
    }

    input {
        File reference_fa
        Int seed
        Int n_mutations
        Int total_depth

        String profile
        String fp_arg
        String fn_arg
        String denoise_vaf_arg
        Float frequency_threshold

        Float min_hf
        Float max_hf
        Float sim_min_hf
        Float sim_max_hf
        String topology
        Float internal_keep

        String validated

        String docker
        RuntimeAttr? runtime_attr_override
    }

    String cell = "seed~{seed}_mut~{n_mutations}_depth~{total_depth}"
    command <<<
        set -euxo pipefail

        # Tools + scripts + models + Himito binary are baked into the image;
        # export the paths the shell scripts look for so they run standalone.
        export HIMITO="/Himito/target/release/Himito"
        export PBSIM_MODEL_DIR="/opt/lineage_sim/pbsim3_models"
        export REF="~{reference_fa}"

        /opt/lineage_sim/run_eval.sh \
            --outdir "~{cell}" \
            --profile "~{profile}" \
            --n-mutations ~{n_mutations} \
            --total-depth ~{total_depth} \
            --seed ~{seed} \
            --ref "$REF" \
            --min-hf ~{min_hf} \
            --max-hf ~{max_hf} \
            --sim-min-hf ~{sim_min_hf} \
            --sim-max-hf ~{sim_max_hf} \
            --topology ~{topology} \
            --internal-keep ~{internal_keep} \
            --frequency-threshold ~{frequency_threshold} \
            ~{fp_arg} \
            ~{fn_arg} \
            ~{denoise_vaf_arg}

        # Prepend the identifying columns run_simulation.sh adds to the combined
        # table. Read the header from the file rather than assuming a column
        # order: score_lineage.py's FIELDS has changed (n_truth_pairs added) and
        # will change again.
        metrics="~{cell}/metrics.tsv"
        header=$(head -n 1 "$metrics")
        row=$(tail -n 1 "$metrics")
        printf 'seed\tn_mutations\ttotal_depth\ttopology\tinternal_keep\t%s\n' "$header" > cell_metrics.tsv
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "~{seed}" "~{n_mutations}" "~{total_depth}" "~{topology}" "~{internal_keep}" "$row" >> cell_metrics.tsv

        # Bundle the full cell (truth/, reads/, himito/) for debugging/inspection.
        tar czf "~{cell}.tar.gz" "~{cell}"
    >>>

    output {
        File cell_metrics = "cell_metrics.tsv"
        File cell_bundle = "~{cell}.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CombineMetrics {

    meta {
        description: "Concatenate per-cell metrics.tsv files into one table, keeping a single header."
    }

    input {
        Array[File] cell_metrics
        String docker
    }

    command <<<
        set -euxo pipefail
        first=1
        : > metrics.tsv
        for f in ~{sep=' ' cell_metrics}; do
            [[ -s "$f" ]] || continue
            if [[ $first -eq 1 ]]; then
                cat "$f" >> metrics.tsv
                first=0
            else
                tail -n +2 "$f" >> metrics.tsv
            fi
        done
        column -t metrics.tsv || true
    >>>

    output {
        File combined = "metrics.tsv"
    }

    runtime {
        cpu:    1
        memory: "2 GiB"
        disks:  "local-disk 10 HDD"
        docker: docker
    }
}
