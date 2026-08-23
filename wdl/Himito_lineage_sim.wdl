version 1.0

# WDL port of run_simulation.sh: a Cartesian sweep of the simulate -> reconstruct
# -> score cycle (run_eval.sh) over seeds x n-mutations x depths. Each grid cell
# runs one full RunEvalCell task; CombineMetrics concatenates the per-cell metrics
# (with identifying seed/n_mutations/total_depth columns) into a single metrics.tsv.

workflow HimitoLineageSim {
    input {
        Array[Int] seeds
        Array[Int] n_mutations
        Array[Int] depths

        String profile = "ont-r10"          # hifi | ont-r10 (hifi needs ccs -> linux/amd64)
        Float pval = 1                       # Himito call -p forwarded to run_eval.sh
        Float fp = 0.001                     # SCITE false-positive rate
        Float fn = 0.05                      # SCITE false-negative rate

        # Tuning knobs. Defaults match run_eval.sh's own defaults, so leaving them
        # unset reproduces a plain `run_eval.sh` invocation. They carry concrete
        # defaults (rather than being bare optionals) because Cromwell fails to
        # look up a no-default optional threaded through nested scatter into a task.
        Float min_hf = 0.01
        Float max_hf = 0.99
        Float sim_min_hf = 0.05
        Float sim_max_hf = 0.99
        Float internal_keep = 0.20

        # Optional reference override; default = rCRS.fasta baked into the image.
        File? reference_fa

        String docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/himito_eval:dev"
        RuntimeAttr? runtime_attr_override
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
                        pval = pval,
                        fp = fp,
                        fn = fn,
                        min_hf = min_hf,
                        max_hf = max_hf,
                        sim_min_hf = sim_min_hf,
                        sim_max_hf = sim_max_hf,
                        internal_keep = internal_keep,
                        reference_fa = reference_fa,
                        docker = docker,
                        runtime_attr_override = runtime_attr_override
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

task RunEvalCell {

    meta {
        description: "One simulate -> reconstruct -> score cycle (run_eval.sh) for a single (seed, n_mutations, total_depth) grid cell."
    }

    parameter_meta {
        seed:          "RNG seed for simulate_tree.py / simulate_reads.sh"
        n_mutations:   "number of heteroplasmic SNVs in the truth tree"
        total_depth:   "total simulated read depth across all clones"
        profile:       "read profile: hifi or ont-r10"
        reference_fa:  "optional reference FASTA; default is rCRS.fasta baked into the image"
    }

    input {
        Int seed
        Int n_mutations
        Int total_depth

        String profile
        Float pval
        Float fp
        Float fn

        Float min_hf
        Float max_hf
        Float sim_min_hf
        Float sim_max_hf
        Float internal_keep

        File? reference_fa

        String docker
        RuntimeAttr? runtime_attr_override
    }

    String cell = "seed~{seed}_mut~{n_mutations}_depth~{total_depth}"
    String ref_arg = if defined(reference_fa) then "~{select_first([reference_fa])}" else "/opt/lineage_sim/rCRS.fasta"

    command <<<
        set -euxo pipefail

        # Tools + scripts + models + Himito binary + rCRS.fasta are baked into the
        # image; export the paths the shell scripts look for so they run standalone.
        export HIMITO="/Himito/target/release/Himito"
        export PBSIM_MODEL_DIR="/opt/lineage_sim/pbsim3_models"
        export REF="~{ref_arg}"

        /opt/lineage_sim/run_eval.sh \
            --outdir "~{cell}" \
            --profile "~{profile}" \
            --n-mutations ~{n_mutations} \
            --total-depth ~{total_depth} \
            --seed ~{seed} \
            --pval ~{pval} \
            --fp ~{fp} \
            --fn ~{fn} \
            --ref "$REF" \
            --min-hf ~{min_hf} \
            --max-hf ~{max_hf} \
            --sim-min-hf ~{sim_min_hf} \
            --sim-max-hf ~{sim_max_hf} \
            --internal-keep ~{internal_keep}

        # Prepend the identifying columns run_simulation.sh adds to the combined table.
        metrics="~{cell}/metrics.tsv"
        header=$(head -n 1 "$metrics")
        row=$(tail -n 1 "$metrics")
        printf 'seed\tn_mutations\ttotal_depth\t%s\n' "$header" > cell_metrics.tsv
        printf '%s\t%s\t%s\t%s\n' "~{seed}" "~{n_mutations}" "~{total_depth}" "$row" >> cell_metrics.tsv

        # Bundle the full cell (truth/, reads/, himito/) for debugging/inspection.
        tar czf "~{cell}.tar.gz" "~{cell}"
    >>>

    output {
        File cell_metrics = "cell_metrics.tsv"
        File cell_bundle = "~{cell}.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
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
