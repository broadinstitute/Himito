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
        File reference_fa

        String profile = "ont-r10"          # hifi | ont-r10 | ont-denoised (hifi needs ccs -> linux/amd64)
        Float pval = 1                       # Himito call -p forwarded to run_eval.sh

        Float? fp
        Float? fn

        Float min_hf = 0.01
        Float max_hf = 0.95
        Float sim_min_hf = 0.01
        Float sim_max_hf = 0.95
        Float internal_keep = 0.20
        

        String docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/himito-lineage-sim:dev"
        RuntimeAttr? runtime_attr_override
    }

    String fp_arg = if defined(fp) then "--fp " + select_first([fp]) else ""
    String fn_arg = if defined(fn) then "--fn " + select_first([fn]) else ""
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, object {}])

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
                        fp_arg = fp_arg,
                        fn_arg = fn_arg,
                        min_hf = min_hf,
                        max_hf = max_hf,
                        sim_min_hf = sim_min_hf,
                        sim_max_hf = sim_max_hf,
                        internal_keep = internal_keep,
                        reference_fa = reference_fa,
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

task RunEvalCell {

    meta {
        description: "One simulate -> reconstruct -> score cycle (run_eval.sh) for a single (seed, n_mutations, total_depth) grid cell."
    }

    parameter_meta {
        seed:          "RNG seed for simulate_tree.py / simulate_reads.sh"
        n_mutations:   "number of heteroplasmic SNVs in the truth tree"
        total_depth:   "total simulated read depth across all clones"
        profile:       "read profile: hifi, ont-r10, or ont-denoised"
        fp_arg:        "optional '--fp <rate>' override; empty uses the profile's Himito default"
        fn_arg:        "optional '--fn <rate>' override; empty uses the profile's Himito default"
        min_hf:        "shared HF floor for Himito call -v and lineage --min-hf"
        reference_fa:  "mitochondrial reference FASTA (required)"
    }

    input {
        File reference_fa
        Int seed
        Int n_mutations
        Int total_depth

        String profile
        Float pval
        String fp_arg
        String fn_arg

        Float min_hf
        Float max_hf
        Float sim_min_hf
        Float sim_max_hf
        Float internal_keep

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
            --pval ~{pval} \
            --ref "$REF" \
            --min-hf ~{min_hf} \
            --max-hf ~{max_hf} \
            --sim-min-hf ~{sim_min_hf} \
            --sim-max-hf ~{sim_max_hf} \
            --internal-keep ~{internal_keep} \
            ~{fp_arg} \
            ~{fn_arg}

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
