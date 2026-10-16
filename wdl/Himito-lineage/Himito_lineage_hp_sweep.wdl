version 1.0

import "https://raw.githubusercontent.com/broadinstitute/Himito/hs_lineage_posterio/wdl/Himito-lineage/Himito_lineage_sim.wdl" as sim

workflow HimitoLineageHpSweep {
    input {
        # Each seed runs one lineage job per (hp_fp, hp_fn) grid cell: 30 seeds x a 5x6 grid is 900 runs.
        Array[Int] seeds
        Int n_mutations
        Int depth
        File reference_fa

        String profile = "ont-denoised"
        Float frequency_threshold = 0.05
        Float? denoise_vaf
        Float? fp
        Float? fn
        Float min_hf = 0.01
        Float max_hf = 0.95
        Float sim_min_hf = 0.05
        Float sim_max_hf = 0.99
        String topology = "random"
        Float internal_keep = 0.20

        String background = "none"
        Int hp_lineage_indels = 0
        String hp_unstable = ""
        Float hp_min_subclone_freq = 0.002
        Boolean collapse_hp = false

        # Multiples of the ont-denoised global rates (0.0001, 0.01), capped at 0.5.
        # Override when --profile is not ont-denoised.
        Array[Float] hp_fp_grid
        Array[Float] hp_fn_grid
        String objective = "hp_intrusion"
        String secondary = "ad_f1"
        Int mcmc_iterations = 2000
        Int mcmc_chains = 2

        String docker = "us.gcr.io/broad-dsp-lrma/hangsuunc/himito-lineage-sim:dev"
    }

    String fp_arg = if defined(fp) then "--fp " + select_first([fp]) else ""
    String fn_arg = if defined(fn) then "--fn " + select_first([fn]) else ""
    String denoise_vaf_arg = if defined(denoise_vaf) then "--denoise-vaf " + select_first([denoise_vaf]) else ""
    String background_arg = if background == "none" then "" else "--background \"" + background + "\""
    String hp_lineage_arg = if hp_lineage_indels > 0 then "--hp-lineage-indels " + hp_lineage_indels else ""
    String hp_unstable_arg = if hp_unstable != "" then "--hp-unstable \"" + hp_unstable + "\"" else ""
    String hp_min_subclone_arg = if hp_unstable != "" then "--hp-min-subclone-freq " + hp_min_subclone_freq else ""

    # One n_mutations value: ValidateConfig still wants the array the sim workflow uses.
    Array[Int] n_mutations_for_config = [n_mutations]

    call sim.ValidateConfig {
        input:
            profile = profile,
            frequency_threshold = frequency_threshold,
            denoise_vaf = denoise_vaf,
            min_hf = min_hf,
            sim_min_hf = sim_min_hf,
            topology = topology,
            internal_keep = internal_keep,
            n_mutations = n_mutations_for_config,
            background = background,
            hp_unstable = hp_unstable,
            hp_lineage_indels = hp_lineage_indels,
            collapse_hp = collapse_hp,
            hp_fp_set = true,
            hp_fn_set = true,
            require_hp_mode = true,
            docker = docker
    }

    # WDL 1.0 has no zip over seeds x rates. Index SimulateCell outputs with
    # range(length(seeds)) from the nested rate scatter.
    scatter (seed_i in range(length(seeds))) {
        call SimulateCell {
            input:
                seed = seeds[seed_i],
                n_mutations = n_mutations,
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
                hp_lineage_indels = hp_lineage_indels,
                hp_unstable = hp_unstable,
                background_arg = background_arg,
                hp_lineage_arg = hp_lineage_arg,
                hp_unstable_arg = hp_unstable_arg,
                hp_min_subclone_arg = hp_min_subclone_arg,
                reference_fa = reference_fa,
                validated = ValidateConfig.ok,
                docker = docker
        }
    }

    scatter (seed_i in range(length(seeds))) {
        scatter (hp_fp_rate in hp_fp_grid) {
            scatter (hp_fn_rate in hp_fn_grid) {
                call LineageScoreCell {
                    input:
                        seed = seeds[seed_i],
                        matrix = SimulateCell.matrix[seed_i],
                        vcf = SimulateCell.vcf[seed_i],
                        truth_tar = SimulateCell.truth_tar[seed_i],
                        hp_fp_rate = hp_fp_rate,
                        hp_fn_rate = hp_fn_rate,
                        collapse_hp = collapse_hp,
                        profile = profile,
                        background = background,
                        min_hf = min_hf,
                        max_hf = max_hf,
                        fp_arg = fp_arg,
                        fn_arg = fn_arg,
                        mcmc_iterations = mcmc_iterations,
                        mcmc_chains = mcmc_chains,
                        docker = docker
                }
            }
        }
    }

    # Three nested scatters; two flattens, same as HimitoLineageSim.
    Array[File] cell_metrics = flatten(flatten(LineageScoreCell.cell_metrics))

    call sim.CombineMetrics {
        input:
            cell_metrics = cell_metrics,
            docker = docker
    }

    call PickBestHpRates {
        input:
            metrics = CombineMetrics.combined,
            objective = objective,
            secondary = secondary,
            docker = docker
    }

    output {
        File metrics = CombineMetrics.combined
        File summary = PickBestHpRates.summary
        File best_hp_rates = PickBestHpRates.best_hp_rates
        File report = PickBestHpRates.report
        Array[File] per_cell_metrics = cell_metrics
    }
}

task SimulateCell {

    meta {
        description: "run_eval.sh through Himito call (--stop-after call). One seed. Lineage is scored later, per HP rate."
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
        Int hp_lineage_indels
        String hp_unstable
        String background_arg
        String hp_lineage_arg
        String hp_unstable_arg
        String hp_min_subclone_arg
        String validated
        String docker
    }

    String cell = "seed~{seed}_mut~{n_mutations}_depth~{total_depth}" + (if hp_lineage_indels > 0 || hp_unstable != "" then "_hp" else "")

    command <<<
        set -euxo pipefail
        echo "config ~{validated}" >&2

        export HIMITO="/Himito/target/release/Himito"
        export PBSIM_MODEL_DIR="/opt/lineage_sim/pbsim3_models"
        export REF="~{reference_fa}"

        # Stops after sim.matrix.csv and sim.vcf. truth/ stays.
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
            --stop-after call \
            ~{fp_arg} \
            ~{fn_arg} \
            ~{denoise_vaf_arg} \
            ~{background_arg} \
            ~{hp_lineage_arg} \
            ~{hp_unstable_arg} \
            ~{hp_min_subclone_arg}

        if [[ -f "~{cell}/truth/background.txt" ]]; then
            cp "~{cell}/truth/background.txt" background.txt
        else
            : > background.txt
        fi
        tar czf truth.tar.gz -C "~{cell}" truth
    >>>

    output {
        File truth_tar = "truth.tar.gz"
        File matrix = "~{cell}/himito/sim.matrix.csv"
        File vcf = "~{cell}/himito/sim.vcf"
        File background_txt = "background.txt"
    }

    runtime {
        cpu:         1
        memory:      "4 GiB"
        disks:       "local-disk 50 SSD"
        preemptible: 2
        docker:      docker
    }
}

task LineageScoreCell {

    meta {
        description: "Himito lineage at one (hp_fp_rate, hp_fn_rate), then score_lineage.py. The global fp/fn stay at the profile default unless fp_arg/fn_arg are set."
    }

    input {
        Int seed
        File matrix
        File vcf
        File truth_tar
        Float hp_fp_rate
        Float hp_fn_rate
        Boolean collapse_hp
        String profile
        String background
        Float min_hf
        Float max_hf
        String fp_arg
        String fn_arg
        Int mcmc_iterations
        Int mcmc_chains
        String docker
    }

    command <<<
        set -euxo pipefail
        export HIMITO="/Himito/target/release/Himito"
        mkdir -p unpacked
        tar xzf "~{truth_tar}" -C unpacked

        PROFILE="~{profile}"
        case "$PROFILE" in
          hifi)          DTYPE=pacbio;       DEFAULT_FP=0.005;  DEFAULT_FN=0.05;;
          ont-denoised)  DTYPE=ont-denoised; DEFAULT_FP=0.005; DEFAULT_FN=0.02;;
          ont-r10)       DTYPE=ont-r10;      DEFAULT_FP=0.005; DEFAULT_FN=0.02;;
          *) echo "profile must be hifi or ont-r10 or ont-denoised" >&2; exit 1;;
        esac

        FP_ARG="~{fp_arg}"
        FN_ARG="~{fn_arg}"
        FP_RATE="$DEFAULT_FP"
        FN_RATE="$DEFAULT_FN"
        [[ -n "$FP_ARG" ]] && FP_RATE="${FP_ARG#--fp }"
        [[ -n "$FN_ARG" ]] && FN_RATE="${FN_ARG#--fn }"

        LINEAGE_ARGS=(
          -m "~{matrix}" -v "~{vcf}"
          --hp-fp-rate ~{hp_fp_rate}
          --hp-fn-rate ~{hp_fn_rate}
          --min-hf ~{min_hf} --max-hf ~{max_hf}
          --mcmc-iterations ~{mcmc_iterations}
          --mcmc-chains ~{mcmc_chains}
          -d "$DTYPE"
          -o sim_lineage
        )
        [[ -n "$FP_ARG" ]] && LINEAGE_ARGS+=(--fp-rate "$FP_RATE")
        [[ -n "$FN_ARG" ]] && LINEAGE_ARGS+=(--fn-rate "$FN_RATE")
        if [[ "~{collapse_hp}" == "true" ]]; then
            LINEAGE_ARGS+=(--collapse-hp)
        fi
        "$HIMITO" lineage "${LINEAGE_ARGS[@]}"

        SCORE_ARGS=()
        [[ -f unpacked/truth/background.txt ]] && SCORE_ARGS+=(--truth-background unpacked/truth/background.txt)
        [[ -f unpacked/truth/hp_tracts.tsv ]] && SCORE_ARGS+=(--truth-hp-tracts unpacked/truth/hp_tracts.tsv)
        [[ -f unpacked/truth/hp_lineage.tsv ]] && SCORE_ARGS+=(--truth-hp-lineage unpacked/truth/hp_lineage.tsv)
        if [[ "~{background}" != "none" ]]; then
            SCORE_ARGS+=(--background "~{background}")
        fi
        HAP_ARGS=()
        if [[ -s sim_lineage.cleaned_matrix_all_variants.csv && -s unpacked/truth/clones.tsv ]]; then
            HAP_ARGS=(--recon-matrix sim_lineage.cleaned_matrix_all_variants.csv
                      --truth-clones unpacked/truth/clones.tsv)
        fi
        python3 /opt/lineage_sim/score_lineage.py \
            --truth-tree unpacked/truth/truth_mutation_tree.tsv \
            --recon-tree sim_lineage.mutation_tree.tsv \
            --truth-variants unpacked/truth/truth_variants.txt \
            --vcf "~{vcf}" \
            --profile "$PROFILE" --fp "$FP_RATE" --fn "$FN_RATE" \
            ${HAP_ARGS[@]+"${HAP_ARGS[@]}"} \
            ${SCORE_ARGS[@]+"${SCORE_ARGS[@]}"} \
            --metrics-out raw_metrics.tsv

        header=$(head -n 1 raw_metrics.tsv)
        row=$(tail -n 1 raw_metrics.tsv)
        printf 'seed\thp_fp_rate\thp_fn_rate\tcollapse_hp\t%s\n' "$header" > cell_metrics.tsv
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "~{seed}" "~{hp_fp_rate}" "~{hp_fn_rate}" "~{collapse_hp}" "$row" >> cell_metrics.tsv
    >>>

    output {
        File cell_metrics = "cell_metrics.tsv"
    }

    runtime {
        cpu:         1
        memory:      "2 GiB"
        disks:       "local-disk 10 HDD"
        preemptible: 3
        docker:      docker
    }
}

task PickBestHpRates {

    meta {
        description: "pick_best_hp_rates.py. A flat grid writes best=null and says to leave the default alone."
    }

    input {
        File metrics
        String objective
        String secondary
        String docker
    }

    command <<<
        set -euxo pipefail
        python3 /opt/lineage_sim/pick_best_hp_rates.py \
            --metrics "~{metrics}" \
            --summary hp_sweep_summary.tsv \
            --json best_hp_rates.json \
            --report hp_sweep_report.md \
            --objective "~{objective}" \
            --secondary "~{secondary}"
    >>>

    output {
        File summary = "hp_sweep_summary.tsv"
        File best_hp_rates = "best_hp_rates.json"
        File report = "hp_sweep_report.md"
    }

    runtime {
        cpu:    1
        memory: "2 GiB"
        disks:  "local-disk 10 HDD"
        docker: docker
    }
}
