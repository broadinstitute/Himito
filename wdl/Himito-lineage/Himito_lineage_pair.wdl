version 1.0

workflow Himito_lineage_eval {
    meta {
        description: "Compare PacBio and ONT lineages across matched target depths using production calling presets."
    }

    parameter_meta {
        desired_coverages: "Requested mean primary alignment depths on the mitochondrial contig. Both platforms are capped at the lower available depth."
        downsample_seed: "Seed for reproducible read-name subsampling."
        call_extra_args: "Optional quick-start overrides. Omit to use production caller presets and permutation testing."
    }

    input {
        File pacbio_bam
        File pacbio_bam_bai
        File ont_bam
        File ont_bam_bai
        File reference_fa
        String prefix
        String sample_id
        String himito_docker = "hangsuunc/himito-lineage-sim:dev"
        String? billing_project
        Int kmer_size = 21
        Float min_hf = 0.01
        Float max_hf = 0.95
        String chromo = "chrM"
        String pacbio_datatype = "pacbio"
        String ont_datatype = "ont-denoised"
        Array[Int]+ desired_coverages = [100, 300, 500, 700, 1000, 2000]
        Int downsample_seed = 7
        String? call_extra_args
        String? lineage_extra_args
        String? merge_extra_args
        Boolean snv_only = true
        Boolean pass_only = true
    }

    call CalculateCoverage as CalculateCoveragePacBio {
        input:
            bam = pacbio_bam,
            bai = pacbio_bam_bai,
            billing_project = billing_project,
            locus = chromo,
            prefix = prefix + ".pb"
    }

    call CalculateCoverage as CalculateCoverageONT {
        input:
            bam = ont_bam,
            bai = ont_bam_bai,
            billing_project = billing_project,
            locus = chromo,
            prefix = prefix + ".ont"
    }

    Float common_source_coverage = if CalculateCoveragePacBio.coverage < CalculateCoverageONT.coverage then CalculateCoveragePacBio.coverage else CalculateCoverageONT.coverage

    scatter (target in desired_coverages) {
        String depth_prefix = prefix + "." + target + "x"
        Float effective_target = if target < common_source_coverage then target * 1.0 else common_source_coverage

        call DownsampleBam as PacBioAtCoverage {
            input:
                bam = CalculateCoveragePacBio.subsetbam,
                bai = CalculateCoveragePacBio.subsetbai,
                locus = chromo,
                prefix = depth_prefix + ".pb",
                source_coverage = CalculateCoveragePacBio.coverage,
                target_coverage = effective_target,
                seed = downsample_seed,
                docker = himito_docker
        }

        call DownsampleBam as OntAtCoverage {
            input:
                bam = CalculateCoverageONT.subsetbam,
                bai = CalculateCoverageONT.subsetbai,
                locus = chromo,
                prefix = depth_prefix + ".ont",
                source_coverage = CalculateCoverageONT.coverage,
                target_coverage = effective_target,
                seed = downsample_seed,
                docker = himito_docker
        }

        call Himito_quickstart as PacbioCall {
            input:
                bam = PacBioAtCoverage.sampled_bam,
                bai = PacBioAtCoverage.sampled_bai,
                reference_fa = reference_fa,
                prefix = depth_prefix + ".pb",
                kmer_size = kmer_size,
                sample_id = sample_id,
                chromo = chromo,
                data_type = pacbio_datatype,
                extra_args = call_extra_args,
                himito_docker = himito_docker
        }

        call Himito_quickstart as OntCall {
            input:
                bam = OntAtCoverage.sampled_bam,
                bai = OntAtCoverage.sampled_bai,
                reference_fa = reference_fa,
                prefix = depth_prefix + ".ont",
                kmer_size = kmer_size,
                sample_id = sample_id,
                chromo = chromo,
                data_type = ont_datatype,
                extra_args = call_extra_args,
                himito_docker = himito_docker
        }

        call identify_shared_vcf as IdentifySharedVCF {
            input:
                pacbio_vcf = PacbioCall.vcf,
                ont_vcf = OntCall.vcf,
                prefix = depth_prefix,
                extra_args = merge_extra_args,
                snv_only = snv_only,
                pass_only = pass_only
        }

        call Himito_lineage as HimitoLineage {
            input:
                shared_vcf = IdentifySharedVCF.vcf,
                shared_vcf_tbi = IdentifySharedVCF.vcf_tbi,
                pacbio_matrix = PacbioCall.read_var_mat,
                ont_matrix = OntCall.read_var_mat,
                prefix = depth_prefix,
                min_hf = min_hf,
                max_hf = max_hf,
                extra_args = lineage_extra_args,
                himito_docker = himito_docker
        }

        call ScorePair {
            input:
                sample_id = sample_id,
                target = target,
                effective_target = effective_target,
                prefix = depth_prefix,
                pb_source_coverage = CalculateCoveragePacBio.coverage,
                ont_source_coverage = CalculateCoverageONT.coverage,
                pb_coverage = PacBioAtCoverage.achieved_coverage,
                ont_coverage = OntAtCoverage.achieved_coverage,
                pb_vcf = PacbioCall.vcf,
                ont_vcf = OntCall.vcf,
                shared_vcf = IdentifySharedVCF.vcf,
                pb_tree = HimitoLineage.pb_mutation_tree,
                ont_tree = HimitoLineage.ont_mutation_tree,
                pb_raw_matrix = PacbioCall.read_var_mat,
                ont_raw_matrix = OntCall.read_var_mat,
                pb_cleaned_matrix = HimitoLineage.pb_read_var_mat,
                ont_cleaned_matrix = HimitoLineage.ont_read_var_mat,
                docker = himito_docker
        }
    }

    call CombineCoverageMetrics {
        input:
            metrics = ScorePair.metrics,
            prefix = prefix,
            docker = himito_docker
    }

    output {
        Array[Int] target_coverages = desired_coverages
        Array[Float] effective_target_coverages = effective_target
        Float pb_source_coverage = CalculateCoveragePacBio.coverage
        Float ont_source_coverage = CalculateCoverageONT.coverage
        Array[Float] pb_achieved_coverage = PacBioAtCoverage.achieved_coverage
        Array[Float] ont_achieved_coverage = OntAtCoverage.achieved_coverage
        File coverage_metrics = CombineCoverageMetrics.combined
        Array[File] per_coverage_metrics = ScorePair.metrics
        Array[File] haplotype_support = ScorePair.support
        Array[File] shared_vcf = IdentifySharedVCF.vcf
        Array[File] shared_vcf_tbi = IdentifySharedVCF.vcf_tbi

        Array[File] pb_vcf = PacbioCall.vcf
        Array[File] pb_raw_matrix = PacbioCall.read_var_mat
        Array[File] pb_mutation_tree = HimitoLineage.pb_mutation_tree
        Array[File] pb_read_lineage = HimitoLineage.pb_read_lineage
        Array[File] pb_var_cooccurrence = HimitoLineage.pb_var_cooccurrence
        Array[File] pb_read_var_mat = HimitoLineage.pb_read_var_mat
        Array[File] pb_haplotype_map = HimitoLineage.pb_haplotype_map
        Array[File] pb_molecule_summary = HimitoLineage.pb_molecule_summary

        Array[File] ont_vcf = OntCall.vcf
        Array[File] ont_raw_matrix = OntCall.read_var_mat
        Array[File] ont_mutation_tree = HimitoLineage.ont_mutation_tree
        Array[File] ont_read_lineage = HimitoLineage.ont_read_lineage
        Array[File] ont_var_cooccurrence = HimitoLineage.ont_var_cooccurrence
        Array[File] ont_read_var_mat = HimitoLineage.ont_read_var_mat
        Array[File] ont_haplotype_map = HimitoLineage.ont_haplotype_map
        Array[File] ont_molecule_summary = HimitoLineage.ont_molecule_summary
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

struct DataTypeParameters {
    Int num_shards
    String map_preset
}

task CalculateCoverage {

    meta {
        description : "Subset a BAM file to a specified locus."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String? billing_project
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euo pipefail
        bam="~{bam}"
        if [[ "$bam" == gs://* ]]; then
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            export GCS_REQUESTER_PAYS_PROJECT="~{default="" billing_project}"
        fi
        samtools view -bhX "$bam" "~{bai}" "~{locus}" > "~{prefix}.bam"
        samtools index "~{prefix}.bam"
        samtools depth -aa -G 0xf04 -r "~{locus}" "~{prefix}.bam" |
            awk '{sum += $3} END {printf "%.8f\n", NR ? sum/NR : 0}' > coverage.txt
    >>>

    output {
        Float coverage = read_float("coverage.txt")
        File subsetbam =  "~{prefix}.bam"
        File subsetbai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/hifiasm:0.25.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Himito_quickstart {
    input {
        File bam
        File bai
        File reference_fa
        String prefix
        Int kmer_size
        String sample_id
        String data_type
        String chromo = "chrM"
        String himito_docker
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        # quick-start replaces the last extension when it names outputs.
        /Himito/target/release/Himito quick-start -i "~{bam}" \
                                                  -c "~{chromo}" \
                                                  -o "~{prefix}.out" \
                                                  -k ~{kmer_size} \
                                                  -r "~{reference_fa}" \
                                                  -s "~{sample_id}" \
                                                  -d ~{data_type} \
                                                  ~{extra_args}

    >>>  

    output {
        File graph = "~{prefix}.methyl.gfa"
        File methyl_bed = "~{prefix}.bed"
        File asm = "~{prefix}.fasta"
        File read_var_mat = "~{prefix}.matrix.csv"
        File read_methyl_mat = "~{prefix}.methylation_per_read.csv" 
        File numts_bam = "~{prefix}.numts.bam"
        File vcf = "~{prefix}.vcf"
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             himito_docker
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

task identify_shared_vcf {
    input {
        File pacbio_vcf
        File ont_vcf
        Boolean snv_only
        Boolean pass_only

        String prefix

        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        # Avoid /tmp space issues on some workers.
        export TMPDIR="${TMPDIR:-$PWD/bcftools_tmp}"
        mkdir -p "$TMPDIR"

        # Normalize to bgzip + tabix (Himito quick-start emits plain .vcf).
        # Optionally restrict to SNPs here (-v snps is a `view` option, not `isec`).
        bcftools view ~{if snv_only then "-v snps" else ""} ~{if pass_only then "-f PASS,." else ""} -Oz -o pacbio.vcf.gz "~{pacbio_vcf}"
        bcftools index -t pacbio.vcf.gz
        bcftools view ~{if snv_only then "-v snps" else ""} ~{if pass_only then "-f PASS,." else ""} -Oz -o ont.vcf.gz "~{ont_vcf}"
        bcftools index -t ont.vcf.gz

        # Sites present in both VCFs (exact CHROM/POS/REF/ALT); keep PacBio records (-w1).
        bcftools isec -c none -n=2 -w1 -Oz -o "~{prefix}.shared.vcf.gz" \
            ~{extra_args} \
            pacbio.vcf.gz ont.vcf.gz
        bcftools index -t "~{prefix}.shared.vcf.gz"
    >>>

    output {
        File vcf = "~{prefix}.shared.vcf.gz"
        File vcf_tbi = "~{prefix}.shared.vcf.gz.tbi"
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
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

task Himito_lineage {
    input {
        File shared_vcf
        File shared_vcf_tbi
        File pacbio_matrix
        File ont_matrix
        String prefix
        Float min_hf
        Float max_hf
        String himito_docker
        String? extra_args
        Boolean collapse_hp = false

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito lineage -m "~{pacbio_matrix}" \
                                              -v "~{shared_vcf}" \
                                              --min-hf ~{min_hf} \
                                              --max-hf ~{max_hf} \
                                              -d "pacbio" \
                                              -o "~{prefix}.pb" \
                                              ~{true="--include-non-snvs --collapse-hp" false="" collapse_hp} \
                                              ~{extra_args}

        /Himito/target/release/Himito lineage -m "~{ont_matrix}" \
                                        -v "~{shared_vcf}" \
                                        --min-hf ~{min_hf} \
                                        --max-hf ~{max_hf} \
                                        -d "ont-denoised" \
                                        -o "~{prefix}.ont" \
                                        ~{true="--include-non-snvs --collapse-hp" false="" collapse_hp} \
                                        ~{extra_args}

    >>>  

    output {
        File pb_mutation_tree = "~{prefix}.pb.mutation_tree.tsv"
        File pb_read_lineage = "~{prefix}.pb.read_lineage.nwk"
        File pb_var_cooccurrence = "~{prefix}.pb.variant_cooccurrence.tsv"
        File pb_read_var_mat = "~{prefix}.pb.cleaned_matrix_all_variants.csv"
        File pb_haplotype_map = "~{prefix}.pb.cleaned_haplotype_map.tsv" 
        File pb_molecule_summary = "~{prefix}.pb.molecule_summary.tsv" 

        File ont_mutation_tree = "~{prefix}.ont.mutation_tree.tsv"
        File ont_read_lineage = "~{prefix}.ont.read_lineage.nwk"
        File ont_var_cooccurrence = "~{prefix}.ont.variant_cooccurrence.tsv"
        File ont_read_var_mat = "~{prefix}.ont.cleaned_matrix_all_variants.csv"
        File ont_haplotype_map = "~{prefix}.ont.cleaned_haplotype_map.tsv" 
        File ont_molecule_summary = "~{prefix}.ont.molecule_summary.tsv" 
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             himito_docker
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

# The target/source fraction and fixed seed follow Downsample_experiment.wdl.
task DownsampleBam {
    input {
        File bam
        File bai
        String locus
        String prefix
        Float source_coverage
        Float target_coverage
        Int seed
        String docker
    }

    Int disk_size = 10 + ceil(3 * size(bam, "GB"))

    command <<<
        set -euo pipefail
        fraction=$(awk -v source="~{source_coverage}" -v target="~{target_coverage}" 'BEGIN {
            if (source <= 0 || target <= 0) {
                print "Source and target coverage must be positive" > "/dev/stderr"
                exit 1
            }
            printf "%.12f", (target >= source ? 1 : target/source)
        }')
        if [[ "$fraction" == 1.000000000000 ]]; then
            cp "~{bam}" "~{prefix}.bam"
            cp "~{bai}" "~{prefix}.bam.bai"
        else
            samtools view -b --subsample-seed "~{seed}" --subsample "$fraction" \
                -o "~{prefix}.bam" "~{bam}"
            samtools index "~{prefix}.bam"
        fi
        samtools depth -aa -G 0xf04 -r "~{locus}" "~{prefix}.bam" |
            awk '{sum += $3} END {
                if (!sum) {
                    print "Downsampled BAM has no callable depth at the locus" > "/dev/stderr"
                    exit 1
                }
                printf "%.8f\n", sum/NR
            }' > achieved_coverage.txt
    >>>

    output {
        File sampled_bam = "~{prefix}.bam"
        File sampled_bai = "~{prefix}.bam.bai"
        Float achieved_coverage = read_float("achieved_coverage.txt")
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " SSD"
        docker: docker
    }
}

task ScorePair {
    input {
        String sample_id
        Int target
        Float effective_target
        String prefix
        Float pb_source_coverage
        Float ont_source_coverage
        Float pb_coverage
        Float ont_coverage
        File pb_vcf
        File ont_vcf
        File shared_vcf
        File pb_tree
        File ont_tree
        File pb_raw_matrix
        File ont_raw_matrix
        File pb_cleaned_matrix
        File ont_cleaned_matrix
        String docker
    }

    command <<<
        set -euo pipefail
        PYTHONPATH=/opt/lineage_sim python - \
            "~{sample_id}" "~{target}" "~{effective_target}" "~{prefix}" \
            "~{pb_source_coverage}" "~{ont_source_coverage}" "~{pb_coverage}" "~{ont_coverage}" \
            "~{pb_vcf}" "~{ont_vcf}" "~{shared_vcf}" \
            "~{pb_tree}" "~{pb_raw_matrix}" "~{pb_cleaned_matrix}" \
            "~{ont_tree}" "~{ont_raw_matrix}" "~{ont_cleaned_matrix}" <<'PYTHON'
        import json
        import re
        import sys
        from score_lineage import detected_variants_with_hf_from_vcf, parse_mutation_tree
        from score_lineage_pair import compare_pair, write_metrics

        sample, target, effective_target, prefix, pb_source, ont_source, pb_depth, ont_depth, *paths = sys.argv[1:]
        depths = dict(sample_id=sample, target_coverage=int(target),
                      effective_target_coverage=float(effective_target),
                      pb_source_coverage=float(pb_source), ont_source_coverage=float(ont_source),
                      pb_achieved_coverage=float(pb_depth), ont_achieved_coverage=float(ont_depth),
                      pb_source_limited=float(pb_source) < int(target),
                      ont_source_limited=float(ont_source) < int(target))
        snv = re.compile(r"m\.\d+[ACGT]>[ACGT]")
        for label, path in zip(("pb", "ont", "shared"), paths[:3]):
            calls = detected_variants_with_hf_from_vcf(path)
            depths[f"{label}_pass_alleles"] = len(calls)
            depths[f"{label}_pass_snvs"] = sum(bool(snv.fullmatch(v)) for v in calls)
        pb_sites = parse_mutation_tree(paths[3]).keys() - {"ROOT"}
        ont_sites = parse_mutation_tree(paths[6]).keys() - {"ROOT"}
        depths.update(pb_tree_sites=len(pb_sites), ont_tree_sites=len(ont_sites),
                      shared_tree_sites=len(pb_sites & ont_sites))
        scores, support = compare_pair(paths[3:6], paths[6:9], 3, 10)
        with open(prefix + ".pair_metrics.tsv", "w") as f:
            write_metrics(f, {**depths, **scores})
        with open(prefix + ".haplotype_support.json", "w") as f:
            json.dump(support, f, indent=2)
            f.write("\n")
        PYTHON
    >>>

    output {
        File metrics = "~{prefix}.pair_metrics.tsv"
        File support = "~{prefix}.haplotype_support.json"
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk 20 SSD"
        docker: docker
    }
}

task CombineCoverageMetrics {
    input {
        Array[File]+ metrics
        String prefix
        String docker
    }

    command <<<
        set -euo pipefail
        python - "~{write_lines(metrics)}" "~{prefix}.coverage_metrics.tsv" <<'PYTHON'
        import csv
        import sys
        from pathlib import Path

        paths = Path(sys.argv[1]).read_text().splitlines()
        fields = None
        with open(sys.argv[2], "w", newline="") as output:
            for path in paths:
                with open(path, newline="") as source:
                    reader = csv.DictReader(source, delimiter="\t")
                    if fields is None:
                        fields = reader.fieldnames
                        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
                        writer.writeheader()
                    rows = list(reader)
                    if reader.fieldnames != fields or len(rows) != 1:
                        raise ValueError(f"Invalid per-coverage metrics: {path}")
                    writer.writerow(rows[0])
        PYTHON
    >>>

    output {
        File combined = "~{prefix}.coverage_metrics.tsv"
    }

    runtime {
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 10 SSD"
        docker: docker
    }
}
