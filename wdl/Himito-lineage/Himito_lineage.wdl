version 1.0

workflow Himito_lineage_eval {
    input {
        File pacbio_bam
        File pacbio_bam_bai
        File reference_fa
        String prefix
        String sample_id
        String? billing_project
        Int kmer_size
        Float min_hf
        Float max_hf
        String chromo
        String? extra_args
    }

    call CalculateCoverage as CalculateCoveragePacBio {
        input:
            bam = pacbio_bam,
            bai = pacbio_bam_bai,
            billing_project = billing_project,
            locus = chromo,
            prefix = sample_id
    }


    call Himito_quickstart as PacbioCall {
        input:
            bam = CalculateCoveragePacBio.subsetbam,
            bai = CalculateCoveragePacBio.subsetbai,
            reference_fa = reference_fa,
            prefix = prefix,
            kmer_size = kmer_size,
            sample_id = sample_id,
            chromo = chromo,
            data_type = "pacbio",
            extra_args = extra_args

    }

    call Himito_lineage as HimitoLineage {
        input:
            shared_vcf = PacbioCall.vcf,
            pacbio_matrix = PacbioCall.read_var_mat,
            prefix = prefix,
            min_hf = min_hf,
            max_hf = max_hf,
            extra_args = extra_args
    }

    output {
        Float chrM_coverage = CalculateCoveragePacBio.coverage
        File graph = PacbioCall.graph
        File methyl_bed = PacbioCall.methyl_bed
        File asm = PacbioCall.asm
        File read_var_mat = PacbioCall.read_var_mat
        File read_methyl_mat = PacbioCall.read_methyl_mat
        File numts_bam = PacbioCall.numts_bam
        File vcf = PacbioCall.vcf

        File pb_mutation_tree = HimitoLineage.pb_mutation_tree
        File pb_read_lineage = HimitoLineage.pb_read_lineage
        File pb_var_cooccurrence = HimitoLineage.pb_var_cooccurrence
        File pb_read_var_mat = HimitoLineage.pb_read_var_mat
        File pb_haplotype_map = HimitoLineage.pb_haplotype_map
        File pb_molecule_summary = HimitoLineage.pb_molecule_summary
    }
}


struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
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
        set -euxo pipefail
        
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        export GCS_REQUESTER_PAYS_PROJECT="~{billing_project}"
        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
        samtools depth -r ~{locus} ~{prefix}.bam | awk '{sum+=$3} END {print sum/NR}' > coverage.txt

    >>>

    output {
        Float coverage = read_float("coverage.txt")
        File subsetbam =  "~{prefix}.bam"
        File subsetbai = " ~{prefix}.bam.bai"
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
        Int maximal_depth = 1000
        String sample_id
        String data_type
        String chromo = "chrM"
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito quick-start -i ~{bam} \
                                                  -c ~{chromo} \
                                                  -o ~{prefix} \
                                                  -k ~{kmer_size} \
                                                  -r ~{reference_fa} \
                                                  -s ~{sample_id} \
                                                  -d ~{data_type} \
                                                  --maximal-mt-depth ~{maximal_depth} \
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
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/himito:dev"
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
        bcftools view -Oz -o pacbio.vcf.gz ~{pacbio_vcf}
        bcftools index -t pacbio.vcf.gz
        bcftools view -Oz -o ont.vcf.gz ~{ont_vcf}
        bcftools index -t ont.vcf.gz

        # Sites present in both VCFs (exact CHROM/POS/REF/ALT); keep PacBio records (-w1).
        # Pass e.g. -f PASS via extra_args to restrict to FILTER=PASS.
        bcftools isec -n=2 -w1 -Oz -o ~{prefix}.shared.vcf.gz \
            ~{extra_args} \
            pacbio.vcf.gz ont.vcf.gz
        bcftools index -t ~{prefix}.shared.vcf.gz
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
        File pacbio_matrix
        String prefix
        Float min_hf
        Float max_hf
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito lineage -m ~{pacbio_matrix} \
                                              -v ~{shared_vcf} \
                                              --min-hf ~{min_hf} \
                                              --max-hf ~{max_hf} \
                                              -d "pacbio" \
                                              -o ~{prefix}.pb \
                                              ~{extra_args}

    >>>  

    output {
        File pb_mutation_tree = "~{prefix}.pb.mutation_tree.tsv"
        File pb_read_lineage = "~{prefix}.pb.read_lineage.nwk"
        File pb_var_cooccurrence = "~{prefix}.pb.variant_cooccurrence.tsv"
        File pb_read_var_mat = "~{prefix}.pb.cleaned_matrix.csv"
        File pb_haplotype_map = "~{prefix}.pb.cleaned_haplotype_map.tsv" 
        File pb_molecule_summary = "~{prefix}.pb.molecule_summary.tsv" 
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/himito:dev"
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