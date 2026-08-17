version 1.0

workflow Modkit_entropy{
    meta{
        description: "a workflow that using Modkit to extract methylated signal from MM/ML tagged bams"
    }

    input{
        File bam
        File bai
        File? bed
        File ref
        File ref_index
        String? locus
        String output_prefix
        Int nthreads
    }
    if (defined(locus)) {
        call SubsetBam{input:
            bam = bam,
            bai = bai,
            locus = select_first([locus]),
            prefix = output_prefix
        }

    }

    call Filter {
        input:
            bam = select_first([SubsetBam.subset_bam, bam]),
            bai = select_first([SubsetBam.subset_bai, bai]),
            prefix = output_prefix
    }

    call Modkit_entropy_task{input:
        bam = select_first([Filter.mt_bam]),
        region_bed = bed,
        outputprefix = output_prefix,
        num_threads = nthreads,
        ref = ref,
        ref_index = ref_index
    }
    
    output{
        File region_files = Modkit_entropy_task.regions
        File entropy_bedgraph = Modkit_entropy_task.entropy_bedgraph
        File log_file = Modkit_entropy_task.log_file
    }
}


task Modkit_entropy_task{
    meta{
        description: "a task that using Modkit to extract methylated entropy from MM/ML tagged bams"
    }
    input{
        File bam
        File? bai
        File? region_bed
        String outputprefix
        File ref
        File ref_index
        Int num_threads

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        mkdir -p output

        # Cromwell often localizes SubsetBam outputs already named *.bam.bai next to
        # the BAM (same inode/path). Blind ln -sf then fails with "are the same file"
        # under set -e, so modkit never runs and glob("output/*") is empty.

        samtools index ~{bam}
        
        if [[ ! -e "~{ref}.fai" ]]; then
            ln -s ~{ref_index} "~{ref}.fai"
        fi

        modkit entropy \
            --in-bam ~{bam} \
            -o output \
            ~{if defined(region_bed) then "--regions " + region_bed else ""} \
            --cpg \
            --ref ~{ref} \
            --threads ~{num_threads} \
            --log-filepath ~{outputprefix}.log

        mv output/regions.bed ~{outputprefix}.regions.bed
        mv output/windows.bedgraph ~{outputprefix}.windows.bedgraph
    >>>

    output{
        File regions = "~{outputprefix}.regions.bed"
        File entropy_bedgraph = "~{outputprefix}.windows.bedgraph"
        File log_file = "~{outputprefix}.log"
    }

    Int disk_size = 100 + ceil(2 * (size(bam, "GiB")))
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       50,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/methylationanalysis:v2"
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


task SubsetBam {

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
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task Filter {
    input {
        File bam
        File bai
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito filter -i ~{bam} -c chrM -m ~{prefix}_mt.bam -n ~{prefix}_numts.bam
    >>>

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    output {
        File mt_bam = "~{prefix}_mt.bam"
        File numts_bam = "~{prefix}_numts.bam"
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1.1.2"
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
