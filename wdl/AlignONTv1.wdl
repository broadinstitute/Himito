version 1.0

workflow AlignONT {

    meta {
        desciption:
        "Given an fastq, align to a reference."
    }


    input {
        File fastq
        String prefix
        String gcs_out_root_dir
        File ref_fa
        String map_preset  # r10 for "lr:hq", r9 for "map-ont"
    }


    call Minimap2_ONT {
        input:
            read_fastq = fastq,
            prefix = prefix,
            ref_fasta = ref_fa,
            gcs_out_root_dir = gcs_out_root_dir,
            map_preset = map_preset
    }

}

task split_fastq {
    input {
        File read_fastq
        Int number_of_files
        Int disk_size
        String output_prefix
        RuntimeAttr? runtime_attr_override
        String disk_type = "SSD" 
    }
    command <<<
        mkdir -p ~{output_prefix}
        seqkit split ~{read_fastq} -p ~{number_of_files} -o ~{output_prefix}
        cd ~{output_prefix}
        ls > file_list.txt

    >>>
    output {
        Array[File] output_fastqs = read_lines("file_list.txt")
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/minimap2:v2.30"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

}


task Minimap2_ONT {
    input {
        File read_fastq
        String prefix
        File ref_fasta
        String gcs_out_root_dir
        String map_preset
        Int cpus = 16
        String disk_type = "SSD"  # options: SSD or LOCAL
        RuntimeAttr? runtime_attr_override
    }
    meta {
        descrpiton: "A wrapper to minimap2 for mapping & aligning (groups of) sequences to a reference. Note that this only works for reads belonging to a single readgroup."
    }

    # we limit the number of CPUs here because the number of input files could be huge (e.g. ONT fastqs)
    # also, the task is mostly CPU-bound (i.e. minimap2)
    Int pd_disk_size = 1 + 10*ceil(size(read_fastq, "GiB") + ceil(size(ref_fasta, "GiB")))
    Int local_disk_size = if(size(read_fastq, "GiB")>150) then 750 else 375
    Int disk_size = if('LOCAL'==disk_type) then local_disk_size else pd_disk_size

    # we limit the number of CPUs here because the number of input files could be huge (e.g. ONT fastqs)
    # also, the task is mostly CPU-bound (i.e. minimap2)
    Int mem = cpus * 5
    Int mm2_threads = cpus - 2


    command <<<
        set -euxo pipefail

        ############
        # parameter setting
        NUM_CPUS=$( cat /proc/cpuinfo | grep '^processor' | tail -n1 | awk '{print $NF+1}' )
        RAM_IN_GB=$( free -g | grep "^Mem" | awk '{print $2}' )
        echo "Memory info:"
        cat /proc/meminfo
        echo ""

        # mimimap2
        MAP_PARAMS="-ayYL --MD --eqx --cs -x ~{map_preset} -t ~{mm2_threads} -K4G               ~{ref_fasta}"


        # samtools sort
        SORT_PARAMS="-@4 -m4G --no-PG -o ~{prefix}.pre.bam"

        ############
        # minimap2
        /minimap2/minimap2-2.30/minimap2 $MAP_PARAMS ~{read_fastq} > tmp.sam


        ############
        # sort BAM
        time \
        samtools sort "${SORT_PARAMS}" tmp.sam -O BAM -o ~{prefix}.pre.bam
        rm tmp.sam

        ############
        # MD tag and index
        time \
        samtools calmd -@"${NUM_CPUS}" -b ~{prefix}.pre.bam ~{ref_fasta} > ~{prefix}.bam
        time \
        samtools index -@"${NUM_CPUS}" ~{prefix}.bam  

        gsutil cp ~{prefix}.bam ~{gcs_out_root_dir}
        gsutil cp ~{prefix}.bam.bai ~{gcs_out_root_dir}      

    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bam_index = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/minimap2:v2.30"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeBams {

    parameter_meta {
        bam_input: {localization_optional: true}
        bai_input: {localization_optional: true}
    }

    input{
        Array[File] bam_input
        Array[File] bai_input
        String prefix
        Int threads_num
        String gcs_out_root_dir
    }

    command <<<
        set -eux

        # we do single-sample phased VCFs localization ourselves
        mkdir -p bams
        time \
        gcloud storage cp ~{sep=" " bam_input} /mnt/disks/cromwell_root/bams/

        time \
        gcloud storage cp ~{sep=" " bai_input} /mnt/disks/cromwell_root/bams/

        # then merge, and safely assume all ssp-VCFs are sorted in the same order, on one chr
        cd bams
        ls *.bam > my_bams.txt

        samtools merge \
            --threads ~{threads_num} \
            -b my_bams.txt \
            -O BAM \
            -o ~{prefix}.bam
        samtools index ~{prefix}.bam

        gsutil cp ~{prefix}.bam ~{gcs_out_root_dir}
        gsutil cp ~{prefix}.bam.bai ~{gcs_out_root_dir}

    >>>

    output{
    }

    runtime {
        cpu: 16
        memory: "32 GiB"
        disks: "local-disk 500 LOCAL"
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/minimap2:v2.30"
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

