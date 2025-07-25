version 1.0

workflow Mitorsaw_call {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        File reference_fai
        String sampleid

    }
    call Mitorsaw {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            reference_fasta = reference_fa,
            reference_fasta_fai = reference_fai,
            prefix = sampleid
    }


    output {
        File vcf = Mitorsaw.vcf
        File vcf_index = Mitorsaw.tbi
        File stats = Mitorsaw.stats
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

task Mitorsaw {

    meta {
        description : "Call Mitochondrial variants using mitorsaw"
    }

    input {
        File bam
        File bai
        File reference_fasta
        File reference_fasta_fai
        String prefix
        Float min_af = 0.01

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail
        /mitorsaw-v0.2.0-x86_64-unknown-linux-gnu/mitorsaw haplotype \
            --reference ~{reference_fasta} \
            --bam ~{bam} \
            --minimum-maf ~{min_af} \
            --output-vcf ~{prefix}.mitorsaw.vcf.gz \
            --output-hap-stats ~{prefix}.mitorsaw.stat

    >>>

    output {
        File vcf = "~{prefix}.mitorsaw.vcf.gz"
        File tbi = "~{prefix}.mitorsaw.vcf.gz.tbi"
        File stats = "~{prefix}.mitorsaw.stat"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "hangsuunc/mitorsaw:v1"
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
