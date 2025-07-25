version 1.0

workflow MosdepthCoverage {

    meta {
        description: "Calculate coverage statistics using mosdepth."
    }
    parameter_meta {
        aligned_bam: "Aligned BAM file."
        aligned_bai: "Aligned BAM index file."
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File aligned_bam
        File aligned_bai
        String region
        String sampleid

        Int bin_length = 1000

        # Runtime parameters
        Int? preemptible_tries = 1
    }

    call SubsetBam as SubsetReads {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                locus = region,
                prefix = sampleid
    }

    call MosDepth {
        input:
            bam = SubsetReads.subset_bam,
            bai = SubsetReads.subset_bai,
            region = region,
            bin_length = bin_length,
            preemptible_tries = preemptible_tries,
    }


    output {
        File cov_stat = MosDepth.coverage_file
    }

}


task MosDepth {

    meta {
        description: "Calculate coverage using mosdepth."
    }
    parameter_meta {
        bam: "Aligned BAM file."
        bai: "Aligned BAM index file."
        bed: "BED file containing regions of interest."
        bin_length: "Length of bins to use for coverage calculation."
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File bam
        File bai
        File? bed
        String region
        Int bin_length

        # Runtime parameters
        Int? preemptible_tries = 3
    }

    String mosdepth_by_region = select_first([bed, bin_length])

    Int disk_size = 2 * ceil(size(bam, "GB") + size(bai, "GB"))
    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.coverage_over_bed"

    command {
        set -euxo pipefail

        # Create symbolic links for bam and bai in the current working directory
        ln -s ~{bam} ./~{basename}.bam
        ln -s ~{bai} ./~{basename}.bai

        mosdepth -t 4 -c ~{region} -x -Q 1 ~{prefix} ./~{basename}.bam
        gunzip ~{prefix}.per-base.bed.gz
    }

    output {
        File coverage_file = "~{prefix}.per-base.bed"
    }

    runtime {
        cpu:                    4
        memory:                 8 + " GiB"
        disks: "local-disk " +  100 + " HDD"
        preemptible:            preemptible_tries
        maxRetries:             1
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-mosdepth:0.3.2"
    }
}

task SubsetBam {

    meta {
        description : "Subset a BAM file to a specified locus, and extract sequence into fa"
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
    }

    input {
        File bam
        File bai
        String locus
        String prefix
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam

        samtools fasta ~{prefix}.bam > ~{prefix}.fasta
        samtools fastq ~{prefix}.bam > ~{prefix}.fastq
    >>>

    output {
        File subset_fasta = "~{prefix}.fasta"
        File subset_fastq = "~{prefix}.fastq"
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }
    
    runtime {
        cpu: 1
        memory: 4 + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }

}