version 1.0

workflow Himito_call {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        String prefix
        String reference_header
        String sampleid
        Int kmer_size

    }
    call Filter {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            prefix = prefix
    }

    call Build {
        input:
            bam = Filter.mt_bam,
            reference = reference_fa,
            prefix = prefix,
            kmer_size = kmer_size,
            sampleid = sampleid
    }

    call Call {
        input:
            graph_gfa = Build.graph,
            reference = reference_fa,
            reference_name = reference_header,
            prefix = prefix,
            kmer_size = kmer_size,
            sampleid=sampleid
    }


    output {
        File graph = Build.graph
        File vcf_file = Call.vcf
    }
}

task Filter {
    input {
        File bam
        File bai
        String prefix
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito filter -i ~{bam} -c chrM -m ~{prefix}_mt.bam -n ~{prefix}_numts.bam
    >>>

    output {
        File mt_bam = "~{prefix}_mt.bam"
        File numts_bam = "~{prefix}_numts.bam"
    }

    runtime {
        docker: "hangsuunc/himito:v1"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 300 SSD"
    }
}

task Build {
    input {
        File bam
        File reference
        String prefix
        String sampleid
        Int kmer_size
    }

    command <<<
        set -euxo pipefail

        /Himito/target/release/Himito build -k ~{kmer_size} -r ~{reference} -i ~{bam} -o ~{sampleid}.~{prefix}.gfa

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.gfa"
    }

    runtime {
        docker: "hangsuunc/himito:v3"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task Call {

    input {
        File graph_gfa
        File reference
        String reference_name
        String prefix
        String sampleid
        Int kmer_size
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito call -g ~{graph_gfa} -r ~{reference_name} -k ~{kmer_size} -s ~{sampleid} -o ~{sampleid}.~{prefix}.vcf

    >>>

    output {
        # File graph = "~{prefix}.annotated.gfa"
        File vcf = "~{sampleid}.~{prefix}.vcf"
    }

    runtime {
        docker: "hangsuunc/himito:v1"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}
