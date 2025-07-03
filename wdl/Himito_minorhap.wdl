version 1.0

workflow Himito_call {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        Int ref_length
        String prefix
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

    call Minorhap {
        input:
            graph_gfa = Build.graph,
            ref_length = ref_length,
            bin_size = 500,
            pad_size = 50,
            min_read_ratio = 0.01,
            prefix = "minor_hap",
            sampleid= sampleid
    }

    output {
        File graph = Build.graph
        File haplotypes = Minorhap.minor_haplotype
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
        docker: "hangsuunc/himito:v3"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
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

task Minorhap {

    input {
        File graph_gfa
        Int ref_length
        Int bin_size
        Int pad_size
        Float min_read_ratio
        String prefix
        String sampleid
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito minorhap -g ~{graph_gfa} -r ~{ref_length} -b ~{bin_size} -p ~{pad_size} -m ~{min_read_ratio} -s ~{sampleid} -o ~{sampleid}.~{prefix}.fasta

    >>>

    output {
        # File graph = "~{prefix}.annotated.gfa"
        File minor_haplotype = "~{sampleid}.~{prefix}.fasta"
    }

    runtime {
        docker: "hangsuunc/himito:v3"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}
