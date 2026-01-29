version 1.0

workflow Himito_build {
    input {
        File fasta
        File reference_fa
        String prefix
        String reference_header
        String sampleid
        String data_type
        Int kmer_size

    }

    call Build {
        input:
            bam = fasta,
            reference = reference_fa,
            prefix = prefix,
            kmer_size = kmer_size,
            sampleid = sampleid
    }




    output {
        File graph = Build.graph
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

        if [[ ~{bam} == *.gz ]]; then
            gunzip -c ~{bam} > ~{prefix}.fasta
            /Himito/target/release/Himito build -k ~{kmer_size} -r ~{reference} -i ~{prefix}.fasta -o ~{sampleid}.~{prefix}.gfa
        else
            /Himito/target/release/Himito build -k ~{kmer_size} -r ~{reference} -i ~{bam} -o ~{sampleid}.~{prefix}.gfa
        fi 
        

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.gfa"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}
