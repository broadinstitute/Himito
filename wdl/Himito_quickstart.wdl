version 1.0

workflow Himito {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        String prefix
        String sample_id
        String data_type
        Int kmer_size
        String chromo
    }

    call QuickStart {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            reference_fa = reference_fa,
            prefix = prefix,
            kmer_size = kmer_size,
            sample_id = sample_id,
            chromo = chromo,
            data_type = data_type,

    }

    output {
        File graph = QuickStart.graph
        File methyl_bed = QuickStart.methyl_bed
        File asm = QuickStart.asm
        File read_var_mat = QuickStart.read_var_mat
        File read_methyl_mat = QuickStart.read_methyl_mat
        File numts_bam = QuickStart.numts_bam
        File vcf = QuickStart.vcf
    }
}

task QuickStart {
    input {
        File bam
        File bai
        File reference_fa
        String prefix
        Int kmer_size
        String sample_id
        String chromo = "chrM"
        String data_type = "pacbio"
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito quick-start -i ~{bam} -c ~{chromo} -o ~{prefix} -k ~{kmer_size} -r ~{reference_fa} -s ~{sample_id} -d ~{data_type}
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

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 300 SSD"
    }
}
