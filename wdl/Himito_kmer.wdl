version 1.0

workflow Himito_kmer {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        File reference_fai
        File truth_vcf
        File truth_tbi
        String prefix
        String sampleid
        String data_type
        Array[Int] kmer_list

    }

    scatter (k in kmer_list) {
        call QuickStart {
            input:
                bam = whole_genome_bam,
                bai = whole_genome_bai,
                reference_fa = reference_fa,
                prefix = prefix + "_" + k,
                kmer_size = k,
                sample_id = sampleid,
                chromo = "chrM",
                data_type = data_type,

        }


        call VCFEval as Himito_Eval {
            input:
                query_vcf = QuickStart.vcf,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                base_vcf = truth_vcf,
                base_vcf_index = truth_tbi,
                query_output_sample_name = sampleid + "_" + k + "_Himito",
        }

    }
    


    output {
        Array[File] fasta_file = QuickStart.asm
        Array[File] Himito_summary_file = Himito_Eval.summary_statistics
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

struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
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
        disks: "local-disk 200 SSD"
    }
}


task VCFEval {
    input {
        # Input VCF Files
        File query_vcf
        File reference_fa
        File reference_fai
        File base_vcf
        File base_vcf_index
        String query_output_sample_name

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") + size(reference_fa, "GB")) + 50,
                                                  "cpu": 8, "memory": 16}
    }

    command <<<
        set -xeuo pipefail

        # Compress and Index vcf files
        bcftools view ~{query_vcf} -O z -o ~{query_output_sample_name}.vcf.gz
        bcftools index -t ~{query_output_sample_name}.vcf.gz

        # split multiallelic sites in the base_vcf
        bcftools norm \
                -f ~{reference_fa} \
                -m -both ~{base_vcf} \
                -O z \
                -o ~{query_output_sample_name}.base.normed.vcf.gz 
        bcftools index -t ~{query_output_sample_name}.base.normed.vcf.gz 
        
        # rtg vcfeval
        rtg format -o rtg_ref ~{reference_fa}
        rtg vcfeval \
            -b ~{query_output_sample_name}.base.normed.vcf.gz  \
            -c ~{query_output_sample_name}.vcf.gz \
            -o reg \
            -t rtg_ref \
            --squash-ploidy \
            --sample ALT,ALT 

        mkdir output_dir
        cp reg/summary.txt output_dir/~{query_output_sample_name}_summary.txt
        
    
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        File summary_statistics = "output_dir/~{query_output_sample_name}_summary.txt"
        
    }
}
