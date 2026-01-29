version 1.0

workflow Himito_call {
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
        Array[Float] max_methylation_threshold
        Int kmer_size

    }

    scatter (t in max_methylation_threshold) {
        call Filter {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            prefix = prefix + "_" + t,
            threshold = t
        }

        call Build {
            input:
                bam = Filter.mt_bam,
                reference = reference_fa,
                prefix = prefix + "_" + t,
                kmer_size = kmer_size,
                sampleid = sampleid
        }

        call Call {
            input:
                graph_gfa = Build.graph,
                reference = reference_fa,
                prefix = prefix + "_" + t,
                kmer_size = kmer_size,
                sampleid=sampleid,
                data_type = data_type
        }


        call VCFEval as Himito_Eval {
            input:
                query_vcf = Call.vcf,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                base_vcf = truth_vcf,
                base_vcf_index = truth_tbi,
                query_output_sample_name = sampleid + "_" + t + "_Himito",
        }

    }
    


    output {
        Array[File] vcf_file = Call.vcf
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

task Filter {
    input {
        File bam
        File bai
        String prefix
        Float threshold
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito filter -i ~{bam} -c chrM -m ~{prefix}_mt.bam -n ~{prefix}_numts.bam -f ~{threshold}
    >>>

    output {
        File mt_bam = "~{prefix}_mt.bam"
        File numts_bam = "~{prefix}_numts.bam"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1"
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
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task Call {

    input {
        File graph_gfa
        File reference
        String prefix
        String sampleid
        String data_type
        Int kmer_size
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito call -g ~{graph_gfa} -r ~{reference} -k ~{kmer_size} -d ~{data_type} -s ~{sampleid} -o ~{sampleid}.~{prefix}.vcf

    >>>

    output {
        # File graph = "~{prefix}.annotated.gfa"
        File vcf = "~{sampleid}.~{prefix}.vcf"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/himito:v1"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
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
