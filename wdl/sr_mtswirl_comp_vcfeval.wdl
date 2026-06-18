version 1.0
workflow VCF_Evaluation {
    input {
        File query_vcf
        File query_vcf_tbi
        File reference_fa
        File reference_fai
        File truth_vcf
        File truth_vcf_tbi
        String query_output_sample_name

        # String query_field
        # String? base_field
        # Float threshold

    }

    call SplitVCFs{
        input:
            query_vcf = query_vcf,
            query_vcf_tbi = query_vcf_tbi,
            base_vcf = truth_vcf,
            base_vcf_index = truth_vcf_tbi,   
    }

    scatter (sample in SplitVCFs.sample_list){
        call VCFEval {
            input:
                query_vcf = query_vcf,
                query_vcf_tbi = query_vcf_tbi,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                base_vcf = truth_vcf,
                base_vcf_index = truth_vcf_tbi,
                sampleid = sample,
                query_output_sample_name = sample + "_" + query_output_sample_name
        }
    }




    

    output {
        Array[File] summary_statistics = VCFEval.summary_statistics
        Array[Array[File]] evaled_vcf = VCFEval.combined_output
        Array[Array[File]] evaled_vcf_index = VCFEval.combined_output_index
    }
}



struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

task SplitVCFs {
    input {
        File query_vcf
        File query_vcf_tbi
        File base_vcf
        File base_vcf_index

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") ) + 50,
                                                  "cpu": 1, "memory": 4}
    }

    command <<<
        set -xeuo pipefail

        bcftools query -l ~{query_vcf} > query_samplelist.txt
        bcftools query -l ~{base_vcf} > base_samplelist.txt
        comm -12 <(sort query_samplelist.txt) <(sort base_samplelist.txt) > shared_samplelist.txt
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        Array[String] sample_list = read_lines("shared_samplelist.txt")
    }
}

task VCFEval {
    input {
        # Input VCF Files
        File query_vcf
        File query_vcf_tbi
        File reference_fa
        File reference_fai
        File base_vcf
        File base_vcf_index
        String sampleid
        String query_output_sample_name

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") + size(reference_fa, "GB")) + 50,
                                                  "cpu": 1, "memory": 4}
    }

    command <<<
        set -xeuo pipefail

        bcftools view -s ~{sampleid} ~{base_vcf} | bcftools view -i 'GT="alt"' -Oz -o ~{sampleid}.base.vcf.gz
        bcftools index -t ~{sampleid}.base.vcf.gz

        bcftools view -s ~{sampleid} ~{query_vcf} | bcftools view -i 'GT="alt"' -Oz -o ~{sampleid}.query.vcf.gz
        bcftools index -t ~{sampleid}.query.vcf.gz
        
        # rtg vcfeval
        rtg format -o rtg_ref ~{reference_fa}
        rtg vcfeval \
            -b ~{sampleid}.base.vcf.gz  \
            -c ~{sampleid}.query.vcf.gz \
            -o reg \
            -t rtg_ref \
            --squash-ploidy \
            --sample ALT,ALT 

        mkdir output_dir
        cp reg/summary.txt output_dir/~{query_output_sample_name}_summary.txt
        cp reg/*.vcf.gz* output_dir/
    
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
        Array[File] combined_output = glob("output_dir/*.vcf.gz")
        Array[File] combined_output_index = glob("output_dir/*.vcf.gz.tbi")
    }
}