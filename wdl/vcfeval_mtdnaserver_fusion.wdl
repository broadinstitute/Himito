version 1.0
workflow VCF_Evaluation {
    input {
        File query_vcf
        File query_vcf_1
        File reference_fa
        File reference_fai
        File truth_vcf
        File truth_vcf_tbi
        String query_output_sample_name

        # String query_field
        # String? base_field
        # Float threshold

    }

    call VCFEval {
        input:
            query_vcf = query_vcf,
            query_vcf_1 = query_vcf_1,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            base_vcf = truth_vcf,
            base_vcf_index = truth_vcf_tbi,
            query_output_sample_name = query_output_sample_name
    }
    

    output {
        File summary_statistics = VCFEval.summary_statistics
        File merged_vcf = VCFEval.merged_vcf
        File merged_vcf_index = VCFEval.merged_vcf_index
    }
}



struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

task VCFEval {
    input {
        # Input VCF Files
        File query_vcf
        File query_vcf_1
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

        bcftools view ~{query_vcf_1} -O z -o ~{query_output_sample_name}_1.vcf.gz
        bcftools index -t ~{query_output_sample_name}_1.vcf.gz

        # Some mtDNA-Server VCFs contain malformed FORMAT/AF values (Number=A mismatch at some sites),
        # which makes `bcftools merge` abort. AF is not used downstream here, so drop it before merge.
        bcftools annotate -x FORMAT/AF ~{query_output_sample_name}.vcf.gz -Oz -o ~{query_output_sample_name}.noAF.vcf.gz
        bcftools index -t ~{query_output_sample_name}.noAF.vcf.gz
        bcftools annotate -x FORMAT/AF ~{query_output_sample_name}_1.vcf.gz -Oz -o ~{query_output_sample_name}_1.noAF.vcf.gz
        bcftools index -t ~{query_output_sample_name}_1.noAF.vcf.gz

        # remove filter
        bcftools merge  ~{query_output_sample_name}.noAF.vcf.gz ~{query_output_sample_name}_1.noAF.vcf.gz -Oz -o ~{query_output_sample_name}_merged.vcf.gz
        bcftools index -t ~{query_output_sample_name}_merged.vcf.gz

        # Normalize FILTER column for downstream comparison: set all records to PASS.
        bcftools view ~{query_output_sample_name}_merged.vcf.gz | \
            awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$7="PASS"; print}' | \
            bgzip -c > ~{query_output_sample_name}_merged.pass.vcf.gz
        bcftools index -t ~{query_output_sample_name}_merged.pass.vcf.gz
        
        # rtg vcfeval
        rtg format -o rtg_ref ~{reference_fa}
        rtg vcfeval \
            -b ~{base_vcf}  \
            -c ~{query_output_sample_name}_merged.pass.vcf.gz \
            -o reg \
            -t rtg_ref \
            --squash-ploidy \
            --sample ALT,ALT 

        mkdir output_dir
        cp reg/summary.txt output_dir/~{query_output_sample_name}_summary.txt
        # cp reg/weighted_roc.tsv.gz output_dir/
        # cp reg/*.vcf.gz* output_dir/
        # cp output_dir/output.vcf.gz output_dir/~{query_output_sample_name}.vcf.gz
        # cp output_dir/output.vcf.gz.tbi output_dir/~{query_output_sample_name}.vcf.gz.tbi

    
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
        File merged_vcf = "~{query_output_sample_name}_merged.vcf.gz"
        File merged_vcf_index = "~{query_output_sample_name}_merged.vcf.gz.tbi"
        # File weighted_roc = "output_dir/weighted_roc.tsv.gz"
        # Array[File] combined_output = glob("output_dir/*.vcf.gz")
        # Array[File] combined_output_index = glob("output_dir/*.vcf.gz.tbi")
    }
}