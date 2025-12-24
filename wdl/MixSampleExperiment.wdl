version 1.0

workflow MixSamples {
    input {
        File first_donor_bam
        File first_donor_bai
        File second_donor_bam
        File second_donor_bai
        File first_donor_vcf
        File first_donor_tbi
        File second_donor_vcf
        File second_donor_tbi

        File reference_fa
        File reference_fai
        File reference_dict
        Array[Int] desiredCoverages
        Array[Float] first_proportion_list
        Int kmer_size = 21

        String sampleid
        String region = "chrM"
        String vcf_score_field_mitograph
        String query_field_mitograph
        String vcf_score_field_mitorsaw
        String query_field_mitorsaw
        String data_type
        Boolean run_mitosaw

    }

    call CalculateCoverage as first_donor_coverage{
        input:
            bam = first_donor_bam,
            bai = first_donor_bai,
            locus = region,
            prefix = sampleid
    }
    call CalculateCoverage as second_donor_coverage{
        input:
            bam = second_donor_bam,
            bai = second_donor_bai,
            locus = region,
            prefix = sampleid
    }

    call merge_vcf {
        input:
        first_donor_vcf = first_donor_vcf,
        first_donor_tbi = first_donor_tbi,
        second_donor_vcf = second_donor_vcf,
        second_donor_tbi = second_donor_tbi,
        reference_fa = reference_fa,
        prefix = sampleid,
    }

    scatter (desiredCoverage in desiredCoverages) {
        scatter (first_proportion in first_proportion_list)  {
            call downsampleBam {input:
                first_input_bam = first_donor_coverage.subsetbam,
                first_input_bam_bai = first_donor_coverage.subsetbai,
                second_input_bam = second_donor_coverage.subsetbam,
                second_input_bam_bai = second_donor_coverage.subsetbai,
                basename = sampleid,
                desiredCoverage = desiredCoverage,
                fraction = first_proportion,
                currentCoverage1 = first_donor_coverage.coverage,
                currentCoverage2 = second_donor_coverage.coverage,
                preemptible_tries = 0
            }

            call mix_sample {
                input:
                    first_donor_bam = downsampleBam.downsampled_bam_1,
                    first_donor_bai = downsampleBam.downsampled_bai_1,
                    second_donor_bam = downsampleBam.downsampled_bam_2,
                    second_donor_bai = downsampleBam.downsampled_bam_2,
                    prefix = sampleid
            }
            if (run_mitosaw) {
                # call mitorsaw
                call Mitorsaw {
                    input:
                        bam = mix_sample.merged_bam,
                        bai = mix_sample.merged_bai,
                        reference_fasta = reference_fa,
                        reference_fasta_fai = reference_fai,
                        prefix = sampleid
                }
                call VCFEval as Mitorsaw_Eval {
                    input:
                        query_vcf = Mitorsaw.vcf,
                        reference_fa = reference_fa,
                        reference_fai = reference_fai,
                        query_output_sample_name = sampleid + "_" + desiredCoverage + "_" + first_proportion,
                        base_vcf = merge_vcf.truth_vcf,
                        base_vcf_index = merge_vcf.truth_tbi,
                        vcf_score_field = vcf_score_field_mitorsaw,
                        query_field = query_field_mitorsaw,
                        threshold = 0,
                        fraction = first_proportion
                }
            }

            # call mitograph assembly
            call QuickStart {
                input:
                    bam = mix_sample.merged_bam,
                    bai = mix_sample.merged_bai,
                    reference_fa = reference_fa,
                    prefix = sampleid,
                    kmer_size = kmer_size,
                    sample_id = sampleid,
                    chromo = region,
                    data_type = data_type,

            }

            call VCFEval as Mitograph_Eval {
                input:
                    query_vcf = QuickStart.vcf,
                    reference_fa = reference_fa,
                    reference_fai = reference_fai,
                    query_output_sample_name = sampleid + "_" + desiredCoverage + "_" + first_proportion,
                    base_vcf = merge_vcf.truth_vcf,
                    base_vcf_index = merge_vcf.truth_tbi,
                    vcf_score_field = vcf_score_field_mitograph,
                    query_field = query_field_mitograph,
                    threshold = 0,
                    fraction = first_proportion
            }


        }

    }
    output {
        Array[Array[File]] mitograph_summary_file = Mitograph_Eval.summary_statistics
        Array[Array[File?]] mitorsaw_summary_file = Mitorsaw_Eval.summary_statistics
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

task CalculateCoverage {

    meta {
        description : "Subset a BAM file to a specified locus."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
        samtools depth -r ~{locus} ~{prefix}.bam | awk '{sum+=$3} END {print sum/NR}' > coverage.txt

    >>>

    output {
        Float coverage = read_float("coverage.txt")
        File subsetbam =  "~{prefix}.bam"
        File subsetbai = " ~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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


task downsampleBam {
  input {
    File first_input_bam
    File first_input_bam_bai
    File second_input_bam
    File second_input_bam_bai
    String basename
    Int desiredCoverage
    Float currentCoverage1
    Float currentCoverage2
    Float fraction
    Int? preemptible_tries
  }

  meta {
    description: "Uses Picard to downsample to desired coverage based on provided estimate of coverage."
  }
  parameter_meta {
  }

    Float second_fraction = 1 - fraction
    Float total_desired_reads = desiredCoverage
    Float first_desired_reads = total_desired_reads * fraction
    Float second_desired_reads = total_desired_reads * second_fraction
    
    Float scalingFactor1 = first_desired_reads / currentCoverage1
    Float scalingFactor2 = second_desired_reads / currentCoverage2

  command <<<
    set -eo pipefail
    echo scalingFactor1
    echo scalingFactor2
    gatk DownsampleSam -I ~{first_input_bam} -O ~{basename}_~{desiredCoverage}x_~{fraction}.bam -R 1 -P ~{scalingFactor1} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true
    gatk DownsampleSam -I ~{second_input_bam} -O ~{basename}_~{desiredCoverage}x_~{second_fraction}.bam -R 1 -P ~{scalingFactor2} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true


  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "8 GB"
    cpu: "2"
    disks: "local-disk 500 HDD"
    docker: "us.gcr.io/broad-gatk/gatk"
  }
  output {
    File downsampled_bam_1 = "~{basename}_~{desiredCoverage}x_~{fraction}.bam"
    File downsampled_bai_1 = "~{basename}_~{desiredCoverage}x_~{fraction}.bai"
    File downsampled_bam_2 = "~{basename}_~{desiredCoverage}x_~{second_fraction}.bam"
    File downsampled_bai_2 = "~{basename}_~{desiredCoverage}x_~{second_fraction}.bai"
  }
}


task mix_sample {

    meta {
        description : "sample reads to a specific proportion and mix the two bam"
    }

    parameter_meta {
        first_donor_bam: {
            description: "first donor bam to subset",
            localization_optional: false
        }
        first_donor_bai:    "index for first donor bam file"
        second_donor_bam: {
            description: "second donor bam to subset",
            localization_optional: false
        }
        second_donor_bai:    "index for second donor bam file"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File first_donor_bam
        File first_donor_bai
        File second_donor_bam
        File second_donor_bai
        String prefix

        RuntimeAttr? runtime_attr_override
    }


    command <<<
        set -euxo pipefail

        samtools merge -f ~{prefix}.merged.bam ~{first_donor_bam} ~{second_donor_bam}

        samtools sort -o ~{prefix}.merged.sorted.bam ~{prefix}.merged.bam
        
        samtools index ~{prefix}.merged.sorted.bam

    >>>

    output {
        File merged_bam = "~{prefix}.merged.sorted.bam"
        File merged_bai = "~{prefix}.merged.sorted.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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


struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

task VCFEval {
    input {
        # Input VCF Files
        File query_vcf
        File reference_fa
        File reference_fai
        String query_output_sample_name
        File base_vcf
        File base_vcf_index
        String vcf_score_field
        String query_field
        Float threshold = 0
        Float fraction

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") + size(reference_fa, "GB")) + 50,
                                                  "cpu": 8, "memory": 16}
    }
    # Float threshold1 = threshold - 0.01
    String query_info = "${query_field}\\>${threshold}" # extract heteroplasmic variants
    

    command <<<
        set -xeuo pipefail

        # Compress and Index vcf files
        bcftools view ~{query_vcf} -O z -o ~{query_vcf}.vcf.gz
        bcftools index -t ~{query_vcf}.vcf.gz

        # extract AF from query vcf file
        bcftools view -i  ~{query_info} ~{query_vcf}.vcf.gz -O z -o ~{query_output_sample_name}.query.~{fraction}.vcf.gz
        bcftools index -t ~{query_output_sample_name}.query.~{fraction}.vcf.gz
        
        # split multiallelic sites in the base_vcf
        bcftools norm \
                -f ~{reference_fa} \
                -m -both ~{base_vcf} \
                -O z \
                -o ~{base_vcf}.normed.vcf.gz 
        bcftools index -t ~{base_vcf}.normed.vcf.gz 
        
        # rtg vcfeval
        rtg format -o rtg_ref ~{reference_fa}
        rtg vcfeval \
            -b ~{base_vcf}.normed.vcf.gz  \
            -c ~{query_output_sample_name}.query.~{fraction}.vcf.gz \
            -o reg \
            -t rtg_ref \
            --squash-ploidy \
            --sample ALT,ALT \
            --vcf-score-field ~{vcf_score_field}

        mkdir output_dir
        cp reg/summary.txt output_dir/~{query_output_sample_name}.~{fraction}.summary.txt
    
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        File summary_statistics = "output_dir/~{query_output_sample_name}.~{fraction}.summary.txt"
        # File weighted_roc = "output_dir/weighted_roc.tsv.gz"
        # Array[File] combined_output = glob("output_dir/*.vcf.gz")
        # Array[File] combined_output_index = glob("output_dir/*.vcf.gz.tbi")
    }
}

task merge_vcf {

    meta {
        description : "merge two homoplasmic vcf as truth set"
    }

    input {
        File first_donor_vcf
        File first_donor_tbi
        File second_donor_vcf
        File second_donor_tbi
        File reference_fa
        String prefix
        

        RuntimeAttr? runtime_attr_override
    }


    command <<<
        set -euxo pipefail

        bcftools norm \
                -f ~{reference_fa} \
                -m -both ~{first_donor_vcf} \
                -O z \
                -o first.normed.vcf.gz 
        bcftools index -t first.normed.vcf.gz

        bcftools norm \
                -f ~{reference_fa} \
                -m -both ~{second_donor_vcf} \
                -O z \
                -o second.normed.vcf.gz 
        bcftools index -t second.normed.vcf.gz

        bcftools merge ~{first_donor_vcf} ~{second_donor_vcf} -O z -o ~{prefix}.merged.vcf.gz
        # bcftools view -H first.normed.vcf.gz
        # bcftools isec -C first.normed.vcf.gz second.normed.vcf.gz -w1 -O z -o ~{prefix}.merged.vcf.gz
        bcftools index -t ~{prefix}.merged.vcf.gz
        

    >>>

    output {
        File truth_vcf = "~{prefix}.merged.vcf.gz"
        File truth_tbi = "~{prefix}.merged.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "broadinstitute/gatk:4.6.1.0"
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

task Mutect2 {

    meta {
        description : "Call Mitochondrial variants using mutect2"
    }

    input {
        File bam
        File bai
        File reference_fasta
        File reference_fasta_fai
        File reference_fasta_dict
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail
        # samtools faidx ~{reference_fasta}
        # gatk CreateSequenceDictionary -R ~{reference_fasta} -O ~{reference_fasta}.dict
        gatk Mutect2 -R ~{reference_fasta} -L chrM --mitochondria-mode -I ~{bam} -O ~{prefix}.mutect2.vcf.gz

    >>>

    output {
        File vcf = "~{prefix}.mutect2.vcf.gz"
        File tbi = "~{prefix}.mutect2.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "broadinstitute/gatk:4.6.1.0"
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
