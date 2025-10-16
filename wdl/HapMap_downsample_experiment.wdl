version 1.0

workflow DownsampleExperiment {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File truth_vcf
        File truth_tbi

        File reference_fa
        File reference_fai
        File reference_dict
        Array[Int] desiredCoverages
        Int kmer_size = 21

        String sampleid
        String data_type
        String region = "chrM"

    }

    call CalculateCoverage as coverage{
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = region,
            prefix = sampleid
    }

    scatter (desiredCoverage in desiredCoverages) {
        
        call downsampleBam {input:
            input_bam = coverage.subsetbam,
            input_bam_bai = coverage.subsetbai,
            basename = sampleid,
            desiredCoverage = desiredCoverage,
            currentCoverage = coverage.coverage,
            preemptible_tries = 0
        }


        call Mitorsaw {
            input:
                bam = downsampleBam.downsampled_bam,
                bai = downsampleBam.downsampled_bai,
                reference_fasta = reference_fa,
                reference_fasta_fai = reference_fai,
                prefix = sampleid + "_" + desiredCoverage
        }

        call Filter {
            input:
                bam = downsampleBam.downsampled_bam,
                bai = downsampleBam.downsampled_bai,
                prefix = sampleid
        }

        call Build {
            input:
                bam = Filter.mt_bam,
                reference = reference_fa,
                prefix = desiredCoverage,
                kmer_size = kmer_size,
                sampleid = sampleid
        }

        call Call {
            input:
                graph_gfa = Build.graph,
                reference_fa = reference_fa,
                prefix = desiredCoverage,
                kmer_size = kmer_size,
                data_type = data_type,
                sampleid=sampleid
        }

        call VCFEval as Himito_Eval {
            input:
                query_vcf = Call.vcf,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                base_vcf = truth_vcf,
                base_vcf_index = truth_tbi,
                query_output_sample_name = sampleid + "_" + desiredCoverage + "_Himito",
        }
    
        call VCFEval as Mitorsaw_Eval {
            input:
                query_vcf = Mitorsaw.vcf,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                base_vcf = truth_vcf,
                base_vcf_index = truth_tbi,
                query_output_sample_name = sampleid + "_" + desiredCoverage + "_Mitorsaw",
        }
        

    }
    output {
        Array[File] Himito_summary_file = Himito_Eval.summary_statistics
        Array[File] mitorsaw_summary_file = Mitorsaw_Eval.summary_statistics
        Array[File] Himito_vcf = Call.vcf
        Array[File] Mitorsaw_vcf = Mitorsaw.vcf
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
        File input_bam
        File input_bam_bai
        String basename
        Int desiredCoverage
        Float currentCoverage
        Int? preemptible_tries
    }

    meta {
        description: "Uses Picard to downsample to desired coverage based on provided estimate of coverage."
    }

    parameter_meta {
    }

    Float scalingFactor = desiredCoverage / currentCoverage


    command <<<
        set -eo pipefail

        gatk DownsampleSam -I ~{input_bam} -O ~{basename}_~{desiredCoverage}x.bam -R 7 -P ~{scalingFactor} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true


    >>>
    runtime {
        preemptible: select_first([preemptible_tries, 5])
        memory: "8 GB"
        cpu: "2"
        disks: "local-disk 500 HDD"
        docker: "us.gcr.io/broad-gatk/gatk"
    }
    output {
        File downsampled_bam = "~{basename}_~{desiredCoverage}x.bam"
        File downsampled_bai = "~{basename}_~{desiredCoverage}x.bai"
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

        /Himito/target/release/Himito build -i ~{bam} -k ~{kmer_size} -r ~{reference} -o ~{sampleid}.~{prefix}.gfa 

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.gfa"
    }

    runtime {
        docker: "hangsuunc/himito:v1"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task Call {

    input {
        File graph_gfa
        File reference_fa
        String prefix
        String sampleid
        String data_type
        Int kmer_size
    }
    
    

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito call -g ~{graph_gfa} -r ~{reference_fa} -k ~{kmer_size} -s ~{sampleid} -d ~{data_type} -o ~{sampleid}.~{prefix}.vcf
        ls

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.annotated.gfa"
        File matrix = "~{sampleid}.~{prefix}.matrix.csv"
        File vcf = "~{sampleid}.~{prefix}.vcf"
    }

    runtime {
        docker: "hangsuunc/himito:v1"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}


struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
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
        cp ~{bam} ~{prefix}.bam
        cp ~{bai} ~{prefix}.bam.bai
        /mitorsaw-v0.2.0-x86_64-unknown-linux-gnu/mitorsaw haplotype \
            --reference ~{reference_fasta} \
            --bam  ~{prefix}.bam \
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
        cpu_cores:          4,
        mem_gb:             16,
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