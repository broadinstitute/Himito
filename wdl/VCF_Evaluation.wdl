version 1.0
workflow VCF_Evaluation {
    input {
        File whole_read_bam
        File whole_read_bai
        File reference_fa
        File reference_fai
        File truth_vcf
        File truth_vcf_tbi
        String prefix
        String sampleid
        Int kmer_size
        Float currentCoverage
        Array[Int] desiredCoverages
        String query_field
        String? base_field
        Float threshold

    }
    
    call SubsetBam as SubsetReads {
            input:
                bam = whole_read_bam,
                bai = whole_read_bai,
                locus = "chrM",
                prefix = sampleid
    }
    
    # downsample Coverage
    scatter (desiredCoverage in desiredCoverages) {
        call downsampleBam {input:
            input_bam = SubsetReads.subset_bam,
            input_bam_bai = SubsetReads.subset_bai,
            basename = sampleid,
            desiredCoverage = desiredCoverage,
            currentCoverage = currentCoverage,
            preemptible_tries = 0
        }

        # call mitograph assembly
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
                prefix = prefix,
                kmer_size = kmer_size,
                sampleid = sampleid
        }

        call Call {
            input:
                graph_gfa = Build.graph,
                reference_fa = reference_fa,
                prefix = desiredCoverage,
                kmer_size = kmer_size,
                sampleid=sampleid
        }

        call VCFEval {
            input:
                query_vcf = Call.vcf,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                query_output_sample_name = desiredCoverage,
                base_vcf = truth_vcf,
                base_vcf_index = truth_vcf_tbi,
                query_field = query_field,
                base_field = base_field,
                threshold = threshold
        }


    }

    output {
        Array[File] summary_statistics = VCFEval.summary_statistics
        Array[File] vcf = Call.vcf
    }
}

task downsampleBam {
  input {
    File input_bam
    File input_bam_bai
    String basename
    Int desiredCoverage
    Float currentCoverage
    Float scalingFactor = desiredCoverage / currentCoverage

    Int? preemptible_tries
  }

  meta {
    description: "Uses Picard to downsample to desired coverage based on provided estimate of coverage."
  }
  parameter_meta {
    basename: "Input is a string specifying the sample name which will be used to locate the file on gs."
    downsampled_bam: "Output is a bam file downsampled to the specified mean coverage."
    downsampled_bai: "Output is the index file for a bam file downsampled to the specified mean coverage."
    desiredCoverage: "Input is an integer of the desired approximate coverage in the output bam file."
  }
  command <<<
    set -eo pipefail
    gatk DownsampleSam -I ~{input_bam} -O ~{basename}_~{desiredCoverage}x.bam -R 1 -P ~{scalingFactor} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true


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
        Int kmer_size
    }
    
    

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito call -g ~{graph_gfa} -r ~{reference_fa} -k ~{kmer_size} -s ~{sampleid} -o ~{sampleid}.~{prefix}.vcf
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

task Methyl {

    input {
        File graph_gfa
        File bam
        String sampleid
        String prefix = "methyl"
        Float min_prob
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito methyl -g ~{graph_gfa} -p ~{min_prob} -b ~{bam} -o ~{sampleid}.~{prefix}.bed

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.methyl.gfa"
        File matrix = "~{sampleid}.~{prefix}.methylation_per_read.csv"
        File bed = "~{sampleid}.~{prefix}.bed"
    }

    runtime {
        docker: "hangsuunc/himito:v1"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task Asm {

    input {
        File graph_gfa
        String prefix
        String sampleid
    }

    command <<<
        set -euxo pipefail
        /Himito/target/release/Himito asm -g ~{graph_gfa} -s ~{sampleid} -o ~{sampleid}.~{prefix}.fasta

    >>>

    output {
        File fasta = "~{sampleid}.~{prefix}.fasta"
    }

    runtime {
        docker: "hangsuunc/himito:v1"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task SubsetBam {

    meta {
        description : "Subset a BAM file to a specified locus, and extract sequence into fa"
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
    }

    input {
        File bam
        File bai
        String locus
        String prefix
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam

        samtools fasta ~{prefix}.bam > ~{prefix}.fasta
        samtools fastq ~{prefix}.bam > ~{prefix}.fastq
    >>>

    output {
        File subset_fasta = "~{prefix}.fasta"
        File subset_fastq = "~{prefix}.fastq"
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }
    
    runtime {
        cpu: 1
        memory: 4 + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
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
        String query_field
        String? base_field
        Float threshold

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") + size(reference_fa, "GB")) + 50,
                                                  "cpu": 8, "memory": 16}
    }
    String query_info = "${query_field}\\>${threshold}"
    String base_info = if defined(base_field) then "${base_field}\\>${threshold}" else ""

    command <<<
        set -xeuo pipefail

        # Compress and Index vcf files
        bcftools view ~{query_vcf} -O z -o ~{query_vcf}.vcf.gz
        bcftools index -t ~{query_vcf}.vcf.gz

        # extract AF from query vcf file
        bcftools view -i  ~{query_info} ~{query_vcf}.vcf.gz -O z -o ~{query_output_sample_name}.query.~{threshold}.vcf.gz
        bcftools index -t ~{query_output_sample_name}.query.~{threshold}.vcf.gz

        # Check if base_field is defined and apply filter if it is
        if [ "~{default='' base_field}" != "" ]; then
            bcftools view -i ~{base_info} ~{base_vcf} -O z -o ~{query_output_sample_name}.base.~{threshold}.vcf.gz
        else
            cp ~{base_vcf} ~{query_output_sample_name}.base.~{threshold}.vcf.gz
        fi
        bcftools index -t ~{query_output_sample_name}.base.~{threshold}.vcf.gz

        
        # split multiallelic sites in the base_vcf
        bcftools norm \
                -f ~{reference_fa} \
                -m -both ~{query_output_sample_name}.base.~{threshold}.vcf.gz \
                -O z \
                -o ~{query_output_sample_name}.base.~{threshold}.normed.vcf.gz 
        bcftools index -t ~{query_output_sample_name}.base.~{threshold}.normed.vcf.gz 
        
        # rtg vcfeval
        rtg format -o rtg_ref ~{reference_fa}
        rtg vcfeval \
            -b ~{query_output_sample_name}.base.~{threshold}.normed.vcf.gz  \
            -c ~{query_output_sample_name}.query.~{threshold}.vcf.gz \
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
        # File weighted_roc = "output_dir/weighted_roc.tsv.gz"
        # Array[File] combined_output = glob("output_dir/*.vcf.gz")
        # Array[File] combined_output_index = glob("output_dir/*.vcf.gz.tbi")
    }
}
