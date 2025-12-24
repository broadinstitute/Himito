version 1.0
workflow CompareFasta {
    input {
        File whole_hap1_bam
        File whole_hap1_bai
        File whole_hap2_bam
        File whole_hap2_bai
        File whole_read_bam
        File whole_read_bai
        File reference_fa
        String prefix
        String sampleid
        Int kmer_size
        Float currentCoverage
        Array[Int] desiredCoverages

    }

    # extract truth assembly
    call SubsetBam as Hap1 {
        input:
            bam = whole_hap1_bam,
            bai = whole_hap1_bai,
            locus = "chrM",
            prefix = sampleid + "hap1"
    }

    call SubsetBam as Hap2 {
        input:
            bam = whole_hap2_bam,
            bai = whole_hap2_bai,
            locus = "chrM",
            prefix = sampleid + "hap2"
    }

    call MergeFasta {
        input:        
            fasta1 = Hap1.subset_fasta,
            fasta2 = Hap2.subset_fasta,
            locus = "chrM",
            prefix = sampleid
    }
    
    call SubsetBam as SubsetReads {
            input:
                bam = whole_read_bam,
                bai = whole_read_bai,
                locus = "chrM",
                prefix = sampleid + "reads"
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
        call QuickStart {
            input:
                bam = downsampleBam.downsampled_bam,
                bai = downsampleBam.downsampled_bai,
                reference_fa = reference_fa,
                prefix = prefix + "_" + desiredCoverage + "x",
                kmer_size = kmer_size,
                sample_id = sampleid,
                chromo = "chrM",
                data_type = "ont"
        }


    }

    output {
        File truth_fasta = MergeFasta.merged_fasta
        Array[File] Himito_fasta = QuickStart.asm
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

task MergeFasta {

    meta {
        description : "cat two fasta file into one"
    }

    parameter_meta {
    }

    input {
        File fasta1
        File fasta2
        String locus
        String prefix
    }

    command <<<
        cat ~{fasta1} ~{fasta2} > ~{prefix}~{locus}.fasta
    >>>

    output {
        File merged_fasta = "~{prefix}~{locus}.fasta"
        
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


task QuickStart {
    input {
        File bam
        File bai
        File reference_fa
        String prefix
        Int kmer_size
        String sample_id
        String chromo
        String data_type
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
