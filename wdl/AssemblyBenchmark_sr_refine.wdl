version 1.0
workflow CompareFasta {
    input {
        File whole_hap1_bam
        File whole_hap1_bai
        File whole_hap2_bam
        File whole_hap2_bai
        File whole_read_bam
        File whole_read_bai
        File whole_sr_bam
        File whole_sr_bai

        File reference_fa
        File reference_gb
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
    
    # subset sr data
    call SubsetBam as SubsetReads_sr {
                input:
                    bam = whole_sr_bam,
                    bai = whole_sr_bai,
                    locus = "chrM",
                    prefix = sampleid + "sr"
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

        call Correct {
            input:
                graph_gfa = Build.graph,
                short_read_fasta = SubsetReads_sr.subset_fasta,
                reference_file = reference_fa,
                prefix = "sr",
                sampleid = sampleid,
                min_support_counts = 10,
                query_length = 99
            
        }


        call Asm as ASM_raw {
            input:
                graph_gfa = Build.graph,
                prefix = prefix + "_Himito_raw_${desiredCoverage}",
                sampleid=sampleid
        }

        call Asm as ASM_corrected {
            input:
                graph_gfa = Correct.corrected_graph,
                prefix = prefix + "_Himito_corrected_${desiredCoverage}",
                sampleid=sampleid
        }

        # call MitoHiFi assembly
        call SubsetBam as SubsetReads_1 {
                input:
                    bam = downsampleBam.downsampled_bam,
                    bai = downsampleBam.downsampled_bai,
                    locus = "chrM",
                    prefix = sampleid + "reads"
        }
        call MitoHifiAsm {
            input:
                reads = SubsetReads_1.subset_fastq,
                reffa = reference_fa,
                refgb = reference_gb,
                prefix = prefix + "_mitohifi${desiredCoverage}",

        }
    }
    
    output {
        File truth_fasta = MergeFasta.merged_fasta
        Array[File] Himito_corrected_fasta = ASM_corrected.fasta
        Array[File] Himito_raw_fasta = ASM_raw.fasta
        Array[File] mitohifi_fasta = MitoHifiAsm.final_fa
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
        docker: "hangsuunc/himito:v2"
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
        docker: "hangsuunc/himito:v2"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task Correct {
    input {
        File graph_gfa
        File short_read_fasta
        File reference_file
        String prefix
        String sampleid
        Int min_support_counts
        Int query_length
    }

    command <<<
        set -euxo pipefail
        msbwt2-build -o sr_msbwt.npy ~{short_read_fasta}
        /Himito/target/release/Himito correct -g ~{graph_gfa} -b sr_msbwt.npy -r ~{reference_file} -o ~{sampleid}.~{prefix}.corrected.gfa -m ~{min_support_counts} -q ~{query_length}
    >>>

    output {
        File corrected_graph = "~{sampleid}.~{prefix}.corrected.gfa"
    }

    runtime {
        docker: "hangsuunc/himito:v2"
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
        docker: "hangsuunc/himito:v2"
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
        docker: "hangsuunc/himito:v2"
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
        docker: "hangsuunc/himito:v2"
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

task MitoHifiAsm{
    input{
        File reads
        File reffa
        File refgb
        String prefix
        Int num_cpus
    }

    #Int disk_size = 50 + ceil(2 * size(reads, "GiB"))

    command <<<
        #set -euxo pipefail
        mitohifi.py -r ~{reads} -f ~{reffa} -g ~{refgb} -t ~{num_cpus} -o 1 
        ls -l
        mv final_mitogenome.fasta ~{prefix}.fasta
        mv contigs_stats.tsv ~{prefix}_contig_stats.tsv
    >>>

    output{
        File final_fa="~{prefix}.fasta"
        File final_stats="~{prefix}_contig_stats.tsv"
        File mapping_file = "coverage_mapping/HiFi-vs-final_mitogenome.sorted.bam"
        

    }
    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk " + "256" + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 64
        preemptible: 2
        maxRetries: 1
        docker: "ghcr.io/marcelauliano/mitohifi:master"
        #docker:"biocontainers/mitohifi:3.0.0_cv1"
    }    
}