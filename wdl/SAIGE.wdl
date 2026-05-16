version 1.0

workflow SAIGE_all{
    meta{
        description: "a workflow for variant association using SAIGE_singlevariant and SAIGE_geneset"
    }
    input{
        File plink_bed_file
        File plink_bim_file
        File plink_fam_file
        File sparseGRM
        File sparseGRM_IDlist
        File phenotype_file
        File homoplasmic_vcf
        File homoplasmic_vcf_csi
        File heteroplasmic_vcf
        File heteroplasmic_vcf_csi
        File GroupFile
        String chromosome
        Array[String] phecode_list
        String trait_type
        String output_prefix
        String saige_docker
        String vcffield
        String memory = "8G"
        Float single_variant_minimal_af = 0.01
        Int single_variant_min_mac = 20
        Float gene_set_minimal_af = 0
        Float gene_set_min_mac = 0.5
        Int cpu = 2
        String disk_size = "local-disk 50 HDD"
        Int preemptible = 1

    }

    scatter (phecode in phecode_list) {
        call RunFitNullGLMM {input:
            plink_bed_file = plink_bed_file,
            plink_bim_file = plink_bim_file,
            plink_fam_file = plink_fam_file,
            sparseGRM = sparseGRM,
            sparseGRM_IDlist = sparseGRM_IDlist,
            phenotype_file = phenotype_file,
            phecode = phecode,
            trait_type =trait_type,
            output_prefix = output_prefix,
            saige_docker = saige_docker,
            memory = memory,
            cpu = cpu,
            disk = disk_size,
            preemptible = preemptible
        }
        
        call RunStep2_singlevariant { input:
            vcf = homoplasmic_vcf,
            vcf_csi = homoplasmic_vcf_csi,
            variance_ratio = RunFitNullGLMM.variance_ratio_txt,
            GMMATmodelFile = RunFitNullGLMM.null_model_rda,
            vcffield = vcffield,
            chromo = chromosome,
            phecode = phecode,
            trait_type = trait_type,
            output_prefix = output_prefix,
            minimal_af = single_variant_minimal_af,
            min_mac = single_variant_min_mac,
            memory = memory,
            saige_docker = saige_docker,
            cpu = cpu,
            disk = disk_size,
            preemptible = preemptible
        }

        call RunStep2_geneset { input:
            vcf = heteroplasmic_vcf,
            vcf_csi = heteroplasmic_vcf_csi,
            variance_ratio = RunFitNullGLMM.variance_ratio_txt,
            GMMATmodelFile = RunFitNullGLMM.null_model_rda,
            vcffield = vcffield,
            chromo = chromosome,
            phecode = phecode,
            output_prefix = output_prefix,
            minimal_af = gene_set_minimal_af,
            min_mac = gene_set_min_mac,
            GroupFile = GroupFile,
            memory = memory,
            saige_docker = saige_docker,
            cpu = cpu,
            disk = disk_size,
            preemptible = preemptible
        }

    }

    

    output {
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


task RunFitNullGLMM {

  input {
    File plink_bed_file
    File plink_bim_file
    File plink_fam_file
    File sparseGRM
    File sparseGRM_IDlist
    File phenotype_file
    String phecode
    String trait_type
    String output_prefix
    String saige_docker

    # Runtime parameters
    String memory
    Int cpu
    String disk
    Int preemptible
  }
  command <<<
    set -euxo pipefail   

    step1_fitNULLGLMM.R \
        --bedFile="~{plink_bed_file}"  \
        --bimFile="~{plink_bim_file}"  \
        --famFile="~{plink_fam_file}"  \
        --useSparseGRMtoFitNULL=TRUE    \
        --sparseGRMFile="~{sparseGRM}" \
        --sparseGRMSampleIDFile="~{sparseGRM_IDlist}"  \
        --phenoFile="~{phenotype_file}" \
        --phenoCol="~{phecode}" \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,mPC1,mPC2,mPC3,mPC4,mPC5,mPC6,mPC7,mPC8,mPC9,mPC10,sex,age,age2,age_sex,age2_sex \
        --qCovarColList=sex \
        --sampleIDColinphenoFile=person_id \
        --traitType=~{trait_type} \
        --isCateVarianceRatio=TRUE      \
        --cateVarRatioMinMACVecExclude=200,500 \
        --cateVarRatioMaxMACVecInclude=500,1000000  \
        --IsOverwriteVarianceRatioFile=TRUE \
        --outputPrefix="~{output_prefix}_step1Out_~{phecode}"
    >>>

  output {
    File null_model_rda = "~{output_prefix}_step1Out_~{phecode}.rda"
    File variance_ratio_txt = "~{output_prefix}_step1Out_~{phecode}.varianceRatio.txt"
  }

  runtime {
    docker: saige_docker
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }

}

task RunStep2_singlevariant {

  input {
    File vcf
    File vcf_csi
    File variance_ratio
    File GMMATmodelFile
    String vcffield
    String chromo
    String phecode
    String trait_type = "binary"
    String output_prefix
    Float minimal_af
    Int min_mac


    # Runtime parameters
    String memory
    String saige_docker
    Int cpu
    String disk
    Int preemptible
  }
  command <<<
    set -euxo pipefail   

    step2_SPAtests.R \
      --vcfFile=~{vcf} \
      --vcfFileIndex=~{vcf_csi} \
      --vcfField=~{vcffield} \
      --AlleleOrder=ref-first \
      --chrom=~{chromo} \
      --SAIGEOutputFile=~{output_prefix}_step2Out_~{phecode}_singlevar \
      --minMAF=~{minimal_af} \
      --minMAC=~{min_mac} \
      --GMMATmodelFile=~{GMMATmodelFile} \
      --varianceRatioFile=~{variance_ratio} \
      --is_Firth_beta=TRUE \
      --LOCO=FALSE \
      --is_output_moreDetails=TRUE

    >>>

  output {
    Array[File] outputfiles = glob("~{output_prefix}_step2Out_~{phecode}*")
  }

  runtime {
    docker: saige_docker
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }

}

task RunStep2_geneset {
   meta {
        desciption:
        "heteroplasmic variant association (all heteroplasmies) analysis using gene set test"
    }

  input {
    File vcf
    File vcf_csi
    File variance_ratio
    File GMMATmodelFile
    String vcffield
    String chromo
    String phecode
    String output_prefix
    Float minimal_af
    Float min_mac
    File GroupFile


    # Runtime parameters
    String memory
    String saige_docker
    Int cpu
    String disk
    Int preemptible
  }
  command <<<
    set -euxo pipefail   
    
    step2_SPAtests.R \
      --vcfFile=~{vcf} \
      --vcfFileIndex=~{vcf_csi} \
      --vcfField=~{vcffield} \
      --AlleleOrder=ref-first \
      --chrom=~{chromo} \
      --SAIGEOutputFile=~{output_prefix}_step2Out_~{phecode}_geneset \
      --minMAF=~{minimal_af} \
      --minMAC=~{min_mac} \
      --GMMATmodelFile=~{GMMATmodelFile} \
      --varianceRatioFile=~{variance_ratio} \
      --is_Firth_beta=TRUE \
      --LOCO=FALSE \
      --groupFile="~{GroupFile}"    \
      --annotation_in_groupTest="non_coding_transcript_exon,start_lost,stop_gained;start_lost,missense,synonymous,frameshift,stop_gained,stop_lost,stop_retained,dloop"        \
      --maxMAF_in_groupTest=0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5 \
      --pCutoffforFirth=0.05 \
      --is_output_markerList_in_groupTest=TRUE \
      --is_output_moreDetails=TRUE

    >>>

  output {
    Array[File] outputfiles = glob("~{output_prefix}_step2Out_~{phecode}*")
  }

  runtime {
    docker: saige_docker
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }

}