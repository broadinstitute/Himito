version 1.0

workflow MergeBam{
    meta{
        description: "a workflow for MHC annotation and gene extraction"
    }
    input{
        Array[File] bams
        String prefix
        String destination_path
        Int batch_size = 100    
    }

    call CreateBatches {
        input:
            bams = bams,
            batch_size = batch_size
    }

    scatter (index in range(length(CreateBatches.bams_batch))) {
        Array[File] splited_bams = read_lines(CreateBatches.bams_batch[index])
        
        call MergeBamFile as MergeBamFile_batch {
            input:
                bams = splited_bams,
                prefix = prefix + "_batch_" + index
        }
    }

    call MergeBamFile as MergeBamFile_Final {
        input:
            bams = select_all(MergeBamFile_batch.merged_bam),
            prefix = prefix
    }

    call FinalizeToFile {
        input: 
            file=MergeBamFile_Final.merged_bam,
            outdir=destination_path
    }


    output {
        File merged_bam = MergeBamFile_Final.merged_bam

    }
}

task CreateBatches {
    input {
        Array[String] bams
        Int batch_size
    }

    command {
        set -euox pipefail

        cat ~{write_lines(bams)} | split -l ~{batch_size} - bams_batch_
        
    }

    output {
        Array[File] bams_batch = glob("bams_batch_*")
    }

    Int disk_size = 100

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsde-methods/slee/kage-lite:pr_29"
    }
}

task MergeBamFile{
    input{
        Array[File] bams
        String prefix
    }
    command <<<
        set -euxo pipefail

        for bam in ~{sep=" " bams}; do
            samtools index "$bam"
        done
       samtools merge -o ~{prefix}.bam ~{sep=" " bams}
    >>>
    
    output{
        File merged_bam="~{prefix}.bam"
        
    }

    Int disk_size = 100

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task FinalizeToFile {

    meta{
        description: "Copies the given file to the specified bucket."
    }

    parameter_meta {
        file: {
            description: "file to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
        name:   "name to set for uploaded file"
    }

    input {
        File file
        String outdir
        String? name

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }



    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euxo pipefail

        gsutil -m cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
