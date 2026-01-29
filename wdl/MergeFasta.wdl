version 1.0

workflow MergeFasta{
    meta{
        description: "a workflow for MHC annotation and gene extraction"
    }
    input{
        Array[File] fastas
        String prefix
        String destination_path
        Int batch_size = 100    
    }

    call CreateBatches {
        input:
            fastas = fastas,
            batch_size = batch_size
    }

    scatter (index in range(length(CreateBatches.fasta_batch))) {
        Array[File] splited_fastas = read_lines(CreateBatches.fasta_batch[index])
        
        call MergeFastaFile as MergeFastaFile_batch {
            input:
                fastas = splited_fastas,
                prefix = prefix + "_batch_" + index
        }
    }

    call MergeFastaFile as MergeFastaFile_Final {
        input:
            fastas = select_all(MergeFastaFile_batch.merged_fasta),
            prefix = prefix
    }

    call FinalizeToFile {
        input: 
            file=MergeFastaFile_Final.merged_fasta,
            outdir=destination_path
    }


    output {
        File merged_fasta = MergeFastaFile_Final.merged_fasta

    }
}

task CreateBatches {
    input {
        Array[String] fastas
        Int batch_size
    }

    command {
        set -euox pipefail

        cat ~{write_lines(fastas)} | split -l ~{batch_size} - fastas_batch_
        
    }

    output {
        Array[File] fasta_batch = glob("fastas_batch_*")
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

task MergeFastaFile{
    input{
        Array[File] fastas
        String prefix
    }
    command <<<
        set -euxo pipefail

        # if fasta files are gzipped, unzip them first
        for fasta in ~{sep=" " fastas}; do
            if [[ "$fasta" == *.gz ]]; then
                gunzip -c "$fasta"
            else
                cat "$fasta"
            fi
        done > ~{prefix}.fasta

        bgzip ~{prefix}.fasta
    >>>
    
    output{
        File merged_fasta="~{prefix}.fasta.gz"
        
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
