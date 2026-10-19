version 1.0

# Epigenetic age from PacBio pb-CpG-tools pileups (combined + haplotype-resolved) with pyaging.
#
#   per sample : extract_clock_cpgs   -> Illumina cg-probe CpGs (GRCh38) with beta + coverage per track
#   cohort     : predict_epigenetic_age -> pyaging clocks on combined / hap1 / hap2 rows + fill QC
#
# Build the image first (see README): docker build -f docker/Dockerfile -t epiage-pyaging:0.5.7 .

workflow pyaging_epigenetic_age {
  meta {
    description: "Estimate DNA methylation age from long-read CpG pileups (pb-CpG-tools or modkit) using pyaging clocks, per sample and per haplotype."
    allowNestedInputs: true
  }

  parameter_meta {
    sample_ids: "Sample identifiers; used as output row labels and file prefixes (must be unique)."
    combined_beds: "pb-CpG-tools <prefix>.combined.bed(.gz), one per sample, same order as sample_ids. Must be on GRCh38 unless probe_manifest is replaced."
    hap1_beds: "Optional pb-CpG-tools <prefix>.hap1.bed(.gz), one per sample (all-or-none)."
    hap2_beds: "Optional pb-CpG-tools <prefix>.hap2.bed(.gz), one per sample (all-or-none)."
    covariates_tsv: "Optional TSV with columns sample_id, age, sex (M/F). Required for clocks that take age/sex as inputs (pcgrimage, grimage2); those are skipped otherwise."
    pileup_format: "'pb-cpg-tools' (default) or 'modkit' (modkit pileup --cpg --combine-strands bedMethyl)."
    probe_manifest: "Optional override of the cg-probe coordinate table (chrom, start = 0-based C, probe_id). Default: GRCh38 table baked into the image."
    min_coverage: "Minimum read coverage for a CpG to count as observed on a track."
    min_probe_fraction: "Fail a sample whose combined track matches fewer than this fraction of manifest probes (build/coordinate mismatch guard)."
    clocks: "Comma-separated pyaging clock names. Clocks not prefetched into the image need offline_weights = false and network access."
    tracks: "Which rows to score: any of combined,hap1,hap2."
    hap_fill: "'combined' fills CpGs not covered on a haplotype with that sample's combined value; 'reference' uses the cohort/reference chain instead."
    cohort_min_samples: "Minimum samples with an observed combined value before a missing CpG is filled with the cohort median (0 disables)."
    offline_weights: "Use only clock weights cached in the image (HF_HUB_OFFLINE=1)."
  }

  input {
    Array[String] sample_ids
    Array[File] combined_beds
    Array[File]? hap1_beds
    Array[File]? hap2_beds
    File? covariates_tsv

    String pileup_format = "pb-cpg-tools"
    File? probe_manifest
    Int min_coverage = 4
    Float min_probe_fraction = 0.5

    String clocks = "horvath2013,skinandblood,hannum,dnamphenoage,pchorvath2013,pcskinandblood,pcphenoage,pcgrimage,dunedinpace"
    String tracks = "combined,hap1,hap2"
    String hap_fill = "combined"
    Int cohort_min_samples = 3
    Boolean offline_weights = true

    String output_prefix = "cohort"
    String docker = "epiage-pyaging:0.5.7"
    Int predict_mem_gb = 32
  }

  scatter (idx in range(length(sample_ids))) {
    if (defined(hap1_beds)) {
      File hap1_bed_i = select_first([hap1_beds])[idx]
    }
    if (defined(hap2_beds)) {
      File hap2_bed_i = select_first([hap2_beds])[idx]
    }

    call extract_clock_cpgs {
      input:
        sample_id = sample_ids[idx],
        combined_bed = combined_beds[idx],
        hap1_bed = hap1_bed_i,
        hap2_bed = hap2_bed_i,
        pileup_format = pileup_format,
        probe_manifest = probe_manifest,
        min_coverage = min_coverage,
        min_probe_fraction = min_probe_fraction,
        docker = docker
    }
  }

  call predict_epigenetic_age {
    input:
      tables = extract_clock_cpgs.clock_cpgs,
      covariates_tsv = covariates_tsv,
      clocks = clocks,
      tracks = tracks,
      hap_fill = hap_fill,
      cohort_min_samples = cohort_min_samples,
      offline_weights = offline_weights,
      output_prefix = output_prefix,
      docker = docker,
      mem_gb = predict_mem_gb
  }

  output {
    File epigenetic_age = predict_epigenetic_age.ages
    File epigenetic_age_qc = predict_epigenetic_age.qc
    File clock_status = predict_epigenetic_age.clock_status
    Array[File] clock_cpg_tables = extract_clock_cpgs.clock_cpgs
    Array[File] extract_qc = extract_clock_cpgs.qc
  }
}

task extract_clock_cpgs {
  input {
    String sample_id
    File combined_bed
    File? hap1_bed
    File? hap2_bed
    String pileup_format
    File? probe_manifest
    Int min_coverage
    Float min_probe_fraction
    String docker
  }

  Int disk_gb = ceil(size(combined_bed, "GB") + size(hap1_bed, "GB") + size(hap2_bed, "GB")) + 20

  command <<<
    set -euo pipefail
    python3 /opt/epiage/scripts/extract_clock_cpgs.py \
      --sample-id "~{sample_id}" \
      --combined "~{combined_bed}" \
      ~{"--hap1 " + hap1_bed} \
      ~{"--hap2 " + hap2_bed} \
      --format "~{pileup_format}" \
      ~{"--manifest " + probe_manifest} \
      --min-coverage ~{min_coverage} \
      --min-probe-fraction ~{min_probe_fraction} \
      --out-prefix "~{sample_id}"
  >>>

  output {
    File clock_cpgs = "~{sample_id}.clock_cpgs.tsv.gz"
    File qc = "~{sample_id}.extract_qc.json"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disk: disk_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 2
    maxRetries: 1
  }
}

task predict_epigenetic_age {
  input {
    Array[File] tables
    File? covariates_tsv
    String clocks
    String tracks
    String hap_fill
    Int cohort_min_samples
    Boolean offline_weights
    String output_prefix
    String docker
    Int mem_gb = 32
    Int threads = 4
  }

  Int disk_gb = ceil(size(tables, "GB")) + 40

  command <<<
    set -euo pipefail
    ~{if offline_weights then "export HF_HUB_OFFLINE=1" else ""}
    export OMP_NUM_THREADS=~{threads} MKL_NUM_THREADS=~{threads}
    python3 /opt/epiage/scripts/predict_epigenetic_age.py \
      --tables ~{sep=" " tables} \
      ~{"--covariates " + covariates_tsv} \
      --clocks "~{clocks}" \
      --tracks "~{tracks}" \
      --hap-fill "~{hap_fill}" \
      --cohort-min-samples ~{cohort_min_samples} \
      --out-prefix "~{output_prefix}"
  >>>

  output {
    File ages = "~{output_prefix}.epigenetic_age.tsv"
    File qc = "~{output_prefix}.epigenetic_age_qc.tsv"
    File clock_status = "~{output_prefix}.clock_status.tsv"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}
