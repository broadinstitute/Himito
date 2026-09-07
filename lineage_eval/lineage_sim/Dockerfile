# Docker image for the Himito lineage simulation sweep (Himito_lineage_sim.wdl).
#
# Bakes in everything run_eval.sh -> {simulate_tree.py, simulate_reads.sh,
# run_himito.sh, score_lineage.py} needs:
#   - the Himito release binary
#   - the read-simulation toolchain (pbsim3, pbccs, minimap2, samtools)
#   - Python + numpy/dendropy
#   - the pbsim3 error models
#   - the lineage_sim scripts and rCRS.fasta
#
# pbccs is linux/amd64-only, so the runtime is always amd64. Do NOT compile
# Rust under QEMU on Apple Silicon: `rustc -vV` SIGSEGVs there. Default is to
# copy a prebuilt linux/amd64 Himito binary. Compile from source only on a
# real amd64 host (--build-arg HIMITO_SOURCE=from-source).
#
# Build context is this directory (lineage_sim/):
#   docker build \
#     -t us.gcr.io/broad-dsp-lrma/hangsuunc/himito-lineage-sim:dev \
#     -f Dockerfile .
# Override the prebuilt binary source if needed:
#   --build-arg HIMITO_IMAGE=us.gcr.io/broad-dsp-lrma/hangsuunc/himito:dev

# ---- optional: compile Himito (linux/amd64 host only; skipped by default) ----
FROM rust:1.84.0 AS himito-from-source
ARG DEBIAN_FRONTEND=noninteractive
ARG HIMITO_REPO=https://github.com/broadinstitute/Himito.git
ARG HIMITO_REF=hs_lineage_
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git ssh build-essential make cmake zlib1g-dev libclang-dev wget util-linux bsdextrautils \
    && rm -rf /var/lib/apt/lists/*
RUN git config --global http.sslverify false \
    && git clone --depth 1 --branch "${HIMITO_REF}" "${HIMITO_REPO}" /Himito \
    && cd /Himito \
    && cargo build --release

# ---- default: reuse a published linux/amd64 Himito (no rustc, QEMU-safe) -----
ARG HIMITO_IMAGE=docker.io/hangsuunc/himito:dev
FROM --platform=linux/amd64 ${HIMITO_IMAGE} AS himito-prebuilt

ARG HIMITO_SOURCE=prebuilt
FROM himito-${HIMITO_SOURCE} AS himito-build

# ---- runtime ----------------------------------------------------------------
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8 AS runtime
USER root
ENV MAMBA_ROOT_PREFIX=/opt/conda

# Simulation/reconstruction toolchain (pbsim3, pbccs, minimap2, samtools, python).
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH

# pbsim3 ships the binary via bioconda but not always its error models; fetch the
# data/ dir the way setup_env.sh does.
RUN git clone --depth 1 https://github.com/yukiteruono/pbsim3.git /tmp/pbsim3_src \
    && mkdir -p /opt/lineage_sim/pbsim3_models \
    && cp /tmp/pbsim3_src/data/*.model /opt/lineage_sim/pbsim3_models/ \
    && rm -rf /tmp/pbsim3_src

# Himito binary + rCRS.fasta from the build stage.
COPY --from=himito-build /Himito/target/release/Himito /Himito/target/release/Himito
COPY rCRS.fasta /opt/lineage_sim/rCRS.fasta

# The lineage_sim scripts (build context = this directory).
COPY run_eval.sh run_himito.sh simulate_reads.sh simulate_tree.py score_lineage.py \
     /opt/lineage_sim/
RUN chmod +x /opt/lineage_sim/*.sh

ENV HIMITO=/Himito/target/release/Himito \
    REF=/opt/lineage_sim/rCRS.fasta \
    PBSIM_MODEL_DIR=/opt/lineage_sim/pbsim3_models

WORKDIR /work
