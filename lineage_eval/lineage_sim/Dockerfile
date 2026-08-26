# Docker image for the Himito lineage simulation sweep (Himito_lineage_sim.wdl).
#
# Bakes in everything run_eval.sh -> {simulate_tree.py, simulate_reads.sh,
# run_himito.sh, score_lineage.py} needs:
#   - the Himito release binary (built from source)
#   - the read-simulation toolchain (pbsim3, pbccs, minimap2, samtools)
#   - Python + numpy/dendropy
#   - the pbsim3 error models
#   - the lineage_sim scripts and rCRS.fasta
#
# ccs (pbccs) has no macOS-ARM64 build, so the HiFi profile needs linux/amd64.
# The stage below is pinned to linux/amd64 so `docker build` works on Apple Silicon.
#
# Build context is this directory (lineage_sim/):
#   docker build --platform=linux/amd64 \
#     -t us.gcr.io/broad-dsp-lrma/hangsuunc/himito-lineage-sim:dev \
#     -f Dockerfile .

# ---- Stage 1: build the Himito release binary --------------------------------
FROM --platform=linux/amd64 rust:1.84.0 AS himito-build
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

# ---- Stage 2: runtime image --------------------------------------------------
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
