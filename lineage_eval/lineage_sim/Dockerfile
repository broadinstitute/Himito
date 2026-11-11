# Linux amd64 image for Himito_lineage_sim.wdl, including HiFi/ccs.
# Build from the repository root using the current checkout:
#   docker buildx build --platform linux/amd64 --provenance=false --sbom=false \
#     --load -t himito-lineage-sim:linux-amd64 \
#     -f lineage_eval/lineage_sim/Dockerfile .
# Default HIMITO_SOURCE=workspace compiles Cargo.lock and src from this checkout.
# Optional modes: --build-arg HIMITO_SOURCE=prebuilt --build-arg HIMITO_IMAGE=...
# or --build-arg HIMITO_SOURCE=from-source --build-arg HIMITO_REF=...

ARG HIMITO_IMAGE=us.gcr.io/broad-dsp-lrma/hangsuunc/himito:dev
ARG HIMITO_SOURCE=workspace

# ---- default: compile the current workspace with a locked dependency set ----
FROM rust:1.90.0-bookworm AS himito-workspace
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential cmake pkg-config zlib1g-dev libbz2-dev liblzma-dev libclang-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /Himito
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release --locked
COPY test_data/rCRS.fasta ./test_data/rCRS.fasta

# Optional remote-source and prebuilt modes remain available explicitly.
FROM rust:1.90.0-bookworm AS himito-from-source
ARG HIMITO_REPO=https://github.com/broadinstitute/mitograph.git
ARG HIMITO_REF=lineage_bench
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git build-essential cmake pkg-config zlib1g-dev libbz2-dev liblzma-dev libclang-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 --branch "${HIMITO_REF}" "${HIMITO_REPO}" /Himito \
    && cd /Himito && cargo build --release --locked

# ---- optional: reuse a published linux/amd64 Himito -----
FROM --platform=linux/amd64 ${HIMITO_IMAGE} AS himito-prebuilt

FROM himito-${HIMITO_SOURCE} AS himito-build

# ---- runtime ----------------------------------------------------------------
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8 AS runtime
USER root
ENV MAMBA_ROOT_PREFIX=/opt/conda

# Simulation/reconstruction toolchain (pbsim3, pbccs, minimap2, samtools, python, pysam).
RUN apt-get update \
    && apt-get install -y --no-install-recommends bsdextrautils ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY lineage_eval/lineage_sim/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH

# pbsim3 ships the binary via bioconda but not always its error models; fetch the
# data/ dir the way setup_env.sh does.
RUN git clone --depth 1 https://github.com/yukiteruono/pbsim3.git /tmp/pbsim3_src \
    && mkdir -p /opt/lineage_sim/pbsim3_models \
    && cp /tmp/pbsim3_src/data/*.model /opt/lineage_sim/pbsim3_models/ \
    && rm -rf /tmp/pbsim3_src

# Binary and reference from the selected workspace/remote/prebuilt stage.
COPY --from=himito-build /Himito/target/release/Himito /Himito/target/release/Himito
COPY --from=himito-build /Himito/test_data/rCRS.fasta /opt/lineage_sim/rCRS.fasta

# The lineage_sim scripts (build context = repository root).
COPY lineage_eval/lineage_sim/*.sh lineage_eval/lineage_sim/*.py /opt/lineage_sim/
RUN chmod +x /opt/lineage_sim/*.sh

ENV HIMITO=/Himito/target/release/Himito \
    REF=/opt/lineage_sim/rCRS.fasta \
    PBSIM_MODEL_DIR=/opt/lineage_sim/pbsim3_models

WORKDIR /work
