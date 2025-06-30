FROM rust:1.84.0
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y --no-install-recommends git ssh build-essential make cmake zlib1g-dev vim wget unzip libclang-dev\
    && rm -rf /var/lib/apt/lists/*

RUN git config --global http.sslverify false
RUN git clone https://github.com/broadinstitute/Himito.git \ 
    && cd Himito \
    && cargo build --release

RUN cargo install msbwt2 \
    && msbwt2-build --help