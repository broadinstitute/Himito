# Installation in restricted network environments

Some users cannot reliably reach GitHub, `crates.io`, Docker Hub, or the default Rust installer CDN. This page documents alternative ways to install **Himito** without depending on those services.

For the default install path, see the [README](../README.md).

## Recommended: Zenodo source archive

Each release is archived on Zenodo with a DOI. Download the zip on any network that can reach Zenodo (often more reliable than GitHub in some regions), then build locally:

```bash
# Latest archived release (check Zenodo for newer versions)
wget https://zenodo.org/records/21522817/files/broadinstitute/Himito-v1.1.2-alpha.zip
unzip broadinstitute/Himito-v1.1.2-alpha.zip
cd Himito-v1.1.2-alpha

# Install system build dependencies (Debian/Ubuntu)
sudo apt-get update
sudo apt-get install -y build-essential cmake pkg-config libssl-dev libclang-dev zlib1g-dev

cargo build --release
./target/release/Himito --help
```

Archive landing page: [https://doi.org/10.5281/zenodo.21522816](https://doi.org/10.5281/zenodo.21522816)

Verify the download (optional):

```bash
# md5 from Zenodo metadata for v1.1.2-alpha
md5sum broadinstitute/Himito-v1.1.2-alpha.zip
# expected: 33e2d3e2a617cd409216563d3c98d420
```

## Prebuilt binary (no compile step)

GitHub Releases include a Linux x86_64 binary. If GitHub is blocked, copy the file from a machine with access, or wait for a binary to be attached to the Zenodo record.

```bash
wget https://github.com/broadinstitute/Himito/releases/download/v1.1.2/Himito
chmod +x Himito
./Himito --help
```

Requires a compatible glibc/Linux x86_64 environment. For other platforms, build from source.

## Rust toolchain without rustup.rs

If `https://sh.rustup.rs` is unreachable, use a mirror or an offline rustup bundle.

**USTC mirror (common in China):**

```bash
export RUSTUP_DIST_SERVER=https://mirrors.ustc.edu.cn/rust-static
export RUSTUP_UPDATE_ROOT=https://mirrors.ustc.edu.cn/rust-static/rustup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
```

**Tsinghua mirror:**

```bash
export RUSTUP_DIST_SERVER=https://mirrors.tuna.tsinghua.edu.cn/rustup
export RUSTUP_UPDATE_ROOT=https://mirrors.tuna.tsinghua.edu.cn/rustup/rustup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
```

Or install Rust from your Linux distribution and ensure `cargo` is ≥ 1.70.

## Cargo / crates.io mirror

Himito pulls Rust crates (including `rust-htslib`, which compiles against system libraries). Point Cargo at a mirror before `cargo build`:

```bash
mkdir -p ~/.cargo
cat >> ~/.cargo/config.toml <<'EOF'
[source.crates-io]
replace-with = "ustc"

[source.ustc]
registry = "sparse+https://mirrors.ustc.edu.cn/crates.io-index/"

[registries.ustc]
index = "sparse+https://mirrors.ustc.edu.cn/crates.io-index/"
EOF
```

Tsinghua alternative: `sparse+https://mirrors.tuna.tsinghua.edu.cn/crates.io-index/`

## Offline / air-gapped build

On a networked machine:

```bash
git clone https://github.com/broadinstitute/Himito.git
cd Himito
cargo fetch
tar czf himito-vendor.tar.gz ~/.cargo/registry ~/.cargo/git target 2>/dev/null || \
  tar czf himito-cargo-registry.tar.gz ~/.cargo/registry ~/.cargo/git
# Also archive the source tree (or use the Zenodo zip)
```

Transfer `Himito/`, `~/.cargo/registry`, and `~/.cargo/git` to the offline host, then:

```bash
cd Himito
cargo build --release --offline
```

If `--offline` fails because a crate is missing, re-run `cargo fetch` on the connected machine with the same `Cargo.lock`.

## Docker without Docker Hub

The [Dockerfile](../docker/Dockerfile) uses `rust:1.84.0` from Docker Hub. Options:

1. **Load a saved image** (prepare on a machine with Docker Hub access):

   ```bash
   docker pull rust:1.84.0
   docker build -f docker/Dockerfile -t himito:latest .
   docker save himito:latest | gzip > himito-docker.tar.gz
   ```

   On the restricted host:

   ```bash
   docker load -i himito-docker.tar.gz
   docker run --rm himito:latest Himito --help
   ```

2. **Use a regional mirror** of Docker Hub (e.g. Aliyun, DaoCloud) and retag `rust:1.84.0` before building.

## Optional: `msbwt2` (short-read graph correction)

The `correct` command uses [msBWT](https://github.com/HudsonAlpha/rust-msbwt). Install with the same Cargo mirror settings:

```bash
cargo install msbwt2
msbwt2-build --help
```

For offline use, run `cargo install msbwt2` on a connected machine and copy `~/.cargo/bin/msbwt2*` to the target host.

## WDL / Terra / Cromwell

Himito WDLs under [`wdl/`](../wdl/) are typically run with Docker images from Google Container Registry. If GCR is unreachable:

- Run Cromwell locally with a locally built Himito binary instead of the cloud Docker task command, or
- Mirror the required image to a regional registry and update the `docker` runtime attribute in your inputs JSON.

## Verify installation

```bash
./target/release/Himito --version 2>/dev/null || ./Himito --help
./target/release/Himito --help
cargo test   # optional, from source checkout
```

## Citation

If you use Himito, please cite the software DOI and the preprint:

- Software: [https://doi.org/10.5281/zenodo.21522816](https://doi.org/10.5281/zenodo.21522816)
- Paper: [https://doi.org/10.1101/2025.11.03.686348](https://doi.org/10.1101/2025.11.03.686348)
