#!/usr/bin/env bash
# Create the himito-eval conda env, fetch pbsim3 error models, and verify tools.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 1. Create env (idempotent).
if ! conda env list | grep -q '^himito-eval '; then
  conda env create -f "$HERE/environment.yml"
else
  echo "env himito-eval already exists; skipping create"
fi

# 2. pbsim3 bioconda package ships the binary but not always the error models.
#    Clone the repo just for its data/ directory.
MODELS_DIR="$HERE/pbsim3_models"
if [[ ! -f "$MODELS_DIR/ERRHMM-SEQUEL.model" ]]; then
  rm -rf "$HERE/.pbsim3_src"
  git clone --depth 1 https://github.com/yukiteruono/pbsim3.git "$HERE/.pbsim3_src"
  mkdir -p "$MODELS_DIR"
  cp "$HERE/.pbsim3_src/data/"*.model "$MODELS_DIR/"
  rm -rf "$HERE/.pbsim3_src"
fi
echo "PBSIM_MODEL_DIR=$MODELS_DIR"

# 3. Verify tools resolve inside the env.
conda run -n himito-eval bash -c '
  set -e
  for t in pbsim ccs minimap2 samtools python; do
    command -v "$t" >/dev/null || { echo "MISSING: $t" >&2; exit 1; }
  done
  python -c "import numpy, dendropy"
  echo "all tools present"
'
ls "$MODELS_DIR"/ERRHMM-SEQUEL.model "$MODELS_DIR"/ERRHMM-ONT-HQ.model
