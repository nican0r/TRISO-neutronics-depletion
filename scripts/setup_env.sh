#!/usr/bin/env bash
# Sets up the micromamba conda environment for this project.
# Run once from the repo root: bash scripts/setup_env.sh
set -euo pipefail

MAMBA="${MAMBA:-/opt/homebrew/opt/micromamba/bin/micromamba}"
MAMBA_ROOT="${MAMBA_ROOT_PREFIX:-$HOME/mamba}"
ENV_NAME="triso-env"

if ! command -v "$MAMBA" &>/dev/null; then
  echo "ERROR: micromamba not found at $MAMBA"
  echo "Install with: brew install micromamba"
  exit 1
fi

echo "Creating/updating conda environment '$ENV_NAME' from environment.yml ..."
# --platform osx-64: OpenMC has no native arm64 package; runs under Rosetta 2 on Apple Silicon.
MAMBA_ROOT_PREFIX="$MAMBA_ROOT" "$MAMBA" env create \
  --platform osx-64 \
  --file environment.yml \
  --yes \
  --name "$ENV_NAME" \
  || MAMBA_ROOT_PREFIX="$MAMBA_ROOT" "$MAMBA" env update \
      --platform osx-64 \
      --file environment.yml \
      --name "$ENV_NAME" \
      --yes

echo ""
echo "Done. Activate with:"
echo "  eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate $ENV_NAME"
echo ""
echo "Then download nuclear data:"
echo "  bash scripts/download_data.sh"
