#!/usr/bin/env bash
# Downloads the NNDC HDF5 cross-section library (ENDF/B-VII.1 processed for OpenMC).
# Source: OpenMC CI script (openmc-dev/openmc tools/ci/download-xs.sh)
# Library: nndc_hdf5, ~1.4 GB compressed
#
# Run from the repo root: bash scripts/download_data.sh
# Sets OPENMC_CROSS_SECTIONS in .envrc after extraction.
set -euo pipefail

DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data"
CROSS_SECTIONS_XML="$DATA_DIR/nndc_hdf5/cross_sections.xml"
DOWNLOAD_URL="https://anl.box.com/shared/static/teaup95cqv8s9nn56hfn7ku8mmelr95p.xz"

mkdir -p "$DATA_DIR"

if [[ -f "$CROSS_SECTIONS_XML" ]]; then
  echo "Nuclear data already present at $CROSS_SECTIONS_XML"
else
  echo "Downloading NNDC HDF5 library (~1.4 GB, this will take a while) ..."
  curl -L --progress-bar "$DOWNLOAD_URL" | tar -C "$DATA_DIR" -xJ
  echo "Extraction complete."
fi

echo ""
echo "NNDC HDF5 (ENDF/B-VII.1) library installed at:"
echo "  $CROSS_SECTIONS_XML"
echo ""
echo "Add the following to your shell environment before running OpenMC:"
echo "  export OPENMC_CROSS_SECTIONS=$CROSS_SECTIONS_XML"
echo ""
echo "Or create a .envrc file (direnv) with that line to auto-load it."

# Write .envrc if not already set
ENVRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/.envrc"
if [[ ! -f "$ENVRC" ]] || ! grep -q "OPENMC_CROSS_SECTIONS" "$ENVRC"; then
  echo "export OPENMC_CROSS_SECTIONS=$CROSS_SECTIONS_XML" >> "$ENVRC"
  echo ".envrc updated with OPENMC_CROSS_SECTIONS."
fi
