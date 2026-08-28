#!/usr/bin/env bash
# Downloads the ENDF/B-VIII.0 HDF5 library for OpenMC (multi-temperature, with kerma).
# After extraction, trims to 900 K and 1200 K only to reduce disk footprint.
#
# Full library temperatures: 250 K, 293.6 K, 600 K, 900 K, 1200 K, 2500 K
# Source: openmc.org official data libraries (hosted on ANL Box)
#
# Run from the repo root: bash scripts/download_endfb80.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$REPO_ROOT/data"
LIB_DIR="$DATA_DIR/endfb80_hdf5"
CROSS_SECTIONS_XML="$LIB_DIR/cross_sections.xml"
DOWNLOAD_URL="https://anl.box.com/shared/static/uhbxlrx7hvxqw27psymfbhi7bx7s6u6a.xz"

mkdir -p "$LIB_DIR"

if [[ -f "$CROSS_SECTIONS_XML" ]]; then
    echo "ENDF/B-VIII.0 library already present at $CROSS_SECTIONS_XML"
else
    echo "Downloading ENDF/B-VIII.0 library (~7–10 GB compressed; allow 20–60 min) ..."
    # --strip-components=1 drops the tarball's top-level directory name so files
    # land directly in $LIB_DIR regardless of what that directory is called inside
    # the archive.
    curl -L --progress-bar "$DOWNLOAD_URL" | tar -C "$LIB_DIR" -xJ --strip-components=1
    echo "Extraction complete."
fi

echo ""
echo "Trimming to 900 K and 1200 K (removing 250 K, 294 K, 600 K, 2500 K) ..."
# Resolve Python from the active conda environment, falling back to the system python3.
PYTHON_BIN="${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python)}"
"$PYTHON_BIN" "$REPO_ROOT/scripts/trim_temperatures.py" "$LIB_DIR" 900 1200
echo "Trimming complete."

echo ""
echo "ENDF/B-VIII.0 (900 K + 1200 K) library installed at:"
echo "  $CROSS_SECTIONS_XML"
echo ""
echo "Updating .envrc to use the new library ..."
ENVRC="$REPO_ROOT/.envrc"
# Replace existing OPENMC_CROSS_SECTIONS line (or append if absent)
if grep -q "OPENMC_CROSS_SECTIONS" "$ENVRC" 2>/dev/null; then
    sed -i.bak "s|export OPENMC_CROSS_SECTIONS=.*|export OPENMC_CROSS_SECTIONS=$CROSS_SECTIONS_XML|" "$ENVRC"
    rm -f "$ENVRC.bak"
else
    echo "export OPENMC_CROSS_SECTIONS=$CROSS_SECTIONS_XML" >> "$ENVRC"
fi
echo ".envrc updated. Run 'direnv allow' to reload, or:"
echo "  export OPENMC_CROSS_SECTIONS=$CROSS_SECTIONS_XML"
