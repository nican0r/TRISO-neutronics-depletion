#!/usr/bin/env bash
# Download a pre-built OpenMC depletion chain file and place it at
# data/chain_endfb71_thermal.xml.
#
# The chain file is required by openmc.deplete.CoupledOperator and is NOT
# bundled with OpenMC — it must be obtained separately.
#
# Chain used: ENDF/B-VII.1 thermal-spectrum chain from the openmc-dev/data
# repository. This is the best pre-built chain available without running
# NJOY to process raw ENDF/B-VIII.0 data. Using the VII.1 chain with
# VIII.0 transport cross sections is standard practice for Stage 0 and
# introduces negligible error in predicted burnup for short irradiation
# times. A native VIII.0 chain can be generated later with
# depletion/generate_endf80_chain.py in the openmc-dev/data repo.
#
# Approximate uncompressed size: ~27 MB.
# Source: https://github.com/openmc-dev/data (master branch)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST_DIR="$REPO_ROOT/data"
DEST_FILE="$DEST_DIR/chain_endfb71_thermal.xml"
URL="https://raw.githubusercontent.com/openmc-dev/data/master/depletion/chain_endfb71_thermal.xml"

mkdir -p "$DEST_DIR"

if [[ -f "$DEST_FILE" ]]; then
    echo "Chain file already present at $DEST_FILE — skipping download."
    exit 0
fi

echo "Downloading ENDF/B-VII.1 thermal depletion chain file..."
echo "  Source: $URL"
echo "  Destination: $DEST_FILE"

curl -fSL --progress-bar "$URL" -o "$DEST_FILE"

echo "Done. Chain file written to: $DEST_FILE"
echo "Nuclide count:"
grep -c '<nuclide ' "$DEST_FILE" || true
