# TRISO Neutronics and Depletion

OpenMC-based neutronics and burnup model for TRISO particle fuel, based on the AGR-1 irradiation test geometry.

---

## What this models

A single AGR-1 TRISO fuel compact — a graphite cylinder (r = 0.62 cm, h = 2.5 cm) packed with ~3,650 TRISO particles at 30% packing fraction. Each particle is a five-layer concentric sphere: UCO fuel kernel (19.75 wt% HALEU) surrounded by a porous carbon buffer, inner PyC, SiC, and outer PyC coating layers. The compact is enclosed by reflective boundary conditions on all faces, approximating an infinite periodic lattice (k-inf geometry).

The model computes:
- Neutron multiplication (k-inf) and flux distribution across TRISO layers
- Fuel depletion to 15% FIMA over ~620 days at 61.6 W/cm³
- Spectral hardening from BOL to EOL (thermal fraction vs. burnup)
- Fuel temperature coefficient (Doppler feedback) at BOL and EOL
- Manufacturing tolerance sensitivity (packing fraction, kernel radius)

---

## How it works

The codebase is structured as a Python package (`src/triso`) with analysis scripts in `scripts/`.

| Module | Role |
|---|---|
| `src/triso/materials.py` | Defines all six material layers as `openmc.Material` objects |
| `src/triso/geometry.py` | Builds the TRISO compact geometry via `pack_spheres` and a rectilinear lattice |
| `src/triso/run.py` | Assembles and runs the eigenvalue model; returns k-eff and tally DataFrames |
| `src/triso/depletion.py` | Configures and runs the CE/CM predictor-corrector burnup calculation |

All scripts in `scripts/` are standalone and import from `src/triso`. Nuclear data (ENDF/B-VIII.0, ~2–3 GB after trimming) and simulation output (HDF5 statepoints, `depletion_results.h5`) are not committed to the repository.

---

## Quickstart

### First-time setup (once per machine)

```bash
# Install micromamba and direnv
brew install micromamba direnv

# Hook direnv into your shell
echo 'eval "$(direnv hook zsh)"' >> ~/.zshrc   # or ~/.bashrc for bash
source ~/.zshrc

# Register the mamba envs directory so 'triso-env' resolves by name
micromamba config append envs_dirs ~/mamba/envs

# Create the conda environment (OpenMC 0.15.3, osx-64 via Rosetta 2)
bash scripts/setup_env.sh

# Activate the environment and install the package in editable mode
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env
pip install -e .

# Download ENDF/B-VIII.0 nuclear data (trimmed to 900 K + 1200 K, ~2–3 GB)
# Allow 20–60 min for the download.
caffeinate bash scripts/download_endfb80.sh

# Download the depletion chain file (~27 MB, required for burnup calculations)
bash scripts/download_chain.sh

# Trust the .envrc so direnv sets OPENMC_CROSS_SECTIONS automatically on cd
direnv allow
```

### Every run

```bash
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env

# Eigenvalue run — k-eff, flux, and tally results (~5–15 min)
caffeinate python -m triso.run

# Depletion run — burnup to 15% FIMA (~1.5–2.2 hr)
# Run from a dedicated output directory to keep statepoints organised
mkdir -p output/depletion && cd output/depletion
caffeinate python -m triso.depletion ../../data/chain_endfb71_thermal.xml
cd ../..

# Physics validation — k-eff range, flux self-shielding checks
python scripts/validate_physics.py

# k-inf vs. burnup plot — saves output/keff_vs_burnup.png
python scripts/plot_depletion.py

# Spectral hardening — saves output/spectrum_hardening.png
python scripts/plot_spectrum.py

# Fuel temperature coefficient / Doppler (~15–30 min, 6 eigenvalue runs)
# Requires a completed depletion run in output/depletion/
caffeinate python scripts/doppler_coefficient.py

# Manufacturing tolerance sensitivity (~25–75 min, 5 eigenvalue runs)
caffeinate python scripts/manufacturing_tolerance.py
```

---

## Implementation details

| Topic | Description |
|---|---|
| [Environment setup](docs/01-environment-setup.md) | conda environment, OpenMC install (Rosetta 2), ENDF/B-VIII.0 nuclear data download and trimming |
| [Neutronics model](docs/02-neutronics-model.md) | Materials (UCO kernel, TRISO coating layers, graphite matrix), geometry (particle packing, rectilinear lattice), eigenvalue run settings, and tallies |
| [Physics validation](docs/03-physics-validation.md) | k-eff range checks, self-shielding flux density analysis, Stage 0 handoff audit, depletion-readiness confirmation |
| [Depletion](docs/04-depletion.md) | Burnup schedule (22 time steps, CE/CM integrator), k-inf history to 15% FIMA, spectral hardening (thermal fraction BOL → EOL) |
| [Reactor physics analysis](docs/05-reactor-physics-analysis.md) | Fuel temperature coefficient (Doppler, BOL and EOL), manufacturing tolerance sensitivity (packing fraction ±2 pp, kernel radius ±5 µm) |

---

## Key results

| Quantity | Value |
|---|---|
| BOL k-inf | 1.102 ± 0.001 |
| EOL k-inf (15% FIMA, 620 EFPD) | 0.925 ± 0.002 |
| k-inf = 1 crossing | ~290 EFPD (~7% FIMA) |
| BOL fuel temperature coefficient | −10.0 pcm/K |
| EOL fuel temperature coefficient | −7.0 pcm/K |
| Thermal fraction drop (BOL → EOL) | ~50% (spectral hardening from fission product accumulation) |
