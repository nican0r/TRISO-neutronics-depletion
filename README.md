# TRISO-neutronics-depletion

OpenMC-based neutronics and depletion model for TRISO particle fuel. Built as a portfolio project demonstrating HTGR fuel performance simulation competencies.

---

## Stage 0 — Environment setup (`project-init`)

### What was implemented

- `environment.yml` — micromamba/conda environment spec pinning OpenMC 0.15.3 and all Python dependencies
- `pyproject.toml` — project metadata and Python dependency list for the `src/triso` package
- `scripts/setup_env.sh` — one-shot script to create the conda environment
- `scripts/download_data.sh` — downloads the NNDC HDF5 cross-section library (~1.4 GB) from ANL Box
- `scripts/validate_env.py` — validates OpenMC import and cross-section library accessibility
- `.envrc` — sets `OPENMC_CROSS_SECTIONS` to the local `data/` path (direnv-compatible)
- `src/triso/__init__.py` — empty package root; model modules added in later steps
- `.gitignore` — excludes `data/` (large binary nuclear data), OpenMC output HDF5 files, and the conda env

### How it works

OpenMC is installed via micromamba from the conda-forge channel. Because there is no native arm64 binary for OpenMC on macOS, the environment uses `--platform osx-64` to install the x86_64 build, which runs transparently under Rosetta 2. Nuclear data (966 nuclide datasets, ENDF/B-VII.1 processed to HDF5 by the NNDC) are downloaded separately to `data/nndc_hdf5/` and are not committed to the repository. The `OPENMC_CROSS_SECTIONS` environment variable points to `data/nndc_hdf5/cross_sections.xml` and must be set before running any simulation.

### Design decisions

- **OpenMC 0.15.3, not 0.16.0** — 0.16.0 released August 2026 but has not yet landed on conda-forge; 0.15.3 is the latest conda-packaged version and is the same used in TRISO literature.
- **osx-64 platform (Rosetta 2)** — no arm64 binary exists on conda-forge; Rosetta is transparent and has no correctional impact on results.
- **NNDC HDF5 library (ENDF/B-VII.1)** — sourced from OpenMC's own CI download script (`tools/ci/download-xs.sh`); this is the most tested library for OpenMC.
- **`data/` gitignored** — 966 HDF5 files total several GB; `download_data.sh` is the authoritative source of truth for reproducing the data directory.
- **micromamba over full conda/miniforge** — single Homebrew binary, no PATH mutation needed; avoids modifying the user's shell init files.
- **TODO (Stage 1+)**: Consider upgrading to ENDF/B-VIII.0 when it becomes straightforwardly available; 235U cross-sections differ by <1% in thermal region but resonance integrals may shift.

---

---

## Stage 1 — Materials (`step-1`)

### What was implemented

- `src/triso/materials.py` — `build_materials()` function returning a `dict[str, openmc.Material]` with six keyed entries: `kernel`, `buffer`, `ipyc`, `sic`, `opyc`, `matrix`.

### How it works

Each TRISO layer is defined as an independent `openmc.Material` object with element atom ratios, density, and temperature. The UCO kernel uses `add_element('U', 1.0, enrichment=19.75)` to specify 19.75 wt% U-235 (HALEU), with C/U = 0.5 and O/U = 0.4 atom ratios (UC₀.₅O₀.₄). All materials are set to 900 K as a single uniform operating temperature. 

The graphite matrix carries the `c_Graphite` S(α,β) thermal scattering kernel (for more realistic modeling vs free-gas approximation); PyC layers do not (see design decisions). Materials are returned as plain Python objects; an `openmc.Model` constructed in later steps will collect them into `openmc.Materials` and handle XML export.

### Design decisions

- **UCO composition UC₀.₅O₀.₄, density 10.5 g/cm³** — AGR-1 as-fabricated targets from INL/EXT-10-19476 (Table 3) and Demkowicz et al., Nucl. Eng. Des. 329 (2018) 102–111.
- **Buffer density 1.0 g/cm³** — AGR-1 target per INL/EXT-10-19476 Table 3; porous layer is ~50% dense relative to solid PyC.
- **IPyC/OPyC density 1.87 g/cm³** — AGR-1 target range 1.85–1.90 g/cm³, midpoint chosen; INL/EXT-10-19476.
- **SiC density 3.20 g/cm³** — AGR-1 target 3.19 g/cm³ rounded to 2 decimal places; INL/EXT-10-19476 Table 3.
- **Graphite matrix density 1.75 g/cm³** — AGR-1 compact matrix target per INL/EXT-10-19476 Table 5.
- **Temperature 900 K** — Single-temperature Stage 0 approximation; HTGR full-power range is 900–1200 K; 900 K is the conservative cool-side estimate. Revisit with multi-temperature model in a later stage.
- **No S(α,β) on PyC layers** — PyC is turbostratic carbon, not crystalline graphite; applying `c_Graphite` would be physically incorrect. Effect is small (thermal scattering in thin PyC layers), marked TODO for Stage 1 review.
- **Enrichment keyword approximation** — `add_element('U', enrichment=19.75)` uses a fixed U-234/U-235 mass ratio of 0.008, valid for centrifuge-enriched product; U-234 ≈ 0.18 wt% of total U, U-236 = 0. Error on k-eff is negligible for Stage 0.
- **`depletable` flag not set** — Left for the depletion step; materials are pure neutronics objects at this stage.

---

## Stage 2 — Geometry (`step-2`)

### What was implemented

- `src/triso/geometry.py` — `build_geometry()` function returning an `openmc.Geometry` containing a single packed TRISO compact.
- `_triso_universe()` — private helper building the five-layer concentric-sphere TRISO particle universe (kernel → buffer → IPyC → SiC → OPyC).
- `build_geometry()` — public entry point: packs TRISO particles into an AGR-1 cylindrical compact at 30% packing fraction, embeds them in a graphite matrix via a rectilinear lattice, and applies reflective boundary conditions on all faces.

### How it works

`openmc.model.pack_spheres` fills the compact cylinder (r = 0.62 cm, h = 2.5 cm) with non-overlapping sphere centers at 30% packing fraction. Each center becomes an `openmc.model.TRISO` object pointing to the shared particle universe. `openmc.model.create_triso_lattice` bins these particles into a rectilinear lattice (16 × 16 × 33 cells, 0.078 cm pitch) with graphite matrix as the background fill. The compact cell is filled with this lattice and enclosed by three reflective surfaces — a cylinder and two end-cap planes — modelling an infinite periodic lattice of identical compacts (k-inf geometry).

### Design decisions

- **AGR-1 nominal particle radii** — kernel 0.0175 cm, buffer outer 0.0275 cm, IPyC 0.0315 cm, SiC 0.0350 cm, OPyC 0.0390 cm; from INL/EXT-10-19476 Table 3.
- **AGR-1 compact dimensions** — r = 0.62 cm, h = 2.50 cm; from INL/EXT-10-19476.
- **30% packing fraction** — design target per step-2 specification; AGR-1 compacts typically range 25–35%.
- **Lattice pitch = 2 × R_OPyC = 0.078 cm** — rule-of-thumb minimum to prevent a particle from spanning two lattice cells, which would cause geometry errors in OpenMC's cell-finding algorithm.
- **All-reflective boundaries** — standard k-inf lattice approximation; no neutron leakage modelled. Marked `# TODO` for vacuum axial boundaries in a later stage.
- **No helium coolant channel** — adding an annular coolant region requires a second outer bounding cylinder and a graphite sleeve, which meaningfully increases geometry complexity for no Stage 0 benefit. Marked `# TODO` in code.
- **Lattice bounding box computed dynamically** — `nx`, `ny`, `nz` are derived from compact dimensions and pitch so the lattice fully covers the compact regardless of future dimension changes.

---

## Quickstart

```bash
# 1. Install micromamba (macOS, Homebrew)
brew install micromamba

# 2. Create the conda environment
bash scripts/setup_env.sh

# 3. Download nuclear data (~1.4 GB)
bash scripts/download_data.sh

# 4. Activate and validate
eval "$(micromamba shell hook --shell bash)"
micromamba activate triso-env
export OPENMC_CROSS_SECTIONS="$(pwd)/data/nndc_hdf5/cross_sections.xml"
python scripts/validate_env.py
```
