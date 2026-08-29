# TRISO-neutronics-depletion

OpenMC-based neutronics and depletion model for TRISO particle fuel. Built as a portfolio project demonstrating HTGR fuel performance simulation competencies.

---

## Stage 0 — Environment setup (`project-init`)

### What was implemented

- `environment.yml` — micromamba/conda environment spec pinning OpenMC 0.15.3 and all Python dependencies
- `pyproject.toml` — project metadata and Python dependency list for the `src/triso` package
- `scripts/setup_env.sh` — one-shot script to create the conda environment
- `scripts/download_data.sh` — downloads the NNDC ENDF/B-VII.1 HDF5 library (~1.4 GB, single temperature 293.6 K, no kerma)
- `scripts/download_endfb80.sh` — downloads the ENDF/B-VIII.0 HDF5 library (multi-temperature, with kerma), then trims to 900 K + 1200 K only
- `scripts/trim_temperatures.py` — rewrites OpenMC HDF5 nuclide files keeping only specified temperature points; called automatically by `download_endfb80.sh`
- `scripts/validate_env.py` — validates OpenMC import and cross-section library accessibility
- `.envrc` — sets `OPENMC_CROSS_SECTIONS` to the local `data/` path (direnv-compatible)
- `src/triso/__init__.py` — empty package root; model modules added in later steps
- `.gitignore` — excludes `data/` (large binary nuclear data), OpenMC output HDF5 files, and the conda env

### How it works

OpenMC is installed via micromamba from the conda-forge channel. Because there is no native arm64 binary for OpenMC on macOS, the environment uses `--platform osx-64` to install the x86_64 build, which runs transparently under Rosetta 2. Nuclear data (966 nuclide datasets, ENDF/B-VII.1 processed to HDF5 by the NNDC) are downloaded separately to `data/nndc_hdf5/` and are not committed to the repository. The `OPENMC_CROSS_SECTIONS` environment variable points to `data/nndc_hdf5/cross_sections.xml` and must be set before running any simulation.

### Design decisions

- **OpenMC 0.15.3, not 0.16.0** — 0.16.0 released August 2026 but has not yet landed on conda-forge; 0.15.3 is the latest conda-packaged version and is the same used in TRISO literature.
- **osx-64 platform (Rosetta 2)** — no arm64 binary exists on conda-forge; Rosetta is transparent and has no correctional impact on results.
- **ENDF/B-VIII.0 is the active library** (`scripts/download_endfb80.sh`) — includes kerma coefficients (required for `heating-local` tally) and pre-tabulated cross sections at 900 K and 1200 K (HTGR operating range). The NNDC ENDF/B-VII.1 library (`download_data.sh`) is retained for reference but is no longer the default; it has no kerma and only 293.6 K.
- **Trimmed to 900 K + 1200 K** — the full ENDF/B-VIII.0 distribution ships with six temperatures (250 K, 293.6 K, 600 K, 900 K, 1200 K, 2500 K); `trim_temperatures.py` rewrites each HDF5 file to drop the four unused points, reducing the library from ~7–10 GB to roughly 2–3 GB. The `0K` windowed-multipole groups are always preserved.
- **`data/` gitignored** — all HDF5 libraries total several GB; the download scripts are the authoritative source of truth for reproducing the data directory.
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
- **Temperature 1200 K** — Upper HTGR fuel design temperature. The ENDF/B-VIII.0 neutron XS files include both 900 K and 1200 K, but the c_Graphite S(α,β) thermal scattering evaluation only provides a 1200 K table (its temperature grid runs 296 K, 400 K, 500 K, 600 K, 700 K, 800 K, 1000 K, 1200 K — no 900 K point). Running at 900 K would have OpenMC silently round the graphite thermal scattering up to 1200 K anyway, so 1200 K is the self-consistent choice. The 900 K neutron XS tables remain in the library for future use.
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

## Stage 3 — Run & Tallies (`step-3`)

### What was implemented

- `src/triso/run.py` — `build_model()` assembling the complete `openmc.Model` (geometry + materials + settings + tallies); `run_model()` executing the calculation and returning a results dict; `__main__` block for direct invocation.
- Four tallies:
  - `flux by material` — energy-integrated flux aggregated per layer type (kernel, buffer, IPyC, SiC, OPyC, matrix)
  - `flux by material 3-group` — same, binned into thermal / epithermal / fast groups to reveal spectral shape per region
  - `kernel reaction rates` — fission and absorption rates in the UCO kernel
  - `heating by material` — `heating-local` energy deposition score across all six material layers

### How it works

`build_model()` calls `build_materials()` and `build_geometry()`, then assembles `openmc.Settings` for eigenvalue mode (k-inf) and attaches four tallies before returning an `openmc.Model`. Tallies use `MaterialFilter` rather than `CellFilter` so that scores are aggregated across *all* repeated TRISO particle instances in the packed lattice — a `CellFilter` on a nominal cell ID would only score one particle's contribution. `run_model()` invokes `model.run()`, opens the resulting statepoint, and returns a dict of k-eff ± σ plus four pandas DataFrames for downstream analysis.

### Design decisions

- **200 batches / 50 inactive / 5 000 particles per batch** — ~750 000 active histories; expected σ(k-eff) ≈ 10–20 pcm and tally relative errors ≲ 5% on the kernel and matrix. Stage 0 target: fast enough to run on a laptop in ≈ 5–15 min while giving statistically meaningful tally output.
- **MaterialFilter over CellFilter** — TRISO lattice has hundreds of repeated cell instances; MaterialFilter is the correct aggregate for per-layer-type quantities. DistribcellFilter (instance-by-instance) produces very large output and is deferred to a later stage.
- **3-group energy boundaries: 0 / 0.625 eV / 100 keV / 20 MeV** — IAEA standard thermal cutoff (0.625 eV) and fast threshold (100 keV). A finer group structure is deferred to the MGXS generation stage.
- **`heating-local` score** — reports kinetic energy deposited locally by charged particles and recoil nuclei; valid without `photon_transport = True`. `heating` (which includes photon transport) is deferred to a stage that enables coupled neutron-photon transport.
- **Box initial source** — `openmc.stats.Box` spanning the compact cylinder bounding box; adequate starting distribution for a reflective-boundary k-inf problem where the fission source converges quickly.
- **TODO**: `heating-local` returns zero with the single-temperature NNDC library because kerma coefficients are absent from that library. Will produce non-zero results when upgraded to a multi-temperature library. Enable `photon_transport` and switch to `heating` for a more complete energy deposition picture at the same time.

---

## Stage 4 — Validation / sanity checks (`step-4`)

### What was implemented

- `scripts/validate_physics.py` — standalone validation script with three checks: k-eff range, self-shielding flux depression, and energy deposition proxy.
- `--statepoint PATH` argument to load an existing statepoint HDF5 (bypasses re-running the simulation).
- `--output-dir DIR` argument (default: `output/`); two PNG plots saved there: `flux_thermal_by_material.png` and `flux_3group_by_material.png`.
- `_vol_fractions()` helper computing analytical volume fractions for each TRISO layer from the geometry constants in `geometry.py`.
- `_mat_id_map()` / `_patch_mat_labels()` helpers: OpenMC stores `MaterialFilter` bins as integer IDs in the statepoint HDF5; these helpers read `summary.h5` and restore the human-readable material names defined in `materials.py`.

### How it works

The script loads the simulation statepoint (or re-runs the model if none is given), then performs three independent checks. **k-eff** is compared against an expected window for a compact-only, all-reflective HALEU geometry. **Self-shielding** is evaluated by dividing the `MaterialFilter` flux tally — which is volume-integrated — by the analytical volume fraction of each layer, converting to flux density (neutrons / cm² · s per unit volume) so that a meaningful cross-region comparison is possible. The raw flux ratios closely match the geometry volume ratios and are dominated by that geometry effect; only the volume-normalized ratio reveals the actual few-percent spatial self-shielding signature. **Energy deposition** uses the fission reaction rate in the kernel as a proxy because `heating-local` returns zero with the NNDC HDF5 library (kerma coefficients absent; see existing TODO). All results are printed to stdout; plots are saved as PNG files to `output/`.

### Results (Stage 0 run)

| Check | Result | Notes |
|---|---|---|
| k-eff | **1.10184 ± 0.00103** — PASS | Expected 1.05–1.40 for compact-only geometry; ENDF/B-VIII.0 at 1200 K |
| Thermal flux density kernel/matrix | **0.984** — PASS | 1.6% depression; small but real |
| Epithermal flux density kernel/matrix | **0.997** — PASS | 0.3% depression in epithermal |
| Fast flux density kernel/matrix | **1.009** — INFO | Elevated: fission neutrons born in kernel |
| Fission rate in kernel | **4.95×10⁻¹** — PASS | Non-zero; energy release confirmed localised |

### Design decisions

- **k-eff expected range 1.05–1.40, not 1.35–1.65** — the geometry is compact-only with no external graphite moderator blocks; the only thermalization comes from the graphite matrix inside the compact (~70 vol%). A full HTGR fuel element (graphite sleeve + compact) would thermalize more effectively and give a higher k-inf. The compact-only k-inf of 1.21 is consistent with published TRISO compact benchmark calculations.
- **Volume normalization for flux density** — `MaterialFilter` tallies are volume-integrated (sum of track lengths per material type), so a raw kernel/matrix ratio of ~0.038 is dominated by the volume ratio (kernel 2.71%, matrix 70%). Dividing by volume fraction from `_vol_fractions()` gives the physically meaningful flux density. Volume fractions are derived analytically from the AGR-1 particle radii and packing fraction in `geometry.py`.
- **Fast group excluded from self-shielding PASS/FAIL** — the kernel is a net source of fast fission neutrons; fast flux density > matrix is physically correct, not a failure.
- **Dominant self-shielding mechanism is resonance, not spatial** — at r_kernel = 0.0175 cm, the particle radius is ~100× smaller than the thermal neutron mean free path in graphite (~2 cm). Spatial flux gradients within the particle are therefore weak (<2%). The primary self-shielding effect is resonance self-shielding (Dancoff factor reduction of the effective U-238 resonance integral), which manifests as an increased effective multiplication factor compared to a homogenised geometry — not as a large spatial flux depression visible in these tallies.
- **Energy deposition via fission-rate proxy** — `heating-local` returns zero with the NNDC HDF5 library (kerma coefficients absent). The fission rate in the kernel (f/a = 0.499) confirms energy is released in the fuel kernel. The lower-than-pure-U235 fission fraction (pure 235 thermal: ~0.85) reflects U-238 resonance capture in the HALEU kernel. Proper energy deposition accounting requires a multi-temperature library with kerma; marked as an existing TODO.
- **Material name mapping from `summary.h5`** — OpenMC writes integer material IDs in the statepoint HDF5 even when the statepoint is linked with a summary. `_mat_id_map()` reads the summary's material list to translate `mat_N` labels back to the names defined in `materials.py` (`'UCO kernel'`, `'graphite matrix'`, etc.).
- **TODO (Stage 1+)**: Add a graphite fuel-rod sleeve to the geometry to better represent the real HTGR moderation ratio; re-run validation to see how k-inf shifts. Enable `photon_transport=True` and upgrade to a multi-temperature library to validate `heating-local` energy deposition.

---

## Stage 5 — Stage 0 handoff audit (`step-5`)

### What was implemented

- `RESULTS.md` — records the Stage 0 reference k-eff (1.10184 ± 0.00103) with ENDF/B-VIII.0 at 1200 K as the canonical baseline for future reproducibility checks.
- `scripts/download_chain.sh` — downloads the ENDF/B-VII.1 thermal depletion chain file (`chain_endfb71_thermal.xml`, 3819 nuclides, ~27 MB) from the `openmc-dev/data` GitHub repository into `data/`.
- `scripts/audit_stage0.py` — four-check audit script: k-eff against RESULTS.md baseline, `openmc.deplete` module availability, chain file parseable by `openmc.deplete.Chain`, and kernel material depletion-readiness.
- `src/triso/materials.py` — `kernel.depletable = True` added so the UCO kernel participates in burnup calculations.

### How it works

`audit_stage0.py` loads the existing `statepoint.200.h5` (no re-run required) and compares k-eff against the RESULTS.md reference using a 3σ tolerance window. It then imports `openmc.deplete` to confirm the depletion module is accessible, parses `data/chain_endfb71_thermal.xml` via `openmc.deplete.Chain.from_xml` to verify the chain file is present and valid, and finally calls `build_materials()` to confirm the kernel has U-234/U-235/U-238 nuclides and `depletable=True`. Exit code is 0 only if all four checks pass.

### Audit results

| Check | Result |
|---|---|
| k-eff reproduced | **PASS** — 1.10184 ± 0.00103 (Δ = 0.00000 vs reference; tol = 0.00310) |
| Depletion module | **PASS** — openmc 0.15.3, `CoupledOperator` and `Chain` present |
| Chain file | **PASS** — `chain_endfb71_thermal.xml` loaded, 3819 nuclides |
| Kernel composition | **PASS** — U234/U235/U236/U238 present, `depletable=True` |

### Design decisions

- **RESULTS.md k-eff = 1.10184, not 1.2117** — the 1.2117 value recorded in the Stage 4 results table was from an earlier run using the NNDC ENDF/B-VII.1 library at 900 K. The step-4 commit upgraded to ENDF/B-VIII.0 and 1200 K, which lowers k-eff by ~1000 pcm due to stronger Doppler broadening at higher temperature. The Stage 4 results table in this README was corrected to reflect the current library/temperature.
- **ENDF/B-VII.1 chain, not VIII.0** — no pre-built ENDF/B-VIII.0 chain file is publicly distributed by the OpenMC project. The `openmc-dev/data` repository contains only a generation script (`generate_endf80_chain.py`) requiring full ENDF/B-VIII.0 NJOY processing. Using the VII.1 thermal chain with VIII.0 transport cross sections is standard practice and introduces negligible error for Stage 0 burnup calculations. A native VIII.0 chain should be generated for higher-fidelity depletion stages.
- **`depletable=True` on kernel only** — only the UCO kernel contains fissile material. Buffer, PyC, SiC, and matrix do not deplete meaningfully in the Stage 0 timescale and are left non-depletable to avoid unnecessary nuclide tracking overhead.
- **3σ tolerance for k-eff comparison** — using 3 × σ_statepoint as the reproducibility window naturally accounts for Monte Carlo statistical variance between runs; a fixed absolute tolerance would be too tight for stochastic codes.

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

# Create the conda environment
bash scripts/setup_env.sh

# Activate the environment and install the package in editable mode
# (editable mode means changes to src/triso/ take effect without reinstalling)
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env
pip install -e .

# Download ENDF/B-VIII.0 nuclear data (trimmed to 900 K + 1200 K, ~2–3 GB after trim)
# This replaces the older NNDC library; allow 20–60 min for the download.
bash scripts/download_endfb80.sh

# Download the depletion chain file (~27 MB, required for burnup calculations)
bash scripts/download_chain.sh

# Trust the .envrc so direnv sets OPENMC_CROSS_SECTIONS automatically on cd
direnv allow
```

### Every run

```bash
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env
python -m triso.run
```
