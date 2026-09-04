# Environment Setup

The simulation runs under **OpenMC 0.15.3**, installed via micromamba from conda-forge into an isolated `triso-env` environment (`environment.yml`). On Apple Silicon, no arm64 binary exists on conda-forge, so the environment is created with `--platform osx-64` and runs transparently under Rosetta 2 — this has no effect on results.

Nuclear data are downloaded separately (not committed to the repository) and live under `data/`. The `OPENMC_CROSS_SECTIONS` environment variable must point to the active library's `cross_sections.xml` before running any simulation. The `.envrc` file sets this automatically when using direnv.

---

## Nuclear data library

The active library is **ENDF/B-VIII.0**, downloaded by `scripts/download_endfb80.sh`. It is preferred over the older NNDC ENDF/B-VII.1 library for two reasons:

- It includes **kerma coefficients**, required to populate the `heating-local` energy deposition tally.
- It ships with cross sections pre-tabulated at **multiple temperatures**, covering the HTGR operating range.

The full ENDF/B-VIII.0 distribution contains six temperature points (250 K, 293.6 K, 600 K, 900 K, 1200 K, 2500 K) and weighs ~7–10 GB. `scripts/trim_temperatures.py` (called automatically by the download script) strips all points except **900 K and 1200 K**, reducing the library to ~2–3 GB. The `0K` energy-grid groups embedded in each HDF5 file are always preserved — these are the base cross-section energy grids needed for interpolation, not a physical 0 K dataset.

Temperatures outside the 900–1200 K tabulated range require `Settings.temperature = {'method': 'interpolation'}` to use sqrt(T) interpolation between the two bounding points.

The NNDC ENDF/B-VII.1 library (`scripts/download_data.sh`) is retained for reference only. It has no kerma and only a single 293.6 K temperature point.

## Depletion chain

Transmutation and decay paths come from `data/chain_endfb71_thermal.xml`, the ENDF/B-VII.1 thermal depletion chain (3,819 nuclides), downloaded by `scripts/download_chain.sh`. No pre-built ENDF/B-VIII.0 chain is publicly distributed by the OpenMC project — generating one requires running full NJOY processing on the raw ENDF files. Using the VII.1 chain with VIII.0 transport cross sections is standard practice and introduces negligible error for the burnup levels modelled here.
