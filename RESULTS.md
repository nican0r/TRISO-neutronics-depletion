# Stage 0 Results

Reference baseline for the Stage 0 handoff audit (step-5).
Produced with `python -m triso.run` using the configuration described below.

## Eigenvalue run

| Metric | Value |
|--------|-------|
| k-eff | **1.10184 ± 0.00103** |
| Cross-section library | ENDF/B-VIII.0 (`data/endfb80_hdf5/cross_sections.xml`) |
| Temperature | 1200 K (uniform — kernel, layers, and matrix) |
| Total batches | 200 |
| Inactive batches | 50 |
| Particles per batch | 5 000 |
| Active particle-histories | 750 000 |
| Statepoint file | `statepoint.200.h5` |

## Geometry

| Metric | Value |
|--------|-------|
| Target packing fraction | 0.30 |
| Compact radius | 0.62 cm |
| Compact height | 2.50 cm |
| Particle outer radius (OPyC) | 0.0390 cm |
| Boundary conditions | All-reflective (k-inf lattice) |

> **Note**: `pack_spheres` is stochastic; the actual packing fraction
> varies ±0.001–0.005 per geometry build. The target (0.30) is the
> design value; exact particle count differs between runs.

## Depletion run (step-6)

Produced with `python -m triso.depletion data/chain_endfb71_thermal.xml` from `output/depletion/`.

| Setting | Value |
|---|---|
| Power density | 61.6 W/cm³ (compact volume) |
| Target burnup | 15% FIMA |
| Total irradiation | 620.2 days (≈ AGR-1 duration) |
| Time-step scheme | 5×1d + 5×10d + 12×47.1d (22 steps) |
| Integrator | CECMIntegrator (CE/CM predictor-corrector, 2 solves/step) |
| Particles per batch | 2 000 |
| Batches per step | 100 (30 inactive) |
| Chain file | `data/chain_endfb71_thermal.xml` (ENDF/B-VII.1, 3819 nuclides) |
| Normalization | `fission-q` |
| Results file | `output/depletion/depletion_results.h5` |

### k-eff vs. burnup

| Step | Time [d] | k-eff | σ |
|---:|---:|---:|---:|
| 0 | 0.0 | 1.10207 | 0.00236 |
| 1 | 1.0 | 1.09526 | 0.00238 |
| 2 | 2.0 | 1.09430 | 0.00231 |
| 3 | 3.0 | 1.09251 | 0.00221 |
| 4 | 4.0 | 1.08972 | 0.00247 |
| 5 | 5.0 | 1.08803 | 0.00197 |
| 6 | 15.0 | 1.08469 | 0.00250 |
| 7 | 25.0 | 1.08013 | 0.00229 |
| 8 | 35.0 | 1.08176 | 0.00223 |
| 9 | 45.0 | 1.07763 | 0.00200 |
| 10 | 55.0 | 1.07382 | 0.00227 |
| 11 | 102.1 | 1.05999 | 0.00214 |
| 12 | 149.2 | 1.04757 | 0.00227 |
| 13 | 196.3 | 1.02931 | 0.00231 |
| 14 | 243.4 | 1.01512 | 0.00222 |
| 15 | 290.5 | 1.00123 | 0.00203 |
| 16 | 337.6 | 0.98676 | 0.00191 |
| 17 | 384.7 | 0.98139 | 0.00204 |
| 18 | 431.8 | 0.96917 | 0.00204 |
| 19 | 478.9 | 0.95444 | 0.00203 |
| 20 | 526.0 | 0.94552 | 0.00187 |
| 21 | 573.1 | 0.93303 | 0.00207 |
| 22 | 620.2 | 0.92479 | 0.00217 |

k-eff crosses unity at approximately **day 290** (end-of-life for the compact in isolation).
k-eff swing: **1.102 → 0.925** over 620 days (Δk ≈ −0.177, ≈ −285 pcm/day average).

---

## Audit tolerance

The handoff audit (`scripts/audit_stage0.py`) accepts a re-run as
matching if |k_new − 1.10184| ≤ 3 σ_new (approximately 3 × 0.001 ≈ 0.003).
