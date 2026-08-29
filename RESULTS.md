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

## Audit tolerance

The handoff audit (`scripts/audit_stage0.py`) accepts a re-run as
matching if |k_new − 1.10184| ≤ 3 σ_new (approximately 3 × 0.001 ≈ 0.003).
