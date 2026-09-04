# Physics Validation

Before burnup calculations begin, the eigenvalue model is checked against three independent physical expectations (`scripts/validate_physics.py`) and a reproducibility audit (`scripts/audit_stage0.py`). The purpose is to confirm that the geometry, materials, and nuclear data are consistent — not to benchmark against a reference code.

---

## Validation checks

**k-eff range.** For a compact-only, all-reflective HALEU geometry at 1200 K with ENDF/B-VIII.0, k-inf is expected in the range 1.05–1.40. The lower bound reflects that the compact has no surrounding graphite moderator blocks — all thermalization comes from the matrix inside the compact (~70 vol%). A full HTGR fuel element with a graphite sleeve would thermalize more effectively and give a higher k-inf.

Result: **1.10184 ± 0.00103 — PASS**

**Flux self-shielding.** The `MaterialFilter` flux tally records volume-integrated track length, so a raw kernel/matrix flux ratio of ~0.038 is dominated by the volume ratio (kernel 2.71%, matrix 70%), not by any physical effect. Dividing by the analytical volume fraction of each layer converts to flux density (neutrons/cm²·s per unit volume), which reveals the actual spatial self-shielding signature.

| Group | Kernel / matrix flux density | Result |
|---|---|---|
| Thermal | 0.984 | PASS — 1.6% depression |
| Epithermal | 0.997 | PASS — 0.3% depression |
| Fast | 1.009 | INFO — elevated; fission neutrons born in kernel |

The fast group is excluded from the PASS/FAIL criteria: the kernel is a net source of fast fission neutrons, so fast flux density exceeding the matrix is physically correct, not a failure.

The dominant self-shielding mechanism here is resonance self-shielding (Dancoff factor reduction of the effective U-238 resonance integral), not spatial flux depression. At r_kernel = 0.0175 cm, the particle radius is ~100× smaller than the thermal neutron mean free path in graphite (~2 cm), so spatial flux gradients within the particle are weak (<2%).

**Energy deposition proxy.** The `heating-local` tally is correctly populated by the ENDF/B-VIII.0 library (which includes kerma coefficients). The fission rate in the kernel (4.95×10⁻¹ per source neutron) confirms that energy release is localised in the fuel — the lower-than-pure-U235 fission fraction (~0.85 for pure thermal U-235) reflects U-238 resonance capture competing with fission in the HALEU kernel.

---

## Reproducibility and depletion readiness

A separate audit (`scripts/audit_stage0.py`) verifies four prerequisites before the burnup calculation is run:

| Check | Result |
|---|---|
| k-eff reproduced within 3σ | **PASS** — 1.10184 ± 0.00103 (Δ = 0.00000 vs reference; tol = 0.00310) |
| Depletion module available | **PASS** — openmc 0.15.3, `CoupledOperator` and `Chain` present |
| Chain file parseable | **PASS** — `chain_endfb71_thermal.xml` loaded, 3819 nuclides |
| Kernel marked depletable | **PASS** — U234/U235/U236/U238 present, `depletable=True` |

The 3σ reproducibility window is chosen deliberately: using 3 × σ_statepoint accounts for the inherent Monte Carlo statistical variance between runs. A fixed absolute tolerance would be too tight for a stochastic code.

The reference k-eff (1.10184 ± 0.00103, ENDF/B-VIII.0 at 1200 K) is recorded in `RESULTS.md` as the canonical Stage 0 baseline.
