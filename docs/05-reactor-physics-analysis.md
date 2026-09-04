# Reactor Physics Analysis

Two sensitivity analyses are run against the baseline eigenvalue model: a fuel temperature coefficient (Doppler feedback) sweep, and a manufacturing tolerance study on packing fraction and kernel radius. Both use the same compact geometry and ENDF/B-VIII.0 nuclear data as the baseline run.

---

## Fuel temperature coefficient (Doppler)

The fuel temperature coefficient (FTC) measures how k-inf changes as the UCO kernel temperature changes, with all other materials held fixed. Isolating the kernel temperature from the moderator temperature ensures the coefficient reflects only the Doppler broadening of fuel resonances, not a mixed fuel+moderator response.

**Method.** Three eigenvalue calculations are run at kernel temperatures of 900 K, 1050 K, and 1200 K, with the buffer, PyC, SiC, and graphite matrix held at 1200 K. The analysis is repeated at both BOL (fresh fuel) and EOL (15% FIMA kernel composition from the depletion run). The FTC in pcm/K is extracted by linear regression of k-inf vs. temperature over the three points.

At 900 K and 1200 K the tabulated cross sections are used directly. At 1050 K, OpenMC uses sqrt(T) interpolation between the two bounding points — the standard treatment and well-conditioned at the midpoint. The trimmed library does not support extrapolation above 1200 K, so the FTC is evaluated over the [900, 1200] K interval only.

**Results.**

| Burnup | T_fuel | k-inf | σ (pcm) |
|---|---|---|---|
| BOL | 900 K | 1.12920 | 236 |
| BOL | 1050 K | 1.11969 | 215 |
| BOL | 1200 K | 1.09918 | 195 |
| EOL | 900 K | 0.94894 | 249 |
| EOL | 1050 K | 0.93662 | 198 |
| EOL | 1200 K | 0.92797 | 212 |

| Burnup | FTC (pcm/K) |
|---|---|
| BOL | −10.0 |
| EOL | −7.0 |

The sign is negative at both burnup points, as expected: Doppler broadening of U-238 resonances increases parasitic absorption as temperature rises, reducing reactivity. The mechanism is a fundamental safety feature of HTGR fuel — it provides immediate negative feedback when the fuel heats up.

The magnitude (−10 to −7 pcm/K) exceeds the published full-core HTGR range of −2 to −6 pcm/K (Kuijper et al., NSE 153 (2006); IAEA-TECDOC-978). This is expected: published values are for compacts surrounded by large graphite moderator blocks that thermalise the spectrum. The compact-only reflective geometry has far less moderator-to-fuel ratio, producing a harder neutron spectrum (see [Depletion — Spectral hardening](04-depletion.md)). In a harder spectrum, more neutrons pass through the U-238 resolved resonance region (1 eV – 100 keV) before thermalising, amplifying the Doppler effect. The compact-only FTC cannot be compared directly to full-core published values without adding the graphite sleeve and moderator blocks to the geometry.

The FTC weakens from BOL (−10.0 pcm/K) to EOL (−7.0 pcm/K) because U-238 — the primary Doppler absorber — partially depletes across the irradiation and is partially replaced by Pu-239. Pu-239 contributes Doppler feedback but has fewer and narrower resonances than U-238 in the resolved region, so the net coefficient decreases in magnitude with burnup.

---

## Manufacturing tolerance sensitivity

This analysis quantifies how much k-inf changes when the two most variable AGR-1 fabrication parameters — packing fraction and kernel diameter — are shifted by their manufacturing specification tolerances.

**Perturbation scheme.** Five BOL eigenvalue calculations are run: a nominal baseline and four perturbations, each varying one parameter at a time.

| Case | Packing fraction | Kernel ⌀ |
|---|---|---|
| Nominal | 30% | 350 µm |
| PF − 2 pp | 28% | 350 µm |
| PF + 2 pp | 32% | 350 µm |
| Kernel − 10 µm | 30% | 340 µm |
| Kernel + 10 µm | 30% | 360 µm |

The ±2 percentage-point packing fraction shift is an absolute change (28% and 32%), consistent with how reactor physics sensitivity studies report perturbations. The ±10 µm kernel diameter shift (±5 µm radius) comes directly from the AGR-1 manufacturing specification in INL/EXT-10-19476 Table 3.

For the kernel perturbations, only the kernel sphere radius changes. All coating-layer outer radii are held at nominal values, so the buffer layer absorbs the delta and the OPyC outer radius — which sets the lattice pitch — remains unchanged. This isolates the fuel-volume effect from any packing geometry effect.

Each run uses the same settings as the baseline eigenvalue run (200 batches, 50 inactive, 5,000 particles per batch), giving σ(k-eff) ≈ 10–20 pcm and Δk statistical uncertainty ≈ 14–28 pcm. This is sufficient to resolve sensitivity differences of a few hundred pcm with better than 10% precision.

The analysis is BOL only. Manufacturing tolerances affect the initial fuel loading; the BOL k-eff sensitivity is the most direct measure of how fabrication variation translates to reactivity uncertainty. How this sensitivity evolves with burnup is deferred to a later stage.
