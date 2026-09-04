# Depletion

The depletion calculation (`src/triso/depletion.py`) irradiates the same compact geometry to **15% FIMA over 620.2 days**, matching the AGR-1 irradiation duration and midpoint burnup target. The UCO kernel is the only depletable material; the coating layers and graphite matrix are held at fixed composition throughout.

---

## Power condition and time-step scheme

The compact is operated at a constant **61.6 W/cm³** (compact volume). This is derived from AGR-1 test parameters: 186 W total compact power, compact volume 3.022 cm³, to achieve 15% FIMA in 620 effective full power days (EFPD). Power is normalised via `normalization_mode='fission-q'`, which uses Q-values from the chain file rather than kerma coefficients — this works with any transport library and does not require photon transport to be enabled.

The 620.2-day irradiation is divided into **22 time steps** with step widths chosen by the timescales of the isotopes that drive the early composition transient:

| Phase | Steps | Width | End day | Reason |
|---|---|---|---|---|
| Xe/Sm equilibration | 5 | 1 day | day 5 | Xe-135 equilibrates in ~2 d (t½ = 9.17 h); Sm-149 in ~10 d |
| FP transient | 5 | 10 days | day 55 | Remaining short-lived fission products stabilise |
| Steady burnup | 12 | 47.1 days | day 620.2 | Composition changes slowly; coarse steps sufficient |

---

## Integration method

The burnup calculation uses the **CE/CM predictor-corrector integrator** (`CECMIntegrator`). Each time step makes two transport solves: a predictor using beginning-of-step reaction rates, and a corrector using midpoint-composition rates. This substantially reduces time-integration error compared to a plain first-order predictor, at the cost of doubling transport work per step.

Each transport solve uses 100 batches (30 inactive) and 2,000 particles per batch — reduced from the 5,000-particle eigenvalue run because depletion transport only needs to normalise power and compute one-group reaction rates. The resulting σ(k-eff) of ~30–50 pcm per step is sufficient for this purpose.

```
22 time steps × 2 solves (CE + CM) = 44 transport calculations
Each solve: 100 batches × 2,000 neutrons
Estimated runtime: ~1.5–2.2 hours on a laptop
```

---

## k-inf history

The compact's reactivity declines steadily across the irradiation:

| Metric | Value |
|---|---|
| BOL k-inf | 1.10207 ± 0.00236 |
| EOL k-inf (620.2 d, 15% FIMA) | 0.92479 ± 0.00217 |
| k-inf = 1 crossing | ~290 EFPD (~7% FIMA) |
| Total reactivity swing | −17,728 pcm |

The decline is near-linear and driven primarily by U-235 depletion. No pronounced Pu-239 plateau is visible at ≤15% FIMA: Pu-239 builds up from U-238 capture and partially offsets U-235 loss, but not enough to flatten the reactivity curve significantly at these burnup levels.

---

## Spectral hardening

A 100-group lethargy spectrum tally (1 meV → 20 MeV, log-spaced) is recorded in the UCO kernel at every depletion transport step. Comparing the probability-density spectrum at BOL, MOL (step 15, ~290 EFPD, 7% FIMA), and EOL reveals how the neutron energy distribution shifts with burnup.

The thermal fraction — the fraction of neutrons below 0.625 eV — drops by approximately **50% from BOL to EOL**. Three processes compete:

- **Fission product accumulation** dominates. Xe-135 (σ_abs ≈ 2.6 × 10⁶ b at thermal) and Sm-149, plus longer-lived fission products, selectively absorb thermal neutrons without contributing to fission. More neutrons are captured before thermalising, so proportionally fewer appear in the thermal peak at EOL.
- **U-235 depletion** acts as a weak opposing effect. U-235 is itself a thermal absorber; as it burns away, there are fewer thermal absorbers in the fuel — a slight softening tendency. Outweighed by fission product accumulation at ≤15% FIMA.
- **Pu-239 buildup** is primarily a reactivity effect at these burnup levels; its contribution to the spectral shape is minor.

The compact-only geometry produces a much harder baseline spectrum than a full HTGR fuel element (BOL thermal fraction ~3% vs. ~20–40% for a full element with graphite sleeve), because there are no surrounding moderator blocks. The hardening signal is the *change* in thermal fraction across burnup, not its absolute value.
