"""TRISO compact depletion: burnup schedule and integration settings.

Power density and irradiation schedule derived from AGR-1 test parameters:
  INL/EXT-10-19476 — AGR-1 Fuel Particle and Compact Characterization Data
  Demkowicz et al., Nucl. Eng. Des. 329 (2018) 102–111

Derivation of POWER_DENSITY (61.6 W/cm³):
  Compact volume = π × 0.62² × 2.5 = 3.022 cm³
  U mass per compact: 3,653 particles × (kernel vol fraction 9.04%)
    × 10.5 g/cm³ × U mass fraction 0.9505 = 0.818 g U
  Avg atomic weight = 0.1975×235 + 0.8025×238 = 237.4 g/mol
  Initial HM atoms = (0.818/237.4) × 6.022e23 = 2.074e21 atoms
  Energy at 15% FIMA = 0.15 × 2.074e21 × 200 MeV × 1.602e-13 J/MeV ≈ 9.97 GJ
  AGR-1 irradiation: 620 EFPD → power = 9.97e9 / (620 × 86400) = 186 W
  Power density = 186 W / 3.022 cm³ = 61.6 W/cm³
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import openmc
import openmc.deplete

import math

from .geometry import (
    _COMPACT_H, _COMPACT_R,
    _R_KERNEL, _R_OPYC, _PACKING_FRACTION,
    build_geometry,
)
from .materials import build_materials

# Analytical kernel volume in the compact [cm³].
# V_compact × packing_fraction × (R_kernel/R_opyc)³ gives the fraction of
# compact volume occupied by kernel material.
# Used to set material.volume for the CoupledOperator (required for converting
# tally reaction rates to atom counts in the Bateman equations).
KERNEL_VOLUME: float = (
    math.pi * _COMPACT_R**2 * _COMPACT_H  # compact cylinder volume
    * _PACKING_FRACTION                   # TRISO particle volume fraction
    * (_R_KERNEL / _R_OPYC) ** 3          # kernel fraction within one particle
)
# ≈ 0.0820 cm³ at nominal AGR-1 geometry and 30% packing fraction.

# ---------------------------------------------------------------------------
# Spectrum tally configuration (shared with scripts/plot_spectrum.py)
# ---------------------------------------------------------------------------
# 100 log-spaced energy bins from 1 meV to 20 MeV — fine enough to resolve
# the thermal peak, epithermal resonance region, and fast fission spectrum.
SPECTRUM_N_GROUPS: int = 100
SPECTRUM_E_BINS: np.ndarray = np.logspace(
    np.log10(1e-3), np.log10(2e7), SPECTRUM_N_GROUPS + 1
)

# ---------------------------------------------------------------------------
# Depletion constants
# ---------------------------------------------------------------------------

# 61.6 W/cm³ at the compact level — derived from AGR-1 (INL/EXT-10-19476):
# 186 W total compact power achieves 15% FIMA in 620 EFPD.
# This is the model power density (entire compact volume), not the kernel only.
POWER_DENSITY: float = 61.6  # W/cm³

# Time-step scheme — total 620.2 days ≈ 620 EFPD (AGR-1 irradiation duration).
# Finer steps during the first ~55 days capture short-lived fission product
# equilibria: Xe-135 (t½ = 9.17 h, equilibrium ~2 d) and Sm-149 (~53 h
# effective, equilibrium ~5–10 d). Uniform coarse steps for steady burnup.
TIMESTEPS: list[float] = (
    [1.0] * 5    # days 0–5:    Xe/Sm equilibration
    + [10.0] * 5   # days 5–55:   short-lived FP transient
    + [47.1] * 12  # days 55–620: steady-state burnup (12 × 47.1 d = 565.2 d)
)

# ---------------------------------------------------------------------------
# Transport settings for depletion steps
# ---------------------------------------------------------------------------
# Reduced vs. the standalone eigenvalue run (200 batches / 5 000 per batch).
# Depletion transport only needs to normalise power and compute one-group
# reaction rates; σ(k-eff) ≈ 30–50 pcm per step is sufficient.
# 2 000 particles × 70 active batches → 140 000 active histories / solve.
# CECMIntegrator: 2 solves/step × 22 steps = 44 solves × ~2.5 min ≈ 1.5–2 hr.
_DEP_BATCHES = 100
_DEP_INACTIVE = 30
_DEP_PARTICLES = 2_000


def build_depletion_model() -> openmc.Model:
    """Build the OpenMC Model for depletion transport steps.

    Identical geometry and materials to the eigenvalue model but with a
    reduced particle count appropriate for depletion normalization transport.
    """
    mats = build_materials()
    # CoupledOperator requires volume on every depletable material to convert
    # volumetric reaction rates (from tallies) into absolute atom counts for
    # the Bateman equations. Analytical estimate; matches the geometry packing.
    mats['kernel'].volume = KERNEL_VOLUME
    geom = build_geometry(mats)

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = _DEP_BATCHES
    settings.inactive = _DEP_INACTIVE
    settings.particles = _DEP_PARTICLES
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            [-_COMPACT_R, -_COMPACT_R, -_COMPACT_H / 2],
            [+_COMPACT_R, +_COMPACT_R, +_COMPACT_H / 2],
        )
    )

    # Fine-group spectrum tally in the kernel — written to every depletion
    # statepoint so that plot_spectrum.py can read it without re-running transport.
    t_spec = openmc.Tally(name='kernel spectrum')
    t_spec.filters = [
        openmc.MaterialFilter([mats['kernel']]),
        openmc.EnergyFilter(SPECTRUM_E_BINS),
    ]
    t_spec.scores = ['flux']

    return openmc.Model(
        geometry=geom,
        materials=openmc.Materials(list(mats.values())),
        settings=settings,
        tallies=openmc.Tallies([t_spec]),
    )


def run_depletion(
    chain_file: str | Path = 'data/chain_endfb71_thermal.xml',
) -> openmc.deplete.Results:
    """Run the burnup calculation and return depletion Results.

    Parameters
    ----------
    chain_file:
        Path to the depletion chain XML (ENDF/B-VII.1 thermal chain,
        downloaded by scripts/download_chain.sh).

    Returns
    -------
    openmc.deplete.Results
        Loaded results for k-eff and composition queries.

    Notes
    -----
    Writes ``depletion_results.h5`` and per-step ``statepoint.N.h5`` files
    to the current working directory.

    Normalization mode is ``'fission-q'``: normalises total power using
    per-nuclide fission Q-values from the chain file. Does not require kerma
    coefficients, so it is consistent across all transport libraries.
    # TODO: switch to 'energy-deposition' once photon transport is enabled
    # and heating-local tally is validated with ENDF/B-VIII.0 kerma data.
    """
    model = build_depletion_model()

    operator = openmc.deplete.CoupledOperator(
        model,
        chain_file=str(chain_file),
        normalization_mode='fission-q',
    )

    # CE/CM predictor-corrector: predicts the end-of-step composition using
    # beginning-of-step reaction rates (CE), then corrects with midpoint rates
    # (CM). Substantially more accurate than plain CE predictor at the cost of
    # 2 transport solves per step instead of 1.
    integrator = openmc.deplete.CECMIntegrator(
        operator,
        TIMESTEPS,
        power_density=POWER_DENSITY,
        timestep_units='d',
    )

    integrator.integrate()

    return openmc.deplete.Results('depletion_results.h5')


if __name__ == '__main__':
    import sys

    chain = sys.argv[1] if len(sys.argv) > 1 else 'data/chain_endfb71_thermal.xml'

    n_steps = len(TIMESTEPS)
    total_days = sum(TIMESTEPS)
    n_solves = 2 * n_steps  # CECMIntegrator: 2 solves per step
    print(f'Depletion schedule: {n_steps} steps, {total_days:.1f} days total')
    print(f'Power density: {POWER_DENSITY} W/cm³')
    print(f'Integrator: CECMIntegrator ({n_solves} transport solves)')
    print(f'Particles/batch: {_DEP_PARTICLES}, batches: {_DEP_BATCHES} '
          f'({_DEP_INACTIVE} inactive)')
    print(f'Chain file: {chain}')
    print()

    results = run_depletion(chain_file=chain)

    time, keff = results.get_keff()
    print(f'{"Step":>5}  {"Time [d]":>10}  {"k-eff":>10}  {"σ":>8}')
    print('-' * 40)
    for i, (t, k) in enumerate(zip(time / 86400, keff)):
        print(f'{i:>5}  {t:>10.1f}  {k[0]:>10.5f}  {k[1]:>8.5f}')
