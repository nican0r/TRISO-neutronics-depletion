"""TRISO compact eigenvalue run: settings, tallies, and result extraction.

Particle count / batch rationale:
  200 total batches, 50 inactive, 5 000 particles/batch
  → ~750 000 active particle-histories.
  Expected σ(k-eff) ≈ 10–20 pcm; key tally relative errors ≲ 5%.
  Chosen as the Stage 0 trade-off between runtime (≈ 5–15 min on a laptop)
  and sufficient tally convergence for spectral characterisation.

Energy group structure (3 groups, eV):
  Thermal     0 → 0.625 eV  (IAEA standard thermal cutoff)
  Epithermal  0.625 eV → 100 keV
  Fast        100 keV → 20 MeV
"""

from __future__ import annotations

from itertools import product as _iproduct
from pathlib import Path

import pandas as pd
import openmc
import openmc.model

from .geometry import _COMPACT_H, _COMPACT_R, build_geometry
from .materials import build_materials

_BATCHES = 200
_INACTIVE = 50
_PARTICLES = 5_000

# Coarse 3-group energy bin edges [eV] — thermal / epithermal / fast.
_ENERGY_BINS = [0.0, 0.625, 1.0e5, 2.0e7]


def _tally_to_df(tally: openmc.Tally) -> pd.DataFrame:
    """Extract tally results into a DataFrame without using get_pandas_dataframe.

    get_pandas_dataframe has a pandas 3.x incompatibility in OpenMC 0.15.x
    (assigns a tuple into a StringArray column → ValueError). This helper
    builds the same structure directly from tally.mean / tally.std_dev.
    """
    filter_names: list[str] = []
    filter_label_lists: list[list] = []
    for f in tally.filters:
        if isinstance(f, openmc.MaterialFilter):
            filter_names.append('material')
            # bins are Material objects when built in-memory, integer IDs when
            # read back from a statepoint — handle both.
            filter_label_lists.append([
                b.name if hasattr(b, 'name') else f'mat_{b}' for b in f.bins
            ])
        elif isinstance(f, openmc.EnergyFilter):
            filter_names.append('energy [eV]')
            bins = f.bins
            if bins.ndim == 2:
                # Statepoint format: (n_bins, 2) array of [low, high] pairs
                labels = [f'[{float(bins[i, 0]):.4g}, {float(bins[i, 1]):.4g})'
                          for i in range(len(bins))]
            else:
                # In-memory format: flat edge array [e0, e1, e2, ...]
                edges = [float(b) for b in bins]
                labels = [f'[{edges[i]:.4g}, {edges[i+1]:.4g})'
                          for i in range(len(edges) - 1)]
            filter_label_lists.append(labels)

    filter_combos = list(_iproduct(*filter_label_lists))
    nuclides = tally.nuclides if tally.nuclides else ['total']
    scores = tally.scores

    records = []
    for fi, combo in enumerate(filter_combos):
        for ni, nuc in enumerate(nuclides):
            for si, score in enumerate(scores):
                record = dict(zip(filter_names, combo))
                record['nuclide'] = nuc
                record['score'] = score
                record['mean'] = float(tally.mean[fi, ni, si])
                record['std_dev'] = float(tally.std_dev[fi, ni, si])
                records.append(record)

    return pd.DataFrame(records)


def build_model() -> openmc.Model:
    """Assemble the full OpenMC model: geometry, materials, settings, tallies."""
    mats = build_materials()
    geom = build_geometry(mats)

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = _BATCHES
    settings.inactive = _INACTIVE
    settings.particles = _PARTICLES
    # Uniform box source slightly larger than the compact cylinder; adequate
    # starting distribution for a reflective-boundary k-inf problem.
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            [-_COMPACT_R, -_COMPACT_R, -_COMPACT_H / 2],
            [+_COMPACT_R, +_COMPACT_R, +_COMPACT_H / 2],
        )
    )

    all_mats = list(mats.values())
    mat_filter = openmc.MaterialFilter(all_mats)
    energy_filter = openmc.EnergyFilter(_ENERGY_BINS)
    kernel_filter = openmc.MaterialFilter([mats['kernel']])

    # Energy-integrated flux aggregated by material layer type.
    # MaterialFilter tallies across *all* instances of each layer in the compact,
    # which is correct for a packed-sphere lattice where cells are repeated.
    t_flux = openmc.Tally(name='flux by material')
    t_flux.filters = [mat_filter]
    t_flux.scores = ['flux']

    # 3-group flux by material — reveals spectral shape per layer.
    t_flux_mg = openmc.Tally(name='flux by material 3-group')
    t_flux_mg.filters = [mat_filter, energy_filter]
    t_flux_mg.scores = ['flux']

    # Fission and absorption reaction rates in the UCO kernel.
    t_rxn = openmc.Tally(name='kernel reaction rates')
    t_rxn.filters = [kernel_filter]
    t_rxn.scores = ['fission', 'absorption']

    # Local energy deposition (heating-local excludes escaping photons;
    # valid without photon_transport=True and appropriate for Stage 0).
    # TODO: heating-local returns zero with the single-temperature NNDC library
    # because it requires kerma (kinetic energy release in matter) coefficients
    # which are not present in that library. Will produce non-zero results when
    # upgraded to a multi-temperature library that includes kerma data.
    t_heat = openmc.Tally(name='heating by material')
    t_heat.filters = [mat_filter]
    t_heat.scores = ['heating-local']

    return openmc.Model(
        geometry=geom,
        materials=openmc.Materials(all_mats),
        settings=settings,
        tallies=openmc.Tallies([t_flux, t_flux_mg, t_rxn, t_heat]),
    )


def run_model(model: openmc.Model | None = None, cwd: str | Path = '.') -> dict:
    """Run the eigenvalue calculation and return a results summary dict.

    Keys
    ----
    k_eff   : (mean, std_dev) tuple of the combined k-eff estimate
    flux    : DataFrame — energy-integrated flux by material
    flux_mg : DataFrame — 3-group flux by material
    rxn     : DataFrame — fission and absorption rates in the kernel
    heat    : DataFrame — heating-local by material
    """
    if model is None:
        model = build_model()

    sp_path = model.run(cwd=cwd)

    # openmc.StatePoint is file that gets written to at end of run
    with openmc.StatePoint(sp_path) as sp:
        k_mean = sp.keff.nominal_value
        k_std = sp.keff.std_dev
        flux_df = _tally_to_df(sp.get_tally(name='flux by material'))
        flux_mg_df = _tally_to_df(sp.get_tally(name='flux by material 3-group'))
        rxn_df = _tally_to_df(sp.get_tally(name='kernel reaction rates'))
        heat_df = _tally_to_df(sp.get_tally(name='heating by material'))

    return {
        'k_eff': (k_mean, k_std),
        'flux': flux_df,
        'flux_mg': flux_mg_df,
        'rxn': rxn_df,
        'heat': heat_df,
    }


if __name__ == '__main__':
    results = run_model()
    k, dk = results['k_eff']
    print(f'k-eff = {k:.5f} ± {dk:.5f}')
    print('\n--- Flux by material (energy-integrated) ---')
    print(results['flux'].to_string())
    print('\n--- Flux by material (3-group) ---')
    print(results['flux_mg'].to_string())
    print('\n--- Kernel reaction rates ---')
    print(results['rxn'].to_string())
    print('\n--- Heating by material ---')
    print(results['heat'].to_string())
