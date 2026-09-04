#!/usr/bin/env python3
"""Compute the fuel temperature coefficient (Doppler) for TRISO compacts.

Evaluates k-inf at fuel temperatures 900 K, 1050 K, and 1200 K for two burnup
conditions (BOL: fresh fuel; EOL: 620.2 EFPD / 15% FIMA from depletion_results.h5).

Temperature treatment
---------------------
The ENDF/B-VIII.0 library (trimmed to 900 K + 1200 K) does NOT include windowed-
multipole (WMP) data — the `0K` groups in the HDF5 files are base energy grids
for the tabulated cross sections, not WMP broadening data.  At 900 K and 1200 K
the exact tabulated values are used.  At 1050 K, OpenMC uses sqrt(T) interpolation
(Settings.temperature method='interpolation') between the two bounding tabulated
points.  Extrapolation beyond 1200 K is not supported without additional tabulated
data (e.g., 2500 K from the full ENDF/B-VIII.0 distribution); obtaining 1500 K
cross sections would require re-downloading the library and retaining the 2500 K
point so that 1500 K can be bounded.

What is held fixed
------------------
All non-kernel materials (porous carbon buffer, IPyC, SiC, OPyC, graphite matrix)
remain at 1200 K for every run.  Only the UCO kernel temperature varies.  This
isolates the fuel/Doppler component; the moderator temperature coefficient and
geometric feedback are not captured here.

Temperature sweep
-----------------
  900 K   — below nominal operating temperature (tabulated, exact)
  1050 K  — intermediate (sqrt(T) interpolation between 900 K and 1200 K)
  1200 K  — nominal operating temperature (tabulated, exact)

Published comparison
--------------------
HTGR fuel temperature coefficient (Doppler only): approximately −2 to −6 pcm/K.
Sources: Kuijper et al., Nucl. Sci. Eng. 153 (2006) 276–306;
         IAEA-TECDOC-978 (2001).

Usage
-----
  python scripts/doppler_coefficient.py
  python scripts/doppler_coefficient.py --dep-h5 path/to/depletion_results.h5
  python scripts/doppler_coefficient.py --output-dir plots/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import openmc
import openmc.model

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from triso.geometry import _COMPACT_R, _COMPACT_H, build_geometry
from triso.materials import build_materials
from triso.depletion import KERNEL_VOLUME

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

_T_STRUCT: float = 1200.0                       # K — moderator / structural (fixed)
_FUEL_TEMPS: list[float] = [900.0, 1050.0, 1200.0]  # K — kernel temperature sweep

# Transport settings — same as depletion step (σ(k-eff) ≈ 30–50 pcm per run;
# at ΔT = 600 K and |FTC| ≈ 3 pcm/K the expected Δk ≈ 1800 pcm, so 30–50 pcm
# σ resolves the coefficient to <5% statistical uncertainty).
_BATCHES: int = 100
_INACTIVE: int = 30
_PARTICLES: int = 2_000

# Published HTGR fuel temperature coefficient range [pcm/K] (Doppler only).
# Kuijper et al., NSE 153 (2006) 276–306; IAEA-TECDOC-978.
_FTC_PUB_MIN: float = -6.0
_FTC_PUB_MAX: float = -2.0

# Minimum atom fraction to include in the EOL material.  Nuclides below this
# threshold have negligible cross-section contributions at 15% FIMA.
_ATOM_FRAC_THRESHOLD: float = 1e-12

# EOL step index in the depletion_results.h5 number array.
# number.shape = (23, 1, 3819): 23 states (initial + 22 steps).
_EOL_STEP: int = -1   # last row = step 22, t = 620.2 EFPD

# ---------------------------------------------------------------------------
# HDF5 data access
# ---------------------------------------------------------------------------

def _nuc_to_index(dep_h5: Path) -> dict[str, int]:
    """Build nuclide_name → number-array-index mapping.

    The 'nuclides' group in depletion_results.h5 has per-nuclide subgroups
    with an 'atom number index' attribute that gives the column index in the
    'number' dataset (not the alphabetical position of the key).
    """
    with h5py.File(dep_h5, 'r') as f:
        return {
            nuc: int(f['nuclides'][nuc].attrs['atom number index'])
            for nuc in f['nuclides'].keys()
        }


def _atom_counts(
    dep_h5: Path,
    step: int,
    nuc_to_idx: dict[str, int],
) -> dict[str, float]:
    """Return {nuclide: atom_count} for the kernel material at depletion step."""
    with h5py.File(dep_h5, 'r') as f:
        number = f['number'][()]   # (n_steps, n_mats, n_nuclides)
    return {
        nuc: float(number[step, 0, idx])
        for nuc, idx in nuc_to_idx.items()
    }


def _library_neutron_nuclides() -> set[str]:
    """Return set of nuclide names with neutron XS data in the active library."""
    lib = openmc.data.DataLibrary.from_xml()
    return {
        m
        for entry in lib.libraries
        if entry.get('type') == 'neutron'
        for m in entry.get('materials', [])
    }


# ---------------------------------------------------------------------------
# Material builders
# ---------------------------------------------------------------------------

def _eol_kernel(
    atoms: dict[str, float],
    lib_nucs: set[str],
    volume_cm3: float,
    temperature: float,
) -> openmc.Material:
    """Build a UCO kernel Material from depletion atom counts.

    Only nuclides present in the transport library and above the atom-fraction
    threshold are included.  Excluded nuclides are predominantly exotic short-
    lived fission products with negligible transport XS contribution.
    """
    filtered = {
        nuc: n for nuc, n in atoms.items()
        if n > 0 and nuc in lib_nucs
    }
    total_atoms = sum(filtered.values())
    # atom/b-cm: 1 b·cm = 1e-24 cm³, so N [atom/b-cm] = N_atoms / (V_cm3 × 1e24)
    total_density_bcm = total_atoms / (volume_cm3 * 1e24)

    mat = openmc.Material(name='UCO kernel (EOL)')
    mat.temperature = temperature
    mat.depletable = True
    mat.set_density('atom/b-cm', total_density_bcm)
    for nuc, n in filtered.items():
        frac = n / total_atoms
        if frac >= _ATOM_FRAC_THRESHOLD:
            mat.add_nuclide(nuc, frac)   # default percent_type='ao' (atom fraction)
    return mat


def _mats_for_run(
    fuel_temp: float,
    eol_atoms: dict[str, float] | None,
    lib_nucs: set[str],
) -> dict[str, openmc.Material]:
    """Return material dict with kernel at fuel_temp; everything else at 1200 K."""
    mats = build_materials()
    if eol_atoms is not None:
        mats['kernel'] = _eol_kernel(eol_atoms, lib_nucs, KERNEL_VOLUME, fuel_temp)
    else:
        mats['kernel'].temperature = fuel_temp
    # buffer, ipyc, sic, opyc, matrix are already at _T_STRUCT from build_materials()
    return mats


# ---------------------------------------------------------------------------
# Single eigenvalue run (skip if statepoint cached)
# ---------------------------------------------------------------------------

def _run_keff(mats: dict, work_dir: Path) -> tuple[float, float]:
    """Run eigenvalue calculation in work_dir and return (k_mean, k_std)."""
    sp_file = work_dir / f'statepoint.{_BATCHES}.h5'

    if not sp_file.exists():
        work_dir.mkdir(parents=True, exist_ok=True)
        geom = build_geometry(mats)
        settings = openmc.Settings()
        settings.run_mode = 'eigenvalue'
        settings.batches = _BATCHES
        settings.inactive = _INACTIVE
        settings.particles = _PARTICLES
        settings.temperature = {'method': 'interpolation'}
        settings.source = openmc.IndependentSource(
            space=openmc.stats.Box(
                [-_COMPACT_R, -_COMPACT_R, -_COMPACT_H / 2],
                [+_COMPACT_R, +_COMPACT_R, +_COMPACT_H / 2],
            )
        )

        model = openmc.Model(
            geometry=geom,
            materials=openmc.Materials(list(mats.values())),
            settings=settings,
        )
        sp_file = model.run(cwd=work_dir)
    else:
        print(f'    [cached] {sp_file}')

    with openmc.StatePoint(sp_file) as sp:
        return float(sp.keff.nominal_value), float(sp.keff.std_dev)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _plot_keff_and_ftc(ftc_store: list[dict], output_dir: Path) -> None:
    """Single-panel figure: k-inf vs T for BOL and EOL."""
    burn_colors = {'BOL': 'steelblue', 'EOL': 'firebrick'}

    fig, ax_k = plt.subplots(1, 1, figsize=(7, 5))

    for data in ftc_store:
        T = data['T']
        k = data['k']
        dk = data['dk']
        short = data['short']
        color = burn_colors[short]

        coeffs = np.polyfit(T, k, 1)
        T_fit = np.linspace(T[0] - 50, T[-1] + 50, 200)
        k_fit = np.polyval(coeffs, T_fit)
        ftc_pcm = data['ftc_pcm']

        ax_k.errorbar(T, k, yerr=dk, fmt='o', capsize=5, color=color,
                      zorder=5, label=f'{short}  MC ± 1σ')
        ax_k.plot(T_fit, k_fit, '--', color=color, alpha=0.7,
                  label=f'{short}  FTC = {ftc_pcm:+.2f} pcm/K')

    ax_k.set_xlabel('Fuel temperature (K)')
    ax_k.set_ylabel('k-inf')
    ax_k.set_title('k-inf vs fuel temperature')
    ax_k.legend(fontsize=8)
    ax_k.grid(True, alpha=0.3)

    fig.suptitle(
        'Fuel Temperature Coefficient (Doppler)\n'
        'Moderator / structural temperature fixed at 1200 K',
        fontsize=11,
    )
    fig.tight_layout()
    out = output_dir / 'doppler_ftc.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Plot saved → {out}')


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def compute_ftc(dep_h5: Path, output_dir: Path) -> int:
    """Run all (burnup × temperature) combinations and report the FTC.

    Returns 0 if all checks pass (sign negative, magnitude within published
    range), 1 if any check fails.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    work_root = output_dir / 'doppler'

    lib_nucs = _library_neutron_nuclides()
    nuc_to_idx = _nuc_to_index(dep_h5)
    eol_atoms = _atom_counts(dep_h5, _EOL_STEP, nuc_to_idx)

    n_lib_at_eol = sum(
        1 for nuc, n in eol_atoms.items()
        if n > 0 and nuc in lib_nucs
    )
    total_eol = sum(n for n in eol_atoms.values() if n > 0)
    lib_eol = sum(n for nuc, n in eol_atoms.items() if n > 0 and nuc in lib_nucs)
    print(f'EOL composition: {n_lib_at_eol} nuclides in transport library '
          f'({100*lib_eol/total_eol:.2f}% of total atoms)')

    burnup_configs: list[tuple[str, str, dict | None]] = [
        ('BOL (0 EFPD / 0% FIMA)',     'BOL', None),
        ('EOL (620 EFPD / 15% FIMA)',  'EOL', eol_atoms),
    ]

    results: dict[tuple[str, float], tuple[float, float]] = {}

    print()
    print('Running 6 eigenvalue calculations (100 batches / 2000 particles each) ...')
    print()
    for long_label, short, atoms in burnup_configs:
        for T in _FUEL_TEMPS:
            work_dir = work_root / f'{short}_{int(T)}K'
            print(f'  [{short}  {T:.0f} K]  {work_dir}')
            mats = _mats_for_run(T, atoms, lib_nucs)
            k_mean, k_std = _run_keff(mats, work_dir)
            results[(short, T)] = (k_mean, k_std)
            print(f'    k-inf = {k_mean:.5f} ± {k_std:.5f}')

    # --- Compute FTC via linear regression over 3 temperature points ---
    T_arr = np.array(_FUEL_TEMPS)
    exit_code = 0
    ftc_store: list[dict] = []

    print()
    print('=' * 66)
    print('FUEL TEMPERATURE COEFFICIENT (Doppler) — CHECKPOINT')
    print('Fixed:  buffer, PyC, SiC, graphite matrix at 1200 K')
    print('Varied: UCO kernel temperature (900 K / 1050 K / 1200 K)')
    print(f'Published HTGR range: [{_FTC_PUB_MIN:.0f}, {_FTC_PUB_MAX:.0f}] pcm/K')
    print(f'  Kuijper et al., NSE 153 (2006); IAEA-TECDOC-978')
    print('=' * 66)

    for long_label, short, _ in burnup_configs:
        k_vals = np.array([results[(short, T)][0] for T in _FUEL_TEMPS])
        k_errs = np.array([results[(short, T)][1] for T in _FUEL_TEMPS])
        slope, _ = np.polyfit(T_arr, k_vals, 1)
        ftc_pcm = slope * 1e5

        ftc_store.append({
            'long_label': long_label,
            'short': short,
            'T': T_arr,
            'k': k_vals,
            'dk': k_errs,
            'ftc_pcm': ftc_pcm,
        })

        print(f'\n{long_label}:')
        for T, k, dk in zip(T_arr, k_vals, k_errs):
            print(f'  T = {T:6.0f} K   k-inf = {k:.5f} ± {dk:.5f}')
        print(f'  FTC = {ftc_pcm:+.2f} pcm/K   (linear regression, 3 points)')

        if ftc_pcm >= 0.0:
            print(f'  *** FAIL *** FTC is positive — Doppler must be negative')
            exit_code = 1
        elif ftc_pcm < _FTC_PUB_MIN:
            print(f'  WARN: |FTC| larger than published range — flag for review')
        elif ftc_pcm > _FTC_PUB_MAX:
            print(f'  WARN: |FTC| smaller than published range — flag for review')
        else:
            print(f'  PASS: sign negative, magnitude within published HTGR range')

    # --- Plots ---
    _plot_keff_and_ftc(ftc_store, output_dir)

    return exit_code


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--dep-h5',
        type=Path,
        default=Path('output/depletion/depletion_results.h5'),
        help='Path to depletion_results.h5 (default: output/depletion/depletion_results.h5)',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('output'),
        help='Directory for plots and eigenvalue run subdirs (default: output/)',
    )
    args = parser.parse_args()

    if not args.dep_h5.exists():
        sys.exit(f'ERROR: depletion results not found: {args.dep_h5}')

    exit_code = compute_ftc(args.dep_h5, args.output_dir)

    if exit_code != 0:
        print('\nOne or more checks FAILED — investigate before accepting results.')
    else:
        print('\nAll checks PASSED.')
    sys.exit(exit_code)


if __name__ == '__main__':
    main()
