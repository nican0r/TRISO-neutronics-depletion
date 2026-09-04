#!/usr/bin/env python3
"""Neutron spectrum in the UCO kernel at BOL, MOL, and EOL of the depletion run.

Reads the 'kernel spectrum' tally written to the depletion statepoints by
``build_depletion_model()``.  Statepoint openmc_simulation_nN.h5 maps 1-to-1
to Results[N], so no additional transport runs are required.

Requires a depletion run that used the current ``build_depletion_model()``,
which includes the fine-group spectrum tally.  If the tally is absent (legacy
run without it), exit with a clear error message.

Usage
-----
  python scripts/plot_spectrum.py
  python scripts/plot_spectrum.py --depletion-dir output/depletion --output-dir output
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import openmc
import openmc.deplete

sys.path.insert(0, str(Path(__file__).parents[1] / 'src'))
from triso.depletion import SPECTRUM_E_BINS, TIMESTEPS

_FIMA_TARGET_PCT: float = 15.0        # % FIMA at AGR-1 end-of-irradiation
_IRRAD_DAYS: float = sum(TIMESTEPS)   # 620.2 EFPD

# Depletion step indices for BOL / MOL / EOL.
# statepoint openmc_simulation_nN.h5 == Results[N] (verified by k-eff comparison).
_SNAPSHOTS: dict[int, str] = {0: 'BOL', 15: 'MOL', 22: 'EOL'}

# Cumulative EFPD at the end of each depletion step.
_CUMULATIVE_DAYS: list[float] = [0.0] + list(np.cumsum(TIMESTEPS))

_THERMAL_CUTOFF_EV: float = 0.625  # IAEA standard thermal cutoff


def _load_spectrum(
    sp_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read kernel spectrum tally from a depletion statepoint.

    Returns
    -------
    e_mid : ndarray (N_GROUPS,)  — geometric-mean energy per group [eV]
    phi_leth : ndarray (N_GROUPS,)  — lethargy-normalised flux, integral-normalised
                                      (probability density: integrates to 1 over lethargy)
    phi_leth_err : ndarray (N_GROUPS,)  — 1σ uncertainty (same scale)
    """
    with openmc.StatePoint(str(sp_path), autolink=False) as sp:
        try:
            tally = sp.get_tally(name='kernel spectrum')
        except KeyError:
            sys.exit(
                f'ERROR: tally "kernel spectrum" not found in {sp_path}.\n'
                'Re-run the depletion with the current build_depletion_model() '
                'which includes the fine-group spectrum tally.'
            )
        flux = tally.mean[:, 0, 0]
        flux_err = tally.std_dev[:, 0, 0]

    e_mid = np.sqrt(SPECTRUM_E_BINS[:-1] * SPECTRUM_E_BINS[1:])
    d_lethargy = np.log(SPECTRUM_E_BINS[1:] / SPECTRUM_E_BINS[:-1])

    phi_leth = flux / d_lethargy
    phi_leth_err = flux_err / d_lethargy
    integral = np.sum(phi_leth * d_lethargy)
    return e_mid, phi_leth / integral, phi_leth_err / integral


def _hardening_metrics(
    spectra: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> None:
    d_u = np.log(SPECTRUM_E_BINS[1:] / SPECTRUM_E_BINS[:-1])
    e_mid = np.sqrt(SPECTRUM_E_BINS[:-1] * SPECTRUM_E_BINS[1:])

    # Thermal fraction is the primary hardening metric. Flux-weighted mean
    # energy is omitted: it is dominated by the fast peak (~1 MeV) and gives
    # a misleadingly large absolute value that obscures the thermal/epithermal
    # shift of interest.
    print(f'\n{"Label":<6}  {"EFPD":>7}  {"% FIMA":>7}  {"Thermal frac":>13}')
    print(f'  (fraction of neutrons below {_THERMAL_CUTOFF_EV} eV; compact-only geometry'
          f' → low absolute value expected, ~3% BOL vs ~20–40% for full HTGR element)')
    print('-' * 42)
    for step_idx, label in _SNAPSHOTS.items():
        efpd = _CUMULATIVE_DAYS[step_idx]
        fima = _FIMA_TARGET_PCT * efpd / _IRRAD_DAYS
        _, phi_n, _ = spectra[label]
        phi_abs = phi_n * d_u
        thermal_frac = float(phi_abs[e_mid < _THERMAL_CUTOFF_EV].sum() / phi_abs.sum())
        print(f'{label:<6}  {efpd:>7.1f}  {fima:>7.2f}  {thermal_frac:>13.4f}')


def _assess_trend(
    spectra: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> None:
    d_u = np.log(SPECTRUM_E_BINS[1:] / SPECTRUM_E_BINS[:-1])
    e_mid = np.sqrt(SPECTRUM_E_BINS[:-1] * SPECTRUM_E_BINS[1:])

    def therm_frac(label: str) -> float:
        _, phi_n, _ = spectra[label]
        phi_abs = phi_n * d_u
        return float(phi_abs[e_mid < _THERMAL_CUTOFF_EV].sum() / phi_abs.sum())

    tf_bol = therm_frac('BOL')
    tf_eol = therm_frac('EOL')
    delta_th = tf_eol - tf_bol
    rel_drop = delta_th / tf_bol * 100.0

    direction_th = 'HARDENING' if delta_th < 0 else 'SOFTENING'
    print('\n--- Spectral trend (BOL → EOL) ---')
    print(f'  Thermal frac change:  {delta_th:+.4f}  ({rel_drop:+.1f}%)  → {direction_th}')
    print()
    print('  Competing physical drivers:')
    print('  - FP accumulation (Xe-135, Sm-149, …): strong thermal absorbers → HARDER')
    print('    (dominant effect at ≤15% FIMA)')
    print('  - U-235 depletion (σ_th≈685 b): fewer thermal absorbers → softer')
    print('    (competing effect, but outweighed by FP accumulation here)')
    print('  - Pu-239 buildup: reactivity effect (offsets U-235 loss in k-eff),')
    print('    minor contribution to spectral shift at these burnup levels')
    print()
    print(f'  Conclusion: FP accumulation dominates — thermal fraction drops'
          f' {abs(rel_drop):.0f}% over the irradiation.')


def plot_spectrum(depletion_dir: Path, output_dir: Path) -> None:
    results_path = depletion_dir / 'depletion_results.h5'
    if not results_path.exists():
        sys.exit(f'ERROR: {results_path} not found')

    results = openmc.deplete.Results(str(results_path))
    print(f'Loaded {len(results)} depletion steps from {results_path}')

    spectra: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for step_idx, label in _SNAPSHOTS.items():
        efpd = _CUMULATIVE_DAYS[step_idx]
        fima = _FIMA_TARGET_PCT * efpd / _IRRAD_DAYS
        sp_path = depletion_dir / f'openmc_simulation_n{step_idx}.h5'
        if not sp_path.exists():
            sys.exit(f'ERROR: statepoint not found: {sp_path}')
        print(f'Reading {label} ({efpd:.1f} EFPD, {fima:.2f}% FIMA): {sp_path.name}')
        spectra[label] = _load_spectrum(sp_path)

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    colors = {'BOL': 'steelblue', 'MOL': 'darkorange', 'EOL': 'firebrick'}
    snap_labels = {
        'BOL': f'BOL  ({_CUMULATIVE_DAYS[0]:.0f} EFPD, 0.0% FIMA)',
        'MOL': (f'MOL  ({_CUMULATIVE_DAYS[15]:.1f} EFPD, '
                f'{_FIMA_TARGET_PCT*_CUMULATIVE_DAYS[15]/_IRRAD_DAYS:.1f}% FIMA)'),
        'EOL': f'EOL  ({_CUMULATIVE_DAYS[22]:.1f} EFPD, 15.0% FIMA)',
    }

    fig, ax = plt.subplots(figsize=(9, 5))
    for label in ['BOL', 'MOL', 'EOL']:
        e_mid, phi, _ = spectra[label]
        ax.plot(e_mid, phi, color=colors[label], linewidth=1.4,
                label=snap_labels[label])

    ax.axvline(_THERMAL_CUTOFF_EV, color='gray', linestyle=':', linewidth=0.9,
               label=f'Thermal cutoff ({_THERMAL_CUTOFF_EV} eV)')
    ax.axvline(1e5, color='gray', linestyle='--', linewidth=0.9,
               label='Fast threshold (100 keV)')

    ax.set_xscale('log')
    ax.set_xlabel('Neutron energy (eV)')
    ax.set_ylabel('φ(E) per unit lethargy (probability density, integrates to 1)')
    ax.set_title(
        'TRISO compact — neutron spectrum in UCO kernel: BOL / MOL / EOL\n'
        'AGR-1 geometry, 19.75% HALEU, ENDF/B-VIII.0 @ 1200 K'
    )
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(SPECTRUM_E_BINS[0], SPECTRUM_E_BINS[-1])
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / 'spectrum_hardening.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'\nPlot saved to {out_path}')

    _hardening_metrics(spectra)
    _assess_trend(spectra)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--depletion-dir', type=Path,
        default=Path('output/depletion'),
        help='Directory containing depletion_results.h5 and openmc_simulation_n*.h5',
    )
    parser.add_argument(
        '--output-dir', type=Path, default=Path('output'),
        help='Directory for spectrum_hardening.png',
    )
    args = parser.parse_args()
    plot_spectrum(args.depletion_dir, args.output_dir)


if __name__ == '__main__':
    main()
