#!/usr/bin/env python3
"""Plot k-inf vs. burnup from a completed TRISO compact depletion run.

Geometry uses all-reflective boundaries (k-inf, no neutron leakage), so the
plotted quantity is k-inf, not k-eff.  This is carried forward from Stage 0.

X-axis: EFPD (effective full-power days) with a secondary % FIMA axis derived
from the linear relationship FIMA% = 15% × (EFPD / 620.2).  This is exact for
constant-power irradiation and matches the AGR-1 target parameters used to set
the power density in src/triso/depletion.py.

Sanity checks performed
-----------------------
1. Monotonic decline overall  — k-inf must not increase between the first and
   last time point.  A monotonically increasing curve would indicate a
   problem with power normalisation or material composition.
2. No large upward discontinuity — no single step should show an increase
   larger than _STEP_JUMP_TOL (300 pcm); a larger jump suggests a corrupted
   step or insufficient particles for power normalization.
3. k-inf crosses unity — confirms meaningful burnup has occurred; failure
   means the compact never reached EOL during the irradiation schedule.

Usage
-----
  # From the project root:
  python scripts/plot_depletion.py

  # Specify a non-default results file:
  python scripts/plot_depletion.py --results output/depletion/depletion_results.h5
  python scripts/plot_depletion.py --output-dir plots/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import openmc.deplete

# ---------------------------------------------------------------------------
# Sanity-check thresholds
# ---------------------------------------------------------------------------
# Maximum single-step k-inf *increase* tolerated (pcm).  Statistical noise at
# ~200–250 pcm per-step (≈1σ) means occasional small rises are expected; 300
# pcm is a conservative 2σ ceiling.
_STEP_JUMP_TOL: int = 300  # pcm

# Target burnup and irradiation duration from AGR-1 parameters
# (INL/EXT-10-19476 / Demkowicz et al. 2018); defines the FIMA secondary axis.
_FIMA_TARGET_PCT: float = 15.0   # % FIMA at end of irradiation
_IRRAD_DAYS: float = 620.2       # total EFPD


def _fima_from_efpd(efpd: np.ndarray) -> np.ndarray:
    """Linear % FIMA approximation from EFPD at constant power."""
    return _FIMA_TARGET_PCT * efpd / _IRRAD_DAYS


def plot_keff_vs_burnup(
    results_path: Path,
    output_dir: Path,
) -> dict[str, bool]:
    """Load results, generate plot, run sanity checks.

    Parameters
    ----------
    results_path:
        Path to ``depletion_results.h5``.
    output_dir:
        Directory for the output PNG.

    Returns
    -------
    dict mapping check name → passed (bool).
    """
    results = openmc.deplete.Results(str(results_path))
    time_s, keff = results.get_keff()

    efpd = time_s / 86400.0
    k_mid = keff[:, 0]   # mean k-inf
    k_sig = keff[:, 1]   # 1σ uncertainty

    fima = _fima_from_efpd(efpd)

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.fill_between(efpd, k_mid - k_sig, k_mid + k_sig,
                    alpha=0.25, color='steelblue', label=r'±1$\sigma$')
    ax.plot(efpd, k_mid, 'o-', color='steelblue', markersize=4,
            linewidth=1.4, label='k-inf (mean)')
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.9, label='k-inf = 1')

    ax.set_xlabel('Irradiation time (EFPD)')
    ax.set_ylabel('k-inf')
    ax.set_title('TRISO compact k-inf vs. burnup\n'
                 'AGR-1 geometry, 19.75% HALEU, ENDF/B-VIII.0 @ 1200 K')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Secondary x-axis: % FIMA
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    efpd_ticks = np.linspace(0, _IRRAD_DAYS, 7)
    ax2.set_xticks(efpd_ticks)
    ax2.set_xticklabels([f'{_fima_from_efpd(t):.1f}%' for t in efpd_ticks])
    ax2.set_xlabel('Burnup (% FIMA)')

    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / 'keff_vs_burnup.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Plot saved to {out_path}')

    # -----------------------------------------------------------------------
    # Sanity checks
    # -----------------------------------------------------------------------
    checks: dict[str, bool] = {}

    # 1. Overall decline: k at end must be below k at start
    overall_decline = float(k_mid[-1]) < float(k_mid[0])
    checks['overall_decline'] = overall_decline

    # 2. No large upward jump between consecutive steps (> _STEP_JUMP_TOL pcm)
    step_deltas_pcm = np.diff(k_mid) * 1e5
    max_jump_pcm = float(np.max(step_deltas_pcm))
    no_large_jump = max_jump_pcm <= _STEP_JUMP_TOL
    checks['no_large_jump'] = no_large_jump

    # 3. k-inf crosses unity during the irradiation
    crosses_unity = bool(np.any(k_mid <= 1.0))
    checks['crosses_unity'] = crosses_unity

    return checks, k_mid, k_sig, efpd


def _print_table(k_mid, k_sig, efpd):
    fima = _fima_from_efpd(efpd)
    print(f'\n{"Step":>5}  {"EFPD":>8}  {"FIMA%":>6}  {"k-inf":>10}  {"σ":>8}')
    print('-' * 48)
    for i, (t, f, k, s) in enumerate(zip(efpd, fima, k_mid, k_sig)):
        print(f'{i:>5}  {t:>8.1f}  {f:>6.2f}  {k:>10.5f}  {s:>8.5f}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--results',
        type=Path,
        default=Path('output/depletion/depletion_results.h5'),
        help='Path to depletion_results.h5 (default: output/depletion/depletion_results.h5)',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('output'),
        help='Directory for keff_vs_burnup.png (default: output/)',
    )
    args = parser.parse_args()

    if not args.results.exists():
        sys.exit(f'ERROR: results file not found: {args.results}')

    checks, k_mid, k_sig, efpd = plot_keff_vs_burnup(args.results, args.output_dir)

    _print_table(k_mid, k_sig, efpd)

    print('\n--- Sanity checks ---')
    all_pass = True
    for name, passed in checks.items():
        status = 'PASS' if passed else 'FAIL'
        if not passed:
            all_pass = False
        print(f'  {name:<22} {status}')

    if not all_pass:
        print('\nOne or more sanity checks FAILED — investigate before proceeding.')
        sys.exit(1)
    else:
        print('\nAll sanity checks PASSED.')
        k_swing = (float(k_mid[0]) - float(k_mid[-1])) * 1e5
        fima_end = _fima_from_efpd(efpd[-1])
        print(f'  k-inf swing: {k_mid[0]:.5f} → {k_mid[-1]:.5f} '
              f'(Δ = {k_swing:.0f} pcm over {efpd[-1]:.1f} EFPD / '
              f'{fima_end:.1f}% FIMA)')


if __name__ == '__main__':
    main()
