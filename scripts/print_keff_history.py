#!/usr/bin/env python3
"""Print per-batch k-eff history from a completed eigenvalue statepoint.

Usage
-----
  python scripts/print_keff_history.py
  python scripts/print_keff_history.py --statepoint path/to/statepoint.200.h5
"""
import argparse
import warnings
from pathlib import Path

import numpy as np
import openmc

warnings.filterwarnings('ignore')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--statepoint', type=Path, default=Path('statepoint.200.h5'),
        help='Path to statepoint HDF5 (default: statepoint.200.h5)',
    )
    args = parser.parse_args()

    if not args.statepoint.exists():
        raise SystemExit(f'ERROR: statepoint not found: {args.statepoint}')

    with openmc.StatePoint(str(args.statepoint)) as sp:
        k_gen = sp.k_generation
        n_inactive = sp.n_inactive

        active_k = k_gen[n_inactive:]
        cumul_mean = np.cumsum(active_k) / np.arange(1, len(active_k) + 1)
        cumul_std = np.array([
            active_k[:i+1].std(ddof=1) / np.sqrt(i+1) if i > 0 else float('nan')
            for i in range(len(active_k))
        ])

        print(f'k-eff = {sp.keff.nominal_value:.5f} ± {sp.keff.std_dev:.5f}'
              f'  ({sp.n_batches} batches, {n_inactive} inactive, '
              f'{sp.n_particles} particles/batch)\n')

        print(f'{"Batch":>6}  {"Phase":>8}  {"k-gen":>10}  {"cumul mean":>12}  {"cumul ±1σ":>10}')
        print('-' * 58)
        for i, k in enumerate(k_gen, start=1):
            phase = 'inactive' if i <= n_inactive else 'active'
            if i > n_inactive:
                ai = i - n_inactive - 1
                c_std = cumul_std[ai]
                std_str = f'{c_std:>10.5f}' if not np.isnan(c_std) else f'{"—":>10}'
                print(f'{i:>6}  {phase:>8}  {k:>10.5f}  {cumul_mean[ai]:>12.5f}  {std_str}')
            else:
                print(f'{i:>6}  {phase:>8}  {k:>10.5f}  {"—":>12}  {"—":>10}')


if __name__ == '__main__':
    main()
