#!/usr/bin/env python3
"""Print 3-group flux by material from a completed eigenvalue statepoint.

Usage
-----
  python scripts/print_flux_3group.py
  python scripts/print_flux_3group.py --statepoint path/to/statepoint.200.h5
"""
import argparse
import sys
from pathlib import Path

import openmc

sys.path.insert(0, str(Path(__file__).parents[1] / 'src'))
from triso.run import _tally_to_df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--statepoint', type=Path, default=Path('statepoint.200.h5'),
        help='Path to statepoint HDF5 (default: statepoint.200.h5)',
    )
    args = parser.parse_args()

    if not args.statepoint.exists():
        raise SystemExit(f'ERROR: statepoint not found: {args.statepoint}')

    su = openmc.Summary(str(args.statepoint.parent / 'summary.h5'))
    id_to_name = {m.id: m.name for m in su.geometry.get_all_materials().values()}

    with openmc.StatePoint(str(args.statepoint)) as sp:
        df = _tally_to_df(sp.get_tally(name='flux by material 3-group'))

    df['material'] = (
        df['material'].str.extract(r'(\d+)')[0]
        .astype(int)
        .apply(lambda x: id_to_name.get(x, f'mat_{x}'))
    )

    groups = [
        ('[0, 0.625)',     'Thermal (0–0.625 eV)'),
        ('[0.625, 1e+05)', 'Epithermal (0.625 eV–100 keV)'),
        ('[1e+05, 2e+07)', 'Fast (100 keV–20 MeV)'),
    ]
    materials = [
        'UCO kernel', 'porous carbon buffer', 'inner PyC',
        'SiC', 'outer PyC', 'graphite matrix',
    ]

    col_mat = max(len(m) for m in materials)

    for g_key, g_label in groups:
        print(f'\n{g_label}')
        top = f'┌{"─"*(col_mat+2)}┬{"─"*11}┬{"─"*9}┐'
        hdr = f'│ {"Material":<{col_mat}} │ {"Mean flux":>9} │ {"±std":>7} │'
        div = f'├{"─"*(col_mat+2)}┼{"─"*11}┼{"─"*9}┤'
        bot = f'└{"─"*(col_mat+2)}┴{"─"*11}┴{"─"*9}┘'
        print(top)
        print(hdr)
        print(div)
        for i, mat in enumerate(materials):
            row = df[(df['material'] == mat) & (df['energy [eV]'] == g_key)]
            mean = row['mean'].values[0]
            std = row['std_dev'].values[0]
            print(f'│ {mat:<{col_mat}} │ {mean:>9.3f} │ {std:>7.5f} │')
            if i < len(materials) - 1:
                print(div)
        print(bot)


if __name__ == '__main__':
    main()
