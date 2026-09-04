#!/usr/bin/env python3
"""Plot 2D cross-sections of the TRISO compact geometry.

Generates two plots:
  - XY slice (radial cross-section through the compact midplane)
  - XZ slice (axial cross-section through the compact centreline)

No simulation run required — only geometry and materials XML are written.

Usage
-----
  python scripts/plot_geometry.py
  python scripts/plot_geometry.py --output-dir output/
"""
import argparse
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'src'))

import openmc

# Geometry plotting doesn't run transport, but OpenMC reads cross_sections.xml
# when expanding enriched elements. Try both libraries; nndc is sufficient here.
_REPO_ROOT = Path(__file__).parents[1]
for _xs_candidate in (
    _REPO_ROOT / 'data' / 'nndc_hdf5' / 'cross_sections.xml',
    _REPO_ROOT / 'data' / 'endfb80_hdf5' / 'cross_sections.xml',
):
    if _xs_candidate.exists():
        openmc.config['cross_sections'] = str(_xs_candidate)
        break

from triso.materials import build_materials
from triso.geometry import build_geometry, _COMPACT_R, _COMPACT_H

# Material display colours (RGB 0–255)
_COLOURS = {
    'UCO kernel':           (255,  80,  80),   # red
    'porous carbon buffer': (200, 140,  60),   # tan/brown
    'inner PyC':            ( 80, 160, 220),   # light blue
    'SiC':                  ( 80, 200, 120),   # green
    'outer PyC':            ( 60, 100, 180),   # dark blue
    'graphite matrix':      (160, 160, 160),   # grey
}


def make_plots(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    mats = build_materials()
    geom = build_geometry(mats)
    model = openmc.Model(
        geometry=geom,
        materials=openmc.Materials(list(mats.values())),
    )

    colours = {mats[k]: v for k, v in {
        'kernel': _COLOURS['UCO kernel'],
        'buffer': _COLOURS['porous carbon buffer'],
        'ipyc':   _COLOURS['inner PyC'],
        'sic':    _COLOURS['SiC'],
        'opyc':   _COLOURS['outer PyC'],
        'matrix': _COLOURS['graphite matrix'],
    }.items()}

    plots = []

    # XY slice — radial view through compact midplane (z=0)
    xy = openmc.Plot(name='xy_midplane')
    xy.basis = 'xy'
    xy.origin = (0.0, 0.0, 0.0)
    xy.width = (2 * _COMPACT_R * 1.02, 2 * _COMPACT_R * 1.02)
    xy.pixels = (1200, 1200)
    xy.color_by = 'material'
    xy.colors = colours
    xy.filename = 'geometry_xy'
    plots.append(xy)

    # XZ slice — axial view through compact centreline (y=0)
    xz = openmc.Plot(name='xz_centreline')
    xz.basis = 'xz'
    xz.origin = (0.0, 0.0, 0.0)
    xz.width = (2 * _COMPACT_R * 1.02, _COMPACT_H * 1.02)
    xz.pixels = (600, int(600 * _COMPACT_H / (2 * _COMPACT_R)))
    xz.color_by = 'material'
    xz.colors = colours
    xz.filename = 'geometry_xz'
    plots.append(xz)

    model.plots = openmc.Plots(plots)

    # Write XMLs and run the plotter in a temp dir to avoid cluttering the cwd
    with tempfile.TemporaryDirectory() as tmp:
        model.export_to_xml(directory=tmp)
        openmc.plot_geometry(cwd=tmp, openmc_exec='/Users/nelsonpereira/mamba/envs/triso-env/bin/openmc')

        # Move PNGs to output_dir
        for name in ('geometry_xy.png', 'geometry_xz.png'):
            src = Path(tmp) / name
            if src.exists():
                dst = output_dir / name
                dst.write_bytes(src.read_bytes())
                print(f'Saved {dst}')
            else:
                print(f'WARNING: {name} not found in {tmp}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output-dir', type=Path, default=Path('output'),
        help='Directory for output PNGs (default: output/)',
    )
    args = parser.parse_args()
    make_plots(args.output_dir)


if __name__ == '__main__':
    main()
