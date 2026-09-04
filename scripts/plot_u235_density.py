#!/usr/bin/env python3
"""Side-by-side cross-section plot showing U-235 number density in TRISO kernels
at the start and end of the depletion irradiation.

All TRISO kernels share a single depletable material, so U-235 density is spatially
uniform within any given timestep. The comparison between the two panels makes the
integrated depletion visible.

Particle centres are regenerated via pack_spheres; the spatial layout is a fresh
random realisation (not the exact one used in the Monte Carlo run), which is
appropriate for a representational cross-section plot.

Usage
-----
  python scripts/plot_u235_density.py
  python scripts/plot_u235_density.py --results output/depletion/depletion_results.h5
  python scripts/plot_u235_density.py --output-dir output/
  python scripts/plot_u235_density.py --basis xz
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np

# Resolve the src/ package path regardless of working directory.
sys.path.insert(0, str(Path(__file__).parents[1] / 'src'))

import openmc
import openmc.deplete

for _xs in (
    Path(__file__).parents[1] / 'data' / 'nndc_hdf5' / 'cross_sections.xml',
    Path(__file__).parents[1] / 'data' / 'endfb80_hdf5' / 'cross_sections.xml',
):
    if _xs.exists():
        openmc.config['cross_sections'] = str(_xs)
        break

from triso.geometry import (
    _COMPACT_R, _COMPACT_H,
    _R_KERNEL, _R_BUFFER, _R_IPYC, _R_SIC, _R_OPYC,
    _PACKING_FRACTION,
)
from triso.depletion import KERNEL_VOLUME

# ---------------------------------------------------------------------------
# Layer colours (RGB 0–1 floats) matching plot_geometry.py conventions.
# Kernel colour is overridden by the U-235 density colormap.
# ---------------------------------------------------------------------------
_C_BUFFER = (200/255, 140/255,  60/255)   # tan/brown
_C_IPYC   = ( 80/255, 160/255, 220/255)   # light blue
_C_SIC    = ( 80/255, 200/255, 120/255)   # green
_C_OPYC   = ( 60/255, 100/255, 180/255)   # dark blue
_C_MATRIX = (160/255, 160/255, 160/255)   # grey

_LAYER_RADII = [_R_OPYC, _R_SIC, _R_IPYC, _R_BUFFER, _R_KERNEL]
_LAYER_COLORS = [_C_OPYC, _C_SIC, _C_IPYC, _C_BUFFER, None]  # None = colormap

# Thickness of the z-slice shown; particles with |cz| < _SLICE_HALF_WIDTH
# appear in the XY cross-section.
_SLICE_HALF_WIDTH = _R_OPYC * 1.5


def _get_particle_centers(basis: str) -> np.ndarray:
    """Return TRISO particle centres for a random compact realisation."""
    from triso.materials import build_materials
    from triso.geometry import build_geometry

    mats = build_materials()
    geom = build_geometry(mats)

    # Re-generate centres using the same geometry objects used in build_geometry.
    cyl = openmc.ZCylinder(r=_COMPACT_R, boundary_type='reflective')
    top = openmc.ZPlane(z0=+_COMPACT_H / 2, boundary_type='reflective')
    bot = openmc.ZPlane(z0=-_COMPACT_H / 2, boundary_type='reflective')
    region = -cyl & -top & +bot

    centers = openmc.model.pack_spheres(
        radius=_R_OPYC,
        region=region,
        pf=_PACKING_FRACTION,
    )
    return np.array(centers)


def _draw_panel(
    ax: plt.Axes,
    centers: np.ndarray,
    u235_density: float,
    kernel_cmap_norm: mcolors.Normalize,
    cmap,
    title: str,
    basis: str,
) -> None:
    """Draw one cross-section panel onto *ax*."""

    # Background (matrix)
    if basis == 'xy':
        background = plt.Circle((0, 0), _COMPACT_R, color=_C_MATRIX, zorder=0)
        ax.add_patch(background)
        ax.set_xlim(-_COMPACT_R * 1.05, _COMPACT_R * 1.05)
        ax.set_ylim(-_COMPACT_R * 1.05, _COMPACT_R * 1.05)
        ax.set_aspect('equal')
        ax.set_xlabel('x (cm)')
        ax.set_ylabel('y (cm)')
        # Compact boundary
        border = plt.Circle((0, 0), _COMPACT_R, fill=False,
                             edgecolor='black', linewidth=0.8, zorder=10)
        ax.add_patch(border)
        slice_axis = 2   # z is the slice axis for XY
        h_axis, v_axis = 0, 1
    else:  # xz
        rect = mpatches.FancyBboxPatch(
            (-_COMPACT_R, -_COMPACT_H / 2), 2 * _COMPACT_R, _COMPACT_H,
            boxstyle='square,pad=0', color=_C_MATRIX, zorder=0,
        )
        ax.add_patch(rect)
        ax.set_xlim(-_COMPACT_R * 1.05, _COMPACT_R * 1.05)
        ax.set_ylim(-_COMPACT_H / 2 * 1.05, _COMPACT_H / 2 * 1.05)
        ax.set_aspect('equal')
        ax.set_xlabel('x (cm)')
        ax.set_ylabel('z (cm)')
        # Compact boundary
        border = mpatches.FancyBboxPatch(
            (-_COMPACT_R, -_COMPACT_H / 2), 2 * _COMPACT_R, _COMPACT_H,
            boxstyle='square,pad=0', fill=False,
            edgecolor='black', linewidth=0.8, zorder=10,
        )
        ax.add_patch(border)
        slice_axis = 1   # y is the slice axis for XZ
        h_axis, v_axis = 0, 2

    kernel_color = cmap(kernel_cmap_norm(u235_density))

    # Sort so that particles closer to slice plane render on top.
    offsets = np.abs(centers[:, slice_axis])
    visible = offsets < _SLICE_HALF_WIDTH
    visible_idx = np.where(visible)[0]
    visible_idx = visible_idx[np.argsort(offsets[visible_idx])[::-1]]

    for idx in visible_idx:
        cx = centers[idx, h_axis]
        cy = centers[idx, v_axis]
        dz = centers[idx, slice_axis]

        # Draw layers from outside in so inner layers render on top.
        for r_layer, color in zip(_LAYER_RADII, _LAYER_COLORS):
            if r_layer <= abs(dz):
                continue  # layer fully above/below slice plane
            r_vis = math.sqrt(r_layer**2 - dz**2)
            c = color if color is not None else kernel_color
            circle = plt.Circle((cx, cy), r_vis, color=c, zorder=5)
            ax.add_patch(circle)

    ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=8)


def plot_u235_density(
    results_path: Path,
    output_dir: Path,
    basis: str = 'xy',
) -> None:
    """Generate side-by-side U-235 density cross-section plot."""

    # -----------------------------------------------------------------------
    # Load depletion data
    # -----------------------------------------------------------------------
    results = openmc.deplete.Results(str(results_path))
    time_s, _ = results.get_keff()
    time_d = time_s / 86400.0

    # Only one depletable material: kernel (mat ID '1')
    _, u235_atoms = results.get_atoms('1', 'U235')

    density_initial = float(u235_atoms[0])  / KERNEL_VOLUME   # atoms/cm³
    density_final   = float(u235_atoms[-1]) / KERNEL_VOLUME

    depletion_pct = (1 - density_final / density_initial) * 100

    print(f'U-235 number density — initial: {density_initial:.4e} atoms/cm³')
    print(f'U-235 number density — final:   {density_final:.4e} atoms/cm³')
    print(f'Depletion: {depletion_pct:.1f}% over {time_d[-1]:.1f} EFPD')

    # -----------------------------------------------------------------------
    # Geometry
    # -----------------------------------------------------------------------
    print('Generating particle centres via pack_spheres …')
    centers = _get_particle_centers(basis)
    print(f'  {len(centers)} particles packed')

    # -----------------------------------------------------------------------
    # Colormap and normalisation
    # -----------------------------------------------------------------------
    cmap = matplotlib.colormaps['YlOrRd']
    pad = (density_initial - density_final) * 0.10
    norm = mcolors.Normalize(
        vmin=density_final - pad,
        vmax=density_initial + pad,
    )

    # -----------------------------------------------------------------------
    # Figure
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(
        1, 2,
        figsize=(12, 5.5),
        gridspec_kw={'wspace': 0.35},
    )

    if basis == 'xy':
        slice_label = 'XY slice (z = 0, radial midplane)'
    else:
        slice_label = 'XZ slice (y = 0, axial centreline)'

    _draw_panel(
        axes[0], centers, density_initial, norm, cmap,
        f'BOL  (t = 0 EFPD)\n'
        f'U-235: {density_initial:.3e} at/cm³',
        basis,
    )
    _draw_panel(
        axes[1], centers, density_final, norm, cmap,
        f'EOL  (t = {time_d[-1]:.0f} EFPD)\n'
        f'U-235: {density_final:.3e} at/cm³',
        basis,
    )

    # Shared colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='vertical',
                        fraction=0.025, pad=0.02, aspect=30)
    cbar.set_label('U-235 number density (atoms/cm³)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Material legend
    legend_handles = [
        mpatches.Patch(color=_C_MATRIX, label='Graphite matrix'),
        mpatches.Patch(color=_C_OPYC,   label='OPyC'),
        mpatches.Patch(color=_C_SIC,    label='SiC'),
        mpatches.Patch(color=_C_IPYC,   label='IPyC'),
        mpatches.Patch(color=_C_BUFFER, label='Porous C buffer'),
        mpatches.Patch(color=cmap(norm(density_initial)), label='UCO kernel (BOL)'),
        mpatches.Patch(color=cmap(norm(density_final)),   label='UCO kernel (EOL)'),
    ]
    fig.legend(
        handles=legend_handles,
        loc='lower center',
        ncol=4,
        fontsize=8,
        framealpha=0.9,
        bbox_to_anchor=(0.45, -0.07),
    )

    fig.suptitle(
        f'TRISO compact U-235 density: BOL vs. EOL\n'
        f'AGR-1 geometry, 19.75% HALEU, {slice_label}',
        fontsize=11,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f'u235_density_{basis}.png'
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_path}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--results',
        type=Path,
        default=Path('output/depletion/depletion_results.h5'),
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('output'),
    )
    parser.add_argument(
        '--basis',
        choices=['xy', 'xz'],
        default='xy',
        help='Cross-section plane: xy (radial) or xz (axial). Default: xy',
    )
    args = parser.parse_args()

    if not args.results.exists():
        sys.exit(f'ERROR: results file not found: {args.results}')

    plot_u235_density(args.results, args.output_dir, args.basis)


if __name__ == '__main__':
    main()
