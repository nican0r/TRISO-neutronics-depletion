#!/usr/bin/env python3
"""Manufacturing tolerance sensitivity — TRISO kernel diameter and packing fraction.

Runs five BOL eigenvalue calculations (nominal + four perturbations) and reports
Δρ in pcm (reactivity difference) to identify which manufacturing parameter has
the larger reactivity effect.  k-eff is dimensionless; Δρ = (k_pert − k_nom) /
(k_pert × k_nom) converts the k-eff shift to reactivity before scaling to pcm.

Baseline
--------
Nominal Step-6 / Step-3 geometry: packing fraction 30%, kernel radius 0.0175 cm
(350 µm diameter, AGR-1 nominal per INL/EXT-10-19476 Table 3).

Perturbations
-------------
  PF − 2 pp :  packing fraction 28%, kernel unchanged
  PF + 2 pp :  packing fraction 32%, kernel unchanged
  Kernel − 5 µm:  kernel radius 0.0170 cm (340 µm ⌀), PF unchanged
  Kernel + 5 µm:  kernel radius 0.0180 cm (360 µm ⌀), PF unchanged

Kernel diameter tolerance of ±10 µm (±5 µm on radius) is cited from
INL/EXT-10-19476 Table 3.  Only the kernel sphere radius is perturbed; all
coating-layer outer radii (buffer 0.0275 cm, IPyC 0.0315 cm, SiC 0.0350 cm,
OPyC 0.0390 cm) are held at their nominal AGR-1 values.  The buffer layer
absorbs the kernel-size delta; the OPyC outer radius is unchanged, so the
packing geometry is unaffected.

Transport settings: 200 batches (50 inactive), 5 000 particles/batch — identical
to the nominal Step-3 eigenvalue run (expected σ(k-eff) ≈ 1–2 × 10⁻⁴, i.e. ~10–20 pcm
in reactivity).

Usage
-----
  caffeinate python scripts/manufacturing_tolerance.py
  caffeinate python scripts/manufacturing_tolerance.py --output-dir plots/
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import openmc
import openmc.model

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from triso.geometry import (
    _COMPACT_R, _COMPACT_H,
    _R_BUFFER, _R_IPYC, _R_SIC, _R_OPYC,
)
from triso.materials import build_materials

# ---------------------------------------------------------------------------
# Nominal geometry constants (Step-6 baseline)
# ---------------------------------------------------------------------------

_PF_NOM: float = 0.30
_R_KERNEL_NOM: float = 0.0175   # cm — 350 µm diameter; AGR-1 nominal, INL/EXT-10-19476 Table 3
_LATTICE_PITCH: float = 2.0 * _R_OPYC   # 0.078 cm; same rule as geometry.py

# Perturbation magnitudes
_PF_DELTA: float = 0.02          # ±2 percentage points
_R_KERNEL_DELTA: float = 0.0005  # ±5 µm radius (= ±10 µm diameter); INL/EXT-10-19476 Table 3

# Transport settings — identical to Step-3 nominal run.
_BATCHES: int = 200
_INACTIVE: int = 50
_PARTICLES: int = 5_000

# ---------------------------------------------------------------------------
# Case definitions
# ---------------------------------------------------------------------------

_CASES: list[dict] = [
    {
        'label': 'Nominal',
        'short': 'nominal',
        'pf': _PF_NOM,
        'r_kernel': _R_KERNEL_NOM,
    },
    {
        'label': 'PF − 2 pp (28%)',
        'short': 'pf_minus2',
        'pf': _PF_NOM - _PF_DELTA,
        'r_kernel': _R_KERNEL_NOM,
    },
    {
        'label': 'PF + 2 pp (32%)',
        'short': 'pf_plus2',
        'pf': _PF_NOM + _PF_DELTA,
        'r_kernel': _R_KERNEL_NOM,
    },
    {
        'label': 'Kernel − 5 µm (340 µm ⌀)',
        'short': 'kern_minus',
        'pf': _PF_NOM,
        'r_kernel': _R_KERNEL_NOM - _R_KERNEL_DELTA,
    },
    {
        'label': 'Kernel + 5 µm (360 µm ⌀)',
        'short': 'kern_plus',
        'pf': _PF_NOM,
        'r_kernel': _R_KERNEL_NOM + _R_KERNEL_DELTA,
    },
]

# ---------------------------------------------------------------------------
# Geometry builder
# ---------------------------------------------------------------------------

def _build_geometry(
    mats: dict,
    pf: float = _PF_NOM,
    r_kernel: float = _R_KERNEL_NOM,
) -> openmc.Geometry:
    """Build a reflective TRISO compact with perturbed packing fraction or kernel radius.

    Only the kernel sphere radius varies; all coating-layer outer radii are fixed
    at their nominal AGR-1 values so that the overall particle outer radius (and
    therefore the lattice pitch and packing geometry) is unchanged.
    """
    s1 = openmc.Sphere(r=r_kernel)
    s2 = openmc.Sphere(r=_R_BUFFER)
    s3 = openmc.Sphere(r=_R_IPYC)
    s4 = openmc.Sphere(r=_R_SIC)
    s5 = openmc.Sphere(r=_R_OPYC)
    triso_univ = openmc.Universe(cells=[
        openmc.Cell(fill=mats['kernel'], region=-s1),
        openmc.Cell(fill=mats['buffer'], region=+s1 & -s2),
        openmc.Cell(fill=mats['ipyc'],   region=+s2 & -s3),
        openmc.Cell(fill=mats['sic'],    region=+s3 & -s4),
        openmc.Cell(fill=mats['opyc'],   region=+s4 & -s5),
    ])

    cyl = openmc.ZCylinder(r=_COMPACT_R, boundary_type='reflective')
    top = openmc.ZPlane(z0=+_COMPACT_H / 2, boundary_type='reflective')
    bot = openmc.ZPlane(z0=-_COMPACT_H / 2, boundary_type='reflective')
    compact_region = -cyl & -top & +bot

    centers = openmc.model.pack_spheres(
        radius=_R_OPYC,
        region=compact_region,
        pf=pf,
    )
    trisos = [
        openmc.model.TRISO(outer_radius=_R_OPYC, fill=triso_univ, center=c)
        for c in centers
    ]

    p = _LATTICE_PITCH
    nx = math.ceil(2 * _COMPACT_R / p)
    ny = nx
    nz = math.ceil(_COMPACT_H / p)

    triso_lat = openmc.model.create_triso_lattice(
        trisos,
        lower_left=(-_COMPACT_R, -_COMPACT_R, -_COMPACT_H / 2),
        pitch=(p, p, p),
        shape=(nx, ny, nz),
        background=mats['matrix'],
    )
    compact_cell = openmc.Cell(fill=triso_lat, region=compact_region)
    return openmc.Geometry(openmc.Universe(cells=[compact_cell]))


# ---------------------------------------------------------------------------
# Eigenvalue runner
# ---------------------------------------------------------------------------

def _run_keff(
    mats: dict,
    geom: openmc.Geometry,
    work_dir: Path,
) -> tuple[float, float]:
    """Run eigenvalue in work_dir (or load cached statepoint) and return (k, σ)."""
    sp_file = work_dir / f'statepoint.{_BATCHES}.h5'
    if not sp_file.exists():
        work_dir.mkdir(parents=True, exist_ok=True)
        settings = openmc.Settings()
        settings.run_mode = 'eigenvalue'
        settings.batches = _BATCHES
        settings.inactive = _INACTIVE
        settings.particles = _PARTICLES
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
# Plotting
# ---------------------------------------------------------------------------

def _plot_sensitivity(results: list[dict], output_dir: Path) -> None:
    """Bar chart of Δρ (pcm) for each perturbed case vs. nominal."""
    k_nom = results[0]['k']
    dk_nom = results[0]['dk']

    labels = [r['label'] for r in results[1:]]
    # Δρ = (k_pert − k_nom) / (k_pert × k_nom); σ via error propagation on 1/k
    drho_pcm = [
        (r['k'] - k_nom) / (r['k'] * k_nom) * 1e5
        for r in results[1:]
    ]
    drho_err = [
        np.sqrt((dk_nom / k_nom ** 2) ** 2 + (r['dk'] / r['k'] ** 2) ** 2) * 1e5
        for r in results[1:]
    ]
    colors = ['#d62728' if d < 0 else '#1f77b4' for d in drho_pcm]

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(labels))
    bars = ax.bar(x, drho_pcm, color=colors, alpha=0.85, zorder=3)
    ax.errorbar(x, drho_pcm, yerr=drho_err, fmt='none', ecolor='black', capsize=5, zorder=4)
    ax.axhline(0.0, color='black', linewidth=0.9, zorder=2)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Δρ  (pcm)')
    ax.set_title(
        'Manufacturing Tolerance Sensitivity — Δρ vs. Nominal\n'
        f'Baseline: PF = {_PF_NOM:.0%},  kernel diameter = {_R_KERNEL_NOM * 2e4:.0f} µm'
    )
    ax.grid(True, axis='y', alpha=0.3, zorder=0)

    for bar, val, err in zip(bars, drho_pcm, drho_err):
        offset = max(err, 50) * (1 if val >= 0 else -1) + (20 if val >= 0 else -30)
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val + offset,
            f'{val:+.0f} pcm',
            ha='center', va='center', fontsize=8,
        )

    fig.tight_layout()
    out = output_dir / 'manufacturing_tolerance.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Plot saved → {out}')


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def run_sensitivity(output_dir: Path) -> int:
    """Run all cases, print the sensitivity table, save plot; return 0 on success."""
    output_dir.mkdir(parents=True, exist_ok=True)
    work_root = output_dir / 'mfg_tolerance'

    results: list[dict] = []
    print(f'Processing {len(_CASES)} cases ({_BATCHES} batches / {_PARTICLES:,} particles each) ...')
    print()

    for case in _CASES:
        diam_um = case['r_kernel'] * 2e4
        work_dir = work_root / case['short']
        sp_file = work_dir / f'statepoint.{_BATCHES}.h5'
        print(
            f"  [{case['label']}]  "
            f"PF = {case['pf']:.0%}  "
            f"kernel ⌀ = {diam_um:.0f} µm  →  {work_dir}"
        )
        if sp_file.exists():
            print(f'    [cached] {sp_file}')
            with openmc.StatePoint(sp_file) as sp:
                k, dk = float(sp.keff.nominal_value), float(sp.keff.std_dev)
        else:
            mats = build_materials()
            geom = _build_geometry(mats, pf=case['pf'], r_kernel=case['r_kernel'])
            k, dk = _run_keff(mats, geom, work_dir)
        results.append({'label': case['label'], 'k': k, 'dk': dk})
        print(f'    k-inf = {k:.5f} ± {dk:.5f}')

    # --- Checkpoint report ---
    k_nom  = results[0]['k']
    dk_nom = results[0]['dk']

    print()
    print('=' * 72)
    print('MANUFACTURING TOLERANCE SENSITIVITY — CHECKPOINT')
    print(f'Baseline: PF = {_PF_NOM:.0%},  kernel ⌀ = {_R_KERNEL_NOM*2e4:.0f} µm')
    print(f'Kernel tolerance: ±{_R_KERNEL_DELTA*1e4:.0f} µm radius '
          f'(±{_R_KERNEL_DELTA*2e4:.0f} µm diameter)  — INL/EXT-10-19476 Table 3')
    print('=' * 72)
    hdr = f"{'Case':<30}  {'k-eff':>9}  {'σ_k ×10⁵':>9}  {'Δρ (pcm)':>10}  {'|Δρ|/σ_Δρ':>10}"
    print(hdr)
    print('-' * len(hdr))

    row = f"{'Nominal':30}  {k_nom:.5f}  {dk_nom*1e5:>9.1f}  {'—':>10}  {'—':>10}"
    print(row)

    max_abs_drho: dict[str, float] = {'pf': 0.0, 'kernel': 0.0}

    for r in results[1:]:
        drho_pcm = (r['k'] - k_nom) / (r['k'] * k_nom) * 1e5
        sigma_drho = np.sqrt((dk_nom / k_nom ** 2) ** 2 + (r['dk'] / r['k'] ** 2) ** 2) * 1e5
        sig_ratio = abs(drho_pcm) / sigma_drho if sigma_drho > 0 else float('inf')
        print(
            f"  {r['label']:<28}  {r['k']:.5f}  {r['dk']*1e5:>9.1f}"
            f"  {drho_pcm:>+10.0f}  {sig_ratio:>9.1f}σ"
        )
        if 'PF' in r['label']:
            max_abs_drho['pf'] = max(max_abs_drho['pf'], abs(drho_pcm))
        else:
            max_abs_drho['kernel'] = max(max_abs_drho['kernel'], abs(drho_pcm))

    print()
    print(f"Max |Δρ|  PF perturbation (±2 pp):               {max_abs_drho['pf']:>6.0f} pcm")
    print(f"Max |Δρ|  kernel diameter perturbation (±10 µm): {max_abs_drho['kernel']:>6.0f} pcm")

    if max_abs_drho['pf'] >= max_abs_drho['kernel']:
        more_sensitive = 'packing fraction (±2 percentage points)'
        less_sensitive = 'kernel diameter (±10 µm)'
    else:
        more_sensitive = 'kernel diameter (±10 µm)'
        less_sensitive = 'packing fraction (±2 percentage points)'

    ratio = (
        max_abs_drho['pf'] / max_abs_drho['kernel']
        if max_abs_drho['kernel'] > 0
        else float('inf')
    )
    print()
    print(f'Reactivity is more sensitive to {more_sensitive}')
    print(f'  (sensitivity ratio: {ratio:.1f}× vs. {less_sensitive})')

    _plot_sensitivity(results, output_dir)
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('output'),
        help='Directory for plots and eigenvalue run subdirs (default: output/)',
    )
    args = parser.parse_args()
    sys.exit(run_sensitivity(args.output_dir))


if __name__ == '__main__':
    main()
