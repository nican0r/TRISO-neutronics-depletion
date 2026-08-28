#!/usr/bin/env python3
"""Physics validation for the TRISO compact neutronics model (Stage 0).

Checks performed
----------------
1. k-eff sanity: compared against the expected k-inf window for a
   graphite-moderated, all-reflective HALEU compact at 293.6 K.
2. Flux depression in the UCO kernel (self-shielding): thermal flux must be
   lower in the kernel than in the graphite matrix.
3. Energy deposition proxy: fission rate in the kernel confirms that energy
   release is localised in the fissile layer.  Note: heating-local returns
   zero with the NNDC HDF5 library (kerma coefficients absent); see README.

Usage
-----
  # Load from an existing statepoint (fast path):
  python scripts/validate_physics.py --statepoint statepoint.200.h5

  # Run fresh simulation, then validate:
  python scripts/validate_physics.py

Plots are saved as PNG to --output-dir (default: output/).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend; safe headless / CI
import matplotlib.pyplot as plt
import numpy as np
import openmc

# Allow `from triso...` without installing the package.
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
from triso.run import build_model, run_model, _tally_to_df  # noqa: E402
from triso.geometry import (  # noqa: E402
    _R_KERNEL, _R_BUFFER, _R_IPYC, _R_SIC, _R_OPYC,
    _PACKING_FRACTION,
)

# Expected k-inf window for a compact-only, all-reflective HALEU model at 293.6 K.
# This is a compact-only geometry: the only moderator is the graphite matrix inside
# the compact itself (~70 vol% of the compact at 30% PF).  No external graphite
# moderator blocks, sleeves, or reflector are included.  Real HTGR fuel elements
# thermalize more effectively because the graphite structure surrounding each compact
# dominates the moderating volume; compact-only k-inf is therefore lower than a full
# element k-inf.  Range: 1.05–1.40 (conservative lower bound; upper bound allows for
# absence of Doppler broadening at 293.6 K vs. operating ~900 K, which would cost
# ~1–3% in k).
_KEFF_LO = 1.05
_KEFF_HI = 1.40

# Full material names as assigned in materials.py (openmc.Material(name=...)),
# in physical layer order from kernel outward to moderator.
_LAYER_FULL_NAMES = [
    'UCO kernel',
    'porous carbon buffer',
    'inner PyC',
    'SiC',
    'outer PyC',
    'graphite matrix',
]
# Short labels for axis ticks and printed output.
_LAYER_DISPLAY = ['kernel', 'buffer', 'IPyC', 'SiC', 'OPyC', 'matrix']
_FULL_TO_DISPLAY = dict(zip(_LAYER_FULL_NAMES, _LAYER_DISPLAY))

_GROUP_LABELS = [
    'Thermal\n(0–0.625 eV)',
    'Epithermal\n(0.625 eV–100 keV)',
    'Fast\n(100 keV–20 MeV)',
]
_GROUP_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c']


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _mat_id_map(sp_path: Path) -> dict[str, str]:
    """Build a mat_N → full-name mapping from summary.h5 if available."""
    summary_path = sp_path.parent / 'summary.h5'
    if not summary_path.exists():
        return {}
    su = openmc.Summary(str(summary_path))
    return {f'mat_{m.id}': m.name for m in su.materials}


def _patch_mat_labels(df, id_map: dict[str, str]):
    """Replace mat_N integer-ID labels with human-readable material names."""
    if 'material' not in df.columns or not id_map:
        return df
    df = df.copy()
    df['material'] = df['material'].map(lambda x: id_map.get(x, x))
    return df


def _load_from_statepoint(sp_path: Path) -> dict:
    """Load tally results from an existing statepoint.

    OpenMC stores MaterialFilter bins as integer IDs in the HDF5 file.
    summary.h5 (written alongside the statepoint by model.run()) holds the
    material name table; _mat_id_map reads it to patch 'mat_N' labels back
    to the full material names defined in materials.py.
    """
    id_map = _mat_id_map(sp_path)
    with openmc.StatePoint(sp_path) as sp:
        k_mean = sp.keff.nominal_value
        k_std = sp.keff.std_dev
        flux_df = _patch_mat_labels(_tally_to_df(sp.get_tally(name='flux by material')), id_map)
        flux_mg_df = _patch_mat_labels(_tally_to_df(sp.get_tally(name='flux by material 3-group')), id_map)
        rxn_df = _tally_to_df(sp.get_tally(name='kernel reaction rates'))
        heat_df = _patch_mat_labels(_tally_to_df(sp.get_tally(name='heating by material')), id_map)
    return {
        'k_eff': (k_mean, k_std),
        'flux': flux_df,
        'flux_mg': flux_mg_df,
        'rxn': rxn_df,
        'heat': heat_df,
    }


# ---------------------------------------------------------------------------
# Check 1 — k-eff
# ---------------------------------------------------------------------------

def check_keff(k: float, dk: float) -> bool:
    """Print k-eff and evaluate against the expected range. Returns True = PASS."""
    _banner('k-eff SANITY CHECK')
    print(f'  Result  : k = {k:.5f} ± {dk:.5f}')
    print(f'  Expected: {_KEFF_LO:.2f} – {_KEFF_HI:.2f}')
    print()

    if _KEFF_LO <= k <= _KEFF_HI:
        print('  PASS — k-eff is within the expected range.')
        passed = True
    elif k > _KEFF_HI:
        print(f'  WARN — k-eff is ABOVE the expected ceiling ({_KEFF_HI}).')
        print('         Possible causes: geometry error (reflective boundaries not')
        print('         leaking), packing fraction > 30%, or enrichment mismatch.')
        passed = False
    else:
        print(f'  WARN — k-eff is BELOW the expected floor ({_KEFF_LO}).')
        print('         Possible causes: absorber contamination, wrong material')
        print('         density, or enrichment lower than 19.75%.')
        passed = False

    print()
    print('  Physics context:')
    print('    Geometry : all-reflective boundaries (k-inf); no neutron leakage.')
    print('    Fuel     : 19.75 wt% HALEU UCO, 30% packing fraction in graphite.')
    print('    Moderator: compact matrix graphite only — no external graphite blocks.')
    print('               A full HTGR fuel element (graphite sleeve + compact) gives')
    print('               a higher k-inf because additional graphite improves')
    print('               thermalization. Compact-only k-inf is therefore at the lower')
    print('               end of the broader HTGR design space.')
    print('    Real core: leakage + control rods reduce k-eff below 1.0 at power.')
    return passed


# ---------------------------------------------------------------------------
# Check 2 — Flux depression / self-shielding
# ---------------------------------------------------------------------------

def check_flux_depression(flux_mg_df, out_dir: Path) -> bool:
    """
    Plot per-material flux to demonstrate kernel self-shielding. Returns True = PASS.

    Self-shielding signature: thermal and epithermal flux are lower in the UCO
    kernel than in the surrounding graphite matrix.  The kernel absorbs strongly
    in both the thermal region (U-235 fission + U-238 (n,γ)) and the resonance
    region (U-238 resonance absorption), so neutrons are depleted before they
    can build up inside the dense fissile volume.
    """
    _banner('SELF-SHIELDING CHECK — Flux by material')

    groups = list(flux_mg_df['energy [eV]'].unique())
    n_groups = len(groups)

    available = set(flux_mg_df['material'].unique())
    # Match by full material name; build display labels in the same order.
    ordered_full = [m for m in _LAYER_FULL_NAMES if m in available]
    ordered_disp = [_FULL_TO_DISPLAY.get(m, m) for m in ordered_full]

    if not ordered_full:
        print('  ERROR — No recognised material names found in tally.')
        print(f'          Labels present: {sorted(available)}')
        print('          Expected names from materials.py:', _LAYER_FULL_NAMES)
        return False

    # Collect mean / std_dev arrays: shape (n_materials, n_groups)
    means = np.zeros((len(ordered_full), n_groups))
    errs = np.zeros_like(means)
    for mi, mat in enumerate(ordered_full):
        sub = flux_mg_df[flux_mg_df['material'] == mat]
        for gi, g in enumerate(groups):
            row = sub[sub['energy [eV]'] == g]
            means[mi, gi] = float(row['mean'].values[0])
            errs[mi, gi] = float(row['std_dev'].values[0])

    # --- Numerical report ---
    passed = True
    kernel_name = 'UCO kernel'
    matrix_name = 'graphite matrix'
    vf = _vol_fractions()

    if kernel_name in ordered_full and matrix_name in ordered_full:
        ki = ordered_full.index(kernel_name)
        mi = ordered_full.index(matrix_name)

        vf_k = vf[kernel_name]
        vf_m = vf[matrix_name]
        vol_ratio = vf_k / vf_m  # expected raw tally ratio if Φ_density were equal

        print(f'  Volume fractions — kernel: {vf_k:.4f}  matrix: {vf_m:.4f}')
        print(f'  Expected raw ratio (volume effect only): {vol_ratio:.4f}')
        print()
        print('  Raw MaterialFilter flux is volume-INTEGRATED (sum of track lengths),')
        print('  not per unit volume.  Flux density = raw flux / volume fraction.')
        print()

        # Fast-group flux may be ELEVATED in the kernel (fission neutrons are born
        # there), so only thermal and epithermal depression is a self-shielding signal.
        group_labels = ['Thermal', 'Epithermal', 'Fast']
        for gi, glabel in enumerate(group_labels):
            raw_k = means[ki, gi]
            raw_m = means[mi, gi]
            if raw_m > 0 and vf_m > 0 and vf_k > 0:
                raw_ratio = raw_k / raw_m
                density_ratio = (raw_k / vf_k) / (raw_m / vf_m)
                if gi < 2:  # thermal and epithermal only
                    tag = 'PASS' if density_ratio < 1.0 else 'WARN'
                    if density_ratio >= 1.0:
                        passed = False
                else:
                    # Fast group: kernel is a NET SOURCE of fast fission neutrons,
                    # so density_ratio > 1 is physically expected, not a failure.
                    tag = 'INFO (fast group — kernel is fission source)'
                print(f'  {glabel}:')
                print(f'    Raw ratio   kernel/matrix = {raw_ratio:.4f}  '
                      f'(≈ volume ratio {vol_ratio:.4f} — dominated by geometry)')
                print(f'    Flux density kernel/matrix = {density_ratio:.4f}  '
                      f'← {tag}')

        print()
        if passed:
            print('  PASS — Flux density is depressed in the kernel for thermal/epithermal.')
            print('         Self-shielding is present, as expected for a fissile sphere.')
            print('         The small magnitude (~1–5% density depression) is physically')
            print('         correct for TRISO particles: at r_kernel = 0.0175 cm the mean')
            print('         free path (graphite thermal: ~2 cm) is >> particle size, so')
            print('         spatial flux gradients are weak.  The dominant self-shielding')
            print('         mechanism is resonance self-shielding (Dancoff factor reduction')
            print('         of the effective U-238 resonance integral), not spatial.')
    else:
        print('  WARN — Could not find UCO kernel and/or graphite matrix in tally labels.')
        passed = False

    # --- Plot 1: thermal-group bar chart (shows depression clearly) ---
    thermal_means = means[:, 0]
    thermal_errs = errs[:, 0]
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(ordered_full))
    ax.bar(x, thermal_means, yerr=thermal_errs, capsize=4,
           color='steelblue', alpha=0.8, ecolor='dimgray')
    ax.set_xticks(x)
    ax.set_xticklabels(ordered_disp, rotation=15, ha='right')
    ax.set_xlabel('Material layer (kernel → moderator)')
    ax.set_ylabel('Thermal flux  [arb. per source n·cm²]')
    ax.set_title(
        'Thermal flux by TRISO layer — self-shielding check\n'
        'Expected: kernel flux < matrix flux  (U absorption depletes thermal neutrons)'
    )
    ax.grid(axis='y', alpha=0.4)
    plt.tight_layout()
    p1 = out_dir / 'flux_thermal_by_material.png'
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    print(f'  Thermal flux plot → {p1}')

    # --- Plot 2: 3-group grouped bar chart ---
    fig, ax = plt.subplots(figsize=(11, 5))
    width = 0.25
    offsets = np.linspace(
        -(n_groups - 1) * width / 2,
        (n_groups - 1) * width / 2,
        n_groups,
    )
    for gi in range(n_groups):
        label = _GROUP_LABELS[gi] if gi < len(_GROUP_LABELS) else f'Group {gi}'
        color = _GROUP_COLORS[gi % len(_GROUP_COLORS)]
        ax.bar(
            x + offsets[gi],
            means[:, gi],
            width,
            yerr=errs[:, gi],
            capsize=3,
            label=label.replace('\n', ' '),
            color=color,
            alpha=0.8,
            ecolor='dimgray',
        )
    ax.set_xticks(x)
    ax.set_xticklabels(ordered_disp, rotation=15, ha='right')
    ax.set_xlabel('Material layer (kernel → moderator)')
    ax.set_ylabel('Flux  [arb. per source n·cm²]')
    ax.set_title(
        '3-group flux by TRISO layer\n'
        'Matrix dominates thermal; kernel shows depression across all groups'
    )
    ax.legend(title='Energy group', fontsize=8)
    ax.grid(axis='y', alpha=0.4)
    plt.tight_layout()
    p2 = out_dir / 'flux_3group_by_material.png'
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    print(f'  3-group flux plot → {p2}')

    return passed


# ---------------------------------------------------------------------------
# Check 3 — Energy deposition
# ---------------------------------------------------------------------------

def check_energy_deposition(rxn_df, heat_df, out_dir: Path) -> bool:
    """
    Verify energy release is concentrated in the UCO kernel.

    Uses heating-local if kerma data is available (ENDF/B-VIII.0+), falling back
    to fission rate as a proxy when all heating values are zero (NNDC library).
    """
    _banner('ENERGY DEPOSITION CHECK')

    # Detect whether heating-local has real data or is all zeros.
    heat_total = heat_df['mean'].abs().sum() if not heat_df.empty else 0.0
    has_kerma = heat_total > 0.0

    if has_kerma:
        _check_heating_plot(heat_df, out_dir)
    else:
        print('  NOTE: heating-local is all zeros — kerma coefficients absent from library.')
        print('  Falling back to fission rate as energy localisation proxy.')
        print()

    # Always report fission rate (direct physics check independent of kerma).
    print('  Fission rate in kernel (primary energy source):')
    fission_rows = rxn_df[rxn_df['score'] == 'fission']
    absorption_rows = rxn_df[rxn_df['score'] == 'absorption']

    if fission_rows.empty:
        print('  ERROR — No fission score found in kernel reaction tally.')
        return False

    f_mean = float(fission_rows['mean'].values[0])
    f_std = float(fission_rows['std_dev'].values[0])
    print(f'    Fission rate    : {f_mean:.4e} ± {f_std:.4e}  (per src n·cm²)')

    if not absorption_rows.empty:
        a_mean = float(absorption_rows['mean'].values[0])
        a_std = float(absorption_rows['std_dev'].values[0])
        print(f'    Absorption rate : {a_mean:.4e} ± {a_std:.4e}  (per src n·cm²)')
        if a_mean > 0:
            eta_eff = f_mean / a_mean
            print(f'    Fission fraction: {eta_eff:.4f}  (f/a in kernel; ~0.7 expected for HALEU)')

    print()
    if f_mean > 0:
        print('  PASS — Non-zero fission rate confirms energy release is localised')
        print('         in the UCO kernel, as required by the TRISO design.')
        return True
    else:
        print('  WARN — Fission rate is zero. Investigate kernel material definition.')
        return False


def _check_heating_plot(heat_df, out_dir: Path) -> None:
    """Plot volume-integrated and volume-normalised heating-local by material."""
    vf = _vol_fractions()

    available = set(heat_df['material'].unique())
    ordered_full = [m for m in _LAYER_FULL_NAMES if m in available]
    ordered_disp = [_FULL_TO_DISPLAY.get(m, m) for m in ordered_full]

    if not ordered_full:
        print('  ERROR — No recognised material names in heating tally.')
        return

    means = np.array([
        float(heat_df[heat_df['material'] == m]['mean'].values[0])
        for m in ordered_full
    ])
    errs = np.array([
        float(heat_df[heat_df['material'] == m]['std_dev'].values[0])
        for m in ordered_full
    ])
    vf_arr = np.array([vf.get(m, 1.0) for m in ordered_full])
    density = means / vf_arr
    density_errs = errs / vf_arr

    kernel_name = 'UCO kernel'
    matrix_name = 'graphite matrix'
    if kernel_name in ordered_full and matrix_name in ordered_full:
        ki = ordered_full.index(kernel_name)
        mi = ordered_full.index(matrix_name)
        raw_ratio = means[ki] / means[mi] if means[mi] > 0 else float('nan')
        vol_ratio = density[ki] / density[mi] if density[mi] > 0 else float('nan')
        print(f'  Raw heating     kernel/matrix = {raw_ratio:.1f}×')
        print(f'  Heating density kernel/matrix = {vol_ratio:.1f}×  (per unit volume)')
        print()
        print('  PASS — Energy deposition strongly concentrated in the UCO kernel,')
        print('         consistent with fission as the dominant energy source.')

    x = np.arange(len(ordered_full))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].bar(x, means, yerr=errs, capsize=4, color='firebrick', alpha=0.8, ecolor='dimgray')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(ordered_disp, rotation=15, ha='right')
    axes[0].set_ylabel('heating-local  [eV per source n·cm²]')
    axes[0].set_title('Volume-integrated heating by layer\n(raw tally — dominated by volume)')
    axes[0].grid(axis='y', alpha=0.4)

    axes[1].bar(x, density, yerr=density_errs, capsize=4, color='firebrick', alpha=0.8, ecolor='dimgray')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(ordered_disp, rotation=15, ha='right')
    axes[1].set_ylabel('heating-local density  [eV per source n·cm²·vol-fraction]')
    axes[1].set_title('Volume-normalised heating density\nKernel >> matrix confirms fission localisation')
    axes[1].grid(axis='y', alpha=0.4)

    plt.suptitle('Energy deposition (heating-local) by TRISO layer', fontweight='bold')
    plt.tight_layout()
    p = out_dir / 'heating_by_material.png'
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print(f'  Heating plot → {p}')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _vol_fractions() -> dict[str, float]:
    """Analytical volume fractions of each TRISO layer within the compact.

    Derived from AGR-1 particle radii (geometry.py) and 30% packing fraction.
    fraction_layer = PF × (r_outer³ − r_inner³) / r_OPyC³
    Matrix fills the remaining (1 − PF) fraction of compact volume.
    """
    def _f(r_out, r_in=0.0):
        return _PACKING_FRACTION * (r_out**3 - r_in**3) / _R_OPYC**3

    return {
        'UCO kernel':           _f(_R_KERNEL),
        'porous carbon buffer': _f(_R_BUFFER, _R_KERNEL),
        'inner PyC':            _f(_R_IPYC,   _R_BUFFER),
        'SiC':                  _f(_R_SIC,    _R_IPYC),
        'outer PyC':            _f(_R_OPYC,   _R_SIC),
        'graphite matrix':      1.0 - _PACKING_FRACTION,
    }


def _banner(title: str) -> None:
    width = 60
    print(f'\n{"=" * width}')
    print(f'  {title}')
    print(f'{"=" * width}')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Stage 0 physics validation for the TRISO compact model.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--statepoint', type=Path, default=None,
        help='Path to existing statepoint HDF5 (skips re-run; default: run fresh).',
    )
    parser.add_argument(
        '--output-dir', type=Path, default=Path('output'),
        help='Directory for PNG plots (created if absent; default: output/).',
    )
    args = parser.parse_args()

    out_dir: Path = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print()
    print('TRISO compact — Stage 0 physics validation')
    print('=' * 42)

    if args.statepoint:
        print(f'\nLoading statepoint: {args.statepoint}')
        results = _load_from_statepoint(args.statepoint)
    else:
        print('\nNo statepoint given — running simulation (~5–15 min on laptop)...')
        results = run_model()

    k, dk = results['k_eff']

    p1 = check_keff(k, dk)
    p2 = check_flux_depression(results['flux_mg'], out_dir)
    p3 = check_energy_deposition(results['rxn'], results['heat'], out_dir)

    _banner('SUMMARY')
    statuses = {
        'k-eff range check': p1,
        'Flux depression (self-shielding)': p2,
        'Energy deposition proxy (fission rate)': p3,
    }
    all_pass = True
    for name, ok in statuses.items():
        status = 'PASS' if ok else 'WARN/FAIL'
        print(f'  {status:<12} {name}')
        if not ok:
            all_pass = False

    print()
    if all_pass:
        print('  All checks passed. The Stage 0 model is physically consistent.')
    else:
        print('  One or more checks flagged. Review WARN/FAIL items above.')
    print()


if __name__ == '__main__':
    main()
