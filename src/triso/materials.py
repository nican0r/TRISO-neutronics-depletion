"""TRISO particle material definitions.

All compositions and densities from AGR-1 fuel specifications:
  Demkowicz et al., Nucl. Eng. Des. 329 (2018) 102–111
  INL/EXT-10-19476 — AGR-1 Fuel Particle and Compact Characterization Data
"""

import openmc

# 293.6 K matches the single temperature point in the NNDC HDF5 library
# downloaded by download_data.sh. Using any other value requires a
# multi-temperature library. Revisit when upgrading the cross-section data.
_T_K = 293.6


def build_materials() -> dict[str, openmc.Material]:
    """Return TRISO layer materials keyed by layer name."""

    # UCO kernel: UC0.5O0.4, 19.75 wt% U-235 (HALEU).
    # Atom ratios C/U = 0.5, O/U = 0.4 from AGR-1 characterization.
    # Density 10.5 g/cm³ (AGR-1 target 10.4–10.6 g/cm³, INL/EXT-10-19476).
    # Enrichment is weight percent U-235 relative to total uranium.
    kernel = openmc.Material(name='UCO kernel')
    kernel.add_element('U', 1.0, enrichment=19.75)
    kernel.add_element('C', 0.5)
    kernel.add_element('O', 0.4)
    kernel.set_density('g/cm3', 10.5)
    kernel.temperature = _T_K

    # Porous carbon buffer: ~50% dense vs. solid PyC.
    # Density 1.0 g/cm³ (AGR-1 target 1.0 g/cm³, INL/EXT-10-19476 Table 3).
    buffer = openmc.Material(name='porous carbon buffer')
    buffer.add_element('C', 1.0)
    buffer.set_density('g/cm3', 1.0)
    buffer.temperature = _T_K

    # Inner pyrolytic carbon (IPyC).
    # Density 1.87 g/cm³ (AGR-1 target 1.85–1.90 g/cm³, INL/EXT-10-19476).
    # TODO: PyC is turbostratic carbon; c_Graphite S(α,β) is an approximation.
    # Stage 0 omits thermal scattering for PyC; revisit in Stage 1.
    ipyc = openmc.Material(name='inner PyC')
    ipyc.add_element('C', 1.0)
    ipyc.set_density('g/cm3', 1.87)
    ipyc.temperature = _T_K

    # Silicon carbide (β-SiC pressure vessel layer).
    # Density 3.20 g/cm³ (AGR-1 target 3.19 g/cm³, INL/EXT-10-19476 Table 3).
    sic = openmc.Material(name='SiC')
    sic.add_element('Si', 1.0)
    sic.add_element('C', 1.0)
    sic.set_density('g/cm3', 3.20)
    sic.temperature = _T_K

    # Outer pyrolytic carbon (OPyC).
    # Same density target as IPyC per AGR-1 spec.
    opyc = openmc.Material(name='outer PyC')
    opyc.add_element('C', 1.0)
    opyc.set_density('g/cm3', 1.87)
    opyc.temperature = _T_K

    # Graphite matrix (compact matrix material — distinct object from PyC layers).
    # Density 1.75 g/cm³ (AGR-1 compact matrix target, INL/EXT-10-19476 Table 5).
    # c_Graphite S(α,β) applied: matrix is nuclear-grade graphite, not PyC.
    matrix = openmc.Material(name='graphite matrix')
    matrix.add_element('C', 1.0)
    matrix.add_s_alpha_beta('c_Graphite')
    matrix.set_density('g/cm3', 1.75)
    matrix.temperature = _T_K

    return {
        'kernel': kernel,
        'buffer': buffer,
        'ipyc': ipyc,
        'sic': sic,
        'opyc': opyc,
        'matrix': matrix,
    }
