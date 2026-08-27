"""TRISO fuel compact geometry.

Particle dimensions from AGR-1 fuel specifications:
  INL/EXT-10-19476 — AGR-1 Fuel Particle and Compact Characterization Data
  Demkowicz et al., Nucl. Eng. Des. 329 (2018) 102–111
"""

import math

import openmc
import openmc.model

from .materials import build_materials

# AGR-1 TRISO particle layer outer radii [cm] — INL/EXT-10-19476 Table 3.
_R_KERNEL = 0.0175  # kernel only
_R_BUFFER = 0.0275  # + buffer (100 µm)
_R_IPYC   = 0.0315  # + IPyC  ( 40 µm)
_R_SIC    = 0.0350  # + SiC   ( 35 µm)
_R_OPYC   = 0.0390  # + OPyC  ( 40 µm)

# AGR-1 compact cylinder dimensions [cm] — INL/EXT-10-19476.
_COMPACT_R = 0.62
_COMPACT_H = 2.50

_PACKING_FRACTION = 0.30          # design target per step-2 specification
_LATTICE_PITCH    = 2.0 * _R_OPYC # 2× outer radius; avoids cell-overlap artefacts


def _triso_universe(mats: dict) -> openmc.Universe:
    s1 = openmc.Sphere(r=_R_KERNEL)
    s2 = openmc.Sphere(r=_R_BUFFER)
    s3 = openmc.Sphere(r=_R_IPYC)
    s4 = openmc.Sphere(r=_R_SIC)
    s5 = openmc.Sphere(r=_R_OPYC)
    cells = [
        openmc.Cell(fill=mats['kernel'], region=-s1),
        openmc.Cell(fill=mats['buffer'], region=+s1 & -s2),
        openmc.Cell(fill=mats['ipyc'],   region=+s2 & -s3),
        openmc.Cell(fill=mats['sic'],    region=+s3 & -s4),
        openmc.Cell(fill=mats['opyc'],   region=+s4 & -s5),
    ]
    return openmc.Universe(cells=cells)


def build_geometry(mats: dict | None = None) -> openmc.Geometry:
    """Build a single TRISO compact with reflective boundary conditions.

    Parameters
    ----------
    mats:
        Material dict from build_materials(). Created internally if omitted.

    Returns
    -------
    openmc.Geometry
        Root geometry containing the TRISO lattice embedded in a graphite
        compact, bounded by reflective surfaces on all faces.
    """
    if mats is None:
        mats = build_materials()

    triso_univ = _triso_universe(mats)

    # All-reflective boundaries model an infinite lattice of identical compacts
    # (standard lattice-physics starting point for k-inf).
    # TODO: switch axial planes to vacuum to approximate a finite fuel column.
    cyl = openmc.ZCylinder(r=_COMPACT_R, boundary_type='reflective')
    top = openmc.ZPlane(z0=+_COMPACT_H / 2, boundary_type='reflective')
    bot = openmc.ZPlane(z0=-_COMPACT_H / 2, boundary_type='reflective')
    compact_region = -cyl & -top & +bot

    centers = openmc.model.pack_spheres(
        radius=_R_OPYC,
        region=compact_region,
        pf=_PACKING_FRACTION,
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

    # TODO: add helium coolant channel surrounding the compact (annular region
    # between compact cylinder and a larger outer cylinder). Omitted in Stage 0
    # because it requires a second bounding cylinder and an outer graphite sleeve,
    # which adds non-trivial geometry complexity without changing the lattice physics.

    compact_cell = openmc.Cell(fill=triso_lat, region=compact_region)
    return openmc.Geometry(openmc.Universe(cells=[compact_cell]))
