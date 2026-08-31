---
name: openmc-skill
description: Use when writing, debugging, or explaining OpenMC Python code — including geometry/material/tally setup, running simulations, fuel depletion (burnup), multi-group cross section generation, source distributions, and reading nuclear data (ACE/ENDF). Also invoke for TRISO fuel modeling, PWR assembly/core models, or any openmc.* API question.
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
  - WebFetch
---

# OpenMC Python API Skill

## What is OpenMC?

OpenMC is an open-source Monte Carlo particle transport code for nuclear reactor simulations. Its Python API (`openmc`) provides a complete interface for building geometry, materials, tallies, and running simulations — all from Python.

## When to use this skill

Use this skill when:
- Building geometry, materials, or tallies for a neutronics simulation
- Running criticality (k-eff) or fixed-source transport calculations
- Performing fuel depletion/burnup analysis
- Generating multi-group cross sections for diffusion or deterministic codes
- Reading/processing nuclear data (ACE, ENDF, NJOY output)
- Defining source distributions or variance reduction parameters

## Key modules

| Module | Purpose |
|--------|---------|
| `openmc` | Core: geometry, materials, tallies, settings, running |
| `openmc.model` | Convenience classes/functions for model assembly |
| `openmc.examples` | Prebuilt reference models (PWR pin, assembly, core) |
| `openmc.deplete` | Fuel depletion (Bateman equation solvers, operators) |
| `openmc.mgxs` | Multi-group cross section generation from tally data |
| `openmc.stats` | Probability distributions for source definitions |
| `openmc.data` | Low-level nuclear data interface (ACE, ENDF, NJOY) |

## Reference files

- [reference/basic-functionality.md](reference/basic-functionality.md) — Materials, geometry, surfaces, tallies, meshes, settings, running, post-processing
- [reference/model-building.md](reference/model-building.md) — `openmc.model`: convenience functions, composite surfaces, TRISO fuel
- [reference/example-models.md](reference/example-models.md) — `openmc.examples`: prebuilt PWR and slab models
- [reference/depletion.md](reference/depletion.md) — `openmc.deplete`: integrators, operators, chain data, solvers
- [reference/mgxs.md](reference/mgxs.md) — `openmc.mgxs`: energy groups, cross section types, libraries
- [reference/statistics.md](reference/statistics.md) — `openmc.stats`: univariate, angular, and spatial distributions
- [reference/nuclear-data.md](reference/nuclear-data.md) — `openmc.data`: neutron/photon data, resonances, ACE/ENDF/NJOY interfaces

## Minimal working example

```python
import openmc

# Material
fuel = openmc.Material()
fuel.add_nuclide('U235', 0.04)
fuel.add_nuclide('U238', 0.96)
fuel.add_element('O', 2.0)
fuel.set_density('g/cm3', 10.5)

water = openmc.Material()
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.set_density('g/cm3', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([fuel, water])

# Geometry
fuel_cyl = openmc.ZCylinder(r=0.4)
outer_cyl = openmc.ZCylinder(r=0.6, boundary_type='reflective')
top = openmc.ZPlane(z0=10.0, boundary_type='reflective')
bot = openmc.ZPlane(z0=-10.0, boundary_type='reflective')

fuel_cell = openmc.Cell(region=-fuel_cyl & +bot & -top, fill=fuel)
clad_cell = openmc.Cell(region=+fuel_cyl & -outer_cyl & +bot & -top, fill=water)

universe = openmc.Universe(cells=[fuel_cell, clad_cell])
geometry = openmc.Geometry(universe)

# Settings
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 20
settings.particles = 1000
settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((0, 0, 0))
)

# Run
openmc.run()
```

## Docs
https://docs.openmc.org/en/stable/pythonapi/index.html
