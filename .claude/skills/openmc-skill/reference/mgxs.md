# OpenMC Multi-Group Cross Section Generation (`openmc.mgxs`)

API docs: https://docs.openmc.org/en/stable/pythonapi/mgxs.html

The `openmc.mgxs` module generates multi-group cross sections (MGXS) from OpenMC continuous-energy Monte Carlo tallies. The output can feed deterministic codes (e.g., OpenMOC, MPACT) or OpenMC's own multi-group mode.

---

## Workflow overview

```python
import openmc
import openmc.mgxs

# 1. Define energy group structure
groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 0.625e-6, 20.0e6])  # 2 groups (eV)

# 2. Create MGXS objects and add tallies to model
total_xs = openmc.mgxs.TotalXS(domain=fuel_cell, groups=groups)
fission_xs = openmc.mgxs.FissionXS(domain=fuel_cell, groups=groups)
scatter_matrix = openmc.mgxs.ScatterMatrixXS(domain=fuel_cell, groups=groups)

# 3. Add required tallies to the model
tallies = openmc.Tallies()
tallies += total_xs.tallies.values()
tallies += fission_xs.tallies.values()
tallies += scatter_matrix.tallies.values()

# 4. Run simulation, then load results
with openmc.StatePoint('statepoint.100.h5') as sp:
    total_xs.load_from_statepoint(sp)
    fission_xs.load_from_statepoint(sp)
    scatter_matrix.load_from_statepoint(sp)

# 5. Export / inspect
print(total_xs.get_pandas_dataframe())
total_xs.export_xs_data(filename='mgxs.h5', xs_type='macro')
```

---

## Energy Groups

### `openmc.mgxs.EnergyGroups(group_edges)`

| Parameter | Description |
|---|---|
| `group_edges` | Array of energy bin edges in eV (length = num_groups + 1) |

Key attributes:
- `.num_groups` — number of energy groups
- `.group_edges` — energy bin boundaries

### Built-in group structures

Access via `openmc.mgxs.GROUP_STRUCTURES` (dict):

| Key | Groups | Description |
|---|---|---|
| `'CASMO-2'` | 2 | CASMO 2-group |
| `'CASMO-8'` | 8 | CASMO 8-group |
| `'CASMO-16'` | 16 | CASMO 16-group |
| `'CASMO-70'` | 70 | CASMO 70-group |
| `'XMAS-172'` | 172 | XMAS 172-group |
| `'SHEM-361'` | 361 | SHEM 361-group |
| `'SCALE-252'` | 252 | SCALE 252-group |
| `'MPACT-51'` | 51 | MPACT 51-group |
| `'ECCO-33'` | 33 | ECCO 33-group |
| `'ECCO-1968'` | 1968 | ECCO 1968-group |

```python
groups = openmc.mgxs.EnergyGroups(
    openmc.mgxs.GROUP_STRUCTURES['CASMO-70']
)
```

### `openmc.mgxs.convert_flux_groups(flux, in_groups, out_groups)`

Convert a flux spectrum between two energy group structures.

---

## MGXS Classes

All MGXS classes share a common interface. Instantiate with `domain` (Cell, Material, Universe, or Mesh) and `groups`.

### Common methods and attributes

| Method/Attribute | Description |
|---|---|
| `.tallies` | Dict of OpenMC `Tally` objects needed to compute this MGXS |
| `.load_from_statepoint(sp)` | Load tally results from a `StatePoint` |
| `.get_xs(groups, subdomains, nuclides, xs_type, order_groups, value)` | Return MGXS array |
| `.get_pandas_dataframe(groups, subdomains, nuclides, xs_type, distribcell_paths)` | Return results as DataFrame |
| `.export_xs_data(filename, domain_type, xs_type, nuclides, row_column)` | Export to HDF5 |
| `.print_xs(subdomains, nuclides, xs_type)` | Print table |
| `.domain` | The spatial domain (Cell, Material, etc.) |
| `.groups` | `EnergyGroups` instance |
| `.energy_groups` | Same as `.groups` |
| `.xs_type` | `'macro'` or `'micro'` |
| `.nuclides` | List of nuclides to score (default `'all'`) |
| `.by_nuclide` | If True, compute per-nuclide cross sections |

### Scalar cross sections (1D: energy group → value)

| Class | Description |
|---|---|
| `openmc.mgxs.TotalXS` | Total cross section Σ_t |
| `openmc.mgxs.TransportXS` | Transport-corrected total Σ_tr |
| `openmc.mgxs.AbsorptionXS` | Absorption Σ_a |
| `openmc.mgxs.CaptureXS` | Radiative capture Σ_c |
| `openmc.mgxs.FissionXS` | Fission Σ_f |
| `openmc.mgxs.KappaFissionXS` | Recoverable fission energy κΣ_f |
| `openmc.mgxs.ScatterXS` | Scattering Σ_s |
| `openmc.mgxs.ReducedAbsorptionXS` | Reduced absorption |
| `openmc.mgxs.Chi` | Fission spectrum χ (fraction of fission neutrons born in each group) |
| `openmc.mgxs.InverseVelocity` | 1/v (used for time-dependent problems) |
| `openmc.mgxs.DiffusionCoefficient` | Diffusion coefficient D |
| `openmc.mgxs.Current` | Partial current |
| `openmc.mgxs.ArbitraryXS(reaction)` | Any MT reaction type |

### Matrix cross sections (2D: incoming × outgoing group)

| Class | Description |
|---|---|
| `openmc.mgxs.ScatterMatrixXS` | Scattering transfer matrix Σ_s(g→g'). Supports Legendre moments via `legendre_order`. |
| `openmc.mgxs.NuFissionMatrixXS` | Fission production matrix νΣ_f χ |
| `openmc.mgxs.MultiplicityMatrixXS` | Scattering multiplicity matrix |
| `openmc.mgxs.ScatterProbabilityMatrix` | Group-to-group scattering probabilities |
| `openmc.mgxs.ArbitraryMatrixXS(reaction)` | Matrix form for any MT reaction |
| `openmc.mgxs.MeshSurfaceMGXS` | MGXS on mesh surfaces |

---

## Multi-Delayed-Group Cross Sections (MDGXS)

For kinetics calculations; require a delayed group structure in addition to energy groups.

| Class | Description |
|---|---|
| `openmc.mgxs.ChiDelayed` | Delayed fission spectrum χ_d(g, delay_group) |
| `openmc.mgxs.DelayedNuFissionXS` | Delayed ν·Σ_f per delayed group |
| `openmc.mgxs.DelayedNuFissionMatrixXS` | Delayed fission production matrix |
| `openmc.mgxs.Beta` | Delayed neutron fraction β per group |
| `openmc.mgxs.DecayRate` | Decay constant λ per delayed precursor group |

```python
chi_d = openmc.mgxs.ChiDelayed(domain=fuel_cell, groups=groups, num_delayed_groups=6)
```

---

## Library

### `openmc.mgxs.Library(geometry, by_nuclide, mgxs_types, energy_groups, num_delayed_groups, correction, xs_type, legendre_order, name)`

Aggregates many MGXS objects across all domains in a geometry into a single library.

```python
lib = openmc.mgxs.Library(geometry)
lib.energy_groups = groups
lib.mgxs_types = ['total', 'absorption', 'nu-fission', 'fission', 'chi', 'scatter matrix']
lib.by_nuclide = False
lib.correction = 'P0'           # transport correction
lib.legendre_order = 3          # scattering anisotropy order
lib.build_library()

# Add tallies to model
tallies = openmc.Tallies()
lib.add_to_tallies_file(tallies, merge=True)

# After running, load results
with openmc.StatePoint('statepoint.h5') as sp:
    lib.load_from_statepoint(sp)

lib.export_to_hdf5('mgxs_library.h5')
```

Key attributes: `all_mgxs` (dict of all MGXS objects organized by domain), `domains` (list of domains).
