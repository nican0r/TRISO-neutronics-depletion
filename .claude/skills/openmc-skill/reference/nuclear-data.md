# OpenMC Nuclear Data Interface (`openmc.data`)

API docs: https://docs.openmc.org/en/stable/pythonapi/data.html

`openmc.data` provides a Python interface to nuclear data files (ACE, ENDF, NJOY). Use it to read, process, and inspect cross sections, decay data, fission yields, thermal scattering, and more.

---

## Quick-reference utility functions

These functions are available directly in `openmc.data` and require no file loading.

| Function | Returns | Example |
|---|---|---|
| `openmc.data.atomic_mass(isotope)` | Mass in amu | `openmc.data.atomic_mass('U235')` |
| `openmc.data.atomic_weight(element)` | Weight in amu | `openmc.data.atomic_weight('U')` |
| `openmc.data.half_life(isotope)` | Half-life in s (or `None`) | `openmc.data.half_life('Cs137')` |
| `openmc.data.decay_constant(isotope)` | Decay constant in s⁻¹ | `openmc.data.decay_constant('I131')` |
| `openmc.data.decay_energy(nuclide)` | Mean decay energy | `openmc.data.decay_energy('Co60')` |
| `openmc.data.decay_photon_energy(nuclide)` | Photon energy distribution from decay | — |
| `openmc.data.isotopes(element)` | List of `(isotope, abundance)` | `openmc.data.isotopes('U')` |
| `openmc.data.gnds_name(nuclide)` | GNDS-convention nuclide name | `openmc.data.gnds_name('U235')` → `'U235'` |
| `openmc.data.zam(nuclide)` | `(Z, A, metastable)` tuple | `openmc.data.zam('Am242_m1')` |
| `openmc.data.water_density(temperature, pressure)` | Water density in g/cm³ | `openmc.data.water_density(600, 15.5e6)` |
| `openmc.data.dose_coefficients(particle, geometry)` | Dose conversion coefficients | — |
| `openmc.data.mass_attenuation_coefficient(element, energy)` | μ/ρ as function of energy | — |
| `openmc.data.mass_energy_absorption_coefficient(element, energy)` | μ_en/ρ as function of energy | — |
| `openmc.data.linearize(function, tolerance)` | Tabulate a callable function | — |
| `openmc.data.thin(x, y, tolerance)` | Remove redundant (x,y) points | — |
| `openmc.data.combine_distributions(distributions, probabilities)` | Mix distributions | — |
| `openmc.data.kalbach_slope(Ea, Za, Aa, Eb, Zb, Ab, Ec, Zc, Ac)` | Kalbach-Mann slope parameter | — |

---

## Data Library

### `openmc.data.DataLibrary`

Manages a collection of cross section data files (nuclear data library).

```python
lib = openmc.data.DataLibrary.from_hdf5('cross_sections.xml')

# Find a specific nuclide
entry = lib.get_by_material('U235')
print(entry['path'])

# Register a new file
lib.register_file('my_nuclide.h5')
lib.export_to_xml('cross_sections_new.xml')
```

---

## Incident Neutron Data

### `openmc.data.IncidentNeutron`

Continuous-energy neutron interaction data for a single nuclide/isotope.

```python
u235 = openmc.data.IncidentNeutron.from_hdf5('U235.h5')
u235 = openmc.data.IncidentNeutron.from_ace('U235.ace')
u235 = openmc.data.IncidentNeutron.from_endf('U235.endf')

# Access reactions by MT number
elastic = u235[2]              # elastic scattering
fission = u235[18]             # fission
capture = u235[102]            # (n,γ)

# Cross section as function of energy
import numpy as np
energies = np.logspace(-5, 7, 1000)   # eV
xs = u235[18].xs['0K'](energies)      # fission XS at 0K

# Resonances
print(u235.resonances)

# Export
u235.export_to_hdf5('U235.h5')
```

### `openmc.data.Reaction`

A single nuclear reaction (identified by MT number).

| Attribute | Description |
|---|---|
| `mt` | MT reaction number |
| `q_value` | Q-value in eV |
| `xs` | Dict of temperature → cross section function |
| `products` | List of `Product` objects |

### `openmc.data.Product`

Secondary particle emitted in a reaction.

| Attribute | Description |
|---|---|
| `particle` | Particle type |
| `emission_mode` | `'prompt'`, `'delayed'`, `'total'` |
| `decay_rate` | Decay rate for delayed products |
| `yield_` | Multiplicity as function of energy |
| `distribution` | `AngleEnergy` distribution |

### `openmc.data.FissionEnergyRelease`

Energy components released by fission reactions (prompt, delayed, recoverable).

---

## Decay Data

### `openmc.data.Decay`

Radioactive decay data for a nuclide.

```python
cs137 = openmc.data.Decay.from_endf('Cs137.endf')
cs137 = openmc.data.Decay.from_hdf5('Cs137.h5')

print(cs137.half_life)           # seconds
print(cs137.modes)               # list of DecayMode
print(cs137.nuclide)             # parent nuclide name

# Decay photon energy distribution
energy_dist = cs137.photon_energy
```

---

## Fission Product Yields

### `openmc.data.FissionProductYields`

Independent and cumulative fission product yields.

```python
u235_fp = openmc.data.FissionProductYields.from_endf('U235.endf')
u235_fp = openmc.data.FissionProductYields.from_hdf5('U235.h5')

# Independent yields at a given energy
yields = u235_fp.independent[0.0253]    # thermal (eV)
print(yields['Xe135'])

# Cumulative yields
cum = u235_fp.cumulative[14.0e6]        # fast
```

---

## Thermal Scattering

### `openmc.data.ThermalScattering`

S(α,β) thermal scattering law data (for materials like H₂O, graphite, ZrH).

```python
h2o = openmc.data.ThermalScattering.from_hdf5('c_H_in_H2O.h5')
h2o = openmc.data.ThermalScattering.from_ace('tsl-HinH2O.ace')

print(h2o.nuclides)       # applicable nuclides
h2o.export_to_hdf5('h2o_tsl.h5')
```

| Class | Description |
|---|---|
| `openmc.data.CoherentElastic` | Coherent elastic scattering from crystalline material |
| `openmc.data.IncoherentElastic` | Incoherent elastic scattering |

---

## Photon Data

### `openmc.data.IncidentPhoton`

Photon interaction data for an element.

```python
pb = openmc.data.IncidentPhoton.from_hdf5('Pb.h5')
pb = openmc.data.IncidentPhoton.from_endf('Pb_photo.endf')

# Coherent scattering cross section
xs = pb[502].xs(energies)
```

### `openmc.data.AtomicRelaxation`

Atomic relaxation data (fluorescence, Auger electrons).

---

## Resonance Data

| Class | Description |
|---|---|
| `openmc.data.Resonances` | Container for resolved and unresolved resonance data |
| `openmc.data.ResonanceRange` | A single resolved resonance energy range |
| `openmc.data.SingleLevelBreitWigner` | SLBW resolved resonance formalism |
| `openmc.data.MultiLevelBreitWigner` | MLBW resolved resonance formalism |
| `openmc.data.ReichMoore` | Reich-Moore resolved resonance formalism |
| `openmc.data.RMatrixLimited` | R-matrix limited formalism |
| `openmc.data.Unresolved` | Unresolved resonance parameters |
| `openmc.data.WindowedMultipole` | Windowed multipole representation for Doppler broadening on-the-fly |
| `openmc.data.ProbabilityTables` | Unresolved resonance probability tables |

---

## 1D Function Classes

Used internally to represent cross sections and other energy-dependent quantities.

| Class | Description |
|---|---|
| `openmc.data.Tabulated1D(x, y, breakpoints, interpolation)` | Tabulated function with interpolation regions |
| `openmc.data.Polynomial(coefficients)` | Power series |
| `openmc.data.Sum(functions)` | Sum of multiple functions |
| `openmc.data.Combination(functions, operators)` | Arbitrary combination |
| `openmc.data.Regions1D(functions, breakpoints)` | Piecewise composition |
| `openmc.data.ResonancesWithBackground` | XS = resonance + background |

All support calling: `f(energy_array)` returns array of values.

---

## Angle-Energy Distributions

| Class | Description |
|---|---|
| `openmc.data.UncorrelatedAngleEnergy(angle, energy)` | Independent angle and energy distributions |
| `openmc.data.CorrelatedAngleEnergy` | Correlated (E', μ) distribution |
| `openmc.data.KalbachMann` | Kalbach-Mann pre-equilibrium model |
| `openmc.data.NBodyPhaseSpace(total_mass, n_particles)` | N-body phase space |
| `openmc.data.LaboratoryAngleEnergy` | Lab-frame angle-energy |
| `openmc.data.MaxwellEnergy(theta)` | Maxwell fission spectrum |
| `openmc.data.Watt(a, b)` | Watt fission spectrum |
| `openmc.data.Evaporation(theta)` | Evaporation spectrum |
| `openmc.data.MadlandNix(e_f_light, e_f_heavy, tm)` | Madland-Nix prompt fission spectrum |
| `openmc.data.ContinuousTabular` | Tabular (E, μ, p) distribution |
| `openmc.data.DiscretePhoton(primary_flag, energy, atomic_weight_ratio)` | Discrete photon emission |
| `openmc.data.LevelInelastic(threshold, mass_ratio)` | Discrete level inelastic scattering |

---

## ACE Format Interface

### `openmc.data.ace.Library(filename)`

```python
lib = openmc.data.ace.Library('endfb8-neutron.ace')
print(lib.tables)            # list of Table objects
t = lib.tables[0]
```

### `openmc.data.ace.Table`

| Attribute | Description |
|---|---|
| `name` | ZAID identifier (e.g., `'92235.80c'`) |
| `atomic_weight_ratio` | AWR |
| `temperature` | Temperature in MeV |
| `pairs` | (NXS, JXS) header arrays |
| `xss` | Cross section data array |

### Functions

| Function | Description |
|---|---|
| `openmc.data.ace.ascii_to_binary(ascii_file, binary_file)` | Convert ACE type-1 ASCII to binary |
| `openmc.data.ace.get_libraries_from_xsdir(path)` | Parse MCNP xsdir → list of ACE file paths |
| `openmc.data.ace.get_libraries_from_xsdata(path)` | Parse Serpent xsdata → list of ACE file paths |

---

## ENDF Format Interface

### `openmc.data.endf.Evaluation(filename)`

```python
ev = openmc.data.endf.Evaluation('U235.endf')
mf2 = ev.section[2, 151]     # access by (MF, MT)
```

### Functions

| Function | Description |
|---|---|
| `openmc.data.endf.get_evaluations(endf_file)` | List all material evaluations in an ENDF file |
| `openmc.data.endf.get_head_record(file_obj)` | Parse HEAD record |
| `openmc.data.endf.get_cont_record(file_obj)` | Parse CONT record |
| `openmc.data.endf.get_tab1_record(file_obj)` | Parse TAB1 record → `Tabulated1D` |
| `openmc.data.endf.get_tab2_record(file_obj)` | Parse TAB2 record |
| `openmc.data.endf.get_text_record(file_obj)` | Parse TEXT record |
| `openmc.data.endf.float_endf(string)` | Convert ENDF-format string to float |

---

## NJOY Interface

Process raw ENDF files into ACE or PENDF format using NJOY.

| Function | Description |
|---|---|
| `openmc.data.njoy.run(commands, tapein, tapeout, input_filename, stdout, njoy_exec)` | Run NJOY with given commands |
| `openmc.data.njoy.make_pendf(endf_file, pendf_file, temperatures, error, stdout)` | Generate linearized/pointwise PENDF |
| `openmc.data.njoy.make_ace(endf_file, temperatures, acer_opts, broadr_opts, thermr_opts, gaspr_opts, heatr_opts, purr_opts, stderr, evaluation)` | Generate ACE file from ENDF |
| `openmc.data.njoy.make_ace_thermal(endf_files, temperatures, ...)` | Generate thermal scattering ACE |

```python
openmc.data.njoy.make_ace(
    endf_file='U235.endf',
    temperatures=[293.6, 600.0, 900.0],
)
```
