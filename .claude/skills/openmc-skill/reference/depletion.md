# OpenMC Depletion (`openmc.deplete`)

API docs: https://docs.openmc.org/en/stable/pythonapi/deplete.html

The `openmc.deplete` module solves the Bateman equations for fuel burnup. It couples OpenMC transport with depletion solvers (CRAM) to evolve material compositions over time.

---

## Minimal depletion workflow

```python
import openmc
import openmc.deplete

# 1. Build model (mark materials as depletable)
fuel.depletable = True

# 2. Choose transport operator
operator = openmc.deplete.CoupledOperator(
    model,
    chain_file='/path/to/chain.xml',
)

# 3. Choose integrator and run
timesteps = [30, 30, 30, 30]       # days
power = 174e6                       # watts (or use power_density)

integrator = openmc.deplete.PredictorIntegrator(
    operator, timesteps, power, timestep_units='d'
)
integrator.integrate()

# 4. Read results
results = openmc.deplete.Results('depletion_results.h5')
time, keff = results.get_keff()
time, atoms = results.get_atoms(fuel, 'Xe135')
```

---

## Transport Operators

### `openmc.deplete.CoupledOperator(model, chain_file, prev_results, diff_burnable_mats, normalization_mode, fission_yield_mode, fission_yield_opts, reaction_rate_opts)`

Transport-coupled operator. Runs OpenMC transport at each depletion step to get reaction rates.

| Parameter | Description |
|---|---|
| `model` | `openmc.model.Model` instance |
| `chain_file` | Path to depletion chain XML file |
| `prev_results` | Path to previous results for restart |
| `diff_burnable_mats` | Split identical materials into unique instances |
| `normalization_mode` | `'energy-deposition'`, `'fission-q'`, or `'source-rate'` |
| `fission_yield_mode` | `'constant'`, `'cutoff'`, `'average'`, or `'tallied'` |

### `openmc.deplete.IndependentOperator(materials, fluxes, micros, chain_file, keff, normalization_mode, fission_yield_mode, prev_results)`

Transport-independent operator: uses precomputed one-group fluxes and cross sections. Suitable for parametric studies or when coupling to another transport code.

```python
micros, flux = openmc.deplete.get_microxs_and_flux(model, domains)
operator = openmc.deplete.IndependentOperator(
    openmc.Materials([fuel]),
    [flux],
    [micros],
    chain_file='chain.xml',
)
```

---

## Integrators (time-stepping algorithms)

All integrators share the same interface:

```python
integrator = openmc.deplete.CFIntegrator(operator, timesteps, power, timestep_units)
integrator.integrate()
```

| Class | Algorithm | Notes |
|---|---|---|
| `PredictorIntegrator` | First-order predictor (Euler forward) | Simple, fast, less accurate |
| `CECMIntegrator` | CE/CM (Predictor-Corrector, Midpoint) | Moderate accuracy |
| `CELIIntegrator` | CE/LI CFQ4 | Good accuracy |
| `CF4Integrator` | CF4 (Commutator-Free, 4-stage) | High accuracy |
| `EPCRK4Integrator` | EPC-RK4 | High accuracy |
| `LEQIIntegrator` | LE/QI CFQ4 | High accuracy |
| `SICELIIntegrator` | SI-CE/LI (Stochastic Implicit CE/LI) | For high-fidelity simulations |
| `SILEQIIntegrator` | SI-LE/QI (Stochastic Implicit LE/QI) | For high-fidelity simulations |

Common parameters:
- `timesteps`: list of time step lengths
- `power` or `power_density` (W/cm³) or `source_rates`
- `timestep_units`: `'s'`, `'min'`, `'h'`, `'d'`, `'MWd/kg'`

---

## Depletion Chain

### `openmc.deplete.Chain`

Full representation of a depletion chain (transmutation paths, decay modes, fission yields).

```python
chain = openmc.deplete.Chain.from_endf(endf_files)
chain = openmc.deplete.Chain.from_xml('chain.xml')
chain.export_to_xml('chain_out.xml')

# Inspect
print(chain.nuclides)          # list of Nuclide names
nuc = chain['U235']
print(nuc.reactions)           # list of ReactionTuple
print(nuc.decay_modes)         # list of DecayTuple
```

### `openmc.deplete.Nuclide`

Single nuclide in the chain.

| Attribute | Description |
|---|---|
| `name` | Nuclide name (e.g. `'U235'`) |
| `half_life` | Half-life in seconds |
| `decay_modes` | List of `DecayTuple` |
| `reactions` | List of `ReactionTuple` |
| `yield_data` | `FissionYieldDistribution` (if fissile) |

### Supporting data classes

| Class | Description |
|---|---|
| `DecayTuple` | `(type, target, branching_ratio)` |
| `ReactionTuple` | `(type, target, Q, branching_ratio)` |
| `FissionYieldDistribution` | Energy-dependent yield distributions for a parent nuclide |
| `FissionYield` | Mapping `{daughter: yield}` at a specific energy |

---

## Results

### `openmc.deplete.Results`

Read and query depletion simulation output.

```python
results = openmc.deplete.Results('depletion_results.h5')

# k-eff over time
time, keff = results.get_keff()                  # time in [s], keff array

# Atom count
time, atoms = results.get_atoms(material, 'Xe135')

# Activity
time, activity = results.get_activity(material, 'Cs137', units='Bq')

# Decay heat
time, heat = results.get_decay_heat(material, units='W')

# Nuclide concentrations at step n
atoms_dict = results[n].get_atoms(material)
```

### `openmc.deplete.StepResult`

Result from a single depletion timestep. Contains `k`, `rates`, `time`, `source_rate`, and `mat_to_ind`.

---

## CRAM Solvers

CRAM (Chebyshev Rational Approximation Method) solves the matrix exponential form of the Bateman equations.

| Class | Description |
|---|---|
| `openmc.deplete.cram.CRAM16` | 16th-order IPF CRAM |
| `openmc.deplete.cram.CRAM48` | 48th-order IPF CRAM (default, more accurate) |
| `openmc.deplete.cram.IPFCramSolver` | Base solver using incomplete partial fractions |

---

## Microscopic cross sections (for IndependentOperator)

### `openmc.deplete.get_microxs_and_flux(model, domains, nuclides, chain_file, run_kwargs)`

Runs a short OpenMC simulation to extract one-group reaction rates and fluxes for depletion domains.

```python
micros, flux = openmc.deplete.get_microxs_and_flux(
    model,
    domains=[fuel_cell],
    nuclides=['U235', 'U238', 'Pu239'],
    chain_file='chain.xml',
)
```

### `openmc.deplete.MicroXS`

Stores microscopic cross sections indexed by `(nuclide, reaction)`.

---

## Continuous Feed and Removal

### `openmc.deplete.transfer_rates.TransferRates(operator, model)`

Define continuous removal of nuclides from one material (and optional transfer to another).

```python
transfer = openmc.deplete.transfer_rates.TransferRates(operator, model)
transfer.set_transfer_rate(material, ['Xe135', 'Kr85'], 1e-5)  # 1/s removal rate
```

### `openmc.deplete.transfer_rates.ExternalSourceRates`

Define external source rates of nuclides into materials.

---

## Rigorous 2-Step (R2S) Method

### `openmc.deplete.R2SManager`

Manages the R2S workflow for shutdown dose rate calculations:
1. Forward neutron transport to get activation rates
2. Depletion to get activated material compositions
3. Photon source generation and transport

```python
r2s = openmc.deplete.R2SManager(model, chain_file='chain.xml', ...)
r2s.run_forward()
r2s.run_depletion()
r2s.run_photon_transport()
```

---

## Module-level variables

| Variable | Description |
|---|---|
| `openmc.deplete.comm` | MPI intercommunicator for OpenMC library calls |
| `openmc.deplete.pool.USE_MULTIPROCESSING` | Enable/disable multiprocessing for Bateman solver |
| `openmc.deplete.pool.NUM_PROCESSES` | Worker count for parallel depletion |
| `openmc.deplete.chain.REACTIONS` | Dict mapping reaction names → MT values |

---

## D1S Method (Direct 1-Step)

| Function | Description |
|---|---|
| `openmc.deplete.d1s.prepare_tallies(model, chain, ...)` | Prepare tallies for D1S shutdown dose calculation. |
| `openmc.deplete.d1s.time_correction_factors(chain, ...)` | Calculate time correction factors. |
| `openmc.deplete.d1s.apply_time_correction(tally, factors)` | Apply time correction to tally results. |
