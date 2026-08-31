# OpenMC Example Models (`openmc.examples`)

API docs: https://docs.openmc.org/en/stable/pythonapi/examples.html

The `openmc.examples` module provides prebuilt, ready-to-run reference models. All functions return an `openmc.model.Model` object.

---

## Simple Models

### `openmc.examples.pwr_pin_cell()`

Create a standard PWR pin-cell model with UO₂ fuel, helium gap, Zircaloy cladding, and borated water.

```python
import openmc.examples

model = openmc.examples.pwr_pin_cell()
model.run()
```

### `openmc.examples.slab_mg()`

Create a 1D slab model using multi-group cross sections.

```python
model = openmc.examples.slab_mg()
model.run()
```

### `openmc.examples.sphere_with_shielded_pocket()`

Continuous-energy deep-shielding model with a far detector pocket. Useful for variance reduction testing.

```python
model = openmc.examples.sphere_with_shielded_pocket()
model.run()
```

### `openmc.examples.random_ray_pin_cell()`

PWR pin-cell example using C5G7 multi-group cross section data, configured for random ray transport.

```python
model = openmc.examples.random_ray_pin_cell()
model.run()
```

### `openmc.examples.random_ray_three_region_cube_with_detectors()`

Three-region cube model with two external tally regions, configured for random ray transport.

```python
model = openmc.examples.random_ray_three_region_cube_with_detectors()
model.run()
```

---

## Reactor Assembly Models

### `openmc.examples.pwr_assembly()`

17×17 PWR fuel assembly with guide tubes, instrument tube, and borated water.

```python
model = openmc.examples.pwr_assembly()
model.run()
```

### `openmc.examples.pwr_core()`

Full-core PWR model built from assemblies.

```python
model = openmc.examples.pwr_core()
model.run()
```

---

## Usage pattern

All example models return an `openmc.model.Model`, so you can inspect and modify them before running:

```python
model = openmc.examples.pwr_pin_cell()

# Inspect settings
print(model.settings.batches)

# Modify settings before running
model.settings.particles = 50000
model.settings.batches = 200

# Add tallies
tally = openmc.Tally(name='fission rate')
tally.filters = [openmc.CellFilter([c.id for c in model.geometry.get_all_cells().values()])]
tally.scores = ['fission']
model.tallies.append(tally)

model.run()
```

---

## Additional Jupyter notebook examples

https://github.com/openmc-dev/openmc/wiki/Example-Jupyter-Notebooks
