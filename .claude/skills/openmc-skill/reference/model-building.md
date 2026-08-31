# OpenMC Model Building (`openmc.model`)

API docs: https://docs.openmc.org/en/stable/pythonapi/model.html

The `openmc.model` module provides convenience functions and classes for assembling complete reactor models more easily than using raw geometry primitives.

---

## Model Container

### `openmc.model.Model`

Bundles geometry, materials, settings, and tallies into a single object that can be exported and run.

```python
import openmc
import openmc.model

model = openmc.model.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies,
)
model.run()
model.export_to_xml()      # writes all XML input files
```

Key methods:
- `model.run(**kwargs)` — runs simulation, accepts same args as `openmc.run()`
- `model.export_to_xml(directory='.')` — exports all input XMLs
- `model.import_properties(filename)` — import material properties from HDF5
- `model.export_properties(filename)` — export material properties

---

## Convenience Functions

### `openmc.model.borated_water(boron_ppm, temperature, pressure)`

Returns an `openmc.Material` with the composition of boron dissolved in water.

```python
water = openmc.model.borated_water(boron_ppm=900, temperature=600.0, pressure=15.5e6)
```

### `openmc.model.pin(surfaces, materials, subdivisions, divide_vols)`

Convenience function for building a fuel pin composed of concentric cylinders.

```python
# surfaces: list of ZCylinder, materials: list of Material (one per region)
pin_universe = openmc.model.pin(
    [fuel_or, clad_ir, clad_or],
    [fuel, void, clad, water]
)
```

### `openmc.model.subdivide(surfaces)`

Create regions separated by a series of surfaces. Returns a list of `Region` objects.

---

## Composite Surfaces

Composite surfaces are high-level geometric shapes built from multiple primitive surfaces. They return a region (intersection of halfspaces) and expose individual surface components.

| Class | Description |
|---|---|
| `openmc.model.RectangularParallelepiped(xmin, xmax, ymin, ymax, zmin, zmax, boundary_type)` | Box bounded by 6 planes. |
| `openmc.model.RectangularPrism(width, height, axis, origin, boundary_type)` | Infinite rectangular prism (4 planar surfaces). |
| `openmc.model.OrthogonalBox` | Arbitrarily oriented orthogonal box. |
| `openmc.model.HexagonalPrism(edge_length, orientation, origin, boundary_type)` | Hexagonal prism (6 planes). `orientation`: `'x'` or `'y'`. |
| `openmc.model.IsogonalOctagon(origin, r1, r2, axis)` | Infinite isogonal octagon prism. |
| `openmc.model.CruciformPrism` | Generalized cruciform prism. |
| `openmc.model.RightCircularCylinder(center_base, height, radius, axis, boundary_type)` | Finite cylinder (cylinder + 2 endcap planes). |
| `openmc.model.ConicalFrustum(z0, z1, r0, r1, axis)` | Frustum (truncated cone). |
| `openmc.model.Vessel(outer_radius, inner_radius, height)` | Cylinder with semi-ellipsoid top and bottom caps. |
| `openmc.model.Polygon(points, basis)` | Polygon formed from closed path of points. |
| `openmc.model.XConeOneSided`, `openmc.model.YConeOneSided`, `openmc.model.ZConeOneSided` | One-sided (open) cones parallel to respective axes. |
| `openmc.model.CylinderSector(r_inner, r_outer, phi0, phi1, axis)` | Cylindrical sector (annular wedge). |

```python
# Example: hex prism containing a fuel pin
hex_prism = openmc.model.HexagonalPrism(edge_length=0.63, orientation='x', boundary_type='reflective')
universe = openmc.Universe(cells=[
    openmc.Cell(region=-fuel_cyl, fill=fuel),
    openmc.Cell(region=+fuel_cyl & -hex_prism, fill=water),
])
```

---

## TRISO Fuel Modeling

TRISO (Tristructural-Isotopic) micro fuel particles consist of a fuel kernel surrounded by buffer, IPyC, SiC, and OPyC layers.

### `openmc.model.TRISO`

Represents a single TRISO particle. Encapsulates the concentric spherical layers.

```python
# Define layer materials: kernel, buffer, IPyC, SiC, OPyC
triso = openmc.model.TRISO(
    outer_radius=0.04,       # outermost layer radius [cm]
    fill=triso_univ,         # Universe describing concentric layers
    center=(0.0, 0.0, 0.0)
)
```

### `openmc.model.pack_spheres(radius, region, pf, initial_packing, seed)`

Generate a random, non-overlapping configuration of sphere centers within a container region.

```python
centers = openmc.model.pack_spheres(
    radius=0.04,
    region=-container_sphere,
    pf=0.35,         # packing fraction
    initial_packing=0.2,
)
```

Returns an array of sphere center coordinates.

### `openmc.model.create_triso_lattice(trisos, lower_left, pitch, shape, background)`

Creates an optimized lattice of TRISO particles for efficient particle tracking.

```python
trisos = [openmc.model.TRISO(outer_radius, fill=triso_univ, center=c) for c in centers]
triso_lat = openmc.model.create_triso_lattice(
    trisos,
    lower_left=(-r, -r, -r),
    pitch=(0.09, 0.09, 0.09),
    shape=(n, n, n),
    background=matrix_material,
)
```

### Full TRISO pebble example

```python
import openmc
import openmc.model

# Layer materials
kernel = openmc.Material(name='UO2 kernel')
kernel.add_nuclide('U235', 0.04); kernel.add_element('O', 2.0)
kernel.set_density('g/cm3', 10.5)

buffer = openmc.Material(name='Buffer')
buffer.add_element('C', 1.0); buffer.set_density('g/cm3', 1.0)

# ... define IPyC, SiC, OPyC materials similarly

# TRISO layer spheres
s1 = openmc.Sphere(r=0.0213)
s2 = openmc.Sphere(r=0.0313)
s3 = openmc.Sphere(r=0.0350)
s4 = openmc.Sphere(r=0.0385)
s5 = openmc.Sphere(r=0.0420)

triso_cells = [
    openmc.Cell(fill=kernel, region=-s1),
    openmc.Cell(fill=buffer, region=+s1 & -s2),
    # ... etc.
]
triso_univ = openmc.Universe(cells=triso_cells)

# Pack spheres into pebble
pebble_sphere = openmc.Sphere(r=3.0)
centers = openmc.model.pack_spheres(radius=0.0420, region=-pebble_sphere, pf=0.35)
trisos = [openmc.model.TRISO(0.0420, triso_univ, c) for c in centers]
lat = openmc.model.create_triso_lattice(trisos, ...)
```
