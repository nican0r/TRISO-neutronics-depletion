# OpenMC Basic Functionality (`openmc`)

API docs: https://docs.openmc.org/en/stable/pythonapi/base.html

---

## Materials

| Class/Function | Description |
|---|---|
| `openmc.Material(material_id, name, temperature)` | A material composed of nuclides/elements. Call `.add_nuclide()`, `.add_element()`, `.set_density()`, `.add_s_alpha_beta()`. |
| `openmc.Materials(materials)` | Collection written to `materials.xml`. Supports `export_to_xml()`. |
| `openmc.plot_xs(this, types, divisor_types, temperature, ...)` | Plot continuous-energy cross sections for a nuclide, element, or material. |

### Material key methods
```python
mat = openmc.Material()
mat.add_nuclide('U235', 0.04)          # by atom/weight fraction
mat.add_element('O', 2.0)
mat.set_density('g/cm3', 10.5)
mat.add_s_alpha_beta('c_H_in_H2O')     # thermal scattering
mat.temperature = 600.0                 # K
mat.depletable = True                   # flag for depletion
```

---

## Surfaces

All surfaces support Boolean CSG operators (`-` for negative halfspace, `+` for positive).

| Class | Form |
|---|---|
| `openmc.Plane(a, b, c, d)` | Ax + By + Cz = D |
| `openmc.XPlane(x0)` | x = x₀ |
| `openmc.YPlane(y0)` | y = y₀ |
| `openmc.ZPlane(z0)` | z = z₀ |
| `openmc.XCylinder(y0, z0, r)` | Cylinder ∥ x-axis |
| `openmc.YCylinder(x0, z0, r)` | Cylinder ∥ y-axis |
| `openmc.ZCylinder(x0, y0, r)` | Cylinder ∥ z-axis |
| `openmc.Sphere(x0, y0, z0, r)` | (x-x₀)²+(y-y₀)²+(z-z₀)² = r² |
| `openmc.XCone(x0, y0, z0, r2)` | Cone ∥ x-axis |
| `openmc.YCone(x0, y0, z0, r2)` | Cone ∥ y-axis |
| `openmc.ZCone(x0, y0, z0, r2)` | Cone ∥ z-axis |
| `openmc.Quadric(a,b,c,d,e,f,g,h,j,k)` | General quadric surface |
| `openmc.XTorus(x0,y0,z0,a,b,c)` | Torus ∥ x-axis |
| `openmc.YTorus(...)` | Torus ∥ y-axis |
| `openmc.ZTorus(...)` | Torus ∥ z-axis |

`boundary_type`: `'transmission'` (default), `'vacuum'`, `'reflective'`, `'periodic'`, `'white'`

### Regions (CSG)

| Class | Operator |
|---|---|
| `openmc.Halfspace` | Returned by `+surface` / `-surface` |
| `openmc.Intersection` | `region_a & region_b` |
| `openmc.Union` | `region_a \| region_b` |
| `openmc.Complement` | `~region` |
| `openmc.BoundingBox` | Axis-aligned bounding box of a region |

---

## Geometry

| Class | Description |
|---|---|
| `openmc.Cell(cell_id, name, fill, region)` | A region of space filled with a material or universe. |
| `openmc.Universe(universe_id, name, cells)` | Collection of cells; can be repeated in lattices. |
| `openmc.DAGMCUniverse(filename)` | References a DAGMC CAD-based geometry file. |
| `openmc.RectLattice(lattice_id, name)` | Rectangular lattice. Set `.pitch`, `.lower_left`, `.universes`. |
| `openmc.HexLattice(lattice_id, name)` | Hexagonal lattice. Set `.pitch`, `.center`, `.universes`, `.orientation`. |
| `openmc.Geometry(root)` | Top-level geometry container. `root` is a Universe or Cell. |

```python
cyl = openmc.ZCylinder(r=0.4)
cell = openmc.Cell(fill=fuel, region=-cyl)
universe = openmc.Universe(cells=[cell])
geometry = openmc.Geometry(universe)
```

---

## Tallies

### Filters

| Filter | Bins on |
|---|---|
| `openmc.CellFilter(bins)` | Cell ID |
| `openmc.MaterialFilter(bins)` | Material ID |
| `openmc.SurfaceFilter(bins)` | Surface ID |
| `openmc.EnergyFilter(bins)` | Incident energy (eV) |
| `openmc.EnergyoutFilter(bins)` | Outgoing energy (eV) |
| `openmc.MeshFilter(mesh)` | Mesh element |
| `openmc.MeshSurfaceFilter(mesh)` | Mesh surface crossings |
| `openmc.ParticleFilter(bins)` | Particle type (`'neutron'`, `'photon'`) |
| `openmc.MuFilter(bins)` | Cosine of scattering angle |
| `openmc.PolarFilter(bins)` | Polar angle |
| `openmc.AzimuthalFilter(bins)` | Azimuthal angle |
| `openmc.DistribcellFilter(cell)` | Repeated cell instances |
| `openmc.DelayedGroupFilter(bins)` | Delayed neutron precursor groups |
| `openmc.LegendreFilter(order)` | Legendre expansion moments |
| `openmc.SphericalHarmonicsFilter(order)` | Spherical harmonic moments |
| `openmc.TimeFilter(bins)` | Particle time |
| `openmc.UniverseFilter(bins)` | Universe ID |
| `openmc.CellBornFilter(bins)` | Cell where particle was born |
| `openmc.CellInstanceFilter(bins)` | Cell instance |
| `openmc.EnergyFunctionFilter(energy, y)` | Multiplies score by f(E) |
| `openmc.ReactionFilter(bins)` | MT reaction number |

### Tally scores (common)

`flux`, `total`, `absorption`, `scatter`, `fission`, `nu-fission`, `kappa-fission`, `current`, `events`, `heating`, `heating-local`, `damage-energy`, `(n,p)`, `(n,a)`, `H1-production`

### Tally classes

```python
tally = openmc.Tally(name='flux tally')
tally.filters = [openmc.CellFilter([cell.id]), openmc.EnergyFilter([0, 0.625e-6, 20.0])]
tally.scores = ['flux']
tally.nuclides = ['U235', 'U238']   # optional: score per nuclide

tallies = openmc.Tallies([tally])
```

| Class | Description |
|---|---|
| `openmc.Tally` | Single tally definition |
| `openmc.Tallies` | Collection, written to `tallies.xml` |
| `openmc.Trigger(trigger_type, threshold)` | Stop simulation on uncertainty criterion (`'std_dev'`, `'rel_err'`, `'variance'`) |
| `openmc.TallyDerivative` | Perturbation derivative applied to a tally |

---

## Meshes

| Class | Description |
|---|---|
| `openmc.RegularMesh()` | Uniform Cartesian mesh. Set `.dimension`, `.lower_left`, `.upper_right` or `.width`. |
| `openmc.RectilinearMesh()` | Non-uniform Cartesian mesh. Set `.x_grid`, `.y_grid`, `.z_grid`. |
| `openmc.CylindricalMesh()` | Cylindrical (r, φ, z). Set `.r_grid`, `.phi_grid`, `.z_grid`. |
| `openmc.SphericalMesh()` | Spherical (r, θ, φ). Set `.r_grid`, `.theta_grid`, `.phi_grid`. |
| `openmc.UnstructuredMesh(filename, library)` | Unstructured mesh from file (LibMesh or MOAB). |

---

## Settings

```python
settings = openmc.Settings()
settings.run_mode = 'eigenvalue'        # or 'fixed source', 'volume', 'plot'
settings.batches = 150
settings.inactive = 30
settings.particles = 10000
settings.seed = 42
settings.temperature = {'method': 'interpolation', 'range': (300, 1200)}
settings.photon_transport = True
settings.source = [src1, src2]
```

Key attributes: `batches`, `inactive`, `particles`, `seed`, `run_mode`, `source`, `temperature`, `energy_mode` (`'continuous-energy'` or `'multi-group'`), `cutoff`, `resonance_scattering`, `weight_windows`, `survival_biasing`, `max_tracks`.

---

## Sources

| Class | Description |
|---|---|
| `openmc.IndependentSource(space, angle, energy, time, strength, particle)` | Phase-space distribution source. `space` is `openmc.stats.Spatial`, `angle` is `openmc.stats.UnitSphere`, `energy` is `openmc.stats.Univariate`. |
| `openmc.FileSource(path, strength)` | Source from HDF5 particle file. |
| `openmc.CompiledSource(library, parameters)` | Source from compiled `.so` shared library. |
| `openmc.MeshSource(mesh, sources)` | Spatially distributed source over mesh elements. |
| `openmc.TokamakSource(...)` | Neutron emission from tokamak plasma (DT fusion). |

```python
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box([-1,-1,-1], [1,1,1]),
    angle=openmc.stats.Isotropic(),
    energy=openmc.stats.Watt(a=0.988e6, b=2.249e-6),
)
```

---

## Running OpenMC

| Function | Description |
|---|---|
| `openmc.run(particles, threads, geometry_debug, restart_file, tracks, output, cwd, openmc_exec, mpi_args, event_based)` | Run transport simulation. |
| `openmc.calculate_volumes(threads, output, cwd, openmc_exec, mpi_args)` | Run stochastic volume calculation. |
| `openmc.plot_geometry(threads, output, cwd, openmc_exec)` | Run in plotting mode. |
| `openmc.plot_inline(plots, openmc_exec)` | Display geometry plots inline (Jupyter). |
| `openmc.search_for_keff(model_builder, target, tol, bracketed_range, guesses, maxiter, print_iterations, run_args)` | Binary search to find parameter value giving target k-eff. |

---

## Post-processing

| Class/Function | Description |
|---|---|
| `openmc.StatePoint(filename)` | Read simulation state; access tally results, k-eff, source. |
| `openmc.Summary(filename)` | Read model summary (geometry, materials). |
| `openmc.Tracks(filename)` | Collection of particle tracks. |
| `openmc.ParticleTrack` | Single particle track. |
| `openmc.read_source_file(path)` | Read HDF5 source file → list of `SourceParticle`. |
| `openmc.write_source_file(source_particles, path)` | Write source particles to HDF5. |
| `openmc.read_collision_track_file(path)` | Read HDF5 or MCPL collision track file. |
| `openmc.voxel_to_vtk(voxel_file, output)` | Convert voxel HDF5 to VTK for visualization. |

### Reading tally results
```python
with openmc.StatePoint('statepoint.100.h5') as sp:
    tally = sp.get_tally(name='flux tally')
    df = tally.get_pandas_dataframe()
    mean = tally.mean          # shape: (filter_bins, nuclides, scores)
    std_dev = tally.std_dev
```

---

## Geometry Plotting

| Class | Description |
|---|---|
| `openmc.SlicePlot` | 2D slice of geometry. Set `.origin`, `.width`, `.basis`, `.pixels`, `.color_by`. |
| `openmc.VoxelPlot` | 3D voxel geometry. Set `.origin`, `.width`, `.pixels`. |
| `openmc.WireframeRayTracePlot` | Wireframe rendering. |
| `openmc.SolidRayTracePlot` | Phong-shaded solid rendering. |
| `openmc.Plots` | Collection of plot definitions, written to `plots.xml`. |

---

## Variance Reduction

| Class/Function | Description |
|---|---|
| `openmc.WeightWindows(mesh, lower_ww_bounds, upper_ww_bounds, particle_type, energy_bounds, survival_ratio, max_lower_bound_ratio, weight_cutoff, max_split)` | Mesh-based weight windows for importance sampling. |
| `openmc.WeightWindowGenerator(mesh, energy_bounds, particle_type, method, max_realizations, update_interval, on_the_fly)` | Automatic weight window generation during simulation. |
| `openmc.wwinp_to_wws(path)` | Load weight windows from MCNP wwinp format. |
| `openmc.hdf5_to_wws(path)` | Load weight windows from HDF5 file. |

---

## Tally Arithmetic (post-processing)

Access via `openmc.arithmetic.*` or through tally operators (`+`, `-`, `*`, `/`).

| Class | Description |
|---|---|
| `CrossScore` | Encapsulates all score combinations across tallies. |
| `CrossFilter` | All filter combinations. |
| `AggregateScore` | Aggregate of scores. |
| `AggregateFilter` | Aggregate of filter bins. |

### Functional expansion reconstruction

| Class/Function | Description |
|---|---|
| `openmc.ZernikeRadial(coeffs, radius)` | Evaluate radial Zernike polynomial from expansion coefficients. |
| `openmc.legendre_from_expcoef(coeffs, domain)` | Return `numpy.polynomial.Legendre` from expansion coefficients. |

---

## Multi-group data (MG mode)

| Class | Description |
|---|---|
| `openmc.XSdata(name, energy_groups)` | A multi-group cross section dataset. |
| `openmc.MGXSLibrary(energy_groups)` | Multi-group cross section library for MG simulation. |
