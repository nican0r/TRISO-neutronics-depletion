# OpenMC Statistics (`openmc.stats`)

API docs: https://docs.openmc.org/en/stable/pythonapi/stats.html

`openmc.stats` provides probability distributions used to define particle source phase-space coordinates: energy, position, and direction.

---

## Usage pattern

Distributions are assigned to `openmc.IndependentSource`:

```python
import openmc

source = openmc.IndependentSource(
    space=openmc.stats.Box([-1,-1,-1], [1,1,1]),
    angle=openmc.stats.Isotropic(),
    energy=openmc.stats.Watt(a=0.988e6, b=2.249e-6),
)
settings.source = source
```

---

## Univariate (energy) distributions

These represent probability distributions of a single variable, most commonly used for energy.

### Continuous distributions

| Class | Description | Key parameters |
|---|---|---|
| `openmc.stats.Uniform(a, b)` | Uniform over [a, b] | `a`, `b` (eV) |
| `openmc.stats.PowerLaw(a, b, n)` | Power-law P∝xⁿ over [a,b] | `a`, `b`, `n` |
| `openmc.stats.Maxwell(theta)` | Maxwellian: P(E)∝√E·exp(-E/θ) | `theta` (eV) |
| `openmc.stats.Watt(a, b)` | Watt fission spectrum: P(E)∝exp(-E/a)·sinh(√(bE)) | `a` (eV), `b` (1/eV) |
| `openmc.stats.Normal(mean, std_dev, low, high)` | Normal (optionally truncated) | `mean`, `std_dev`, `low`, `high` |
| `openmc.stats.Tabular(x, p, interpolation)` | Piecewise continuous (tabulated). `interpolation`: `'histogram'` or `'linear-linear'`. | `x`, `p` |
| `openmc.stats.Legendre(coefficients)` | PDF given by Legendre polynomial expansion | `coefficients` |
| `openmc.stats.Mixture(probability, distribution)` | Mixture of distributions | `probability` (list), `distribution` (list) |
| `openmc.stats.DecaySpectrum(material)` | Energy distribution from decay photons of a nuclide mixture | `material` |

### Discrete distributions

| Class/Function | Description |
|---|---|
| `openmc.stats.Discrete(x, p)` | Probability mass function (list of values and probabilities) |
| `openmc.stats.delta_function(x)` | Convenience: `Discrete([x], [1.0])` |

### Convenience functions

| Function | Description |
|---|---|
| `openmc.stats.fusion_neutron_spectrum(T_ion, mean_energy)` | Gaussian energy distribution for DT/DD fusion neutron emission. |
| `openmc.stats.muir(e0, m_rat, kt)` | Muir Gaussian energy spectrum for fusion plasmas. |

### Watt parameters for common fissile nuclides

| Nuclide | a (MeV) | b (MeV⁻¹) |
|---|---|---|
| U-235 thermal | 0.988 | 2.249 |
| Pu-239 thermal | 0.966 | 2.842 |
| U-233 thermal | 0.977 | 2.546 |

---

## Angular distributions

Used for the `angle` parameter of `IndependentSource`.

| Class | Description |
|---|---|
| `openmc.stats.Isotropic()` | Uniform on unit sphere |
| `openmc.stats.Monodirectional(reference_uvw)` | Single direction. `reference_uvw`: unit vector tuple. |
| `openmc.stats.PolarAzimuthal(mu, phi, reference_uvw)` | Separate distributions for polar (μ=cos θ) and azimuthal (φ) angles. `mu` and `phi` are `Univariate` distributions. |

```python
# Beam source: monodirectional +z
angle = openmc.stats.Monodirectional(reference_uvw=(0.0, 0.0, 1.0))

# Forward-peaked source using tabular mu distribution
mu = openmc.stats.Tabular([-1, 0, 1], [0.1, 0.3, 0.6])
phi = openmc.stats.Uniform(0, 2*3.14159)
angle = openmc.stats.PolarAzimuthal(mu=mu, phi=phi)
```

---

## Spatial distributions

Used for the `space` parameter of `IndependentSource`.

### Cartesian / general

| Class | Description | Parameters |
|---|---|---|
| `openmc.stats.Point(xyz)` | Delta function at a point | `xyz`: (x, y, z) |
| `openmc.stats.Box(lower_left, upper_right, only_fissionable)` | Uniform in rectangular box | `lower_left`, `upper_right` (3-tuples) |
| `openmc.stats.CartesianIndependent(x, y, z)` | Independent distributions for x, y, z | `x`, `y`, `z`: `Univariate` each |

### Cylindrical / spherical

| Class | Description | Parameters |
|---|---|---|
| `openmc.stats.CylindricalIndependent(r, phi, z, origin)` | Cylindrical coordinates with independent r, φ, z | `r`, `phi`, `z`: `Univariate` |
| `openmc.stats.SphericalIndependent(r, cos_theta, phi, origin)` | Spherical coordinates with independent r, θ, φ | `r`, `cos_theta`, `phi`: `Univariate` |
| `openmc.stats.spherical_uniform(r_inner, r_outer, cos_theta_min, cos_theta_max, phi_min, phi_max, origin)` | Uniform in a spherical shell sector | — |
| `openmc.stats.cylindrical_uniform(r_inner, r_outer, phi_min, phi_max, z_min, z_max, origin)` | Uniform in a cylindrical shell sector | — |

### Mesh-based

| Class | Description |
|---|---|
| `openmc.stats.MeshSpatial(mesh, strengths, volume_normalized)` | Uniform within each mesh element, weighted by `strengths`. |
| `openmc.stats.PointCloud(filename)` | Sample from a file of (x, y, z, weight) points. |

```python
# Volumetric source on a mesh (e.g., from activation)
mesh = openmc.RegularMesh()
mesh.lower_left, mesh.upper_right = [-10,-10,-10], [10,10,10]
mesh.dimension = [10, 10, 10]
strengths = np.ones(mesh.dimension)   # equal strength per cell

space = openmc.stats.MeshSpatial(mesh, strengths, volume_normalized=True)
```

---

## Base classes (for type checking or subclassing)

| Class | Description |
|---|---|
| `openmc.stats.Univariate` | Base class for all 1D distributions |
| `openmc.stats.UnitSphere` | Base class for angular distributions on unit sphere |
| `openmc.stats.Spatial` | Base class for 3D spatial distributions |
