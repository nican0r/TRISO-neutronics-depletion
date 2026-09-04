# Neutronics Model

The neutronics model (`src/triso/`) represents a single AGR-1 TRISO fuel compact treated as an infinite periodic lattice — reflective boundaries on all faces, no neutron leakage. The eigenvalue computed is therefore k-inf. The model is composed of three parts: materials, geometry, and the eigenvalue calculation with tallies.

---

## Materials

Each of the six material regions in the TRISO compact is an independent OpenMC material (`src/triso/materials.py`). All materials are set to **1200 K**, the upper HTGR fuel design temperature and the only choice that is self-consistent with the available thermal scattering data (the `c_Graphite` S(α,β) table provides a 1200 K point but not a 900 K one; running at 900 K would cause OpenMC to silently round graphite thermal scattering up to 1200 K).

| Region | Composition | Density | Source |
|---|---|---|---|
| UCO kernel | UC₀.₅O₀.₄, 19.75 wt% HALEU | 10.5 g/cm³ | INL/EXT-10-19476 Table 3; Demkowicz et al. 2018 |
| Buffer | Porous carbon | 1.0 g/cm³ | INL/EXT-10-19476 Table 3 (~50% dense) |
| IPyC / OPyC | Pyrolytic carbon | 1.87 g/cm³ | INL/EXT-10-19476 (midpoint of 1.85–1.90 g/cm³) |
| SiC | Silicon carbide | 3.20 g/cm³ | INL/EXT-10-19476 Table 3 |
| Graphite matrix | Graphite + c_Graphite S(α,β) | 1.75 g/cm³ | INL/EXT-10-19476 Table 5 |

The `c_Graphite` thermal scattering kernel is applied to the matrix only. PyC layers are turbostratic carbon, not crystalline graphite — applying `c_Graphite` to them would be physically incorrect.

The uranium enrichment is specified via `add_element('U', enrichment=19.75)`, which uses a fixed U-234/U-235 mass ratio of 0.008 (valid for centrifuge-enriched material). U-234 contributes ~0.18 wt% of total uranium; the effect on k-eff is negligible. Only the UCO kernel is marked `depletable=True` — the coating layers and matrix do not contain fissile material and are held at fixed composition throughout all burnup calculations.

---

## Geometry

The geometry (`src/triso/geometry.py`) has two levels: the individual TRISO particle and the packed compact.

**TRISO particle:** Each particle is a set of five concentric spheres sharing a common center, with outer radii taken directly from AGR-1 as-fabricated targets (INL/EXT-10-19476 Table 3):

| Layer | Outer radius |
|---|---|
| Kernel | 0.0175 cm |
| Buffer | 0.0275 cm |
| IPyC | 0.0315 cm |
| SiC | 0.0350 cm |
| OPyC | 0.0390 cm |

**Compact:** `openmc.model.pack_spheres` fills the AGR-1 compact cylinder (r = 0.62 cm, h = 2.5 cm) with non-overlapping sphere centers at **30% packing fraction** (the AGR-1 design target; typical range 25–35%). Each center becomes a TRISO object pointing to the shared particle universe. The particles are then embedded in a rectilinear lattice with a pitch of 0.078 cm (2 × R_OPyC) — the minimum pitch that prevents any particle from spanning two lattice cells, which would cause geometry errors in OpenMC's cell-finding algorithm. The graphite matrix fills all lattice cells not occupied by a particle.

Three reflective surfaces — a bounding cylinder and two end-cap planes — enclose the compact. This models an infinite periodic array of identical compacts and removes any leakage term from the eigenvalue, giving k-inf. A helium coolant channel is not included in this geometry; adding one requires a second bounding cylinder and graphite sleeve and is deferred.

---

## Eigenvalue calculation and tallies

The eigenvalue run (`src/triso/run.py`) uses **200 batches, 50 inactive, 5,000 particles per batch** (~750,000 active histories). This gives σ(k-eff) ≈ 10–20 pcm and tally relative errors ≲ 5% on the kernel and matrix — sufficient for physics validation on a laptop in ~5–15 min.

Four tallies are attached, all using `MaterialFilter` rather than `CellFilter`. The TRISO lattice contains hundreds of repeated cell instances; a `CellFilter` on a nominal cell ID would score only one particle's contribution. `MaterialFilter` correctly aggregates across all instances of each material type.

| Tally | Score | Purpose |
|---|---|---|
| `flux by material` | flux | Energy-integrated flux per layer type |
| `flux by material 3-group` | flux | Thermal / epithermal / fast flux per region (boundaries: 0.625 eV, 100 keV) |
| `kernel reaction rates` | fission, absorption | Reaction rates in the UCO kernel |
| `heating by material` | heating-local | Local energy deposition per layer (charged particles + recoil nuclei) |

The 3-group boundaries follow the IAEA standard: 0.625 eV thermal cutoff, 100 keV fast threshold. The `heating-local` score is valid without photon transport enabled; a full `heating` score (including photon contribution) requires coupled neutron-photon transport and is deferred.
