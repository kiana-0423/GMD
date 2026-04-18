# GMD — Good Molecular Dynamics  `v0.4`

GMD is a C++20 molecular dynamics engine built around a small set of composable runtime abstractions:

- `ForceProvider` — force and energy evaluation (LJ, bonded, Ewald/PME, ML)
- `CompositeForceProvider` — sums multiple providers; total energy is automatically accumulated
- `NeighborBuilder` — Verlet-list construction with skin-distance rebuild check
- `Integrator` — velocity-Verlet time stepping with optional thermostat / barostat hooks
- `Simulation` — orchestrator that wires everything together

---

## What's New in v0.4

| Feature | Files |
|---|---|
| **Bonded force provider** — harmonic bonds, harmonic angles, periodic dihedrals, harmonic impropers | `include/gmd/force/bonded_params.hpp` `include/gmd/force/bonded_force_provider.hpp` `src/force/bonded_force_provider.cpp` |
| **Topology** — per-molecule connectivity (`BondTerm`, `AngleTerm`, `DihedralTerm`, `ImproperTerm`) | `include/gmd/system/topology.hpp` |
| **Molecular force field loader** (`load_molecular_ff`) — atom types with mass/ε/σ/q, explicit pair overrides, bonded parameter tables, inline comment support | `include/gmd/io/config_loader.hpp` `src/io/config_loader.cpp` |
| **Topology loader** (`load_topology`) — reads `.top` files with bond/angle/dihedral/improper sections | same files above |
| **Ethane demo** — complete working example with `.ff`, `.top`, `.xyz`, `.run` and a standalone build script | `examples/ethane_demo/` `run_ethane_demo.sh` |

---

## Current Scope

Implemented:

- periodic boundary conditions and minimum-image convention
- cell-list Verlet neighbor lists with skin-distance rebuild check
- shifted Lennard-Jones pair interactions (multi-element, Lorentz-Berthelot mixing, explicit pair overrides)
- **harmonic bonds** (`V = k(r − r₀)²`)
- **harmonic angles** (`V = k(θ − θ₀)²`)
- **periodic proper dihedrals** (`V = k[1 + cos(nφ − δ)]`)
- **harmonic improper dihedrals** (`V = k(φ − φ₀)²`)
- Velocity Verlet integration
- velocity-rescaling thermostat
- Nosé-Hoover thermostat
- Berendsen barostat
- long-range Coulomb via Ewald and PME
- extended XYZ trajectory output and energy logging
- inline `run.in` LJ force-field definitions
- external `.ff` files (LJ and molecular)
- external `.top` topology files

Not yet implemented:

- 1-2 / 1-3 nonbonded exclusion lists for molecular systems
- SHAKE / RATTLE bond constraints
- CUDA compute backend
- usable ML force-field runtime
- checkpoint / restart
- Python bindings
- full unit / regression test coverage

---

## Requirements

- C++ compiler with C++20 support (GCC 11+, Clang 14+, MSVC 2019+)
- CMake 3.24+

---

## Build

```bash
cd /path/to/GMD
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The main executable is produced at `build/gmd`.

Available CMake options:

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` / `Debug` |
| `GMD_ENABLE_CUDA` | `OFF` | Enable CUDA language (runtime still CPU-only) |
| `GMD_BUILD_PYTHON` | `OFF` | Python bindings stub |

---

## Quick Start — Lennard-Jones Smoke Case

```bash
./build/gmd tests/smoke_lj.xyz tests/smoke_lj.run
```

Output files:

- `output.xyz` — extended XYZ trajectory
- `output.log` — tabular energy / temperature log

---

## Ethane Demo — Molecular Force Field

A fully self-contained example in `examples/ethane_demo/` demonstrating the bonded force field:

```bash
cd /path/to/GMD
bash run_ethane_demo.sh
```

The script compiles all sources with `c++` (no CMake required), copies the input files into `build_demo/`, runs the simulation, and prints a per-step energy table.

The main CLI also supports molecular runs directly:

```bash
./build/gmd \
  examples/ethane_demo/ethane.xyz \
  examples/ethane_demo/ethane.run \
  examples/ethane_demo/ethane.ff \
  examples/ethane_demo/ethane.top
```

Expected output (truncated):

```
=== Loading inputs ===
  FF file       : ethane.ff
  Atom types    : 2
  Bond types    : 2
  Angle types   : 3
  Dihedral types: 1
  Bonds         : 7
  Angles        : 12
  Dihedrals     : 9

=== Initial state ===
  Bonded PE  =   0.276 kcal/mol
  LJ PE      = ...
  Total PE   = ...

=== Running MD  (500 steps) ===
Step    Time(fs)    PE(kcal/mol)    KE(kcal/mol)    E_tot(kcal/mol)  T(K)
------------------------------------------------------------------------
0       0.000       ...
10      10.000      ...
...
500     500.000     4039.238        ...

=== Done ===
  Frames written : 11
  Final PE       : 175.155249 eV  (4039.238 kcal/mol)
```

---

## Input Files

GMD accepts up to four plain-text input files. Lines starting with `#` and inline comments (everything after `#` on a line) are ignored.

### `xyz.in` — initial structure

```text
<N>
<Lx> <Ly> <Lz>
<atom-line> ...
```

Atom-line formats (`type` is 1-based; element symbol also accepted):

| Columns | Content |
|---|---|
| `type x y z` | mass from FF or built-in element table |
| `type x y z mass` | explicit mass |
| `type x y z vx vy vz` | velocities from file |
| `type x y z vx vy vz mass` | both |
| replace `type` with `Ar`, `C`, `H`, … | same variants |

### `run.in` — simulation parameters

Core directives:

| Directive | Meaning |
|---|---|
| `run <N>` | number of MD steps |
| `time_step <dt>` | time step [fs] |
| `velocity <T>` | target / initial temperature [K] |
| `velocity_init random\|input` | randomize or read from xyz |
| `velocity_seed <s>` | RNG seed |
| `remove_com_velocity true\|false` | remove center-of-mass drift |

Inline LJ force field:

```ini
force_field lj
cutoff 8.5
type 1 Ar epsilon 0.01032 sigma 3.405 charge 0.0
type 2 Ne epsilon 0.00312 sigma 2.749
```

Thermostat:

```ini
thermostat nose_hoover
thermostat_tau 100.0    # coupling time [fs]
```

or `thermostat velocity_rescaling`.

Barostat:

```ini
barostat berendsen
pressure 1.0
barostat_tau 2000.0
compressibility 4.5e-5
```

Long-range Coulomb (Ewald or PME):

```ini
coulomb ewald
ewald_alpha 0.3
ewald_kmax 7
ewald_cutoff 10.0
```

```ini
coulomb pme
pme_alpha 0.3
pme_cutoff 10.0
pme_order 4
pme_grid 32 32 32
```

### `molecule.ff` — molecular force field

Used with `ConfigLoader::load_molecular_ff()` for bonded systems.  
All `kcal/mol` and degree values are **converted automatically** to internal `eV / rad` units.

```ini
force_field  molecular

lj_cutoff  10.0

# Atom types: mass [amu], epsilon [kcal/mol], sigma [Å], charge [e]
type 1  C  mass 12.011  epsilon 0.0660  sigma 3.500  charge -0.180
type 2  H  mass  1.008  epsilon 0.0300  sigma 2.500  charge  0.060

# Explicit cross-pair override (overrides Lorentz-Berthelot mixing)
pair 1 2  epsilon 0.0447  sigma 3.000

# Bond types: k [kcal/mol/Å²], r0 [Å]
bond_type 1  k 310.0  r0 1.526   # C-C
bond_type 2  k 340.0  r0 1.090   # C-H

# Angle types: k [kcal/mol/rad²], theta0 [deg]
angle_type 1  k 40.0  theta0 112.7   # C-C-C
angle_type 2  k 50.0  theta0 110.7   # H-C-H
angle_type 3  k 50.0  theta0 110.7   # H-C-C

# Dihedral types: k [kcal/mol], n (periodicity), delta [deg]
dihedral_type 1  k 0.150  n 3  delta 0.0

# Improper types: k [kcal/mol/rad²], phi0 [deg]
improper_type 1  k 10.5  phi0 0.0
```

### `molecule.top` — topology

Used with `ConfigLoader::load_topology()`. All indices are 1-based.

```text
bonds  7
1 2  bond_type 2   # C1-H1a
1 5  bond_type 1   # C1-C2

angles  3
2 1 5  angle_type 3   # j = vertex atom

dihedrals  9
2 1 5 6  dihedral_type 1

impropers  0
```

---

## Combining Force Providers

`CompositeForceProvider` accumulates energy and forces from an ordered list of providers:

```cpp
gmd::ConfigLoader loader;

auto mff  = loader.load_molecular_ff("molecule.ff");
auto topo = loader.load_topology("molecule.top");

auto bonded = std::make_shared<gmd::BondedForceProvider>(topo);
for (auto& bp : mff.bond_types)     bonded->add_bond_type(bp);
for (auto& ap : mff.angle_types)    bonded->add_angle_type(ap);
for (auto& dp : mff.dihedral_types) bonded->add_dihedral_type(dp);
for (auto& ip : mff.improper_types) bonded->add_improper_type(ip);

auto lj = std::make_shared<gmd::ClassicalForceProvider>(mff.lj);

auto composite = std::make_shared<gmd::CompositeForceProvider>();
composite->add(lj);      // non-bonded LJ
composite->add(bonded);  // bonds + angles + dihedrals + impropers

simulation.set_force_provider(composite);
// Total PE = LJ + bond + angle + dihedral + improper
```

---

## Bonded Potential Functions

| Term | Potential | Parameters |
|---|---|---|
| Bond | $V = k(r - r_0)^2$ | `k` [eV/Å²], `r0` [Å] |
| Angle | $V = k(\theta - \theta_0)^2$ | `k` [eV/rad²], `theta0` [rad] |
| Dihedral | $V = k[1 + \cos(n\varphi - \delta)]$ | `k` [eV], `n` (int), `delta` [rad] |
| Improper | $V = k(\varphi - \varphi_0)^2$ | `k` [eV/rad²], `phi0` [rad] |

Forces are derived analytically. Dihedral and improper forces use the Blondel–Karplus gradient distribution to all four atoms.

Unit conversion helpers available in `bonded_force_provider.hpp`:

```cpp
gmd::kcal_per_mol_to_eV   // = 4.336410e-2
gmd::deg_to_rad           // = π / 180
```

---

## Testing

```bash
cd build
ctest --output-on-failure
```

Currently runs the bundled smoke/integration test based on
`tests/smoke_lj.xyz` and `tests/smoke_lj.run`. The rest of `tests/` is still
reserved for future unit and regression cases.

---

## Project Layout

```
GMD/
├── app/gmd_main.cpp               CLI entry point (LJ / Ewald / PME workflow)
├── cmake/                         Build helper scripts
├── examples/
│   └── ethane_demo/               Complete molecular FF example
│       ├── ethane_demo.cpp        Standalone demo driver
│       ├── ethane.ff              Molecular force field (OPLS-AA subset)
│       ├── ethane.top             Bond/angle/dihedral topology
│       ├── ethane.xyz             Initial coordinates (C2H6, 8 atoms)
│       └── ethane.run             MD parameters (500 steps, 300 K, NHC)
├── include/gmd/
│   ├── boundary/                  PBC & minimum-image
│   ├── core/                      Simulation orchestrator
│   ├── force/
│   │   ├── force_provider.hpp         abstract interface
│   │   ├── composite_force_provider.hpp
│   │   ├── classical_force_provider.hpp   (LJ)
│   │   ├── bonded_params.hpp             (BondParams, AngleParams, …)
│   │   ├── bonded_force_provider.hpp     (bonds / angles / dihedrals)
│   │   ├── ewald_force_provider.hpp
│   │   └── pme_force_provider.hpp
│   ├── integrator/                VV + thermostat + barostat
│   ├── io/                        ConfigLoader (xyz/run/ff/top), TrajectoryWriter
│   ├── ml/                        ML force provider (placeholder)
│   ├── neighbor/                  VerletNeighborBuilder
│   ├── runtime/                   RuntimeContext
│   └── system/                    System, Box, Topology
├── src/                           Implementation (mirrors include/gmd/)
├── tests/                         Unit / regression tests (reserved)
├── run_ethane_demo.sh             One-command build + run for ethane demo
└── CMakeLists.txt
```

---

## Architecture Overview

```
ConfigLoader ─────────────────────────────────────────────────┐
  load_xyz()          → System (coords, masses, types)         │
  load_run()          → RunConfig                              │
  load_molecular_ff() → MolecularForceFieldConfig              ├──► Simulation
  load_topology()     → Topology                               │       │
                                                               │  NeighborBuilder::rebuild()
ForceProvider (interface)                                      │  ForceProvider::compute()
  ClassicalForceProvider   (LJ non-bonded)                     │  Integrator::step()
  BondedForceProvider      (bonds/angles/dihedrals/impropers)  │       │
  EwaldForceProvider       (k-space Coulomb)                   │  TrajectoryWriter::write_frame()
  PMEForceProvider         (mesh Ewald)                        │
  CompositeForceProvider   (sums all above) ───────────────────┘
```

---

## Known Limitations

- No 1-2 / 1-3 nonbonded exclusion lists — adjacent bonded atoms still interact via LJ
- No SHAKE / RATTLE bond constraints
- No GPU execution (CUDA option present but CPU-only)
- No working ML inference backend
- No checkpoint / restart
- Berendsen barostat requires virial from every active force term; current Ewald / PME providers still disable it
- Limited automated test coverage

---

## Contributing

Contributions are welcome. If you change behavior:

- Keep public headers in `include/gmd/` and implementations in `src/`
- Add or expand tests in `tests/` where possible
- Update this README when input formats or runtime behavior change
- Keep `ForceProvider`, `Integrator`, and `NeighborBuilder` as the primary extension points

---

## License

This project is licensed under the terms in [LICENSE](LICENSE).

## Citation

If you use GMD in research, you can cite it as software:

```bibtex
@software{gmd2026,
  title={GMD: Good Molecular Dynamics},
  author={Contributors},
  year={2026},
  url={https://github.com/kiana-0423/GMD}
}
```
