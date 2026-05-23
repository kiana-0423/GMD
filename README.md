# GMD — Good Molecular Dynamics  `v2.2`

Release history is tracked in [CHANGELOG.md](CHANGELOG.md). The current documented release is `v2.2` from `2026-05-23`.

GMD is a C++20 molecular dynamics engine built around a small set of composable runtime abstractions:

- `ForceProvider` — force and energy evaluation (LJ, bonded, Ewald/PME, ML)
- `CompositeForceProvider` — sums multiple providers; total energy is automatically accumulated
- `NeighborBuilder` — Verlet-list construction with skin-distance rebuild check
- `Integrator` — velocity-Verlet time stepping with optional thermostat / barostat hooks
- `Simulation` — orchestrator that wires everything together

---

## What's New in v2.2

| Feature | Files |
|---|---|
| **3D MPI domain decomposition** — full 3D process grids with face/edge/corner ghost exchange, 3D periodic wraparound, `MPI_Alltoallv` reverse force accumulation, `MPI_Allgatherv` atom redistribution, and per-dimension periodicity control | `include/gmd/parallel/domain_decomposition.hpp` `include/gmd/parallel/mpi_communicator.hpp` `src/parallel/` |
| **`--proc-grid Px Py Pz` CLI flag** and **`mpi_grid` run.in directive** — user-selectable 3D process grids validated against MPI world size | `app/gmd_main.cpp` `include/gmd/io/config_loader.hpp` |
| **Comprehensive MPI test suite** — 17+ CTest targets covering 1D periodic ghost exchange/force/migration, 3D face/edge/corner ghost exchange, 3D reverse force, 8-rank LJ/Ewald/PME consistency, and 2-rank smoke tests | `tests/mpi_periodic_1d.cpp` `tests/mpi_ghost_exchange_3d.cpp` `tests/mpi_domain_decomposition_3d.cpp` `CMakeLists.txt` |
| **3D ghost image flags in Verlet lists** — `S[dim]` computed for all three dimensions, fixing neighbor images in 2D/3D process grids | `src/system/verlet_neighbor_builder.cpp` |
| **Periodic boundary migration test data** — atoms placed to exercise wraparound ghost exchange and redistribution | `tests/smoke_mpi_periodic_lj.{xyz,run}` `tests/smoke_mpi_boundary_migration.{xyz,run}` |

---

## What's New in v2.1

| Feature | Files |
|---|---|
| **MPI spatial domain decomposition** — balanced or user-selected 3D process grids, ghost-atom coordinate exchange, reverse force accumulation, migration, and allreduce helpers | `include/gmd/parallel/domain_decomposition.hpp` `include/gmd/parallel/mpi_communicator.hpp` `include/gmd/parallel/mpi_environment.hpp` `src/parallel/` |
| **Replicated PME under MPI** — PME charge meshes are summed across ranks, then each rank runs the same full 3D FFT; pencil metadata exists but distributed mesh/FFT communication is not implemented | `include/gmd/force/pme_force_provider.hpp` `include/gmd/parallel/pme_parallel.hpp` `src/force/pme_force_provider.cpp` |
| **`GMD_ENABLE_MPI` CMake option** — opt-in MPI build with `find_package(MPI REQUIRED)`, `MPI::MPI_CXX` linkage, and `GMD_ENABLE_MPI` compile definition | `CMakeLists.txt` `cmake/MPIOptions.cmake` |
| **MPI smoke tests** — multi-rank LJ/Ewald/PME runs and serial-vs-parallel logged energy consistency checks | `tests/smoke_mpi_lj.xyz` `tests/compare_energy_logs.cpp` `cmake/RunMpiLJConsistency.cmake` `cmake/RunMpiConsistency.cmake` |
| **CLI `--np N` flag** — validates expected process count against MPI world size | `app/gmd_main.cpp` |
| **Rank-aware I/O** — only rank 0 writes trajectory/log files; global coordinate gather for output frames | `app/gmd_main.cpp` |

---

## What's New in v2.0

| Feature | Files |
|---|---|
| **ML force provider integration** via TorchScript (LibTorch) — loads `.pt` models, reads `local_cutoff`, calls `forward(species, positions, edge_index, edge_shift)` | `include/gmd/force/torchscript_adapter.hpp` `include/gmd/force/ml_force_provider.hpp` `src/force/` |
| **`edge_index` and `edge_shift` tensors** — `NeighborList::image_flags` stores per-pair integer shift vectors; `VerletNeighborBuilder` computes them | `include/gmd/system/system.hpp` `src/system/verlet_neighbor_builder.cpp` |
| **Atomic numbers (`Z`)** — populated per atom during XYZ loading via built-in element table; forwarded to ML models | `include/gmd/system/system.hpp` `src/io/config_loader.cpp` |
| **`GMD_ENABLE_TORCH` CMake option** — opt-in LibTorch linkage; `force_field ml` and `model_path` directives in run files | `CMakeLists.txt` `app/gmd_main.cpp` |
| **PME smoke test** — validates Particle-Mesh Ewald path | `tests/smoke_pme.run` |

---

## What's New in v1.1

| Feature | Files |
|---|---|
| **Selectable LJ mixing rules** — `mixing_rule` / `mixing` now works for inline LJ in `run.in`, standalone LJ `.ff`, and molecular `.ff` files | `include/gmd/io/config_loader.hpp` `src/io/config_loader.cpp` |
| **Three supported combining rules** — `lorentz_berthelot` (default), `geometric`, and `waldman_hagler` for cross-type LJ pairs without explicit overrides | same files above |
| **Rule aliases and flexible spelling** — accepts `lb`, `geom`, `wh`, plus hyphenated names such as `lorentz-berthelot` and `waldman-hagler` | same files above |
| **Clearer non-bonded documentation** — examples and input format docs now describe `mixing_rule` and note that explicit `pair` entries override the configured fallback rule | `README.md` `examples/ethane_demo/ethane.ff` |

---

## What's New in v1.0

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
- cell-list Verlet neighbor lists with skin-distance rebuild check and per-pair image shift vectors
- shifted Lennard-Jones pair interactions (multi-element, selectable mixing rules, explicit pair overrides that take priority over mixing)
- **harmonic bonds** (`V = k(r − r₀)²`)
- **harmonic angles** (`V = k(θ − θ₀)²`)
- **periodic proper dihedrals** (`V = k[1 + cos(nφ − δ)]`)
- **harmonic improper dihedrals** (`V = k(φ − φ₀)²`)
- Velocity Verlet integration
- velocity-rescaling thermostat
- Nosé-Hoover thermostat (VVNH splitting)
- Berendsen barostat (weak-coupling, requires virial)
- **Monte Carlo barostat** (isotropic NPT, Metropolis criterion, no virial required, adaptive step-size)
- long-range Coulomb via Ewald and PME (self-contained 3D FFT, B-spline orders 4/6)
- **ML force provider** — TorchScript backend (`GMD_ENABLE_TORCH=ON`), SE3-GNN compatible, reads `local_cutoff` from model
- **MPI spatial domain decomposition** — balanced or user-selected 3D process grids, face/edge/corner ghost-atom exchange, reverse force accumulation via `MPI_Alltoallv`, atom redistribution via `MPI_Allgatherv`, and allreduce collectives (`GMD_ENABLE_MPI=ON`)
- replicated-mesh PME under MPI with pencil-decomposition metadata reserved for a future distributed FFT
- extended XYZ trajectory output and energy logging (rank-0 I/O with global coordinate gather in MPI mode)
- inline `run.in` LJ force-field definitions
- external `.ff` files (LJ and molecular; molecular runs default to bonded-only)
- external `.top` topology files
- atomic numbers (`Z`) populated from element symbols during XYZ loading

Not yet implemented:

- 1-2 / 1-3 nonbonded exclusion lists for molecular systems
- SHAKE / RATTLE bond constraints
- CUDA compute backend (CMake option present but CPU-only)
- checkpoint / restart
- Python bindings
- full unit / regression test coverage
- distributed PME mesh and FFT communication (current MPI PME path replicates the full mesh per rank)

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
| `GMD_ENABLE_TORCH` | `OFF` | Enable TorchScript ML force provider (requires LibTorch) |
| `GMD_ENABLE_MPI` | `OFF` | Enable MPI spatial domain decomposition (requires MPI) |

MPI build example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGMD_ENABLE_MPI=ON
cmake --build build --parallel
```

Torch + MPI build example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DGMD_ENABLE_TORCH=ON -DCMAKE_PREFIX_PATH=/path/to/libtorch \
      -DGMD_ENABLE_MPI=ON
cmake --build build --parallel
```

---

## Quick Start — Lennard-Jones Smoke Case

```bash
./build/gmd tests/smoke_lj.xyz tests/smoke_lj.run
```

Output files:

- `output.xyz` — extended XYZ trajectory
- `output.log` — tabular energy / temperature log

---

## MPI Parallel Run

Build with MPI enabled:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGMD_ENABLE_MPI=ON
cmake --build build --parallel
```

Run with multiple processes:

```bash
mpirun -np 4 ./build/gmd tests/smoke_mpi_lj.xyz tests/smoke_mpi_lj.run
```

The `--np N` CLI flag validates that the MPI world size matches expectations:

```bash
mpirun -np 4 ./build/gmd --np 4 input.xyz run.in
```

In MPI mode:
- Atoms are distributed across ranks via a balanced 3D process grid (computed by `MPI_Dims_create` or a balanced factorisation); `--proc-grid Px Py Pz` or `mpi_grid Px Py Pz` can request a compatible grid
- Atoms that cross domain boundaries are redistributed to their owning rank after each drift step
- Ghost atoms are exchanged across face, edge, and corner neighbors each step based on the neighbor-list cutoff + skin
- Forces on ghost atoms are reverse-accumulated back to their home ranks via `MPI_Alltoallv`
- Only rank 0 writes trajectory and energy log files
- Global coordinates are gathered for output frames
- Ewald and PME run with the 3D domain decomposition. PME currently sums rank-local charge assignment onto a replicated full mesh and performs a full FFT on every rank; it is not distributed PME.
- TorchScript ML force fields and the Monte Carlo barostat are rejected in multi-rank runs because their distributed ownership/trial-move contracts are not implemented.

---

## Ethane Demo — Molecular Force Field

A fully self-contained example in `examples/ethane_demo/` demonstrating the bonded force field:

```bash
cd /path/to/GMD
bash run_ethane_demo.sh
```

The script compiles all sources with `c++` (no CMake required), copies the input files into `build_demo/`, runs the simulation, and prints a per-step energy table.

The main CLI also supports molecular runs directly. By default, external
molecular FF runs use **bonded-only** mode because 1-2/1-3 non-bonded
exclusions are not implemented yet:

```bash
./build/gmd \
  examples/ethane_demo/ethane.xyz \
  examples/ethane_demo/ethane.run \
  examples/ethane_demo/ethane.ff \
  examples/ethane_demo/ethane.top
```

To explicitly enable molecular LJ anyway, add this to `run.in`:

```ini
molecular_nonbonded lj_unsafe
```

That opt-in path is intended only for diagnostics and smoke testing until
exclusion lists are implemented.

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

=== Building force providers ===
  BondedForceProvider ready
  LJ non-bonded: disabled (no 1-2/1-3 exclusions implemented)

=== Initial state (bonded only) ===
  Bond PE     = 0.011972 eV  (0.276080 kcal/mol)

=== Running MD  (500 steps) ===
Step    Time(fs)    PE(kcal/mol)    KE(kcal/mol)    E_tot(kcal/mol)  T(K)
------------------------------------------------------------------------------
0       0.000       0.276           6.260           6.536           300.000
10      10.000      1.603           5.057           6.660           242.346
...
500     500.000     4.119           3.478           7.597           166.691

=== Done ===
  Frames written : 11
  Final PE       : 0.178612 eV  (4.118959 kcal/mol)
```

> **Note:** The ethane demo runs bonded interactions only (bonds + angles +
> dihedrals). LJ is intentionally disabled because GMD does not yet implement
> 1-2/1-3 non-bonded exclusion lists; including LJ would put bonded C-H atoms
> (at 1.09 Å) deep inside the repulsive wall (σ_CH = 3.0 Å), which is
> unphysical. When exclusions are added this restriction will be lifted.

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
| `molecular_nonbonded none\|lj_unsafe` | external molecular FF mode; default `none` |

Inline LJ force field:

```ini
force_field lj
cutoff 8.5
mixing_rule lorentz_berthelot   # or geometric / waldman_hagler
type 1 Ar epsilon 0.01032 sigma 3.405 charge 0.0
type 2 Ne epsilon 0.00312 sigma 2.749
```

Supported mixing-rule spellings:

- `lorentz_berthelot`, `lorentz-berthelot`, `lb`
- `geometric`, `geom`
- `waldman_hagler`, `waldman-hagler`, `wh`

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

or `barostat monte_carlo`:

```ini
barostat monte_carlo
pressure 1.0
mc_frequency 25
mc_volume_step 0.01
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

ML force field (requires `GMD_ENABLE_TORCH=ON`):

```ini
force_field ml
model_path /path/to/model.pt
```

### `molecule.ff` — molecular force field

Used with `ConfigLoader::load_molecular_ff()` for bonded systems.  
All `kcal/mol` and degree values are **converted automatically** to internal `eV / rad` units.

```ini
force_field  molecular

mixing_rule  lorentz_berthelot   # optional; default is lorentz_berthelot
lj_cutoff  10.0

# Atom types: mass [amu], epsilon [kcal/mol], sigma [Å], charge [e]
type 1  C  mass 12.011  epsilon 0.0660  sigma 3.500  charge -0.180
type 2  H  mass  1.008  epsilon 0.0300  sigma 2.500  charge  0.060

# Explicit cross-pair override (overrides mixing_rule fallback)
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

// CLI default for molecular FFs: bonded-only (LJ excluded — no 1-2/1-3 exclusions yet).
simulation.set_force_provider(bonded);
// Total PE = bond + angle + dihedral + improper

// Explicit unsafe opt-in for diagnostics before exclusions are available:
// auto lj = std::make_shared<gmd::ClassicalForceProvider>(mff.lj);
// auto composite = std::make_shared<gmd::CompositeForceProvider>();
// composite->add(lj); composite->add(bonded);
// simulation.set_force_provider(composite);
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
# Non-MPI tests
cd build
ctest --output-on-failure

# MPI tests (requires GMD_ENABLE_MPI=ON)
ctest --output-on-failure -L mpi

# All tests
ctest --output-on-failure
```

When `GMD_ENABLE_MPI=ON`, CTest additionally registers MPI smoke tests for LJ,
Ewald, PME, and the minimal bonded molecular path, plus serial-vs-MPI LJ,
Ewald, and replicated-PME consistency checks. The consistency checks validate
the full logged trajectory row-by-row as well as the final energy drift.


| Test name | Exercises |
|---|---|
| `gmd_smoke_inline_lj` | inline LJ force field, NVT |
| `gmd_smoke_ewald` | Ewald electrostatics, NVT |
| `gmd_smoke_pme` | PME electrostatics, NVT |
| `gmd_smoke_mc_barostat` | MC barostat NPT, Nosé-Hoover |
| `gmd_smoke_molecular` | molecular FF (`.ff` + `.top`), bonded forces, explicit unsafe LJ opt-in / pair override path |
| `gmd_mpi_periodic_1d_ghost_exchange_2proc` | 1D periodic wraparound ghost exchange (rank 0 ↔ rank 1) |
| `gmd_mpi_periodic_1d_force_consistency_2proc` | distributed vs serial LJ force/energy across a periodic x boundary |
| `gmd_mpi_periodic_1d_migration_2proc` | atom migration across periodic x high→low and low→high boundaries |
| `gmd_mpi_domain_decomposition_3d_grid_4proc` | 3D process grid creation and coordinate mapping (4 ranks) |
| `gmd_mpi_domain_decomposition_3d_migration_8proc` | 3D multi-axis atom migration with periodic wrapping (8 ranks) |
| `gmd_mpi_ghost_exchange_3d_face_8proc` | 3D face ghost exchange (8 ranks) |
| `gmd_mpi_ghost_exchange_3d_edge_8proc` | 3D edge ghost exchange (8 ranks) |
| `gmd_mpi_ghost_exchange_3d_corner_8proc` | 3D corner ghost exchange (8 ranks) |
| `gmd_mpi_ghost_exchange_3d_periodic_8proc` | 3D periodic corner wraparound ghost exchange (8 ranks) |
| `gmd_mpi_reverse_force_3d_face_8proc` | 3D face reverse force accumulation (8 ranks) |
| `gmd_mpi_reverse_force_3d_edge_8proc` | 3D edge reverse force accumulation (8 ranks) |
| `gmd_mpi_reverse_force_3d_corner_8proc` | 3D corner reverse force accumulation (8 ranks) |
| `gmd_mpi_lj_3d_periodic_consistency_8proc` | 8-rank periodic LJ energy/force vs serial consistency |
| `gmd_smoke_mpi_lj_2proc` | 2-process MPI LJ NVT run (verifies execution does not crash) |
| `gmd_smoke_mpi_ewald_2proc` | 2-process MPI Ewald smoke run on a true cross-rank charged system |
| `gmd_smoke_mpi_pme_2proc` | 2-process MPI PME smoke run on a true cross-rank charged system |
| `gmd_smoke_mpi_molecular_2proc` | 2-process MPI bonded molecular smoke run across a rank boundary |
| `gmd_smoke_mpi_lj_consistency` | serial vs 2-process MPI energy consistency check (verifies numerical equivalence) |
| `gmd_smoke_mpi_ewald_consistency_8proc` | serial vs 8-process Ewald consistency on a 3D process grid |
| `gmd_smoke_mpi_pme_consistency_8proc` | serial vs 8-process replicated-PME consistency on a 3D process grid, including empty local domains |

---

## Project Layout

```
GMD/
├── app/gmd_main.cpp               CLI entry point (LJ / Ewald / PME / ML / MPI workflow)
├── cmake/                         Build helper scripts (CompilerOptions, CUDAOptions, MPIOptions, …)
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
│   ├── cuda/                      CUDA backend (reserved)
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
│   ├── ml/                        ML force provider + TorchScript adapter
│   ├── neighbor/                  VerletNeighborBuilder
│   ├── parallel/                  MPI domain decomposition + communication
│   │   ├── domain_decomposition.hpp
│   │   ├── mpi_communicator.hpp
│   │   ├── mpi_environment.hpp
│   │   └── pme_parallel.hpp
│   ├── runtime/                   RuntimeContext
│   ├── system/                    System, Box, Topology
│   └── utils/                     Utilities (reserved)
├── src/                           Implementation (mirrors include/gmd/)
├── tests/                         Smoke / integration / MPI unit / consistency tests
│   ├── mpi_periodic_1d.cpp            1D periodic MPI unit tests
│   ├── mpi_ghost_exchange_3d.cpp      3D ghost exchange + reverse force unit tests
│   ├── mpi_domain_decomposition_3d.cpp 3D grid + migration unit tests
│   ├── compare_energy_logs.cpp        Serial-vs-MPI energy log comparator
│   └── smoke_*                        Smoke test inputs and run files
├── run_ethane_demo.sh             One-command build + run for ethane demo
└── CMakeLists.txt
```

---

## Architecture Overview

```
ConfigLoader ─────────────────────────────────────────────────┐
  load_xyz()          → System (coords, masses, types, Z)      │
  load_run()          → RunConfig                              │
  load_molecular_ff() → MolecularForceFieldConfig              ├──► Simulation
  load_topology()     → Topology                               │       │
                                                               │  DomainDecomposition (MPI)
ForceProvider (interface)                                      │  MpiCommunicator (MPI)
  ClassicalForceProvider   (LJ non-bonded)                     │  NeighborBuilder::rebuild()
  BondedForceProvider      (bonds/angles/dihedrals/impropers)  │  ForceProvider::compute()
  EwaldForceProvider       (k-space Coulomb)                   │  Integrator::step()
  PMEForceProvider         (mesh Ewald)                        │       │
  MLForceProvider          (TorchScript)                       │  TrajectoryWriter::write_frame()
  CompositeForceProvider   (sums all above) ───────────────────┘
```

---

## Known Limitations

- No 1-2 / 1-3 nonbonded exclusion lists — molecular CLI therefore defaults to bonded-only, and `molecular_nonbonded lj_unsafe` remains explicitly unsafe
- No SHAKE / RATTLE bond constraints
- No GPU execution (CUDA option present but CPU-only)
- No checkpoint / restart
- Berendsen barostat requires virial from every active force term; current `Ewald` and `PME` paths provide it via a coordinate-virial approximation; barostat pressure is computed with MPI-allreduced kinetic energy and virial
- PME mesh storage and FFT work are replicated across MPI ranks; pencil-decomposition metadata exists, but distributed PME mesh ownership and FFT transpose communication are not active
- ML force provider requires `GMD_ENABLE_TORCH=ON` and a compatible TorchScript model; MPI domain decomposition is rejected because local-plus-ghost model energy ownership and message-passing halo depth are not defined
- MPI NVE/NVT is supported for implemented classical short-range, Ewald, replicated-PME, and current bonded paths; Berendsen NPT is also supported, but the Monte Carlo barostat is serial-only because MPI trial-volume coordination is not implemented
- The MPI LJ NVE regression tracks serial logged energies and drift within test tolerances (`gmd_smoke_mpi_lj_consistency`)
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
