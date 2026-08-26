# GMD — Good Molecular Dynamics  `v2.4`

Release history is tracked in [CHANGELOG.md](CHANGELOG.md). The current documented release is `v2.4` from `2026-05-24`.

GMD is a C++20 molecular dynamics engine built around a small set of composable runtime abstractions:

- `ForceProvider` — force and energy evaluation (LJ, bonded, Ewald/PME, ML)
- `CompositeForceProvider` — sums multiple providers; total energy is automatically accumulated
- `NeighborBuilder` — Verlet-list construction with skin-distance rebuild check
- `Integrator` — velocity-Verlet time stepping with optional thermostat / barostat hooks
- `Simulation` — orchestrator that wires everything together

---

## What's New in v2.4

| Feature | Files |
|---|---|
| **Validation-progress release status** — README and validation docs now distinguish analytic validation, regression baselines, tested workflows, and prototype interfaces | `README.md` `validation/README.md` |
| **Special-pair / 1-4 scaling validation** — topology-derived 1-2/1-3 exclusions and 1-4 LJ/Coulomb scaling are documented as analytic-validation covered | `include/gmd/system/special_pair_map.hpp` `validation/static_special_pairs/` |
| **SHAKE/RATTLE and checkpoint/restart status corrected** — both are documented as tested features, with MPI constraint and checkpoint scalability limits called out | `README.md` `tests/restart_continuity.py` `tests/constraint_solver_tests.cpp` |
| **PME external validation plan** — replicated PME remains tested but external-validation pending; LAMMPS PPPM/OpenMM PME reference work is planned separately | `validation/pme_external/README.md` |
| **Release evidence tooling** — collect serial/MPI configure, build, CTest, validation summaries, and environment info without touching existing build dirs | `scripts/collect_release_logs.sh` `docs/release_logs/README.md` |
| **Manual/nightly full MPI CI** — full MPI CTest workflow separate from ordinary PR CI | `.github/workflows/full-mpi.yml` |

### Correctness fixes (post-v2.4)

These change simulation results for the configurations they affect.

| Fix | Files |
|---|---|
| **Degrees of freedom now account for constraints and the COM setting** — one shared `compute_degrees_of_freedom()` replaces the hardcoded `3N - 3` in both thermostats and in trajectory temperature output. A constrained run previously reported a temperature that was too low, so the thermostat drove the system hotter than its target | `include/gmd/integrator/thermostat.hpp` `src/integrator/thermostat.cpp` `src/integrator/{nose_hoover,velocity_rescaling}_thermostat.cpp` `src/integrator/velocity_verlet_integrator.cpp` `app/gmd_main.cpp` |
| **Neighbor-list rebuild decided collectively under MPI** *(correctness safeguard, not a result-changing fix)* — local rebuild flags are combined with a logical OR so every rank rebuilds during the same force evaluation, or none does. Under the current domain-decomposition path ghost cleanup already clears the list every step, so this is expected to be a no-op there and may not alter existing decomposed trajectories; it prevents divergent control flow wherever a valid list survives between steps, and protects future implementations | `src/core/simulation.cpp` `src/parallel/mpi_communicator.cpp` |
| **Forces refreshed after a barostat volume change** — `VelocityVerletIntegrator::step()` no longer returns holding forces computed for the pre-scaling geometry; `Simulation::step()` skips the refresh when the cell did not change | `src/integrator/velocity_verlet_integrator.cpp` `src/core/simulation.cpp` |
| **Physically derived virial for Ewald, PME and bonded terms** — see *Virial and pressure* below | `src/force/{ewald,pme,bonded}_force_provider.cpp` `include/gmd/force/special_pair_coulomb.hpp` |
| **Box dimensions validated** — zero, negative and non-finite edge lengths are rejected where a box is installed, instead of reaching periodic wrapping or the neighbor builder | `include/gmd/system/box.hpp` `include/gmd/system/system.hpp` `src/system/verlet_neighbor_builder.cpp` |
| **Proper dihedral and improper forces corrected** — the terminal-atom force had a flipped sign and the middle-atom projection used the coefficients for the opposite `b1` convention. Any run with proper dihedrals or impropers was integrating incorrect torsional forces. Found by component-wise virial validation: a torsion angle is invariant under isotropic scaling, so the error is invisible in `tr(W)` | `src/force/bonded_force_provider.cpp` |
| **Nose-Hoover restart validates, never overwrites, the DOF** — a checkpoint's `dof` is checked against the run's authoritative count and an incompatible restart is rejected rather than silently continued with a mismatched thermostat mass | `src/integrator/nose_hoover_thermostat.cpp` `include/gmd/integrator/thermostat.hpp` |
| **Constraint lists normalised** — `(i,j)` and `(j,i)` collapse to one constraint, exact duplicates are dropped, conflicting target distances and self-constraints are rejected, so the DOF subtraction counts distinct constraints | `src/integrator/constraint_solver.cpp` `include/gmd/integrator/constraint_solver.hpp` |

---

## Virial and pressure

The virial tensor `W` feeds the barostats through `P = (2*KE + tr(W)) / 3V`, and
each force term now contributes it in the form that term actually requires:

| Term | Virial |
|---|---|
| LJ, Ewald/PME real space, special-pair Coulomb correction | pair virial `r_ij ⊗ F_ij` from the minimum-image separation |
| Bonds, angles, dihedrals, impropers | per-interaction `Σ_a (r_a - r_ref) ⊗ F_a`, with positions taken relative to one atom of the interaction |
| Ewald/PME reciprocal space | `Σ_k E_k [δ_ab - 2(1/4α² + 1/k²) k_a k_b]` |
| Ewald/PME net-charge correction | isotropic `U_net · δ_ab` (the term scales as `1/V`) |
| Ewald/PME self-energy | none — it does not depend on the cell |
| SHAKE/RATTLE constraints | `Σ_c (2Λ_c/dt)·(r_c ⊗ r_c)`, the endpoint force recovered from the converged RATTLE multipliers |

### Constraint (SHAKE/RATTLE) virial

**The splitting.** Constrained dynamics uses the standard velocity-Verlet
SHAKE/RATTLE splitting, in which the constraint force enters **twice** per step:

```
(1)  v_i += (dt/2m_i)·F_i(t)                    first half-kick
(2)  v_i += w_i·Σ_c s_ic·Γ_c·r_c(t)             SHAKE impulse      } solved
(3)  r_i += dt·v_i                              drift              } together
(4)  v_i += (dt/2m_i)·F_i(t+dt)                 second half-kick
(5)  v_i += w_i·Σ_c s_ic·Λ_c·r_c(t+dt)          RATTLE impulse
```

(2) is the constraint force at time level `t`, paired with `F(t)`; (5) is the
same force at `t+dt`, paired with `F(t+dt)`. They are different quantities.

**Which one the pressure needs, and the factor.** The providers are evaluated at
`r(t+dt)`, so the provider virial belongs to `t+dt` and its partner is the
RATTLE impulse. Matching (5) against (4) term by term,

```
(dt/2m_i)·G_i(t+dt) = w_i·Λ_c·r_c(t+dt)   ⇒   G_i(t+dt) = (2/dt)·Σ_c s_ic·Λ_c·r_c
W_constraint(t+dt) = Σ_c (2Λ_c/dt)·(r_c ⊗ r_c)
```

The factor is **`2/dt`**, and it is exact — this is an endpoint quantity, not a
step average, so **both halves of the reported virial live at the same time
level and there is no `O(dt)` mixing**. Pinned from the outside rather than
assumed: a freely rotating rigid dimer's endpoint constraint force must equal
the centripetal force `μω²d`, and it does, to `O(dt²)`. `1/dt` — the value that
would be right if SHAKE omitted its velocity impulse (2), leaving RATTLE to
carry the whole step — undershoots it by exactly two.

**Central and symmetric by construction.** RATTLE does not move atoms, so
`r_c(t+dt)` is fixed for the whole solve and summing the *scalar* multipliers
gives a pair force exactly along `r_c` and exactly antisymmetric — for coupled
constraints such as a rigid triangle as much as for isolated ones. Nothing is
projected and no residual is discarded.

**Sign.** `λ = -(r_c·v_c)/((w_i+w_j)|r_c|²)`, so a separating pair gives `λ < 0`
and an attractive restoring force with a negative trace — matching the negative
virial an attractive Lennard-Jones pair produces under the same convention.

**Why the kinetic term does not already cover it.** `2K` and the configurational
constraint virial are different quantities, and projecting the velocities does
not remove the need for the second. For a rigid rotor of reduced mass `μ`,

```
2K = μω²d²        tr W_constraint = -μω²d²        2K + tr W = 0
```

which is the physically required answer — a rigid body's frozen internal
coordinate contributes nothing to the pressure. Drop the constraint virial and
`2K` is left uncancelled, reporting a spurious `2K/3V`.

### Reported pressure vs current-geometry virial

These are two different things and the engine stores them separately.

- **`System::last_virial()`** is the *current geometry's* tensor: provider plus,
  when one belongs to that geometry, the constraint term. After a barostat
  rescale the forces are re-evaluated for a geometry no dynamics integrated, and
  that tensor has no constraint partner — it is then marked **invalid** rather
  than passed off as a complete pressure virial.
- **`System::step_thermodynamics()`** is the state *of the step that just
  finished* — pressure, virial, `2K`, volume and potential energy — captured
  together at the end of `finish_step()`, before any barostat runs. Its pressure
  is bit-for-bit the number the barostat itself consumed, and a later force
  evaluation deliberately does not touch the record.

The log row is that record: **`PE`, `KE`, `E_total`, `T`, `P` and `V` all come
from one state**, so `P = (2K + tr W)/3V` holds among exactly those numbers and
`E_total` is the energy of the configuration that produced them. Keeping the
completed step and the current geometry in one slot would force a rescale to
destroy the pressure the step actually measured, which is why constrained NPT
previously had no pressure to report.

The one thing that is *not* from that state is the **`.xyz` configuration**:
those are the current coordinates, which under a barostat are the rescaled ones
the next step starts from. Snapshotting coordinates every step to match would
cost an `N x 3` copy per step for a difference no thermodynamic quantity in the
row depends on. Under a barostat, therefore, the row's `V` lags the `.xyz`
coordinates by one rescale, by design.

**Validity is explicit, never a sentinel.** A numerical zero is an ordinary
pressure. When no complete pressure exists — the initial frame of a constrained
run, before any step has produced a RATTLE multiplier — the log writes `nan` and
sets the `P_valid` column to `0`.

**MPI output.** The MPI output path writes a separate `System` holding the
gathered global coordinates. It takes part in no dynamics, so every non-atomic
field the writer reads is brought across from the distributed `System` by
`System::copy_frame_state_from()` — box, potential energy, completed-step
record, SHAKE/RATTLE diagnostics and the virial state used as the initial-frame
fallback. Nothing there is reduced: each field is either replicated (box,
constraint virial) or already global (provider virial, and the completed-step
record built from globally reduced quantities), so every rank already holds the
identical value and reducing would multiply it. That identity is asserted
directly in `tests/mpi_constraint_virial.cpp` rather than assumed, and
`tests/constrained_pressure_reporting.py` compares whole logs from serial,
`np=2` and `np=4` runs of the real executable column by column.

**MPI constraint solving.** `ConstraintSolver` allgathers every owned atom and works from a
constraint list replicated against global atom tags, so every rank converges to
the same multipliers and the tensor each rank computes is *already global*. It
is therefore **not** reduced — an allreduce would multiply it by the rank count.
`tests/mpi_constraint_virial.cpp` runs the same fixture at 1, 2 and 4 ranks,
which catches both a dropped and a rank-multiplied contribution.

**Validated scope.** Against the centripetal force of a rigid rotor — a
reference independent of the implementation — for a dimer and for a rigid
triangle of three coupled constraints: the magnitude of the recovered endpoint
force, the sign, the trace, the identity `2K + tr W = 0`, and second-order
convergence of both errors in `dt`. Plus exact tensor symmetry, exact centrality
of every pair force, the two constraint half-impulses being measurably distinct,
energy conservation of a free rotor, per-constraint additivity, translation
invariance, invariance across a periodic boundary, rank invariance, and restart
continuity of the reported pressure. Constrained pressure is complete and
analytically validated for its **trace and diagonal** on orthorhombic cells. It
has **no independent external reference**: LAMMPS does not expose a separately
extractable constraint virial with a proven-equivalent definition, and the
quantity is dynamical, so it cannot be compared from a static configuration.

None of these depend on the coordinate origin or on how atoms happen to be
wrapped into the cell. In particular the reciprocal sum is **not** computed as
`Σ_i r_i ⊗ F_i`: the reciprocal energy depends on the cell explicitly, through
both the `1/V` prefactor and `k = 2πn/L`, so that expression is not the virial
at all there.

All of this is checked against the definition rather than against a closed form,
in two ways:

- **Isotropic.** `tr(W)` against `-dU/ds` under `r → s·r, L → s·L`.
- **Per-axis.** Each diagonal component `W_aa` against `-dU/dε` under a normal
  strain on that axis alone (`L_a → (1+ε)L_a`, fractional coordinates fixed),
  on deliberately non-cubic cells so no error hides behind cubic symmetry.

`tests/virial_finite_difference_tests.cpp` runs both for LJ, bonded (intact and
wrapped across a boundary), Ewald and PME, neutral and net-charged.
`tests/mpi_dof_virial.cpp` repeats both under 4-way domain decomposition, where a
virial reduced the wrong number of times fails by a factor of the rank count.

**Scope of the claim.** `Box` stores three edge lengths, so the engine
represents orthorhombic cells only and no shear strain can be applied.
**Only the three diagonal components are validated. The off-diagonal components
are not.** They are checked for the symmetry `W_ab == W_ba`, which is necessary
but not sufficient. Validating them would need a triclinic box representation
and a shear deformation, neither of which exists here.

The per-axis check earns its keep: it is what exposed the dihedral force bug
listed above, which the isotropic trace could not see because a torsion angle is
unchanged by isotropic scaling and so contributes zero to `tr(W)` either way.

---

## What's New in v2.3

| Feature | Files |
|---|---|
| **Directory restructuring** — merged 12 subdirectories into 6 for better cohesion; no code logic changed | `include/gmd/` `src/` `CMakeLists.txt` |
| **Include path cleanup** — old `gmd/boundary`, `gmd/neighbor`, `gmd/runtime`, and `gmd/ml` paths moved under `system`, `core`, and `force` | `include/gmd/` `src/` |

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
| **PME execution-mode interface** — replicated PME remains the numerical backend; `pme_mode distributed` is a distributed-PME interface/prototype path for workflow testing and future FFTW-MPI integration, not a completed distributed mesh/FFT implementation | `include/gmd/force/pme_force_provider.hpp` `src/force/pme_force_provider.cpp` `CMakeLists.txt` |
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
- shifted Lennard-Jones pair interactions (multi-element, selectable mixing rules, explicit pair overrides that take priority over mixing); the static LJ cluster has an analytic validation reference
- **harmonic bonds** (`V = k(r − r₀)²`)
- **harmonic angles** (`V = k(θ − θ₀)²`)
- **periodic proper dihedrals** (`V = k[1 + cos(nφ − δ)]`)
- **harmonic improper dihedrals** (`V = k(φ − φ₀)²`)
- topology-derived **special pairs**: 1-2/1-3 LJ and Coulomb exclusions plus configurable 1-4 scaling; the four-atom chain validation uses an analytic reference
- Velocity Verlet integration
- velocity-rescaling thermostat
- Nosé-Hoover thermostat (VVNH splitting)
- Berendsen barostat (weak-coupling, requires virial)
- **Monte Carlo barostat** (isotropic NPT, Metropolis criterion, no virial required, adaptive step-size)
- long-range Coulomb via Ewald and replicated PME (self-contained 3D FFT, B-spline orders 4/6); Ewald has an analytic static validation reference, while replicated PME is tested against a provisional regression baseline and still needs external validation
- **ML force provider** — TorchScript backend (`GMD_ENABLE_TORCH=ON`), SE3-GNN compatible, reads `local_cutoff` from model
- **MPI spatial domain decomposition** — balanced or user-selected 3D process grids, face/edge/corner ghost-atom exchange, reverse force accumulation via `MPI_Alltoallv`, atom redistribution via `MPI_Allgatherv`, and allreduce collectives (`GMD_ENABLE_MPI=ON`); covered by 2/4/8-rank tests
- replicated-mesh PME under MPI; the `pme_mode distributed` path currently exercises an interface/prototype while still using the replicated numerical backend
- SHAKE / RATTLE bond constraints; serial constraints are tested, and MPI constraints currently use a correctness-first global-gather projection rather than a scalable owner-based solver
- checkpoint / restart through the CLI; restart continuity is tested by comparing continuous 100-step runs against 50-step + checkpoint + restart + 50-step runs in serial and MPI
- extended XYZ trajectory output and energy logging (rank-0 I/O with global coordinate gather in MPI mode)
- inline `run.in` LJ force-field definitions
- external `.ff` files (LJ and molecular; molecular runs apply topology special pairs by default)
- external `.top` topology files
- atomic numbers (`Z`) populated from element symbols during XYZ loading

Not yet implemented:

- CUDA compute backend (CMake option present but CPU-only)
- Python bindings
- full unit / regression test coverage
- independent external validation for replicated PME (LAMMPS PPPM or OpenMM PME reference pending)
- true distributed PME/FFT: distributed charge assignment, grid decomposition, MPI FFT communication, and force interpolation from a distributed reciprocal grid

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
- `output.log` — tabular thermodynamic log with `PE`, `KE`, `E_total`, `T`, `P`, `V`

The write cadence defaults to every 100 steps and can be overridden with:

```ini
output_interval 10
```

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

The main CLI also supports molecular runs directly. External molecular FF
runs enable LJ with topology-derived special pairs by default: bonded 1-2
and angled 1-3 pairs are excluded, while dihedral endpoint 1-4 pairs use
configurable scaling.

```bash
./build/gmd \
  examples/ethane_demo/ethane.xyz \
  examples/ethane_demo/ethane.run \
  examples/ethane_demo/ethane.ff \
  examples/ethane_demo/ethane.top
```

To make that behavior explicit, or to customize OPLS-style 1-4 factors, add:

```ini
molecular_nonbonded special
lj_scale_14 0.5
coul_scale_14 0.833333333333
```

`molecular_nonbonded none` remains available for bonded-only diagnostics.

### SHAKE / RATTLE Constraints

Molecular runs can constrain selected bond lengths with SHAKE after the
Velocity-Verlet position update and RATTLE after the final velocity update:

```ini
constraints on
constraint_tolerance 1e-6
constraint_max_iterations 100
rattle on
constrain_bond_type 1
```

`constrain_bond_type N` uses the `r0` value from `bond_type N` as the target
distance. Multiple constrained bond types can be listed with repeated
`constrain_bond_type` lines or one `constrain_bond_types 1 2 ...` line.
Topology files may also include explicit constraints:

```text
constraints 2
  1 2 0.9572
  1 3 0.9572
```

The thermodynamic log appends SHAKE/RATTLE iteration counts and maximum
constraint errors for every written frame. In MPI runs, constraints are applied
with a correctness-first global gather over owned atoms keyed by stable atom
tags, so constraints may cross domain boundaries. This is intentionally simple;
it is not yet the scalable local constraint-communication algorithm needed for
very large constrained biomolecular systems.

### Checkpoint / Restart

Long simulations can write a restartable text checkpoint that is separate from
the normal `.xyz` trajectory and `.log` thermodynamic output:

```ini
write_checkpoint_every 1000
checkpoint_file state.gmdchk
```

The checkpoint starts with `GMD_CHECKPOINT 1` and stores the current step,
simulation time, box, atom tags/types/molecule ids/masses/charges, positions,
velocities, topology entries, force-field/topology file references, boundary
summary, run/config summary, thermostat state, barostat state, and RNG state
where the active module owns one. Neighbor lists, ghost atoms, PME grids,
forces, and domain decomposition are intentionally rebuilt after restart.

To continue from a checkpoint, use a run input containing:

```ini
restart_from state.gmdchk
run 50000
```

`run` means "additional steps" after the checkpoint step. Restart writes a
fresh `output.xyz` and `output.log` with continuous step/time values beginning
at the checkpoint step; it does not append to an existing trajectory. MPI
restart reads the same global checkpoint on every rank, then rebuilds domain
decomposition for the current MPI size. This supports changing the process
count for the currently replicated text checkpoint format, assuming the same
force-field/topology inputs are available.

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
  LJ non-bonded: enabled (1-2/1-3 excluded, 1-4 scaled)

=== Initial state ===
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

> **Note:** The standalone ethane demo driver may choose its provider set
> independently. The main CLI molecular workflow now enables non-bonded LJ
> safely through its topology-derived special-pair table.

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
| `molecular_nonbonded special\|none` | molecular LJ mode; default `special` (`lj_unsafe` accepted as a legacy alias) |
| `lj_scale_12`, `coul_scale_12` | 1-2 pair factors; default `0.0`, `0.0` |
| `lj_scale_13`, `coul_scale_13` | 1-3 pair factors; default `0.0`, `0.0` |
| `lj_scale_14`, `coul_scale_14` | 1-4 pair factors; default `0.5`, `0.833333333333` |

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
pme_mode replicated
pme_benchmark off
```

PME modes:

```ini
pme_mode replicated    # default, full PME mesh and FFT on every MPI rank
pme_mode distributed   # interface/prototype path; still uses replicated PME numerics
pme_mode auto          # reserved for future automatic backend choice
pme_benchmark on       # print charge/FFT/Green/interpolation timings
```

The distributed-PME prototype interface is guarded by:

```bash
cmake -S . -B build-mpi-pme \
  -DGMD_ENABLE_MPI=ON \
  -DGMD_ENABLE_DISTRIBUTED_PME=ON
```

`GMD_ENABLE_DISTRIBUTED_PME=ON` performs FFTW-MPI dependency checks
(`fftw3-mpi.h`, `libfftw3_mpi`, `libfftw3`) and fails at configure time if
they are missing. The default build keeps replicated PME and does not require
FFTW. In the current implementation, requesting `pme_mode distributed` selects
the distributed-PME interface/prototype path and timing hooks, but the
numerical work is still performed by the replicated PME backend and a warning
is printed. It is useful for exercising input/configuration, CMake dependency
checks, tests, and workflow compatibility; it does not reduce PME mesh memory
per rank and should not be used as evidence of scalable distributed-PME
performance.

Future true distributed PME still needs distributed charge assignment,
distributed grid ownership, an MPI FFT backend using slab or pencil
decomposition, reciprocal-space Green's-function application on distributed
mesh data, and force interpolation back from the distributed reciprocal grid
to domain-decomposed atoms.

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

system.set_special_pair_map(std::make_shared<gmd::SpecialPairMap>(*topo));
auto lj = std::make_shared<gmd::ClassicalForceProvider>(mff.lj);
auto composite = std::make_shared<gmd::CompositeForceProvider>();
composite->add(lj); composite->add(bonded);
simulation.set_force_provider(composite);
// Total PE = scaled/excluded LJ + bond + angle + dihedral + improper
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
| `gmd_smoke_molecular` | molecular FF (`.ff` + `.top`), bonded forces and special-pair LJ path |
| `gmd_special_pair_rules` | four-atom chain: hand-checked LJ exclusion/scaling plus Ewald/PME Coulomb corrections |
| `gmd_constraint_solver` | SHAKE/RATTLE two-atom, water geometry, and constrained NVE stability checks |
| `gmd_checkpoint_restart` | serial 100-step continuity vs 50+checkpoint+restart+50, state round-trip, bad-file errors |
| `gmd_pme_modes` | fixed charged system: replicated PME vs distributed-PME interface/prototype selection path energy/force equivalence, both using replicated numerics |
| `gmd_validation_static_lj_cluster` | validation suite: fixed-coordinate LJ cluster vs stored baseline |
| `gmd_validation_static_special_pairs` | validation suite: analytic 1-2/1-3 exclusion and 1-4 scaling reference, including a modified-scale variant |
| `gmd_validation_static_coulomb` | validation suite: charged Ewald analytic single-point reference plus provisional PME regression baseline |
| `gmd_mpi_periodic_1d_ghost_exchange_2proc` | 1D periodic wraparound ghost exchange (rank 0 ↔ rank 1) |
| `gmd_mpi_periodic_1d_force_consistency_2proc` | distributed vs serial LJ force/energy across a periodic x boundary |
| `gmd_mpi_periodic_1d_migration_2proc` | atom migration across periodic x high→low and low→high boundaries |
| `gmd_mpi_special_pair_consistency_2proc` | cross-rank four-atom chain: tag-based LJ/Ewald special-pair behavior |
| `gmd_mpi_constraint_solver_2proc` | cross-rank SHAKE/RATTLE constraint projection by global atom tag |
| `gmd_mpi_checkpoint_restart_2proc` | MPI ranks read the same global checkpoint and recover atom state consistently |
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

## Validation

GMD now ships a lightweight validation suite under `validation/`:

- `static_*` cases are short, deterministic, and registered in CTest.
- `nve_*`, `nvt_*`, `npt_*`, and `diffusion_*` are optional longer workflows.
- Static single-point checks use `gmd_validate`, which writes JSON containing component energies and per-atom forces.
- Dynamics cases post-process `output.log` and `output.xyz` into `summary.json`, `thermo_series.csv`, and, where relevant, `msd_series.csv`.

Examples:

```bash
# short validation cases
python3 validation/run_validation.py \
  --case short \
  --gmd build/gmd \
  --gmd-validate build/gmd_validate \
  --work-root build/validation

# optional long workflow validation
cmake --build build --target gmd_validation_long
```

Each case directory contains:

- GMD input files
- a reference result JSON
- where available, a LAMMPS/OpenMM reproduction note or input stub
- a Python analysis script
- an explicit tolerance file
- a case README documenting intent, settings, and limitations

Reference provenance is stated per case. Entries marked `analytic_reference_*`
are independent analytic checks. Entries marked `provisional_gmd_baseline` are
repeatable regression baselines only; they still need an external software
result to become scientific cross-code validation.

Current validation status:

- **LJ static cluster**: analytic reference, suitable for short CI.
- **Special-pair / 1-4 scaling static chain**: analytic shifted-LJ + Ewald special-pair reference, including a modified-scale variant.
- **Ewald static Coulomb**: analytic periodic Ewald reference.
- **Replicated PME static Coulomb**: tested as a regression baseline only; external LAMMPS PPPM or OpenMM PME validation is pending.
- **Distributed PME mode**: interface/prototype selection path only; it uses replicated PME numerics and is not evidence of scalable distributed FFT performance.
- **SHAKE/RATTLE**: tested for serial projection, water geometry, constrained NVE stability, and a cross-rank MPI projection; the MPI algorithm is correctness-first global gather.
- **Checkpoint/restart**: tested through the real `gmd` CLI restart-continuity workflow in serial and MPI.

---

## Continuous Integration

GitHub Actions runs a minimal CI workflow in `.github/workflows/ci.yml`:

- serial `RelWithDebInfo` configure/build with Ninja
- full serial CTest suite, including short static validation cases
- MPI `RelWithDebInfo` configure/build with OpenMPI
- MPI test registration check via `ctest -N`
- selected 2-rank MPI smoke/unit tests, including ghost exchange, migration,
  special pairs, constraints, checkpoint/restart, restart continuity, LJ,
  Ewald, PME, and molecular smoke tests

Long validation workflows (`nve_*`, `nvt_*`, `npt_*`, `diffusion_*`) are not
part of every CI run. They remain optional/manual validation because they are
slower and many currently use `provisional_gmd_baseline` regression references:

```bash
cmake --build build --target gmd_validation_long
```

A separate manual/nightly workflow can run the full MPI CTest suite. It is kept
out of ordinary pull-request CI because the 4/8-rank tests are slower and more
sensitive to runner MPI availability.

Local CI reproduction:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --parallel
ctest --test-dir build --output-on-failure --verbose

cmake -S . -B build-mpi -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DGMD_ENABLE_MPI=ON
cmake --build build-mpi --parallel
ctest --test-dir build-mpi -N
ctest --test-dir build-mpi --output-on-failure --verbose \
  -R "gmd_mpi_periodic_1d_ghost_exchange_2proc|gmd_mpi_periodic_1d_force_consistency_2proc|gmd_mpi_periodic_1d_migration_2proc|gmd_mpi_special_pair_consistency_2proc|gmd_mpi_constraint_solver_2proc|gmd_mpi_checkpoint_restart_2proc|gmd_mpi_restart_continuity_cli_2proc|gmd_smoke_mpi_lj_2proc|gmd_smoke_mpi_ewald_2proc|gmd_smoke_mpi_pme_2proc|gmd_smoke_mpi_molecular_2proc"
```

Release evidence can be collected without touching an existing `build/`
directory:

```bash
scripts/collect_release_logs.sh /tmp/gmd-v2.4-release-logs
```

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

- SHAKE / RATTLE constraints use a correctness-first MPI global-gather projection; scalable local constraint communication and SETTLE/LINCS are not implemented
- Constraint forces contribute to the virial as an endpoint quantity at `t+dt` (see *Virial and pressure*), so constrained pressure is complete and single-time-level. It has **no independent external reference**, only an analytical one
- Constraint **independence is assumed, not verified**. The solver normalises its list to distinct pairs and rejects duplicates, conflicts and self-constraints, but a redundant closed topology (for example all pairs among five or more atoms) is accepted and every distinct constraint is counted, which over-subtracts degrees of freedom. Deciding this in general needs the rank of the constraint Jacobian, which is configuration dependent. Configure independent constraints
- The **initial frame of a constrained run** has no completed step behind it and no RATTLE multiplier for its force evaluation, so it has no complete pressure. It is reported as `nan` with `P_valid 0`, never as a number
- The SHAKE reference-gradient linearisation has no solution when a bond turns through ~90° in a single step (`r(t+dt)·r(t) → 0`). That is diagnosed as a hard error naming the pair and asking for a smaller time step, rather than being allowed to produce a wild correction
- Off-diagonal virial components are unvalidated; see *Virial and pressure*
- A Nose-Hoover checkpoint can only be restarted into a run with the same degrees of freedom. Changing the constraint set, the centre-of-mass removal setting or the atom count is rejected, as is a checkpoint predating constraint-aware DOF accounting (the old `3N-3` rule) whenever the two disagree. There is no migration path: the thermostat mass `Q` and friction variable `xi` belong to the DOF they were generated under. Restart from the input instead
- No GPU execution (CUDA option present but CPU-only)
- Checkpoint/restart uses a readable replicated text file; large-scale binary/parallel checkpoint I/O is not implemented
- Berendsen barostat requires virial from every active force term; barostat pressure is computed with MPI-allreduced kinetic energy and virial. Every provider now reports a physically derived virial rather than a coordinate approximation (see *Virial and pressure* below), and each contribution is reduced exactly once across the communicator
- `pme_mode distributed` currently uses the replicated PME numerical backend. It is an interface/prototype for workflow compatibility, dependency checks, and timing hooks; it does not reduce PME grid memory per rank, does not perform distributed FFT communication, and should not be used as evidence of scalable distributed-PME performance
- Ewald/PME special-pair Coulomb scaling is applied with analytical `(scale - 1) q_i q_j/r` corrections; PME retains its normal mesh discretization error and still needs broader accuracy regression coverage
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
