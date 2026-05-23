# Changelog

All notable user-facing changes in GMD are documented here.

## [v2.3] - 2026-05-23

### Changed

- **Directory restructuring** — merged 12 subdirectories into 6 for better cohesion:
  - `boundary/` → `system/` (PBC and minimum-image are system geometry properties)
  - `neighbor/` → `system/` (neighbor lists are part of the system data structure)
  - `runtime/` → `core/` (runtime context is a simulation-core utility)
  - `ml/` → `force/` (ML force providers are force-field implementations)
  - Removed empty `cuda/` and `utils/` placeholder directories.
- All `#include` paths updated accordingly (e.g. `gmd/boundary/` → `gmd/system/`).
- `CMakeLists.txt` and `run_ethane_demo.sh` source file lists updated.
- No code logic or architecture changed — only file locations and include paths.

## [v2.2] - 2026-05-23

### Added

- **3D MPI domain decomposition** — the 1D x-axis decomposition introduced in v2.1 is now a full 3D process grid.
  - `DomainDecomposition` supports arbitrary `{Px, Py, Pz}` process grids with face, edge, and corner neighbor communication.
  - `choose_processor_grid()` uses `MPI_Dims_create` when MPI is available, falling back to a balanced factorisation for non-MPI builds.
  - `create_1d_decomposition()` remains available for backward compatibility and simple x-only splitting.
  - `neighbor_rank(offset)` correctly wraps periodic neighbors and returns `no_rank` (-1) for non-periodic boundaries.
  - `owner_rank(pos)` wraps periodic coordinates before determining the owning process.
- **`--proc-grid Px Py Pz` CLI flag** — explicitly specify a custom 3D process grid at runtime.
- **`mpi_grid Px Py Pz` run.in directive** — configure the process grid from the run file.
- **Comprehensive MPI unit tests** (17+ CTest targets when `GMD_ENABLE_MPI=ON`):
  - `gmd_mpi_periodic_1d_*` — 1D periodic ghost exchange, force consistency, and atom migration (2 ranks).
  - `gmd_mpi_domain_decomposition_3d_*` — 3D grid creation, coordinate mapping, and multi-axis atom migration (4/8 ranks).
  - `gmd_mpi_ghost_exchange_3d_*` — 3D face/edge/corner ghost exchange, periodic corner wrapping, and reverse force accumulation (8 ranks).
  - `gmd_mpi_lj_3d_periodic_consistency` — 8-rank periodic LJ energy/force consistency vs serial (8 ranks).
  - `gmd_smoke_mpi_ewald_2proc` / `gmd_smoke_mpi_pme_2proc` — 2-rank Ewald/PME smoke tests on cross-rank charged systems.
  - `gmd_smoke_mpi_molecular_2proc` — 2-rank bonded molecular smoke test across a rank boundary.
  - `gmd_smoke_mpi_ewald_consistency_8proc` / `gmd_smoke_mpi_pme_consistency_8proc` — serial vs 8-rank Ewald/PME consistency on 3D process grids.
- **Periodic 1D MPI boundary migration test data** — `tests/smoke_mpi_periodic_lj.{xyz,run}` and `tests/smoke_mpi_boundary_migration.{xyz,run}` with atoms placed to exercise wraparound ghost exchange and redistribution.

### Changed

- **`MpiCommunicator` ghost exchange now uses 3D neighbor offsets** — `neighbor_offsets()` iterates up to 26 face/edge/corner neighbors, automatically skipping offsets along axes where `proc_grid[dim] == 1`. `periodic_shift()` wraps coordinates for periodic boundaries in all three dimensions.
- **`MpiCommunicator` reverse force accumulation uses `MPI_Alltoallv`** — forces on ghost atoms are packed per home rank and exchanged in a single all-to-all collective, replacing the previous pair-wise `MPI_Sendrecv` scheme.
- **`MpiCommunicator` atom redistribution uses `MPI_Allgatherv`** — all atom states are gathered globally, periodic coordinates are normalized, and `owner_rank()` determines the new owning process.
- **`VerletNeighborBuilder` ghost image flags extended to 3D** — `S[dim]` is now computed for all three dimensions (previously only x), correctly handling ghost atoms that are periodically shifted along y or z in 3D decompositions.
- **`DomainDecomposition` periodic state stored as `std::array<bool, 3>`** — replaces the previous single `periodic_x` flag, enabling per-dimension periodicity control.
- **`validate_rank_grid()` added** — verifies that the MPI rank matches the domain process grid coordinate at the start of every collective operation.
- **`MpiCommunicator` exposes `pack_send_buffer()` and `unpack_recv_buffer()`** as private helpers with clear 3D offset semantics.

### Fixed

- Ghost image shift vectors (`S[dim]`) in `VerletNeighborBuilder` now correctly computed for all three spatial dimensions, fixing incorrect neighbor images in 2D and 3D process grids.
- `MPI_Allreduce` collectives in force providers (`ClassicalForceProvider`, `EwaldForceProvider`, `PMEForceProvider`, `BondedForceProvider`) are now placed outside early-return paths so ranks with zero atoms still participate.
- `compute_twice_ke()` and thermostat initializers use globally allreduced atom counts for consistent `dof_` across all ranks.
- `VelocityInitializer` COM velocity removal and temperature rescaling use `MPI_Allreduce` for consistent scaling factors.
- Double-allreduce of kinetic energy removed from `write_global_frame` in `gmd_main.cpp`.

---

## [v2.1] - 2026-05-22

### Added

- **MPI spatial domain decomposition** for distributed parallel simulations (initial 1D x-axis release, upgraded to full 3D in v2.2).
  - `DomainDecomposition` (`include/gmd/parallel/domain_decomposition.hpp`, `src/parallel/domain_decomposition.cpp`) provides 1D x-axis box splitting with configurable ghost width derived from cutoff + skin.
  - `MpiCommunicator` (`include/gmd/parallel/mpi_communicator.hpp`, `src/parallel/mpi_communicator.cpp`) encapsulates allreduce (scalar/vector), broadcast, barrier, ghost-atom coordinate exchange, and reverse force accumulation.
  - `MpiEnvironment` (`include/gmd/parallel/mpi_environment.hpp`, `src/parallel/mpi_environment.cpp`) provides RAII-style `MPI_Init`/`MPI_Finalize` with non-MPI fallback.
- **PME pencil decomposition** infrastructure.
  - `PmeParallelDecomposition` (`include/gmd/parallel/pme_parallel.hpp`, `src/parallel/pme_parallel.cpp`) defines a 2D process grid (Py × Pz) with transpose collectives (`x→y`, `y→z`, `z→y`, `y→x`). Currently the PME mesh is replicated; distributed FFT will be activated in a future release.
- **`GMD_ENABLE_MPI` CMake option** — opt-in MPI build with `find_package(MPI REQUIRED)`, `MPI::MPI_CXX` linkage, and `GMD_ENABLE_MPI` compile definition. See `cmake/MPIOptions.cmake` for the `gmd_configure_mpi_target` helper.
- **MPI smoke tests** — two new CTest targets (only registered when `GMD_ENABLE_MPI=ON`):
  - `gmd_smoke_mpi_lj_2proc`: runs a 2-process LJ simulation.
  - `gmd_smoke_mpi_lj_consistency`: compares serial and 2-process energy logs to verify numerical consistency.
  - New test collateral: `tests/smoke_mpi_lj.xyz`, `tests/smoke_mpi_lj.run`, `tests/compare_energy_logs.cpp`, `cmake/RunMpiLJConsistency.cmake`.
- **CLI `--np N` flag** — validates that the MPI world size matches the expected process count.
- **Rank-aware I/O** — only rank 0 writes trajectory and energy log files; global coordinates are gathered via `allreduce_vector` for output frames.
- **System MPI extensions** — `System` gains `atom_tags_`, `atom_owners_`, `num_local_atoms_` with corresponding accessors (`atom_tag()`, `atom_owner()`, `num_local_atoms()`, `mutable_atom_tags()`, `mutable_atom_owners()`). `resize()` accepts an optional `local_count` parameter.

### Changed (MPI correctness fixes)

- **`compute_twice_ke()`** (`src/integrator/thermostat.cpp`) — added `MPI_Allreduce` so all ranks return the same global kinetic energy. This ensures that temperature-dependent operations (thermostat scaling, barostat pressure) are consistent across ranks.
- **Nose-Hoover thermostat `initialize()`** (`src/integrator/nose_hoover_thermostat.cpp`) — `dof_` is now computed from the globally allreduced atom count so that the friction coefficient xi evolves identically on every rank.
- **Velocity rescaling thermostat `initialize()`** (`src/integrator/velocity_rescaling_thermostat.cpp`) — same global `dof_` fix.
- **`VelocityInitializer`** (`src/system/initializer.cpp`):
  - `remove_center_of_mass_velocity()` — added `MPI_Allreduce` for total mass and COM momentum so that all ranks subtract the same COM velocity.
  - `rescale_temperature()` — added `MPI_Allreduce` for kinetic energy and atom count so that the scaling factor is consistent.
  - `initialize()` — in MPI mode no longer early-returns when `atom_count() == 0`, because doing so would skip collective MPI calls and cause `MPI_ERR_TRUNCATE`.
- **`ClassicalForceProvider::compute()`** (`src/force/classical_force_provider.cpp`) — virial `MPI_Allreduce` moved outside the early-return path so that ranks with zero atoms still participate in the collective call.
- **`BondedForceProvider::compute()`** (`src/force/bonded_force_provider.cpp`) — added virial `MPI_Allreduce` for global bonded virial, also placed unconditionally after the computation block.
- **`EwaldForceProvider`** (`src/force/ewald_force_provider.cpp`):
  - `compute_reciprocal()` — structure factor S(k) built from local atoms only, then allreduced to obtain the global value; forces computed on local atoms only.
  - `compute_self_correction()` — q²_sum and Q_net computed from local atoms only, then allreduced.
  - `compute_real_space()` — added pair deduplication via atom_tag comparison for ghost-atom boundary pairs (same pattern as `ClassicalForceProvider`); outer i-loop limited to local atoms.
  - `compute()` — refactored `has_charges` early-exit into a coordinated check using `MPI_LOR` allreduce, ensuring all ranks reach the same decision.
- **`write_global_frame`** (`app/gmd_main.cpp`) — removed redundant `allreduce_scalar` wrapper around `compute_twice_ke`, which already performs its own allreduce. Previously this double-allreduced the kinetic energy in MPI builds.

### Changed

- `Simulation` gains `set_mpi_communicator()` and `set_domain_decomposition()`; `step()` now performs ghost exchange before force evaluation and reverse force accumulation afterward, with allreduce on potential energy.
- `ClassicalForceProvider::compute()` respects `num_local_atoms()` — only local atoms are iterated as primary index `i`; ghost atoms `j` still participate in pair evaluation.
- `VerletNeighborBuilder::rebuild()` builds neighbor lists for local atoms only; ghost atoms are included as secondary partners.
- `BondedForceProvider::compute()` uses owner-rank resolution for bonded terms spanning process boundaries.
- `app/gmd_main.cpp` initializes `MpiEnvironment`, distributes atoms via `keep_rank_local_atoms()`, sets up `DomainDecomposition`, and uses rank-gated logging throughout.

### Compatibility Notes

- All existing serial workflows are unaffected. Omitting `-DGMD_ENABLE_MPI=ON` produces an identical binary to previous releases.
- MPI and non-MPI builds share the same input file format; no input changes are required to run in parallel.
- The `--np N` flag is optional and only validated when MPI is active.

---

## [v2.0] - 2026-04-20

### Added

- **ML force provider integration** via TorchScript (LibTorch) backend.
  - New `TorchScriptModelRuntimeAdapter` (`include/gmd/force/torchscript_adapter.hpp`, `src/force/torchscript_adapter.cpp`) loads a `.pt` model exported by `gmd_se3gnn`, reads the `local_cutoff` attribute, and calls `forward(species, positions, edge_index, edge_shift)`.
  - `MLForceProvider` now exposes `float cutoff() const noexcept` (delegated to the adapter) so the correct `VerletNeighborBuilder` cutoff is set automatically.
  - Run files accept two new directives: `force_field ml` and `model_path /path/to/model.pt`.
- **`edge_index` and `edge_shift` tensors** for periodic-boundary ML models.
  - `NeighborList` stores per-pair integer image shift vectors (`image_flags [E,3]`) alongside the existing CSR neighbor array.
  - `VerletNeighborBuilder::rebuild()` computes and stores these shift vectors; `TorchScriptModelRuntimeAdapter` expands CSR half-pairs into a directed full graph and computes Cartesian `edge_shift = S * box.lengths` (Å).
- **Atomic numbers** (`Z`) populated per atom during XYZ loading.
  - `System` gains an `atomic_numbers_` array with `atomic_numbers()` / `mutable_atomic_numbers()` accessors.
  - `ConfigLoader::load_xyz` resolves element symbols to atomic numbers via a built-in `kElementAtomicNumbers` table and writes them into `System`.
  - `ModelEvaluationRequest` now carries an `atomic_numbers` span and a `neighbor_list` pointer, both forwarded by `MLForceProvider::compute()`.
- **`GMD_ENABLE_TORCH` CMake option**: when `ON`, `find_package(Torch REQUIRED)` is called, `src/force/torchscript_adapter.cpp` is compiled into `gmd_core`, and `${TORCH_LIBRARIES}` is linked. Build with `-DGMD_ENABLE_TORCH=ON -DCMAKE_PREFIX_PATH=/path/to/libtorch`.

### Changed

- `ModelRuntimeAdapter` interface extended with a default-impl `virtual float cutoff() const noexcept` and `ModelEvaluationRequest` fields for `atomic_numbers` and `neighbor_list`.
- `RunConfig` gains `force_field_type` and `ml_model_path` fields; the `force_field` parser in `ConfigLoader::load_run` now accepts `ml` in addition to `lj`.

### Compatibility Notes

- All existing `lj`, `molecular`, and inline force-field inputs are unaffected.
- The Torch backend is fully opt-in: omitting `-DGMD_ENABLE_TORCH=ON` produces an identical binary to previous releases, and a clear runtime error is emitted if a run file requests `force_field ml` on a non-Torch build.
- `image_flags` is a new field on `NeighborList`; code that constructs or clears `NeighborList` directly must call `clear()` (which now also clears `image_flags`).

## [v1.1] - 2026-04-20

### Added

- Configurable LJ cross-type mixing via `mixing_rule` / `mixing` for inline `run.in` force fields, standalone LJ `.ff` files, and molecular `.ff` files.
- Support for three built-in combining rules: `lorentz_berthelot` (default), `geometric`, and `waldman_hagler`.
- Flexible rule aliases and spellings, including `lb`, `geom`, `wh`, `lorentz-berthelot`, and `waldman-hagler`.

### Changed

- Explicit `pair` entries remain the highest-priority source for cross-type LJ parameters; the configured mixing rule is only used as a fallback.
- Non-bonded examples and documentation now describe `mixing_rule` directly, including the override behavior for explicit cross-pair definitions.

### Compatibility Notes

- Existing inputs remain compatible because the default mixing rule is still `lorentz_berthelot`.
- This release does not change the current molecular non-bonded safety model: external molecular force fields still default to bonded-only mode unless `molecular_nonbonded lj_unsafe` is requested explicitly.

## [v1.0] - 2026-04-19

### Added

- Bonded molecular mechanics support with harmonic bonds, harmonic angles, periodic dihedrals, and harmonic impropers.
- Topology parsing for `.top` files with bond, angle, dihedral, and improper sections.
- Molecular force-field loading for `.ff` files with per-type mass, epsilon, sigma, charge, bonded parameter tables, and explicit pair overrides.
- A complete ethane demo with `.xyz`, `.run`, `.ff`, and `.top` inputs plus a standalone build-and-run script.

### Notes

- `v1.0` established the first end-to-end molecular workflow in GMD on top of the existing LJ, neighbor-list, and integration infrastructure.
