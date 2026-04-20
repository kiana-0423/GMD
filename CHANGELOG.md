# Changelog

All notable user-facing changes in GMD are documented here.

## [v2.0] - 2026-04-20

### Added

- **ML force provider integration** via TorchScript (LibTorch) backend.
  - New `TorchScriptModelRuntimeAdapter` (`include/gmd/ml/torchscript_adapter.hpp`, `src/ml/torchscript_adapter.cpp`) loads a `.pt` model exported by `gmd_se3gnn`, reads the `local_cutoff` attribute, and calls `forward(species, positions, edge_index, edge_shift)`.
  - `MLForceProvider` now exposes `float cutoff() const noexcept` (delegated to the adapter) so the correct `VerletNeighborBuilder` cutoff is set automatically.
  - Run files accept two new directives: `force_field ml` and `model_path /path/to/model.pt`.
- **`edge_index` and `edge_shift` tensors** for periodic-boundary ML models.
  - `NeighborList` stores per-pair integer image shift vectors (`image_flags [E,3]`) alongside the existing CSR neighbor array.
  - `VerletNeighborBuilder::rebuild()` computes and stores these shift vectors; `TorchScriptModelRuntimeAdapter` expands CSR half-pairs into a directed full graph and computes Cartesian `edge_shift = S * box.lengths` (Å).
- **Atomic numbers** (`Z`) populated per atom during XYZ loading.
  - `System` gains an `atomic_numbers_` array with `atomic_numbers()` / `mutable_atomic_numbers()` accessors.
  - `ConfigLoader::load_xyz` resolves element symbols to atomic numbers via a built-in `kElementAtomicNumbers` table and writes them into `System`.
  - `ModelEvaluationRequest` now carries an `atomic_numbers` span and a `neighbor_list` pointer, both forwarded by `MLForceProvider::compute()`.
- **`GMD_ENABLE_TORCH` CMake option**: when `ON`, `find_package(Torch REQUIRED)` is called, `src/ml/torchscript_adapter.cpp` is compiled into `gmd_core`, and `${TORCH_LIBRARIES}` is linked. Build with `-DGMD_ENABLE_TORCH=ON -DCMAKE_PREFIX_PATH=/path/to/libtorch`.

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
