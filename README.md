# GMD

GMD is a molecular dynamics software project built around a C++ core and CUDA backends, with optional Python bindings for workflow orchestration and machine learning integration.

## Project Goal

GMD aims to provide a stable and extensible molecular dynamics backend with the following priorities:

- C++ as the simulation core language
- CUDA as the high-performance compute backend
- Clear separation between simulation control, force evaluation, and model runtime
- Support for both classical force fields and machine learning force fields
- Python as an optional outer interface rather than a runtime dependency

## Architecture Principles

The project should follow these rules from the beginning:

1. The simulation core runs independently in C++.
2. CUDA kernels are scheduled and owned by the C++ runtime.
3. Python is only used for external orchestration, scripting, and ecosystem integration.
4. Neighbor construction, memory management, and integration stay inside the engine.
5. Force evaluation uses a unified backend interface so that classical and ML force fields share the same execution path.

The core execution path is:

```text
CLI / Optional Python API
    |
    v
C++ Simulation Engine
    |
    v
ForceProvider Interface
    |
    v
CUDA Kernels / ML Runtime Adapter
```

## Recommended Project Layout

The following structure is a practical starting point for a C++-first implementation:

```text
GMD/
├── CMakeLists.txt
├── README.md
├── cmake/
│   ├── Dependencies.cmake
│   ├── CompilerOptions.cmake
│   └── CUDAOptions.cmake
├── docs/
│   ├── architecture/
│   ├── algorithms/
│   └── design-notes/
├── include/
│   └── gmd/
│       ├── core/
│       ├── system/
│       ├── force/
│       ├── neighbor/
│       ├── integrator/
│       ├── runtime/
│       ├── io/
│       ├── ml/
│       ├── cuda/
│       └── utils/
├── src/
│   ├── core/
│   ├── system/
│   ├── force/
│   ├── neighbor/
│   ├── integrator/
│   ├── runtime/
│   ├── io/
│   ├── ml/
│   ├── cuda/
│   └── utils/
├── kernels/
│   ├── neighbor/
│   ├── pair/
│   ├── bonded/
│   ├── integrate/
│   ├── reduce/
│   └── features/
├── app/
│   └── gmd_main.cpp
├── configs/
│   ├── examples/
│   └── testcases/
├── examples/
│   ├── lj_fluid/
│   ├── water_box/
│   └── mlff_demo/
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── regression/
│   └── performance/
├── benchmarks/
│   ├── single_gpu/
│   └── ml_inference/
├── third_party/
└── python/
    └── gmd/
```

## Module Design

The engine should be split into the following code modules.

### 1. `core`

This is the top-level simulation control layer.

Main responsibilities:

- Simulation lifecycle management
- Time-step loop control
- Module scheduling
- Global configuration dispatch
- Runtime state transitions

Suggested classes:

- `Simulation`
- `SimulationContext`
- `StepController`
- `ExecutionPolicy`

Typical functions:

- initialize simulation state
- run a fixed number of steps
- perform one simulation step
- trigger module updates by stage
- coordinate checkpoint and restart

### 2. `system`

This module stores the physical system definition and simulation state.

Main responsibilities:

- atom and topology metadata
- particle state storage
- periodic box definition
- host-side and device-side views
- restartable system snapshots

Suggested classes:

- `System`
- `Topology`
- `ParticleData`
- `Box`
- `DeviceState`

Typical functions:

- allocate positions, velocities, and forces
- load topology and initial coordinates
- synchronize host/device views when needed
- apply periodic boundary metadata
- expose read-only views for force kernels

### 3. `neighbor`

This module controls short-range spatial search and neighbor maintenance.

Main responsibilities:

- cell list construction
- Verlet neighbor list construction
- rebuild criteria checks
- cutoff and skin handling
- particle displacement tracking

Suggested classes:

- `NeighborBuilder`
- `CellList`
- `NeighborList`
- `RebuildPolicy`

Typical functions:

- bin particles into spatial cells
- build pair candidate lists
- rebuild neighbors when displacement exceeds threshold
- provide neighbor views to force providers

### 4. `force`

This module defines the unified force backend abstraction.

Main responsibilities:

- force backend interface definition
- classical force field implementation
- dispatch of bonded and nonbonded terms
- energy, force, and virial accumulation

Suggested classes:

- `ForceProvider`
- `ClassicalForceProvider`
- `PairForce`
- `BondedForce`
- `ForceAccumulator`

Typical functions:

- evaluate forces from current state
- zero and accumulate force buffers
- compute total energy and optional virial
- combine multiple force terms into one result

### 5. `integrator`

This module advances the system in time.

Main responsibilities:

- time integration
- thermostat and barostat hooks
- constraint application hooks
- position and velocity updates

Suggested classes:

- `Integrator`
- `VelocityVerletIntegrator`
- `Thermostat`
- `Barostat`
- `ConstraintSolver`

Typical functions:

- half-step and full-step updates
- update positions and velocities
- apply temperature or pressure control
- wrap coordinates under periodic boundary conditions

### 6. `runtime`

This module owns backend execution resources.

Main responsibilities:

- CUDA stream and event management
- execution order and synchronization policy
- memory pools and temporary buffers
- backend capability checks

Suggested classes:

- `RuntimeContext`
- `CudaContext`
- `StreamPool`
- `BufferManager`

Typical functions:

- create and destroy execution resources
- provide stream-aware launch context
- stage temporary buffers for kernels
- isolate backend-specific errors

### 7. `io`

This module handles input, output, logging, and checkpoints.

Main responsibilities:

- configuration parsing
- structure and topology loading
- trajectory output
- checkpoint save and restore
- logging and diagnostics

Suggested classes:

- `ConfigLoader`
- `TrajectoryWriter`
- `CheckpointManager`
- `Logger`

Typical functions:

- parse run configuration files
- write coordinates, energies, and metadata
- save and load restart files
- emit structured runtime logs

### 8. `ml`

This module adapts machine learning force fields to the engine.

Main responsibilities:

- feature construction for ML force fields
- runtime adapter abstraction
- model metadata and capability checks
- optional uncertainty and active learning hooks

Suggested classes:

- `MLForceProvider`
- `FeatureBuilder`
- `ModelRuntimeAdapter`
- `TorchScriptProvider`
- `ONNXProvider`
- `TensorRTProvider`

Typical functions:

- build model inputs from positions and neighbors
- execute model inference on device
- return energy, forces, and optional stress
- validate model cutoff and supported atom types

### 9. `cuda`

This module contains CUDA-facing wrappers and backend-specific support code.

Main responsibilities:

- launch wrappers for kernels
- device utility code
- backend-specialized data layouts
- error translation from CUDA runtime

Suggested classes or components:

- `CudaLaunch`
- `DeviceBuffer`
- `KernelDispatch`

Typical functions:

- prepare launch parameters
- dispatch kernels on selected streams
- check CUDA runtime status
- provide reusable device helpers

### 10. `utils`

This module contains shared infrastructure that should not leak into domain logic.

Main responsibilities:

- error handling
- timers and profiling helpers
- type aliases and math helpers
- string and file utilities

Suggested classes:

- `Status`
- `Profiler`
- `Timer`
- `MathUtils`

## CUDA Kernel Groups

To keep compute code organized, device kernels should be grouped by behavior rather than by file size.

### `kernels/neighbor`

- build cell occupancy
- reorder particles by cell if needed
- generate neighbor list indices
- mark neighbor rebuild conditions

### `kernels/pair`

- Lennard-Jones interactions
- short-range Coulomb interactions
- pairwise ML local environment preprocessing

### `kernels/bonded`

- bond forces
- angle forces
- dihedral and improper forces

### `kernels/integrate`

- velocity half-step update
- position full-step update
- periodic wrapping
- optional constraint support kernels

### `kernels/reduce`

- total energy reduction
- temperature reduction
- virial and pressure reduction

### `kernels/features`

- graph edge construction
- radial basis feature generation
- angular feature support
- atom-wise feature packing for ML models

## Main Code Paths Required In This Project

The project should explicitly implement the following end-to-end code paths.

### A. Simulation startup path

- read configuration
- build system and topology
- allocate host and device buffers
- initialize runtime resources
- initialize selected force backend

### B. Single-step simulation path

- check whether neighbor rebuild is required
- rebuild neighbor structures if needed
- clear force buffers
- evaluate force provider
- integrate positions and velocities
- reduce observables if requested
- write output on configured intervals

### C. ML force field path

- construct local environment from neighbor data
- build model features on device
- run inference using selected ML runtime
- convert model outputs into engine force buffers
- feed results back to the integrator

### D. Checkpoint and restart path

- serialize current system state
- store runtime metadata and step counter
- restore device-visible state from checkpoint
- resume time integration safely

### E. Validation path

- compare forces with reference values
- monitor energy drift in NVE runs
- verify neighbor list correctness
- benchmark throughput and GPU occupancy

## Suggested Milestones

The implementation should proceed in this order:

1. Build the C++ core, runtime layer, and CLI.
2. Add CUDA neighbor list, pair force, and Velocity Verlet support.
3. Make the engine run a simple Lennard-Jones fluid case.
4. Add checkpoint, logging, and regression tests.
5. Add the `ForceProvider` abstraction for ML backends.
6. Integrate one stable ML runtime backend in C++.
7. Add optional Python bindings with `pybind11`.

## Recommended First Deliverable

The first stable version of GMD should support:

- single GPU execution
- periodic boundary conditions
- Verlet neighbor lists
- Lennard-Jones nonbonded interactions
- Velocity Verlet integration
- trajectory output
- checkpoint and restart
- one pluggable force backend interface

This keeps the first version narrow enough to implement, while preserving the correct architecture for later ML force field integration.
