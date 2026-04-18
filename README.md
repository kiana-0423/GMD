# GMD

GMD is a molecular dynamics software project built around a C++ core and CUDA backends, with optional Python bindings for workflow orchestration and machine learning integration.

## Table of Contents

- [Project Goal](#project-goal)
- [Requirements](#requirements)
- [Building](#building)
- [Quick Start](#quick-start)
- [Architecture Principles](#architecture-principles)
- [Project Structure](#project-structure)
  - [Folder Descriptions](#folder-descriptions)
- [Module Design](#module-design)
- [CUDA Kernel Groups](#cuda-kernel-groups)
- [Main Code Paths](#main-code-paths-required-in-this-project)
- [Suggested Milestones](#suggested-milestones)
- [Development Workflow](#development-workflow)
- [Contributing](#contributing)
- [Project Status](#project-status)
- [Development Log](#development-log)
- [License](#license)
- [Support and Contact](#support-and-contact)

## Project Goal

GMD aims to provide a stable and extensible molecular dynamics backend with the following priorities:

- C++ as the simulation core language
- CUDA as the high-performance compute backend
- Clear separation between simulation control, force evaluation, and model runtime
- Support for both classical force fields and machine learning force fields
- Python as an optional outer interface rather than a runtime dependency

## Requirements

- **C++ Compiler**: GCC 11+, Clang 14+, or MSVC 2019+ (C++20 standard required)
- **CMake**: Version 3.24 or higher
- **CUDA Toolkit** (optional): Version 11.0+ for GPU acceleration
  - NVIDIA GPU with compute capability 7.0+ recommended for optimal performance
- **Python** (optional): Python 3.8+ for Python bindings (requires `pybind11`)

## Building

### Basic Build (CPU only)

```bash
cd /path/to/GMD
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel $(nproc)
```

The compiled executable will be at `build/gmd`.

### Build with CUDA Support

```bash
cd /path/to/GMD
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DGMD_ENABLE_CUDA=ON
cmake --build . --parallel $(nproc)
```

### Build with Python Bindings

```bash
cd /path/to/GMD
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DGMD_BUILD_PYTHON=ON
cmake --build . --parallel $(nproc)
```

### Build with Both CUDA and Python

```bash
cd /path/to/GMD
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DGMD_ENABLE_CUDA=ON -DGMD_BUILD_PYTHON=ON
cmake --build . --parallel $(nproc)
```

### Build Options

- `-DCMAKE_BUILD_TYPE=Release|Debug`: Set optimization level (Release is recommended for simulations)
- `-DGMD_ENABLE_CUDA=ON|OFF`: Enable/disable CUDA support (default: OFF)
- `-DGMD_BUILD_PYTHON=ON|OFF`: Enable/disable Python bindings (default: OFF)
- `-DCMAKE_CUDA_ARCHITECTURES=80;86`: Specify GPU compute capability (default: auto-detect)

## Quick Start

### Running a Simple Example

After building, you can run one of the included examples:

```bash
cd build
./gmd ../examples/lj_fluid/config.txt
```

Available examples:
- `examples/lj_fluid/`: Lennard-Jones fluid simulation
- `examples/water_box/`: Water system with classical force field
- `examples/mlff_demo/`: Machine learning force field demonstration

### Configuration Files

Configuration files are located in `configs/examples/` and `configs/testcases/`. Create a new configuration file or modify existing ones to customize your simulation:

```text
[Simulation]
timesteps = 1000
dt = 0.001

[System]
initial_config = path/to/structure.gro

[ForceField]
type = classical  # or ml_model
```

### Testing

Run the test suite to verify installation:

```bash
cd build
ctest --parallel $(nproc)
```

Run specific test categories:
```bash
ctest -R unit       # Unit tests only
ctest -R integration # Integration tests only
ctest -R performance # Performance benchmarks
```

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

## Project Structure

The following structure organizes the GMD project across multiple layers:

```text
GMD/
├── CMakeLists.txt              # Main build configuration
├── README.md                   # This file
├── LICENSE                     # License information
├── .git/                       # Version control
├── .gitignore                  # Git ignore rules
│
├── cmake/                      # CMake configuration modules
│   ├── Dependencies.cmake      # External dependency configuration
│   ├── CompilerOptions.cmake   # Compiler flags and standards
│   └── CUDAOptions.cmake       # CUDA toolchain configuration
│
├── docs/                       # Technical documentation
│   ├── architecture/           # System architecture diagrams and explanations
│   ├── algorithms/             # Algorithm descriptions and implementations
│   └── design-notes/           # Design decisions and rationale
│
├── include/gmd/                # Public header files (interfaces)
│   ├── core/                   # Simulation control interface
│   ├── system/                 # System state and topology interface
│   ├── force/                  # Force field backend abstraction
│   ├── neighbor/               # Neighbor list construction interface
│   ├── integrator/             # Time integration interface
│   ├── runtime/                # Runtime context and resource management
│   ├── io/                     # I/O and configuration interface
│   ├── ml/                     # Machine learning force field interface
│   ├── cuda/                   # CUDA kernel dispatch interface
│   └── utils/                  # Utility functions and helpers
│
├── src/                        # Implementation files
│   ├── core/                   # Simulation lifecycle and step control
│   ├── system/                 # System state and topology implementation
│   ├── force/                  # Classical force field implementation
│   ├── neighbor/               # Neighbor list construction and maintenance
│   ├── integrator/             # Time integration algorithms
│   ├── runtime/                # CUDA context and resource management
│   ├── io/                     # Configuration loading, trajectory output, checkpoints
│   ├── ml/                     # ML force provider and model adapters
│   ├── cuda/                   # CUDA kernel wrappers and utilities
│   └── utils/                  # Error handling, profiling, math utilities
│
├── kernels/                    # CUDA kernel implementations
│   ├── neighbor/               # Cell list and neighbor list kernels
│   ├── pair/                   # Pairwise interaction kernels (LJ, Coulomb, etc.)
│   ├── bonded/                 # Bond, angle, and dihedral kernels
│   ├── integrate/              # Velocity and position update kernels
│   ├── reduce/                 # Energy, temperature, virial reduction kernels
│   └── features/               # ML feature extraction kernels
│
├── app/                        # Main application
│   └── gmd_main.cpp            # Command-line entry point
│
├── configs/                    # Configuration templates and examples
│   ├── examples/               # Example configuration files
│   └── testcases/              # Test case configurations
│
├── examples/                   # Sample simulation projects
│   ├── lj_fluid/               # Lennard-Jones fluid in periodic box
│   ├── water_box/              # Water system example with classical FF
│   └── mlff_demo/              # Machine learning force field demonstration
│
├── tests/                      # Test suite
│   ├── unit/                   # Unit tests for individual modules
│   ├── integration/            # Integration tests for module interactions
│   ├── regression/             # Regression tests for correctness validation
│   └── performance/            # Performance benchmarks
│
├── benchmarks/                 # Performance measurement tools
│   ├── single_gpu/             # Single GPU performance tests
│   └── ml_inference/           # ML model inference performance tests
│
├── third_party/                # External libraries and dependencies
│   └── (various)               # Third-party library integrations
│
├── python/                     # Python bindings and interface
│   └── gmd/                    # Python package for GMD
│
└── build/                      # Build output directory (generated)
    ├── CMakeFiles/             # CMake generated files
    ├── gmd                     # Compiled executable
    ├── CMakeLists.txt          # Generated
    └── Makefile                # Generated build rules
```

## Folder Descriptions

### Root Level Files

- **CMakeLists.txt**: Main build configuration file that defines compilation targets, dependencies, and build options
- **README.md**: Project overview, architecture principles, and development guidelines
- **LICENSE**: Project licensing information

### `cmake/` - Build System Configuration

Modular CMake configuration files for managing dependencies and compiler settings:

- **Dependencies.cmake**: Configures external libraries (CUDA, PyTorch, ONNX, etc.)
- **CompilerOptions.cmake**: Sets C++ standard, optimization flags, and compiler warnings
- **CUDAOptions.cmake**: CUDA toolkit version, compute capability targets, and GPU memory settings

### `docs/` - Technical Documentation

- **architecture/**: System design, component interactions, and module dependencies
- **algorithms/**: Detailed descriptions of numerical algorithms (integration, neighbor lists, force calculation)
- **design-notes/**: Rationale behind design decisions and alternative approaches considered

### `include/gmd/` - Public API Headers

Interface definitions for all modules. Organized by domain:

- **core/**: `Simulation`, `StepController` - main simulation loop control
- **system/**: `System`, `ParticleData`, `Box` - system state and metadata
- **force/**: `ForceProvider`, `ClassicalForceProvider` - unified force field abstraction
- **neighbor/**: `NeighborBuilder`, `NeighborList` - spatial search and list management
- **integrator/**: `Integrator`, `VelocityVerletIntegrator` - time stepping algorithms
- **runtime/**: `RuntimeContext`, `CudaContext` - execution resource management
- **io/**: `ConfigLoader`, `TrajectoryWriter` - input/output and file handling
- **ml/**: `MLForceProvider`, `ModelRuntimeAdapter` - ML model integration
- **cuda/**: `KernelDispatch` - CUDA kernel launching interface
- **utils/**: `Status`, `Profiler`, `Timer` - shared utilities

### `src/` - Implementation Files

Parallel structure to `include/` containing the concrete implementations of all public interfaces.

### `kernels/` - CUDA Kernel Code

High-performance device code organized by functionality:

- **neighbor/**: Cell list construction, neighbor list generation
- **pair/**: Lennard-Jones, Coulomb, and other pairwise forces
- **bonded/**: Bond stretches, angle bending, dihedral rotations
- **integrate/**: Velocity Verlet position/velocity updates
- **reduce/**: Parallel reductions for energy, temperature, virial calculations
- **features/**: ML input feature extraction from atomic positions and neighbors

### `app/` - Application Entry Points

- **gmd_main.cpp**: Command-line interface and main execution entry point

### `configs/` - Configuration Data

- **examples/**: Template configuration files showing all available options
- **testcases/**: Pre-configured test scenarios for regression testing

### `examples/` - Sample Projects

Runnable example simulations demonstrating different features:

- **lj_fluid/**: Basic Lennard-Jones fluid system in periodic boundary conditions
- **water_box/**: Water molecules with classical force field parameters
- **mlff_demo/**: Integration of machine learning force field models

### `tests/` - Testing Framework

- **unit/**: Tests for individual module functionality in isolation
- **integration/**: Tests for interactions between multiple modules
- **regression/**: Golden-reference tests to catch unintended behavioral changes
- **performance/**: Timing and throughput benchmarks

### `benchmarks/` - Performance Analysis

- **single_gpu/**: Execution time and memory benchmarks on single GPU
- **ml_inference/**: ML model inference throughput and latency measurements

### `third_party/` - External Dependencies

Integrated external libraries like CUDA headers, JSON parsers, linear algebra libraries, etc.

### `python/` - Python Bindings

Python package providing high-level interface to C++ simulation engine:

- Workflow orchestration
- Post-processing and analysis
- Interactive parameter exploration

### `build/` - Build Artifacts (Generated)

Contains all compiled objects, executables, and intermediate build files. Typically ignored in version control.

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

1. ✓ Build the C++ core, runtime layer, and CLI.
2. ✓ Implement CPU neighbor list (cell list + Verlet skin), LJ pair force, and Velocity Verlet integrator.
3. ✓ Run a simple Lennard-Jones fluid case with correct physics (shifted potential, MIC, Newton III).
4. ✓ Add thermostat (velocity rescaling, Nosé-Hoover) and barostat (Berendsen).
5. ✓ Add multi-element LJ support (Lorentz-Berthelot mixing rules) and trajectory output.
6. Add CUDA neighbor list, pair force, and integration kernels.
7. Add checkpoint, logging, and regression tests.
8. Add the `ForceProvider` abstraction for ML backends.
9. Integrate one stable ML runtime backend (TorchScript or ONNX).
10. Add optional Python bindings with `pybind11`.

## Recommended First Deliverable

The first stable version of GMD should support:

- ✓ periodic boundary conditions
- ✓ Verlet neighbor lists (cell list, O(N))
- ✓ Lennard-Jones nonbonded interactions (single and multi-element)
- ✓ Velocity Verlet integration
- ✓ thermostat (velocity rescaling, Nosé-Hoover) and barostat (Berendsen)
- ✓ trajectory output (extended XYZ + energy log)
- ✓ one pluggable force backend interface (`ForceProvider`)
- ☐ single GPU execution (CUDA backend not yet implemented)
- ☐ checkpoint and restart

## Development Workflow

### Code Organization Principles

1. **Header-Implementation Separation**: All public interfaces reside in `include/gmd/` with implementations in `src/`
2. **Module Independence**: Each module should have minimal coupling with others through well-defined interfaces
3. **CUDA Separation**: Device code is isolated in `kernels/` with C++ wrappers in `src/cuda/`
4. **Test Proximity**: Tests live alongside code with clear naming conventions

### Coding Standards

- C++20 standard required
- Use `const` correctness and modern C++ idioms
- Prefer `std::span` and `std::vector` over raw pointers
- Document public APIs with inline comments
- Use CUDA error checking with custom assertions

### Adding New Features

1. Define the interface in `include/gmd/{module}/`
2. Implement core logic in `src/{module}/`
3. Add CUDA kernels to `kernels/` if GPU acceleration needed
4. Write tests in `tests/{unit|integration}/`
5. Update `CMakeLists.txt` with new source files
6. Document design choices in `docs/design-notes/`

## Contributing

We welcome contributions! Please follow these guidelines:

1. Fork the repository and create a feature branch
2. Ensure your code follows the coding standards above
3. Write tests for new functionality
4. Run the full test suite before submitting a pull request
5. Update documentation as needed
6. Include a clear commit message describing your changes

For major changes, please open an issue first to discuss the proposed changes.

## Project Status

**Current Version**: 0.2.0 (Active Development)

### Implemented ✓

| Module | Status | Key Files |
|--------|--------|-----------|
| Data / System | ✓ Complete | `box.hpp`, `system.hpp` |
| Boundary Conditions | ✓ Complete | `periodic_boundary`, `minimum_image` |
| Neighbor List | ✓ Complete | `verlet_neighbor_builder` (cell list + Verlet skin) |
| Classical Force Field (LJ) | ✓ Complete | `classical_force_provider` — shifted potential, MIC, Newton III, O(N) neighbor-list path |
| Multi-element LJ | ✓ Complete | Lorentz-Berthelot mixing rules; `element` directive in `.ff` files |
| Integrator (Velocity Verlet) | ✓ Complete | `velocity_verlet_integrator` |
| Thermostat | ✓ Complete | `VelocityRescalingThermostat`, `NoseHooverThermostat` (VVNH splitting) |
| Barostat | ✓ Complete | `BerendsenBarostat` (kinetic + virial pressure estimate) |
| Velocity Initializer | ✓ Complete | Maxwell-Boltzmann sampling, COM velocity removal, temperature scaling |
| Input Parsing | ✓ Complete | `ConfigLoader` — xyz (5/8 col), run config, `.ff` force field files |
| Trajectory Output | ✓ Complete | `TrajectoryWriter` — extended XYZ + energy/temperature log |
| Main Loop | ✓ Complete | `Simulation` — initialize, step, run; neighbor rebuild check; t=0 force bug fixed |

### Not Yet Implemented

- **Bonded interactions**: bond stretches, angle bending, dihedral/improper torsions
- **Constraints**: SHAKE/RATTLE for rigid bonds
- **CUDA backend**: all computation currently runs on CPU
- **ML force fields**: `MLForceProvider` / `ModelRuntimeAdapter` interfaces exist but are not implemented
- **Checkpoint / restart**: no save/restore of simulation state
- **Python bindings**: pybind11 layer not yet built
- **Tests**: unit and regression test suite not yet populated

---

## Development Log

> 以下为各模块实现进度的详细记录（中文）。

### 数据模块 ✓
对应文件：`box.hpp`、`system.hpp`、`config_loader.hpp/.cpp`。

`Box` 保存盒长和半盒长；`System` 保存质量、原子类型、坐标、速度、受力、势能、邻居表；`ConfigLoader::load_xyz` 把原子数、坐标、质量、盒子尺寸、原子类型一次性读进 `System`。

### 初始化模块 ✓
对应文件：`initializer.hpp/.cpp`。

- 随机速度初始化（Maxwell-Boltzmann 高斯采样）
- 输入速度模式（从 xyz 8 列格式读取 vx/vy/vz）
- 去掉质心速度（可开关）
- 按目标温度缩放（含自由度修正，去质心后使用 3N-3）

### 边界条件模块 ✓
对应文件：`boundary/periodic_boundary`、`boundary/minimum_image`。

- PBC：`wrap_coordinate`、`wrap_position`、`wrap_positions`
- MIC：`apply_minimum_image_component`、`apply_minimum_image`
- 均已接入 `VelocityVerletIntegrator`（坐标回卷）和 `ClassicalForceProvider::compute`（位移 MIC 处理）

### 邻居表模块 ✓
对应文件：`neighbor_builder.hpp`（抽象接口）、`verlet_neighbor_builder.hpp/.cpp`。

- **建表（`rebuild()`）**：O(N) cell list 算法——模拟盒划分为边长 ≥ r_list 的格子，linked list 填格，只检查 3×3×3 个邻居格，半对 CSR 存储（j > i，Newton III）
- **重建判断（`needs_rebuild()`）**：任意原子位移超过 r_skin/2 即触发重建
- `ForceRequest` 透传 `NeighborList`，力场自动走 O(N) 路径

### 力场模块 ✓
对应文件：`force_provider.hpp`、`classical_force_provider.hpp/.cpp`，以及 `config_loader`（读取接口）。

**LJ 物理实现：**
- Lennard-Jones (12-6)，含截断和 shifted potential（截断处能量连续）
- MIC 接入，Newton III 配对，O(N) 邻居表路径 / O(N²) 回退

**多元素支持（Lorentz-Berthelot）：**
- `LJForceFieldConfig` 包含 `elements` 向量和 `element_name_to_type` 映射
- `pair_params(type_i, type_j)` 自动计算交叉参数（几何平均 ε，算术平均 σ）
- `ClassicalForceProvider` 以 `pair_table_`（type×type 预计算 `PairCache`）索引查表
- `.ff` 文件格式：支持 `element Ar epsilon X sigma Y` 多元素指令；旧格式向下兼容
- 示例文件：`configs/examples/lj_argon.ff`（Ar），`configs/examples/lj_ne_ar.ff`（Ne/Ar 混合）

### 积分模块 ✓
对应文件：`velocity_verlet_integrator.hpp/.cpp`、`thermostat.hpp/.cpp`、`velocity_rescaling_thermostat`、`nose_hoover_thermostat`、`barostat.hpp`、`berendsen_barostat.hpp/.cpp`。

**Velocity Verlet 积分器：**
步序：thermostat pre-half-kick → 第一半步力+位移 → 重新计算力 → 第二半步速度 → thermostat post-half-kick + full apply → barostat

**Thermostat：**
- `VelocityRescalingThermostat`：每步将全体速度缩放至目标温度，自由度 3N-3
- `NoseHooverThermostat`：VVNH 劈裂，摩擦变量 ξ 更新，Q 惰性初始化（`Q = dof·kB·T·τ²`）

**Barostat：**
- `BerendsenBarostat`：按 `mu = cbrt(1 − β·dt/τ_P·(P_target − P))` 等比缩放盒长和坐标；压强由动能+维里估算

**接口（在 `VelocityVerletIntegrator` 上）：**
`set_thermostat()`、`set_target_temperature()`、`set_barostat()`、`set_target_pressure()`

### 输入模块 ✓
对应文件：`config_loader.hpp/.cpp`。

- `xyz` 5 列（type x y z mass）或 8 列（+vx vy vz）格式；第 0 列类型写入 `System::atom_types_`
- `run` 文件：步数、时间步、温度、`velocity_init`、`velocity_seed`、`remove_com_velocity`
- `.ff` 力场参数文件（单/多元素）

### 输出模块 ✓
对应文件：`trajectory_writer.hpp/.cpp`。

- **XYZ 轨迹文件**（`*.xyz`）：extended XYZ 格式，注释行含 step/time/PE/KE/T
- **能量日志文件**（`*.log`）：`step  time[fs]  PE[eV]  KE[eV]  E_total[eV]  T[K]`
- `write_frame_if()` 按步长间隔过滤；析构时自动关闭并刷新

### 主循环模块 ✓
对应文件：`simulation.hpp/.cpp`，入口 `gmd_main.cpp`。

- `initialize()`：速度初始化 → 邻居表首次建立 → t=0 初始力计算
- `step()`：自动检查 `needs_rebuild()` 并在必要时重建邻居表
- `run()`：调用 `step()` 指定次数

### 已解决的问题
- ~~初始力未在时间循环前计算~~（已修复）
- ~~截断势未做 shifted potential~~（已修复）
- ~~受力计算 O(N²) 双循环无邻居表~~（已修复：cell list O(N)）
- ~~无 thermostat / barostat / 多元素力场~~（已修复）
- ~~输出模块未实现~~（已修复）
- ~~未存储原子类型、解析错误信息不含行号~~（已修复）

### 尚未实现
- 键合相互作用（bond / angle / dihedral）
- 约束（SHAKE / RATTLE）
- CUDA 后端
- ML 力场运行时
- Checkpoint / restart
- Python bindings
- 单元测试和回归测试

## License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.

## Citation

If you use GMD in your research, please cite:

```bibtex
@software{gmd2024,
  title={GMD: Molecular Dynamics with CUDA and Machine Learning},
  author={Contributors},
  year={2024},
  url={https://github.com/kiana-0423/GMD}
}
```

## Support and Contact

For questions, bug reports, or feature requests, please:

- Open an issue on GitHub
- Check existing documentation in `docs/`
- Review design notes in `docs/design-notes/`

## Acknowledgments

GMD builds upon established molecular dynamics methodologies and benefits from the broader open-source scientific computing community.
