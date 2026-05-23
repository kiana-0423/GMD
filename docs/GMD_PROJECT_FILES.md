# GMD 项目完整文件清单与说明

> **项目:** GMD (Generalized Molecular Dynamics) v2.2  
> **语言:** C++20 | **构建:** CMake | **许可证:** MIT  
> **生成日期:** 2026-05-23（目录重组后）

---

## 目录

1. [项目结构总览](#1-项目结构总览)
2. [根目录文件](#2-根目录文件)
3. [app/ — 应用程序入口](#3-app--应用程序入口)
4. [include/gmd/ — 公共 API 头文件](#4-includegmd--公共-api-头文件)
5. [src/ — 源代码实现](#5-src--源代码实现)
6. [tests/ — 测试文件](#6-tests--测试文件)
7. [examples/ — 示例](#7-examples--示例)
8. [cmake/ — CMake 模块](#8-cmake--cmake-模块)
9. [build/ & build-mpi-check/ — 构建目录](#9-build--build-mpi-check--构建目录)
10. [目录重组说明](#10-目录重组说明)
11. [统计汇总](#11-统计汇总)

---

## 1. 项目结构总览

> **v2.2 目录重组:** 原 12 个子目录合并为 6 个，详见[第 10 节](#10-目录重组说明)。

```
GMD/
├── CMakeLists.txt                 # 顶层构建配置
├── README.md                      # 项目说明
├── CHANGELOG.md                   # 版本更新日志
├── LICENSE                        # MIT 许可证
├── run_ethane_demo.sh             # 乙烷示例快速构建脚本
├── app/
│   └── gmd_main.cpp               # 主程序入口 (~700行)
├── include/gmd/                   # 公共 API 头文件 (34个, 6个子目录)
│   ├── core/          (2)         #   模拟核心 + 运行时上下文
│   ├── force/         (10)        #   力场计算 + ML 力场
│   ├── integrator/    (8)         #   积分器/热浴/压浴
│   ├── io/            (2)         #   输入输出
│   ├── parallel/      (4)         #   MPI 并行化
│   └── system/        (8)         #   系统数据结构 + 边界条件 + 近邻列表
├── src/                           # 源代码实现 (27个, 6个子目录)
│   ├── core/          (2)
│   ├── force/         (7)
│   ├── integrator/    (7)
│   ├── io/            (2)
│   ├── parallel/      (4)
│   └── system/        (5)
├── tests/                         # 测试文件
│   ├── *.cpp                      #   单元测试程序 (4)
│   ├── smoke_*.run                #   冒烟测试运行配置 (10+)
│   ├── smoke_*.xyz                #   冒烟测试坐标文件 (7+)
│   ├── smoke_*.ff                 #   冒烟测试力场参数 (1)
│   └── smoke_*.top                #   冒烟测试拓扑文件 (1)
├── examples/ethane_demo/          # 乙烷分子模拟示例
│   ├── ethane_demo.cpp
│   ├── ethane.xyz / .run / .ff / .top
├── cmake/                         # CMake 辅助模块 (5)
│   ├── CompilerOptions.cmake
│   ├── CUDAOptions.cmake
│   ├── MPIOptions.cmake
│   ├── RunMpiConsistency.cmake
│   └── RunMpiLJConsistency.cmake
├── build/                         # 构建输出目录
├── build-mpi-check/               # MPI 构建输出目录
└── docs/                          # 文档
    ├── GMD_SOURCE_CODE_ARCHITECTURE.md
    └── GMD_PROJECT_FILES.md       # ← 本文件
```

---

## 2. 根目录文件

### `CMakeLists.txt`
**顶层 CMake 构建配置文件**

- 要求 C++20 标准
- 定义构建选项: `GMD_ENABLE_CUDA`, `GMD_ENABLE_PYTHON`, `GMD_ENABLE_TORCH`, `GMD_ENABLE_MPI`
- 构建目标: **`gmd_core`** 静态库 + **`gmd`** 可执行文件
- 注册 **40+ 个 CTest 测试目标**
- 构建时自动复制测试数据文件到 `build/tests/`

### `README.md`
项目说明文档，包含 v2.2 版本亮点、构建指南、快速开始和 MPI 并行运行说明。

### `CHANGELOG.md`
详细版本更新日志（v2.2 → v1.0）。

### `LICENSE`
MIT 开源许可证，Copyright © 2026 Ēlýsion。

### `run_ethane_demo.sh`
乙烷示例快速构建与运行脚本，无需 CMake，直接用 C++20 编译器编译所有源文件。

---

## 3. app/ — 应用程序入口

### `app/gmd_main.cpp`
**主模拟器程序（~700 行）**

核心函数:

| 函数 | 功能 |
|------|------|
| `parse_command_line()` | 解析命令行参数：`xyz run [ff] [top]` + `--np N` `--proc-grid Px Py Pz` |
| `detect_force_field_file_kind()` | 自动检测 LJ 或分子力场类型 |
| `keep_rank_local_atoms()` | MPI 域分解后过滤本进程拥有的原子 |

**主工作流:**
1. 初始化 MPI 环境 → 加载配置 → 加载力场/拓扑 → 加载坐标
2. 选择力场提供者（LJ / 分子 / ML / 默认 Ar）+ 库仑方法（Ewald / PME）
3. 创建近邻列表 → 设置 MPI 域分解 → 配置积分器 + 热浴/压浴
4. 初始化速度 → MD 循环 → 输出轨迹/日志（仅 rank-0 I/O）

**输出:** `output.xyz`（轨迹）+ `output.log`（能量日志）

---

## 4. include/gmd/ — 公共 API 头文件

> 共 **34 个头文件**，按 6 个模块组织。

### 4.1 core/ — 模拟核心 (2 个)

| 文件 | 说明 |
|------|------|
| `simulation.hpp` | 主 Simulation 类，MD 模拟总控制器（PIMPL 模式） |
| `runtime_context.hpp` | `RuntimeContext` 类，提供 rank/size 等运行时信息查询 |

### 4.2 force/ — 力场计算 (10 个)

| 文件 | 说明 |
|------|------|
| `force_provider.hpp` | 抽象力场接口：`compute()` / `initialize()` / `finalize()` |
| `classical_force_provider.hpp` | Lennard-Jones 12-6 势，支持势能截断偏移 |
| `bonded_force_provider.hpp` | 分子内力场：键伸缩、键角弯曲、二面角、异常二面角 |
| `bonded_params.hpp` | 参数结构体：`BondTerm` / `AngleTerm` / `DihedralTerm` / `ImproperTerm` |
| `ewald_force_provider.hpp` | Ewald 求和方法计算长程库仑力 |
| `pme_force_provider.hpp` | Particle-Mesh Ewald（3D FFT 加速） |
| `composite_force_provider.hpp` | 组合多个力场提供者（如 LJ + PME） |
| `ml_force_provider.hpp` | 包装 ML 模型的力场提供者 |
| `model_runtime_adapter.hpp` | 抽象 ML 推理后端接口 |
| `torchscript_adapter.hpp` | LibTorch TorchScript 具体实现 |

### 4.3 integrator/ — 积分器与控温控压 (8 个)

| 文件 | 说明 |
|------|------|
| `integrator.hpp` | 抽象积分器接口 |
| `velocity_verlet_integrator.hpp` | 二阶 Velocity-Verlet 算法 |
| `thermostat.hpp` | 抽象恒温器基类 |
| `velocity_rescaling_thermostat.hpp` | 简单速度重缩放恒温器 |
| `nose_hoover_thermostat.hpp` | Nosé-Hoover 扩展系统恒温器（严格 NVT） |
| `barostat.hpp` | 抽象恒压器基类 |
| `berendsen_barostat.hpp` | Berendsen 弱耦合恒压器 |
| `mc_barostat.hpp` | 蒙特卡洛 NPT 恒压器 |

### 4.4 io/ — 输入输出 (2 个)

| 文件 | 说明 |
|------|------|
| `config_loader.hpp` | 解析 `.run` / `.ff` / `.xyz` / `.top` 配置文件 |
| `trajectory_writer.hpp` | 写入 `.xyz` 轨迹文件和 `.log` 能量日志 |

### 4.5 parallel/ — MPI 并行化 (4 个)

> 条件编译: `GMD_ENABLE_MPI`

| 文件 | 说明 |
|------|------|
| `mpi_environment.hpp` | RAII 风格的 MPI 初始化/终结管理 |
| `mpi_communicator.hpp` | MPI 通信：Allreduce、Ghost 交换、力反向累加、原子重分配 |
| `domain_decomposition.hpp` | 3D 空间域分解，含 `DomainInfo` 结构体 |
| `pme_parallel.hpp` | PME 铅笔分解（分布式 FFT 基础设施） |

### 4.6 system/ — 系统数据结构 (8 个)

| 文件 | 说明 |
|------|------|
| `system.hpp` | `System` 类：每原子状态（质量、电荷、坐标、速度、力、近邻列表） |
| `box.hpp` | `Box` 类：模拟盒子几何（边长、半边长） |
| `topology.hpp` | `Topology` 类：分子连接性（键、角、二面角、异常二面角） |
| `initializer.hpp` | `VelocityInitializer` 类：Maxwell-Boltzmann 速度初始化 |
| `periodic_boundary.hpp` | 将原子坐标包裹回主元胞（PBC wrapping） |
| `minimum_image.hpp` | 计算最近镜像位移向量（minimum image convention） |
| `neighbor_builder.hpp` | 抽象近邻列表构建器接口 |
| `verlet_neighbor_builder.hpp` | 基于元胞列表的 Verlet 近邻列表（O(N) + skin 距离 + 3D 镜像标志） |

---

## 5. src/ — 源代码实现

> 共 **27 个 .cpp 实现文件**，与 `include/gmd/` 中的头文件对应。

### 5.1 core/ (2 个)

| 文件 | 实现内容 |
|------|---------|
| `simulation.cpp` | Simulation 内部实现类（PIMPL），MD 主循环逻辑 |
| `runtime_context.cpp` | RuntimeContext 实现 |

### 5.2 force/ (7 个)

| 文件 | 实现内容 |
|------|---------|
| `classical_force_provider.cpp` | LJ 12-6 力计算、Lorentz-Berthelot 混合规则 |
| `bonded_force_provider.cpp` | 键/角/二面角/异常二面角的力与能量计算 |
| `ewald_force_provider.cpp` | Ewald 求和：实空间 + 倒空间 + 自能项 |
| `pme_force_provider.cpp` | PME：B-样条插值、3D FFT 网格电荷分配与力插值 |
| `composite_force_provider.cpp` | 组合多个 ForceProvider 的代理实现 |
| `ml_force_provider.cpp` | ML 力场包装器实现 |
| `model_runtime_adapter.cpp` | 抽象 ML 后端的默认实现 |
| `torchscript_adapter.cpp` | LibTorch 模型加载与推理（条件编译） |

### 5.3 integrator/ (7 个)

| 文件 | 实现内容 |
|------|---------|
| `velocity_verlet_integrator.cpp` | Velocity-Verlet 三步积分算法 |
| `thermostat.cpp` | Thermostat 基类实现 |
| `velocity_rescaling_thermostat.cpp` | 速度重缩放算法 |
| `nose_hoover_thermostat.cpp` | Nosé-Hoover 扩展系统动力学 |
| `barostat.cpp` | Barostat 基类实现 |
| `berendsen_barostat.cpp` | Berendsen 弱耦合压力控制 |
| `mc_barostat.cpp` | 蒙特卡洛体积移动 + 自适应步长 |

### 5.4 io/ (2 个)

| 文件 | 实现内容 |
|------|---------|
| `config_loader.cpp` | 解析 `.run` / `.ff` / `.xyz` / `.top` |
| `trajectory_writer.cpp` | XYZ 轨迹格式输出、能量日志输出 |

### 5.5 parallel/ (4 个)

| 文件 | 实现内容 |
|------|---------|
| `mpi_environment.cpp` | MPI_Init / MPI_Finalize RAII 包装 |
| `mpi_communicator.cpp` | 点对点 Ghost 交换、力累加、集合通信 |
| `domain_decomposition.cpp` | 3D 笛卡尔网格划分、原子分配与迁移 |
| `pme_parallel.cpp` | PME 分布式 FFT 铅笔分解框架 |

### 5.6 system/ (5 个)

| 文件 | 实现内容 |
|------|---------|
| `initializer.cpp` | Maxwell-Boltzmann 速度初始化、质心速度归零 |
| `minimum_image.cpp` | 最小镜像约定的具体计算 |
| `periodic_boundary.cpp` | PBC 坐标包裹的实现 |
| `verlet_neighbor_builder.cpp` | 元胞列表构建、Verlet 缓冲半径、增量重建判断 |

---

## 6. tests/ — 测试文件

### 6.1 单元测试 C++ 程序 (4 个)

| 文件 | 说明 |
|------|------|
| `mpi_periodic_1d.cpp` | 2 进程 1D 周期性测试：ghost 交换、LJ 力一致性、原子跨域迁移 |
| `mpi_domain_decomposition_3d.cpp` | 4/8 进程 3D 网格测试：网格创建与所有权、原子迁移 |
| `mpi_ghost_exchange_3d.cpp` | 8 进程 2×2×2 网格测试：面/边/角 ghost 交换、力反向累加 |
| `compare_energy_logs.cpp` | 工具程序：比较两个 `.log` 文件的势能/温度/能量漂移 |

### 6.2 冒烟测试配置文件

| 类别 | 文件 |
|------|------|
| 运行配置 (.run) | `smoke_lj.run`, `smoke_ewald.run`, `smoke_pme.run`, `smoke_mc_barostat.run`, `smoke_molecular.run`, `smoke_mpi_lj.run`, `smoke_mpi_ewald_consistency.run`, `smoke_mpi_pme_consistency.run`, `smoke_mpi_periodic_lj.run`, `smoke_mpi_boundary_migration.run` |
| 坐标文件 (.xyz) | `smoke_lj.xyz`, `smoke_ewald.xyz`, `smoke_mc_barostat.xyz`, `smoke_molecular.xyz`, `smoke_mpi_lj.xyz`, `smoke_mpi_ewald.xyz`, `smoke_mpi_periodic_lj.xyz`, `smoke_mpi_boundary_migration.xyz`, `smoke_mpi_molecular.xyz` |
| 力场/拓扑 | `smoke_molecular.ff`, `smoke_molecular.top` |

---

## 7. examples/ — 示例

### 7.1 examples/ethane_demo/ — 乙烷分子模拟

完整的乙烷（C₂H₆）OPLS-AA 力场模拟示例。

| 文件 | 说明 |
|------|------|
| `ethane_demo.cpp` | 程序化模拟设置代码（可独立编译运行） |
| `ethane.xyz` | 8 原子初始结构，盒子 20×20×20 Å |
| `ethane.run` | 运行参数：dt=1.0 fs，500 步，T=300K，Nosé-Hoover |
| `ethane.ff` | OPLS-AA 力场参数（CT 碳、HC 氢、键、角、二面角） |
| `ethane.top` | 拓扑连接：7 键、12 角、9 二面角 |

---

## 8. cmake/ — CMake 模块

| 文件 | 说明 |
|------|------|
| `CompilerOptions.cmake` | `gmd_apply_default_warnings()`：MSVC `/W4`，GCC/Clang `-Wall -Wextra -Wpedantic` |
| `CUDAOptions.cmake` | `gmd_configure_cuda_target()`：CUDA C++17 标准（预留） |
| `MPIOptions.cmake` | `gmd_configure_mpi_target()`：链接 `MPI::MPI_CXX`，定义 `GMD_ENABLE_MPI` |
| `RunMpiConsistency.cmake` | MPI 一致性测试框架：串行参考 → MPI 多进程 → 比较能量日志 |
| `RunMpiLJConsistency.cmake` | LJ 力场专用 MPI 一致性测试框架 |

---

## 9. build/ & build-mpi-check/ — 构建目录

### build/
标准构建输出目录：`CMakeCache.txt`, `Makefile`, `gmd` 可执行文件, `libgmd_core.a` 静态库, `tests/` 测试数据。

### build-mpi-check/
MPI 专用构建配置目录，含独立的 `CMakeCache.txt` 和 `DartConfiguration.tcl`。

---

## 10. 目录重组说明

> **v2.2 重组 (2026-05-23):** 原 12 个子目录合并为 6 个。

| 原目录 | 合并到 | 原因 |
|--------|--------|------|
| `boundary/` (2文件) | `system/` | PBC/最小镜像是系统几何属性 |
| `neighbor/` (2文件) | `system/` | 近邻列表是系统数据结构的一部分 |
| `runtime/` (1文件) | `core/` | 运行时上下文辅助模拟核心 |
| `ml/` (3文件) | `force/` | ML 力场是力场的一种实现方式 |
| `cuda/` (空) | 删除 | 空预留目录 |
| `utils/` (空) | 删除 | 空预留目录 |

**代码和架构完全不变**，仅移动文件并更新 `#include` 路径。

---

## 11. 统计汇总

| 类别 | 数量 |
|------|------|
| 公共头文件 (`include/gmd/`) | 34 |
| 源代码文件 (`src/`) | 27 |
| 子目录数 (include + src) | 6 + 6 |
| 应用程序入口 (`app/`) | 1 |
| 测试 C++ 程序 | 4 |
| 冒烟测试配置文件 | 15+ |
| CMake 模块 | 5 |
| 示例文件 | 5 |
| CTest 注册测试目标 | 40+ |

### 模块文件分布

```
force/          ██████████████████ 17 (10h + 7cpp)
system/         █████████████ 13 (8h + 5cpp)
integrator/     ███████████████ 15 (8h + 7cpp)
parallel/       ████████ 8 (4h + 4cpp)
core/           ████ 4 (2h + 2cpp)
io/             ████ 4 (2h + 2cpp)
```

---

## 附录: 条件编译宏

| 宏 | 控制范围 | 涉及文件 |
|----|---------|---------|
| `GMD_ENABLE_MPI` | MPI 并行化 | `parallel/` 全部, `app/gmd_main.cpp`, `core/simulation.cpp` |
| `GMD_ENABLE_TORCH` | ML 力场 | `force/torchscript_adapter.cpp`, `force/ml_force_provider.*` |
| `GMD_ENABLE_CUDA` | GPU 加速 | 预留 |
| `BUILD_TESTING` | 测试注册 | `CMakeLists.txt`（注册 40+ CTest 目标） |

---

> **文档生成:** 2026-05-23 | **基于:** GMD v2.2 目录重组后完整项目源码分析
