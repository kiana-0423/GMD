# File Description

This file is from https://github.com/brucefan1983/Molecular-Dynamics-Simulation.git

Main instruction file for the Molecular Dynamics Simulation project. It contains the build configuration and instructions for compiling the project using CMake. The file specifies the source files, dependencies, and build targets for the project.

Including the following modules:

- 数据模块:对应config_loader.hpp 和 box.hpp,负责加载配置文件和定义模拟盒子,并在system模块进行读取；
- 初始化模块:initializer.hpp 和 initializer.cpp,主要负责速度初始化（随机高斯采样/输入速度、去质心速度、温度缩放、动能计算）。
- 边界条件模块:periodic_boundary.hpp/.cpp 和 minimum_image.hpp/.cpp。前者负责周期边界回卷（PBC），后者负责最小镜像位移处理（MIC）。
- 力场模块
- 积分模块
- 输入模块:config_loader.hpp 和 config_loader.cpp,负责读取 xyz/run 输入（支持 xyz 的 5 列 type x y z mass 或 8 列 type x y z vx vy vz mass，支持 velocity_init、velocity_seed、remove_com_velocity 等初始化参数）
- 输出模块:负责输出模拟结果,如坐标、能量等,并提供接口供用户调用。
- 主循环模块:simulation.hpp 和 simulation.cpp（入口在 gmd_main.cpp）,负责主循环的实现，调用其他模块完成模拟过程。

日志：

- 数据模块：基本完成  
  对应文件是 box.hpp、system.hpp、config_loader.hpp 和 config_loader.cpp。  
  现在 `Box` 已经能保存盒长和半盒长，`System` 已经能保存质量、坐标、速度、受力、势能，`ConfigLoader::load_xyz` 能把原子数、坐标、质量、盒子尺寸读进 `System`。  
  这一块是当前最完整的模块之一。

- 初始化模块：基本完成（速度初始化）  
  对应文件是 initializer.hpp 和 initializer.cpp。  
  当前已实现：
  - 随机速度初始化（Maxwell-Boltzmann 高斯采样）
  - 输入速度模式（从 xyz 8 列格式读取 vx/vy/vz）
  - 去掉质心速度（可开关）
  - 按目标温度缩放（含自由度修正，去质心后使用 3N-3）
  - 计算动能  
  说明：原子位置初始化仍由 config_loader.cpp 从输入文件加载，因此该模块应定义为“速度初始化模块”。

- 边界条件模块：已拆分为独立文件，完成并已接入力场  
  对应文件是 boundary/periodic_boundary.hpp、boundary/periodic_boundary.cpp、boundary/minimum_image.hpp 和 boundary/minimum_image.cpp。  
  当前已实现：
  - PBC：`wrap_coordinate`、`wrap_position`、`wrap_positions`
  - MIC：`apply_minimum_image_component`、`apply_minimum_image`  
  当前接线状态：
  - `VelocityVerletIntegrator` 已使用 `wrap_position` 做坐标周期回卷
  - `apply_minimum_image` 已接入 `ClassicalForceProvider::compute`，每对原子的位移向量均经过 MIC 处理

- 力场模块：LJ 力场物理实现已完成，力场读取接口已完成，邻居表已接入  
  对应文件是 force_provider.hpp、classical_force_provider.hpp、classical_force_provider.cpp，以及 config_loader.hpp/config_loader.cpp（读取接口）。  
  当前已实现：
  - `ClassicalForceProvider::compute` 实现了完整的 Lennard-Jones (12-6)，含截断（cutoff）和势能平移（shifted potential，截断处能量连续）
  - 使用 `apply_minimum_image` 处理周期边界，MIC 已真正接入受力计算
  - 牛顿第三定律配对更新，每对原子仅计算一次
  - 当 `ForceRequest::neighbor_list` 有效时走 O(N) 邻居表路径，否则回退 O(N²) 双循环
  - 构造函数支持显式参数或从 `LJForceFieldConfig` 直接构造
  - 参数访问器：`epsilon()`、`sigma()`、`cutoff()`、`energy_shift()`
  - `set_params()` 支持运行时替换参数并重算缓存
  - `ConfigLoader::load_force_field()` 读取 `.ff` 格式配置文件（`force_field lj`、`epsilon`、`sigma`、`cutoff`），含合法性校验
  - 提供示例配置文件 `configs/examples/lj_argon.ff`（氩原子 LJ 参数）
  多元素支持（本次新增）：
  - `LJForceFieldConfig` 重构为含 `elements` 向量和 `element_name_to_type` 映射
  - `pair_params(type_i, type_j)` 按 Lorentz-Berthelot 规则（几何平均 ε，算术平均 σ）自动计算交叉对参数
  - `ClassicalForceProvider` 重构为 `pair_table_`（type×type 预计算 `PairCache`），`compute()` 按原子类型索引查表
  - `.ff` 文件新增 `element <name> epsilon <val> sigma <val>` 多元素指令；单元素旧格式向下兼容
  - 提供示例 `configs/examples/lj_ne_ar.ff`（Ne/Ar 二元混合）

- 积分模块：已完成，支持 thermostat 和 barostat 接入  
  对应文件是 integrator.hpp、velocity_verlet_integrator.hpp/.cpp、thermostat.hpp/.cpp、velocity_rescaling_thermostat.hpp/.cpp、nose_hoover_thermostat.hpp/.cpp、barostat.hpp、berendsen_barostat.hpp/.cpp。  
  当前已实现：
  - 抽象积分器接口
  - `VelocityVerletIntegrator` 具体实现（坐标+速度更新、周期回卷、邻居表透传）
  - **Thermostat 抽象接口**：`apply()`（全步，用于速度缩放类）和 `apply_half_kick()`（半步，用于扩展系统类）
  - **VelocityRescalingThermostat**：每步将全体速度缩放至目标温度，自由度修正为 3N−3
  - **NoseHooverThermostat**：VVNH 劈裂——两次 `apply_half_kick()` 包夹 VV 步，摩擦变量 ξ 更新；Q 首次调用时惰性初始化（`Q = dof·kB·T·τ²`）
  - **Barostat 抽象接口**：`apply(system, dt, P_target, virial_trace)`
  - **BerendsenBarostat**：每步末尾按 `mu = cbrt(1 − β·dt/τ_P·(P_target − P))` 等比缩放盒长和坐标；压强由动能+维里估算
  - `VelocityVerletIntegrator` 新增 `set_thermostat()`、`set_target_temperature()`、`set_barostat()`、`set_target_pressure()` 接口
  - 步序：thermostat pre-half-kick → 第一半步力+位移 → 重新计算力 → 第二半步速度 → thermostat post-half-kick + full apply → barostat
  - 辅助函数：`compute_twice_ke()`、`temperature_from_twice_ke()`、`kBoltzmann = 8.617333262e-5 eV/K`

- 邻居表模块：已完成（cell list + Verlet skin）  
  对应文件是 neighbor_builder.hpp（抽象接口）、verlet_neighbor_builder.hpp 和 verlet_neighbor_builder.cpp。  
  当前已实现：
  - `VerletNeighborBuilder`：继承 `NeighborBuilder` 抽象接口
  - **建表（`rebuild()`）**：O(N) cell list 算法——将模拟盒划分为边长 ≥ r_list 的格子，linked list 填格，只检查 3×3×3=27 个邻居格，半对存储（j > i，利用 Newton III），MIC 距离筛选，CSR 压缩存储到 `System::NeighborList`
  - **重建判断（`needs_rebuild()`）**：O(N) 最大位移准则——任意原子位移超过 r_skin/2 即触发重建
  - **`initialize()`**：调用第一次 `rebuild()`
  - `NeighborList` 以 CSR 格式存入 `System`（`counts`、`offsets`、`neighbors`、`ref_coordinates`、`valid`）
  - `ForceRequest` 增加 `neighbor_list` 字段，积分器和主循环均已透传
  参数：`r_cut`（与力场 cutoff 一致）、`r_skin`（典型 1.5 Å）

- 输入模块：已完成并扩展初始化相关配置及力场读取  
  对应文件是 config_loader.hpp 和 config_loader.cpp。  
  目前能读取：
  - `xyz` 文件中的原子数、坐标、质量、盒子尺寸
  - `xyz` 8 列格式中的输入速度 vx/vy/vz
  - `run` 文件中的步数、时间步长、温度
  - `run` 中的初始化控制项：`velocity_init`、`velocity_seed`、`remove_com_velocity`
  - `.ff` 力场参数文件，支持单/多元素格式（`load_force_field()`）  
  健壮性改进：`xyz` 第 0 列整数类型写入 `System::atom_types_`；原子行错误含行号；质量非正数报告具体原子序号。

- 输出模块：已完成基础实现  
  对应文件是 trajectory_writer.hpp 和 trajectory_writer.cpp。  
  当前已实现：
  - **XYZ 轨迹文件**（`*.xyz`）：extended XYZ 格式，每帧含原子数行、注释行（step/time/PE/KE/T），随后每原子一行 `type  x  y  z`
  - **能量日志文件**（`*.log`）：空格分隔表格，列为 `step  time[fs]  PE[eV]  KE[eV]  E_total[eV]  T[K]`
  - `write_frame()` 支持外部传入 `twice_ke`（避免重复遍历速度）或自动内部计算
  - `write_frame_if()` 按步长间隔过滤，`frame_count()` 追踪已写帧数
  - 析构时自动关闭并刷新文件

- 主循环模块：已完成基本骨架，初始力 bug 已修复，邻居表已接入  
  主循环并不在 system.hpp。`System` 现在只是状态容器。  
  真正的主循环模块是 simulation.hpp 和 simulation.cpp，程序装配入口在 gmd_main.cpp。  
  当前这部分已经实现了：
  - 装配 `System`
  - 配置 `ForceProvider`
  - 配置 `Integrator`
  - 配置 `VelocityInitializer`
  - 配置 `NeighborBuilder`（可选）
  - `initialize()`：速度初始化 → 邻居表首次建立 → t=0 初始力计算，确保 Velocity Verlet 第一个半步使用真实力；初始力计算时透传邻居表
  - `step()`：自动检查 `needs_rebuild()` 并在必要时重建邻居表
  - `run()`  
  主循环骨架已完成，配合已实现的 LJ 力场和 cell-list 邻居表，现在能真正驱动物理正确、O(N) 复杂度的原子动力学。

Problems to be solved:

- ~~最关键的问题是初始力没有在进入时间循环前先算一次~~（已修复：`Simulation::initialize()` 末尾现在会执行一次初始力计算）
- ~~截断势没有做 shifted potential 处理~~（已修复：`ClassicalForceProvider` 在构造时预计算 `energy_shift_`，各对势能均减去 `V(r_cut)`）
- ~~没有邻居表，受力计算是 O(N²) 双循环~~（已修复：`VerletNeighborBuilder` 以 cell list 建表，`ClassicalForceProvider` 走 O(N) 邻居表路径）
- ~~没有 thermostat、barostat、不同元素力场~~（已修复：VelocityRescaling/Nosé-Hoover thermostat、Berendsen barostat、多元素 LJ 均已实现并接入积分器）
- ~~输出模块（TrajectoryWriter）尚未实现~~（已修复：XYZ 轨迹 + 能量日志均已实现，支持按间隔输出）
- ~~输入解析脆弱：未存储原子类型、错误信息缺行号~~（已修复：`atom_types_` 现在正确从 xyz 第 0 列读取；错误消息含原子序号）
- 仍缺少：约束（SHAKE/RATTLE）、键角/二面角力场，因此仍非通用 MD 程序。
