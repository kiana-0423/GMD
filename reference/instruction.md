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

- 边界条件模块：已拆分为独立文件，部分完成  
  对应文件是 boundary/periodic_boundary.hpp、boundary/periodic_boundary.cpp、boundary/minimum_image.hpp 和 boundary/minimum_image.cpp。  
  当前已实现：
  - PBC：`wrap_coordinate`、`wrap_position`、`wrap_positions`
  - MIC：`apply_minimum_image_component`、`apply_minimum_image`  
  当前接线状态：
  - `VelocityVerletIntegrator` 已使用 `wrap_position` 做坐标周期回卷
  - `minimum_image` 逻辑已提供接口，但尚未接入 `ClassicalForceProvider` 的真实受力计算

- 力场模块：接口完成，实际物理实现未完成  
  对应文件是 force_provider.hpp、classical_force_provider.hpp 和 classical_force_provider.cpp。  
  当前已经有统一的 `ForceProvider` 抽象，也已经把 `ClassicalForceProvider` 接进了主流程。  
  但 `ClassicalForceProvider::compute` 目前只返回零势能和零力，不是实际的 LJ 或其他真实势函数。  
  所以它是“框架接通了，但物理力场还没实现”。

- 积分模块：基本完成  
  对应文件是 integrator.hpp、velocity_verlet_integrator.hpp 和 velocity_verlet_integrator.cpp。  
  现在已经有：
  - 抽象积分器接口
  - `VelocityVerletIntegrator` 具体实现
  - 与 `ForceProvider` 的实际联动
  - 坐标和速度更新
  - 步进中的周期回卷  
  从工程结构上，这个模块已经不是空壳了，已经能真正驱动 `Simulation::step()`。  
  但因为力场还是零力占位，所以目前积分模块跑起来更像“自由粒子演化 + 周期边界”。

- 输入模块：已完成并扩展初始化相关配置  
  对应文件是 config_loader.hpp 和 config_loader.cpp。  
  目前能读取：
  - `xyz` 文件中的原子数、坐标、质量、盒子尺寸
  - `xyz` 8 列格式中的输入速度 vx/vy/vz
  - `run` 文件中的步数、时间步长、温度
  - `run` 中的初始化控制项：`velocity_init`、`velocity_seed`、`remove_com_velocity`  
  并且包含基础合法性检查。

- 输出模块：基本未完成  
  对应接口 trajectory_writer.hpp 目前还是前向声明，没有实现。  
  现在程序唯一的输出只是 gmd_main.cpp 里最后打印一行启动/运行摘要，不具备：
  - 坐标轨迹输出
  - 能量输出文件
  - 温度/热力学量输出  
  所以输出模块目前应判断为“未完成”。

- 主循环模块：已完成基本骨架，但文件归属应修正  
  主循环并不在 system.hpp。`System` 现在只是状态容器。  
  真正的主循环模块是 simulation.hpp 和 simulation.cpp，程序装配入口在 gmd_main.cpp。  
  当前这部分已经实现了：
  - 装配 `System`
  - 配置 `ForceProvider`
  - 配置 `Integrator`
  - 配置 `VelocityInitializer`
  - `initialize()`
  - `step()`
  - `run()`  
  所以主循环骨架已经完成，而且已经能实际跑若干步。

Problems to be solved:

- 最关键的问题是初始力没有在进入时间循环前先算一次。现在的顺序是先初始化速度，然后直接进入循环做积分，再算力。这会导致第一步的前半步速度更新使用的是全零力，而不是t=0时刻的真实力。严格的Velocity Verlet通常应先做一次初始 force evaluation，再进入 half-kick / drift / force / half-kick。
- 没有邻居表，受力计算是 O(N²) 双循环，体系一大就不实用。这也是它和正式框架中 NeighborBuilder 抽象的明显差别。
- 截断势没有做 shifted potential 或 switched force 处理，截止半径处能量和力都有不连续，长时间能量守恒会比较差。
- 没有 thermostat、barostat、约束、轨迹输出、不同元素力场、键角二面角等，所以它更像“稀有气体 LJ 演示器”，不是通用 MD 程序。
- 输入解析也比较脆弱，比如默认直接访问第二个token，没有完整健壮性检查，精度转换里 getDouble 里还先用了 float 临时变量。
