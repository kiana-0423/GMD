# GMD Validation Suite

本目录用于建立可重复运行的科学验证基线，目标是把两类事情分开：

- `static_*`：小体系、固定坐标、强数值校验，适合进 CI。
- `nve_* / nvt_* / npt_* / diffusion_*`：较长工作流验证，默认作为可选 validation，不放进每次 CI。

统一约定：

- 静态 case 通过 `gmd_validate` 生成 `result.json`。
- 动力学 case 通过 `gmd` 生成 `output.log` 和 `output.xyz`，再由 `analyze.py` 归档为 `summary.json` 与 `series.csv`。
- 容差全部放在每个 case 的 `tolerance.json`。
- `source` 标注 `analytic_reference_*` 表示 independent analytic validation，不依赖 GMD 输出。
- `source` 标注 `provisional_gmd_baseline` 表示 regression baseline：该 case 已具备重复运行与自动比较能力，但不能作为 independent scientific validation。

当前 reference 成熟度：

- `static_lj_cluster`：解析 LJ pair-sum reference。
- `static_special_pairs`：解析 shifted-LJ + Ewald + special-pair scaling reference，并覆盖修改 1-4 scale 的变体。
- `static_coulomb/reference_ewald.json`：解析周期 Ewald reference。
- `static_bonded_reference`：bond / angle / proper dihedral / improper 的 LAMMPS 外部 reference（LAMMPS 22 Jul 2025 - Update 5）。单位、functional form、improper 的 sign convention mapping 与 atom ordering 全部记录在该 case 的 `README.md`；`generate_reference.py` 只在 reference 需要重新生成时手动运行，validation 本身不依赖 LAMMPS。
- `static_coulomb/reference_pme.json`：仍为 provisional GMD PME regression baseline，等待 LAMMPS PPPM 或 OpenMM PME 外部参考。**该 baseline 已于 2026-08-27 重新生成**：修复三处 PME reciprocal force 缺陷后，旧 baseline 记录的 force 缺失约 94%（reciprocal force 恒为零）。重新生成后的 force 与同 case 的解析 Ewald reference 一致到 8.2e-4，即 grid 16^3 / order 4 / alpha 0.3 下应有的 mesh error；旧值与之相差 8.9e-2。energy 未变。
- `pme_external/`：**已完成的 PME 外部参考**。OpenMM 8.6.0 PME 为主参考，LAMMPS 22 Jul 2025 - Update 5 的 PPPM 与 exact Ewald 为第二引擎。Coulomb-only、非立方 18x22x26 A、12 原子、精确电中性、无对称性的 fixture。alpha / grid / cutoff / 边界条件 / 无 exclusion 全部精确对齐；B-spline order 无法与 OpenMM 对齐（其固定为 5，无 API 暴露）；LAMMPS PPPM 用的是 optimised Green's function，本身就是另一种 mesh 近似。三个代码对同一物理常数的取舍不同：GMD 用 CODATA 2022 的 14.3996454686836（`gmd::kCoulombConstant`），LAMMPS 用 14.399645（自身六位小数，低 3.255e-08），OpenMM 折算为 14.399645478（高 6.765e-10）；两个引擎的常数在生成时实测而非引用文档。由于每一项都精确携带一个 k_e 因子，按 k_e^GMD / k_e^engine 线性缩放是精确修正。**在 GMD 的常数被修正之前**该因子还要补偿 GMD 自身 3.16e-06 的截断，那个偏差会超过 grid 128 的 mesh error 并被误读为收敛下限；现在只剩引擎自身的舍入。三个引擎在 16/32/64/128 网格上单调收敛到同一个 exact Ewald 极限。reference settings（grid 64^3，GMD order 6 vs OpenMM order 5）下 energy 差 1.13e-06 eV、force 最大分量差 1.53e-06 eV/A；容差取 |GMD-exact| + |OpenMM-exact| 三角不等式界的两倍，而非把观测值向上取整。**virial 张量全部九个分量对 LAMMPS 验证**（exact Ewald 1.6e-07 eV，PPPM 8.6e-07 eV）；OpenMM 完全不暴露 virial。重新生成 reference 需要两个外部引擎，CI 比较不需要。
- 长时间 NVE/NVT/NPT/diffusion cases：仍为 provisional workflow/regression baselines。

约束体系初速度的 DOF 修正与 Berendsen serial/MPI 强制比较（2026-08-30）：

- **缺陷**：约束体系请求 300 K，实际从 **514 K** 起步（单个刚性水分子）。
  `VelocityInitializer` 按 `2K = (3N−3)·k_B·T` 缩放，随后 integrator 才把速度投影到约束
  切空间——把刚刚放进去的动能又拿掉一部分——而温度报告用的是权威的 `3N − rank − 3`。
  两个自由度计数、两个阶段，彼此都不知道对方；initializer 既没有 constraint solver
  也没有 rank，本来也无从得知。
- 两个误差**不会相消**。投影移除的是约束模式上的能量，平均约占 `rank/(3N−3)`，与两个
  计数之差几乎相同——于是误差的**均值近零而涨落不然**，而每次运行只抽一次样：
  1 分子 +71.3%、3 分子 +22.9%、10 分子 +4.8%、40 分子 +1.6%。
- **修正**：`Simulation::initialize()` 改为**先** integrator、**后** velocity initializer。
  initializer 需要的一切只在 integrator 跑完后才存在：位置被 SHAKE 投影到约束流形上，
  `require_independent()` 在**投影后的几何**上接受该约束集（rank 是动力学真正起步的那个
  构型的性质，而非输入构型的），权威的 `3N − rank − 3` 随之确定。
- **顺序**：采样 → 移除全局质心速度 → RATTLE 投影 → 用权威 DOF 缩放。**一遍即可**，
  原因是具体的：质心移除保持切空间性质（距离约束的 Jacobian 行为 `(+r_ij, −r_ij)`，
  故整体平移落在 `J` 的零空间中）；投影保持零动量（修正量为 `+λr_ij/m_i` 与
  `−λr_ij/m_j`，动量变化相消）；标量缩放两者都保持。同一个零空间事实也说明平动模式与
  约束模式不重叠，因此同时减去 3 和 rank **不会重复扣除**。
  实现仍写成有上限的循环，**验证**两个全局标量性质并在不满足时明确报错，而不是默默接受
  一个不自洽的速度场。投影本身调用 integrator 自己的 `apply_velocity_constraints()`，
  约束数学仍留在 constraint solver 内。
- **保证**：切空间残差 `|r·v| < 1e-11`；质心动量为归约舍入量级；`2K = dof·k_B·T` 达到
  3.3e-16；DOF 与温度报告/thermostat/压强所用的是同一个；np=1/2/4、跨 rank 约束、
  空 rank 下按 tag 一致；存储顺序置换后按 tag 差异 8.3e-17。
- **baseline 不变**，且是验证过的：validation 中没有任何 case 使用约束，十个 case 全部
  逐位复现（abs_error 恰为 0）。约束体系本身的轨迹确实改变——那正是目的。
- **restart 不受影响**：restart 根本不安装 velocity initializer。

- **Berendsen serial/MPI 比较由"仅报告"改为强制断言**。原先只报告不断言，是因为写它时
  `velocity_init random` 确实依赖 rank。抽样改为按全局 tag 之后该理由已不存在：
  np=1/2/4 的能量日志**逐列逐帧完全相同**；最终每原子状态在 np=1 逐位相同，
  np=2 / np=4 分别为 1.7e-14 / 2.1e-14（容差 1e-12，约五十倍）。
  结构性检查（帧数、有限性）先于任何差异度量执行，因此 NaN 或缺帧会被如实报出而不是
  混进误差里；诊断信息给出日志差异的列名与步数、状态差异的 atom tag 与分量。

随机初速度的 rank 无关性修正（2026-08-30）：

- **缺陷**：`velocity_init random` 依赖 rank。`VelocityInitializer` 用单个
  `std::mt19937` 按**本地存储顺序**推进，因此一个物理原子拿到哪几个随机数由它在数组中的
  位置决定。区域分解下每个 rank 都从相同的种子状态出发，并把最初几个抽样给了**自己的**
  第一个本地原子——在不同 rank、不同 rank 数下那是不同的物理原子。随后对目标温度的全局
  rescale 又把差异掩盖：总动能无论如何都被拉到目标，所以 step 0 的温度与势能完全相同，
  而底层速度场不同。MPI 运行因此从不复现同种子的串行轨迹。
- **修正**：抽样改为身份的纯函数，`value = f(seed, stream, 全局 atom tag, 分量)`，
  原子之间不携带任何状态，也不依赖遍历顺序。实现见
  `include/gmd/core/keyed_random.hpp`：SplitMix64（Steele/Lea/Flood, OOPSLA 2014，
  即 `java.util.SplittableRandom` 所用的 finalizer），两轮 multiply-xorshift，
  常数有据可查，纯整数运算。
  - 均匀数取 `((bits >> 12) + 0.5) / 2^52`，严格落在开区间 `(0,1)`。**丢弃 12 位而非 11 位**：
    用 53 位时最大值 `(2^53−1)+0.5` 不可表示（该量级间距为 1.0），按 half-to-even 进位到
    `2^53`，商恰为 `1.0`，于是 `log(u)` 恰为 0。测试专门断言 53 位版本确实会进位到 1.0。
  - 无取模，故无 modulo bias。
  - x/y/z **各自**用自己的 key 抽取自己的一对均匀数，而不是共用一对取 cos/sin——后者会让
    两个输出严格相关。
  - `RandomStream` 分离流，将来新增随机功能不会挪动这里的数值。
- **tag 校验**：tag 成为承载语义的量，故使用前先校验。负 tag 与重复 tag 均被拒绝；重复可能
  **跨 rank**（各 rank 内部唯一、两个 rank 共用一个），只有全局检查能发现。所有 rank 抛出
  同一条诊断：只有部分 rank 抛出会让其余 rank 独自阻塞在下一个 collective，把错误输入变成
  死锁而非报错。**不会**静默回退到数组下标——那正是缺陷本身。
- **ghost 原子不参与初始化**。ghost 携带宿主的 tag，按 tag 抽样反而会得到正确的值；真正错误
  的是归约：该原子的动量与动能会被计入质心和与 rescale 两次。采样与三处归约、以及 tag 校验
  都只遍历 owned 原子，因此 halo 中重复的 tag 不会被误判为重复。
- **保证**：按 tag 比较，np=1/2/4 以及本地存储顺序反转下均一致，误差 **2.8e-17**；抽样本身
  **逐位相同**。残差来自归约顺序（rank 内按存储顺序累加、跨 rank 由 `MPI_Allreduce` 合并，
  两者在不同 rank 数下都不同），经由唯一的质心平移与唯一的 scale 因子传播到每个原子。
  修正前等价性 fixture 的 24 个原子**全部**不同，最差分量 4.7e-01——改善十六个数量级。
  端到端：50 步轨迹在 serial/np=1/2/4 下日志**完全相同**，最终状态在 np=1 逐位相同、
  np≥2 内 8.9e-16（约 28 ulp）。
- **串行结果改变**，且**刻意不保留**旧的串行随机流：保留它就等于保留"遍历顺序即随机身份"，
  那就是缺陷。**分布未变**：目标温度、`sqrt(k_B T/m)` 宽度、质心移除、3N−3 rescale 全部不变，
  且整个速度场已对照该规范的独立重实现验证到 2.8e-17。
- **restart 不受影响**：restart 根本不安装 velocity initializer，checkpoint 中的速度被恢复
  而非重新采样，因此无需序列化任何 RNG 状态。串行运行跨 checkpoint 拆分后与连续运行逐位相同。
- **baseline**：五个 dynamics reference 全部重新生成。`diffusion_lj_fluid` 是全仓库唯一被推出
  原容差的指标（1.796 → 3.273 Å²/ps）；该容差是**固定种子下的可复现性**边界，而非 D 的测量
  精度：同一 fixture 换五个种子给出 3.273、2.664、2.260、2.224、2.205，离散度约 48%，新值只是
  其中一次普通抽样。**所有 static energy / force / virial reference 均未改变**——它们都不赋速度。

内部时间单位审计与 Berendsen 轨迹验证（2026-08-30）：

- **内部时间单位不是自由选择**。GMD 以 `v += (F/m)·dt`、`r += v·dt` 积分，`F` 为 eV/Å、
  `m` 为 amu，要求 `F/m` 是本单位制下的加速度即唯一确定
  `T = Å·√(amu/eV) = √(m_u/e)·1e5 fs = 10.1805057178711931...` fs。
  与 `k_B`、压强换算不同，**这个常数带真实不确定度**：`m_u` 的 3.1e-10 相对不确定度经开方
  减半为 1.57e-10。速度单位为其倒数 `√(eV/amu) = 0.0982269474...` Å/fs。
- **三个缺陷**：
  1. `kInternalTimeUnitsPerFs = 1.018051e+1` 为七位有效数字取整，比精确值高
     **4.206204e-07**，约为 CODATA 不确定度的 2700 倍；且其**名称与用法互为倒数**
     （它被用作 fs 时间步的除数）。现为 `gmd::kFemtosecondsPerInternalTime` 与
     `gmd::kInternalTimePerFemtosecond`，由 `static_assert` 保证互为精确倒数。
  2. **Nosé–Hoover 的 `tau` 从未换算**：`Q = dof·k_B·T·tau²` 用的是 fs 数值，而 `ξ`
     以内部单位的 `dt` 积分。**请求 tau = 100 fs 实际弛豫时间为 1018.05 fs**（恰好差 `T`），
     `Q` 大了 `T² = 103.6427` 倍。
  3. **Berendsen 的 `tau_P` 同样未换算**（`mu³ = 1 − beta·(dt/tau)·ΔP`）：
     **请求 2000 fs 实际耦合为 20361 fs**。
- 两处 tau 缺陷在各自文件内**不可见**：每个量都自洽，且两个对象都直接由 `RunConfig` 构造。
  `RunConfig` 现按 `time_step_fs`/`time_step` 的既有模式，新增
  `thermostat_tau_fs`/`thermostat_tau` 与 `barostat_tau_fs`/`barostat_tau`。
  **未改动任何 thermostat/barostat 方程。**
- **原本就正确的部分**：日志 `time[fs]` 列（由调用方按 `step × time_step_fs` 计算）、
  checkpoint 的 `time_fs`、以及扩散系数 Å²/ps 背后的 fs→ps 换算。
  `Simulation` 传给 force provider 的 `force_time` 用的是内部单位，那是另一个量，
  不进入任何日志或 checkpoint。
- **baseline 影响**：`nvt`（T mean 123.81→120.41 K、stddev 15.66→**22.00** K，
  stddev 是全仓库唯一被此修正推出原容差的指标）、`npt`（T mean 126.44→119.92 K、
  stddev 17.04→21.17 K、pressure mean 45.99→35.56 bar，均在容差内）、
  `diffusion`（1.79598→1.79603 Å²/ps，仅 2.8e-05，因该 run 未配置 thermostat）已重新生成；
  `nve` **逐位不变**（无 thermostat/barostat，且其 drift 指标分辨不出 4.2e-07）。
  **所有 static energy / force / virial baseline 与时间无关，未重新生成也未变化。**
  耦合变强使温度涨落上升、同时均值更贴近目标，是该修正的预期特征而非漂移。

- **新增 `berendsen_npt_lj`**：仓库中唯一在轨迹层面覆盖 Berendsen barostat 的 case
  （`npt_lj_fluid` 用的是 MC barostat）。**属 dynamics regression，不是外部科学参考**。
  它是符号敏感的成对运行，且两个已修正的缺陷都由构造捕获，阈值由**实测缺陷代码**确定：
  恢复 bar/eV/Å³ 混比会让两个 run 都压缩（1.0129 之外的方向断言直接失败）；
  恢复未换算的 tau 会把响应从最小 0.0490 压到最大 0.0129，`min_volume_response = 0.03`
  正落在两者之间。另含体积/能量/温度/压强有限性、体积不越界、以及 checkpoint 前后
  最终体积**逐位相同**的连续性检查。该 case 已注册为 CTest 测试。

- **已知限制**：`velocity_init random` 的 MPI 运行**无法复现同种子的串行轨迹**。
  `VelocityInitializer` 用单个 generator 按**本地**原子下标顺序抽样，因此每个 rank 都把
  generator 的前三个抽样给了自己的第一个原子——在区域分解下那是与串行不同的物理原子；
  随后对目标温度的 rescale 又把差异掩盖（step 0 的温度与势能完全相同）。
  这不是 barostat 或受力的问题：同一 fixture 以纯 NVE 运行同样发散，而以
  `velocity 0.0` 运行则 50 步**逐位相同**。修复需改为按全局 tag 抽样，会改变所有随机初速度
  轨迹，故此处仅记录。因此 `berendsen_npt_lj` 只报告而不断言 serial/MPI 差异；
  可断言的 rank 无关性在 `tests/mpi_berendsen_barostat.cpp` 中，
  它在每个 rank 上按全局 tag 构造相同的速度场后再比较耦合因子（np=1/2/4，np=4 含空 rank）。

压强单位换算修正（2026-08-30）：

- bar ⇄ eV/Å³ 的换算因子曾以**两个独立字面量**存在：`src/io/trajectory_writer.cpp`
  与 `include/gmd/integrator/mc_barostat.hpp` 各自的 `6.2415091e-7`，两者都比精确值
  高 **4.091837e-09**。现统一为 `gmd::kEVPerAngstromCubedToBar` 与
  `gmd::kBarToEVPerAngstromCubed`。与 k_e、k_B 不同，这不是实测量而是**纯单位恒等式**，
  四个成分（bar = 100000 Pa、Pa = 1 J/m³、Å = 1e-10 m、eV = 1.602176634e-19 J）全部
  按定义精确，故 1 bar = 500/801088317 eV/Å³，没有可供取整的不确定度。逆方向是**有限
  小数**：1 eV/Å³ = 1602176.634 bar（精确），因此以该方向为字面量，正方向由其取倒数
  导出——这样两个方向各自都是对应精确有理数的最近 double，且在 double 运算下互为精确
  倒数（由三条 `static_assert` 保证）。
- **Berendsen barostat 曾把 bar 与 eV/Å³ 直接相减**。它按 `(2K + tr W) / 3V` 计算瞬时
  压强（eV/Å³），却直接减去以 bar 传入的目标压强，两侧都未换算。于是请求 1 bar 实际等于
  请求 1 eV/Å³，即 1602176.634 bar。这不只是尺度错误：该差值决定耦合的**符号**，因此对
  任何常规目标压强，盒子的缩放方向都与真实压强无关。`beta` 本就是 bar⁻¹（默认 4.5e-5 为
  液态水值），故比较改在 bar 下进行，换算的是瞬时压强。算法、符号约定与 virial 处理均未改动。
- **力、能量与 virial 完全不受影响**：压强在此只是被报告和被控制的量，从不作为力的输入。
  所有 static baseline 按构造不变。
- **只有 `npt_lj_fluid` 一个 baseline 变化，且轨迹本身没有变**。日志中除压强外的每一列
  （step / time / PE / KE / E_total / T / V）在全部 201 帧上**逐位相同**，说明 MC barostat
  接受了完全相同的体积移动序列；该 case 的 temperature 指标因此也逐位不变。压强指标的相对
  变化（+4.33e-09 / +4.75e-09）**不等于**常数本身的 4.0918e-09 比值，原因是日志精度而非
  物理：压强约 46 bar 而以六位小数打印，一个打印单位是 1e-6 bar，换算只移动约 1.9e-7 bar，
  指标是被量化后取的平均。逐帧看效应完全有界：201 帧中 62 帧发生变化，每一帧都恰好变化
  一个末位打印单位，51 帧上移（压强全为正）、11 帧下移（压强全为负），无例外；由计数预测的
  净偏移 (51−11)×1e-6/201 = 1.9900497512e-07 bar 与实测均值偏移 1.9900497250e-07 bar
  相符到 2.6e-15。
- `nve` / `nvt` / `diffusion` 重新运行后**逐位不变**（每个指标 abs_error 恰为 0.0）：三者
  都没有 barostat，换算只到达它们不测量的 P[bar] 列。
- **restart**：checkpoint 存的是内部单位 eV/Å³，换算**不被序列化**，因此旧 checkpoint 会以
  修正后的换算恢复；存储的数值不变，由它报告出的 bar 值移动 4.09e-09。既有日志与轨迹中的
  `P[bar]` 列由旧因子写出，与新值在该量级以下不可比。
- 容差未变，且是刻意如此：其推导不依赖该换算，且两个压强指标只移动约 2e-07 bar，容差为
  2000 bar；重新生成只是因为原值本可逐位复现而现在不再如此。

温度相关 baseline 与 Boltzmann 常数修正（2026-08-30）：

- 生产代码曾同时存在**四个** Boltzmann 常数：`src/system/initializer.cpp` 的
  `8.617343e-5`、`thermostat.hpp` 与 `mc_barostat.hpp` 的 `8.617333262e-5`、
  以及 `examples/ethane_demo` 的 `8.617333e-5`。现统一为
  `gmd::kBoltzmannConstantEVPerKelvin = 8.617333262145177e-5`（由 2019 SI 精确定义
  推导：k_B = 1.380649e-23 J/K 与 e = 1.602176634e-19 C 均为精确值，故
  k_B[eV/K] = 1380649/16021766340 是**精确有理数**，NIST 记作
  "8.617 333 262... x 10^-5 eV/K, exact"）。
- **可观测的缺陷**：速度初始化与温度报告用了不同常数，初始化到 300 K 的体系
  报告 300.000339014 K；修正后为 300 K（2.2e-16）。
- **仅影响随机初速度的动力学 case**。`nvt` / `npt` / `diffusion` 三个 baseline
  已重新生成（相对变化 1.1e-07 ~ 2.6e-04，全部仍在原容差内）；`nve` 的
  energy drift **逐位不变**并已在文件中给出证明——该指标是两个 `%.6f` 打印能量之差，
  而修正只移动约 3e-08 eV，低于最后一位打印精度两个数量级。
- **所有 static case（LJ / Coulomb / special pairs / bonded / pme_external）完全不受影响**：
  没有任何 Coulomb、LJ 或 bonded 量依赖 k_B。

当前 release-facing 状态：

- LJ：已通过解析 reference 验证。
- Ewald：已通过解析 periodic Ewald reference 验证。
- replicated PME：**已有独立外部参考**（OpenMM PME + LAMMPS PPPM，见 `pme_external/`）。`static_coulomb/reference_pme.json` 仍是 GMD self-baseline，作为 regression 保留，不再是 PME 的唯一证据。
- `pme_mode distributed`：仅为 interface/prototype，数值 backend 仍是 replicated PME。
- SHAKE/RATTLE：串行与一个 cross-rank MPI correctness-first global-gather 路径已有测试；不是 scalable distributed constraint solver。
- 约束动力学：已改为标准 SHAKE/RATTLE splitting（修正了投影方向与缺失的半步速度冲量），含约束轨迹改变；约束 virial 作为 `t+dt` 端点量进入 reported pressure，验证参考为刚性转子的向心力解析解，不依赖实现公式；详见下节。**没有外部引擎参考。**
- checkpoint/restart：真实 `gmd` CLI restart-continuity 测试覆盖 serial、MPI np=2、MPI np=4。
- virial：**每一个 virial source 都已逐分量（全部九个分量）对独立 reference 验证**，见 README.md 的 *Virial validation coverage* 表。`Box` 只能表示 orthorhombic cell、无法施加 shear strain，因此 `tests/virial_finite_difference_tests.cpp` 只验证 trace 与三个对角分量；off-diagonal 由三类不需要引擎 shear 的 reference 覆盖：解析 pair identity、独立 bonded force moment，以及一个 test-only、写在一般 3x3 cell 上的 reciprocal energy（**可以** shear），对其做 general strain 数值微分即得全部九个分量。rotation covariance 是必要条件而非 shear derivative 的替代品。**静电 virial 现已有外部引擎参考**：`pme_external` 把 Ewald 与 PME 的全部六个独立分量（张量对称性被显式测量而非假定，故覆盖全部九个）对 LAMMPS 比较，分别到 1.6e-07 eV 与 8.6e-07 eV。**其余 virial source 仍无外部参考**：LJ、bonded、special-pair、constraint 都只有解析/数值独立参考——LAMMPS 只导出合并后的压强张量，单独拆出某一项需要未经证明的等价假设；OpenMM 根本不暴露 virial。

## 约束动力学与约束 virial 的验证范围

**积分器已改为标准的 constrained velocity-Verlet（SHAKE/RATTLE）splitting**，修正了两处缺陷，
因此**所有含约束的轨迹都会改变**：

- SHAKE 之前沿*漂移后*的键矢量 `r_c(t+dt)` 投影，而 `sigma_c = |r_c|^2 - d_c^2` 的梯度应取在步首
  构型，即沿 `r_c(t)`。两者都落在约束流形上，但落点不同，只有前者才对应约束运动方程。沿漂移键
  投影等价于对键长做缩放，非 symplectic，会持续泄漏能量：自由刚性转子在 `omega*dt = 0.04` 下
  4000 步损失约 **86%** 动能，修正后守恒到 1e-13。
- SHAKE 缺少与位置修正配套的半步速度冲量 `v += dr/dt`。

一步之内约束力出现**两次**（w = 1/m，s_ic 为 ±1）：

```
(1) v_i += (dt/2m_i) F_i(t)                 (2) v_i += w_i sum_c s_ic Gamma_c r_c(t)      <- SHAKE
(3) r_i += dt v_i                           (4) v_i += (dt/2m_i) F_i(t+dt)
(5) v_i += w_i sum_c s_ic Lambda_c r_c(t+dt)                                              <- RATTLE
```

(2) 是 `t` 时刻的约束力，与 `F(t)` 配对；(5) 是 `t+dt` 时刻的，与 `F(t+dt)` 配对。provider virial
在 `r(t+dt)` 上求值，故取 (5)。把 (5) 与 (4) 逐项对应：

```
(dt/2m_i) G_i(t+dt) = w_i Lambda_c r_c(t+dt)
W_constraint(t+dt)  = sum_c (2 Lambda_c / dt) (r_c (x) r_c)
```

**系数 2/dt 是精确的**，且这是**端点量而非步平均量**——报告的 virial 两部分处在同一时间层，
不存在 `O(dt)` 的时间层混用。系数由外部物理结果反推确认而非假定：刚性转子的端点约束力必须等于
向心力 `mu omega^2 d`，用 `1/dt` 会正好小一倍。

**中心性与对称性**由构造保证：RATTLE 不移动原子，`r_c(t+dt)` 在整个求解过程中固定，逐约束累加
**标量**乘子得到的成对力严格沿 `r_c` 且严格反对称，耦合约束（刚性三角形）同样如此，没有任何投影，
也不丢弃任何残差。

### 报告压强与当前几何 virial 是两回事

- `System::last_virial()` 是**当前几何**的张量。barostat 缩放盒子后重新计算的 provider virial
  没有对应的约束乘子，此时被标记为 **invalid**，而不是冒充完整的压强 virial。
- `System::step_pressure()` 是**刚结束那一步**的压强，在 barostat 运行之前捕获，与 barostat
  自身消费的数值逐位相同；之后的任何一次力求值都不会改动它。报告的压强取自这里。

日志一行内的 `PE / KE / E_total / T / P / V` 现在全部取自**同一个状态**（完成步，在 barostat
缩放之前一并捕获），因此互相自洽：`P = (2K + tr W)/3V` 在这几个数之间严格成立。唯一不取自该状态的是
`.xyz` 里的坐标，它是当前构型（barostat 运行时即缩放后、下一步的起点），故在 barostat 运行时该行的
`V` 相对 `.xyz` 坐标滞后一次缩放。

**MPI 输出路径**另建一个只装全局坐标的 System，因此写帧前必须把所有非原子字段从分布式 System
拷过去（`System::copy_frame_state_from()`）：盒子、势能、完成步记录、SHAKE/RATTLE 诊断、
以及初始帧回退所用的 virial 状态。这些**都不做 reduce**——它们要么是复制量（盒子、约束 virial），
要么已是全局量（provider virial、由全局归约量构成的完成步记录）；该"每个 rank 都相同"的假设由
`tests/mpi_constraint_virial.cpp` 显式断言，并由 `tests/constrained_pressure_reporting.py`
在 serial / np=2 / np=4 之间逐列比较整份日志来端到端验证。

**有效性是显式的，绝不用哨兵值。** 数值 0 是一个正常的压强，因此日志新增 `P_valid` 列，
无法给出完整压强时写 `nan`。目前唯一会出现的情形是**含约束运行的初始帧**：还没有任何一步完成，
初始力求值也没有 RATTLE 乘子可配对。无约束运行每一帧（含初始帧）都给出有效压强。

### 验证参考

`tests/constraint_virial_tests.cpp`、`tests/mpi_constraint_virial.cpp`、
`tests/constrained_pressure_reporting.py`：

- **刚性转子的向心力**，与实现公式无关的独立参考。无外力的刚性转子完全由约束力维系，故
  `G_i = -m_i omega^2 (r_i - R_com)`，`tr W = -I omega^2 = -2K`。据此检查恢复出的**端点**力的
  大小、符号、迹，以及压强恒等式 `2K + tr W = 0`；两者的误差都随 `dt` 二阶收敛（实测：dimer 在
  `dt = 0.4 → 0.05` 上相对误差 `+4.0e-4 → +6.25e-6`，每次减半降为四分之一）。
- **两个约束半冲量分别测量**并证明确实不同：各自沿自己时间层的键矢量，方向相差恰好为该步转过的
  角度 `omega*dt`；virial 取的是 RATTLE 那一个。这是 "端点" 一词有意义的直接证据。
- **约束 NVE 能量**：自由刚性转子 4000 步相对漂移 ~1e-13（各 dt）；真实力场的四原子链在固定物理
  时长上漂移随 `dt^2` 下降，且不高于同一体系去掉约束后的对照值。
- 所有 fixture 的初始状态都**同时**位于约束流形上：`|r_ij| = d` 且 `r_ij·v_ij = 0`；步后同时检查
  位置约束与切空间速度约束。
- 耦合约束：三个共享原子的距离约束（刚性三角形），检查张量严格对称、每个成对力严格沿键、
  逐约束可加性，以及同样的 `2K + tr W → 0`。
- 周期边界不变性、整体平移不变性（相对容差，因为张量由 O(30 Å) 的坐标差构成）。
- MPI rank 不变性：同一 fixture 在 1/2/4 rank 下逐分量完全一致；约束 virial 已是全局量，
  **不得再做 allreduce**。
- **约束 Berendsen NPT（serial 与 MPI np=2 / np=4，均走真实 `gmd` 可执行文件）**：连续数十步都发生
  盒子缩放，压强列不得退化为一列 0，`P_valid` 必须全为 1；MPI 与 serial 的整份日志逐列比较，差异不得
  超过日志本身的打印精度（`%.6f`，即最后一位 1 个单位）；
  并且第 1 步的压强必须与"除去 barostat 外完全相同"的运行**逐位相等**——barostat 在第 1 步末尾才
  首次动作，此时完成步压强已经捕获。
- restart 连续性：checkpoint（格式版本 2）显式携带约束 virial 及其**时间层名称**、当前几何
  provider virial、以及完成步压强，`tests/nose_hoover_restart_continuity.py` 比较重启前后逐步的
  pressure。

明确**不**声称：

- **没有外部引擎参考。** LAMMPS 不以可证明等价的定义单独导出约束 virial，且该量是动力学量，
  无法从静态构型比较。此处只声称解析参考。
- **约束 virial 的全部九个分量**已验证：倾斜刚性 dimer 对 `W_ab = -mu omega^2 d^2 u_a u_b`、
  倾斜刚性三角形对二阶矩张量 `W = -omega^2 sum_i m_i a_i (x) a_i`（均由刚体力学独立写出，不调用被测
  累加函数），逐分量二阶收敛；另有旋转协变性 `W -> R W R^T` 检查。仍未验证的是**力场 provider 自身**
  的 off-diagonal virial——盒子只存三个边长，无法施加剪切形变做有限差分。
- 约束独立性现已**验证并强制**：按连通分量计算质量加权约束 Jacobian `J_M = J M^(-1/2)` 的秩
  （单边 Jacobi SVD，秩容差 `max(rows, cols) * eps * sigma_max`），秩小于约束数即在动力学开始前
  报错拒绝。仅在**初始构型**验证：秩与构型有关，轨迹中途退化不会被再次检查。
- 单步内键矢量转过 ~90° 时 SHAKE 的参考梯度线性化无解，此时报硬错误并要求减小步长。

## 生产路径上的约束动力学验证 case（2026-09-16）

在此之前，约束体系的覆盖全部是 unit test 与 `tests/constrained_pressure_reporting.py`；
`validation/` 下**没有任何 case 使用约束**。现补入两个走真实 `gmd` CLI 的 validation case。

### 为什么残差必须从 checkpoint 重算，而不能读日志

`output.log` 的残差列是 `%.6f`：**1e-13 与 1e-7 都打印成 `0.000000`**。想验证约束满足到
1e-10，从日志里根本做不到。checkpoint 则是 17 位有效数字，且携带 tag、质量、坐标、速度、
约束 virial 及其时间层、provider virial、完成步压强。因此两个 case 都在
`validation/constrained_common.py` 里**重新计算**

```
位置残差  ||r_ij| - d_ij|            [A]
速度残差  |r_ij . v_ij| / |r_ij|     [A / internal time unit]
```

这是对生产路径**实际产出**的独立测量，而不是把 solver 自己的 converged 标志读回来。

**单位提醒**：GMD 的速度以内部时间单位存储（`A sqrt(amu/eV)` = 10.180505717871194 fs），
不是 fs。`output.log` 表头把该列标成 `rattle_error[A/fs]`，**这个标注差了 10.18 倍**；
数值本身是内部单位的。

### `constrained_nve_water` — 自由刚性水分子

单个刚性水分子（2 个 O-H 约束 + 1 个等价于固定 H-O-H 角的 H-H 约束），质量不等、一般性旋转、
每个距离初始偏离目标 **+0.25%**、初速度带**非切向**分量——后两者是为了让初始 SHAKE 与 RATTLE
投影真的有活干，否则"投影后距离等于目标"是同义反复。

bonded 项力常数**全为 0**，只用于产生分子内 non-bonded exclusion；分子在盒子里独处，因此
**完全无受力**。这是刻意的：自由刚体的运动解析已知，所以这个 case 的多数检查对的是力学而非
上一次的 GMD 输出——能量、线动量、对质心的角动量严格守恒，质心走直线，以及最锋利的一条：

> 约束 virial 恰好抵消**转动**动能，于是报告压强退化为质心的理想气体压强 `M |v_com|^2 / 3V`。

实测 `2K_rot + tr W = 3.6e-08`（该量本身为 7.7e-03），压强与解析值相对差 **2.4e-07**。
该残差随 `dt` 二阶下降（dt=1.0/0.5/0.25 → 1.2e-07 / 3.6e-08 / 1.8e-08），dt=0.125 起触及
~2e-08 的地板。这一条同时钉住了约束 virial 的**量级（2/dt 系数）、符号与时间层**。

### `constrained_nvt_cluster` — 四个相互作用的刚性分子 + Nose-Hoover

12 原子、12 约束、约束图 4 个连通分量。与上一个 case 的关键区别：分子之间**真的有 LJ 相互作用**，
因此 provider virial 非零，完成步 virial 是**两个独立产生的张量在同一时间层上的和**。
"总 virial == provider + constraint" 这一条正是**约束项取错时间层**会失败的检查，而对着一个
全零张量根本无法做出。

权威自由度 `3N - constraints - 3 = 36 - 12 - 3 = **21**`（无约束应为 33）。这个数是承重的：
Nose-Hoover 的 `Q = dof kB T tau^2` 由它算出，用错会热浴到错误的动能。case 直接断言该值，
并断言初始温度**精确**落在 300 K——只有"先投影后缩放"的顺序与约束自由度两者都对时才会如此。

MPI 布局是特意排的：分子 C 跨 `x = Lx/2`、分子 D 跨 `y = Ly/2`，np=2 与 np=4 下都有约束跨 rank；
`x > Lx/2, y > Ly/2` 象限**留空**，np=4 时有一个 rank 不拥有任何原子。

**MPI 容差是推导出来的，不是猜的**：domain decomposition 改变了力求和顺序，400 步后 MPI 与
serial 的约束 virial 差 `1.03e-11`（相对最大分量 5.8e-10）；而这条容差要抓的失效——张量被按
rank 数重复 reduce——会让 `zz` 偏 `9.6e-03`，**高九个数量级**。np=2 与 np=4 的差异**完全相同**，
这本身就说明它是求和重排而非 rank 缩放。

**明确不声称**：400 步确定性轨迹是 **dynamics / regression 覆盖，不是正则系综的统计验证**。
温度边界是动力学边界，不是涨落定理的预言。

## 严格容差下的 SHAKE/RATTLE 收敛审计（`tests/constraint_convergence_audit.py`）

针对"刚性三角形偏离目标几何约 0.3% 时在 `constraint_tolerance 1e-13` 下耗尽迭代"的报告。

**报告的归因不对**：0.3% 的 slack 是偶然的。在 slack / 步长 / 容差 / 迭代上限 / 质量比 / 取向
上做扫描，失败在 slack 上**均匀分布**，而**完全集中在几何上**。形状良好的三角形在 0.3% slack 下
收敛到 1e-13 毫无困难。

**solver 的数值没有缺陷**，由两个**不使用生产 solver**的独立参考证实：

- RATTLE 在固定几何下对乘子是**线性**的，因此这里用**精确有理数**（`fractions.Fraction`）
  直接消元求闭式解。生产迭代与之一致，且**与精确解的距离按 1/h 标度**（h 为三角形的高）：
  `err*h` 在 h 变化 25 倍的范围内恒定在 2.5 倍以内。这正是接近秩亏的 Jacobian 的条件数放大，
  是**正确**的迭代在病态系统上应有的表现。
- SHAKE 是二次的，用 **60 位十进制** Newton 迭代求解，其自身残差地板远低于 float64。

两种机制会让严格容差不可达，**且补救方向相反**：

| | 成因 | 表现 | 补救 |
|---|---|---|---|
| **类别 2** | 迭代上限不足 | 三角形趋于共线时所需 sweep 数约按 `1/h^2` 增长；收敛仍在继续 | 提高 `constraint_max_iterations`（扫描中 500 次失败的 case 在 4000 次全部收敛） |
| **类别 1** | 浮点地板 | `|r_i - r_j|` 是**坐标之差**，分辨率约为 `eps * max|坐标|` 而非 `eps * 键长` | 放松 `constraint_tolerance`，或让坐标靠近原点 |

类别 1 的实测：**同一个** fixture 在距原点 10 A 处收敛到 1e-13，在 20000 A 处**无论多少次迭代
都不可达**。地板随坐标量级线性增长（坐标 ~1e3 A 时约 8.5e-14，~1e4 A 时约 1.1e-12）。

类别 3/4/5/6 均被排除：两个 solver 各自的残差范数在量纲上自洽；迭代确实到达精确解；精确有理消元
从未遇到奇异矩阵，故没有 fixture 是不可行或秩亏的；而在下述修正之前，失败**没有**被正确诊断。

**关于单调性**：本审计早期曾把收敛描述为单调，这是错的，测试也不这样断言。Gauss-Seidel 在
max 范数下不是下降法——扫描一个约束会扰动其余两个，"当前最差"在 sweep 之间来回切换，因此最差残差
单个 sweep 内最多可上升约 1.7 倍，而求解完全正常。测试断言的是**不发散**：单 sweep 增长有界，
且包络下降十个数量级。

### 唯一的生产改动：失败诊断

审计没有发现数值缺陷，因此**没有改动任何收敛判据、容差或残差界**，失败仍然是硬失败。改的是这句话：

```
Error: SHAKE failed to converge within 500 iterations; max bond error = 0.000000
```

`std::to_string` 是**六位定点**，于是严格容差失败能携带的任何残差（1e-13、1e-11、1e-8）
**全部打印成 `0.000000`**——报告"误差恰好为零"的同时拒绝收敛。本文件其实早就为秩分析诊断过
同一个问题（`format_number()` 就在几百行之上），只是 SHAKE / RATTLE 的失败路径一直没切过去；
tolerance-equivalent duplicate 诊断也把容差本身打成了 `0.000000`。

比不可读更糟的是：它无法区分上表那两种**补救方向相反**的情形。现在失败会记录最优残差与它停止
改善的位置，并说明属于哪一类；停滞的 SHAKE 还会给出坐标量级与它蕴含的分辨率。RATTLE **故意不**
给这个地板估计——它的残差 `|r.v|/|r|` 是**速度**，地板由速度尺度而非离原点的距离决定，在那里给出
SHAKE 的数字会是一个看着合理却不适用的值。

诊断所需的一切都取自已经复制到各 rank 的原子数组，因此每个 rank 构造出相同的消息，失败保持 collective。

### 负控制（negative controls）

在一次性 worktree 里把生产代码按下表逐条弄坏、重建、跑 validation，**要求 validation 失败**。
控制通过（validation 仍然绿）说明该检查其实没在测它声称的东西。所有变异跑完后全部还原，
并比对整棵树的 checksum；**变异本身从不提交**。

| 变异 | 结果 | 被哪个检查抓到 |
|---|---|---|
| 跳过速度投影（RATTLE） | CAUGHT | `max_velocity_tangency_residual` 等 |
| 跳过位置投影（SHAKE） | CAUGHT | `initial_target_distances_after_projection`, `energy_drift_per_atom_ps` |
| 约束自由度退回 `3N−3` | CAUGHT | `authoritative_constrained_dof` |
| 移除约束 virial | CAUGHT | `constraint_virial_cancels_rotational_ke`, `pressure_matches_analytic_ideal_gas_com` |
| 约束 virial 系数错一倍（`1/dt` 而非 `2/dt`，即时间层取错时会出现的值） | CAUGHT | 同上两条 |
| MPI 下 virial 再按 rank 数乘一遍 | CAUGHT | `mpi_agreement` |
| 漏掉一条跨 rank 约束 | CAUGHT | `authoritative_constrained_dof`, 两个残差检查, `mpi_agreement` |
| wrapped 坐标用错 image | CAUGHT | `wrapped_unwrapped_equivalence` |
| solver 在超出容差时报告成功 | CAUGHT | `per_frame_log_bounds`, `max_velocity_tangency_residual` 等 |
| 把 NaN 折成 0 再算误差 | CAUGHT | `_finite()` 在形成任何误差度量**之前**就拒绝 |
| restart 同时丢掉约束 virial **与**完成步压强 | CAUGHT | `restart_restores_completed_step_pressure` |
| restart 只丢掉约束 virial | **EXPECTED MISS** | — |
| restart 只丢掉完成步压强 | **EXPECTED MISS** | — |

最后两条**不是缺陷，是冗余**：restart 会恢复约束 virial，也会在 checkpoint 坐标上重新求值
provider virial 并单独恢复完成步压强记录；**任何一条单独存在就足以重建被报告的那个压强**，
而下一步的 RATTLE 无论如何都会重新算出约束 virial。因此单独丢一条在 CLI 的任何输出上都不可观测；
**两条同时丢就会被抓到**。约束 virial 张量本身在库层面由 `tests/mpi_constraint_virial.cpp` 覆盖。

### 本次工作中发现、但不属于本次改动范围的两个既有缺陷

两个都**与约束求解器无关**（都能在 `constraints off` 下复现），是在构造审计 fixture 时撞到的，
没有在本次提交中修复：

1. **`VerletNeighborBuilder::rebuild` 越界写**。盒子很大且坐标远离原点时 cell 索引越界，
   `EXC_BAD_ACCESS`。复现：3 原子、盒子 25000 Å、坐标约 (5137, 6411, 7229)、`lj_cutoff 0.95`。
   `constraints off` 同样崩溃。

2. **只含 `constraints` 段、不含任何 bonded 项的 topology 在 MPI 下死锁**，条件是某个 rank
   不拥有任何原子。加入哪怕一条 bond 即可避免。`constraints off` 同样死锁，所以与约束求解无关。
   本仓库的 fixture 都声明了（力常数为 0 的）bond 用于 exclusion，因此不会触发；
   `tests/constraint_convergence_audit.py` 里显式注明了这一点，以免被当成"随手加的"。

## 约束秩分析：fallback 语义与 near-tie 选择规则（2026-09-17）

### 分析管线里的两个判断，以及它们为什么要交叉核对

`analyze_independence()` 对每个连通分量分别回答两个**不同**的问题，用的是两种**不同**的方法：

| 问题 | 方法 | 阈值 |
|---|---|---|
| **有多少**行冗余 | 质量加权 Jacobian `J_M = J M^(-1/2)` 的奇异值（单边 Jacobi） | `tol = max(rows, cols) * eps * sigma_max` |
| **哪些**行冗余 | 列主元 modified Gram–Schmidt（`select_independent_columns`） | 同一个 `tol` |

第二个问题本质上更不可靠（主元残差范数只是 `sigma` 的上界），因此**以奇异值为准**，并把两者
**交叉核对**：当 MGS 接受的列数 `kept.size()` 不等于奇异值秩 `component.rank` 时，这个分量
**无法说出**哪几行是冗余的，于是退化为点名**整个分量**，并把
`redundant_constraints_identified` 置为 false。这就是 fallback。

**fallback 不改变接受/拒绝。** `require_independent()` 对**任何**秩亏集合都抛异常；fallback
只影响诊断里点名的是哪几条约束。一个不可信的秩**永远不会**被用来算自由度——被拒绝的初始化
不会留下速度场、thermostat 状态或完成的 rank report，`test_fallback_still_rejects_and_leaves_no_state()`
直接断言这一点。

### 是否存在自然的 fallback fixture：**有**

`tests/constraint_rank_fallback_tests.cpp` 用的是真实管线，不是手搓 report。

四个原子 + 全部六条两两距离约束：`3*4 - 6 = 6` 个内部自由度对六条约束，一般构型下**独立**。
让四个原子**共面**，四点六距离的 Cayley–Menger 行列式为零，六条距离之间出现一个关系，秩掉到 5。
于是**扫描离面高度 h** 会让 `sigma_min` **连续地**穿过 rank tolerance，而两种方法正是在这个穿越
区间里给出不同答案。

选这个形状而不是"近共线链"是刻意的：任何三个互相约束且接近共线的原子，会被
`reject_degenerate_target_triangles()` 在**构造期**就拒掉（它只看 target 距离），根本到不了秩分析。
四个近共面原子不含这种三元组。

**实测区间**约为 `h ∈ [2.6e-15, 4.4e-15]`，即距离共面几个 1e-15。本机上 31 个秩亏高度中有 20 个
命中 fallback。测试**扫描该区间并要求至少一次命中**，而不是钉死某一个 h——钉死一个值等于赌它
在每个编译器上都落在窄带内。

**诚实的范围说明**：这是**真实生产路径**的 fixture，不是注入的故障；但它也**不是**任何人做物理
建模会得到的构型。到达这个分支本来就需要这样的输入。文档不把它描述成"自然发生的物理 fallback"。

不变性已验证：约束顺序、原子存储顺序、旋转、平移、质量比全部不改变 fallback 结论。
MPI 下 np=1/2/4 一致，且有一个 fixture 被限制在单个 octant 内，使 2/4 rank 时存在**不拥有任何原子**
的 rank，它仍必须给出相同结论。

### near-tie 规则的精确定义

```
tie_relative_tolerance = sqrt(eps)                 = 1.4901161193847656e-08
cutoff                 = best_norm * (1 - sqrt(eps))
```

每一步主元选择中，所有**未被取走**且残差范数 `>= cutoff` 的列都算**并列**，并列者取
**canonical key 最小**的那个——`(min tag, max tag)`，即物理约束本身的属性，而不是它在输入列表里的位置。

- **eps**：`std::numeric_limits<double>::epsilon()` = 2.22e-16；`sqrt(eps)` 是"同一个量的两条不同
  算术路径算出来的结果不可区分"的常规判据。
- **零残差**：`best_norm > tolerance` 不成立时循环直接结束，不会走到 tie 比较；`cutoff` 也就不会
  出现 `0 * (1-x)` 的退化比较。
- **为什么绝对 `1e-300` 不合适**：`1 - 1e-300` 在 double 里**就是 1.0**，所以该判据实际等价于
  "只有 bit 完全相同才算并列"。物理等价的两列是经由不同算术路径得到范数的，只在舍入级别一致，
  于是绝对判据把它们判为不同，选择权就交回给了输入顺序/取向/单位——见下表。
- **承诺的确定性**：在**同一个浮点环境**内，结论只依赖物理约束集合，不依赖约束书写顺序、原子存储
  顺序、全局 tag 置换、刚体旋转、平移、均匀坐标缩放、质量比。
- **不承诺**：跨**根本不同**的浮点环境（不同的 FMA 收缩、不同的 `sqrt`/`cos` 实现、x87 扩展精度）
  逐位相同的选择。窗口宽达 `sqrt(eps)`，比这些差异大约八个数量级，因此实践中足够；但这是**余量论证**，
  不是逐位保证。

### near-tie 规则的**必要性**是怎么证明的

旧的对称 fixture 无法证明：它们是**轴对齐**的，并列范数**逐位相同**，两种判据给出同样答案。
`test_axis_aligned_tie_does_not_discriminate()` 把这条局限**钉在测试里**，以免日后被误当成覆盖。

有区分力的 fixture 是**一般性旋转**后的正方形 + 不等质量。旋转破坏逐位相同，但不动物理对称性。
判决步上四列的实测范数：

| 约束 | 残差范数 | 相对差 |
|---|---|---|
| (0,1) | 1.80277563773199456e+00 | 3.695e-16 |
| (0,3) | 1.80277563773199523e+00 | leader |
| (1,2) | 1.80277563773199456e+00 | 3.695e-16 |
| (2,3) | 1.80277563773199501e+00 | 1.232e-16 |

相对差**非零**（所以 1e-300 判为不同），且比 `sqrt(eps)` 小**八个数量级**（所以相对判据判为并列）。

在一次性 worktree 里把判据换回 `1e-300` 后，**14 项检查失败**，每一项都点名被选中的**物理约束**：

| 变化 | 当前（相对） | 旧（绝对 1e-300） |
|---|---|---|
| 基准 | (2,3) | **(1,3)** |
| 约束顺序/存储顺序置换 | (2,3) | **(1,3)** |
| 坐标 × 0.1 | (2,3) | **(0,1)** |
| 坐标 × 10 | (2,3) | **(1,2)** |
| 质量比 1.0 | (2,3) | **(1,2)** |

**均匀坐标缩放**是其中最锋利的一条：所有残差范数与 rank tolerance 一起缩放，相对判据**按构造**不变，
绝对判据则会移动。

### 反方向：窗口不能吞掉真实差异

四个**不规则共面四边形**的残差范数相差远大于 `sqrt(eps)`，主元顺序完全由大小决定。它们点名的分别是
**(1,2)、(0,3)、(0,1)、(0,2)**——都**不是**最大的 canonical key（tie 规则倾向于留下的那个），
且在**两种判据下答案相同**。若窗口宽到能吞掉真实差异，它们就会改口去点名 canonical-key 答案。

### 负控制（negative controls）

在一次性 worktree 里逐条弄坏生产代码、重建、跑相关测试，**要求测试失败**；全部还原后比对整棵树
checksum。变异从不提交。

| 变异 | 结果 |
|---|---|
| fallback 被伪造成"两种方法一致" | CAUGHT |
| fallback 的 dependent 索引被清空 | CAUGHT |
| report 级聚合用 OR 而非 AND | CAUGHT |
| canonical key 退回"按位置取第一个" | CAUGHT |
| near-tie 窗口退回绝对 `1e-300` | CAUGHT |
| near-tie 窗口放宽到 0.5（吞掉真实差异） | CAUGHT |
| MPI 各 rank 做出不同的 fallback 判断 | CAUGHT |
| Full MPI workflow 恢复过时的 `-R` 过滤器 | CAUGHT（registration guard） |

**其中一条负控制第一轮没抓到，并因此补强了测试**：把 report 级聚合从 AND 改成 OR 时，最初的
"一个 fallback 分量 + 一个独立分量"fixture **抓不到**——因为没有"已识别的秩亏分量"可供两种聚合
产生分歧，OR 和 AND 给出相同答案。补上"**一个已识别秩亏分量（5 原子全配对笼，10 约束 / 秩 9）
+ 一个 fallback 分量**"的 fixture 后才真正区分开。这正是负控制该起的作用：它暴露的是测试的
盲区，不是生产代码的缺陷。

运行方式：

```bash
python3 validation/run_validation.py \
  --case short \
  --gmd build/gmd \
  --gmd-validate build/gmd_validate \
  --work-root build/validation
```

或运行完整套件：

```bash
python3 validation/run_validation.py \
  --case all \
  --gmd build/gmd \
  --gmd-validate build/gmd_validate \
  --work-root build/validation
```

PME 外部参考的 fixture、对齐范围、收敛表、容差推导与遗留限制记录在 `validation/pme_external/README.md`。
