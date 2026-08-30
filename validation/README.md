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
