# berendsen_npt_lj

目的：

- 在真实轨迹上验证 **Berendsen barostat**。此前仓库里没有任何 case 跑过它：
  `npt_lj_fluid` 用的是 Monte Carlo barostat，Berendsen 路径只有直接单元测试。
- 该 case 是 **dynamics regression validation，不是外部科学参考**。所有数值均由
  GMD 自身产生，没有运行任何外部引擎，也没有任何值来自理论推导。不要将其引用为
  Berendsen barostat 的 cross-code validation。

## 结构：符号敏感的成对运行

同一压缩后的 32 原子 LJ fixture，两个 run 只有目标压强不同：

| run | 目标压强 | 期望 |
|---|---|---|
| `expand.run` | 200 bar（**低于** fixture 运行期间约 2300 bar 的平均压强） | 盒子必须**变大** |
| `compress.run` | 4000 bar（**高于**平均压强） | 盒子必须**变小** |

`analyze.py` 中的方向与幅度断言不是 regression 值，不依赖 `reference.json`：
它们对任何正确的 barostat 都成立，也是这个 case 存在的理由。

## 它能抓住的两个生产缺陷

1. **bar 与 `eV/Å³` 混比**（已在 `fix: convert the Berendsen barostat's pressure
   comparison to bar` 中修正）。把 200 bar 当作 200 eV/Å³ 就是 3.2e8 bar，高于任何
   瞬时压强，于是**两个 run 都会压缩**。实测：expand 比值 0.9911、compress 0.8355，
   `expand_direction` 与 `directions_differ` 直接失败。

2. **relaxation time 未做单位换算**（已在 `fix: convert thermostat and barostat
   relaxation times to internal units` 中修正）。方向仍然正确，但耦合弱约十倍。
   实测：expand 比值 1.0129、compress 0.9962，最大响应仅 0.0129；正确代码最小响应
   为 0.0490。`min_volume_response = 0.03` 正落在两者之间。

## 其它检查

- 全程体积有限、为正、且不越出初始体积的 `[0.25, 4.0]` 倍（排除失控）。
- 全程 PE / KE / E_total / T / P 有限，温度非负。
- **restart 连续性**：`restart_first.run` 跑 250 步并写 checkpoint，
  `restart_second.run` 续跑 250 步，最终体积必须与连续 500 步运行**逐位相同**。
- **serial/MPI 一致性**：传入 `--mpiexec` 时重跑 expand 并比较最终体积；
  Berendsen 的缩放因子来自全局归约量，本身与 rank 数无关。

## 容差依据

全部由实测确定，记录在 `reference.json` 的 `tolerance_basis` 中：温度与压强的
mean 容差取三倍均值标准误（32 原子体系涨落极大），volume ratio 容差取正确值与
缺陷值差距的约五分之一，`min_volume_response` 取两者实测响应之间。

固定参数：

- fixture：`lj32_compressed.xyz`（`npt_lj_fluid` 的 lj32 按 0.74 等比压缩）
- `time_step 2.0`，`run 500`，`cutoff 6.0`
- `thermostat nose_hoover`，`thermostat_tau 100.0`
- `barostat berendsen`，`barostat_tau 1000.0`，`compressibility 4.5e-5`
- `velocity 300.0`，`velocity_seed 20260830`（确定性）
