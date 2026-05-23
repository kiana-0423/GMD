# static_coulomb

目的：

- 验证带电小体系的 LJ + Coulomb 单点能量与力。
- 分别覆盖 `Ewald` 与 `PME`。

当前基线来源：

- `reference_ewald.json` 由 `validation/analytic_references.py` 生成，
  是独立 Python 周期 Ewald 解析参考，不使用 GMD 输出。
- `reference_pme.json` 仍是 `provisional_gmd_baseline`。PME mesh 结果需要
  LAMMPS PPPM 或 OpenMM PME 外部参考后才能升级为 independent validation。
- 目录内提供了匹配参数的 `lammps_ewald.in` 和 `lammps_pme.in`，便于后续替换为外部参考软件结果。

固定参数：

- 截断 `8.0 A`
- `ewald_alpha = 0.3`, `ewald_kmax = 3`
- `pme_alpha = 0.3`, `pme_order = 4`, `pme_grid = 16 16 16`

解析 Ewald 公式：

- `U_real = k_e sum q_i q_j erfc(alpha r_ij) / r_ij`
- `U_recip = 1/2 sum_{k!=0} k_e 4 pi exp(-k^2/4 alpha^2) |S(k)|^2 / (V k^2)`
- `U_self = -k_e alpha / sqrt(pi) sum q_i^2`
- `U_background = -k_e pi Q^2 / (2 V alpha^2)`
