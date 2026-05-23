# static_special_pairs

目的：

- 验证分子 `1-2 / 1-3` exclusion 与 `1-4` scaling。
- 同时覆盖 LJ 与 Coulomb，且输出分量能量和力误差。

GMD 设置：

- 四原子直链，固定坐标
- `lj_scale_14 = 0.5`
- `coul_scale_14 = 0.833333333333`
- Coulomb 使用 Ewald real + reciprocal + self

参考说明：

- `reference.json` 由 `validation/analytic_references.py` 生成，来源为独立
  Python 解析参考，而不是 GMD 输出。
- 解析参考包含 shifted LJ、周期 Ewald real/reciprocal/self/background
  项，以及 special-pair 直接修正 `(scale - 1) k_e q_i q_j / r`。
- 同一个 case 还验证修改 `lj_scale_14 = 0.25`、
  `coul_scale_14 = 0.5` 后，1-4 energy/force 按解析公式变化。
- `lammps.in` 对应 `special_bonds` 设定，便于后续导入外部软件结果。
