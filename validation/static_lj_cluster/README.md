# static_lj_cluster

目的：

- 验证固定坐标下的 LJ 总能量与原子力。
- 输出 `lj/bonded/coulomb/total` 分量能量，以及 force RMS / max error。

GMD 设置：

- `cutoff = 8.5 A`
- `pair_modify shift yes` 等价于 GMD 当前截断平移 LJ
- 不积分，只做单点力评估

参考说明：

- `reference.json` 设计为解析两体求和基线。
- `lammps.in` 给出对应的 LAMMPS 复现实验输入草案。

当前状态：

- 框架已接入自动分析。
- `reference.json` 将由本次基线运行填充为确定值。
