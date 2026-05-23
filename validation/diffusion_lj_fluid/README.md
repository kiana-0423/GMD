# diffusion_lj_fluid

目的：

- 验证 LJ fluid 的 MSD 曲线与扩散系数提取工作流。

方法：

- 从 `output.xyz` 读取坐标，按固定盒长做最小镜像反缠绕。
- 对 MSD 后半段做线性拟合，得到 `D = slope / 6`。

说明：

- 当前基线仍是 `provisional_gmd_baseline`，主要先验证 workflow 稳定性与输出格式。
