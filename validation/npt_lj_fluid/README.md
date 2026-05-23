# npt_lj_fluid

目的：

- 验证带温压耦合工作流下的温度与压力统计量。

统计说明：

- 当前 GMD NPT case 使用 `Monte Carlo` barostat。
- 因此压力统计来自每个输出步的瞬时 virial 压力，而不是 barostat 内部 acceptance 历史。
- `reference.json` 先记录 `provisional_gmd_baseline`，后续再替换成外部软件参考。
