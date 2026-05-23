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
- `static_coulomb/reference_pme.json`：仍为 provisional GMD PME regression baseline，等待 LAMMPS PPPM 或 OpenMM PME 外部参考。
- `pme_external/`：PME 外部验证设计与待办事项；当前不包含已完成 reference。
- 长时间 NVE/NVT/NPT/diffusion cases：仍为 provisional workflow/regression baselines。

当前 release-facing 状态：

- LJ：已通过解析 reference 验证。
- Ewald：已通过解析 periodic Ewald reference 验证。
- replicated PME：可运行并进入测试，但仍只有 provisional regression baseline。
- `pme_mode distributed`：仅为 interface/prototype，数值 backend 仍是 replicated PME。
- SHAKE/RATTLE：串行与一个 cross-rank MPI correctness-first global-gather 路径已有测试；不是 scalable distributed constraint solver。
- checkpoint/restart：真实 `gmd` CLI restart-continuity 测试覆盖 serial、MPI np=2、MPI np=4。

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

PME 外部参考验证的准备工作记录在 `validation/pme_external/README.md`。
