# GMD Next 开发边界

- `next/` 是独立项目。默认只修改此目录；不要将新目标接入父目录的
  `CMakeLists.txt`，也不要调整现有引擎来迁就新架构。
- 首个硬件后端为 NVIDIA CUDA。CPU 代码承担数学参考职责。
- 新代码使用 `gmd_next` 命名空间；不包含父项目的 `gmd/*` 头文件，
  不链接 `gmd_core`，不依赖旧 `System` / `Simulation` / `ForceProvider`。
- 数学依据为 `../paper/md_graph_tensor_preprint.tex`；工程契约见
  `docs/operator_contracts.md`，实现状态见 `docs/roadmap.md`。
- 标明“已实现”“设计中”“待 GPU 验证”的差别。CUDA 环境探针不代表
  已实现 GPU 力计算；没有硬件结果时，不声称 GPU 测试或性能验证通过。
- 保持物理单位、位移方向、截断、排除、缩放、维里及自由度约定明确。
  CPU/GPU 结果相同不能替代独立物理检查。
- 新的 GPU 热路径不得通过旧主机对象逐步交换状态；缓冲区生命周期、
  stream/event 顺序和错误状态必须显式定义。
- 验证使用独立构建目录。只运行与修改相关的检查；不要自动更新父项目
  的物理基线，也不要为让测试通过而放宽科学误差阈值。
