# GMD Next：图—张量算子与 CUDA 执行架构

本目录在 GMD 仓库内独立搭建下一代引擎。数学依据是
[图—张量预印本](../paper/md_graph_tensor_preprint.tex)，首个硬件后端确定为
**NVIDIA CUDA**。目录、构建目标、命名空间和执行状态均与当前引擎隔离。

当前交付是 **P0 架构与参考骨架** 加 **P1.1 主机契约**：可运行 CPU 数学
检查；`GMDNext::contracts` 提供稳定 ID、状态版本、有效性、单位换算、
观测请求及模型/运行校验（仅 CPU 验证）；可选择构建 CUDA Runtime / CUB
环境探针与 P1.2 的 `GMDNext::cuda_runtime`（上下文、完成状态、设备缓冲区、
workspace；**尚未在 CUDA 机上编译或运行**）。GPU 力算子、GPU 邻居构建和 MD 时间推进尚未实现。CPU 参考层中的 `std::vector` 不定义未来设备状态的存储接口。

## 入口

| 文档 | 内容 |
|---|---|
| [程序实施规划](docs/program_plan.md) | 首个交付范围、文件/构建职责、接口、任务依赖和验收顺序 |
| [架构](docs/architecture.md) | 层次、所有权、设备常驻状态、执行依赖、与旧引擎及现有输出接口的边界 |
| [输出层设计](docs/output_pipeline.md) | 轻量记录、共享序列化、独立频率、缓冲区生命周期与兼容契约 |
| [算子契约](docs/operator_contracts.md) | 数学定义、计数、周期镜像、有效性、单位、错误处理 |
| [CUDA 后端](docs/cuda_backend.md) | Runtime / CUB / 自定义 kernel / cuFFT 分工与构建 |
| [实施计划](docs/roadmap.md) | 当前状态、分阶段交付及验收条件 |
| [验证方案](docs/validation.md) | 独立物理参考、CPU/GPU 对照、正确性及性能测量 |
| [本次验证记录](docs/verification.md) | 实际执行结果、环境与尚未验证的 CUDA 部分 |

## 独立构建

在仓库根目录运行，无需 CUDA 或第三方测试框架：

```sh
cmake -S next -B next/build -DCMAKE_BUILD_TYPE=Release
cmake --build next/build --parallel
ctest --test-dir next/build --output-on-failure
./next/build/gmd_next_reference_demo
```

在具备 CUDA Toolkit 与 NVIDIA GPU 的开发机上：

```sh
cmake -S next -B next/build-cuda \
  -DCMAKE_BUILD_TYPE=Release \
  -DGMD_NEXT_ENABLE_CUDA=ON \
  -DGMD_NEXT_TEST_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=native
cmake --build next/build-cuda --parallel
ctest --test-dir next/build-cuda --output-on-failure
```

`native` 要求配置时能识别 GPU。交叉编译或无 GPU 的构建机需按部署设备
显式设置架构编号，并关闭 `GMD_NEXT_TEST_CUDA`。最低配置要求为 CMake
3.24、C++20；CUDA 探针要求 Toolkit 12.0 或更新版本。最低版本是构建约束，
不是已完成的兼容性测试矩阵。

父项目的 `cmake -S . ...` 不会包含这里的目标；这里也不读取父项目的
CMake 模块或头文件。可将整个 `next/` 复制到其他位置独立构建。
论文链接和未来旧引擎对照路径仅用于开发文档与外围验证。

## 当前目录

```text
next/
├── CMakeLists.txt             独立项目与配置开关
├── docs/                     设计与验收契约
├── include/gmd_next/         CPU 数学参考 API（reference.hpp）
│   ├── core/                 ID、版本、有效性/错误、单位、观测请求
│   └── model/                LJ 模型、盒、主机状态、运行规格与校验
├── src/                      参考实现与主机契约实现
├── app/                      小系统参考演示
├── backends/cuda/            CUDA 探针与运行时资源（include/、runtime/、tests/）
└── tests/                    数学与物理契约检查
```

后续 runtime、storage、relations、operators、potentials、integrators、
communication 和 io 模块按里程碑引入，职责已在架构文档定义；不以空实现
或静默 CPU 回退伪装 GPU 支持。

后续输出采用“订阅驱动的按需捕获 → 轻量只读记录 → 独立序列化层”的路线，
保留现有文件格式与旧 API 兼容入口；旧 `System` 桥接只作过渡对照。
日志、轨迹和 checkpoint 分别采样，新核心不依赖旧引擎。目前仅完成设计，
同步接入、续跑验收及后续异步优化见[实施计划](docs/roadmap.md)。
