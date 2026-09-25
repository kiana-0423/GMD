# CUDA 实现入口

- `probe.cu`：检查 CUDA Runtime、CUB FP64 归约和可见 GPU。它不是力计算后端。
- `include/gmd_next/cuda/runtime.hpp`、`runtime/`：P1.2 执行上下文、完成状态、
  设备缓冲区与 workspace（`GMDNext::cuda_runtime`）。公共头不含 CUDA 类型；
  `runtime/native.hpp` 仅供后端内部 `.cpp/.cu` 使用。
- `tests/runtime_tests.cu`：GPU 检查，仅 `GMD_NEXT_TEST_CUDA=ON` 时注册。

P1.2 代码尚未经 nvcc 编译或在 GPU 上运行；在目标机通过前不视为完成。探针采用 CUDA C++17，与参考层的 C++20 API 没有依赖；
新 CUDA 目标使用 CUDA C++20，以共享 contracts 头文件；不从探针推断完整
引擎兼容性。

构建命令见 [主 README](../../README.md)，库与 kernel 分工见
[CUDA 设计](../../docs/cuda_backend.md)。CUB 使用 Toolkit 附带的头文件；
没有联网下载或额外安装逻辑。cuFFT 只在 PME 实现阶段引入。
