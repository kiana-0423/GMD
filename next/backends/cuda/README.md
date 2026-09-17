# CUDA 实现入口

当前只有 `probe.cu`：检查 CUDA Runtime、CUB FP64 归约和可见 GPU。
它不是力计算后端。探针采用 CUDA C++17，与参考层的 C++20 API 没有依赖；
实际算子的语言要求在 P1 引入时明确，不从探针推断完整引擎兼容性。

构建命令见 [主 README](../../README.md)，库与 kernel 分工见
[CUDA 设计](../../docs/cuda_backend.md)。CUB 使用 Toolkit 附带的头文件；
没有联网下载或额外安装逻辑。cuFFT 只在 PME 实现阶段引入。
