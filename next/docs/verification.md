# 验证记录

## P0 搭建

日期：2026-09-17。环境：Apple arm64，AppleClang 21.0.0，CMake 4.4.2。

| 检查 | 结果 |
|---|---|
| 独立 CPU Release 配置与编译 | 通过 |
| `gmd_next.operators` | 通过：Gather/重复索引装配、伴随、空关系及输入形状 |
| `gmd_next.lj` | 通过：解析 LJ、力有限差分、九分量维里形变导数、半/全表、排列、缩放、PBC、cutoff 及错误输入 |
| CPU 示例 | `U=0, F0.x=-24, F1.x=24, Wxx=24`，与解析双粒子值一致 |
| 将源码复制到临时独立目录后 configure/build/test | 通过，2/2；无父项目头文件或库依赖 |
| 明确请求 CUDA，但本机缺少 nvcc | 按设计清楚失败，没有静默回退 CPU |
| 搭建前记录的 300 个现有文件（含论文）SHA-256 比较 | 全部保持不变 |
| 新源码 include 边界及本地文档链接 | 通过 |
| CUDA 探针编译、GPU 执行、GPU 算子正确性/性能 | **未验证**；本机无 CUDA 编译器，MD GPU 算子尚未实现 |

默认 macOS SDK 搜索首先选到了与当前 linker 不兼容的 CommandLineTools
SDK，普通编译器探针报 `unknown architecture arm64e.x1`。此次验证仅对
新项目的构建缓存显式指定已安装的 Xcode SDK 后通过：

```sh
cmake -S next -B next/build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_SYSROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
cmake --build next/build --parallel
ctest --test-dir next/build --output-on-failure
```

没有修改系统 SDK 选择、父项目构建配置或旧测试基线。旧引擎测试未在本次
重复运行；隔离证据来自未变的现有文件、独立构建图及脱离父目录的构建。
本记录不声称预印本原始 Python 实验已复现。

## P1.1 主机契约

日期：2026-09-21。环境：Apple arm64，AppleClang 21.0.0，CMake 4.4.2；
与 P0 相同，需显式指定 Xcode SDK。构建目录 `next/build-p11`（Release）与
`next/build-p11-asan`（Debug，`-fsanitize=address,undefined`）。

| 检查 | 结果 |
|---|---|
| Release 配置与编译，`-Wall -Wextra -Wpedantic` | 通过，0 条警告 |
| `gmd_next.operators`、`gmd_next.lj`（已有参考检查） | 通过；reference 头文件、源码、测试与演示无改动 |
| `gmd_next.contracts` | 通过：ID/局部索引容量、版本匹配与失效域、三态有效性、单位换算、观测请求合并、记录的缺失/无效/版本检查 |
| `gmd_next.model` | 通过：生产/参考两种范围的接受与拒绝、规范化排序、去质心与 DOF、oracle 映射 |
| 单位换算对照 | 由 SI 定义独立推导：`1e5*sqrt(m_u/e)` fs、`e/1e-30/1e5` bar、`1.380649e-23/e` eV/K，相对误差 ≤1e-12 |
| oracle 对照 | 手算 ε=σ=1、r=1、rc=2.5：potential_shift `U=0.016316891136, F=±24, Wxx=24`；force_shift `F=24.0389994774528` |
| Debug + ASan/UBSan | 4/4 通过，无运行时报告 |
| 变异检查（临时副本） | 6 个人为缺陷（rc+skin 边界、DOF、fs 方向、未设版本匹配、请求能量默认 0、重复 ID 检查）均使测试失败；其中“默认 0”初次未被发现，已补测试 |
| 将 `next/` 复制到临时目录独立构建测试 | 通过，4/4；无父项目头文件或库依赖 |
| 明确请求 CUDA，但本机缺少 nvcc | 按设计清楚失败 |
| **CUDA 资源、设备状态、GPU 力计算、性能** | **未实现、未验证**；本次无 GPU 代码 |

这些都是主机侧检查，不说明 GPU 正确性或性能。

```sh
cmake -S next -B next/build-p11 -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_SYSROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
cmake --build next/build-p11 --parallel
ctest --test-dir next/build-p11 --output-on-failure
```

## P1.2 CUDA 资源

日期：2026-09-21。本机仍无 `nvcc`、`nvidia-smi` 与可用 GPU；Docker 守护进程未
运行，未使用容器。

| 检查 | 结果 |
|---|---|
| `gmd_next.runtime`（主机簿记） | 通过：扩容几何/上限/溢出、view 代际与尺寸/设备/缓冲区检查、释放、完成票据、故障位解码 |
| 变异检查 | 5 个人为缺陷（尺寸变化不失效、忽略未知故障位、扩容不钳制、允许两个未决完成、无效设备相等）均被发现 |
| Release、ASan/UBSan、独立复制构建 | 均 5/5 通过，0 条警告 |
| `runtime/context.cpp`、`buffer.cpp` 类型检查 | 用 NVIDIA 的 CUDA 12.4 头文件 wheel（`nvidia-cuda-runtime-cu12`、`nvidia-cuda-nvcc-cu12`，仅解压到临时目录），AppleClang `-fsyntax-only -Wall -Wextra -Wpedantic`：无诊断 |
| `tests/runtime_tests.cu` 前端解析 | clang CUDA 模式，主机端与 `sm_60` 设备端均无诊断；因临时工具链版本识别问题，需一个声明 `cudaConfigureCall` 的 shim 才能解析 `<<<>>>`，并以空 `curand_mtgp32_kernel.h` 占位 |
| 公共头 `gmd_next/cuda/runtime.hpp` | 无任何 CUDA 头文件时可编译；公共头中无 CUDA 类型 |
| **nvcc 编译、链接 cudart、GPU 运行 `gmd_next.cuda_runtime`** | **未执行**；P1.2 退出条件未满足 |

语法检查不能代替 nvcc 编译，更不能代替 GPU 运行。以下行为只由尚未运行的
GPU 测试覆盖：`cudaMallocAsync` OOM 返回码及恢复、`cudaEventQuery` 未完成时
不设置 last error 的假设、流序释放与上下文析构顺序、故障字回读与清零、
启动配置错误使上下文失效。

目标机验收命令：

```sh
cmake -S next -B next/build-cuda -DCMAKE_BUILD_TYPE=Release \
  -DGMD_NEXT_ENABLE_CUDA=ON -DGMD_NEXT_TEST_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=native
cmake --build next/build-cuda --parallel
ctest --test-dir next/build-cuda --output-on-failure
```
