# P0 搭建验证记录

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
