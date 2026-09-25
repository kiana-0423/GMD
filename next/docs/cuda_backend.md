# CUDA 后端设计

首后端：CUDA。P0 已提供可选工具链配置与 Runtime/CUB 归约探针；
下表的 MD 映射是后续实现计划。

## 库与 kernel 分工

| 工作 | 初始实现选择 | 边界 |
|---|---|---|
| device、内存、stream/event、拷贝 | CUDA Runtime | 后端内部拥有；上层显式提交与等待 |
| 原子/边局部函数、几何和积分 | 自定义 CUDA kernel | 可以融合 Gather、几何、导数，不构造稠密矩阵 |
| 全局能量、位移界、维里归约 | CUB DeviceReduce / BlockReduce | 显式 stream 与可复用 workspace |
| cell/邻居分组与偏移 | CUB 排序、scan；配合自定义 kernel | 构图需单独验证容量、重复和漏边 |
| 力装配 | 首版半表 FP64 原子累加；保留分段归约/全表实现路径 | 首版不宣称确定性；依硬件能力检查所需原子操作 |
| 分段归约 | 有正确分组/offset 后使用 CUB 相应原语 | 无序 scatter 不能直接当作连续分段输入 |
| PME FFT（后续） | cuFFT | charge spread、Green 函数、归一化、gather、自能等另行实现 |
| 多 GPU 通信（后续） | 独立 transport 接口 | 先定义 owned/ghost，再选择 CUDA-aware MPI 或 staging |

不预先引入 cuBLAS/cuSPARSE：目前二体和超边主要是索引访问、局部计算和
归约。若后续模型确有稠密或稀疏线性代数需求，再添加相应库。

## stream 与 workspace

P1.2 已按下文编写 `GMDNext::cuda_runtime`（未经 nvcc 编译、未在 GPU 运行，
细节见[实施规划第 9 节](program_plan.md#9-p12-实现记录)）。缓冲区使用流序
`cudaMallocAsync`/`cudaFreeAsync`，因此单 stream 下释放自动排在已提交工作之后；
设备要求 sm_60+（FP64 `atomicAdd`）并支持内存池，否则创建上下文时拒绝。

第一版只使用显式单 stream。CUB 采用查询 scratch 大小、分配、提交的接口；
容量变化时重新查询。将 scratch 视为写资源，跨 stream 复用要有事件依赖。
正常步重用分配；错误不能靠每次 `cudaDeviceSynchronize` 隐藏依赖问题。

未来 cuFFT plan 绑定明确的 stream 与 workspace，保存布局和归一化约定。
CUDA Graphs 等到固定容量的单步路径通过验证后再引入，构图扩容和拓扑变化
必须更新或重新捕获执行计划。

## 输出捕获边界（设计中）

初期仅在输出点等待相关事件，按请求字段下载；不引入每步全设备同步。
订阅在力计算前进入计划，动力学必需的观测量始终保留，输出请求可增加采样。
正常计算的设备状态布局不由 writer 指定。

后续异步路径在 compute stream 将所需数据打包到独立输出槽，copy stream
等待捕获事件后下载到有界 pinned 主机池；writer 只读取下载完成的槽。
下一步只能在捕获后改写原状态，设备槽在下载后、主机槽在消费后才可复用。
打包带来的设备读写与实际重叠收益需要和同步路径比较，不能假定异步总是更快。
完整生命周期与验收见[输出层设计](output_pipeline.md#5-捕获下载与写入的生命周期)。

## 配置与探针

默认 `GMD_NEXT_ENABLE_CUDA=OFF` 允许无 CUDA 的机器检查参考层。
这不改变首后端选择，也不表示 CPU 是运行时自动回退。

开启后检查 CUDA 编译器并查找 Toolkit，编译 `gmd_next_cuda_probe`。
探针在一个显式 stream 上使用 CUB 对四个 FP64 数求和，验证结果并输出
GPU 名称、SM 版本、driver/runtime 版本。缺少设备、库调用失败或结果错误
均返回非零。只有 `GMD_NEXT_TEST_CUDA=ON` 才将它注册为 GPU 测试。

P0 探针使用少量临时分配并在结束前同步，只检查开发环境，不代表未来
热路径的内存管理或性能实现。当前没有 GPU 力计算目标，不能据探针通过
声称 GPU MD 已实现。

选择具体 Toolkit、驱动和 `CMAKE_CUDA_ARCHITECTURES` 后，记录到运行报告。
首次 GPU 验证前，这些组合均属于待验证状态。

## 依据

以下链接用于 API 与资源语义核查；不作为本项目的性能证据。

- [NVIDIA CUDA Best Practices：主机/设备数据传输](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/#data-transfer-between-host-and-device)
- [CUB DeviceReduce API](https://nvidia.github.io/cccl/cub/api/structcub_1_1DeviceReduce.html)
- [CUB device-scope 与 workspace](https://nvidia.github.io/cccl/cub/developer/device_scope.html)
- [cuFFT 文档](https://docs.nvidia.com/cuda/cufft/index.html)

CCCL 在线文档随版本变化。当前探针采用显式 scratch/stream 接口；发布前
以实际固定的 Toolkit 所带头文件和 GPU 测试作为兼容性证据。
