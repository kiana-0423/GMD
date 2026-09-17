#include <cuda_runtime.h>
#include <cub/device/device_reduce.cuh>

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void check(cudaError_t status, const char* operation) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
    }
}

// Only an environment probe. Production storage/workspace ownership is P1.
struct Resources {
    cudaStream_t stream = nullptr;
    double* input = nullptr;
    double* output = nullptr;
    void* scratch = nullptr;

    Resources() = default;
    Resources(const Resources&) = delete;
    Resources& operator=(const Resources&) = delete;
    ~Resources() {
        if (stream) cudaStreamSynchronize(stream);
        if (scratch) cudaFree(scratch);
        if (output) cudaFree(output);
        if (input) cudaFree(input);
        if (stream) cudaStreamDestroy(stream);
    }
};

}  // namespace

int main() {
    try {
        int count = 0;
        check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
        if (count == 0) {
            throw std::runtime_error("CUDA probe requires a visible NVIDIA GPU");
        }
        check(cudaSetDevice(0), "cudaSetDevice");
        cudaDeviceProp properties{};
        check(cudaGetDeviceProperties(&properties, 0), "cudaGetDeviceProperties");
        int driver = 0;
        int runtime = 0;
        check(cudaDriverGetVersion(&driver), "cudaDriverGetVersion");
        check(cudaRuntimeGetVersion(&runtime), "cudaRuntimeGetVersion");

        const std::array<double, 4> input{1.0, -2.0, 3.0, 4.0};
        double output = 0.0;
        Resources resources;
        check(cudaStreamCreateWithFlags(&resources.stream, cudaStreamNonBlocking), "create stream");
        check(cudaMalloc(reinterpret_cast<void**>(&resources.input), sizeof(input)), "allocate input");
        check(cudaMalloc(reinterpret_cast<void**>(&resources.output), sizeof(output)), "allocate output");
        check(cudaMemcpyAsync(resources.input, input.data(), sizeof(input),
                              cudaMemcpyHostToDevice, resources.stream), "upload input");

        std::size_t scratch_bytes = 0;
        check(cub::DeviceReduce::Sum(nullptr, scratch_bytes, resources.input, resources.output,
                                     static_cast<int>(input.size()), resources.stream), "query CUB scratch");
        check(cudaMalloc(&resources.scratch, std::max(scratch_bytes, std::size_t{1})), "allocate scratch");
        check(cub::DeviceReduce::Sum(resources.scratch, scratch_bytes, resources.input, resources.output,
                                     static_cast<int>(input.size()), resources.stream), "CUB sum");
        check(cudaMemcpyAsync(&output, resources.output, sizeof(output),
                              cudaMemcpyDeviceToHost, resources.stream), "download result");
        check(cudaStreamSynchronize(resources.stream), "complete probe");
        if (output != 6.0) {
            throw std::runtime_error("CUB reduction returned an incorrect result");
        }
        std::cout << "CUDA/CUB environment probe passed; GPU MD operators are not implemented.\n"
                  << "device=" << properties.name << " sm=" << properties.major << properties.minor
                  << " driver=" << driver << " runtime=" << runtime << " sum=" << output << '\n';
    } catch (const std::exception& error) {
        std::cerr << "CUDA probe failed: " << error.what() << '\n';
        return 1;
    }
}
