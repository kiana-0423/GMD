// GPU checks for the P1.2 runtime. Registered only with GMD_NEXT_TEST_CUDA=ON;
// without a device this is unverified, not passed.

#include "runtime/native.hpp"

#include "test_support.hpp"

#include <array>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

namespace {

using namespace gmd_next;
using core::ContractError;
using core::ErrorCode;

__global__ void raise_fault(std::uint32_t* word, std::uint32_t bits) {
    if (blockIdx.x == 0 && threadIdx.x == 0) atomicOr(word, bits);
}

__global__ void scale(double* values, std::size_t count, double factor) {
    const std::size_t i = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x;
    if (i < count) values[i] *= factor;
}

bool has_code(const core::ValidationReport& report, ErrorCode code) {
    return report.contains(code);
}

void device_selection() {
    throws_code<ContractError>(ErrorCode::device_unavailable,
                               [] { cuda::Context::create(runtime::DeviceId{-1}); });
    int count = 0;
    cudaGetDeviceCount(&count);
    throws_code<ContractError>(ErrorCode::device_unavailable,
                               [&] { cuda::Context::create(runtime::DeviceId{count}); });
}

void round_trip_and_kernel(cuda::Context& context) {
    cuda::DeviceBuffer<double> buffer(context);
    const std::vector<double> input{1.0, -2.5, 3.25, 0.0, 7.0};
    buffer.resize(input.size(), runtime::Contents::discard);
    buffer.upload(input);
    auto view = buffer.view();
    scale<<<1, 32, 0, cuda::native_stream(context)>>>(view.data(), view.size(), 2.0);
    cuda::check_launch(context, "scale");
    std::vector<double> output(input.size());
    buffer.download(output);
    auto completion = context.record_completion();
    expect(completion.wait().succeeded(), "Round trip completion failed");
    for (std::size_t i = 0; i < input.size(); ++i) near(output[i], 2.0 * input[i], 0.0, 0.0);

    throws_code<ContractError>(ErrorCode::inconsistent_size,
                               [&] { buffer.upload(std::span(input).first(2)); });
}

void growth_invalidation_and_reuse(cuda::Context& context) {
    cuda::DeviceBuffer<double> buffer(context);
    const std::vector<double> input{4.0, 5.0, 6.0};
    buffer.resize(input.size(), runtime::Contents::discard);
    buffer.upload(input);
    const auto small = buffer.view();

    buffer.resize(1000, runtime::Contents::preserve);
    expect(has_code(buffer.check(small), ErrorCode::stale_view), "Growth must stale old views");
    const auto large = buffer.view();
    expect(large.data() != small.data(), "Growth past capacity must reallocate");
    std::vector<double> grown(1000);
    buffer.download(grown);
    expect(context.record_completion().wait().succeeded(), "Preserve completion failed");
    for (std::size_t i = 0; i < input.size(); ++i) near(grown[i], input[i], 0.0, 0.0);

    // Within capacity nothing is allocated: same address, same capacity.
    const auto capacity = buffer.capacity();
    buffer.resize(800, runtime::Contents::preserve);
    expect(buffer.capacity() == capacity && buffer.view().data() == large.data(),
           "Shrinking must reuse the allocation");
    expect(has_code(buffer.check(large), ErrorCode::stale_view), "A size change stales views");
    const auto kept = buffer.view();
    buffer.resize(800, runtime::Contents::preserve);
    expect(buffer.check(kept).ok(), "A no-op resize must keep views valid");

    // A view claiming another device is rejected before any launch uses it.
    const runtime::DeviceSpan<double> foreign(kept.data(), kept.size(),
                                              runtime::DeviceId{context.device().ordinal() + 1},
                                              kept.identity());
    expect(has_code(buffer.check(foreign), ErrorCode::device_mismatch), "Device mismatch missed");
    expect(has_code(context.check_device(foreign.device(), "view"), ErrorCode::device_mismatch),
           "Context accepted a foreign view");

    cuda::DeviceBuffer<double> moved(std::move(buffer));
    expect(moved.check(kept).ok(), "Moving a buffer must keep its views valid");
    throws_code<ContractError>(ErrorCode::malformed_input,
                               [&] { buffer.resize(1, runtime::Contents::discard); });

    cuda::Workspace workspace(context);
    workspace.reserve(256);
    const auto scratch = workspace.view();
    workspace.reserve(64);
    expect(workspace.size_bytes() == 256 && workspace.check(scratch).ok(),
           "A smaller workspace request must not reallocate");
}

void allocation_failure_is_recoverable(cuda::Context& context) {
    cuda::DeviceBuffer<double> buffer(context);
    buffer.resize(4, runtime::Contents::discard);
    const auto before = buffer.view();
    const std::size_t impossible = context.properties().global_memory_bytes;  // x8 bytes
    throws_code<ContractError>(ErrorCode::allocation_failed,
                               [&] { buffer.resize(impossible, runtime::Contents::preserve); });
    expect(buffer.size() == 4 && buffer.check(before).ok(), "Failed growth changed the buffer");
    expect(!context.has_failed(), "Out of memory must not poison the context");
    expect(context.record_completion().wait().succeeded(), "Context unusable after OOM");
}

void device_faults(cuda::Context& context) {
    auto word = context.fault_word();
    raise_fault<<<1, 1, 0, cuda::native_stream(context)>>>(
        word.data(), runtime::bit(runtime::DeviceFault::capacity_exceeded));
    cuda::check_launch(context, "raise_fault");
    auto faulted = context.record_completion();
    const auto& status = faulted.wait();
    expect(status.state == runtime::CompletionState::failed, "A device fault must fail completion");
    expect(status.failures.contains(ErrorCode::capacity_overflow), "Fault bit was misdecoded");
    // Recoverable: the word was cleared and the context keeps working.
    expect(!context.has_failed(), "A device fault must not poison the context");
    expect(context.record_completion().wait().succeeded(), "Fault word was not reset");

    auto pending = context.record_completion();
    throws_code<ContractError>(ErrorCode::resource_busy, [&] { context.record_completion(); });
    expect(pending.wait().succeeded(), "Outstanding completion failed");
}

void launch_failure_poisons_context() {
    auto context = cuda::Context::create(runtime::DeviceId{0});
    cuda::DeviceBuffer<double> buffer(context);
    buffer.resize(8, runtime::Contents::discard);
    // 4096 threads per block exceeds every CUDA device limit: a launch error.
    scale<<<1, 4096, 0, cuda::native_stream(context)>>>(buffer.view().data(), 8, 1.0);
    throws_code<ContractError>(ErrorCode::execution_failed,
                               [&] { cuda::check_launch(context, "scale"); });
    expect(context.has_failed(), "A launch error must poison the context");
    throws_code<ContractError>(ErrorCode::execution_failed, [&] { context.record_completion(); });
    throws_code<ContractError>(ErrorCode::execution_failed,
                               [&] { buffer.resize(16, runtime::Contents::discard); });
    // Buffers and context are released at scope exit without throwing.
}

}  // namespace

int main() {
    try {
        device_selection();
        auto context = cuda::Context::create(runtime::DeviceId{0});
        const auto& p = context.properties();
        std::cout << "device=" << p.name << " sm=" << p.major << p.minor
                  << " driver=" << p.driver_version << " runtime=" << p.runtime_version << '\n';
        round_trip_and_kernel(context);
        growth_invalidation_and_reuse(context);
        allocation_failure_is_recoverable(context);
        device_faults(context);
        launch_failure_poisons_context();
        // A fresh context after a failed one works.
        auto again = cuda::Context::create(runtime::DeviceId{0});
        expect(again.record_completion().wait().succeeded(), "Replacement context failed");
        std::cout << "CUDA runtime contracts passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
