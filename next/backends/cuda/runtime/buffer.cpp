#include "runtime/native.hpp"

#include <algorithm>
#include <string>
#include <utility>

namespace gmd_next::cuda {

using core::ErrorCode;

RawBuffer::RawBuffer(const Context& context, std::size_t element_bytes)
    : state_(context.state()), tracker_(element_bytes) {}

RawBuffer::RawBuffer(RawBuffer&& other) noexcept
    : state_(std::move(other.state_)), tracker_(other.tracker_), data_(other.data_) {
    other.state_.reset();
    other.data_ = nullptr;
    // The moved-from object keeps the identity but not the state, so any view
    // checked against it is stale rather than silently valid.
    other.tracker_.release();
}

RawBuffer& RawBuffer::operator=(RawBuffer&& other) noexcept {
    if (this != &other) {
        free_storage();
        state_ = std::move(other.state_);
        tracker_ = other.tracker_;
        data_ = other.data_;
        other.state_.reset();
        other.data_ = nullptr;
        other.tracker_.release();
    }
    return *this;
}

RawBuffer::~RawBuffer() { free_storage(); }

void RawBuffer::free_storage() noexcept {
    if (!state_ || !data_) return;
    // Stream-ordered: the free happens after all work already enqueued on this
    // context's stream, so pending kernels and copies keep valid memory.
    cudaSetDevice(state_->device.ordinal());
    if (cudaFreeAsync(data_, state_->stream) != cudaSuccess) {
        state_->fail({ErrorCode::execution_failed, "cudaFreeAsync", "buffer release failed"});
        cudaGetLastError();
    }
    data_ = nullptr;
    tracker_.release();
}

detail::ContextState& RawBuffer::live() const {
    if (!state_) {
        throw core::ContractError({ErrorCode::malformed_input, "buffer",
                                   "operation on a moved-from buffer"});
    }
    state_->ensure_usable();
    return *state_;
}

runtime::DeviceId RawBuffer::device() const {
    return state_ ? state_->device : runtime::DeviceId{};
}

void RawBuffer::resize(std::size_t count, runtime::Contents contents) {
    auto& state = live();
    const auto plan = tracker_.plan(count);
    if (!plan.reallocate) {
        tracker_.commit(plan);
        return;
    }
    state.activate();
    const std::size_t bytes = plan.capacity * tracker_.element_bytes();
    void* fresh = nullptr;
    const cudaError_t allocated = cudaMallocAsync(&fresh, bytes, state.stream);
    if (allocated == cudaErrorMemoryAllocation) {
        // Out of memory is recoverable: clear it and leave the buffer as it was.
        cudaGetLastError();
        throw core::ContractError({ErrorCode::allocation_failed, "buffer",
                                   std::to_string(bytes) + " bytes could not be allocated"});
    }
    state.check(allocated, "cudaMallocAsync");

    if (contents == runtime::Contents::preserve && data_ && tracker_.size() > 0) {
        const std::size_t kept = std::min(tracker_.size(), count) * tracker_.element_bytes();
        const cudaError_t copied =
            cudaMemcpyAsync(fresh, data_, kept, cudaMemcpyDeviceToDevice, state.stream);
        if (copied != cudaSuccess) {
            cudaFreeAsync(fresh, state.stream);
            state.check(copied, "cudaMemcpyAsync(preserve)");
        }
    }
    void* old = data_;
    data_ = fresh;
    tracker_.commit(plan);
    if (old) state.check(cudaFreeAsync(old, state.stream), "cudaFreeAsync");
}

void RawBuffer::copy_from_host(const void* host, std::size_t count) {
    auto& state = live();
    if (count != tracker_.size()) {
        throw core::ContractError({ErrorCode::inconsistent_size, "upload",
                                   std::to_string(count) + " host elements for a buffer of " +
                                       std::to_string(tracker_.size())});
    }
    if (count == 0) return;
    state.activate();
    state.check(cudaMemcpyAsync(data_, host, count * tracker_.element_bytes(),
                                cudaMemcpyHostToDevice, state.stream),
                "cudaMemcpyAsync(upload)");
}

void RawBuffer::copy_to_host(void* host, std::size_t count) const {
    auto& state = live();
    if (count != tracker_.size()) {
        throw core::ContractError({ErrorCode::inconsistent_size, "download",
                                   std::to_string(count) + " host elements for a buffer of " +
                                       std::to_string(tracker_.size())});
    }
    if (count == 0) return;
    state.activate();
    state.check(cudaMemcpyAsync(host, data_, count * tracker_.element_bytes(),
                                cudaMemcpyDeviceToHost, state.stream),
                "cudaMemcpyAsync(download)");
}

core::ValidationReport RawBuffer::check(runtime::BufferIdentity identity, std::size_t view_size,
                                        runtime::DeviceId view_device) const {
    return tracker_.check(identity, view_size, view_device, device());
}

}  // namespace gmd_next::cuda
