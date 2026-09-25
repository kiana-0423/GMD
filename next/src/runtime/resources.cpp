#include "gmd_next/runtime/completion.hpp"
#include "gmd_next/runtime/device_view.hpp"

#include <atomic>
#include <limits>
#include <string>

namespace gmd_next::runtime {
namespace {

using core::ErrorCode;

constexpr std::size_t kMaxSize = std::numeric_limits<std::size_t>::max();

std::uint64_t next_buffer_id() {
    static std::atomic<std::uint64_t> counter{0};
    return ++counter;
}

}  // namespace

std::size_t grown_capacity(std::size_t current, std::size_t required, std::size_t element_bytes) {
    if (element_bytes == 0) {
        throw core::ContractError({ErrorCode::malformed_input, "element_bytes",
                                   "element size must be positive"});
    }
    const std::size_t max_elements = kMaxSize / element_bytes;
    if (required > max_elements) {
        throw core::ContractError({ErrorCode::capacity_overflow, "capacity",
                                   std::to_string(required) + " elements of " +
                                       std::to_string(element_bytes) +
                                       " bytes exceed the addressable size"});
    }
    if (required <= current) return current;
    // current + current/2 cannot overflow unless current > 2/3 of max; clamp.
    const std::size_t geometric =
        current > max_elements - current / 2 ? max_elements : current + current / 2;
    return required > geometric ? required : geometric;
}

CapacityTracker::CapacityTracker(std::size_t element_bytes)
    : element_bytes_(element_bytes), buffer_(next_buffer_id()) {
    if (element_bytes == 0) {
        throw core::ContractError({ErrorCode::malformed_input, "element_bytes",
                                   "element size must be positive"});
    }
}

ResizePlan CapacityTracker::plan(std::size_t count) const {
    const std::size_t capacity = grown_capacity(capacity_, count, element_bytes_);
    return ResizePlan{capacity != capacity_, capacity, count};
}

void CapacityTracker::commit(const ResizePlan& plan) {
    if (plan.size > plan.capacity || (!plan.reallocate && plan.capacity != capacity_)) {
        throw core::ContractError({ErrorCode::malformed_input, "resize_plan",
                                   "plan does not describe this buffer"});
    }
    if (plan.reallocate || plan.size != size_) ++generation_;
    capacity_ = plan.capacity;
    size_ = plan.size;
}

void CapacityTracker::release() {
    size_ = 0;
    capacity_ = 0;
    ++generation_;
}

core::ValidationReport CapacityTracker::check(BufferIdentity identity, std::size_t view_size,
                                              DeviceId view_device, DeviceId owner_device) const {
    core::ValidationReport report = check_same_device(owner_device, view_device, "view.device");
    if (identity.buffer != buffer_) {
        report.add(ErrorCode::stale_view, "view.buffer", "view was taken from another buffer");
    } else if (identity.generation != generation_) {
        report.add(ErrorCode::stale_view, "view.generation",
                   "buffer was resized or reallocated after the view was taken (view " +
                       std::to_string(identity.generation) + ", buffer " +
                       std::to_string(generation_) + ")");
    }
    if (view_size != size_) {
        report.add(ErrorCode::inconsistent_size, "view.size",
                   "view holds " + std::to_string(view_size) + " elements, buffer " +
                       std::to_string(size_));
    }
    return report;
}

core::ValidationReport check_same_device(DeviceId resource, DeviceId view, std::string_view field) {
    core::ValidationReport report;
    if (!resource.is_valid() || !view.is_valid() || resource != view) {
        report.add(ErrorCode::device_mismatch, std::string(field),
                   "device " + std::to_string(view.ordinal()) + " used with a resource on device " +
                       std::to_string(resource.ordinal()));
    }
    return report;
}

std::string_view name(CompletionState state) {
    switch (state) {
    case CompletionState::pending: return "pending";
    case CompletionState::succeeded: return "succeeded";
    case CompletionState::failed: return "failed";
    }
    return "unknown";
}

core::ValidationReport decode_device_faults(std::uint32_t word) {
    core::ValidationReport report;
    if (word & bit(DeviceFault::invalid_argument)) {
        report.add(ErrorCode::malformed_input, "device_status", "a kernel rejected its input");
    }
    if (word & bit(DeviceFault::capacity_exceeded)) {
        report.add(ErrorCode::capacity_overflow, "device_status",
                   "a kernel needed more output capacity than was reserved");
    }
    if (word & bit(DeviceFault::non_finite_result)) {
        report.add(ErrorCode::non_finite_value, "device_status",
                   "a kernel produced a non-finite value");
    }
    if (const auto unknown = word & ~kKnownDeviceFaults) {
        report.add(ErrorCode::execution_failed, "device_status",
                   "unknown device fault bits " + std::to_string(unknown));
    }
    return report;
}

std::uint64_t CompletionLedger::begin() {
    if (outstanding_) {
        throw core::ContractError({ErrorCode::resource_busy, "completion",
                                   "ticket " + std::to_string(*outstanding_) +
                                       " must be resolved before new work is recorded"});
    }
    outstanding_ = next_++;
    return *outstanding_;
}

void CompletionLedger::resolve(std::uint64_t ticket) {
    if (!is_outstanding(ticket)) {
        throw core::ContractError({ErrorCode::malformed_input, "completion",
                                   "ticket " + std::to_string(ticket) + " is not outstanding"});
    }
    outstanding_.reset();
}

}  // namespace gmd_next::runtime
