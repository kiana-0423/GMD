#pragma once

#include "gmd_next/core/status.hpp"

#include <compare>
#include <cstddef>
#include <cstdint>

namespace gmd_next::runtime {

// A device ordinal. No CUDA type appears in this header, so host code and the
// contracts target can reason about device ownership without a toolkit.
class DeviceId {
public:
    constexpr DeviceId() = default;
    constexpr explicit DeviceId(int ordinal) : ordinal_(ordinal) {}
    constexpr int ordinal() const { return ordinal_; }
    constexpr bool is_valid() const { return ordinal_ >= 0; }
    friend constexpr auto operator<=>(const DeviceId&, const DeviceId&) = default;

private:
    int ordinal_ = -1;
};

// Which buffer a view was taken from, and in which of its states. Any
// reallocation or change of logical size moves the generation, which ends
// every view taken before it.
struct BufferIdentity {
    std::uint64_t buffer = 0;      // 0 means "no buffer"
    std::uint64_t generation = 0;
    friend constexpr bool operator==(const BufferIdentity&, const BufferIdentity&) = default;
};

// Non-owning view of device memory. data() is a device address: it may be
// passed to a kernel or a stream-ordered copy, never dereferenced on the host,
// which is why this type has no element access and is not a std::span.
template<class T>
class DeviceSpan {
public:
    constexpr DeviceSpan() = default;
    constexpr DeviceSpan(T* data, std::size_t size, DeviceId device, BufferIdentity identity)
        : data_(data), size_(size), device_(device), identity_(identity) {}

    constexpr T* data() const { return data_; }
    constexpr std::size_t size() const { return size_; }
    constexpr std::size_t size_bytes() const { return size_ * sizeof(T); }
    constexpr bool empty() const { return size_ == 0; }
    constexpr DeviceId device() const { return device_; }
    constexpr BufferIdentity identity() const { return identity_; }

private:
    T* data_ = nullptr;
    std::size_t size_ = 0;
    DeviceId device_{};
    BufferIdentity identity_{};
};

// What happens to existing elements when a buffer must be reallocated.
//   discard   old contents are undefined afterwards; no device copy is made
//   preserve  the first min(old size, new size) elements are copied over
enum class Contents : std::uint8_t { discard, preserve };

// Geometric growth (x1.5, at least the request) so that a slowly growing
// neighbour count does not reallocate every rebuild. Throws capacity_overflow
// when the byte count would not fit in size_t.
std::size_t grown_capacity(std::size_t current, std::size_t required, std::size_t element_bytes);

// Outcome of planning a resize, applied by the owner only after the device
// work it implies has been enqueued successfully.
struct ResizePlan {
    bool reallocate = false;
    std::size_t capacity = 0;
    std::size_t size = 0;
};

// Host-side bookkeeping of one device buffer: logical size, capacity and the
// generation that views are checked against. It never touches device memory.
class CapacityTracker {
public:
    explicit CapacityTracker(std::size_t element_bytes);

    std::size_t size() const { return size_; }
    std::size_t capacity() const { return capacity_; }
    std::size_t element_bytes() const { return element_bytes_; }
    BufferIdentity identity() const { return {buffer_, generation_}; }

    // Reallocation only when the capacity is too small; shrinking keeps it.
    ResizePlan plan(std::size_t count) const;
    void commit(const ResizePlan& plan);
    // Storage returned: size and capacity drop to zero and old views go stale.
    void release();

    // A view is usable only for the exact buffer state it was taken from.
    core::ValidationReport check(BufferIdentity identity, std::size_t view_size,
                                 DeviceId view_device, DeviceId owner_device) const;

private:
    std::size_t element_bytes_;
    std::size_t size_ = 0;
    std::size_t capacity_ = 0;
    std::uint64_t buffer_;
    std::uint64_t generation_ = 1;
};

// Resources bound to one device must not receive views of another.
core::ValidationReport check_same_device(DeviceId resource, DeviceId view, std::string_view field);

}  // namespace gmd_next::runtime
