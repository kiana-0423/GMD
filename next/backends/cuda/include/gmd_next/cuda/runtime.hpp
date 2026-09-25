#pragma once

#include "gmd_next/core/status.hpp"
#include "gmd_next/runtime/completion.hpp"
#include "gmd_next/runtime/device_view.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>

// CUDA execution resources. No CUDA type appears here; kernels inside the
// backend reach the stream through runtime/native.hpp. Every class is bound to
// one Context and is not thread-safe: one host thread drives one context.
namespace gmd_next::cuda {

namespace detail {
struct ContextState;  // defined in backends/cuda/runtime/native.hpp
}

struct DeviceProperties {
    std::string name;
    int major = 0;
    int minor = 0;
    std::size_t global_memory_bytes = 0;
    int driver_version = 0;
    int runtime_version = 0;
};

class Completion;

// One device, one non-blocking stream, one completion event and one device
// status word that kernels OR DeviceFault bits into. All work is ordered on
// that single stream, which is what makes stream-ordered allocation and free
// safe for buffers used only through this context.
class Context {
public:
    // Throws ContractError:
    //   device_unavailable         no driver/device, or ordinal out of range
    //   unsupported_configuration  below sm_60 (no native FP64 atomicAdd) or
    //                              no stream-ordered memory pools
    static Context create(runtime::DeviceId device);

    Context(Context&&) noexcept;
    Context& operator=(Context&&) noexcept;
    ~Context();

    runtime::DeviceId device() const;
    const DeviceProperties& properties() const;

    // A runtime error (not a device fault) poisons the context: every later
    // submission throws execution_failed. Recovery means a new context.
    bool has_failed() const;
    void ensure_usable() const;

    // Enqueues the status word read-back and reset, then the event. Only one
    // completion may be outstanding (resource_busy otherwise).
    Completion record_completion();

    runtime::DeviceSpan<std::uint32_t> fault_word() const;
    core::ValidationReport check_device(runtime::DeviceId view_device, std::string_view field) const;

    const std::shared_ptr<detail::ContextState>& state() const { return state_; }

private:
    explicit Context(std::shared_ptr<detail::ContextState> state);
    std::shared_ptr<detail::ContextState> state_;
};

// Outcome of all work enqueued on the context before record_completion().
// Device faults fail the completion but leave the context usable, so a
// capacity overflow can be retried; runtime errors also poison the context.
class Completion {
public:
    Completion(Completion&&) noexcept;
    Completion& operator=(Completion&&) noexcept;
    // Waits if still pending, so a fault can never be dropped unobserved.
    ~Completion();

    const runtime::CompletionStatus& query();
    const runtime::CompletionStatus& wait();
    const runtime::CompletionStatus& status() const { return status_; }

private:
    friend class Context;
    Completion(std::shared_ptr<detail::ContextState> state, std::uint64_t ticket);
    void resolve(int cuda_error);
    void abandon() noexcept;

    std::shared_ptr<detail::ContextState> state_;
    std::uint64_t ticket_ = 0;
    runtime::CompletionStatus status_;
};

// Untyped, move-only device allocation with stream-ordered alloc/free.
// Capacity is reserved geometrically and never shrinks; resizing within the
// capacity allocates nothing. Any resize that changes the size or
// reallocates makes earlier views stale.
class RawBuffer {
public:
    RawBuffer(const Context& context, std::size_t element_bytes);
    RawBuffer(RawBuffer&& other) noexcept;
    RawBuffer& operator=(RawBuffer&& other) noexcept;
    RawBuffer(const RawBuffer&) = delete;
    RawBuffer& operator=(const RawBuffer&) = delete;
    ~RawBuffer();

    // Throws allocation_failed and leaves the buffer unchanged when the device
    // is out of memory; the context stays usable in that case.
    void resize(std::size_t count, runtime::Contents contents);

    std::size_t size() const { return tracker_.size(); }
    std::size_t capacity() const { return tracker_.capacity(); }
    runtime::BufferIdentity identity() const { return tracker_.identity(); }
    runtime::DeviceId device() const;
    void* data() const { return data_; }

    // Stream-ordered copies of exactly size() elements. The host memory must
    // stay valid, and for downloads unread, until a later completion succeeds.
    void copy_from_host(const void* host, std::size_t count);
    void copy_to_host(void* host, std::size_t count) const;

    core::ValidationReport check(runtime::BufferIdentity identity, std::size_t view_size,
                                 runtime::DeviceId view_device) const;

private:
    detail::ContextState& live() const;
    void free_storage() noexcept;

    std::shared_ptr<detail::ContextState> state_;
    runtime::CapacityTracker tracker_;
    void* data_ = nullptr;
};

template<class T>
class DeviceBuffer {
    static_assert(std::is_trivially_copyable_v<T>, "device buffers hold trivially copyable data");

public:
    explicit DeviceBuffer(const Context& context) : raw_(context, sizeof(T)) {}

    void resize(std::size_t count, runtime::Contents contents) { raw_.resize(count, contents); }
    std::size_t size() const { return raw_.size(); }
    std::size_t capacity() const { return raw_.capacity(); }

    runtime::DeviceSpan<T> view() {
        return {static_cast<T*>(raw_.data()), raw_.size(), raw_.device(), raw_.identity()};
    }
    runtime::DeviceSpan<const T> view() const {
        return {static_cast<const T*>(raw_.data()), raw_.size(), raw_.device(), raw_.identity()};
    }
    core::ValidationReport check(const runtime::DeviceSpan<const T>& view) const {
        return raw_.check(view.identity(), view.size(), view.device());
    }
    core::ValidationReport check(const runtime::DeviceSpan<T>& view) const {
        return raw_.check(view.identity(), view.size(), view.device());
    }

    void upload(std::span<const T> host) { raw_.copy_from_host(host.data(), host.size()); }
    void download(std::span<T> host) const { raw_.copy_to_host(host.data(), host.size()); }

private:
    RawBuffer raw_;
};

// Scratch for library calls such as CUB. Contents are never preserved and the
// size only grows, so a size query followed by reserve() allocates at most
// once per new maximum. A workspace is ordered by its context's stream; using
// it from any other stream needs an explicit event dependency (not provided).
class Workspace {
public:
    explicit Workspace(const Context& context) : raw_(context, 1) {}

    void reserve(std::size_t bytes) {
        if (bytes > raw_.size()) raw_.resize(bytes, runtime::Contents::discard);
    }
    std::size_t size_bytes() const { return raw_.size(); }
    runtime::DeviceSpan<std::byte> view() {
        return {static_cast<std::byte*>(raw_.data()), raw_.size(), raw_.device(), raw_.identity()};
    }
    core::ValidationReport check(const runtime::DeviceSpan<std::byte>& view) const {
        return raw_.check(view.identity(), view.size(), view.device());
    }

private:
    RawBuffer raw_;
};

}  // namespace gmd_next::cuda
