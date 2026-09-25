#pragma once

// Backend-internal: the only header that exposes CUDA types. Include it from
// backend .cpp/.cu files, never from public headers.

#include "gmd_next/cuda/runtime.hpp"

#include <cuda_runtime.h>

#include <optional>

namespace gmd_next::cuda::detail {

struct ContextState {
    runtime::DeviceId device;
    DeviceProperties properties;
    cudaStream_t stream = nullptr;
    cudaEvent_t event = nullptr;
    std::uint32_t* device_status = nullptr;  // written by kernels
    std::uint32_t* host_status = nullptr;    // pinned read-back slot
    runtime::CompletionLedger ledger;
    std::optional<core::Diagnostic> failure;

    ContextState() = default;
    ContextState(const ContextState&) = delete;
    ContextState& operator=(const ContextState&) = delete;
    ~ContextState();

    // The runtime API is stateful per host thread; every entry point selects
    // this context's device before touching it.
    void activate();
    // The first runtime error sticks; later ones do not overwrite it.
    void fail(core::Diagnostic diagnostic);
    void ensure_usable() const;
    // Throws execution_failed for any error and poisons the context.
    void check(cudaError_t status, const char* operation);
};

}  // namespace gmd_next::cuda::detail

namespace gmd_next::cuda {

cudaStream_t native_stream(const Context& context);
// Call right after a kernel launch: turns a launch error into execution_failed.
void check_launch(const Context& context, const char* kernel);

}  // namespace gmd_next::cuda
