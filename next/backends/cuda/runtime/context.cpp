#include "runtime/native.hpp"

#include <string>
#include <utility>

namespace gmd_next::cuda {
namespace detail {

using core::ErrorCode;

namespace {

core::Diagnostic runtime_error(cudaError_t status, const char* operation) {
    return {ErrorCode::execution_failed, operation, cudaGetErrorString(status)};
}

}  // namespace

ContextState::~ContextState() {
    // Failure cleanup: best effort, never throws. Synchronizing first keeps
    // stream-ordered frees and pending copies from outliving their memory.
    if (device.is_valid()) cudaSetDevice(device.ordinal());
    if (stream) cudaStreamSynchronize(stream);
    if (host_status) cudaFreeHost(host_status);
    if (device_status) cudaFree(device_status);
    if (event) cudaEventDestroy(event);
    if (stream) cudaStreamDestroy(stream);
}

void ContextState::activate() {
    const cudaError_t status = cudaSetDevice(device.ordinal());
    if (status != cudaSuccess) check(status, "cudaSetDevice");
}

void ContextState::fail(core::Diagnostic diagnostic) {
    if (!failure) failure = std::move(diagnostic);
}

void ContextState::ensure_usable() const {
    if (failure) {
        throw core::ContractError({ErrorCode::execution_failed, "context",
                                   "context failed earlier: " + failure->message()});
    }
}

void ContextState::check(cudaError_t status, const char* operation) {
    if (status == cudaSuccess) return;
    auto diagnostic = runtime_error(status, operation);
    fail(diagnostic);
    throw core::ContractError(std::move(diagnostic));
}

}  // namespace detail

using core::ErrorCode;

Context Context::create(runtime::DeviceId device) {
    int count = 0;
    const cudaError_t counted = cudaGetDeviceCount(&count);
    if (counted != cudaSuccess) {
        cudaGetLastError();
        throw core::ContractError({ErrorCode::device_unavailable, "device",
                                   std::string("no usable CUDA device: ") +
                                       cudaGetErrorString(counted)});
    }
    if (!device.is_valid() || device.ordinal() >= count) {
        throw core::ContractError({ErrorCode::device_unavailable, "device",
                                   "ordinal " + std::to_string(device.ordinal()) + " of " +
                                       std::to_string(count) + " visible devices"});
    }

    // Constructed first so that any failure below releases what was created.
    auto state = std::make_shared<detail::ContextState>();
    state->device = device;
    state->activate();

    cudaDeviceProp properties{};
    state->check(cudaGetDeviceProperties(&properties, device.ordinal()), "cudaGetDeviceProperties");
    state->properties.name = properties.name;
    state->properties.major = properties.major;
    state->properties.minor = properties.minor;
    state->properties.global_memory_bytes = properties.totalGlobalMem;
    state->check(cudaDriverGetVersion(&state->properties.driver_version), "cudaDriverGetVersion");
    state->check(cudaRuntimeGetVersion(&state->properties.runtime_version), "cudaRuntimeGetVersion");

    if (properties.major < 6) {
        throw core::ContractError({ErrorCode::unsupported_configuration, "device",
                                   "sm_" + std::to_string(properties.major) +
                                       std::to_string(properties.minor) +
                                       " lacks native FP64 atomicAdd; sm_60 or newer is required"});
    }
    int pools = 0;
    state->check(cudaDeviceGetAttribute(&pools, cudaDevAttrMemoryPoolsSupported, device.ordinal()),
                 "cudaDeviceGetAttribute");
    if (pools == 0) {
        throw core::ContractError({ErrorCode::unsupported_configuration, "device",
                                   "stream-ordered memory pools are not supported"});
    }

    state->check(cudaStreamCreateWithFlags(&state->stream, cudaStreamNonBlocking),
                 "cudaStreamCreateWithFlags");
    state->check(cudaEventCreateWithFlags(&state->event, cudaEventDisableTiming),
                 "cudaEventCreateWithFlags");
    state->check(cudaMalloc(reinterpret_cast<void**>(&state->device_status), sizeof(std::uint32_t)),
                 "cudaMalloc(status)");
    state->check(cudaMallocHost(reinterpret_cast<void**>(&state->host_status), sizeof(std::uint32_t)),
                 "cudaMallocHost(status)");
    *state->host_status = 0;
    state->check(cudaMemsetAsync(state->device_status, 0, sizeof(std::uint32_t), state->stream),
                 "cudaMemsetAsync(status)");
    state->check(cudaStreamSynchronize(state->stream), "cudaStreamSynchronize");
    return Context{std::move(state)};
}

Context::Context(std::shared_ptr<detail::ContextState> state) : state_(std::move(state)) {}
Context::Context(Context&&) noexcept = default;
Context& Context::operator=(Context&&) noexcept = default;
Context::~Context() = default;

runtime::DeviceId Context::device() const { return state_->device; }
const DeviceProperties& Context::properties() const { return state_->properties; }
bool Context::has_failed() const { return state_->failure.has_value(); }
void Context::ensure_usable() const { state_->ensure_usable(); }

runtime::DeviceSpan<std::uint32_t> Context::fault_word() const {
    return {state_->device_status, 1, state_->device, runtime::BufferIdentity{}};
}

core::ValidationReport Context::check_device(runtime::DeviceId view_device,
                                             std::string_view field) const {
    return runtime::check_same_device(state_->device, view_device, field);
}

Completion Context::record_completion() {
    auto& state = *state_;
    state.ensure_usable();
    state.activate();
    const auto ticket = state.ledger.begin();
    try {
        // Read back, then clear, so each completion reports only the faults
        // raised by work enqueued since the previous one.
        state.check(cudaMemcpyAsync(state.host_status, state.device_status, sizeof(std::uint32_t),
                                    cudaMemcpyDeviceToHost, state.stream),
                    "cudaMemcpyAsync(status)");
        state.check(cudaMemsetAsync(state.device_status, 0, sizeof(std::uint32_t), state.stream),
                    "cudaMemsetAsync(status)");
        state.check(cudaEventRecord(state.event, state.stream), "cudaEventRecord");
    } catch (...) {
        state.ledger.resolve(ticket);
        throw;
    }
    return Completion{state_, ticket};
}

Completion::Completion(std::shared_ptr<detail::ContextState> state, std::uint64_t ticket)
    : state_(std::move(state)), ticket_(ticket) {}

Completion::Completion(Completion&& other) noexcept
    : state_(std::move(other.state_)), ticket_(other.ticket_), status_(std::move(other.status_)) {
    other.state_.reset();
}

Completion& Completion::operator=(Completion&& other) noexcept {
    if (this != &other) {
        abandon();
        state_ = std::move(other.state_);
        ticket_ = other.ticket_;
        status_ = std::move(other.status_);
        other.state_.reset();
    }
    return *this;
}

Completion::~Completion() { abandon(); }

void Completion::abandon() noexcept {
    if (!state_ || status_.state != runtime::CompletionState::pending) return;
    try {
        wait();
    } catch (...) {
        // wait() already recorded any runtime failure on the context.
    }
}

const runtime::CompletionStatus& Completion::query() {
    if (!state_ || status_.state != runtime::CompletionState::pending) return status_;
    state_->activate();
    const cudaError_t status = cudaEventQuery(state_->event);
    if (status == cudaErrorNotReady) return status_;
    resolve(static_cast<int>(status));
    return status_;
}

const runtime::CompletionStatus& Completion::wait() {
    if (!state_ || status_.state != runtime::CompletionState::pending) return status_;
    state_->activate();
    resolve(static_cast<int>(cudaEventSynchronize(state_->event)));
    return status_;
}

void Completion::resolve(int cuda_error) {
    auto& state = *state_;
    state.ledger.resolve(ticket_);
    const auto status = static_cast<cudaError_t>(cuda_error);
    if (status != cudaSuccess) {
        core::Diagnostic diagnostic{ErrorCode::execution_failed, "completion",
                                    cudaGetErrorString(status)};
        state.fail(diagnostic);
        status_.failures.add(diagnostic.code, diagnostic.field, diagnostic.detail);
        status_.state = runtime::CompletionState::failed;
        return;
    }
    // The pinned slot is complete once the event is: it was copied before it.
    status_.failures = runtime::decode_device_faults(*state.host_status);
    status_.state = status_.failures.ok() ? runtime::CompletionState::succeeded
                                          : runtime::CompletionState::failed;
}

cudaStream_t native_stream(const Context& context) { return context.state()->stream; }

void check_launch(const Context& context, const char* kernel) {
    context.state()->check(cudaGetLastError(), kernel);
}

}  // namespace gmd_next::cuda
