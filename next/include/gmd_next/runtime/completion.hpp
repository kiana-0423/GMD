#pragma once

#include "gmd_next/core/status.hpp"

#include <cstdint>
#include <optional>
#include <string_view>

namespace gmd_next::runtime {

// Submitting work is not finishing it. A completion is pending until the
// device reports, and succeeds only if neither the runtime nor the device
// status word reported a problem.
enum class CompletionState : std::uint8_t { pending, succeeded, failed };

std::string_view name(CompletionState state);

struct CompletionStatus {
    CompletionState state = CompletionState::pending;
    core::ValidationReport failures;  // empty unless state == failed

    bool succeeded() const { return state == CompletionState::succeeded; }
};

// Bits kernels may OR into the context's device status word. Each maps to one
// ErrorCode; bits outside kKnownDeviceFaults are reported, never ignored.
enum class DeviceFault : std::uint32_t {
    invalid_argument = 1u << 0,   // an index or parameter a kernel had to reject
    capacity_exceeded = 1u << 1,  // output needed more room than was reserved
    non_finite_result = 1u << 2,  // a computed value was NaN or infinite
};

inline constexpr std::uint32_t kKnownDeviceFaults = 0x7u;

constexpr std::uint32_t bit(DeviceFault fault) {
    return static_cast<std::uint32_t>(fault);
}

core::ValidationReport decode_device_faults(std::uint32_t word);

// The first backend reads the device status word back through one pinned slot,
// so only one completion may be outstanding per context. Tickets make a stale
// or duplicated resolution an explicit error instead of a lost fault.
class CompletionLedger {
public:
    // Throws resource_busy while a previous ticket is unresolved.
    std::uint64_t begin();
    bool is_outstanding(std::uint64_t ticket) const {
        return outstanding_ && *outstanding_ == ticket;
    }
    // Throws malformed_input for a ticket that is not the outstanding one.
    void resolve(std::uint64_t ticket);

private:
    std::uint64_t next_ = 1;
    std::optional<std::uint64_t> outstanding_;
};

}  // namespace gmd_next::runtime
