#pragma once

#include "gmd_next/core/numeric.hpp"

#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace gmd_next::core {

// Every rejection carries one of these codes. Adding a code is a contract
// change: callers switch on it, so no configuration may fail without one.
enum class ErrorCode : std::uint8_t {
    none,
    malformed_input,            // structurally wrong value, e.g. a negative count
    non_finite_value,           // NaN or infinity where a number is required
    value_out_of_range,         // finite but outside the admissible interval
    inconsistent_size,          // arrays that must describe the same atoms do not
    duplicate_identity,         // the same stable id appears twice
    unsupported_configuration,  // well-formed, but outside the implemented scope
    unit_mismatch,              // the unit system does not define the quantity
    version_mismatch,           // a recorded state version is stale
    capacity_overflow,          // the request exceeds an index or buffer capacity
    not_requested,              // the quantity was never asked for
    not_evaluated,              // requested, but no value was produced
    device_unavailable,         // no such device, or it cannot be used at all
    device_mismatch,            // a view or resource belongs to another device
    allocation_failed,          // the device could not provide the memory
    execution_failed,           // a launch, copy or completion reported an error
    stale_view,                 // the buffer was reallocated after the view was taken
    resource_busy,              // a single-use resource is still owned by pending work
};

std::string_view name(ErrorCode code);

struct Diagnostic {
    ErrorCode code = ErrorCode::none;
    std::string field;   // the input, option or observable that was rejected
    std::string detail;  // what was wrong with it
    std::string message() const;
};

class ContractError : public std::runtime_error {
public:
    explicit ContractError(Diagnostic diagnostic);
    ErrorCode code() const noexcept { return diagnostic_.code; }
    const Diagnostic& diagnostic() const noexcept { return diagnostic_; }

private:
    Diagnostic diagnostic_;
};

// Collects every problem found in one input instead of stopping at the first,
// so a caller can report a whole malformed configuration at once.
class ValidationReport {
public:
    void add(ErrorCode code, std::string field, std::string detail);
    void merge(const ValidationReport& other);
    bool ok() const { return diagnostics_.empty(); }
    bool contains(ErrorCode code) const;
    const Diagnostic* find(ErrorCode code) const;
    std::span<const Diagnostic> diagnostics() const { return diagnostics_; }
    std::string summary() const;
    // Throws ContractError built from the first diagnostic; no-op when ok().
    void require_ok() const;

private:
    std::vector<Diagnostic> diagnostics_;
};

// A result is valid, unavailable or invalid, and never silently a zero.
//   valid       a finite value was produced for this state
//   unavailable nothing was computed: not requested, or not sampled here
//   invalid     computation was attempted and the result must not be consumed
enum class Availability : std::uint8_t { valid, unavailable, invalid };

std::string_view name(Availability availability);

// A payload with no host-side numbers, used for device-resident results whose
// presence and reason are the only parts a host record carries.
inline bool all_finite(std::monostate) { return true; }

// Carrier for a physical quantity plus its availability. Unavailable and
// invalid quantities keep the reason so output records can state why a column
// is missing rather than writing a plausible number.
template<class T>
class Quantity {
public:
    Quantity() : error_(ErrorCode::not_requested), reason_("not requested") {}

    static Quantity from_value(T value) {
        Quantity quantity;
        if (!all_finite(value)) {
            return failed(ErrorCode::non_finite_value, "value is not finite");
        }
        quantity.value_ = std::move(value);
        quantity.availability_ = Availability::valid;
        quantity.error_ = ErrorCode::none;
        quantity.reason_.clear();
        return quantity;
    }

    static Quantity not_requested() { return Quantity{}; }

    // Requested but not produced here, for example sampled at another stage.
    static Quantity missing(std::string reason) {
        Quantity quantity;
        quantity.error_ = ErrorCode::not_evaluated;
        quantity.reason_ = std::move(reason);
        return quantity;
    }

    static Quantity failed(ErrorCode code, std::string reason) {
        Quantity quantity;
        quantity.availability_ = Availability::invalid;
        quantity.error_ = code;
        quantity.reason_ = std::move(reason);
        return quantity;
    }

    Availability availability() const { return availability_; }
    bool is_valid() const { return availability_ == Availability::valid; }
    ErrorCode error() const { return error_; }
    std::string_view reason() const { return reason_; }

    // Throws instead of handing out the default-constructed payload.
    const T& value() const {
        if (availability_ != Availability::valid) {
            throw ContractError({error_, "quantity", reason_});
        }
        return value_;
    }

    T value_or(T fallback) const { return is_valid() ? value_ : fallback; }

private:
    T value_{};
    Availability availability_ = Availability::unavailable;
    ErrorCode error_ = ErrorCode::not_requested;
    std::string reason_;
};

using Presence = Quantity<std::monostate>;
// Availability of a payload that is not held on the host, such as a device
// force array. Only the presence and its reason travel in host records.
inline Presence present() { return Presence::from_value(std::monostate{}); }

// File and log convention: a missing scalar is written as NaN, because zero is
// a valid physical value and must stay distinguishable from "not sampled".
double value_or_nan(const Quantity<double>& quantity);

}  // namespace gmd_next::core
