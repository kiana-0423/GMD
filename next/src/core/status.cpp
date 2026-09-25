#include "gmd_next/core/status.hpp"

#include <algorithm>
#include <limits>

namespace gmd_next::core {

std::string_view name(ErrorCode code) {
    switch (code) {
    case ErrorCode::none: return "none";
    case ErrorCode::malformed_input: return "malformed_input";
    case ErrorCode::non_finite_value: return "non_finite_value";
    case ErrorCode::value_out_of_range: return "value_out_of_range";
    case ErrorCode::inconsistent_size: return "inconsistent_size";
    case ErrorCode::duplicate_identity: return "duplicate_identity";
    case ErrorCode::unsupported_configuration: return "unsupported_configuration";
    case ErrorCode::unit_mismatch: return "unit_mismatch";
    case ErrorCode::version_mismatch: return "version_mismatch";
    case ErrorCode::capacity_overflow: return "capacity_overflow";
    case ErrorCode::not_requested: return "not_requested";
    case ErrorCode::not_evaluated: return "not_evaluated";
    case ErrorCode::device_unavailable: return "device_unavailable";
    case ErrorCode::device_mismatch: return "device_mismatch";
    case ErrorCode::allocation_failed: return "allocation_failed";
    case ErrorCode::execution_failed: return "execution_failed";
    case ErrorCode::stale_view: return "stale_view";
    case ErrorCode::resource_busy: return "resource_busy";
    }
    return "unknown";
}

std::string_view name(Availability availability) {
    switch (availability) {
    case Availability::valid: return "valid";
    case Availability::unavailable: return "unavailable";
    case Availability::invalid: return "invalid";
    }
    return "unknown";
}

std::string Diagnostic::message() const {
    std::string text(name(code));
    text += " [";
    text += field;
    text += "]: ";
    text += detail;
    return text;
}

ContractError::ContractError(Diagnostic diagnostic)
    : std::runtime_error(diagnostic.message()), diagnostic_(std::move(diagnostic)) {}

void ValidationReport::add(ErrorCode code, std::string field, std::string detail) {
    diagnostics_.push_back(Diagnostic{code, std::move(field), std::move(detail)});
}

void ValidationReport::merge(const ValidationReport& other) {
    diagnostics_.insert(diagnostics_.end(), other.diagnostics_.begin(), other.diagnostics_.end());
}

bool ValidationReport::contains(ErrorCode code) const {
    return find(code) != nullptr;
}

const Diagnostic* ValidationReport::find(ErrorCode code) const {
    const auto found = std::find_if(diagnostics_.begin(), diagnostics_.end(),
                                    [code](const Diagnostic& d) { return d.code == code; });
    return found == diagnostics_.end() ? nullptr : &*found;
}

std::string ValidationReport::summary() const {
    if (diagnostics_.empty()) return "no problems";
    std::string text;
    for (const auto& diagnostic : diagnostics_) {
        if (!text.empty()) text += "; ";
        text += diagnostic.message();
    }
    return text;
}

void ValidationReport::require_ok() const {
    if (diagnostics_.empty()) return;
    Diagnostic first = diagnostics_.front();
    if (diagnostics_.size() > 1) {
        first.detail += " (and " + std::to_string(diagnostics_.size() - 1) + " more)";
    }
    throw ContractError(std::move(first));
}

double value_or_nan(const Quantity<double>& quantity) {
    return quantity.value_or(std::numeric_limits<double>::quiet_NaN());
}

}  // namespace gmd_next::core
