#pragma once

#include "gmd_next/core/status.hpp"
#include "test_support.hpp"

#include <string>
#include <string_view>

// Checks that a report rejected exactly the named field with the named code.
inline void rejects(const gmd_next::core::ValidationReport& report,
                    gmd_next::core::ErrorCode code, std::string_view field) {
    for (const auto& diagnostic : report.diagnostics()) {
        if (diagnostic.code == code && diagnostic.field == field) return;
    }
    throw std::runtime_error("Expected " + std::string(name(code)) + " on " + std::string(field) +
                             ", received: " + report.summary());
}

inline void accepts(const gmd_next::core::ValidationReport& report) {
    expect(report.ok(), "Expected a valid configuration, received: " + report.summary());
}
