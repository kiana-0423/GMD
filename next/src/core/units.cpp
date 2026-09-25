#include "gmd_next/core/units.hpp"

#include "gmd_next/core/numeric.hpp"

#include <string>

namespace gmd_next::core {
namespace {

void require_finite(double value, std::string_view field) {
    if (!all_finite(value)) {
        throw ContractError({ErrorCode::non_finite_value, std::string(field),
                             "conversion input is not finite"});
    }
}

// reduced units carry no scale to seconds, pascal or kelvin, so a conversion
// into a physical unit has no defined answer and is rejected instead.
void require_metal(UnitSystem system, std::string_view field, std::string_view quantity) {
    if (system != UnitSystem::metal) {
        throw ContractError({ErrorCode::unit_mismatch, std::string(field),
                             std::string(name(system)) + " units do not define " +
                                 std::string(quantity)});
    }
}

}  // namespace

std::string_view name(UnitSystem system) {
    switch (system) {
    case UnitSystem::metal: return "metal";
    case UnitSystem::reduced: return "reduced";
    }
    return "unknown";
}

std::string_view name(TimeUnit unit) {
    switch (unit) {
    case TimeUnit::femtoseconds: return "femtoseconds";
    case TimeUnit::internal: return "internal";
    }
    return "unknown";
}

double to_internal_time(double value, TimeUnit unit, UnitSystem system) {
    require_finite(value, "time");
    switch (unit) {
    case TimeUnit::internal:
        return value;
    case TimeUnit::femtoseconds:
        require_metal(system, "time", "femtoseconds");
        return value * units::kInternalTimePerFemtosecond;
    }
    throw ContractError({ErrorCode::malformed_input, "time", "unknown time unit"});
}

double from_internal_time(double internal_value, TimeUnit unit, UnitSystem system) {
    require_finite(internal_value, "time");
    switch (unit) {
    case TimeUnit::internal:
        return internal_value;
    case TimeUnit::femtoseconds:
        require_metal(system, "time", "femtoseconds");
        return internal_value * units::kFemtosecondsPerInternalTime;
    }
    throw ContractError({ErrorCode::malformed_input, "time", "unknown time unit"});
}

double pressure_to_bar(double ev_per_cubic_angstrom, UnitSystem system) {
    require_finite(ev_per_cubic_angstrom, "pressure");
    require_metal(system, "pressure", "bar");
    return ev_per_cubic_angstrom * units::kEvPerCubicAngstromToBar;
}

double pressure_from_bar(double bar, UnitSystem system) {
    require_finite(bar, "pressure");
    require_metal(system, "pressure", "bar");
    return bar * units::kBarToEvPerCubicAngstrom;
}

ValidationReport validate_step_stamp(const StepStamp& stamp, std::string_view field) {
    ValidationReport report;
    const std::string prefix(field);
    if (stamp.step < 0) {
        report.add(ErrorCode::value_out_of_range, prefix + ".step",
                   "step index must be nonnegative, received " + std::to_string(stamp.step));
    }
    if (!all_finite(stamp.time_internal)) {
        report.add(ErrorCode::non_finite_value, prefix + ".time", "time must be finite");
    } else if (stamp.time_internal < 0.0) {
        report.add(ErrorCode::value_out_of_range, prefix + ".time",
                   "time must be nonnegative, received " + std::to_string(stamp.time_internal));
    }
    return report;
}

}  // namespace gmd_next::core
