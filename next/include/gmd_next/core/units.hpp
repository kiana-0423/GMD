#pragma once

#include "gmd_next/core/status.hpp"

#include <cstdint>
#include <string_view>

namespace gmd_next::core {

// metal:   length Angstrom, energy eV, mass amu, charge e, temperature K.
//          Integrating v += (F/m) dt with F in eV/A and m in amu fixes the
//          internal time unit to A * sqrt(amu/eV); fs is an interface unit.
// reduced: dimensionless Lennard-Jones units. No scale relates them to fs, bar
//          or K, so those conversions are rejected rather than assumed.
enum class UnitSystem : std::uint8_t { metal, reduced };
enum class TimeUnit : std::uint8_t { femtoseconds, internal };

std::string_view name(UnitSystem system);
std::string_view name(TimeUnit unit);

namespace units {

// Internal time unit in femtoseconds, T = 1e-10 * sqrt(m_u / e) / 1e-15 fs with
//   e   = 1.602176634e-19 J per eV   exact, SI 2019
//   m_u = 1.66053906892(52)e-27 kg   CODATA 2022, relative 3.1e-10
// so T = 10.1805057178711931... fs. The literal is that rounded to a double;
// it does not claim precision beyond the 1.57e-10 relative uncertainty of m_u.
inline constexpr double kFemtosecondsPerInternalTime = 10.180505717871194;
inline constexpr double kInternalTimePerFemtosecond = 1.0 / kFemtosecondsPerInternalTime;

// Pressure at the interface boundary. 1 eV/A^3 = 1e-30 m^3 per A^3 and
// 1.602176634e-19 J per eV over 1e5 Pa per bar = 1602176.634 bar, exactly.
// The terminating direction is the literal so both directions round once.
inline constexpr double kEvPerCubicAngstromToBar = 1602176.634;
inline constexpr double kBarToEvPerCubicAngstrom = 1.0 / kEvPerCubicAngstromToBar;

// k_B = 1.380649e-23 J/K over 1.602176634e-19 J/eV; both are exact since the
// 2019 SI redefinition. Used by later temperature records, not by this layer.
inline constexpr double kBoltzmannEvPerKelvin = 8.617333262145177e-5;

static_assert(kFemtosecondsPerInternalTime * kInternalTimePerFemtosecond == 1.0,
              "the internal-time conversions must be exact reciprocals");
static_assert(kEvPerCubicAngstromToBar * kBarToEvPerCubicAngstrom == 1.0,
              "the pressure conversions must be exact reciprocals");
static_assert(kEvPerCubicAngstromToBar == 801088317.0 / 500.0,
              "eV/A^3 -> bar must be the nearest double to 801088317/500");

}  // namespace units

// Conversions to and from the integrator's time unit. Non-finite input and unit
// systems that do not define the requested unit are rejected, never guessed.
double to_internal_time(double value, TimeUnit unit, UnitSystem system);
double from_internal_time(double internal_value, TimeUnit unit, UnitSystem system);

// Pressure and virial live in eV/A^3 internally; bar is an interface unit.
double pressure_to_bar(double ev_per_cubic_angstrom, UnitSystem system);
double pressure_from_bar(double bar, UnitSystem system);

// A point on the trajectory: the step index and its time in internal units.
struct StepStamp {
    std::int64_t step = 0;
    double time_internal = 0.0;

    friend bool operator==(const StepStamp&, const StepStamp&) = default;
};

ValidationReport validate_step_stamp(const StepStamp& stamp, std::string_view field);

}  // namespace gmd_next::core
