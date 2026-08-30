#pragma once

#include <cstdint>
#include <string_view>

#include "gmd/integrator/barostat.hpp"

namespace gmd {

// Berendsen barostat: weakly couples the simulation box to a pressure bath.
//
// Each step, coordinates and box lengths are scaled by a factor mu:
//   mu = (1 - beta * dt / tau_P * (P_target - P_current))^(1/3)
//
// P_current is estimated from the kinetic + virial pressure:
//   P = (2 * KE + virial_trace) / (3 * V)   [eV/A^3, then converted to bar]
//
// Parameters:
//   tau_P     — barostat relaxation time [same time unit as dt]
//   beta      — isothermal compressibility [1/bar] (default: the liquid-water
//               value, 4.5e-5 bar^-1)
//
// Pressure units: target_pressure arrives in bar, as the run input's `pressure`
// field. P_current is computed in GMD's internal eV/A^3 and converted to bar
// before the comparison, so the difference beta_ multiplies is in bar and beta_
// is a genuine bar^-1. This used to subtract the two without converting either.
class BerendsenBarostat final : public Barostat {
public:
    explicit BerendsenBarostat(double tau_P = 2000.0,
                               double beta  = 4.5e-5) noexcept
        : tau_P_(tau_P), beta_(beta) {}

    ~BerendsenBarostat() override = default;

    std::string_view name() const noexcept override { return "berendsen"; }

    void apply(System& system,
               ForceProvider& provider,
               RuntimeContext& runtime,
               std::uint64_t step,
               double dt,
               double temperature,
               double target_pressure,
               double virial_trace) override;

    void set_tau_P(double t) noexcept { tau_P_ = t; }
    void set_beta(double b)  noexcept { beta_  = b; }

private:
    double tau_P_;
    double beta_;
};

}  // namespace gmd
