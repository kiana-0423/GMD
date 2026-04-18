#pragma once

#include "gmd/integrator/barostat.hpp"

namespace gmd {

// Berendsen barostat: weakly couples the simulation box to a pressure bath.
//
// Each step, coordinates and box lengths are scaled by a factor mu:
//   mu = (1 - beta * dt / tau_P * (P_target - P_current))^(1/3)
//
// P_current is estimated from the kinetic + virial pressure:
//   P = (N_dof * kB * T + virial_trace/3) / (3 * V)
//   (the kinetic term uses the instantaneous temperature passed in)
//
// Parameters:
//   tau_P     — barostat relaxation time [same time unit as dt]
//   beta      — isothermal compressibility [1/pressure unit] (default: liquid-water value ~4.5e-5 bar^-1, but units are driver-specific)
class BerendsenBarostat final : public Barostat {
public:
    explicit BerendsenBarostat(double tau_P = 2000.0,
                               double beta  = 4.5e-5) noexcept
        : tau_P_(tau_P), beta_(beta) {}

    ~BerendsenBarostat() override = default;

    void apply(System& system, double dt, double target_pressure,
               double virial_trace) override;

    double tau_P() const noexcept { return tau_P_; }
    double beta()  const noexcept { return beta_;  }

    void set_tau_P(double t) noexcept { tau_P_ = t; }
    void set_beta(double b)  noexcept { beta_  = b; }

private:
    double tau_P_;
    double beta_;
};

}  // namespace gmd
