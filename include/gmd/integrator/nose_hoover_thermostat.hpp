#pragma once

#include <cstddef>
#include <string_view>

#include "gmd/integrator/thermostat.hpp"

namespace gmd {

// Nosé-Hoover chain thermostat (single-chain, i.e., one auxiliary variable).
//
// The extended-system equations of motion are:
//
//   dp_i/dt  = F_i - xi * p_i
//   dxi/dt   = (G - Q*xi^2) / Q,   G = sum_i p_i^2/m_i - dof*kB*T
//
// where xi is the thermostat friction variable and Q is the thermostat mass
// (effective coupling strength).  This generates the canonical (NVT) ensemble.
//
// Integration is done with the Nosé-Hoover / velocity-Verlet (VVNH) splitting:
//   half-kick for particles AND half-update for xi (apply_half_kick)
//   drift (coordinates)
//   force evaluation
//   half-kick for particles AND half-update for xi (apply_half_kick)
//
// Coupling strength Q = dof * kB * T_target * tau^2.
// Typical relaxation time tau: 0.05–0.2 ps.
class NoseHooverThermostat final : public Thermostat {
public:
    // tau [ps-equivalent time units as used in the simulation, typically fs].
    explicit NoseHooverThermostat(double tau = 100.0) noexcept;
    ~NoseHooverThermostat() override = default;

    std::string_view name() const noexcept override { return "nose_hoover"; }

    void initialize(const System& system) noexcept override;

    // Full-step apply: a second half-kick is already handled by apply_half_kick
    // inside the integrator.  For VVNH this is effectively a no-op here because
    // both half-kicks are driven via apply_half_kick().
    void apply(System& /*system*/, double /*dt*/, double /*T*/) override {}

    // Half-kick: evolve xi by dt/2, then scale all velocities by exp(-xi*dt/2).
    void apply_half_kick(System& system,
                         double dt_half,
                         double target_temperature) noexcept override;

private:
    double tau_  = 100.0;   // relaxation time (same units as dt)
    double xi_   = 0.0;     // friction variable
    double Q_    = 0.0;     // thermostat mass, computed in initialize()
    std::size_t dof_ = 0;
    double current_temperature_ = 0.0;
};

}  // namespace gmd
