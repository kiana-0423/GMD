#pragma once

#include <string_view>

#include "gmd/integrator/thermostat.hpp"

namespace gmd {

// Simple velocity-rescaling thermostat.
//
// At the end of each step all velocities are uniformly scaled so that the
// instantaneous kinetic temperature equals T_target exactly:
//
//   lambda = sqrt(T_target / T_current)
//   v_i   *= lambda   for all i
//
// This approach guarantees the target temperature every step but does NOT
// sample the canonical (NVT) ensemble correctly — use Nosé-Hoover for that.
// It is however an excellent choice for equilibration runs.
//
// Degrees of freedom are supplied by the integrator via
// set_degrees_of_freedom(); see compute_degrees_of_freedom() in thermostat.hpp.
class VelocityRescalingThermostat final : public Thermostat {
public:
    VelocityRescalingThermostat() = default;
    ~VelocityRescalingThermostat() override = default;

    std::string_view name() const noexcept override { return "velocity_rescaling"; }

    void initialize(const System& system) noexcept override;

    // Scale all velocities to match target_temperature exactly.
    void apply(System& system, double dt, double target_temperature) override;

private:
    // dof_ is inherited from Thermostat.
    double current_temperature_ = 0.0; // for diagnostics
};

}  // namespace gmd
