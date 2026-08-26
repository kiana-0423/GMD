#include "gmd/integrator/velocity_rescaling_thermostat.hpp"

#include <cmath>

#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

void VelocityRescalingThermostat::initialize(const System& system) noexcept {
    // Default assumption: unconstrained system with the COM velocity removed.
    // VelocityVerletIntegrator::initialize() overrides this immediately with
    // the authoritative count once the constraint solver and the COM-removal
    // setting are known.
    set_default_degrees_of_freedom(
        compute_degrees_of_freedom(system, DegreesOfFreedomConfig{}));
}

void VelocityRescalingThermostat::apply(System& system,
                                        double /*dt*/,
                                        double target_temperature) {
    if (target_temperature <= 0.0 || dof_ == 0) return;

    const double twice_ke = compute_twice_ke(system);
    current_temperature_  = temperature_from_twice_ke(twice_ke, dof_);

    if (current_temperature_ < 1e-12) return;  // avoid division by zero

    const double lambda = std::sqrt(target_temperature / current_temperature_);
    auto velocities = system.mutable_velocities();
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        velocities[i][0] *= lambda;
        velocities[i][1] *= lambda;
        velocities[i][2] *= lambda;
    }
    current_temperature_ = target_temperature;
}

}  // namespace gmd
