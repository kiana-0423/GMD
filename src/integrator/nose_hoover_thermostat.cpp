#include "gmd/integrator/nose_hoover_thermostat.hpp"

#include <cmath>

#include "gmd/system/system.hpp"

namespace gmd {

NoseHooverThermostat::NoseHooverThermostat(double tau) noexcept
    : tau_(tau) {}

void NoseHooverThermostat::initialize(const System& system) noexcept {
    const std::size_t n = system.atom_count();
    dof_ = (n >= 2) ? 3 * n - 3 : 3 * n;
    xi_  = 0.0;
    // Q will be set on first apply_half_kick once target_temperature is known.
    Q_ = 0.0;
}

void NoseHooverThermostat::apply_half_kick(System& system,
                                           double dt_half,
                                           double target_temperature) noexcept {
    if (dof_ == 0 || target_temperature <= 0.0) return;

    // Lazily initialise Q when target temperature is first known.
    if (Q_ <= 0.0) {
        Q_ = static_cast<double>(dof_) * kBoltzmann * target_temperature * tau_ * tau_;
    }

    const double twice_ke    = compute_twice_ke(system);
    current_temperature_     = temperature_from_twice_ke(twice_ke, dof_);
    const double G           = twice_ke - static_cast<double>(dof_) * kBoltzmann * target_temperature;

    // Half-step update of friction variable xi (velocity Verlet for xi).
    xi_ += (G / Q_) * dt_half;

    // Rescale velocities: v_i *= exp(-xi * dt_half).
    const double scale = std::exp(-xi_ * dt_half);
    auto velocities    = system.mutable_velocities();
    for (std::size_t i = 0; i < system.atom_count(); ++i) {
        velocities[i][0] *= scale;
        velocities[i][1] *= scale;
        velocities[i][2] *= scale;
    }
}

}  // namespace gmd
