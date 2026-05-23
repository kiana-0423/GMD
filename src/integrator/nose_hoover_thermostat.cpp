#include "gmd/integrator/nose_hoover_thermostat.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

NoseHooverThermostat::NoseHooverThermostat(double tau) noexcept
    : tau_(tau) {}

void NoseHooverThermostat::initialize(const System& system) noexcept {
    const std::size_t n = global_atom_count(system);
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
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        velocities[i][0] *= scale;
        velocities[i][1] *= scale;
        velocities[i][2] *= scale;
    }
}

std::string NoseHooverThermostat::checkpoint_state() const {
    std::ostringstream out;
    out.precision(17);
    out << "tau " << tau_
        << " xi " << xi_
        << " Q " << Q_
        << " dof " << dof_
        << " current_temperature " << current_temperature_;
    return out.str();
}

void NoseHooverThermostat::load_checkpoint_state(const std::string& state) {
    if (state.empty() || state == "stateless") {
        return;
    }
    std::istringstream input(state);
    std::string key;
    if (!(input >> key) || key != "tau" || !(input >> tau_) ||
        !(input >> key) || key != "xi" || !(input >> xi_) ||
        !(input >> key) || key != "Q" || !(input >> Q_) ||
        !(input >> key) || key != "dof" || !(input >> dof_) ||
        !(input >> key) || key != "current_temperature" ||
        !(input >> current_temperature_)) {
        throw std::runtime_error("Invalid Nose-Hoover thermostat checkpoint state");
    }
}

}  // namespace gmd
