#include "gmd/integrator/thermostat.hpp"

#include "gmd/system/system.hpp"

namespace gmd {

double compute_twice_ke(const System& system) noexcept {
    const auto masses    = system.masses();
    const auto velocities = system.velocities();
    double twice_ke = 0.0;
    for (std::size_t i = 0; i < system.atom_count(); ++i) {
        const double m = masses[i];
        const auto&  v = velocities[i];
        twice_ke += m * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    return twice_ke;
}

double temperature_from_twice_ke(double twice_ke, std::size_t dof) noexcept {
    if (dof == 0) return 0.0;
    return twice_ke / (static_cast<double>(dof) * kBoltzmann);
}

}  // namespace gmd
