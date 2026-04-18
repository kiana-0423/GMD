#include "gmd/integrator/berendsen_barostat.hpp"

#include <cmath>

#include "gmd/integrator/thermostat.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace gmd {

void BerendsenBarostat::apply(System& system,
                               ForceProvider& /*provider*/,
                               RuntimeContext& /*runtime*/,
                               std::uint64_t  /*step*/,
                               double dt,
                               double /*temperature*/,
                               double target_pressure,
                               double virial_trace) {
    const std::size_t n = system.atom_count();
    if (n == 0) return;

    const Box& box = system.box();
    const double volume = box.lengths[0] * box.lengths[1] * box.lengths[2];
    if (volume < 1e-30) return;

    // Instantaneous kinetic pressure contribution using current kE.
    const double twice_ke = compute_twice_ke(system);
    // P_kin = (2 * KE) / (3 * V)  (each direction contributes kT to pressure)
    const double p_kin = twice_ke / (3.0 * volume);

    // Total pressure from virial theorem: P = (2*KE + W) / (3*V)
    const double p_current = (twice_ke + virial_trace) / (3.0 * volume);

    // Berendsen coupling factor: scale volume.
    // mu³ = 1 - (beta * dt / tau_P) * (P_target - P_current)
    const double mu3 = 1.0 - beta_ * (dt / tau_P_) * (target_pressure - p_current);
    if (mu3 <= 0.0) return;  // guard against unphysical scaling
    const double mu = std::cbrt(mu3);

    // Scale box lengths.
    const std::array<double, 3> new_lengths = {
        box.lengths[0] * mu,
        box.lengths[1] * mu,
        box.lengths[2] * mu
    };
    system.mutable_box().set_lengths(new_lengths);

    // Scale all coordinates and velocities by mu.
    auto coords = system.mutable_coordinates();
    for (auto& c : coords) {
        c[0] *= mu;
        c[1] *= mu;
        c[2] *= mu;
    }

    // Velocities are NOT scaled in the standard Berendsen barostat;
    // only coordinates move with the box.
}

}  // namespace gmd
