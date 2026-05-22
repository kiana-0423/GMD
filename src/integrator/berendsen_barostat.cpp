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
    const std::size_t n = global_atom_count(system);
    if (n == 0) return;

    const Box& box = system.box();
    const double volume = box.lengths[0] * box.lengths[1] * box.lengths[2];
    if (volume < 1e-30) return;

    const double twice_ke = compute_twice_ke(system);
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

    // Scale only locally owned coordinates; ghost atoms are refreshed later.
    auto coords = system.mutable_coordinates();
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        coords[atom_index][0] *= mu;
        coords[atom_index][1] *= mu;
        coords[atom_index][2] *= mu;
    }
    system.mutable_neighbor_list().valid = false;

    // Velocities are NOT scaled in the standard Berendsen barostat;
    // only coordinates move with the box.
}

}  // namespace gmd
