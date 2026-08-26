#include "gmd/integrator/nose_hoover_thermostat.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

NoseHooverThermostat::NoseHooverThermostat(double tau) noexcept
    : tau_(tau) {}

void NoseHooverThermostat::initialize(const System& system) noexcept {
    // Default assumption: unconstrained system with the COM velocity removed.
    // VelocityVerletIntegrator::initialize() overrides this immediately with
    // the authoritative count once the constraint solver and the COM-removal
    // setting are known.
    set_default_degrees_of_freedom(
        compute_degrees_of_freedom(system, DegreesOfFreedomConfig{}));
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

// Restart policy for the Nose-Hoover thermostat.
//
// The degrees of freedom are recomputed from the current system and run
// configuration before a restart reaches this point, and that value -- not the
// one in the checkpoint -- is authoritative. The checkpoint's `dof` is treated
// as a compatibility stamp: it is validated, never installed.
//
// Restoring only dof_ from the checkpoint would be unsound, because the
// thermostat mass Q = dof * kB * T_target * tau^2 and the friction variable xi
// were both produced under the checkpoint's DOF. A different DOF means a
// different extended system, so xi and Q carry no meaning across the change and
// there is no rescaling that makes the continued trajectory the same run.
// Mismatches are therefore rejected rather than migrated.
//
// Legacy checkpoints: the on-disk format is unchanged, so a checkpoint written
// before degrees of freedom accounted for constraints and the COM setting is
// simply one whose `dof` holds the old 3N-3 value. It is accepted transparently
// when that happens to equal the authoritative count (an unconstrained run with
// the COM velocity removed) and rejected with an explanatory error otherwise.
// Such a run must be restarted from its input rather than continued.
void NoseHooverThermostat::load_checkpoint_state(const std::string& state) {
    if (state.empty() || state == "stateless") {
        return;
    }

    double checkpoint_tau = 0.0;
    double checkpoint_xi = 0.0;
    double checkpoint_q = 0.0;
    long long checkpoint_dof = 0;
    double checkpoint_temperature = 0.0;

    std::istringstream input(state);
    std::string key;
    // Field order and names are unchanged, so checkpoints written by earlier
    // versions parse here exactly as they always did.
    if (!(input >> key) || key != "tau" || !(input >> checkpoint_tau) ||
        !(input >> key) || key != "xi" || !(input >> checkpoint_xi) ||
        !(input >> key) || key != "Q" || !(input >> checkpoint_q) ||
        !(input >> key) || key != "dof" || !(input >> checkpoint_dof) ||
        !(input >> key) || key != "current_temperature" ||
        !(input >> checkpoint_temperature)) {
        throw std::runtime_error("Invalid Nose-Hoover thermostat checkpoint state");
    }

    if (checkpoint_dof <= 0) {
        throw std::runtime_error(
            "Nose-Hoover checkpoint records " + std::to_string(checkpoint_dof) +
            " degrees of freedom; a restart needs a positive count");
    }

    if (!degrees_of_freedom_is_authoritative()) {
        throw std::runtime_error(
            "Nose-Hoover checkpoint cannot be restored before the run's degrees "
            "of freedom are known: initialize the integrator (which installs the "
            "authoritative count from the atom count, the centre-of-mass setting "
            "and the active constraints) before loading thermostat state");
    }
    if (dof_ == 0) {
        throw std::runtime_error(
            "Nose-Hoover restart has zero degrees of freedom after removing "
            "centre-of-mass motion and constraints; temperature is undefined");
    }

    if (static_cast<std::size_t>(checkpoint_dof) != dof_) {
        throw std::runtime_error(
            "Nose-Hoover checkpoint was written with " +
            std::to_string(checkpoint_dof) + " degrees of freedom but this run has " +
            std::to_string(dof_) +
            ". The thermostat mass Q and friction variable xi were generated under "
            "the checkpoint's value, so the run cannot be continued. Likely causes: "
            "the constraint configuration changed, the centre-of-mass removal "
            "setting changed, the atom count changed, or the checkpoint predates "
            "degrees of freedom accounting for constraints (the old 3N-3 rule). "
            "Restart from the input instead of the checkpoint.");
    }

    // Compatible: restore the extended-system state verbatim so the continued
    // trajectory is deterministic. dof_ is left at the authoritative value,
    // which the check above proved equal to the checkpoint's.
    tau_ = checkpoint_tau;
    xi_ = checkpoint_xi;
    Q_ = checkpoint_q;
    current_temperature_ = checkpoint_temperature;
}

}  // namespace gmd
