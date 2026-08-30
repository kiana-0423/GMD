#pragma once

#include <cstddef>
#include <string>
#include <string_view>

#include "gmd/core/physical_constants.hpp"

namespace gmd {

class System;

// Abstract interface for temperature coupling schemes.
// A Thermostat modifies the velocities (and possibly auxiliary variables)
// of a System to drive it towards a target temperature.
//
// Thermostats are applied by the integrator at the end of each step
// after the second half-kick.  Nosé-Hoover additionally needs to inject
// its friction terms inside the half-kicks; it overrides apply_half_kick().
class Thermostat {
public:
    virtual ~Thermostat() = default;

    virtual std::string_view name() const noexcept = 0;

    // Called once before the first step with the system at t=0.
    virtual void initialize(const System&) noexcept {}

    // Applied at the end of a full Velocity Verlet step.
    // For simple rescaling schemes (Velocity Rescaling, Berendsen) this is
    // the only hook that needs to be implemented.
    virtual void apply(System& system, double dt, double target_temperature) = 0;

    // Optional hook called in both half-kicks for extended-system thermostats
    // (e.g., Nosé-Hoover).  Default implementation is a no-op so that simple
    // thermostats do not need to override it.
    virtual void apply_half_kick(System&,
                                 double,
                                 double) noexcept {}

    // Text state used by checkpoint/restart. Stateless thermostats can keep
    // the defaults; stateful implementations should include enough data to
    // continue a deterministic run.
    virtual std::string checkpoint_state() const { return "stateless"; }
    virtual void load_checkpoint_state(const std::string& /*state*/) {}

    // Degrees of freedom used to convert kinetic energy into temperature.
    //
    // initialize() installs a provisional default that assumes an unconstrained
    // system with the centre-of-mass velocity removed. The integrator then calls
    // set_degrees_of_freedom() with the authoritative count -- the one that also
    // accounts for constraints and for runs that keep the COM velocity -- which
    // additionally marks the value as authoritative.
    //
    // The distinction matters at restart: load_checkpoint_state() validates the
    // checkpoint against the authoritative count and must be able to tell that
    // the count it is checking against really came from the run configuration.
    void set_degrees_of_freedom(std::size_t dof) noexcept {
        dof_ = dof;
        dof_is_authoritative_ = true;
    }
    std::size_t degrees_of_freedom() const noexcept { return dof_; }
    bool degrees_of_freedom_is_authoritative() const noexcept {
        return dof_is_authoritative_;
    }

protected:
    // Used by initialize() to install the provisional value; deliberately does
    // not mark it authoritative.
    void set_default_degrees_of_freedom(std::size_t dof) noexcept {
        dof_ = dof;
        dof_is_authoritative_ = false;
    }

    std::size_t dof_ = 0;

private:
    bool dof_is_authoritative_ = false;
};

// --- Shared kinetic-energy helpers used by multiple thermostats ---------------

// Returns 2 * KE = sum_i  m_i * v_i^2  (mass-weighted velocity squared sum).
double compute_twice_ke(const System& system) noexcept;
std::size_t global_atom_count(const System& system) noexcept;

// Returns the instantaneous temperature for a system with `dof` degrees of
// freedom.  kB is in units consistent with the rest of the code (eV/K).
double temperature_from_twice_ke(double twice_ke, std::size_t dof) noexcept;

// --- Shared degrees-of-freedom calculation ----------------------------------
//
// Every temperature consumer (thermostats, trajectory output, diagnostics)
// must agree on the same DOF count, so the formula lives here and nowhere
// else.  Starting from 3N:
//
//   - subtract 3 when the centre-of-mass velocity is removed (needs N >= 2,
//     otherwise the three translational modes are all the system has),
//   - subtract one per distinct active holonomic constraint.
//
// `constraint_count` is a *global* count of constraints, and it is the number
// of degrees of freedom they remove because a dependent set is rejected before
// dynamics start: VelocityVerletIntegrator::initialize() calls
// ConstraintSolver::require_independent(), which throws unless the count equals
// the rank of the mass-weighted constraint Jacobian. ConstraintSolver stores
// constraints against stable global atom tags and replicates the same list on
// every rank, so this count must not be reduced again under MPI.
struct DegreesOfFreedomConfig {
    bool remove_center_of_mass_velocity = true;
    std::size_t constraint_count = 0;
};

// Global-atom-count overload. Saturates at zero rather than wrapping around:
// an over-constrained or degenerate system reports 0 DOF, and callers treat
// that as "temperature is undefined" instead of dividing by a huge number.
std::size_t compute_degrees_of_freedom(std::size_t global_atom_count,
                                       const DegreesOfFreedomConfig& config) noexcept;

// Convenience overload that resolves the global atom count from `system`
// (MPI-aware via global_atom_count()).
std::size_t compute_degrees_of_freedom(const System& system,
                                       const DegreesOfFreedomConfig& config) noexcept;

// Boltzmann constant [eV/K]. The name is kept for the call sites that
// already use it; the value has exactly one definition, in
// gmd/core/physical_constants.hpp, and this must never become a second
// literal.
inline constexpr double kBoltzmann = kBoltzmannConstantEVPerKelvin;

}  // namespace gmd
