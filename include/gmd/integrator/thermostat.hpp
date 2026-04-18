#pragma once

#include <cstddef>
#include <string_view>

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
    virtual void initialize(const System& system) noexcept {}

    // Applied at the end of a full Velocity Verlet step.
    // For simple rescaling schemes (Velocity Rescaling, Berendsen) this is
    // the only hook that needs to be implemented.
    virtual void apply(System& system, double dt, double target_temperature) = 0;

    // Optional hook called in both half-kicks for extended-system thermostats
    // (e.g., Nosé-Hoover).  Default implementation is a no-op so that simple
    // thermostats do not need to override it.
    virtual void apply_half_kick(System& system,
                                 double dt_half,
                                 double target_temperature) noexcept {}
};

// --- Shared kinetic-energy helpers used by multiple thermostats ---------------

// Returns 2 * KE = sum_i  m_i * v_i^2  (mass-weighted velocity squared sum).
double compute_twice_ke(const System& system) noexcept;

// Returns the instantaneous temperature for a system with `dof` degrees of
// freedom.  kB is in units consistent with the rest of the code (eV/K).
double temperature_from_twice_ke(double twice_ke, std::size_t dof) noexcept;

// Boltzmann constant [eV/K].
inline constexpr double kBoltzmann = 8.617333262e-5;

}  // namespace gmd
