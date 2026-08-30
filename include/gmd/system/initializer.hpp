#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>

namespace gmd {

class System;

enum class VelocityInitMode : std::uint8_t {
    Random,
    FromInput,
};

// Initializes per-atom velocities after coordinates and masses have been loaded.
//
// Random velocities are a deterministic function of the atom's stable global
// TAG, not of its position in local storage. The same physical atom therefore
// receives the same velocity whether the run is serial or decomposed over any
// number of ranks, and whatever order the atoms happen to be stored in. See
// gmd/core/keyed_random.hpp.
// What the initializer needs to know about constraints, supplied by whoever
// owns the constraint solver.
//
// The initializer cannot work this out for itself: it holds no solver, and the
// authoritative rank is only established once the positions have been projected
// onto the constraint manifold. Simulation::initialize() therefore runs the
// integrator's initialize() FIRST -- which projects the positions, accepts the
// rank and installs the thermostat's degrees of freedom -- and hands the result
// down here.
//
// Both fields are optional. A run with no constraints leaves them empty and the
// initializer falls back to 3N-3 (or 3N), which for an unconstrained system is
// the same number the authoritative count would give.
struct VelocityConstraintContext {
    // The authoritative count: 3N - rank - 3, the same value temperature
    // reporting, the thermostats and the pressure use. Zero means "not
    // supplied", not "no degrees of freedom".
    std::size_t degrees_of_freedom = 0;

    // The authoritative velocity projection -- RATTLE, reached through the
    // integrator rather than reimplemented here. Empty means unconstrained.
    std::function<void(System&)> project_velocities;

    bool has_constraints() const noexcept {
        return static_cast<bool>(project_velocities);
    }
};

class VelocityInitializer {
public:
    explicit VelocityInitializer(std::uint32_t seed = 5489u) noexcept;

    void initialize(System& system,
                    double target_temperature,
                    VelocityInitMode mode = VelocityInitMode::Random,
                    bool remove_center_of_mass_velocity = true) const;

    // Constraint-aware form. See VelocityConstraintContext, and the note on the
    // sequence in src/system/initializer.cpp.
    void initialize(System& system,
                    double target_temperature,
                    VelocityInitMode mode,
                    bool remove_center_of_mass_velocity,
                    const VelocityConstraintContext& constraints) const;

    // Kinetic energy of this rank's OWN atoms, in eV. Not reduced: callers that
    // want the whole system's must reduce it themselves, as rescale_temperature
    // does. Ghost atoms are excluded -- they are copies of atoms another rank
    // owns, and counting them would double their contribution.
    double kinetic_energy(const System& system) const noexcept;

    std::uint32_t seed() const noexcept { return seed_; }

private:
    void sample_random_velocities(System& system, double target_temperature) const;
    void remove_center_of_mass_velocity(System& system) const;
    // `override_dof` of zero means "use the unconstrained 3N-3 / 3N rule".
    void rescale_temperature(System& system,
                             double target_temperature,
                             bool center_of_mass_removed,
                             std::size_t override_dof = 0) const;
    // Global 2K, reduced across ranks.
    double global_twice_kinetic_energy(const System& system) const;
    // Global |sum m v|, reduced across ranks.
    double global_momentum_magnitude(const System& system) const;
    // Collective under MPI. Throws on every rank if any global atom tag is
    // duplicated or negative, because the tag is the random draw's identity.
    void validate_atom_tags(const System& system) const;

    // No generator state: every draw is a pure function of
    // (seed, stream, atom tag, component). See gmd/core/keyed_random.hpp for
    // why, and for what is and is not reproducible.
    std::uint32_t seed_;
};

}  // namespace gmd