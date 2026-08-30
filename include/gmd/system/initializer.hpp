#pragma once

#include <cstdint>

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
class VelocityInitializer {
public:
    explicit VelocityInitializer(std::uint32_t seed = 5489u) noexcept;

    void initialize(System& system,
                    double target_temperature,
                    VelocityInitMode mode = VelocityInitMode::Random,
                    bool remove_center_of_mass_velocity = true) const;

    // Kinetic energy of this rank's OWN atoms, in eV. Not reduced: callers that
    // want the whole system's must reduce it themselves, as rescale_temperature
    // does. Ghost atoms are excluded -- they are copies of atoms another rank
    // owns, and counting them would double their contribution.
    double kinetic_energy(const System& system) const noexcept;

    std::uint32_t seed() const noexcept { return seed_; }

private:
    void sample_random_velocities(System& system, double target_temperature) const;
    void remove_center_of_mass_velocity(System& system) const;
    void rescale_temperature(System& system,
                             double target_temperature,
                             bool center_of_mass_removed) const;
    // Collective under MPI. Throws on every rank if any global atom tag is
    // duplicated or negative, because the tag is the random draw's identity.
    void validate_atom_tags(const System& system) const;

    // No generator state: every draw is a pure function of
    // (seed, stream, atom tag, component). See gmd/core/keyed_random.hpp for
    // why, and for what is and is not reproducible.
    std::uint32_t seed_;
};

}  // namespace gmd