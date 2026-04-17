#pragma once

#include <cstdint>
#include <random>

namespace gmd {

class System;

enum class VelocityInitMode : std::uint8_t {
    Random,
    FromInput,
};

// Initializes per-atom velocities after coordinates and masses have been loaded.
class VelocityInitializer {
public:
    explicit VelocityInitializer(std::uint32_t seed = 5489u) noexcept;

    void initialize(System& system,
                    double target_temperature,
                    VelocityInitMode mode = VelocityInitMode::Random,
                    bool remove_center_of_mass_velocity = true) const;
    double kinetic_energy(const System& system) const noexcept;

private:
    void sample_random_velocities(System& system, double target_temperature) const;
    void remove_center_of_mass_velocity(System& system) const;
    void rescale_temperature(System& system,
                             double target_temperature,
                             bool center_of_mass_removed) const;

    mutable std::mt19937 generator_;
};

}  // namespace gmd