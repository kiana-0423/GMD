#include "gmd/system/initializer.hpp"

#include <cmath>
#include <stdexcept>

#include "gmd/system/system.hpp"

namespace gmd {

namespace {

constexpr double kBoltzmannConstant = 8.617343e-5;

}  // namespace

VelocityInitializer::VelocityInitializer(std::uint32_t seed) noexcept
    : generator_(seed) {}

double VelocityInitializer::kinetic_energy(const System& system) const noexcept {
    const auto masses = system.masses();
    const auto velocities = system.velocities();

    double kinetic_energy = 0.0;
    for (std::size_t atom_index = 0; atom_index < velocities.size(); ++atom_index) {
        const auto& velocity = velocities[atom_index];
        const double velocity_squared = velocity[0] * velocity[0] +
                                        velocity[1] * velocity[1] +
                                        velocity[2] * velocity[2];
        kinetic_energy += masses[atom_index] * velocity_squared;
    }
    return kinetic_energy * 0.5;
}

void VelocityInitializer::sample_random_velocities(System& system, double target_temperature) const {
    auto masses = system.masses();
    auto velocities = system.mutable_velocities();

    for (std::size_t atom_index = 0; atom_index < system.atom_count(); ++atom_index) {
        if (masses[atom_index] <= 0.0) {
            throw std::runtime_error("All particle masses must be positive before velocity initialization");
        }

        if (target_temperature == 0.0) {
            velocities[atom_index] = {0.0, 0.0, 0.0};
            continue;
        }

        const double sigma = std::sqrt(kBoltzmannConstant * target_temperature / masses[atom_index]);
        std::normal_distribution<double> distribution(0.0, sigma);
        velocities[atom_index] = {distribution(generator_), distribution(generator_), distribution(generator_)};
    }
}

void VelocityInitializer::remove_center_of_mass_velocity(System& system) const {
    auto masses = system.masses();
    auto velocities = system.mutable_velocities();
    double total_mass = 0.0;
    System::Vec3 center_of_mass_velocity = {0.0, 0.0, 0.0};
    for (std::size_t atom_index = 0; atom_index < velocities.size(); ++atom_index) {
        const auto mass = masses[atom_index];
        if (mass <= 0.0) {
            throw std::runtime_error("All particle masses must be positive before velocity initialization");
        }
        total_mass += mass;
        const auto& velocity = velocities[atom_index];
        center_of_mass_velocity[0] += masses[atom_index] * velocity[0];
        center_of_mass_velocity[1] += masses[atom_index] * velocity[1];
        center_of_mass_velocity[2] += masses[atom_index] * velocity[2];
    }

    if (total_mass <= 0.0) {
        throw std::runtime_error("Total system mass must be positive before velocity initialization");
    }

    center_of_mass_velocity[0] /= total_mass;
    center_of_mass_velocity[1] /= total_mass;
    center_of_mass_velocity[2] /= total_mass;

    for (auto& velocity : velocities) {
        velocity[0] -= center_of_mass_velocity[0];
        velocity[1] -= center_of_mass_velocity[1];
        velocity[2] -= center_of_mass_velocity[2];
    }
}

void VelocityInitializer::rescale_temperature(System& system,
                                              double target_temperature,
                                              bool center_of_mass_removed) const {
    if (target_temperature == 0.0) {
        auto velocities = system.mutable_velocities();
        for (auto& velocity : velocities) {
            velocity = {0.0, 0.0, 0.0};
        }
        return;
    }

    const double current_kinetic_energy = kinetic_energy(system);
    if (current_kinetic_energy <= 0.0) {
        throw std::runtime_error("Velocity initialization produced zero kinetic energy");
    }

    const auto atom_count = static_cast<double>(system.atom_count());
    const auto dof = center_of_mass_removed ? 3.0 * atom_count - 3.0 : 3.0 * atom_count;
    if (dof <= 0.0) {
        throw std::runtime_error("Not enough degrees of freedom to define a temperature");
    }

    const double current_temperature = current_kinetic_energy * 2.0 /
                                       (dof * kBoltzmannConstant);
    const double scale_factor = std::sqrt(target_temperature / current_temperature);

    auto velocities = system.mutable_velocities();
    for (auto& velocity : velocities) {
        velocity[0] *= scale_factor;
        velocity[1] *= scale_factor;
        velocity[2] *= scale_factor;
    }
}

void VelocityInitializer::initialize(System& system,
                                     double target_temperature,
                                     VelocityInitMode mode,
                                     bool remove_center_of_mass_velocity_flag) const {
    if (target_temperature < 0.0) {
        throw std::runtime_error("Target temperature must be non-negative");
    }
    if (system.atom_count() == 0) {
        return;
    }

    if (mode == VelocityInitMode::Random) {
        sample_random_velocities(system, target_temperature);
    }

    if (remove_center_of_mass_velocity_flag) {
        remove_center_of_mass_velocity(system);
    }

    rescale_temperature(system, target_temperature, remove_center_of_mass_velocity_flag);
}

}  // namespace gmd