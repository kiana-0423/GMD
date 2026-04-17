#pragma once

#include <cstdint>
#include <memory>

#include "gmd/system/initializer.hpp"

namespace gmd {

class ForceProvider;
class Integrator;
class NeighborBuilder;
class RuntimeContext;
class System;

// Coordinates the main MD execution loop and owns the runtime-facing pipeline.
class Simulation {
public:
    Simulation() noexcept;
    explicit Simulation(System* system) noexcept;
    ~Simulation();

    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept;
    Simulation& operator=(Simulation&&) noexcept;

    void set_system(System* system) noexcept;
    void set_force_provider(std::shared_ptr<ForceProvider> provider) noexcept;
    void set_neighbor_builder(std::shared_ptr<NeighborBuilder> builder) noexcept;
    void set_integrator(std::shared_ptr<Integrator> integrator) noexcept;
    void set_velocity_initializer(std::shared_ptr<VelocityInitializer> initializer) noexcept;
    void set_velocity_init_mode(VelocityInitMode mode) noexcept;
    void set_remove_center_of_mass_velocity(bool enabled) noexcept;
    void set_initial_temperature(double temperature) noexcept;
    void set_time_step(double time_step) noexcept;

    // Returns true when the simulation has the minimum components required to step.
    bool ready() const noexcept;

    // Initializes registered components against the selected runtime backend.
    void initialize(RuntimeContext& runtime);
    // Advances the system by one integration step.
    void step(RuntimeContext& runtime);
    // Repeats step() for the requested number of iterations.
    void run(RuntimeContext& runtime, std::uint64_t steps);

    const System* system() const noexcept;
    System* mutable_system() noexcept;

private:
    // Keeps implementation details out of the public ABI while headers are stabilizing.
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace gmd
