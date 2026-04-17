#include "gmd/core/simulation.hpp"

#include <stdexcept>
#include <utility>

#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/integrator.hpp"
#include "gmd/neighbor/neighbor_builder.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

namespace gmd {

class Simulation::Impl {
public:
    System* system = nullptr;
    std::shared_ptr<ForceProvider> force_provider;
    std::shared_ptr<NeighborBuilder> neighbor_builder;
    std::shared_ptr<Integrator> integrator;
    std::shared_ptr<VelocityInitializer> velocity_initializer;
    VelocityInitMode velocity_init_mode = VelocityInitMode::Random;
    bool remove_center_of_mass_velocity = true;
    double initial_temperature = 0.0;
    double time_step = 0.0;
    std::uint64_t step = 0;
};

Simulation::Simulation() noexcept
    : impl_(std::make_unique<Impl>()) {}

Simulation::Simulation(System* system) noexcept
    : Simulation() {
    set_system(system);
}

Simulation::~Simulation() = default;

Simulation::Simulation(Simulation&&) noexcept = default;

Simulation& Simulation::operator=(Simulation&&) noexcept = default;

void Simulation::set_system(System* system) noexcept {
    impl_->system = system;
}

void Simulation::set_force_provider(std::shared_ptr<ForceProvider> provider) noexcept {
    impl_->force_provider = std::move(provider);
}

void Simulation::set_neighbor_builder(std::shared_ptr<NeighborBuilder> builder) noexcept {
    impl_->neighbor_builder = std::move(builder);
}

void Simulation::set_integrator(std::shared_ptr<Integrator> integrator) noexcept {
    impl_->integrator = std::move(integrator);
}

void Simulation::set_velocity_initializer(std::shared_ptr<VelocityInitializer> initializer) noexcept {
    impl_->velocity_initializer = std::move(initializer);
}

void Simulation::set_velocity_init_mode(VelocityInitMode mode) noexcept {
    impl_->velocity_init_mode = mode;
}

void Simulation::set_remove_center_of_mass_velocity(bool enabled) noexcept {
    impl_->remove_center_of_mass_velocity = enabled;
}

void Simulation::set_initial_temperature(double temperature) noexcept {
    impl_->initial_temperature = temperature;
}

void Simulation::set_time_step(double time_step) noexcept {
    impl_->time_step = time_step;
}

bool Simulation::ready() const noexcept {
    return impl_->system != nullptr && impl_->force_provider != nullptr && impl_->integrator != nullptr;
}

void Simulation::initialize(RuntimeContext& runtime) {
    if (impl_->system == nullptr) {
        throw std::runtime_error("Simulation requires a System before initialization");
    }

    if (impl_->velocity_initializer != nullptr) {
        impl_->velocity_initializer->initialize(*impl_->system,
                                                impl_->initial_temperature,
                                                impl_->velocity_init_mode,
                                                impl_->remove_center_of_mass_velocity);
    }

    if (impl_->neighbor_builder != nullptr) {
        impl_->neighbor_builder->initialize(*impl_->system, runtime);
    }
    if (impl_->force_provider != nullptr) {
        impl_->force_provider->initialize(runtime);
    }
    if (impl_->integrator != nullptr) {
        impl_->integrator->initialize(*impl_->system, runtime);
    }
    impl_->step = 0;
}

void Simulation::step(RuntimeContext& runtime) {
    if (!ready()) {
        throw std::runtime_error("Simulation is not ready to step");
    }

    if (impl_->neighbor_builder != nullptr && impl_->neighbor_builder->needs_rebuild(*impl_->system, impl_->step)) {
        impl_->neighbor_builder->rebuild(*impl_->system, runtime, nullptr);
    }

    const IntegratorStepContext step_context{.step = impl_->step, .dt = impl_->time_step};
    impl_->integrator->step(*impl_->system, *impl_->force_provider, step_context, runtime);
    ++impl_->step;
}

void Simulation::run(RuntimeContext& runtime, std::uint64_t steps) {
    for (std::uint64_t iteration = 0; iteration < steps; ++iteration) {
        step(runtime);
    }
}

const System* Simulation::system() const noexcept {
    return impl_->system;
}

System* Simulation::mutable_system() noexcept {
    return impl_->system;
}

}  // namespace gmd