#include "gmd/core/simulation.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>

#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/integrator.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/neighbor/neighbor_builder.hpp"
#include "gmd/neighbor/verlet_neighbor_builder.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
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
    std::shared_ptr<MpiCommunicator> mpi_comm;
    std::shared_ptr<DomainDecomposition> domain_decomposition;
    std::shared_ptr<VelocityInitializer> velocity_initializer;
    VelocityInitMode velocity_init_mode = VelocityInitMode::Random;
    bool remove_center_of_mass_velocity = true;
    double initial_temperature = 0.0;
    double time_step = 0.0;
    std::uint64_t step = 0;

    void sync_domain_box() {
        if (domain_decomposition != nullptr && system != nullptr) {
            domain_decomposition->refresh(system->box());
        }
    }

    void redistribute_owned_atoms() {
        if (system == nullptr || mpi_comm == nullptr || domain_decomposition == nullptr) {
            return;
        }

        sync_domain_box();
        mpi_comm->redistribute_atoms(*system, *domain_decomposition);
    }

    void prepare_force_evaluation(RuntimeContext& runtime, std::uint64_t force_step) {
        if (system == nullptr) {
            throw std::runtime_error("Simulation force evaluation requires a System");
        }

        if (mpi_comm != nullptr && domain_decomposition != nullptr) {
            sync_domain_box();
            mpi_comm->exchange_ghost_coordinates(*system, *domain_decomposition);
        }

        if (neighbor_builder != nullptr &&
            (!system->neighbor_list().valid || neighbor_builder->needs_rebuild(*system, force_step))) {
            neighbor_builder->rebuild(*system, runtime, nullptr);
        }
    }

    ForceResult evaluate_force(std::uint64_t force_step,
                               double force_time,
                               RuntimeContext& runtime) {
        prepare_force_evaluation(runtime, force_step);

        ForceResult result;
        const auto coordinates = system->coordinates();
        const NeighborList& neighbor_list = system->neighbor_list();
        ForceRequest request{
            .system = system,
            .box = &system->box(),
            .step = force_step,
            .time = force_time,
            .coordinates = std::span<const Coordinate3D>(coordinates.data(), coordinates.size()),
            .neighbor_list = neighbor_list.valid ? &neighbor_list : nullptr,
        };
        force_provider->compute(request, result, runtime);
        if (!result.success) {
            throw std::runtime_error("Force provider reported an unsuccessful force evaluation");
        }

        auto system_forces = system->mutable_forces();
        const auto copy_count = std::min(system_forces.size(), result.forces.size());
        for (std::size_t index = 0; index < copy_count; ++index) {
            system_forces[index] = result.forces[index];
        }
        for (std::size_t index = copy_count; index < system_forces.size(); ++index) {
            system_forces[index] = {0.0, 0.0, 0.0};
        }
        system->set_potential_energy(result.potential_energy);

        if (mpi_comm != nullptr && domain_decomposition != nullptr) {
            mpi_comm->reverse_accumulate_ghost_forces(*system, *domain_decomposition);
        }
        if (mpi_comm != nullptr) {
            system->set_potential_energy(
                mpi_comm->allreduce_scalar(system->potential_energy()));
        }

        return result;
    }
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
    if (impl_->domain_decomposition != nullptr) {
        auto verlet_builder =
            std::dynamic_pointer_cast<VerletNeighborBuilder>(impl_->neighbor_builder);
        if (verlet_builder != nullptr) {
            verlet_builder->set_domain_decomposition(impl_->domain_decomposition);
        }
    }
}

void Simulation::set_integrator(std::shared_ptr<Integrator> integrator) noexcept {
    impl_->integrator = std::move(integrator);
}

void Simulation::set_mpi_communicator(std::shared_ptr<MpiCommunicator> comm) noexcept {
    impl_->mpi_comm = std::move(comm);
}

void Simulation::set_domain_decomposition(std::shared_ptr<DomainDecomposition> dd) noexcept {
    impl_->domain_decomposition = std::move(dd);
    auto verlet_builder =
        std::dynamic_pointer_cast<VerletNeighborBuilder>(impl_->neighbor_builder);
    if (verlet_builder != nullptr) {
        verlet_builder->set_domain_decomposition(impl_->domain_decomposition);
    }
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

    if (impl_->force_provider != nullptr) {
        impl_->force_provider->initialize(runtime);
    }
    if (impl_->integrator != nullptr) {
        impl_->integrator->initialize(*impl_->system, runtime);
    }

    // Compute initial forces at t=0 so the first half-kick uses the current state.
    if (impl_->force_provider != nullptr) {
        impl_->evaluate_force(0, 0.0, runtime);
    }

    impl_->step = 0;
}

void Simulation::step(RuntimeContext& runtime) {
    if (!ready()) {
        throw std::runtime_error("Simulation is not ready to step");
    }

    const IntegratorStepContext step_context{.step = impl_->step, .dt = impl_->time_step};
    auto velocity_verlet =
        std::dynamic_pointer_cast<VelocityVerletIntegrator>(impl_->integrator);
    if (velocity_verlet == nullptr) {
        if (impl_->mpi_comm != nullptr &&
            impl_->domain_decomposition != nullptr &&
            impl_->mpi_comm->size() > 1) {
            throw std::runtime_error(
                "MPI execution currently requires VelocityVerletIntegrator");
        }

        if (impl_->neighbor_builder != nullptr &&
            (!impl_->system->neighbor_list().valid ||
             impl_->neighbor_builder->needs_rebuild(*impl_->system, impl_->step))) {
            impl_->neighbor_builder->rebuild(*impl_->system, runtime, nullptr);
        }
        impl_->integrator->step(*impl_->system, *impl_->force_provider, step_context, runtime);
        ++impl_->step;
        return;
    }

    velocity_verlet->begin_step(*impl_->system, step_context);
    impl_->redistribute_owned_atoms();

    const double next_time = static_cast<double>(impl_->step + 1) * impl_->time_step;
    const ForceResult next_force = impl_->evaluate_force(impl_->step + 1, next_time, runtime);
    velocity_verlet->finish_step(*impl_->system,
                                 step_context,
                                 next_force.virial_valid,
                                 next_force.virial);

    if (velocity_verlet->has_barostat()) {
        velocity_verlet->apply_barostat(*impl_->system,
                                        *impl_->force_provider,
                                        runtime,
                                        step_context);
        impl_->redistribute_owned_atoms();
        impl_->system->mutable_neighbor_list().valid = false;
        impl_->evaluate_force(impl_->step + 1, next_time, runtime);
    }

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
