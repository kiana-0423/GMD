#include "gmd/integrator/velocity_verlet_integrator.hpp"

#include <algorithm>
#include <stdexcept>

#include "gmd/boundary/periodic_boundary.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/system.hpp"

namespace gmd {

namespace {

ForceResult evaluate_force(System& system,
                           ForceProvider& force_provider,
                           std::uint64_t step,
                           double time,
                           RuntimeContext& runtime) {
    ForceResult result;
    const auto coordinates = system.coordinates();
    ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .step = step,
        .time = time,
        .coordinates = std::span<const Coordinate3D>(coordinates.data(), coordinates.size()),
    };
    force_provider.compute(request, result, runtime);
    if (!result.success) {
        throw std::runtime_error("Force provider reported an unsuccessful force evaluation");
    }
    return result;
}

void copy_forces_to_system(System& system, const ForceResult& result) {
    auto system_forces = system.mutable_forces();
    const auto copy_count = std::min(system_forces.size(), result.forces.size());
    for (std::size_t index = 0; index < copy_count; ++index) {
        system_forces[index] = result.forces[index];
    }
    for (std::size_t index = copy_count; index < system_forces.size(); ++index) {
        system_forces[index] = {0.0, 0.0, 0.0};
    }
    system.set_potential_energy(result.potential_energy);
}

}  // namespace

VelocityVerletIntegrator::VelocityVerletIntegrator(double dt) noexcept
    : dt_(dt) {}

std::string_view VelocityVerletIntegrator::name() const noexcept {
    return "velocity_verlet";
}

void VelocityVerletIntegrator::initialize(System& system, RuntimeContext& runtime) {
    (void)runtime;
    auto forces = system.mutable_forces();
    for (auto& force : forces) {
        force = {0.0, 0.0, 0.0};
    }
}

void VelocityVerletIntegrator::step(System& system,
                                    ForceProvider& force_provider,
                                    const IntegratorStepContext& ctx,
                                    RuntimeContext& runtime) {
    const double dt = ctx.dt > 0.0 ? ctx.dt : dt_;
    if (dt <= 0.0) {
        throw std::runtime_error("VelocityVerletIntegrator requires a positive time step");
    }

    ForceResult current_force = evaluate_force(system, force_provider, ctx.step, ctx.step * dt, runtime);
    copy_forces_to_system(system, current_force);

    const auto masses = system.masses();
    auto coordinates = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();

    for (std::size_t atom_index = 0; atom_index < system.atom_count(); ++atom_index) {
        if (masses[atom_index] <= 0.0) {
            throw std::runtime_error("VelocityVerletIntegrator requires strictly positive masses");
        }

        const double inverse_mass = 1.0 / masses[atom_index];
        for (std::size_t dim = 0; dim < 3; ++dim) {
            velocities[atom_index][dim] += 0.5 * current_force.forces[atom_index][dim] * inverse_mass * dt;
            coordinates[atom_index][dim] += velocities[atom_index][dim] * dt;
        }
        wrap_position(coordinates[atom_index], system.box());
    }

    ForceResult next_force = evaluate_force(system, force_provider, ctx.step + 1, (ctx.step + 1) * dt, runtime);
    copy_forces_to_system(system, next_force);

    for (std::size_t atom_index = 0; atom_index < system.atom_count(); ++atom_index) {
        const double inverse_mass = 1.0 / masses[atom_index];
        for (std::size_t dim = 0; dim < 3; ++dim) {
            velocities[atom_index][dim] += 0.5 * next_force.forces[atom_index][dim] * inverse_mass * dt;
        }
    }
}

double VelocityVerletIntegrator::dt() const noexcept {
    return dt_;
}

void VelocityVerletIntegrator::set_dt(double dt) noexcept {
    dt_ = dt;
}

}  // namespace gmd