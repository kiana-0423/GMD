#include "gmd/integrator/velocity_verlet_integrator.hpp"

#include <algorithm>
#include <stdexcept>

#include "gmd/boundary/periodic_boundary.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/barostat.hpp"
#include "gmd/integrator/thermostat.hpp"
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
    const NeighborList& nl = system.neighbor_list();
    ForceRequest request{
        .system        = &system,
        .box           = &system.box(),
        .step          = step,
        .time          = time,
        .coordinates   = std::span<const Coordinate3D>(coordinates.data(), coordinates.size()),
        .neighbor_list = nl.valid ? &nl : nullptr,
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
    last_virial_trace_ = 0.0;
    last_virial_valid_ = false;
    auto forces = system.mutable_forces();
    for (auto& force : forces) {
        force = {0.0, 0.0, 0.0};
    }
    if (thermostat_) {
        thermostat_->initialize(system);
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

    // --- Thermostat pre-kick (Nosé-Hoover first half-kick or no-op) ---
    if (thermostat_) {
        thermostat_->apply_half_kick(system, 0.5 * dt, target_temperature_);
    }

    ForceResult current_force = evaluate_force(system, force_provider, ctx.step, ctx.step * dt, runtime);
    copy_forces_to_system(system, current_force);

    const auto masses = system.masses();
    auto coordinates = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();

    // First half-kick (v += 0.5*a*dt) + full position update.
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

    // Cache virial trace for barostat (uses full virial tensor if available).
    last_virial_valid_ = next_force.virial_valid;
    if (last_virial_valid_) {
        last_virial_trace_ = next_force.virial[0] + next_force.virial[4] + next_force.virial[8];
    }

    // Second half-kick.
    for (std::size_t atom_index = 0; atom_index < system.atom_count(); ++atom_index) {
        const double inverse_mass = 1.0 / masses[atom_index];
        for (std::size_t dim = 0; dim < 3; ++dim) {
            velocities[atom_index][dim] += 0.5 * next_force.forces[atom_index][dim] * inverse_mass * dt;
        }
    }

    // --- Thermostat post-kick (Nosé-Hoover second half-kick, or full rescaling step) ---
    if (thermostat_) {
        thermostat_->apply_half_kick(system, 0.5 * dt, target_temperature_);
        thermostat_->apply(system, dt, target_temperature_);
    }

    // --- Barostat (end of step) ---
    // MCBarostat works unconditionally (no virial required).
    // BerendsenBarostat is guarded: skip when virial is unavailable so that
    // pressure is not driven by a zero-virial (incorrect) estimate.
    if (barostat_) {
        const bool can_run = !barostat_->requires_virial() || last_virial_valid_;
        if (can_run) {
            barostat_->apply(system, force_provider, runtime,
                             ctx.step, dt, target_temperature_,
                             target_pressure_, last_virial_trace_);
        }
    }
}

double VelocityVerletIntegrator::dt() const noexcept {
    return dt_;
}

void VelocityVerletIntegrator::set_dt(double dt) noexcept {
    dt_ = dt;
}

void VelocityVerletIntegrator::set_thermostat(std::shared_ptr<Thermostat> thermostat) noexcept {
    thermostat_ = std::move(thermostat);
}

void VelocityVerletIntegrator::set_target_temperature(double temperature) noexcept {
    target_temperature_ = temperature;
}

void VelocityVerletIntegrator::set_barostat(std::shared_ptr<Barostat> barostat) noexcept {
    barostat_ = std::move(barostat);
}

void VelocityVerletIntegrator::set_target_pressure(double pressure) noexcept {
    target_pressure_ = pressure;
}

void VelocityVerletIntegrator::set_last_virial_trace(double virial_trace) noexcept {
    last_virial_trace_ = virial_trace;
}

}  // namespace gmd
