#include "gmd/integrator/velocity_verlet_integrator.hpp"

#include <algorithm>
#include <stdexcept>

#include "gmd/system/box.hpp"
#include "gmd/system/periodic_boundary.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/barostat.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/core/runtime_context.hpp"
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

    // Put the initial state ON the constraint manifold, in position and in
    // velocity, before any dynamics run. A state that satisfies neither is not a
    // constrained state, and projecting it on the first step would make that
    // step's multipliers a one-off correction of an invalid state rather than a
    // constraint force. A restarted run is already on the manifold, so this is a
    // no-op there and restart continuity is unaffected.
    if (has_constraints()) {
        apply_position_constraints(system);
        apply_velocity_constraints(system);
    }

    auto forces = system.mutable_forces();
    for (auto& force : forces) {
        force = {0.0, 0.0, 0.0};
    }
    if (thermostat_) {
        thermostat_->initialize(system);
        // initialize() only knows the atom count, so it installs the
        // unconstrained/COM-removed default. Override it with the count that
        // also accounts for constraints and for the run's COM-removal setting.
        thermostat_->set_degrees_of_freedom(degrees_of_freedom(system));
    }
}

void VelocityVerletIntegrator::set_remove_center_of_mass_velocity(bool enabled) noexcept {
    remove_center_of_mass_velocity_ = enabled;
}

std::size_t VelocityVerletIntegrator::constraint_count() const noexcept {
    // ConstraintSolver normalises its list to distinct pairs and replicates the
    // same tag-based list on every rank, so this is already a global count and
    // must not be reduced again. It counts *distinct* constraints; independence
    // is assumed rather than proven (see the ConstraintSolver class comment).
    return has_constraints() ? constraints_->active_constraint_count() : 0;
}

std::size_t VelocityVerletIntegrator::degrees_of_freedom(const System& system) const noexcept {
    return compute_degrees_of_freedom(
        system,
        DegreesOfFreedomConfig{remove_center_of_mass_velocity_, constraint_count()});
}

void VelocityVerletIntegrator::step(System& system,
                                    ForceProvider& force_provider,
                                    const IntegratorStepContext& ctx,
                                    RuntimeContext& runtime) {
    const double force_time = (ctx.step + 1) * (ctx.dt > 0.0 ? ctx.dt : dt_);

    begin_step(system, ctx);
    ForceResult next_force = evaluate_force(system,
                                            force_provider,
                                            ctx.step + 1,
                                            force_time,
                                            runtime);
    copy_forces_to_system(system, next_force);
    finish_step(system, ctx, next_force.virial_valid, next_force.virial);

    // A barostat rescales the box and the coordinates *after* the forces above
    // were computed, which would leave the system holding forces for the
    // pre-scaling geometry. Detect an accepted volume change and refresh.
    const Box box_before = system.box();
    apply_barostat(system, force_provider, runtime, ctx);
    if (!box_lengths_equal(box_before, system.box())) {
        refresh_after_barostat(system, force_provider, runtime, ctx.step + 1, force_time);
    }
}

void VelocityVerletIntegrator::refresh_after_barostat(System& system,
                                                      ForceProvider& force_provider,
                                                      RuntimeContext& runtime,
                                                      std::uint64_t force_step,
                                                      double force_time) {
    // Constraints act on the rescaled coordinates, then the stale neighbor list
    // is dropped so the provider rebuilds against the new box, and only then are
    // forces recomputed for the geometry the caller will actually see.
    apply_position_constraints(system);
    system.mutable_neighbor_list().valid = false;

    ForceResult rescaled = evaluate_force(system,
                                          force_provider,
                                          force_step,
                                          force_time,
                                          runtime);
    copy_forces_to_system(system, rescaled);
    system.set_last_virial(rescaled.virial, rescaled.virial_valid);
    last_virial_valid_ = rescaled.virial_valid;
    if (last_virial_valid_) {
        last_virial_trace_ = rescaled.virial[0] + rescaled.virial[4] + rescaled.virial[8];
    }
}

void VelocityVerletIntegrator::begin_step(System& system,
                                          const IntegratorStepContext& ctx) {
    const double dt = ctx.dt > 0.0 ? ctx.dt : dt_;
    if (dt <= 0.0) {
        throw std::runtime_error("VelocityVerletIntegrator requires a positive time step");
    }

    // --- Thermostat pre-kick (Nosé-Hoover first half-kick or no-op) ---
    if (thermostat_) {
        thermostat_->apply_half_kick(system, 0.5 * dt, target_temperature_);
    }

    // The constraint geometry the step starts from. Standard SHAKE corrects
    // along these gradients, so they must be taken BEFORE the drift moves the
    // atoms. Collective under MPI, and replicated, like the projection itself.
    ConstraintReference reference;
    const bool constrained = has_constraints();
    if (constrained) {
        reference = constraints_->capture_reference(system);
    }

    const auto masses = system.masses();
    const auto forces = system.forces();
    auto coordinates = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();

    // First half-kick (v += 0.5*a*dt) + full position update.
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        if (masses[atom_index] <= 0.0) {
            throw std::runtime_error("VelocityVerletIntegrator requires strictly positive masses");
        }

        const double inverse_mass = 1.0 / masses[atom_index];
        for (std::size_t dim = 0; dim < 3; ++dim) {
            velocities[atom_index][dim] += 0.5 * forces[atom_index][dim] * inverse_mass * dt;
            coordinates[atom_index][dim] += velocities[atom_index][dim] * dt;
        }
        wrap_position(coordinates[atom_index], system.box());
    }

    // Steps (2) and (3) of the SHAKE/RATTLE splitting: solve for the constraint
    // multipliers at time level t and apply BOTH the position correction and its
    // matching half-step velocity impulse dr/dt.
    if (constrained) {
        system.set_last_shake_stats(constraints_->apply_shake(system, reference, dt));
    }
}

void VelocityVerletIntegrator::finish_step(System& system,
                                           const IntegratorStepContext& ctx,
                                           bool virial_valid,
                                           const std::array<double, 9>& virial) {
    const double dt = ctx.dt > 0.0 ? ctx.dt : dt_;
    if (dt <= 0.0) {
        throw std::runtime_error("VelocityVerletIntegrator requires a positive time step");
    }

    // Cache virial trace for barostat (uses full virial tensor if available).
    last_virial_valid_ = virial_valid;
    if (last_virial_valid_) {
        last_virial_trace_ = virial[0] + virial[4] + virial[8];
    }

    const auto masses = system.masses();
    const auto forces = system.forces();
    auto velocities = system.mutable_velocities();

    // Second half-kick.
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        const double inverse_mass = 1.0 / masses[atom_index];
        for (std::size_t dim = 0; dim < 3; ++dim) {
            velocities[atom_index][dim] += 0.5 * forces[atom_index][dim] * inverse_mass * dt;
        }
    }

    // --- Thermostat post-kick (Nosé-Hoover second half-kick, or full rescaling step) ---
    if (thermostat_) {
        thermostat_->apply_half_kick(system, 0.5 * dt, target_temperature_);
        thermostat_->apply(system, dt, target_temperature_);
    }

    apply_velocity_constraints(system);
}

void VelocityVerletIntegrator::apply_barostat(System& system,
                                              ForceProvider& force_provider,
                                              RuntimeContext& runtime,
                                              const IntegratorStepContext& ctx) {
    const double dt = ctx.dt > 0.0 ? ctx.dt : dt_;
    if (dt <= 0.0) {
        throw std::runtime_error("VelocityVerletIntegrator requires a positive time step");
    }

    if (barostat_) {
        const bool can_run = !barostat_->requires_virial() || last_virial_valid_;
        if (can_run) {
            barostat_->apply(system, force_provider, runtime,
                             ctx.step, dt, target_temperature_,
                             target_pressure_, last_virial_trace_);
        }
    }
}

void VelocityVerletIntegrator::apply_position_constraints(System& system) {
    if (constraints_ != nullptr && constraints_->enabled()) {
        system.set_last_shake_stats(constraints_->apply_shake(system));
    }
}

void VelocityVerletIntegrator::apply_velocity_constraints(System& system) {
    if (constraints_ != nullptr && constraints_->enabled()) {
        system.set_last_rattle_stats(constraints_->apply_rattle(system));
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

void VelocityVerletIntegrator::set_constraint_solver(
        std::shared_ptr<ConstraintSolver> constraints) noexcept {
    constraints_ = std::move(constraints);
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
