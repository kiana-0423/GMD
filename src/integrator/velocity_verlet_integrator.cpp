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
    // Tell the System whether a constraint term is required in the reported
    // virial. With constraints on, the state starts Unavailable: no step has run
    // yet, so no constraint multiplier belongs to the initial provider virial and
    // the initial pressure must not be reported as complete. A restart attaches
    // the checkpointed value afterwards (see the restart handling in gmd_main).
    system.set_constraints_active(has_constraints());

    // Put the initial state ON the constraint manifold, in position and in
    // velocity, before any dynamics run. Without this the first RATTLE of the
    // first step has to absorb the whole initial violation, and the multiplier it
    // converges to is a one-off correction of an invalid state rather than a
    // constraint force -- which would show up as a spurious pressure spike on
    // step one. dt is deliberately omitted: this projection closes no step and
    // must not produce a virial. A restarted run is already on the manifold, so
    // this is a no-op there and restart continuity is unaffected.
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
    // The provider half of the reported virial, at r(t+dt). Its constraint
    // partner is not known yet: it comes from the RATTLE inside finish_step,
    // which runs on these same coordinates.
    system.set_provider_virial(next_force.virial, next_force.virial_valid);
    finish_step(system, ctx);

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
    //
    // Both projections are geometric: rescaling the cell moves the atoms, so the
    // bonds are off their targets and the velocities are no longer tangent to
    // them. dt is deliberately omitted from both -- these corrections close no
    // step, so no constraint virial is recovered from them.
    apply_position_constraints(system);
    apply_velocity_constraints(system);
    system.mutable_neighbor_list().valid = false;

    ForceResult rescaled = evaluate_force(system,
                                          force_provider,
                                          force_step,
                                          force_time,
                                          runtime);
    copy_forces_to_system(system, rescaled);
    // The provider virial now belongs to the rescaled geometry while the step's
    // constraint multipliers belong to the geometry before the rescale. There is
    // no contemporaneous constraint term for this state, and set_provider_virial
    // drops the stale one rather than combining across the rescale. With
    // constraints active the resulting virial is therefore reported INVALID: a
    // provider-only tensor is not a complete pressure virial when constraints
    // act. The next step's RATTLE restores a valid one.
    system.set_provider_virial(rescaled.virial, rescaled.virial_valid);
    last_virial_valid_ = system.last_virial_valid();
    if (last_virial_valid_) {
        const auto& combined = system.last_virial();
        last_virial_trace_ = combined[0] + combined[4] + combined[8];
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
    // matching half-step velocity impulse dr/dt. The multipliers here pair with
    // F(t); the ones the reported virial needs come from RATTLE at t+dt.
    if (constrained) {
        system.set_last_shake_stats(constraints_->apply_shake(system, reference, dt));
    }
}

void VelocityVerletIntegrator::finish_step(System& system,
                                           const IntegratorStepContext& ctx) {
    const double dt = ctx.dt > 0.0 ? ctx.dt : dt_;
    if (dt <= 0.0) {
        throw std::runtime_error("VelocityVerletIntegrator requires a positive time step");
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

    // Dynamical projection closing the step. Its multipliers are the constraint
    // partners of the forces evaluated at these coordinates, so this is where the
    // step's constraint virial comes from and where the combined tensor is
    // completed. Nothing may read the reported virial between the force
    // evaluation and this call.
    apply_velocity_constraints(system, dt);

    // Cache the COMBINED trace for the barostat, read back off the System so the
    // barostat, the trajectory log and the checkpoint cannot disagree.
    last_virial_valid_ = system.last_virial_valid();
    if (last_virial_valid_) {
        const auto& combined = system.last_virial();
        last_virial_trace_ = combined[0] + combined[4] + combined[8];
    }

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
    if (constraints_ == nullptr || !constraints_->enabled()) {
        return;
    }
    system.set_last_shake_stats(constraints_->apply_shake(system));
}

void VelocityVerletIntegrator::apply_velocity_constraints(System& system, double dt) {
    if (constraints_ == nullptr || !constraints_->enabled()) {
        // No constraints in this run: the provider virial is complete on its own.
        system.set_constraints_active(false);
        return;
    }
    system.set_constraints_active(true);

    if (dt > 0.0) {
        ConstraintVirialResult constraint_virial;
        system.set_last_rattle_stats(constraints_->apply_rattle(system, dt, constraint_virial));
        if (constraint_virial.valid) {
            system.set_constraint_virial(constraint_virial.virial);
        } else {
            // RATTLE disabled: constraints act but no multiplier is available, so
            // the reported virial stays incomplete rather than silently short.
            system.mark_constraint_virial_unavailable();
        }
        return;
    }

    system.set_last_rattle_stats(constraints_->apply_rattle(system));
    system.mark_constraint_virial_unavailable();
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
