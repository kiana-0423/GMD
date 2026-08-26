#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/integrator.hpp"

namespace gmd {

class Thermostat;
class Barostat;

// Default second-order integrator for deterministic molecular dynamics.
// Optionally supports:
//   - Velocity-Verlet Nosé-Hoover (VVNH): thermostat half-kicks wrap the VV steps.
//   - Velocity-rescaling thermostat: applied once per step after the second half-kick.
//   - Berendsen barostat: applied once per step at the end.
class VelocityVerletIntegrator final : public Integrator {
public:
    explicit VelocityVerletIntegrator(double dt) noexcept;
    ~VelocityVerletIntegrator() override = default;

    std::string_view name() const noexcept override;

    void initialize(System& system, RuntimeContext& runtime) override;
    void step(System& system,
              ForceProvider& force_provider,
              const IntegratorStepContext& ctx,
              RuntimeContext& runtime) override;
    void begin_step(System& system, const IntegratorStepContext& ctx);
    // Second half-kick, thermostat, then RATTLE.
    //
    // RATTLE is where the step's constraint virial comes from, so the combined
    // virial is assembled and cached only after it has run. The provider virial
    // must already have been installed on the System (System::set_provider_virial)
    // before this is called; the tensor this caches for the barostat is the
    // combined one, read back off the System.
    void finish_step(System& system, const IntegratorStepContext& ctx);
    void apply_barostat(System& system,
                        ForceProvider& force_provider,
                        RuntimeContext& runtime,
                        const IntegratorStepContext& ctx);
    // Geometric projection of positions onto the constraint manifold. Not a
    // step of the dynamics -- for initial and barostat-rescaled states only.
    // The dynamical SHAKE lives inside begin_step(), which has the reference
    // geometry the projection must correct along.
    void apply_position_constraints(System& system);

    // Velocity projection (RATTLE).
    //
    // `dt > 0` marks the dynamical projection that closes a step, so the
    // ENDPOINT constraint virial is recovered from the converged multipliers and
    // attached to the provider virial already installed on the System. `dt == 0`
    // marks a projection that is not closing a step -- the one that puts an
    // initial or rescaled state onto the constraint manifold -- whose
    // multipliers are a one-off correction of an invalid state rather than a
    // constraint force; that form leaves the constraint virial unavailable.
    void apply_velocity_constraints(System& system, double dt = 0.0);

    // Returns the configured default time step for this integrator instance.
    double dt() const noexcept;
    void set_dt(double dt) noexcept;

    // --- Thermostat ---
    void set_thermostat(std::shared_ptr<Thermostat> thermostat) noexcept;
    void set_target_temperature(double temperature) noexcept;
    double target_temperature() const noexcept { return target_temperature_; }
    void set_constraint_solver(std::shared_ptr<ConstraintSolver> constraints) noexcept;
    bool has_constraints() const noexcept { return constraints_ != nullptr && constraints_->enabled(); }

    // --- Degrees of freedom ---
    // Whether the run removes the centre-of-mass velocity. This only affects
    // the DOF count; the removal itself is performed by VelocityInitializer.
    // Simulation forwards its own setting here during initialize().
    void set_remove_center_of_mass_velocity(bool enabled) noexcept;
    bool remove_center_of_mass_velocity() const noexcept {
        return remove_center_of_mass_velocity_;
    }

    // Number of distinct active holonomic constraints (0 when disabled).
    // Independence is assumed, not verified; see ConstraintSolver.
    std::size_t constraint_count() const noexcept;

    // Authoritative DOF count for this run: 3N, less 3 for COM removal, less
    // one per distinct constraint. Every temperature consumer must use this.
    std::size_t degrees_of_freedom(const System& system) const noexcept;

    // --- Barostat ---
    void set_barostat(std::shared_ptr<Barostat> barostat) noexcept;
    void set_target_pressure(double pressure) noexcept;
    double target_pressure() const noexcept { return target_pressure_; }
    bool has_barostat() const noexcept { return barostat_ != nullptr; }

    // Last virial trace computed during force evaluation (used for barostat).
    // Only valid when virial_valid is true in the last ForceResult.
    void set_last_virial_trace(double virial_trace) noexcept;

private:
    // Re-establish forces, virial and neighbor-list state after a barostat has
    // rescaled the box and coordinates.
    void refresh_after_barostat(System& system,
                                ForceProvider& force_provider,
                                RuntimeContext& runtime,
                                std::uint64_t force_step,
                                double force_time);

    double dt_ = 0.0;
    double target_temperature_ = 300.0;   // [K]
    double target_pressure_    = 1.0;     // user-defined pressure units
    double last_virial_trace_  = 0.0;
    bool last_virial_valid_    = false;
    bool remove_center_of_mass_velocity_ = true;

    std::shared_ptr<Thermostat> thermostat_;
    std::shared_ptr<Barostat>   barostat_;
    std::shared_ptr<ConstraintSolver> constraints_;
};

}  // namespace gmd
