// Regression tests for the Nose-Hoover restart / degrees-of-freedom contract.
//
// The integrator computes the authoritative DOF from the global atom count, the
// centre-of-mass removal setting and the active constraints. A checkpoint also
// carries a `dof` field, and restoring it blindly would leave the thermostat
// integrating with one DOF while trajectory output reported temperature with
// another -- and would install a thermostat mass Q and friction variable xi
// that belong to a different extended system.
//
// The policy under test: the checkpoint's DOF is validated, never installed.
// Compatible checkpoints restore xi, Q, tau and the cached temperature exactly;
// incompatible ones are rejected with an explanatory error.

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

namespace {

int failures = 0;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[nh restart] " << message << '\n';
        ++failures;
    }
}

template <typename Body>
std::string capture_error(Body body) {
    try {
        body();
    } catch (const std::exception& error) {
        return error.what();
    }
    return {};
}

constexpr std::size_t kAtoms = 10;

gmd::System make_system() {
    gmd::System system;
    system.resize(kAtoms, kAtoms);
    gmd::Box box;
    box.set_lengths({20.0, 20.0, 20.0});
    system.set_box(box);
    for (std::size_t i = 0; i < kAtoms; ++i) {
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {1.5 * static_cast<double>(i), 5.0, 5.0};
        system.mutable_velocities()[i] = {0.002 * static_cast<double>(i + 1), -0.001, 0.0015};
    }
    return system;
}

// Chain of `count` constraints over consecutive atoms.
std::shared_ptr<gmd::ConstraintSolver> make_chain_constraints(std::size_t count) {
    std::vector<gmd::BondConstraint> constraints;
    for (std::size_t i = 0; i < count; ++i) {
        constraints.push_back(gmd::BondConstraint{static_cast<int>(i),
                                                  static_cast<int>(i + 1), 1.5});
    }
    return std::make_shared<gmd::ConstraintSolver>(std::move(constraints),
                                                   gmd::ConstraintSettings{});
}

// Builds an integrator + thermostat pair wired the way Simulation wires them.
struct Rig {
    std::shared_ptr<gmd::NoseHooverThermostat> thermostat;
    gmd::VelocityVerletIntegrator integrator{1.0};
};

std::unique_ptr<Rig> make_rig(gmd::System& system,
                              bool remove_com,
                              std::size_t constraint_count) {
    auto rig = std::make_unique<Rig>();
    rig->thermostat = std::make_shared<gmd::NoseHooverThermostat>(50.0);
    rig->integrator.set_thermostat(rig->thermostat);
    rig->integrator.set_target_temperature(300.0);
    rig->integrator.set_remove_center_of_mass_velocity(remove_com);
    if (constraint_count > 0) {
        rig->integrator.set_constraint_solver(make_chain_constraints(constraint_count));
    }
    gmd::RuntimeContext runtime;
    rig->integrator.initialize(system, runtime);
    return rig;
}

// Advances the thermostat so xi, Q and the cached temperature are non-trivial.
void spin_up(Rig& rig, gmd::System& system) {
    for (int i = 0; i < 4; ++i) {
        rig.thermostat->apply_half_kick(system, 0.5, 300.0);
    }
}

double field(const std::string& state, const std::string& name) {
    std::istringstream input(state);
    std::string key;
    double value = 0.0;
    while (input >> key) {
        if (key == name) {
            input >> value;
            return value;
        }
        std::string skip;
        input >> skip;
    }
    return std::nan("");
}

// --- Compatible restarts ---------------------------------------------------

void test_compatible_unconstrained_restart() {
    gmd::System system = make_system();
    auto original = make_rig(system, /*remove_com=*/true, /*constraints=*/0);
    check(original->thermostat->degrees_of_freedom() == 3 * kAtoms - 3,
          "unconstrained COM-removed DOF should be 3N-3");
    spin_up(*original, system);
    const std::string state = original->thermostat->checkpoint_state();

    gmd::System restart_system = make_system();
    auto restored = make_rig(restart_system, true, 0);
    const std::string error = capture_error([&] {
        restored->thermostat->load_checkpoint_state(state);
    });
    check(error.empty(), "compatible unconstrained restart was rejected: " + error);

    check(restored->thermostat->degrees_of_freedom() == 3 * kAtoms - 3,
          "restart must keep the authoritative DOF");
    check(restored->thermostat->checkpoint_state() == state,
          "compatible restart must reproduce the checkpoint state exactly");
}

void test_compatible_constrained_restart() {
    gmd::System system = make_system();
    auto original = make_rig(system, /*remove_com=*/true, /*constraints=*/9);
    check(original->thermostat->degrees_of_freedom() == 3 * kAtoms - 3 - 9,
          "constrained COM-removed DOF should be 3N-3-Nc");
    spin_up(*original, system);
    const std::string state = original->thermostat->checkpoint_state();

    gmd::System restart_system = make_system();
    auto restored = make_rig(restart_system, true, 9);
    const std::string error = capture_error([&] {
        restored->thermostat->load_checkpoint_state(state);
    });
    check(error.empty(), "compatible constrained restart was rejected: " + error);
    check(restored->thermostat->checkpoint_state() == state,
          "compatible constrained restart must reproduce the checkpoint state");
}

// --- Continuity of the extended-system state -------------------------------

void test_restart_continuity() {
    gmd::System system = make_system();
    auto original = make_rig(system, true, 9);
    spin_up(*original, system);

    const std::string state = original->thermostat->checkpoint_state();
    const double xi_before = field(state, "xi");
    const double q_before = field(state, "Q");
    const double temperature_before = field(state, "current_temperature");
    const double tau_before = field(state, "tau");

    check(std::abs(xi_before) > 0.0, "xi should be non-zero after spin-up");
    check(q_before > 0.0, "Q should be set after spin-up");
    check(temperature_before > 0.0, "temperature should be cached after spin-up");

    // The restarted run resumes from the system state as it stood at the
    // checkpoint, which is what read_checkpoint() restores alongside the
    // thermostat state. Starting from a fresh system would compare two
    // different trajectories rather than testing continuity.
    gmd::System restart_system = system;
    auto restored = make_rig(restart_system, true, 9);
    restored->thermostat->load_checkpoint_state(state);

    const std::string after = restored->thermostat->checkpoint_state();
    check(field(after, "xi") == xi_before, "xi did not survive the restart");
    check(field(after, "Q") == q_before, "Q did not survive the restart");
    check(field(after, "tau") == tau_before, "tau did not survive the restart");
    check(field(after, "current_temperature") == temperature_before,
          "cached temperature did not survive the restart");

    // The restarted thermostat must now produce the same trajectory as the
    // original would have: same velocities after the same further half-kicks.
    for (int i = 0; i < 3; ++i) {
        original->thermostat->apply_half_kick(system, 0.5, 300.0);
        restored->thermostat->apply_half_kick(restart_system, 0.5, 300.0);
    }

    double worst = 0.0;
    for (std::size_t i = 0; i < kAtoms; ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            worst = std::max(worst, std::abs(system.velocities()[i][d] -
                                             restart_system.velocities()[i][d]));
        }
    }
    check(worst == 0.0,
          "restarted trajectory diverged from the continued one (max velocity "
          "difference " + std::to_string(worst) + ")");
}

// --- Rejected restarts -----------------------------------------------------

// Writes a checkpoint state with an arbitrary dof, standing in for a checkpoint
// produced under a different configuration or by an older version.
std::string synthetic_state(double tau, double xi, double q, long long dof,
                            double temperature) {
    std::ostringstream out;
    out.precision(17);
    out << "tau " << tau << " xi " << xi << " Q " << q << " dof " << dof
        << " current_temperature " << temperature;
    return out.str();
}

void test_constraint_change_is_rejected() {
    gmd::System system = make_system();
    auto original = make_rig(system, true, 9);
    spin_up(*original, system);
    const std::string state = original->thermostat->checkpoint_state();

    // Same system, same COM setting, but a different number of constraints.
    gmd::System restart_system = make_system();
    auto restored = make_rig(restart_system, true, 4);
    const std::string error = capture_error([&] {
        restored->thermostat->load_checkpoint_state(state);
    });
    check(!error.empty(),
          "a checkpoint whose DOF disagrees with the constraint configuration "
          "must be rejected");
    check(error.find("degrees of freedom") != std::string::npos,
          "constraint-mismatch error should mention degrees of freedom, got: " + error);
    check(error.find("18") != std::string::npos && error.find("23") != std::string::npos,
          "error should report both the checkpoint and the current DOF, got: " + error);
}

void test_com_setting_change_is_rejected() {
    gmd::System system = make_system();
    auto original = make_rig(system, /*remove_com=*/true, 0);
    spin_up(*original, system);
    const std::string state = original->thermostat->checkpoint_state();

    gmd::System restart_system = make_system();
    auto restored = make_rig(restart_system, /*remove_com=*/false, 0);
    const std::string error = capture_error([&] {
        restored->thermostat->load_checkpoint_state(state);
    });
    check(!error.empty(),
          "changing the centre-of-mass removal setting must be rejected");
    check(error.find("centre-of-mass") != std::string::npos,
          "COM-mismatch error should mention the centre-of-mass setting, got: " + error);
}

void test_legacy_checkpoint_policy() {
    gmd::System system = make_system();

    // A legacy checkpoint is one whose dof holds the old 3N-3 value. Where that
    // happens to equal the authoritative count -- unconstrained, COM removed --
    // it is accepted transparently.
    {
        const std::string legacy =
            synthetic_state(50.0, 0.01, 1.0e-3, 3 * kAtoms - 3, 295.0);
        auto restored = make_rig(system, true, 0);
        const std::string error = capture_error([&] {
            restored->thermostat->load_checkpoint_state(legacy);
        });
        check(error.empty(),
              "a legacy 3N-3 checkpoint matching the current DOF must be accepted: " +
                  error);
        check(restored->thermostat->degrees_of_freedom() == 3 * kAtoms - 3,
              "accepted legacy checkpoint must keep the authoritative DOF");
    }

    // With constraints active the old rule no longer agrees, and the run cannot
    // be continued: Q and xi were generated for a different extended system.
    {
        const std::string legacy =
            synthetic_state(50.0, 0.01, 1.0e-3, 3 * kAtoms - 3, 295.0);
        auto restored = make_rig(system, true, 9);
        const std::string error = capture_error([&] {
            restored->thermostat->load_checkpoint_state(legacy);
        });
        check(!error.empty(),
              "a legacy 3N-3 checkpoint must be rejected once constraints are active");
        check(error.find("3N-3") != std::string::npos,
              "legacy rejection should name the old 3N-3 rule as a likely cause, got: " +
                  error);
    }
}

void test_invalid_dof_is_rejected() {
    gmd::System system = make_system();

    for (long long bad_dof : {0LL, -5LL}) {
        auto restored = make_rig(system, true, 0);
        const std::string error = capture_error([&] {
            restored->thermostat->load_checkpoint_state(
                synthetic_state(50.0, 0.01, 1.0e-3, bad_dof, 295.0));
        });
        check(!error.empty(),
              "a checkpoint recording dof=" + std::to_string(bad_dof) +
                  " must be rejected");
        check(error.find("positive") != std::string::npos,
              "invalid-dof error should ask for a positive count, got: " + error);
    }
}

void test_restore_before_dof_is_known_is_rejected() {
    gmd::System system = make_system();
    gmd::NoseHooverThermostat bare(50.0);
    bare.initialize(system);   // installs only the provisional default

    const std::string error = capture_error([&] {
        bare.load_checkpoint_state(
            synthetic_state(50.0, 0.01, 1.0e-3, 3 * kAtoms - 3, 295.0));
    });
    check(!error.empty(),
          "restoring before the authoritative DOF is installed must be rejected, "
          "even when the checkpoint value happens to match the provisional one");
    check(error.find("authoritative") != std::string::npos,
          "ordering error should explain the authoritative-DOF requirement, got: " +
              error);
}

void test_malformed_state_still_rejected() {
    gmd::System system = make_system();
    auto restored = make_rig(system, true, 0);
    check(!capture_error([&] {
              restored->thermostat->load_checkpoint_state("tau 1.0 xi 2.0");
          }).empty(),
          "a truncated thermostat state must still be rejected");

    // Empty and "stateless" remain no-ops for backward compatibility.
    check(capture_error([&] {
              restored->thermostat->load_checkpoint_state("stateless");
          }).empty(),
          "\"stateless\" must remain accepted as a no-op");
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif

    test_compatible_unconstrained_restart();
    test_compatible_constrained_restart();
    test_restart_continuity();
    test_constraint_change_is_rejected();
    test_com_setting_change_is_rejected();
    test_legacy_checkpoint_policy();
    test_invalid_dof_is_rejected();
    test_restore_before_dof_is_known_is_rejected();
    test_malformed_state_still_rejected();

    if (failures != 0) {
        std::cerr << "[nh restart] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[nh restart] all checks passed\n";
    return 0;
}
