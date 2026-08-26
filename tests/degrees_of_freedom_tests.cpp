// Regression tests for the shared degrees-of-freedom calculation.
//
// Before this suite existed, every temperature consumer hardcoded 3N-3: the
// centre-of-mass setting was ignored and SHAKE/RATTLE constraints were never
// subtracted, so a constrained run reported a temperature that was too low and
// the thermostat compensated by heating the system past its target.

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/integrator/velocity_rescaling_thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

namespace {

int failures = 0;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[dof] " << message << '\n';
        ++failures;
    }
}

void check_eq(std::size_t actual, std::size_t expected, const std::string& message) {
    if (actual != expected) {
        std::cerr << "[dof] " << message << ": expected " << expected
                  << ", got " << actual << '\n';
        ++failures;
    }
}

gmd::System make_system(std::size_t atom_count) {
    gmd::System system;
    system.resize(atom_count, atom_count);
    gmd::Box box;
    box.set_lengths({20.0, 20.0, 20.0});
    system.set_box(box);
    for (std::size_t i = 0; i < atom_count; ++i) {
        system.mutable_masses()[i] = 1.0;
        system.mutable_coordinates()[i] = {static_cast<double>(i), 0.0, 0.0};
    }
    return system;
}

// A chain of `atom_count` atoms with every consecutive bond constrained.
std::shared_ptr<gmd::ConstraintSolver> make_constraints(std::size_t atom_count) {
    gmd::Topology topology;
    for (std::size_t i = 0; i + 1 < atom_count; ++i) {
        topology.bonds.push_back({static_cast<int>(i), static_cast<int>(i + 1), 0});
    }
    // bond.type_idx is 0 for every bond above, so type 0 is the constrained one.
    auto constraints = gmd::constraints_from_bond_types(topology, {0}, {1.0});
    return std::make_shared<gmd::ConstraintSolver>(std::move(constraints),
                                                   gmd::ConstraintSettings{});
}

// --- The formula itself ---------------------------------------------------

void test_formula() {
    using gmd::compute_degrees_of_freedom;
    using Config = gmd::DegreesOfFreedomConfig;

    // Unconstrained, COM removed: 3N - 3.
    check_eq(compute_degrees_of_freedom(10, Config{true, 0}), 27,
             "10 atoms, COM removed");

    // Unconstrained, COM kept: the full 3N.
    check_eq(compute_degrees_of_freedom(10, Config{false, 0}), 30,
             "10 atoms, COM kept");

    // Constrained: one DOF per constraint on top of the COM removal.
    check_eq(compute_degrees_of_freedom(10, Config{true, 9}), 18,
             "10 atoms, COM removed, 9 constraints");
    check_eq(compute_degrees_of_freedom(10, Config{false, 9}), 21,
             "10 atoms, COM kept, 9 constraints");

    // A single atom keeps all three translations: there is no relative motion
    // for the COM subtraction to be meaningful against.
    check_eq(compute_degrees_of_freedom(1, Config{true, 0}), 3, "1 atom, COM removed");
    check_eq(compute_degrees_of_freedom(1, Config{false, 0}), 3, "1 atom, COM kept");

    // Two atoms with COM removal: 6 - 3 = 3.
    check_eq(compute_degrees_of_freedom(2, Config{true, 0}), 3, "2 atoms, COM removed");
    // ... and a constraint between them leaves 2 rotational DOF.
    check_eq(compute_degrees_of_freedom(2, Config{true, 1}), 2,
             "2 atoms, COM removed, 1 constraint");

    // Degenerate and over-constrained systems saturate at zero rather than
    // wrapping around through unsigned underflow.
    check_eq(compute_degrees_of_freedom(0, Config{true, 0}), 0, "0 atoms");
    check_eq(compute_degrees_of_freedom(0, Config{false, 0}), 0, "0 atoms, COM kept");
    check_eq(compute_degrees_of_freedom(2, Config{true, 3}), 0,
             "2 atoms, COM removed, 3 constraints (exactly zero)");
    check_eq(compute_degrees_of_freedom(2, Config{true, 99}), 0,
             "2 atoms over-constrained must saturate at 0, not underflow");
    check_eq(compute_degrees_of_freedom(1, Config{true, 5}), 0,
             "1 atom over-constrained must saturate at 0, not underflow");
}

// --- The integrator's authoritative count ---------------------------------

void test_integrator_dof() {
    gmd::System system = make_system(10);

    {   // Unconstrained, COM removed (the default).
        gmd::VelocityVerletIntegrator integrator(1.0);
        integrator.set_remove_center_of_mass_velocity(true);
        check_eq(integrator.degrees_of_freedom(system), 27,
                 "integrator: unconstrained, COM removed");
        check_eq(integrator.constraint_count(), 0, "integrator: no constraints");
    }

    {   // Unconstrained, COM kept.
        gmd::VelocityVerletIntegrator integrator(1.0);
        integrator.set_remove_center_of_mass_velocity(false);
        check_eq(integrator.degrees_of_freedom(system), 30,
                 "integrator: unconstrained, COM kept");
    }

    {   // Constrained chain: 9 bonds between 10 atoms.
        gmd::VelocityVerletIntegrator integrator(1.0);
        integrator.set_remove_center_of_mass_velocity(true);
        integrator.set_constraint_solver(make_constraints(10));
        check_eq(integrator.constraint_count(), 9, "integrator: 9 chain constraints");
        check_eq(integrator.degrees_of_freedom(system), 18,
                 "integrator: constrained, COM removed");

        integrator.set_remove_center_of_mass_velocity(false);
        check_eq(integrator.degrees_of_freedom(system), 21,
                 "integrator: constrained, COM kept");
    }
}

// --- The thermostats must receive that same count -------------------------

void test_thermostat_receives_dof() {
    gmd::RuntimeContext runtime;
    gmd::System system = make_system(10);

    {   // Velocity rescaling, constrained, COM removed.
        auto thermostat = std::make_shared<gmd::VelocityRescalingThermostat>();
        gmd::VelocityVerletIntegrator integrator(1.0);
        integrator.set_thermostat(thermostat);
        integrator.set_remove_center_of_mass_velocity(true);
        integrator.set_constraint_solver(make_constraints(10));
        integrator.initialize(system, runtime);

        check_eq(thermostat->degrees_of_freedom(), 18,
                 "velocity rescaling thermostat DOF after initialize");
        check_eq(thermostat->degrees_of_freedom(),
                 integrator.degrees_of_freedom(system),
                 "velocity rescaling thermostat agrees with integrator");
    }

    {   // Nose-Hoover, unconstrained, COM kept.
        auto thermostat = std::make_shared<gmd::NoseHooverThermostat>(100.0);
        gmd::VelocityVerletIntegrator integrator(1.0);
        integrator.set_thermostat(thermostat);
        integrator.set_remove_center_of_mass_velocity(false);
        integrator.initialize(system, runtime);

        check_eq(thermostat->degrees_of_freedom(), 30,
                 "Nose-Hoover thermostat DOF after initialize (COM kept)");
        check_eq(thermostat->degrees_of_freedom(),
                 integrator.degrees_of_freedom(system),
                 "Nose-Hoover thermostat agrees with integrator");
    }

    {   // Nose-Hoover, constrained: this is the case the old 3N-3 formula got
        // wrong, and the one where the error changes the sampled temperature.
        auto thermostat = std::make_shared<gmd::NoseHooverThermostat>(100.0);
        gmd::VelocityVerletIntegrator integrator(1.0);
        integrator.set_thermostat(thermostat);
        integrator.set_remove_center_of_mass_velocity(true);
        integrator.set_constraint_solver(make_constraints(10));
        integrator.initialize(system, runtime);

        check_eq(thermostat->degrees_of_freedom(), 18,
                 "Nose-Hoover thermostat DOF with constraints");
        check(thermostat->degrees_of_freedom() != 27,
              "Nose-Hoover thermostat must not fall back to the old 3N-3 value");
    }
}

// The reported temperature scales as 1/dof, so an incorrect DOF count is
// directly an incorrect temperature. Pin that relationship down.
void test_temperature_uses_dof() {
    const double twice_ke = 1.0;
    const double t_27 = gmd::temperature_from_twice_ke(twice_ke, 27);
    const double t_18 = gmd::temperature_from_twice_ke(twice_ke, 18);

    check(t_18 > t_27, "fewer DOF must report a higher temperature");
    check(std::abs(t_18 / t_27 - 27.0 / 18.0) < 1.0e-12,
          "temperature must scale exactly as 1/dof");

    // Zero DOF is explicitly defined as "temperature undefined", reported as
    // zero rather than a division by zero.
    check(gmd::temperature_from_twice_ke(twice_ke, 0) == 0.0,
          "zero DOF must report zero temperature, not a division by zero");
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    // An MPI-enabled build still runs this as a single-process test, but the
    // force providers issue collectives that require an initialised MPI.
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif
    test_formula();
    test_integrator_dof();
    test_thermostat_receives_dof();
    test_temperature_uses_dof();

    if (failures != 0) {
        std::cerr << "[dof] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[dof] all checks passed\n";
    return 0;
}
