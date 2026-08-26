// Regression tests for constraint-list normalisation and validation.
//
// The degrees-of-freedom calculation subtracts one per entry in the constraint
// list, which is only meaningful if that list holds distinct constraints. These
// tests pin down what ConstraintSolver accepts, what it de-duplicates and what
// it rejects -- and document the one case it cannot decide: redundancy in a
// closed constraint topology.

#include <cstddef>
#include <iostream>
#include <memory>
#include <limits>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

namespace {

int failures = 0;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[constraints] " << message << '\n';
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

gmd::ConstraintSolver make_solver(std::vector<gmd::BondConstraint> constraints) {
    return gmd::ConstraintSolver(std::move(constraints), gmd::ConstraintSettings{});
}

gmd::System make_system(std::size_t atom_count) {
    gmd::System system;
    system.resize(atom_count, atom_count);
    gmd::Box box;
    box.set_lengths({30.0, 30.0, 30.0});
    system.set_box(box);
    for (std::size_t i = 0; i < atom_count; ++i) {
        system.mutable_masses()[i] = 1.0;
    }
    return system;
}

// --- Normalisation ---------------------------------------------------------

void test_reversed_duplicate_is_one_constraint() {
    auto solver = make_solver({{0, 1, 1.5}, {1, 0, 1.5}});
    check(solver.active_constraint_count() == 1,
          "(i,j) and (j,i) must count as one constraint, got " +
              std::to_string(solver.active_constraint_count()));
    check(solver.dropped_duplicate_count() == 1,
          "the reversed repeat should be reported as a dropped duplicate");
    check(solver.normalization_diagnostics().exact_duplicates == 1,
          "a reversed repeat with an identical target is an EXACT duplicate");
    check(solver.normalization_diagnostics().tolerance_equivalent_duplicates == 0,
          "a bit-identical target must not be classed as tolerance-equivalent");
    check(solver.normalization_diagnostics().discarded_targets.empty(),
          "an identical target is not a discarded target");

    // And it must be stored in normalised (min, max) order.
    check(solver.constraints().size() == 1 &&
              solver.constraints()[0].i == 0 && solver.constraints()[0].j == 1,
          "constraint should be stored as (min, max)");
}

void test_bit_identical_duplicates() {
    auto solver = make_solver({{2, 5, 1.0}, {2, 5, 1.0}, {2, 5, 1.0}});
    check(solver.active_constraint_count() == 1,
          "bit-identical duplicates must collapse to one constraint, got " +
              std::to_string(solver.active_constraint_count()));
    check(solver.normalization_diagnostics().exact_duplicates == 2,
          "two of the three entries should be counted as exact duplicates");
    check(solver.normalization_diagnostics().tolerance_equivalent_duplicates == 0,
          "bit-identical repeats must not be counted as tolerance-equivalent");
}

// --- Tolerance-equivalent duplicates and the tolerance boundary -----------

gmd::ConstraintSolver make_solver_with_tolerance(
        std::vector<gmd::BondConstraint> constraints, double tolerance) {
    gmd::ConstraintSettings settings;
    settings.tolerance = tolerance;
    return gmd::ConstraintSolver(std::move(constraints), settings);
}

void test_tolerance_equivalent_duplicate() {
    constexpr double tolerance = 1.0e-6;
    auto solver = make_solver_with_tolerance(
        {{0, 1, 1.5}, {0, 1, 1.5 + 1.0e-9}}, tolerance);

    check(solver.active_constraint_count() == 1,
          "a repeat within tolerance must collapse to one constraint");
    check(solver.normalization_diagnostics().tolerance_equivalent_duplicates == 1,
          "a non-identical target within tolerance is a tolerance-equivalent duplicate");
    check(solver.normalization_diagnostics().exact_duplicates == 0,
          "a non-identical target must not be counted as an exact duplicate");

    // First value wins, deterministically.
    check(solver.constraints()[0].target_distance == 1.5,
          "the first target distance must be the one kept");

    // The diagnostic must name both distances and the tolerance.
    const auto& discarded = solver.normalization_diagnostics().discarded_targets;
    check(discarded.size() == 1, "one discarded target should be recorded");
    if (discarded.size() == 1) {
        check(discarded[0].find("1.500000") != std::string::npos,
              "diagnostic should report the kept distance, got: " + discarded[0]);
        check(discarded[0].find("discarded") != std::string::npos,
              "diagnostic should report the discarded distance, got: " + discarded[0]);
        check(discarded[0].find("tolerance") != std::string::npos,
              "diagnostic should report the tolerance, got: " + discarded[0]);
    }
}

// The classification boundary is |d1 - d2| <= tolerance: exactly at the
// tolerance is accepted, just beyond it is a conflict.
void test_tolerance_boundary() {
    constexpr double tolerance = 1.0e-6;

    // Exactly at the boundary: accepted as tolerance-equivalent.
    {
        const double other = 1.5 + tolerance;
        const std::string error = capture_error([&] {
            make_solver_with_tolerance({{0, 1, 1.5}, {0, 1, other}}, tolerance);
        });
        check(error.empty(),
              "a difference exactly equal to the tolerance must be accepted, got: " + error);

        auto solver = make_solver_with_tolerance({{0, 1, 1.5}, {0, 1, other}}, tolerance);
        check(solver.normalization_diagnostics().tolerance_equivalent_duplicates == 1,
              "a difference exactly at the tolerance is tolerance-equivalent");
        check(solver.constraints()[0].target_distance == 1.5,
              "first value must still win at the boundary");
    }

    // Just outside: rejected as a conflict.
    {
        const double other = 1.5 + tolerance * 4.0;
        const std::string error = capture_error([&] {
            make_solver_with_tolerance({{0, 1, 1.5}, {0, 1, other}}, tolerance);
        });
        check(!error.empty(),
              "a difference beyond the tolerance must be rejected as a conflict");
        check(error.find("exceeds the constraint tolerance") != std::string::npos,
              "conflict error should name the tolerance, got: " + error);
    }
}

// Whichever orientation each pair arrives in, the kept value and the counts
// must be the same.
void test_orientation_independence() {
    constexpr double tolerance = 1.0e-6;
    const double near = 1.5 + 1.0e-9;

    auto forward = make_solver_with_tolerance({{0, 1, 1.5}, {0, 1, near}}, tolerance);
    auto reversed = make_solver_with_tolerance({{1, 0, 1.5}, {0, 1, near}}, tolerance);
    auto both_reversed = make_solver_with_tolerance({{1, 0, 1.5}, {1, 0, near}}, tolerance);

    for (const auto* solver : {&forward, &reversed, &both_reversed}) {
        check(solver->active_constraint_count() == 1,
              "orientation must not change the constraint count");
        check(solver->constraints()[0].i == 0 && solver->constraints()[0].j == 1,
              "orientation must not change the stored (min, max) pair");
        check(solver->constraints()[0].target_distance == 1.5,
              "orientation must not change which target distance is kept");
        check(solver->normalization_diagnostics().tolerance_equivalent_duplicates == 1,
              "orientation must not change the duplicate classification");
    }
}

void test_conflicting_duplicate_is_rejected() {
    const std::string error = capture_error([] {
        make_solver({{0, 1, 1.5}, {0, 1, 2.5}});
    });
    check(!error.empty(),
          "the same pair with two different target distances must be rejected");
    check(error.find("Conflicting") != std::string::npos,
          "conflict error should say so, got: " + error);

    // Reversed order is the same pair, so it must conflict too.
    check(!capture_error([] { make_solver({{0, 1, 1.5}, {1, 0, 2.5}}); }).empty(),
          "a reversed conflicting duplicate must also be rejected");

    // A difference below the solver tolerance is not a conflict: SHAKE could
    // not tell the two targets apart anyway.
    gmd::ConstraintSettings settings;
    settings.tolerance = 1.0e-6;
    const std::string tiny = capture_error([&] {
        gmd::ConstraintSolver(
            {gmd::BondConstraint{0, 1, 1.5}, gmd::BondConstraint{0, 1, 1.5 + 1.0e-9}},
            settings);
    });
    check(tiny.empty(),
          "targets closer than the constraint tolerance must not be a conflict: " + tiny);
}

void test_self_constraint_is_rejected() {
    const std::string error = capture_error([] { make_solver({{3, 3, 1.0}}); });
    check(!error.empty(), "a self-constraint must be rejected");
    check(error.find("itself") != std::string::npos,
          "self-constraint error should say so, got: " + error);
}

void test_invalid_tags_and_distances_are_rejected() {
    check(!capture_error([] { make_solver({{-1, 2, 1.0}}); }).empty(),
          "a negative atom tag must be rejected");
    check(!capture_error([] { make_solver({{0, 1, 0.0}}); }).empty(),
          "a zero target distance must be rejected");
    check(!capture_error([] { make_solver({{0, 1, -1.0}}); }).empty(),
          "a negative target distance must be rejected");
    check(!capture_error([] {
              make_solver({{0, 1, std::numeric_limits<double>::quiet_NaN()}});
          }).empty(),
          "a NaN target distance must be rejected");
    check(!capture_error([] {
              make_solver({{0, 1, std::numeric_limits<double>::infinity()}});
          }).empty(),
          "an infinite target distance must be rejected");
    check(!capture_error([] {
              make_solver({{0, 1, -std::numeric_limits<double>::infinity()}});
          }).empty(),
          "a negative-infinite target distance must be rejected");
}

// --- Valid topologies ------------------------------------------------------

void test_independent_chain() {
    // 5 atoms, 4 consecutive bonds: independent, nothing to drop.
    auto solver = make_solver({{0, 1, 1.0}, {1, 2, 1.0}, {2, 3, 1.0}, {3, 4, 1.0}});
    check(solver.active_constraint_count() == 4,
          "an independent chain of 4 constraints should count 4, got " +
              std::to_string(solver.active_constraint_count()));
    check(solver.dropped_duplicate_count() == 0,
          "an independent chain should drop nothing");

    // And the DOF calculation follows it: 3*5 - 3 - 4 = 8.
    gmd::System system = make_system(5);
    gmd::VelocityVerletIntegrator integrator(1.0);
    integrator.set_remove_center_of_mass_velocity(true);
    integrator.set_constraint_solver(
        std::make_shared<gmd::ConstraintSolver>(std::move(solver)));
    check(integrator.constraint_count() == 4, "integrator should see 4 constraints");
    check(integrator.degrees_of_freedom(system) == 8,
          "5 atoms with 4 constraints and COM removal should leave 8 DOF, got " +
              std::to_string(integrator.degrees_of_freedom(system)));
}

// A rigid triangle: 3 atoms, 3 distance constraints. These *are* independent
// (9 coordinates - 3 constraints = 6, the rigid-body DOF of a planar triangle),
// so this closed topology is both accepted and correctly counted.
void test_closed_but_independent_triangle() {
    auto solver = make_solver({{0, 1, 1.0}, {1, 2, 1.0}, {0, 2, 1.0}});
    check(solver.active_constraint_count() == 3,
          "a rigid triangle has 3 independent constraints, got " +
              std::to_string(solver.active_constraint_count()));

    gmd::System system = make_system(3);
    gmd::VelocityVerletIntegrator integrator(1.0);
    integrator.set_remove_center_of_mass_velocity(false);
    integrator.set_constraint_solver(
        std::make_shared<gmd::ConstraintSolver>(std::move(solver)));
    check(integrator.degrees_of_freedom(system) == 6,
          "a free rigid triangle should retain 6 DOF, got " +
              std::to_string(integrator.degrees_of_freedom(system)));
}

// A redundant closed topology: all 10 pairs among 5 atoms. A rigid body has
// only 15 - 6 = 9 removable degrees of freedom, so one of the 10 constraints is
// dependent on the others.
//
// This documents the SUPPORTED BEHAVIOUR, which is not the physically ideal
// one: the set contains no duplicates, so all 10 are accepted and counted, and
// the DOF comes out one too low. Detecting this needs the rank of the
// constraint Jacobian, which is configuration dependent and is not attempted.
// The requirement that configured constraints be independent is documented on
// ConstraintSolver.
void test_redundant_closed_topology_is_accepted_and_overcounts() {
    std::vector<gmd::BondConstraint> all_pairs;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            all_pairs.push_back(gmd::BondConstraint{i, j, 1.0});
        }
    }
    check(all_pairs.size() == 10, "the test should build 10 distinct pairs");

    auto solver = make_solver(all_pairs);
    check(solver.active_constraint_count() == 10,
          "all 10 distinct pairs are accepted (redundancy is not detected), got " +
              std::to_string(solver.active_constraint_count()));
    check(solver.dropped_duplicate_count() == 0,
          "no duplicates exist in an all-pairs set");

    gmd::System system = make_system(5);
    gmd::VelocityVerletIntegrator integrator(1.0);
    integrator.set_remove_center_of_mass_velocity(false);
    integrator.set_constraint_solver(
        std::make_shared<gmd::ConstraintSolver>(std::move(solver)));

    // Documented behaviour: 15 - 10 = 5, whereas a rigid body actually retains
    // 6. Pinning this down means the day independence detection is added, this
    // test fails loudly and gets updated deliberately.
    check(integrator.degrees_of_freedom(system) == 5,
          "documented (over-subtracting) behaviour for a redundant set is 5 DOF, got " +
              std::to_string(integrator.degrees_of_freedom(system)));
}

// --- Deduplication must flow through to the DOF ----------------------------

void test_duplicates_do_not_inflate_dof() {
    gmd::System system = make_system(4);

    // Three consecutive bonds, but each supplied twice and one reversed.
    auto solver = std::make_shared<gmd::ConstraintSolver>(
        std::vector<gmd::BondConstraint>{
            {0, 1, 1.0}, {1, 0, 1.0},
            {1, 2, 1.0}, {1, 2, 1.0},
            {2, 3, 1.0}, {3, 2, 1.0},
        },
        gmd::ConstraintSettings{});

    check(solver->active_constraint_count() == 3,
          "six entries describing three bonds must count 3, got " +
              std::to_string(solver->active_constraint_count()));

    gmd::VelocityVerletIntegrator integrator(1.0);
    integrator.set_remove_center_of_mass_velocity(true);
    integrator.set_constraint_solver(solver);
    // 3*4 - 3 - 3 = 6, not 3*4 - 3 - 6 = 3.
    check(integrator.degrees_of_freedom(system) == 6,
          "duplicate constraint entries must not inflate the DOF subtraction, got " +
              std::to_string(integrator.degrees_of_freedom(system)));
}

// Constraints built from a topology go through the same normalisation, so a
// bond also listed as an explicit constraint is not counted twice.
void test_topology_derived_duplicates_collapse() {
    gmd::Topology topology;
    topology.bonds = {{0, 1, 0}, {1, 2, 0}};
    topology.constraints = {gmd::BondConstraint{1, 0, 1.0}};  // same as bond 0-1

    auto constraints = gmd::constraints_from_bond_types(topology, {0}, {1.0});
    check(constraints.size() == 3,
          "the raw list should hold the explicit constraint plus both bonds");

    auto solver = make_solver(constraints);
    check(solver.active_constraint_count() == 2,
          "0-1 listed both explicitly and as a bond must collapse to one, leaving 2, got " +
              std::to_string(solver.active_constraint_count()));
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif

    test_reversed_duplicate_is_one_constraint();
    test_bit_identical_duplicates();
    test_tolerance_equivalent_duplicate();
    test_tolerance_boundary();
    test_orientation_independence();
    test_conflicting_duplicate_is_rejected();
    test_self_constraint_is_rejected();
    test_invalid_tags_and_distances_are_rejected();
    test_independent_chain();
    test_closed_but_independent_triangle();
    test_redundant_closed_topology_is_accepted_and_overcounts();
    test_duplicates_do_not_inflate_dof();
    test_topology_derived_duplicates_collapse();

    if (failures != 0) {
        std::cerr << "[constraints] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[constraints] all checks passed\n";
    return 0;
}
