// Constraint independence: does a constraint set remove as many degrees of
// freedom as it has constraints?
//
// THE CRITERION. A holonomic constraint sigma_c = |r_c|^2 - d_c^2 removes one
// degree of freedom only if its gradient is linearly independent of the others,
// compared in the metric the dynamics uses. The constrained equations of motion
// involve J M^-1 J^T, so the relevant matrix is the mass-weighted Jacobian
// J_M = J M^(-1/2), and the number of degrees of freedom removed is rank(J_M).
//
// Distinct is not independent. Every pair among four atoms is six distinct
// constraints over a body with 3*4 - 6 = 6 internal degrees of freedom, so it is
// exactly rigid; add a fifth atom and the all-pairs set has ten constraints
// against nine. Three collinear atoms have three distinct pair distances of
// which only two are free. None of that is visible to a duplicate check.
//
// These tests pin both directions: independent sets must be accepted with the
// rank equal to the count, and each way of being dependent must be rejected with
// a diagnostic that names the atoms.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

using Vec3 = gmd::System::Vec3;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[constraint independence] " << message << '\n';
        ++failures;
    }
}

gmd::ConstraintSettings settings() {
    gmd::ConstraintSettings s;
    s.tolerance = 1.0e-10;
    s.max_iterations = 500;
    return s;
}

// A system with the given atom positions, unit-ish but distinct masses so that
// the mass weighting is exercised rather than cancelling out.
gmd::System make_system(const std::vector<Vec3>& positions) {
    gmd::System system;
    system.resize(positions.size(), positions.size());
    gmd::Box box;
    box.set_lengths({40.0, 44.0, 48.0});
    system.set_box(box);
    auto masses = system.mutable_masses();
    auto tags = system.mutable_atom_tags();
    auto coordinates = system.mutable_coordinates();
    for (std::size_t i = 0; i < positions.size(); ++i) {
        masses[i] = 1.0 + 0.5 * static_cast<double>(i);   // 1.0, 1.5, 2.0, ...
        tags[i] = static_cast<int>(i);
        coordinates[i] = positions[i];
    }
    return system;
}

double distance(const gmd::System& system, int a, int b) {
    const auto c = system.coordinates();
    double sum = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        const double delta = c[static_cast<std::size_t>(a)][d] - c[static_cast<std::size_t>(b)][d];
        sum += delta * delta;
    }
    return std::sqrt(sum);
}

// Constraints holding the pairs at whatever distance they currently sit,
// so the fixture starts on the constraint manifold.
std::vector<gmd::BondConstraint> hold(const gmd::System& system,
                                      const std::vector<std::array<int, 2>>& pairs) {
    std::vector<gmd::BondConstraint> constraints;
    for (const auto& pair : pairs) {
        constraints.push_back({pair[0], pair[1], distance(system, pair[0], pair[1])});
    }
    return constraints;
}

// Runs the analysis and returns the report; `expect_accepted` decides whether
// require_independent() is expected to pass or throw.
gmd::ConstraintRankReport analyze(const gmd::System& system,
                                  const std::vector<gmd::BondConstraint>& constraints,
                                  const std::string& label,
                                  bool expect_accepted,
                                  const std::vector<std::string>& expected_fragments = {}) {
    gmd::ConstraintSolver solver(constraints, settings());
    const auto report = solver.analyze_independence(system);

    std::string thrown;
    try {
        solver.require_independent(system);
    } catch (const std::exception& error) {
        thrown = error.what();
    }

    if (expect_accepted) {
        check(report.independent, label + ": must be accepted as independent");
        check(thrown.empty(), label + ": require_independent must not throw, got: " + thrown);
        check(report.rank == report.constraint_count,
              label + ": rank " + std::to_string(report.rank) + " must equal the constraint "
              "count " + std::to_string(report.constraint_count));
    } else {
        check(!report.independent, label + ": must be rejected as dependent");
        check(!thrown.empty(), label + ": require_independent must throw");
        check(report.rank < report.constraint_count,
              label + ": rank " + std::to_string(report.rank) + " must be below the "
              "constraint count " + std::to_string(report.constraint_count));
        check(!report.problems.empty(), label + ": must report at least one problem");
        for (const auto& fragment : expected_fragments) {
            const bool in_message = thrown.find(fragment) != std::string::npos;
            check(in_message, label + ": diagnostic must mention \"" + fragment +
                              "\"; got:\n" + thrown);
        }
        std::cout << "  rejected " << label << ":\n";
        for (const auto& problem : report.problems) {
            std::cout << "      " << problem << '\n';
        }
    }
    return report;
}

// --- Independent sets ------------------------------------------------------

void test_independent_chain() {
    gmd::System system = make_system({{10.0, 10.0, 10.0},
                                      {11.5, 10.2, 10.3},
                                      {12.6, 11.3, 10.9},
                                      {13.9, 11.4, 11.8}});
    const auto constraints = hold(system, {{0, 1}, {1, 2}, {2, 3}});
    const auto report = analyze(system, constraints, "independent 4-atom chain", true);
    check(report.components.size() == 1, "a bonded chain is one connected component");
    check(report.components.front().atom_tags.size() == 4,
          "the component must contain all four atoms");
    check(std::isfinite(report.components.front().condition_number),
          "a well-separated chain must be well conditioned");
}

void test_independent_triangle() {
    gmd::System system = make_system({{20.0, 20.0, 20.0},
                                      {21.4, 20.0, 20.0},
                                      {20.5, 21.1, 20.0}});
    const auto constraints = hold(system, {{0, 1}, {1, 2}, {0, 2}});
    const auto report = analyze(system, constraints, "independent rigid triangle", true);
    // 3 atoms, 9 coordinates, 6 rigid-body motions -> exactly 3 internal DOF.
    check(report.rank == 3, "a non-degenerate triangle has three independent constraints");
}

// Two molecules that share no atom. Their Jacobian rows have disjoint support,
// so the rank is additive and the analysis must split them.
void test_disjoint_components() {
    gmd::System system = make_system({{10.0, 10.0, 10.0}, {11.5, 10.0, 10.0},
                                      {30.0, 30.0, 30.0}, {31.5, 30.2, 30.1}});
    const auto constraints = hold(system, {{0, 1}, {2, 3}});
    const auto report = analyze(system, constraints, "two disjoint dimers", true);
    check(report.components.size() == 2,
          "two molecules sharing no atom must be two components, got " +
              std::to_string(report.components.size()));
    check(report.rank == 2, "and contribute one independent constraint each");
}

// --- Duplicates and conflicts, caught during normalisation -----------------

void test_exact_and_reversed_duplicates() {
    gmd::System system = make_system({{10.0, 10.0, 10.0}, {11.5, 10.0, 10.0}});
    const double d = distance(system, 0, 1);

    gmd::ConstraintSolver exact({{0, 1, d}, {0, 1, d}}, settings());
    check(exact.active_constraint_count() == 1,
          "an exact duplicate must collapse to one constraint");
    check(exact.normalization_diagnostics().exact_duplicates == 1,
          "and be reported as an exact duplicate");
    check(exact.normalization_diagnostics().reversed_duplicates == 0,
          "a same-order duplicate is not a reversed one");
    check(exact.analyze_independence(system).independent,
          "after collapsing, the remaining single constraint is independent");

    gmd::ConstraintSolver reversed({{0, 1, d}, {1, 0, d}}, settings());
    check(reversed.active_constraint_count() == 1,
          "a reversed duplicate is the same constraint and must collapse too");
    check(reversed.normalization_diagnostics().reversed_duplicates == 1,
          "and must be reported as reversed, since a topology listing a bond twice "
          "in opposite orders is usually a generation bug");
    check(reversed.normalization_diagnostics().exact_duplicates == 1,
          "a reversed repeat with an identical target is also an exact duplicate");
}

void test_tolerance_equivalent_duplicate() {
    gmd::ConstraintSettings loose = settings();
    loose.tolerance = 1.0e-6;
    gmd::ConstraintSolver solver({{0, 1, 1.5}, {0, 1, 1.5 + 1.0e-9}}, loose);
    check(solver.active_constraint_count() == 1,
          "targets closer than the solver tolerance are one constraint");
    check(solver.normalization_diagnostics().tolerance_equivalent_duplicates == 1,
          "and must be counted as a tolerance-equivalent duplicate");
    check(solver.normalization_diagnostics().discarded_targets.size() == 1,
          "with the discarded target recorded");
}

void test_conflicting_targets() {
    bool threw = false;
    std::string message;
    try {
        gmd::ConstraintSolver solver({{0, 1, 1.5}, {0, 1, 1.9}}, settings());
    } catch (const std::exception& error) {
        threw = true;
        message = error.what();
    }
    check(threw, "two different targets for the same pair must be rejected");
    check(message.find("Conflicting") != std::string::npos,
          "the diagnostic must say the constraints conflict, got: " + message);
    check(message.find("(0, 1)") != std::string::npos,
          "and must name the atom pair, got: " + message);
    std::cout << "  rejected conflicting targets:\n      " << message << '\n';
}

// --- Dependent sets, caught by the rank -----------------------------------

// Every pair among five atoms: ten constraints over a body with at most
// 3*5 - 6 = 9 internal degrees of freedom.
void test_redundant_closed_network() {
    gmd::System system = make_system({{20.0, 20.0, 20.0},
                                      {21.5, 20.1, 20.2},
                                      {20.6, 21.4, 20.3},
                                      {20.4, 20.2, 21.6},
                                      {21.7, 21.3, 21.5}});
    std::vector<std::array<int, 2>> pairs;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) pairs.push_back({i, j});
    }
    const auto report = analyze(system, hold(system, pairs),
                                "all-pairs cage over five atoms", false,
                                {"over-constrained", "rank 9"});
    check(report.constraint_count == 10, "the fixture must supply ten constraints");
    check(report.rank == 9, "and the rigid body can only remove nine degrees of freedom");
    check(report.redundant_count() == 1, "leaving exactly one redundant constraint");
}

// Three collinear atoms holding their CURRENT distances. Those distances satisfy
// the triangle inequality only to round-off, so the implied triangle height is
// zero to the solver's resolution and the TARGET check rejects the set before
// any geometry is involved. That is the right place for it: no configuration
// these targets admit is non-degenerate, so no amount of projecting would help.
void test_over_constrained_collinear() {
    gmd::System system = make_system({{10.0, 10.0, 10.0},
                                      {11.5, 10.0, 10.0},
                                      {13.7, 10.0, 10.0}});
    const auto constraints = hold(system, {{0, 1}, {1, 2}, {0, 2}});
    bool threw = false;
    std::string message;
    try {
        gmd::ConstraintSolver solver(constraints, settings());
    } catch (const std::exception& error) {
        threw = true;
        message = error.what();
    }
    check(threw, "three collinear pair distances must be rejected");
    check(message.find("Degenerate constraint targets") != std::string::npos,
          "as degenerate targets, got: " + message);
    check(message.find("COLLINEAR") != std::string::npos,
          "and the diagnostic must say why, got: " + message);
    check(message.find("triangle height") != std::string::npos,
          "and must report the implied triangle height it judged, got: " + message);
    std::cout << "  rejected collinear triple targets:\n      " << message << "\n";
}

// Four collinear atoms with the chain plus its closing bond. There is NO
// constrained triple here -- (0,2) and (1,3) are not constrained -- so the
// triangle-inequality check on the targets cannot see this at all, and only the
// rank of the mass-weighted Jacobian does. On a straight line every gradient
// points along the same direction, so the four of them span three dimensions.
void test_geometric_dependence_without_a_constrained_triple() {
    gmd::System system = make_system({{10.0, 10.0, 10.0},
                                      {11.3, 10.0, 10.0},
                                      {12.9, 10.0, 10.0},
                                      {14.6, 10.0, 10.0}});
    const auto report = analyze(system, hold(system, {{0, 1}, {1, 2}, {2, 3}, {0, 3}}),
                                "four collinear atoms, chain plus closing bond", false,
                                {"rank 3", "geometric"});
    check(report.rank == 3,
          "four collinear pair distances have only three independent gradients");
    check(report.redundant_identified(),
          "and the redundant row must be identified");
}

// A bond of zero length. grad sigma = 2 r_ij vanishes, so the constraint cannot
// remove anything and SHAKE has nothing to project along.
void test_degenerate_zero_length_gradient() {
    gmd::System system = make_system({{10.0, 10.0, 10.0},
                                      {10.0, 10.0, 10.0},
                                      {12.0, 10.0, 10.0}});
    const std::vector<gmd::BondConstraint> constraints{{0, 1, 1.5}, {1, 2, 2.0}};
    analyze(system, constraints, "coincident atoms, zero-length gradient", false,
            {"degenerate gradient", "(0, 1)"});
}

// Independent, but close to the configuration where it would not be. Four atoms
// with all six pair distances is exactly rigid -- 3*4 - 6 = 6 -- but only while
// they are not coplanar: flatten them and reflection through the plane becomes a
// zero mode, dropping the rank to five. Here the fourth atom stands 1e-8 A off
// the plane of the other three.
//
// This must be ACCEPTED -- rejecting a valid geometry because it is awkward would
// be worse than the problem -- and flagged. Note that every one of the six
// triangles is perfectly well formed, so the target check cannot see this: the
// degeneracy is out of plane, and only the rank of the Jacobian exposes it.
void test_near_singular_but_valid() {
    const double height = 1.0e-8;
    gmd::System system = make_system({{20.0, 20.0, 20.0},
                                      {21.4, 20.0, 20.0},
                                      {20.5, 21.1, 20.0},
                                      {20.633333333333, 20.366666666667, 20.0 + height}});
    std::vector<std::array<int, 2>> pairs;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) pairs.push_back({i, j});
    }
    const auto report = analyze(system, hold(system, pairs),
                                "near-coplanar four-atom cage", true);
    check(report.rank == 6,
          "four non-coplanar atoms with all six distances are exactly rigid, rank 6; "
          "got " + std::to_string(report.rank));
    check(!report.warnings.empty(),
          "and being 1e-8 A from coplanar must be flagged as ill-conditioned rather "
          "than silently accepted");
    if (!report.warnings.empty()) {
        std::cout << "  warned (accepted) near-coplanar cage:\n      "
                  << report.warnings.front() << '\n';
    }
    check(report.components.front().condition_number > 1.0e7,
          "the condition number must actually be large, got " +
              std::to_string(report.components.front().condition_number));
}

// The rank tolerance is RELATIVE to the largest singular value, so the verdict
// cannot depend on the unit the coordinates happen to be expressed in. Scaling
// every length by 1000 scales every singular value by 1000 and must change
// nothing.
void test_rank_tolerance_is_scale_aware() {
    const std::vector<Vec3> base{{20.0, 20.0, 20.0}, {21.4, 20.0, 20.0}, {20.5, 21.1, 20.0}};
    gmd::System small = make_system(base);
    const auto small_report =
        gmd::ConstraintSolver(hold(small, {{0, 1}, {1, 2}, {0, 2}}), settings())
            .analyze_independence(small);

    std::vector<Vec3> scaled;
    for (const auto& p : base) scaled.push_back({p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0});
    gmd::System large;
    large.resize(3, 3);
    gmd::Box box;
    box.set_lengths({100000.0, 100000.0, 100000.0});
    large.set_box(box);
    for (std::size_t i = 0; i < 3; ++i) {
        large.mutable_masses()[i] = 1.0 + 0.5 * static_cast<double>(i);
        large.mutable_atom_tags()[i] = static_cast<int>(i);
        large.mutable_coordinates()[i] = scaled[i];
    }
    const auto large_report =
        gmd::ConstraintSolver(hold(large, {{0, 1}, {1, 2}, {0, 2}}), settings())
            .analyze_independence(large);

    check(small_report.rank == large_report.rank,
          "scaling every coordinate by 1000 must not change the rank verdict");
    check(large_report.components.front().rank_tolerance >
              small_report.components.front().rank_tolerance * 100.0,
          "the tolerance must scale with the matrix, not be a fixed absolute number");
    const double small_condition = small_report.components.front().condition_number;
    const double large_condition = large_report.components.front().condition_number;
    check(std::abs(small_condition - large_condition) < 1.0e-6 * small_condition,
          "and the condition number, being a ratio, must be unchanged: " +
              std::to_string(small_condition) + " vs " + std::to_string(large_condition));
}

// --- Degrees of freedom ----------------------------------------------------

void test_degrees_of_freedom_uses_the_verified_rank() {
    gmd::System system = make_system({{10.0, 10.0, 10.0},
                                      {11.5, 10.2, 10.3},
                                      {12.6, 11.3, 10.9},
                                      {13.9, 11.4, 11.8}});
    auto solver = std::make_shared<gmd::ConstraintSolver>(
        hold(system, {{0, 1}, {1, 2}, {2, 3}}), settings());
    gmd::VelocityVerletIntegrator integrator(0.5);
    integrator.set_constraint_solver(solver);
    integrator.set_remove_center_of_mass_velocity(true);
    gmd::RuntimeContext runtime;
    integrator.initialize(system, runtime);

    // 3N - 3 (COM) - 3 (constraints) = 12 - 3 - 3 = 6.
    check(integrator.degrees_of_freedom(system) == 6,
          "DOF must be 3N - 3 - rank, got " +
              std::to_string(integrator.degrees_of_freedom(system)));
    check(integrator.constraint_count() == 3, "and the constraint count must be the rank");
    check(integrator.constraint_rank_report().independent,
          "initialize() must have recorded the independence check");
    check(integrator.constraint_rank_report().rank == 3, "with rank 3");
}

// A dependent set must stop the run before any dynamics, not quietly
// over-subtract degrees of freedom. Four collinear atoms with the chain plus its
// closing bond: no constrained triple, so only the rank of the projected geometry
// can catch it, and it has to catch it inside initialize().
void test_initialize_rejects_a_dependent_set() {
    gmd::System system = make_system({{10.0, 10.0, 10.0},
                                      {11.3, 10.0, 10.0},
                                      {12.9, 10.0, 10.0},
                                      {14.6, 10.0, 10.0}});
    auto solver = std::make_shared<gmd::ConstraintSolver>(
        hold(system, {{0, 1}, {1, 2}, {2, 3}, {0, 3}}), settings());
    gmd::VelocityVerletIntegrator integrator(0.5);
    integrator.set_constraint_solver(solver);
    gmd::RuntimeContext runtime;

    bool threw = false;
    std::string message;
    try {
        integrator.initialize(system, runtime);
    } catch (const std::exception& error) {
        threw = true;
        message = error.what();
    }
    check(threw, "initialize() must reject a dependent constraint set");
    check(message.find("not independent") != std::string::npos,
          "with a diagnostic naming the problem, got: " + message);
    check(message.find("(0,3)") != std::string::npos ||
              message.find("(0, 3)") != std::string::npos,
          "and naming the constrained atom pairs, got: " + message);
    check(integrator.constraint_rank_report().constraint_count == 0,
          "and must leave no accepted report behind for degrees_of_freedom() to use");
}


// --- The geometry that decides is the PROJECTED one ------------------------

// Full rank as supplied, degenerate as targeted. The three atoms start as a
// proper triangle, but their TARGET distances satisfy d02 = d01 + d12, so the
// only configuration satisfying all three is collinear -- where the three
// gradients span two dimensions. A rank check on the supplied geometry passes
// this; the run must still be rejected.
void test_full_rank_input_projecting_onto_a_degenerate_target() {
    gmd::System system = make_system({{10.0, 10.0, 10.0},
                                      {11.0, 10.3, 10.0},
                                      {12.0, 10.0, 10.0}});
    const std::vector<gmd::BondConstraint> as_supplied{
        {0, 1, distance(system, 0, 1)},
        {1, 2, distance(system, 1, 2)},
        {0, 2, distance(system, 0, 2)},
    };
    const auto supplied_report =
        gmd::ConstraintSolver(as_supplied, settings()).analyze_independence(system);
    check(supplied_report.independent && supplied_report.rank == 3,
          "the fixture must be full rank AS SUPPLIED, or it does not test anything: "
          "rank " + std::to_string(supplied_report.rank));

    const std::vector<gmd::BondConstraint> degenerate_targets{
        {0, 1, 1.0}, {1, 2, 1.0}, {0, 2, 2.0}};
    bool threw = false;
    std::string message;
    try {
        auto solver = std::make_shared<gmd::ConstraintSolver>(degenerate_targets, settings());
        gmd::VelocityVerletIntegrator integrator(0.5);
        integrator.set_constraint_solver(solver);
        gmd::RuntimeContext runtime;
        integrator.initialize(system, runtime);
    } catch (const std::exception& error) {
        threw = true;
        message = error.what();
    }
    check(threw,
          "a degenerate target must be rejected even though the supplied coordinates "
          "are full rank");
    check(message.find("Degenerate constraint targets") != std::string::npos ||
              message.find("not independent") != std::string::npos ||
              message.find("DEGENERATE TARGET") != std::string::npos,
          "and the diagnostic must identify the targets, not blame the solver; got:\n" +
              message);
    check(message.find("COLLINEAR") != std::string::npos ||
              message.find("collinear") != std::string::npos,
          "and should say why the targets are degenerate; got:\n" + message);
    std::cout << "  rejected degenerate-target triangle:\n      " << message << "\n";
}

// The authoritative report describes the PROJECTED geometry, not the supplied
// one. Here the two genuinely differ: the atoms are supplied nearly collinear
// (ill-conditioned) while their targets describe a healthy triangle, so the
// projection improves the conditioning by orders of magnitude. If the analysis
// still ran on the input, the accepted report would carry a warning it should
// not.
void test_authoritative_report_describes_the_projected_geometry() {
    gmd::System system = make_system({{20.0, 20.0, 20.0},
                                      {21.4, 20.0, 20.0},
                                      {22.8, 20.0 + 2.0e-8, 20.0}});
    const std::vector<gmd::BondConstraint> constraints{
        {0, 1, 1.4}, {1, 2, 1.4}, {0, 2, 1.6}};   // a healthy, well-spread triangle

    const auto supplied =
        gmd::ConstraintSolver(constraints, settings()).analyze_independence(system);
    check(supplied.independent, "the supplied geometry is full rank, just awkward");
    check(!supplied.warnings.empty(),
          "and must be flagged ill-conditioned as supplied, or the two geometries do "
          "not differ enough for this test to mean anything");
    const double supplied_condition =
        supplied.components.empty() ? 0.0 : supplied.components.front().condition_number;

    auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, settings());
    gmd::VelocityVerletIntegrator integrator(0.5);
    integrator.set_constraint_solver(solver);
    gmd::RuntimeContext runtime;
    integrator.initialize(system, runtime);

    const auto& authoritative = integrator.constraint_rank_report();
    check(authoritative.independent, "the projected geometry must be accepted");
    check(!authoritative.components.empty(), "with a component recorded");
    const double projected_condition =
        authoritative.components.empty() ? 0.0
                                         : authoritative.components.front().condition_number;
    check(authoritative.warnings.empty(),
          "and must carry no ill-conditioning warning, because the geometry it "
          "describes is the projected one, not the supplied one");
    check(projected_condition < supplied_condition / 1.0e3,
          "the projection must have improved the conditioning by orders of magnitude, "
          "which is what makes the two reports distinguishable: supplied " +
              std::to_string(supplied_condition) + " vs projected " +
              std::to_string(projected_condition));
    std::cout << "  authoritative report is post-projection: condition number "
              << supplied_condition << " as supplied -> " << projected_condition
              << " once projected\n";
}

// Moderately off the manifold but perfectly valid: the targets are 8% away from
// the supplied distances, so the projection has real work to do, and the result
// is an ordinary non-degenerate triangle.
void test_off_manifold_but_valid_target_projects_and_stays_independent() {
    gmd::System system = make_system({{20.0, 20.0, 20.0},
                                      {21.4, 20.0, 20.0},
                                      {20.5, 21.1, 20.0}});
    std::vector<gmd::BondConstraint> constraints;
    for (const auto& pair : std::vector<std::array<int, 2>>{{0, 1}, {1, 2}, {0, 2}}) {
        constraints.push_back({pair[0], pair[1], 1.08 * distance(system, pair[0], pair[1])});
    }
    const double before = distance(system, 0, 1);

    auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, settings());
    gmd::VelocityVerletIntegrator integrator(0.5);
    integrator.set_constraint_solver(solver);
    gmd::RuntimeContext runtime;
    integrator.initialize(system, runtime);

    check(integrator.constraint_rank_report().independent,
          "a valid target that merely needs projecting must be accepted");
    check(integrator.constraint_rank_report().rank == 3, "with rank 3");
    for (const auto& constraint : constraints) {
        const double achieved = distance(system, constraint.i, constraint.j);
        check(std::abs(achieved - constraint.target_distance) < 1.0e-9,
              "and the projection must actually have reached the target: pair (" +
                  std::to_string(constraint.i) + ", " + std::to_string(constraint.j) +
                  ") is at " + std::to_string(achieved) + " against a target of " +
                  std::to_string(constraint.target_distance));
    }
    check(std::abs(distance(system, 0, 1) - before) > 0.05,
          "the fixture must have started meaningfully off the manifold, or the "
          "projection was a no-op");
    const auto coordinates = system.coordinates();
    const auto velocities = system.velocities();
    for (const auto& constraint : constraints) {
        double dot = 0.0;
        for (std::size_t d = 0; d < 3; ++d) {
            dot += (coordinates[static_cast<std::size_t>(constraint.i)][d] -
                    coordinates[static_cast<std::size_t>(constraint.j)][d]) *
                   (velocities[static_cast<std::size_t>(constraint.i)][d] -
                    velocities[static_cast<std::size_t>(constraint.j)][d]);
        }
        check(std::abs(dot) < 1.0e-10,
              "velocities must be projected onto the tangent space once the set is "
              "accepted, r_ij . v_ij = " + std::to_string(dot));
    }
}

// --- The SVD kernel, against independently known singular values -----------

// A matrix with EXACTLY the requested singular values, built as A = U S V^T from
// orthonormal factors and returned as columns. Nothing here uses the routine
// under test, so the singular values are known independently rather than
// measured.
std::vector<std::vector<double>> matrix_with_singular_values(
        std::size_t rows, const std::vector<double>& sigma, unsigned seed) {
    const std::size_t n = sigma.size();
    unsigned state = seed * 2654435761u + 1u;
    auto next = [&state]() {
        state = state * 1664525u + 1013904223u;
        return static_cast<double>(state % 20000u) / 10000.0 - 1.0;
    };
    std::vector<std::vector<double>> u;
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<double> v(rows);
        for (double& c : v) c = next();
        for (int pass = 0; pass < 2; ++pass) {          // classical GS, twice
            for (const auto& b : u) {
                double projection = 0.0;
                for (std::size_t k = 0; k < rows; ++k) projection += v[k] * b[k];
                for (std::size_t k = 0; k < rows; ++k) v[k] -= projection * b[k];
            }
        }
        double norm = 0.0;
        for (double c : v) norm += c * c;
        norm = std::sqrt(norm);
        for (double& c : v) c /= norm;
        u.push_back(v);
    }
    std::vector<std::vector<double>> V(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) V[i][i] = 1.0;
    for (std::size_t p = 0; p + 1 < n; ++p) {
        for (std::size_t q = p + 1; q < n; ++q) {
            const double angle = 0.31 * static_cast<double>(p + 1) +
                                 0.17 * static_cast<double>(q + 1) + 0.05 * seed;
            const double c = std::cos(angle);
            const double sn = std::sin(angle);
            for (std::size_t k = 0; k < n; ++k) {
                const double a = V[k][p];
                const double b = V[k][q];
                V[k][p] = c * a - sn * b;
                V[k][q] = sn * a + c * b;
            }
        }
    }
    std::vector<std::vector<double>> columns(n, std::vector<double>(rows, 0.0));
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double weight = sigma[i] * V[j][i];
            for (std::size_t k = 0; k < rows; ++k) columns[j][k] += weight * u[i][k];
        }
    }
    return columns;
}

// A second, independent computation: cyclic Jacobi EIGENvalues of the Gram
// matrix A^T A, carried in long double. Different algorithm (two-sided, on the
// normal equations) at higher precision, so agreement is evidence rather than a
// tautology.
//
// Squaring costs half the digits, so this reference is only meaningful for
// well-conditioned matrices; the exactly- and nearly-deficient cases are checked
// against the analytically constructed singular values instead.
std::vector<double> gram_reference_singular_values(
        const std::vector<std::vector<double>>& columns) {
    const std::size_t n = columns.size();
    const std::size_t rows = columns.front().size();
    std::vector<std::vector<long double>> g(n, std::vector<long double>(n, 0.0L));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            long double sum = 0.0L;
            for (std::size_t k = 0; k < rows; ++k) {
                sum += static_cast<long double>(columns[i][k]) *
                       static_cast<long double>(columns[j][k]);
            }
            g[i][j] = sum;
        }
    }
    for (int sweep = 0; sweep < 100; ++sweep) {
        long double off = 0.0L;
        for (std::size_t p = 0; p + 1 < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) off += g[p][q] * g[p][q];
        }
        if (off <= 1.0e-40L) break;
        for (std::size_t p = 0; p + 1 < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                if (g[p][q] == 0.0L) continue;
                const long double theta = (g[q][q] - g[p][p]) / (2.0L * g[p][q]);
                const long double t = (theta >= 0.0L ? 1.0L : -1.0L) /
                                      (std::fabs(theta) + std::sqrt(1.0L + theta * theta));
                const long double c = 1.0L / std::sqrt(1.0L + t * t);
                const long double sn = c * t;
                for (std::size_t k = 0; k < n; ++k) {
                    const long double gkp = g[k][p];
                    const long double gkq = g[k][q];
                    g[k][p] = c * gkp - sn * gkq;
                    g[k][q] = sn * gkp + c * gkq;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const long double gpk = g[p][k];
                    const long double gqk = g[q][k];
                    g[p][k] = c * gpk - sn * gqk;
                    g[q][k] = sn * gpk + c * gqk;
                }
            }
        }
    }
    std::vector<double> values;
    for (std::size_t i = 0; i < n; ++i) {
        values.push_back(static_cast<double>(std::sqrt(std::fabs(g[i][i]))));
    }
    std::sort(values.begin(), values.end(), std::greater<double>());
    return values;
}

void test_singular_values_against_constructed_matrices() {
    struct Case {
        const char* label;
        std::vector<double> sigma;
        double tolerance;      // relative to sigma_max
        bool gram_reference;   // is the Gram cross-check meaningful here?
    };
    const std::vector<Case> cases{
        {"full rank, well conditioned", {5.0, 3.0, 2.0, 1.0}, 1.0e-13, true},
        {"full rank, spread over six orders", {1.0e3, 1.0e1, 1.0e-1, 1.0e-3}, 1.0e-12, false},
        {"exactly rank deficient", {4.0, 2.0, 1.0, 0.0}, 1.0e-13, true},
        {"two exactly zero", {7.0, 3.0, 0.0, 0.0}, 1.0e-13, true},
        {"near deficient at 1e-12", {1.0, 0.5, 0.25, 1.0e-12}, 1.0e-11, false},
        {"near deficient at 1e-9", {2.0, 1.0, 1.0e-9}, 1.0e-12, false},
        {"larger component, ten columns",
         {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5}, 1.0e-13, true},
    };

    std::cout << "  singular values against independently constructed matrices:\n";
    for (const auto& test_case : cases) {
        for (unsigned seed = 1; seed <= 3; ++seed) {
            const auto columns =
                matrix_with_singular_values(3 * test_case.sigma.size(), test_case.sigma, seed);
            const auto measured = gmd::jacobi_singular_values(columns, test_case.label);
            check(measured.size() == test_case.sigma.size(),
                  std::string(test_case.label) + ": wrong number of singular values");

            std::vector<double> expected = test_case.sigma;
            std::sort(expected.begin(), expected.end(), std::greater<double>());
            double worst = 0.0;
            for (std::size_t i = 0; i < expected.size(); ++i) {
                worst = std::max(worst, std::abs(measured[i] - expected[i]) / expected.front());
            }
            check(worst <= test_case.tolerance,
                  std::string(test_case.label) + " (seed " + std::to_string(seed) +
                      "): singular values differ from the constructed ones by " +
                      std::to_string(worst) + " relative to sigma_max, above " +
                      std::to_string(test_case.tolerance));

            double gram_worst = -1.0;
            if (test_case.gram_reference) {
                const auto reference = gram_reference_singular_values(columns);
                gram_worst = 0.0;
                for (std::size_t i = 0; i < expected.size(); ++i) {
                    gram_worst = std::max(gram_worst,
                                          std::abs(measured[i] - reference[i]) /
                                              expected.front());
                }
                // The Gram route squares the condition number, so it is held to a
                // looser bound than the analytic construction: sqrt(eps) relative.
                check(gram_worst <= 1.0e-8,
                      std::string(test_case.label) +
                          ": disagrees with the long-double Gram-matrix eigenvalue "
                          "reference by " + std::to_string(gram_worst));
            }
            if (seed == 1) {
                std::cout << "      " << test_case.label << ": max relative error "
                          << worst;
                if (gram_worst >= 0.0) std::cout << ", vs Gram reference " << gram_worst;
                std::cout << '\n';
            }
        }
    }
}

// Permuting the columns must not change the singular values.
void test_singular_values_are_permutation_invariant() {
    const std::vector<double> sigma{6.0, 3.0, 1.5, 1.0e-10, 0.0};
    const auto columns = matrix_with_singular_values(20, sigma, 7);
    const auto reference = gmd::jacobi_singular_values(columns, "unpermuted");

    std::vector<std::size_t> order(columns.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    for (int attempt = 0; attempt < 5; ++attempt) {
        std::rotate(order.begin(), order.begin() + 1, order.end());
        if (attempt == 4) std::reverse(order.begin(), order.end());
        std::vector<std::vector<double>> permuted;
        for (std::size_t index : order) permuted.push_back(columns[index]);
        const auto measured = gmd::jacobi_singular_values(permuted, "permuted");
        for (std::size_t i = 0; i < reference.size(); ++i) {
            check(std::abs(measured[i] - reference[i]) <= 1.0e-13 * reference.front(),
                  "permuting the columns changed singular value " + std::to_string(i) +
                      ": " + std::to_string(reference[i]) + " -> " +
                      std::to_string(measured[i]));
        }
    }
}

// A rank cannot be decided from input that is not a number.
void test_non_finite_input_is_rejected() {
    std::vector<std::vector<double>> columns{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    columns[1][2] = std::numeric_limits<double>::quiet_NaN();
    bool threw = false;
    std::string message;
    try {
        gmd::jacobi_singular_values(columns, "component {3, 4}");
    } catch (const std::exception& error) {
        threw = true;
        message = error.what();
    }
    check(threw, "a non-finite Jacobian entry must be rejected, not silently ranked");
    check(message.find("component {3, 4}") != std::string::npos,
          "and the message must name the affected component, got: " + message);
}

// --- Stress: the analysis on larger and badly scaled but legitimate systems -

void test_stress_larger_and_ill_scaled_components() {
    std::vector<Vec3> positions;
    for (int i = 0; i < 10; ++i) {
        positions.push_back({10.0 + 1.3 * i, 10.0 + 0.4 * std::sin(0.7 * i),
                             10.0 + 0.3 * std::cos(1.1 * i)});
    }
    gmd::System system;
    system.resize(positions.size(), positions.size());
    gmd::Box box;
    box.set_lengths({60.0, 60.0, 60.0});
    system.set_box(box);
    for (std::size_t i = 0; i < positions.size(); ++i) {
        // Masses spanning four orders, so the mass-weighted columns differ in
        // scale by a factor of a hundred.
        system.mutable_masses()[i] = (i % 2 == 0) ? 1.0e-2 : 1.0e2;
        system.mutable_atom_tags()[i] = static_cast<int>(i);
        system.mutable_coordinates()[i] = positions[i];
    }
    std::vector<std::array<int, 2>> pairs;
    for (int i = 0; i + 1 < 10; ++i) pairs.push_back({i, i + 1});
    const auto report = analyze(system, hold(system, pairs),
                                "ten-atom chain, masses spanning 10^4", true);
    check(report.rank == 9, "a nine-bond chain has nine independent constraints");

    auto ring_pairs = pairs;
    ring_pairs.push_back({9, 0});
    const auto ring = analyze(system, hold(system, ring_pairs),
                              "closed ten-atom ring, masses spanning 10^4", true);
    check(ring.rank == 10, "a ten-membered ring has ten independent constraints");
}

// The verdict must not depend on the order the constraints were listed in, and
// neither must the set of rows named redundant.
void test_analysis_is_invariant_to_constraint_order() {
    gmd::System system = make_system({{20.0, 20.0, 20.0},
                                      {21.5, 20.1, 20.2},
                                      {20.6, 21.4, 20.3},
                                      {20.4, 20.2, 21.6},
                                      {21.7, 21.3, 21.5}});
    std::vector<std::array<int, 2>> pairs;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) pairs.push_back({i, j});
    }

    const auto base_constraints = hold(system, pairs);
    const auto base =
        gmd::ConstraintSolver(base_constraints, settings()).analyze_independence(system);
    check(base.redundant_identified(),
          "the pivoted factorisation must agree with the singular values about how "
          "many rows are redundant");
    check(base.dependent_constraints().size() == base.redundant_count(),
          "and must name exactly that many, got " +
              std::to_string(base.dependent_constraints().size()) + " for " +
              std::to_string(base.redundant_count()) + " redundant");

    auto pairs_of = [&](const gmd::ConstraintRankReport& report,
                        const std::vector<gmd::BondConstraint>& list) {
        std::vector<std::array<int, 2>> out;
        for (std::size_t index : report.dependent_constraints()) {
            out.push_back({list[index].i, list[index].j});
        }
        std::sort(out.begin(), out.end());
        return out;
    };
    const auto base_pairs = pairs_of(base, base_constraints);

    auto order = pairs;
    for (int attempt = 0; attempt < 4; ++attempt) {
        std::rotate(order.begin(), order.begin() + 3, order.end());
        if (attempt == 3) std::reverse(order.begin(), order.end());
        const auto shuffled_constraints = hold(system, order);
        const auto shuffled = gmd::ConstraintSolver(shuffled_constraints, settings())
                                  .analyze_independence(system);
        check(shuffled.rank == base.rank,
              "reordering the constraint list changed the rank: " +
                  std::to_string(base.rank) + " -> " + std::to_string(shuffled.rank));
        check(shuffled.redundant_identified() == base.redundant_identified(),
              "and changed whether the redundant rows could be identified");
        check(pairs_of(shuffled, shuffled_constraints) == base_pairs,
              "and named a different set of atom pairs as redundant");
    }
    std::cout << "  order invariance: rank " << base.rank << " of "
              << base.constraint_count << ", redundant pair(s)";
    for (const auto& pair : base_pairs) {
        std::cout << " (" << pair[0] << "," << pair[1] << ")";
    }
    std::cout << ", unchanged under 4 reorderings\n";
}

// --- Reversed duplicates, all four orientations ---------------------------

void test_reversed_duplicate_orientations() {
    struct Case {
        const char* label;
        int first_i, first_j, second_i, second_j;
        std::size_t expected_reversed;
    };
    const std::vector<Case> cases{
        {"(0,1) then (1,0)", 0, 1, 1, 0, 1},
        {"(1,0) then (0,1)", 1, 0, 0, 1, 1},
        {"(0,1) then (0,1)", 0, 1, 0, 1, 0},
        {"(1,0) then (1,0)", 1, 0, 1, 0, 0},
    };

    for (const auto& test_case : cases) {
        {   // Identical targets: an exact duplicate.
            gmd::ConstraintSolver solver({{test_case.first_i, test_case.first_j, 1.5},
                                          {test_case.second_i, test_case.second_j, 1.5}},
                                         settings());
            const auto& diagnostics = solver.normalization_diagnostics();
            check(solver.active_constraint_count() == 1,
                  std::string(test_case.label) + ", exact: must collapse to one constraint");
            check(diagnostics.exact_duplicates == 1,
                  std::string(test_case.label) + ", exact: must count one exact duplicate");
            check(diagnostics.tolerance_equivalent_duplicates == 0,
                  std::string(test_case.label) + ", exact: must not count a "
                  "tolerance-equivalent duplicate");
            check(diagnostics.reversed_duplicates == test_case.expected_reversed,
                  std::string(test_case.label) + ", exact: expected " +
                      std::to_string(test_case.expected_reversed) +
                      " reversed duplicate(s), got " +
                      std::to_string(diagnostics.reversed_duplicates));
        }
        {   // Targets differing by less than the solver tolerance.
            gmd::ConstraintSettings loose = settings();
            loose.tolerance = 1.0e-6;
            gmd::ConstraintSolver solver(
                {{test_case.first_i, test_case.first_j, 1.5},
                 {test_case.second_i, test_case.second_j, 1.5 + 1.0e-9}},
                loose);
            const auto& diagnostics = solver.normalization_diagnostics();
            check(solver.active_constraint_count() == 1,
                  std::string(test_case.label) + ", tolerance-equivalent: must collapse");
            check(diagnostics.tolerance_equivalent_duplicates == 1,
                  std::string(test_case.label) +
                      ", tolerance-equivalent: must count one such duplicate");
            check(diagnostics.exact_duplicates == 0,
                  std::string(test_case.label) +
                      ", tolerance-equivalent: must not count an exact duplicate");
            check(diagnostics.reversed_duplicates == test_case.expected_reversed,
                  std::string(test_case.label) + ", tolerance-equivalent: expected " +
                      std::to_string(test_case.expected_reversed) +
                      " reversed duplicate(s), got " +
                      std::to_string(diagnostics.reversed_duplicates));
        }
    }
}


// --- Per-component identification state ------------------------------------

// The report-level view is DERIVED from the components, so it cannot describe a
// mixed result as fully identified. This checks that directly by assembling a
// report by hand: one rank-deficient component whose redundant rows were
// identified, one that fell back to naming its whole membership.
void test_mixed_identification_cannot_be_reported_as_identified() {
    gmd::ConstraintRankReport report;
    report.constraint_count = 7;
    report.rank = 5;

    gmd::ConstraintComponent identified;
    identified.atom_tags = {0, 1, 2};
    identified.constraint_indices = {0, 1, 2};
    identified.rank = 2;
    identified.redundant_constraints_identified = true;
    identified.redundant_constraint_indices = {2};

    gmd::ConstraintComponent fell_back;
    fell_back.atom_tags = {5, 6, 7, 8};
    fell_back.constraint_indices = {3, 4, 5, 6};
    fell_back.rank = 3;
    fell_back.redundant_constraints_identified = false;
    fell_back.redundant_constraint_indices = {3, 4, 5, 6};   // the whole component

    report.components = {identified, fell_back};

    check(!report.redundant_identified(),
          "a report with one identified and one fallen-back component must NOT claim "
          "to be fully identified");
    check(report.dependent_constraints().size() == 5,
          "and its constraint list must be the concatenation of the per-component "
          "lists, 1 + 4 = 5, got " +
              std::to_string(report.dependent_constraints().size()));

    // Order the other way round: the answer cannot depend on which component the
    // failure happened to be.
    report.components = {fell_back, identified};
    check(!report.redundant_identified(),
          "and must still not claim identification when the fallen-back component "
          "comes first");

    // With both identified it is true, so the check above is not vacuous.
    report.components = {identified, identified};
    check(report.redundant_identified(),
          "two identified components must report identification");

    // Vacuously true when nothing is deficient.
    gmd::ConstraintComponent healthy;
    healthy.constraint_indices = {0, 1};
    healthy.rank = 2;
    gmd::ConstraintRankReport clean;
    clean.components = {healthy};
    check(clean.redundant_identified(),
          "an independent set has nothing to identify, so the flag is vacuously true");
    check(clean.dependent_constraints().empty(),
          "and names no redundant constraints");
}

// Two disconnected rank-deficient components in one system. Each must carry its
// own identification state and its own redundant rows, and the report-level view
// must summarise both.
void test_two_disconnected_deficient_components() {
    // Component A: five atoms, every pair -- ten constraints, rank 9.
    // Component B: four collinear atoms, chain plus closing bond -- rank 3.
    std::vector<Vec3> positions{{10.0, 10.0, 10.0}, {11.5, 10.1, 10.2},
                                {10.6, 11.4, 10.3}, {10.4, 10.2, 11.6},
                                {11.7, 11.3, 11.5},
                                {30.0, 30.0, 30.0}, {31.3, 30.0, 30.0},
                                {32.9, 30.0, 30.0}, {34.6, 30.0, 30.0}};
    gmd::System system = make_system(positions);

    std::vector<std::array<int, 2>> pairs;
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) pairs.push_back({i, j});
    }
    for (const auto& pair : std::vector<std::array<int, 2>>{{5, 6}, {6, 7}, {7, 8}, {5, 8}}) {
        pairs.push_back(pair);
    }

    const auto constraints = hold(system, pairs);
    const auto report =
        gmd::ConstraintSolver(constraints, settings()).analyze_independence(system);

    check(!report.independent, "both components are deficient, so the set is not");
    check(report.components.size() == 2,
          "two molecules sharing no atom must be two components, got " +
              std::to_string(report.components.size()));
    check(report.constraint_count == 14 && report.rank == 12,
          "10 + 4 constraints of rank 9 + 3; got " + std::to_string(report.constraint_count) +
              " and " + std::to_string(report.rank));

    std::size_t named = 0;
    for (const auto& component : report.components) {
        check(!component.independent(), "each component here is rank deficient");
        // The invariant that must hold in BOTH branches.
        if (component.redundant_constraints_identified) {
            check(component.redundant_constraint_indices.size() ==
                      component.constraint_indices.size() - component.rank,
                  "an identified component must name exactly count - rank rows, got " +
                      std::to_string(component.redundant_constraint_indices.size()));
        } else {
            check(component.redundant_constraint_indices == component.constraint_indices,
                  "a fallen-back component must name its entire membership");
        }
        named += component.redundant_constraint_indices.size();
    }
    check(report.dependent_constraints().size() == named,
          "the report-level list must be the concatenation of the per-component lists");
    check(report.redundant_identified() ==
              (report.components[0].redundant_constraints_identified &&
               report.components[1].redundant_constraints_identified),
          "and the report-level flag must be the conjunction of the per-component ones");

    std::cout << "  two disconnected deficient components: ranks";
    for (const auto& component : report.components) {
        std::cout << ' ' << component.rank << '/' << component.constraint_indices.size()
                  << (component.redundant_constraints_identified ? " (identified)"
                                                                 : " (fallback)");
    }
    std::cout << ", " << report.dependent_constraints().size() << " redundant row(s) named\n";
}

// --- Canonical, order-independent tie-breaking ----------------------------

// Symmetric fixtures produce genuinely tied pivot candidates. Breaking those
// ties by column index would make the answer depend on the order the constraints
// were listed in; breaking them by the canonical atom pair does not. Every
// permutation of the input must name the same PHYSICAL constraint redundant.
void test_symmetric_fixtures_name_the_same_pair_under_every_permutation() {
    struct Fixture {
        const char* label;
        std::vector<Vec3> positions;
        std::vector<std::array<int, 2>> pairs;
        double mass;   // one mass for every atom, so nothing breaks the symmetry
        std::string label_storage;   // used for the rotated variants
    };

    std::vector<Fixture> fixtures;
    // Four equally spaced collinear atoms, equal masses: (0,1) and (2,3) are
    // mirror images of each other and tie exactly.
    fixtures.push_back({"four equally spaced collinear atoms, equal masses",
                        {{10.0, 10.0, 10.0}, {11.5, 10.0, 10.0},
                         {13.0, 10.0, 10.0}, {14.5, 10.0, 10.0}},
                        {{0, 1}, {1, 2}, {2, 3}, {0, 3}},
                        12.0, {}});
    // A regular pentagon with every pair constrained: ten constraints over a
    // planar body, rank 7, and the fivefold symmetry ties many candidates.
    {
        std::vector<Vec3> ring;
        for (int k = 0; k < 5; ++k) {
            const double angle = 2.0 * M_PI * static_cast<double>(k) / 5.0;
            ring.push_back({20.0 + 1.6 * std::cos(angle), 20.0 + 1.6 * std::sin(angle), 20.0});
        }
        std::vector<std::array<int, 2>> all_pairs;
        for (int i = 0; i < 5; ++i) {
            for (int j = i + 1; j < 5; ++j) all_pairs.push_back({i, j});
        }
        fixtures.push_back({"regular pentagon, every pair constrained, equal masses",
                            ring, all_pairs, 12.0, {}});
    }

    // Each fixture is run twice: axis-aligned, where the symmetric columns reach
    // bit-identical residual norms, and rotated into a general orientation,
    // where they reach the same value by different arithmetic and agree only to
    // round-off. The second is what makes the RELATIVE tie criterion necessary:
    // an absolute comparison would call those columns distinct and hand the
    // choice back to input order.
    const Vec3 tilt_axis = [] {
        Vec3 a{1.0, -2.0, 3.0};
        const double n = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
        for (double& c : a) c /= n;
        return a;
    }();
    auto rotate = [&](const Vec3& v, const Vec3& centre) {
        const double angle = 0.7;
        const double c = std::cos(angle);
        const double sn = std::sin(angle);
        const Vec3 r{v[0] - centre[0], v[1] - centre[1], v[2] - centre[2]};
        const double kd = tilt_axis[0]*r[0] + tilt_axis[1]*r[1] + tilt_axis[2]*r[2];
        const Vec3 k{tilt_axis[1]*r[2] - tilt_axis[2]*r[1],
                     tilt_axis[2]*r[0] - tilt_axis[0]*r[2],
                     tilt_axis[0]*r[1] - tilt_axis[1]*r[0]};
        Vec3 out{};
        for (std::size_t d = 0; d < 3; ++d) {
            out[d] = centre[d] + r[d] * c + k[d] * sn + tilt_axis[d] * kd * (1.0 - c);
        }
        return out;
    };

    std::vector<Fixture> all = fixtures;
    for (const auto& fixture : fixtures) {
        Fixture tilted = fixture;
        tilted.label_storage = std::string(fixture.label) + " (rotated into a general "
                                                            "orientation)";
        const Vec3 centre{fixture.positions.front()[0], fixture.positions.front()[1],
                          fixture.positions.front()[2]};
        for (auto& position : tilted.positions) position = rotate(position, centre);
        all.push_back(tilted);
    }

    for (auto& fixture_ref : all) {
        const Fixture& fixture = fixture_ref;
        const std::string label = fixture.label_storage.empty()
                                      ? std::string(fixture.label)
                                      : fixture.label_storage;
        gmd::System system;
        system.resize(fixture.positions.size(), fixture.positions.size());
        gmd::Box box;
        box.set_lengths({60.0, 60.0, 60.0});
        system.set_box(box);
        for (std::size_t i = 0; i < fixture.positions.size(); ++i) {
            system.mutable_masses()[i] = fixture.mass;
            system.mutable_atom_tags()[i] = static_cast<int>(i);
            system.mutable_coordinates()[i] = fixture.positions[i];
        }

        auto redundant_pairs = [&](const std::vector<std::array<int, 2>>& order) {
            const auto list = hold(system, order);
            const auto report =
                gmd::ConstraintSolver(list, settings()).analyze_independence(system);
            std::vector<std::array<int, 2>> out;
            for (std::size_t index : report.dependent_constraints()) {
                out.push_back({std::min(list[index].i, list[index].j),
                               std::max(list[index].i, list[index].j)});
            }
            std::sort(out.begin(), out.end());
            return std::make_pair(report.rank, out);
        };

        auto order = fixture.pairs;
        std::sort(order.begin(), order.end());
        const auto reference = redundant_pairs(order);
        check(!reference.second.empty(),
              label + ": must be rank deficient, or there is nothing to name");

        // Every permutation, up to a cap for the larger fixture.
        std::size_t permutations = 0;
        constexpr std::size_t limit = 5000;
        do {
            const auto measured = redundant_pairs(order);
            check(measured.first == reference.first,
                  label + ": rank changed under a permutation, " +
                      std::to_string(reference.first) + " -> " +
                      std::to_string(measured.first));
            check(measured.second == reference.second,
                  label +
                      ": a permutation of the input named a different physical "
                      "constraint redundant, which means the tie-break is still "
                      "positional");
            ++permutations;
        } while (permutations < limit && std::next_permutation(order.begin(), order.end()));

        std::cout << "  " << label << ": rank " << reference.first << ", redundant";
        for (const auto& pair : reference.second) {
            std::cout << " (" << pair[0] << "," << pair[1] << ")";
        }
        std::cout << ", identical across " << permutations << " permutations\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif

    test_independent_chain();
    test_independent_triangle();
    test_disjoint_components();
    test_exact_and_reversed_duplicates();
    test_reversed_duplicate_orientations();
    test_tolerance_equivalent_duplicate();
    test_conflicting_targets();
    test_redundant_closed_network();
    test_over_constrained_collinear();
    test_geometric_dependence_without_a_constrained_triple();
    test_degenerate_zero_length_gradient();
    test_near_singular_but_valid();
    test_rank_tolerance_is_scale_aware();
    test_degrees_of_freedom_uses_the_verified_rank();
    test_initialize_rejects_a_dependent_set();
    test_full_rank_input_projecting_onto_a_degenerate_target();
    test_authoritative_report_describes_the_projected_geometry();
    test_off_manifold_but_valid_target_projects_and_stays_independent();
    test_singular_values_against_constructed_matrices();
    test_singular_values_are_permutation_invariant();
    test_non_finite_input_is_rejected();
    test_stress_larger_and_ill_scaled_components();
    test_analysis_is_invariant_to_constraint_order();
    test_mixed_identification_cannot_be_reported_as_identified();
    test_two_disconnected_deficient_components();
    test_symmetric_fixtures_name_the_same_pair_under_every_permutation();

    if (failures != 0) {
        std::cerr << "[constraint independence] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[constraint independence] all checks passed\n";
    return 0;
}
