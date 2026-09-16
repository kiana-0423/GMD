// The constraint-rank FALLBACK branch, exercised end to end through the real
// analysis pipeline.
//
// WHAT THE FALLBACK IS. analyze_independence() decides two different things
// about a rank-deficient component, by two different methods:
//
//   HOW MANY rows are redundant -- from the singular values of the
//       mass-weighted Jacobian, counted against
//       tol = max(rows, cols) * eps * sigma_max;
//   WHICH rows they are -- from a pivoted modified Gram-Schmidt selection,
//       which accepts a column when its residual after projecting out the basis
//       so far still clears the same tol.
//
// The second is the less reliable question, so the first is treated as
// authoritative and the two are CROSS-CHECKED: when the number of columns MGS
// accepted does not equal the singular-value rank, the component cannot say
// which of its rows are the redundant ones, and it falls back to naming its
// ENTIRE membership. `redundant_constraints_identified` records which happened.
//
// WHY THIS FILE EXISTS. That branch was previously covered only by assembling a
// ConstraintRankReport by hand at the API level, which cannot show that the
// production pipeline ever reaches it, nor that the report it then produces is
// self-consistent. The fixture below drives it through
// ConstraintSolver::analyze_independence() on real coordinates.
//
// THE FIXTURE, and why it is shaped this way. Four atoms with all six pairwise
// distances constrained are exactly rigid: 6 constraints against 3*4 - 6 = 6
// internal degrees of freedom, so the set is INDEPENDENT at a generic geometry.
// Make the four atoms COPLANAR and the Cayley-Menger determinant relating six
// distances of four coplanar points vanishes, giving one relation among them:
// the rank drops to 5 and the set becomes dependent.
//
// Sweeping the out-of-plane height h therefore carries the smallest singular
// value CONTINUOUSLY from O(1) down through the rank tolerance, and it is in
// that crossing region that the two methods disagree -- MGS's residual norms
// over-estimate sigma_min, so there is a band where MGS still accepts six
// columns while the singular values already say five.
//
// This shape is chosen deliberately over a near-collinear one: any three
// mutually constrained atoms that are close to collinear are rejected far
// earlier, by reject_degenerate_target_triangles(), which works on the target
// distances alone. Four near-coplanar atoms contain no such triple.
//
// HONEST SCOPE. The band is narrow -- roughly h in [2.6e-15, 4.4e-15] for this
// fixture, i.e. a configuration a few parts in 1e15 away from coplanar. This is
// a real production-path fixture and not an injected fault, but it is not a
// configuration anyone would arrive at from physical modelling; it is what
// reaching this branch actually requires. The test therefore sweeps the band
// rather than betting on one h, so it does not depend on a single value
// surviving a change of compiler or library.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;
using Vec3 = gmd::System::Vec3;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[constraint rank fallback] " << message << '\n';
        ++failures;
    }
}

gmd::ConstraintSettings settings() {
    gmd::ConstraintSettings s;
    s.tolerance = 1.0e-12;
    s.max_iterations = 500;
    return s;
}

// The base quadrilateral. Deliberately irregular in the plane, so the only
// degeneracy in play is the out-of-plane one the sweep controls.
std::vector<Vec3> quad(double height) {
    return {{0.0, 0.0, 0.0},
            {1.13, 0.0, 0.0},
            {0.21, 1.07, 0.0},
            {0.79, 0.62, height}};
}

const std::vector<std::array<int, 2>> kAllSixPairs = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

std::vector<Vec3> rotate(std::vector<Vec3> points, double a, double b, double c) {
    for (auto& q : points) {
        const double x = q[0], y = q[1], z = q[2];
        const double x1 = x * std::cos(a) - y * std::sin(a);
        const double y1 = x * std::sin(a) + y * std::cos(a);
        const double y2 = y1 * std::cos(b) - z * std::sin(b);
        const double z2 = y1 * std::sin(b) + z * std::cos(b);
        const double z3 = z2 * std::cos(c) - x1 * std::sin(c);
        const double x3 = z2 * std::sin(c) + x1 * std::cos(c);
        q = {x3, y2, z3};
    }
    return points;
}

std::vector<Vec3> translate(std::vector<Vec3> points, const Vec3& by) {
    for (auto& q : points) {
        for (std::size_t d = 0; d < 3; ++d) q[d] += by[d];
    }
    return points;
}

double distance(const std::vector<Vec3>& p, int a, int b) {
    double sum = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        const double delta = p[static_cast<std::size_t>(a)][d] - p[static_cast<std::size_t>(b)][d];
        sum += delta * delta;
    }
    return std::sqrt(sum);
}

// `tag_of_slot[k]` is the global tag stored in slot k, so the storage order and
// the tag order can be permuted independently.
gmd::System make_system(const std::vector<Vec3>& positions,
                        const std::vector<int>& tag_of_slot,
                        double mass_ratio) {
    gmd::System system;
    system.resize(positions.size(), positions.size());
    gmd::Box box;
    box.set_lengths({400.0, 440.0, 480.0});
    system.set_box(box);
    auto masses = system.mutable_masses();
    auto tags = system.mutable_atom_tags();
    auto coordinates = system.mutable_coordinates();
    for (std::size_t k = 0; k < positions.size(); ++k) {
        tags[k] = tag_of_slot[k];
        coordinates[k] = positions[k];
        // Mass keyed to the TAG, not the slot, so a storage permutation moves
        // the atom and its mass together.
        masses[k] = std::pow(mass_ratio, static_cast<double>(tag_of_slot[k] % 3));
    }
    return system;
}

// Constraints holding each pair at the distance it currently sits at, expressed
// against global tags. `slot_of_tag` inverts the storage permutation.
std::vector<gmd::BondConstraint> hold(const std::vector<Vec3>& positions,
                                      const std::vector<std::array<int, 2>>& pairs,
                                      const std::vector<int>& slot_of_tag) {
    std::vector<gmd::BondConstraint> constraints;
    for (const auto& pair : pairs) {
        constraints.push_back({pair[0], pair[1],
                               distance(positions, slot_of_tag[static_cast<std::size_t>(pair[0])],
                                        slot_of_tag[static_cast<std::size_t>(pair[1])])});
    }
    return constraints;
}

struct Outcome {
    bool analysed = false;
    bool fell_back = false;
    bool independent = true;
    std::size_t rank = 0;
    std::size_t columns = 0;
    std::size_t dependent_count = 0;
    bool report_identified = true;
    double smallest_singular_value = 0.0;
    double rank_tolerance = 0.0;
};

// Runs the real analysis and summarises the one deficient component.
Outcome analyse(const std::vector<Vec3>& positions,
                const std::vector<std::array<int, 2>>& pairs,
                const std::vector<int>& tag_of_slot,
                double mass_ratio) {
    std::vector<int> slot_of_tag(tag_of_slot.size());
    for (std::size_t k = 0; k < tag_of_slot.size(); ++k) {
        slot_of_tag[static_cast<std::size_t>(tag_of_slot[k])] = static_cast<int>(k);
    }
    const gmd::ConstraintSolver solver(hold(positions, pairs, slot_of_tag), settings());
    const auto system = make_system(positions, tag_of_slot, mass_ratio);
    const auto report = solver.analyze_independence(system);

    Outcome out;
    out.report_identified = report.redundant_identified();
    for (const auto& component : report.components) {
        if (component.independent()) continue;
        out.analysed = true;
        out.independent = false;
        out.rank = component.rank;
        out.columns = component.constraint_indices.size();
        out.dependent_count = component.redundant_constraint_indices.size();
        out.fell_back = !component.redundant_constraints_identified;
        out.smallest_singular_value = component.smallest_singular_value;
        out.rank_tolerance = component.rank_tolerance;
    }
    if (!out.analysed) out.independent = true;
    return out;
}

// The invariant every deficient component must satisfy, whichever branch it
// took. Asserted everywhere the analysis runs, not only on the fallback path.
void check_component_invariant(const Outcome& o, const std::string& where) {
    if (o.independent) return;
    if (o.fell_back) {
        check(o.dependent_count == o.columns,
              where + ": a component that fell back must name its WHOLE membership, "
                      "so dependent count " + std::to_string(o.dependent_count) +
                  " must equal the column count " + std::to_string(o.columns));
        check(!o.report_identified,
              where + ": a report containing a fallen-back component must not claim "
                      "identification");
    } else {
        check(o.dependent_count == o.columns - o.rank,
              where + ": an identified component must name exactly count - rank = " +
                  std::to_string(o.columns - o.rank) + " redundant constraint(s), got " +
                  std::to_string(o.dependent_count));
    }
}

// ---------------------------------------------------------------------------
// 1. The fallback is reachable through the production path.
// ---------------------------------------------------------------------------
// The sweep is the assertion. A single hand-picked height would be a bet that
// one floating-point value lands in a narrow band on every compiler; sweeping
// the band and requiring at least one hit is the same claim without the bet.
void test_fallback_is_reached_end_to_end() {
    int hits = 0;
    int deficient = 0;
    double lowest_hit = 0.0, highest_hit = 0.0;

    for (int step = 0; step <= 40; ++step) {
        // 2.0e-15 .. 5.0e-15, logarithmically.
        const double height = 2.0e-15 * std::pow(5.0 / 2.0, static_cast<double>(step) / 40.0);
        const auto outcome = analyse(quad(height), kAllSixPairs, {0, 1, 2, 3}, 1.5);
        if (outcome.independent) continue;
        ++deficient;
        check_component_invariant(outcome, "sweep h=" + std::to_string(height));
        if (outcome.fell_back) {
            ++hits;
            if (lowest_hit == 0.0) lowest_hit = height;
            highest_hit = height;
            check(outcome.rank == 5 && outcome.columns == 6,
                  "the fallback fixture must be 6 constraints at rank 5, got " +
                      std::to_string(outcome.columns) + " at rank " +
                      std::to_string(outcome.rank));
        }
    }

    check(deficient > 0,
          "the sweep must produce rank-deficient components; none appeared, so the "
          "fixture no longer straddles the coplanar degeneracy");
    check(hits > 0,
          "no height in [2e-15, 5e-15] made the pivoted rank selection disagree with "
          "the singular values, so the fallback branch was never reached through the "
          "production path. If this fires after a change to the rank tolerance or to "
          "the selection, re-derive the band with the sweep in the file comment "
          "rather than widening it blindly");
    if (hits > 0) {
        std::cout << "[constraint rank fallback] production-path fallback reached at "
                  << hits << " of " << deficient << " deficient heights, h in ["
                  << lowest_hit << ", " << highest_hit << "]\n";
    }
}

// ---------------------------------------------------------------------------
// 2. The fallback decision does not depend on how the input was written down.
// ---------------------------------------------------------------------------
void test_fallback_is_invariant() {
    // A height inside the band, used for the invariance matrix. Picked from the
    // middle of the measured band so that round-off differences between
    // platforms stay inside it.
    const double height = 3.311e-15;
    const auto base = analyse(quad(height), kAllSixPairs, {0, 1, 2, 3}, 1.5);
    check(base.fell_back,
          "the invariance fixture must itself reach the fallback; if this fires the "
          "band has moved and the height above needs re-deriving from the sweep");
    check_component_invariant(base, "invariance base");

    const auto same = [&](const Outcome& o, const std::string& what) {
        check(o.independent == base.independent && o.fell_back == base.fell_back &&
                  o.rank == base.rank && o.columns == base.columns &&
                  o.dependent_count == base.dependent_count,
              what + " changed the fallback outcome: expected fell_back=" +
                  std::to_string(base.fell_back) + " rank=" + std::to_string(base.rank) +
                  " dependent=" + std::to_string(base.dependent_count) + ", got fell_back=" +
                  std::to_string(o.fell_back) + " rank=" + std::to_string(o.rank) +
                  " dependent=" + std::to_string(o.dependent_count));
        check_component_invariant(o, what);
    };

    auto reversed = kAllSixPairs;
    std::reverse(reversed.begin(), reversed.end());
    same(analyse(quad(height), reversed, {0, 1, 2, 3}, 1.5), "reversing the constraint order");

    auto rotated_list = kAllSixPairs;
    std::rotate(rotated_list.begin(), rotated_list.begin() + 3, rotated_list.end());
    same(analyse(quad(height), rotated_list, {0, 1, 2, 3}, 1.5), "rotating the constraint list");

    // Storage permutation: the same atoms and tags, held in different slots.
    const auto points = quad(height);
    same(analyse({points[3], points[1], points[0], points[2]}, kAllSixPairs, {3, 1, 0, 2}, 1.5),
         "permuting the atom storage order");
    same(analyse({points[2], points[3], points[1], points[0]}, kAllSixPairs, {2, 3, 1, 0}, 1.5),
         "permuting the atom storage order the other way");

    same(analyse(rotate(quad(height), 0.37, 0.61, 0.23), kAllSixPairs, {0, 1, 2, 3}, 1.5),
         "rotating the fixture");
    same(analyse(translate(quad(height), {37.25, -11.5, 4.75}), kAllSixPairs, {0, 1, 2, 3}, 1.5),
         "translating the fixture");

    for (double mass_ratio : {1.0, 16.0, 200.0}) {
        same(analyse(quad(height), kAllSixPairs, {0, 1, 2, 3}, mass_ratio),
             "mass ratio " + std::to_string(mass_ratio));
    }
}

// ---------------------------------------------------------------------------
// 3. Two disconnected components, only one of which falls back.
// ---------------------------------------------------------------------------
// The report-level flag is the AND over components, so a mixed result must not
// be reported as identified. This is the same claim the hand-assembled report
// test makes, made instead by the real pipeline.
void test_mixed_components_through_the_pipeline() {
    const double height = 3.311e-15;
    auto points = quad(height);
    // A second, well-conditioned rigid triangle far away: independent, so it
    // contributes no redundancy of its own.
    points.push_back({120.0, 0.0, 0.0});
    points.push_back({121.3, 0.0, 0.0});
    points.push_back({120.6, 1.1, 0.0});

    auto pairs = kAllSixPairs;
    pairs.push_back({4, 5});
    pairs.push_back({4, 6});
    pairs.push_back({5, 6});

    const gmd::ConstraintSolver solver(
        hold(points, pairs, {0, 1, 2, 3, 4, 5, 6}), settings());
    const auto system = make_system(points, {0, 1, 2, 3, 4, 5, 6}, 1.5);
    const auto report = solver.analyze_independence(system);

    check(report.components.size() == 2,
          "the fixture must split into exactly two constraint components, got " +
              std::to_string(report.components.size()));

    std::size_t deficient = 0, independent = 0, fell_back = 0;
    for (const auto& component : report.components) {
        if (component.independent()) {
            ++independent;
            check(component.redundant_constraint_indices.empty(),
                  "an independent component must name no redundant constraints");
        } else {
            ++deficient;
            if (!component.redundant_constraints_identified) ++fell_back;
        }
    }
    check(independent == 1 && deficient == 1,
          "expected one independent and one deficient component, got " +
              std::to_string(independent) + " and " + std::to_string(deficient));
    check(fell_back == 1,
          "the near-coplanar component must still fall back when analysed alongside "
          "an independent one");
    check(!report.redundant_identified(),
          "a report with a fallen-back component must not claim identification, even "
          "when its other component is perfectly healthy");
    // The independent component must not have contributed any dependent indices,
    // so the report-level list is exactly the fallen-back component's membership.
    check(report.dependent_constraints().size() == 6,
          "the report's dependent list must be the fallen-back component's whole "
          "membership and nothing else, expected 6, got " +
              std::to_string(report.dependent_constraints().size()));
}

// ---------------------------------------------------------------------------
// 3b. One component IDENTIFIED-deficient, one FALLEN BACK.
// ---------------------------------------------------------------------------
// The sharper aggregation case, and the one that actually distinguishes AND
// from OR. With a fallen-back component next to an INDEPENDENT one, an OR over
// "identified" and an AND agree -- there is no identified-deficient component to
// disagree about. It takes a component that IS cleanly identified, alongside one
// that fell back, for the two to differ.
//
// Five atoms with all ten pairwise distances constrained supply that: 10
// constraints against 3*5 - 6 = 9 internal degrees of freedom, so the set is
// structurally deficient by exactly one at any generic geometry, and its
// smallest singular value sits far enough below the rank tolerance that the two
// rank methods agree cleanly -- it is identified, not a fallback.
void test_identified_and_fallen_back_components_together() {
    const double height = 3.311e-15;
    auto points = quad(height);                       // falls back, 6 constraints
    // A five-atom cage, far away so it is a separate component.
    points.push_back({200.0, 0.0, 0.0});
    points.push_back({201.13, 0.0, 0.0});
    points.push_back({200.21, 1.07, 0.0});
    points.push_back({200.97, 0.91, 0.0});
    points.push_back({200.53, 0.47, 0.83});

    auto pairs = kAllSixPairs;
    for (int a = 4; a < 9; ++a) {
        for (int b = a + 1; b < 9; ++b) pairs.push_back({a, b});
    }

    std::vector<int> tags(9);
    std::iota(tags.begin(), tags.end(), 0);
    const gmd::ConstraintSolver solver(hold(points, pairs, tags), settings());
    const auto system = make_system(points, tags, 1.5);
    const auto report = solver.analyze_independence(system);

    std::size_t identified_deficient = 0, fallen_back = 0;
    std::size_t dependent_total = 0;
    for (const auto& component : report.components) {
        if (component.independent()) continue;
        dependent_total += component.redundant_constraint_indices.size();
        if (component.redundant_constraints_identified) {
            ++identified_deficient;
            check(component.redundant_constraint_indices.size() ==
                      component.constraint_indices.size() - component.rank,
                  "the identified component must name exactly count - rank redundant "
                  "constraint(s)");
        } else {
            ++fallen_back;
        }
    }

    check(identified_deficient == 1 && fallen_back == 1,
          "this fixture must produce exactly one IDENTIFIED-deficient component and "
          "one FALLEN-BACK component -- that combination is what distinguishes an AND "
          "aggregation from an OR -- got " + std::to_string(identified_deficient) +
              " identified and " + std::to_string(fallen_back) + " fallen back");
    check(!report.redundant_identified(),
          "a report holding one identified component AND one that fell back must NOT "
          "claim identification; an OR over the components would wrongly say it does");
    // 6 from the fallen-back component's whole membership, 1 from the cage.
    check(dependent_total == 7,
          "the dependent list must be the fallen-back component's whole membership (6) "
          "plus the cage's single identified redundant constraint (1), expected 7, got " +
              std::to_string(dependent_total));
}

// ---------------------------------------------------------------------------
// 4. A fallen-back component still REJECTS the run, and leaves nothing behind.
// ---------------------------------------------------------------------------
// The fallback is about which rows are named, not about whether the set is
// usable: an untrusted identification must never be turned into an accepted run
// or into a degree-of-freedom count.
void test_fallback_still_rejects_and_leaves_no_state() {
    const double height = 3.311e-15;
    const auto points = quad(height);
    const auto constraints = hold(points, kAllSixPairs, {0, 1, 2, 3});

    auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, settings());
    auto system = make_system(points, {0, 1, 2, 3}, 1.5);

    bool threw = false;
    std::string what;
    try {
        solver->require_independent(system, "the test geometry");
    } catch (const std::invalid_argument& error) {
        threw = true;
        what = error.what();
    }
    check(threw, "require_independent() must reject a rank-deficient set even when the "
                 "redundant rows could not be identified");
    check(what.find("not independent") != std::string::npos,
          "the rejection must say the set is not independent; got: " + what);
    check(what.find("rank 5") != std::string::npos,
          "the rejection must report the rank it found; got: " + what);

    // Through the integrator: initialize() must throw, and must not leave a
    // thermostat initialised or a velocity field behind.
    auto thermostat = std::make_shared<gmd::NoseHooverThermostat>(100.0);
    gmd::VelocityVerletIntegrator integrator(1.0);
    integrator.set_constraint_solver(solver);
    integrator.set_thermostat(thermostat);

    auto velocities = system.mutable_velocities();
    for (auto& v : velocities) v = {0.0, 0.0, 0.0};

    gmd::RuntimeContext runtime;
    bool integrator_threw = false;
    try {
        integrator.initialize(system, runtime);
    } catch (const std::invalid_argument&) {
        integrator_threw = true;
    }
    check(integrator_threw,
          "the integrator must refuse to initialize on a set whose rank could not be "
          "trusted");

    double speed = 0.0;
    for (const auto& v : system.velocities()) {
        for (std::size_t d = 0; d < 3; ++d) speed += std::abs(v[d]);
    }
    check(speed == 0.0,
          "a rejected initialization must leave no velocity field behind, found total "
          "absolute velocity " + std::to_string(speed));
    check(integrator.constraint_rank_report().components.empty(),
          "a rejected initialization must not leave a completed rank report behind");
}

}  // namespace

int main() {
    test_fallback_is_reached_end_to_end();
    test_fallback_is_invariant();
    test_mixed_components_through_the_pipeline();
    test_identified_and_fallen_back_components_together();
    test_fallback_still_rejects_and_leaves_no_state();

    if (failures == 0) {
        std::cout << "[constraint rank fallback] all checks passed\n";
    } else {
        std::cerr << "[constraint rank fallback] " << failures << " check(s) failed\n";
    }
    return failures == 0 ? 0 : 1;
}
