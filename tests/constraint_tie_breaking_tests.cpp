// Rank-column tie-breaking: is the RELATIVE near-tie window necessary, and is
// the result invariant under everything that should not matter?
//
// THE RULE UNDER TEST. select_independent_columns() picks, at each step, the
// column whose residual is largest after projecting out the basis so far. When
// several columns are effectively tied it must not break the tie by position in
// the input, because that would make the answer depend on the order the
// constraints happened to be listed in. It instead takes the smallest CANONICAL
// ATOM PAIR (min tag, max tag), which is a property of the physical constraint.
//
// A tie is decided RELATIVELY:
//
//     cutoff = best_norm * (1 - sqrt(eps))          sqrt(eps) = 1.4901e-08
//
// and every live column at or above the cutoff is an equally good pivot.
//
// WHY RELATIVE, AND WHY THIS FILE EXISTS. Two columns that are physically
// equivalent reach their residual norms by different arithmetic paths, so they
// agree only to within round-off accumulated over the projections. The previous
// coverage used SYMMETRIC, AXIS-ALIGNED fixtures whose tied norms come out
// BIT-IDENTICAL; for those, an absolute criterion and a relative one agree, so
// they cannot show the relative window is needed. test_axis_aligned_tie_does_
// not_discriminate() below pins that limitation in place so it is not mistaken
// for coverage.
//
// The fixture that does discriminate is a GENERICALLY ROTATED square with
// UNEQUAL masses. The rotation destroys bit-identity without touching the
// physical symmetry. Measured at the deciding pivot step, four columns sit
// within 1.2e-16 to 3.7e-16 of the leader:
//
//     (0,1)  1.80277563773199456e+00     relative gap 3.695e-16
//     (0,3)  1.80277563773199523e+00     leader
//     (1,2)  1.80277563773199456e+00     relative gap 3.695e-16
//     (2,3)  1.80277563773199501e+00     relative gap 1.232e-16
//
// Those gaps are NOT zero -- so an absolute 1e-300 criterion calls the four
// columns distinct and hands the choice to raw arithmetic -- and they are some
// eight orders of magnitude below sqrt(eps), so the relative rule calls them
// tied and the canonical key decides. The numbers above are recorded for the
// reader; the assertions below are on the PHYSICAL CONSTRAINT the production
// path selects, not on any recomputed norm.
//
// NECESSITY is demonstrated outside this file, by restoring the old absolute
// 1e-300 criterion in a disposable worktree: the rotated fixture then names
// (1,3) instead of (2,3), and uniform coordinate scaling moves it again to
// (1,2) and (0,1). See validation/README.md for the recorded results.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;
using Vec3 = gmd::System::Vec3;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[constraint tie breaking] " << message << '\n';
        ++failures;
    }
}

gmd::ConstraintSettings settings() {
    gmd::ConstraintSettings s;
    s.tolerance = 1.0e-12;
    s.max_iterations = 500;
    return s;
}

const std::vector<std::array<int, 2>> kAllSixPairs = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

double distance(const std::vector<Vec3>& p, int a, int b) {
    double sum = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        const double delta = p[static_cast<std::size_t>(a)][d] - p[static_cast<std::size_t>(b)][d];
        sum += delta * delta;
    }
    return std::sqrt(sum);
}

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

std::vector<Vec3> scale(std::vector<Vec3> points, double factor) {
    for (auto& q : points) {
        for (std::size_t d = 0; d < 3; ++d) q[d] *= factor;
    }
    return points;
}

// The physical identity of the constraints named redundant, as "(i,j)" pairs in
// canonical order. This is what the assertions compare -- naming the physical
// constraint rather than an index into a list whose order the test is varying.
std::string redundant_pairs(const std::vector<Vec3>& positions,
                            const std::vector<std::array<int, 2>>& pairs,
                            const std::vector<int>& tag_of_slot,
                            double mass_ratio) {
    std::vector<int> slot_of_tag(tag_of_slot.size());
    for (std::size_t k = 0; k < tag_of_slot.size(); ++k) {
        slot_of_tag[static_cast<std::size_t>(tag_of_slot[k])] = static_cast<int>(k);
    }

    gmd::System system;
    system.resize(positions.size(), positions.size());
    gmd::Box box;
    box.set_lengths({4000.0, 4400.0, 4800.0});
    system.set_box(box);
    auto masses = system.mutable_masses();
    auto tags = system.mutable_atom_tags();
    auto coordinates = system.mutable_coordinates();
    for (std::size_t k = 0; k < positions.size(); ++k) {
        tags[k] = tag_of_slot[k];
        coordinates[k] = positions[k];
        masses[k] = std::pow(mass_ratio, static_cast<double>(tag_of_slot[k] % 2));
    }

    std::vector<gmd::BondConstraint> constraints;
    for (const auto& pair : pairs) {
        constraints.push_back({pair[0], pair[1],
                               distance(positions, slot_of_tag[static_cast<std::size_t>(pair[0])],
                                        slot_of_tag[static_cast<std::size_t>(pair[1])])});
    }

    const gmd::ConstraintSolver solver(constraints, settings());
    const auto report = solver.analyze_independence(system);

    std::string out;
    for (const auto& component : report.components) {
        if (component.independent()) continue;
        if (!component.redundant_constraints_identified) out += "FALLBACK:";
        for (const auto index : component.redundant_constraint_indices) {
            const auto& c = solver.constraints()[index];
            out += "(" + std::to_string(c.i) + "," + std::to_string(c.j) + ")";
        }
    }
    return out.empty() ? "(independent)" : out;
}

// A unit square. Coplanar, so the Cayley-Menger relation among the six
// distances of four coplanar points makes the set rank 5 of 6; the square's
// symmetry makes several pivot candidates physically equivalent.
const std::vector<Vec3> kSquare = {{0.0, 0.0, 0.0},
                                   {1.0, 0.0, 0.0},
                                   {1.0, 1.0, 0.0},
                                   {0.0, 1.0, 0.0}};

// ---------------------------------------------------------------------------
// 1. The discriminating fixture: rotated, unequal masses.
// ---------------------------------------------------------------------------
void test_near_tie_selection_is_invariant() {
    const auto rotated = rotate(kSquare, 0.37, 0.61, 0.23);
    const std::string expected = "(2,3)";

    const auto expect = [&](const std::string& got, const std::string& what) {
        check(got == expected,
              what + " selected a different physical constraint: expected the "
                     "redundant constraint to be " + expected + ", got " + got +
                  ". Selection must depend only on the physical constraint set, not "
                  "on how it was written down, oriented or scaled");
    };

    expect(redundant_pairs(rotated, kAllSixPairs, {0, 1, 2, 3}, 16.0), "the base fixture");

    auto reversed = kAllSixPairs;
    std::reverse(reversed.begin(), reversed.end());
    expect(redundant_pairs(rotated, reversed, {0, 1, 2, 3}, 16.0), "reversing the constraint order");

    auto rotated_list = kAllSixPairs;
    std::rotate(rotated_list.begin(), rotated_list.begin() + 2, rotated_list.end());
    expect(redundant_pairs(rotated, rotated_list, {0, 1, 2, 3}, 16.0), "rotating the constraint list");

    auto swapped = kAllSixPairs;
    std::swap(swapped[0], swapped[5]);
    expect(redundant_pairs(rotated, swapped, {0, 1, 2, 3}, 16.0),
           "swapping the first and last constraint");

    // Storage permutations: the same tags in different slots.
    expect(redundant_pairs({rotated[2], rotated[0], rotated[3], rotated[1]}, kAllSixPairs,
                           {2, 0, 3, 1}, 16.0),
           "permuting the atom storage order");
    expect(redundant_pairs({rotated[3], rotated[2], rotated[1], rotated[0]}, kAllSixPairs,
                           {3, 2, 1, 0}, 16.0),
           "reversing the atom storage order");

    // Uniform coordinate scaling. This is the sharpest of the set: the rank
    // tolerance and every residual norm scale together, so a RELATIVE rule is
    // invariant while an absolute one is not.
    for (double factor : {0.1, 10.0, 1000.0}) {
        expect(redundant_pairs(scale(rotated, factor), kAllSixPairs, {0, 1, 2, 3}, 16.0),
               "scaling the coordinates by " + std::to_string(factor));
    }

    for (double mass_ratio : {1.0, 1.5, 16.0, 200.0}) {
        expect(redundant_pairs(rotated, kAllSixPairs, {0, 1, 2, 3}, mass_ratio),
               "mass ratio " + std::to_string(mass_ratio));
    }

    // A different generic rotation must reach the same physical answer.
    expect(redundant_pairs(rotate(kSquare, 1.1, 0.2, 0.9), kAllSixPairs, {0, 1, 2, 3}, 16.0),
           "a different generic rotation");
    expect(redundant_pairs(rotate(kSquare, 2.4, 1.7, 0.4), kAllSixPairs, {0, 1, 2, 3}, 16.0),
           "a third generic rotation");
}

// ---------------------------------------------------------------------------
// 2. The axis-aligned fixture is NOT a discriminating test.
// ---------------------------------------------------------------------------
// Pinned deliberately. Its tied norms are bit-identical, so it agrees under any
// tie criterion and proves nothing about the relative window. Keeping the fact
// asserted stops it being mistaken for coverage of that rule.
void test_axis_aligned_tie_does_not_discriminate() {
    const auto axis_aligned = redundant_pairs(kSquare, kAllSixPairs, {0, 1, 2, 3}, 16.0);
    check(axis_aligned == "(2,3)",
          "the axis-aligned square should still name (2,3); got " + axis_aligned);
    // It agrees with the rotated fixture here, which is exactly why it cannot
    // tell the two criteria apart -- see the file comment.
    const auto rotated =
        redundant_pairs(rotate(kSquare, 0.37, 0.61, 0.23), kAllSixPairs, {0, 1, 2, 3}, 16.0);
    check(axis_aligned == rotated,
          "axis-aligned and rotated fixtures must agree under the relative rule");
}

// ---------------------------------------------------------------------------
// 3. Genuinely different candidates must NOT be collapsed into a tie.
// ---------------------------------------------------------------------------
// The negative direction. These coplanar quadrilaterals are irregular, so their
// residual norms differ by far more than sqrt(eps) and the pivot order is
// decided by magnitude alone. The constraint each one names is therefore NOT
// the largest canonical key -- which is what the tie rule would tend to leave
// over -- and it is the same under either tie criterion. If the relative window
// ever grew wide enough to swallow real differences, these would start naming
// the canonical-key answer instead.
void test_distinct_candidates_are_not_treated_as_tied() {
    struct Case {
        const char* name;
        std::vector<Vec3> points;
        const char* expected;
    };
    const std::vector<Case> cases = {
        {"kite",      {{0,0,0},{2.7,0,0},{1.35,0.42,0},{1.35,-3.1,0}},   "(1,2)"},
        {"sliver",    {{0,0,0},{3.4,0,0},{3.3,0.11,0},{0.9,0.07,0}},     "(0,3)"},
        {"tiny-edge", {{0,0,0},{0.31,0,0},{2.9,1.4,0},{1.1,3.3,0}},      "(0,1)"},
        {"skew",      {{0,0,0},{4.3,0.2,0},{1.1,0.9,0},{3.2,2.8,0}},     "(0,2)"},
    };

    int not_largest_key = 0;
    for (const auto& c : cases) {
        for (double mass_ratio : {1.0, 16.0}) {
            const auto got = redundant_pairs(c.points, kAllSixPairs, {0, 1, 2, 3}, mass_ratio);
            check(got == c.expected,
                  std::string("the ") + c.name + " quadrilateral (mass ratio " +
                      std::to_string(mass_ratio) + ") must name " + c.expected +
                      ", decided by residual MAGNITUDE rather than by the canonical key; "
                      "got " + got + ". A tie window wide enough to swallow this "
                      "difference would name the canonical-key answer instead");
        }
        if (std::string(c.expected) != "(2,3)") ++not_largest_key;
    }
    check(not_largest_key == static_cast<int>(cases.size()),
          "every negative fixture must name something other than the largest canonical "
          "key, otherwise it cannot distinguish a magnitude decision from a tie");
}

}  // namespace

int main() {
    test_near_tie_selection_is_invariant();
    test_axis_aligned_tie_does_not_discriminate();
    test_distinct_candidates_are_not_treated_as_tied();

    if (failures == 0) {
        std::cout << "[constraint tie breaking] all checks passed\n";
    } else {
        std::cerr << "[constraint tie breaking] " << failures << " check(s) failed\n";
    }
    return failures == 0 ? 0 : 1;
}
