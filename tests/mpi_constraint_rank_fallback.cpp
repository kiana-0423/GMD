// The rank FALLBACK and the near-tie SELECTION must be collective decisions.
//
// analyze_independence() gathers all owned atoms and works from a constraint
// list held against global atom tags, so every rank builds the same
// mass-weighted Jacobian, computes the same singular values, runs the same
// pivoted selection and reaches the same verdict. Nothing about that is
// per-rank. This test checks it rather than assuming it, for the two decisions
// that are new here:
//
//   * whether the pivoted rank calculation and the singular values DISAGREED,
//     which is the fallback; and
//   * WHICH physical constraint a near-tie selected.
//
// A rank that decided either on its own would show up as a disagreement in the
// allreduce below, and -- because require_independent() throws -- as some ranks
// throwing while others did not, which is the shape of a deadlock.
//
// The fallback fixture is swept rather than pinned: the band of out-of-plane
// heights that makes the two rank methods disagree is narrow, and translating a
// fixture into a domain-decomposed cell changes the round-off in the bond
// vectors. Sweeping and requiring at least one hit asserts the same thing
// without depending on one height surviving that translation.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

#include <mpi.h>

namespace {

constexpr double kBox = 20.0;
using Vec3 = gmd::System::Vec3;

void check(bool condition, const std::string& message, int rank, int& failures) {
    if (!condition) {
        std::cerr << "[mpi constraint rank fallback][rank " << rank << "] " << message << '\n';
        ++failures;
    }
}

gmd::ConstraintSettings settings() {
    gmd::ConstraintSettings s;
    s.tolerance = 1.0e-12;
    s.max_iterations = 500;
    return s;
}

double distance(const std::vector<Vec3>& p, int a, int b) {
    double sum = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        const double delta = p[static_cast<std::size_t>(a)][d] - p[static_cast<std::size_t>(b)][d];
        sum += delta * delta;
    }
    return std::sqrt(sum);
}

// This rank's share of a globally tagged fixture.
gmd::System local_share(const std::vector<Vec3>& positions,
                        const gmd::DomainDecomposition& dd,
                        const gmd::Box& box, int rank, double mass_ratio) {
    std::vector<int> owned;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (dd.owner_rank(box, positions[i]) == rank) owned.push_back(static_cast<int>(i));
    }
    gmd::System system;
    system.resize(owned.size(), owned.size());
    system.set_box(box);
    for (std::size_t k = 0; k < owned.size(); ++k) {
        const auto index = static_cast<std::size_t>(owned[k]);
        system.mutable_masses()[k] = std::pow(mass_ratio, static_cast<double>(owned[k] % 2));
        system.mutable_atom_tags()[k] = owned[k];
        system.mutable_atom_owners()[k] = rank;
        system.mutable_coordinates()[k] = positions[index];
    }
    return system;
}

std::vector<gmd::BondConstraint> hold(const std::vector<Vec3>& positions,
                                      const std::vector<std::array<int, 2>>& pairs) {
    std::vector<gmd::BondConstraint> constraints;
    for (const auto& pair : pairs) {
        constraints.push_back({pair[0], pair[1], distance(positions, pair[0], pair[1])});
    }
    return constraints;
}

const std::vector<std::array<int, 2>> kAllSixPairs = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

// Four near-coplanar atoms, centred so they straddle the domain boundary at
// x = 10 (np=2 and np=4) and y = 10 (np=4).
std::vector<Vec3> straddling_quad(double height) {
    return {{9.43, 9.62, 10.0},
            {10.56, 9.62, 10.0},
            {9.64, 10.69, 10.0},
            {10.22, 10.24, 10.0 + height}};
}

// Everything inside one octant, so at np=2 and np=4 at least one rank owns
// nothing at all and still has to reach the same verdict.
std::vector<Vec3> single_domain_quad(double height) {
    return {{2.0, 2.0, 2.0},
            {3.13, 2.0, 2.0},
            {2.21, 3.07, 2.0},
            {2.79, 2.62, 2.0 + height}};
}

// Encodes one rank's verdict as integers so it can be compared across ranks by
// an allreduce of min and max: equal min and max means every rank agreed.
struct Verdict {
    int deficient = 0;
    int fell_back = 0;
    int rank_value = 0;
    int dependent_count = 0;
    int report_identified = 0;
    int selected_i = -1;     // canonical pair of the FIRST redundant constraint
    int selected_j = -1;
};

Verdict analyse(const std::vector<Vec3>& positions,
                const std::vector<std::array<int, 2>>& pairs,
                const gmd::DomainDecomposition& dd, const gmd::Box& box,
                int rank, double mass_ratio) {
    const gmd::ConstraintSolver solver(hold(positions, pairs), settings());
    const auto system = local_share(positions, dd, box, rank, mass_ratio);
    const auto report = solver.analyze_independence(system);

    Verdict v;
    v.report_identified = report.redundant_identified() ? 1 : 0;
    for (const auto& component : report.components) {
        if (component.independent()) continue;
        v.deficient = 1;
        v.fell_back = component.redundant_constraints_identified ? 0 : 1;
        v.rank_value = static_cast<int>(component.rank);
        v.dependent_count = static_cast<int>(component.redundant_constraint_indices.size());
        if (!component.redundant_constraint_indices.empty()) {
            const auto& c = solver.constraints()[component.redundant_constraint_indices.front()];
            v.selected_i = c.i;
            v.selected_j = c.j;
        }
    }
    return v;
}

// Every field must be identical on every rank.
void require_agreement(const Verdict& v, const std::string& label, int rank, int size,
                       int& failures) {
    const std::array<int, 7> mine = {v.deficient, v.fell_back, v.rank_value,
                                     v.dependent_count, v.report_identified,
                                     v.selected_i, v.selected_j};
    std::array<int, 7> lo{}, hi{};
    MPI_Allreduce(mine.data(), lo.data(), 7, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(mine.data(), hi.data(), 7, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    static const char* names[7] = {"deficient", "fell_back", "rank", "dependent_count",
                                   "report_identified", "selected_i", "selected_j"};
    for (std::size_t k = 0; k < 7; ++k) {
        check(lo[k] == hi[k],
              label + ": ranks disagreed about " + names[k] + " (min " +
                  std::to_string(lo[k]) + ", max " + std::to_string(hi[k]) +
                  ") on " + std::to_string(size) +
                  " rank(s); the rank analysis must be a collective decision",
              rank, failures);
    }
}

int run(int rank, int size) {
    int failures = 0;

    gmd::Box box;
    box.set_lengths({kBox, kBox, kBox});
    gmd::DomainDecomposition dd;
    dd.create_decomposition(box, size, rank);

    // ---- 1. the fallback, swept over the band, straddling a rank boundary ----
    int hits = 0;
    for (int step = 0; step <= 24; ++step) {
        const double height = 2.0e-15 * std::pow(5.0 / 2.0, static_cast<double>(step) / 24.0);
        const auto positions = straddling_quad(height);
        const auto verdict = analyse(positions, kAllSixPairs, dd, box, rank, 1.5);
        require_agreement(verdict, "straddling fallback sweep h=" + std::to_string(height),
                          rank, size, failures);
        if (verdict.fell_back == 1) {
            ++hits;
            check(verdict.dependent_count == 6,
                  "a fallen-back component must name its whole membership, got " +
                      std::to_string(verdict.dependent_count),
                  rank, failures);
            check(verdict.report_identified == 0,
                  "a report containing a fallen-back component must not claim "
                  "identification", rank, failures);
        }
    }
    check(hits > 0,
          "no height in the swept band reached the fallback under MPI; the band is "
          "derived in tests/constraint_rank_fallback_tests.cpp and may need "
          "re-deriving rather than widening", rank, failures);

    // ---- 2. the same, with every atom inside one domain (empty ranks) ----
    int empty_hits = 0;
    for (int step = 0; step <= 24; ++step) {
        const double height = 2.0e-15 * std::pow(5.0 / 2.0, static_cast<double>(step) / 24.0);
        const auto positions = single_domain_quad(height);
        const auto verdict = analyse(positions, kAllSixPairs, dd, box, rank, 1.5);
        require_agreement(verdict, "empty-rank fallback sweep h=" + std::to_string(height),
                          rank, size, failures);
        if (verdict.fell_back == 1) ++empty_hits;
    }
    check(empty_hits > 0,
          "no height reached the fallback with the fixture confined to one domain; a "
          "rank owning no atoms must still reach the same verdict", rank, failures);

    // A rank that owns nothing must still have participated: with size > 1 the
    // single-domain fixture leaves at least one rank empty by construction.
    const auto empty_system = local_share(single_domain_quad(3.311e-15), dd, box, rank, 1.5);
    int local_empty = empty_system.num_local_atoms() == 0 ? 1 : 0;
    int empty_ranks = 0;
    MPI_Allreduce(&local_empty, &empty_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (size > 1) {
        check(empty_ranks > 0,
              "the single-domain fixture must leave at least one rank owning no atoms, "
              "so the empty-rank path is genuinely exercised; found " +
                  std::to_string(empty_ranks), rank, failures);
    }

    // ---- 3. near-tie selection must name the same physical constraint ----
    const std::vector<Vec3> square = {{9.5, 9.5, 10.0},
                                      {10.5, 9.5, 10.0},
                                      {10.5, 10.5, 10.0},
                                      {9.5, 10.5, 10.0}};
    auto rotate = [](std::vector<Vec3> p, double a, double b, double c) {
        for (auto& q : p) {
            const double x = q[0] - 10.0, y = q[1] - 10.0, z = q[2] - 10.0;
            const double x1 = x * std::cos(a) - y * std::sin(a);
            const double y1 = x * std::sin(a) + y * std::cos(a);
            const double y2 = y1 * std::cos(b) - z * std::sin(b);
            const double z2 = y1 * std::sin(b) + z * std::cos(b);
            const double z3 = z2 * std::cos(c) - x1 * std::sin(c);
            const double x3 = z2 * std::sin(c) + x1 * std::cos(c);
            q = {x3 + 10.0, y2 + 10.0, z3 + 10.0};
        }
        return p;
    };
    for (double mass_ratio : {1.0, 16.0}) {
        const auto verdict =
            analyse(rotate(square, 0.37, 0.61, 0.23), kAllSixPairs, dd, box, rank, mass_ratio);
        require_agreement(verdict,
                          "near-tie selection, mass ratio " + std::to_string(mass_ratio),
                          rank, size, failures);
        check(verdict.deficient == 1,
              "the rotated square must be rank deficient (coplanar, 6 constraints at "
              "rank 5)", rank, failures);
    }

    // ---- 4. rejection is collective: every rank must throw ----
    {
        const auto positions = straddling_quad(3.311e-15);
        const gmd::ConstraintSolver solver(hold(positions, kAllSixPairs), settings());
        const auto system = local_share(positions, dd, box, rank, 1.5);
        int threw = 0;
        try {
            solver.require_independent(system, "the MPI test geometry");
        } catch (const std::invalid_argument&) {
            threw = 1;
        }
        int lo = 0, hi = 0;
        MPI_Allreduce(&threw, &lo, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&threw, &hi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        check(lo == hi,
              "some ranks rejected the deficient set and others did not (min " +
                  std::to_string(lo) + ", max " + std::to_string(hi) +
                  "); a non-collective rejection is how this deadlocks in production",
              rank, failures);
        check(lo == 1, "every rank must reject a rank-deficient constraint set",
              rank, failures);
    }

    return failures;
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int failures = run(rank, size);
    int total_failures = 0;
    MPI_Allreduce(&failures, &total_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (total_failures == 0) {
            std::cout << "[mpi constraint rank fallback] all checks passed on " << size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi constraint rank fallback] " << total_failures
                      << " check(s) failed on " << size << " rank(s)\n";
        }
    }
    return total_failures == 0 ? 0 : 1;
}
