// Constraint independence must reach the same verdict at every rank count.
//
// The analysis gathers all owned atoms and works from a constraint list
// replicated against global atom tags, so every rank builds the same
// mass-weighted Jacobian and computes the same rank. That is the claim; this
// test checks it rather than assuming it, for an independent set, a redundant
// closed network and a geometrically degenerate one, at 1, 2 and 4 ranks.
//
// The fixtures straddle a domain boundary so the constraints genuinely span
// owners: a rank-local analysis would see a different graph on each rank.

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
        std::cerr << "[mpi constraint independence][rank " << rank << "] " << message << '\n';
        ++failures;
    }
}

gmd::ConstraintSettings settings() {
    gmd::ConstraintSettings s;
    s.tolerance = 1.0e-10;
    s.max_iterations = 500;
    return s;
}

// This rank's share of a globally tagged fixture centred on x = 10, which is a
// domain boundary for both 2 and 4 ranks.
gmd::System local_share(const std::vector<Vec3>& positions,
                        const gmd::DomainDecomposition& dd,
                        const gmd::Box& box, int rank) {
    std::vector<int> owned;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (dd.owner_rank(box, positions[i]) == rank) owned.push_back(static_cast<int>(i));
    }
    gmd::System system;
    system.resize(owned.size(), owned.size());
    system.set_box(box);
    for (std::size_t k = 0; k < owned.size(); ++k) {
        const auto index = static_cast<std::size_t>(owned[k]);
        system.mutable_masses()[k] = 1.0 + 0.5 * static_cast<double>(index);
        system.mutable_atom_tags()[k] = owned[k];
        system.mutable_atom_owners()[k] = rank;
        system.mutable_coordinates()[k] = positions[index];
    }
    return system;
}

double distance(const std::vector<Vec3>& positions, int a, int b) {
    double sum = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        const double delta = positions[static_cast<std::size_t>(a)][d] -
                             positions[static_cast<std::size_t>(b)][d];
        sum += delta * delta;
    }
    return std::sqrt(sum);
}

std::vector<gmd::BondConstraint> hold(const std::vector<Vec3>& positions,
                                      const std::vector<std::array<int, 2>>& pairs) {
    std::vector<gmd::BondConstraint> constraints;
    for (const auto& pair : pairs) {
        constraints.push_back({pair[0], pair[1], distance(positions, pair[0], pair[1])});
    }
    return constraints;
}

// The verdict, the total rank and the per-component ranks must be identical on
// every rank, and must match what a serial run would conclude.
void check_case(const std::string& label,
                const std::vector<Vec3>& positions,
                const std::vector<std::array<int, 2>>& pairs,
                bool expect_independent,
                std::size_t expected_rank,
                int rank, int size, int& failures) {
    gmd::Box box;
    box.set_lengths({kBox, kBox, kBox});
    gmd::DomainDecomposition dd;
    dd.create_decomposition(box, size, rank, 4.0, 1.0, {true, true, true});
    gmd::System system = local_share(positions, dd, box, rank);

    // The premise: with more than one rank the fixture must actually be split.
    int local_atoms = static_cast<int>(system.num_local_atoms());
    int total = 0;
    int most = 0;
    MPI_Allreduce(&local_atoms, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_atoms, &most, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    check(total == static_cast<int>(positions.size()),
          label + ": the fixture must have " + std::to_string(positions.size()) +
              " atoms in total, got " + std::to_string(total), rank, failures);
    if (size > 1) {
        check(most < static_cast<int>(positions.size()),
              label + ": the constrained atoms must be split across ranks, but one rank "
                      "holds all " + std::to_string(most), rank, failures);
    }

    gmd::ConstraintSolver solver(hold(positions, pairs), settings());
    const auto report = solver.analyze_independence(system);

    check(report.independent == expect_independent,
          label + ": independence verdict is " + std::to_string(report.independent) +
              " on " + std::to_string(size) + " rank(s), expected " +
              std::to_string(expect_independent), rank, failures);
    check(report.rank == expected_rank,
          label + ": rank is " + std::to_string(report.rank) + ", expected " +
              std::to_string(expected_rank), rank, failures);

    // Every rank must agree, exactly.
    long long local_rank_value = static_cast<long long>(report.rank);
    long long lo = 0;
    long long hi = 0;
    MPI_Allreduce(&local_rank_value, &lo, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_rank_value, &hi, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    check(lo == hi, label + ": ranks disagree on the Jacobian rank (min " +
                        std::to_string(lo) + ", max " + std::to_string(hi) + ")",
          rank, failures);

    int local_verdict = report.independent ? 1 : 0;
    int verdict_lo = 0;
    int verdict_hi = 0;
    MPI_Allreduce(&local_verdict, &verdict_lo, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_verdict, &verdict_hi, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    check(verdict_lo == verdict_hi,
          label + ": ranks disagree on the independence verdict", rank, failures);

    // Component decomposition must be identical too: a rank-local view of the
    // constraint graph would split it differently.
    long long components = static_cast<long long>(report.components.size());
    long long clo = 0;
    long long chi = 0;
    MPI_Allreduce(&components, &clo, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&components, &chi, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    check(clo == chi, label + ": ranks disagree on the number of connected components",
          rank, failures);

    // The singular values decide the rank, so they must agree bitwise as well.
    for (const auto& component : report.components) {
        double value = component.smallest_singular_value;
        double slo = 0.0;
        double shi = 0.0;
        MPI_Allreduce(&value, &slo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&value, &shi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        check(slo == shi,
              label + ": ranks disagree on a smallest singular value (min " +
                  std::to_string(slo) + ", max " + std::to_string(shi) + ")",
              rank, failures);
    }
}

int run(int rank, int size) {
    int failures = 0;

    // Independent: a four-atom chain straddling x = 10.
    check_case("independent chain",
               {{9.2, 10.0, 10.0}, {10.4, 10.3, 10.2}, {11.5, 11.2, 10.8}, {12.7, 11.4, 11.7}},
               {{0, 1}, {1, 2}, {2, 3}}, true, 3, rank, size, failures);

    // Redundant: every pair among five atoms, ten constraints over a body with
    // at most nine internal degrees of freedom.
    {
        const std::vector<Vec3> positions{{9.3, 10.0, 10.0}, {10.5, 10.1, 10.2},
                                          {9.8, 11.3, 10.3}, {9.6, 10.2, 11.4},
                                          {10.7, 11.2, 11.3}};
        std::vector<std::array<int, 2>> pairs;
        for (int i = 0; i < 5; ++i) {
            for (int j = i + 1; j < 5; ++j) pairs.push_back({i, j});
        }
        check_case("all-pairs cage over five atoms", positions, pairs, false, 9,
                   rank, size, failures);
    }

    // Geometrically degenerate, and with no constrained triple, so only the rank
    // of the Jacobian can see it: four collinear atoms, the chain plus its
    // closing bond. On a straight line all four gradients point the same way and
    // span three dimensions, not four.
    check_case("four collinear atoms, chain plus closing bond",
               {{8.9, 10.0, 10.0}, {10.2, 10.0, 10.0},
                {11.6, 10.0, 10.0}, {12.9, 10.0, 10.0}},
               {{0, 1}, {1, 2}, {2, 3}, {0, 3}}, false, 3, rank, size, failures);

    // Two molecules sharing no atom, on opposite sides of the boundary.
    check_case("two disjoint dimers",
               {{9.0, 10.0, 10.0}, {10.4, 10.0, 10.0},
                {15.0, 15.0, 15.0}, {16.4, 15.2, 15.1}},
               {{0, 1}, {2, 3}}, true, 2, rank, size, failures);

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
            std::cout << "[mpi constraint independence] all checks passed on " << size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi constraint independence] " << total_failures
                      << " check(s) failed on " << size << " rank(s)\n";
        }
    }
    return total_failures == 0 ? 0 : 1;
}
