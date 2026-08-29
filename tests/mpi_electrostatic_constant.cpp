// The electrostatic constant under MPI domain decomposition.
//
// tests/electrostatic_constant_tests.cpp proves that every serial
// electrostatic path carries exactly one factor of k_e. That proof does not
// carry over to MPI on its own, because the providers allreduce their own
// energies and virials: a constant applied once per rank instead of once per
// pair, or a self/background term that every rank contributes in full, would
// show up as a k_e that scales with the rank count while staying perfectly
// self-consistent inside any single rank.
//
// So this file measures the constant the same way -- provider result divided
// by a k_e = 1 reference for the whole global system -- and requires the
// answer to be the same number at 1, 2 and 4 ranks, and the same number every
// rank reports. The reference is the global sum, computed identically on every
// rank from the full configuration, which is what makes rank-count dependence
// visible rather than cancelling out.

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;
int global_rank = 0;
int global_size = 1;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[mpi k_e][rank " << global_rank << "] " << message << '\n';
        ++failures;
    }
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

// Must track tests/electrostatic_constant_tests.cpp; both are updated by any
// commit that changes the production constant.
constexpr double kExpectedProduction = 14.3996454686836;
constexpr double kAlpha = 0.32;
constexpr double kCutoff = 8.0;
constexpr int kKmax = 10;

// Eight charges, exactly neutral, irregular positions in a non-cubic box, so
// that a decomposition splitting them differently at 2 and 4 ranks still has
// pairs straddling every boundary.
//
// Twenty of the twenty-eight pairs are inside the 8 A cutoff, minimum
// separation 2.94 A, so the real-space erfc term carries real weight. The
// first version of this fixture was spread across the whole box with EVERY
// pair outside the cutoff, which meant the real-space sum contributed nothing
// and a defect confined to it could not be seen here at all. That is asserted
// rather than assumed, in test_real_space_actually_contributes().
const std::array<std::array<double, 3>, 8> kPositions = {{
    {2.13, 3.41, 4.77},  {4.67, 4.20, 6.02},  {7.31, 2.86, 3.95},
    {3.95, 8.12, 7.44},  {9.84, 6.53, 9.17},  {6.28, 10.35, 4.61},
    {11.42, 9.78, 6.30}, {8.06, 5.17, 11.83}}};
// Multiples of 1/8, so the sum is exactly zero in binary floating point.
const std::array<double, 8> kCharges = {0.75, -0.5, 0.625, -0.5,
                                        -0.375, 0.25, -0.625, 0.375};
const std::array<double, 3> kLengths = {17.0, 21.0, 25.0};

// The same k_e = 1 Ewald sum as the serial audit, over the whole global system.
double reference_energy_unit_constant(int kmax, bool include_real_space = true) {
    const std::size_t n = kCharges.size();
    long double energy = 0.0L;
    const long double alpha = kAlpha;
    const long double two_alpha_over_root_pi =
        2.0L * alpha / std::sqrt(std::numbers::pi_v<long double>);
    (void)two_alpha_over_root_pi;

    for (std::size_t i = 0; i + 1 < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            long double dr[3];
            for (std::size_t d = 0; d < 3; ++d) {
                dr[d] = kPositions[i][d] - kPositions[j][d];
                dr[d] -= kLengths[d] * std::round(dr[d] / kLengths[d]);
            }
            const long double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if (!include_real_space) continue;
            if (r2 >= static_cast<long double>(kCutoff) * kCutoff) continue;
            const long double r = std::sqrt(r2);
            energy += kCharges[i] * kCharges[j] * std::erfc(alpha * r) / r;
        }
    }

    const long double volume = kLengths[0] * kLengths[1] * kLengths[2];
    const long double inv_four_alpha_sq = 1.0L / (4.0L * alpha * alpha);
    const long double two_pi = 2.0L * std::numbers::pi_v<long double>;
    for (int nx = -kmax; nx <= kmax; ++nx) {
      for (int ny = -kmax; ny <= kmax; ++ny) {
        for (int nz = -kmax; nz <= kmax; ++nz) {
          if (nx == 0 && ny == 0 && nz == 0) continue;
          const long double k[3] = {two_pi * nx / kLengths[0],
                                    two_pi * ny / kLengths[1],
                                    two_pi * nz / kLengths[2]};
          const long double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
          const long double g = 4.0L * std::numbers::pi_v<long double>
                              / (volume * k2) * std::exp(-k2 * inv_four_alpha_sq);
          long double sre = 0.0L, sim = 0.0L;
          for (std::size_t i = 0; i < n; ++i) {
              const long double phase = k[0]*kPositions[i][0] + k[1]*kPositions[i][1]
                                      + k[2]*kPositions[i][2];
              sre += kCharges[i] * std::cos(phase);
              sim += kCharges[i] * std::sin(phase);
          }
          energy += 0.5L * g * (sre * sre + sim * sim);
        }
      }
    }

    long double q2 = 0.0L;
    for (double q : kCharges) q2 += static_cast<long double>(q) * q;
    energy -= alpha / std::sqrt(std::numbers::pi_v<long double>) * q2;
    // The fixture is exactly neutral, so there is no background term.
    return static_cast<double>(energy);
}

double reference_energy_unit_constant_no_real_space(int kmax) {
    return reference_energy_unit_constant(kmax, /*include_real_space=*/false);
}

// Contiguous blocks of tags per rank: at 2 ranks the split falls between tags
// 3 and 4, at 4 ranks between 1/2, 3/4 and 5/6, so different pairs straddle a
// boundary in each configuration.
// A local domain holds only its own atoms, so the provider must be given the
// full configuration to evaluate against; the ownership fields are what stop
// each pair being counted more than once. This mirrors what the production
// path does after a ghost exchange.
gmd::System replicated_system_with_local_ownership() {
    gmd::System system;
    gmd::Box box;
    box.set_lengths({kLengths[0], kLengths[1], kLengths[2]});
    // Local atoms first, then the rest as ghosts, which is the layout
    // num_local_atoms() describes.
    std::vector<int> order;
    for (int tag = 0; tag < 8; ++tag)
        if (tag % global_size == global_rank) order.push_back(tag);
    for (int tag = 0; tag < 8; ++tag)
        if (tag % global_size != global_rank) order.push_back(tag);

    std::size_t local_count = 0;
    for (int tag = 0; tag < 8; ++tag)
        if (tag % global_size == global_rank) ++local_count;

    system.resize(8, local_count);
    system.set_box(box);
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto tag = static_cast<std::size_t>(order[i]);
        system.mutable_coordinates()[i] = kPositions[tag];
        system.mutable_charges()[i] = kCharges[tag];
        system.mutable_masses()[i] = 1.0;
        system.mutable_atom_tags()[i] = order[i];
        system.mutable_atom_owners()[i] =
            (order[i] % global_size == global_rank) ? global_rank
                                                    : (order[i] % global_size);
    }
    return system;
}

// The providers return a RANK-LOCAL PARTIAL energy: the reciprocal, self and
// net-charge terms are each divided by the rank count inside the provider so
// that the allreduce Simulation performs afterwards reconstitutes them exactly
// once. Measuring the constant therefore means reproducing that reduction. It
// is also what makes this test sharp: drop the division by mpi_size() and the
// reduced energy is rank-count times too large, which shows up here directly
// as a k_e that scales with np.
double measure(gmd::ForceProvider& provider, const gmd::System& system) {
    gmd::RuntimeContext runtime;
    const auto coordinates = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .coordinates = std::span<const gmd::Coordinate3D>(coordinates.data(),
                                                          coordinates.size()),
        .neighbor_list = nullptr,
    };
    gmd::ForceResult result;
    provider.initialize(runtime);
    provider.compute(request, result, runtime);
    provider.finalize(runtime);

    double global_energy = 0.0;
    MPI_Allreduce(&result.potential_energy, &global_energy, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    return global_energy;
}

// A term that contributes nothing cannot be checked by dividing it out, so
// the fixture has to actually exercise the real-space sum.
void test_real_space_actually_contributes() {
    const double full = reference_energy_unit_constant(kKmax);
    // The reciprocal, self and background terms alone: same call with every
    // pair pushed outside the cutoff.
    const double reciprocal_only = reference_energy_unit_constant_no_real_space(kKmax);
    const double share = std::fabs((full - reciprocal_only) / full);
    if (global_rank == 0) {
        std::cout << "  real-space share of the fixture's energy: "
                  << (100.0 * share) << " %\n";
    }
    check(share > 0.01,
          "the real-space term is only " + number(100.0 * share) +
              "% of this fixture's energy, so a defect confined to it would not "
              "be measurable here");
}

void test_constant_is_independent_of_rank_count() {
    const gmd::System system = replicated_system_with_local_ownership();
    const double reference = reference_energy_unit_constant(kKmax);

    gmd::EwaldForceProvider ewald(kAlpha, kKmax, kCutoff);
    const double ewald_energy = measure(ewald, system);
    const double ewald_ke = ewald_energy / reference;

    std::cout << std::setprecision(17);
    if (global_rank == 0) {
        std::cout << "  np=" << global_size << "  Ewald k_e = " << ewald_ke << '\n';
    }

    check(std::fabs(ewald_ke / kExpectedProduction - 1.0) <= 1.0e-11,
          "Ewald measures k_e = " + number(ewald_ke) + " at " +
              std::to_string(global_size) + " rank(s), expected " +
              number(kExpectedProduction) +
              ". A value scaling with the rank count means a term is applied "
              "once per rank instead of once per pair");

    // Every rank must agree bit for bit: the providers allreduce, so a rank
    // that disagreed would mean the reduction is not covering every term.
    double minimum = 0.0, maximum = 0.0;
    MPI_Allreduce(&ewald_ke, &minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&ewald_ke, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    check(minimum == maximum,
          "ranks disagree about k_e: min " + number(minimum) + ", max " +
              number(maximum));
}

void test_pme_matches_ewald_constant() {
    const gmd::System system = replicated_system_with_local_ownership();

    gmd::EwaldForceProvider ewald(kAlpha, 16, kCutoff);
    const double ewald_energy = measure(ewald, system);

    gmd::PMEForceProvider pme(kAlpha, kCutoff, 6, {64, 64, 64});
    const double pme_energy = measure(pme, system);

    // At order 6 on a 64^3 mesh the residual is mesh error, orders of
    // magnitude below the 3.16e-06 change this audit is about; what would fail
    // here is PME carrying a different constant from Ewald.
    const double relative = std::fabs(pme_energy / ewald_energy - 1.0);
    if (global_rank == 0) {
        std::cout << "  np=" << global_size << "  PME/Ewald - 1 = " << relative << '\n';
    }
    check(relative <= 1.0e-6,
          "PME and Ewald energies differ by " + number(relative) +
              " relative at " + std::to_string(global_size) +
              " rank(s); far beyond mesh error, so the two paths are not "
              "using the same Coulomb constant");
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    test_real_space_actually_contributes();
    test_constant_is_independent_of_rank_count();
    test_pme_matches_ewald_constant();

    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_rank == 0) {
        if (total == 0) {
            std::cout << "[mpi k_e] all checks passed on " << global_size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi k_e] " << total << " check(s) failed on "
                      << global_size << " rank(s)\n";
        }
    }
    return total == 0 ? 0 : 1;
}
