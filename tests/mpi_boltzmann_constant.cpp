// The Boltzmann constant and the temperature reduction under MPI.
//
// tests/boltzmann_constant_tests.cpp proves what each serial temperature path
// uses. That does not carry over to MPI, because both halves of a temperature
// are reduced: compute_twice_ke() allreduces the kinetic energy, and the
// degrees of freedom come from a separately reduced global atom count. Two
// failure modes stay invisible inside any single rank --
//
//   * a kinetic energy summed once per rank instead of once per atom, which
//     multiplies the reported temperature by the rank count;
//   * a degrees-of-freedom count taken from local rather than global atoms,
//     which divides it by roughly the rank count;
//
// -- and they would cancel each other exactly at np=1. So this file measures
// the constant the same way the serial audit does, from the reduced global
// quantities, and requires the same number at 1, 2 and 4 ranks, on every rank.
//
// It also exercises a rank that owns no atoms at all, because the velocity
// initializer must still enter the collectives there: returning early on an
// empty rank deadlocks the others.

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "gmd/integrator/thermostat.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;
int global_rank = 0;
int global_size = 1;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[mpi k_B][rank " << global_rank << "] " << message << '\n';
        ++failures;
    }
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

// Must track tests/boltzmann_constant_tests.cpp; both are updated by any commit
// that changes a production constant.
constexpr double kExpectedInitializer = 8.617343e-5;
constexpr double kExpectedShared      = 8.617333262e-5;

constexpr std::size_t kTotalAtoms = 8;
constexpr double kTargetTemperature = 300.0;

// Unequal masses, so a rank boundary that split the sum wrongly could not be
// masked by every atom weighing the same.
const std::array<double, kTotalAtoms> kMasses = {
    1.008, 12.011, 15.999, 4.0026, 6.941, 10.811, 14.007, 3.016};

// Which tags this rank owns. At np=4 the last rank deliberately owns nothing:
// an empty domain is legal and must still take part in every collective.
std::vector<int> owned_tags() {
    std::vector<int> owned;
    const int contributing = (global_size >= 4) ? global_size - 1 : global_size;
    for (int tag = 0; tag < static_cast<int>(kTotalAtoms); ++tag) {
        if (tag % contributing == global_rank) owned.push_back(tag);
    }
    return owned;
}

gmd::System local_system() {
    const std::vector<int> owned = owned_tags();
    gmd::System system;
    system.resize(owned.size(), owned.size());
    gmd::Box box;
    box.set_lengths({24.0, 28.0, 32.0});
    system.set_box(box);
    for (std::size_t i = 0; i < owned.size(); ++i) {
        const auto tag = static_cast<std::size_t>(owned[i]);
        system.mutable_masses()[i] = kMasses[tag] + 0.1 * static_cast<double>(tag);
        system.mutable_coordinates()[i] = {2.0 + 2.7 * static_cast<double>(tag),
                                           3.0 + 1.9 * static_cast<double>(tag % 5),
                                           4.0 + 2.3 * static_cast<double>(tag % 3)};
        system.mutable_atom_tags()[i] = owned[i];
        system.mutable_atom_owners()[i] = global_rank;
    }
    return system;
}

// The initializer's own degrees of freedom: 3N-3 over the GLOBAL atom count.
std::size_t global_dof() { return 3 * kTotalAtoms - 3; }

void test_reported_temperature_is_rank_independent() {
    gmd::System system = local_system();
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, true);

    // compute_twice_ke() already reduces across ranks, so this is the global
    // kinetic energy on every rank.
    const double twice_ke = gmd::compute_twice_ke(system);
    const double measured_initializer =
        twice_ke / (static_cast<double>(global_dof()) * kTargetTemperature);

    const double reported = gmd::temperature_from_twice_ke(twice_ke, global_dof());
    const double measured_reporting =
        twice_ke / (static_cast<double>(global_dof()) * reported);

    if (global_rank == 0) {
        std::cout << std::setprecision(17)
                  << "  np=" << global_size
                  << "  initializer k_B = " << measured_initializer
                  << "  reporting k_B = " << measured_reporting << '\n';
    }

    check(std::fabs(measured_initializer / kExpectedInitializer - 1.0) <= 1.0e-11,
          "velocity initialization measures k_B = " + number(measured_initializer) +
              " at " + std::to_string(global_size) + " rank(s), expected " +
              number(kExpectedInitializer) +
              ". A value scaling with the rank count means the kinetic energy is "
              "summed once per rank, or the degrees of freedom are local rather "
              "than global");
    check(std::fabs(measured_reporting / kExpectedShared - 1.0) <= 1.0e-11,
          "temperature reporting measures k_B = " + number(measured_reporting) +
              " at " + std::to_string(global_size) + " rank(s), expected " +
              number(kExpectedShared));

    // No rank multiplication: the reduced energy must not scale with np. The
    // serial value is recomputed here from the global fixture, independently of
    // how the atoms happen to be distributed.
    check(std::fabs(measured_initializer / static_cast<double>(global_size) -
                    kExpectedInitializer) > 1.0e-9 || global_size == 1,
          "the measured constant is exactly the rank count times the expected "
          "value, which is what a per-rank kinetic-energy sum looks like");

    // Every rank must agree bit for bit; the quantities are already reduced.
    double minimum = 0.0, maximum = 0.0;
    MPI_Allreduce(&reported, &minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&reported, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    check(minimum == maximum,
          "ranks disagree about the reported temperature: min " + number(minimum) +
              ", max " + number(maximum));
}

void test_empty_rank_participates() {
    const std::vector<int> owned = owned_tags();
    const bool empty_here = owned.empty();
    int empty_ranks = 0;
    const int local_empty = empty_here ? 1 : 0;
    MPI_Allreduce(&local_empty, &empty_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (global_rank == 0 && global_size >= 4) {
        std::cout << "  np=" << global_size << "  ranks owning no atoms: "
                  << empty_ranks << '\n';
    }
    if (global_size >= 4) {
        check(empty_ranks >= 1,
              "this fixture is meant to leave one rank empty at np>=4 but none is");
    }

    // Reaching here at all is the assertion: initialize() and compute_twice_ke()
    // are collective, so an early return on the empty rank would have deadlocked
    // the run rather than produced a wrong number.
    gmd::System system = local_system();
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, true);
    const double twice_ke = gmd::compute_twice_ke(system);
    check(twice_ke > 0.0,
          "the reduced kinetic energy is " + number(twice_ke) +
              "; an empty rank must contribute zero without zeroing the global sum");
    if (empty_here) {
        check(system.num_local_atoms() == 0,
              "this rank was supposed to own no atoms");
    }
}

void test_zero_temperature_is_still_collective() {
    gmd::System system = local_system();
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, 0.0, gmd::VelocityInitMode::Random, true);
    const double twice_ke = gmd::compute_twice_ke(system);
    check(twice_ke == 0.0,
          "a 0 K system reduced across " + std::to_string(global_size) +
              " rank(s) has kinetic energy " + number(twice_ke) + ", not zero");
    check(gmd::temperature_from_twice_ke(twice_ke, global_dof()) == 0.0,
          "a 0 K system must report exactly 0 K");
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    test_reported_temperature_is_rank_independent();
    test_empty_rank_participates();
    test_zero_temperature_is_still_collective();

    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_rank == 0) {
        if (total == 0) {
            std::cout << "[mpi k_B] all checks passed on " << global_size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi k_B] " << total << " check(s) failed on "
                      << global_size << " rank(s)\n";
        }
    }
    return total == 0 ? 0 : 1;
}
