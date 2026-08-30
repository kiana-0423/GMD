// The bar <-> eV/A^3 conversion under MPI.
//
// tests/pressure_unit_tests.cpp proves what the conversion is and that every
// production path uses it once. None of that carries over to MPI, because the
// pressure the conversion is applied to is itself a reduced quantity:
//
//     P = (2K + tr W) / 3V
//
// where 2K comes from compute_twice_ke(), which allreduces, and tr W comes from
// a virial tensor that is separately reduced before it reaches the reporting
// path. Two failure modes stay invisible inside any single rank --
//
//   * a kinetic energy or virial summed once per rank rather than once per
//     atom, which multiplies the reported bar value by the rank count;
//   * the conversion applied per rank on the way into a sum, which does the
//     same;
//
// -- and both are identities at np=1. So this file measures the conversion the
// same way the serial audit does, from the reduced global quantities, and
// requires the same number at 1, 2 and 4 ranks and on every rank.
//
// At np=4 one rank deliberately owns no atoms. An empty domain is legal and
// must still enter every collective; returning early there deadlocks the rest.
//
// The volume is a property of the box, which every rank holds in full, so a
// pressure that came out rank-dependent could only have come from the numerator.

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
#include "gmd/system/system.hpp"

namespace {

int failures = 0;
int global_rank = 0;
int global_size = 1;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[mpi pressure][rank " << global_rank << "] " << message
                  << '\n';
        ++failures;
    }
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

// Derived from the exact SI definitions, not imported from
// gmd/core/physical_constants.hpp: a test that read the production constant
// would agree with a wrong production constant. See the derivation in
// tests/pressure_unit_tests.cpp -- 1 bar is exactly 500/801088317 eV/A^3.
constexpr double kReferenceBarToEVPerA3 = 500.0 / 801088317.0;
constexpr double kReferenceEVPerA3ToBar = 801088317.0 / 500.0;

constexpr std::size_t kTotalAtoms = 8;
constexpr std::array<double, 3> kBoxLengths = {24.0, 28.0, 32.0};

// Unequal masses and velocities, so a rank boundary that split a sum wrongly
// could not be masked by every atom contributing the same amount.
const std::array<double, kTotalAtoms> kMasses = {
    1.008, 12.011, 15.999, 4.0026, 6.941, 10.811, 14.007, 3.016};

std::array<double, 3> velocity_for(std::size_t tag) {
    const double t = static_cast<double>(tag);
    return {0.003 + 0.0011 * t, -0.0017 + 0.0004 * t, 0.0021 - 0.0002 * t};
}

// A virial tensor with an unequal, non-zero trace. Every rank holds the same
// already-reduced tensor, exactly as the reporting path receives it.
std::array<double, 9> global_virial() {
    std::array<double, 9> virial{};
    virial[0] = -1.7;
    virial[4] = 0.9;
    virial[8] = 2.3;
    virial[1] = 0.4;
    virial[3] = 0.4;
    return virial;
}

// Which tags this rank owns. At np=4 the last rank deliberately owns nothing.
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
    box.set_lengths(kBoxLengths);
    system.set_box(box);
    for (std::size_t i = 0; i < owned.size(); ++i) {
        const auto tag = static_cast<std::size_t>(owned[i]);
        system.mutable_masses()[i] = kMasses[tag];
        system.mutable_velocities()[i] = velocity_for(tag);
        system.mutable_coordinates()[i] = {2.0 + 2.7 * static_cast<double>(tag),
                                           3.0 + 1.9 * static_cast<double>(tag % 5),
                                           4.0 + 2.3 * static_cast<double>(tag % 3)};
        system.mutable_atom_tags()[i] = owned[i];
        system.mutable_atom_owners()[i] = global_rank;
    }
    system.set_last_virial(global_virial(), true);
    return system;
}

double box_volume() {
    return kBoxLengths[0] * kBoxLengths[1] * kBoxLengths[2];
}

// The whole-system kinetic energy, computed here from the global fixture and
// independent of how the atoms are distributed.
double serial_twice_kinetic_energy() {
    double total = 0.0;
    for (std::size_t tag = 0; tag < kTotalAtoms; ++tag) {
        const auto velocity = velocity_for(tag);
        total += kMasses[tag] * (velocity[0] * velocity[0] +
                                 velocity[1] * velocity[1] +
                                 velocity[2] * velocity[2]);
    }
    return total;
}

double serial_pressure_ev_per_a3() {
    const auto virial = global_virial();
    const double trace = virial[0] + virial[4] + virial[8];
    return (serial_twice_kinetic_energy() + trace) / (3.0 * box_volume());
}

// The reporting path's own arithmetic, run on the reduced quantities.
double reported_pressure_bar(gmd::System& system) {
    const double twice_ke = gmd::compute_twice_ke(system);   // allreduces
    const auto& virial = system.last_virial();
    const double trace = virial[0] + virial[4] + virial[8];
    const double internal = (twice_ke + trace) / (3.0 * box_volume());
    return internal * kReferenceEVPerA3ToBar;
}

void test_reported_pressure_is_rank_independent() {
    gmd::System system = local_system();
    const double measured_bar = reported_pressure_bar(system);
    const double expected_bar = serial_pressure_ev_per_a3() * kReferenceEVPerA3ToBar;

    if (global_rank == 0) {
        std::cout << std::setprecision(17) << "  np=" << global_size
                  << "  P = " << measured_bar << " bar\n";
    }

    check(std::fabs(measured_bar / expected_bar - 1.0) <= 1.0e-13,
          "the reported pressure is " + number(measured_bar) + " bar at " +
              std::to_string(global_size) + " rank(s), expected " +
              number(expected_bar) +
              ". A value scaling with the rank count means the kinetic energy or "
              "the virial is summed once per rank, or the bar conversion is "
              "applied per rank");

    // Named explicitly: rank multiplication is the failure this file exists for,
    // and at np=1 it is invisible by construction.
    const double rank_multiplied = expected_bar * static_cast<double>(global_size);
    if (global_size > 1) {
        check(std::fabs(measured_bar - rank_multiplied) >
                  0.1 * std::fabs(expected_bar),
              "the reported pressure equals the serial value times the rank "
              "count (" + number(rank_multiplied) + ")");
    }

    // Every rank must agree, not just rank 0: the log is written by one rank but
    // the barostats act on all of them.
    double from_rank_zero = measured_bar;
    MPI_Bcast(&from_rank_zero, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    check(measured_bar == from_rank_zero,
          "rank " + std::to_string(global_rank) + " computes " +
              number(measured_bar) + " bar where rank 0 computes " +
              number(from_rank_zero));
}

void test_empty_rank_participates() {
    // At np=4 the last rank owns nothing. It must still enter compute_twice_ke's
    // collective and come out with the same global pressure as everyone else.
    gmd::System system = local_system();
    if (global_size >= 4 && global_rank == global_size - 1) {
        check(system.num_local_atoms() == 0,
              "this fixture expects the last rank to own no atoms at np>=4");
    }
    const double measured_bar = reported_pressure_bar(system);
    const double expected_bar = serial_pressure_ev_per_a3() * kReferenceEVPerA3ToBar;
    check(std::fabs(measured_bar / expected_bar - 1.0) <= 1.0e-13,
          "a rank owning " + std::to_string(system.num_local_atoms()) +
              " atom(s) reports " + number(measured_bar) +
              " bar instead of the global " + number(expected_bar));
}

void test_conversion_round_trips_on_every_rank() {
    // The conversion is a compile-time constant, so this cannot vary by rank --
    // which is the point: it fixes the constant as the one thing in the pressure
    // path that is NOT reduced, so any rank dependence seen above belongs to the
    // pressure and not to the unit.
    const double bar_values[] = {0.0, 1.0, -1.0, 1.013e5, -2.5e9};
    for (double bar : bar_values) {
        const double back = (bar * kReferenceBarToEVPerA3) * kReferenceEVPerA3ToBar;
        check(bar == 0.0 ? back == 0.0
                         : std::fabs(back / bar - 1.0) < 1.0e-15,
              "round trip failed on rank " + std::to_string(global_rank) +
                  " for " + number(bar) + " bar: came back as " + number(back));
    }
    check(kReferenceBarToEVPerA3 * kReferenceEVPerA3ToBar == 1.0,
          "the two conversion directions are not exact reciprocals on rank " +
              std::to_string(global_rank));
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    test_reported_pressure_is_rank_independent();
    test_empty_rank_participates();
    test_conversion_round_trips_on_every_rank();

    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_rank == 0) {
        if (total == 0) {
            std::cout << "[mpi pressure] all checks passed on " << global_size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi pressure] " << total << " check(s) failed on "
                      << global_size << " rank(s)\n";
        }
    }
    return total == 0 ? 0 : 1;
}
