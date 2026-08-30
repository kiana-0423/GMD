// The Berendsen barostat's coupling factor under MPI.
//
// validation/berendsen_npt_lj exercises the barostat over a trajectory, but it
// cannot compare rank counts: a run started from `velocity_init random` does
// not reproduce across them, because VelocityInitializer draws from one
// generator sequentially by LOCAL atom index. Each rank therefore hands its own
// first atom the generator's first three draws, which under decomposition is a
// different physical atom than in a serial run, and the rescale to the target
// temperature hides the difference in every reported scalar. That is a
// pre-existing limitation of the initializer and has nothing to do with the
// barostat -- the same fixture run with `velocity 0.0` stays bit-identical
// across ranks for fifty steps.
//
// So this file removes the initializer from the question. Every rank builds the
// same velocity field by construction, from the global atom tag rather than
// from a generator, and the barostat's own arithmetic is what gets compared.
//
// What could go wrong here and nowhere else:
//
//   * BerendsenBarostat::apply() calls global_atom_count() and
//     compute_twice_ke(), both of which reduce. A kinetic energy summed once
//     per rank instead of once per atom would scale the instantaneous pressure
//     with the rank count and change mu.
//   * It scales only locally owned coordinates, so a rank that owns none must
//     still enter every collective inside the reduction; returning early there
//     deadlocks the others. At np=4 one rank deliberately owns nothing.
//   * The resulting box must be identical on every rank, not merely on rank 0,
//     because every rank integrates against it on the next step.

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <mpi.h>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/berendsen_barostat.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;
int global_rank = 0;
int global_size = 1;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[mpi berendsen][rank " << global_rank << "] " << message << '\n';
        ++failures;
    }
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

class ZeroForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "zero_force"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request, gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.forces.assign(request.coordinates.size(), {0.0, 0.0, 0.0});
        result.potential_energy = 0.0;
        result.virial_valid = false;
        result.success = true;
    }
};

constexpr std::size_t kTotalAtoms = 8;
constexpr double kBoxLength = 20.0;
constexpr double kBeta = 4.5e-5;
constexpr double kTau = 500.0;
constexpr double kTimeStep = 1.0;
constexpr double kTargetPressureBar = 3000.0;

// Unequal masses, so a mis-split sum could not be masked by uniformity.
const std::array<double, kTotalAtoms> kMasses = {
    1.008, 12.011, 15.999, 4.0026, 6.941, 10.811, 14.007, 3.016};

// Velocity as a pure function of the GLOBAL tag: the whole point is that no
// generator, and therefore no atom ordering, enters.
std::array<double, 3> velocity_for(std::size_t tag) {
    const double t = static_cast<double>(tag);
    return {0.004 + 0.0013 * t, -0.0021 + 0.0005 * t, 0.0032 - 0.0003 * t};
}

std::array<double, 9> global_virial() {
    std::array<double, 9> virial{};
    virial[0] = -1.9;
    virial[4] = 1.1;
    virial[8] = 2.7;
    virial[1] = 0.5;
    virial[3] = 0.5;
    return virial;
}

// At np>=4 the last rank deliberately owns nothing.
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
    box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    system.set_box(box);
    for (std::size_t i = 0; i < owned.size(); ++i) {
        const auto tag = static_cast<std::size_t>(owned[i]);
        system.mutable_masses()[i] = kMasses[tag];
        system.mutable_velocities()[i] = velocity_for(tag);
        system.mutable_coordinates()[i] = {2.0 + 2.1 * static_cast<double>(tag),
                                           3.0 + 1.3 * static_cast<double>(tag % 5),
                                           4.0 + 1.7 * static_cast<double>(tag % 3)};
        system.mutable_atom_tags()[i] = owned[i];
        system.mutable_atom_owners()[i] = global_rank;
    }
    system.set_last_virial(global_virial(), true);
    return system;
}

// The whole-system kinetic energy, computed from the global fixture and
// independent of how the atoms happen to be distributed.
double serial_twice_kinetic_energy() {
    double total = 0.0;
    for (std::size_t tag = 0; tag < kTotalAtoms; ++tag) {
        const auto v = velocity_for(tag);
        total += kMasses[tag] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    return total;
}

// Applies the barostat once and returns the resulting box length.
double coupled_box_length() {
    gmd::System system = local_system();
    ZeroForceProvider provider;
    gmd::RuntimeContext runtime;
    gmd::BerendsenBarostat barostat(kTau, kBeta);
    const auto virial = global_virial();
    barostat.apply(system, provider, runtime, 0, kTimeStep, 300.0,
                   kTargetPressureBar, virial[0] + virial[4] + virial[8]);
    return system.box().lengths[0];
}

void test_coupling_is_rank_independent() {
    const double length = coupled_box_length();

    // The expected value, computed here from the global fixture using the
    // reference conversions rather than any production symbol:
    //   mu^3 = 1 - beta * (dt/tau) * (P_target - P_current[bar])
    const auto virial = global_virial();
    const double trace = virial[0] + virial[4] + virial[8];
    const double volume = kBoxLength * kBoxLength * kBoxLength;
    const double pressure_internal =
        (serial_twice_kinetic_energy() + trace) / (3.0 * volume);
    // 1 eV/A^3 = 801088317/500 bar exactly; see tests/pressure_unit_tests.cpp.
    const double pressure_bar = pressure_internal * (801088317.0 / 500.0);
    const double mu_cubed =
        1.0 - kBeta * (kTimeStep / kTau) * (kTargetPressureBar - pressure_bar);
    const double expected = kBoxLength * std::cbrt(mu_cubed);

    if (global_rank == 0) {
        std::cout << std::setprecision(17) << "  np=" << global_size
                  << "  box " << length << " A (expected " << expected << ")\n";
    }

    check(std::fabs(length / expected - 1.0) <= 1.0e-13,
          "the coupled box length is " + number(length) + " A at " +
              std::to_string(global_size) + " rank(s), expected " +
              number(expected) +
              ". A value that moves with the rank count means the kinetic energy "
              "or the atom count behind the instantaneous pressure is summed "
              "once per rank rather than once per atom");

    // Every rank must hold the same cell, not just rank 0.
    double from_rank_zero = length;
    MPI_Bcast(&from_rank_zero, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    check(length == from_rank_zero,
          "rank " + std::to_string(global_rank) + " ends with box " +
              number(length) + " where rank 0 has " + number(from_rank_zero));
}

void test_empty_rank_participates() {
    // At np>=4 the last rank owns nothing. It must still enter the reductions
    // inside apply() and come out with the same cell as everyone else.
    gmd::System system = local_system();
    if (global_size >= 4 && global_rank == global_size - 1) {
        check(system.num_local_atoms() == 0,
              "this fixture expects the last rank to own no atoms at np>=4");
    }
    const double length = coupled_box_length();
    double from_rank_zero = length;
    MPI_Bcast(&from_rank_zero, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    check(length == from_rank_zero,
          "a rank owning " + std::to_string(system.num_local_atoms()) +
              " atom(s) produced box " + number(length) + " instead of " +
              number(from_rank_zero));
}

void test_direction_is_rank_independent() {
    // The sign of the response is the property validation/berendsen_npt_lj
    // asserts over a trajectory. It must not depend on the decomposition either.
    gmd::System low = local_system();
    gmd::System high = local_system();
    ZeroForceProvider provider;
    gmd::RuntimeContext runtime;
    const auto virial = global_virial();
    const double trace = virial[0] + virial[4] + virial[8];

    gmd::BerendsenBarostat expand(kTau, kBeta);
    expand.apply(low, provider, runtime, 0, kTimeStep, 300.0, -50000.0, trace);
    gmd::BerendsenBarostat compress(kTau, kBeta);
    compress.apply(high, provider, runtime, 0, kTimeStep, 300.0, 50000.0, trace);

    check(low.box().lengths[0] > kBoxLength,
          "a target far below the instantaneous pressure must expand the cell; "
          "at " + std::to_string(global_size) + " rank(s) it gave " +
              number(low.box().lengths[0]));
    check(high.box().lengths[0] < kBoxLength,
          "a target far above the instantaneous pressure must compress the cell; "
          "at " + std::to_string(global_size) + " rank(s) it gave " +
              number(high.box().lengths[0]));
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    test_coupling_is_rank_independent();
    test_empty_rank_participates();
    test_direction_is_rank_independent();

    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_rank == 0) {
        if (total == 0) {
            std::cout << "[mpi berendsen] all checks passed on " << global_size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi berendsen] " << total << " check(s) failed on "
                      << global_size << " rank(s)\n";
        }
    }
    return total == 0 ? 0 : 1;
}
