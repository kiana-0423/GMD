// Emits the velocity field produced by random initialization, keyed by atom
// tag, so that runs at different rank counts can be compared.
//
// The comparison cannot be done inside a single MPI run. VelocityInitializer's
// centre-of-mass removal and temperature rescale are collective, so a rank
// cannot quietly initialize a second, whole-system copy on the side to compare
// against: it would either deadlock or fold the extra copy into everyone's
// reductions. So this binary produces one artifact per invocation and
// tests/velocity_init_equivalence.py runs it at np=1, 2 and 4 and compares.
//
// The same reasoning is why the field is written keyed by TAG. Comparing by
// array index would compare atom i of one decomposition against a different
// physical atom i of another, which is exactly the confusion the defect lives
// in.
//
// Storage arrangements. --reverse-storage reverses each rank's local order
// while leaving the physical system, the tags and the ownership untouched. A
// generator advanced in local order produces a different field; one keyed on
// the tag does not. That gives a second arrangement at the same rank count, so
// order dependence and rank dependence can be told apart.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <mpi.h>

#include "gmd/system/box.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

namespace {

int global_rank = 0;
int global_size = 1;

constexpr std::size_t kAtomCount = 24;
constexpr double kTargetTemperature = 300.0;
constexpr std::uint32_t kSeed = 20260830u;

double mass_for(std::size_t tag) {
    static const std::array<double, 8> kMasses = {
        1.008, 12.011, 15.999, 4.0026, 6.941, 10.811, 14.007, 32.065};
    return kMasses[tag % kMasses.size()];
}

std::array<double, 3> position_for(std::size_t tag) {
    const double t = static_cast<double>(tag);
    return {2.0 + 1.7 * t, 3.0 + 1.1 * std::fmod(t, 7.0), 4.0 + 1.3 * std::fmod(t, 5.0)};
}

// Which tags this rank owns. At np>=4 the last rank deliberately owns none: an
// empty domain is legal and must still enter every collective.
std::vector<int> owned_tags() {
    std::vector<int> owned;
    const int contributing = (global_size >= 4) ? global_size - 1 : global_size;
    for (int tag = 0; tag < static_cast<int>(kAtomCount); ++tag) {
        if (tag % contributing == global_rank) owned.push_back(tag);
    }
    return owned;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    std::string output;
    bool reverse_storage = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            output = argv[++i];
        } else if (std::strcmp(argv[i], "--reverse-storage") == 0) {
            reverse_storage = true;
        }
    }
    if (output.empty()) {
        if (global_rank == 0) std::cerr << "usage: --out <path> [--reverse-storage]\n";
        MPI_Finalize();
        return 2;
    }

    std::vector<int> owned = owned_tags();
    if (reverse_storage) std::reverse(owned.begin(), owned.end());

    gmd::System system;
    system.resize(owned.size(), owned.size());
    gmd::Box box;
    box.set_lengths({60.0, 60.0, 60.0});
    system.set_box(box);
    for (std::size_t slot = 0; slot < owned.size(); ++slot) {
        const auto tag = static_cast<std::size_t>(owned[slot]);
        system.mutable_masses()[slot] = mass_for(tag);
        system.mutable_coordinates()[slot] = position_for(tag);
        system.mutable_atom_tags()[slot] = owned[slot];
        system.mutable_atom_owners()[slot] = global_rank;
    }

    gmd::VelocityInitializer initializer(kSeed);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, true);

    // Gather (tag, vx, vy, vz) to rank 0.
    const int local_count = static_cast<int>(system.num_local_atoms());
    std::vector<double> local(static_cast<std::size_t>(local_count) * 4);
    for (int i = 0; i < local_count; ++i) {
        const auto v = system.velocities()[static_cast<std::size_t>(i)];
        local[static_cast<std::size_t>(i) * 4 + 0] =
            static_cast<double>(system.atom_tag(static_cast<std::size_t>(i)));
        local[static_cast<std::size_t>(i) * 4 + 1] = v[0];
        local[static_cast<std::size_t>(i) * 4 + 2] = v[1];
        local[static_cast<std::size_t>(i) * 4 + 3] = v[2];
    }

    std::vector<int> counts(static_cast<std::size_t>(global_size), 0);
    const int local_doubles = local_count * 4;
    MPI_Gather(&local_doubles, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> displacements(static_cast<std::size_t>(global_size), 0);
    int total = 0;
    if (global_rank == 0) {
        for (int r = 0; r < global_size; ++r) {
            displacements[static_cast<std::size_t>(r)] = total;
            total += counts[static_cast<std::size_t>(r)];
        }
    }
    std::vector<double> gathered(global_rank == 0 ? static_cast<std::size_t>(total) : 0);
    MPI_Gatherv(local.data(), local_doubles, MPI_DOUBLE,
                gathered.data(), counts.data(), displacements.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Global momentum and kinetic energy, reduced the same way the initializer
    // reduces them, so the artifact records what the whole system carries.
    double local_momentum[3] = {0.0, 0.0, 0.0};
    double local_twice_ke = 0.0;
    for (int i = 0; i < local_count; ++i) {
        const auto v = system.velocities()[static_cast<std::size_t>(i)];
        const double m = system.masses()[static_cast<std::size_t>(i)];
        for (int d = 0; d < 3; ++d) local_momentum[d] += m * v[d];
        local_twice_ke += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    double momentum[3] = {0.0, 0.0, 0.0};
    double twice_ke = 0.0;
    MPI_Allreduce(local_momentum, momentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_twice_ke, &twice_ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (global_rank == 0) {
        std::vector<std::array<double, 4>> rows;
        rows.reserve(static_cast<std::size_t>(total) / 4);
        for (std::size_t i = 0; i + 3 < gathered.size(); i += 4) {
            rows.push_back({gathered[i], gathered[i + 1], gathered[i + 2], gathered[i + 3]});
        }
        std::sort(rows.begin(), rows.end(),
                  [](const auto& a, const auto& b) { return a[0] < b[0]; });

        std::ofstream out(output);
        out << std::setprecision(17);
        out << "# np " << global_size << " reverse_storage "
            << (reverse_storage ? 1 : 0) << '\n';
        out << "# atoms " << rows.size() << '\n';
        out << "momentum " << momentum[0] << ' ' << momentum[1] << ' ' << momentum[2] << '\n';
        out << "twice_kinetic_energy " << twice_ke << '\n';
        for (const auto& row : rows) {
            out << static_cast<long long>(row[0]) << ' '
                << row[1] << ' ' << row[2] << ' ' << row[3] << '\n';
        }
    }

    MPI_Finalize();
    return 0;
}
