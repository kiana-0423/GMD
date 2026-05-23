#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/io/checkpoint.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/core/runtime_context.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

constexpr double box_length = 10.0;
constexpr double position_tolerance = 1e-12;
constexpr double force_tolerance = 1e-10;

gmd::Box make_box() {
    gmd::Box box;
    box.set_lengths({box_length, box_length, box_length});
    return box;
}

gmd::System make_local_atom(const gmd::Box& box, double x, int tag, int owner) {
    gmd::System system;
    system.resize(1, 1);
    system.set_box(box);
    system.mutable_masses()[0] = 1.0;
    system.mutable_charges()[0] = 0.0;
    system.mutable_coordinates()[0] = {x, 5.0, 5.0};
    system.mutable_atom_tags()[0] = tag;
    system.mutable_atom_owners()[0] = owner;
    system.mutable_molecule_ids()[0] = 1000 + tag;
    return system;
}

gmd::DomainDecomposition make_decomposition(const gmd::Box& box,
                                            int size,
                                            int rank,
                                            bool periodic_x) {
    gmd::DomainDecomposition decomposition;
    decomposition.create_1d_decomposition(box, size, rank, 1.0, 0.0, periodic_x);
    return decomposition;
}

bool close(double lhs, double rhs, double tolerance) {
    return std::abs(lhs - rhs) <= tolerance;
}

void check(bool condition, const std::string& message, int rank, int& failures) {
    if (condition) {
        return;
    }

    std::cerr << "[mpi periodic 1d rank " << rank << "] " << message << '\n';
    ++failures;
}

gmd::ForceResult compute_lj(gmd::System& system, gmd::RuntimeContext& runtime) {
    gmd::ClassicalForceProvider provider(0.2, 0.5, 1.5);
    const auto coordinates = system.coordinates();
    const gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .coordinates =
            std::span<const gmd::Coordinate3D>(coordinates.data(), coordinates.size()),
    };

    gmd::ForceResult result;
    provider.compute(request, result, runtime);
    return result;
}

void copy_forces(gmd::System& system, const gmd::ForceResult& result) {
    auto forces = system.mutable_forces();
    for (std::size_t atom_index = 0; atom_index < forces.size(); ++atom_index) {
        forces[atom_index] = result.forces[atom_index];
    }
}

void test_periodic_ghost_exchange(const gmd::MpiCommunicator& communicator,
                                  int rank,
                                  int size,
                                  int& failures) {
    const gmd::Box box = make_box();
    gmd::System system = make_local_atom(box, rank == 0 ? 0.2 : 9.8, rank, rank);
    gmd::DomainDecomposition periodic = make_decomposition(box, size, rank, true);

    communicator.exchange_ghost_coordinates(system, periodic);

    check(system.num_ghost_atoms() == 1,
          "periodic edge exchange should add exactly one wraparound ghost",
          rank,
          failures);
    if (system.num_ghost_atoms() == 1) {
        const std::size_t ghost = system.num_local_atoms();
        check(system.atom_tag(ghost) == (rank == 0 ? 1 : 0),
              "wraparound ghost has the wrong tag",
              rank,
              failures);
        check(system.atom_owner(ghost) == (rank == 0 ? 1 : 0),
              "wraparound ghost has the wrong owner rank",
              rank,
              failures);
        check(close(system.coordinates()[ghost][0], rank == 0 ? -0.2 : 10.2,
                    position_tolerance),
              "wraparound ghost is not shifted into the local halo image",
              rank,
              failures);
    }

    system = make_local_atom(box, rank == 0 ? 0.2 : 9.8, rank, rank);
    gmd::DomainDecomposition nonperiodic = make_decomposition(box, size, rank, false);
    communicator.exchange_ghost_coordinates(system, nonperiodic);
    check(system.num_ghost_atoms() == 0,
          "nonperiodic x endpoints must not exchange wraparound ghosts",
          rank,
          failures);
}

void test_periodic_force_consistency(const gmd::MpiCommunicator& communicator,
                                     gmd::RuntimeContext& runtime,
                                     int rank,
                                     int size,
                                     int& failures) {
    const gmd::Box box = make_box();
    gmd::System distributed = make_local_atom(box, rank == 0 ? 0.2 : 9.8, rank, rank);
    const gmd::DomainDecomposition decomposition = make_decomposition(box, size, rank, true);
    communicator.exchange_ghost_coordinates(distributed, decomposition);

    const gmd::ForceResult distributed_result = compute_lj(distributed, runtime);
    copy_forces(distributed, distributed_result);
    communicator.reverse_accumulate_ghost_forces(distributed, decomposition);

    std::vector<double> local_forces(6, 0.0);
    for (std::size_t atom_index = 0; atom_index < distributed.num_local_atoms(); ++atom_index) {
        const int tag = distributed.atom_tag(atom_index);
        for (std::size_t dim = 0; dim < 3; ++dim) {
            local_forces[static_cast<std::size_t>(3 * tag) + dim] =
                distributed.forces()[atom_index][dim];
        }
    }

    std::vector<double> global_forces;
    communicator.allreduce_vector(local_forces, global_forces);
    const double global_energy =
        communicator.allreduce_scalar(distributed_result.potential_energy);

    gmd::System serial;
    serial.resize(2, 2);
    serial.set_box(box);
    serial.mutable_masses()[0] = 1.0;
    serial.mutable_masses()[1] = 1.0;
    serial.mutable_coordinates()[0] = {0.2, 5.0, 5.0};
    serial.mutable_coordinates()[1] = {9.8, 5.0, 5.0};
    serial.mutable_atom_tags()[0] = 0;
    serial.mutable_atom_tags()[1] = 1;
    const gmd::ForceResult serial_result = compute_lj(serial, runtime);

    check(close(global_energy, serial_result.potential_energy, force_tolerance),
          "periodic rank-edge LJ energy differs from serial",
          rank,
          failures);
    for (std::size_t atom_index = 0; atom_index < serial_result.forces.size(); ++atom_index) {
        for (std::size_t dim = 0; dim < 3; ++dim) {
            const std::size_t force_index = 3 * atom_index + dim;
            check(close(global_forces[force_index],
                        serial_result.forces[atom_index][dim],
                        force_tolerance),
                  "periodic rank-edge LJ force differs from serial",
                  rank,
                  failures);
        }
    }
}

void test_periodic_migration(const gmd::MpiCommunicator& communicator,
                             int rank,
                             int size,
                             int& failures) {
    const gmd::Box box = make_box();
    gmd::System system = make_local_atom(box,
                                         rank == 0 ? -0.1 : 10.1,
                                         rank == 0 ? 10 : 20,
                                         rank);
    const gmd::DomainDecomposition decomposition = make_decomposition(box, size, rank, true);

    communicator.redistribute_atoms(system, decomposition);

    check(system.num_local_atoms() == 1 && system.num_ghost_atoms() == 0,
          "periodic edge migration should keep one owned atom on each rank",
          rank,
          failures);
    if (system.num_local_atoms() == 1) {
        check(system.atom_tag(0) == (rank == 0 ? 20 : 10),
              "periodic edge migration delivered the wrong atom",
              rank,
              failures);
        check(system.molecule_ids()[0] == 1000 + (rank == 0 ? 20 : 10),
              "periodic edge migration did not preserve molecule id",
              rank,
              failures);
        check(system.atom_owner(0) == rank,
              "periodic edge migration did not reset the local owner rank",
              rank,
              failures);
        check(close(system.coordinates()[0][0], rank == 0 ? 0.1 : 9.9,
                    position_tolerance),
              "periodic edge migration did not wrap the x coordinate",
              rank,
              failures);

        const auto checkpoint_path =
            std::filesystem::current_path() /
            ("mpi_periodic_migration_rank_" + std::to_string(rank) + ".gmdchk");
        gmd::CheckpointMetadata metadata;
        metadata.step = 1;
        metadata.time_fs = 0.5;
        metadata.config_summary = "mpi periodic migration molecule-id test";
        gmd::CheckpointData checkpoint{metadata, &system, nullptr};
        gmd::write_checkpoint(checkpoint_path, checkpoint);

        gmd::System checkpoint_system;
        gmd::read_checkpoint(checkpoint_path, checkpoint_system, nullptr);
        check(checkpoint_system.atom_count() == 1 &&
                  checkpoint_system.molecule_ids()[0] == system.molecule_ids()[0],
              "checkpoint after migration did not retain molecule id",
              rank,
              failures);
        std::filesystem::remove(checkpoint_path);
    }
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment mpi_environment(argc, argv);
#endif

    gmd::MpiCommunicator communicator;
    gmd::RuntimeContext runtime;
    const int rank = communicator.rank();
    const int size = communicator.size();
    int failures = 0;

    check(size == 2, "test requires exactly two MPI ranks", rank, failures);
    if (size == 2) {
        const std::string test_case = argc > 1 ? argv[1] : "all";
        if (test_case == "ghost" || test_case == "all") {
            test_periodic_ghost_exchange(communicator, rank, size, failures);
        }
        if (test_case == "force" || test_case == "all") {
            test_periodic_force_consistency(communicator, runtime, rank, size, failures);
        }
        if (test_case == "migration" || test_case == "all") {
            test_periodic_migration(communicator, rank, size, failures);
        }
        check(test_case == "ghost" ||
                  test_case == "force" ||
                  test_case == "migration" ||
                  test_case == "all",
              "unknown test case " + test_case,
              rank,
              failures);
    }

    const double global_failures = communicator.allreduce_scalar(static_cast<double>(failures));
    return global_failures == 0.0 ? 0 : 1;
}
