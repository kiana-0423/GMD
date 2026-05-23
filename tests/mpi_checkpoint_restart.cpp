#include <cmath>
#include <filesystem>
#include <iostream>

#include "gmd/io/checkpoint.hpp"
#include "gmd/system/system.hpp"

#include <mpi.h>

namespace {

gmd::System make_global_system() {
    gmd::System system;
    system.resize(4, 4);
    gmd::Box box;
    box.set_lengths({16.0, 17.0, 18.0});
    system.set_box(box);
    auto masses = system.mutable_masses();
    auto charges = system.mutable_charges();
    auto types = system.mutable_atom_types();
    auto molecules = system.mutable_molecule_ids();
    auto coords = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();
    for (std::size_t i = 0; i < system.atom_count(); ++i) {
        masses[i] = 1.0 + static_cast<double>(i);
        charges[i] = -0.05 * static_cast<double>(i);
        types[i] = static_cast<int>(i % 2);
        molecules[i] = static_cast<int>(i / 2);
        coords[i] = {static_cast<double>(i),
                     1.0 + static_cast<double>(i),
                     2.0 + static_cast<double>(i)};
        velocities[i] = {0.1 * static_cast<double>(i + 1),
                         0.2 * static_cast<double>(i + 1),
                         0.3 * static_cast<double>(i + 1)};
    }
    return system;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const auto path = std::filesystem::current_path() / "mpi_checkpoint_restart.gmdchk";
    if (rank == 0) {
        auto system = make_global_system();
        gmd::CheckpointMetadata metadata;
        metadata.step = 12;
        metadata.time_fs = 6.0;
        metadata.config_summary = "mpi checkpoint test";
        gmd::CheckpointData checkpoint{metadata, &system, nullptr};
        gmd::write_checkpoint(path, checkpoint);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    gmd::System restarted;
    const auto metadata = gmd::read_checkpoint(path, restarted, nullptr);
    int ok = 1;
    ok = ok && metadata.step == 12;
    ok = ok && restarted.atom_count() == 4;
    ok = ok && std::abs(restarted.coordinates()[3][2] - 5.0) < 1.0e-12;
    ok = ok && std::abs(restarted.velocities()[2][1] - 0.6) < 1.0e-12;

    int global_ok = 0;
    MPI_Allreduce(&ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Finalize();
    if (!global_ok) {
        if (rank == 0) {
            std::cerr << "MPI checkpoint restart test failed\n";
        }
        return 1;
    }
    return 0;
}
