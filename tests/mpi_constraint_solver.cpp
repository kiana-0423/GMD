#include <cmath>
#include <iostream>
#include <vector>

#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/system.hpp"

#include <mpi.h>

namespace {

double local_distance_error(const gmd::System& system) {
    const auto coords = system.coordinates();
    if (system.num_local_atoms() != 1) {
        return 0.0;
    }
    return coords[0][0];
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment env(argc, argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        if (rank == 0) {
            std::cerr << "mpi_constraint_solver requires exactly 2 ranks\n";
        }
        return 1;
    }

    gmd::System system;
    system.resize(1, 1);
    gmd::Box box;
    box.set_lengths({20.0, 20.0, 20.0});
    system.set_box(box);
    system.mutable_masses()[0] = 1.0;
    system.mutable_atom_tags()[0] = rank;
    system.mutable_atom_owners()[0] = rank;
    system.mutable_coordinates()[0] = rank == 0
        ? gmd::System::Vec3{9.4, 10.0, 10.0}
        : gmd::System::Vec3{10.8, 10.0, 10.0};
    system.mutable_velocities()[0] = rank == 0
        ? gmd::System::Vec3{0.2, 0.0, 0.0}
        : gmd::System::Vec3{-0.1, 0.0, 0.0};

    gmd::ConstraintSolver solver({gmd::BondConstraint{0, 1, 1.0}},
                                 gmd::ConstraintSettings{1.0e-10, 100, true});
    solver.apply_shake(system);
    solver.apply_rattle(system);

    double local_x = local_distance_error(system);
    double local_vx = system.velocities()[0][0];
    std::vector<double> xs(2, 0.0);
    std::vector<double> vxs(2, 0.0);
    MPI_Allgather(&local_x, 1, MPI_DOUBLE, xs.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&local_vx, 1, MPI_DOUBLE, vxs.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    const double length = std::abs(xs[1] - xs[0]);
    const double rdotv = (xs[0] - xs[1]) * (vxs[0] - vxs[1]);

    if (rank == 0 && std::abs(length - 1.0) > 1.0e-9) {
        std::cerr << "MPI SHAKE cross-rank length error: " << std::abs(length - 1.0) << "\n";
        return 1;
    }
    if (rank == 0 && std::abs(rdotv) > 1.0e-9) {
        std::cerr << "MPI RATTLE cross-rank velocity error: " << std::abs(rdotv) << "\n";
        return 1;
    }
    return 0;
}
