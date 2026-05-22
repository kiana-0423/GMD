#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

constexpr double tolerance = 1e-12;

struct TaggedPosition {
    int tag;
    gmd::System::Vec3 position;
};

gmd::Box make_box() {
    gmd::Box box;
    box.set_lengths({8.0, 8.0, 8.0});
    return box;
}

gmd::System make_system(const gmd::Box& box,
                        const std::vector<TaggedPosition>& atoms,
                        int owner) {
    gmd::System system;
    system.resize(atoms.size(), atoms.size());
    system.set_box(box);
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        system.mutable_masses()[index] = 1.0;
        system.mutable_coordinates()[index] = atoms[index].position;
        system.mutable_atom_tags()[index] = atoms[index].tag;
        system.mutable_atom_owners()[index] = owner;
    }
    return system;
}

bool close(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= tolerance;
}

void check(bool condition, const std::string& message, int rank, int& failures) {
    if (condition) {
        return;
    }

    std::cerr << "[mpi domain 3d rank " << rank << "] " << message << '\n';
    ++failures;
}

int grid_size(const std::array<int, 3>& grid) {
    return grid[0] * grid[1] * grid[2];
}

gmd::System::Vec3 cell_center(const std::array<int, 3>& coord) {
    return {
        2.0 + 4.0 * static_cast<double>(coord[0]),
        2.0 + 4.0 * static_cast<double>(coord[1]),
        2.0 + 4.0 * static_cast<double>(coord[2]),
    };
}

std::size_t find_tag(const gmd::System& system, int tag) {
    for (std::size_t index = 0; index < system.num_local_atoms(); ++index) {
        if (system.atom_tag(index) == tag) {
            return index;
        }
    }
    return system.num_local_atoms();
}

void test_four_rank_grid(int rank, int size, int& failures) {
    const gmd::Box box = make_box();
    gmd::DomainDecomposition automatic;
    automatic.create_decomposition(box, size, rank, 0.25, 0.25, {true, true, true});
    const auto auto_grid = automatic.info().proc_grid;
    check(grid_size(auto_grid) == size,
          "automatic processor grid product must match MPI size",
          rank,
          failures);
    check(auto_grid[0] == 1 || auto_grid[1] == 1 || auto_grid[2] == 1,
          "four ranks should leave one processor-grid axis at extent one",
          rank,
          failures);
    check(automatic.proc_coord_from_rank(rank) == automatic.info().proc_coord,
          "automatic rank-to-coordinate mapping is inconsistent",
          rank,
          failures);

    gmd::DomainDecomposition periodic;
    periodic.create_decomposition(box, {2, 2, 1}, rank, 0.25, 0.5, {true, true, true});
    const auto coord = periodic.info().proc_coord;
    check(periodic.rank_from_proc_coord(coord) == rank,
          "explicit coordinate-to-rank mapping is inconsistent",
          rank,
          failures);
    check(periodic.info().proc_grid[2] == 1,
          "explicit processor grid should retain its unit z axis",
          rank,
          failures);
    check(close(periodic.info().owned_lo[0], 4.0 * static_cast<double>(coord[0])) &&
              close(periodic.info().owned_hi[0],
                    4.0 * static_cast<double>(coord[0] + 1)),
          "x ownership bounds are wrong",
          rank,
          failures);
    check(close(periodic.info().owned_lo[1], 4.0 * static_cast<double>(coord[1])) &&
              close(periodic.info().owned_hi[1],
                    4.0 * static_cast<double>(coord[1] + 1)),
          "y ownership bounds are wrong",
          rank,
          failures);
    check(close(periodic.info().owned_lo[2], 0.0) &&
              close(periodic.info().owned_hi[2], 8.0),
          "unit z-axis ownership bounds are wrong",
          rank,
          failures);
    check(periodic.owner_rank(box, {6.0, 6.0, 2.0}) == 3,
          "owner_rank should map xyz positions through the 3D process grid",
          rank,
          failures);
    check(periodic.owner_rank(box, {-0.1, 8.1, 8.1}) == 1,
          "periodic xyz ownership should wrap each dimension",
          rank,
          failures);
    check(periodic.neighbor_rank({0, 0, 1}) == rank,
          "periodic neighbor on a unit axis should be the local rank",
          rank,
          failures);

    gmd::DomainDecomposition nonperiodic;
    nonperiodic.create_decomposition(box, {2, 2, 1}, rank, 0.0, 0.0, {false, false, false});
    check(nonperiodic.owner_rank(box, {-1.0, 9.0, 2.0}) == 2,
          "nonperiodic ownership should clamp at x/y domain edges",
          rank,
          failures);
    if (nonperiodic.info().proc_coord[0] == 0) {
        check(nonperiodic.neighbor_rank({-1, 0, 0}) ==
                  gmd::DomainDecomposition::no_rank,
              "nonperiodic low-x face should not have a neighbor",
              rank,
              failures);
    }
}

void test_xyz_migration(const gmd::MpiCommunicator& communicator,
                        int rank,
                        int& failures) {
    const gmd::Box box = make_box();
    gmd::DomainDecomposition decomposition;
    decomposition.create_decomposition(box, {2, 2, 2}, rank, 0.0, 0.0, {true, true, true});

    std::vector<TaggedPosition> atoms;
    const auto coord = decomposition.info().proc_coord;
    for (int dim = 0; dim < 3; ++dim) {
        auto target = coord;
        target[static_cast<std::size_t>(dim)] =
            1 - target[static_cast<std::size_t>(dim)];
        atoms.push_back(TaggedPosition{
            .tag = dim * 100 + rank,
            .position = cell_center(target),
        });
    }

    gmd::System system = make_system(box, atoms, rank);
    communicator.redistribute_atoms(system, decomposition);

    check(system.num_local_atoms() == 3,
          "xyz migration should deliver one atom from each face direction",
          rank,
          failures);
    for (int dim = 0; dim < 3; ++dim) {
        auto source = coord;
        source[static_cast<std::size_t>(dim)] =
            1 - source[static_cast<std::size_t>(dim)];
        const int source_rank = decomposition.rank_from_proc_coord(source);
        check(find_tag(system, dim * 100 + source_rank) != system.num_local_atoms(),
              "xyz migration delivered the wrong owner set",
              rank,
              failures);
    }
}

void test_xyz_periodic_wrap(const gmd::MpiCommunicator& communicator,
                            int rank,
                            int& failures) {
    const gmd::Box box = make_box();
    gmd::DomainDecomposition decomposition;
    decomposition.create_decomposition(box, {2, 2, 2}, rank, 0.0, 0.0, {true, true, true});

    const auto coord = decomposition.info().proc_coord;
    std::vector<TaggedPosition> atoms;
    for (int dim = 0; dim < 3; ++dim) {
        auto position = cell_center(coord);
        const int side = coord[static_cast<std::size_t>(dim)];
        position[static_cast<std::size_t>(dim)] = side == 0 ? -0.1 : 8.1;
        atoms.push_back(TaggedPosition{
            .tag = 1000 * dim + 10 * rank + side,
            .position = position,
        });
    }

    gmd::System system = make_system(box, atoms, rank);
    communicator.redistribute_atoms(system, decomposition);

    check(system.num_local_atoms() == 3,
          "periodic xyz wrap should deliver one wrapped atom per axis",
          rank,
          failures);
    for (int dim = 0; dim < 3; ++dim) {
        auto source = coord;
        source[static_cast<std::size_t>(dim)] =
            1 - source[static_cast<std::size_t>(dim)];
        const int source_rank = decomposition.rank_from_proc_coord(source);
        const int source_side = source[static_cast<std::size_t>(dim)];
        const int expected_tag = 1000 * dim + 10 * source_rank + source_side;
        const std::size_t local_index = find_tag(system, expected_tag);
        check(local_index != system.num_local_atoms(),
              "periodic xyz wrap delivered the wrong tag",
              rank,
              failures);
        if (local_index != system.num_local_atoms()) {
            const double expected_coordinate = coord[static_cast<std::size_t>(dim)] == 0
                ? 0.1
                : 7.9;
            check(close(system.coordinates()[local_index][static_cast<std::size_t>(dim)],
                        expected_coordinate),
                  "periodic xyz wrap did not normalize the coordinate",
                  rank,
                  failures);
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment mpi_environment(argc, argv);
#endif

    gmd::MpiCommunicator communicator;
    const int rank = communicator.rank();
    const int size = communicator.size();
    const std::string test_case = argc > 1 ? argv[1] : "";
    int failures = 0;

    if (test_case == "grid4") {
        check(size == 4, "grid4 requires exactly four MPI ranks", rank, failures);
        if (size == 4) {
            test_four_rank_grid(rank, size, failures);
        }
    } else if (test_case == "migrate8") {
        check(size == 8, "migrate8 requires exactly eight MPI ranks", rank, failures);
        if (size == 8) {
            test_xyz_migration(communicator, rank, failures);
            test_xyz_periodic_wrap(communicator, rank, failures);
        }
    } else {
        check(false, "expected grid4 or migrate8 test case", rank, failures);
    }

    const double global_failures = communicator.allreduce_scalar(static_cast<double>(failures));
    return global_failures == 0.0 ? 0 : 1;
}
