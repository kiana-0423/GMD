#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/core/runtime_context.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

constexpr double box_length = 8.0;
constexpr double halo_width = 0.5;
constexpr double coordinate_tolerance = 1e-12;
constexpr double force_tolerance = 1e-10;

using Vec3 = gmd::System::Vec3;

struct TaggedAtom {
    int tag;
    Vec3 position;
};

enum class HaloKind {
    face,
    edge,
    corner,
};

gmd::Box make_box() {
    gmd::Box box;
    box.set_lengths({box_length, box_length, box_length});
    return box;
}

gmd::DomainDecomposition make_decomposition(const gmd::Box& box,
                                            int rank,
                                            bool periodic) {
    gmd::DomainDecomposition decomposition;
    decomposition.create_decomposition(box,
                                       {2, 2, 2},
                                       rank,
                                       halo_width,
                                       0.0,
                                       {periodic, periodic, periodic});
    return decomposition;
}

gmd::System make_system(const gmd::Box& box,
                        const std::vector<TaggedAtom>& atoms,
                        int owner) {
    gmd::System system;
    system.resize(atoms.size(), atoms.size());
    system.set_box(box);
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        system.mutable_masses()[index] = 1.0;
        system.mutable_charges()[index] = 0.25;
        system.mutable_atom_types()[index] = 0;
        system.mutable_atom_tags()[index] = atoms[index].tag;
        system.mutable_atom_owners()[index] = owner;
        system.mutable_coordinates()[index] = atoms[index].position;
    }
    return system;
}

bool close(double lhs, double rhs, double tolerance = coordinate_tolerance) {
    return std::abs(lhs - rhs) <= tolerance;
}

void check(bool condition, const std::string& message, int rank, int& failures) {
    if (condition) {
        return;
    }

    std::cerr << "[mpi ghost 3d rank " << rank << "] " << message << '\n';
    ++failures;
}

int source_rank(HaloKind kind) {
    switch (kind) {
    case HaloKind::face:
        return 1;
    case HaloKind::edge:
        return 3;
    case HaloKind::corner:
        return 7;
    }

    return -1;
}

int source_tag(HaloKind kind) {
    return 100 + source_rank(kind);
}

Vec3 source_position(HaloKind kind) {
    switch (kind) {
    case HaloKind::face:
        return {4.1, 2.0, 2.0};
    case HaloKind::edge:
        return {4.1, 4.1, 2.0};
    case HaloKind::corner:
        return {4.1, 4.1, 4.1};
    }

    return {0.0, 0.0, 0.0};
}

Vec3 returned_force(HaloKind kind) {
    switch (kind) {
    case HaloKind::face:
        return {1.0, 2.0, 3.0};
    case HaloKind::edge:
        return {4.0, 5.0, 6.0};
    case HaloKind::corner:
        return {7.0, 8.0, 9.0};
    }

    return {0.0, 0.0, 0.0};
}

std::size_t count_ghost_tag(const gmd::System& system, int tag) {
    std::size_t count = 0;
    for (std::size_t index = system.num_local_atoms(); index < system.atom_count(); ++index) {
        if (system.atom_tag(index) == tag) {
            ++count;
        }
    }
    return count;
}

std::size_t find_ghost_tag(const gmd::System& system, int tag) {
    for (std::size_t index = system.num_local_atoms(); index < system.atom_count(); ++index) {
        if (system.atom_tag(index) == tag) {
            return index;
        }
    }
    return system.atom_count();
}

void check_position(const Vec3& actual,
                    const Vec3& expected,
                    const std::string& label,
                    int rank,
                    int& failures) {
    for (std::size_t dim = 0; dim < actual.size(); ++dim) {
        check(close(actual[dim], expected[dim]),
              label + " has the wrong coordinate",
              rank,
              failures);
    }
}

void test_direct_ghost_exchange(const gmd::MpiCommunicator& communicator,
                                HaloKind kind,
                                int rank,
                                int& failures) {
    const gmd::Box box = make_box();
    const int sender = source_rank(kind);
    const std::vector<TaggedAtom> atoms = rank == sender
        ? std::vector<TaggedAtom>{{source_tag(kind), source_position(kind)}}
        : std::vector<TaggedAtom>{};
    gmd::System system = make_system(box, atoms, rank);
    const gmd::DomainDecomposition decomposition = make_decomposition(box, rank, false);

    communicator.exchange_ghost_coordinates(system, decomposition);

    if (rank != 0) {
        return;
    }

    check(system.num_ghost_atoms() == 1,
          "direct face/edge/corner path should deliver one ghost to rank zero",
          rank,
          failures);
    check(count_ghost_tag(system, source_tag(kind)) == 1,
          "ghost exchange should not duplicate an owner/tag pair",
          rank,
          failures);
    const std::size_t ghost = find_ghost_tag(system, source_tag(kind));
    check(ghost != system.atom_count(), "expected 3D ghost was not received", rank, failures);
    if (ghost != system.atom_count()) {
        check(system.atom_owner(ghost) == sender,
              "3D ghost has the wrong home rank",
              rank,
              failures);
        check_position(system.coordinates()[ghost],
                       source_position(kind),
                       "3D ghost",
                       rank,
                       failures);
    }
}

void test_periodic_corner_exchange(const gmd::MpiCommunicator& communicator,
                                   int rank,
                                   int& failures) {
    const gmd::Box box = make_box();
    const int tag = 707;
    const std::vector<TaggedAtom> atoms = rank == 7
        ? std::vector<TaggedAtom>{{tag, {7.9, 7.9, 7.9}}}
        : std::vector<TaggedAtom>{};
    gmd::System system = make_system(box, atoms, rank);
    const gmd::DomainDecomposition decomposition = make_decomposition(box, rank, true);

    communicator.exchange_ghost_coordinates(system, decomposition);

    if (rank != 0) {
        return;
    }

    check(system.num_ghost_atoms() == 1,
          "periodic corner exchange should deliver one wrapped ghost",
          rank,
          failures);
    const std::size_t ghost = find_ghost_tag(system, tag);
    check(ghost != system.atom_count(), "periodic corner ghost was not received", rank, failures);
    if (ghost != system.atom_count()) {
        check(system.atom_owner(ghost) == 7,
              "periodic corner ghost has the wrong owner",
              rank,
              failures);
        check_position(system.coordinates()[ghost],
                       {-0.1, -0.1, -0.1},
                       "periodic corner ghost",
                       rank,
                       failures);
    }
}

void test_reverse_force(const gmd::MpiCommunicator& communicator,
                        HaloKind kind,
                        int rank,
                        int& failures) {
    const gmd::Box box = make_box();
    const int owner = source_rank(kind);
    const std::vector<TaggedAtom> atoms = rank == owner
        ? std::vector<TaggedAtom>{{source_tag(kind), source_position(kind)}}
        : std::vector<TaggedAtom>{};
    gmd::System system = make_system(box, atoms, rank);
    const gmd::DomainDecomposition decomposition = make_decomposition(box, rank, false);

    communicator.exchange_ghost_coordinates(system, decomposition);
    if (rank == 0) {
        const std::size_t ghost = find_ghost_tag(system, source_tag(kind));
        check(ghost != system.atom_count(),
              "reverse-force setup did not receive its source ghost",
              rank,
              failures);
        if (ghost != system.atom_count()) {
            system.mutable_forces()[ghost] = returned_force(kind);
        }
    }

    communicator.reverse_accumulate_ghost_forces(system, decomposition);

    check(system.num_ghost_atoms() == 0,
          "reverse force accumulation should clear ghost atoms",
          rank,
          failures);
    if (rank != owner) {
        return;
    }

    check(system.num_local_atoms() == 1,
          "reverse force owner should retain its local source atom",
          rank,
          failures);
    if (system.num_local_atoms() == 1) {
        check_position(system.forces()[0],
                       returned_force(kind),
                       "reverse force",
                       rank,
                       failures);
    }
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

void test_periodic_lj_consistency(const gmd::MpiCommunicator& communicator,
                                  gmd::RuntimeContext& runtime,
                                  int rank,
                                  int& failures) {
    const gmd::Box box = make_box();
    std::vector<TaggedAtom> atoms;
    if (rank == 0) {
        atoms.push_back({0, {0.2, 0.2, 0.2}});
    } else if (rank == 7) {
        atoms.push_back({1, {7.8, 7.8, 7.8}});
    }

    gmd::System distributed = make_system(box, atoms, rank);
    gmd::DomainDecomposition decomposition;
    decomposition.create_decomposition(box, {2, 2, 2}, rank, 1.5, 0.0, {true, true, true});
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

    gmd::System serial = make_system(box,
                                     {{0, {0.2, 0.2, 0.2}},
                                      {1, {7.8, 7.8, 7.8}}},
                                     0);
    const gmd::ForceResult serial_result = compute_lj(serial, runtime);

    check(close(global_energy, serial_result.potential_energy, force_tolerance),
          "8-rank periodic LJ energy differs from serial",
          rank,
          failures);
    for (std::size_t atom_index = 0; atom_index < serial_result.forces.size(); ++atom_index) {
        for (std::size_t dim = 0; dim < 3; ++dim) {
            const std::size_t index = 3 * atom_index + dim;
            check(close(global_forces[index],
                        serial_result.forces[atom_index][dim],
                        force_tolerance),
                  "8-rank periodic LJ force differs from serial",
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
    gmd::RuntimeContext runtime;
    const int rank = communicator.rank();
    const int size = communicator.size();
    const std::string test_case = argc > 1 ? argv[1] : "";
    int failures = 0;

    check(size == 8, "test requires exactly eight MPI ranks", rank, failures);
    if (size == 8) {
        if (test_case == "face") {
            test_direct_ghost_exchange(communicator, HaloKind::face, rank, failures);
        } else if (test_case == "edge") {
            test_direct_ghost_exchange(communicator, HaloKind::edge, rank, failures);
        } else if (test_case == "corner") {
            test_direct_ghost_exchange(communicator, HaloKind::corner, rank, failures);
        } else if (test_case == "periodic") {
            test_periodic_corner_exchange(communicator, rank, failures);
        } else if (test_case == "reverse-face") {
            test_reverse_force(communicator, HaloKind::face, rank, failures);
        } else if (test_case == "reverse-edge") {
            test_reverse_force(communicator, HaloKind::edge, rank, failures);
        } else if (test_case == "reverse-corner") {
            test_reverse_force(communicator, HaloKind::corner, rank, failures);
        } else if (test_case == "lj") {
            test_periodic_lj_consistency(communicator, runtime, rank, failures);
        } else {
            check(false, "unknown test case " + test_case, rank, failures);
        }
    }

    const double global_failures = communicator.allreduce_scalar(static_cast<double>(failures));
    return global_failures == 0.0 ? 0 : 1;
}
