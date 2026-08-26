#include "gmd/parallel/mpi_communicator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {
namespace {

#ifdef GMD_ENABLE_MPI
bool mpi_is_available() {
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    return is_initialized != 0 && is_finalized == 0;
}
#endif

constexpr int ghost_record_width = 7;
constexpr int reverse_force_record_width = 13;
#ifdef GMD_ENABLE_MPI
constexpr int atom_state_record_width = 12;
#endif

#ifdef GMD_ENABLE_MPI
int mpi_neighbor(int rank) {
    return rank == DomainDecomposition::no_rank ? MPI_PROC_NULL : rank;
}

void validate_rank_grid(const DomainDecomposition& dd, int mpi_rank, int mpi_size) {
    const DomainInfo& domain = dd.info();
    int nprocs = 1;
    for (int extent : domain.proc_grid) {
        nprocs *= extent;
    }
    if (nprocs != mpi_size ||
        dd.rank_from_proc_coord(domain.proc_coord) != mpi_rank) {
        throw std::invalid_argument("MPI rank does not match the domain process grid");
    }
}

using NeighborOffset = std::array<int, 3>;

std::vector<NeighborOffset> neighbor_offsets(const DomainInfo& domain) {
    std::vector<NeighborOffset> offsets;
    offsets.reserve(26);
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                const NeighborOffset offset{dx, dy, dz};
                if (offset == NeighborOffset{0, 0, 0}) {
                    continue;
                }

                bool uses_unit_axis = false;
                for (std::size_t dim = 0; dim < offset.size(); ++dim) {
                    if (offset[dim] != 0 && domain.proc_grid[dim] == 1) {
                        uses_unit_axis = true;
                        break;
                    }
                }
                if (!uses_unit_axis) {
                    offsets.push_back(offset);
                }
            }
        }
    }

    return offsets;
}

NeighborOffset opposite(const NeighborOffset& offset) {
    return {-offset[0], -offset[1], -offset[2]};
}

System::Vec3 periodic_shift(const DomainDecomposition& dd,
                            const Box& box,
                            const NeighborOffset& offset) {
    const DomainInfo& domain = dd.info();
    System::Vec3 shift{0.0, 0.0, 0.0};
    for (std::size_t dim = 0; dim < offset.size(); ++dim) {
        if (!domain.periodic[dim]) {
            continue;
        }
        if (offset[dim] < 0 && domain.proc_coord[dim] == 0) {
            shift[dim] = box.lengths[dim];
        } else if (offset[dim] > 0 &&
                   domain.proc_coord[dim] == domain.proc_grid[dim] - 1) {
            shift[dim] = -box.lengths[dim];
        }
    }

    return shift;
}
#endif

}  // namespace

double MpiCommunicator::allreduce_scalar(double local_value) const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        double global_value = 0.0;
        MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return global_value;
    }
#endif

    return local_value;
}

void MpiCommunicator::allreduce_vector(const std::vector<double>& local,
                                       std::vector<double>& global) const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        global.resize(local.size());
        MPI_Allreduce(local.data(),
                      global.data(),
                      static_cast<int>(local.size()),
                      MPI_DOUBLE,
                      MPI_SUM,
                      MPI_COMM_WORLD);
        return;
    }
#endif

    global = local;
}

bool MpiCommunicator::allreduce_logical_or(bool local_value) const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        int local_flag = local_value ? 1 : 0;
        int global_flag = 0;
        MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        return global_flag != 0;
    }
#endif

    return local_value;
}

void MpiCommunicator::broadcast_box(Box& box, int root) const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        MPI_Bcast(box.lengths.data(),
                  static_cast<int>(box.lengths.size()),
                  MPI_DOUBLE,
                  root,
                  MPI_COMM_WORLD);
        box.set_lengths(box.lengths);
        return;
    }
#else
    (void)root;
#endif

    box.set_lengths(box.lengths);
}

void MpiCommunicator::barrier() const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
}

void MpiCommunicator::exchange_ghost_coordinates(System& system,
                                                 const DomainDecomposition& dd) const {
    system.clear_ghost_atoms();

#ifdef GMD_ENABLE_MPI
    if (!mpi_is_available() || size() <= 1) {
        return;
    }

    const int my_rank = rank();
    validate_rank_grid(dd, my_rank, size());

    const DomainInfo& domain = dd.info();
    const Box& box = system.box();
    for (const auto& offset : neighbor_offsets(domain)) {
        const int send_neighbor = mpi_neighbor(dd.neighbor_rank(offset));
        const int recv_neighbor = mpi_neighbor(dd.neighbor_rank(opposite(offset)));
        std::vector<double> send_buffer = pack_send_buffer(system, dd, offset);
        const auto shift = periodic_shift(dd, box, offset);
        for (std::size_t record = 0; record < send_buffer.size();
             record += ghost_record_width) {
            for (std::size_t dim = 0; dim < shift.size(); ++dim) {
                send_buffer[record + dim] += shift[dim];
            }
        }
        if (send_buffer.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("Ghost exchange buffer is too large for MPI_Sendrecv");
        }

        const int send_count = static_cast<int>(send_buffer.size());
        int recv_count = 0;
        MPI_Sendrecv(&send_count,
                     1,
                     MPI_INT,
                     send_neighbor,
                     110,
                     &recv_count,
                     1,
                     MPI_INT,
                     recv_neighbor,
                     110,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        if (recv_count < 0 || recv_count % ghost_record_width != 0) {
            throw std::runtime_error("Ghost exchange received a malformed buffer size");
        }

        std::vector<double> recv_buffer(static_cast<std::size_t>(recv_count));
        MPI_Sendrecv(send_buffer.data(),
                     send_count,
                     MPI_DOUBLE,
                     send_neighbor,
                     111,
                     recv_buffer.data(),
                     recv_count,
                     MPI_DOUBLE,
                     recv_neighbor,
                     111,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        if (recv_neighbor != MPI_PROC_NULL) {
            unpack_recv_buffer(system, recv_buffer, recv_neighbor);
        }
    }
#else
    (void)dd;
#endif
}

void MpiCommunicator::reverse_accumulate_ghost_forces(
        System& system,
        const DomainDecomposition& dd) const {
    system.clear_reverse_ghost_virial();

#ifdef GMD_ENABLE_MPI
    if (!mpi_is_available() || size() <= 1) {
        system.clear_ghost_atoms();
        return;
    }

    const int my_rank = rank();
    const int mpi_size = size();
    validate_rank_grid(dd, my_rank, mpi_size);

    std::vector<std::vector<double>> per_rank(static_cast<std::size_t>(mpi_size));
    std::vector<int> send_counts(static_cast<std::size_t>(mpi_size), 0);
    for (int home_rank = 0; home_rank < mpi_size; ++home_rank) {
        auto packed = pack_reverse_force_buffer(system, home_rank);
        if (packed.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("Reverse force exchange buffer is too large for MPI");
        }
        send_counts[static_cast<std::size_t>(home_rank)] =
            static_cast<int>(packed.size());
        per_rank[static_cast<std::size_t>(home_rank)] = std::move(packed);
    }

    std::vector<int> recv_counts(static_cast<std::size_t>(mpi_size), 0);
    MPI_Alltoall(send_counts.data(),
                 1,
                 MPI_INT,
                 recv_counts.data(),
                 1,
                 MPI_INT,
                 MPI_COMM_WORLD);

    std::vector<int> send_displs(static_cast<std::size_t>(mpi_size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(mpi_size), 0);
    int send_total = 0;
    int recv_total = 0;
    for (int other_rank = 0; other_rank < mpi_size; ++other_rank) {
        send_displs[static_cast<std::size_t>(other_rank)] = send_total;
        recv_displs[static_cast<std::size_t>(other_rank)] = recv_total;
        send_total += send_counts[static_cast<std::size_t>(other_rank)];
        recv_total += recv_counts[static_cast<std::size_t>(other_rank)];
    }
    if (recv_total < 0 || recv_total % reverse_force_record_width != 0) {
        throw std::runtime_error("Reverse force exchange received a malformed buffer size");
    }

    std::vector<double> send_buffer(static_cast<std::size_t>(send_total));
    for (int other_rank = 0; other_rank < mpi_size; ++other_rank) {
        const auto& packed = per_rank[static_cast<std::size_t>(other_rank)];
        const auto displacement =
            static_cast<std::size_t>(send_displs[static_cast<std::size_t>(other_rank)]);
        std::copy(packed.begin(), packed.end(), send_buffer.begin() + displacement);
    }

    std::vector<double> recv_buffer(static_cast<std::size_t>(recv_total));
    MPI_Alltoallv(send_buffer.data(),
                  send_counts.data(),
                  send_displs.data(),
                  MPI_DOUBLE,
                  recv_buffer.data(),
                  recv_counts.data(),
                  recv_displs.data(),
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);
    for (int recv_count : recv_counts) {
        if (recv_count < 0 || recv_count % reverse_force_record_width != 0) {
            throw std::runtime_error("Reverse force exchange received a malformed buffer size");
        }
    }
    unpack_reverse_force_buffer(system, recv_buffer);
#else
    (void)dd;
#endif

    system.clear_ghost_atoms();
}

void MpiCommunicator::redistribute_atoms(System& system,
                                         const DomainDecomposition& dd) const {
#ifdef GMD_ENABLE_MPI
    if (!mpi_is_available() || size() <= 1) {
        return;
    }

    const int mpi_size = size();
    const int mpi_rank = rank();
    validate_rank_grid(dd, mpi_rank, mpi_size);

    const auto periodic = dd.periodic();
    const Box box = system.box();

    std::vector<double> local_state;
    local_state.reserve(system.num_local_atoms() * atom_state_record_width);

    const auto masses = system.masses();
    const auto charges = system.charges();
    const auto atom_types = system.atom_types();
    const auto molecule_ids = system.molecule_ids();
    const auto atomic_numbers = system.atomic_numbers();
    const auto coordinates = system.coordinates();
    const auto velocities = system.velocities();
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        local_state.push_back(masses[atom_index]);
        local_state.push_back(charges[atom_index]);
        local_state.push_back(static_cast<double>(atom_types[atom_index]));
        local_state.push_back(static_cast<double>(molecule_ids[atom_index]));
        local_state.push_back(static_cast<double>(atomic_numbers[atom_index]));
        local_state.push_back(static_cast<double>(system.atom_tag(atom_index)));
        local_state.push_back(coordinates[atom_index][0]);
        local_state.push_back(coordinates[atom_index][1]);
        local_state.push_back(coordinates[atom_index][2]);
        local_state.push_back(velocities[atom_index][0]);
        local_state.push_back(velocities[atom_index][1]);
        local_state.push_back(velocities[atom_index][2]);
    }

    const int send_count = static_cast<int>(local_state.size());
    std::vector<int> recv_counts(static_cast<std::size_t>(mpi_size), 0);
    MPI_Allgather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> displs(static_cast<std::size_t>(mpi_size), 0);
    int total_count = 0;
    for (int index = 0; index < mpi_size; ++index) {
        displs[static_cast<std::size_t>(index)] = total_count;
        total_count += recv_counts[static_cast<std::size_t>(index)];
    }
    if (total_count < 0 || total_count % atom_state_record_width != 0) {
        throw std::runtime_error("MPI atom redistribution received a malformed atom-state buffer");
    }

    std::vector<double> global_state(static_cast<std::size_t>(total_count));
    MPI_Allgatherv(local_state.data(),
                   send_count,
                   MPI_DOUBLE,
                   global_state.data(),
                   recv_counts.data(),
                   displs.data(),
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);

    struct AtomState {
        double mass;
        double charge;
        int atom_type;
        int molecule_id;
        int atomic_number;
        int tag;
        System::Vec3 coordinate;
        System::Vec3 velocity;
    };

    std::vector<AtomState> local_atoms;
    local_atoms.reserve(global_state.size() / atom_state_record_width);
    for (std::size_t offset = 0; offset < global_state.size(); offset += atom_state_record_width) {
        AtomState atom{
            .mass = global_state[offset],
            .charge = global_state[offset + 1],
            .atom_type = static_cast<int>(global_state[offset + 2]),
            .molecule_id = static_cast<int>(global_state[offset + 3]),
            .atomic_number = static_cast<int>(global_state[offset + 4]),
            .tag = static_cast<int>(global_state[offset + 5]),
            .coordinate = {
                global_state[offset + 6],
                global_state[offset + 7],
                global_state[offset + 8],
            },
            .velocity = {
                global_state[offset + 9],
                global_state[offset + 10],
                global_state[offset + 11],
            },
        };

        for (std::size_t dim = 0; dim < periodic.size(); ++dim) {
            if (!periodic[dim]) {
                continue;
            }

            atom.coordinate[dim] = std::fmod(atom.coordinate[dim], box.lengths[dim]);
            if (atom.coordinate[dim] < 0.0) {
                atom.coordinate[dim] += box.lengths[dim];
            }
        }

        if (dd.owner_rank(system.box(), atom.coordinate) != mpi_rank) {
            continue;
        }
        local_atoms.push_back(atom);
    }

    system.resize(local_atoms.size(), local_atoms.size());
    system.set_box(box);

    auto local_masses = system.mutable_masses();
    auto local_charges = system.mutable_charges();
    auto local_atom_types = system.mutable_atom_types();
    auto local_molecule_ids = system.mutable_molecule_ids();
    auto local_atomic_numbers = system.mutable_atomic_numbers();
    auto local_coordinates = system.mutable_coordinates();
    auto local_velocities = system.mutable_velocities();
    auto local_tags = system.mutable_atom_tags();
    auto local_owners = system.mutable_atom_owners();
    for (std::size_t atom_index = 0; atom_index < local_atoms.size(); ++atom_index) {
        const AtomState& atom = local_atoms[atom_index];
        local_masses[atom_index] = atom.mass;
        local_charges[atom_index] = atom.charge;
        local_atom_types[atom_index] = atom.atom_type;
        local_molecule_ids[atom_index] = atom.molecule_id;
        local_atomic_numbers[atom_index] = atom.atomic_number;
        local_coordinates[atom_index] = atom.coordinate;
        local_velocities[atom_index] = atom.velocity;
        local_tags[atom_index] = atom.tag;
        local_owners[atom_index] = mpi_rank;
    }
#else
    (void)system;
    (void)dd;
#endif
}

std::vector<double> MpiCommunicator::pack_send_buffer(
        const System& system,
        const DomainDecomposition& dd,
        const std::array<int, 3>& offset) const {
    std::vector<double> send_buffer;
    send_buffer.reserve(system.num_local_atoms() * ghost_record_width);

    const DomainInfo& domain = dd.info();
    const double ghost_width = dd.ghost_width();
    const auto coordinates = system.coordinates();
    const auto masses = system.masses();
    const auto charges = system.charges();
    const auto atom_types = system.atom_types();
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        const auto& position = coordinates[atom_index];
        bool in_halo = true;
        for (std::size_t dim = 0; dim < offset.size(); ++dim) {
            if (offset[dim] < 0) {
                in_halo = in_halo &&
                    position[dim] >= domain.owned_lo[dim] &&
                    position[dim] < domain.owned_lo[dim] + ghost_width;
            } else if (offset[dim] > 0) {
                in_halo = in_halo &&
                    position[dim] >= domain.owned_hi[dim] - ghost_width &&
                    position[dim] < domain.owned_hi[dim];
            }
        }
        if (!in_halo) {
            continue;
        }

        send_buffer.push_back(position[0]);
        send_buffer.push_back(position[1]);
        send_buffer.push_back(position[2]);
        send_buffer.push_back(masses[atom_index]);
        send_buffer.push_back(charges[atom_index]);
        send_buffer.push_back(static_cast<double>(system.atom_tag(atom_index)));
        send_buffer.push_back(static_cast<double>(atom_types[atom_index]));
    }

    return send_buffer;
}

void MpiCommunicator::unpack_recv_buffer(System& system,
                                         const std::vector<double>& recv_buffer,
                                         int owner) const {
    if (recv_buffer.size() % ghost_record_width != 0) {
        throw std::invalid_argument("Ghost receive buffer does not contain complete atoms");
    }

    for (std::size_t offset = 0; offset < recv_buffer.size(); offset += ghost_record_width) {
        const int tag = static_cast<int>(recv_buffer[offset + 5]);
        bool already_present = false;
        for (std::size_t atom_index = system.num_local_atoms();
             atom_index < system.atom_count();
             ++atom_index) {
            if (system.atom_tag(atom_index) == tag &&
                system.atom_owner(atom_index) == owner) {
                already_present = true;
                break;
            }
        }
        if (already_present) {
            continue;
        }

        const System::Vec3 position{
            recv_buffer[offset],
            recv_buffer[offset + 1],
            recv_buffer[offset + 2]
        };
        system.add_ghost_atom(recv_buffer[offset + 3],
                              recv_buffer[offset + 4],
                              position,
                              tag,
                              owner);
        system.mutable_atom_types()[system.atom_count() - 1] =
            static_cast<int>(recv_buffer[offset + 6]);
    }
}

std::vector<double> MpiCommunicator::pack_reverse_force_buffer(const System& system,
                                                               int home_rank) const {
    std::vector<double> send_buffer;
    send_buffer.reserve(system.num_ghost_atoms() * reverse_force_record_width);

    const auto forces = system.forces();
    for (std::size_t atom_index = system.num_local_atoms();
         atom_index < system.atom_count();
         ++atom_index) {
        if (system.atom_owner(atom_index) != home_rank) {
            continue;
        }

        send_buffer.push_back(static_cast<double>(system.atom_tag(atom_index)));
        send_buffer.push_back(forces[atom_index][0]);
        send_buffer.push_back(forces[atom_index][1]);
        send_buffer.push_back(forces[atom_index][2]);

        // Pair virials are reduced at the force-provider level today. Keep the
        // reverse packet layout ready for providers that attach ghost virials.
        for (std::size_t index = 0; index < 9; ++index) {
            send_buffer.push_back(0.0);
        }
    }

    return send_buffer;
}

void MpiCommunicator::unpack_reverse_force_buffer(System& system,
                                                  const std::vector<double>& recv_buffer) const {
    if (recv_buffer.size() % reverse_force_record_width != 0) {
        throw std::invalid_argument("Reverse force buffer does not contain complete atoms");
    }

    auto forces = system.mutable_forces();
    for (std::size_t offset = 0; offset < recv_buffer.size(); offset += reverse_force_record_width) {
        const int tag = static_cast<int>(recv_buffer[offset]);
        std::size_t local_index = system.num_local_atoms();
        for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
            if (system.atom_tag(atom_index) == tag) {
                local_index = atom_index;
                break;
            }
        }
        if (local_index == system.num_local_atoms()) {
            throw std::runtime_error("Reverse force exchange received an unknown local atom tag");
        }

        forces[local_index][0] += recv_buffer[offset + 1];
        forces[local_index][1] += recv_buffer[offset + 2];
        forces[local_index][2] += recv_buffer[offset + 3];

        std::array<double, 9> virial{};
        for (std::size_t index = 0; index < virial.size(); ++index) {
            virial[index] = recv_buffer[offset + 4 + index];
        }
        system.accumulate_reverse_ghost_virial(virial);
    }
}

int MpiCommunicator::rank() const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        int mpi_rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        return mpi_rank;
    }
#endif

    return 0;
}

int MpiCommunicator::size() const {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        int mpi_size = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        return mpi_size;
    }
#endif

    return 1;
}

}  // namespace gmd
