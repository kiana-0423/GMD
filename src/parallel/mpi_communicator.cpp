#include "gmd/parallel/mpi_communicator.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

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
constexpr int atom_state_record_width = 11;

#ifdef GMD_ENABLE_MPI
void validate_1d_rank_grid(const DomainDecomposition& dd, int mpi_rank, int mpi_size) {
    const DomainInfo& domain = dd.info();
    const int nprocs = domain.proc_grid[0];
    const int my_x_coord = domain.proc_coord[0];
    if (nprocs != mpi_size ||
        my_x_coord < 0 ||
        my_x_coord >= nprocs ||
        mpi_rank != my_x_coord) {
        throw std::invalid_argument("MPI exchange requires a 1D rank-ordered domain grid");
    }
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
    validate_1d_rank_grid(dd, my_rank, size());

    const DomainInfo& domain = dd.info();
    const int nprocs = domain.proc_grid[0];
    const int my_x_coord = domain.proc_coord[0];
 const bool periodic = domain.periodic_x;
 const double ghost_width = dd.ghost_width();
 const double domain_width = system.box().lengths[0] / static_cast<double>(nprocs);
 const double owned_lo_x = domain_width * static_cast<double>(my_x_coord);
 const double owned_hi_x = my_x_coord == nprocs - 1
 ? system.box().lengths[0]
 : domain_width * static_cast<double>(my_x_coord + 1);

 auto exchange_boundary = [&](int neighbor,
 double send_lo_x,
 double send_hi_x,
 double shift_x) {
 std::vector<double> send_buffer =
 pack_send_buffer(system, send_lo_x, send_hi_x);

 // Apply periodic shift to outgoing ghost positions so that the
 // receiving rank sees them at the correct (nearby) coordinate.
 if (std::abs(shift_x) > 0.0) {
 for (std::size_t i = 0; i < send_buffer.size(); i += ghost_record_width) {
 send_buffer[i] += shift_x; // x-coordinate is the first field
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
 neighbor,
 110,
 &recv_count,
 1,
 MPI_INT,
 neighbor,
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
 neighbor,
 111,
 recv_buffer.data(),
 recv_count,
 MPI_DOUBLE,
 neighbor,
 111,
 MPI_COMM_WORLD,
 MPI_STATUS_IGNORE);
 unpack_recv_buffer(system, recv_buffer, neighbor);
 };

 // Left neighbor (rank-1). With periodic x, rank 0's left neighbor is rank nprocs-1.
 const bool has_left = (my_x_coord > 0) || periodic;
 if (has_left) {
 const int left_neighbor = (my_x_coord > 0) ? my_rank - 1 : nprocs - 1;
 // Send atoms in [owned_lo_x, owned_lo_x + ghost_width) to the left neighbor.
 // Those atoms become the left neighbor's right-side ghosts.
 // When crossing the periodic boundary (rank 0 -> rank nprocs-1),
 // shift positions by -Lx so they appear near the right edge of the box.
 const double shift = (my_x_coord == 0) ? -system.box().lengths[0] : 0.0;
 exchange_boundary(left_neighbor, owned_lo_x, owned_lo_x + ghost_width, shift);
 }

 // Right neighbor (rank+1). With periodic x, rank nprocs-1's right neighbor is rank 0.
 const bool has_right = (my_x_coord + 1 < nprocs) || periodic;
 if (has_right) {
 const int right_neighbor = (my_x_coord + 1 < nprocs) ? my_rank + 1 : 0;
 const double shift = (my_x_coord == nprocs - 1) ? system.box().lengths[0] : 0.0;
 exchange_boundary(right_neighbor, owned_hi_x - ghost_width, owned_hi_x, shift);
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
 validate_1d_rank_grid(dd, my_rank, size());

 const DomainInfo& domain = dd.info();
 const int nprocs = domain.proc_grid[0];
 const int my_x_coord = domain.proc_coord[0];
 const bool periodic = domain.periodic_x;

 auto exchange_forces = [&](int neighbor) {
 const std::vector<double> send_buffer = pack_reverse_force_buffer(system, neighbor);
 if (send_buffer.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
 throw std::overflow_error("Reverse force exchange buffer is too large for MPI_Sendrecv");
 }

 const int send_count = static_cast<int>(send_buffer.size());
 int recv_count = 0;
 MPI_Sendrecv(&send_count,
 1,
 MPI_INT,
 neighbor,
 120,
 &recv_count,
 1,
 MPI_INT,
 neighbor,
 120,
 MPI_COMM_WORLD,
 MPI_STATUS_IGNORE);
 if (recv_count < 0 || recv_count % reverse_force_record_width != 0) {
 throw std::runtime_error("Reverse force exchange received a malformed buffer size");
 }

 std::vector<double> recv_buffer(static_cast<std::size_t>(recv_count));
 MPI_Sendrecv(send_buffer.data(),
 send_count,
 MPI_DOUBLE,
 neighbor,
 121,
 recv_buffer.data(),
 recv_count,
 MPI_DOUBLE,
 neighbor,
 121,
 MPI_COMM_WORLD,
 MPI_STATUS_IGNORE);
 unpack_reverse_force_buffer(system, recv_buffer);
 };

 // Left neighbor. With periodic x, rank 0 sends reverse forces to rank nprocs-1.
 const bool has_left = (my_x_coord > 0) || periodic;
 if (has_left) {
 const int left_neighbor = (my_x_coord > 0) ? my_rank - 1 : nprocs - 1;
 exchange_forces(left_neighbor);
 }

 // Right neighbor. With periodic x, rank nprocs-1 sends reverse forces to rank 0.
 const bool has_right = (my_x_coord + 1 < nprocs) || periodic;
 if (has_right) {
 const int right_neighbor = (my_x_coord + 1 < nprocs) ? my_rank + 1 : 0;
 exchange_forces(right_neighbor);
 }
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
 validate_1d_rank_grid(dd, mpi_rank, mpi_size);

 const DomainInfo& domain = dd.info();
 const int nprocs = domain.proc_grid[0];
 const int my_x_coord = domain.proc_coord[0];
 const bool periodic = domain.periodic_x;
 const double Lx = system.box().lengths[0];
 const double domain_width = Lx / static_cast<double>(nprocs);

 std::vector<double> local_state;
 local_state.reserve(system.num_local_atoms() * atom_state_record_width);

 const auto masses = system.masses();
 const auto charges = system.charges();
 const auto atom_types = system.atom_types();
 const auto atomic_numbers = system.atomic_numbers();
 const auto coordinates = system.coordinates();
 const auto velocities = system.velocities();
 for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
 local_state.push_back(masses[atom_index]);
 local_state.push_back(charges[atom_index]);
 local_state.push_back(static_cast<double>(atom_types[atom_index]));
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
 .atomic_number = static_cast<int>(global_state[offset + 3]),
 .tag = static_cast<int>(global_state[offset + 4]),
 .coordinate = {
 global_state[offset + 5],
 global_state[offset + 6],
 global_state[offset + 7],
 },
 .velocity = {
 global_state[offset + 8],
 global_state[offset + 9],
 global_state[offset + 10],
 },
 };

 // Apply periodic x wrapping before determining the owner rank.
 if (periodic) {
 if (atom.coordinate[0] < 0.0) {
 atom.coordinate[0] += Lx;
 } else if (atom.coordinate[0] >= Lx) {
 atom.coordinate[0] -= Lx;
 }
 }

        if (dd.owner_rank(system.box(), atom.coordinate) != mpi_rank) {
            continue;
        }
        local_atoms.push_back(atom);
    }

    const Box box = system.box();
    system.resize(local_atoms.size(), local_atoms.size());
    system.set_box(box);

    auto local_masses = system.mutable_masses();
    auto local_charges = system.mutable_charges();
    auto local_atom_types = system.mutable_atom_types();
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

std::vector<double> MpiCommunicator::pack_send_buffer(const System& system,
                                                      double boundary_lo_x,
                                                      double boundary_hi_x) const {
    std::vector<double> send_buffer;
    send_buffer.reserve(system.num_local_atoms() * ghost_record_width);

    const auto coordinates = system.coordinates();
    const auto masses = system.masses();
    const auto charges = system.charges();
    const auto atom_types = system.atom_types();
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        const auto& position = coordinates[atom_index];
        if (position[0] < boundary_lo_x || position[0] >= boundary_hi_x) {
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
        const System::Vec3 position{
            recv_buffer[offset],
            recv_buffer[offset + 1],
            recv_buffer[offset + 2]
        };
        system.add_ghost_atom(recv_buffer[offset + 3],
                              recv_buffer[offset + 4],
                              position,
                              static_cast<int>(recv_buffer[offset + 5]),
                              owner);
        system.mutable_atom_types()[system.atom_count() - 1] =
            static_cast<int>(recv_buffer[offset + 6]);
    }
}

std::vector<double> MpiCommunicator::pack_reverse_force_buffer(const System& system,
                                                               int home_rank) const {
    std::vector<double> send_buffer;
    send_buffer.reserve(system.num_ghost_atoms() * reverse_force_record_width);

    const auto coordinates = system.coordinates();
    const auto forces = system.forces();
    for (std::size_t atom_index = system.num_local_atoms();
         atom_index < system.atom_count();
         ++atom_index) {
        if (system.atom_owner(atom_index) != home_rank) {
            continue;
        }

        const auto& position = coordinates[atom_index];
        const auto& force = forces[atom_index];
        send_buffer.push_back(static_cast<double>(system.atom_tag(atom_index)));
        send_buffer.push_back(force[0]);
        send_buffer.push_back(force[1]);
        send_buffer.push_back(force[2]);
        for (double coordinate : position) {
            for (double component : force) {
                send_buffer.push_back(coordinate * component);
            }
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
