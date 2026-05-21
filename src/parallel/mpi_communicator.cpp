#include "gmd/parallel/mpi_communicator.hpp"

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
    const double ghost_width = dd.ghost_width();
    const double domain_width = system.box().lengths[0] / static_cast<double>(nprocs);
    const double owned_lo_x = domain_width * static_cast<double>(my_x_coord);
    const double owned_hi_x = my_x_coord == nprocs - 1
        ? system.box().lengths[0]
        : domain_width * static_cast<double>(my_x_coord + 1);

    auto exchange_boundary = [&](int neighbor,
                                 double send_lo_x,
                                 double send_hi_x) {
        const std::vector<double> send_buffer =
            pack_send_buffer(system, send_lo_x, send_hi_x);
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

    if (my_x_coord > 0) {
        exchange_boundary(my_rank - 1, owned_lo_x, owned_lo_x + ghost_width);
    }
    if (my_x_coord + 1 < nprocs) {
        exchange_boundary(my_rank + 1, owned_hi_x - ghost_width, owned_hi_x);
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

    if (my_x_coord > 0) {
        exchange_forces(my_rank - 1);
    }
    if (my_x_coord + 1 < nprocs) {
        exchange_forces(my_rank + 1);
    }
#else
    (void)dd;
#endif

    system.clear_ghost_atoms();
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
