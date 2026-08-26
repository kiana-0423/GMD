#pragma once

#include <array>
#include <vector>

namespace gmd {

struct Box;
class DomainDecomposition;
class System;

class MpiCommunicator {
public:
    double allreduce_scalar(double local_value) const;
    void allreduce_vector(const std::vector<double>& local,
                          std::vector<double>& global) const;

    // Logical OR across every rank. Used for decisions that must be taken
    // identically on all ranks (e.g. whether to rebuild the neighbor list),
    // so that control flow cannot diverge and desynchronise later collectives.
    bool allreduce_logical_or(bool local_value) const;

    void broadcast_box(Box& box, int root) const;
    void barrier() const;

    void exchange_ghost_coordinates(System& system, const DomainDecomposition& dd) const;
    void reverse_accumulate_ghost_forces(System& system, const DomainDecomposition& dd) const;
    void redistribute_atoms(System& system, const DomainDecomposition& dd) const;

    int rank() const;
    int size() const;

private:
    std::vector<double> pack_send_buffer(const System& system,
                                         const DomainDecomposition& dd,
                                         const std::array<int, 3>& offset) const;
    void unpack_recv_buffer(System& system,
                            const std::vector<double>& recv_buffer,
                            int owner) const;
    std::vector<double> pack_reverse_force_buffer(const System& system, int home_rank) const;
    void unpack_reverse_force_buffer(System& system,
                                     const std::vector<double>& recv_buffer) const;
};

}  // namespace gmd
