#pragma once

#include <array>

namespace gmd {

struct Box;

struct DomainInfo {
    // Halo-expanded subdomain bounds.
    std::array<double, 3> lo{0.0, 0.0, 0.0};
    std::array<double, 3> hi{0.0, 0.0, 0.0};
    // Ownership bounds without ghost width.
    std::array<double, 3> owned_lo{0.0, 0.0, 0.0};
    std::array<double, 3> owned_hi{0.0, 0.0, 0.0};
    std::array<int, 3> proc_grid{1, 1, 1};
    std::array<int, 3> proc_coord{0, 0, 0};
    std::array<bool, 3> periodic{true, true, true};
};

class DomainDecomposition {
public:
    static constexpr int no_rank = -1;

    static std::array<int, 3> choose_processor_grid(int nprocs);

    // Split the box uniformly over a balanced 3D process grid.
    void create_decomposition(const Box& box, int nprocs, int my_rank);

    void create_decomposition(const Box& box,
                              int nprocs,
                              int my_rank,
                              double cutoff,
                              double skin,
                              const std::array<bool, 3>& periodic);

    // Split over an explicitly supplied 3D process grid.
    void create_decomposition(const Box& box,
                              const std::array<int, 3>& proc_grid,
                              int my_rank,
                              double cutoff,
                              double skin,
                              const std::array<bool, 3>& periodic);

    // Split the box uniformly along x using the currently configured ghost width.
    void create_1d_decomposition(const Box& box, int nprocs, int my_rank);

    // Configure the ghost width from the neighbor-list radius before splitting.
    void create_1d_decomposition(const Box& box,
                                 int nprocs,
                                 int my_rank,
                                 double cutoff,
                                 double skin);

    // Configure an explicit x-periodicity mode for the 1D rank ring.
    void create_1d_decomposition(const Box& box,
                                 int nprocs,
                                 int my_rank,
                                 double cutoff,
                                 double skin,
                                 bool periodic_x);

    bool is_local(const std::array<double, 3>& pos) const;
    int owner_rank(const Box& box, const std::array<double, 3>& pos) const;
    int rank_from_proc_coord(const std::array<int, 3>& proc_coord) const;
    std::array<int, 3> proc_coord_from_rank(int rank) const;
    int neighbor_rank(const std::array<int, 3>& offset) const;
    void refresh(const Box& box);
    const DomainInfo& info() const;
    double ghost_width() const noexcept;
    bool periodic_x() const noexcept;
    const std::array<bool, 3>& periodic() const noexcept;

private:
    void configure(const Box& box,
                   const std::array<int, 3>& proc_grid,
                   int my_rank);

    DomainInfo info_{};
    double ghost_width_{0.0};
};

}  // namespace gmd
