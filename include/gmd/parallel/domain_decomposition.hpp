#pragma once

#include <array>

namespace gmd {

struct Box;

struct DomainInfo {
    std::array<double, 3> lo{0.0, 0.0, 0.0};
    std::array<double, 3> hi{0.0, 0.0, 0.0};
    std::array<int, 3> proc_grid{1, 1, 1};
    std::array<int, 3> proc_coord{0, 0, 0};
};

class DomainDecomposition {
public:
    // Split the box uniformly along x using the currently configured ghost width.
    void create_1d_decomposition(const Box& box, int nprocs, int my_rank);

    // Configure the ghost width from the neighbor-list radius before splitting.
    void create_1d_decomposition(const Box& box,
                                 int nprocs,
                                 int my_rank,
                                 double cutoff,
                                 double skin);

    bool is_local(const std::array<double, 3>& pos) const;
    const DomainInfo& info() const;
    double ghost_width() const noexcept;

private:
    DomainInfo info_{};
    double ghost_width_{0.0};
};

}  // namespace gmd
