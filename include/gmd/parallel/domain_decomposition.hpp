#pragma once

#include <array>

namespace gmd {

struct Box;

struct DomainInfo {
    std::array<double, 3> lo{0.0, 0.0, 0.0};
    std::array<double, 3> hi{0.0, 0.0, 0.0};
    std::array<int, 3> proc_grid{1, 1, 1};
    std::array<int, 3> proc_coord{0, 0, 0};
 bool periodic_x{true}; // Whether the x-axis uses periodic boundary conditions.
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

 // Create a 1D decomposition with explicit periodic-x flag.
 void create_1d_decomposition(const Box& box,
 int nprocs,
 int my_rank,
 double cutoff,
 double skin,
 bool periodic_x);

 bool is_local(const std::array<double, 3>& pos) const;
 int owner_rank(const Box& box, const std::array<double, 3>& pos) const;
 void refresh(const Box& box);
 const DomainInfo& info() const;
 double ghost_width() const noexcept;
 bool periodic_x() const noexcept;

}  // namespace gmd
