#pragma once

#include <array>
#include <cstdint>
#include <string_view>
#include <vector>

#include "gmd/neighbor/neighbor_builder.hpp"

namespace gmd {

// Verlet neighbor list built with a cell-list algorithm.
//
// Strategy
// --------
//   rebuild()  — O(N):  partition atoms into cells of side >= r_list,
//                        then only check atoms in the 3x3x3 neighbourhood
//                        of each cell's cube.  Half-pair storage (j > i)
//                        with Newton's-third-law symmetry.
//   needs_rebuild() — O(N): maximum-displacement criterion.
//                        Triggers rebuild when any atom has moved more
//                        than r_skin/2 since the last rebuild.
//
// Parameters (all in Angstrom)
//   r_cut   — force cutoff radius (must match ForceProvider)
//   r_skin  — extra "skin" distance so lists remain valid for several steps
//              typical value: 1.0–2.0 Å
//
// The resulting NeighborList stored in System uses CSR format:
//   atom i's half-pair neighbors: system.neighbor_list().neighbors[
//       offsets[i] .. offsets[i] + counts[i] )
class VerletNeighborBuilder final : public NeighborBuilder {
public:
    // r_cut [Å]: force cutoff;  r_skin [Å]: Verlet skin thickness.
    explicit VerletNeighborBuilder(double r_cut, double r_skin = 1.5) noexcept;
    ~VerletNeighborBuilder() override = default;

    std::string_view name() const noexcept override;

    void initialize(System& system, RuntimeContext& runtime) override;
    bool needs_rebuild(const System& system, std::uint64_t step) const override;
    void rebuild(System& system, RuntimeContext& runtime, NeighborBuildStats* stats) override;

    double r_cut()  const noexcept { return r_cut_; }
    double r_skin() const noexcept { return r_skin_; }
    double r_list() const noexcept { return r_cut_ + r_skin_; }

private:
    double r_cut_;
    double r_skin_;

    // --- cell list internal storage ---
    // Number of cells along each axis.
    std::array<int, 3> num_cells_{1, 1, 1};
    // cell_size_[d] = box_length[d] / num_cells_[d]
    std::array<double, 3> cell_size_{0.0, 0.0, 0.0};
    // Head-of-chain array: head_[cell_idx] = first atom index in cell (-1 if empty).
    std::vector<int> head_;
    // Linked-list: next_[i] = next atom in the same cell (-1 if last).
    std::vector<int> next_;

    // Compute the flat cell index for a fractional position.
    int cell_index(int cx, int cy, int cz) const noexcept;
    // Map a coordinate component to a cell index (with PBC wrap).
    int coord_to_cell(double x, double cell_size, int num_cells) const noexcept;
};

}  // namespace gmd
