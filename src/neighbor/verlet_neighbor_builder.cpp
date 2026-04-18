#include "gmd/neighbor/verlet_neighbor_builder.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "gmd/boundary/minimum_image.hpp"
#include "gmd/system/system.hpp"

namespace gmd {

VerletNeighborBuilder::VerletNeighborBuilder(double r_cut, double r_skin) noexcept
    : r_cut_(r_cut), r_skin_(r_skin) {}

std::string_view VerletNeighborBuilder::name() const noexcept {
    return "verlet_neighbor_builder";
}

// ---------------------------------------------------------------------------
// Helper: flat cell index from (cx, cy, cz)
// ---------------------------------------------------------------------------
int VerletNeighborBuilder::cell_index(int cx, int cy, int cz) const noexcept {
    return cx + num_cells_[0] * (cy + num_cells_[1] * cz);
}

// ---------------------------------------------------------------------------
// Helper: map one coordinate component → cell index (PBC)
// ---------------------------------------------------------------------------
int VerletNeighborBuilder::coord_to_cell(double x, double cell_sz, int nc) const noexcept {
    int c = static_cast<int>(std::floor(x / cell_sz));
    // Clamp with wrap-around for atoms exactly on the boundary.
    if (c < 0)   c += nc;
    if (c >= nc) c -= nc;
    return c;
}

// ---------------------------------------------------------------------------
// initialize: first rebuild on the initial coordinates.
// ---------------------------------------------------------------------------
void VerletNeighborBuilder::initialize(System& system, RuntimeContext& runtime) {
    rebuild(system, runtime, nullptr);
}

// ---------------------------------------------------------------------------
// needs_rebuild: maximum-displacement criterion.
// Rebuild when any atom has moved more than r_skin/2 since last rebuild.
// ---------------------------------------------------------------------------
bool VerletNeighborBuilder::needs_rebuild(const System& system, std::uint64_t /*step*/) const {
    const auto& nl = system.neighbor_list();
    if (!nl.valid) return true;

    const double trigger_sq = (r_skin_ * 0.5) * (r_skin_ * 0.5);
    const auto coords = system.coordinates();
    const Box& box = system.box();
    const std::size_t n = coords.size();

    for (std::size_t i = 0; i < n; ++i) {
        std::array<double, 3> dr = {
            coords[i][0] - nl.ref_coordinates[i][0],
            coords[i][1] - nl.ref_coordinates[i][1],
            coords[i][2] - nl.ref_coordinates[i][2]
        };
        apply_minimum_image(dr, box);
        const double dx = dr[0];
        const double dy = dr[1];
        const double dz = dr[2];
        if (dx*dx + dy*dy + dz*dz > trigger_sq) {
            return true;
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
// rebuild: cell-list O(N) construction followed by half-pair Verlet list.
// ---------------------------------------------------------------------------
void VerletNeighborBuilder::rebuild(System& system,
                                    RuntimeContext& /*runtime*/,
                                    NeighborBuildStats* stats) {
    const auto coords  = system.coordinates();
    const std::size_t n = coords.size();
    const Box& box     = system.box();
    NeighborList& nl   = system.mutable_neighbor_list();

    const double r_list_sq = r_list() * r_list();

    // -----------------------------------------------------------------------
    // 1. Compute cell grid dimensions.
    //    Each cell side must be >= r_list so that only the 3x3x3 shell of
    //    neighbouring cells needs to be checked.
    // -----------------------------------------------------------------------
    for (int d = 0; d < 3; ++d) {
        num_cells_[d] = std::max(1, static_cast<int>(std::floor(box.lengths[d] / r_list())));
        cell_size_[d] = box.lengths[d] / num_cells_[d];
    }
    const int total_cells = num_cells_[0] * num_cells_[1] * num_cells_[2];

    // -----------------------------------------------------------------------
    // 2. Build linked-list cell structure.
    // -----------------------------------------------------------------------
    head_.assign(total_cells, -1);
    next_.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        // Wrap coordinate into [0, L) before converting to cell index.
        double xi = coords[i][0] - std::floor(coords[i][0] / box.lengths[0]) * box.lengths[0];
        double yi = coords[i][1] - std::floor(coords[i][1] / box.lengths[1]) * box.lengths[1];
        double zi = coords[i][2] - std::floor(coords[i][2] / box.lengths[2]) * box.lengths[2];

        const int cx = coord_to_cell(xi, cell_size_[0], num_cells_[0]);
        const int cy = coord_to_cell(yi, cell_size_[1], num_cells_[1]);
        const int cz = coord_to_cell(zi, cell_size_[2], num_cells_[2]);
        const int cidx = cell_index(cx, cy, cz);

        next_[i]     = head_[cidx];
        head_[cidx]  = static_cast<int>(i);
    }

    // -----------------------------------------------------------------------
    // 3. Build half-pair neighbor list using the cell structure.
    //    For each atom i, only visit atoms j in the 3x3x3 neighbourhood of
    //    i's cell AND with j > i (Newton III half-pairs).
    // -----------------------------------------------------------------------
    nl.counts.assign(n, 0);
    nl.offsets.assign(n + 1, 0);
    nl.neighbors.clear();
    nl.ref_coordinates.resize(n);

    // Temporary per-atom neighbor vectors to avoid overallocation.
    std::vector<std::vector<int>> tmp(n);

    for (std::size_t i = 0; i < n; ++i) {
        double xi = coords[i][0] - std::floor(coords[i][0] / box.lengths[0]) * box.lengths[0];
        double yi = coords[i][1] - std::floor(coords[i][1] / box.lengths[1]) * box.lengths[1];
        double zi = coords[i][2] - std::floor(coords[i][2] / box.lengths[2]) * box.lengths[2];

        const int cx = coord_to_cell(xi, cell_size_[0], num_cells_[0]);
        const int cy = coord_to_cell(yi, cell_size_[1], num_cells_[1]);
        const int cz = coord_to_cell(zi, cell_size_[2], num_cells_[2]);

        // Visit all 27 neighbouring cells (including own cell).
        for (int ddz = -1; ddz <= 1; ++ddz) {
            for (int ddy = -1; ddy <= 1; ++ddy) {
                for (int ddx = -1; ddx <= 1; ++ddx) {
                    // PBC wrap cell indices.
                    int nx = (cx + ddx + num_cells_[0]) % num_cells_[0];
                    int ny = (cy + ddy + num_cells_[1]) % num_cells_[1];
                    int nz = (cz + ddz + num_cells_[2]) % num_cells_[2];
                    int ncidx = cell_index(nx, ny, nz);

                    // Walk the linked list of the neighbouring cell.
                    int j = head_[ncidx];
                    while (j != -1) {
                        if (static_cast<std::size_t>(j) > i) {
                            // Half-pair: only store (i, j) with j > i.
                            std::array<double, 3> dr = {
                                coords[i][0] - coords[j][0],
                                coords[i][1] - coords[j][1],
                                coords[i][2] - coords[j][2]
                            };
                            apply_minimum_image(dr, box);
                            const double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                            if (r2 < r_list_sq) {
                                tmp[i].push_back(j);
                            }
                        }
                        j = next_[j];
                    }
                }
            }
        }
        nl.counts[i] = static_cast<int>(tmp[i].size());
    }

    // -----------------------------------------------------------------------
    // 4. Compact into CSR flat storage.
    // -----------------------------------------------------------------------
    nl.offsets[0] = 0;
    for (std::size_t i = 0; i < n; ++i) {
        nl.offsets[i + 1] = nl.offsets[i] + nl.counts[i];
    }
    nl.neighbors.resize(static_cast<std::size_t>(nl.offsets[n]));
    for (std::size_t i = 0; i < n; ++i) {
        for (int k = 0; k < nl.counts[i]; ++k) {
            nl.neighbors[nl.offsets[i] + k] = tmp[i][k];
        }
    }

    // -----------------------------------------------------------------------
    // 5. Save reference coordinates for next needs_rebuild() check.
    // -----------------------------------------------------------------------
    for (std::size_t i = 0; i < n; ++i) {
        nl.ref_coordinates[i] = coords[i];
    }
    nl.valid = true;

    if (stats != nullptr) {
        stats->rebuilt    = true;
        stats->pair_count = static_cast<std::uint64_t>(nl.neighbors.size());
    }
}

}  // namespace gmd
