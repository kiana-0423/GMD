#include "gmd/parallel/domain_decomposition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include "gmd/system/box.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {
namespace {

void validate_box(const Box& box) {
    for (double length : box.lengths) {
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("Domain decomposition requires positive box lengths");
        }
    }
}

int process_count(const std::array<int, 3>& grid) {
    int count = 1;
    for (int extent : grid) {
        if (extent <= 0) {
            throw std::invalid_argument(
                "Domain decomposition process-grid extents must be positive");
        }
        if (count > std::numeric_limits<int>::max() / extent) {
            throw std::invalid_argument("Domain decomposition process grid is too large");
        }
        count *= extent;
    }
    return count;
}

std::array<int, 3> factor_balanced_grid(int nprocs) {
    std::array<int, 3> best{nprocs, 1, 1};
    int best_spread = nprocs - 1;
    int best_largest = nprocs;

    for (int px = 1; px <= nprocs; ++px) {
        if (nprocs % px != 0) {
            continue;
        }

        const int yz_count = nprocs / px;
        for (int py = 1; py <= yz_count; ++py) {
            if (yz_count % py != 0) {
                continue;
            }

            std::array<int, 3> candidate{px, py, yz_count / py};
            std::sort(candidate.begin(), candidate.end(), std::greater<int>{});
            const int spread = candidate[0] - candidate[2];
            if (spread < best_spread ||
                (spread == best_spread && candidate[0] < best_largest)) {
                best = candidate;
                best_spread = spread;
                best_largest = candidate[0];
            }
        }
    }

    return best;
}

#ifdef GMD_ENABLE_MPI
bool mpi_is_available() noexcept {
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    return is_initialized != 0 && is_finalized == 0;
}
#endif

}  // namespace

std::array<int, 3> DomainDecomposition::choose_processor_grid(int nprocs) {
    if (nprocs <= 0) {
        throw std::invalid_argument("Domain decomposition requires at least one process");
    }

#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        std::array<int, 3> dims{0, 0, 0};
        MPI_Dims_create(nprocs, static_cast<int>(dims.size()), dims.data());
        return dims;
    }
#endif

    return factor_balanced_grid(nprocs);
}

void DomainDecomposition::create_decomposition(const Box& box,
                                               int nprocs,
                                               int my_rank) {
    configure(box, choose_processor_grid(nprocs), my_rank);
}

void DomainDecomposition::create_decomposition(const Box& box,
                                               int nprocs,
                                               int my_rank,
                                               double cutoff,
                                               double skin,
                                               const std::array<bool, 3>& periodic) {
    create_decomposition(box,
                         choose_processor_grid(nprocs),
                         my_rank,
                         cutoff,
                         skin,
                         periodic);
}

void DomainDecomposition::create_decomposition(
        const Box& box,
        const std::array<int, 3>& proc_grid,
        int my_rank,
        double cutoff,
        double skin,
        const std::array<bool, 3>& periodic) {
    if (!std::isfinite(cutoff) || !std::isfinite(skin) || cutoff < 0.0 || skin < 0.0) {
        throw std::invalid_argument("Domain decomposition ghost distances must be non-negative");
    }

    ghost_width_ = cutoff + skin;
    info_.periodic = periodic;
    configure(box, proc_grid, my_rank);
}

void DomainDecomposition::configure(const Box& box,
                                    const std::array<int, 3>& proc_grid,
                                    int my_rank) {
    validate_box(box);
    const int nprocs = process_count(proc_grid);
    if (my_rank < 0 || my_rank >= nprocs) {
        throw std::invalid_argument("Domain decomposition rank is outside the process grid");
    }

    info_.proc_grid = proc_grid;
    info_.proc_coord = proc_coord_from_rank(my_rank);
    for (std::size_t dim = 0; dim < info_.proc_grid.size(); ++dim) {
        const double domain_width =
            box.lengths[dim] / static_cast<double>(info_.proc_grid[dim]);
        const int coord = info_.proc_coord[dim];
        info_.owned_lo[dim] = domain_width * static_cast<double>(coord);
        info_.owned_hi[dim] = coord == info_.proc_grid[dim] - 1
            ? box.lengths[dim]
            : domain_width * static_cast<double>(coord + 1);
        info_.lo[dim] = info_.owned_lo[dim] - ghost_width_;
        info_.hi[dim] = info_.owned_hi[dim] + ghost_width_;
    }
}

void DomainDecomposition::create_1d_decomposition(const Box& box,
                                                  int nprocs,
                                                  int my_rank) {
    configure(box, {nprocs, 1, 1}, my_rank);
}

void DomainDecomposition::create_1d_decomposition(const Box& box,
                                                  int nprocs,
                                                  int my_rank,
                                                  double cutoff,
                                                  double skin,
                                                  bool periodic_x) {
    if (!std::isfinite(cutoff) || !std::isfinite(skin) || cutoff < 0.0 || skin < 0.0) {
        throw std::invalid_argument("Domain decomposition ghost distances must be non-negative");
    }

    ghost_width_ = cutoff + skin;
    info_.periodic[0] = periodic_x;
    create_1d_decomposition(box, nprocs, my_rank);
}

void DomainDecomposition::create_1d_decomposition(const Box& box,
                                                  int nprocs,
                                                  int my_rank,
                                                  double cutoff,
                                                  double skin) {
    create_1d_decomposition(box, nprocs, my_rank, cutoff, skin, info_.periodic[0]);
}

bool DomainDecomposition::is_local(const std::array<double, 3>& pos) const {
    for (int d = 0; d < 3; ++d) {
        if (pos[d] < info_.lo[d] || pos[d] >= info_.hi[d]) {
            return false;
        }
    }

    return true;
}

int DomainDecomposition::owner_rank(const Box& box,
                                    const std::array<double, 3>& pos) const {
    validate_box(box);
    process_count(info_.proc_grid);

    std::array<int, 3> owner_coord{};
    for (std::size_t dim = 0; dim < owner_coord.size(); ++dim) {
        double coordinate = pos[dim];
        if (info_.periodic[dim]) {
            coordinate = std::fmod(coordinate, box.lengths[dim]);
            if (coordinate < 0.0) {
                coordinate += box.lengths[dim];
            }
        } else {
            coordinate =
                std::clamp(coordinate, 0.0, std::nextafter(box.lengths[dim], 0.0));
        }

        const double width =
            box.lengths[dim] / static_cast<double>(info_.proc_grid[dim]);
        const int coord = static_cast<int>(coordinate / width);
        owner_coord[dim] = std::min(coord, info_.proc_grid[dim] - 1);
    }

    return rank_from_proc_coord(owner_coord);
}

int DomainDecomposition::rank_from_proc_coord(
        const std::array<int, 3>& proc_coord) const {
    process_count(info_.proc_grid);
    for (std::size_t dim = 0; dim < proc_coord.size(); ++dim) {
        if (proc_coord[dim] < 0 || proc_coord[dim] >= info_.proc_grid[dim]) {
            throw std::invalid_argument("Process coordinate is outside the domain grid");
        }
    }

    return proc_coord[0] +
           info_.proc_grid[0] * (proc_coord[1] + info_.proc_grid[1] * proc_coord[2]);
}

std::array<int, 3> DomainDecomposition::proc_coord_from_rank(int rank) const {
    const int nprocs = process_count(info_.proc_grid);
    if (rank < 0 || rank >= nprocs) {
        throw std::invalid_argument("Rank is outside the domain grid");
    }

    std::array<int, 3> coord{};
    coord[0] = rank % info_.proc_grid[0];
    const int yz_rank = rank / info_.proc_grid[0];
    coord[1] = yz_rank % info_.proc_grid[1];
    coord[2] = yz_rank / info_.proc_grid[1];
    return coord;
}

int DomainDecomposition::neighbor_rank(const std::array<int, 3>& offset) const {
    std::array<int, 3> neighbor = info_.proc_coord;
    for (std::size_t dim = 0; dim < neighbor.size(); ++dim) {
        neighbor[dim] += offset[dim];
        if (neighbor[dim] >= 0 && neighbor[dim] < info_.proc_grid[dim]) {
            continue;
        }
        if (!info_.periodic[dim]) {
            return no_rank;
        }

        neighbor[dim] %= info_.proc_grid[dim];
        if (neighbor[dim] < 0) {
            neighbor[dim] += info_.proc_grid[dim];
        }
    }

    return rank_from_proc_coord(neighbor);
}

void DomainDecomposition::refresh(const Box& box) {
    configure(box, info_.proc_grid, rank_from_proc_coord(info_.proc_coord));
}

const DomainInfo& DomainDecomposition::info() const {
    return info_;
}

double DomainDecomposition::ghost_width() const noexcept {
    return ghost_width_;
}

bool DomainDecomposition::periodic_x() const noexcept {
    return info_.periodic[0];
}

const std::array<bool, 3>& DomainDecomposition::periodic() const noexcept {
    return info_.periodic;
}

}  // namespace gmd
