#include "gmd/parallel/pme_parallel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {
namespace {

int current_rank() noexcept {
#ifdef GMD_ENABLE_MPI
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    if (is_initialized != 0 && is_finalized == 0) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
    }
#endif
    return 0;
}

std::array<int, 2> split_range(int count, int parts, int coord) {
    const int base = count / parts;
    const int remainder = count % parts;
    const int size = base + (coord < remainder ? 1 : 0);
    const int start = coord * base + std::min(coord, remainder);
    return {start, size};
}

}  // namespace

void PmeParallelDecomposition::setup(int nprocs, const std::array<int, 3>& grid) {
    if (nprocs <= 0) {
        throw std::invalid_argument("PME decomposition requires a positive process count");
    }
    for (int extent : grid) {
        if (extent <= 0) {
            throw std::invalid_argument("PME decomposition requires a positive FFT grid");
        }
    }

    grid_ = grid;

    int best_py = 1;
    int best_pz = nprocs;
    int best_balance = std::abs(best_pz - best_py);
    for (int py = 1; py * py <= nprocs; ++py) {
        if (nprocs % py != 0) {
            continue;
        }

        const int pz = nprocs / py;
        const int balance = std::abs(pz - py);
        if (balance < best_balance) {
            best_py = py;
            best_pz = pz;
            best_balance = balance;
        }
    }

    proc_grid_ = {best_py, best_pz};
    const int rank = current_rank();
    if (rank < 0 || rank >= nprocs) {
        throw std::runtime_error("PME rank does not fit the configured process grid");
    }
    proc_coord_ = {rank / proc_grid_[1], rank % proc_grid_[1]};

    const auto y_range = split_range(grid_[1], proc_grid_[0], proc_coord_[0]);
    const auto z_range = split_range(grid_[2], proc_grid_[1], proc_coord_[1]);
    local_start_ = {0, y_range[0], z_range[0]};
    local_size_ = {grid_[0], y_range[1], z_range[1]};
}

std::array<int, 3> PmeParallelDecomposition::local_grid_start() const noexcept {
    return local_start_;
}

std::array<int, 3> PmeParallelDecomposition::local_grid_size() const noexcept {
    return local_size_;
}

void PmeParallelDecomposition::transpose_x_to_y(
        std::vector<std::complex<double>>& data) const {
    validate_replicated_mesh(data);
}

void PmeParallelDecomposition::transpose_y_to_z(
        std::vector<std::complex<double>>& data) const {
    validate_replicated_mesh(data);
}

void PmeParallelDecomposition::transpose_z_to_y(
        std::vector<std::complex<double>>& data) const {
    validate_replicated_mesh(data);
}

void PmeParallelDecomposition::transpose_y_to_x(
        std::vector<std::complex<double>>& data) const {
    validate_replicated_mesh(data);
}

void PmeParallelDecomposition::validate_replicated_mesh(
        const std::vector<std::complex<double>>& data) const {
    const std::size_t expected_size = static_cast<std::size_t>(grid_[0])
                                    * static_cast<std::size_t>(grid_[1])
                                    * static_cast<std::size_t>(grid_[2]);
    if (data.size() != expected_size) {
        throw std::invalid_argument(
            "PME transpose currently expects a replicated full mesh");
    }
}

}  // namespace gmd
