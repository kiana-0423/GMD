#include "gmd/parallel/domain_decomposition.hpp"

#include <cmath>
#include <stdexcept>

#include "gmd/system/box.hpp"

namespace gmd {
namespace {

void validate_box(const Box& box) {
    for (double length : box.lengths) {
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("Domain decomposition requires positive box lengths");
        }
    }
}

}  // namespace

void DomainDecomposition::create_1d_decomposition(const Box& box,
                                                  int nprocs,
                                                  int my_rank) {
    validate_box(box);
    if (nprocs <= 0) {
        throw std::invalid_argument("Domain decomposition requires at least one process");
    }
    if (my_rank < 0 || my_rank >= nprocs) {
        throw std::invalid_argument("Domain decomposition rank is outside the process grid");
    }

    const double domain_width = box.lengths[0] / static_cast<double>(nprocs);
    const double owned_lo_x = domain_width * static_cast<double>(my_rank);
    const double owned_hi_x = my_rank == nprocs - 1
        ? box.lengths[0]
        : domain_width * static_cast<double>(my_rank + 1);

    info_.lo = {owned_lo_x - ghost_width_, 0.0, 0.0};
    info_.hi = {owned_hi_x + ghost_width_, box.lengths[1], box.lengths[2]};
    info_.proc_grid = {nprocs, 1, 1};
    info_.proc_coord = {my_rank, 0, 0};
}

void DomainDecomposition::create_1d_decomposition(const Box& box,
                                                  int nprocs,
                                                  int my_rank,
                                                  double cutoff,
                                                  double skin) {
    if (!std::isfinite(cutoff) || !std::isfinite(skin) || cutoff < 0.0 || skin < 0.0) {
        throw std::invalid_argument("Domain decomposition ghost distances must be non-negative");
    }

    ghost_width_ = cutoff + skin;
    create_1d_decomposition(box, nprocs, my_rank);
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
    const int nprocs = info_.proc_grid[0];
    if (nprocs <= 0) {
        throw std::runtime_error("Domain decomposition is not initialized");
    }

    double x = std::fmod(pos[0], box.lengths[0]);
    if (x < 0.0) {
        x += box.lengths[0];
    }

    const double domain_width = box.lengths[0] / static_cast<double>(nprocs);
    const int owner = static_cast<int>(x / domain_width);
    return std::min(owner, nprocs - 1);
}

void DomainDecomposition::refresh(const Box& box) {
    create_1d_decomposition(box, info_.proc_grid[0], info_.proc_coord[0]);
}

const DomainInfo& DomainDecomposition::info() const {
    return info_;
}

double DomainDecomposition::ghost_width() const noexcept {
    return ghost_width_;
}

}  // namespace gmd
