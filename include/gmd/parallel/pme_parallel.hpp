#pragma once

#include <array>
#include <complex>
#include <vector>

namespace gmd {

// Process-grid metadata for a PME Py x Pz pencil decomposition.
//
// Phase 4 currently keeps a replicated PME mesh on every rank. The range
// metadata below still describes the x-pencil each rank will own once PME uses
// distributed local FFTs and MPI transpose collectives.
class PmeParallelDecomposition {
public:
    void setup(int nprocs, const std::array<int, 3>& grid);

    std::array<int, 3> local_grid_start() const noexcept;
    std::array<int, 3> local_grid_size() const noexcept;

    void transpose_x_to_y(std::vector<std::complex<double>>& data) const;
    void transpose_y_to_z(std::vector<std::complex<double>>& data) const;
    void transpose_z_to_y(std::vector<std::complex<double>>& data) const;
    void transpose_y_to_x(std::vector<std::complex<double>>& data) const;

private:
    std::array<int, 3> grid_ = {0, 0, 0};
    std::array<int, 3> local_start_ = {0, 0, 0};
    std::array<int, 3> local_size_ = {0, 0, 0};
    std::array<int, 2> proc_grid_ = {1, 1};
    std::array<int, 2> proc_coord_ = {0, 0};

    void validate_replicated_mesh(const std::vector<std::complex<double>>& data) const;
};

}  // namespace gmd
