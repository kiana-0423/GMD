#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace gmd_next::reference {

// Host-only mathematical oracle, not the future DeviceState/view ABI.
using Vec3 = std::array<double, 3>;
using Tensor3 = std::array<double, 9>;  // row-major

// Flattened ordered interaction slots. Repeated indices are intentional.
std::vector<Vec3> gather(std::span<const Vec3> values,
                         std::span<const std::size_t> indices);
// Returns a fresh zero-initialized array, then sums all contributions.
std::vector<Vec3> assemble_sum(std::size_t atom_count,
                              std::span<const std::size_t> indices,
                              std::span<const Vec3> contributions);

struct Pair {
    std::size_t source;
    std::size_t target;
    double scale = 1.0;
};

enum class CutoffMode { none, potential_shift, force_shift };

struct LennardJones {
    double epsilon;
    double sigma;
    CutoffMode cutoff_mode = CutoffMode::none;
    double cutoff = 0.0;  // must be zero for none; otherwise finite and > 0
};

struct OrthorhombicCell {
    Vec3 lengths;
    std::array<bool, 3> periodic{true, true, true};
};

struct PairResult {
    double energy = 0.0;
    std::vector<Vec3> forces;
    Tensor3 virial{};
};

// Half-list evaluator: exactly one entry per unordered pair; either direction.
// The caller owns completeness of the explicit interaction set. No neighbor
// builder, MPI, time integrator, type mixing, or implicit host/device transfers.
// Displacement is target - source; source receives +g, target receives -g.
PairResult evaluate_lj(std::span<const Vec3> positions,
                       std::span<const Pair> pairs,
                       const LennardJones& model,
                       std::optional<OrthorhombicCell> cell = std::nullopt);

}  // namespace gmd_next::reference
