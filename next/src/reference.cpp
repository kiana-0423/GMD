#include "gmd_next/reference.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

namespace gmd_next::reference {
namespace {

void require_finite(double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Reference input must be finite");
    }
}

void validate_model(const LennardJones& model) {
    require_finite(model.epsilon);
    require_finite(model.sigma);
    require_finite(model.cutoff);
    if (model.epsilon <= 0.0 || model.sigma <= 0.0) {
        throw std::invalid_argument("LJ epsilon and sigma must be positive");
    }
    switch (model.cutoff_mode) {
    case CutoffMode::none:
        if (model.cutoff != 0.0) {
            throw std::invalid_argument("Uncut LJ requires cutoff == 0");
        }
        break;
    case CutoffMode::potential_shift:
    case CutoffMode::force_shift:
        if (model.cutoff <= 0.0) {
            throw std::invalid_argument("Truncated LJ requires a positive cutoff");
        }
        break;
    default:
        throw std::invalid_argument("Unknown LJ cutoff mode");
    }
}

struct RadialValue {
    double energy;
    double derivative;
};

RadialValue radial_lj(double r, const LennardJones& model) {
    const double s = model.sigma / r;
    const double s2 = s * s;
    const double s6 = s2 * s2 * s2;
    const double s12 = s6 * s6;
    const RadialValue result{
        4.0 * model.epsilon * (s12 - s6),
        (24.0 * model.epsilon / r) * (s6 - 2.0 * s12),
    };
    if (!std::isfinite(result.energy) || !std::isfinite(result.derivative)) {
        throw std::overflow_error("LJ radial evaluation overflowed");
    }
    return result;
}

}  // namespace

std::vector<Vec3> gather(std::span<const Vec3> values,
                         std::span<const std::size_t> indices) {
    std::vector<Vec3> result;
    result.reserve(indices.size());
    for (const auto index : indices) {
        if (index >= values.size()) {
            throw std::out_of_range("Gather index exceeds input size");
        }
        result.push_back(values[index]);
    }
    return result;
}

std::vector<Vec3> assemble_sum(std::size_t atom_count,
                              std::span<const std::size_t> indices,
                              std::span<const Vec3> contributions) {
    if (indices.size() != contributions.size()) {
        throw std::invalid_argument("Assembly indices and contributions differ in size");
    }
    std::vector<Vec3> result(atom_count, Vec3{});
    for (std::size_t slot = 0; slot < indices.size(); ++slot) {
        const auto atom = indices[slot];
        if (atom >= atom_count) {
            throw std::out_of_range("Assembly index exceeds atom count");
        }
        for (std::size_t axis = 0; axis < 3; ++axis) {
            result[atom][axis] += contributions[slot][axis];
        }
    }
    return result;
}

PairResult evaluate_lj(std::span<const Vec3> positions,
                       std::span<const Pair> pairs,
                       const LennardJones& model,
                       std::optional<OrthorhombicCell> cell) {
    validate_model(model);
    for (const auto& position : positions) {
        for (const double value : position) {
            require_finite(value);
        }
    }
    if (cell) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto length = cell->lengths[axis];
            require_finite(length);
            if (length <= 0.0) {
                throw std::invalid_argument("Cell lengths must be positive");
            }
            if (cell->periodic[axis] &&
                (model.cutoff_mode == CutoffMode::none || model.cutoff >= 0.5 * length)) {
                throw std::invalid_argument("Periodic reference requires cutoff < half box length");
            }
        }
    }

    std::set<std::pair<std::size_t, std::size_t>> seen;
    for (const auto& pair : pairs) {
        if (pair.source >= positions.size() || pair.target >= positions.size()) {
            throw std::out_of_range("Pair endpoint exceeds position count");
        }
        if (pair.source == pair.target) {
            throw std::invalid_argument("Self pairs are invalid");
        }
        if (!seen.emplace(std::min(pair.source, pair.target),
                          std::max(pair.source, pair.target)).second) {
            throw std::invalid_argument("Half list contains a duplicate unordered pair");
        }
        require_finite(pair.scale);
        if (pair.scale < 0.0) {
            throw std::invalid_argument("Pair scale must be nonnegative");
        }
    }

    PairResult result;
    result.forces.resize(positions.size(), Vec3{});
    const RadialValue at_cutoff = model.cutoff_mode == CutoffMode::none
        ? RadialValue{0.0, 0.0} : radial_lj(model.cutoff, model);
    for (const auto& pair : pairs) {
        // Exclusions must be applied before singular geometry operations.
        if (pair.scale == 0.0) {
            continue;
        }
        Vec3 d{};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            d[axis] = positions[pair.target][axis] - positions[pair.source][axis];
            if (cell && cell->periodic[axis]) {
                // Odd tie rule preserves opposite directions at half-box ties.
                d[axis] -= cell->lengths[axis] * std::round(d[axis] / cell->lengths[axis]);
            }
        }
        const double r = std::hypot(d[0], d[1], d[2]);
        if (!std::isfinite(r)) {
            throw std::overflow_error("Pair displacement overflowed");
        }
        if (r == 0.0) {
            throw std::domain_error("Active LJ pair has coincident positions");
        }
        if (model.cutoff_mode != CutoffMode::none && r >= model.cutoff) {
            continue;
        }
        auto radial = radial_lj(r, model);
        radial.energy -= at_cutoff.energy;
        if (model.cutoff_mode == CutoffMode::force_shift) {
            radial.energy -= (r - model.cutoff) * at_cutoff.derivative;
            radial.derivative -= at_cutoff.derivative;
        }
        result.energy += pair.scale * radial.energy;
        Vec3 g{};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            g[axis] = pair.scale * radial.derivative * (d[axis] / r);
            result.forces[pair.source][axis] += g[axis];
            result.forces[pair.target][axis] -= g[axis];
        }
        for (std::size_t a = 0; a < 3; ++a) {
            for (std::size_t b = 0; b < 3; ++b) {
                result.virial[3 * a + b] -= d[a] * g[b];
            }
        }
    }
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!finite(result.energy) ||
        !std::all_of(result.virial.begin(), result.virial.end(), finite)) {
        throw std::overflow_error("LJ energy or virial accumulation overflowed");
    }
    for (const auto& force : result.forces) {
        if (!std::all_of(force.begin(), force.end(), finite)) {
            throw std::overflow_error("LJ force accumulation overflowed");
        }
    }
    return result;
}

}  // namespace gmd_next::reference
