#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "gmd/system/topology.hpp"

namespace gmd {

class System;

struct ConstraintSettings {
    double tolerance = 1.0e-6;
    int max_iterations = 100;
    bool enable_rattle = true;
};

struct ConstraintProjectionStats {
    bool enabled = false;
    bool converged = true;
    int iterations = 0;
    double max_error = 0.0;
    std::string stage;
};

class ConstraintSolver {
public:
    ConstraintSolver() = default;
    ConstraintSolver(std::vector<BondConstraint> constraints,
                     ConstraintSettings settings);

    bool enabled() const noexcept { return !constraints_.empty(); }
    const ConstraintSettings& settings() const noexcept { return settings_; }
    const std::vector<BondConstraint>& constraints() const noexcept { return constraints_; }

    ConstraintProjectionStats apply_shake(System& system) const;
    ConstraintProjectionStats apply_rattle(System& system) const;

private:
    std::vector<BondConstraint> constraints_;
    ConstraintSettings settings_;
};

std::vector<BondConstraint> constraints_from_bond_types(
    const Topology& topology,
    const std::vector<int>& constrained_bond_types,
    const std::vector<double>& bond_type_distances);

}  // namespace gmd
