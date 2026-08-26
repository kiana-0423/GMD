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

// Diagnostics produced while normalising a constraint list.
//
// A repeated atom pair falls into exactly one of three categories, decided by
// how far apart the two target distances are:
//
//   |d1 - d2| == 0                    exact duplicate
//   0 < |d1 - d2| <= tolerance        tolerance-equivalent duplicate
//   |d1 - d2| >  tolerance            conflict (rejected outright)
//
// The middle category is a judgement call, so it is counted separately and each
// occurrence is recorded in `discarded_targets` with both distances and the
// tolerance that was applied. The first value encountered wins.
struct ConstraintNormalizationDiagnostics {
    std::size_t exact_duplicates = 0;
    std::size_t tolerance_equivalent_duplicates = 0;

    // One entry per tolerance-equivalent duplicate whose target differed from
    // the value that was kept. Human-readable; the application decides how to
    // surface these (gmd prints them once at start-up).
    std::vector<std::string> discarded_targets;

    bool empty() const noexcept {
        return exact_duplicates == 0 && tolerance_equivalent_duplicates == 0;
    }
};

// Holonomic bond-length constraints solved with SHAKE/RATTLE.
//
// Constraint list normalisation (performed by the constructor):
//   - atom pairs are stored as (min, max), so (i, j) and (j, i) are the same
//     constraint and are normalised identically regardless of the orientation
//     they were supplied in;
//   - repeated pairs are collapsed to one entry, keeping the FIRST target
//     distance encountered, which makes the result deterministic and
//     independent of pair orientation;
//   - a repeat whose target differs by more than the constraint tolerance is
//     rejected as a conflict, since no geometry satisfies both;
//   - negative atom tags, self-constraints (i == j) and non-positive or
//     non-finite target distances are rejected.
//
// A repeat whose target differs by more than zero but at most the tolerance is
// a "tolerance-equivalent duplicate": SHAKE converges to within `tolerance`, so
// the two targets are not distinguishable by the solver. It is accepted, the
// first value is kept, and it is reported through
// normalization_diagnostics().
//
// INDEPENDENCE IS ASSUMED, NOT VERIFIED. After normalisation the list contains
// distinct constraints, but distinct is not the same as independent: a closed
// topology such as an all-pairs cage over five or more atoms contains more
// distance constraints than the rigid body has removable degrees of freedom.
// Deciding that in general means computing the rank of the 3N x M constraint
// Jacobian, which is configuration dependent (the rank can drop at particular
// geometries) and is not attempted here. A redundant set is accepted and every
// distinct constraint is counted, so degrees_of_freedom() over-subtracts and
// the reported temperature comes out high. Configure independent constraints.
class ConstraintSolver {
public:
    ConstraintSolver() = default;
    ConstraintSolver(std::vector<BondConstraint> constraints,
                     ConstraintSettings settings);

    bool enabled() const noexcept { return !constraints_.empty(); }
    const ConstraintSettings& settings() const noexcept { return settings_; }

    // The normalised list: distinct pairs, each stored as (min, max).
    const std::vector<BondConstraint>& constraints() const noexcept { return constraints_; }

    // Number of distinct active constraints, i.e. the number of degrees of
    // freedom removed *assuming the configured constraints are independent*.
    std::size_t active_constraint_count() const noexcept { return constraints_.size(); }

    // Total repeats dropped, of either kind.
    std::size_t dropped_duplicate_count() const noexcept {
        return diagnostics_.exact_duplicates +
               diagnostics_.tolerance_equivalent_duplicates;
    }

    // Per-category counts and the human-readable record of any non-identical
    // target that was discarded. ConstraintSolver deliberately does not print
    // anything itself: it is a low-level library class and the project has no
    // logging abstraction, so diagnostics are exposed to the caller instead.
    // `gmd` reports them on stdout once, right after the constraint summary.
    const ConstraintNormalizationDiagnostics& normalization_diagnostics() const noexcept {
        return diagnostics_;
    }

    ConstraintProjectionStats apply_shake(System& system) const;
    ConstraintProjectionStats apply_rattle(System& system) const;

private:
    std::vector<BondConstraint> constraints_;
    ConstraintSettings settings_;
    ConstraintNormalizationDiagnostics diagnostics_;
};

std::vector<BondConstraint> constraints_from_bond_types(
    const Topology& topology,
    const std::vector<int>& constrained_bond_types,
    const std::vector<double>& bond_type_distances);

}  // namespace gmd
