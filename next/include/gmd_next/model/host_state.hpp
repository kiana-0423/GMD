#pragma once

#include "gmd_next/core/identity.hpp"
#include "gmd_next/core/numeric.hpp"
#include "gmd_next/core/units.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace gmd_next::model {

// The single atom type the first production scope supports.
inline constexpr std::int32_t kSingleAtomType = 0;

// Input state owned by the host. It is the authority until upload; after that
// the device state advances and this snapshot is not updated per step.
struct HostState {
    std::vector<core::AtomId> ids;        // stable, nonnegative, no duplicates
    std::vector<core::Vec3> positions;    // Angstrom, not required to be wrapped
    std::vector<core::Vec3> velocities;   // Angstrom per internal time unit;
                                          // may be empty only for static use
    std::vector<double> masses;           // amu, strictly positive
    std::vector<std::int32_t> types;      // empty, or one entry per atom, all equal
    std::int64_t step = 0;                // step index this state belongs to
    double time = 0.0;                    // expressed in time_unit
    core::TimeUnit time_unit = core::TimeUnit::femtoseconds;

    std::size_t atom_count() const { return ids.size(); }
};

// A validated state in a fixed, reproducible order: ascending stable id. The
// order is what output rows and device slots are mapped from, so it may not
// depend on how the input file happened to list the atoms.
class NormalizedState {
public:
    // Validates in the production scope and normalizes. Throws ContractError
    // when the state, the unit system or the requested normalization is invalid.
    NormalizedState(const HostState& state, core::UnitSystem units,
                    bool remove_center_of_mass_motion);

    std::size_t atom_count() const { return ids_.size(); }
    std::span<const core::AtomId> ids() const { return ids_; }
    std::span<const core::Vec3> positions() const { return positions_; }
    std::span<const core::Vec3> velocities() const { return velocities_; }
    std::span<const double> masses() const { return masses_; }
    std::int32_t type() const { return kSingleAtomType; }
    // Normalized slot -> index in the input arrays, so a caller can report a
    // rejection against the row the user actually wrote.
    std::span<const std::size_t> input_order() const { return input_order_; }

    core::StepStamp origin() const { return origin_; }
    double total_mass() const { return total_mass_; }
    std::int64_t degrees_of_freedom() const { return degrees_of_freedom_; }
    bool center_of_mass_motion_removed() const { return center_of_mass_removed_; }
    // The velocity that was subtracted; zero when nothing was removed.
    core::Vec3 removed_center_of_mass_velocity() const { return removed_velocity_; }

private:
    std::vector<core::AtomId> ids_;
    std::vector<core::Vec3> positions_;
    std::vector<core::Vec3> velocities_;
    std::vector<double> masses_;
    std::vector<std::size_t> input_order_;
    core::StepStamp origin_;
    double total_mass_ = 0.0;
    std::int64_t degrees_of_freedom_ = 0;
    bool center_of_mass_removed_ = false;
    core::Vec3 removed_velocity_{0.0, 0.0, 0.0};
};

}  // namespace gmd_next::model
