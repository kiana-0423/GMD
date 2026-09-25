#include "gmd_next/model/validate.hpp"

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

namespace gmd_next::model {
namespace {

using core::ErrorCode;
using core::ValidationReport;

std::string with_index(std::string_view field, std::size_t index) {
    return std::string(field) + "[" + std::to_string(index) + "]";
}

bool check_finite(ValidationReport& report, double value, std::string field) {
    if (core::all_finite(value)) return true;
    report.add(ErrorCode::non_finite_value, std::move(field), "value must be finite");
    return false;
}

bool check_finite(ValidationReport& report, const core::Vec3& value, std::string field) {
    if (core::all_finite(value)) return true;
    report.add(ErrorCode::non_finite_value, std::move(field), "every component must be finite");
    return false;
}

void check_array_size(ValidationReport& report, std::size_t actual, std::size_t expected,
                      std::string field) {
    if (actual == expected) return;
    report.add(ErrorCode::inconsistent_size, std::move(field),
               "holds " + std::to_string(actual) + " entries for " + std::to_string(expected) +
                   " atoms");
}

void validate_cell_geometry(ValidationReport& report, const CellSpec& cell) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (!check_finite(report, cell.lengths[axis], with_index("cell.lengths", axis))) continue;
        if (cell.lengths[axis] <= 0.0) {
            report.add(ErrorCode::value_out_of_range, with_index("cell.lengths", axis),
                       "edge length must be positive");
        }
    }
}

}  // namespace

std::string_view name(ValidationScope scope) {
    switch (scope) {
    case ValidationScope::reference_static: return "reference_static";
    case ValidationScope::production_r1: return "production_r1";
    }
    return "unknown";
}

core::ValidationReport validate_lj_model(const LjModel& model,
                                         const std::optional<CellSpec>& cell,
                                         ValidationScope scope) {
    ValidationReport report;
    const bool epsilon_finite = check_finite(report, model.epsilon, "model.epsilon");
    const bool sigma_finite = check_finite(report, model.sigma, "model.sigma");
    const bool cutoff_finite = check_finite(report, model.cutoff, "model.cutoff");
    const bool skin_finite = check_finite(report, model.skin, "model.skin");

    if (epsilon_finite && model.epsilon <= 0.0) {
        report.add(ErrorCode::value_out_of_range, "model.epsilon", "epsilon must be positive");
    }
    if (sigma_finite && model.sigma <= 0.0) {
        report.add(ErrorCode::value_out_of_range, "model.sigma", "sigma must be positive");
    }

    if (cutoff_finite) {
        if (model.cutoff_mode == CutoffMode::none) {
            // An uncut model has no radius to consume, so a stray cutoff value
            // is a configuration error rather than an ignored field.
            if (model.cutoff != 0.0) {
                report.add(ErrorCode::malformed_input, "model.cutoff",
                           "cutoff must be exactly zero when cutoff_mode is none");
            }
        } else if (model.cutoff <= 0.0) {
            report.add(ErrorCode::value_out_of_range, "model.cutoff",
                       std::string(name(model.cutoff_mode)) + " requires a positive cutoff");
        }
    }

    if (skin_finite && model.skin < 0.0) {
        report.add(ErrorCode::value_out_of_range, "model.skin", "skin must not be negative");
    }

    if (scope == ValidationScope::reference_static) {
        // The static oracle evaluates an explicit half list and builds nothing,
        // so it cannot honour a list margin.
        if (skin_finite && model.skin != 0.0) {
            report.add(ErrorCode::unsupported_configuration, "model.skin",
                       "static reference evaluation builds no neighbour list; skin must be zero");
        }
    } else {
        if (model.cutoff_mode == CutoffMode::none) {
            report.add(ErrorCode::unsupported_configuration, "model.cutoff_mode",
                       "the production scope requires a finite cutoff");
        }
        if (skin_finite && model.skin <= 0.0) {
            report.add(ErrorCode::value_out_of_range, "model.skin",
                       "the production scope requires skin > 0");
        }
    }

    if (!cell) {
        if (scope == ValidationScope::production_r1) {
            report.add(ErrorCode::unsupported_configuration, "cell",
                       "the production scope requires a periodic orthorhombic cell");
        }
        return report;
    }

    validate_cell_geometry(report, *cell);
    const auto shortest = shortest_periodic_length(*cell);

    if (scope == ValidationScope::reference_static) {
        // Minimum image on one nearest copy: every periodic axis needs a finite
        // cutoff below half its edge.
        if (shortest && core::all_finite(*shortest) && *shortest > 0.0) {
            if (model.cutoff_mode == CutoffMode::none) {
                report.add(ErrorCode::unsupported_configuration, "model.cutoff_mode",
                           "a periodic reference cell requires a finite cutoff");
            } else if (cutoff_finite && model.cutoff >= 0.5 * *shortest) {
                report.add(ErrorCode::value_out_of_range, "model.cutoff",
                           "cutoff must stay below half the shortest periodic edge");
            }
        }
        return report;
    }

    if (!is_fully_periodic(*cell)) {
        report.add(ErrorCode::unsupported_configuration, "cell.periodic",
                   "the production scope requires all three axes to be periodic");
    }
    if (shortest && core::all_finite(*shortest) && *shortest > 0.0 && cutoff_finite &&
        skin_finite) {
        const double radius = list_radius(model);
        if (radius >= 0.5 * *shortest) {
            report.add(ErrorCode::value_out_of_range, "model.skin",
                       "cutoff + skin = " + std::to_string(radius) +
                           " must stay below half the shortest periodic edge " +
                           std::to_string(0.5 * *shortest));
        }
    }
    return report;
}

core::ValidationReport validate_host_state(const HostState& state, ValidationScope scope) {
    ValidationReport report;
    const std::size_t count = state.atom_count();
    const bool production = scope == ValidationScope::production_r1;

    if (production && count == 0) {
        report.add(ErrorCode::unsupported_configuration, "state.ids",
                   "a production run needs at least one atom");
    }
    if (!core::fits_local_index_space(count)) {
        report.add(ErrorCode::capacity_overflow, "state.ids",
                   std::to_string(count) + " atoms exceed the 32-bit device index space");
    }

    check_array_size(report, state.positions.size(), count, "state.positions");
    check_array_size(report, state.masses.size(), count, "state.masses");
    if (production || !state.velocities.empty()) {
        check_array_size(report, state.velocities.size(), count, "state.velocities");
    }
    if (!state.types.empty()) {
        check_array_size(report, state.types.size(), count, "state.types");
    }

    for (std::size_t atom = 0; atom < count; ++atom) {
        if (!state.ids[atom].is_valid()) {
            report.add(ErrorCode::malformed_input, with_index("state.ids", atom),
                       "atom id must be nonnegative");
        }
    }
    // Duplicate stable ids would make output order and restart mapping ambiguous.
    std::vector<core::AtomId> sorted(state.ids);
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t slot = 1; slot < sorted.size(); ++slot) {
        if (sorted[slot] == sorted[slot - 1]) {
            report.add(ErrorCode::duplicate_identity, "state.ids",
                       "atom id " + std::to_string(sorted[slot].value()) + " appears more than once");
            break;
        }
    }

    for (std::size_t atom = 0; atom < state.positions.size(); ++atom) {
        check_finite(report, state.positions[atom], with_index("state.positions", atom));
    }
    for (std::size_t atom = 0; atom < state.velocities.size(); ++atom) {
        check_finite(report, state.velocities[atom], with_index("state.velocities", atom));
    }
    for (std::size_t atom = 0; atom < state.masses.size(); ++atom) {
        if (!check_finite(report, state.masses[atom], with_index("state.masses", atom))) continue;
        if (state.masses[atom] <= 0.0) {
            report.add(ErrorCode::value_out_of_range, with_index("state.masses", atom),
                       "mass must be strictly positive");
        }
    }
    for (std::size_t atom = 0; atom < state.types.size(); ++atom) {
        if (state.types[atom] != kSingleAtomType) {
            report.add(ErrorCode::unsupported_configuration, with_index("state.types", atom),
                       "only the single Lennard-Jones type " + std::to_string(kSingleAtomType) +
                           " is supported");
            break;
        }
    }

    if (state.step < 0) {
        report.add(ErrorCode::value_out_of_range, "state.step",
                   "step index must be nonnegative, received " + std::to_string(state.step));
    }
    if (check_finite(report, state.time, "state.time") && state.time < 0.0) {
        report.add(ErrorCode::value_out_of_range, "state.time", "time must be nonnegative");
    }
    return report;
}

core::ValidationReport validate_run_spec(const RunSpec& spec) {
    ValidationReport report;
    if (spec.units != core::UnitSystem::metal) {
        report.add(ErrorCode::unsupported_configuration, "spec.units",
                   std::string(name(spec.units)) +
                       " units are not implemented; the production scope uses metal units");
    }
    if (spec.ensemble != Ensemble::nve) {
        report.add(ErrorCode::unsupported_configuration, "spec.ensemble",
                   std::string(name(spec.ensemble)) +
                       " is not implemented; the production scope integrates nve");
    }
    if (spec.timestep_unit == core::TimeUnit::femtoseconds &&
        spec.units != core::UnitSystem::metal) {
        report.add(ErrorCode::unit_mismatch, "spec.timestep_unit",
                   std::string(name(spec.units)) + " units have no scale to femtoseconds");
    }
    if (check_finite(report, spec.timestep, "spec.timestep") && spec.timestep <= 0.0) {
        report.add(ErrorCode::value_out_of_range, "spec.timestep",
                   "the time step must be strictly positive");
    }
    if (spec.steps < 0) {
        report.add(ErrorCode::value_out_of_range, "spec.steps",
                   "step count must be nonnegative, received " + std::to_string(spec.steps));
    }
    report.merge(validate_lj_model(spec.model, spec.cell, ValidationScope::production_r1));
    return report;
}

core::ValidationReport validate_run(const RunSpec& spec, const HostState& state) {
    ValidationReport report = validate_run_spec(spec);
    report.merge(validate_host_state(state, ValidationScope::production_r1));
    if (state.time_unit == core::TimeUnit::femtoseconds && spec.units != core::UnitSystem::metal) {
        report.add(ErrorCode::unit_mismatch, "state.time_unit",
                   std::string(name(spec.units)) + " units have no scale to femtoseconds");
    }
    if (spec.remove_initial_center_of_mass_motion && state.atom_count() < 2) {
        report.add(ErrorCode::unsupported_configuration, "spec.remove_initial_center_of_mass_motion",
                   "removing net momentum leaves no degrees of freedom below two atoms");
    }
    return report;
}

NormalizedState::NormalizedState(const HostState& state, core::UnitSystem units,
                                 bool remove_center_of_mass_motion) {
    validate_host_state(state, ValidationScope::production_r1).require_ok();
    if (remove_center_of_mass_motion && state.atom_count() < 2) {
        throw core::ContractError({ErrorCode::unsupported_configuration,
                                   "remove_center_of_mass_motion",
                                   "removing net momentum leaves no degrees of freedom below "
                                   "two atoms"});
    }

    const std::size_t count = state.atom_count();
    input_order_.resize(count);
    std::iota(input_order_.begin(), input_order_.end(), std::size_t{0});
    // Ascending stable id is the one order that survives device permutation,
    // migration and restart, so output rows are mapped from it.
    std::sort(input_order_.begin(), input_order_.end(),
              [&state](std::size_t lhs, std::size_t rhs) {
                  return state.ids[lhs] < state.ids[rhs];
              });

    ids_.reserve(count);
    positions_.reserve(count);
    velocities_.reserve(count);
    masses_.reserve(count);
    for (const std::size_t source : input_order_) {
        ids_.push_back(state.ids[source]);
        positions_.push_back(state.positions[source]);
        velocities_.push_back(state.velocities.empty() ? core::Vec3{0.0, 0.0, 0.0}
                                                       : state.velocities[source]);
        masses_.push_back(state.masses[source]);
    }

    total_mass_ = std::accumulate(masses_.begin(), masses_.end(), 0.0);
    degrees_of_freedom_ = 3 * static_cast<std::int64_t>(count);
    if (remove_center_of_mass_motion) {
        core::Vec3 momentum{0.0, 0.0, 0.0};
        for (std::size_t atom = 0; atom < count; ++atom) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                momentum[axis] += masses_[atom] * velocities_[atom][axis];
            }
        }
        for (std::size_t axis = 0; axis < 3; ++axis) {
            removed_velocity_[axis] = momentum[axis] / total_mass_;
        }
        for (auto& velocity : velocities_) {
            for (std::size_t axis = 0; axis < 3; ++axis) velocity[axis] -= removed_velocity_[axis];
        }
        // Three constraints removed with the net momentum.
        degrees_of_freedom_ -= 3;
        center_of_mass_removed_ = true;
    }
    if (!core::all_finite(total_mass_) || !core::all_finite(removed_velocity_)) {
        throw core::ContractError({ErrorCode::non_finite_value, "state.masses",
                                   "mass or momentum accumulation overflowed"});
    }

    origin_ = core::StepStamp{state.step,
                              core::to_internal_time(state.time, state.time_unit, units)};
    core::validate_step_stamp(origin_, "state").require_ok();
}

namespace {

// Validates before the normalized member is constructed, so a rejected run
// never reaches a partially built ValidatedRun.
const HostState& checked(const RunSpec& spec, const HostState& state) {
    validate_run(spec, state).require_ok();
    return state;
}

}  // namespace

ValidatedRun::ValidatedRun(const RunSpec& spec, const HostState& state)
    : spec_(spec),
      state_(checked(spec, state), spec.units, spec.remove_initial_center_of_mass_motion),
      timestep_internal_(core::to_internal_time(spec.timestep, spec.timestep_unit, spec.units)),
      initial_versions_(core::StateVersions{}.stamp()) {}

core::ObservableRequest ValidatedRun::force_request() const {
    // NVE needs only the force. Energy and virial are computed because output
    // asked for them, which is what keeps them out of a no-output run.
    return core::ObservableRequest::for_dynamics_step(
        spec_.output_observables, core::SamplingStage::after_force_evaluation);
}

}  // namespace gmd_next::model
