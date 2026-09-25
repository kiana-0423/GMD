#pragma once

#include "gmd_next/core/observables.hpp"
#include "gmd_next/core/status.hpp"
#include "gmd_next/core/version.hpp"
#include "gmd_next/model/cell.hpp"
#include "gmd_next/model/host_state.hpp"
#include "gmd_next/model/lj_model.hpp"
#include "gmd_next/model/run_spec.hpp"

#include <cstdint>
#include <optional>
#include <string_view>

namespace gmd_next::model {

// Two scopes, deliberately kept apart.
//
// reference_static  what the host CPU oracle in gmd_next::reference evaluates:
//                   one explicit half list, no neighbour build, no time
//                   advance. Uncut and partially periodic models stay legal
//                   here, and a nonzero skin is rejected because no list is
//                   built that could consume it.
// production_r1     the first deliverable: single GPU, FP64, metal units,
//                   single-type Lennard-Jones with a finite cutoff, fully
//                   periodic fixed orthorhombic cell, skin > 0 and
//                   rc + skin < min(L)/2, velocity Verlet NVE.
//
// Widening the production scope is a milestone, not a validation relaxation.
enum class ValidationScope : std::uint8_t { reference_static, production_r1 };

std::string_view name(ValidationScope scope);

// Cell is optional in the reference scope only; production requires one.
core::ValidationReport validate_lj_model(const LjModel& model,
                                         const std::optional<CellSpec>& cell,
                                         ValidationScope scope);
core::ValidationReport validate_host_state(const HostState& state, ValidationScope scope);
// Model, cell, units, ensemble, time step and step count of a production run.
core::ValidationReport validate_run_spec(const RunSpec& spec);
// Everything above plus the checks that need both sides, reported together.
core::ValidationReport validate_run(const RunSpec& spec, const HostState& state);

// A frozen, normalized run. Constructing one is the only way to reach the
// normalized state, so nothing downstream can consume unvalidated input.
class ValidatedRun {
public:
    // Throws ContractError carrying the first diagnostic of validate_run.
    ValidatedRun(const RunSpec& spec, const HostState& state);

    const RunSpec& spec() const { return spec_; }
    const NormalizedState& state() const { return state_; }
    // Time step in the integrator's internal unit, converted exactly once.
    double timestep_internal() const { return timestep_internal_; }
    std::int64_t degrees_of_freedom() const { return state_.degrees_of_freedom(); }
    // Coordinates, cell, model and permutation all start at their first valid
    // version; the runtime that owns the state advances them from here.
    core::VersionStamp initial_versions() const { return initial_versions_; }

    // What one force evaluation must produce: the force the integrator needs,
    // plus whichever of energy and virial were requested for output.
    core::ObservableRequest force_request() const;

private:
    RunSpec spec_;
    NormalizedState state_;
    double timestep_internal_ = 0.0;
    core::VersionStamp initial_versions_;
};

}  // namespace gmd_next::model
