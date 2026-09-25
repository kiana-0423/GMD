#include "gmd_next/core/observables.hpp"
#include "gmd_next/core/units.hpp"
#include "gmd_next/core/version.hpp"
#include "gmd_next/model/validate.hpp"
#include "gmd_next/reference.hpp"

#include "contract_support.hpp"
#include "test_support.hpp"

#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

namespace {
using namespace gmd_next;
using core::ErrorCode;
using model::ValidationScope;

constexpr double kNan = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// Argon-like parameters in metal units, inside the production scope:
// rc + skin = 9.5 A < min(L)/2 = 15 A.
model::RunSpec production_spec() {
    model::RunSpec spec;
    spec.timestep = 1.0;
    spec.steps = 10;
    spec.cell = model::CellSpec{{30.0, 32.0, 34.0}};
    spec.model = model::LjModel{0.0103, 3.4, model::CutoffMode::potential_shift, 8.5, 1.0};
    return spec;
}

// Deliberately listed out of id order, with a nonzero net momentum.
model::HostState production_state() {
    model::HostState state;
    state.ids = {core::AtomId{7}, core::AtomId{2}, core::AtomId{5}};
    state.positions = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    state.velocities = {{0.01, 0.0, -0.02}, {0.03, 0.01, 0.0}, {-0.01, 0.02, 0.01}};
    state.masses = {39.948, 20.0, 12.5};
    return state;
}

void production_run_is_normalized() {
    const auto spec = production_spec();
    const auto state = production_state();
    accepts(model::validate_run(spec, state));
    const model::ValidatedRun run(spec, state);

    // Ascending stable id, with the input rows carried along consistently.
    const auto& normalized = run.state();
    expect(normalized.ids()[0] == core::AtomId{2} && normalized.ids()[1] == core::AtomId{5} &&
               normalized.ids()[2] == core::AtomId{7},
           "Normalized state must be ordered by stable id");
    for (std::size_t slot = 0; slot < normalized.atom_count(); ++slot) {
        const auto source = normalized.input_order()[slot];
        expect(normalized.ids()[slot] == state.ids[source], "Input order lost an id");
        expect(normalized.positions()[slot] == state.positions[source], "Positions were not permuted");
        expect(normalized.velocities()[slot] == state.velocities[source], "Velocities were not permuted");
        expect(normalized.masses()[slot] == state.masses[source], "Masses were not permuted");
    }

    // 1 fs in A*sqrt(amu/eV), derived from SI values rather than the header.
    const double femtoseconds_per_internal = 1.0e5 * std::sqrt(1.66053906892e-27 / 1.602176634e-19);
    near(run.timestep_internal(), 1.0 / femtoseconds_per_internal, 0.0, 1e-12);
    expect(run.degrees_of_freedom() == 9, "Without momentum removal DOF is 3N");
    expect(!normalized.center_of_mass_motion_removed(), "Nothing was asked to be removed");
    expect(run.initial_versions().is_complete() &&
               run.initial_versions().matches(core::StateVersions{}.stamp()),
           "A validated run must start from the first valid version of every domain");

    // No output: the evaluation produces only what NVE needs.
    const auto quiet = run.force_request();
    expect(quiet.required() == core::ObservableSet{core::Observable::forces},
           "A no-output run must not compute energy or virial");
    core::ObservableRecord record(quiet, run.initial_versions(), normalized.origin());
    record.store_forces();
    accepts(record.check_complete());
    expect(record.potential_energy().availability() == core::Availability::unavailable,
           "An unrequested energy must be unavailable, not zero");

    auto logging_spec = spec;
    logging_spec.output_observables = {core::Observable::potential_energy, core::Observable::virial};
    const auto logging = model::ValidatedRun(logging_spec, state).force_request();
    expect(logging.required().size() == 3, "Output must add energy and virial to the force");
    expect(logging.is_output_only(core::Observable::potential_energy) &&
               logging.is_output_only(core::Observable::virial),
           "Energy and virial are output-only under NVE");

    // A restart state converts its time at the boundary, exactly once.
    auto restart = state;
    restart.step = 40;
    restart.time = 20.0;
    const model::ValidatedRun resumed(spec, restart);
    expect(resumed.state().origin().step == 40, "Origin lost its step");
    near(resumed.state().origin().time_internal, 20.0 / femtoseconds_per_internal, 0.0, 1e-12);
}

void center_of_mass_and_dof() {
    auto spec = production_spec();
    spec.remove_initial_center_of_mass_motion = true;
    const auto state = production_state();
    const model::ValidatedRun run(spec, state);

    core::Vec3 momentum{};
    double mass = 0.0;
    for (std::size_t atom = 0; atom < state.atom_count(); ++atom) {
        mass += state.masses[atom];
        for (std::size_t axis = 0; axis < 3; ++axis) {
            momentum[axis] += state.masses[atom] * state.velocities[atom][axis];
        }
    }
    const auto& normalized = run.state();
    expect(normalized.center_of_mass_motion_removed(), "Removal was requested");
    expect(run.degrees_of_freedom() == 6, "Removing net momentum removes three DOF");
    near(normalized.total_mass(), mass);
    for (std::size_t axis = 0; axis < 3; ++axis) {
        near(normalized.removed_center_of_mass_velocity()[axis], momentum[axis] / mass);
        double residual = 0.0;
        for (std::size_t slot = 0; slot < normalized.atom_count(); ++slot) {
            residual += normalized.masses()[slot] * normalized.velocities()[slot][axis];
        }
        near(residual, 0.0, 1e-14, 0.0);
    }

    auto single = state;
    single.ids.resize(1);
    single.positions.resize(1);
    single.velocities.resize(1);
    single.masses.resize(1);
    rejects(model::validate_run(spec, single), ErrorCode::unsupported_configuration,
            "spec.remove_initial_center_of_mass_motion");
    throws_code<core::ContractError>(ErrorCode::unsupported_configuration,
                                     [&] { model::ValidatedRun(spec, single); });
    spec.remove_initial_center_of_mass_motion = false;
    expect(model::ValidatedRun(spec, single).degrees_of_freedom() == 3, "One atom has 3 DOF");
}

void run_spec_rejections() {
    const auto state = production_state();
    const auto check = [&](auto edit, ErrorCode code, const char* field) {
        auto spec = production_spec();
        edit(spec);
        rejects(model::validate_run(spec, state), code, field);
        throws_code<core::ContractError>(code, [&] { model::ValidatedRun(spec, state); });
    };
    using Spec = model::RunSpec;
    check([](Spec& s) { s.timestep = 0.0; }, ErrorCode::value_out_of_range, "spec.timestep");
    check([](Spec& s) { s.timestep = -0.5; }, ErrorCode::value_out_of_range, "spec.timestep");
    check([](Spec& s) { s.timestep = kNan; }, ErrorCode::non_finite_value, "spec.timestep");
    check([](Spec& s) { s.steps = -1; }, ErrorCode::value_out_of_range, "spec.steps");
    check([](Spec& s) { s.ensemble = model::Ensemble::nvt; },
          ErrorCode::unsupported_configuration, "spec.ensemble");
    check([](Spec& s) { s.ensemble = model::Ensemble::npt; },
          ErrorCode::unsupported_configuration, "spec.ensemble");
    check([](Spec& s) { s.units = core::UnitSystem::reduced; },
          ErrorCode::unsupported_configuration, "spec.units");
    check([](Spec& s) { s.model.epsilon = 0.0; }, ErrorCode::value_out_of_range, "model.epsilon");
    check([](Spec& s) { s.model.sigma = -1.0; }, ErrorCode::value_out_of_range, "model.sigma");
    check([](Spec& s) { s.model.epsilon = kInf; }, ErrorCode::non_finite_value, "model.epsilon");
    check([](Spec& s) { s.model.skin = 0.0; }, ErrorCode::value_out_of_range, "model.skin");
    check([](Spec& s) { s.model.skin = -0.1; }, ErrorCode::value_out_of_range, "model.skin");
    check([](Spec& s) { s.model.cutoff = 0.0; }, ErrorCode::value_out_of_range, "model.cutoff");
    check([](Spec& s) {
        s.model.cutoff_mode = model::CutoffMode::none;
        s.model.cutoff = 0.0;
    }, ErrorCode::unsupported_configuration, "model.cutoff_mode");
    check([](Spec& s) { s.cell.periodic[2] = false; },
          ErrorCode::unsupported_configuration, "cell.periodic");
    check([](Spec& s) { s.cell.lengths[1] = 0.0; }, ErrorCode::value_out_of_range,
          "cell.lengths[1]");
    check([](Spec& s) { s.cell.lengths[0] = kNan; }, ErrorCode::non_finite_value,
          "cell.lengths[0]");

    // cutoff + skin against half of the shortest edge (30 A here): the list
    // radius 15 is exactly at the limit and must fail, just below it passes.
    check([](Spec& s) { s.model.skin = 6.5; }, ErrorCode::value_out_of_range, "model.skin");
    check([](Spec& s) { s.model.cutoff = 14.9; }, ErrorCode::value_out_of_range, "model.skin");
    auto limit = production_spec();
    limit.model.skin = 6.5 - 1e-9;
    accepts(model::validate_run(limit, state));
    // A cutoff that fits alone but not with the skin is still rejected.
    limit.model.cutoff = 14.0;
    limit.model.skin = 1.0;
    rejects(model::validate_run(limit, state), ErrorCode::value_out_of_range, "model.skin");

    // Every problem in one input is reported together.
    auto broken = production_spec();
    broken.timestep = 0.0;
    broken.model.skin = 0.0;
    broken.ensemble = model::Ensemble::npt;
    const auto report = model::validate_run(broken, state);
    expect(report.diagnostics().size() == 3, "Expected three diagnostics: " + report.summary());
}

void host_state_rejections() {
    const auto spec = production_spec();
    const auto check = [&](auto edit, ErrorCode code, const char* field) {
        auto state = production_state();
        edit(state);
        rejects(model::validate_run(spec, state), code, field);
        throws_code<core::ContractError>(code, [&] { model::ValidatedRun(spec, state); });
    };
    using State = model::HostState;
    check([](State& s) { s.masses[1] = 0.0; }, ErrorCode::value_out_of_range, "state.masses[1]");
    check([](State& s) { s.masses[0] = -4.0; }, ErrorCode::value_out_of_range, "state.masses[0]");
    check([](State& s) { s.masses[2] = kNan; }, ErrorCode::non_finite_value, "state.masses[2]");
    check([](State& s) { s.positions[1][2] = kNan; }, ErrorCode::non_finite_value,
          "state.positions[1]");
    check([](State& s) { s.velocities[0][0] = kInf; }, ErrorCode::non_finite_value,
          "state.velocities[0]");
    check([](State& s) { s.positions.pop_back(); }, ErrorCode::inconsistent_size, "state.positions");
    check([](State& s) { s.masses.push_back(1.0); }, ErrorCode::inconsistent_size, "state.masses");
    check([](State& s) { s.velocities.clear(); }, ErrorCode::inconsistent_size, "state.velocities");
    check([](State& s) { s.types = {0, 0}; }, ErrorCode::inconsistent_size, "state.types");
    check([](State& s) { s.types = {0, 1, 0}; }, ErrorCode::unsupported_configuration,
          "state.types[1]");
    check([](State& s) { s.ids[2] = core::AtomId{7}; }, ErrorCode::duplicate_identity, "state.ids");
    check([](State& s) { s.ids[0] = core::AtomId{-3}; }, ErrorCode::malformed_input, "state.ids[0]");
    check([](State& s) { s.step = -2; }, ErrorCode::value_out_of_range, "state.step");
    check([](State& s) { s.time = kNan; }, ErrorCode::non_finite_value, "state.time");
    check([](State& s) { s = model::HostState{}; }, ErrorCode::unsupported_configuration,
          "state.ids");

    // Uniform types equal to the supported one are accepted explicitly.
    auto typed = production_state();
    typed.types = {0, 0, 0};
    accepts(model::validate_run(spec, typed));

    // A time in fs cannot be converted in reduced units; the state says so itself.
    auto reduced = production_spec();
    reduced.units = core::UnitSystem::reduced;
    reduced.timestep_unit = core::TimeUnit::internal;
    rejects(model::validate_run(reduced, production_state()), ErrorCode::unit_mismatch,
            "state.time_unit");
}

void reference_scope_is_separate() {
    const model::LjModel uncut{1.0, 1.0};
    const model::LjModel cut{1.0, 1.0, model::CutoffMode::force_shift, 2.5};
    const model::CellSpec cell{{10.0, 10.0, 10.0}};

    // Static cases the oracle already evaluates remain legal in the reference scope.
    accepts(model::validate_lj_model(uncut, std::nullopt, ValidationScope::reference_static));
    accepts(model::validate_lj_model(cut, cell, ValidationScope::reference_static));
    model::CellSpec open_z{{10.0, 10.0, 3.0}, {true, true, false}};
    accepts(model::validate_lj_model(cut, open_z, ValidationScope::reference_static));
    model::CellSpec open{{10.0, 10.0, 10.0}, {false, false, false}};
    accepts(model::validate_lj_model(uncut, open, ValidationScope::reference_static));

    // ...and are still outside the production scope.
    rejects(model::validate_lj_model(uncut, std::nullopt, ValidationScope::production_r1),
            ErrorCode::unsupported_configuration, "model.cutoff_mode");
    rejects(model::validate_lj_model(cut, std::nullopt, ValidationScope::production_r1),
            ErrorCode::unsupported_configuration, "cell");
    rejects(model::validate_lj_model(cut, open_z, ValidationScope::production_r1),
            ErrorCode::unsupported_configuration, "cell.periodic");

    // Reference-scope rules match what the oracle enforces.
    rejects(model::validate_lj_model(uncut, cell, ValidationScope::reference_static),
            ErrorCode::unsupported_configuration, "model.cutoff_mode");
    const model::LjModel too_long{1.0, 1.0, model::CutoffMode::potential_shift, 5.0};
    rejects(model::validate_lj_model(too_long, cell, ValidationScope::reference_static),
            ErrorCode::value_out_of_range, "model.cutoff");
    rejects(model::validate_lj_model({1.0, 1.0, model::CutoffMode::none, 2.5}, std::nullopt,
                                     ValidationScope::reference_static),
            ErrorCode::malformed_input, "model.cutoff");
    // A skin is list configuration; the static path builds no list and says so.
    rejects(model::validate_lj_model({1.0, 1.0, model::CutoffMode::potential_shift, 2.5, 0.3}, cell,
                                     ValidationScope::reference_static),
            ErrorCode::unsupported_configuration, "model.skin");

    // A static state needs neither velocities nor atoms.
    model::HostState positions_only;
    positions_only.ids = {core::AtomId{0}, core::AtomId{1}};
    positions_only.positions = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    positions_only.masses = {1.0, 1.0};
    accepts(model::validate_host_state(positions_only, ValidationScope::reference_static));
    accepts(model::validate_host_state(model::HostState{}, ValidationScope::reference_static));
    rejects(model::validate_host_state(positions_only, ValidationScope::production_r1),
            ErrorCode::inconsistent_size, "state.velocities");
    // Duplicates and bad masses are rejected in both scopes.
    positions_only.ids[1] = core::AtomId{0};
    rejects(model::validate_host_state(positions_only, ValidationScope::reference_static),
            ErrorCode::duplicate_identity, "state.ids");
}

// The contract is mapped onto the oracle by hand here, so neither library
// depends on the other and a mapping error shows up as a wrong number.
reference::LennardJones to_reference(const model::LjModel& lj) {
    const auto mode = lj.cutoff_mode == model::CutoffMode::potential_shift
                          ? reference::CutoffMode::potential_shift
                          : lj.cutoff_mode == model::CutoffMode::force_shift
                                ? reference::CutoffMode::force_shift
                                : reference::CutoffMode::none;
    return reference::LennardJones{lj.epsilon, lj.sigma, mode, lj.cutoff};
}

void oracle_agreement() {
    auto spec = production_spec();
    spec.model = model::LjModel{1.0, 1.0, model::CutoffMode::potential_shift, 2.5, 0.5};
    model::HostState state;
    state.ids = {core::AtomId{1}, core::AtomId{0}};
    state.positions = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    state.velocities = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    state.masses = {1.0, 1.0};
    const model::ValidatedRun run(spec, state);
    const auto positions = run.state().positions();
    expect(positions[0] == core::Vec3{1.0, 0.0, 0.0}, "Id 0 must occupy the first slot");

    const std::array<reference::Pair, 1> pairs{{{0, 1}}};
    const reference::OrthorhombicCell cell{spec.cell.lengths, spec.cell.periodic};
    const std::vector<reference::Vec3> x(positions.begin(), positions.end());

    // Hand-derived at r = 1, rc = 2.5, epsilon = sigma = 1:
    //   2.5^-6 = 64/15625 = 0.004096,  2.5^-12 = 1.6777216e-5
    //   phi(2.5)  = 4 (2.5^-12 - 2.5^-6)            = -0.016316891136
    //   phi'(1)   = 24 (1 - 2)                      = -24
    //   phi'(2.5) = (24/2.5)(2.5^-6 - 2 * 2.5^-12)  =  0.0389994774528
    // Slot 0 sits at +x of slot 1 and is pushed further along +x.
    const auto shifted = reference::evaluate_lj(x, pairs, to_reference(spec.model), cell);
    near(shifted.energy, 0.016316891136);
    near(shifted.forces[0][0], 24.0);
    near(shifted.forces[1][0], -24.0);
    near(shifted.virial[0], 24.0);

    spec.model.cutoff_mode = model::CutoffMode::force_shift;
    const auto forced = reference::evaluate_lj(x, pairs, to_reference(spec.model), cell);
    near(forced.forces[0][0], 24.0389994774528);
    // u = phi(1) - phi(rc) - (1 - rc) phi'(rc)
    near(forced.energy, 0.016316891136 + 1.5 * 0.0389994774528);
}

}  // namespace

int main() {
    try {
        production_run_is_normalized();
        center_of_mass_and_dof();
        run_spec_rejections();
        host_state_rejections();
        reference_scope_is_separate();
        oracle_agreement();
        std::cout << "Run/state validation, normalization, scopes and oracle mapping passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
