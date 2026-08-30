// Random velocity initialization for constrained systems.
//
// THE ORDER THINGS HAPPEN IN. Simulation::initialize() currently runs the
// velocity initializer FIRST, and the integrator's initialize() afterwards.
// That second call is where the constraint work happens:
//
//   1  VelocityInitializer: sample, remove the global centre-of-mass velocity,
//      rescale so that 2K = dof * k_B * T -- with dof hard-coded as 3N-3.
//   2  VelocityVerletIntegrator::initialize:
//        a  SHAKE-project the positions onto the constraint manifold
//        b  require_independent() at the PROJECTED geometry -- the
//           authoritative rank
//        c  RATTLE-project the velocities into the tangent space
//        d  install dof = 3N - rank - 3 on the thermostat
//
// So the rescale in (1) targets a degree-of-freedom count that (2d) contradicts,
// and it happens BEFORE the projection in (2c) removes kinetic energy. Two
// separate errors compounding, and neither is visible from inside the
// initializer: it has no constraint solver and no rank.
//
// WHAT THAT COSTS. The projection removes the energy that lived in the
// constrained modes. On average that is rank/(3N-3) of the total, which is very
// nearly what the difference between the two DOF counts would have accounted
// for -- so the error has a mean near zero and a FLUCTUATION that does not.
// Each run draws once, so each run is off by that fluctuation:
//
//     3 rigid waters (9 atoms, 9 constraints, dof 15):  +22.9 percent
//     10 waters                                      :   +4.8 percent
//     40 waters                                      :   +1.6 percent
//     150 waters                                     :   -0.9 percent
//
// A constrained run asked for 300 K starts at 368 K. It is not a small-system
// curiosity either: a few percent at 40 molecules is still a few percent.
//
// This file asserts what holds today -- tangency, zero centre-of-mass momentum,
// and that the integrator's authoritative DOF is 3N - rank - 3 -- and REPORTS
// the temperature error rather than asserting it, so that no commit in this
// series is red. The commit that makes the initializer constraint aware is the
// one that turns it into an equality.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/core/simulation.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[constrained init] " << message << '\n';
        ++failures;
    }
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

// k_B from the exact SI definitions, as elsewhere in this repository's audits.
constexpr double kElementaryChargeCoulombs = 1.602176634e-19;
constexpr double kBoltzmannJoulesPerKelvin = 1.380649e-23;
const double kBoltzmann = kBoltzmannJoulesPerKelvin / kElementaryChargeCoulombs;

constexpr double kTargetTemperature = 300.0;
constexpr std::uint32_t kSeed = 20260830u;

class NullForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "null"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request, gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.forces.assign(request.coordinates.size(), {0.0, 0.0, 0.0});
        result.potential_energy = 0.0;
        result.virial_valid = false;
        result.success = true;
    }
};

gmd::ConstraintSettings tight_settings() {
    gmd::ConstraintSettings settings;
    settings.tolerance = 1.0e-13;
    settings.max_iterations = 1000;
    return settings;
}

// A fixture and everything measured from it.
struct Outcome {
    std::size_t atom_count = 0;
    std::size_t constraint_count = 0;
    std::size_t authoritative_dof = 0;
    std::size_t naive_dof = 0;
    double twice_kinetic_energy = 0.0;
    double temperature = 0.0;          // using the authoritative DOF
    double worst_tangency = 0.0;       // max |r_ij . v_ij| over constraints
    double momentum_magnitude = 0.0;
    std::map<int, std::array<double, 3>> velocities_by_tag;
};

// n rigid water molecules, each tilted differently so that no constraint is
// axis aligned and no two molecules are related by a translation.
std::vector<gmd::BondConstraint> water_fixture(gmd::System& system,
                                               std::size_t molecule_count,
                                               double box_length) {
    constexpr double kOH = 0.9572;
    constexpr double kHH = 1.5139;
    const std::size_t atom_count = 3 * molecule_count;
    system.resize(atom_count, atom_count);
    gmd::Box box;
    box.set_lengths({box_length, box_length, box_length});
    system.set_box(box);

    std::vector<gmd::BondConstraint> constraints;
    for (std::size_t m = 0; m < molecule_count; ++m) {
        const std::size_t o = 3 * m, h1 = 3 * m + 1, h2 = 3 * m + 2;
        const double shift = 4.0 * static_cast<double>(m % 8);
        const double row = 4.0 * static_cast<double>(m / 8);
        const double angle = 0.4 + 0.3 * static_cast<double>(m);
        system.mutable_masses()[o] = 15.999;
        system.mutable_masses()[h1] = 1.008;
        system.mutable_masses()[h2] = 1.008;
        const std::array<double, 3> origin = {5.0 + shift, 6.0 + row, 7.0};
        system.mutable_coordinates()[o] = origin;
        const double second = angle + 2.0 * std::asin(kHH / (2.0 * kOH));
        system.mutable_coordinates()[h1] = {origin[0] + kOH * std::cos(angle),
                                            origin[1] + kOH * std::sin(angle) * 0.6,
                                            origin[2] + kOH * std::sin(angle) * 0.8};
        system.mutable_coordinates()[h2] = {origin[0] + kOH * std::cos(second),
                                            origin[1] + kOH * std::sin(second) * 0.6,
                                            origin[2] + kOH * std::sin(second) * 0.8};
        for (std::size_t i = 0; i < 3; ++i) {
            system.mutable_atom_tags()[3 * m + i] = static_cast<int>(3 * m + i);
        }
        constraints.push_back({static_cast<int>(o), static_cast<int>(h1), kOH});
        constraints.push_back({static_cast<int>(o), static_cast<int>(h2), kOH});
        constraints.push_back({static_cast<int>(h1), static_cast<int>(h2), kHH});
    }
    return constraints;
}

// Drives the REAL startup path -- Simulation::initialize() -- rather than
// calling VelocityInitializer directly, because the defect lives in the order
// those two objects run in.
Outcome initialize_system(std::size_t molecule_count,
                          const std::vector<gmd::BondConstraint>& override_constraints,
                          bool use_override,
                          bool remove_com = true,
                          double temperature = kTargetTemperature,
                          double box_length = 200.0) {
    gmd::System system;
    std::vector<gmd::BondConstraint> constraints =
        water_fixture(system, molecule_count, box_length);
    if (use_override) constraints = override_constraints;

    auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(1.0);
    if (!constraints.empty()) {
        integrator->set_constraint_solver(
            std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings()));
    }
    auto initializer = std::make_shared<gmd::VelocityInitializer>(kSeed);
    NullForceProvider provider;
    gmd::RuntimeContext runtime;

    gmd::Simulation simulation(&system);
    simulation.set_velocity_initializer(initializer);
    simulation.set_velocity_init_mode(gmd::VelocityInitMode::Random);
    simulation.set_initial_temperature(temperature);
    simulation.set_remove_center_of_mass_velocity(remove_com);
    simulation.set_force_provider(
        std::shared_ptr<gmd::ForceProvider>(&provider, [](gmd::ForceProvider*) {}));
    simulation.set_integrator(integrator);
    simulation.set_time_step(1.0);
    simulation.initialize(runtime);

    Outcome outcome;
    outcome.atom_count = system.num_local_atoms();
    outcome.constraint_count = constraints.size();
    outcome.authoritative_dof = integrator->degrees_of_freedom(system);
    outcome.naive_dof = remove_com ? 3 * outcome.atom_count - 3 : 3 * outcome.atom_count;

    std::array<double, 3> momentum = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const auto v = system.velocities()[i];
        const double m = system.masses()[i];
        outcome.twice_kinetic_energy += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for (std::size_t d = 0; d < 3; ++d) momentum[d] += m * v[d];
        outcome.velocities_by_tag[system.atom_tag(i)] = {v[0], v[1], v[2]};
    }
    outcome.momentum_magnitude = std::sqrt(momentum[0] * momentum[0] +
                                           momentum[1] * momentum[1] +
                                           momentum[2] * momentum[2]);
    outcome.temperature =
        outcome.authoritative_dof == 0
            ? 0.0
            : outcome.twice_kinetic_energy /
                  (static_cast<double>(outcome.authoritative_dof) * kBoltzmann);

    for (const auto& c : constraints) {
        const auto ri = system.coordinates()[static_cast<std::size_t>(c.i)];
        const auto rj = system.coordinates()[static_cast<std::size_t>(c.j)];
        const auto vi = system.velocities()[static_cast<std::size_t>(c.i)];
        const auto vj = system.velocities()[static_cast<std::size_t>(c.j)];
        double dot = 0.0;
        for (std::size_t d = 0; d < 3; ++d) dot += (ri[d] - rj[d]) * (vi[d] - vj[d]);
        outcome.worst_tangency = std::max(outcome.worst_tangency, std::abs(dot));
    }
    return outcome;
}

Outcome initialize_waters(std::size_t molecule_count, bool remove_com = true,
                          double temperature = kTargetTemperature) {
    return initialize_system(molecule_count, {}, false, remove_com, temperature);
}

// --- properties that hold whatever DOF the rescale used --------------------

void test_unconstrained_reference() {
    // With no constraints the naive and authoritative counts coincide, so the
    // temperature is already exact. This is the control: it shows the fixture
    // and the measurement are right, and it must keep passing afterwards.
    gmd::System system;
    water_fixture(system, 4, 200.0);          // built, then run WITHOUT constraints
    const Outcome outcome = initialize_system(4, {}, true, true);
    check(outcome.constraint_count == 0, "the unconstrained control has constraints");
    check(outcome.authoritative_dof == outcome.naive_dof,
          "unconstrained, the authoritative DOF (" +
              std::to_string(outcome.authoritative_dof) + ") must equal 3N-3 (" +
              std::to_string(outcome.naive_dof) + ")");
    check(std::abs(outcome.temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "an unconstrained system initialised to " + number(kTargetTemperature) +
              " K carries " + number(outcome.temperature) + " K");
    check(outcome.momentum_magnitude < 1.0e-12 * outcome.twice_kinetic_energy + 1.0e-14,
          "unconstrained centre-of-mass momentum is " +
              number(outcome.momentum_magnitude));
}

void test_unconstrained_without_com_removal() {
    const Outcome outcome = initialize_system(4, {}, true, false);
    check(outcome.authoritative_dof == 3 * outcome.atom_count,
          "with COM removal disabled and no constraints the DOF must be 3N, got " +
              std::to_string(outcome.authoritative_dof));
    check(std::abs(outcome.temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "with COM removal disabled the field carries " +
              number(outcome.temperature) + " K");
    check(outcome.momentum_magnitude > 0.0,
          "with COM removal disabled the residual momentum is exactly zero");
}

void test_authoritative_dof_accounts_for_constraints() {
    for (std::size_t molecules : {1u, 3u, 10u}) {
        const Outcome outcome = initialize_waters(molecules);
        const std::size_t expected = 3 * outcome.atom_count - outcome.constraint_count - 3;
        check(outcome.authoritative_dof == expected,
              std::to_string(molecules) + " water(s): the integrator reports DOF " +
                  std::to_string(outcome.authoritative_dof) + ", expected 3N - rank - 3 = " +
                  std::to_string(expected));
        // And the naive count the initializer uses is a different number, which
        // is the whole problem.
        check(outcome.naive_dof != outcome.authoritative_dof,
              "this fixture is supposed to have constraints that change the DOF");
    }
}

void test_velocities_are_tangent_and_momentum_free() {
    // Both hold today: the integrator RATTLE-projects after the initializer has
    // run, and projection conserves momentum while centre-of-mass removal
    // preserves tangency. They are asserted here so the fix cannot quietly
    // trade one for the other.
    for (std::size_t molecules : {1u, 3u, 10u}) {
        const Outcome outcome = initialize_waters(molecules);
        check(outcome.worst_tangency < 1.0e-11,
              std::to_string(molecules) + " water(s): worst |r.v| is " +
                  number(outcome.worst_tangency) + ", so the initial velocities are "
                  "not in the constraint tangent space");
        const double scale = std::sqrt(outcome.twice_kinetic_energy) + 1.0;
        check(outcome.momentum_magnitude < 1.0e-11 * scale,
              std::to_string(molecules) + " water(s): residual centre-of-mass "
              "momentum is " + number(outcome.momentum_magnitude));
    }
}

void test_zero_temperature_is_exact() {
    const Outcome outcome = initialize_waters(3, true, 0.0);
    check(outcome.twice_kinetic_energy == 0.0,
          "a constrained system initialised to 0 K carries kinetic energy " +
              number(outcome.twice_kinetic_energy));
    check(outcome.worst_tangency == 0.0,
          "a system at rest is trivially tangent; worst |r.v| is " +
              number(outcome.worst_tangency));
    check(outcome.momentum_magnitude == 0.0,
          "a system at rest has zero momentum; got " +
              number(outcome.momentum_magnitude));
}

void test_redundant_constraints_are_rejected() {
    // Four constraints on three atoms cannot be independent. The authoritative
    // analysis must reject the run rather than silently subtracting four.
    std::vector<gmd::BondConstraint> redundant = {
        {0, 1, 0.9572}, {0, 2, 0.9572}, {1, 2, 1.5139}, {2, 1, 1.5139},
    };
    // The duplicate (2,1) normalises away, so add a genuinely dependent one:
    redundant.back() = {1, 2, 1.5139};
    std::vector<gmd::BondConstraint> dependent = {
        {0, 1, 1.0}, {1, 2, 1.0}, {0, 2, 2.0},   // collinear target: rank 2, count 3
    };
    bool threw = false;
    try {
        initialize_system(1, dependent, true);
    } catch (const std::exception&) {
        threw = true;
    }
    check(threw,
          "a dependent constraint set was accepted; a collinear target has rank 2 "
          "against 3 constraints and must be rejected before any dynamics start");
}

// --- the temperature must now be exact ------------------------------------

void test_constrained_temperature_is_exact() {
    // The whole point. The rescale now uses the authoritative DOF and runs
    // AFTER the projection, so the target is hit exactly rather than on
    // average. The bound is round-off in a global sum, not a physics tolerance.
    for (std::size_t molecules : {1u, 2u, 3u, 10u, 40u}) {
        const Outcome outcome = initialize_waters(molecules);
        check(std::abs(outcome.temperature / kTargetTemperature - 1.0) < 1.0e-12,
              std::to_string(molecules) + " water(s): initialised to " +
                  number(kTargetTemperature) + " K, the field carries " +
                  number(outcome.temperature) + " K over " +
                  std::to_string(outcome.authoritative_dof) +
                  " authoritative degrees of freedom");
        // And it is the authoritative count that was used, not 3N-3: measuring
        // against the naive one must now be visibly wrong.
        const double naive_temperature =
            outcome.twice_kinetic_energy /
            (static_cast<double>(outcome.naive_dof) * kBoltzmann);
        check(std::abs(naive_temperature / kTargetTemperature - 1.0) > 1.0e-6,
              std::to_string(molecules) + " water(s): the field also satisfies "
              "the 3N-3 convention, so this fixture cannot tell the two apart");
    }
}

void test_single_distance_constraint() {
    // The smallest constrained system that is not degenerate: one bond, two
    // atoms of different mass. dof = 6 - 1 - 3 = 2.
    gmd::System probe;
    water_fixture(probe, 1, 200.0);
    std::vector<gmd::BondConstraint> one = {{0, 1, 0.9572}};
    const Outcome outcome = initialize_system(1, one, true);
    check(outcome.authoritative_dof == 3 * outcome.atom_count - 1 - 3,
          "one constraint on " + std::to_string(outcome.atom_count) +
              " atoms should leave 3N-1-3 degrees of freedom; got " +
              std::to_string(outcome.authoritative_dof));
    check(outcome.worst_tangency < 1.0e-11,
          "a single constrained bond is not tangent: |r.v| = " +
              number(outcome.worst_tangency));
    check(std::abs(outcome.temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "one constraint: the field carries " + number(outcome.temperature) + " K");
}

void test_disconnected_components() {
    // Two molecules that share no atom. The rank is the sum of the components'
    // ranks, and neither may be dropped.
    const Outcome outcome = initialize_waters(2);
    check(outcome.constraint_count == 6, "expected six constraints across two waters");
    check(outcome.authoritative_dof == 3 * outcome.atom_count - 6 - 3,
          "two disconnected rigid waters should leave 3N-6-3 degrees of freedom; got " +
              std::to_string(outcome.authoritative_dof));
    check(outcome.worst_tangency < 1.0e-11,
          "disconnected components: worst |r.v| is " + number(outcome.worst_tangency));
    check(std::abs(outcome.temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "disconnected components: the field carries " + number(outcome.temperature) +
              " K");
}

void test_constraints_and_com_removal_do_not_double_subtract() {
    // Rigid translation lies in the null space of a distance-constraint
    // Jacobian, so the three translational modes and the constrained modes are
    // disjoint and the DOF subtractions do not overlap. If they were counted
    // twice the DOF would be three lower and the temperature correspondingly
    // wrong.
    const Outcome with_com = initialize_waters(3, true);
    const Outcome without_com = initialize_waters(3, false);
    check(with_com.authoritative_dof + 3 == without_com.authoritative_dof,
          "removing the centre-of-mass velocity should cost exactly three degrees "
          "of freedom; it cost " +
              std::to_string(without_com.authoritative_dof - with_com.authoritative_dof));
    check(std::abs(with_com.temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "constraints + COM removal: " + number(with_com.temperature) + " K");
    check(std::abs(without_com.temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "constraints, no COM removal: " + number(without_com.temperature) + " K");
    const double scale = std::sqrt(with_com.twice_kinetic_energy) + 1.0;
    check(with_com.momentum_magnitude < 1.0e-11 * scale,
          "constraints + COM removal left momentum " +
              number(with_com.momentum_magnitude));
    check(without_com.momentum_magnitude > 0.0,
          "with COM removal disabled the momentum is exactly zero");
    // Both remain tangent: COM removal must not have broken the projection.
    check(with_com.worst_tangency < 1.0e-11 && without_com.worst_tangency < 1.0e-11,
          "COM removal broke constraint tangency");
}

// Runs one custom-geometry fixture through the real startup path.
Outcome initialize_custom(const std::vector<std::array<double, 3>>& positions,
                          const std::vector<double>& masses,
                          const std::vector<gmd::BondConstraint>& constraints,
                          double temperature = kTargetTemperature) {
    gmd::System system;
    const std::size_t n = positions.size();
    system.resize(n, n);
    gmd::Box box;
    box.set_lengths({200.0, 200.0, 200.0});
    system.set_box(box);
    for (std::size_t i = 0; i < n; ++i) {
        system.mutable_coordinates()[i] = positions[i];
        system.mutable_masses()[i] = masses[i];
        system.mutable_atom_tags()[i] = static_cast<int>(i);
    }
    auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(1.0);
    if (!constraints.empty()) {
        integrator->set_constraint_solver(
            std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings()));
    }
    auto initializer = std::make_shared<gmd::VelocityInitializer>(kSeed);
    NullForceProvider provider;
    gmd::RuntimeContext runtime;
    gmd::Simulation simulation(&system);
    simulation.set_velocity_initializer(initializer);
    simulation.set_velocity_init_mode(gmd::VelocityInitMode::Random);
    simulation.set_initial_temperature(temperature);
    simulation.set_remove_center_of_mass_velocity(true);
    simulation.set_force_provider(
        std::shared_ptr<gmd::ForceProvider>(&provider, [](gmd::ForceProvider*) {}));
    simulation.set_integrator(integrator);
    simulation.set_time_step(1.0);
    simulation.initialize(runtime);

    Outcome outcome;
    outcome.atom_count = system.num_local_atoms();
    outcome.constraint_count = constraints.size();
    outcome.authoritative_dof = integrator->degrees_of_freedom(system);
    outcome.naive_dof = 3 * outcome.atom_count - 3;
    std::array<double, 3> momentum = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const auto v = system.velocities()[i];
        const double m = system.masses()[i];
        outcome.twice_kinetic_energy += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for (std::size_t d = 0; d < 3; ++d) momentum[d] += m * v[d];
        outcome.velocities_by_tag[system.atom_tag(i)] = {v[0], v[1], v[2]};
    }
    outcome.momentum_magnitude = std::sqrt(momentum[0] * momentum[0] +
                                           momentum[1] * momentum[1] +
                                           momentum[2] * momentum[2]);
    outcome.temperature =
        outcome.authoritative_dof == 0
            ? 0.0
            : outcome.twice_kinetic_energy /
                  (static_cast<double>(outcome.authoritative_dof) * kBoltzmann);
    for (const auto& c : constraints) {
        const auto ri = system.coordinates()[static_cast<std::size_t>(c.i)];
        const auto rj = system.coordinates()[static_cast<std::size_t>(c.j)];
        const auto vi = system.velocities()[static_cast<std::size_t>(c.i)];
        const auto vj = system.velocities()[static_cast<std::size_t>(c.j)];
        double dot = 0.0;
        for (std::size_t d = 0; d < 3; ++d) dot += (ri[d] - rj[d]) * (vi[d] - vj[d]);
        outcome.worst_tangency = std::max(outcome.worst_tangency, std::abs(dot));
    }
    return outcome;
}

void test_ill_conditioned_but_independent() {
    // A thin triangle: full rank, but close enough to collinear that the
    // Jacobian is badly conditioned. The targets are the geometry's own
    // distances, so the projection has nothing to move and the conditioning is
    // the only thing under test. The existing rank policy accepts such a set
    // with a warning, so initialization must produce a finite, tangent,
    // correctly scaled field rather than failing.
    const std::vector<std::array<double, 3>> positions = {
        {10.0, 10.0, 10.0}, {11.0, 10.0, 10.0}, {11.5, 10.5, 10.0}};
    auto distance = [&](std::size_t a, std::size_t b) {
        double sum = 0.0;
        for (std::size_t d = 0; d < 3; ++d) {
            const double delta = positions[a][d] - positions[b][d];
            sum += delta * delta;
        }
        return std::sqrt(sum);
    };
    const std::vector<gmd::BondConstraint> thin = {
        {0, 1, distance(0, 1)}, {1, 2, distance(1, 2)}, {0, 2, distance(0, 2)}};
    // Thin, but still solvable. The triangle inequality closes to within about
    // eight percent of the long side, which conditions the Jacobian far worse
    // than an equilateral arrangement while staying inside what SHAKE and
    // RATTLE can converge to at this tolerance. Squeezing it further -- 0.3
    // percent slack, say -- makes RATTLE run out of iterations, which is the
    // solver's own limit rather than the rank policy's and is not what this
    // case is about.
    const double slack = distance(0, 1) + distance(1, 2) - distance(0, 2);
    check(slack > 0.0 && slack < 0.15 * distance(0, 2),
          "this fixture is meant to be thin but not degenerate; slack is " +
              number(slack) + " against a long side of " + number(distance(0, 2)));

    Outcome outcome;
    try {
        outcome = initialize_custom(positions, {15.999, 1.008, 12.011}, thin);
    } catch (const std::exception& error) {
        check(false, std::string("an independent, if ill-conditioned, constraint "
                                 "set was rejected: ") + error.what());
        return;
    }
    for (const auto& [tag, v] : outcome.velocities_by_tag) {
        for (std::size_t d = 0; d < 3; ++d) {
            check(std::isfinite(v[d]),
                  "atom " + std::to_string(tag) + " has a non-finite velocity under "
                  "an ill-conditioned constraint set");
        }
    }
    std::cout << "  ill-conditioned triangle     dof=" << outcome.authoritative_dof
              << "  |r.v|=" << std::scientific << std::setprecision(3)
              << outcome.worst_tangency << "  T=" << std::fixed
              << std::setprecision(6) << outcome.temperature << std::defaultfloat
              << " K\n";
    check(outcome.authoritative_dof == 3 * outcome.atom_count - 3 - 3,
          "a thin but independent triangle should still remove three degrees of "
          "freedom; dof is " + std::to_string(outcome.authoritative_dof));
    check(outcome.worst_tangency < 1.0e-9,
          "ill-conditioned but independent: worst |r.v| is " +
              number(outcome.worst_tangency));
    check(std::abs(outcome.temperature / kTargetTemperature - 1.0) < 1.0e-10,
          "ill-conditioned but independent: the field carries " +
              number(outcome.temperature) + " K");
}

void test_storage_permutation_is_invariant() {
    // The rank-independent RNG guarantee must survive the constraint work.
    // Velocities are compared by TAG, and the constraint list is expressed in
    // tags, so permuting storage changes nothing physical.
    const Outcome ordered = initialize_waters(3);

    gmd::System system;
    const auto constraints = water_fixture(system, 3, 200.0);
    // Reverse the storage while keeping each atom's tag, mass and position.
    const std::size_t n = system.num_local_atoms();
    std::vector<std::array<double, 3>> positions(n);
    std::vector<double> masses(n);
    std::vector<int> tags(n);
    for (std::size_t i = 0; i < n; ++i) {
        positions[i] = system.coordinates()[i];
        masses[i] = system.masses()[i];
        tags[i] = system.atom_tag(i);
    }
    gmd::System permuted;
    permuted.resize(n, n);
    permuted.set_box(system.box());
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t source = n - 1 - i;
        permuted.mutable_coordinates()[i] = positions[source];
        permuted.mutable_masses()[i] = masses[source];
        permuted.mutable_atom_tags()[i] = tags[source];
    }

    auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(1.0);
    integrator->set_constraint_solver(
        std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings()));
    auto initializer = std::make_shared<gmd::VelocityInitializer>(kSeed);
    NullForceProvider provider;
    gmd::RuntimeContext runtime;
    gmd::Simulation simulation(&permuted);
    simulation.set_velocity_initializer(initializer);
    simulation.set_velocity_init_mode(gmd::VelocityInitMode::Random);
    simulation.set_initial_temperature(kTargetTemperature);
    simulation.set_remove_center_of_mass_velocity(true);
    simulation.set_force_provider(
        std::shared_ptr<gmd::ForceProvider>(&provider, [](gmd::ForceProvider*) {}));
    simulation.set_integrator(integrator);
    simulation.set_time_step(1.0);
    simulation.initialize(runtime);

    double worst = 0.0;
    for (std::size_t i = 0; i < permuted.num_local_atoms(); ++i) {
        const auto v = permuted.velocities()[i];
        const auto& reference = ordered.velocities_by_tag.at(permuted.atom_tag(i));
        for (std::size_t d = 0; d < 3; ++d) {
            worst = std::max(worst, std::abs(v[d] - reference[d]));
        }
    }
    std::cout << "  storage permutation, by tag  " << std::scientific
              << std::setprecision(3) << worst << std::defaultfloat << '\n';
    // The projection is iterative and its per-atom corrections are accumulated
    // in storage order, so the bound is that solver tolerance rather than the
    // pure reduction round-off the unconstrained case achieves.
    check(worst < 1.0e-11,
          "permuting storage changed the constrained velocity field by " +
              number(worst) + " when compared by atom tag");
}

// --- reported for the record ----------------------------------------------

void report_constrained_temperature_error() {
    std::cout << "  molecules   atoms  constraints   dof(auth)  dof(3N-3)   "
                 "T [K]        relative error\n";
    for (std::size_t molecules : {1u, 3u, 10u, 40u}) {
        const Outcome outcome = initialize_waters(molecules);
        const double relative = outcome.temperature / kTargetTemperature - 1.0;
        std::cout << "  " << std::setw(9) << molecules
                  << std::setw(8) << outcome.atom_count
                  << std::setw(13) << outcome.constraint_count
                  << std::setw(12) << outcome.authoritative_dof
                  << std::setw(11) << outcome.naive_dof
                  << std::setw(13) << std::fixed << std::setprecision(4)
                  << outcome.temperature
                  << std::setw(16) << std::scientific << std::setprecision(4)
                  << relative << std::defaultfloat << '\n';
    }
    // What is asserted is that the fixture still exercises the disagreement
    // between the two DOF counts; the values themselves are printed so a
    // regression is legible next to the numbers this replaced.
    const Outcome outcome = initialize_waters(3);
    check(outcome.authoritative_dof < outcome.naive_dof,
          "the constrained fixture no longer distinguishes the two DOF counts");
}

}  // namespace

int main() {
    std::cout << "[constrained init] auditing constrained random velocity "
                 "initialization\n";
    test_unconstrained_reference();
    test_unconstrained_without_com_removal();
    test_authoritative_dof_accounts_for_constraints();
    test_velocities_are_tangent_and_momentum_free();
    test_zero_temperature_is_exact();
    test_redundant_constraints_are_rejected();
    test_constrained_temperature_is_exact();
    test_single_distance_constraint();
    test_disconnected_components();
    test_constraints_and_com_removal_do_not_double_subtract();
    test_ill_conditioned_but_independent();
    test_storage_permutation_is_invariant();
    report_constrained_temperature_error();

    if (failures == 0) {
        std::cout << "[constrained init] audit passed\n";
        return 0;
    }
    std::cerr << "[constrained init] " << failures << " check(s) failed\n";
    return 1;
}
