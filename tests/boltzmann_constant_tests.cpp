// Audit of the Boltzmann constant across every production temperature path.
//
// GMD works in eV, Angstrom and kelvin, so k_B is the eV-per-kelvin conversion
// that turns a kinetic energy into a temperature. Five production paths use it:
// Maxwell-Boltzmann velocity sampling, the rescale-to-target that follows it,
// instantaneous temperature reporting, the Nose-Hoover thermostat mass and
// friction force, and the Monte Carlo barostat's Metropolis exponent. They must
// all use the same number, or a system initialised to 300 K does not report
// 300 K.
//
// HOW THE CONSTANT IS MEASURED RATHER THAN ASSUMED
//
// Nothing here reads a production constant directly: the initializer's lives in
// an anonymous namespace in a .cpp, and the barostat's is a private class
// member. Each path is instead driven with known inputs and the constant
// recovered from its output:
//
//   initializer       initialize() ends in a rescale that forces
//                     2K / (dof * k_B) == T, so k_B = 2K / (dof * T)
//   reporting         temperature_from_twice_ke() inverts the same relation
//   Nose-Hoover       the thermostat mass is Q = dof * k_B * T * tau^2, and Q
//                     is exposed at 17 digits through checkpoint_state()
//   rescaling         the thermostat drives 2K / (dof * k_B) to the target
//
// A path that omitted k_B would measure 1; one that applied it twice would
// measure k_B^2; one that inverted the conversion would measure 1/k_B. Each is
// a named failure below rather than a generic tolerance.
//
// THE AUTHORITATIVE VALUE
//
// Unusually, k_B in eV/K is EXACT and has no uncertainty at all. Since the 2019
// SI redefinition both of its ingredients are exact by definition:
//
//   k_B = 1.380649e-23 J/K          https://physics.nist.gov/cgi-bin/cuu/Value?k
//   e   = 1.602176634e-19 C         https://physics.nist.gov/cgi-bin/cuu/Value?e
//
// and 1 eV is e joules exactly, so
//
//   k_B [eV/K] = 1380649 / 16021766340 = 8.6173332621451774336...e-5
//
// which is a non-terminating rational: every decimal literal is a truncation of
// it. NIST tabulates it as "8.617 333 262... x 10^-5 eV/K, exact" -- the
// ellipsis is their notation for exactly this.
// (https://physics.nist.gov/cgi-bin/cuu/Value?kev, CODATA 2022.)
//
// The reference below is written as those two SI literals and divided, rather
// than as one pre-computed number, so the derivation is visible and checkable
// instead of copied.

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/mc_barostat.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/integrator/velocity_rescaling_thermostat.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[k_B audit] " << message << '\n';
        ++failures;
    }
}

std::string number(double value, int digits = 17) {
    std::ostringstream out;
    out << std::setprecision(digits) << value;
    return out.str();
}

// --- the authoritative value, derived here from the exact SI definitions ----

constexpr double kBoltzmannJoulesPerKelvin = 1.380649e-23;      // exact, SI 2019
constexpr double kElementaryChargeCoulombs = 1.602176634e-19;   // exact, SI 2019
constexpr double kReferenceBoltzmann =
    kBoltzmannJoulesPerKelvin / kElementaryChargeCoulombs;

// What every production path is expected to use: gmd::kBoltzmannConstantEVPerKelvin.
// One pin, because there is now one constant. Both names are kept so the
// assertions below still read as "initialization agrees with reporting", which
// is the property that was broken.
constexpr double kExpectedShared      = 8.617333262145177e-5;
constexpr double kExpectedInitializer = kExpectedShared;

// The rounding policy: the production literal must reproduce the exact rational
// k_B/e to within double-precision representation. There is no physical
// uncertainty to hide behind here -- both ingredients are exact by SI
// definition -- so any looser literal is an arbitrary truncation. For scale,
// the values this replaced were 1.7e-11 and 1.13e-06 low.
constexpr double kRoundingPolicyRelative = 1.0e-15;

static_assert(kExpectedShared / kReferenceBoltzmann - 1.0 < kRoundingPolicyRelative &&
              1.0 - kExpectedShared / kReferenceBoltzmann < kRoundingPolicyRelative,
              "the pinned production constant no longer matches k_B/e derived "
              "from the exact SI definitions to the documented rounding policy");

// Round-off budget. Every measurement below is a ratio of quantities the
// engine computes in double from the same inputs, so agreement is limited by
// summation order, not by the constant. 1e-12 is several orders above the
// observed spread and six below the 1.13e-06 split under audit.
constexpr double kRelativeTolerance = 1.0e-12;

bool close(double a, double b, double tolerance = kRelativeTolerance) {
    return std::fabs(a / b - 1.0) <= tolerance;
}

// --- fixtures ---------------------------------------------------------------
//
// Deliberately unequal masses and an irregular geometry: a fixture with equal
// masses cannot distinguish a per-atom mass factor from a global one, and one
// with symmetric velocities can hide a sign or axis error in the kinetic sum.

gmd::System make_system(std::size_t atom_count) {
    gmd::System system;
    system.resize(atom_count, atom_count);
    gmd::Box box;
    box.set_lengths({24.0, 28.0, 32.0});
    system.set_box(box);
    // Masses spanning a factor of ~12, no two alike.
    const std::array<double, 8> masses = {1.008, 12.011, 15.999, 4.0026,
                                          6.941, 10.811, 14.007, 3.016};
    for (std::size_t i = 0; i < atom_count; ++i) {
        system.mutable_masses()[i] = masses[i % masses.size()] + 0.1 * static_cast<double>(i);
        system.mutable_coordinates()[i] = {2.0 + 2.7 * static_cast<double>(i),
                                           3.0 + 1.9 * static_cast<double>(i % 5),
                                           4.0 + 2.3 * static_cast<double>(i % 3)};
        system.mutable_atom_tags()[i] = static_cast<int>(i);
        system.mutable_atom_owners()[i] = 0;
    }
    return system;
}

// The initializer computes its own degrees of freedom as 3N-3 when the centre
// of mass has been removed. Reproduced here rather than imported so the two
// cannot drift silently.
std::size_t reference_dof(std::size_t atom_count, bool com_removed) {
    return com_removed ? 3 * atom_count - 3 : 3 * atom_count;
}

double twice_kinetic_energy(const gmd::System& system) {
    const auto masses = system.masses();
    const auto velocities = system.velocities();
    double total = 0.0;
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const auto& v = velocities[i];
        total += masses[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    return total;
}

// --- measurements -----------------------------------------------------------

constexpr double kTargetTemperature = 300.0;
constexpr std::size_t kAtomCount = 8;

// initialize() always ends in rescale_temperature(), which forces
// 2K / (dof * k_B) == T exactly. So k_B falls straight out.
double measure_initializer_constant(double target = kTargetTemperature,
                                    bool com_removed = true) {
    gmd::System system = make_system(kAtomCount);
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, target, gmd::VelocityInitMode::Random, com_removed);
    const double twice_ke = twice_kinetic_energy(system);
    return twice_ke / (static_cast<double>(reference_dof(kAtomCount, com_removed)) * target);
}

// temperature_from_twice_ke() is the inverse relation, so driving it with a
// known 2K and dof recovers the constant it uses.
double measure_reporting_constant() {
    gmd::System system = make_system(kAtomCount);
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, true);
    const double twice_ke = gmd::compute_twice_ke(system);
    const std::size_t dof = reference_dof(kAtomCount, true);
    const double reported = gmd::temperature_from_twice_ke(twice_ke, dof);
    return twice_ke / (static_cast<double>(dof) * reported);
}

double parse_checkpoint_field(const std::string& state, const std::string& field) {
    std::istringstream input(state);
    std::string key;
    double value = 0.0;
    while (input >> key) {
        if (!(input >> value)) break;
        if (key == field) return value;
    }
    throw std::runtime_error("field '" + field + "' not found in: " + state);
}

// Q = dof * k_B * T * tau^2, exposed at 17 digits through checkpoint_state().
double measure_nose_hoover_constant(double target = kTargetTemperature,
                                    double tau = 40.0,
                                    std::size_t atom_count = kAtomCount) {
    gmd::System system = make_system(atom_count);
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, target, gmd::VelocityInitMode::Random, true);

    gmd::NoseHooverThermostat thermostat(tau);
    thermostat.initialize(system);
    const std::size_t dof = thermostat.degrees_of_freedom();
    // One half-kick is what makes the thermostat compute Q lazily.
    thermostat.apply_half_kick(system, 0.5, target);
    const double q = parse_checkpoint_field(thermostat.checkpoint_state(), "Q");
    return q / (static_cast<double>(dof) * target * tau * tau);
}

// The rescaling thermostat drives 2K / (dof * k_B) to the target, so after one
// application the constant follows from the kinetic energy it left behind.
double measure_velocity_rescaling_constant(double target = kTargetTemperature) {
    gmd::System system = make_system(kAtomCount);
    gmd::VelocityInitializer initializer(20260830u);
    // Start away from the target so the thermostat actually has work to do.
    initializer.initialize(system, target * 2.5, gmd::VelocityInitMode::Random, true);

    gmd::VelocityRescalingThermostat thermostat;
    thermostat.initialize(system);
    const std::size_t dof = thermostat.degrees_of_freedom();
    thermostat.apply(system, 1.0, target);
    const double twice_ke = gmd::compute_twice_ke(system);
    return twice_ke / (static_cast<double>(dof) * target);
}

// --- tests ------------------------------------------------------------------

void test_production_matches_the_authoritative_value() {
    const double relative = kExpectedShared / kReferenceBoltzmann - 1.0;
    std::cout << "\n  k_B from exact SI definitions = " << number(kReferenceBoltzmann)
              << " eV/K\n";
    std::cout << "  production                   = " << number(kExpectedShared)
              << "        relative deviation " << number(relative, 4) << '\n';
    check(std::fabs(relative) <= kRoundingPolicyRelative,
          "the production constant deviates from k_B/e by " + number(relative, 4) +
              ", outside the documented rounding policy of " +
              number(kRoundingPolicyRelative, 3));
}

void test_reference_matches_the_authoritative_value() {
    // NIST tabulates k_B as 8.617 333 262... x 10^-5 eV/K, exact. The quoted
    // digits are a truncation of the exact rational, so the derived reference
    // must agree with them to the last digit quoted and no further.
    constexpr double quoted = 8.617333262e-5;
    std::cout << "\n  k_B from exact SI definitions = " << number(kReferenceBoltzmann)
              << " eV/K\n";
    std::cout << "  NIST quoted digits            = " << number(quoted, 10)
              << "...   agreement " << number(kReferenceBoltzmann / quoted - 1.0, 4)
              << '\n';
    check(std::fabs(kReferenceBoltzmann / quoted - 1.0) < 5.0e-11,
          "the reference derived from k_B/e disagrees with NIST's quoted digits by " +
              number(kReferenceBoltzmann / quoted - 1.0, 4));
}

void test_known_velocities_give_exact_temperature() {
    // Two atoms, hand-chosen masses and velocities, all three components
    // different and non-zero, so an axis swap or a dropped term cannot pass.
    gmd::System system;
    system.resize(2, 2);
    gmd::Box box;
    box.set_lengths({20.0, 20.0, 20.0});
    system.set_box(box);
    system.mutable_masses()[0] = 2.0;
    system.mutable_masses()[1] = 5.0;
    system.mutable_velocities()[0] = {0.25, -0.5, 0.75};
    system.mutable_velocities()[1] = {-0.125, 0.375, -0.625};

    // 2K = m0 v0^2 + m1 v1^2, computed here independently.
    const double expected_twice_ke =
        2.0 * (0.25 * 0.25 + 0.5 * 0.5 + 0.75 * 0.75) +
        5.0 * (0.125 * 0.125 + 0.375 * 0.375 + 0.625 * 0.625);
    const double twice_ke = gmd::compute_twice_ke(system);
    check(close(twice_ke, expected_twice_ke, 1.0e-15),
          "compute_twice_ke() returned " + number(twice_ke) + ", expected " +
              number(expected_twice_ke));

    for (const std::size_t dof : {std::size_t{3}, std::size_t{5}, std::size_t{6}}) {
        const double reported = gmd::temperature_from_twice_ke(twice_ke, dof);
        const double expected = expected_twice_ke /
                                (static_cast<double>(dof) * kExpectedShared);
        check(close(reported, expected, 1.0e-13),
              "temperature at dof=" + std::to_string(dof) + " is " + number(reported) +
                  " K, expected " + number(expected) + " K");
    }
    check(gmd::temperature_from_twice_ke(twice_ke, 0) == 0.0,
          "zero degrees of freedom must report 0 K rather than dividing by zero");
}

void test_paths_are_mutually_consistent() {
    const double initializer = measure_initializer_constant();
    const double reporting   = measure_reporting_constant();
    const double nose_hoover = measure_nose_hoover_constant();
    const double rescaling   = measure_velocity_rescaling_constant();

    std::cout << "\n  Constant recovered from each production path\n";
    std::cout << "    velocity initialization       " << number(initializer) << '\n';
    std::cout << "    temperature reporting         " << number(reporting) << '\n';
    std::cout << "    Nose-Hoover thermostat mass   " << number(nose_hoover) << '\n';
    std::cout << "    velocity-rescaling thermostat " << number(rescaling) << '\n';

    check(close(initializer, kExpectedInitializer),
          "velocity initialization measures k_B = " + number(initializer) +
              ", expected " + number(kExpectedInitializer));
    for (const auto& [label, value] :
         {std::pair{"temperature reporting", reporting},
          std::pair{"Nose-Hoover thermostat mass", nose_hoover},
          std::pair{"velocity-rescaling thermostat", rescaling}}) {
        check(close(value, kExpectedShared),
              std::string(label) + " measures k_B = " + number(value) +
                  ", expected " + number(kExpectedShared));
    }

    // The property that was broken: initialization and reporting must agree.
    // Asserted directly rather than via each matching the pin, so that two
    // paths drifting together could not pass.
    check(close(initializer, reporting),
          "velocity initialization measures k_B = " + number(initializer) +
              " but temperature reporting measures " + number(reporting) +
              " (relative " + number(initializer / reporting - 1.0, 4) +
              "); a system initialised to a target temperature would not report "
              "that temperature back");
    for (const auto& [label, value] :
         {std::pair{"Nose-Hoover thermostat mass", nose_hoover},
          std::pair{"velocity-rescaling thermostat", rescaling}}) {
        check(close(value, initializer),
              std::string(label) + " measures k_B = " + number(value) +
                  " but velocity initialization measures " + number(initializer) +
                  "; two production paths are using different constants");
    }
}

void test_initialize_then_report_round_trip() {
    // The observable consequence of the split: ask for 300 K, read back the
    // temperature the engine itself reports.
    gmd::System system = make_system(kAtomCount);
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, true);
    const double reported = gmd::temperature_from_twice_ke(
        gmd::compute_twice_ke(system), reference_dof(kAtomCount, true));
    const double error = reported / kTargetTemperature - 1.0;
    std::cout << "\n  Initialised to " << number(kTargetTemperature, 6)
              << " K, engine reports " << number(reported, 12)
              << " K   (relative " << number(error, 4) << ")\n";
    check(std::fabs(error) <= kRelativeTolerance,
          "a system initialised to " + number(kTargetTemperature, 6) +
              " K reports " + number(reported, 12) + " K, a relative error of " +
              number(error, 4) + ". Initialization and reporting are not using "
              "the same Boltzmann constant.");
}

void test_zero_temperature_and_degrees_of_freedom() {
    gmd::System system = make_system(kAtomCount);
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, 0.0, gmd::VelocityInitMode::Random, true);
    for (std::size_t i = 0; i < system.atom_count(); ++i) {
        const auto& v = system.velocities()[i];
        check(v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0,
              "atom " + std::to_string(i) +
                  " has a non-zero velocity after initialising to 0 K");
    }
    check(gmd::compute_twice_ke(system) == 0.0,
          "a 0 K system must have exactly zero kinetic energy");

    // Centre-of-mass removal must actually remove it, and must change the
    // degrees of freedom the temperature is divided by.
    gmd::System moving = make_system(kAtomCount);
    gmd::VelocityInitializer other(4242u);
    other.initialize(moving, kTargetTemperature, gmd::VelocityInitMode::Random, true);
    const auto masses = moving.masses();
    const auto velocities = moving.velocities();
    double total_mass = 0.0;
    std::array<double, 3> momentum = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < moving.atom_count(); ++i) {
        total_mass += masses[i];
        for (std::size_t d = 0; d < 3; ++d) momentum[d] += masses[i] * velocities[i][d];
    }
    for (std::size_t d = 0; d < 3; ++d) {
        check(std::fabs(momentum[d] / total_mass) < 1.0e-12,
              "centre-of-mass velocity component " + std::to_string(d) +
                  " is " + number(momentum[d] / total_mass) + ", not removed");
    }

    gmd::DegreesOfFreedomConfig with_com;
    with_com.remove_center_of_mass_velocity = true;
    gmd::DegreesOfFreedomConfig without_com;
    without_com.remove_center_of_mass_velocity = false;
    check(gmd::compute_degrees_of_freedom(kAtomCount, with_com) ==
              reference_dof(kAtomCount, true),
          "compute_degrees_of_freedom disagrees with 3N-3 for a COM-removed system");
    check(gmd::compute_degrees_of_freedom(kAtomCount, without_com) ==
              reference_dof(kAtomCount, false),
          "compute_degrees_of_freedom disagrees with 3N for a free system");
}

void test_initializer_hits_its_target_scale() {
    // The rescale must be exact for every target and both COM settings: the
    // recovered constant may not drift with either.
    for (const double target : {1.0, 77.0, 300.0, 1500.0}) {
        for (const bool com : {true, false}) {
            const double measured = measure_initializer_constant(target, com);
            check(close(measured, kExpectedInitializer),
                  "initialising to " + number(target, 6) + " K with COM removal " +
                      (com ? "on" : "off") + " measures k_B = " + number(measured) +
                      ", so the rescale is not exact");
        }
    }
}

void test_nose_hoover_mass_scales_correctly() {
    // Q = dof * k_B * T * tau^2. Vary each factor independently; the recovered
    // constant must not move.
    for (const double target : {50.0, 300.0, 900.0}) {
        const double measured = measure_nose_hoover_constant(target, 40.0, kAtomCount);
        check(close(measured, kExpectedShared),
              "Nose-Hoover Q at T=" + number(target, 6) + " K measures k_B = " +
                  number(measured) + "; Q is not proportional to T");
    }
    for (const double tau : {10.0, 40.0, 250.0}) {
        const double measured = measure_nose_hoover_constant(300.0, tau, kAtomCount);
        check(close(measured, kExpectedShared),
              "Nose-Hoover Q at tau=" + number(tau, 6) + " measures k_B = " +
                  number(measured) + "; Q is not proportional to tau^2");
    }
    for (const std::size_t atoms : {std::size_t{4}, std::size_t{8}, std::size_t{16}}) {
        const double measured = measure_nose_hoover_constant(300.0, 40.0, atoms);
        check(close(measured, kExpectedShared),
              "Nose-Hoover Q with " + std::to_string(atoms) +
                  " atoms measures k_B = " + number(measured) +
                  "; Q is not proportional to the degrees of freedom");
    }
}

void test_no_missing_inverted_or_doubled_factor() {
    const double paths[] = {measure_initializer_constant(),
                            measure_reporting_constant(),
                            measure_nose_hoover_constant(),
                            measure_velocity_rescaling_constant()};
    const char* labels[] = {"velocity initialization", "temperature reporting",
                            "Nose-Hoover thermostat mass",
                            "velocity-rescaling thermostat"};
    for (std::size_t i = 0; i < 4; ++i) {
        const double v = paths[i];
        check(std::fabs(v - 1.0) > 0.5,
              std::string(labels[i]) +
                  " measures 1, i.e. k_B is missing from that path entirely");
        check(std::fabs(v / (kExpectedShared * kExpectedShared) - 1.0) > 1.0e-3,
              std::string(labels[i]) + " measures k_B squared (" + number(v) +
                  "), i.e. the constant is applied twice");
        check(std::fabs(v * kExpectedShared - 1.0) > 1.0e-3,
              std::string(labels[i]) + " measures 1/k_B (" + number(v) +
                  "), i.e. the temperature conversion is inverted");
    }
}

// --- restart -----------------------------------------------------------------

void test_checkpoints_carry_the_constant_indirectly() {
    // The constant is not serialized anywhere by name. But the Nose-Hoover
    // thermostat mass IS serialized, and Q = dof * k_B * T * tau^2, so a
    // checkpoint written under one constant carries that constant's Q into any
    // run that resumes it. load_checkpoint_state() restores Q verbatim -- by
    // design, so the continued trajectory is deterministic -- which means a
    // resumed run keeps the OLD thermostat mass while every other path uses the
    // new constant. That is a real restart consequence, and it is asserted here
    // rather than assumed.
    constexpr double tau = 40.0;
    gmd::System system = make_system(kAtomCount);
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, true);

    gmd::NoseHooverThermostat writer(tau);
    writer.initialize(system);
    writer.apply_half_kick(system, 0.5, kTargetTemperature);
    const std::string state = writer.checkpoint_state();
    const double written_q = parse_checkpoint_field(state, "Q");

    check(state.find("Q ") != std::string::npos,
          "the Nose-Hoover checkpoint no longer records Q; this test's premise "
          "about restart behaviour needs revisiting");
    check(state.find("8.617") == std::string::npos,
          "a Boltzmann constant literal appears in the checkpoint state string; "
          "the constant is supposed to enter only through Q");

    gmd::NoseHooverThermostat reader(tau);
    reader.initialize(system);
    // load_checkpoint_state() refuses a provisional degrees-of-freedom count;
    // the integrator installs the authoritative one, so do the same here.
    reader.set_degrees_of_freedom(writer.degrees_of_freedom());
    reader.load_checkpoint_state(state);
    const double restored_q = parse_checkpoint_field(reader.checkpoint_state(), "Q");
    check(restored_q == written_q,
          "Q was not restored verbatim: wrote " + number(written_q) +
              ", read back " + number(restored_q));

    // And a checkpoint written under a different constant stays different: the
    // reader does not recompute Q from the current constant.
    const double foreign_q = written_q * 1.5;
    std::ostringstream foreign;
    foreign.precision(17);
    foreign << "tau " << tau << " xi " << parse_checkpoint_field(state, "xi")
            << " Q " << foreign_q
            << " dof " << writer.degrees_of_freedom()
            << " current_temperature " << kTargetTemperature;
    gmd::NoseHooverThermostat resumed(tau);
    resumed.initialize(system);
    resumed.set_degrees_of_freedom(writer.degrees_of_freedom());
    resumed.load_checkpoint_state(foreign.str());
    check(parse_checkpoint_field(resumed.checkpoint_state(), "Q") == foreign_q,
          "a checkpoint's Q was overwritten on load; resuming would silently "
          "change the extended-system dynamics");
}

// --- MC barostat ------------------------------------------------------------
//
// The barostat's constant cannot be recovered the way the others can: its
// Metropolis exponent is -w/kT with w = dU + P_ext dV - N kT ln(V'/V), and the
// accept/reject outcome also depends on an unobservable uniform draw, so one
// measurement cannot separate k_B from that draw.
//
// What IS observable, exactly and deterministically, is the structure of the
// exponent. With dU held at zero the exponent is
//
//     -P_ext dV / (k_B T) + N ln(V'/V)
//
// so it depends on P_ext and T only through the ratio P_ext / T. Running the
// same seeded barostat at (T, P) and at (cT, cP) must therefore produce a
// bit-identical trajectory of accepted volumes, while (T, P) and (T, 2P) must
// not. That is the sign-and-scaling check.
//
// It is not, on its own, sensitive to a small change in the constant: the
// outcome is a discrete accept or reject, and a 1e-11 shift in the exponent
// flips nothing. So there is a second test. For a fixed seed the first trial
// move draws a fixed dV and a fixed uniform u, and the decision flips at the
// temperature where
//
//     ln u = -P_ext dV / (k_B T*) + N ln(V'/V)
//
// Everything on the right except k_B T* is fixed by the seed and the geometry,
// so T* is EXACTLY inversely proportional to the barostat's Boltzmann constant
// and to nothing else. Bisecting for T* therefore measures that constant, to
// whatever precision the bisection reaches. Restoring the previous
// 8.617333262e-5 moves T* from 553693.28283268539 to 553693.28284201352, a
// relative shift of 1.684e-11 -- exactly the constants' own difference, which
// is the proof that T* tracks it. The tolerance below catches that with
// seventeen times' margin, and the value is identical between -O0 and -O2
// builds.

class ConstantEnergyProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "constant_energy"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request, gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        // Volume-independent so dU is exactly zero for every trial move, which
        // is what leaves P_ext/T as the only temperature-dependent handle.
        result.success = true;
        result.potential_energy = 0.0;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        result.virial_valid = true;
    }
};

std::vector<double> run_barostat(double temperature, double pressure_bar,
                                 int attempts = 60) {
    gmd::System system = make_system(kAtomCount);
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, temperature, gmd::VelocityInitMode::Random, true);
    system.set_potential_energy(0.0);

    ConstantEnergyProvider provider;
    gmd::RuntimeContext runtime;
    gmd::MCBarostat barostat(/*frequency=*/1, /*max_delta_ln_V=*/0.02,
                             /*adjust_interval=*/1000000, /*seed=*/777u);
    std::vector<double> volumes;
    for (int step = 0; step < attempts; ++step) {
        barostat.apply(system, provider, runtime, static_cast<std::uint64_t>(step),
                       1.0, temperature, pressure_bar, 0.0);
        const auto& lengths = system.box().lengths;
        volumes.push_back(lengths[0] * lengths[1] * lengths[2]);
    }
    return volumes;
}

// Bisects for the temperature at which the first trial move's accept/reject
// decision flips. See the block comment above: T* is exactly inversely
// proportional to the barostat's Boltzmann constant.
bool first_move_changes_volume(double temperature, double pressure_bar) {
    gmd::System system = make_system(kAtomCount);
    gmd::VelocityInitializer initializer(20260830u);
    initializer.initialize(system, 300.0, gmd::VelocityInitMode::Random, true);
    system.set_potential_energy(0.0);
    const auto& before = system.box().lengths;
    const double volume_before = before[0] * before[1] * before[2];

    ConstantEnergyProvider provider;
    gmd::RuntimeContext runtime;
    gmd::MCBarostat barostat(1, 0.02, 1000000, 777u);
    barostat.apply(system, provider, runtime, 0, 1.0, temperature, pressure_bar, 0.0);
    const auto& after = system.box().lengths;
    return after[0] * after[1] * after[2] != volume_before;
}

double bisect_critical_temperature(double pressure_bar) {
    double low = 1.0, high = 1.0e7;
    const bool high_state = first_move_changes_volume(high, pressure_bar);
    if (first_move_changes_volume(low, pressure_bar) == high_state) return 0.0;
    for (int i = 0; i < 200; ++i) {
        const double middle = 0.5 * (low + high);
        if (first_move_changes_volume(middle, pressure_bar) == high_state) {
            high = middle;
        } else {
            low = middle;
        }
    }
    return 0.5 * (low + high);
}

void test_mc_barostat_constant_matches() {
    // A regression pin, and deliberately so: T* cannot be predicted without the
    // uniform draw the barostat never exposes. What makes it a measurement of
    // the constant rather than an arbitrary number is that T* is exactly
    // inversely proportional to it, so any relative change in the barostat's
    // k_B appears as the same relative change here.
    constexpr double kCriticalTemperature = 553693.28283268539;
    constexpr double kTolerance = 1.0e-12;   // 17x below the 1.684e-11 signal

    const double measured = bisect_critical_temperature(1.0e6);
    std::cout << "\n  MC barostat accept/reject flips at T* = "
              << number(measured) << " K\n";
    check(measured > 0.0,
          "no accept/reject flip was found in [1, 1e7] K, so this fixture no "
          "longer measures the barostat's constant at all");
    check(std::fabs(measured / kCriticalTemperature - 1.0) <= kTolerance,
          "the MC barostat's accept/reject boundary is at T* = " +
              number(measured) + " K, not " + number(kCriticalTemperature) +
              " K (relative " + number(measured / kCriticalTemperature - 1.0, 4) +
              "). T* is exactly inversely proportional to the barostat's "
              "Boltzmann constant, so this says the barostat is no longer using "
              "the same constant as everything else. Restoring the previous "
              "8.617333262e-5 puts T* at 553693.28284201352.");
}

void test_mc_barostat_beta_structure() {
    const auto base = run_barostat(300.0, 500.0);

    // Same P/T ratio: the exponent is identical, so the accepted-volume
    // trajectory must be bit-identical.
    for (const double c : {0.5, 2.0, 7.0}) {
        const auto scaled = run_barostat(300.0 * c, 500.0 * c);
        bool identical = scaled.size() == base.size();
        for (std::size_t i = 0; identical && i < base.size(); ++i) {
            identical = (scaled[i] == base[i]);
        }
        check(identical,
              "scaling temperature and pressure together by " + number(c, 3) +
                  " changed the accepted-volume trajectory; the Metropolis "
                  "exponent does not depend on P_ext and T only through "
                  "P_ext / (k_B T)");
    }

    // Different P/T ratio: it must NOT be identical, or the check above would
    // be satisfied by an exponent that ignores pressure and temperature alike.
    const auto doubled = run_barostat(300.0, 1000.0);
    bool differs = doubled.size() != base.size();
    for (std::size_t i = 0; !differs && i < base.size(); ++i) {
        differs = (doubled[i] != base[i]);
    }
    check(differs,
          "doubling the target pressure at fixed temperature left the "
          "accepted-volume trajectory unchanged, so the acceptance test is not "
          "using P_ext / (k_B T) at all");

    // Sign: at a high enough pressure the box must end up compressed relative
    // to a low-pressure run, since the P dV work term penalises expansion.
    const auto low = run_barostat(300.0, 1.0);
    const auto high = run_barostat(300.0, 200000.0);
    check(high.back() < low.back(),
          "raising the target pressure from 1 bar to 200000 bar did not compress "
          "the box (final volumes " + number(high.back(), 8) + " vs " +
              number(low.back(), 8) + "); the sign of the P dV term is wrong");
}

}  // namespace

int main() {
    std::cout << "Boltzmann constant audit\n";
    test_reference_matches_the_authoritative_value();
    test_production_matches_the_authoritative_value();
    test_known_velocities_give_exact_temperature();
    test_paths_are_mutually_consistent();
    test_initialize_then_report_round_trip();
    test_initializer_hits_its_target_scale();
    test_zero_temperature_and_degrees_of_freedom();
    test_nose_hoover_mass_scales_correctly();
    test_no_missing_inverted_or_doubled_factor();
    test_checkpoints_carry_the_constant_indirectly();
    test_mc_barostat_beta_structure();
    test_mc_barostat_constant_matches();

    if (failures != 0) {
        std::cerr << "\nBoltzmann constant audit failed: " << failures << '\n';
        return 1;
    }
    std::cout << "\nBoltzmann constant audit passed\n";
    return 0;
}
