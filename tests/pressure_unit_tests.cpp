// Audit of the conversion between GMD's internal pressure unit and bar.
//
// GMD computes pressure in its own unit system as P = (2K + tr W) / 3V, which
// with energies in eV and lengths in Angstrom makes every internal pressure an
// energy density in eV/A^3. Users never see that number. They set a target
// pressure in bar and they read a P[bar] column out of the log, so exactly one
// conversion stands between the two, and it is used in both directions:
//
//   trajectory writer   P[bar]  = P_internal / (bar -> eV/A^3)
//   Monte Carlo barostat P_ext  = P_target[bar] * (bar -> eV/A^3)
//   Berendsen barostat   compares a target against an instantaneous pressure
//
// Nothing in the codebase forces those sites to agree. Each holds its own
// literal, so this file measures the conversion back out of each path rather
// than reading any production symbol, and requires the measurements to match
// each other.
//
// HOW EACH PATH IS MEASURED
//
//   reporting     a System is given a completed-step pressure of known internal
//                 magnitude, one frame is written, and the bar column is read
//                 back: conversion = P_internal / P_bar.
//
//   MC barostat   the conversion sits inside a Metropolis exponent, so it is
//                 recovered from the target pressure at which a fixed-seed trial
//                 move flips from accepted to rejected. See the derivation above
//                 measure_mc_barostat_conversion(): a second run at a different
//                 atom count cancels the random draw exactly, leaving a closed
//                 form for the conversion with no unknown in it.
//
//   Berendsen     the coupling factor mu is observable from the box, and
//                 mu^3 = 1 - beta*(dt/tau)*(P_target - P_current) inverts to the
//                 instantaneous pressure the barostat believed it had. Comparing
//                 that against the true internal pressure says which unit the
//                 comparison was performed in.
//
// THE AUTHORITATIVE VALUE
//
// This conversion is a pure unit identity with no measured quantity in it. All
// four ingredients are exact by definition:
//
//   1 bar = 100000 Pa            (definition of the bar)
//   1 Pa  = 1 J/m^3              (definition of the pascal)
//   1 A   = 1e-10 m              (definition of the angstrom)
//   1 eV  = 1.602176634e-19 J    (exact since the 2019 SI redefinition,
//                                 https://physics.nist.gov/cgi-bin/cuu/Value?e)
//
// so
//
//   1 bar = 1e5 J/m^3 = 1e5 * 1e-30 J/A^3 = 1e-25 J/A^3
//         = 1e-25 / 1.602176634e-19 eV/A^3
//         = 1e-6 / 1.602176634 eV/A^3
//         = 1000 / 1602176634  =  500 / 801088317 eV/A^3
//         = 6.24150907446076260777624098...e-7 eV/A^3
//
// The reverse direction is the reciprocal, 801088317 / 500, which unlike the
// forward direction TERMINATES:
//
//   1 eV/A^3 = 1602176.634 bar exactly.
//
// The reference below is written as those two integers so that the derivation is
// visible and checkable rather than copied, and cross-checked against the naive
// floating-point route through the SI literals.

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/berendsen_barostat.hpp"
#include "gmd/integrator/mc_barostat.hpp"
#include "gmd/io/checkpoint.hpp"
#include "gmd/io/trajectory_writer.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[pressure units] " << message << '\n';
        ++failures;
    }
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

double relative_difference(double measured, double expected) {
    if (expected == 0.0) return std::abs(measured);
    return std::abs(measured - expected) / std::abs(expected);
}

// --- the independent reference -------------------------------------------
//
// 500 / 801088317 and its reciprocal. Both numerators and denominators are
// integers well inside the exactly-representable range, so each quotient is a
// single correctly-rounded division of two exact doubles: the results are the
// nearest doubles to the exact rationals, which is the best any literal could
// do. Deriving them here rather than importing gmd::kBarToEVPerAngstromCubed is
// deliberate -- a test that imported the production constant would agree with a
// wrong production constant.
constexpr double kReferenceBarToEVPerA3 = 500.0 / 801088317.0;
constexpr double kReferenceEVPerA3ToBar = 801088317.0 / 500.0;

// The value this correction supersedes: the forward conversion rounded to eight
// significant figures, 4.091837e-09 relative high.
constexpr double kSupersededValue = 6.2415091e-7;

// The same numbers by the naive route through the SI literals. This is not the
// definition -- it rounds three times rather than once, and lands one ulp high
// on the forward direction -- but it is an independent path to the same value
// and it is what catches a transcription error in the integers above.
constexpr double kElementaryChargeCoulombs = 1.602176634e-19;  // exact, SI
constexpr double kBoltzmannJoulesPerKelvin = 1.380649e-23;     // exact, SI
// Derived here rather than imported for the same reason as the pressure
// reference: tests/boltzmann_constant_tests.cpp audits the production k_B, and
// this file must not depend on that audit having passed.
constexpr double kReferenceBoltzmannEVPerKelvin =
    kBoltzmannJoulesPerKelvin / kElementaryChargeCoulombs;
constexpr double kBarInPascals             = 1.0e5;            // exact
constexpr double kCubicAngstromInCubicM    = 1.0e-30;          // exact

void test_reference_derivation_is_self_consistent() {
    const double naive =
        kBarInPascals * kCubicAngstromInCubicM / kElementaryChargeCoulombs;
    // One ulp at this magnitude is ~1.6e-16 relative; allow two.
    check(relative_difference(naive, kReferenceBarToEVPerA3) < 4.0e-16,
          "the rational reference " + number(kReferenceBarToEVPerA3) +
              " disagrees with the SI-literal route " + number(naive));

    // The reciprocal relation must hold to the last bit, or the two directions
    // are not the same conversion.
    check(kReferenceBarToEVPerA3 * kReferenceEVPerA3ToBar == 1.0,
          "the two reference directions are not exact reciprocals: product " +
              number(kReferenceBarToEVPerA3 * kReferenceEVPerA3ToBar));
    check(1.0 / kReferenceEVPerA3ToBar == kReferenceBarToEVPerA3,
          "1 / (eV/A^3 -> bar) does not reproduce (bar -> eV/A^3)");

    // The terminating decimal, stated as such.
    check(kReferenceEVPerA3ToBar == 1602176.634,
          "1 eV/A^3 must be exactly 1602176.634 bar; got " +
              number(kReferenceEVPerA3ToBar));
}

// --- path 1: pressure reporting -------------------------------------------

const char* kLogStem = "gmd_pressure_unit_probe";

// Reads the P[bar] column out of the single frame written by write_probe_frame.
double read_pressure_bar(const std::filesystem::path& log_path) {
    std::ifstream input(log_path);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        double step = 0.0, time = 0.0, pe = 0.0, ke = 0.0, total = 0.0;
        double temperature = 0.0, pressure_bar = 0.0;
        if (!(fields >> step >> time >> pe >> ke >> total >> temperature >>
              pressure_bar)) {
            continue;
        }
        return pressure_bar;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// Writes one frame whose completed-step pressure is exactly `internal_pressure`
// and returns the bar value the log reported.
double reported_bar_for_internal_pressure(double internal_pressure) {
    gmd::System system;
    system.resize(2, 2);
    gmd::Box box;
    box.set_lengths({10.0, 10.0, 10.0});
    system.set_box(box);
    system.mutable_masses()[0] = 12.0;
    system.mutable_masses()[1] = 12.0;
    system.mutable_coordinates()[0] = {1.0, 1.0, 1.0};
    system.mutable_coordinates()[1] = {4.0, 2.0, 3.0};

    // The completed-step record is the writer's preferred source and carries the
    // pressure directly, so the reported number is a pure conversion of this one
    // field and nothing else.
    gmd::System::StepThermodynamics completed;
    completed.valid = true;
    completed.pressure = internal_pressure;
    completed.volume = 1000.0;
    completed.twice_kinetic_energy = 0.0;
    completed.potential_energy = 0.0;
    system.set_step_thermodynamics(completed);

    const std::filesystem::path stem =
        std::filesystem::temp_directory_path() / kLogStem;
    gmd::TrajectoryWriter writer;
    writer.open(stem);
    writer.write_frame(system, 0, 0.0, 0.0, 3);
    writer.close();

    const std::filesystem::path log_path = stem.string() + ".log";
    const double bar = read_pressure_bar(log_path);
    std::filesystem::remove(log_path);
    std::filesystem::remove(stem.string() + ".xyz");
    return bar;
}

// conversion = P_internal / P_bar.
double measure_reporting_conversion(double internal_pressure) {
    const double bar = reported_bar_for_internal_pressure(internal_pressure);
    return internal_pressure / bar;
}

void test_reporting_carries_exactly_one_conversion() {
    // A large internal pressure so that the log's six fixed decimals leave far
    // more significant digits than double precision needs: at 1e6 eV/A^3 the
    // printed bar value is ~1.6e12 and one printed unit is 6e-19 relative.
    const double measured = measure_reporting_conversion(1.0e6);
    check(std::isfinite(measured),
          "the reporting path produced no finite conversion");

    // Structural, not a value pin: a path that omitted the conversion would
    // measure 1, one that inverted it would measure 1/k ~ 1.6e6, one that
    // applied it twice would measure k^2 ~ 3.9e-13. All three are orders of
    // magnitude away from a conversion near 6e-7.
    check(measured > 1.0e-8 && measured < 1.0e-5,
          "the reported bar value does not correspond to a single bar-to-eV/A^3 "
          "conversion; measured " + number(measured));
    std::cout << "  reporting conversion      " << number(measured) << '\n';
}

void test_reporting_sign_and_special_values() {
    // Zero is a pressure, not a sentinel: it must convert to exactly zero and
    // still be reported as valid.
    const double zero_bar = reported_bar_for_internal_pressure(0.0);
    check(zero_bar == 0.0,
          "zero internal pressure must report exactly 0 bar, got " +
              number(zero_bar));

    // A negative pressure -- a system under tension -- must stay negative. A
    // conversion that took a magnitude, or that was subtracted rather than
    // divided, would show up here and nowhere else.
    const double positive = reported_bar_for_internal_pressure(1.0e6);
    const double negative = reported_bar_for_internal_pressure(-1.0e6);
    check(negative < 0.0,
          "negative internal pressure must report negative bar, got " +
              number(negative));
    check(relative_difference(negative, -positive) < 1.0e-15,
          "the conversion is not odd: +P reports " + number(positive) +
              " but -P reports " + number(negative));

    // Same conversion across eleven orders of magnitude, which is what rules out
    // an additive offset masquerading as a scale factor.
    const double small = measure_reporting_conversion(1.0e-3);
    const double large = measure_reporting_conversion(1.0e8);
    check(relative_difference(small, large) < 1.0e-9,
          "the reporting conversion is not a pure scale factor: " +
              number(small) + " at small pressure vs " + number(large) +
              " at large pressure");
}

void test_round_trip_through_the_reference() {
    // bar -> eV/A^3 -> bar for a spread of magnitudes and both signs.
    const double values[] = {0.0, 1.0, -1.0, 1.013e5, -2.5e9, 1.0e-12, 1.0e15};
    for (double bar : values) {
        const double internal = bar * kReferenceBarToEVPerA3;
        const double back = internal * kReferenceEVPerA3ToBar;
        check(relative_difference(back, bar) < 1.0e-15,
              "round trip failed for " + number(bar) + " bar: came back as " +
                  number(back));
    }
}

// --- path 2: the Monte Carlo barostat -------------------------------------

// A provider whose energy never changes, so that dU = 0 in the Metropolis
// exponent and the acceptance depends only on the pressure-work term and the
// phase-space term.
class ConstantEnergyProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "constant_energy"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        const std::size_t n = request.coordinates.size();
        result.forces.assign(n, {0.0, 0.0, 0.0});
        result.potential_energy = 0.0;
        result.virial_valid = false;
        result.success = true;
    }
};

gmd::System make_barostat_system(std::size_t atom_count, double box_length) {
    gmd::System system;
    system.resize(atom_count, atom_count);
    gmd::Box box;
    box.set_lengths({box_length, box_length, box_length});
    system.set_box(box);
    for (std::size_t i = 0; i < atom_count; ++i) {
        system.mutable_masses()[i] = 12.0;
        const double t = static_cast<double>(i);
        system.mutable_coordinates()[i] = {1.0 + 0.37 * t, 2.0 + 0.11 * t,
                                           3.0 + 0.23 * t};
    }
    system.set_potential_energy(0.0);
    return system;
}

constexpr std::uint32_t kBarostatSeed = 12345;
constexpr double kBarostatMaxDeltaLnV = 0.01;
constexpr double kBarostatTemperature = 300.0;

// Runs one fixed-seed trial move and reports whether the cell changed, i.e.
// whether the move was accepted.
bool first_move_is_accepted(std::size_t atom_count, double box_length,
                            double target_pressure_bar) {
    gmd::System system = make_barostat_system(atom_count, box_length);
    ConstantEnergyProvider provider;
    gmd::RuntimeContext runtime;
    gmd::MCBarostat barostat(1, kBarostatMaxDeltaLnV, 1000000, kBarostatSeed);
    const double before = system.box().lengths[0];
    barostat.apply(system, provider, runtime, 0, 1.0, kBarostatTemperature,
                   target_pressure_bar, 0.0);
    return system.box().lengths[0] != before;
}

// The volume ratio the fixed seed proposes, read off an accepted move. It is a
// property of the seed alone, not of the pressure.
double proposed_ln_volume_ratio(std::size_t atom_count, double box_length) {
    gmd::System system = make_barostat_system(atom_count, box_length);
    ConstantEnergyProvider provider;
    gmd::RuntimeContext runtime;
    gmd::MCBarostat barostat(1, kBarostatMaxDeltaLnV, 1000000, kBarostatSeed);
    const double before = system.box().lengths[0];
    // A hugely negative target makes the pressure-work term favour the move, so
    // w <= 0 and the move is accepted deterministically.
    barostat.apply(system, provider, runtime, 0, 1.0, kBarostatTemperature,
                   -1.0e12, 0.0);
    const double after = system.box().lengths[0];
    if (after == before) {
        // Rejected. With dU = 0 and a target that large and negative, the only
        // way w can be positive is for the pressure-work term to enter with the
        // wrong sign -- which is a real finding, and one worth reporting here
        // rather than letting it become a division by a zero volume change and
        // a downstream NaN.
        check(false,
              "the MC barostat rejected a volume move at a target pressure of "
              "-1e12 bar. With no energy change, the pressure-work term is the "
              "only thing that can oppose it, so it is entering the Metropolis "
              "weight with the wrong sign");
        return 0.0;
    }
    return 3.0 * std::log(after / before);
}

// Bisects for the target pressure at which the fixed-seed first move flips from
// accepted to rejected.
double critical_target_pressure(std::size_t atom_count, double box_length,
                                double low, double high) {
    // Precondition: accepted at `low`, rejected at `high`.
    for (int iteration = 0; iteration < 200; ++iteration) {
        const double mid = 0.5 * (low + high);
        if (mid == low || mid == high) break;
        if (first_move_is_accepted(atom_count, box_length, mid)) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return 0.5 * (low + high);
}

// Recovers the bar -> eV/A^3 conversion from the barostat's own arithmetic.
//
// With dU = 0 the Metropolis weight is
//
//     w = P_target * c * dV - N * kT * ln(V'/V)
//
// where c is the conversion under test. The move is accepted when w <= 0 or
// when a uniform draw u satisfies u < exp(-w/kT); at the flip point of a
// fixed-seed run the second condition is an equality, so
//
//     P* * c * dV - N * kT * L = -kT * ln u        with L = ln(V'/V).
//
// u is unknown, but it depends only on the seed. Running the identical bisection
// at a second atom count N2 -- same box, same seed, same proposed move, same
// draw, and dU still zero -- gives a second equation with the same u, and
// subtracting eliminates it entirely:
//
//     (P1* - P2*) * c * dV = kT * L * (N1 - N2)
//
//     c = kT * L * (N1 - N2) / ((P1* - P2*) * dV)
//
// Nothing random survives. kT here is the barostat's own k_B * T, which is
// audited separately in tests/boltzmann_constant_tests.cpp; taking it from
// gmd::kBoltzmannConstantEVPerKelvin would make this measurement depend on that
// constant being right, so it is instead cancelled the same way u is -- see
// measure_mc_barostat_conversion(), which forms a ratio in which kT drops out.
struct BarostatFlip {
    double critical_pressure;
    double volume_change;
    double ln_ratio;
};

BarostatFlip barostat_flip(std::size_t atom_count, double box_length) {
    const double ln_ratio = proposed_ln_volume_ratio(atom_count, box_length);
    const double volume = box_length * box_length * box_length;
    const double volume_change = volume * (std::exp(ln_ratio) - 1.0);

    // Bracket: at zero target pressure the phase-space term alone favours an
    // expansion, so the move is accepted; a large enough target rejects it.
    double high = 1.0;
    while (first_move_is_accepted(atom_count, box_length, high) && high < 1.0e30) {
        high *= 4.0;
    }
    const double critical =
        critical_target_pressure(atom_count, box_length, 0.0, high);
    return BarostatFlip{critical, volume_change, ln_ratio};
}

// Returns c / (k_B), i.e. the conversion divided by the Boltzmann constant, so
// that no assumption about k_B enters. The caller compares two such ratios, or
// multiplies by an independently derived k_B.
double measure_mc_barostat_conversion_over_kT() {
    constexpr std::size_t kFirstAtomCount = 4;
    constexpr std::size_t kSecondAtomCount = 12;
    constexpr double kBoxLength = 20.0;

    const BarostatFlip first = barostat_flip(kFirstAtomCount, kBoxLength);
    const BarostatFlip second = barostat_flip(kSecondAtomCount, kBoxLength);

    // Same seed, same box: the proposed move must be identical, or the
    // subtraction below is not cancelling like with like.
    check(first.ln_ratio == second.ln_ratio,
          "the two barostat runs did not propose the same volume move: " +
              number(first.ln_ratio) + " vs " + number(second.ln_ratio));

    const double atom_difference =
        static_cast<double>(kFirstAtomCount) - static_cast<double>(kSecondAtomCount);
    const double pressure_difference =
        first.critical_pressure - second.critical_pressure;
    // c / (k_B * T) = L * (N1 - N2) / ((P1* - P2*) * dV)
    return first.ln_ratio * atom_difference /
           (pressure_difference * first.volume_change);
}

void test_mc_barostat_pressure_work_uses_one_conversion() {
    const double ratio = measure_mc_barostat_conversion_over_kT();
    check(std::isfinite(ratio) && ratio > 0.0,
          "the MC barostat conversion measurement did not converge: " +
              number(ratio));

    const double measured =
        ratio * kReferenceBoltzmannEVPerKelvin * kBarostatTemperature;

    std::cout << "  MC barostat conversion    " << number(measured) << '\n';

    // Structural bound, as for the reporting path: near 6e-7 rather than 1,
    // 1.6e6 or 4e-13.
    check(measured > 1.0e-8 && measured < 1.0e-5,
          "the MC barostat's pressure-work term does not carry a single "
          "bar-to-eV/A^3 conversion; measured " + number(measured));
}

void test_mc_barostat_pressure_work_has_the_right_sign() {
    // Raising the external pressure must make expansion harder, never easier.
    // If P_ext*dV entered with the wrong sign, the ordering below would invert.
    constexpr double kBoxLength = 20.0;
    constexpr std::size_t kAtoms = 4;
    const BarostatFlip flip = barostat_flip(kAtoms, kBoxLength);
    check(flip.volume_change > 0.0,
          "this fixture assumes the seeded move is an expansion; it is not");
    check(flip.critical_pressure > 0.0,
          "there is no positive target pressure at which the expansion is "
          "rejected, so the pressure-work term is not opposing expansion");
    check(first_move_is_accepted(kAtoms, kBoxLength,
                                 flip.critical_pressure * 0.5),
          "halving the target pressure should keep an expansion accepted");
    check(!first_move_is_accepted(kAtoms, kBoxLength,
                                  flip.critical_pressure * 2.0),
          "doubling the target pressure should reject the expansion");
}

// --- path 3: the Berendsen barostat ---------------------------------------
//
// The Berendsen barostat holds no conversion constant of its own. It compares a
// target against an instantaneous pressure, and the only question is which unit
// that comparison happens in. The coupling factor
//
//     mu^3 = 1 - beta * (dt / tau) * (P_target - P_current)
//
// is fully observable -- mu is the ratio of the box lengths before and after --
// so the relation inverts to the instantaneous pressure the barostat believed it
// had:
//
//     P_believed = P_target - (1 - mu^3) * tau / (beta * dt)
//
// Comparing P_believed against the true internal pressure, expressed in bar,
// says whether the two sides of the subtraction were in the same unit. They were
// not: P_believed came back equal to the raw eV/A^3 number, so a run asking for
// 1 bar was in fact asking for 1 eV/A^3, which is 1602176.634 bar.

class NullForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "null"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.forces.assign(request.coordinates.size(), {0.0, 0.0, 0.0});
        result.potential_energy = 0.0;
        result.virial_valid = false;
        result.success = true;
    }
};

struct BerendsenProbe {
    double internal_pressure;   // [eV/A^3], computed here from the same inputs
    double believed_pressure;   // [bar or eV/A^3 -- that is what is under test]
};

BerendsenProbe probe_berendsen(double target_pressure_bar) {
    constexpr double kBoxLength = 12.0;
    constexpr double kTimeStep = 1.0;
    constexpr double kTau = 10.0;
    constexpr double kBeta = 4.5e-5;
    constexpr std::size_t kAtoms = 4;

    gmd::System system;
    system.resize(kAtoms, kAtoms);
    gmd::Box box;
    box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    system.set_box(box);
    for (std::size_t i = 0; i < kAtoms; ++i) {
        const double t = static_cast<double>(i);
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {1.0 + 2.0 * t, 2.0 + 0.5 * t,
                                           3.0 + 0.25 * t};
        system.mutable_velocities()[i] = {0.004, -0.002, 0.003};
    }
    // A virial with an unequal, non-zero trace, so that the pressure is not just
    // the kinetic term and a dropped virial would show.
    std::array<double, 9> virial{};
    virial[0] = 0.5;
    virial[4] = 0.3;
    virial[8] = 0.2;
    system.set_last_virial(virial, true);

    double twice_ke = 0.0;
    for (std::size_t i = 0; i < kAtoms; ++i) {
        const auto velocity = system.velocities()[i];
        twice_ke += system.masses()[i] *
                    (velocity[0] * velocity[0] + velocity[1] * velocity[1] +
                     velocity[2] * velocity[2]);
    }
    const double volume = kBoxLength * kBoxLength * kBoxLength;
    const double trace = virial[0] + virial[4] + virial[8];
    const double internal = (twice_ke + trace) / (3.0 * volume);

    NullForceProvider provider;
    gmd::RuntimeContext runtime;
    gmd::BerendsenBarostat barostat(kTau, kBeta);
    barostat.apply(system, provider, runtime, 0, kTimeStep, 300.0,
                   target_pressure_bar, trace);

    const double mu = system.box().lengths[0] / kBoxLength;
    const double mu_cubed = mu * mu * mu;
    const double believed =
        target_pressure_bar - (1.0 - mu_cubed) * kTau / (kBeta * kTimeStep);
    return BerendsenProbe{internal, believed};
}

void test_berendsen_compares_in_one_unit() {
    const BerendsenProbe probe = probe_berendsen(1.0);
    const double internal_in_bar = probe.internal_pressure * kReferenceEVPerA3ToBar;

    std::cout << "  Berendsen P_current       " << number(probe.believed_pressure)
              << " bar (true " << number(internal_in_bar) << " bar)\n";

    // The recovery goes through a cube root and a subtraction of two nearly
    // equal numbers, so it is a finite-precision inversion rather than an exact
    // one; 1e-6 is far tighter than the 1.6e6 factor it has to be able to see.
    check(relative_difference(probe.believed_pressure, internal_in_bar) < 1.0e-6,
          "the Berendsen barostat compared a target in bar against an "
          "instantaneous pressure of " + number(probe.believed_pressure) +
              ", but the true pressure in bar is " + number(internal_in_bar));

    // Named explicitly: the failure mode this replaces was the raw internal
    // number reaching the subtraction unconverted.
    check(relative_difference(probe.believed_pressure,
                              probe.internal_pressure) > 1.0e3,
          "the Berendsen barostat is still comparing the target against a raw "
          "eV/A^3 pressure (" + number(probe.internal_pressure) + ")");
}

void test_berendsen_and_reporting_agree() {
    // Both barostats and the log must mean the same thing by "1 bar". This is
    // the check that ties the Berendsen path, which holds no constant of its
    // own, to the one the other two share.
    const BerendsenProbe probe = probe_berendsen(1.0);
    const double reporting = measure_reporting_conversion(1.0e6);
    const double implied = probe.internal_pressure / probe.believed_pressure;
    check(relative_difference(implied, reporting) < 1.0e-6,
          "the Berendsen barostat's pressure unit implies a conversion of " +
              number(implied) + ", but the log reports with " +
              number(reporting));
}

void test_berendsen_responds_to_the_target_in_bar() {
    // A target far above the instantaneous pressure must compress the cell, and
    // one far below must expand it. Under the unconverted comparison a target of
    // a few thousand bar was numerically enormous next to a pressure of ~1e-4,
    // so every ordinary target compressed and the sign carried no information.
    const BerendsenProbe reference = probe_berendsen(0.0);
    const double true_bar = reference.internal_pressure * kReferenceEVPerA3ToBar;

    constexpr double kBoxLength = 12.0;
    auto box_after = [](double target) {
        // Re-runs the same fixture and returns the resulting box length.
        gmd::System system;
        system.resize(4, 4);
        gmd::Box box;
        box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
        system.set_box(box);
        for (std::size_t i = 0; i < 4; ++i) {
            const double t = static_cast<double>(i);
            system.mutable_masses()[i] = 12.0;
            system.mutable_coordinates()[i] = {1.0 + 2.0 * t, 2.0 + 0.5 * t,
                                               3.0 + 0.25 * t};
            system.mutable_velocities()[i] = {0.004, -0.002, 0.003};
        }
        std::array<double, 9> virial{};
        virial[0] = 0.5;
        virial[4] = 0.3;
        virial[8] = 0.2;
        system.set_last_virial(virial, true);
        NullForceProvider provider;
        gmd::RuntimeContext runtime;
        gmd::BerendsenBarostat barostat(10.0, 4.5e-5);
        barostat.apply(system, provider, runtime, 0, 1.0, 300.0, target,
                       virial[0] + virial[4] + virial[8]);
        return system.box().lengths[0];
    };

    check(box_after(true_bar * 10.0) < kBoxLength,
          "a target ten times the instantaneous pressure must compress the box");
    check(box_after(true_bar * 0.1) > kBoxLength,
          "a target a tenth of the instantaneous pressure must expand the box");
    // And the crossing is at the instantaneous pressure itself, not somewhere
    // 1.6e6 away from it.
    check(std::abs(box_after(true_bar) - kBoxLength) < 1.0e-12 * kBoxLength,
          "a target equal to the instantaneous pressure must leave the box "
          "unchanged; it moved to " + number(box_after(true_bar)));
}

// --- path 4: restart ------------------------------------------------------
//
// The checkpoint stores step_pressure in eV/A^3, GMD's internal unit, and not in
// bar. That is the right choice -- it is the number the integrator and the
// barostats actually hold -- and it has a consequence worth pinning: the
// conversion is NOT serialized, so a checkpoint written before this correction
// resumes with the corrected conversion. The stored pressure is unchanged, and
// the first frame after the restart reports it in bar slightly differently than
// the last frame before it did, by the constant's own 4.09e-09.
//
// What must hold is that the round trip through the file is exact and that the
// units on both sides of it are the same one.

void test_checkpoint_stores_internal_units_and_round_trips() {
    const double internal_pressure = -3.7182956231e-4;   // [eV/A^3]

    gmd::System system;
    system.resize(2, 2);
    gmd::Box box;
    box.set_lengths({10.0, 11.0, 12.0});
    system.set_box(box);
    system.mutable_masses()[0] = 12.0;
    system.mutable_masses()[1] = 15.999;

    gmd::CheckpointData checkpoint;
    checkpoint.system = &system;
    checkpoint.metadata.step = 17;
    checkpoint.metadata.time_fs = 34.0;
    checkpoint.metadata.step_pressure_valid = true;
    checkpoint.metadata.step_pressure = internal_pressure;
    checkpoint.metadata.step_pressure_volume = 10.0 * 11.0 * 12.0;
    checkpoint.metadata.step_pressure_twice_ke = 0.25;
    checkpoint.metadata.step_pressure_potential_energy = -1.5;

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "gmd_pressure_unit_probe.chk";
    gmd::write_checkpoint(path, checkpoint);

    gmd::System restored;
    const gmd::CheckpointMetadata metadata = gmd::read_checkpoint(path, restored);
    std::filesystem::remove(path);

    check(metadata.step_pressure_valid,
          "the restored checkpoint lost its pressure validity flag");

    // Bit-exact, not approximately: the checkpoint writes enough digits to
    // round-trip a double, and a pressure that came back merely close would mean
    // a restarted run does not continue the trajectory it left.
    check(metadata.step_pressure == internal_pressure,
          "the checkpointed pressure did not round-trip exactly: wrote " +
              number(internal_pressure) + ", read " +
              number(metadata.step_pressure));

    // The stored number is the internal one. If a conversion had been applied on
    // the way in, the value would be the bar figure, 1.6e6 times larger.
    check(relative_difference(metadata.step_pressure, internal_pressure) < 1.0e-15,
          "the checkpoint appears to store pressure in bar rather than eV/A^3");

    // And a frame written from the restored state reports the same bar value a
    // frame written before the checkpoint would have, which is the actual
    // continuity requirement.
    gmd::System::StepThermodynamics completed;
    completed.valid = true;
    completed.pressure = metadata.step_pressure;
    completed.volume = metadata.step_pressure_volume;
    completed.twice_kinetic_energy = metadata.step_pressure_twice_ke;
    completed.potential_energy = metadata.step_pressure_potential_energy;
    restored.set_step_thermodynamics(completed);

    const double before = reported_bar_for_internal_pressure(internal_pressure);
    const std::filesystem::path stem =
        std::filesystem::temp_directory_path() / kLogStem;
    gmd::TrajectoryWriter writer;
    writer.open(stem);
    writer.write_frame(restored, 17, 34.0, completed.twice_kinetic_energy, 3);
    writer.close();
    const double after = read_pressure_bar(stem.string() + ".log");
    std::filesystem::remove(stem.string() + ".log");
    std::filesystem::remove(stem.string() + ".xyz");

    check(before == after,
          "the pressure reported after a restart (" + number(after) +
              " bar) differs from the one reported before it (" +
              number(before) + " bar)");
}

// --- the two production paths must agree ----------------------------------

void test_reporting_and_barostat_agree() {
    const double reporting = measure_reporting_conversion(1.0e6);
    const double ratio = measure_mc_barostat_conversion_over_kT();
    const double barostat =
        ratio * kReferenceBoltzmannEVPerKelvin * kBarostatTemperature;

    // These two constants live in different translation units with nothing
    // linking them, so this is the check that would catch one being edited
    // alone. The bisection resolves the barostat's value to roughly 1e-12
    // relative; the bound is loose enough for that and far tighter than any
    // rounding difference between two independently written literals.
    check(relative_difference(reporting, barostat) < 1.0e-9,
          "the reporting path and the MC barostat disagree about the "
          "bar-to-eV/A^3 conversion: reporting " + number(reporting) +
              " vs barostat " + number(barostat));
}

// --- the production value must be the authoritative one --------------------
//
// The structural checks above hold for any self-consistent conversion, correct
// or not: they all passed against the superseded 6.2415091e-7, which is exactly
// what they did in the commit that introduced them. These are the assertions
// that require the right number, and they are what the correction is for.

void test_reporting_matches_the_authoritative_value() {
    const double measured = measure_reporting_conversion(1.0e6);
    // The log carries the value at full double precision here -- at an internal
    // pressure of 1e6 eV/A^3 the printed bar value has nineteen significant
    // digits before its six fixed decimals begin -- so the only slack needed is
    // the division that recovers the constant.
    check(relative_difference(measured, kReferenceBarToEVPerA3) < 1.0e-15,
          "the reported pressure does not use the exact bar conversion: "
          "measured " + number(measured) + ", reference " +
              number(kReferenceBarToEVPerA3));

    // The superseded value, named, so that a revert is reported as a revert
    // rather than as an anonymous tolerance failure.
    check(relative_difference(measured, kSupersededValue) > 1.0e-9,
          "the reported pressure still uses the superseded 6.2415091e-7, which "
          "is 4.091837e-09 relative above the exact value");
}

void test_mc_barostat_matches_the_authoritative_value() {
    const double ratio = measure_mc_barostat_conversion_over_kT();
    const double measured = ratio * kReferenceBoltzmannEVPerKelvin *
                            kBarostatTemperature;

    // Looser than the reporting bound, because this number comes from a
    // bisection on a discrete accept/reject outcome rather than from a printed
    // value. The bisection itself runs to double resolution; the limit is the
    // conditioning of the two-run subtraction, measured at a few parts in 1e15.
    // 1e-12 leaves room for that and is still three orders of magnitude tighter
    // than the 4.09e-09 deviation it has to be able to see.
    check(relative_difference(measured, kReferenceBarToEVPerA3) < 1.0e-12,
          "the MC barostat's pressure-work term does not use the exact bar "
          "conversion: measured " + number(measured) + ", reference " +
              number(kReferenceBarToEVPerA3));

    check(relative_difference(measured, kSupersededValue) > 1.0e-9,
          "the MC barostat still uses the superseded 6.2415091e-7");
}

void report_measured_values() {
    const double measured = measure_reporting_conversion(1.0e6);
    const double relative = (measured - kReferenceBarToEVPerA3) /
                            kReferenceBarToEVPerA3;
    std::cout << "  reference (500/801088317) " << number(kReferenceBarToEVPerA3)
              << '\n'
              << "  production deviation      " << std::scientific
              << std::setprecision(6) << relative << std::defaultfloat << '\n';
}

}  // namespace

int main() {
    std::cout << "[pressure units] auditing the bar <-> eV/A^3 conversion\n";
    test_reference_derivation_is_self_consistent();
    test_reporting_carries_exactly_one_conversion();
    test_reporting_sign_and_special_values();
    test_round_trip_through_the_reference();
    test_mc_barostat_pressure_work_uses_one_conversion();
    test_mc_barostat_pressure_work_has_the_right_sign();
    test_reporting_and_barostat_agree();
    test_checkpoint_stores_internal_units_and_round_trips();
    test_berendsen_compares_in_one_unit();
    test_berendsen_and_reporting_agree();
    test_berendsen_responds_to_the_target_in_bar();
    test_reporting_matches_the_authoritative_value();
    test_mc_barostat_matches_the_authoritative_value();
    report_measured_values();

    if (failures == 0) {
        std::cout << "[pressure units] audit passed\n";
        return 0;
    }
    std::cerr << "[pressure units] " << failures << " check(s) failed\n";
    return 1;
}
