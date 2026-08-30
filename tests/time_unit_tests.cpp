// Audit of GMD's internal time unit.
//
// GMD works in eV, Angstrom and atomic mass units, and it integrates
//
//     v += (F / m) * dt        r += v * dt
//
// with F in eV/A and m in amu. Nothing in those two lines mentions seconds, so
// the time unit is not free: it is fixed by requiring that F/m, which carries
// units of eV/(A*amu), be an acceleration in A per (time unit) squared:
//
//     eV / (A * amu) = A / T^2   =>   T = A * sqrt(amu / eV)
//
// Users never see T. They give a timestep in femtoseconds and read a time
// column in femtoseconds, so exactly one conversion stands between them, and it
// is the only place the value appears.
//
// WHAT THIS FILE MEASURES RATHER THAN ASSUMES
//
//   ballistic     a free particle carrying one internal velocity unit covers a
//                 measurable distance per femtosecond; that distance IS the
//                 velocity unit, and its reciprocal is T.
//
//   harmonic      a harmonic oscillator of known mass and force constant has an
//                 analytically exact period, 2*pi*sqrt(m/k) in internal units.
//                 Measured in femtoseconds and divided by that, the answer is T
//                 again -- this time from a real physical observable rather than
//                 from a displacement.
//
//   thermostat    the Nose-Hoover mass is Q = dof*k_B*T_target*tau^2, and for an
//                 ideal gas the linearised thermostat oscillates the temperature
//                 with period 2*pi*tau/sqrt(2). Both are exposed, so the tau the
//                 thermostat actually used is recoverable in femtoseconds and can
//                 be compared against the tau that was asked for.
//
//   barostat      the Berendsen coupling factor mu is observable from the box, so
//                 mu^3 = 1 - beta*(dt/tau)*(P_target - P) inverts to the tau the
//                 barostat actually used.
//
// THE AUTHORITATIVE VALUE
//
// Unlike the Boltzmann constant and the pressure conversion, this one is not
// exact: it inherits the uncertainty of the atomic mass constant. Its
// ingredients are
//
//     e   = 1.602176634e-19 J per eV     exact, SI 2019
//           https://physics.nist.gov/cgi-bin/cuu/Value?evj
//     m_u = 1.66053906892(52)e-27 kg     CODATA 2022, relative 3.1e-10
//           https://physics.nist.gov/cgi-bin/cuu/Value?ukg
//     1 A = 1e-10 m                      exact, by definition
//     1 fs = 1e-15 s                     exact, by definition
//
// so
//
//     T = 1e-10 * sqrt(m_u / e) s = sqrt(m_u / e) * 1e5 fs
//       = 10.1805057178711931077510010336...  fs
//
// with a relative uncertainty of 1.57e-10, half that of m_u because of the
// square root. The reference below is written as that expression rather than as
// a number, so the derivation is visible and checkable instead of copied.

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <filesystem>
#include <fstream>
#include <vector>

#include "gmd/core/physical_constants.hpp"
#include "gmd/core/runtime_context.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/berendsen_barostat.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/io/trajectory_writer.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[time units] " << message << '\n';
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
// Derived here from the SI definitions rather than imported: a test that read
// the production constant would agree with a wrong production constant.
//
// T[fs] = 1e-10 m/A * sqrt(m_u / e) / 1e-15 s/fs = sqrt(m_u / e) * 1e5.
// Grouping it that way matters: routing the 1e-10 and the 1e-15 through
// separate multiplications rounds twice and lands one ulp low. Written as a
// single square root scaled by 1e5, the result is the nearest double to the
// exact value.
constexpr double kElementaryChargeCoulombs = 1.602176634e-19;   // exact, SI
constexpr double kAtomicMassConstantKg     = 1.66053906892e-27; // CODATA 2022

const double kReferenceInternalTimeUnitFs =
    std::sqrt(kAtomicMassConstantKg / kElementaryChargeCoulombs) * 1.0e5;

// The velocity unit is its reciprocal: one internal velocity unit is this many
// Angstrom per femtosecond. Numerically sqrt(eV/amu) expressed in A/fs.
const double kReferenceVelocityUnitAPerFs = 1.0 / kReferenceInternalTimeUnitFs;

// What the production constant was before this audit: the same quantity rounded
// to seven significant figures.
constexpr double kSupersededValue = 1.018051e+1;

// Converts a femtosecond timestep the way config_loader does, so that every
// fixture here is driven exactly as a real run would be. Takes the conversion
// as an argument rather than reading production, so the audit can drive the
// integrator with its OWN reference and ask what the integrator then does.
double internal_dt(double dt_fs, double time_unit_fs) {
    return dt_fs / time_unit_fs;
}

void test_reference_derivation_is_self_consistent() {
    // The naive route: separate factors for the metre-to-Angstrom and
    // second-to-femtosecond conversions. Not the definition -- it rounds more
    // than once -- but an independent path to the same value, and what would
    // catch a transcription error in the two SI literals.
    const double naive =
        1.0e-10 * std::sqrt(kAtomicMassConstantKg / kElementaryChargeCoulombs) / 1.0e-15;
    check(relative_difference(naive, kReferenceInternalTimeUnitFs) < 4.0e-16,
          "the grouped reference " + number(kReferenceInternalTimeUnitFs) +
              " disagrees with the separated-factor route " + number(naive));

    // Time unit and velocity unit must be exact reciprocals, or they are not
    // descriptions of the same unit system.
    check(kReferenceInternalTimeUnitFs * kReferenceVelocityUnitAPerFs == 1.0,
          "the time and velocity units are not exact reciprocals: product " +
              number(kReferenceInternalTimeUnitFs * kReferenceVelocityUnitAPerFs));

    // Sanity on the magnitude: a vibration period of a light atom is tens of
    // femtoseconds, so the unit must be of order 10 fs, not 1 or 100.
    check(kReferenceInternalTimeUnitFs > 10.0 && kReferenceInternalTimeUnitFs < 11.0,
          "the derived internal time unit is not of order 10 fs: " +
              number(kReferenceInternalTimeUnitFs));
}

// --- force providers -------------------------------------------------------

class ZeroForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "zero_force"; }
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

// A fixed force on every atom, in eV/A.
class ConstantForceProvider final : public gmd::ForceProvider {
public:
    explicit ConstantForceProvider(std::array<double, 3> force) : force_(force) {}
    std::string_view name() const noexcept override { return "constant_force"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request, gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.forces.assign(request.coordinates.size(), force_);
        result.potential_energy = 0.0;
        result.virial_valid = false;
        result.success = true;
    }
private:
    std::array<double, 3> force_;
};

// F = -k (r - r0), an isotropic harmonic well in eV/A with k in eV/A^2.
class HarmonicProvider final : public gmd::ForceProvider {
public:
    HarmonicProvider(double k, std::array<double, 3> centre) : k_(k), centre_(centre) {}
    std::string_view name() const noexcept override { return "harmonic"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request, gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        const std::size_t n = request.coordinates.size();
        result.forces.assign(n, {0.0, 0.0, 0.0});
        double energy = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t d = 0; d < 3; ++d) {
                const double displacement = request.coordinates[i][d] - centre_[d];
                result.forces[i][d] = -k_ * displacement;
                energy += 0.5 * k_ * displacement * displacement;
            }
        }
        result.potential_energy = energy;
        result.virial_valid = false;
        result.success = true;
    }
private:
    double k_;
    std::array<double, 3> centre_;
};

gmd::System one_atom(double mass, std::array<double, 3> position,
                     std::array<double, 3> velocity) {
    gmd::System system;
    system.resize(1, 1);
    gmd::Box box;
    // Large, and every fixture below sits near its middle. The integrator wraps
    // coordinates into the cell each step, so a particle placed at a negative
    // coordinate -- or one that oscillates through zero -- reappears at the far
    // face and every displacement measured from it is nonsense. Keeping the
    // motion far from both faces is what makes the free-flight and oscillator
    // identities readable.
    box.set_lengths({1000.0, 1000.0, 1000.0});
    system.set_box(box);
    system.mutable_masses()[0] = mass;
    system.mutable_coordinates()[0] = position;
    system.mutable_velocities()[0] = velocity;
    return system;
}

// Evaluates forces once and installs them on the System.
//
// Simulation::initialize() does this before the first step, and it is not
// optional: begin_step() opens with a half-kick that reads System::forces(), so
// a run that starts with them zeroed loses half of the first step's impulse.
// Velocity Verlet is exact for a constant force only when it is seeded this way.
void seed_forces(gmd::System& system, gmd::ForceProvider& provider,
                 gmd::RuntimeContext& runtime) {
    const auto coordinates = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .step = 0,
        .time = 0.0,
        .coordinates = std::span<const gmd::Coordinate3D>(coordinates.data(),
                                                          coordinates.size()),
        .neighbor_list = nullptr,
    };
    gmd::ForceResult result;
    provider.compute(request, result, runtime);
    auto forces = system.mutable_forces();
    for (std::size_t i = 0; i < forces.size() && i < result.forces.size(); ++i) {
        forces[i] = result.forces[i];
    }
}

// Runs `steps` integrator steps at the given internal timestep.
void run(gmd::System& system, gmd::ForceProvider& provider, double dt, int steps) {
    gmd::RuntimeContext runtime;
    gmd::VelocityVerletIntegrator integrator(dt);
    integrator.initialize(system, runtime);
    seed_forces(system, provider, runtime);
    for (int i = 0; i < steps; ++i) {
        const gmd::IntegratorStepContext ctx{.step = static_cast<std::uint64_t>(i),
                                             .dt = dt};
        integrator.step(system, provider, ctx, runtime);
    }
}

// --- ballistic: the velocity unit, measured -------------------------------

void test_ballistic_motion_measures_the_velocity_unit() {
    // One atom, one internal velocity unit along x, no force. After t
    // femtoseconds it has moved (1 velocity unit) * t, so the displacement per
    // femtosecond IS the velocity unit in A/fs.
    const double dt_fs = 1.0;
    const int steps = 100;
    const double dt = internal_dt(dt_fs, kReferenceInternalTimeUnitFs);

    gmd::System system = one_atom(1.0, {64.0, 80.0, 96.0}, {1.0, 0.0, 0.0});
    ZeroForceProvider provider;
    run(system, provider, dt, steps);

    const double displacement = system.coordinates()[0][0] - 64.0;
    const double measured_velocity_unit = displacement / (dt_fs * steps);

    std::cout << "  velocity unit  " << number(measured_velocity_unit)
              << " A/fs  (reference " << number(kReferenceVelocityUnitAPerFs) << ")\n";

    // Driven with the reference conversion the integrator must reproduce it
    // exactly -- this is the free-flight identity r = v t and nothing else, so
    // any deviation is an integrator defect rather than a unit question.
    // Not a physics tolerance: the coordinate is of order 64 A while the
    // displacement being measured is ~10 A, so a hundred accumulated additions
    // leave a few hundred ulp of the coordinate in the difference, which is
    // ~1e-14 relative. The identity itself is exact.
    check(relative_difference(measured_velocity_unit, kReferenceVelocityUnitAPerFs) < 1.0e-12,
          "a free particle carrying one internal velocity unit moved " +
              number(measured_velocity_unit) + " A/fs, not " +
              number(kReferenceVelocityUnitAPerFs));

    // Motion stays on the axis it started on.
    check(system.coordinates()[0][1] == 80.0 && system.coordinates()[0][2] == 96.0,
          "a particle moving along x acquired transverse displacement");
}

// --- constant force: t and t^2 scaling ------------------------------------

void test_constant_force_scaling() {
    // Under a constant force, exactly: v(t) = v0 + a t and
    // r(t) = r0 + v0 t + 1/2 a t^2. Velocity Verlet reproduces both to machine
    // precision for a constant force, so these are equalities and not
    // approximations -- which makes this a test of the integrator's arithmetic
    // rather than of its accuracy.
    const double mass = 12.011;
    const std::array<double, 3> force = {0.35, -0.7, 0.125};
    const double dt_fs = 0.5;
    const double dt = internal_dt(dt_fs, kReferenceInternalTimeUnitFs);
    const std::array<double, 3> r0 = {64.0, 80.25, 96.75};
    const std::array<double, 3> v0 = {0.02, 0.011, -0.004};

    for (int steps : {1, 2, 10, 200}) {
        gmd::System system = one_atom(mass, r0, v0);
        ConstantForceProvider provider(force);
        run(system, provider, dt, steps);

        const double t = dt * steps;   // internal units throughout
        for (std::size_t d = 0; d < 3; ++d) {
            const double acceleration = force[d] / mass;
            const double expected_v = v0[d] + acceleration * t;
            const double expected_r = r0[d] + v0[d] * t + 0.5 * acceleration * t * t;
            const double got_v = system.velocities()[0][d];
            const double got_r = system.coordinates()[0][d];
            check(relative_difference(got_v, expected_v) < 1.0e-13,
                  "after " + std::to_string(steps) + " steps, velocity component " +
                      std::to_string(d) + " is " + number(got_v) + ", not " +
                      number(expected_v));
            check(relative_difference(got_r, expected_r) < 1.0e-12,
                  "after " + std::to_string(steps) + " steps, position component " +
                      std::to_string(d) + " is " + number(got_r) + ", not " +
                      number(expected_r));
        }
    }

    // Reversing the force reverses the acceleration and nothing else: the
    // displacement about the free-flight line must be exactly opposite.
    const int steps = 50;
    const double t = dt * steps;
    gmd::System forward = one_atom(mass, r0, v0);
    gmd::System reverse = one_atom(mass, r0, v0);
    ConstantForceProvider push(force);
    ConstantForceProvider pull({-force[0], -force[1], -force[2]});
    run(forward, push, dt, steps);
    run(reverse, pull, dt, steps);
    for (std::size_t d = 0; d < 3; ++d) {
        const double free_flight = r0[d] + v0[d] * t;
        const double a = forward.coordinates()[0][d] - free_flight;
        const double b = reverse.coordinates()[0][d] - free_flight;
        // Same round-off budget as above: the displacement about the
        // free-flight line is ~1 A recovered from coordinates of order 100.
        check(relative_difference(a, -b) < 1.0e-10,
              "reversing the force did not reverse the displacement in component " +
                  std::to_string(d) + ": " + number(a) + " vs " + number(b));
    }
}

// --- harmonic oscillator: the time unit from a physical observable ---------

// Measures the oscillation period in femtoseconds by linear interpolation of
// upward zero crossings of x(t) - centre.
double measure_period_fs(double mass, double k, double dt_fs, double time_unit_fs,
                         int steps) {
    const double dt = internal_dt(dt_fs, time_unit_fs);
    const std::array<double, 3> centre = {128.0, 128.0, 128.0};
    gmd::System system = one_atom(mass, {129.0, 128.0, 128.0}, {0.0, 0.0, 0.0});
    HarmonicProvider provider(k, centre);
    gmd::RuntimeContext runtime;
    gmd::VelocityVerletIntegrator integrator(dt);
    integrator.initialize(system, runtime);
    seed_forces(system, provider, runtime);

    std::vector<double> crossings;
    double previous = system.coordinates()[0][0] - centre[0];
    for (int i = 0; i < steps; ++i) {
        const gmd::IntegratorStepContext ctx{.step = static_cast<std::uint64_t>(i),
                                             .dt = dt};
        integrator.step(system, provider, ctx, runtime);
        const double current = system.coordinates()[0][0] - centre[0];
        if (previous < 0.0 && current >= 0.0) {
            // Linear interpolation between the two straddling samples, so the
            // resolution is not the timestep.
            const double fraction = -previous / (current - previous);
            crossings.push_back((static_cast<double>(i) + fraction) * dt_fs);
        }
        previous = current;
    }
    if (crossings.size() < 2) return 0.0;
    return (crossings.back() - crossings.front()) /
           static_cast<double>(crossings.size() - 1);
}

void test_harmonic_period_measures_the_time_unit() {
    // A harmonic oscillator's period is 2*pi*sqrt(m/k) in whatever time unit the
    // integrator works in. Measure it in femtoseconds and the ratio is the time
    // unit itself -- an observable that depends on mass, force and time together,
    // so it tests the whole unit system rather than one conversion.
    const double mass = 12.011;     // amu
    const double k = 1.0;           // eV/A^2
    const double dt_fs = 0.02;
    const double period_internal = 2.0 * M_PI * std::sqrt(mass / k);
    const double expected_fs = period_internal * kReferenceInternalTimeUnitFs;

    const double measured = measure_period_fs(mass, k, dt_fs, kReferenceInternalTimeUnitFs,
                                              200000);
    check(measured > 0.0, "no oscillation period could be measured");

    std::cout << "  harmonic period " << number(measured) << " fs  (analytic "
              << number(expected_fs) << ")\n";

    // Velocity Verlet's period error is O(dt^2): at dt_fs = 0.02 against a
    // ~221 fs period the fractional error is about (omega*dt)^2/24 ~ 1.4e-8.
    check(relative_difference(measured, expected_fs) < 1.0e-6,
          "the harmonic period is " + number(measured) + " fs where mass " +
              number(mass) + " amu and k " + number(k) +
              " eV/A^2 require " + number(expected_fs) + " fs");
}

void test_harmonic_convergence_is_second_order() {
    // Halving the timestep must quarter the period error. This is what
    // distinguishes "the time unit is right" from "the integrator happens to be
    // accurate at one step size".
    const double mass = 12.011;
    const double k = 1.0;
    const double expected = 2.0 * M_PI * std::sqrt(mass / k) * kReferenceInternalTimeUnitFs;

    double previous_error = 0.0;
    for (double dt_fs : {0.32, 0.16, 0.08}) {
        const int steps = static_cast<int>(40.0 * expected / dt_fs);
        const double measured = measure_period_fs(mass, k, dt_fs,
                                                  kReferenceInternalTimeUnitFs, steps);
        const double error = std::abs(measured - expected);
        if (previous_error > 0.0) {
            const double ratio = previous_error / error;
            check(ratio > 3.5 && ratio < 4.5,
                  "halving the timestep changed the period error by a factor of " +
                      number(ratio) + ", not the ~4 a second-order integrator "
                      "requires");
        }
        previous_error = error;
    }
}

// --- the reported time column ---------------------------------------------
//
// The log header says time[fs], and app/gmd_main.cpp supplies step *
// time_step_fs, so the femtosecond value is computed by the caller and the
// writer must pass it through untouched. A conversion applied a second time
// inside the writer would be invisible to every test that only reads back what
// it wrote, so this one states the number it expects.
//
// Note that Simulation passes force_time in INTERNAL units to the force
// provider's ForceRequest, which is a different quantity from this column and
// is not what any log or checkpoint reports.

void test_trajectory_time_column_is_passed_through_in_femtoseconds() {
    gmd::System system;
    system.resize(1, 1);
    gmd::Box box;
    box.set_lengths({20.0, 20.0, 20.0});
    system.set_box(box);
    system.mutable_masses()[0] = 12.011;
    system.mutable_coordinates()[0] = {1.0, 2.0, 3.0};

    const double dt_fs = 2.5;
    const std::filesystem::path stem =
        std::filesystem::temp_directory_path() / "gmd_time_unit_probe";
    gmd::TrajectoryWriter writer;
    writer.open(stem);
    for (std::uint64_t step = 0; step < 4; ++step) {
        writer.write_frame(system, step, static_cast<double>(step) * dt_fs, 0.0, 3);
    }
    writer.close();

    std::ifstream input(stem.string() + ".log");
    std::string line;
    std::uint64_t expected_step = 0;
    int checked = 0;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        double step = 0.0, time = 0.0;
        fields >> step >> time;
        const double expected = static_cast<double>(expected_step) * dt_fs;
        check(std::abs(time - expected) < 1.0e-9,
              "frame " + std::to_string(expected_step) + " reports time " +
                  number(time) + " where step * dt_fs is " + number(expected));
        ++expected_step;
        ++checked;
    }
    std::filesystem::remove(stem.string() + ".log");
    std::filesystem::remove(stem.string() + ".xyz");
    check(checked == 4, "expected four logged frames, found " + std::to_string(checked));
}

// --- thermostat and barostat time scales, measured ------------------------

// Recovers the tau the Nose-Hoover thermostat actually used, in femtoseconds,
// from the thermostat mass it built. Q = dof * k_B * T * tau^2 is exposed at
// full precision through checkpoint_state(), so this inverts exactly.
double measure_nose_hoover_tau_fs(double requested_tau, double time_unit_fs) {
    constexpr std::size_t kAtoms = 8;
    constexpr double kTarget = 300.0;
    gmd::System system;
    system.resize(kAtoms, kAtoms);
    gmd::Box box;
    box.set_lengths({50.0, 50.0, 50.0});
    system.set_box(box);
    for (std::size_t i = 0; i < kAtoms; ++i) {
        system.mutable_masses()[i] = 12.011;
        const double t = static_cast<double>(i);
        system.mutable_coordinates()[i] = {2.0 + 3.0 * t, 3.0 + 1.5 * t, 4.0 + 2.0 * t};
        system.mutable_velocities()[i] = {0.01, -0.005, 0.007};
    }
    auto thermostat = std::make_shared<gmd::NoseHooverThermostat>(requested_tau);
    ZeroForceProvider provider;
    gmd::RuntimeContext runtime;
    const double dt = internal_dt(1.0, time_unit_fs);
    gmd::VelocityVerletIntegrator integrator(dt);
    integrator.set_thermostat(thermostat);
    integrator.set_target_temperature(kTarget);
    integrator.initialize(system, runtime);
    const gmd::IntegratorStepContext ctx{.step = 0, .dt = dt};
    integrator.step(system, provider, ctx, runtime);

    // Parse Q and dof back out of the checkpoint string.
    std::istringstream state(thermostat->checkpoint_state());
    std::string key;
    double tau = 0.0, xi = 0.0, q = 0.0, temperature = 0.0;
    std::size_t dof = 0;
    state >> key >> tau >> key >> xi >> key >> q >> key >> dof >> key >> temperature;

    // k_B derived independently, as elsewhere in this repository's audits.
    constexpr double kBoltzmannJoulesPerKelvin = 1.380649e-23;
    const double boltzmann = kBoltzmannJoulesPerKelvin / kElementaryChargeCoulombs;
    // Q = dof * k_B * T * tau_internal^2  =>  tau_internal = sqrt(Q / (dof k_B T))
    const double tau_internal =
        std::sqrt(q / (static_cast<double>(dof) * boltzmann * kTarget));
    return tau_internal * time_unit_fs;
}

// Recovers the tau the Berendsen barostat actually used, in femtoseconds, by
// inverting its observable coupling factor.
double measure_berendsen_tau_fs(double requested_tau, double time_unit_fs) {
    constexpr double kBoxLength = 12.0;
    constexpr double kBeta = 4.5e-5;
    constexpr double kTargetPressureBar = 5000.0;
    const double dt_fs = 1.0;
    const double dt = internal_dt(dt_fs, time_unit_fs);

    gmd::System system;
    system.resize(4, 4);
    gmd::Box box;
    box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    system.set_box(box);
    for (std::size_t i = 0; i < 4; ++i) {
        const double t = static_cast<double>(i);
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {1.0 + 2.0 * t, 2.0 + 0.5 * t, 3.0 + 0.25 * t};
        system.mutable_velocities()[i] = {0.004, -0.002, 0.003};
    }
    std::array<double, 9> virial{};
    virial[0] = 0.5;
    virial[4] = 0.3;
    virial[8] = 0.2;
    system.set_last_virial(virial, true);

    double twice_ke = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        const auto v = system.velocities()[i];
        twice_ke += system.masses()[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    const double trace = virial[0] + virial[4] + virial[8];
    const double volume = kBoxLength * kBoxLength * kBoxLength;
    const double pressure_internal = (twice_ke + trace) / (3.0 * volume);
    // The barostat compares in bar; 1 eV/A^3 = 1602176.634 bar exactly. Derived
    // here rather than imported, as in tests/pressure_unit_tests.cpp.
    const double pressure_bar = pressure_internal * (801088317.0 / 500.0);

    gmd::BerendsenBarostat barostat(requested_tau, kBeta);
    ZeroForceProvider provider;
    gmd::RuntimeContext runtime;
    barostat.apply(system, provider, runtime, 0, dt, 300.0, kTargetPressureBar, trace);

    const double mu = system.box().lengths[0] / kBoxLength;
    const double mu_cubed = mu * mu * mu;
    // mu^3 = 1 - beta * (dt / tau_used) * (P_target - P), everything else known.
    const double tau_used_internal =
        kBeta * dt * (kTargetPressureBar - pressure_bar) / (1.0 - mu_cubed);
    return tau_used_internal * time_unit_fs;
}

void report_relaxation_time_scales() {
    // Measured and reported here, not asserted. The audit commit establishes the
    // measurement; the commit that converts these times is the one that pins
    // them, so that no commit in this series is red.
    const double requested = 100.0;
    const double nose_hoover =
        measure_nose_hoover_tau_fs(requested, kReferenceInternalTimeUnitFs);
    const double berendsen =
        measure_berendsen_tau_fs(2000.0, kReferenceInternalTimeUnitFs);
    std::cout << "  Nose-Hoover tau requested " << number(requested)
              << " fs, effective " << number(nose_hoover) << " fs  (ratio "
              << number(nose_hoover / requested) << ")\n";
    std::cout << "  Berendsen   tau requested " << number(2000.0)
              << " fs, effective " << number(berendsen) << " fs  (ratio "
              << number(berendsen / 2000.0) << ")\n";
}

// --- the production conversion must be the authoritative one ---------------
//
// The checks above drive the integrator with this file's own reference, so they
// hold whatever production uses. This one asks what production actually does
// with a femtosecond timestep, by running the same free-flight measurement
// through the real config-loader conversion.

void test_production_conversion_matches_the_reference() {
    // config_loader computes time_step = time_step_fs * kInternalTimePerFemtosecond.
    // Reproducing that here and re-measuring the velocity unit says whether a
    // run driven by a real input file sees the right time base.
    const double dt_fs = 1.0;
    const int steps = 100;
    const double production_dt = dt_fs * gmd::kInternalTimePerFemtosecond;

    gmd::System system = one_atom(1.0, {64.0, 80.0, 96.0}, {1.0, 0.0, 0.0});
    ZeroForceProvider provider;
    run(system, provider, production_dt, steps);
    const double measured =
        (system.coordinates()[0][0] - 64.0) / (dt_fs * steps);

    check(relative_difference(measured, kReferenceVelocityUnitAPerFs) < 1.0e-12,
          "a run using the production femtosecond conversion measures a velocity "
          "unit of " + number(measured) + " A/fs, not " +
              number(kReferenceVelocityUnitAPerFs));

    // The superseded value, named, so a revert reads as a revert rather than as
    // an anonymous tolerance failure.
    const double superseded_unit = 1.0 / kSupersededValue;
    check(relative_difference(measured, superseded_unit) > 1.0e-7,
          "the production conversion still uses the superseded 1.018051e+1, "
          "which is 4.206204e-07 relative high");

    // And the two directions production exposes must be exact reciprocals, as
    // this file's own pair is.
    check(gmd::kFemtosecondsPerInternalTime * gmd::kInternalTimePerFemtosecond == 1.0,
          "the production time conversions are not exact reciprocals");
    check(relative_difference(gmd::kFemtosecondsPerInternalTime,
                              kReferenceInternalTimeUnitFs) < 1.0e-15,
          "the production internal time unit is " +
              number(gmd::kFemtosecondsPerInternalTime) + " fs, not " +
              number(kReferenceInternalTimeUnitFs));
}

void report_production_constant() {
    std::cout << "  reference T   " << number(kReferenceInternalTimeUnitFs)
              << " fs\n"
              << "  production    " << number(gmd::kFemtosecondsPerInternalTime)
              << " fs\n"
              << "  superseded    " << number(kSupersededValue) << " fs  (relative "
              << std::scientific << std::setprecision(6)
              << (kSupersededValue - kReferenceInternalTimeUnitFs) /
                     kReferenceInternalTimeUnitFs
              << std::defaultfloat << ")\n";
}

}  // namespace

int main() {
    std::cout << "[time units] auditing the internal time unit T = A*sqrt(amu/eV)\n";
    test_reference_derivation_is_self_consistent();
    test_ballistic_motion_measures_the_velocity_unit();
    test_constant_force_scaling();
    test_harmonic_period_measures_the_time_unit();
    test_harmonic_convergence_is_second_order();
    test_trajectory_time_column_is_passed_through_in_femtoseconds();
    test_production_conversion_matches_the_reference();
    report_relaxation_time_scales();
    report_production_constant();

    if (failures == 0) {
        std::cout << "[time units] audit passed\n";
        return 0;
    }
    std::cerr << "[time units] " << failures << " check(s) failed\n";
    return 1;
}
