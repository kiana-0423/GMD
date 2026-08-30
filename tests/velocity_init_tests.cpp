// Random velocity initialization: identity, determinism and order dependence.
//
// VelocityInitializer draws Maxwell-Boltzmann velocities for every atom. The
// question this file exists to ask is what decides WHICH draw a given physical
// atom receives.
//
// Today the answer is its position in the local storage array: one std::mt19937
// is advanced in a loop over atom_count(), so the first atom stored gets the
// generator's first draws. Under MPI each rank starts from the same seeded
// state and applies it to its own first atom, which is a different physical
// atom on every rank and at every rank count. The global temperature rescale
// then hides it -- total kinetic energy is forced to the target either way, so
// the reported temperature matches while the per-atom field does not.
//
// The tests below are in two groups.
//
//   PROPERTIES that hold whatever decides the draw: determinism for a fixed
//   seed, seed sensitivity, finiteness, non-triviality, the temperature and
//   centre-of-mass constraints, the zero-temperature case, and the 1/sqrt(mass)
//   scaling of the Gaussian width. These are asserted.
//
//   CHARACTERIZATION of the defect: the same physical system with its storage
//   permuted, compared by atom tag. This is measured and reported here rather
//   than asserted, because it fails today by construction; the commit that
//   makes the draw a function of the tag is the one that turns it into an
//   assertion, so that no commit in this series is red.
//
// Velocities are always compared BY TAG. Comparing by array index would make
// a permutation test vacuous -- it would compare atom i of one arrangement with
// a different physical atom i of the other.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "gmd/core/physical_constants.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[velocity init] " << message << '\n';
        ++failures;
    }
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

// k_B derived from the exact SI definitions rather than imported, as elsewhere
// in this repository's audits: a test that read the production constant would
// agree with a wrong production constant.
constexpr double kElementaryChargeCoulombs = 1.602176634e-19;  // exact, SI
constexpr double kBoltzmannJoulesPerKelvin = 1.380649e-23;     // exact, SI
const double kReferenceBoltzmannEVPerKelvin =
    kBoltzmannJoulesPerKelvin / kElementaryChargeCoulombs;

constexpr std::size_t kAtomCount = 24;

// Deliberately unequal masses spanning an order of magnitude, so that the
// 1/sqrt(mass) scaling has something to measure and a mass mix-up cannot be
// masked by uniformity.
double mass_for(std::size_t tag) {
    static const std::array<double, 8> kMasses = {
        1.008, 12.011, 15.999, 4.0026, 6.941, 10.811, 14.007, 32.065};
    return kMasses[tag % kMasses.size()];
}

std::array<double, 3> position_for(std::size_t tag) {
    const double t = static_cast<double>(tag);
    return {2.0 + 1.7 * t, 3.0 + 1.1 * std::fmod(t, 7.0), 4.0 + 1.3 * std::fmod(t, 5.0)};
}

// Builds the same physical system every time, with atoms stored in the given
// order. `order[i]` is the tag of the atom stored at slot i.
gmd::System system_with_order(const std::vector<int>& order) {
    gmd::System system;
    system.resize(order.size(), order.size());
    gmd::Box box;
    box.set_lengths({60.0, 60.0, 60.0});
    system.set_box(box);
    for (std::size_t slot = 0; slot < order.size(); ++slot) {
        const auto tag = static_cast<std::size_t>(order[slot]);
        system.mutable_masses()[slot] = mass_for(tag);
        system.mutable_coordinates()[slot] = position_for(tag);
        system.mutable_atom_tags()[slot] = order[slot];
    }
    return system;
}

std::vector<int> identity_order(std::size_t count) {
    std::vector<int> order(count);
    std::iota(order.begin(), order.end(), 0);
    return order;
}

// A fixed, non-trivial permutation: reversal composed with a rotation, so no
// atom keeps its slot and the mapping is not a simple shift.
std::vector<int> permuted_order(std::size_t count) {
    std::vector<int> order = identity_order(count);
    std::reverse(order.begin(), order.end());
    std::rotate(order.begin(), order.begin() + 5, order.end());
    return order;
}

using VelocityByTag = std::map<int, std::array<double, 3>>;

VelocityByTag velocities_by_tag(const gmd::System& system) {
    VelocityByTag field;
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const auto v = system.velocities()[i];
        field[system.atom_tag(i)] = {v[0], v[1], v[2]};
    }
    return field;
}

VelocityByTag initialize(const std::vector<int>& order, std::uint32_t seed,
                         double temperature, bool remove_com = true) {
    gmd::System system = system_with_order(order);
    gmd::VelocityInitializer initializer(seed);
    initializer.initialize(system, temperature, gmd::VelocityInitMode::Random, remove_com);
    return velocities_by_tag(system);
}

// Largest absolute component difference between two fields, matched by tag.
double worst_difference(const VelocityByTag& a, const VelocityByTag& b) {
    double worst = 0.0;
    for (const auto& [tag, va] : a) {
        const auto found = b.find(tag);
        if (found == b.end()) return std::numeric_limits<double>::infinity();
        for (std::size_t d = 0; d < 3; ++d) {
            worst = std::max(worst, std::abs(va[d] - found->second[d]));
        }
    }
    return worst;
}

double velocity_scale(const VelocityByTag& field) {
    double sum = 0.0;
    for (const auto& [tag, v] : field) {
        sum += v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    }
    return std::sqrt(sum / static_cast<double>(field.size()));
}

constexpr double kTargetTemperature = 300.0;
constexpr std::uint32_t kSeed = 20260830u;

// --- properties that must hold whatever decides the draw ------------------

void test_serial_determinism() {
    const auto order = identity_order(kAtomCount);
    const auto first = initialize(order, kSeed, kTargetTemperature);
    const auto second = initialize(order, kSeed, kTargetTemperature);

    check(worst_difference(first, second) == 0.0,
          "two initializations with the same seed and the same storage order "
          "produced different velocity fields");

    // A different seed must produce a genuinely different field, not a shifted
    // or rescaled one.
    const auto other_seed = initialize(order, kSeed + 1u, kTargetTemperature);
    const double difference = worst_difference(first, other_seed);
    check(difference > 1.0e-6,
          "changing the seed barely changed the velocity field (worst component "
          "difference " + number(difference) + ")");
}

void test_output_is_finite_and_nontrivial() {
    const auto field = initialize(identity_order(kAtomCount), kSeed, kTargetTemperature);
    check(field.size() == kAtomCount,
          "expected " + std::to_string(kAtomCount) + " atoms, got " +
              std::to_string(field.size()));

    std::size_t nonzero = 0;
    for (const auto& [tag, v] : field) {
        for (std::size_t d = 0; d < 3; ++d) {
            check(std::isfinite(v[d]),
                  "atom " + std::to_string(tag) + " component " + std::to_string(d) +
                      " is not finite");
            if (v[d] != 0.0) ++nonzero;
        }
    }
    check(nonzero == 3 * kAtomCount,
          "a nonzero-temperature initialization produced exactly-zero components; "
          + std::to_string(3 * kAtomCount - nonzero) + " of " +
          std::to_string(3 * kAtomCount) + " were zero");

    // Not all the same: a constant field would satisfy every check above.
    const auto& first = field.begin()->second;
    bool varies = false;
    for (const auto& [tag, v] : field) {
        if (v[0] != first[0] || v[1] != first[1] || v[2] != first[2]) varies = true;
    }
    check(varies, "every atom received the same velocity");
}

void test_temperature_and_centre_of_mass() {
    gmd::System system = system_with_order(identity_order(kAtomCount));
    gmd::VelocityInitializer initializer(kSeed);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, true);

    // The initializer's own convention: 3N-3 when the centre-of-mass velocity
    // has been removed, over the GLOBAL atom count.
    double twice_ke = 0.0;
    std::array<double, 3> momentum = {0.0, 0.0, 0.0};
    double total_mass = 0.0;
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const auto v = system.velocities()[i];
        const double m = system.masses()[i];
        twice_ke += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for (std::size_t d = 0; d < 3; ++d) momentum[d] += m * v[d];
        total_mass += m;
    }
    const double dof = 3.0 * static_cast<double>(kAtomCount) - 3.0;
    const double temperature = twice_ke / (dof * kReferenceBoltzmannEVPerKelvin);
    check(std::abs(temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "initialised to " + number(kTargetTemperature) + " K, the field carries " +
              number(temperature) + " K");

    // Momentum is removed before the rescale, and the rescale is a single
    // scalar, so it stays removed.
    const double scale = velocity_scale(velocities_by_tag(system)) * total_mass;
    for (std::size_t d = 0; d < 3; ++d) {
        check(std::abs(momentum[d]) < 1.0e-12 * scale,
              "centre-of-mass momentum component " + std::to_string(d) + " is " +
                  number(momentum[d]) + ", not zero");
    }
}

void test_no_centre_of_mass_removal_when_disabled() {
    gmd::System system = system_with_order(identity_order(kAtomCount));
    gmd::VelocityInitializer initializer(kSeed);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, false);

    std::array<double, 3> momentum = {0.0, 0.0, 0.0};
    double twice_ke = 0.0;
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const auto v = system.velocities()[i];
        const double m = system.masses()[i];
        for (std::size_t d = 0; d < 3; ++d) momentum[d] += m * v[d];
        twice_ke += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    const double magnitude = std::sqrt(momentum[0] * momentum[0] +
                                       momentum[1] * momentum[1] +
                                       momentum[2] * momentum[2]);
    check(magnitude > 0.0,
          "with centre-of-mass removal disabled the residual momentum is exactly "
          "zero, which means it was removed anyway");

    // And the degrees of freedom change to 3N, which is a different rescale.
    const double dof = 3.0 * static_cast<double>(kAtomCount);
    const double temperature = twice_ke / (dof * kReferenceBoltzmannEVPerKelvin);
    check(std::abs(temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "with removal disabled the field should carry " + number(kTargetTemperature) +
              " K over 3N degrees of freedom; it carries " + number(temperature));
}

void test_zero_temperature() {
    gmd::System system = system_with_order(identity_order(kAtomCount));
    gmd::VelocityInitializer initializer(kSeed);
    initializer.initialize(system, 0.0, gmd::VelocityInitMode::Random, true);
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const auto v = system.velocities()[i];
        check(v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0,
              "atom " + std::to_string(system.atom_tag(i)) +
                  " is not exactly at rest at 0 K: " + number(v[0]) + ", " +
                  number(v[1]) + ", " + number(v[2]));
    }
}

void test_gaussian_width_scales_with_inverse_root_mass() {
    // The rescale is one global scalar and the centre-of-mass shift is one
    // global vector, so the RATIO of the per-mass velocity spreads survives
    // both. sigma = sqrt(k_B T / m), so atoms of mass m1 and m2 must have
    // spreads in the ratio sqrt(m2/m1).
    //
    // A large sample is used because this is a statistical property of the
    // draw. It is a weak test by construction; the exact per-atom check that
    // replaces it arrives with the keyed generator, which makes each atom's
    // pre-rescale value predictable.
    constexpr std::size_t kSampleAtoms = 4000;
    std::vector<int> order = identity_order(kSampleAtoms);
    gmd::System system;
    system.resize(kSampleAtoms, kSampleAtoms);
    gmd::Box box;
    box.set_lengths({400.0, 400.0, 400.0});
    system.set_box(box);
    const double light = 1.008;
    const double heavy = 32.065;
    for (std::size_t i = 0; i < kSampleAtoms; ++i) {
        system.mutable_masses()[i] = (i % 2 == 0) ? light : heavy;
        system.mutable_coordinates()[i] = {0.05 * static_cast<double>(i), 1.0, 2.0};
        system.mutable_atom_tags()[i] = static_cast<int>(i);
    }
    gmd::VelocityInitializer initializer(kSeed);
    initializer.initialize(system, kTargetTemperature, gmd::VelocityInitMode::Random, true);

    double light_sum = 0.0, heavy_sum = 0.0;
    std::size_t light_n = 0, heavy_n = 0;
    for (std::size_t i = 0; i < kSampleAtoms; ++i) {
        const auto v = system.velocities()[i];
        const double square = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if (i % 2 == 0) { light_sum += square; ++light_n; }
        else            { heavy_sum += square; ++heavy_n; }
    }
    const double measured = std::sqrt((light_sum / static_cast<double>(light_n)) /
                                      (heavy_sum / static_cast<double>(heavy_n)));
    const double expected = std::sqrt(heavy / light);
    std::cout << "  sigma ratio light/heavy  " << number(measured)
              << "  (sqrt(m_heavy/m_light) = " << number(expected) << ")\n";
    // 2000 atoms per mass gives a relative standard error on the ratio of
    // roughly 1/sqrt(2*3*2000) ~ 0.9 percent; 5 percent is a comfortable bound
    // that still rejects a missing or inverted mass factor by a wide margin.
    check(std::abs(measured / expected - 1.0) < 0.05,
          "the Gaussian width does not scale as 1/sqrt(mass): measured ratio " +
              number(measured) + " against " + number(expected));
}

// --- characterization of the defect ---------------------------------------

void report_storage_order_dependence() {
    const auto ordered = initialize(identity_order(kAtomCount), kSeed, kTargetTemperature);
    const auto permuted = initialize(permuted_order(kAtomCount), kSeed, kTargetTemperature);

    std::size_t differing = 0;
    for (const auto& [tag, v] : ordered) {
        const auto& w = permuted.at(tag);
        if (v[0] != w[0] || v[1] != w[1] || v[2] != w[2]) ++differing;
    }
    const double worst = worst_difference(ordered, permuted);

    std::cout << "  storage-order dependence  " << differing << " of "
              << kAtomCount << " atoms differ, worst component "
              << std::scientific << std::setprecision(3) << worst
              << std::defaultfloat << '\n';

    // Reported, not asserted. Every atom is present under both arrangements and
    // both fields are physically valid; what differs is which draw each atom
    // received. The commit that keys the draw on the tag turns this into an
    // equality assertion.
    check(ordered.size() == permuted.size(),
          "permuting storage changed how many atoms were initialized");
    for (const auto& [tag, v] : ordered) {
        check(permuted.count(tag) == 1,
              "atom tag " + std::to_string(tag) +
                  " disappeared when storage was permuted");
    }
}

}  // namespace

int main() {
    std::cout << "[velocity init] auditing random velocity initialization\n";
    test_serial_determinism();
    test_output_is_finite_and_nontrivial();
    test_temperature_and_centre_of_mass();
    test_no_centre_of_mass_removal_when_disabled();
    test_zero_temperature();
    test_gaussian_width_scales_with_inverse_root_mass();
    report_storage_order_dependence();

    if (failures == 0) {
        std::cout << "[velocity init] audit passed\n";
        return 0;
    }
    std::cerr << "[velocity init] " << failures << " check(s) failed\n";
    return 1;
}
