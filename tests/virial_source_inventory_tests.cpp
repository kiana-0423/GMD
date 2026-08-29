// Independent validation of every virial source in GMD, component by component.
//
// The existing tests/virial_finite_difference_tests.cpp validates tr(W) and the
// three diagonal components against numerical strain derivatives of the
// engine's own energy. It cannot validate the off-diagonal components, because
// `gmd::Box` holds three edge lengths and no shear strain is expressible.
//
// This file closes that gap without pretending the engine can shear. Each
// source is held against a reference built from its own definition:
//
//   * pair terms (LJ, Ewald/PME real space, special-pair Coulomb) against the
//     exact analytical pair identity W_ab = r_a F_b, with the force recomputed
//     here from the potential and cross-checked against a central difference
//     of that potential;
//   * bonded terms against sum_i (r_i - r_ref)_a F_i,b, where the per-atom
//     forces are first verified against a central difference of an
//     independently written energy function, and the relative coordinates are
//     built here by walking the interaction chain. The production helper
//     accum_interaction_virial() is deliberately not used;
//   * the Ewald reciprocal term against BOTH an independently derived analytic
//     tensor and the numerical strain derivative of a test-only reciprocal
//     energy written for a general 3x3 cell -- which, unlike `gmd::Box`, can be
//     sheared, so all nine components are reached;
//   * the net-charge correction against its analytically expected isotropic
//     tensor, and the self term against the claim that it contributes nothing.
//
// EVIDENCE LABELS. Every check below is tagged in its label with the class of
// evidence it rests on, so the coverage table in validation/README.md can be
// read back off the test output:
//
//   [analytic-pair]   exact analytical pair reference
//   [bonded-moment]   independent bonded force moment
//   [ewald-derived]   independent Ewald reciprocal derivation
//   [rotation]        rotation covariance
//   [symmetry]        regression-only symmetry / invariance
//
// The MLForceProvider contract lives in tests/ml_virial_contract_tests.cpp,
// because it is a claim about a provider that reports NO virial and belongs
// with the fix that made that contract explicit rather than with the tensors.
//
// ROTATION COVARIANCE IS NOT A SHEAR DERIVATIVE. W -> R W R^T is a necessary
// condition on any Cartesian rank-2 tensor and it does constrain off-diagonal
// components, but it is not equivalent to differentiating the energy with
// respect to a shear strain: it tests how the tensor transforms, not that it is
// the derivative of anything. It is also only applicable where the periodic
// cell plays no role, since rotating a configuration inside a fixed
// orthorhombic box does not rotate the lattice. Where the lattice matters --
// the reciprocal-space terms -- the strain-derivative reference in
// tests/virial_reference.hpp is what carries the claim.

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/bonded_force_provider.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/composite_force_provider.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/io/config_loader.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/special_pair_map.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

#include "virial_reference.hpp"

namespace vr = virial_ref;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[virial inventory] " << message << '\n';
        ++failures;
    }
}

// Worst observed relative error per labelled comparison, printed as a table at
// the end. The point is not decoration: a reviewer reading "all nine
// components validated" needs to see how close each one actually came, and a
// check whose error is exactly zero for every fixture usually means the
// comparison is not doing what it claims.
struct ErrorRecord {
    std::string label;
    double worst_error = 0.0;
    std::size_t worst_component = 0;
    double tolerance = 0.0;
};
std::vector<ErrorRecord> error_records;

void record_error(const std::string& label, double error, std::size_t component,
                  double tolerance) {
    for (auto& record : error_records) {
        if (record.label == label) {
            if (error > record.worst_error) {
                record.worst_error = error;
                record.worst_component = component;
            }
            return;
        }
    }
    error_records.push_back({label, error, component, tolerance});
}

void print_error_summary() {
    std::cout << "\nmaximum relative error per validated tensor "
              << "(relative to max |W| of the reference)\n";
    std::cout << "  worst   component  tolerance   check\n";
    static const char* names[9] = {"xx", "xy", "xz", "yx", "yy", "yz",
                                   "zx", "zy", "zz"};
    for (const auto& record : error_records) {
        std::ostringstream line;
        line << "  " << std::scientific << std::setprecision(2) << record.worst_error
             << "   W_" << names[record.worst_component] << "       "
             << std::scientific << std::setprecision(0) << record.tolerance
             << "      " << record.label;
        std::cout << line.str() << '\n';
    }
}

// --- evaluation helpers ---------------------------------------------------

gmd::ForceResult evaluate(gmd::ForceProvider& provider, gmd::System& system) {
    gmd::RuntimeContext runtime;
    gmd::ForceResult result;
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
    provider.compute(request, result, runtime);
    return result;
}

const char* component_name(std::size_t a, std::size_t b) {
    static const char* names[9] = {"xx", "xy", "xz", "yx", "yy", "yz",
                                   "zx", "zy", "zz"};
    return names[a * 3 + b];
}

// Compares all nine components against a reference tensor.
//
// `minimum_relative_magnitude` is the guard that stops a blank or collapsed
// tensor from passing: every component of the REFERENCE must be at least this
// fraction of the largest one. A fixture that cannot meet it is the wrong
// fixture for a nine-component claim.
void check_tensor(const std::array<double, 9>& actual,
                  const vr::Mat3& expected,
                  const std::string& label,
                  double tolerance,
                  double minimum_relative_magnitude = 1.0e-3) {
    double scale = 0.0;
    for (const auto value : expected) scale = std::max(scale, std::fabs(static_cast<double>(value)));
    check(scale > 0.0, label + ": reference tensor is entirely zero");
    if (scale <= 0.0) return;

    for (std::size_t a = 0; a < 3; ++a) {
        for (std::size_t b = 0; b < 3; ++b) {
            const double reference = static_cast<double>(expected[a * 3 + b]);
            const double measured = actual[a * 3 + b];
            const double error = std::fabs(measured - reference) / scale;
            record_error(label, error, a * 3 + b, tolerance);

            check(error < tolerance,
                  label + ": W_" + component_name(a, b) + " = " +
                      std::to_string(measured) + ", expected " +
                      std::to_string(reference) + " (relative error " +
                      std::to_string(error) + ")");

            check(std::fabs(reference) >= minimum_relative_magnitude * scale,
                  label + ": W_" + component_name(a, b) +
                      " is too small for this fixture to test it (|" +
                      std::to_string(reference) + "| vs scale " +
                      std::to_string(scale) + "); a nine-component claim needs "
                      "every component to carry signal");

            // Sign is stated separately: a component that is right in
            // magnitude and wrong in sign is a different defect from one that
            // is merely inaccurate, and the relative check above reports it
            // only through the size of the difference.
            if (std::fabs(reference) > 1.0e-3 * scale) {
                check((measured > 0.0) == (reference > 0.0),
                      label + ": W_" + component_name(a, b) +
                          " has the wrong sign (" + std::to_string(measured) +
                          " vs " + std::to_string(reference) + ")");
            }
        }
    }
}

double max_component_difference(const std::array<double, 9>& lhs,
                                const std::array<double, 9>& rhs) {
    double worst = 0.0;
    for (std::size_t i = 0; i < 9; ++i) worst = std::max(worst, std::fabs(lhs[i] - rhs[i]));
    return worst;
}

double tensor_scale(const std::array<double, 9>& tensor) {
    double scale = 0.0;
    for (const auto value : tensor) scale = std::max(scale, std::fabs(value));
    return scale;
}

std::array<double, 9> subtract(const std::array<double, 9>& lhs,
                               const std::array<double, 9>& rhs) {
    std::array<double, 9> out{};
    for (std::size_t i = 0; i < 9; ++i) out[i] = lhs[i] - rhs[i];
    return out;
}

vr::Mat3 to_reference(const std::array<double, 9>& tensor) {
    vr::Mat3 out{};
    for (std::size_t i = 0; i < 9; ++i) out[i] = static_cast<vr::Real>(tensor[i]);
    return out;
}

// --- system construction --------------------------------------------------

gmd::System make_system(const std::array<double, 3>& lengths,
                        const std::vector<vr::Vec3>& coordinates,
                        const std::vector<double>& charges = {}) {
    gmd::System system;
    system.resize(coordinates.size(), coordinates.size());
    gmd::Box box;
    box.set_lengths(lengths);
    system.set_box(box);
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {static_cast<double>(coordinates[i][0]),
                                           static_cast<double>(coordinates[i][1]),
                                           static_cast<double>(coordinates[i][2])};
        if (i < charges.size()) system.mutable_charges()[i] = charges[i];
    }
    return system;
}

std::vector<vr::Vec3> coordinates_of(const gmd::System& system) {
    std::vector<vr::Vec3> out;
    const auto coordinates = system.coordinates();
    out.reserve(coordinates.size());
    for (const auto& coordinate : coordinates) {
        out.push_back({static_cast<vr::Real>(coordinate[0]),
                       static_cast<vr::Real>(coordinate[1]),
                       static_cast<vr::Real>(coordinate[2])});
    }
    return out;
}

// Minimum image against an orthorhombic cell, written here rather than taken
// from the engine so a defect in the engine's wrapping cannot cancel itself.
vr::Vec3 minimum_image(vr::Vec3 delta, const std::array<double, 3>& lengths) {
    for (std::size_t d = 0; d < 3; ++d) {
        const auto length = static_cast<vr::Real>(lengths[d]);
        while (delta[d] > 0.5L * length) delta[d] -= length;
        while (delta[d] < -0.5L * length) delta[d] += length;
    }
    return delta;
}

// ===========================================================================
// A. Pair interactions: Lennard-Jones
// ===========================================================================

struct LennardJonesPair {
    vr::Real energy;
    vr::Real force_factor;  // F_i = force_factor * dr
};

// V(r) = 4 eps [(sig/r)^12 - (sig/r)^6], grouped differently from the engine's
// eps4 * (12 s12 - 6 s6) / r^2 so that a transcription slip cannot be shared.
LennardJonesPair lennard_jones(vr::Real r_squared, vr::Real epsilon, vr::Real sigma) {
    const vr::Real ratio_squared = sigma * sigma / r_squared;
    const vr::Real ratio_six = ratio_squared * ratio_squared * ratio_squared;
    const vr::Real ratio_twelve = ratio_six * ratio_six;
    return {4.0L * epsilon * (ratio_twelve - ratio_six),
            24.0L * epsilon * (2.0L * ratio_twelve - ratio_six) / r_squared};
}

// F = -dV/dr, checked against a central difference of the same V. This is the
// step that makes the analytical pair reference trustworthy on its own terms.
void check_lennard_jones_force_against_derivative() {
    const vr::Real epsilon = 0.0104L;
    const vr::Real sigma = 3.4L;
    for (vr::Real r : {2.9L, 3.4L, 3.8L, 4.5L, 6.0L}) {
        const vr::Real step = 1.0e-6L;
        const vr::Real up = lennard_jones((r + step) * (r + step), epsilon, sigma).energy;
        const vr::Real down = lennard_jones((r - step) * (r - step), epsilon, sigma).energy;
        const vr::Real numerical = -(up - down) / (2.0L * step);
        const vr::Real analytic = lennard_jones(r * r, epsilon, sigma).force_factor * r;
        const auto error = static_cast<double>(std::fabs(numerical - analytic)
                                               / std::max(std::fabs(analytic), 1.0e-12L));
        check(error < 1.0e-6,
              "[analytic-pair] LJ reference force disagrees with -dV/dr at r = " +
                  std::to_string(static_cast<double>(r)) + " (relative error " +
                  std::to_string(error) + ")");
    }
}

// The exact pair tensor: W_ab = r_a F_b summed over evaluated pairs.
vr::Mat3 lennard_jones_reference(const gmd::System& system,
                                 vr::Real epsilon,
                                 vr::Real sigma,
                                 vr::Real cutoff,
                                 const std::function<vr::Real(std::size_t, std::size_t)>& scale_of,
                                 const std::function<void(std::size_t, std::size_t,
                                                          vr::Real&, vr::Real&)>& params_of = nullptr) {
    const auto coordinates = coordinates_of(system);
    const auto lengths = system.box().lengths;
    vr::Mat3 virial{};
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        for (std::size_t j = i + 1; j < coordinates.size(); ++j) {
            const vr::Vec3 delta =
                minimum_image(vr::sub(coordinates[i], coordinates[j]), lengths);
            const vr::Real r_squared = vr::dot(delta, delta);
            if (r_squared >= cutoff * cutoff) continue;
            const vr::Real scale = scale_of ? scale_of(i, j) : 1.0L;
            if (scale == 0.0L) continue;

            vr::Real pair_epsilon = epsilon;
            vr::Real pair_sigma = sigma;
            if (params_of) params_of(i, j, pair_epsilon, pair_sigma);

            const vr::Real factor =
                scale * lennard_jones(r_squared, pair_epsilon, pair_sigma).force_factor;
            const vr::Vec3 force = {factor * delta[0], factor * delta[1], factor * delta[2]};
            vr::accumulate_outer(virial, delta, force);
        }
    }
    return virial;
}

// A pair whose separation has all three components non-zero and comfortably
// unequal, so no component of the tensor can hide behind a symmetry.
std::vector<vr::Vec3> tilted_pair(vr::Real ox = 0.0L, vr::Real oy = 0.0L, vr::Real oz = 0.0L) {
    return {{8.0L + ox, 8.0L + oy, 8.0L + oz},
            {10.13L + ox, 9.47L + oy, 10.71L + oz}};
}

gmd::ClassicalForceProvider make_single_type_lj(double epsilon, double sigma, double cutoff) {
    return gmd::ClassicalForceProvider(epsilon, sigma, cutoff);
}

void test_lennard_jones_pair_tensor() {
    const double epsilon = 0.0104;
    const double sigma = 3.4;
    const double cutoff = 8.5;

    // Attractive branch: r > 2^(1/6) sigma, so W_xx must be negative.
    // Repulsive branch: r < 2^(1/6) sigma, so W_xx must be positive. Running
    // both pins the sign convention down in a way a single distance cannot.
    struct Case {
        std::string name;
        std::vector<vr::Vec3> coordinates;
        bool expect_attractive;
    };
    // r_min = 2^(1/6) sigma = 3.816 for sigma = 3.4. The two separations below
    // are 4.60 and 1.83, which straddle it comfortably; the generic
    // tilted_pair() fixture sits at 3.75 and is therefore (mildly) repulsive,
    // which is why it is not reused here.
    const std::vector<Case> cases = {
        {"attractive", {{8.0L, 8.0L, 8.0L}, {10.615L, 9.804L, 11.327L}}, true},
        {"repulsive", {{8.0L, 8.0L, 8.0L}, {9.05L, 8.72L, 9.31L}}, false},
    };

    for (const auto& test_case : cases) {
        gmd::System system = make_system({20.0, 20.0, 20.0}, test_case.coordinates);
        auto provider = make_single_type_lj(epsilon, sigma, cutoff);
        const gmd::ForceResult result = evaluate(provider, system);
        check(result.virial_valid, "[analytic-pair] LJ " + test_case.name +
                                       ": virial must be reported valid");

        const vr::Mat3 reference = lennard_jones_reference(
            system, epsilon, sigma, cutoff, nullptr);
        check_tensor(result.virial, reference,
                     "[analytic-pair] LJ " + test_case.name, 1.0e-12);

        const bool attractive = result.virial[0] < 0.0;
        check(attractive == test_case.expect_attractive,
              "[analytic-pair] LJ " + test_case.name +
                  ": W_xx sign does not match the expected branch of the potential");
    }
}

void test_lennard_jones_mixed_types_and_override() {
    const double cutoff = 8.5;

    gmd::LJForceFieldConfig config;
    config.cutoff = cutoff;
    config.mixing_rule = "lorentz_berthelot";
    config.elements = {
        {"Ar", 39.948, 0.0104, 3.40, 0.0},
        {"Ne", 20.180, 0.0031, 2.75, 0.0},
    };

    // Three atoms: types 0, 1, 0. Pairs (0,1) and (1,2) are cross-type and
    // (0,2) is same-type, so one evaluation exercises both table paths.
    const std::vector<vr::Vec3> coordinates = {
        {8.0L, 8.0L, 8.0L}, {10.31L, 9.44L, 10.62L}, {6.13L, 10.77L, 9.35L}};

    auto run = [&](const gmd::LJForceFieldConfig& active, const std::string& label) {
        gmd::System system = make_system({20.0, 20.0, 20.0}, coordinates);
        system.mutable_atom_types()[0] = 0;
        system.mutable_atom_types()[1] = 1;
        system.mutable_atom_types()[2] = 0;

        gmd::ClassicalForceProvider provider(active);
        const gmd::ForceResult result = evaluate(provider, system);

        const std::vector<int> types = {0, 1, 0};
        const vr::Mat3 reference = lennard_jones_reference(
            system, 0.0L, 0.0L, static_cast<vr::Real>(cutoff), nullptr,
            [&](std::size_t i, std::size_t j, vr::Real& pair_epsilon, vr::Real& pair_sigma) {
                double epsilon = 0.0;
                double sigma = 0.0;
                active.pair_params(types[i], types[j], epsilon, sigma);
                pair_epsilon = static_cast<vr::Real>(epsilon);
                pair_sigma = static_cast<vr::Real>(sigma);
            });
        check_tensor(result.virial, reference, label, 1.0e-12);
        return result.virial;
    };

    const auto mixed = run(config, "[analytic-pair] LJ Lorentz-Berthelot mixing");

    // An explicit override for the cross pair must change the tensor, and must
    // change it to the value the override implies. Asserting only the second
    // half would pass if the override were silently ignored.
    gmd::LJForceFieldConfig overridden = config;
    overridden.pair_overrides[{0, 1}] = gmd::ExplicitPairOverride{0.0210, 4.05};
    const auto override_result = run(overridden, "[analytic-pair] LJ explicit pair override");

    check(max_component_difference(mixed, override_result) >
              1.0e-6 * tensor_scale(mixed),
          "[analytic-pair] LJ explicit pair override: the override did not change "
          "the tensor, so this fixture proves nothing about it");
}

void test_lennard_jones_special_pair_scaling() {
    const double epsilon = 0.0104;
    const double sigma = 3.4;
    const double cutoff = 8.5;

    // A four-atom chain: bonds make (0,1) (1,2) (2,3) 1-2 pairs, angles make
    // (0,2) (1,3) 1-3 pairs, and the dihedral makes (0,3) the 1-4 pair. With
    // 1-2 and 1-3 fully excluded, the entire LJ virial is the 1-4 pair, so the
    // scale factor is read directly off the tensor.
    auto topology = std::make_shared<gmd::Topology>();
    topology->bonds = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}};
    topology->angles = {{0, 1, 2, 0}, {1, 2, 3, 0}};
    topology->dihedrals = {{0, 1, 2, 3, 0}};

    const std::vector<vr::Vec3> coordinates = {
        {7.00L, 7.00L, 7.00L}, {8.42L, 7.63L, 7.91L},
        {9.11L, 9.02L, 8.34L}, {10.37L, 9.31L, 9.72L}};

    vr::Mat3 unscaled_reference{};
    std::array<double, 9> unscaled_measured{};

    for (const double lj_scale : {0.0, 0.5, 1.0}) {
        gmd::SpecialPairScaleConfig scales;
        scales.pair_12 = {0.0, 0.0};
        scales.pair_13 = {0.0, 0.0};
        scales.pair_14 = {lj_scale, 0.0};
        auto special_pairs = std::make_shared<gmd::SpecialPairMap>(*topology, scales);

        gmd::System system = make_system({24.0, 24.0, 24.0}, coordinates);
        system.set_special_pair_map(special_pairs);

        auto provider = make_single_type_lj(epsilon, sigma, cutoff);
        const gmd::ForceResult result = evaluate(provider, system);

        const vr::Mat3 reference = lennard_jones_reference(
            system, epsilon, sigma, cutoff,
            [&](std::size_t i, std::size_t j) {
                return static_cast<vr::Real>(
                    special_pairs->scale_for(static_cast<int>(i), static_cast<int>(j)).lj);
            });

        const std::string label =
            "[analytic-pair] LJ 1-4 scale " + std::to_string(lj_scale);

        if (lj_scale == 0.0) {
            // Full exclusion: the tensor must be exactly zero, not merely
            // small. check_tensor() cannot express that, so state it directly.
            check(tensor_scale(result.virial) == 0.0,
                  label + ": an excluded 1-4 pair must contribute exactly zero, got scale " +
                      std::to_string(tensor_scale(result.virial)));
            continue;
        }

        check_tensor(result.virial, reference, label, 1.0e-12);

        if (lj_scale == 1.0) {
            unscaled_reference = reference;
            unscaled_measured = result.virial;
        }
    }

    // Linearity in the scale factor: the half-scaled tensor must be exactly
    // half the unscaled one, component by component.
    gmd::SpecialPairScaleConfig half_scales;
    half_scales.pair_12 = {0.0, 0.0};
    half_scales.pair_13 = {0.0, 0.0};
    half_scales.pair_14 = {0.5, 0.0};
    auto half_map = std::make_shared<gmd::SpecialPairMap>(*topology, half_scales);
    gmd::System half_system = make_system({24.0, 24.0, 24.0}, coordinates);
    half_system.set_special_pair_map(half_map);
    auto provider = make_single_type_lj(epsilon, sigma, cutoff);
    const auto half = evaluate(provider, half_system).virial;

    double worst = 0.0;
    for (std::size_t i = 0; i < 9; ++i) {
        worst = std::max(worst, std::fabs(half[i] - 0.5 * unscaled_measured[i]));
    }
    check(worst < 1.0e-15 * tensor_scale(unscaled_measured),
          "[analytic-pair] LJ 1-4 scaling is not linear in the scale factor (worst "
          "component deviation " + std::to_string(worst) + ")");
    (void)unscaled_reference;
}

// Every face and the corner. A pair that straddles a boundary must produce the
// same tensor as the same pair sitting in the middle of the cell.
void test_pair_virial_across_periodic_boundaries() {
    const double epsilon = 0.0104;
    const double sigma = 3.4;
    const double cutoff = 8.5;
    const std::array<double, 3> lengths = {20.0, 22.0, 18.0};

    gmd::System interior = make_system(lengths, tilted_pair());
    auto provider = make_single_type_lj(epsilon, sigma, cutoff);
    const auto expected = evaluate(provider, interior).virial;
    check(tensor_scale(expected) > 0.0,
          "[analytic-pair] PBC baseline pair tensor is zero");

    struct Placement {
        std::string name;
        vr::Vec3 first;
        vr::Vec3 second;
    };
    // Each placement puts the first atom just inside one face (or corner) and
    // the second just outside, so the minimum image is the only thing that can
    // reconstruct the pair.
    const std::vector<Placement> placements = {
        {"x face", {19.40L, 8.00L, 8.00L}, {1.53L, 9.47L, 10.71L}},
        {"y face", {8.00L, 21.50L, 8.00L}, {10.13L, 0.97L, 10.71L}},
        {"z face", {8.00L, 8.00L, 17.30L}, {10.13L, 9.47L, 2.01L}},
        {"corner", {19.40L, 21.50L, 17.30L}, {1.53L, 0.97L, 2.01L}},
    };

    for (const auto& placement : placements) {
        gmd::System system = make_system(lengths, {placement.first, placement.second});
        auto boundary_provider = make_single_type_lj(epsilon, sigma, cutoff);
        const auto measured = evaluate(boundary_provider, system).virial;

        const vr::Mat3 reference = lennard_jones_reference(
            system, epsilon, sigma, cutoff, nullptr);
        check_tensor(measured, reference,
                     "[analytic-pair] LJ across the " + placement.name, 1.0e-12);

        // And it must reproduce the interior placement, which is the actual
        // physical claim: wrapping is not allowed to change the tensor.
        const double difference = max_component_difference(measured, expected);
        check(difference < 1.0e-12 * tensor_scale(expected),
              "[symmetry] LJ tensor across the " + placement.name +
                  " differs from the interior placement by " +
                  std::to_string(difference));
    }
}

void test_pair_virial_is_translation_invariant() {
    const double epsilon = 0.0104;
    const double sigma = 3.4;
    const double cutoff = 8.5;
    const std::array<double, 3> lengths = {20.0, 22.0, 18.0};

    gmd::System base = make_system(lengths, tilted_pair());
    auto provider = make_single_type_lj(epsilon, sigma, cutoff);
    const auto expected = evaluate(provider, base).virial;

    // A translation that is not a lattice vector, followed by wrapping into
    // the cell. Both atoms move; the separation does not.
    const std::array<double, 3> shift = {13.7, 17.3, 9.1};
    std::vector<vr::Vec3> shifted;
    for (const auto& coordinate : tilted_pair()) {
        vr::Vec3 moved{};
        for (std::size_t d = 0; d < 3; ++d) {
            vr::Real value = coordinate[d] + static_cast<vr::Real>(shift[d]);
            const auto length = static_cast<vr::Real>(lengths[d]);
            while (value >= length) value -= length;
            while (value < 0.0L) value += length;
            moved[d] = value;
        }
        shifted.push_back(moved);
    }

    gmd::System translated = make_system(lengths, shifted);
    auto translated_provider = make_single_type_lj(epsilon, sigma, cutoff);
    const auto measured = evaluate(translated_provider, translated).virial;

    const double difference = max_component_difference(measured, expected);
    check(difference < 1.0e-12 * tensor_scale(expected),
          "[symmetry] LJ tensor is not invariant under a whole-box translation "
          "and re-wrap (worst component difference " + std::to_string(difference) + ")");
}

// Rotating a compact cluster inside a cubic box, far from every face, leaves
// the periodic lattice irrelevant, so the tensor must obey W -> R W R^T.
void test_pair_virial_rotation_covariance() {
    const double epsilon = 0.0104;
    const double sigma = 3.4;
    const double cutoff = 8.5;
    const std::array<double, 3> lengths = {40.0, 40.0, 40.0};

    const std::vector<vr::Vec3> cluster = {
        {20.00L, 20.00L, 20.00L}, {23.11L, 20.42L, 21.37L},
        {20.63L, 23.24L, 21.06L}, {21.42L, 21.07L, 23.31L}};

    gmd::System base = make_system(lengths, cluster);
    auto provider = make_single_type_lj(epsilon, sigma, cutoff);
    const auto measured = evaluate(provider, base).virial;
    check(tensor_scale(measured) > 0.0,
          "[rotation] LJ rotation fixture produced a zero tensor");

    const vr::Mat3 rotation =
        vr::rotation_matrix({0.37L, -0.62L, 0.79L}, 0.9137L);

    std::vector<vr::Vec3> rotated;
    const vr::Vec3 centre = {20.0L, 20.0L, 20.0L};
    for (const auto& coordinate : cluster) {
        const vr::Vec3 relative = vr::sub(coordinate, centre);
        const vr::Vec3 turned = vr::mat_vec(rotation, relative);
        rotated.push_back({turned[0] + centre[0], turned[1] + centre[1],
                           turned[2] + centre[2]});
    }

    gmd::System rotated_system = make_system(lengths, rotated);
    auto rotated_provider = make_single_type_lj(epsilon, sigma, cutoff);
    const auto rotated_measured = evaluate(rotated_provider, rotated_system).virial;

    const vr::Mat3 expected = vr::rotate_tensor(to_reference(measured), rotation);
    check_tensor(rotated_measured, expected,
                 "[rotation] LJ cluster covariance W -> R W R^T", 1.0e-10);
}

}  // namespace

// ===========================================================================
// B. Bonded interactions
// ===========================================================================

namespace {

// Independent geometry. The torsion is written with the IUPAC form
//   phi = atan2(|b2| b1.(b2 x b3), (b1 x b2).(b2 x b3))
// which is algebraically the same convention the engine reaches through
// m = b1 x b2, n = b2 x b3 and sin phi = (m x n).b2 / (|m||n||b2|), but arrives
// at it by a different route.
vr::Vec3 cross(const vr::Vec3& a, const vr::Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

vr::Real torsion_angle(const vr::Vec3& ri, const vr::Vec3& rj,
                       const vr::Vec3& rk, const vr::Vec3& rl) {
    const vr::Vec3 b1 = vr::sub(rj, ri);
    const vr::Vec3 b2 = vr::sub(rk, rj);
    const vr::Vec3 b3 = vr::sub(rl, rk);
    const vr::Vec3 m = cross(b1, b2);
    const vr::Vec3 n = cross(b2, b3);
    return std::atan2(vr::norm(b2) * vr::dot(b1, n), vr::dot(m, n));
}

vr::Real bend_angle(const vr::Vec3& ri, const vr::Vec3& rj, const vr::Vec3& rk) {
    const vr::Vec3 a = vr::sub(ri, rj);
    const vr::Vec3 b = vr::sub(rk, rj);
    vr::Real cosine = vr::dot(a, b) / (vr::norm(a) * vr::norm(b));
    if (cosine > 1.0L) cosine = 1.0L;
    if (cosine < -1.0L) cosine = -1.0L;
    return std::acos(cosine);
}

// Which bonded term a fixture carries, and its independently written energy.
enum class BondedTerm { Bond, Angle, Dihedral, Improper };

struct BondedFixture {
    std::string name;
    BondedTerm term;
    std::vector<vr::Vec3> coordinates;
    std::vector<int> atoms;   // the atoms the single interaction connects
    // parameters
    double k = 0.0;
    double r0 = 0.0;
    double theta0 = 0.0;
    int periodicity = 1;
    double phase = 0.0;
    double phi0 = 0.0;
};

vr::Real bonded_energy(const BondedFixture& fixture,
                       const std::vector<vr::Vec3>& coordinates) {
    const auto& atoms = fixture.atoms;
    switch (fixture.term) {
        case BondedTerm::Bond: {
            const vr::Real r = vr::norm(vr::sub(coordinates[atoms[1]], coordinates[atoms[0]]));
            const vr::Real deviation = r - static_cast<vr::Real>(fixture.r0);
            return static_cast<vr::Real>(fixture.k) * deviation * deviation;
        }
        case BondedTerm::Angle: {
            const vr::Real theta = bend_angle(coordinates[atoms[0]], coordinates[atoms[1]],
                                              coordinates[atoms[2]]);
            const vr::Real deviation = theta - static_cast<vr::Real>(fixture.theta0);
            return static_cast<vr::Real>(fixture.k) * deviation * deviation;
        }
        case BondedTerm::Dihedral: {
            const vr::Real phi = torsion_angle(coordinates[atoms[0]], coordinates[atoms[1]],
                                               coordinates[atoms[2]], coordinates[atoms[3]]);
            return static_cast<vr::Real>(fixture.k)
                 * (1.0L + std::cos(static_cast<vr::Real>(fixture.periodicity) * phi
                                    - static_cast<vr::Real>(fixture.phase)));
        }
        case BondedTerm::Improper: {
            const vr::Real phi = torsion_angle(coordinates[atoms[0]], coordinates[atoms[1]],
                                               coordinates[atoms[2]], coordinates[atoms[3]]);
            const vr::Real deviation = phi - static_cast<vr::Real>(fixture.phi0);
            return static_cast<vr::Real>(fixture.k) * deviation * deviation;
        }
    }
    return 0.0L;
}

std::shared_ptr<gmd::Topology> topology_for(const BondedFixture& fixture) {
    auto topology = std::make_shared<gmd::Topology>();
    const auto& a = fixture.atoms;
    switch (fixture.term) {
        case BondedTerm::Bond: topology->bonds = {{a[0], a[1], 0}}; break;
        case BondedTerm::Angle: topology->angles = {{a[0], a[1], a[2], 0}}; break;
        case BondedTerm::Dihedral: topology->dihedrals = {{a[0], a[1], a[2], a[3], 0}}; break;
        case BondedTerm::Improper: topology->impropers = {{a[0], a[1], a[2], a[3], 0}}; break;
    }
    return topology;
}

std::shared_ptr<gmd::BondedForceProvider> provider_for(const BondedFixture& fixture) {
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology_for(fixture));
    switch (fixture.term) {
        case BondedTerm::Bond:
            provider->add_bond_type({fixture.k, fixture.r0});
            break;
        case BondedTerm::Angle:
            provider->add_angle_type({fixture.k, fixture.theta0});
            break;
        case BondedTerm::Dihedral:
            provider->add_dihedral_type({fixture.k, fixture.periodicity, fixture.phase});
            break;
        case BondedTerm::Improper:
            provider->add_improper_type({fixture.k, fixture.phi0});
            break;
    }
    return provider;
}

// Non-planar, non-axis-aligned fixtures. A planar torsion or an angle lying in
// a coordinate plane silently zeroes tensor components, which would make a
// nine-component claim vacuous.
std::vector<BondedFixture> bonded_fixtures() {
    const std::vector<vr::Vec3> chain = {
        {9.00L, 9.00L, 9.00L},
        {10.43L, 9.61L, 9.92L},
        {11.17L, 11.02L, 9.38L},
        {12.61L, 11.34L, 10.46L}};

    std::vector<BondedFixture> fixtures;
    fixtures.push_back({"harmonic bond", BondedTerm::Bond, chain, {0, 1},
                        /*k=*/3.1, /*r0=*/1.35});
    fixtures.push_back({"harmonic angle", BondedTerm::Angle, chain, {0, 1, 2},
                        /*k=*/2.4, 0.0, /*theta0=*/1.9106});
    // A non-zero phase makes the torsion sign-sensitive: with delta = 0 the
    // energy is even in phi and a sign error in the angle convention cancels.
    BondedFixture dihedral{"proper dihedral", BondedTerm::Dihedral, chain, {0, 1, 2, 3},
                           /*k=*/0.62};
    dihedral.periodicity = 3;
    dihedral.phase = 0.7853981633974483;  // pi/4
    fixtures.push_back(dihedral);

    BondedFixture improper{"harmonic improper", BondedTerm::Improper, chain, {0, 1, 2, 3},
                           /*k=*/1.9};
    improper.phi0 = -0.35;  // non-zero, and of the opposite sign to phi
    fixtures.push_back(improper);
    return fixtures;
}

// Step 1 of the bonded reference: the provider's per-atom forces must be
// -dV/dr of the independently written energy. Everything downstream depends on
// this, so it is checked on its own and reported on its own.
std::vector<vr::Vec3> verify_bonded_forces(const BondedFixture& fixture,
                                           const gmd::ForceResult& result,
                                           const std::vector<vr::Vec3>& coordinates) {
    std::vector<vr::Vec3> forces(coordinates.size(), vr::Vec3{0.0L, 0.0L, 0.0L});
    const vr::Real step = 2.0e-6L;

    double worst = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            std::vector<vr::Vec3> up = coordinates;
            std::vector<vr::Vec3> down = coordinates;
            up[i][d] += step;
            down[i][d] -= step;
            const vr::Real derivative =
                (bonded_energy(fixture, up) - bonded_energy(fixture, down)) / (2.0L * step);
            forces[i][d] = -derivative;
            worst = std::max(worst, std::fabs(static_cast<double>(forces[i][d])
                                              - result.forces[i][d]));
            scale = std::max(scale, std::fabs(static_cast<double>(forces[i][d])));
        }
    }

    check(scale > 1.0e-6,
          "[bonded-moment] " + fixture.name +
              ": fixture produces no force, so it cannot validate a virial");
    check(worst < 1.0e-6 * std::max(scale, 1.0e-12),
          "[bonded-moment] " + fixture.name +
              ": provider forces disagree with -dV/dr of the independent energy "
              "(worst component " + std::to_string(worst) + ", scale " +
              std::to_string(scale) + ")");

    // Return the PROVIDER's forces for the moment sum. They have just been
    // certified against an independent derivative, and using them avoids
    // feeding finite-difference noise into the tensor comparison.
    std::vector<vr::Vec3> certified(coordinates.size());
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        certified[i] = {static_cast<vr::Real>(result.forces[i][0]),
                        static_cast<vr::Real>(result.forces[i][1]),
                        static_cast<vr::Real>(result.forces[i][2])};
    }
    return certified;
}

// Step 2: W_ab = sum_i (r_i - r_ref)_a F_i,b over the interaction's own atoms,
// with the relative coordinates built by walking the chain under the minimum
// image. accum_interaction_virial() is not used, and neither is any engine
// helper: the walk is written out here.
vr::Mat3 bonded_moment_reference(const BondedFixture& fixture,
                                 const std::vector<vr::Vec3>& coordinates,
                                 const std::vector<vr::Vec3>& forces,
                                 const std::array<double, 3>& lengths,
                                 std::size_t reference_atom) {
    const auto& atoms = fixture.atoms;
    std::vector<vr::Vec3> relative(atoms.size());
    relative[reference_atom] = {0.0L, 0.0L, 0.0L};

    // Walk outward from the reference atom, accumulating minimum-image steps,
    // so the unwrapped interaction-local geometry is reconstructed even when
    // the molecule straddles a face.
    for (std::size_t offset = 1; offset < atoms.size(); ++offset) {
        if (reference_atom + offset < atoms.size()) {
            const std::size_t current = reference_atom + offset;
            const vr::Vec3 step = minimum_image(
                vr::sub(coordinates[atoms[current]], coordinates[atoms[current - 1]]), lengths);
            relative[current] = {relative[current - 1][0] + step[0],
                                 relative[current - 1][1] + step[1],
                                 relative[current - 1][2] + step[2]};
        }
    }
    for (std::size_t index = reference_atom; index-- > 0;) {
        const vr::Vec3 step = minimum_image(
            vr::sub(coordinates[atoms[index]], coordinates[atoms[index + 1]]), lengths);
        relative[index] = {relative[index + 1][0] + step[0],
                           relative[index + 1][1] + step[1],
                           relative[index + 1][2] + step[2]};
    }

    vr::Mat3 virial{};
    for (std::size_t a = 0; a < atoms.size(); ++a) {
        vr::accumulate_outer(virial, relative[a], forces[static_cast<std::size_t>(atoms[a])]);
    }
    return virial;
}

void test_bonded_term_tensors() {
    const std::array<double, 3> lengths = {24.0, 26.0, 22.0};

    for (const auto& fixture : bonded_fixtures()) {
        gmd::System system = make_system(lengths, fixture.coordinates);
        auto provider = provider_for(fixture);
        gmd::RuntimeContext runtime;
        provider->initialize(runtime);

        const gmd::ForceResult result = evaluate(*provider, system);
        check(result.virial_valid,
              "[bonded-moment] " + fixture.name + ": virial must be reported valid");

        const auto coordinates = coordinates_of(system);

        // The independently written energy must agree with the provider's, or
        // the two are not describing the same interaction and nothing below
        // means anything.
        const auto energy_error =
            static_cast<double>(std::fabs(bonded_energy(fixture, coordinates)
                                          - static_cast<vr::Real>(result.potential_energy)));
        check(energy_error < 1.0e-12 * std::max(std::fabs(result.potential_energy), 1.0e-12),
              "[bonded-moment] " + fixture.name +
                  ": independent energy disagrees with the provider (" +
                  std::to_string(energy_error) + ")");

        const auto forces = verify_bonded_forces(fixture, result, coordinates);
        const vr::Mat3 reference =
            bonded_moment_reference(fixture, coordinates, forces, lengths, 0);
        check_tensor(result.virial, reference,
                     "[bonded-moment] " + fixture.name, 1.0e-11);

        // The moment sum is independent of which atom is called the origin,
        // because the interaction's forces sum to zero. If it is not, the
        // forces do not sum to zero and the tensor is origin-dependent.
        for (std::size_t origin = 1; origin < fixture.atoms.size(); ++origin) {
            const vr::Mat3 shifted =
                bonded_moment_reference(fixture, coordinates, forces, lengths, origin);
            double worst = 0.0;
            double scale = 0.0;
            for (std::size_t i = 0; i < 9; ++i) {
                worst = std::max(worst, std::fabs(static_cast<double>(shifted[i] - reference[i])));
                scale = std::max(scale, std::fabs(static_cast<double>(reference[i])));
            }
            check(worst < 1.0e-10 * scale,
                  "[bonded-moment] " + fixture.name +
                      ": moment sum depends on the reference atom (origin " +
                      std::to_string(origin) + ", worst " + std::to_string(worst) + ")");
        }
    }
}

// The same molecule, once intact in the middle of the cell and once translated
// so that it straddles a face and its atoms wrap. The tensor must not move.
void test_bonded_wrapped_equals_intact() {
    const std::array<double, 3> lengths = {24.0, 26.0, 22.0};

    for (const auto& fixture : bonded_fixtures()) {
        gmd::System intact = make_system(lengths, fixture.coordinates);
        auto provider = provider_for(fixture);
        gmd::RuntimeContext runtime;
        provider->initialize(runtime);
        const auto expected = evaluate(*provider, intact).virial;
        check(tensor_scale(expected) > 0.0,
              "[symmetry] " + fixture.name + ": intact tensor is zero");

        // Shift so the chain crosses the x and z faces, then wrap every atom.
        std::vector<vr::Vec3> wrapped;
        const std::array<double, 3> shift = {14.6, 0.0, 12.9};
        for (const auto& coordinate : fixture.coordinates) {
            vr::Vec3 moved{};
            for (std::size_t d = 0; d < 3; ++d) {
                vr::Real value = coordinate[d] + static_cast<vr::Real>(shift[d]);
                const auto length = static_cast<vr::Real>(lengths[d]);
                while (value >= length) value -= length;
                while (value < 0.0L) value += length;
                moved[d] = value;
            }
            wrapped.push_back(moved);
        }

        gmd::System crossing = make_system(lengths, wrapped);
        auto crossing_provider = provider_for(fixture);
        crossing_provider->initialize(runtime);
        const auto measured = evaluate(*crossing_provider, crossing).virial;

        const double difference = max_component_difference(measured, expected);
        check(difference < 1.0e-11 * tensor_scale(expected),
              "[symmetry] " + fixture.name +
                  ": wrapped molecule gives a different tensor from the intact one "
                  "(worst component difference " + std::to_string(difference) + ")");
    }
}

void test_bonded_rotation_covariance() {
    const std::array<double, 3> lengths = {40.0, 40.0, 40.0};
    const vr::Mat3 rotation = vr::rotation_matrix({-0.51L, 0.44L, 0.74L}, 1.2731L);
    const vr::Vec3 centre = {20.0L, 20.0L, 20.0L};

    for (const auto& fixture : bonded_fixtures()) {
        // Re-centre the fixture in the larger cubic box so that rotation keeps
        // it well away from every face and the lattice stays irrelevant.
        std::vector<vr::Vec3> base;
        for (const auto& coordinate : fixture.coordinates) {
            base.push_back({coordinate[0] + 10.0L, coordinate[1] + 9.0L,
                            coordinate[2] + 10.0L});
        }

        gmd::System system = make_system(lengths, base);
        auto provider = provider_for(fixture);
        gmd::RuntimeContext runtime;
        provider->initialize(runtime);
        const auto measured = evaluate(*provider, system).virial;

        std::vector<vr::Vec3> rotated;
        for (const auto& coordinate : base) {
            const vr::Vec3 turned = vr::mat_vec(rotation, vr::sub(coordinate, centre));
            rotated.push_back({turned[0] + centre[0], turned[1] + centre[1],
                               turned[2] + centre[2]});
        }

        gmd::System rotated_system = make_system(lengths, rotated);
        auto rotated_provider = provider_for(fixture);
        rotated_provider->initialize(runtime);
        const auto rotated_measured = evaluate(*rotated_provider, rotated_system).virial;

        const vr::Mat3 expected = vr::rotate_tensor(to_reference(measured), rotation);
        check_tensor(rotated_measured, expected,
                     "[rotation] " + fixture.name + " covariance W -> R W R^T",
                     1.0e-9);
    }
}

}  // namespace

// ===========================================================================
// C. Ewald: real space, reciprocal space, self, net charge, special pairs
// ===========================================================================

namespace {

// A tilted, deliberately irregular charge set. Positions are chosen so that
// every off-diagonal component of every term carries signal.
std::vector<vr::Vec3> charge_positions() {
    return {{3.10L, 4.20L, 2.40L}, {7.31L, 5.13L, 6.24L}, {9.42L, 12.31L, 3.11L},
            {4.23L, 9.94L, 8.72L}, {11.13L, 3.34L, 5.51L}, {6.62L, 12.24L, 9.13L}};
}

std::vector<double> neutral_charges() {
    return {0.7, -0.5, 0.9, -0.4, 0.3, -1.0};
}

std::vector<double> charged_charges() {
    // Sum = +0.6, so the net-charge correction is active.
    return {0.7, -0.5, 0.9, -0.4, 0.3, -0.4};
}

constexpr std::array<double, 3> kEwaldBox = {14.0, 17.0, 11.0};
constexpr double kAlpha = 0.32;
constexpr int kKmax = 6;
// Below every interatomic separation in the fixture, so no real-space pair is
// evaluated and the result is reciprocal + self + net charge alone.
constexpr double kNoRealSpace = 0.5;

std::vector<vr::Real> to_real(const std::vector<double>& values) {
    std::vector<vr::Real> out;
    out.reserve(values.size());
    for (const double value : values) out.push_back(static_cast<vr::Real>(value));
    return out;
}

gmd::System make_charge_system(const std::vector<double>& charges) {
    return make_system(kEwaldBox, charge_positions(), charges);
}

// The Ewald real-space pair force, derived here from
//   V(r) = k_e q_i q_j erfc(alpha r) / r
//   -dV/dr = k_e q_i q_j [ erfc(alpha r)/r^2 + (2 alpha/sqrt(pi)) exp(-(alpha r)^2)/r ]
vr::Mat3 ewald_real_space_reference(const gmd::System& system,
                                    const std::vector<double>& charges,
                                    double alpha, double cutoff) {
    const auto coordinates = coordinates_of(system);
    vr::Mat3 virial{};
    const vr::Real two_alpha_over_root_pi =
        2.0L * static_cast<vr::Real>(alpha) / std::sqrt(std::numbers::pi_v<vr::Real>);

    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        for (std::size_t j = i + 1; j < coordinates.size(); ++j) {
            const vr::Vec3 delta =
                minimum_image(vr::sub(coordinates[i], coordinates[j]), kEwaldBox);
            const vr::Real r_squared = vr::dot(delta, delta);
            if (r_squared >= static_cast<vr::Real>(cutoff) * static_cast<vr::Real>(cutoff)) {
                continue;
            }
            const vr::Real r = std::sqrt(r_squared);
            const vr::Real ar = static_cast<vr::Real>(alpha) * r;
            const vr::Real factor =
                vr::kCoulomb * static_cast<vr::Real>(charges[i]) * static_cast<vr::Real>(charges[j])
                * (std::erfc(ar) / r + two_alpha_over_root_pi * std::exp(-ar * ar))
                / r_squared;
            const vr::Vec3 force = {factor * delta[0], factor * delta[1], factor * delta[2]};
            vr::accumulate_outer(virial, delta, force);
        }
    }
    return virial;
}

// Real space is isolated by difference: the reciprocal, self, net-charge and
// special-pair terms do not depend on the real-space cutoff, so subtracting a
// run with a cutoff below every separation leaves exactly the pair sum.
void test_ewald_real_space_tensor() {
    const double cutoff = 6.5;
    gmd::System system = make_charge_system(neutral_charges());

    gmd::EwaldForceProvider with_pairs(kAlpha, kKmax, cutoff);
    gmd::EwaldForceProvider without_pairs(kAlpha, kKmax, kNoRealSpace);
    const auto full = evaluate(with_pairs, system).virial;
    const auto reciprocal_only = evaluate(without_pairs, system).virial;
    const auto real_space = subtract(full, reciprocal_only);

    const vr::Mat3 reference =
        ewald_real_space_reference(system, neutral_charges(), kAlpha, cutoff);
    check_tensor(real_space, reference, "[analytic-pair] Ewald real space", 1.0e-10);
}

void test_ewald_reciprocal_tensor() {
    for (const bool neutral : {true, false}) {
        const auto charges = neutral ? neutral_charges() : charged_charges();
        const std::string suffix = neutral ? " (neutral)" : " (net charged)";
        gmd::System system = make_charge_system(charges);

        gmd::EwaldForceProvider provider(kAlpha, kKmax, kNoRealSpace);
        const gmd::ForceResult result = evaluate(provider, system);

        auto configuration = vr::make_orthorhombic(
            {static_cast<vr::Real>(kEwaldBox[0]), static_cast<vr::Real>(kEwaldBox[1]),
             static_cast<vr::Real>(kEwaldBox[2])},
            charge_positions(), to_real(charges));

        // Energy first: if the two implementations do not agree on the energy
        // they are not summing the same k-vectors and the tensor comparison
        // would be meaningless.
        const vr::Real reference_energy =
            vr::ewald_reciprocal_energy(configuration, kAlpha, kKmax)
            + vr::ewald_self_energy(to_real(charges), kAlpha)
            + vr::ewald_net_charge_energy(to_real(charges), configuration.volume(), kAlpha);
        const auto energy_error =
            static_cast<double>(std::fabs(reference_energy
                                          - static_cast<vr::Real>(result.potential_energy)));
        check(energy_error < 1.0e-12 * std::fabs(result.potential_energy),
              "[ewald-derived] reciprocal energy" + suffix + " disagrees with the "
              "independent sum (" + std::to_string(energy_error) + ")");

        // Forces, from an independently derived d|S|^2/dr_i.
        const auto reference_forces =
            vr::ewald_reciprocal_forces(configuration, kAlpha, kKmax);
        double worst_force = 0.0;
        double force_scale = 0.0;
        for (std::size_t i = 0; i < reference_forces.size(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) {
                worst_force = std::max(worst_force,
                                       std::fabs(result.forces[i][d]
                                                 - static_cast<double>(reference_forces[i][d])));
                force_scale = std::max(force_scale,
                                       std::fabs(static_cast<double>(reference_forces[i][d])));
            }
        }
        check(worst_force < 1.0e-12 * force_scale,
              "[ewald-derived] reciprocal forces" + suffix +
                  " disagree with the independent derivation (worst " +
                  std::to_string(worst_force) + ")");

        // The tensor, against the independent analytic derivation...
        vr::Mat3 analytic = vr::ewald_reciprocal_virial_analytic(configuration, kAlpha, kKmax);
        const vr::Mat3 net_charge = vr::ewald_net_charge_virial(
            to_real(charges), configuration.volume(), kAlpha);
        for (std::size_t i = 0; i < 9; ++i) analytic[i] += net_charge[i];
        check_tensor(result.virial, analytic,
                     "[ewald-derived] reciprocal tensor vs analytic derivation" + suffix,
                     1.0e-11);

        // ...and against the numerical strain derivative of the same energy,
        // which shares no algebra with either implementation and which reaches
        // the off-diagonal components through a shear the engine cannot apply.
        const auto charge_values = to_real(charges);
        const vr::Mat3 strain = vr::strain_derivative_virial(
            [&](const vr::Mat3& deformation) {
                const auto deformed = vr::deform(configuration, deformation);
                return vr::ewald_reciprocal_energy(deformed, kAlpha, kKmax)
                     + vr::ewald_net_charge_energy(charge_values, deformed.volume(), kAlpha)
                     + vr::ewald_self_energy(charge_values, kAlpha);
            });
        check_tensor(result.virial, strain,
                     "[ewald-derived] reciprocal tensor vs shear strain derivative" + suffix,
                     1.0e-8);
    }
}

// The self term carries energy but no cell dependence, so it must contribute
// exactly nothing to the tensor. Demonstrated by showing that the energy does
// contain it while the tensor equals the reciprocal reference without it.
void test_ewald_self_term_has_no_virial() {
    const auto charges = neutral_charges();
    gmd::System system = make_charge_system(charges);
    gmd::EwaldForceProvider provider(kAlpha, kKmax, kNoRealSpace);
    const gmd::ForceResult result = evaluate(provider, system);

    const auto configuration = vr::make_orthorhombic(
        {static_cast<vr::Real>(kEwaldBox[0]), static_cast<vr::Real>(kEwaldBox[1]),
         static_cast<vr::Real>(kEwaldBox[2])},
        charge_positions(), to_real(charges));

    const vr::Real self_energy = vr::ewald_self_energy(to_real(charges), kAlpha);
    check(std::fabs(static_cast<double>(self_energy)) > 1.0e-3,
          "[ewald-derived] self term is negligible in this fixture, so 'it "
          "contributes no virial' would be untestable");

    const vr::Real reciprocal_energy =
        vr::ewald_reciprocal_energy(configuration, kAlpha, kKmax);
    check(std::fabs(static_cast<double>(result.potential_energy
                                        - static_cast<double>(reciprocal_energy + self_energy)))
              < 1.0e-12 * std::fabs(result.potential_energy),
          "[ewald-derived] the reported energy does not contain the self term, so "
          "this fixture is not testing what it claims");

    const vr::Mat3 reciprocal_only =
        vr::ewald_reciprocal_virial_analytic(configuration, kAlpha, kKmax);
    check_tensor(result.virial, reciprocal_only,
                 "[ewald-derived] self term contributes zero virial", 1.0e-11);
}

void test_ewald_net_charge_correction() {
    const auto charges = charged_charges();
    gmd::System system = make_charge_system(charges);

    gmd::EwaldForceProvider provider(kAlpha, kKmax, kNoRealSpace);
    const auto measured = evaluate(provider, system).virial;

    const auto configuration = vr::make_orthorhombic(
        {static_cast<vr::Real>(kEwaldBox[0]), static_cast<vr::Real>(kEwaldBox[1]),
         static_cast<vr::Real>(kEwaldBox[2])},
        charge_positions(), to_real(charges));

    const vr::Mat3 reciprocal =
        vr::ewald_reciprocal_virial_analytic(configuration, kAlpha, kKmax);
    std::array<double, 9> isolated{};
    for (std::size_t i = 0; i < 9; ++i) {
        isolated[i] = measured[i] - static_cast<double>(reciprocal[i]);
    }

    const vr::Real expected_energy =
        vr::ewald_net_charge_energy(to_real(charges), configuration.volume(), kAlpha);
    check(std::fabs(static_cast<double>(expected_energy)) > 1.0e-4,
          "[ewald-derived] net-charge correction is negligible in this fixture");

    // The analytically expected tensor: isotropic, with U_net on the diagonal.
    const double diagonal = static_cast<double>(expected_energy);
    const double scale = std::fabs(diagonal);
    for (std::size_t a = 0; a < 3; ++a) {
        for (std::size_t b = 0; b < 3; ++b) {
            const double expected = (a == b) ? diagonal : 0.0;
            const double error = std::fabs(isolated[a * 3 + b] - expected);
            check(error < 1.0e-9 * scale,
                  std::string("[ewald-derived] net-charge tensor W_") +
                      component_name(a, b) + " = " + std::to_string(isolated[a * 3 + b]) +
                      ", expected " + std::to_string(expected));
        }
    }

    // A neutral system must produce no such term at all: without this the
    // check above would pass even if the correction were applied always.
    gmd::System neutral = make_charge_system(neutral_charges());
    gmd::EwaldForceProvider neutral_provider(kAlpha, kKmax, kNoRealSpace);
    const auto neutral_measured = evaluate(neutral_provider, neutral).virial;
    const auto neutral_configuration = vr::make_orthorhombic(
        {static_cast<vr::Real>(kEwaldBox[0]), static_cast<vr::Real>(kEwaldBox[1]),
         static_cast<vr::Real>(kEwaldBox[2])},
        charge_positions(), to_real(neutral_charges()));
    const vr::Mat3 neutral_reciprocal =
        vr::ewald_reciprocal_virial_analytic(neutral_configuration, kAlpha, kKmax);
    double residue = 0.0;
    for (std::size_t i = 0; i < 9; ++i) {
        residue = std::max(residue, std::fabs(neutral_measured[i]
                                              - static_cast<double>(neutral_reciprocal[i])));
    }
    check(residue < 1.0e-11 * tensor_scale(neutral_measured),
          "[ewald-derived] a neutral system still receives a net-charge tensor (residue " +
              std::to_string(residue) + ")");
}

// The special-pair Coulomb correction is (scale - 1) times the bare Coulomb
// interaction of the pair, and its virial is the pair tensor of that force.
void test_special_pair_coulomb_tensor() {
    auto topology = std::make_shared<gmd::Topology>();
    topology->bonds = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}};
    topology->angles = {{0, 1, 2, 0}, {1, 2, 3, 0}};
    topology->dihedrals = {{0, 1, 2, 3, 0}};

    const std::vector<vr::Vec3> coordinates = {
        {5.00L, 5.00L, 4.00L}, {6.42L, 5.63L, 4.91L},
        {7.11L, 7.02L, 5.34L}, {8.37L, 7.31L, 6.72L}};
    const std::vector<double> charges = {0.6, -0.4, 0.5, -0.7};
    const double coulomb_scale = 0.8333333333333333;

    gmd::SpecialPairScaleConfig scales;
    scales.pair_12 = {0.0, 0.0};
    scales.pair_13 = {0.0, 0.0};
    scales.pair_14 = {0.5, coulomb_scale};
    auto special_pairs = std::make_shared<gmd::SpecialPairMap>(*topology, scales);

    const std::array<double, 3> lengths = {20.0, 22.0, 18.0};

    gmd::System without = make_system(lengths, coordinates, charges);
    gmd::System with = make_system(lengths, coordinates, charges);
    with.set_special_pair_map(special_pairs);

    gmd::EwaldForceProvider plain(kAlpha, kKmax, 8.0);
    gmd::EwaldForceProvider corrected(kAlpha, kKmax, 8.0);
    const auto baseline = evaluate(plain, without).virial;
    const auto adjusted = evaluate(corrected, with).virial;
    const auto correction = subtract(adjusted, baseline);

    // Reference: the 1-2 and 1-3 pairs are fully excluded (scale 0, so
    // delta = -1) and the 1-4 pair is scaled, each contributing
    // (scale - 1) k_e q_i q_j / r^3 * dr (x) dr.
    vr::Mat3 reference{};
    for (const auto& pair : special_pairs->entries()) {
        const vr::Real delta = static_cast<vr::Real>(pair.scale.coulomb) - 1.0L;
        if (delta == 0.0L) continue;
        const auto a = static_cast<std::size_t>(pair.atom_tag_a);
        const auto b = static_cast<std::size_t>(pair.atom_tag_b);
        const vr::Vec3 separation =
            minimum_image(vr::sub(coordinates[a], coordinates[b]), lengths);
        const vr::Real r_squared = vr::dot(separation, separation);
        const vr::Real r = std::sqrt(r_squared);
        const vr::Real factor = delta * vr::kCoulomb
                              * static_cast<vr::Real>(charges[a])
                              * static_cast<vr::Real>(charges[b]) / (r_squared * r);
        const vr::Vec3 force = {factor * separation[0], factor * separation[1],
                                factor * separation[2]};
        vr::accumulate_outer(reference, separation, force);
    }

    check_tensor(correction, reference,
                 "[analytic-pair] special-pair Coulomb correction", 1.0e-10);

    // The correction must be a real change, not a rounding artefact.
    check(tensor_scale(correction) > 1.0e-6 * tensor_scale(baseline),
          "[analytic-pair] special-pair correction is negligible relative to the "
          "uncorrected tensor, so this fixture cannot detect its omission");
}

// The reciprocal-space term is the one place where sum_i r_i (x) F_i is not
// the virial. Demonstrate that explicitly: if the two agreed, the analytic
// tensor would be unnecessary and the test above would prove nothing.
void test_reciprocal_is_not_the_force_moment() {
    const auto charges = neutral_charges();
    gmd::System system = make_charge_system(charges);
    gmd::EwaldForceProvider provider(kAlpha, kKmax, kNoRealSpace);
    const gmd::ForceResult result = evaluate(provider, system);

    const auto coordinates = coordinates_of(system);
    std::array<double, 9> moment{};
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        for (std::size_t a = 0; a < 3; ++a) {
            for (std::size_t b = 0; b < 3; ++b) {
                moment[a * 3 + b] +=
                    static_cast<double>(coordinates[i][a]) * result.forces[i][b];
            }
        }
    }

    const double difference = max_component_difference(moment, result.virial);
    check(difference > 0.05 * tensor_scale(result.virial),
          "[ewald-derived] sum_i r_i (x) F_i is indistinguishable from the reported "
          "reciprocal tensor in this fixture, so the fixture cannot show that the "
          "cell-derivative form is required (difference " +
              std::to_string(difference) + ")");
}

// A 90-degree axis permutation maps an orthorhombic cell onto another
// orthorhombic cell exactly, so it is the one rotation the reciprocal term can
// be tested against inside the engine's own box representation. A general
// rotation would rotate the charges but not the lattice.
void test_reciprocal_axis_permutation_covariance() {
    const auto charges = neutral_charges();
    gmd::System system = make_charge_system(charges);
    gmd::EwaldForceProvider provider(kAlpha, kKmax, kNoRealSpace);
    const auto measured = evaluate(provider, system).virial;

    // (x, y, z) -> (y, z, x), a proper rotation by 120 degrees about [1,1,1].
    std::vector<vr::Vec3> permuted;
    for (const auto& coordinate : charge_positions()) {
        permuted.push_back({coordinate[2], coordinate[0], coordinate[1]});
    }
    const std::array<double, 3> permuted_box = {kEwaldBox[2], kEwaldBox[0], kEwaldBox[1]};

    gmd::System rotated = make_system(permuted_box, permuted, charges);
    gmd::EwaldForceProvider rotated_provider(kAlpha, kKmax, kNoRealSpace);
    const auto rotated_measured = evaluate(rotated_provider, rotated).virial;

    // R maps e_x -> e_y, e_y -> e_z, e_z -> e_x.
    const vr::Mat3 rotation = {0.0L, 0.0L, 1.0L, 1.0L, 0.0L, 0.0L, 0.0L, 1.0L, 0.0L};
    const vr::Mat3 expected = vr::rotate_tensor(to_reference(measured), rotation);
    check_tensor(rotated_measured, expected,
                 "[rotation] Ewald reciprocal axis-permutation covariance", 1.0e-12);
}

}  // namespace

// ===========================================================================
// D. CompositeForceProvider
// ===========================================================================

namespace {

// A provider with a fixed, known answer, so composition can be checked exactly
// rather than through the noise of a real force field.
class StubProvider final : public gmd::ForceProvider {
public:
    StubProvider(std::array<double, 9> virial, bool valid, double energy)
        : virial_(virial), valid_(valid), energy_(energy) {}

    std::string_view name() const noexcept override { return "stub"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}

    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        const std::size_t n = request.coordinates.size();
        result.success = true;
        result.potential_energy = energy_;
        result.forces.assign(n, gmd::Force3D{0.0, 0.0, 0.0});
        result.virial = virial_;
        result.virial_valid = valid_;
    }

private:
    std::array<double, 9> virial_;
    bool valid_;
    double energy_;
};

void test_composite_adds_components_exactly() {
    const std::array<double, 9> first = {1.5, -2.25, 3.125, -0.5, 4.75, -6.375,
                                         7.0625, -8.5, 9.25};
    const std::array<double, 9> second = {-0.25, 1.125, -2.5, 3.75, -4.0625, 5.5,
                                          -6.75, 7.875, -8.125};

    auto composite = std::make_shared<gmd::CompositeForceProvider>();
    composite->add(std::make_shared<StubProvider>(first, true, 1.0));
    composite->add(std::make_shared<StubProvider>(second, true, 2.0));

    gmd::System system = make_system({20.0, 20.0, 20.0}, tilted_pair());
    const gmd::ForceResult result = evaluate(*composite, system);

    check(result.virial_valid,
          "[symmetry] composite: two valid children must give a valid virial");
    for (std::size_t i = 0; i < 9; ++i) {
        // Exact equality, not a tolerance: the values above are dyadic, the
        // composite accumulates in provider order, and the sum is therefore
        // representable. Anything other than exact addition is a defect.
        check(result.virial[i] == first[i] + second[i],
              "[symmetry] composite: component " + std::to_string(i) +
                  " is not the exact sum of its children (" +
                  std::to_string(result.virial[i]) + " vs " +
                  std::to_string(first[i] + second[i]) + ")");
    }
    check(result.potential_energy == 3.0,
          "[symmetry] composite: energies must add");
}

void test_composite_propagates_invalid_child() {
    const std::array<double, 9> valid = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    const std::array<double, 9> ignored = {100.0, 0.0, 0.0, 0.0, 100.0, 0.0,
                                           0.0, 0.0, 100.0};

    for (const bool invalid_first : {true, false}) {
        auto composite = std::make_shared<gmd::CompositeForceProvider>();
        if (invalid_first) {
            composite->add(std::make_shared<StubProvider>(ignored, false, 1.0));
            composite->add(std::make_shared<StubProvider>(valid, true, 2.0));
        } else {
            composite->add(std::make_shared<StubProvider>(valid, true, 2.0));
            composite->add(std::make_shared<StubProvider>(ignored, false, 1.0));
        }

        gmd::System system = make_system({20.0, 20.0, 20.0}, tilted_pair());
        const gmd::ForceResult result = evaluate(*composite, system);

        // Order must not matter: an invalid child anywhere invalidates the sum.
        check(!result.virial_valid,
              std::string("[symmetry] composite: a child without a virial must ") +
                  "invalidate the composite (invalid child " +
                  (invalid_first ? "first" : "second") + ")");
        check(result.success,
              "[symmetry] composite: an absent virial is not a failed evaluation");
        check(result.potential_energy == 3.0,
              "[symmetry] composite: energy must still add when a virial is absent");
    }
}

// A composite of real providers must equal the sum of the same providers run
// separately -- the check that composition does not re-reduce or double-count.
void test_composite_matches_separate_providers() {
    gmd::System system = make_charge_system(neutral_charges());

    auto lennard_jones_provider =
        std::make_shared<gmd::ClassicalForceProvider>(0.0104, 3.4, 6.5);
    auto ewald_provider =
        std::make_shared<gmd::EwaldForceProvider>(kAlpha, kKmax, 6.5);

    const auto lj = evaluate(*lennard_jones_provider, system).virial;
    const auto ewald = evaluate(*ewald_provider, system).virial;

    auto composite = std::make_shared<gmd::CompositeForceProvider>();
    composite->add(std::make_shared<gmd::ClassicalForceProvider>(0.0104, 3.4, 6.5));
    composite->add(std::make_shared<gmd::EwaldForceProvider>(kAlpha, kKmax, 6.5));
    const auto combined = evaluate(*composite, system).virial;

    std::array<double, 9> expected{};
    for (std::size_t i = 0; i < 9; ++i) expected[i] = lj[i] + ewald[i];

    const double difference = max_component_difference(combined, expected);
    check(difference < 1.0e-12 * tensor_scale(expected),
          "[symmetry] composite of LJ + Ewald does not equal the sum of the two run "
          "separately (worst component difference " + std::to_string(difference) + ")");
    check(tensor_scale(lj) > 0.0 && tensor_scale(ewald) > 0.0,
          "[symmetry] composite fixture: both children must contribute for this to "
          "test addition");
}

}  // namespace

int main() {
    check_lennard_jones_force_against_derivative();
    test_lennard_jones_pair_tensor();
    test_lennard_jones_mixed_types_and_override();
    test_lennard_jones_special_pair_scaling();
    test_pair_virial_across_periodic_boundaries();
    test_pair_virial_is_translation_invariant();
    test_pair_virial_rotation_covariance();

    test_bonded_term_tensors();
    test_bonded_wrapped_equals_intact();
    test_bonded_rotation_covariance();

    test_ewald_real_space_tensor();
    test_ewald_reciprocal_tensor();
    test_ewald_self_term_has_no_virial();
    test_ewald_net_charge_correction();
    test_special_pair_coulomb_tensor();
    test_reciprocal_is_not_the_force_moment();
    test_reciprocal_axis_permutation_covariance();

    test_composite_adds_components_exactly();
    test_composite_propagates_invalid_child();
    test_composite_matches_separate_providers();

    print_error_summary();

    if (failures != 0) {
        std::cerr << "virial source inventory tests failed: " << failures << '\n';
        return 1;
    }
    std::cout << "\nvirial source inventory tests passed\n";
    return 0;
}
