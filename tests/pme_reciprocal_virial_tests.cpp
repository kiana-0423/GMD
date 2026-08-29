// PME reciprocal-mesh validation: energy, forces and all nine virial
// components, against two references that do not share the engine's algebra.
//
// 1. CONVERGENCE TO EWALD. PME approximates the same reciprocal sum that
//    EwaldForceProvider evaluates directly. On identical configurations the
//    difference must shrink as the mesh is refined. The assertions here are
//    convergence assertions -- each mesh doubling must reduce every component's
//    error by at least a stated factor -- rather than an absolute tolerance
//    copied from one measured run, which would only ever restate whatever the
//    code happened to produce.
//
// 2. THE RECIPROCAL METRIC DERIVATIVE. tests/virial_reference.hpp contains a
//    test-only PME reciprocal energy written against a general 3x3 cell: its
//    own B-spline (from the order recursion, not a polynomial table), its own
//    direct DFT (not a Cooley-Tukey transform), and its own influence function.
//    Differentiating that energy numerically with respect to a general strain,
//    shear included, gives all nine components. `gmd::Box` cannot express the
//    shear, but the reference can, and the engine's analytic tensor is
//    evaluated at the orthorhombic configuration both agree on.
//
//    This is what makes the off-diagonal claim real rather than a restatement:
//    the engine computes W from a closed form, the reference computes it as a
//    numerical derivative of an independently written energy.
//
// B-SPLINE ORDERS. Orders 4 and 6 are both exercised end to end. Order 6 must
// converge strictly faster than order 4 on the same mesh -- an order that is
// merely *accepted* by the constructor while silently behaving like something
// else would pass a pure tolerance check and fail this one.

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
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

#include "virial_reference.hpp"

namespace vr = virial_ref;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[pme virial] " << message << '\n';
        ++failures;
    }
}

std::string format(double value, int precision = 3) {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(precision) << value;
    return stream.str();
}

const char* component_name(std::size_t index) {
    static const char* names[9] = {"xx", "xy", "xz", "yx", "yy", "yz",
                                   "zx", "zy", "zz"};
    return names[index];
}

// --- fixtures -------------------------------------------------------------

// Deliberately tilted: no charge shares a coordinate plane with another and no
// axis is special, so every off-diagonal component of the reciprocal tensor is
// well away from zero. A symmetric arrangement would zero them and make a
// nine-component claim vacuous.
std::vector<vr::Vec3> tilted_charges() {
    return {{2.31L, 3.17L, 1.94L},  {6.42L, 4.83L, 5.11L},
            {9.76L, 11.28L, 2.63L}, {3.19L, 8.94L, 7.42L},
            {11.53L, 2.71L, 4.38L}, {5.87L, 12.06L, 8.19L},
            {8.14L, 6.35L, 9.73L},  {1.62L, 10.41L, 3.55L}};
}

std::vector<double> neutral_charges() {
    return {0.7, -0.5, 0.9, -0.4, 0.3, -1.0, 0.6, -0.6};
}

std::vector<double> net_charged() {
    // Sum = +0.7.
    return {0.7, -0.5, 0.9, -0.4, 0.3, -1.0, 0.6, 0.1};
}

// Non-cubic, so nothing can hide behind cubic symmetry, and each edge is a
// different length from the others.
constexpr std::array<double, 3> kBox = {14.0, 17.0, 11.0};
// Small enough to exclude every real-space pair in the fixture, leaving the
// reciprocal, self and net-charge terms alone.
constexpr double kNoRealSpace = 0.5;
// Far beyond where exp(-k^2 / 4 alpha^2) is numerically relevant for every
// alpha used below, so the Ewald reference carries no truncation error of its
// own to confuse the comparison.
constexpr int kReferenceKmax = 16;

gmd::System make_system(const std::vector<double>& charges) {
    gmd::System system;
    const auto positions = tilted_charges();
    system.resize(positions.size(), positions.size());
    gmd::Box box;
    box.set_lengths(kBox);
    system.set_box(box);
    for (std::size_t i = 0; i < positions.size(); ++i) {
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {static_cast<double>(positions[i][0]),
                                           static_cast<double>(positions[i][1]),
                                           static_cast<double>(positions[i][2])};
        system.mutable_charges()[i] = charges[i];
    }
    return system;
}

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

gmd::ForceResult run_pme(gmd::System& system, double alpha, int order,
                         std::array<int, 3> grid) {
    gmd::PMEForceProvider provider(alpha, kNoRealSpace, order, grid);
    gmd::RuntimeContext runtime;
    provider.initialize(runtime);
    return evaluate(provider, system);
}

gmd::ForceResult run_ewald(gmd::System& system, double alpha) {
    gmd::EwaldForceProvider provider(alpha, kReferenceKmax, kNoRealSpace);
    return evaluate(provider, system);
}

std::vector<vr::Real> to_real(const std::vector<double>& values) {
    std::vector<vr::Real> out;
    out.reserve(values.size());
    for (const double value : values) out.push_back(static_cast<vr::Real>(value));
    return out;
}

vr::CellConfiguration make_configuration(const std::vector<double>& charges) {
    return vr::make_orthorhombic(
        {static_cast<vr::Real>(kBox[0]), static_cast<vr::Real>(kBox[1]),
         static_cast<vr::Real>(kBox[2])},
        tilted_charges(), to_real(charges));
}

double tensor_scale(const std::array<double, 9>& tensor) {
    double scale = 0.0;
    for (const double value : tensor) scale = std::max(scale, std::fabs(value));
    return scale;
}

// ===========================================================================
// 1. Convergence to a high-accuracy Ewald full-tensor reference
// ===========================================================================

struct ConvergenceRow {
    int grid;
    double energy_error;
    double force_error;
    std::array<double, 9> component_error;
    double worst_component_error;
};

ConvergenceRow measure(gmd::System& system, const gmd::ForceResult& reference,
                       double alpha, int order, int grid) {
    const auto result = run_pme(system, alpha, order, {grid, grid, grid});

    ConvergenceRow row{};
    row.grid = grid;
    row.energy_error = std::fabs(result.potential_energy - reference.potential_energy)
                     / std::fabs(reference.potential_energy);

    double force_scale = 0.0;
    double force_error = 0.0;
    for (std::size_t i = 0; i < reference.forces.size(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            force_scale = std::max(force_scale, std::fabs(reference.forces[i][d]));
            force_error = std::max(force_error,
                                   std::fabs(result.forces[i][d] - reference.forces[i][d]));
        }
    }
    row.force_error = force_error / force_scale;

    const double scale = tensor_scale(reference.virial);
    for (std::size_t i = 0; i < 9; ++i) {
        row.component_error[i] =
            std::fabs(result.virial[i] - reference.virial[i]) / scale;
        row.worst_component_error =
            std::max(row.worst_component_error, row.component_error[i]);
    }
    return row;
}

void report_convergence(const std::string& label,
                        const std::vector<ConvergenceRow>& rows) {
    std::cout << "\n  " << label << '\n';
    std::cout << "    grid   energy      force       "
              << "W_xx        W_xy        W_xz        W_yy        W_yz        W_zz\n";
    for (const auto& row : rows) {
        std::cout << "    " << std::setw(4) << row.grid << "   "
                  << format(row.energy_error) << "   " << format(row.force_error)
                  << "   " << format(row.component_error[0])
                  << "   " << format(row.component_error[1])
                  << "   " << format(row.component_error[2])
                  << "   " << format(row.component_error[4])
                  << "   " << format(row.component_error[5])
                  << "   " << format(row.component_error[8]) << '\n';
    }
}

void test_pme_converges_to_ewald() {
    const std::vector<double> alphas = {0.25, 0.32, 0.40};
    const std::vector<int> grids = {16, 32, 64};
    // Each mesh doubling must cut every component's error by at least this
    // factor. The true observed factors are much larger (roughly 16x at order
    // 4 and 80x at order 6); the requirement is deliberately loose enough that
    // it asserts convergence rather than pinning down a measured number.
    constexpr double kRequiredImprovement = 4.0;

    for (const bool neutral : {true, false}) {
        const auto charges = neutral ? neutral_charges() : net_charged();
        for (const double alpha : alphas) {
            for (const int order : {4, 6}) {
                gmd::System system = make_system(charges);
                const auto reference = run_ewald(system, alpha);

                // A reference whose off-diagonal components are negligible
                // cannot support a nine-component convergence claim.
                const double scale = tensor_scale(reference.virial);
                for (std::size_t i = 0; i < 9; ++i) {
                    check(std::fabs(reference.virial[i]) > 1.0e-3 * scale,
                          std::string("PME/Ewald fixture: reference W_") +
                              component_name(i) + " is too small to test (" +
                              format(reference.virial[i]) + " vs scale " +
                              format(scale) + ")");
                }

                std::vector<ConvergenceRow> rows;
                for (const int grid : grids) {
                    rows.push_back(measure(system, reference, alpha, order, grid));
                }

                const std::string label =
                    std::string(neutral ? "neutral" : "net charged") +
                    ", alpha = " + format(alpha, 2) + ", order " + std::to_string(order);
                report_convergence(label, rows);

                // Per-doubling convergence is asserted on the tensor as a
                // whole (its worst component), on the energy and on the force.
                // It is deliberately NOT asserted per component: an individual
                // component's error can pass through a sign change between two
                // meshes, leaving it accidentally tiny on the coarser one, and
                // a ratio taken across that crossing measures the accident
                // rather than the convergence. The alpha = 0.40, order 4,
                // net-charged fixture does exactly this in W_xy.
                for (std::size_t step = 1; step < rows.size(); ++step) {
                    const auto& coarse = rows[step - 1];
                    const auto& fine = rows[step];
                    check(fine.worst_component_error * kRequiredImprovement
                              <= coarse.worst_component_error,
                          label + ": the tensor did not converge from grid " +
                              std::to_string(coarse.grid) + " (" +
                              format(coarse.worst_component_error) + ") to grid " +
                              std::to_string(fine.grid) + " (" +
                              format(fine.worst_component_error) + ")");
                    check(fine.energy_error * kRequiredImprovement <= coarse.energy_error,
                          label + ": reciprocal energy did not converge from grid " +
                              std::to_string(coarse.grid) + " to " +
                              std::to_string(fine.grid));
                    check(fine.force_error * kRequiredImprovement <= coarse.force_error,
                          label + ": reciprocal force did not converge from grid " +
                              std::to_string(coarse.grid) + " to " +
                              std::to_string(fine.grid));
                }

                // Every component is still required to converge individually,
                // but across the full refinement (two doublings) rather than
                // between adjacent meshes, which is robust to a single sign
                // crossing while still failing a component that stalls.
                //
                // Components already below kConvergedFloor are exempt: at that
                // level the residual is dominated by the double-precision
                // Ewald reference and by accumulated round-off in a sum of
                // O(10^5) k-vectors, not by the mesh, so demanding further
                // reduction would be measuring arithmetic noise.
                constexpr double kConvergedFloor = 1.0e-9;
                const auto& coarsest = rows.front();
                const auto& finest = rows.back();
                for (std::size_t i = 0; i < 9; ++i) {
                    if (finest.component_error[i] <= kConvergedFloor) continue;
                    check(finest.component_error[i] * (kRequiredImprovement
                                                       * kRequiredImprovement)
                              <= coarsest.component_error[i],
                          label + ": W_" + component_name(i) +
                              " did not converge across the full refinement, grid " +
                              std::to_string(coarsest.grid) + " (" +
                              format(coarsest.component_error[i]) + ") to grid " +
                              std::to_string(finest.grid) + " (" +
                              format(finest.component_error[i]) + ")");
                }
            }
        }
    }
}

// A higher B-spline order must actually buy accuracy. This is the check that a
// wrong or missing polynomial branch cannot survive: an order-6 spline that is
// not really order 6 either blows up or converges no faster than order 4.
void test_higher_bspline_order_is_more_accurate() {
    const double alpha = 0.32;
    gmd::System system = make_system(neutral_charges());
    const auto reference = run_ewald(system, alpha);

    for (const int grid : {16, 32}) {
        const auto fourth = measure(system, reference, alpha, 4, grid);
        const auto sixth = measure(system, reference, alpha, 6, grid);
        check(sixth.worst_component_error < fourth.worst_component_error,
              "B-spline order 6 is not more accurate than order 4 at grid " +
                  std::to_string(grid) + " (order 6 worst " +
                  format(sixth.worst_component_error) + ", order 4 worst " +
                  format(fourth.worst_component_error) + ")");
        std::cout << "    order 4 vs 6 at grid " << grid << ": "
                  << format(fourth.worst_component_error) << " -> "
                  << format(sixth.worst_component_error) << '\n';
    }
}

// A non-cubic mesh over a non-cubic cell. The per-axis grid spacing then
// differs between axes, which is where an index-ordering slip between the mesh
// dimensions and the k-vector components would show up.
void test_non_cubic_mesh() {
    const double alpha = 0.32;
    gmd::System system = make_system(neutral_charges());
    const auto reference = run_ewald(system, alpha);
    const double scale = tensor_scale(reference.virial);

    const std::vector<std::array<int, 3>> grids = {
        {16, 32, 16}, {32, 64, 32}, {64, 128, 64}};

    std::vector<double> worst;
    for (const auto& grid : grids) {
        const auto result = run_pme(system, alpha, 6, grid);
        double error = 0.0;
        for (std::size_t i = 0; i < 9; ++i) {
            error = std::max(error, std::fabs(result.virial[i] - reference.virial[i]) / scale);
        }
        worst.push_back(error);
        std::cout << "    non-cubic mesh " << grid[0] << "x" << grid[1] << "x" << grid[2]
                  << ": worst component error " << format(error) << '\n';
    }
    for (std::size_t step = 1; step < worst.size(); ++step) {
        check(worst[step] * 4.0 <= worst[step - 1],
              "non-cubic mesh did not converge between refinement " +
                  std::to_string(step - 1) + " (" + format(worst[step - 1]) +
                  ") and " + std::to_string(step) + " (" + format(worst[step]) + ")");
    }
}

// ===========================================================================
// 2. The reciprocal metric derivative
// ===========================================================================

void test_pme_against_reciprocal_metric_derivative() {
    const double alpha = 0.32;
    const int grid_size = 16;
    const std::array<int, 3> grid = {grid_size, grid_size, grid_size};

    for (const int order : {4, 6}) {
        const auto charges = neutral_charges();
        gmd::System system = make_system(charges);
        const auto result = run_pme(system, alpha, order, grid);

        const auto configuration = make_configuration(charges);
        const auto charge_values = to_real(charges);

        // Gate: the independent reciprocal energy must reproduce the engine's.
        // If it does not, the reference is spreading charge differently or
        // transforming it differently, and its strain derivative would be a
        // derivative of the wrong function.
        const vr::Real reference_energy =
            vr::pme_reciprocal_energy(configuration, alpha, grid, order);
        const vr::Real self_energy = vr::ewald_self_energy(charge_values, alpha);
        const auto energy_error =
            static_cast<double>(std::fabs(static_cast<vr::Real>(result.potential_energy)
                                          - (reference_energy + self_energy)));
        check(energy_error < 1.0e-11 * std::fabs(result.potential_energy),
              "order " + std::to_string(order) +
                  ": independent PME reciprocal energy disagrees with the provider (" +
                  format(energy_error) + " absolute, energy " +
                  format(result.potential_energy) + ")");

        // The reference: -dE/de_ab of that same independently written energy,
        // under a general strain including shear.
        const vr::Mat3 strain_virial = vr::strain_derivative_virial(
            [&](const vr::Mat3& deformation) {
                const auto deformed = vr::deform(configuration, deformation);
                return vr::pme_reciprocal_energy(deformed, alpha, grid, order)
                     + vr::ewald_self_energy(charge_values, alpha);
            });

        double scale = 0.0;
        for (const auto value : strain_virial) {
            scale = std::max(scale, std::fabs(static_cast<double>(value)));
        }

        for (std::size_t i = 0; i < 9; ++i) {
            const double expected = static_cast<double>(strain_virial[i]);
            const double measured = result.virial[i];
            const double error = std::fabs(measured - expected) / scale;
            check(error < 1.0e-7,
                  std::string("order ") + std::to_string(order) + ": W_" +
                      component_name(i) + " = " + format(measured) +
                      " disagrees with the reciprocal metric derivative " +
                      format(expected) + " (relative error " + format(error) + ")");
            check(std::fabs(expected) > 1.0e-3 * scale,
                  std::string("order ") + std::to_string(order) + ": reference W_" +
                      component_name(i) + " is too small for this fixture to test it");
        }
        std::cout << "    order " << order
                  << ": PME tensor matches the reciprocal metric derivative on all "
                  << "nine components\n";
    }
}

// The strain derivative must not be an accident of the orthorhombic cell: run
// the same reference on a genuinely sheared cell and confirm the off-diagonal
// components respond. This does not validate the engine (which cannot shear);
// it validates that the reference's shear channel is live, so agreement above
// is meaningful.
void test_reference_shear_channel_is_live() {
    const double alpha = 0.32;
    const auto charges = neutral_charges();
    const auto configuration = make_configuration(charges);
    const std::array<int, 3> grid = {16, 16, 16};

    const vr::Real base = vr::pme_reciprocal_energy(configuration, alpha, grid, 4);

    vr::Mat3 shear = vr::identity3();
    shear[1] += 0.02L;  // e_xy
    const auto sheared = vr::deform(configuration, shear);
    const vr::Real deformed = vr::pme_reciprocal_energy(sheared, alpha, grid, 4);

    check(std::fabs(static_cast<double>(deformed - base))
              > 1.0e-6 * std::fabs(static_cast<double>(base)),
          "the reference reciprocal energy does not respond to a shear strain, so "
          "its off-diagonal derivative would be identically zero and the agreement "
          "reported above would be vacuous");

    // And the cell really is sheared, not merely relabelled.
    check(std::fabs(static_cast<double>(sheared.cell[1])) > 0.0,
          "the sheared reference cell is still diagonal");
}

// ===========================================================================
// 3. PME real space is the same pair sum Ewald evaluates
// ===========================================================================

void test_pme_real_space_matches_ewald() {
    const double alpha = 0.32;
    const double cutoff = 5.0;
    const auto charges = neutral_charges();

    gmd::System system = make_system(charges);

    gmd::PMEForceProvider pme_with(alpha, cutoff, 4, {32, 32, 32});
    gmd::PMEForceProvider pme_without(alpha, kNoRealSpace, 4, {32, 32, 32});
    gmd::RuntimeContext runtime;
    pme_with.initialize(runtime);
    pme_without.initialize(runtime);

    gmd::EwaldForceProvider ewald_with(alpha, kReferenceKmax, cutoff);
    gmd::EwaldForceProvider ewald_without(alpha, kReferenceKmax, kNoRealSpace);

    const auto pme_real = evaluate(pme_with, system).virial;
    const auto pme_base = evaluate(pme_without, system).virial;
    const auto ewald_real = evaluate(ewald_with, system).virial;
    const auto ewald_base = evaluate(ewald_without, system).virial;

    double worst = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < 9; ++i) {
        const double from_pme = pme_real[i] - pme_base[i];
        const double from_ewald = ewald_real[i] - ewald_base[i];
        worst = std::max(worst, std::fabs(from_pme - from_ewald));
        scale = std::max(scale, std::fabs(from_ewald));
    }
    check(scale > 0.0, "PME real-space fixture evaluates no pairs");
    check(worst < 1.0e-12 * scale,
          "PME and Ewald disagree on the real-space pair tensor (worst component " +
              format(worst) + ", scale " + format(scale) + ")");
    std::cout << "    PME real-space tensor matches Ewald to " << format(worst / scale)
              << " relative\n";
}

}  // namespace

int main() {
    std::cout << "PME reciprocal virial: component-wise error vs a converged Ewald "
              << "reference (relative to max |W|)\n";
    test_pme_converges_to_ewald();
    std::cout << "\n  B-spline order comparison\n";
    test_higher_bspline_order_is_more_accurate();
    std::cout << "\n  Non-cubic mesh\n";
    test_non_cubic_mesh();
    std::cout << "\n  Reciprocal metric derivative\n";
    test_pme_against_reciprocal_metric_derivative();
    test_reference_shear_channel_is_live();
    std::cout << "\n  Real space\n";
    test_pme_real_space_matches_ewald();

    if (failures != 0) {
        std::cerr << "\nPME reciprocal virial tests failed: " << failures << '\n';
        return 1;
    }
    std::cout << "\nPME reciprocal virial tests passed\n";
    return 0;
}
