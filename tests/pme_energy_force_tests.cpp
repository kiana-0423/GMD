// Direct regression tests for the production PME energy and force path.
//
// Two independent claims, in this order, because the second is worthless
// without the first:
//
//   1. THE TRANSFORM CONVENTION. The reciprocal potential and its coordinate
//      gradient are recomputed here by explicit direct DFTs with a stated
//      normalisation (tests/pme_reference.hpp), on a small mesh, and compared
//      against the production FFT path. This pins the K1*K2*K3 factor in the
//      force interpolation: remove it and every force is N times too small,
//      apply it twice and every force is N times too large. Either fails here
//      immediately, on an absolute comparison, with no appeal to convergence.
//
//   2. ENERGY-FORCE CONSISTENCY. F_i,a = -dU/dr_i,a for the provider's own
//      energy, by central finite differences.
//
//      Mesh-to-Ewald convergence, which lives in
//      tests/pme_reciprocal_virial_tests.cpp, cannot substitute for either.
//      Convergence compares PME against a different method and passes as long
//      as the two agree in the limit; it says nothing about whether the force
//      the engine hands the integrator is the gradient of the energy it
//      reports. The historical defects were exactly that: the reciprocal
//      energy was correct at order 4 the whole time, while the reciprocal force
//      was first zero and then N times too small.
//
// PARAMETER STABILITY. `PMEForceProvider::resolve_params()` fills in alpha and
// the real-space cutoff when either is non-positive, from the box. Every case
// below passes strictly positive values for both, so resolution is a no-op and
// nothing can change between the plus and minus evaluations of a difference.
// One provider instance is reused across each derivative, and the box, grid and
// spline order are fixed for its lifetime. `test_parameters_are_stable()`
// checks that repeated evaluation is bitwise reproducible, which is what would
// break first if any of that stopped being true.

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

#include "pme_reference.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[pme energy-force] " << message << '\n';
        ++failures;
    }
}

std::string format(double value, int precision = 3) {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(precision) << value;
    return stream.str();
}

// --- fixtures -------------------------------------------------------------

// Non-cubic, and every edge a different length.
constexpr std::array<double, 3> kBox = {16.0, 24.0, 12.0};

// Fully asymmetric: no two atoms share a coordinate, nothing sits on a plane
// of symmetry, and no force component of any atom comes out near zero. A
// symmetric arrangement would let a whole component of the gradient go
// untested.
std::vector<std::array<double, 3>> asymmetric_positions() {
    return {{2.37, 3.11, 1.83}, {6.41, 7.93, 5.17}, {11.62, 14.28, 3.44},
            {4.19, 19.07, 8.71}, {13.53, 5.62, 10.26}, {8.84, 11.35, 6.98},
            {1.26, 21.44, 4.05}, {14.77, 16.81, 9.32}};
}

// One atom within 0.07 A of each face, one just inside a corner, and one in
// the interior. Under a finite difference these step across the boundary, so
// the wrapping arithmetic in the charge spreading is exercised by the
// derivative itself rather than only by a static evaluation.
//
// The atoms are spread along the faces they sit on rather than stacked, so the
// closest approach in the fixture is about 4.9 A. Placing two of them against
// opposite faces at the same transverse position would put them 0.1 A apart
// through the boundary, where the Coulomb force is four orders of magnitude
// larger than anything else in the cell -- which does not test wrapping, it
// just swamps every other component.
std::vector<std::array<double, 3>> boundary_positions() {
    return {{0.05, 4.00, 3.00},   {15.95, 12.00, 9.00},
            {6.00, 0.06, 8.00},   {10.00, 23.94, 2.00},
            {3.00, 16.00, 0.04},  {12.00, 8.00, 11.96},
            {0.07, 0.05, 0.06},   {8.30, 14.70, 5.90}};
}

std::vector<double> neutral_charges() {
    return {0.7, -0.5, 0.9, -0.4, 0.3, -1.0, 0.6, -0.6};
}

std::vector<double> net_charges() {
    // Sum = +0.7.
    return {0.7, -0.5, 0.9, -0.4, 0.3, -1.0, 0.6, 0.1};
}

gmd::System make_system(const std::vector<std::array<double, 3>>& positions,
                        const std::vector<double>& charges) {
    gmd::System system;
    system.resize(positions.size(), positions.size());
    gmd::Box box;
    box.set_lengths(kBox);
    system.set_box(box);
    for (std::size_t i = 0; i < positions.size(); ++i) {
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = positions[i];
        system.mutable_charges()[i] = charges[i];
    }
    return system;
}

gmd::ForceResult evaluate(gmd::PMEForceProvider& provider, gmd::System& system) {
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

// Below every interatomic separation in both fixtures, so no real-space pair
// is evaluated and the result is the reciprocal mesh plus the position-
// independent self and net-charge terms.
constexpr double kNoRealSpace = 0.5;

// ===========================================================================
// 1. The transform convention
// ===========================================================================

void test_fft_normalisation_against_direct_dft() {
    // A small mesh, so a direct DFT is cheap and the comparison is exact
    // rather than asymptotic.
    const std::array<int, 3> grid = {8, 8, 8};
    const double alpha = 0.35;

    std::cout << "\n  FFT normalisation vs explicit direct DFT (grid 8x8x8)\n";
    std::cout << "    order  charges     dE_recip     dE_parseval   max|dF|      rel|dF|\n";

    for (const int order : {4, 6}) {
        for (const bool neutral : {true, false}) {
            const auto charges = neutral ? neutral_charges() : net_charges();
            const auto positions = asymmetric_positions();

            gmd::System system = make_system(positions, charges);
            gmd::PMEForceProvider provider(alpha, kNoRealSpace, order, grid);
            gmd::RuntimeContext runtime;
            provider.initialize(runtime);
            const auto result = evaluate(provider, system);

            std::vector<std::array<pme_ref::Real, 3>> reference_positions;
            for (const auto& p : positions) {
                reference_positions.push_back({static_cast<pme_ref::Real>(p[0]),
                                               static_cast<pme_ref::Real>(p[1]),
                                               static_cast<pme_ref::Real>(p[2])});
            }
            std::vector<pme_ref::Real> reference_charges;
            for (const double q : charges) {
                reference_charges.push_back(static_cast<pme_ref::Real>(q));
            }
            const std::array<pme_ref::Real, 3> lengths = {
                static_cast<pme_ref::Real>(kBox[0]), static_cast<pme_ref::Real>(kBox[1]),
                static_cast<pme_ref::Real>(kBox[2])};

            const auto reference = pme_ref::reciprocal_reference(
                reference_positions, reference_charges, lengths, alpha, grid, order);

            // Two independent expressions for the same reciprocal energy: the
            // k-space form the engine uses, and the real-space Parseval form,
            // which is only equal to it under the stated transform convention.
            const pme_ref::Real volume = lengths[0] * lengths[1] * lengths[2];
            const pme_ref::Real corrections =
                pme_ref::self_energy(reference_charges, alpha)
                + pme_ref::net_charge_energy(reference_charges, volume, alpha);

            const double energy_reciprocal = std::fabs(
                result.potential_energy
                - static_cast<double>(reference.energy_reciprocal_space + corrections));
            const double energy_parseval = std::fabs(
                result.potential_energy
                - static_cast<double>(reference.energy_real_space + corrections));

            check(energy_reciprocal < 1.0e-11 * std::fabs(result.potential_energy),
                  "order " + std::to_string(order) +
                      ": reciprocal-space energy disagrees with the direct DFT by " +
                      format(energy_reciprocal));
            check(energy_parseval < 1.0e-10 * std::fabs(result.potential_energy),
                  "order " + std::to_string(order) +
                      ": the real-space Parseval form of the energy disagrees with the "
                      "provider by " + format(energy_parseval) +
                      ", which means the forward and inverse transform normalisations "
                      "are not consistent with each other");

            // The force. This is the assertion that pins K1*K2*K3: the
            // reference's inverse transform carries no 1/N, so its gradient
            // needs no compensating factor, and a provider that omits or
            // doubles one is off by exactly N = 512 here.
            double worst = 0.0;
            double scale = 0.0;
            for (std::size_t i = 0; i < positions.size(); ++i) {
                for (std::size_t d = 0; d < 3; ++d) {
                    const auto expected = static_cast<double>(reference.forces[i][d]);
                    worst = std::max(worst, std::fabs(result.forces[i][d] - expected));
                    scale = std::max(scale, std::fabs(expected));
                }
            }
            check(scale > 1.0e-4,
                  "order " + std::to_string(order) +
                      ": the reference reciprocal force is negligible, so this fixture "
                      "cannot detect a normalisation error");
            check(worst < 1.0e-10 * scale,
                  "order " + std::to_string(order) +
                      ": reciprocal force disagrees with the direct-DFT gradient by " +
                      format(worst) + " (scale " + format(scale) + ", ratio " +
                      format(worst / scale) + "). A missing or doubled mesh-point-count "
                      "factor would show up here as a ratio near " +
                      std::to_string(grid[0] * grid[1] * grid[2]) + ".");

            std::cout << "    " << std::setw(5) << order << "  "
                      << (neutral ? "neutral   " : "net +0.7  ") << "  "
                      << format(energy_reciprocal) << "   " << format(energy_parseval)
                      << "   " << format(worst) << "   " << format(worst / scale) << '\n';
        }
    }
}

// ===========================================================================
// 2. Energy-force consistency
// ===========================================================================

struct GradientReport {
    std::string label;
    double max_absolute = 0.0;
    double max_relative = 0.0;
    double force_scale = 0.0;
};

// Central differences of the provider's own energy. `provider` is used for
// every evaluation, so nothing about the parameters can drift between the two
// sides of a difference.
GradientReport check_energy_force_consistency(
        const std::string& label,
        gmd::PMEForceProvider& provider,
        const std::vector<std::array<double, 3>>& positions,
        const std::vector<double>& charges,
        double tolerance,
        bool require_every_component_nonzero,
        double step = 2.0e-5) {
    gmd::System system = make_system(positions, charges);
    const auto result = evaluate(provider, system);

    GradientReport report;
    report.label = label;

    for (std::size_t i = 0; i < positions.size(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            auto displaced = positions;
            displaced[i][d] += step;
            gmd::System up = make_system(displaced, charges);
            displaced[i][d] -= 2.0 * step;
            gmd::System down = make_system(displaced, charges);

            const double energy_up = evaluate(provider, up).potential_energy;
            const double energy_down = evaluate(provider, down).potential_energy;
            const double numerical = -(energy_up - energy_down) / (2.0 * step);
            const double analytic = result.forces[i][d];

            report.force_scale = std::max(report.force_scale, std::fabs(analytic));
            report.max_absolute = std::max(report.max_absolute,
                                           std::fabs(analytic - numerical));
        }
    }
    report.max_relative = report.max_absolute / std::max(report.force_scale, 1.0e-30);

    check(report.force_scale > 1.0e-4,
          label + ": every force component is negligible, so this case cannot "
                  "validate a gradient (largest |F| = " + format(report.force_scale) + ")");
    check(report.max_relative < tolerance,
          label + ": F != -dU/dr, worst absolute " + format(report.max_absolute) +
              ", relative " + format(report.max_relative) + ", tolerance " +
              format(tolerance));

    // For the asymmetric fixture, every component of every force must carry
    // signal, or "the gradient is right" would be a claim about a subset of the
    // components. The boundary fixture is not held to this: it is placed to
    // exercise wrapping, and where an atom sits on a face its transverse
    // components are legitimately small.
    if (require_every_component_nonzero) {
        double smallest = 1.0e300;
        for (const auto& force : result.forces) {
            for (const double component : force) {
                smallest = std::min(smallest, std::fabs(component));
            }
        }
        check(smallest > 1.0e-4 * report.force_scale,
              label + ": some force component is essentially zero (" + format(smallest) +
                  " against a scale of " + format(report.force_scale) +
                  "), so the configuration is not fully asymmetric");
    }

    return report;
}

std::vector<GradientReport> reports;

void test_energy_force_consistency() {
    // The reciprocal mesh alone. The real-space erfc term is truncated
    // abruptly at the cutoff, so a pair that crosses it during a finite
    // difference puts a step in the energy; excluding real space entirely
    // removes that from the comparison, and the case below reinstates it with
    // a cutoff no pair sits near.
    for (const int order : {4, 6}) {
        for (const int mesh : {16, 32}) {
            for (const double alpha : {0.25, 0.35, 0.45}) {
                for (const bool neutral : {true, false}) {
                    const auto charges = neutral ? neutral_charges() : net_charges();
                    gmd::PMEForceProvider provider(alpha, kNoRealSpace, order,
                                                   {mesh, mesh, mesh});
                    gmd::RuntimeContext runtime;
                    provider.initialize(runtime);

                    const std::string label =
                        "reciprocal, order " + std::to_string(order) + ", mesh " +
                        std::to_string(mesh) + ", alpha " + format(alpha, 2) +
                        (neutral ? ", neutral" : ", net +0.7");
                    reports.push_back(check_energy_force_consistency(
                        label, provider, asymmetric_positions(), charges, 2.0e-6, true));
                }
            }
        }
    }

    // Atoms on and across the periodic faces and a corner.
    for (const int order : {4, 6}) {
        for (const int mesh : {16, 32}) {
            gmd::PMEForceProvider provider(0.35, kNoRealSpace, order, {mesh, mesh, mesh});
            gmd::RuntimeContext runtime;
            provider.initialize(runtime);
            const std::string label = "boundary atoms, order " + std::to_string(order) +
                                      ", mesh " + std::to_string(mesh);
            reports.push_back(check_energy_force_consistency(
                label, provider, boundary_positions(), neutral_charges(), 2.0e-6, false));
        }
    }
}

// With the real-space term switched on. The cutoff is checked against every
// pair separation first: an abrupt truncation makes the energy discontinuous
// wherever a pair sits at the cutoff, and a finite difference straddling that
// measures the step rather than the gradient.
void test_energy_force_consistency_with_real_space() {
    const double cutoff = 5.5;
    const auto positions = asymmetric_positions();

    double closest_approach = 1.0e300;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        for (std::size_t j = i + 1; j < positions.size(); ++j) {
            double r_squared = 0.0;
            for (std::size_t d = 0; d < 3; ++d) {
                double delta = positions[i][d] - positions[j][d];
                while (delta > 0.5 * kBox[d]) delta -= kBox[d];
                while (delta < -0.5 * kBox[d]) delta += kBox[d];
                r_squared += delta * delta;
            }
            closest_approach = std::min(closest_approach,
                                        std::fabs(std::sqrt(r_squared) - cutoff));
        }
    }
    check(closest_approach > 0.05,
          "real-space fixture: some pair sits within " + format(closest_approach) +
              " of the cutoff, so a finite difference would straddle the truncation "
              "discontinuity rather than measure the gradient");

    for (const int order : {4, 6}) {
        gmd::PMEForceProvider provider(0.35, cutoff, order, {32, 32, 32});
        gmd::RuntimeContext runtime;
        provider.initialize(runtime);
        reports.push_back(check_energy_force_consistency(
            "real space + reciprocal, order " + std::to_string(order),
            provider, positions, neutral_charges(), 2.0e-6, true));
    }
}

// ===========================================================================
// 3. Invariances and parameter stability
// ===========================================================================

// TWO DIFFERENT CLAIMS, kept apart because only one of them is exact.
//
//   Lattice translation is an exact symmetry. Moving every atom by a whole
//   number of box vectors -- leaving the coordinates OUTSIDE the cell, so the
//   provider's own wrapping has to do the work -- describes the identical
//   periodic system. Fractional coordinates land on the same mesh points, and
//   the energy and every force must come back unchanged to round-off.
//
//   Arbitrary translation is NOT an exact symmetry of PME. The exact Ewald sum
//   is translation invariant, but the mesh is not: shifting the system moves
//   the charges relative to the grid and changes the B-spline aliasing error.
//   Asserting exact invariance here would be asserting something false -- it is
//   what a first draft of this test did assert, and PME failed it by 9.7e-5 at
//   order 4, which is simply the size of the mesh error. What can honestly be
//   required is that the deviation is bounded by that error and shrinks when
//   the mesh is refined.
void test_lattice_translation_is_exact() {
    // Different multiples on different axes, and negative ones, so the result
    // is not reachable by any single wrap.
    const std::array<double, 3> lattice_shift = {2.0 * kBox[0], -1.0 * kBox[1],
                                                 3.0 * kBox[2]};

    for (const int order : {4, 6}) {
        gmd::PMEForceProvider provider(0.35, kNoRealSpace, order, {32, 32, 32});
        gmd::RuntimeContext runtime;
        provider.initialize(runtime);

        const auto positions = asymmetric_positions();
        gmd::System base = make_system(positions, neutral_charges());
        const auto reference = evaluate(provider, base);

        auto shifted = positions;
        for (auto& position : shifted) {
            for (std::size_t d = 0; d < 3; ++d) position[d] += lattice_shift[d];
        }
        // Deliberately NOT wrapped back: the coordinates handed to the provider
        // are outside the cell, which is the path this test exists to cover.
        bool outside = false;
        for (const auto& position : shifted) {
            for (std::size_t d = 0; d < 3; ++d) {
                if (position[d] < 0.0 || position[d] >= kBox[d]) outside = true;
            }
        }
        check(outside, "order " + std::to_string(order) +
                           ": the lattice-shifted fixture is still inside the cell");

        gmd::System moved = make_system(shifted, neutral_charges());
        const auto result = evaluate(provider, moved);

        const double energy_error =
            std::fabs(result.potential_energy - reference.potential_energy);
        check(energy_error < 1.0e-11 * std::fabs(reference.potential_energy),
              "order " + std::to_string(order) +
                  ": energy changed under a whole-lattice translation by " +
                  format(energy_error));

        double worst = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < positions.size(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) {
                scale = std::max(scale, std::fabs(reference.forces[i][d]));
                worst = std::max(worst,
                                 std::fabs(result.forces[i][d] - reference.forces[i][d]));
            }
        }
        check(worst < 1.0e-10 * scale,
              "order " + std::to_string(order) +
                  ": forces changed under a whole-lattice translation by " +
                  format(worst) + " (scale " + format(scale) + ")");
        std::cout << "    order " << order << " lattice translation: dE = "
                  << format(energy_error) << ", max dF = " << format(worst) << '\n';
    }
}

void test_arbitrary_translation_is_bounded_by_mesh_error() {
    const std::array<double, 3> shift = {7.3, 17.9, 5.1};

    for (const int order : {4, 6}) {
        std::vector<double> deviations;
        for (const int mesh : {16, 32, 64}) {
            gmd::PMEForceProvider provider(0.35, kNoRealSpace, order, {mesh, mesh, mesh});
            gmd::RuntimeContext runtime;
            provider.initialize(runtime);

            const auto positions = asymmetric_positions();
            gmd::System base = make_system(positions, neutral_charges());
            const auto reference = evaluate(provider, base);

            auto translated = positions;
            bool wrapped = false;
            for (auto& position : translated) {
                for (std::size_t d = 0; d < 3; ++d) {
                    position[d] += shift[d];
                    while (position[d] >= kBox[d]) {
                        position[d] -= kBox[d];
                        wrapped = true;
                    }
                    while (position[d] < 0.0) position[d] += kBox[d];
                }
            }
            check(wrapped, "the arbitrary-translation fixture wraps no atom");

            gmd::System moved = make_system(translated, neutral_charges());
            const auto result = evaluate(provider, moved);

            double worst = 0.0;
            double scale = 0.0;
            for (std::size_t i = 0; i < positions.size(); ++i) {
                for (std::size_t d = 0; d < 3; ++d) {
                    scale = std::max(scale, std::fabs(reference.forces[i][d]));
                    worst = std::max(worst, std::fabs(result.forces[i][d]
                                                      - reference.forces[i][d]));
                }
            }
            deviations.push_back(worst / scale);
            std::cout << "    order " << order << " arbitrary translation, mesh "
                      << std::setw(3) << mesh << ": relative force deviation "
                      << format(worst / scale) << '\n';
        }

        for (std::size_t step = 1; step < deviations.size(); ++step) {
            check(deviations[step] * 2.0 <= deviations[step - 1],
                  "order " + std::to_string(order) +
                      ": the translation deviation does not shrink with mesh "
                      "refinement (" + format(deviations[step - 1]) + " -> " +
                      format(deviations[step]) + "), so it is not mesh error");
        }
    }
}

// Repeated evaluation must be bitwise identical. If alpha or the cutoff were
// being re-resolved from the box, or the influence function rebuilt with
// different parameters, the second evaluation would differ -- and every finite
// difference above would be differencing two slightly different functions.
void test_parameters_are_stable() {
    gmd::PMEForceProvider provider(0.35, kNoRealSpace, 6, {32, 32, 32});
    gmd::RuntimeContext runtime;
    provider.initialize(runtime);

    gmd::System system = make_system(asymmetric_positions(), neutral_charges());
    const auto first = evaluate(provider, system);
    const auto second = evaluate(provider, system);
    const auto third = evaluate(provider, system);

    check(first.potential_energy == second.potential_energy
              && second.potential_energy == third.potential_energy,
          "repeated evaluation is not bitwise reproducible, so a finite difference "
          "would compare two different functions");
    for (std::size_t i = 0; i < first.forces.size(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            check(first.forces[i][d] == third.forces[i][d],
                  "repeated evaluation produced different forces for atom " +
                      std::to_string(i));
        }
    }
}

void print_reports() {
    std::cout << "\n  energy-force consistency: F vs -dU/dr, central differences\n";
    std::cout << "    max |dF|    rel        |F| scale   case\n";
    for (const auto& report : reports) {
        std::cout << "    " << format(report.max_absolute) << "   "
                  << format(report.max_relative) << "   "
                  << format(report.force_scale) << "   " << report.label << '\n';
    }
}

}  // namespace

int main() {
    std::cout << "PME production energy and force path\n";
    test_fft_normalisation_against_direct_dft();
    test_parameters_are_stable();
    test_energy_force_consistency();
    test_energy_force_consistency_with_real_space();
    std::cout << "\n  Invariance\n";
    test_lattice_translation_is_exact();
    test_arbitrary_translation_is_bounded_by_mesh_error();
    print_reports();

    if (failures != 0) {
        std::cerr << "\nPME energy-force tests failed: " << failures << '\n';
        return 1;
    }
    std::cout << "\nPME energy-force tests passed\n";
    return 0;
}
