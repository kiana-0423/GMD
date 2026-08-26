// Finite-difference validation of the virial.
//
// Two families of check, both comparing a provider's reported virial against
// the numerical derivative of that same provider's energy.
//
// 1. Isotropic. Under r_i -> s * r_i, L -> s * L,
//
//        dU/ds |_(s=1)  =  -sum_i F_i . r_i  =  -tr(W)
//
//    which is exactly the combination the barostats consume through
//    P = (2*KE + tr(W)) / 3V.
//
// 2. Per-axis (component-wise). Under a normal strain along one axis only --
//    L_a -> (1+e) L_a and r_ia -> (1+e) r_ia, fractional coordinates fixed --
//
//        dU/de |_(e=0)  =  -sum_i F_ia r_ia  =  -W_aa
//
//    so each diagonal component can be pinned down on its own. This catches
//    errors that a trace check cannot: a formula that distributes the right
//    total across the wrong axes passes (1) and fails (2). Non-cubic boxes are
//    used throughout so that no error can hide behind cubic symmetry.
//
// SCOPE OF THE CLAIM. `Box` stores three edge lengths and nothing else, so the
// engine represents orthorhombic cells only and there is no way to apply a
// shear strain. **Only the three diagonal components W_xx, W_yy, W_zz are
// validated here. The off-diagonal components are not validated.** They are
// only checked for the symmetry W_ab == W_ba, which is a necessary but not a
// sufficient condition -- a formula with both off-diagonals wrong in the same
// way would still pass. Validating them needs a triclinic box representation
// and a shear-strain deformation, neither of which exists in this engine.
//
// This matters for the reciprocal-space Ewald/PME term, where the naive
// sum_i r_i (x) F_i is not the virial at all, because the energy also depends
// on the cell explicitly through the 1/V prefactor and through k = 2*pi*n/L.

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/force/bonded_force_provider.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

namespace {

int failures = 0;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[virial fd] " << message << '\n';
        ++failures;
    }
}

// Evaluates a provider against the system as it currently stands.
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

// Returns a copy of `system` with coordinates and box scaled isotropically.
gmd::System scaled_copy(const gmd::System& system, double s) {
    gmd::System copy = system;
    gmd::Box box;
    box.set_lengths({system.box().lengths[0] * s,
                     system.box().lengths[1] * s,
                     system.box().lengths[2] * s});
    copy.set_box(box);
    auto coordinates = copy.mutable_coordinates();
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            coordinates[i][d] *= s;
        }
    }
    return copy;
}

double energy_at_scale(gmd::ForceProvider& provider, const gmd::System& system, double s) {
    gmd::System copy = scaled_copy(system, s);
    return evaluate(provider, copy).potential_energy;
}

// Returns a copy of `system` strained along one axis only: that edge length and
// the matching coordinate component are multiplied by (1 + e), the other two
// axes are untouched. Fractional coordinates are therefore unchanged, which is
// what makes dU/de the corresponding diagonal virial component.
gmd::System strained_copy(const gmd::System& system, std::size_t axis, double e) {
    gmd::System copy = system;
    auto lengths = system.box().lengths;
    lengths[axis] *= (1.0 + e);
    gmd::Box box;
    box.set_lengths(lengths);
    copy.set_box(box);

    auto coordinates = copy.mutable_coordinates();
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        coordinates[i][axis] *= (1.0 + e);
    }
    return copy;
}

double energy_at_strain(gmd::ForceProvider& provider,
                        const gmd::System& system,
                        std::size_t axis,
                        double e) {
    gmd::System copy = strained_copy(system, axis, e);
    return evaluate(provider, copy).potential_energy;
}

// Central difference of U with respect to the isotropic scale factor.
double numerical_dU_ds(gmd::ForceProvider& provider,
                       const gmd::System& system,
                       double eps) {
    const double up = energy_at_scale(provider, system, 1.0 + eps);
    const double down = energy_at_scale(provider, system, 1.0 - eps);
    return (up - down) / (2.0 * eps);
}

double trace(const std::array<double, 9>& virial) {
    return virial[0] + virial[4] + virial[8];
}

// The core assertion: tr(W) must equal -dU/ds.
//
// `tolerance` is relative to the larger of the two magnitudes, because the
// absolute scale varies by orders of magnitude between LJ, bonded and Coulomb
// terms. `eps` trades truncation error against round-off.
void check_virial_matches_finite_difference(gmd::ForceProvider& provider,
                                            gmd::System& system,
                                            const std::string& label,
                                            double tolerance = 1.0e-6,
                                            double eps = 1.0e-6) {
    const gmd::ForceResult result = evaluate(provider, system);
    check(result.success, label + ": force evaluation must succeed");
    check(result.virial_valid, label + ": virial must be reported valid");
    if (!result.virial_valid) {
        return;
    }

    const double analytic = trace(result.virial);
    const double numerical = -numerical_dU_ds(provider, system, eps);

    const double scale = std::max({std::abs(analytic), std::abs(numerical), 1.0e-12});
    const double relative_error = std::abs(analytic - numerical) / scale;

    check(relative_error < tolerance,
          label + ": virial trace does not match -dU/ds (analytic " +
              std::to_string(analytic) + ", finite difference " +
              std::to_string(numerical) + ", relative error " +
              std::to_string(relative_error) + ")");
}

// Validates the three diagonal components independently, each against its own
// per-axis finite difference. Index ordering is asserted explicitly: the tensor
// is row-major, so W_aa lives at virial[a*3 + a].
void check_virial_diagonal(gmd::ForceProvider& provider,
                           gmd::System& system,
                           const std::string& label,
                           double tolerance = 1.0e-5,
                           double eps = 1.0e-6) {
    const gmd::ForceResult result = evaluate(provider, system);
    check(result.success, label + ": force evaluation must succeed");
    check(result.virial_valid, label + ": virial must be reported valid");
    if (!result.virial_valid) {
        return;
    }

    static const char* axis_name[3] = {"xx", "yy", "zz"};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double analytic = result.virial[axis * 3 + axis];

        const double up = energy_at_strain(provider, system, axis, eps);
        const double down = energy_at_strain(provider, system, axis, -eps);
        const double numerical = -(up - down) / (2.0 * eps);

        const double scale =
            std::max({std::abs(analytic), std::abs(numerical), 1.0e-12});
        const double relative_error = std::abs(analytic - numerical) / scale;

        check(relative_error < tolerance,
              label + ": W_" + axis_name[axis] + " does not match -dU/de (analytic " +
                  std::to_string(analytic) + ", finite difference " +
                  std::to_string(numerical) + ", relative error " +
                  std::to_string(relative_error) + ")");

        // Sign must agree too. A component that matched in magnitude but not in
        // sign would still be reported above only through the magnitude of the
        // difference, so state it separately.
        if (std::abs(numerical) > 1.0e-9) {
            check((analytic > 0.0) == (numerical > 0.0),
                  label + ": W_" + axis_name[axis] + " has the wrong sign (analytic " +
                      std::to_string(analytic) + ", finite difference " +
                      std::to_string(numerical) + ")");
        }
    }

    // The diagonal must also add up to the trace the isotropic check uses.
    const double diagonal_sum =
        result.virial[0] + result.virial[4] + result.virial[8];
    check(std::abs(diagonal_sum - trace(result.virial)) < 1.0e-12,
          label + ": diagonal components must sum to the trace");

    // Off-diagonal components are NOT validated here -- see the file header.
    // Symmetry is a necessary condition and cheap to assert, so check it, but
    // it does not establish that the magnitudes are right.
    const double off_scale = std::max({std::abs(result.virial[1]),
                                       std::abs(result.virial[2]),
                                       std::abs(result.virial[5]), 1.0e-12});
    check(std::abs(result.virial[1] - result.virial[3]) < 1.0e-8 * off_scale ||
              std::abs(result.virial[1] - result.virial[3]) < 1.0e-12,
          label + ": virial must be symmetric in xy/yx");
    check(std::abs(result.virial[2] - result.virial[6]) < 1.0e-8 * off_scale ||
              std::abs(result.virial[2] - result.virial[6]) < 1.0e-12,
          label + ": virial must be symmetric in xz/zx");
    check(std::abs(result.virial[5] - result.virial[7]) < 1.0e-8 * off_scale ||
              std::abs(result.virial[5] - result.virial[7]) < 1.0e-12,
          label + ": virial must be symmetric in yz/zy");
}

// --- Systems --------------------------------------------------------------

constexpr double kBoxLength = 16.0;

// A neutral, irregular cluster of charges well inside the box.
gmd::System make_charged_cluster() {
    gmd::System system;
    const std::size_t n = 8;
    system.resize(n, n);
    gmd::Box box;
    box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    system.set_box(box);

    const double positions[8][3] = {
        {4.0, 4.0, 4.0}, {7.1, 4.2, 4.0}, {4.0, 7.3, 4.1}, {7.2, 7.0, 4.0},
        {4.1, 4.0, 7.2}, {7.0, 4.1, 7.0}, {4.0, 7.1, 7.3}, {7.3, 7.2, 7.1},
    };
    for (std::size_t i = 0; i < n; ++i) {
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {positions[i][0], positions[i][1], positions[i][2]};
        // Alternating charges keep the cell neutral.
        system.mutable_charges()[i] = (i % 2 == 0) ? 0.6 : -0.6;
    }
    return system;
}

// A four-atom chain carrying a bond, an angle and a dihedral.
gmd::System make_bonded_chain() {
    gmd::System system;
    system.resize(4, 4);
    gmd::Box box;
    box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    system.set_box(box);
    system.mutable_coordinates()[0] = {6.0, 6.0, 6.0};
    system.mutable_coordinates()[1] = {7.5, 6.1, 6.0};
    system.mutable_coordinates()[2] = {8.2, 7.4, 6.3};
    system.mutable_coordinates()[3] = {9.6, 7.2, 7.1};
    for (std::size_t i = 0; i < 4; ++i) {
        system.mutable_masses()[i] = 12.0;
    }
    return system;
}

std::shared_ptr<gmd::Topology> make_chain_topology() {
    auto topology = std::make_shared<gmd::Topology>();
    topology->bonds = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}};
    topology->angles = {{0, 1, 2, 0}, {1, 2, 3, 0}};
    topology->dihedrals = {{0, 1, 2, 3, 0}};
    return topology;
}

// --- Tests ----------------------------------------------------------------

void test_lennard_jones() {
    gmd::System system = make_charged_cluster();
    // Pair virial: this provider already accumulates r_ij (x) F_ij per pair, so
    // it serves as the positive control for the whole harness.
    gmd::ClassicalForceProvider provider(0.25, 2.5, 6.0);
    check_virial_matches_finite_difference(provider, system, "lennard-jones");
}

void test_bonded() {
    gmd::System system = make_bonded_chain();
    auto topology = make_chain_topology();
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_bond_type({300.0, 1.5});
    provider->add_angle_type({60.0, 112.0 * gmd::deg_to_rad});
    provider->add_dihedral_type({1.5, 3, 0.0});
    check_virial_matches_finite_difference(*provider, system, "bonded");
}

// The same chain, translated so that it straddles the x boundary. The bonded
// energy is unchanged (it is built from minimum-image separations), so the
// virial must be unchanged too. An implementation that forms the virial from
// absolute wrapped coordinates instead of interaction vectors gets this wrong.
gmd::System make_wrapped_bonded_chain() {
    gmd::System system = make_bonded_chain();
    auto coordinates = system.mutable_coordinates();
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        coordinates[i][0] -= 7.0;   // moves atoms 0-1 below zero ...
        if (coordinates[i][0] < 0.0) {
            coordinates[i][0] += kBoxLength;   // ... and wraps them to the far side
        }
    }
    return system;
}

void test_bonded_across_periodic_boundary() {
    gmd::System system = make_wrapped_bonded_chain();
    auto topology = make_chain_topology();
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_bond_type({300.0, 1.5});
    provider->add_angle_type({60.0, 112.0 * gmd::deg_to_rad});
    provider->add_dihedral_type({1.5, 3, 0.0});
    check_virial_matches_finite_difference(*provider, system,
                                           "bonded across periodic boundary");

    // Belt and braces: the wrapped molecule must report the same virial as the
    // intact one, since wrapping is a pure relabelling of the same geometry.
    gmd::System intact = make_bonded_chain();
    const double wrapped_trace = trace(evaluate(*provider, system).virial);
    const double intact_trace = trace(evaluate(*provider, intact).virial);
    check(std::abs(wrapped_trace - intact_trace) < 1.0e-9,
          "bonded virial changed when the molecule was wrapped across the "
          "boundary (intact " + std::to_string(intact_trace) + ", wrapped " +
              std::to_string(wrapped_trace) + "); the virial must not depend on "
              "how atoms are wrapped");
}

void test_ewald() {
    gmd::System system = make_charged_cluster();
    // Explicit alpha / kmax / cutoff so that resolve_params() does not re-derive
    // them from the box: the splitting parameter must stay fixed while the cell
    // is scaled, otherwise the finite difference measures a different sum.
    gmd::EwaldForceProvider provider(0.35, 6, 7.0);
    check_virial_matches_finite_difference(provider, system, "ewald");
}

void test_pme() {
    gmd::System system = make_charged_cluster();
    gmd::PMEForceProvider provider(0.35, 7.0, 6, {32, 32, 32});
    check_virial_matches_finite_difference(provider, system, "pme");
}

// ---------------------------------------------------------------------------
// Component-wise checks, on deliberately non-cubic cells.
// ---------------------------------------------------------------------------

// Edge lengths are all different so that an error which happens to be
// symmetric under x/y/z exchange cannot cancel. min(L)/2 = 7.0 stays above
// every cutoff used below.
constexpr double kLx = 14.0;
constexpr double kLy = 17.0;
constexpr double kLz = 20.0;

gmd::System make_orthorhombic_cluster(bool net_charged) {
    gmd::System system;
    const std::size_t n = 8;
    system.resize(n, n);
    gmd::Box box;
    box.set_lengths({kLx, kLy, kLz});
    system.set_box(box);

    // Positions are asymmetric between axes on purpose.
    const double positions[8][3] = {
        {3.0, 4.0, 5.0},  {6.2, 4.4, 5.1},  {3.1, 8.3, 5.0},  {6.4, 8.1, 5.2},
        {3.0, 4.1, 9.3},  {6.1, 4.2, 9.0},  {3.2, 8.0, 9.4},  {6.3, 8.2, 9.1},
    };
    for (std::size_t i = 0; i < n; ++i) {
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {positions[i][0], positions[i][1], positions[i][2]};
        system.mutable_charges()[i] = (i % 2 == 0) ? 0.6 : -0.6;
    }
    if (net_charged) {
        // Break neutrality so the 1/V net-charge correction is exercised. That
        // term contributes an isotropic U_net * delta_ab, which a per-axis
        // strain probes one component at a time.
        system.mutable_charges()[0] = 1.4;
    }
    return system;
}

// The bonded chain from above, translated so that it straddles the x boundary.
gmd::System make_orthorhombic_bonded_chain(bool wrapped) {
    gmd::System system;
    system.resize(4, 4);
    gmd::Box box;
    box.set_lengths({kLx, kLy, kLz});
    system.set_box(box);
    system.mutable_coordinates()[0] = {6.0, 6.0, 6.0};
    system.mutable_coordinates()[1] = {7.5, 6.1, 6.4};
    system.mutable_coordinates()[2] = {8.2, 7.4, 7.3};
    system.mutable_coordinates()[3] = {9.6, 7.2, 8.1};
    for (std::size_t i = 0; i < 4; ++i) {
        system.mutable_masses()[i] = 12.0;
    }
    if (wrapped) {
        auto coordinates = system.mutable_coordinates();
        for (std::size_t i = 0; i < coordinates.size(); ++i) {
            coordinates[i][0] -= 7.0;
            if (coordinates[i][0] < 0.0) {
                coordinates[i][0] += kLx;
            }
        }
    }
    return system;
}

std::shared_ptr<gmd::BondedForceProvider> make_bonded_provider() {
    auto provider = std::make_shared<gmd::BondedForceProvider>(make_chain_topology());
    provider->add_bond_type({300.0, 1.5});
    provider->add_angle_type({60.0, 112.0 * gmd::deg_to_rad});
    provider->add_dihedral_type({1.5, 3, 0.0});
    return provider;
}

void test_components_lennard_jones() {
    gmd::System system = make_orthorhombic_cluster(false);
    gmd::ClassicalForceProvider provider(0.25, 2.5, 6.0);
    check_virial_diagonal(provider, system, "lennard-jones (orthorhombic)");
}

void test_components_bonded() {
    gmd::System intact = make_orthorhombic_bonded_chain(false);
    auto provider = make_bonded_provider();
    check_virial_diagonal(*provider, intact, "bonded (orthorhombic)");

    // The case that motivated the per-interaction virial: a molecule split
    // across the x boundary. Only W_xx is affected by the wrapping, so a
    // per-axis check localises the error where a trace check would smear it.
    gmd::System wrapped = make_orthorhombic_bonded_chain(true);
    check_virial_diagonal(*provider, wrapped, "bonded across x boundary");

    // Wrapping is a relabelling, so every component must be unchanged by it.
    const gmd::ForceResult a = evaluate(*provider, intact);
    const gmd::ForceResult b = evaluate(*provider, wrapped);
    static const char* axis_name[3] = {"xx", "yy", "zz"};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        check(std::abs(a.virial[axis * 3 + axis] - b.virial[axis * 3 + axis]) < 1.0e-9,
              std::string("bonded W_") + axis_name[axis] +
                  " changed when the molecule was wrapped across the boundary");
    }
}

void test_components_ewald() {
    gmd::System neutral = make_orthorhombic_cluster(false);
    gmd::EwaldForceProvider neutral_provider(0.35, 6, 6.0);
    check_virial_diagonal(neutral_provider, neutral, "ewald (neutral, orthorhombic)");

    gmd::System charged = make_orthorhombic_cluster(true);
    gmd::EwaldForceProvider charged_provider(0.35, 6, 6.0);
    check_virial_diagonal(charged_provider, charged, "ewald (net charged)");
}

void test_components_pme() {
    gmd::System neutral = make_orthorhombic_cluster(false);
    gmd::PMEForceProvider neutral_provider(0.35, 6.0, 6, {32, 32, 32});
    check_virial_diagonal(neutral_provider, neutral, "pme (neutral, orthorhombic)");

    gmd::System charged = make_orthorhombic_cluster(true);
    gmd::PMEForceProvider charged_provider(0.35, 6.0, 6, {32, 32, 32});
    check_virial_diagonal(charged_provider, charged, "pme (net charged)");
}

// A cubic box can hide an axis mix-up; assert the fixtures really are not cubic
// so the component checks above keep their power.
void test_fixtures_are_non_cubic() {
    gmd::System system = make_orthorhombic_cluster(false);
    const auto& lengths = system.box().lengths;
    check(lengths[0] != lengths[1] && lengths[1] != lengths[2] &&
              lengths[0] != lengths[2],
          "component-wise fixtures must use three distinct edge lengths");
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    // An MPI-enabled build still runs this as a single-process test, but the
    // force providers issue collectives that require an initialised MPI.
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif
    test_lennard_jones();
    test_bonded();
    test_bonded_across_periodic_boundary();
    test_ewald();
    test_pme();

    test_fixtures_are_non_cubic();
    test_components_lennard_jones();
    test_components_bonded();
    test_components_ewald();
    test_components_pme();

    if (failures != 0) {
        std::cerr << "[virial fd] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[virial fd] all checks passed\n";
    return 0;
}
