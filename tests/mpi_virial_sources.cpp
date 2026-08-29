// MPI validation of every virial source, all nine components, at np = 1, 2, 4.
//
// WHAT THIS TEST COMPARES AGAINST. Not a serial run of the same code. Every
// provider is held against the independent references in
// tests/virial_reference.hpp -- analytical pair sums, an independently derived
// reciprocal tensor, and a force moment built from finite differences of an
// independently written bonded energy. Those references are pure test-side
// arithmetic over the whole fixture: they contain no MPI, they do not change
// with the rank count, and they are the same numbers the serial test
// tests/virial_source_inventory_tests.cpp checks against. A provider that
// agrees with them at np = 1, 2 and 4 therefore agrees with the serial result
// and with the physics, and the two claims are not entangled.
//
// THE COUNTING FAILURES THIS IS BUILT TO CATCH:
//
//   * a tensor reduced twice, or reduced when it was already global, comes out
//     `size` times too large -- caught by the reference comparison, and stated
//     separately by an explicit per-component guard against N x reference;
//   * a tensor never reduced comes out roughly 1/size too small;
//   * a local/ghost pair counted on both of the ranks that can see it comes
//     out too large by less than a clean factor, which is why the tolerance
//     here is tight rather than generous;
//   * a special-pair correction, or a bonded term spanning a rank boundary,
//     applied by every rank that owns any of its atoms rather than by one.
//
// Every component also carries a non-zero guard: a fixture whose W_xy is
// negligible cannot detect anything about W_xy, and a blank tensor must not be
// able to pass by agreeing with a blank reference.

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include <mpi.h>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/bonded_force_provider.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/composite_force_provider.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/special_pair_map.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

#include "virial_reference.hpp"

namespace vr = virial_ref;

namespace {

int failures = 0;
int global_rank = 0;
int global_size = 1;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[mpi virial sources][rank " << global_rank << "] " << message
                  << '\n';
        ++failures;
    }
}

const char* component_name(std::size_t index) {
    static const char* names[9] = {"xx", "xy", "xz", "yx", "yy", "yz",
                                   "zx", "zy", "zz"};
    return names[index];
}

// --- the fixture ----------------------------------------------------------

constexpr std::array<double, 3> kBox = {18.0, 21.0, 15.0};
constexpr double kCutoff = 6.0;
constexpr double kGhostWidth = 7.0;   // > cutoff, so every in-range pair is visible
constexpr double kAlpha = 0.32;
// The k-space truncation the Ewald PROVIDER is run with, and therefore also
// the truncation its own reference must use, so the two agree exactly.
constexpr int kKmax = 6;
// A converged truncation, used only where PME is the thing under test. PME
// sums the mesh out to |n| = K/2, so comparing it against a kmax = 6 Ewald sum
// would measure Ewald's truncation error (~4e-4 relative here) rather than
// PME's mesh error. exp(-k^2 / 4 alpha^2) is below 1e-17 at this kmax for
// every axis of the fixture box.
constexpr int kConvergedKmax = 16;
constexpr double kNoRealSpace = 0.5;

// Twelve charges spread over the whole cell so that every rank count produces
// a genuinely different partition, and tilted so that no tensor component is
// small. Several pairs are within the cutoff across the periodic faces.
std::vector<vr::Vec3> fixture_positions() {
    return {{2.31L, 3.17L, 1.94L},   {6.42L, 4.83L, 5.11L},
            {9.76L, 11.28L, 2.63L},  {3.19L, 8.94L, 7.42L},
            {14.53L, 2.71L, 4.38L},  {5.87L, 16.06L, 8.19L},
            {11.14L, 6.35L, 12.73L}, {1.62L, 14.41L, 3.55L},
            {16.05L, 9.62L, 10.17L}, {8.43L, 18.77L, 13.26L},
            {13.28L, 13.44L, 6.81L}, {4.76L, 10.19L, 11.53L}};
}

std::vector<double> fixture_charges() {
    return {0.7, -0.5, 0.9, -0.4, 0.3, -1.0, 0.6, -0.6, 0.5, -0.35, 0.45, -0.6};
}

std::array<double, 3> as_double(const vr::Vec3& v) {
    return {static_cast<double>(v[0]), static_cast<double>(v[1]),
            static_cast<double>(v[2])};
}

// This rank's share, owned atoms only. Ghosts are added afterwards by
// exchange_ghost_coordinates(), which is the path a real run takes.
gmd::System make_local_system(const gmd::DomainDecomposition& decomposition,
                              bool with_charges) {
    const auto positions = fixture_positions();
    const auto charges = fixture_charges();

    gmd::Box box;
    box.set_lengths(kBox);

    std::vector<int> owned;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (decomposition.owner_rank(box, as_double(positions[i])) == global_rank) {
            owned.push_back(static_cast<int>(i));
        }
    }

    gmd::System system;
    system.resize(owned.size(), owned.size());
    system.set_box(box);
    for (std::size_t local = 0; local < owned.size(); ++local) {
        const auto tag = static_cast<std::size_t>(owned[local]);
        system.mutable_masses()[local] = 12.0;
        system.mutable_atom_types()[local] = 0;
        system.mutable_atom_tags()[local] = owned[local];
        system.mutable_atom_owners()[local] = global_rank;
        system.mutable_coordinates()[local] = as_double(positions[tag]);
        if (with_charges) system.mutable_charges()[local] = charges[tag];
    }
    return system;
}

gmd::DomainDecomposition make_decomposition() {
    gmd::Box box;
    box.set_lengths(kBox);
    gmd::DomainDecomposition decomposition;
    decomposition.create_decomposition(box, global_size, global_rank, kGhostWidth, 0.0,
                                       {true, true, true});
    return decomposition;
}

std::array<double, 9> evaluate_with_ghosts(gmd::ForceProvider& provider,
                                           gmd::System& system,
                                           const gmd::DomainDecomposition& decomposition,
                                           const gmd::MpiCommunicator& communicator,
                                           gmd::RuntimeContext& runtime) {
    communicator.exchange_ghost_coordinates(system, decomposition);
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
    check(result.virial_valid, "provider reported an invalid virial under MPI");
    return result.virial;
}

// --- shared assertions ----------------------------------------------------

double tensor_scale(const std::array<double, 9>& tensor) {
    double scale = 0.0;
    for (const double value : tensor) scale = std::max(scale, std::fabs(value));
    return scale;
}

// Bitwise agreement across the communicator. Comparing min against max is
// stronger than comparing every rank to rank 0 and needs one collective.
void check_identical_on_every_rank(const std::array<double, 9>& tensor,
                                   const std::string& label) {
    std::array<double, 9> minimum{};
    std::array<double, 9> maximum{};
    MPI_Allreduce(tensor.data(), minimum.data(), 9, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(tensor.data(), maximum.data(), 9, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (std::size_t i = 0; i < 9; ++i) {
        check(minimum[i] == maximum[i],
              label + ": W_" + component_name(i) +
                  " is not bitwise identical on every rank (min " +
                  std::to_string(minimum[i]) + ", max " + std::to_string(maximum[i]) + ")");
    }
}

// The core comparison: all nine components against an independent reference,
// with a non-zero guard and an explicit statement that the value is not the
// reference multiplied by the rank count.
void check_against_reference(const std::array<double, 9>& measured,
                             const vr::Mat3& reference,
                             const std::string& label,
                             double tolerance) {
    double scale = 0.0;
    for (const auto value : reference) {
        scale = std::max(scale, std::fabs(static_cast<double>(value)));
    }
    check(scale > 0.0, label + ": the reference tensor is zero");
    if (scale <= 0.0) return;

    for (std::size_t i = 0; i < 9; ++i) {
        const double expected = static_cast<double>(reference[i]);
        const double error = std::fabs(measured[i] - expected) / scale;
        check(error < tolerance,
              label + ": W_" + component_name(i) + " = " + std::to_string(measured[i]) +
                  ", expected " + std::to_string(expected) + " (relative error " +
                  std::to_string(error) + ")");
        check(std::fabs(expected) > 1.0e-3 * scale,
              label + ": reference W_" + component_name(i) +
                  " is too small for this fixture to test it");
    }

    // Say the rank-multiplication failure out loud rather than relying on the
    // tolerance above to imply it. At np = 1 the two coincide, so the guard is
    // only meaningful for size > 1 and is skipped otherwise instead of being
    // silently trivial.
    if (global_size > 1) {
        for (std::size_t i = 0; i < 9; ++i) {
            const double multiplied = static_cast<double>(reference[i]) * global_size;
            check(std::fabs(measured[i] - multiplied) > 0.5 * std::fabs(multiplied - static_cast<double>(reference[i])),
                  label + ": W_" + component_name(i) +
                      " is consistent with the reference multiplied by the rank count (" +
                      std::to_string(measured[i]) + " vs " + std::to_string(multiplied) + ")");
        }
    }
}

// --- references -----------------------------------------------------------

vr::Vec3 minimum_image(vr::Vec3 delta) {
    for (std::size_t d = 0; d < 3; ++d) {
        const auto length = static_cast<vr::Real>(kBox[d]);
        while (delta[d] > 0.5L * length) delta[d] -= length;
        while (delta[d] < -0.5L * length) delta[d] += length;
    }
    return delta;
}

vr::Mat3 lennard_jones_reference(double epsilon, double sigma, double cutoff) {
    const auto positions = fixture_positions();
    vr::Mat3 virial{};
    for (std::size_t i = 0; i < positions.size(); ++i) {
        for (std::size_t j = i + 1; j < positions.size(); ++j) {
            const vr::Vec3 delta = minimum_image(vr::sub(positions[i], positions[j]));
            const vr::Real r_squared = vr::dot(delta, delta);
            if (r_squared >= static_cast<vr::Real>(cutoff) * static_cast<vr::Real>(cutoff)) {
                continue;
            }
            const vr::Real ratio_squared =
                static_cast<vr::Real>(sigma) * static_cast<vr::Real>(sigma) / r_squared;
            const vr::Real ratio_six = ratio_squared * ratio_squared * ratio_squared;
            const vr::Real ratio_twelve = ratio_six * ratio_six;
            const vr::Real factor = 24.0L * static_cast<vr::Real>(epsilon)
                                  * (2.0L * ratio_twelve - ratio_six) / r_squared;
            const vr::Vec3 force = {factor * delta[0], factor * delta[1], factor * delta[2]};
            vr::accumulate_outer(virial, delta, force);
        }
    }
    return virial;
}

vr::Mat3 ewald_real_space_reference(double cutoff) {
    const auto positions = fixture_positions();
    const auto charges = fixture_charges();
    const vr::Real two_alpha_over_root_pi =
        2.0L * static_cast<vr::Real>(kAlpha) / std::sqrt(std::numbers::pi_v<vr::Real>);

    vr::Mat3 virial{};
    for (std::size_t i = 0; i < positions.size(); ++i) {
        for (std::size_t j = i + 1; j < positions.size(); ++j) {
            const vr::Vec3 delta = minimum_image(vr::sub(positions[i], positions[j]));
            const vr::Real r_squared = vr::dot(delta, delta);
            if (r_squared >= static_cast<vr::Real>(cutoff) * static_cast<vr::Real>(cutoff)) {
                continue;
            }
            const vr::Real r = std::sqrt(r_squared);
            const vr::Real ar = static_cast<vr::Real>(kAlpha) * r;
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

vr::CellConfiguration reference_configuration() {
    std::vector<vr::Real> charges;
    for (const double q : fixture_charges()) charges.push_back(static_cast<vr::Real>(q));
    return vr::make_orthorhombic(
        {static_cast<vr::Real>(kBox[0]), static_cast<vr::Real>(kBox[1]),
         static_cast<vr::Real>(kBox[2])},
        fixture_positions(), charges);
}

// ===========================================================================
// 1. Lennard-Jones pair virial
// ===========================================================================

void test_lennard_jones(const gmd::MpiCommunicator& communicator,
                        gmd::RuntimeContext& runtime) {
    const double epsilon = 0.0104;
    const double sigma = 3.4;

    const auto decomposition = make_decomposition();
    gmd::System system = make_local_system(decomposition, false);
    gmd::ClassicalForceProvider provider(epsilon, sigma, kCutoff);

    const auto measured =
        evaluate_with_ghosts(provider, system, decomposition, communicator, runtime);

    check_identical_on_every_rank(measured, "LJ");
    check_against_reference(measured, lennard_jones_reference(epsilon, sigma, kCutoff),
                            "LJ pair tensor", 1.0e-11);
}

// ===========================================================================
// 2. Ewald real space, reciprocal space and the net-charge term
// ===========================================================================

void test_ewald(const gmd::MpiCommunicator& communicator,
                gmd::RuntimeContext& runtime) {
    const auto decomposition = make_decomposition();

    // Real space is isolated by difference against a run whose cutoff excludes
    // every pair, exactly as in the serial test.
    gmd::System full_system = make_local_system(decomposition, true);
    gmd::System bare_system = make_local_system(decomposition, true);
    gmd::EwaldForceProvider with_pairs(kAlpha, kKmax, kCutoff);
    gmd::EwaldForceProvider without_pairs(kAlpha, kKmax, kNoRealSpace);

    const auto full =
        evaluate_with_ghosts(with_pairs, full_system, decomposition, communicator, runtime);
    const auto reciprocal =
        evaluate_with_ghosts(without_pairs, bare_system, decomposition, communicator, runtime);

    check_identical_on_every_rank(full, "Ewald total");
    check_identical_on_every_rank(reciprocal, "Ewald reciprocal");

    std::array<double, 9> real_space{};
    for (std::size_t i = 0; i < 9; ++i) real_space[i] = full[i] - reciprocal[i];
    check_against_reference(real_space, ewald_real_space_reference(kCutoff),
                            "Ewald real-space tensor", 1.0e-10);

    // The fixture carries net charge, so the reciprocal-only result is the
    // reciprocal tensor plus the isotropic net-charge correction.
    const auto configuration = reference_configuration();
    std::vector<vr::Real> charges;
    for (const double q : fixture_charges()) charges.push_back(static_cast<vr::Real>(q));

    vr::Mat3 expected =
        vr::ewald_reciprocal_virial_analytic(configuration, kAlpha, kKmax);
    const vr::Mat3 net_charge =
        vr::ewald_net_charge_virial(charges, configuration.volume(), kAlpha);
    for (std::size_t i = 0; i < 9; ++i) expected[i] += net_charge[i];

    check_against_reference(reciprocal, expected,
                            "Ewald reciprocal tensor", 1.0e-10);
}

// ===========================================================================
// 3. PME
// ===========================================================================

// PME replicates the mesh on every rank and stores one equal share of the
// reciprocal result, so the sum over ranks reconstructs it exactly once. That
// is a different counting scheme from the pair terms, and it has its own way
// of going wrong, so it gets its own check.
void test_pme(const gmd::MpiCommunicator& communicator,
              gmd::RuntimeContext& runtime) {
    const auto decomposition = make_decomposition();
    gmd::System system = make_local_system(decomposition, true);

    gmd::PMEForceProvider provider(kAlpha, kNoRealSpace, 6, {64, 64, 64});
    provider.initialize(runtime);
    const auto measured =
        evaluate_with_ghosts(provider, system, decomposition, communicator, runtime);

    check_identical_on_every_rank(measured, "PME reciprocal");

    // Compared against the exact Ewald reciprocal tensor, so the tolerance is
    // the PME mesh error at this grid and order (~3e-7 here), not machine
    // precision. The serial test pins that error down properly by refining the
    // mesh; here the point is only that the rank count does not change it, and
    // every counting failure this test exists to catch is off by tens of
    // percent or more -- orders of magnitude clear of the mesh error.
    const auto configuration = reference_configuration();
    std::vector<vr::Real> charges;
    for (const double q : fixture_charges()) charges.push_back(static_cast<vr::Real>(q));

    vr::Mat3 expected =
        vr::ewald_reciprocal_virial_analytic(configuration, kAlpha, kConvergedKmax);
    const vr::Mat3 net_charge =
        vr::ewald_net_charge_virial(charges, configuration.volume(), kAlpha);
    for (std::size_t i = 0; i < 9; ++i) expected[i] += net_charge[i];

    check_against_reference(measured, expected, "PME reciprocal tensor", 1.0e-5);
}

// ===========================================================================
// 4. Special-pair Coulomb corrections across rank boundaries
// ===========================================================================

std::shared_ptr<gmd::Topology> chain_topology() {
    auto topology = std::make_shared<gmd::Topology>();
    // A chain over atoms whose positions put them in different domains at
    // np = 2 and np = 4, so the corrections and the bonded terms below have to
    // survive a rank boundary.
    //
    // The last atom is 11, the highest tag in the fixture, and that is not
    // cosmetic: BondedForceProvider sizes its global tag table from the
    // largest index the topology mentions and rejects any local tag beyond it,
    // so a chain ending at 10 aborts as soon as a rank owns atom 11.
    topology->bonds = {{0, 4, 0}, {4, 8, 0}, {8, 11, 0}};
    topology->angles = {{0, 4, 8, 0}, {4, 8, 11, 0}};
    topology->dihedrals = {{0, 4, 8, 11, 0}};
    return topology;
}

void test_special_pair_corrections(const gmd::MpiCommunicator& communicator,
                                   gmd::RuntimeContext& runtime) {
    auto topology = chain_topology();
    gmd::SpecialPairScaleConfig scales;
    scales.pair_12 = {0.0, 0.0};
    scales.pair_13 = {0.0, 0.0};
    scales.pair_14 = {0.5, 0.8333333333333333};
    auto special_pairs = std::make_shared<gmd::SpecialPairMap>(*topology, scales);

    const auto decomposition = make_decomposition();

    gmd::System without = make_local_system(decomposition, true);
    gmd::System with = make_local_system(decomposition, true);
    with.set_special_pair_map(special_pairs);

    gmd::EwaldForceProvider plain(kAlpha, kKmax, kCutoff);
    gmd::EwaldForceProvider corrected(kAlpha, kKmax, kCutoff);

    const auto baseline =
        evaluate_with_ghosts(plain, without, decomposition, communicator, runtime);
    const auto adjusted =
        evaluate_with_ghosts(corrected, with, decomposition, communicator, runtime);

    std::array<double, 9> correction{};
    for (std::size_t i = 0; i < 9; ++i) correction[i] = adjusted[i] - baseline[i];

    check_identical_on_every_rank(correction, "special-pair correction");

    // Reference: each corrected pair contributes exactly once, whichever rank
    // happens to own either endpoint.
    const auto positions = fixture_positions();
    const auto charges = fixture_charges();
    vr::Mat3 reference{};
    for (const auto& pair : special_pairs->entries()) {
        const vr::Real delta = static_cast<vr::Real>(pair.scale.coulomb) - 1.0L;
        if (delta == 0.0L) continue;
        const auto a = static_cast<std::size_t>(pair.atom_tag_a);
        const auto b = static_cast<std::size_t>(pair.atom_tag_b);
        const vr::Vec3 separation = minimum_image(vr::sub(positions[a], positions[b]));
        const vr::Real r_squared = vr::dot(separation, separation);
        const vr::Real r = std::sqrt(r_squared);
        const vr::Real factor = delta * vr::kCoulomb * static_cast<vr::Real>(charges[a])
                              * static_cast<vr::Real>(charges[b]) / (r_squared * r);
        const vr::Vec3 force = {factor * separation[0], factor * separation[1],
                                factor * separation[2]};
        vr::accumulate_outer(reference, separation, force);
    }

    check_against_reference(correction, reference,
                            "special-pair Coulomb correction", 1.0e-10);
}

// ===========================================================================
// 5. Bonded interactions spanning rank boundaries
// ===========================================================================

// Independently written bonded energy for the chain above, evaluated on the
// unwrapped chain geometry. Used only through finite differences, so it shares
// nothing with the provider beyond the definition of the interaction.
vr::Real chain_energy(const std::vector<vr::Vec3>& unwrapped,
                      double bond_k, double bond_r0,
                      double angle_k, double angle_theta0,
                      double dihedral_k, int periodicity, double phase) {
    auto bond = [&](std::size_t i, std::size_t j) {
        const vr::Real r = vr::norm(vr::sub(unwrapped[j], unwrapped[i]));
        const vr::Real deviation = r - static_cast<vr::Real>(bond_r0);
        return static_cast<vr::Real>(bond_k) * deviation * deviation;
    };
    auto angle = [&](std::size_t i, std::size_t j, std::size_t k) {
        const vr::Vec3 a = vr::sub(unwrapped[i], unwrapped[j]);
        const vr::Vec3 b = vr::sub(unwrapped[k], unwrapped[j]);
        vr::Real cosine = vr::dot(a, b) / (vr::norm(a) * vr::norm(b));
        if (cosine > 1.0L) cosine = 1.0L;
        if (cosine < -1.0L) cosine = -1.0L;
        const vr::Real deviation = std::acos(cosine) - static_cast<vr::Real>(angle_theta0);
        return static_cast<vr::Real>(angle_k) * deviation * deviation;
    };
    auto cross = [](const vr::Vec3& a, const vr::Vec3& b) {
        return vr::Vec3{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0]};
    };
    auto torsion = [&](std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
        const vr::Vec3 b1 = vr::sub(unwrapped[j], unwrapped[i]);
        const vr::Vec3 b2 = vr::sub(unwrapped[k], unwrapped[j]);
        const vr::Vec3 b3 = vr::sub(unwrapped[l], unwrapped[k]);
        const vr::Vec3 m = cross(b1, b2);
        const vr::Vec3 n = cross(b2, b3);
        const vr::Real phi = std::atan2(vr::norm(b2) * vr::dot(b1, n), vr::dot(m, n));
        return static_cast<vr::Real>(dihedral_k)
             * (1.0L + std::cos(static_cast<vr::Real>(periodicity) * phi
                                - static_cast<vr::Real>(phase)));
    };

    return bond(0, 1) + bond(1, 2) + bond(2, 3)
         + angle(0, 1, 2) + angle(1, 2, 3)
         + torsion(0, 1, 2, 3);
}

void test_bonded_across_ranks(const gmd::MpiCommunicator& communicator,
                              gmd::RuntimeContext& runtime) {
    const double bond_k = 3.1;
    const double bond_r0 = 4.5;
    const double angle_k = 2.4;
    const double angle_theta0 = 1.9106;
    const double dihedral_k = 0.62;
    const int periodicity = 3;
    const double phase = 0.7853981633974483;

    auto topology = chain_topology();
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_bond_type({bond_k, bond_r0});
    provider->add_angle_type({angle_k, angle_theta0});
    provider->add_dihedral_type({dihedral_k, periodicity, phase});
    provider->initialize(runtime);

    const auto decomposition = make_decomposition();
    gmd::System system = make_local_system(decomposition, false);
    const auto measured =
        evaluate_with_ghosts(*provider, system, decomposition, communicator, runtime);

    check_identical_on_every_rank(measured, "bonded");

    // Reference: unwrap the chain by walking minimum-image steps, take
    // finite-difference forces of the independent energy, and form the moment
    // sum about the first chain atom.
    const auto positions = fixture_positions();
    const std::array<std::size_t, 4> chain = {0, 4, 8, 11};

    std::vector<vr::Vec3> unwrapped(4);
    unwrapped[0] = positions[chain[0]];
    for (std::size_t index = 1; index < 4; ++index) {
        const vr::Vec3 step =
            minimum_image(vr::sub(positions[chain[index]], positions[chain[index - 1]]));
        unwrapped[index] = {unwrapped[index - 1][0] + step[0],
                            unwrapped[index - 1][1] + step[1],
                            unwrapped[index - 1][2] + step[2]};
    }

    std::vector<vr::Vec3> forces(4, vr::Vec3{0.0L, 0.0L, 0.0L});
    const vr::Real step = 2.0e-6L;
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            std::vector<vr::Vec3> up = unwrapped;
            std::vector<vr::Vec3> down = unwrapped;
            up[i][d] += step;
            down[i][d] -= step;
            const vr::Real derivative =
                (chain_energy(up, bond_k, bond_r0, angle_k, angle_theta0, dihedral_k,
                              periodicity, phase)
                 - chain_energy(down, bond_k, bond_r0, angle_k, angle_theta0, dihedral_k,
                                periodicity, phase))
                / (2.0L * step);
            forces[i][d] = -derivative;
        }
    }

    vr::Mat3 reference{};
    for (std::size_t i = 0; i < 4; ++i) {
        vr::accumulate_outer(reference, vr::sub(unwrapped[i], unwrapped[0]), forces[i]);
    }

    // The tolerance is set by the finite-difference forces, not by the tensor
    // arithmetic. A term counted twice fails by 100 %, which is six orders of
    // magnitude clear of this.
    check_against_reference(measured, reference,
                            "bonded tensor across rank boundaries", 2.0e-6);
}

// ===========================================================================
// 6. CompositeForceProvider must not reduce an already-global tensor
// ===========================================================================

void test_composite_does_not_reduce_again(const gmd::MpiCommunicator& communicator,
                                          gmd::RuntimeContext& runtime) {
    const auto decomposition = make_decomposition();

    gmd::System lj_system = make_local_system(decomposition, true);
    gmd::System ewald_system = make_local_system(decomposition, true);
    gmd::System composite_system = make_local_system(decomposition, true);

    gmd::ClassicalForceProvider lennard_jones(0.0104, 3.4, kCutoff);
    gmd::EwaldForceProvider ewald(kAlpha, kKmax, kCutoff);

    const auto lj_only =
        evaluate_with_ghosts(lennard_jones, lj_system, decomposition, communicator, runtime);
    const auto ewald_only =
        evaluate_with_ghosts(ewald, ewald_system, decomposition, communicator, runtime);

    auto composite = std::make_shared<gmd::CompositeForceProvider>();
    composite->add(std::make_shared<gmd::ClassicalForceProvider>(0.0104, 3.4, kCutoff));
    composite->add(std::make_shared<gmd::EwaldForceProvider>(kAlpha, kKmax, kCutoff));
    const auto combined = evaluate_with_ghosts(*composite, composite_system, decomposition,
                                               communicator, runtime);

    check_identical_on_every_rank(combined, "composite");

    std::array<double, 9> expected{};
    for (std::size_t i = 0; i < 9; ++i) expected[i] = lj_only[i] + ewald_only[i];

    double scale = tensor_scale(expected);
    check(scale > 0.0, "composite fixture produced a zero tensor");
    for (std::size_t i = 0; i < 9; ++i) {
        const double error = std::fabs(combined[i] - expected[i]) / scale;
        check(error < 1.0e-12,
              std::string("composite: W_") + component_name(i) +
                  " is not the sum of its children under MPI (" +
                  std::to_string(combined[i]) + " vs " + std::to_string(expected[i]) +
                  ", relative error " + std::to_string(error) + ")");
        check(std::fabs(expected[i]) > 1.0e-3 * scale,
              std::string("composite: expected W_") + component_name(i) +
                  " is too small for this fixture to test it");
    }

    // Each child already allreduces its own tensor. If the composite reduced
    // the sum again the result would be `size` times too large; say so
    // explicitly rather than leaving it implied by the tolerance.
    if (global_size > 1) {
        for (std::size_t i = 0; i < 9; ++i) {
            check(std::fabs(combined[i] - expected[i] * global_size)
                      > 0.5 * std::fabs(expected[i] * (global_size - 1)),
                  std::string("composite: W_") + component_name(i) +
                      " is consistent with a second reduction over the communicator");
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    gmd::RuntimeContext runtime;
    gmd::MpiCommunicator communicator;

    test_lennard_jones(communicator, runtime);
    test_ewald(communicator, runtime);
    test_pme(communicator, runtime);
    test_special_pair_corrections(communicator, runtime);
    test_bonded_across_ranks(communicator, runtime);
    test_composite_does_not_reduce_again(communicator, runtime);

    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_rank == 0) {
        if (total == 0) {
            std::cout << "[mpi virial sources] all checks passed on " << global_size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi virial sources] " << total << " check(s) failed on "
                      << global_size << " rank(s)\n";
        }
    }
    return total == 0 ? 0 : 1;
}
