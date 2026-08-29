#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <span>
#include <string>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/system/special_pair_map.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

namespace {

// CODATA 2022 k_e = E_h * a0, matching gmd::kCoulombConstant. Declared
// independently so this reference is not a restatement of production.
constexpr double kCoulomb = 14.3996454686836;
constexpr double tolerance = 1.0e-10;

void check(bool value, const std::string& message, int& failures) {
    if (!value) {
        std::cerr << "[special pair] " << message << '\n';
        ++failures;
    }
}

bool close(double lhs, double rhs, double tol = tolerance) {
    return std::abs(lhs - rhs) <= tol;
}

gmd::Topology chain_topology() {
    gmd::Topology topology;
    topology.bonds = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}};
    topology.angles = {{0, 1, 2, 0}, {1, 2, 3, 0}};
    topology.dihedrals = {{0, 1, 2, 3, 0}};
    return topology;
}

gmd::System make_chain(const gmd::SpecialPairScaleConfig& scales) {
    gmd::System system;
    system.resize(4, 4);
    gmd::Box box;
    box.set_lengths({30.0, 30.0, 30.0});
    system.set_box(box);
    system.mutable_coordinates()[0] = {5.0, 5.0, 5.0};
    system.mutable_coordinates()[1] = {6.0, 5.0, 5.0};
    system.mutable_coordinates()[2] = {7.0, 5.0, 5.0};
    system.mutable_coordinates()[3] = {8.0, 5.0, 5.0};
    system.mutable_charges()[0] = 1.0;
    system.mutable_charges()[1] = -1.0;
    system.mutable_charges()[2] = 1.0;
    system.mutable_charges()[3] = -1.0;
    system.set_special_pair_map(
        std::make_shared<gmd::SpecialPairMap>(chain_topology(), scales));
    return system;
}

gmd::ForceRequest request_for(const gmd::System& system) {
    const auto coords = system.coordinates();
    return {
        .system = &system,
        .box = &system.box(),
        .coordinates = std::span<const gmd::Coordinate3D>(coords.data(), coords.size()),
    };
}

gmd::ForceResult compute_lj(gmd::System& system, gmd::RuntimeContext& runtime) {
    gmd::ClassicalForceProvider provider(0.2, 1.0, 9.0);
    gmd::ForceResult result;
    provider.compute(request_for(system), result, runtime);
    return result;
}

gmd::ForceResult compute_ewald(gmd::System& system, gmd::RuntimeContext& runtime) {
    gmd::EwaldForceProvider provider(0.3, 7, 9.0);
    gmd::ForceResult result;
    provider.compute(request_for(system), result, runtime);
    return result;
}

gmd::ForceResult compute_pme(gmd::System& system, gmd::RuntimeContext& runtime) {
    gmd::PMEForceProvider provider(0.3, 9.0, 4, {16, 16, 16});
    gmd::ForceResult result;
    provider.compute(request_for(system), result, runtime);
    return result;
}

double lj_energy(double r) {
    constexpr double eps4 = 0.8;
    const double inv_r2 = 1.0 / (r * r);
    const double s6 = inv_r2 * inv_r2 * inv_r2;
    return eps4 * (s6 * s6 - s6);
}

double lj_force_factor(double r) {
    constexpr double eps4 = 0.8;
    const double inv_r2 = 1.0 / (r * r);
    const double s6 = inv_r2 * inv_r2 * inv_r2;
    return eps4 * (12.0 * s6 * s6 - 6.0 * s6) * inv_r2;
}

void test_lj_exclusions_and_scaling(gmd::RuntimeContext& runtime, int& failures) {
    gmd::SpecialPairScaleConfig scales;
    auto chain = make_chain(scales);
    const auto result = compute_lj(chain, runtime);

    const double expected_energy = 0.5 * (lj_energy(3.0) - lj_energy(9.0));
    const double expected_fx_atom0 = 0.5 * lj_force_factor(3.0) * -3.0;
    check(close(result.potential_energy, expected_energy),
          "only scaled 1-4 LJ energy should remain", failures);
    check(close(result.forces[0][0], expected_fx_atom0) &&
              close(result.forces[3][0], -expected_fx_atom0),
          "scaled 1-4 LJ force is incorrect", failures);
    check(close(result.forces[1][0], 0.0) && close(result.forces[2][0], 0.0),
          "excluded 1-2 and 1-3 LJ pairs contributed force", failures);

    scales.pair_14.lj = 0.25;
    auto quarter_chain = make_chain(scales);
    const auto quarter_result = compute_lj(quarter_chain, runtime);
    check(close(quarter_result.potential_energy, result.potential_energy * 0.5),
          "changing lj_scale_14 should scale energy linearly", failures);
    check(close(quarter_result.forces[0][0], result.forces[0][0] * 0.5),
          "changing lj_scale_14 should scale force linearly", failures);
}

template <typename Compute>
void test_coulomb_correction(const char* name,
                             Compute compute,
                             gmd::RuntimeContext& runtime,
                             int& failures) {
    gmd::SpecialPairScaleConfig full;
    full.pair_12.coulomb = 1.0;
    full.pair_13.coulomb = 1.0;
    full.pair_14.coulomb = 1.0;
    auto full_chain = make_chain(full);
    const auto baseline = compute(full_chain, runtime);

    gmd::SpecialPairScaleConfig special;
    auto special_chain = make_chain(special);
    const auto scaled = compute(special_chain, runtime);

    const std::array<double, 4> q = {1.0, -1.0, 1.0, -1.0};
    std::array<double, 4> expected_fx = {0.0, 0.0, 0.0, 0.0};
    double expected_energy_delta = 0.0;
    const std::array<std::array<int, 3>, 6> pairs = {{
        {0, 1, 12}, {1, 2, 12}, {2, 3, 12},
        {0, 2, 13}, {1, 3, 13}, {0, 3, 14},
    }};
    for (const auto& pair : pairs) {
        const int i = pair[0], j = pair[1];
        const double scale = pair[2] == 14 ? special.pair_14.coulomb : 0.0;
        const double delta = scale - 1.0;
        const double r = static_cast<double>(j - i);
        const double prefactor = delta * kCoulomb * q[i] * q[j];
        expected_energy_delta += prefactor / r;
        const double dr = -r;
        expected_fx[static_cast<std::size_t>(i)] += prefactor * dr / (r * r * r);
        expected_fx[static_cast<std::size_t>(j)] -= prefactor * dr / (r * r * r);
    }

    check(close(scaled.potential_energy - baseline.potential_energy,
                expected_energy_delta, 2.0e-9),
          std::string(name) + " special-pair energy correction is incorrect", failures);
    for (std::size_t i = 0; i < expected_fx.size(); ++i) {
        check(close(scaled.forces[i][0] - baseline.forces[i][0],
                    expected_fx[i], 2.0e-9),
              std::string(name) + " special-pair force correction is incorrect", failures);
    }
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif
    gmd::RuntimeContext runtime;
    int failures = 0;
    test_lj_exclusions_and_scaling(runtime, failures);
    test_coulomb_correction("Ewald", compute_ewald, runtime, failures);
    test_coulomb_correction("PME", compute_pme, runtime, failures);
    return failures == 0 ? 0 : 1;
}
