#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/special_pair_map.hpp"
#include "gmd/system/system.hpp"

namespace {

// CODATA 2022 k_e = E_h * a0, matching gmd::kCoulombConstant. Declared
// independently so this reference is not a restatement of production.
constexpr double kCoulomb = 14.3996454686836;
constexpr double tolerance = 2.0e-9;
const std::array<double, 4> positions = {4.0, 4.4, 5.0, 5.4};
const std::array<double, 4> charges = {1.0, -1.0, 1.0, -1.0};

bool close(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= tolerance;
}

void check(bool condition, const std::string& message, int rank, int& failures) {
    if (!condition) {
        std::cerr << "[mpi special pair rank " << rank << "] " << message << '\n';
        ++failures;
    }
}

gmd::Topology topology() {
    gmd::Topology top;
    top.bonds = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}};
    top.angles = {{0, 1, 2, 0}, {1, 2, 3, 0}};
    top.dihedrals = {{0, 1, 2, 3, 0}};
    return top;
}

gmd::Box box() {
    gmd::Box value;
    value.set_lengths({10.0, 10.0, 10.0});
    return value;
}

gmd::DomainDecomposition decomposition(const gmd::Box& value,
                                       int rank,
                                       double ghost_width = 2.0) {
    gmd::DomainDecomposition dd;
    dd.create_1d_decomposition(value, 2, rank, ghost_width, 0.0, true);
    return dd;
}

gmd::System local_chain(int rank, const gmd::SpecialPairScaleConfig& scales) {
    gmd::System system;
    system.resize(2, 2);
    system.set_box(box());
    for (int local = 0; local < 2; ++local) {
        const int tag = 2 * rank + local;
        system.mutable_masses()[static_cast<std::size_t>(local)] = 1.0;
        system.mutable_charges()[static_cast<std::size_t>(local)] =
            charges[static_cast<std::size_t>(tag)];
        system.mutable_coordinates()[static_cast<std::size_t>(local)] =
            {positions[static_cast<std::size_t>(tag)], 5.0, 5.0};
        system.mutable_atom_tags()[static_cast<std::size_t>(local)] = tag;
        system.mutable_atom_owners()[static_cast<std::size_t>(local)] = rank;
    }
    system.set_special_pair_map(
        std::make_shared<gmd::SpecialPairMap>(topology(), scales));
    return system;
}

struct GlobalResult {
    double energy = 0.0;
    std::vector<double> forces;
};

template <typename Provider>
GlobalResult evaluate(Provider& provider,
                      gmd::System& system,
                      const gmd::DomainDecomposition& dd,
                      const gmd::MpiCommunicator& comm,
                      gmd::RuntimeContext& runtime) {
    comm.exchange_ghost_coordinates(system, dd);
    const auto coords = system.coordinates();
    const gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .coordinates = std::span<const gmd::Coordinate3D>(coords.data(), coords.size()),
    };
    gmd::ForceResult result;
    provider.compute(request, result, runtime);
    auto system_forces = system.mutable_forces();
    for (std::size_t i = 0; i < result.forces.size(); ++i) {
        system_forces[i] = result.forces[i];
    }
    comm.reverse_accumulate_ghost_forces(system, dd);

    std::vector<double> local(12, 0.0);
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const std::size_t offset = static_cast<std::size_t>(system.atom_tag(i)) * 3;
        local[offset] = system.forces()[i][0];
        local[offset + 1] = system.forces()[i][1];
        local[offset + 2] = system.forces()[i][2];
    }
    GlobalResult global;
    global.energy = comm.allreduce_scalar(result.potential_energy);
    comm.allreduce_vector(local, global.forces);
    return global;
}

void test_lj(const gmd::MpiCommunicator& comm,
             gmd::RuntimeContext& runtime,
             int rank,
             int& failures) {
    gmd::SpecialPairScaleConfig scales;
    auto system = local_chain(rank, scales);
    auto dd = decomposition(system.box(), rank);
    gmd::ClassicalForceProvider provider(0.2, 0.5, 2.0);
    const auto result = evaluate(provider, system, dd, comm, runtime);

    const double r = positions[3] - positions[0];
    const double s2 = 0.25 / (r * r);
    const double s6 = s2 * s2 * s2;
    const double rc_s2 = 0.25 / 4.0;
    const double rc_s6 = rc_s2 * rc_s2 * rc_s2;
    const double expected_energy =
        0.5 * 0.8 * ((s6*s6 - s6) - (rc_s6*rc_s6 - rc_s6));
    const double ff = 0.5 * 0.8 * (12.0*s6*s6 - 6.0*s6) / (r*r);
    check(close(result.energy, expected_energy),
          "cross-rank LJ exclusion/scaling energy differs from hand result", rank, failures);
    check(close(result.forces[0], ff * -r) &&
              close(result.forces[9], ff * r),
          "cross-rank LJ scaled force differs from hand result", rank, failures);
}

void test_ewald_correction(const gmd::MpiCommunicator& comm,
                           gmd::RuntimeContext& runtime,
                           int rank,
                           int& failures) {
    gmd::SpecialPairScaleConfig full;
    full.pair_12.coulomb = 1.0;
    full.pair_13.coulomb = 1.0;
    full.pair_14.coulomb = 1.0;
    auto full_system = local_chain(rank, full);
    // The intentionally narrow halo omits cross-domain topology partners.
    // Coulomb special-pair correction must therefore use global tags/position
    // exchange rather than depend on a ghost being present.
    auto full_dd = decomposition(full_system.box(), rank, 0.05);
    gmd::EwaldForceProvider full_provider(0.3, 7, 2.0);
    const auto baseline = evaluate(full_provider, full_system, full_dd, comm, runtime);

    gmd::SpecialPairScaleConfig scaled_values;
    auto scaled_system = local_chain(rank, scaled_values);
    auto scaled_dd = decomposition(scaled_system.box(), rank, 0.05);
    gmd::EwaldForceProvider scaled_provider(0.3, 7, 2.0);
    const auto scaled = evaluate(scaled_provider, scaled_system, scaled_dd, comm, runtime);

    const std::array<std::array<int, 3>, 6> pairs = {{
        {0, 1, 12}, {1, 2, 12}, {2, 3, 12},
        {0, 2, 13}, {1, 3, 13}, {0, 3, 14},
    }};
    double expected_energy = 0.0;
    std::array<double, 4> expected_force = {0.0, 0.0, 0.0, 0.0};
    for (const auto& pair : pairs) {
        const int i = pair[0], j = pair[1];
        const double scale = pair[2] == 14 ? scaled_values.pair_14.coulomb : 0.0;
        const double r = positions[static_cast<std::size_t>(j)] -
                         positions[static_cast<std::size_t>(i)];
        const double prefactor =
            (scale - 1.0) * kCoulomb *
            charges[static_cast<std::size_t>(i)] *
            charges[static_cast<std::size_t>(j)];
        expected_energy += prefactor / r;
        expected_force[static_cast<std::size_t>(i)] += prefactor * -r / (r*r*r);
        expected_force[static_cast<std::size_t>(j)] -= prefactor * -r / (r*r*r);
    }
    check(close(scaled.energy - baseline.energy, expected_energy),
          "cross-rank Ewald special correction energy is incorrect", rank, failures);
    for (std::size_t tag = 0; tag < expected_force.size(); ++tag) {
        check(close(scaled.forces[tag * 3] - baseline.forces[tag * 3],
                    expected_force[tag]),
              "cross-rank Ewald special correction force is incorrect", rank, failures);
    }
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);
    gmd::MpiCommunicator comm;
    gmd::RuntimeContext runtime;
    const int rank = comm.rank();
    int failures = 0;
    check(comm.size() == 2, "test requires two ranks", rank, failures);
    if (comm.size() == 2) {
        test_lj(comm, runtime, rank, failures);
        test_ewald_correction(comm, runtime, rank, failures);
    }
    return comm.allreduce_scalar(static_cast<double>(failures)) == 0.0 ? 0 : 1;
}
