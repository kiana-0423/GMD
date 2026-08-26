// MPI regression tests for degrees-of-freedom accounting and virial reduction.
//
// "dof"    — the DOF count must be a global quantity: identical on every rank
//   and equal to the value a single-process run would use. It is built from
//   global_atom_count() (an allreduce over locally owned atoms) and from the
//   constraint list, which ConstraintSolver replicates against global atom tags
//   and which therefore must NOT be reduced a second time.
//
// "virial" — each force provider reduces its own virial across the
//   communicator, so what compute() returns is already the global tensor. This
//   test checks that it is reduced exactly once, by holding it against two
//   finite-difference identities: tr(W) = -dU/ds under isotropic scaling, and
//   W_aa = -dU/de under a normal strain along axis a alone. A virial reduced
//   twice, or not at all, fails by a factor of the rank count; one that
//   double-counts ghost contributions fails by a smaller margin.
//
//   Only the diagonal components are checked. `Box` is orthorhombic, so no
//   shear strain can be applied and the off-diagonal components are NOT
//   validated here or anywhere else.

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

#include <mpi.h>

namespace {

constexpr double kBoxLength = 16.0;

using Vec3 = gmd::System::Vec3;

void check(bool condition, const std::string& message, int rank, int& failures) {
    if (!condition) {
        std::cerr << "[mpi dof/virial][rank " << rank << "] " << message << '\n';
        ++failures;
    }
}

std::vector<Vec3> cluster_positions() {
    return {
        {4.0, 4.0, 4.0}, {7.1, 4.2, 4.0}, {4.0, 7.3, 4.1}, {7.2, 7.0, 4.0},
        {4.1, 4.0, 7.2}, {7.0, 4.1, 7.0}, {4.0, 7.1, 7.3}, {7.3, 7.2, 7.1},
        {11.0, 4.0, 4.2}, {12.4, 7.1, 4.0}, {11.2, 4.1, 7.0}, {12.1, 7.3, 7.2},
    };
}

// Builds this rank's share of the cluster, with each axis scaled by its own
// factor. An isotropic scaling passes the same factor three times; a per-axis
// strain passes (1+e) on one axis and 1.0 on the others.
gmd::System make_decomposed_system(const gmd::DomainDecomposition& dd,
                                   int rank,
                                   const std::array<double, 3>& s) {
    gmd::Box box;
    box.set_lengths({kBoxLength * s[0], kBoxLength * s[1], kBoxLength * s[2]});

    // Ownership is decided on the unscaled geometry so that the same atoms stay
    // on the same rank as `s` varies; otherwise the finite difference would
    // compare two different decompositions.
    gmd::Box reference_box;
    reference_box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    const auto positions = cluster_positions();

    std::vector<int> owned;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (dd.owner_rank(reference_box, positions[i]) == rank) {
            owned.push_back(static_cast<int>(i));
        }
    }

    gmd::System system;
    system.resize(owned.size(), owned.size());
    system.set_box(box);
    for (std::size_t k = 0; k < owned.size(); ++k) {
        const auto tag = static_cast<std::size_t>(owned[k]);
        system.mutable_masses()[k] = 12.0;
        system.mutable_atom_types()[k] = 0;
        system.mutable_atom_tags()[k] = owned[k];
        system.mutable_atom_owners()[k] = rank;
        system.mutable_charges()[k] = (tag % 2 == 0) ? 0.6 : -0.6;
        system.mutable_coordinates()[k] = {positions[tag][0] * s[0],
                                           positions[tag][1] * s[1],
                                           positions[tag][2] * s[2]};
    }
    return system;
}

// ---------------------------------------------------------------------------
// Test 1: DOF is global and identical on every rank.
// ---------------------------------------------------------------------------
int test_dof_consistency(int rank, int size) {
    int failures = 0;

    gmd::Box reference_box;
    reference_box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    gmd::DomainDecomposition dd;
    dd.create_decomposition(reference_box, size, rank, 4.0, 1.0, {true, true, true});

    gmd::System system = make_decomposed_system(dd, rank, {1.0, 1.0, 1.0});
    const std::size_t global_atoms = cluster_positions().size();

    // The atoms really are split up; otherwise this test proves nothing.
    int local_atoms = static_cast<int>(system.num_local_atoms());
    int summed_atoms = 0;
    MPI_Allreduce(&local_atoms, &summed_atoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    check(summed_atoms == static_cast<int>(global_atoms),
          "local atom counts must sum to the global count, got " +
              std::to_string(summed_atoms),
          rank, failures);
    check(gmd::global_atom_count(system) == global_atoms,
          "global_atom_count() must return the global count on every rank",
          rank, failures);

    // A chain of constraints over the global tags, replicated on every rank
    // exactly as ConstraintSolver expects.
    gmd::Topology topology;
    for (std::size_t i = 0; i + 1 < global_atoms; ++i) {
        topology.bonds.push_back({static_cast<int>(i), static_cast<int>(i + 1), 0});
    }
    auto constraints = std::make_shared<gmd::ConstraintSolver>(
        gmd::constraints_from_bond_types(topology, {0}, {1.0}),
        gmd::ConstraintSettings{});

    struct Case {
        bool remove_com;
        bool constrained;
        std::size_t expected;
    };
    // 12 atoms -> 3N = 36; 11 chain constraints.
    const Case cases[] = {
        {true,  false, 33},
        {false, false, 36},
        {true,  true,  22},
        {false, true,  25},
    };

    for (const Case& c : cases) {
        gmd::VelocityVerletIntegrator integrator(1.0);
        integrator.set_remove_center_of_mass_velocity(c.remove_com);
        if (c.constrained) {
            integrator.set_constraint_solver(constraints);
        }

        const std::size_t dof = integrator.degrees_of_freedom(system);
        const std::string label = std::string("COM ") +
            (c.remove_com ? "removed" : "kept") +
            (c.constrained ? ", constrained" : ", unconstrained");

        check(dof == c.expected,
              label + ": expected DOF " + std::to_string(c.expected) + ", got " +
                  std::to_string(dof),
              rank, failures);

        // And every rank must agree, bit for bit.
        long long local_dof = static_cast<long long>(dof);
        long long min_dof = 0;
        long long max_dof = 0;
        MPI_Allreduce(&local_dof, &min_dof, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_dof, &max_dof, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        check(min_dof == max_dof,
              label + ": ranks disagree on DOF (min " + std::to_string(min_dof) +
                  ", max " + std::to_string(max_dof) + ")",
              rank, failures);
    }

    return failures;
}

// ---------------------------------------------------------------------------
// Test 2: the virial is reduced exactly once.
// ---------------------------------------------------------------------------

// Global potential energy for a provider at scale `s`. compute() returns the
// rank-local energy (Simulation is what normally reduces it), so the reduction
// is done explicitly here. The virial, by contrast, is already global.
double global_energy_and_virial(gmd::ForceProvider& provider,
                                const gmd::DomainDecomposition& dd_template,
                                int rank,
                                int size,
                                const std::array<double, 3>& s,
                                std::array<double, 9>* virial_out) {
    gmd::Box reference_box;
    reference_box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    gmd::DomainDecomposition dd;
    dd.create_decomposition(reference_box, size, rank, 4.0, 1.0, {true, true, true});
    (void)dd_template;

    gmd::System system = make_decomposed_system(dd, rank, s);
    dd.refresh(system.box());

    gmd::MpiCommunicator comm;
    comm.exchange_ghost_coordinates(system, dd);

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

    if (virial_out != nullptr) {
        *virial_out = result.virial;
    }

    double global = 0.0;
    MPI_Allreduce(&result.potential_energy, &global, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    return global;
}

void check_provider_virial(gmd::ForceProvider& provider,
                           const std::string& label,
                           int rank,
                           int size,
                           int& failures) {
    gmd::DomainDecomposition dd_template;
    constexpr double eps = 1.0e-6;
    static const char* axis_name[3] = {"xx", "yy", "zz"};

    std::array<double, 9> virial{};
    global_energy_and_virial(provider, dd_template, rank, size, {1.0, 1.0, 1.0}, &virial);

    // --- Isotropic: tr(W) = -dU/ds ---
    {
        const double analytic = virial[0] + virial[4] + virial[8];
        const double up = global_energy_and_virial(provider, dd_template, rank, size,
                                                   {1.0 + eps, 1.0 + eps, 1.0 + eps},
                                                   nullptr);
        const double down = global_energy_and_virial(provider, dd_template, rank, size,
                                                     {1.0 - eps, 1.0 - eps, 1.0 - eps},
                                                     nullptr);
        const double numerical = -(up - down) / (2.0 * eps);
        const double scale = std::max({std::abs(analytic), std::abs(numerical), 1.0e-12});
        const double relative_error = std::abs(analytic - numerical) / scale;

        check(relative_error < 1.0e-5,
              label + ": virial trace does not match -dU/ds under domain "
                      "decomposition (analytic " + std::to_string(analytic) +
                  ", finite difference " + std::to_string(numerical) +
                  ", relative error " + std::to_string(relative_error) +
                  "). A factor near the rank count means the virial was reduced the "
                  "wrong number of times.",
              rank, failures);
    }

    // --- Per-axis: W_aa = -dU/de, each diagonal component on its own ---
    for (std::size_t axis = 0; axis < 3; ++axis) {
        std::array<double, 3> up_scale = {1.0, 1.0, 1.0};
        std::array<double, 3> down_scale = {1.0, 1.0, 1.0};
        up_scale[axis] = 1.0 + eps;
        down_scale[axis] = 1.0 - eps;

        const double up = global_energy_and_virial(provider, dd_template, rank, size,
                                                   up_scale, nullptr);
        const double down = global_energy_and_virial(provider, dd_template, rank, size,
                                                     down_scale, nullptr);
        const double numerical = -(up - down) / (2.0 * eps);
        const double analytic = virial[axis * 3 + axis];

        const double scale = std::max({std::abs(analytic), std::abs(numerical), 1.0e-12});
        const double relative_error = std::abs(analytic - numerical) / scale;

        check(relative_error < 1.0e-5,
              label + ": W_" + axis_name[axis] +
                  " does not match -dU/de under domain decomposition (analytic " +
                  std::to_string(analytic) + ", finite difference " +
                  std::to_string(numerical) + ", relative error " +
                  std::to_string(relative_error) + ")",
              rank, failures);
    }

    // Every rank must hold the same reduced tensor, component by component.
    for (std::size_t component = 0; component < virial.size(); ++component) {
        double local = virial[component];
        double lo = 0.0;
        double hi = 0.0;
        MPI_Allreduce(&local, &lo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local, &hi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        check(std::abs(hi - lo) < 1.0e-9,
              label + ": ranks disagree on reduced virial component " +
                  std::to_string(component) + " (min " + std::to_string(lo) +
                  ", max " + std::to_string(hi) + ")",
              rank, failures);
    }
}

int test_virial_reduction(int rank, int size) {
    int failures = 0;

    {
        gmd::ClassicalForceProvider provider(0.25, 2.5, 4.0);
        check_provider_virial(provider, "lennard-jones", rank, size, failures);
    }
    {
        gmd::EwaldForceProvider provider(0.35, 6, 7.0);
        check_provider_virial(provider, "ewald", rank, size, failures);
    }
    {
        gmd::PMEForceProvider provider(0.35, 7.0, 6, {32, 32, 32});
        check_provider_virial(provider, "pme", rank, size, failures);
    }

    return failures;
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const std::string mode = argc > 1 ? argv[1] : "dof";

    int failures = 0;
    if (mode == "dof") {
        failures = test_dof_consistency(rank, size);
    } else if (mode == "virial") {
        failures = test_virial_reduction(rank, size);
    } else {
        if (rank == 0) {
            std::cerr << "[mpi dof/virial] unknown mode: " << mode << '\n';
        }
        failures = 1;
    }

    int total_failures = 0;
    MPI_Allreduce(&failures, &total_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        if (total_failures == 0) {
            std::cout << "[mpi dof/virial] " << mode << ": all checks passed\n";
        } else {
            std::cerr << "[mpi dof/virial] " << mode << ": " << total_failures
                      << " check(s) failed\n";
        }
    }

    return total_failures == 0 ? 0 : 1;
}
