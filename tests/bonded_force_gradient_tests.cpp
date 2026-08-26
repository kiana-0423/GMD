// Finite-difference validation of bonded forces against their own energy.
//
// For every term, F_ia must equal -dU/dr_ia. This is independent of the virial:
// a force error shows up here directly, whereas in the virial it can hide.
// The proper-dihedral kernel is the reason this file exists. Its force had two
// defects that no existing check caught:
//
//   - the overall sign of the terminal-atom forces was flipped, and
//   - the middle-atom projection used the coefficients belonging to the
//     opposite convention b1 = r_i - r_j, while b1 is built as r_j - r_i.
//
// Neither is visible in tr(W): a torsion angle is invariant under isotropic
// scaling, so the dihedral virial trace is zero either way. It took a
// per-component check to expose it, and this file pins the forces themselves.

#include <algorithm>
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
#include "gmd/force/force_provider.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

namespace {

int failures = 0;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[bonded grad] " << message << '\n';
        ++failures;
    }
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

double energy_with_offset(gmd::ForceProvider& provider,
                          const gmd::System& system,
                          std::size_t atom,
                          std::size_t dim,
                          double delta) {
    gmd::System copy = system;
    copy.mutable_coordinates()[atom][dim] += delta;
    return evaluate(provider, copy).potential_energy;
}

// F_ia == -dU/dr_ia for every atom and every component.
//
// `tolerance` is RELATIVE to the largest force in the fixture, not absolute.
// These fixtures span two orders of magnitude in force (a few eV/A for a lone
// torsion, ~170 eV/A with a stiff bond active), so a single absolute bound is
// either meaningless for one case or unreachably tight for the other.
void check_forces_match_gradient(gmd::ForceProvider& provider,
                                 gmd::System& system,
                                 const std::string& label,
                                 double tolerance = 1.0e-8) {
    const gmd::ForceResult result = evaluate(provider, system);
    check(result.success, label + ": force evaluation must succeed");

    double worst = 0.0;
    double largest_force = 0.0;

    // Central-difference step. The error is round-off dominated below ~1e-6 and
    // truncation dominated above ~1e-4; measured on the combined fixture
    // (U ~ 296 eV, max|F| ~ 170 eV/A) the max deviation runs
    //
    //   h=1e-8: 7.0e-6   h=1e-7: 9.4e-7   h=1e-6: 1.1e-7
    //   h=1e-5: 1.3e-8   h=1e-4: 1.1e-6   h=1e-3: 1.1e-4
    //
    // so 1e-5 sits at the minimum. An earlier 1e-7 left only a few percent of
    // margin against the bound and flipped between compilers -- it passed on
    // clang/macOS and failed on gcc/Linux -- which is a property of the step
    // size, not of the forces (the relative error there was still 5.5e-9).
    constexpr double h = 1.0e-5;

    for (std::size_t atom = 0; atom < system.num_local_atoms(); ++atom) {
        for (std::size_t dim = 0; dim < 3; ++dim) {
            const double up = energy_with_offset(provider, system, atom, dim, h);
            const double down = energy_with_offset(provider, system, atom, dim, -h);
            const double numerical = -(up - down) / (2.0 * h);
            const double analytic = result.forces[atom][dim];

            worst = std::max(worst, std::abs(analytic - numerical));
            largest_force = std::max(largest_force, std::abs(numerical));

            // Sign agreement is called out separately: a flipped force is the
            // exact defect this file was written for.
            if (std::abs(numerical) > 1.0e-6) {
                check((analytic > 0.0) == (numerical > 0.0),
                      label + ": force on atom " + std::to_string(atom) +
                          " component " + std::to_string(dim) +
                          " has the wrong sign (analytic " + std::to_string(analytic) +
                          ", -dU/dx " + std::to_string(numerical) + ")");
            }
        }
    }

    const double limit = tolerance * std::max(1.0, largest_force);
    check(worst <= limit,
          label + ": forces do not match -dU/dx (max component difference " +
              std::to_string(worst) + ", limit " + std::to_string(limit) +
              " = " + std::to_string(tolerance) + " x max|F| " +
              std::to_string(largest_force) + ")");

    // Guard against a vacuous pass on a configuration that produces no force.
    check(largest_force > 1.0e-3,
          label + ": the test geometry must produce non-trivial forces, largest was " +
              std::to_string(largest_force));

    // Bonded terms are internal, so they must not push the molecule as a whole.
    double net[3] = {0.0, 0.0, 0.0};
    for (std::size_t atom = 0; atom < system.num_local_atoms(); ++atom) {
        for (std::size_t dim = 0; dim < 3; ++dim) {
            net[dim] += result.forces[atom][dim];
        }
    }
    check(std::abs(net[0]) < 1.0e-9 && std::abs(net[1]) < 1.0e-9 &&
              std::abs(net[2]) < 1.0e-9,
          label + ": bonded forces must sum to zero");
}

gmd::System make_chain(std::size_t atom_count) {
    gmd::System system;
    system.resize(atom_count, atom_count);
    gmd::Box box;
    // Non-cubic, so an axis mix-up cannot hide.
    box.set_lengths({14.0, 17.0, 20.0});
    system.set_box(box);

    const double positions[4][3] = {
        {6.0, 6.0, 6.0}, {7.5, 6.1, 6.4}, {8.2, 7.4, 7.3}, {9.6, 7.2, 8.1},
    };
    for (std::size_t i = 0; i < atom_count && i < 4; ++i) {
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {positions[i][0], positions[i][1], positions[i][2]};
    }
    return system;
}

// Same geometry, translated so the molecule straddles the x boundary.
gmd::System wrapped_chain(std::size_t atom_count) {
    gmd::System system = make_chain(atom_count);
    auto coordinates = system.mutable_coordinates();
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        coordinates[i][0] -= 7.0;
        if (coordinates[i][0] < 0.0) {
            coordinates[i][0] += system.box().lengths[0];
        }
    }
    return system;
}

void test_bonds() {
    auto topology = std::make_shared<gmd::Topology>();
    topology->bonds = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}};
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_bond_type({300.0, 1.5});

    gmd::System system = make_chain(4);
    check_forces_match_gradient(*provider, system, "bonds");

    gmd::System wrapped = wrapped_chain(4);
    check_forces_match_gradient(*provider, wrapped, "bonds across x boundary");
}

void test_angles() {
    auto topology = std::make_shared<gmd::Topology>();
    topology->angles = {{0, 1, 2, 0}, {1, 2, 3, 0}};
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_angle_type({60.0, 112.0 * gmd::deg_to_rad});

    gmd::System system = make_chain(4);
    check_forces_match_gradient(*provider, system, "angles");

    gmd::System wrapped = wrapped_chain(4);
    check_forces_match_gradient(*provider, wrapped, "angles across x boundary");
}

void test_dihedrals() {
    auto topology = std::make_shared<gmd::Topology>();
    topology->dihedrals = {{0, 1, 2, 3, 0}};
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_dihedral_type({1.5, 3, 0.0});

    gmd::System system = make_chain(4);
    check_forces_match_gradient(*provider, system, "proper dihedral");

    gmd::System wrapped = wrapped_chain(4);
    check_forces_match_gradient(*provider, wrapped, "proper dihedral across x boundary");
}

// A phase shift and a different multiplicity exercise a different point on the
// torsion profile, where dV/dphi has the opposite sign.
void test_dihedrals_shifted_phase() {
    auto topology = std::make_shared<gmd::Topology>();
    topology->dihedrals = {{0, 1, 2, 3, 0}};
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_dihedral_type({2.0, 2, 1.0});

    gmd::System system = make_chain(4);
    check_forces_match_gradient(*provider, system, "proper dihedral (n=2, delta=1)");
}

void test_impropers() {
    auto topology = std::make_shared<gmd::Topology>();
    topology->impropers = {{0, 1, 2, 3, 0}};
    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_improper_type({40.0, 10.0 * gmd::deg_to_rad});

    gmd::System system = make_chain(4);
    check_forces_match_gradient(*provider, system, "improper dihedral");
}

// All four term types active together, which is how a real force field runs.
void test_combined() {
    auto topology = std::make_shared<gmd::Topology>();
    topology->bonds = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}};
    topology->angles = {{0, 1, 2, 0}, {1, 2, 3, 0}};
    topology->dihedrals = {{0, 1, 2, 3, 0}};
    topology->impropers = {{0, 1, 2, 3, 0}};

    auto provider = std::make_shared<gmd::BondedForceProvider>(topology);
    provider->add_bond_type({300.0, 1.5});
    provider->add_angle_type({60.0, 112.0 * gmd::deg_to_rad});
    provider->add_dihedral_type({1.5, 3, 0.0});
    provider->add_improper_type({40.0, 10.0 * gmd::deg_to_rad});

    gmd::System system = make_chain(4);
    check_forces_match_gradient(*provider, system, "combined bonded terms");
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif

    test_bonds();
    test_angles();
    test_dihedrals();
    test_dihedrals_shifted_phase();
    test_impropers();
    test_combined();

    if (failures != 0) {
        std::cerr << "[bonded grad] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[bonded grad] all checks passed\n";
    return 0;
}
