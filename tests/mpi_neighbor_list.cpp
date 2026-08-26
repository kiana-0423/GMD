// MPI regression tests for neighbor-list correctness across ranks.
//
// "rebuild" — VerletNeighborBuilder::needs_rebuild() only inspects locally
//   owned atoms, so when a single rank's atom moves past skin/2 that rank alone
//   used to decide to rebuild while the others carried on with a list that no
//   longer covered their ghost copies of it. Simulation now reduces the local
//   flags with a logical OR, making the rebuild all-or-nothing.
//
//   Note on coverage: with a DomainDecomposition attached, the list is already
//   rebuilt every step for an unrelated reason -- reverse_accumulate_ghost_forces()
//   ends in clear_ghost_atoms(), which clears the neighbor list, so the
//   `!valid` term short-circuits needs_rebuild() before it is ever consulted.
//   To exercise the collective decision itself, this test runs Simulation with a
//   communicator but no domain decomposition, which is the configuration in
//   which a valid list survives from one step to the next.
//
// "image" — image_flags[k] must satisfy, for ghosts exactly as for local atoms,
//
//       r_j + S * L - r_i  ==  minimum_image(r_j - r_i)
//
//   where the ghost's stored coordinate is already the image adjacent to this
//   rank's subdomain. That is the convention the TorchScript adapter feeds to
//   the model as edge_shift.

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/core/simulation.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/minimum_image.hpp"
#include "gmd/system/neighbor_builder.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/verlet_neighbor_builder.hpp"

#include <mpi.h>

namespace {

constexpr double kBoxLength = 12.0;
constexpr double kCutoff = 2.5;
constexpr double kSkin = 1.0;

using Vec3 = gmd::System::Vec3;

void check(bool condition, const std::string& message, int rank, int& failures) {
    if (!condition) {
        std::cerr << "[mpi neighbor][rank " << rank << "] " << message << '\n';
        ++failures;
    }
}

gmd::Box make_box() {
    gmd::Box box;
    box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    return box;
}

// A regular lattice of atoms spread across the box, so that a 1D split gives
// every rank a populated subdomain with ghosts on both faces.
std::vector<Vec3> lattice_positions() {
    std::vector<Vec3> positions;
    for (int ix = 0; ix < 6; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                positions.push_back({1.0 + 2.0 * ix, 3.0 + 5.0 * iy, 3.0 + 5.0 * iz});
            }
        }
    }
    return positions;
}

// Builds the system holding only the atoms this rank owns, tagged globally.
gmd::System make_local_system(const gmd::Box& box,
                              const std::vector<Vec3>& all_positions,
                              const gmd::DomainDecomposition& dd,
                              int rank) {
    std::vector<int> owned;
    for (std::size_t i = 0; i < all_positions.size(); ++i) {
        if (dd.owner_rank(box, all_positions[i]) == rank) {
            owned.push_back(static_cast<int>(i));
        }
    }

    gmd::System system;
    system.resize(owned.size(), owned.size());
    system.set_box(box);
    for (std::size_t k = 0; k < owned.size(); ++k) {
        const auto tag = static_cast<std::size_t>(owned[k]);
        system.mutable_masses()[k] = 1.0;
        system.mutable_charges()[k] = 0.0;
        system.mutable_atom_types()[k] = 0;
        system.mutable_atom_tags()[k] = owned[k];
        system.mutable_atom_owners()[k] = rank;
        system.mutable_coordinates()[k] = all_positions[tag];
    }
    return system;
}

// ---------------------------------------------------------------------------
// Test 1: the rebuild decision must be collective.
// ---------------------------------------------------------------------------

// Produces no forces, so begin_step() leaves the coordinates untouched and the
// only displacement in play is the one this test injects by hand.
class ZeroForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "zero"; }
    void initialize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.success = true;
        result.potential_energy = 0.0;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        result.virial_valid = true;
    }
    void finalize(gmd::RuntimeContext&) override {}
};

// Counts rebuild() calls, forwarding everything to a real VerletNeighborBuilder
// (which is final, hence composition rather than inheritance).
//
// Counting is the right observable here: comparing reference coordinates would
// not work, because rebuild() refreshes them from the current positions and on a
// rank whose atoms did not move that leaves them bit-identical to the previous
// snapshot -- a rank that rebuilt would look exactly like one that did not.
class CountingNeighborBuilder final : public gmd::NeighborBuilder {
public:
    CountingNeighborBuilder(double r_cut, double r_skin) noexcept
        : inner_(r_cut, r_skin) {}

    std::string_view name() const noexcept override { return "counting_neighbor_builder"; }

    void initialize(gmd::System& system, gmd::RuntimeContext& runtime) override {
        ++rebuild_count;
        inner_.initialize(system, runtime);
    }

    bool needs_rebuild(const gmd::System& system, std::uint64_t step) const override {
        return inner_.needs_rebuild(system, step);
    }

    void rebuild(gmd::System& system,
                 gmd::RuntimeContext& runtime,
                 gmd::NeighborBuildStats* stats) override {
        ++rebuild_count;
        inner_.rebuild(system, runtime, stats);
    }

    int rebuild_count = 0;

private:
    gmd::VerletNeighborBuilder inner_;
};

int test_rebuild_is_collective(int rank, int /*size*/) {
    int failures = 0;

    const gmd::Box box = make_box();
    const auto positions = lattice_positions();

    // Every rank holds the same replicated system, all atoms local. Without a
    // DomainDecomposition there is no ghost exchange and therefore no
    // clear_ghost_atoms(), so a valid neighbor list survives across the step and
    // the rebuild decision is genuinely taken by needs_rebuild().
    gmd::System system;
    system.resize(positions.size(), positions.size());
    system.set_box(box);
    for (std::size_t i = 0; i < positions.size(); ++i) {
        system.mutable_masses()[i] = 1.0;
        system.mutable_atom_tags()[i] = static_cast<int>(i);
        system.mutable_atom_owners()[i] = rank;
        system.mutable_coordinates()[i] = positions[i];
    }

    auto builder = std::make_shared<CountingNeighborBuilder>(kCutoff, kSkin);
    auto comm = std::make_shared<gmd::MpiCommunicator>();
    auto provider = std::make_shared<ZeroForceProvider>();
    auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(1.0);

    gmd::Simulation simulation(&system);
    simulation.set_force_provider(provider);
    simulation.set_neighbor_builder(builder);
    simulation.set_integrator(integrator);
    simulation.set_mpi_communicator(comm);
    simulation.set_time_step(1.0);

    gmd::RuntimeContext runtime;
    simulation.initialize(runtime);

    check(system.neighbor_list().valid,
          "neighbor list must be valid after initialize", rank, failures);

    check(!system.neighbor_list().ref_coordinates.empty(),
          "initialize must have installed reference coordinates", rank, failures);

    const int rebuilds_after_initialize = builder->rebuild_count;
    check(rebuilds_after_initialize > 0,
          "initialize must have built the neighbor list at least once", rank, failures);

    // Displace one atom on rank 0 only, by clearly more than the skin/2 trigger.
    if (rank == 0) {
        system.mutable_coordinates()[0][1] += kSkin;
    }

    // Premise of the test: only rank 0 detects the displacement locally. If this
    // stops holding, the test no longer exercises the divergence it was written
    // for, so assert it rather than assuming it.
    const bool local_flag = builder->needs_rebuild(system, 1);
    check(local_flag == (rank == 0),
          std::string("expected local needs_rebuild() to be ") +
              (rank == 0 ? "true on rank 0" : "false off rank 0") +
              " but it was " + (local_flag ? "true" : "false"),
          rank, failures);

    int local_flag_count = local_flag ? 1 : 0;
    int total_flag_count = 0;
    MPI_Allreduce(&local_flag_count, &total_flag_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    check(total_flag_count == 1,
          "exactly one rank should have detected the displacement locally, got " +
              std::to_string(total_flag_count),
          rank, failures);

    // The primitive the fix is built on.
    check(comm->allreduce_logical_or(local_flag),
          "allreduce_logical_or must be true on every rank when any rank is true",
          rank, failures);
    check(!comm->allreduce_logical_or(false),
          "allreduce_logical_or must be false when no rank is true", rank, failures);

    // Drive one step through Simulation, which is where the collective decision
    // is taken. The zero-force provider means nothing else moves the atoms.
    simulation.step(runtime);

    // Every rank must have rebuilt during that step, including the three that
    // saw no local motion at all. Before the collective reduction, only rank 0
    // would have rebuilt here and the others would have carried a stale list.
    check(builder->rebuild_count > rebuilds_after_initialize,
          "this rank did not rebuild its neighbor list even though another rank "
          "moved an atom past skin/2 (stale-list divergence)",
          rank, failures);

    // And the converse: with nothing moving, no rank rebuilds. This is what
    // keeps the fix from degenerating into "rebuild every step".
    const int rebuilds_before_quiet_step = builder->rebuild_count;
    simulation.step(runtime);
    check(builder->rebuild_count == rebuilds_before_quiet_step,
          "no rank should rebuild on a step where no atom moved past skin/2",
          rank, failures);

    return failures;
}

// ---------------------------------------------------------------------------
// Test 2: ghost pairs must honour the image-flag convention.
// ---------------------------------------------------------------------------
int test_ghost_image_flags(int rank, int size) {
    int failures = 0;

    const gmd::Box box = make_box();
    const auto positions = lattice_positions();

    auto dd = std::make_shared<gmd::DomainDecomposition>();
    dd->create_decomposition(box, size, rank, kCutoff, kSkin, {true, true, true});

    gmd::System system = make_local_system(box, positions, *dd, rank);

    gmd::MpiCommunicator comm;
    dd->refresh(system.box());
    comm.exchange_ghost_coordinates(system, *dd);

    check(system.num_ghost_atoms() > 0,
          "expected this rank to receive ghost atoms; the test needs them to be "
          "meaningful",
          rank, failures);

    gmd::VerletNeighborBuilder builder(kCutoff, kSkin);
    builder.set_domain_decomposition(dd);
    gmd::RuntimeContext runtime;
    builder.rebuild(system, runtime, nullptr);

    const auto& nl = system.neighbor_list();
    const auto coords = system.coordinates();

    std::size_t ghost_pairs_checked = 0;
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const int start = nl.offsets[i];
        for (int k = 0; k < nl.counts[i]; ++k) {
            const auto idx = static_cast<std::size_t>(start + k);
            const auto j = static_cast<std::size_t>(nl.neighbors[idx]);
            const std::array<int, 3>& S = nl.image_flags[idx];

            std::array<double, 3> expected = {
                coords[j][0] - coords[i][0],
                coords[j][1] - coords[i][1],
                coords[j][2] - coords[i][2],
            };
            gmd::apply_minimum_image(expected, box);

            for (std::size_t d = 0; d < 3; ++d) {
                const double actual = coords[j][d]
                                    + static_cast<double>(S[d]) * box.lengths[d]
                                    - coords[i][d];
                check(std::abs(actual - expected[d]) < 1.0e-9,
                      std::string(system.is_local_atom(j) ? "local" : "ghost") +
                          " pair image flag mismatch on axis " + std::to_string(d) +
                          ": r_j + S*L - r_i = " + std::to_string(actual) +
                          ", minimum image = " + std::to_string(expected[d]),
                      rank, failures);
            }
            if (!system.is_local_atom(j)) {
                ++ghost_pairs_checked;
            }
        }
    }

    check(ghost_pairs_checked > 0,
          "expected at least one local/ghost pair in the neighbor list",
          rank, failures);

    return failures;
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const std::string mode = argc > 1 ? argv[1] : "rebuild";

    int failures = 0;
    if (mode == "rebuild") {
        failures = test_rebuild_is_collective(rank, size);
    } else if (mode == "image") {
        failures = test_ghost_image_flags(rank, size);
    } else {
        if (rank == 0) {
            std::cerr << "[mpi neighbor] unknown mode: " << mode << '\n';
        }
        failures = 1;
    }

    int total_failures = 0;
    MPI_Allreduce(&failures, &total_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        if (total_failures == 0) {
            std::cout << "[mpi neighbor] " << mode << ": all checks passed\n";
        } else {
            std::cerr << "[mpi neighbor] " << mode << ": " << total_failures
                      << " check(s) failed\n";
        }
    }

    return total_failures == 0 ? 0 : 1;
}
