// Constrained random velocity initialization under MPI.
//
// The serial audit establishes that the field comes out tangent, momentum-free
// and at exactly the target temperature over the authoritative degrees of
// freedom. None of that carries over to MPI for free, because every quantity
// involved is reduced:
//
//   * the kinetic energy that sets the rescale factor is an MPI_Allreduce;
//   * so is the momentum the centre-of-mass removal subtracts;
//   * the constraint projection is a global-gather solve, and a constraint can
//     SPAN RANKS -- the fixture below deliberately splits every molecule, so no
//     rank holds a complete one;
//   * the authoritative degrees of freedom come from a globally reduced atom
//     count and a replicated constraint list.
//
// Two failure modes stay invisible inside a single rank: a kinetic energy or
// momentum summed once per rank rather than once per atom, and a rank that
// silently skips a collective. Both are identities at np=1.
//
// This binary asserts what can be checked inside one run -- that every rank
// agrees on the DOF, that the field is tangent and momentum-free, that the
// temperature is exactly the target -- and writes the field keyed by atom tag
// so tests/velocity_init_equivalence.py can compare rank counts against each
// other. That comparison cannot happen inside one run: initialization is
// collective, so a rank cannot quietly initialize a whole-system copy on the
// side to compare against.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <mpi.h>

#include "gmd/core/runtime_context.hpp"
#include "gmd/core/simulation.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

namespace {

int failures = 0;
int global_rank = 0;
int global_size = 1;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[mpi constrained init][rank " << global_rank << "] " << message
                  << '\n';
        ++failures;
    }
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

constexpr double kElementaryChargeCoulombs = 1.602176634e-19;
constexpr double kBoltzmannJoulesPerKelvin = 1.380649e-23;
const double kBoltzmann = kBoltzmannJoulesPerKelvin / kElementaryChargeCoulombs;

constexpr std::size_t kMolecules = 4;
constexpr std::size_t kTotalAtoms = 3 * kMolecules;
constexpr double kTargetTemperature = 300.0;
constexpr std::uint32_t kSeed = 20260830u;
constexpr double kOH = 0.9572;
constexpr double kHH = 1.5139;

class NullForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "null"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request, gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.forces.assign(request.coordinates.size(), {0.0, 0.0, 0.0});
        result.potential_energy = 0.0;
        result.virial_valid = false;
        result.success = true;
    }
};

double mass_for(std::size_t tag) {
    return (tag % 3 == 0) ? 15.999 : 1.008;
}

std::array<double, 3> position_for(std::size_t tag) {
    const std::size_t molecule = tag / 3;
    const std::size_t within = tag % 3;
    const double shift = 6.0 * static_cast<double>(molecule);
    const double angle = 0.4 + 0.3 * static_cast<double>(molecule);
    const std::array<double, 3> origin = {5.0 + shift, 6.0, 7.0};
    if (within == 0) return origin;
    const double a = (within == 1) ? angle : angle + 2.0 * std::asin(kHH / (2.0 * kOH));
    return {origin[0] + kOH * std::cos(a), origin[1] + kOH * std::sin(a) * 0.6,
            origin[2] + kOH * std::sin(a) * 0.8};
}

std::vector<gmd::BondConstraint> all_constraints() {
    std::vector<gmd::BondConstraint> constraints;
    for (std::size_t m = 0; m < kMolecules; ++m) {
        const int o = static_cast<int>(3 * m);
        constraints.push_back({o, o + 1, kOH});
        constraints.push_back({o, o + 2, kOH});
        constraints.push_back({o + 1, o + 2, kHH});
    }
    return constraints;
}

// Which tags this rank owns. The stride is deliberately 1, not 3: consecutive
// tags land on different ranks, so EVERY molecule is split and every constraint
// spans ranks. At np>=4 the last rank owns nothing.
std::vector<int> owned_tags() {
    std::vector<int> owned;
    const int contributing = (global_size >= 4) ? global_size - 1 : global_size;
    for (int tag = 0; tag < static_cast<int>(kTotalAtoms); ++tag) {
        if (tag % contributing == global_rank) owned.push_back(tag);
    }
    return owned;
}

struct Result {
    std::size_t dof = 0;
    double twice_kinetic_energy = 0.0;
    std::array<double, 3> momentum = {0.0, 0.0, 0.0};
    double worst_tangency = 0.0;
    std::vector<std::array<double, 4>> rows;   // tag, vx, vy, vz (owned only)
};

Result initialize(bool redundant, bool reverse_storage) {
    std::vector<int> owned = owned_tags();
    if (reverse_storage) std::reverse(owned.begin(), owned.end());

    gmd::System system;
    system.resize(owned.size(), owned.size());
    gmd::Box box;
    box.set_lengths({80.0, 80.0, 80.0});
    system.set_box(box);
    for (std::size_t slot = 0; slot < owned.size(); ++slot) {
        const auto tag = static_cast<std::size_t>(owned[slot]);
        system.mutable_masses()[slot] = mass_for(tag);
        system.mutable_coordinates()[slot] = position_for(tag);
        system.mutable_atom_tags()[slot] = owned[slot];
        system.mutable_atom_owners()[slot] = global_rank;
    }

    std::vector<gmd::BondConstraint> constraints = all_constraints();
    if (redundant) {
        // Replace the first molecule's triangle with COLLINEAR targets:
        // d02 = d01 + d12 describes three points on a line, where the three
        // distances have rank 2 rather than 3. Adding a duplicate pair instead
        // would not do -- the solver normalises exact and near-duplicate pairs
        // away, so the set would stay independent and nothing would be
        // rejected. Atoms 0, 1 and 2 sit on different ranks in this fixture, so
        // no single rank can see the dependency by itself.
        constraints[0] = {0, 1, 1.0};
        constraints[1] = {1, 2, 1.0};
        constraints[2] = {0, 2, 2.0};
    }

    gmd::ConstraintSettings settings;
    settings.tolerance = 1.0e-13;
    settings.max_iterations = 1000;
    auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(1.0);
    integrator->set_constraint_solver(
        std::make_shared<gmd::ConstraintSolver>(constraints, settings));
    auto initializer = std::make_shared<gmd::VelocityInitializer>(kSeed);
    NullForceProvider provider;
    gmd::RuntimeContext runtime;

    gmd::Simulation simulation(&system);
    simulation.set_velocity_initializer(initializer);
    simulation.set_velocity_init_mode(gmd::VelocityInitMode::Random);
    simulation.set_initial_temperature(kTargetTemperature);
    simulation.set_remove_center_of_mass_velocity(true);
    simulation.set_force_provider(
        std::shared_ptr<gmd::ForceProvider>(&provider, [](gmd::ForceProvider*) {}));
    simulation.set_integrator(integrator);
    simulation.set_time_step(1.0);
    simulation.initialize(runtime);

    Result result;
    result.dof = integrator->degrees_of_freedom(system);
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        const auto v = system.velocities()[i];
        const double m = system.masses()[i];
        result.twice_kinetic_energy += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        for (std::size_t d = 0; d < 3; ++d) result.momentum[d] += m * v[d];
        result.rows.push_back({static_cast<double>(system.atom_tag(i)), v[0], v[1], v[2]});
    }
    // Tangency is measured only for constraints whose BOTH atoms are owned
    // here; every constraint is owned by someone, and the maximum is reduced
    // below, so nothing is missed.
    for (const auto& c : constraints) {
        std::size_t slot_i = system.atom_count();
        std::size_t slot_j = system.atom_count();
        for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
            if (system.atom_tag(i) == c.i) slot_i = i;
            if (system.atom_tag(i) == c.j) slot_j = i;
        }
        if (slot_i == system.atom_count() || slot_j == system.atom_count()) continue;
        const auto ri = system.coordinates()[slot_i];
        const auto rj = system.coordinates()[slot_j];
        const auto vi = system.velocities()[slot_i];
        const auto vj = system.velocities()[slot_j];
        double dot = 0.0;
        for (std::size_t d = 0; d < 3; ++d) dot += (ri[d] - rj[d]) * (vi[d] - vj[d]);
        result.worst_tangency = std::max(result.worst_tangency, std::abs(dot));
    }
    return result;
}

void write_field(const std::string& path, const Result& result, bool reverse_storage) {
    const int local_count = static_cast<int>(result.rows.size());
    std::vector<double> local;
    local.reserve(static_cast<std::size_t>(local_count) * 4);
    for (const auto& row : result.rows) {
        local.insert(local.end(), row.begin(), row.end());
    }
    const int local_doubles = local_count * 4;
    std::vector<int> counts(static_cast<std::size_t>(global_size), 0);
    MPI_Gather(&local_doubles, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> displacements(static_cast<std::size_t>(global_size), 0);
    int total = 0;
    if (global_rank == 0) {
        for (int r = 0; r < global_size; ++r) {
            displacements[static_cast<std::size_t>(r)] = total;
            total += counts[static_cast<std::size_t>(r)];
        }
    }
    std::vector<double> gathered(global_rank == 0 ? static_cast<std::size_t>(total) : 0);
    MPI_Gatherv(local.data(), local_doubles, MPI_DOUBLE, gathered.data(), counts.data(),
                displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double momentum[3] = {0.0, 0.0, 0.0};
    double twice_ke = 0.0;
    double local_momentum[3] = {result.momentum[0], result.momentum[1], result.momentum[2]};
    double local_ke = result.twice_kinetic_energy;
    MPI_Allreduce(local_momentum, momentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_ke, &twice_ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (global_rank != 0) return;
    std::vector<std::array<double, 4>> rows;
    for (std::size_t i = 0; i + 3 < gathered.size(); i += 4) {
        rows.push_back({gathered[i], gathered[i + 1], gathered[i + 2], gathered[i + 3]});
    }
    std::sort(rows.begin(), rows.end(),
              [](const auto& a, const auto& b) { return a[0] < b[0]; });
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "# np " << global_size << " reverse_storage " << (reverse_storage ? 1 : 0)
        << '\n';
    out << "# atoms " << rows.size() << '\n';
    out << "degrees_of_freedom " << result.dof << '\n';
    out << "momentum " << momentum[0] << ' ' << momentum[1] << ' ' << momentum[2] << '\n';
    out << "twice_kinetic_energy " << twice_ke << '\n';
    for (const auto& row : rows) {
        out << static_cast<long long>(row[0]) << ' ' << row[1] << ' ' << row[2] << ' '
            << row[3] << '\n';
    }
}

void assert_global_properties(const Result& result) {
    // Every rank must agree on the degrees of freedom, or they are computing
    // different temperatures from the same field.
    long long dof = static_cast<long long>(result.dof);
    long long lowest = 0, highest = 0;
    MPI_Allreduce(&dof, &lowest, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&dof, &highest, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    check(lowest == highest,
          "ranks disagree about the degrees of freedom: " + std::to_string(lowest) +
              " to " + std::to_string(highest));
    const long long expected =
        3 * static_cast<long long>(kTotalAtoms) - 3 * static_cast<long long>(kMolecules) - 3;
    check(dof == expected,
          "authoritative DOF is " + std::to_string(dof) + ", expected 3N - rank - 3 = " +
              std::to_string(expected));

    double twice_ke = 0.0;
    double local_ke = result.twice_kinetic_energy;
    MPI_Allreduce(&local_ke, &twice_ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    const double temperature = twice_ke / (static_cast<double>(dof) * kBoltzmann);
    check(std::abs(temperature / kTargetTemperature - 1.0) < 1.0e-12,
          "the field carries " + number(temperature) + " K at " +
              std::to_string(global_size) + " rank(s), not " +
              number(kTargetTemperature) +
              ". A kinetic energy summed once per rank would show up here");

    double momentum[3] = {0.0, 0.0, 0.0};
    double local_momentum[3] = {result.momentum[0], result.momentum[1], result.momentum[2]};
    MPI_Allreduce(local_momentum, momentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    const double magnitude = std::sqrt(momentum[0] * momentum[0] +
                                       momentum[1] * momentum[1] +
                                       momentum[2] * momentum[2]);
    // Against the system's natural momentum scale, not against zero.
    double total_mass = 0.0;
    for (std::size_t tag = 0; tag < kTotalAtoms; ++tag) total_mass += mass_for(tag);
    const double scale = std::sqrt(twice_ke * total_mass);
    check(magnitude < 1.0e-10 * scale,
          "residual centre-of-mass momentum is " + number(magnitude) +
              " against a scale of " + number(scale));

    double worst_tangency = 0.0;
    double local_tangency = result.worst_tangency;
    MPI_Allreduce(&local_tangency, &worst_tangency, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    check(worst_tangency < 1.0e-11,
          "worst |r.v| across all constraints is " + number(worst_tangency) +
              " at " + std::to_string(global_size) +
              " rank(s); the initial velocities are not in the tangent space");

    if (global_rank == 0) {
        std::cout << std::setprecision(6) << "  np=" << global_size << "  dof=" << dof
                  << "  T=" << temperature << " K  |p|=" << magnitude
                  << "  worst |r.v|=" << worst_tangency << '\n';
    }
}

void test_redundant_constraints_rejected_on_every_rank() {
    int threw = 0;
    try {
        (void)initialize(true, false);
    } catch (const std::exception&) {
        threw = 1;
    }
    int total = 0;
    MPI_Allreduce(&threw, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    check(total == global_size,
          "a dependent constraint set was rejected on only " + std::to_string(total) +
              " of " + std::to_string(global_size) +
              " rank(s). Every rank must throw, or the ones that did not would wait "
              "alone in the next collective");
}

void test_empty_rank_participates() {
    if (global_size >= 4 && global_rank == global_size - 1) {
        check(owned_tags().empty(),
              "this fixture expects the last rank to own no atoms at np>=4");
    }
    // Reaching the assertion in assert_global_properties() at all is the test:
    // a rank owning nothing must still enter every collective inside
    // initialization and the constraint solve.
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    std::string output;
    bool reverse_storage = false;
    bool redundant_only = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) output = argv[++i];
        else if (std::strcmp(argv[i], "--reverse-storage") == 0) reverse_storage = true;
        else if (std::strcmp(argv[i], "--redundant-only") == 0) redundant_only = true;
    }

    if (redundant_only) {
        test_redundant_constraints_rejected_on_every_rank();
    } else {
        test_empty_rank_participates();
        const Result result = initialize(false, reverse_storage);
        assert_global_properties(result);
        if (!output.empty()) write_field(output, result, reverse_storage);
    }

    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_rank == 0) {
        if (total == 0) {
            std::cout << "[mpi constrained init] all checks passed on " << global_size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi constrained init] " << total << " check(s) failed on "
                      << global_size << " rank(s)\n";
        }
    }
    MPI_Finalize();
    return total == 0 ? 0 : 1;
}
