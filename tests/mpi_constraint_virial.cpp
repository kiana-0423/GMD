// MPI regression test for the constraint virial: the answer must not depend on
// the rank count.
//
// The fixture is the rigid rotating dimer from tests/constraint_virial_tests.cpp
// -- a pair on the constraint manifold in both position and velocity, held
// together only by its constraint -- positioned so the two atoms straddle a
// domain boundary for every rank count exercised. Its ENDPOINT constraint force
// at t+dt is the centripetal force, which is a closed-form number:
//
//     |G| -> mu omega^2 d,     tr W -> -mu omega^2 d^2 = -2K
//
// Running the same fixture at 1, 2 and 4 ranks and comparing against that
// constant catches both failure modes directly: a contribution dropped on
// non-owning ranks shows up as zero or a fraction, and one that was reduced when
// it should not have been shows up multiplied by the rank count.
//
// OWNERSHIP RULE. There is no per-rank split to make. ConstraintSolver
// allgathers every owned atom (load_atoms) and works from a constraint list
// replicated against global atom tags, so every rank iterates the same
// constraints over the same coordinates and velocities and converges to the same
// multipliers. The tensor each rank computes is therefore ALREADY the global one
// and must NOT be reduced again. Ranks owning no atoms still enter the allgather
// inside load_atoms, so no rank takes a different collective path.

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/integrator.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

#include <mpi.h>

namespace {

constexpr double kBox = 20.0;
constexpr double kBond = 1.5;
constexpr double kOmega = 0.1;
constexpr double kTimeStep = 0.05;
constexpr double kMass0 = 1.0;
constexpr double kMass1 = 3.0;

using Vec3 = gmd::System::Vec3;

void check(bool condition, const std::string& message, int rank, int& failures) {
    if (!condition) {
        std::cerr << "[mpi constraint virial][rank " << rank << "] " << message << '\n';
        ++failures;
    }
}

void check_close(double actual, double expected, double tolerance,
                 const std::string& message, int rank, int& failures) {
    if (!(std::abs(actual - expected) <= tolerance)) {
        std::cerr << "[mpi constraint virial][rank " << rank << "] " << message
                  << ": expected " << expected << ", got " << actual
                  << " (|diff| " << std::abs(actual - expected) << ")\n";
        ++failures;
    }
}

gmd::ConstraintSettings tight_settings() {
    gmd::ConstraintSettings settings;
    settings.tolerance = 1.0e-13;
    settings.max_iterations = 1000;
    return settings;
}

struct Atom {
    Vec3 position;
    Vec3 velocity;
    double mass;
};

// Rigid rotor in the xy-plane, centre of mass at rest on x = 10, which is a
// domain boundary for both 2 and 4 ranks on a 1D split of a 20 A cell.
std::vector<Atom> rotating_dimer() {
    const double total = kMass0 + kMass1;
    const double x0 = (kMass1 / total) * kBond;
    const double x1 = -(kMass0 / total) * kBond;
    return {
        {{10.0 + x0, 10.0, 10.0}, {0.0, kOmega * x0, 0.0}, kMass0},
        {{10.0 + x1, 10.0, 10.0}, {0.0, kOmega * x1, 0.0}, kMass1},
    };
}

// This rank's share of the fixture, tagged globally.
gmd::System local_share(const gmd::DomainDecomposition& dd, const gmd::Box& box, int rank) {
    const auto atoms = rotating_dimer();
    std::vector<int> owned;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (dd.owner_rank(box, atoms[i].position) == rank) {
            owned.push_back(static_cast<int>(i));
        }
    }

    gmd::System system;
    system.resize(owned.size(), owned.size());
    system.set_box(box);
    for (std::size_t k = 0; k < owned.size(); ++k) {
        const auto index = static_cast<std::size_t>(owned[k]);
        system.mutable_masses()[k] = atoms[index].mass;
        system.mutable_atom_tags()[k] = owned[k];
        system.mutable_atom_owners()[k] = rank;
        system.mutable_coordinates()[k] = atoms[index].position;
        system.mutable_velocities()[k] = atoms[index].velocity;
    }
    return system;
}

// --- Tilted rotor: every tensor component, at every rank count -------------

Vec3 normalised(Vec3 v) {
    const double n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    for (double& c : v) c /= n;
    return v;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}

// Rodrigues rotation of `v` by `angle` about the unit axis `axis`.
Vec3 rotate_about(const Vec3& axis, double angle, const Vec3& v) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const Vec3 k = cross(axis, v);
    const double kd = axis[0]*v[0] + axis[1]*v[1] + axis[2]*v[2];
    Vec3 out{};
    for (std::size_t d = 0; d < 3; ++d) {
        out[d] = v[d] * c + k[d] * s + axis[d] * kd * (1.0 - c);
    }
    return out;
}

const char* component_name(std::size_t index) {
    static const char* names[9] = {"W_xx","W_xy","W_xz","W_yx","W_yy","W_yz",
                                   "W_zx","W_zy","W_zz"};
    return names[index];
}

// (2, 3, 6) has length exactly 7, so every component of the bond direction --
// and so every entry of u (x) u -- is non-zero and distinct. Centred on x = 10,
// a domain boundary for 2 and 4 ranks, so the pair still straddles two owners.
const Vec3 kTiltDirection = normalised({2.0, 3.0, 6.0});
const Vec3 kTiltAxis = normalised(cross(kTiltDirection, {0.0, 0.0, 1.0}));

std::vector<Atom> tilted_dimer() {
    const double total = kMass0 + kMass1;
    const double s0 = (kMass1 / total) * kBond;
    const double s1 = -(kMass0 / total) * kBond;
    const Vec3 tangent = cross(kTiltAxis, kTiltDirection);
    Atom a{}, b{};
    for (std::size_t d = 0; d < 3; ++d) {
        const double centre = 10.0;
        a.position[d] = centre + s0 * kTiltDirection[d];
        b.position[d] = centre + s1 * kTiltDirection[d];
        a.velocity[d] = kOmega * s0 * tangent[d];
        b.velocity[d] = kOmega * s1 * tangent[d];
    }
    a.mass = kMass0;
    b.mass = kMass1;
    return {a, b};
}

// Every one of the nine components must equal the analytical rigid-rotor value
// and must be identical on every rank. A contribution reduced when it is already
// global shows up as a factor of the rank count in EVERY component, which is
// checked component by component rather than only in the trace.
void check_tilted_rotor_components(int rank, int size, int& failures) {
    gmd::Box box;
    box.set_lengths({kBox, kBox, kBox});
    gmd::DomainDecomposition dd;
    dd.create_decomposition(box, size, rank, 4.0, 1.0, {true, true, true});

    const auto atoms = tilted_dimer();
    std::vector<int> owned;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (dd.owner_rank(box, atoms[i].position) == rank) owned.push_back(static_cast<int>(i));
    }
    gmd::System system;
    system.resize(owned.size(), owned.size());
    system.set_box(box);
    for (std::size_t k = 0; k < owned.size(); ++k) {
        const auto index = static_cast<std::size_t>(owned[k]);
        system.mutable_masses()[k] = atoms[index].mass;
        system.mutable_atom_tags()[k] = owned[k];
        system.mutable_atom_owners()[k] = rank;
        system.mutable_coordinates()[k] = atoms[index].position;
        system.mutable_velocities()[k] = atoms[index].velocity;
    }

    int local_atoms = static_cast<int>(system.num_local_atoms());
    int max_on_one_rank = 0;
    MPI_Allreduce(&local_atoms, &max_on_one_rank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (size > 1) {
        check(max_on_one_rank == 1,
              "the tilted pair must span two ranks; one rank holds " +
                  std::to_string(max_on_one_rank), rank, failures);
    }

    gmd::ConstraintSolver solver({gmd::BondConstraint{0, 1, kBond}}, tight_settings());
    const gmd::ConstraintReference reference = solver.capture_reference(system);
    {
        auto coordinates = system.mutable_coordinates();
        const auto velocities = system.velocities();
        for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) coordinates[i][d] += velocities[i][d] * kTimeStep;
        }
    }
    solver.apply_shake(system, reference, kTimeStep);
    gmd::ConstraintVirialResult result;
    solver.apply_rattle(system, kTimeStep, result);
    check(result.valid, "the tilted rotor must produce a valid constraint virial",
          rank, failures);

    // Analytical endpoint reference: W_ab = -mu omega^2 d^2 u_a u_b, with u the
    // bond direction after the step, which for a free rigid rotor is the initial
    // one rotated by omega*dt about the spin axis.
    const double mu = (kMass0 * kMass1) / (kMass0 + kMass1);
    const Vec3 u = rotate_about(kTiltAxis, kOmega * kTimeStep, kTiltDirection);
    std::array<double, 9> expected{};
    for (std::size_t a = 0; a < 3; ++a) {
        for (std::size_t b = 0; b < 3; ++b) {
            expected[a * 3 + b] = -mu * kOmega * kOmega * kBond * kBond * u[a] * u[b];
        }
    }
    const double scale =
        std::abs(expected[0] + expected[4] + expected[8]);
    // The endpoint estimator's discretisation error is (omega dt)^2 / 4 relative,
    // which is 6.25e-6 here. The bound below is set from that, not fitted.
    const double bound = scale * 1.0e-4;

    for (std::size_t index = 0; index < 9; ++index) {
        check(std::abs(expected[index]) > 0.02 * scale,
              std::string("reference ") + component_name(index) +
                  " must be non-trivial, or this proves nothing", rank, failures);
        check_close(result.virial[index], expected[index], bound,
                    std::string(component_name(index)) +
                        " must match the analytical tilted-rotor value at any rank count",
                    rank, failures);

        // Rank invariance, component by component and exact: every rank runs the
        // same arithmetic over the same allgathered atoms.
        double local = result.virial[index];
        double lo = 0.0;
        double hi = 0.0;
        MPI_Allreduce(&local, &lo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local, &hi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        check(hi == lo,
              std::string("ranks disagree on ") + component_name(index) +
                  " (min " + std::to_string(lo) + ", max " + std::to_string(hi) + ")",
              rank, failures);

        // And name the rank-multiplication failure mode explicitly.
        const double ratio = result.virial[index] / expected[index];
        check(std::abs(ratio - 1.0) < 1.0e-3,
              std::string(component_name(index)) + " is " + std::to_string(ratio) +
                  "x its analytical value; ~" + std::to_string(size) +
                  "x means the already-global constraint virial was reduced across "
                  "ranks, and ~0 means non-owning ranks dropped it",
              rank, failures);
    }

    if (rank == 0) {
        std::cout << "[mpi constraint virial] tilted rotor: all nine components match "
                     "the analytical value on " << size << " rank(s)\n";
    }
}

// A provider with a fixed, deliberately asymmetric virial. Real providers
// allreduce their own virial before returning it, so every rank seeing the same
// tensor is what production looks like.
class FixedVirialForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "fixed_virial"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.potential_energy = 0.0;
        result.virial = {1.0, 0.25, -0.5, 0.25, 2.0, 0.75, -0.5, 0.75, 3.0};
        result.virial_valid = true;
        result.success = true;
    }
};

// The completed-step record the trajectory writer reports must be IDENTICAL on
// every rank, because the MPI output path copies it from whichever rank happens
// to be writing without reducing it. That is an assumption, so it is checked
// rather than left implicit: every field is either replicated (the box, the
// constraint virial, which every rank computes from the same allgathered atoms)
// or already global (the provider virial, which providers allreduce, and 2K,
// which compute_twice_ke allreduces). Bitwise equality is the correct
// expectation, so bitwise equality is what is asserted.
void check_step_thermodynamics_agree_across_ranks(int rank, int size, int& failures) {
    gmd::Box box;
    box.set_lengths({kBox, kBox, kBox});
    gmd::DomainDecomposition dd;
    dd.create_decomposition(box, size, rank, 4.0, 1.0, {true, true, true});

    gmd::System system = local_share(dd, box, rank);
    auto solver = std::make_shared<gmd::ConstraintSolver>(
        std::vector<gmd::BondConstraint>{{0, 1, kBond}}, tight_settings());
    gmd::VelocityVerletIntegrator integrator(kTimeStep);
    integrator.set_constraint_solver(solver);
    FixedVirialForceProvider provider;
    gmd::RuntimeContext runtime;
    integrator.initialize(system, runtime);
    const gmd::IntegratorStepContext ctx{.step = 0, .dt = kTimeStep};
    integrator.step(system, provider, ctx, runtime);

    const auto& completed = system.step_thermodynamics();
    check(completed.valid, "a completed constrained step must record its thermodynamics",
          rank, failures);

    auto agrees = [&](double value, const std::string& what) {
        double lo = 0.0;
        double hi = 0.0;
        MPI_Allreduce(&value, &lo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&value, &hi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        check(hi == lo,
              "ranks disagree on " + what + " (min " + std::to_string(lo) + ", max " +
                  std::to_string(hi) + "); the MPI output path copies this from one rank "
                  "without reducing it, so it has to be identical everywhere",
              rank, failures);
    };
    agrees(completed.pressure, "the completed-step pressure");
    agrees(completed.twice_kinetic_energy, "the completed-step 2K");
    agrees(completed.volume, "the completed-step volume");
    agrees(completed.potential_energy, "the completed-step potential energy");
    for (std::size_t k = 0; k < 9; ++k) {
        agrees(completed.virial[k],
               "completed-step virial component " + std::to_string(k));
    }

    // Non-vacuity: the record must carry a real constraint contribution, or the
    // agreement above would be an agreement about nothing.
    const double provider_trace = 1.0 + 2.0 + 3.0;
    check(std::abs(completed.virial[0] + completed.virial[4] + completed.virial[8] -
                   provider_trace) > 1.0e-6,
          "the completed-step virial must differ from the provider virial alone, or "
          "the constraint term never reached it", rank, failures);
    check(completed.twice_kinetic_energy > 1.0e-9,
          "the fixture must carry kinetic energy", rank, failures);
}

int run(int rank, int size) {
    int failures = 0;

    gmd::Box box;
    box.set_lengths({kBox, kBox, kBox});
    gmd::DomainDecomposition dd;
    dd.create_decomposition(box, size, rank, 4.0, 1.0, {true, true, true});

    gmd::System system = local_share(dd, box, rank);

    // Confirm the premise: the two constrained atoms really are on different
    // ranks whenever there is more than one, otherwise the test proves nothing.
    int local_atoms = static_cast<int>(system.num_local_atoms());
    int total_atoms = 0;
    int max_on_one_rank = 0;
    MPI_Allreduce(&local_atoms, &total_atoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_atoms, &max_on_one_rank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    check(total_atoms == 2, "the fixture must have two atoms in total, got " +
                                std::to_string(total_atoms), rank, failures);
    if (size > 1) {
        check(max_on_one_rank == 1,
              "the constrained pair must span two ranks; one rank holds " +
                  std::to_string(max_on_one_rank),
              rank, failures);
    }

    // One step of the standard SHAKE/RATTLE splitting with no forces acting:
    // capture the reference geometry, drift, SHAKE (position + velocity impulse),
    // then RATTLE, whose multipliers give the endpoint constraint force.
    gmd::ConstraintSolver solver({gmd::BondConstraint{0, 1, kBond}}, tight_settings());
    const gmd::ConstraintReference reference = solver.capture_reference(system);
    {
        auto coordinates = system.mutable_coordinates();
        const auto velocities = system.velocities();
        for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) {
                coordinates[i][d] += velocities[i][d] * kTimeStep;
            }
        }
    }

    solver.apply_shake(system, reference, kTimeStep);
    gmd::ConstraintVirialResult result;
    const auto stats = solver.apply_rattle(system, kTimeStep, result);

    check(stats.converged, "RATTLE must converge", rank, failures);
    check(result.valid, "the constraint virial must be reported valid", rank, failures);
    check(result.contributions.size() == 1,
          "one contribution per constraint regardless of rank count", rank, failures);

    // Centripetal reference: mu = m0 m1 / (m0 + m1).
    const double mu = (kMass0 * kMass1) / (kMass0 + kMass1);
    const double expected_trace = -mu * kOmega * kOmega * kBond * kBond;
    // The discretisation error is second order in dt and identical on every rank;
    // the tolerance below is set from that, not from the observed spread.
    const double discretisation = std::abs(expected_trace) * (kOmega * kTimeStep) *
                                  (kOmega * kTimeStep);

    check_close(result.trace(), expected_trace, discretisation,
                "endpoint tr W must equal -mu omega^2 d^2 at any rank count", rank,
                failures);

    // Spell the two failure modes out, so a regression names itself.
    const double ratio = result.trace() / expected_trace;
    check(std::abs(ratio - 1.0) < 1.0e-3,
          "constraint virial is " + std::to_string(ratio) +
              "x the expected value; ~" + std::to_string(size) +
              "x means it was reduced across ranks when it is already global, "
              "and ~0 means non-owning ranks dropped it",
          rank, failures);

    // Every rank must hold the identical tensor, component by component. This is
    // the actual rank-invariance claim, and it is exact: every rank runs the same
    // arithmetic over the same allgathered data.
    for (std::size_t k = 0; k < 9; ++k) {
        double local = result.virial[k];
        double lo = 0.0;
        double hi = 0.0;
        MPI_Allreduce(&local, &lo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local, &hi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        check(hi == lo,
              "ranks disagree on constraint virial component " + std::to_string(k) +
                  " (min " + std::to_string(lo) + ", max " + std::to_string(hi) + ")",
              rank, failures);
    }

    // The pressure identity must survive decomposition too: 2K is summed over
    // owned atoms across ranks, tr W is already global.
    double local_twice_ke = 0.0;
    {
        const auto masses = system.masses();
        const auto velocities = system.velocities();
        for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) {
                local_twice_ke += masses[i] * velocities[i][d] * velocities[i][d];
            }
        }
    }
    double twice_ke = 0.0;
    MPI_Allreduce(&local_twice_ke, &twice_ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    check(twice_ke > 1.0e-9, "the fixture must carry kinetic energy", rank, failures);
    check(std::abs(twice_ke + result.trace()) < 1.0e-3 * twice_ke,
          "2K + tr W must vanish for a rigid rotor at any rank count; residual "
          "fraction " + std::to_string((twice_ke + result.trace()) / twice_ke),
          rank, failures);

    check_tilted_rotor_components(rank, size, failures);
    check_step_thermodynamics_agree_across_ranks(rank, size, failures);

    return failures;
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int failures = run(rank, size);

    int total_failures = 0;
    MPI_Allreduce(&failures, &total_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (total_failures == 0) {
            std::cout << "[mpi constraint virial] all checks passed on " << size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi constraint virial] " << total_failures
                      << " check(s) failed on " << size << " rank(s)\n";
        }
    }
    return total_failures == 0 ? 0 : 1;
}
