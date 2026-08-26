// Tests for the standard constrained velocity-Verlet (SHAKE/RATTLE) splitting.
//
// The reference is a rigid rotor, not the implementation. A dimer rotating with
// no external force is held together entirely by its constraint, so two things
// follow that can be written down without reference to any solver:
//
//   1. Energy. Constraint forces do no work and there is no potential, so the
//      kinetic energy is exactly constant. This is the sharpest test of the
//      SHAKE projection itself: correcting along the drifted bond instead of the
//      reference gradient still lands on the constraint manifold, but is not
//      symplectic and bleeds energy secularly.
//
//   2. Two half-impulses. The splitting delivers the constraint force twice per
//      step -- through SHAKE, along r_c(t), paired with F(t), and through
//      RATTLE, along r_c(t+dt), paired with F(t+dt). Both are measured here as
//      the velocity change each projection actually produces, and they lie along
//      different bonds, separated by the angle the step turned through.
//
// Both fixtures start ON the constraint manifold in BOTH position and velocity
// (|r_ij| = d and r_ij . v_ij = 0). A state that violates either is not a
// constrained state, and projecting it measures the size of the initial error
// rather than any physical constraint force.

//      of the SHAKE projection itself: correcting along the drifted bond instead
//      of the reference gradient still lands on the manifold, but is not
//      symplectic and bleeds energy secularly.
//
// Every fixture starts ON the constraint manifold in BOTH position and velocity
// (|r_ij| = d and r_ij . v_ij = 0). A state that violates either is not a
// constrained state, and the projection of such a state measures the size of the
// initial error rather than any physical constraint force.
//
// CONVENTIONS (read off the codebase, not assumed):
//   W = sum_i r_i (x) F_i,  W_ab = virial[a*3 + b]  (ClassicalForceProvider)
//   P = (2K + tr W) / 3V                            (BerendsenBarostat, writer)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/barostat.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/integrator.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

using Vec3 = gmd::System::Vec3;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[constraint virial] " << message << '\n';
        ++failures;
    }
}

void check_close(double actual, double expected, double tolerance,
                 const std::string& message) {
    if (!(std::abs(actual - expected) <= tolerance)) {
        std::cerr << "[constraint virial] " << message << ": expected " << std::setprecision(12)
                  << expected << ", got " << actual << " (|diff| "
                  << std::abs(actual - expected) << ", tolerance " << tolerance << ")\n";
        ++failures;
    }
}

gmd::ConstraintSettings tight_settings() {
    gmd::ConstraintSettings settings;
    settings.tolerance = 1.0e-13;
    settings.max_iterations = 1000;
    return settings;
}

gmd::Box cubic_box(double length) {
    gmd::Box box;
    box.set_lengths({length, length, length});
    return box;
}

// A force provider that contributes nothing. Used so that the only force acting
// on a fixture is the constraint force, which is what makes the centripetal
// reference exact.
class NullForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "null"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.potential_energy = 0.0;
        result.virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        result.virial_valid = true;
        result.success = true;
    }
};

// --- Physically valid constrained fixtures ---------------------------------

// A rigid dimer rotating in the xy-plane about its own centre of mass, which is
// placed at `centre` and is at rest. Positions satisfy |r_ij| = d exactly and
// velocities satisfy r_ij . v_ij = 0 exactly, so the state is on the constraint
// manifold in both position and velocity before any solver runs.
gmd::System make_rotating_dimer(double m0, double m1, double d, double omega,
                                const gmd::Box& box, const Vec3& centre) {
    const double total = m0 + m1;
    gmd::System system;
    system.resize(2, 2);
    system.set_box(box);
    auto masses = system.mutable_masses();
    auto tags = system.mutable_atom_tags();
    auto coordinates = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();

    masses[0] = m0;
    masses[1] = m1;
    tags[0] = 0;
    tags[1] = 1;
    // Offsets along x, split so the centre of mass sits at `centre`.
    const double x0 = (m1 / total) * d;
    const double x1 = -(m0 / total) * d;
    coordinates[0] = {centre[0] + x0, centre[1], centre[2]};
    coordinates[1] = {centre[0] + x1, centre[1], centre[2]};
    // v = omega z_hat x (r - centre): purely tangential, zero net momentum.
    velocities[0] = {0.0, omega * x0, 0.0};
    velocities[1] = {0.0, omega * x1, 0.0};
    return system;
}

double distance(const gmd::System& system, std::size_t a, std::size_t b) {
    const auto coordinates = system.coordinates();
    double sum = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        const double delta = coordinates[a][d] - coordinates[b][d];
        sum += delta * delta;
    }
    return std::sqrt(sum);
}

double twice_kinetic_energy(const gmd::System& system) {
    const auto masses = system.masses();
    const auto velocities = system.velocities();
    double sum = 0.0;
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            sum += masses[i] * velocities[i][d] * velocities[i][d];
        }
    }
    return sum;
}

// Confirms a fixture really is a valid constrained state before it is used.
void require_on_manifold(const gmd::System& system,
                         const std::vector<gmd::BondConstraint>& constraints,
                         const std::string& label) {
    const auto coordinates = system.coordinates();
    const auto velocities = system.velocities();
    for (const auto& constraint : constraints) {
        const auto i = static_cast<std::size_t>(constraint.i);
        const auto j = static_cast<std::size_t>(constraint.j);
        double dot = 0.0;
        for (std::size_t d = 0; d < 3; ++d) {
            dot += (coordinates[i][d] - coordinates[j][d]) *
                   (velocities[i][d] - velocities[j][d]);
        }
        check_close(distance(system, i, j), constraint.target_distance, 1.0e-13,
                    label + ": fixture positions must start on the constraint manifold");
        check(std::abs(dot) < 1.0e-13,
              label + ": fixture velocities must start tangent to the constraint "
                      "surface, r_ij . v_ij = " + std::to_string(dot));
    }
}

// ---------------------------------------------------------------------------
// 1b. The two constraint half-impulses are distinct quantities.
// ---------------------------------------------------------------------------
//
// This is what makes "endpoint" mean something. The splitting delivers the
// constraint force twice per step: through SHAKE, along r_c(t), paired with
// F(t); and through RATTLE, along r_c(t+dt), paired with F(t+dt). Both are
// measured here as the actual velocity change each projection produces, and the
// test shows they point along different bonds, separated by exactly the angle
// the bond turned through during the step. The virial must use the second.
void test_two_constraint_half_impulses_are_distinct() {
    constexpr double m0 = 1.0;
    constexpr double m1 = 3.0;
    constexpr double bond = 1.5;
    constexpr double omega = 0.1;
    constexpr double dt = 0.4;              // coarse on purpose: the two differ by omega*dt
    const double mu = (m0 * m1) / (m0 + m1);

    const gmd::Box box = cubic_box(30.0);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};
    gmd::System system = make_rotating_dimer(m0, m1, bond, omega, box, {15.0, 15.0, 15.0});
    require_on_manifold(system, constraints, "half-impulse dimer");

    gmd::ConstraintSolver solver(constraints, tight_settings());

    auto relative_velocity = [&]() {
        const auto v = system.velocities();
        return Vec3{v[0][0] - v[1][0], v[0][1] - v[1][1], v[0][2] - v[1][2]};
    };
    auto bond_vector = [&]() {
        const auto r = system.coordinates();
        return Vec3{r[0][0] - r[1][0], r[0][1] - r[1][1], r[0][2] - r[1][2]};
    };
    auto normalised = [](Vec3 v) {
        double n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        for (double& c : v) c /= n;
        return v;
    };
    // Impulse on atom 0 from a relative-velocity change: J_0 = mu * d(v_0 - v_1).
    auto impulse_from = [&](const Vec3& before, const Vec3& after) {
        return Vec3{mu * (after[0] - before[0]),
                    mu * (after[1] - before[1]),
                    mu * (after[2] - before[2])};
    };

    const Vec3 reference_bond = bond_vector();

    // Steps (1)-(3): no forces, so the first half-kick is a no-op; drift, then
    // the SHAKE projection with its velocity impulse.
    const gmd::ConstraintReference reference = solver.capture_reference(system);
    {
        auto coordinates = system.mutable_coordinates();
        const auto velocities = system.velocities();
        for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) coordinates[i][d] += velocities[i][d] * dt;
        }
    }
    const Vec3 before_shake = relative_velocity();
    solver.apply_shake(system, reference, dt);
    const Vec3 after_shake = relative_velocity();
    const Vec3 shake_impulse = impulse_from(before_shake, after_shake);

    const Vec3 endpoint_bond = bond_vector();

    // Steps (4)-(5): no forces, so the second half-kick is a no-op; RATTLE.
    gmd::ConstraintVirialResult result;
    solver.apply_rattle(system, dt, result);
    const Vec3 after_rattle = relative_velocity();
    const Vec3 rattle_impulse = impulse_from(after_shake, after_rattle);

    auto magnitude = [](const Vec3& v) {
        return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    };
    const double shake_force = 2.0 * magnitude(shake_impulse) / dt;
    const double endpoint_force = 2.0 * magnitude(rattle_impulse) / dt;
    const double step_mean_force =
        magnitude(Vec3{shake_impulse[0] + rattle_impulse[0],
                       shake_impulse[1] + rattle_impulse[1],
                       shake_impulse[2] + rattle_impulse[2]}) / dt;

    std::cout << "[constraint virial] the step's two constraint half-impulses at dt="
              << dt << ":\n"
              << "    SHAKE  (t,    pairs with F(t))     |G| = "
              << std::setprecision(10) << shake_force << '\n'
              << "    RATTLE (t+dt, pairs with F(t+dt))  |G| = " << endpoint_force
              << "   <- the virial uses this one\n"
              << "    step mean (their sum over dt)      |G| = " << step_mean_force << '\n';

    check(shake_force > 1.0e-6 && endpoint_force > 1.0e-6,
          "both half-impulses must be non-trivial, or the comparison is vacuous");

    // Each impulse lies along the bond of ITS OWN time level.
    const Vec3 r0 = normalised(reference_bond);
    const Vec3 r1 = normalised(endpoint_bond);
    const Vec3 s_hat = normalised(shake_impulse);
    const Vec3 r_hat = normalised(rattle_impulse);
    auto dot3 = [](const Vec3& a, const Vec3& b) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    check_close(std::abs(dot3(s_hat, r0)), 1.0, 1.0e-12,
                "the SHAKE impulse must lie along the REFERENCE bond r_c(t)");
    check_close(std::abs(dot3(r_hat, r1)), 1.0, 1.0e-12,
                "the RATTLE impulse must lie along the ENDPOINT bond r_c(t+dt)");

    // ...and those bonds differ by the angle the step turned through, which is
    // exactly why the two time levels are not interchangeable.
    const double turned = std::acos(std::min(1.0, std::abs(dot3(r0, r1))));
    std::cout << "    bond turned through " << std::setprecision(4) << turned
              << " rad during the step (omega*dt = " << omega * dt << ")\n";
    check(turned > 0.5 * omega * dt,
          "the fixture must actually rotate, or endpoint and reference coincide");
    check(std::abs(dot3(s_hat, r_hat)) < 1.0 - 1.0e-6,
          "the two half-impulses must not be parallel, or the test proves nothing");

    // And the virial is built from the endpoint one.
    check(result.contributions.size() == 1, "one contribution");
    const auto& contribution = result.contributions.front();
    for (std::size_t k = 0; k < 3; ++k) {
        check_close(contribution.force[k], 2.0 * rattle_impulse[k] / dt,
                    1.0e-12 * std::max(1.0, endpoint_force),
                    "the recorded force must be 2/dt times the RATTLE impulse "
                    "(component " + std::to_string(k) + ")");
    }
}

// ---------------------------------------------------------------------------
// 1c. Constrained NVE: a free rigid rotor must not lose energy.
// ---------------------------------------------------------------------------
//
// Constraint forces do no work and there is no potential, so the kinetic energy
// of a free rigid rotor is exactly constant. This is the test that pins the
// SHAKE projection direction: displacing along the drifted bond rather than
// along the reference gradient also lands on the manifold, but is not symplectic
// and drains energy secularly -- measurably, ~86% over 4000 steps at the coarse
// timestep used here.
void test_free_rotor_conserves_energy() {
    constexpr double m0 = 1.0;
    constexpr double m1 = 3.0;
    constexpr double bond = 1.5;
    constexpr double omega = 0.1;
    constexpr int steps = 4000;

    const gmd::Box box = cubic_box(60.0);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};

    std::cout << "[constraint virial] constrained NVE free rotor, " << steps << " steps\n";

    for (const double dt : {0.4, 0.2, 0.1}) {
        gmd::System system =
            make_rotating_dimer(m0, m1, bond, omega, box, {30.0, 30.0, 30.0});
        require_on_manifold(system, constraints, "NVE rotor");

        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        gmd::RuntimeContext runtime;
        integrator.initialize(system, runtime);

        const double initial = twice_kinetic_energy(system);
        double worst_excursion = 0.0;
        double worst_bond = 0.0;
        for (int step = 0; step < steps; ++step) {
            const gmd::IntegratorStepContext ctx{
                .step = static_cast<std::uint64_t>(step), .dt = dt};
            integrator.step(system, provider, ctx, runtime);
            worst_excursion = std::max(
                worst_excursion, std::abs(twice_kinetic_energy(system) - initial) / initial);
            worst_bond = std::max(worst_bond, std::abs(distance(system, 0, 1) - bond));
        }
        const double final_drift = (twice_kinetic_energy(system) - initial) / initial;

        std::cout << "    dt=" << std::setw(5) << dt
                  << "  final relative energy drift " << std::setprecision(3) << final_drift
                  << "  max excursion " << worst_excursion
                  << "  max bond error " << worst_bond << '\n';

        // Exact conservation up to round-off accumulated over 4000 steps. The
        // bound is not fitted: it is the drifting-projection failure mode being
        // absent, which would show up here at the 1e-1 level.
        check(std::abs(final_drift) < 1.0e-9,
              "a free rigid rotor must conserve energy exactly; relative drift over " +
                  std::to_string(steps) + " steps at dt=" + std::to_string(dt) + " was " +
                  std::to_string(final_drift) +
                  " (a large negative value means SHAKE is projecting along the drifted "
                  "bond instead of the reference gradient)");
        check(worst_excursion < 1.0e-9,
              "and must not excurse either; worst was " + std::to_string(worst_excursion));
        check(worst_bond < 1.0e-11,
              "the constraint must hold throughout; worst bond error " +
                  std::to_string(worst_bond));
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

    test_two_constraint_half_impulses_are_distinct();
    test_free_rotor_conserves_energy();

    if (failures != 0) {
        std::cerr << "[constraint virial] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[constraint virial] all checks passed\n";
    return 0;
}
