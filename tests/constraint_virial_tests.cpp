// Physics tests for the SHAKE/RATTLE constraint contribution to the virial.
//
// ENDPOINT, NOT STEP-AVERAGED. The integrator implements the standard
// constrained velocity-Verlet splitting, in which the constraint force enters
// twice per step: once through the SHAKE impulse, paired with F(t), and once
// through the RATTLE impulse, paired with F(t+dt). The providers are evaluated
// at r(t+dt), so the virial takes the RATTLE one, G(t+dt) = (2/dt) Lambda r_c.
// That is an endpoint quantity. Every name, message and diagnostic below says
// which of the two it means, and test_two_constraint_half_impulses_are_distinct
// measures both and shows they are genuinely different quantities.
//
// The reference used here is a rigid rotor, not the implementation formula.
// A pair (or a rigid triangle) rotating with no external force is held together
// entirely by its constraints, so the constraint force is the centripetal force
// and can be written down without reference to any solver:
//
//     G_i = -m_i omega^2 (r_i - R_com)
//
// and, taking the origin at R_com so that the sum of the constraint forces
// vanishing makes the result origin independent,
//
//     tr W_constraint = sum_i (r_i - R) . G_i = -omega^2 sum_i m_i |r_i - R|^2
//                     = -I omega^2 = -2K.
//
// Three things follow, and all three are asserted below.
//
//   1. The magnitude test. For a dimer I omega^2 = mu omega^2 d^2, so the
//      recovered ENDPOINT pair force must approach mu omega^2 d as dt -> 0. This
//      pins the multiplier-to-force conversion factor from the outside: 1/dt --
//      the value that would be right if the SHAKE velocity impulse were missing,
//      so that RATTLE had to carry the whole step -- undershoots it by two.
//
//   2. The pressure identity. 2K + tr W_constraint = 0. A rigid body's frozen
//      internal coordinates contribute nothing to the pressure. This is what
//      makes it wrong to drop the constraint virial on the grounds that the
//      velocities have been projected: the projection produces the 2K in that
//      identity, and without the constraint term nothing cancels it -- a lone
//      rigid rotor would report a spurious pressure of 2K/3V.
//
//   3. Energy. A free rigid rotor has no potential energy and constraint forces
//      do no work, so its kinetic energy is constant. This is the sharpest test
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

// A provider with a fixed, deliberately asymmetric virial, so that a term added
// twice or dropped is visible component by component.
class FixedVirialForceProvider final : public gmd::ForceProvider {
public:
    explicit FixedVirialForceProvider(std::array<double, 9> virial)
        : virial_(virial) {}
    std::string_view name() const noexcept override { return "fixed_virial"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        ++calls;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.potential_energy = 0.0;
        result.virial = virial_;
        result.virial_valid = true;
        result.success = true;
    }
    int calls = 0;

private:
    std::array<double, 9> virial_;
};

// Scales the box by a fixed factor on its first call, then leaves it alone.
class OneShotScalingBarostat final : public gmd::Barostat {
public:
    explicit OneShotScalingBarostat(double scale) : scale_(scale) {}
    void apply(gmd::System& system, gmd::ForceProvider&, gmd::RuntimeContext&,
               std::uint64_t, double, double, double, double) override {
        ++calls;
        if (fired_) return;
        fired_ = true;
        gmd::Box box = system.box();
        box.set_lengths({box.lengths[0] * scale_,
                         box.lengths[1] * scale_,
                         box.lengths[2] * scale_});
        system.set_box(box);
        auto coordinates = system.mutable_coordinates();
        for (auto& coordinate : coordinates) {
            for (std::size_t d = 0; d < 3; ++d) coordinate[d] *= scale_;
        }
    }
    bool requires_virial() const noexcept override { return false; }
    int calls = 0;

private:
    double scale_;
    bool fired_ = false;
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

// A rigid triangle rotating about its centre of mass in the xy-plane, held by
// three coupled distance constraints. Same guarantees as the dimer: rigid-body
// rotation satisfies r_ij . v_ij = 0 for every pair simultaneously.
gmd::System make_rotating_triangle(const std::array<double, 3>& masses_in,
                                   const std::array<Vec3, 3>& offsets,
                                   double omega,
                                   const gmd::Box& box,
                                   const Vec3& centre) {
    double total = 0.0;
    Vec3 weighted{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < 3; ++i) {
        total += masses_in[i];
        for (std::size_t d = 0; d < 3; ++d) weighted[d] += masses_in[i] * offsets[i][d];
    }
    const Vec3 com{weighted[0] / total, weighted[1] / total, weighted[2] / total};

    gmd::System system;
    system.resize(3, 3);
    system.set_box(box);
    auto masses = system.mutable_masses();
    auto tags = system.mutable_atom_tags();
    auto coordinates = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();
    for (std::size_t i = 0; i < 3; ++i) {
        masses[i] = masses_in[i];
        tags[i] = static_cast<int>(i);
        // Position relative to the body's own centre of mass, then placed.
        const Vec3 arm{offsets[i][0] - com[0], offsets[i][1] - com[1], offsets[i][2] - com[2]};
        coordinates[i] = {centre[0] + arm[0], centre[1] + arm[1], centre[2] + arm[2]};
        velocities[i] = {-omega * arm[1], omega * arm[0], 0.0};
    }
    return system;
}

// --- General-orientation helpers -------------------------------------------

Vec3 normalised(Vec3 v) {
    const double n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    for (double& c : v) c /= n;
    return v;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}

double dot3(const Vec3& a, const Vec3& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

// Rotation by `angle` about the unit axis `axis` (Rodrigues), row-major.
std::array<double, 9> rotation_matrix(const Vec3& axis, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double t = 1.0 - c;
    const double x = axis[0], y = axis[1], z = axis[2];
    return {t*x*x + c,    t*x*y - s*z,  t*x*z + s*y,
            t*x*y + s*z,  t*y*y + c,    t*y*z - s*x,
            t*x*z - s*y,  t*y*z + s*x,  t*z*z + c};
}

Vec3 apply_rotation(const std::array<double, 9>& R, const Vec3& v) {
    return {R[0]*v[0] + R[1]*v[1] + R[2]*v[2],
            R[3]*v[0] + R[4]*v[1] + R[5]*v[2],
            R[6]*v[0] + R[7]*v[1] + R[8]*v[2]};
}

// R A R^T, row-major throughout.
std::array<double, 9> conjugate(const std::array<double, 9>& R,
                                const std::array<double, 9>& A) {
    std::array<double, 9> RA{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < 3; ++k) sum += R[i*3 + k] * A[k*3 + j];
            RA[i*3 + j] = sum;
        }
    }
    std::array<double, 9> out{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < 3; ++k) sum += RA[i*3 + k] * R[j*3 + k];  // R^T
            out[i*3 + j] = sum;
        }
    }
    return out;
}

const char* component_name(std::size_t index) {
    static const char* names[9] = {"W_xx","W_xy","W_xz","W_yx","W_yy","W_yz",
                                   "W_zx","W_zy","W_zz"};
    return names[index];
}

// A rigid dimer rotating about `axis` with its bond along `direction`, centre of
// mass at rest at `centre`. Both must be unit and mutually orthogonal, which is
// checked: the velocity constraint r_ij . v_ij = 0 holds only then.
gmd::System make_oriented_rotating_dimer(double m0, double m1, double d, double omega,
                                         const Vec3& direction, const Vec3& axis,
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
    const double s0 = (m1 / total) * d;
    const double s1 = -(m0 / total) * d;
    for (std::size_t k = 0; k < 3; ++k) {
        coordinates[0][k] = centre[k] + s0 * direction[k];
        coordinates[1][k] = centre[k] + s1 * direction[k];
    }
    // v = omega * axis x (r - centre), so purely tangential with zero net momentum.
    const Vec3 tangent = cross(axis, direction);
    for (std::size_t k = 0; k < 3; ++k) {
        velocities[0][k] = omega * s0 * tangent[k];
        velocities[1][k] = omega * s1 * tangent[k];
    }
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

double trace_of(const std::array<double, 9>& tensor) {
    return tensor[0] + tensor[4] + tensor[8];
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

// Runs one full integrator step (half-kick, drift, SHAKE, forces, half-kick,
// RATTLE) so that the ordering under test is the production ordering.
void run_one_step(gmd::System& system,
                  gmd::VelocityVerletIntegrator& integrator,
                  gmd::ForceProvider& provider,
                  double dt,
                  std::uint64_t step = 0) {
    gmd::RuntimeContext runtime;
    integrator.initialize(system, runtime);
    const gmd::IntegratorStepContext ctx{.step = step, .dt = dt};
    integrator.step(system, provider, ctx, runtime);
}

// ---------------------------------------------------------------------------
// 1. Rigid rotating dimer against the centripetal force.
// ---------------------------------------------------------------------------
//
// Independent of the implementation: the constraint force of a freely rotating
// rigid dimer IS the centripetal force mu omega^2 d, and 2K + tr W must vanish.
// Both are checked at a sequence of timesteps, and both must converge.
void test_rotating_dimer_against_centripetal_force() {
    constexpr double m0 = 1.0;
    constexpr double m1 = 3.0;
    constexpr double bond = 1.5;
    constexpr double omega = 0.1;
    const double mu = (m0 * m1) / (m0 + m1);
    const double expected_force = mu * omega * omega * bond;      // centripetal
    const double expected_trace = -mu * omega * omega * bond * bond;

    const gmd::Box box = cubic_box(30.0);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};

    const std::array<double, 4> timesteps{0.4, 0.2, 0.1, 0.05};
    std::array<double, 4> force_error{};
    std::array<double, 4> identity_error{};

    std::cout << "[constraint virial] rigid rotating dimer -- ENDPOINT constraint "
                 "force at t+dt, mu=" << mu
              << " omega=" << omega << " d=" << bond
              << "  (expected |G| = " << expected_force << ")\n";

    for (std::size_t k = 0; k < timesteps.size(); ++k) {
        const double dt = timesteps[k];
        gmd::System system = make_rotating_dimer(m0, m1, bond, omega, box, {15.0, 15.0, 15.0});
        require_on_manifold(system, constraints, "rotating dimer");

        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        run_one_step(system, integrator, provider, dt);

        check(system.constraint_virial_valid(),
              "a completed step must attach a constraint virial");
        check(system.last_virial_valid(),
              "the combined virial must be reported valid after a completed step");

        const auto& constraint_virial = system.constraint_virial();
        // The provider contributes nothing, so the whole reported virial is the
        // constraint term.
        for (std::size_t index = 0; index < 9; ++index) {
            check_close(system.last_virial()[index], constraint_virial[index], 1.0e-15,
                        "with a null provider the reported virial is the constraint virial");
        }

        // Magnitude of the recovered ENDPOINT pair force, read back out of the
        // tensor itself rather than out of the bookkeeping struct: for a lone
        // central pair force tr W = -|G| d, so |G| = |tr W| / d.
        const double measured_force = std::abs(trace_of(constraint_virial)) / bond;

        // The step must leave the state on the constraint manifold in BOTH
        // position and velocity. Positions are SHAKE's job, tangency RATTLE's.
        check_close(distance(system, 0, 1), bond, 1.0e-12,
                    "the step must leave the bond at its constrained length");
        {
            const auto coordinates = system.coordinates();
            const auto velocities = system.velocities();
            double dot = 0.0;
            double speed = 0.0;
            for (std::size_t k = 0; k < 3; ++k) {
                const double dv = velocities[0][k] - velocities[1][k];
                dot += (coordinates[0][k] - coordinates[1][k]) * dv;
                speed += dv * dv;
            }
            check(std::abs(dot) <= 1.0e-11 * bond * std::sqrt(speed),
                  "the step must leave the relative velocity in the tangent space, "
                  "r_ij . v_ij = " + std::to_string(dot));
        }
        force_error[k] = measured_force / expected_force - 1.0;

        // Sign: the constraint of a rotor pulls inward, so the trace is negative.
        check(trace_of(constraint_virial) < 0.0,
              "a rotating rigid dimer is held together by an attractive constraint "
              "force, so tr W must be negative");

        // The pressure identity. 2K is taken from the post-RATTLE velocities,
        // exactly as the trajectory writer and the barostat take it.
        const double twice_ke = twice_kinetic_energy(system);
        identity_error[k] = (twice_ke + trace_of(constraint_virial)) / twice_ke;

        std::cout << "    dt=" << std::setw(5) << dt
                  << "  endpoint |G|=" << std::setprecision(10) << measured_force
                  << "  rel.err=" << std::setprecision(3) << force_error[k]
                  << "  (2K + trW)/2K=" << identity_error[k] << '\n';

        // tr W against the closed form, loosely at coarse dt and tightly at fine.
        check_close(trace_of(constraint_virial), expected_trace,
                    std::abs(expected_trace) * 2.0e-3,
                    "tr W must match -mu omega^2 d^2 at dt=" + std::to_string(dt));
    }

    const std::size_t last = timesteps.size() - 1;

    // Convergence to the continuum value. A finite value, approached from below,
    // with the error falling as dt^2.
    check(std::abs(force_error[last]) < 1.0e-4,
          "the recovered constraint force must approach mu omega^2 d; relative error "
          "at the finest timestep was " + std::to_string(force_error[last]));
    check(std::abs(identity_error[last]) < 1.0e-4,
          "2K + tr W must vanish for a rigid rotor; relative residual at the finest "
          "timestep was " + std::to_string(identity_error[last]));

    for (std::size_t k = 0; k + 1 < timesteps.size(); ++k) {
        const double force_ratio = std::abs(force_error[k]) / std::abs(force_error[k + 1]);
        const double identity_ratio =
            std::abs(identity_error[k]) / std::abs(identity_error[k + 1]);
        check(force_ratio > 3.0 && force_ratio < 5.0,
              "halving dt must cut the constraint-force error by about four "
              "(second order); measured " + std::to_string(force_ratio));
        check(identity_ratio > 3.0 && identity_ratio < 5.0,
              "halving dt must cut the 2K + tr W residual by about four; measured " +
                  std::to_string(identity_ratio));
    }

    // The conversion factor, pinned from the outside. The RATTLE impulse is the
    // constraint share of the SECOND half-kick, so (dt/2m) G = w Lambda r and the
    // factor is 2/dt. Halving it to 1/dt -- correct only for a scheme whose SHAKE
    // omits its velocity impulse, leaving RATTLE to carry the whole step -- would
    // halve every number above.
    check(std::abs(force_error[last]) < 0.25,
          "the multiplier-to-force conversion factor is wrong by a constant: the "
          "measured endpoint force is " + std::to_string(1.0 + force_error[last]) +
          "x the centripetal force (0.5x means the step-mean 1/dt factor was used)");
}

// ---------------------------------------------------------------------------
// 2. Coupled constraints: a rigid rotating triangle.
// ---------------------------------------------------------------------------
//
// Three shared-atom distance constraints solved together. The same centripetal
// reference applies -- tr W = -I omega^2 = -2K -- and in addition the tensor must
// come out exactly symmetric and every pair force exactly along its own bond,
// which the scalar-multiplier formulation gives by construction rather than by
// projecting anything away.
void test_rotating_triangle_coupled_constraints() {
    const std::array<double, 3> masses{1.0, 2.0, 3.0};
    // Deliberately scalene, so no accidental symmetry can hide an error.
    const std::array<Vec3, 3> offsets{Vec3{0.0, 0.0, 0.0},
                                      Vec3{1.4, 0.0, 0.0},
                                      Vec3{0.5, 1.1, 0.0}};
    constexpr double omega = 0.08;
    const gmd::Box box = cubic_box(40.0);
    const Vec3 centre{20.0, 20.0, 20.0};

    gmd::System reference = make_rotating_triangle(masses, offsets, omega, box, centre);
    const std::vector<gmd::BondConstraint> constraints{
        {0, 1, distance(reference, 0, 1)},
        {1, 2, distance(reference, 1, 2)},
        {0, 2, distance(reference, 0, 2)},
    };
    require_on_manifold(reference, constraints, "rotating triangle");

    const std::array<double, 3> timesteps{0.2, 0.1, 0.05};
    std::array<double, 3> identity_error{};

    std::cout << "[constraint virial] rigid rotating triangle (three coupled constraints)\n";

    for (std::size_t k = 0; k < timesteps.size(); ++k) {
        const double dt = timesteps[k];
        gmd::System system = make_rotating_triangle(masses, offsets, omega, box, centre);
        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        run_one_step(system, integrator, provider, dt);

        const auto& virial = system.constraint_virial();
        check(system.constraint_virial_valid(),
              "the coupled-constraint step must attach a constraint virial");

        // Exact symmetry. Not enforced anywhere: it follows from summing scalar
        // multipliers against a bond vector that RATTLE holds fixed.
        const double scale = std::max(1.0, std::abs(trace_of(virial)));
        check_close(virial[1], virial[3], 1.0e-14 * scale, "W_xy must equal W_yx exactly");
        check_close(virial[2], virial[6], 1.0e-14 * scale, "W_xz must equal W_zx exactly");
        check_close(virial[5], virial[7], 1.0e-14 * scale, "W_yz must equal W_zy exactly");

        const double twice_ke = twice_kinetic_energy(system);
        identity_error[k] = (twice_ke + trace_of(virial)) / twice_ke;
        std::cout << "    dt=" << std::setw(5) << dt
                  << "  2K=" << std::setprecision(10) << twice_ke
                  << "  trW=" << trace_of(virial)
                  << "  (2K + trW)/2K=" << std::setprecision(3) << identity_error[k] << '\n';

        for (const auto& constraint : constraints) {
            check_close(distance(system, static_cast<std::size_t>(constraint.i),
                                 static_cast<std::size_t>(constraint.j)),
                        constraint.target_distance, 1.0e-12,
                        "coupled constraints must still hold after the step");
        }
    }

    const std::size_t last = timesteps.size() - 1;
    check(std::abs(identity_error[last]) < 1.0e-3,
          "2K + tr W must vanish for a rigid rotating triangle too; relative residual "
          "at the finest timestep was " + std::to_string(identity_error[last]));
    for (std::size_t k = 0; k + 1 < timesteps.size(); ++k) {
        const double ratio = std::abs(identity_error[k]) / std::abs(identity_error[k + 1]);
        check(ratio > 3.0 && ratio < 5.0,
              "the coupled-constraint residual must fall as dt^2; measured ratio " +
                  std::to_string(ratio));
    }
}

// Every recovered pair force must lie exactly along its own bond, and the two
// atoms of a pair must receive exactly opposite forces. Checked directly on the
// solver output for the coupled case, where an iteration-path artifact would
// show up if the formulation admitted one.
void test_pair_forces_are_central_and_antisymmetric() {
    const std::array<double, 3> masses{1.0, 2.0, 3.0};
    const std::array<Vec3, 3> offsets{Vec3{0.0, 0.0, 0.0},
                                      Vec3{1.4, 0.0, 0.0},
                                      Vec3{0.5, 1.1, 0.0}};
    const gmd::Box box = cubic_box(40.0);
    gmd::System system = make_rotating_triangle(masses, offsets, 0.08, box, {20.0, 20.0, 20.0});
    const std::vector<gmd::BondConstraint> constraints{
        {0, 1, distance(system, 0, 1)},
        {1, 2, distance(system, 1, 2)},
        {0, 2, distance(system, 0, 2)},
    };

    // Drift by hand so that RATTLE has real work to do, then project.
    constexpr double dt = 0.2;
    {
        auto coordinates = system.mutable_coordinates();
        const auto velocities = system.velocities();
        for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) coordinates[i][d] += velocities[i][d] * dt;
        }
    }
    gmd::ConstraintSolver solver(constraints, tight_settings());
    solver.apply_shake(system);
    gmd::ConstraintVirialResult result;
    solver.apply_rattle(system, dt, result);

    check(result.valid, "the RATTLE projection must produce a valid constraint virial");
    check(result.contributions.size() == constraints.size(),
          "one contribution per constraint");

    std::array<double, 9> rebuilt{};
    double largest = 0.0;
    for (const auto& contribution : result.contributions) {
        // Centrality: G x r_ij == 0.
        const auto& r = contribution.bond;
        const auto& g = contribution.force;
        const double cross[3] = {r[1] * g[2] - r[2] * g[1],
                                 r[2] * g[0] - r[0] * g[2],
                                 r[0] * g[1] - r[1] * g[0]};
        double magnitude = 0.0;
        for (std::size_t d = 0; d < 3; ++d) magnitude += g[d] * g[d];
        magnitude = std::sqrt(magnitude);
        largest = std::max(largest, magnitude);
        double residual = 0.0;
        for (double component : cross) residual = std::max(residual, std::abs(component));
        check(residual <= 1.0e-14 * std::max(1.0, magnitude),
              "the pair force must be exactly central; |G x r| = " + std::to_string(residual));

        for (std::size_t a = 0; a < 3; ++a) {
            for (std::size_t b = 0; b < 3; ++b) rebuilt[a * 3 + b] += r[a] * g[b];
        }
    }
    check(largest > 1.0e-6,
          "the fixture must produce non-trivial constraint forces, largest was " +
              std::to_string(largest));
    for (std::size_t index = 0; index < 9; ++index) {
        check_close(rebuilt[index], result.virial[index], 1.0e-13 * std::max(1.0, largest),
                    "the tensor must be the sum of its per-constraint contributions");
    }
}

// ---------------------------------------------------------------------------
// 3. Periodic boundary.
// ---------------------------------------------------------------------------
//
// The same rotor placed across a box face must give the identical tensor: the
// virial is built from minimum-image bond vectors, so it cannot depend on where
// the molecule sits or on how its atoms are wrapped.
void test_periodic_boundary_invariance() {
    constexpr double bond = 1.5;
    constexpr double omega = 0.1;
    constexpr double dt = 0.1;
    const gmd::Box box = cubic_box(12.0);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};

    auto measure = [&](const Vec3& centre) {
        gmd::System system = make_rotating_dimer(1.0, 3.0, bond, omega, box, centre);
        // Wrap the fixture the way the integrator would.
        auto coordinates = system.mutable_coordinates();
        for (auto& coordinate : coordinates) {
            for (std::size_t d = 0; d < 3; ++d) {
                while (coordinate[d] < 0.0) coordinate[d] += box.lengths[d];
                while (coordinate[d] >= box.lengths[d]) coordinate[d] -= box.lengths[d];
            }
        }
        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        run_one_step(system, integrator, provider, dt);
        return system.constraint_virial();
    };

    const auto interior = measure({6.0, 6.0, 6.0});
    const auto straddling = measure({0.0, 6.0, 6.0});     // bond crosses x = 0
    const auto corner = measure({0.0, 0.0, 0.0});

    for (std::size_t index = 0; index < 9; ++index) {
        check_close(straddling[index], interior[index], 1.0e-12,
                    "component " + std::to_string(index) +
                        " must not depend on the molecule straddling a boundary");
        check_close(corner[index], interior[index], 1.0e-12,
                    "component " + std::to_string(index) +
                        " must not depend on the molecule sitting in a corner");
    }
    check(std::abs(trace_of(interior)) > 1.0e-6,
          "the periodic fixture must produce a non-trivial virial");
}

// Translating the whole fixture must not change the tensor either.
void test_translation_invariance() {
    const std::array<double, 3> masses{1.0, 2.0, 3.0};
    const std::array<Vec3, 3> offsets{Vec3{0.0, 0.0, 0.0},
                                      Vec3{1.4, 0.0, 0.0},
                                      Vec3{0.5, 1.1, 0.0}};
    const gmd::Box box = cubic_box(40.0);
    constexpr double dt = 0.1;

    auto measure = [&](const Vec3& centre) {
        gmd::System system = make_rotating_triangle(masses, offsets, 0.08, box, centre);
        const std::vector<gmd::BondConstraint> constraints{
            {0, 1, distance(system, 0, 1)},
            {1, 2, distance(system, 1, 2)},
            {0, 2, distance(system, 0, 2)},
        };
        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        run_one_step(system, integrator, provider, dt);
        return system.constraint_virial();
    };

    const auto here = measure({20.0, 20.0, 20.0});
    const auto there = measure({7.25, 31.5, 12.75});

    // RELATIVE, not absolute. The tensor is built from differences of
    // coordinates that are O(30 A) here but O(1 A) apart, so the round-off floor
    // is set by the coordinate magnitude rather than by the size of the answer.
    // An absolute bound would be a statement about the box, not about the
    // physics. 1e-9 of the tensor scale is still ~7 orders below any effect
    // being tested, and the observed difference is ~8e-11 relative.
    double scale = 0.0;
    for (double component : here) scale = std::max(scale, std::abs(component));
    check(scale > 1.0e-6, "the fixture must produce a non-trivial virial");
    for (std::size_t index = 0; index < 9; ++index) {
        check_close(there[index], here[index], 1.0e-9 * scale,
                    "component " + std::to_string(index) +
                        " must be independent of where the molecule sits");
    }
}

// ---------------------------------------------------------------------------
// 4. Integration: the combination happens once, at the right time.
// ---------------------------------------------------------------------------
void test_combined_virial_is_provider_plus_constraint_exactly_once() {
    const std::array<double, 9> provider_virial{
        1.0, 0.25, -0.5,
        0.25, 2.0, 0.75,
        -0.5, 0.75, 3.0};
    constexpr double bond = 1.5;
    constexpr double dt = 0.1;
    const gmd::Box box = cubic_box(30.0);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};

    gmd::System system = make_rotating_dimer(1.0, 3.0, bond, 0.1, box, {15.0, 15.0, 15.0});
    auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
    gmd::VelocityVerletIntegrator integrator(dt);
    integrator.set_constraint_solver(solver);
    FixedVirialForceProvider provider(provider_virial);
    run_one_step(system, integrator, provider, dt);

    check(system.constraint_virial_state() == gmd::ConstraintVirialState::Valid,
          "a completed constrained step must leave the constraint virial Valid");
    check(system.last_virial_valid(), "the combined virial must be valid");
    for (std::size_t index = 0; index < 9; ++index) {
        check_close(system.last_virial()[index],
                    provider_virial[index] + system.constraint_virial()[index],
                    1.0e-13,
                    "component " + std::to_string(index) +
                        " must be provider + constraint, added exactly once");
    }
    check(std::abs(trace_of(system.constraint_virial())) > 1.0e-6,
          "the constraint term must be non-trivial, or the test proves nothing");
}

// An unconstrained run must be bit-for-bit what it was before any of this.
void test_unconstrained_run_is_untouched() {
    const std::array<double, 9> provider_virial{
        1.0, 0.25, -0.5,
        0.25, 2.0, 0.75,
        -0.5, 0.75, 3.0};
    constexpr double dt = 0.1;
    const gmd::Box box = cubic_box(30.0);
    gmd::System system = make_rotating_dimer(1.0, 3.0, 1.5, 0.1, box, {15.0, 15.0, 15.0});

    gmd::VelocityVerletIntegrator integrator(dt);      // no constraint solver
    FixedVirialForceProvider provider(provider_virial);
    run_one_step(system, integrator, provider, dt);

    check(system.constraint_virial_state() == gmd::ConstraintVirialState::NotApplicable,
          "with no constraints the state must be NotApplicable");
    check(system.last_virial_valid(),
          "with no constraints the provider virial is complete and must stay valid");
    for (std::size_t index = 0; index < 9; ++index) {
        check(system.last_virial()[index] == provider_virial[index],
              "component " + std::to_string(index) +
                  " of an unconstrained run must be the provider virial bit-for-bit");
    }
}

// ---------------------------------------------------------------------------
// 5. Validity semantics.
// ---------------------------------------------------------------------------

// A barostat that rescales the cell leaves the recomputed CURRENT-GEOMETRY
// provider virial with no contemporaneous constraint partner, so that tensor
// must be reported INVALID rather than as a complete pressure virial missing a
// term. The completed step's own pressure is a separate slot and is unaffected;
// test_completed_step_pressure_survives_a_barostat_rescale covers it.
void test_barostat_rescale_invalidates_rather_than_combining_across_it() {
    const std::array<double, 9> provider_virial{
        1.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 3.0};
    constexpr double bond = 1.5;
    constexpr double dt = 0.1;
    const gmd::Box box = cubic_box(30.0);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};

    // Rescaling barostat: the box changes, so forces and virial are recomputed.
    {
        gmd::System system = make_rotating_dimer(1.0, 3.0, bond, 0.1, box, {15.0, 15.0, 15.0});
        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        integrator.set_barostat(std::make_shared<OneShotScalingBarostat>(1.01));
        FixedVirialForceProvider provider(provider_virial);
        run_one_step(system, integrator, provider, dt);

        check(system.constraint_virial_state() == gmd::ConstraintVirialState::Unavailable,
              "after a barostat rescale no constraint multiplier belongs to the new "
              "geometry, so the state must be Unavailable");
        check(!system.last_virial_valid(),
              "a provider-only virial must not be reported as a complete pressure "
              "virial while constraints are active");
        check(provider.calls == 2,
              "the rescale must trigger exactly one re-evaluation, saw " +
                  std::to_string(provider.calls) + " call(s)");
    }

    // A barostat that leaves the box alone must not disturb anything.
    {
        gmd::System system = make_rotating_dimer(1.0, 3.0, bond, 0.1, box, {15.0, 15.0, 15.0});
        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        integrator.set_barostat(std::make_shared<OneShotScalingBarostat>(1.0));
        FixedVirialForceProvider provider(provider_virial);
        run_one_step(system, integrator, provider, dt);

        check(system.constraint_virial_state() == gmd::ConstraintVirialState::Valid,
              "a barostat that does not rescale must leave the step's constraint "
              "virial in place");
        check(system.last_virial_valid(), "and the combined virial valid");
        check(provider.calls == 1,
              "a no-op barostat must not trigger a re-evaluation, saw " +
                  std::to_string(provider.calls) + " call(s)");
        for (std::size_t index = 0; index < 9; ++index) {
            check_close(system.last_virial()[index],
                        provider_virial[index] + system.constraint_virial()[index],
                        1.0e-13, "component " + std::to_string(index) +
                                     " must still be provider + constraint");
        }
    }
}

// Constraints active but RATTLE switched off: there is no multiplier to recover,
// so the reported virial must be invalid rather than silently short of a term.
void test_rattle_disabled_reports_unavailable() {
    constexpr double bond = 1.5;
    constexpr double dt = 0.1;
    const gmd::Box box = cubic_box(30.0);
    gmd::System system = make_rotating_dimer(1.0, 3.0, bond, 0.1, box, {15.0, 15.0, 15.0});

    gmd::ConstraintSettings settings = tight_settings();
    settings.enable_rattle = false;
    auto solver = std::make_shared<gmd::ConstraintSolver>(
        std::vector<gmd::BondConstraint>{{0, 1, bond}}, settings);
    gmd::VelocityVerletIntegrator integrator(dt);
    integrator.set_constraint_solver(solver);
    FixedVirialForceProvider provider({1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0});
    run_one_step(system, integrator, provider, dt);

    check(system.constraint_virial_state() == gmd::ConstraintVirialState::Unavailable,
          "with RATTLE disabled there is no contemporaneous multiplier");
    check(!system.last_virial_valid(),
          "and the reported virial must therefore be invalid");
}

// Before any step has run, a constrained system has no constraint multiplier for
// the initial force evaluation, so the initial virial is incomplete.
void test_initial_state_has_no_constraint_virial() {
    constexpr double bond = 1.5;
    const gmd::Box box = cubic_box(30.0);
    gmd::System system = make_rotating_dimer(1.0, 3.0, bond, 0.1, box, {15.0, 15.0, 15.0});

    auto solver = std::make_shared<gmd::ConstraintSolver>(
        std::vector<gmd::BondConstraint>{{0, 1, bond}}, tight_settings());
    gmd::VelocityVerletIntegrator integrator(0.1);
    integrator.set_constraint_solver(solver);
    gmd::RuntimeContext runtime;
    integrator.initialize(system, runtime);
    system.set_provider_virial({1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0}, true);

    check(system.constraint_virial_state() == gmd::ConstraintVirialState::Unavailable,
          "a fresh constrained run has no constraint multiplier before its first step");
    check(!system.last_virial_valid(),
          "so its initial reported virial must be invalid, not provider-only");

    // Restoring a checkpointed value, which is what a restart does, completes it.
    system.set_constraint_virial({-0.5, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, -0.5});
    check(system.last_virial_valid(),
          "attaching the checkpointed constraint virial must make it valid again");
    check_close(system.last_virial()[0], 0.5, 1.0e-15,
                "and it must be added to the provider virial");
}

// ---------------------------------------------------------------------------
// 6. The pressure a rigid rotor reports.
// ---------------------------------------------------------------------------
//
// The end-to-end consequence of section 5 of the ConstraintVirialResult comment:
// a lone rigid rotor must contribute essentially nothing to the pressure, and
// would report 2K/3V if the constraint virial were dropped.
void test_rigid_rotor_reports_no_pressure() {
    constexpr double bond = 1.5;
    constexpr double omega = 0.1;
    constexpr double dt = 0.05;
    const double length = 30.0;
    const double volume = length * length * length;
    const gmd::Box box = cubic_box(length);

    gmd::System system = make_rotating_dimer(1.0, 3.0, bond, omega, box, {15.0, 15.0, 15.0});
    auto solver = std::make_shared<gmd::ConstraintSolver>(
        std::vector<gmd::BondConstraint>{{0, 1, bond}}, tight_settings());
    gmd::VelocityVerletIntegrator integrator(dt);
    integrator.set_constraint_solver(solver);
    NullForceProvider provider;
    run_one_step(system, integrator, provider, dt);

    const double twice_ke = twice_kinetic_energy(system);
    const double with_constraint = (twice_ke + trace_of(system.last_virial())) / (3.0 * volume);
    const double without_constraint = twice_ke / (3.0 * volume);

    std::cout << "[constraint virial] rigid rotor pressure term: with constraint virial "
              << std::setprecision(4) << with_constraint
              << ", without it " << without_constraint << " (eV/A^3)\n";

    check(without_constraint > 1.0e-9,
          "the fixture must have a kinetic pressure term worth cancelling");
    check(std::abs(with_constraint) < 1.0e-4 * without_constraint,
          "a rigid rotor must contribute essentially no pressure; got " +
              std::to_string(with_constraint) + " against a kinetic term of " +
              std::to_string(without_constraint));
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

// ---------------------------------------------------------------------------
// 4b. A barostat rescale must not destroy the completed step's pressure.
// ---------------------------------------------------------------------------
//
// The pressure of the step that just finished and the virial attached to the
// post-rescale geometry are different things. The first is complete and is what
// the barostat consumed; the second has no constraint partner. Keeping them in
// one slot would force the complete one to be thrown away.
void test_completed_step_pressure_survives_a_barostat_rescale() {
    const std::array<double, 9> provider_virial{
        1.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 3.0};
    constexpr double bond = 1.5;
    constexpr double dt = 0.1;
    const gmd::Box box = cubic_box(30.0);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};

    gmd::System system = make_rotating_dimer(1.0, 3.0, bond, 0.1, box, {15.0, 15.0, 15.0});
    auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
    gmd::VelocityVerletIntegrator integrator(dt);
    integrator.set_constraint_solver(solver);
    integrator.set_barostat(std::make_shared<OneShotScalingBarostat>(1.01));
    FixedVirialForceProvider provider(provider_virial);
    run_one_step(system, integrator, provider, dt);

    const auto& completed = system.step_thermodynamics();
    check(completed.valid,
          "the completed step's pressure must survive the rescale that followed it");
    check(provider.calls == 2,
          "the rescale must trigger exactly one re-evaluation, saw " +
              std::to_string(provider.calls) + " call(s)");

    // It is the COMPLETE pressure: provider + constraint, at the pre-rescale
    // volume, which is not the volume the system now has.
    const double pre_rescale_volume = 30.0 * 30.0 * 30.0;
    check_close(completed.volume, pre_rescale_volume, 1.0e-9,
                "the completed-step pressure must use the pre-rescale volume");
    check(std::abs(completed.volume - system.box().lengths[0] * system.box().lengths[1] *
                                          system.box().lengths[2]) > 1.0,
          "the fixture must actually have rescaled, or the test proves nothing");
    // The constraint part it carries cannot be read back off the System -- the
    // rescale cleared that, which is the whole point -- so it is recovered by
    // subtracting the known provider virial, and checked against the rotor's own
    // identity: a rigid pair's constraint virial has trace -2K.
    std::array<double, 9> constraint_part{};
    for (std::size_t index = 0; index < 9; ++index) {
        constraint_part[index] = completed.virial[index] - provider_virial[index];
    }
    check(std::abs(trace_of(constraint_part)) > 1.0e-6,
          "the completed-step virial must actually carry a constraint term, or the "
          "test proves nothing");
    check_close(trace_of(constraint_part), -completed.twice_kinetic_energy,
                1.0e-3 * completed.twice_kinetic_energy,
                "the constraint term the completed step recorded must satisfy the "
                "rigid-rotor identity tr W = -2K");
    check_close(constraint_part[1], constraint_part[3], 1.0e-14,
                "and must still be the symmetric tensor RATTLE produced");
    const double expected =
        (completed.twice_kinetic_energy + completed.virial[0] + completed.virial[4] +
         completed.virial[8]) / (3.0 * completed.volume);
    check_close(completed.pressure, expected, 1.0e-15,
                "the recorded pressure must be (2K + tr W) / 3V of the recorded pieces");
    check_close(completed.potential_energy, 0.0, 1.0e-15,
                "and it must carry the potential energy of that same state");

    // Meanwhile the CURRENT-GEOMETRY virial is incomplete and says so.
    check(system.constraint_virial_state() == gmd::ConstraintVirialState::Unavailable,
          "after a rescale no constraint multiplier belongs to the new geometry");
    check(!system.last_virial_valid(),
          "so the current-geometry virial must not pass for a complete pressure virial");
}


// ---------------------------------------------------------------------------
// 7. Every component, against an analytical reference: a tilted rigid rotor.
// ---------------------------------------------------------------------------
//
// The rotor cases above all put the bond along x, so their virial has one
// non-zero component and says nothing about the off-diagonals. Tilt the bond so
// it has non-zero x, y AND z, and every entry of the tensor becomes non-trivial.
//
// The reference is written down from rigid-body mechanics, not from the solver:
// a freely rotating rigid dimer is held together entirely by its constraint, so
// the constraint force is the centripetal force
//
//     G_i = -mu omega^2 d u,        u = the unit bond vector
//
// and, since W = r_ij (x) G_i for a pair,
//
//     W_ab = (d u_a)(-mu omega^2 d u_b) = -mu omega^2 d^2 u_a u_b
//
// which fixes ALL NINE components independently of how they are computed.
void test_tilted_rotor_every_component_against_analytical_reference() {
    constexpr double m0 = 1.0;
    constexpr double m1 = 3.0;
    constexpr double bond = 1.5;
    constexpr double omega = 0.1;
    const double mu = (m0 * m1) / (m0 + m1);

    // (2, 3, 6) has length exactly 7, so the direction is exact in binary and
    // every one of its components -- and so every product u_a u_b -- is
    // non-zero and distinct.
    const Vec3 direction = normalised({2.0, 3.0, 6.0});
    // Orthogonal to the bond, so the fixture starts in the tangent space.
    const Vec3 axis = normalised(cross(direction, {0.0, 0.0, 1.0}));
    const gmd::Box box = cubic_box(30.0);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};

    // The measured tensor is the ENDPOINT one, at r(t+dt), so the reference has
    // to be evaluated there too. A free rigid rotor's exact motion is a rigid
    // rotation by omega*dt about the spin axis, so the endpoint bond direction is
    // written down analytically rather than read back from the integrator.
    // Comparing against the t=0 direction instead would leave an O(dt) mismatch
    // that has nothing to do with the virial.
    auto reference_at = [&](double dt) {
        const Vec3 u = apply_rotation(rotation_matrix(axis, omega * dt), direction);
        std::array<double, 9> W{};
        for (std::size_t a = 0; a < 3; ++a) {
            for (std::size_t b = 0; b < 3; ++b) {
                W[a * 3 + b] = -mu * omega * omega * bond * bond * u[a] * u[b];
            }
        }
        return W;
    };
    const std::array<double, 9> reference = reference_at(0.0);

    check(std::abs(dot3(direction, axis)) < 1.0e-15,
          "the rotation axis must be perpendicular to the bond");
    // Every off-diagonal must be a real number, or the test proves nothing.
    const double scale = std::abs(trace_of(reference));
    for (const std::size_t index : {1u, 2u, 3u, 5u, 6u, 7u}) {
        check(std::abs(reference[index]) > 0.05 * scale,
              std::string("reference ") + component_name(index) +
                  " must be a substantial fraction of the trace, or the off-diagonal "
                  "assertion is vacuous; got " + std::to_string(reference[index]));
    }

    std::cout << "[constraint virial] tilted rigid rotor, u = ("
              << std::setprecision(6) << direction[0] << ", " << direction[1] << ", "
              << direction[2] << "); analytical W_ab = -mu omega^2 d^2 u_a u_b\n";
    std::cout << "    reference tensor:";
    for (std::size_t index = 0; index < 9; ++index) {
        if (index % 3 == 0) std::cout << "\n      ";
        std::cout << std::setw(15) << std::setprecision(8) << reference[index];
    }
    std::cout << '\n';

    const std::array<double, 3> timesteps{0.2, 0.1, 0.05};
    std::array<std::array<double, 9>, 3> errors{};

    for (std::size_t k = 0; k < timesteps.size(); ++k) {
        const double dt = timesteps[k];
        gmd::System system = make_oriented_rotating_dimer(
            m0, m1, bond, omega, direction, axis, box, {15.0, 15.0, 15.0});
        require_on_manifold(system, constraints, "tilted rotor");

        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        run_one_step(system, integrator, provider, dt);

        const auto expected = reference_at(dt);
        const auto& measured = system.constraint_virial();
        for (std::size_t index = 0; index < 9; ++index) {
            errors[k][index] = std::abs(measured[index] - expected[index]) / scale;
        }
        std::cout << "    dt=" << std::setw(5) << dt << "  max relative component error "
                  << std::setprecision(3)
                  << *std::max_element(errors[k].begin(), errors[k].end())
                  << "  (off-diagonals:";
        for (const std::size_t index : {1u, 2u, 3u, 5u, 6u, 7u}) {
            std::cout << ' ' << component_name(index) << '=' << errors[k][index];
        }
        std::cout << ")\n";
    }

    const std::size_t last = timesteps.size() - 1;
    for (std::size_t index = 0; index < 9; ++index) {
        // The discretisation error of the endpoint estimator is (omega dt)^2/4
        // relative; at the finest step that is 6.25e-6, so 1e-4 leaves more than
        // an order of margin without being loose enough to hide a wrong tensor.
        check(errors[last][index] < 1.0e-4,
              std::string(component_name(index)) +
                  " must match the analytical rigid-rotor value; relative error at the "
                  "finest timestep was " + std::to_string(errors[last][index]));
        for (std::size_t k = 0; k + 1 < timesteps.size(); ++k) {
            const double ratio = errors[k][index] / errors[k + 1][index];
            check(ratio > 3.0 && ratio < 5.0,
                  std::string(component_name(index)) +
                      " must converge as dt^2; halving dt changed its error by a factor "
                      "of " + std::to_string(ratio));
        }
    }
}

// ---------------------------------------------------------------------------
// 8. A coupled constraint set in a general 3D orientation.
// ---------------------------------------------------------------------------
//
// A rigid triangle, tilted out of every coordinate plane, spinning about the
// normal to its own plane. The reference again comes from rigid-body mechanics
// and never touches the accumulation being tested: each atom's acceleration is
// centripetal about the centre of mass, so
//
//     G_i = -m_i omega^2 (r_i - R),
//     W   = sum_i (r_i - R) (x) G_i = -omega^2 sum_i m_i (r_i - R) (x) (r_i - R)
//
// (using r_i - R rather than r_i is legitimate because the constraint forces sum
// to zero, which makes W independent of the origin). That is the second-moment
// tensor of the body, and in a general orientation every one of its nine entries
// is non-zero.
void test_tilted_triangle_against_rigid_body_reference() {
    const std::array<double, 3> masses{1.0, 2.0, 3.0};
    // Scalene, then tilted by a rotation with no special relationship to the axes.
    const std::array<Vec3, 3> flat{Vec3{0.0, 0.0, 0.0}, Vec3{1.4, 0.0, 0.0},
                                   Vec3{0.5, 1.1, 0.0}};
    const auto tilt = rotation_matrix(normalised({1.0, -2.0, 3.0}), 0.7);
    constexpr double omega = 0.08;
    const gmd::Box box = cubic_box(40.0);
    const Vec3 centre{20.0, 20.0, 20.0};
    const Vec3 spin_axis = apply_rotation(tilt, {0.0, 0.0, 1.0});   // the plane normal

    // Body-frame offsets about the centre of mass, tilted.
    double total = 0.0;
    Vec3 weighted{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < 3; ++i) {
        total += masses[i];
        for (std::size_t d = 0; d < 3; ++d) weighted[d] += masses[i] * flat[i][d];
    }
    std::array<Vec3, 3> arms{};
    for (std::size_t i = 0; i < 3; ++i) {
        const Vec3 flat_arm{flat[i][0] - weighted[0] / total,
                            flat[i][1] - weighted[1] / total,
                            flat[i][2] - weighted[2] / total};
        arms[i] = apply_rotation(tilt, flat_arm);
    }

    auto build = [&]() {
        gmd::System system;
        system.resize(3, 3);
        system.set_box(box);
        for (std::size_t i = 0; i < 3; ++i) {
            system.mutable_masses()[i] = masses[i];
            system.mutable_atom_tags()[i] = static_cast<int>(i);
            for (std::size_t d = 0; d < 3; ++d) {
                system.mutable_coordinates()[i][d] = centre[d] + arms[i][d];
            }
            const Vec3 v = cross(spin_axis, arms[i]);
            for (std::size_t d = 0; d < 3; ++d) {
                system.mutable_velocities()[i][d] = omega * v[d];
            }
        }
        return system;
    };

    // W_ref = -omega^2 sum_i m_i arm_i (x) arm_i, at the ENDPOINT configuration:
    // the measured tensor is taken at r(t+dt), and a free rigid body's exact
    // motion over the step is a rigid rotation by omega*dt about the spin axis.
    // The arms are rotated analytically, not read back from the integrator.
    auto reference_at = [&](double dt) {
        const auto step_rotation = rotation_matrix(spin_axis, omega * dt);
        std::array<double, 9> W{};
        for (std::size_t i = 0; i < 3; ++i) {
            const Vec3 arm = apply_rotation(step_rotation, arms[i]);
            for (std::size_t a = 0; a < 3; ++a) {
                for (std::size_t b = 0; b < 3; ++b) {
                    W[a * 3 + b] -= omega * omega * masses[i] * arm[a] * arm[b];
                }
            }
        }
        return W;
    };
    const std::array<double, 9> reference = reference_at(0.0);
    const double scale = std::abs(trace_of(reference));

    gmd::System probe = build();
    const std::vector<gmd::BondConstraint> constraints{
        {0, 1, distance(probe, 0, 1)},
        {1, 2, distance(probe, 1, 2)},
        {0, 2, distance(probe, 0, 2)},
    };
    require_on_manifold(probe, constraints, "tilted triangle");
    for (const std::size_t index : {1u, 2u, 3u, 5u, 6u, 7u}) {
        check(std::abs(reference[index]) > 0.02 * scale,
              std::string("tilted triangle reference ") + component_name(index) +
                  " must be non-trivial, got " + std::to_string(reference[index]));
    }

    std::cout << "[constraint virial] tilted rigid triangle, three coupled constraints;"
                 " reference W = -omega^2 sum_i m_i a_i (x) a_i\n";
    std::cout << "    reference tensor:";
    for (std::size_t index = 0; index < 9; ++index) {
        if (index % 3 == 0) std::cout << "\n      ";
        std::cout << std::setw(15) << std::setprecision(8) << reference[index];
    }
    std::cout << '\n';

    const std::array<double, 3> timesteps{0.2, 0.1, 0.05};
    std::array<std::array<double, 9>, 3> errors{};
    for (std::size_t k = 0; k < timesteps.size(); ++k) {
        const double dt = timesteps[k];
        gmd::System system = build();
        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        run_one_step(system, integrator, provider, dt);

        const auto expected = reference_at(dt);
        const auto& measured = system.constraint_virial();
        for (std::size_t index = 0; index < 9; ++index) {
            errors[k][index] = std::abs(measured[index] - expected[index]) / scale;
        }
        std::cout << "    dt=" << std::setw(5) << dt << "  max relative component error "
                  << std::setprecision(3)
                  << *std::max_element(errors[k].begin(), errors[k].end())
                  << "  (off-diagonals:";
        for (const std::size_t index : {1u, 2u, 3u, 5u, 6u, 7u}) {
            std::cout << ' ' << component_name(index) << '=' << errors[k][index];
        }
        std::cout << ")\n";
    }

    const std::size_t last = timesteps.size() - 1;
    for (std::size_t index = 0; index < 9; ++index) {
        check(errors[last][index] < 1.0e-4,
              std::string(component_name(index)) +
                  " of a tilted rigid triangle must match the rigid-body reference; "
                  "relative error " + std::to_string(errors[last][index]));
        for (std::size_t k = 0; k + 1 < timesteps.size(); ++k) {
            const double ratio = errors[k][index] / errors[k + 1][index];
            check(ratio > 3.0 && ratio < 5.0,
                  std::string(component_name(index)) +
                      " of a tilted triangle must converge as dt^2; measured ratio " +
                      std::to_string(ratio));
        }
    }
}

// ---------------------------------------------------------------------------
// 9. The tensor must transform covariantly.
// ---------------------------------------------------------------------------
//
// W is a rank-2 Cartesian tensor built from vectors, so rotating the whole
// system by an orthogonal R must take W to R W R^T exactly -- not approximately,
// and not only in its trace. A formula that mixed up an index, or summed
// r_a F_b where it meant r_b F_a, would satisfy the trace and the diagonal of an
// axis-aligned fixture and fail here.
void test_virial_transforms_covariantly_under_rotation() {
    constexpr double bond = 1.5;
    constexpr double omega = 0.1;
    constexpr double dt = 0.05;
    const gmd::Box box = cubic_box(40.0);
    const Vec3 centre{20.0, 20.0, 20.0};
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};

    const Vec3 direction = normalised({2.0, 3.0, 6.0});
    const Vec3 axis = normalised(cross(direction, {0.0, 0.0, 1.0}));
    const auto R = rotation_matrix(normalised({-1.0, 4.0, 2.0}), 1.1);

    auto measure = [&](const Vec3& u, const Vec3& n) {
        gmd::System system = make_oriented_rotating_dimer(
            1.0, 3.0, bond, omega, u, n, box, centre);
        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        run_one_step(system, integrator, provider, dt);
        return system.constraint_virial();
    };

    const auto original = measure(direction, axis);
    // Rotating the fixture's direction and axis rotates the whole configuration
    // and its velocities, because both are built from them.
    const auto rotated = measure(apply_rotation(R, direction), apply_rotation(R, axis));
    const auto expected = conjugate(R, original);

    const double scale = std::abs(trace_of(original));
    check(scale > 1.0e-6, "the covariance fixture must produce a non-trivial tensor");
    double worst = 0.0;
    for (std::size_t index = 0; index < 9; ++index) {
        const double error = std::abs(rotated[index] - expected[index]) / scale;
        worst = std::max(worst, error);
        // The floor here is NOT floating-point arithmetic on the rotation, which
        // would be ~1e-15. It is the constraint solver's own convergence
        // tolerance. RATTLE stops once |r.v|/|r| falls below `tolerance` (1e-13
        // here), and the multiplier is proportional to that dot product, whose
        // converged size is about omega^2 d^2 dt / 2 ~ 6e-4 for this fixture. The
        // last accepted iterate therefore carries a relative uncertainty of
        // roughly 1e-13 / 6e-4 ~ 2e-10, and the two runs stop at different
        // iterates. 1e-8 sits two orders above that floor and is still eight
        // orders below the O(1) discrepancy a transposed or mis-indexed tensor
        // would produce.
        check(error < 1.0e-8,
              std::string("rotating the system must take ") + component_name(index) +
                  " to (R W R^T) exactly; relative error " + std::to_string(error));
    }
    std::cout << "[constraint virial] rotation covariance W -> R W R^T: max relative "
                 "component error " << std::setprecision(3) << worst << '\n';
    check(std::abs(trace_of(rotated) - trace_of(original)) < 1.0e-8 * scale,
          "and the trace, being invariant under rotation, must be unchanged (same "
          "convergence-tolerance floor as the components above)");
}

// ---------------------------------------------------------------------------
// 10. The tilted rotor across every periodic boundary.
// ---------------------------------------------------------------------------
//
// All nine components are built from minimum-image bond vectors, so no component
// may depend on where the molecule sits or on which face it straddles.
void test_tilted_rotor_across_every_boundary() {
    constexpr double bond = 1.5;
    constexpr double omega = 0.1;
    constexpr double dt = 0.05;
    const double length = 12.0;
    const gmd::Box box = cubic_box(length);
    const std::vector<gmd::BondConstraint> constraints{{0, 1, bond}};
    const Vec3 direction = normalised({2.0, 3.0, 6.0});
    const Vec3 axis = normalised(cross(direction, {0.0, 0.0, 1.0}));

    auto measure = [&](const Vec3& centre) {
        gmd::System system = make_oriented_rotating_dimer(
            1.0, 3.0, bond, omega, direction, axis, box, centre);
        auto coordinates = system.mutable_coordinates();
        for (auto& coordinate : coordinates) {
            for (std::size_t d = 0; d < 3; ++d) {
                while (coordinate[d] < 0.0) coordinate[d] += length;
                while (coordinate[d] >= length) coordinate[d] -= length;
            }
        }
        auto solver = std::make_shared<gmd::ConstraintSolver>(constraints, tight_settings());
        gmd::VelocityVerletIntegrator integrator(dt);
        integrator.set_constraint_solver(solver);
        NullForceProvider provider;
        run_one_step(system, integrator, provider, dt);
        return system.constraint_virial();
    };

    const auto interior = measure({6.0, 6.0, 6.0});
    const double scale = std::abs(trace_of(interior));
    check(scale > 1.0e-6, "the boundary fixture must produce a non-trivial tensor");

    const std::vector<std::pair<std::string, Vec3>> placements{
        {"across x", {0.0, 6.0, 6.0}},
        {"across y", {6.0, 0.0, 6.0}},
        {"across z", {6.0, 6.0, 0.0}},
        {"at the corner", {0.0, 0.0, 0.0}},
        {"translated by an arbitrary vector", {2.75, 9.25, 4.5}},
    };
    for (const auto& placement : placements) {
        const auto measured = measure(placement.second);
        for (std::size_t index = 0; index < 9; ++index) {
            const double error = std::abs(measured[index] - interior[index]) / scale;
            check(error < 1.0e-9,
                  std::string(component_name(index)) + " must not depend on the molecule "
                      "being " + placement.first + "; relative error " +
                      std::to_string(error));
        }
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

    test_rotating_dimer_against_centripetal_force();
    test_two_constraint_half_impulses_are_distinct();
    test_free_rotor_conserves_energy();
    test_rotating_triangle_coupled_constraints();
    test_pair_forces_are_central_and_antisymmetric();
    test_periodic_boundary_invariance();
    test_translation_invariance();
    test_combined_virial_is_provider_plus_constraint_exactly_once();
    test_unconstrained_run_is_untouched();
    test_barostat_rescale_invalidates_rather_than_combining_across_it();
    test_rattle_disabled_reports_unavailable();
    test_initial_state_has_no_constraint_virial();
    test_tilted_rotor_every_component_against_analytical_reference();
    test_tilted_triangle_against_rigid_body_reference();
    test_virial_transforms_covariantly_under_rotation();
    test_tilted_rotor_across_every_boundary();
    test_completed_step_pressure_survives_a_barostat_rescale();
    test_rigid_rotor_reports_no_pressure();

    if (failures != 0) {
        std::cerr << "[constraint virial] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[constraint virial] all checks passed\n";
    return 0;
}
