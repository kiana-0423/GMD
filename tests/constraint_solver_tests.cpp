#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>
#include <span>
#include <string>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/system.hpp"

namespace {

class ZeroForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "zero_force"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.success = true;
        result.potential_energy = 0.0;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.virial_valid = true;
    }
};

class HarmonicBondForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "harmonic_bond"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.success = true;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.virial_valid = true;
        const auto& ri = request.coordinates[0];
        const auto& rj = request.coordinates[1];
        const std::array<double, 3> dr{
            ri[0] - rj[0],
            ri[1] - rj[1],
            ri[2] - rj[2],
        };
        const double r = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
        const double delta = r - 1.0;
        result.potential_energy = k_ * delta * delta;
        const double factor = r > 1.0e-12 ? -2.0 * k_ * delta / r : 0.0;
        for (std::size_t dim = 0; dim < 3; ++dim) {
            result.forces[0][dim] += factor * dr[dim];
            result.forces[1][dim] -= factor * dr[dim];
        }
    }

private:
    double k_ = 10.0;
};

void check(bool condition, const std::string& message, int& failures) {
    if (!condition) {
        std::cerr << "[constraints] " << message << "\n";
        ++failures;
    }
}

double distance(const gmd::System& system, std::size_t i, std::size_t j) {
    const auto coords = system.coordinates();
    const double dx = coords[i][0] - coords[j][0];
    const double dy = coords[i][1] - coords[j][1];
    const double dz = coords[i][2] - coords[j][2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double velocity_dot_bond(const gmd::System& system, std::size_t i, std::size_t j) {
    const auto coords = system.coordinates();
    const auto velocities = system.velocities();
    const std::array<double, 3> r{
        coords[i][0] - coords[j][0],
        coords[i][1] - coords[j][1],
        coords[i][2] - coords[j][2],
    };
    const std::array<double, 3> v{
        velocities[i][0] - velocities[j][0],
        velocities[i][1] - velocities[j][1],
        velocities[i][2] - velocities[j][2],
    };
    return r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
}

gmd::System make_system(std::size_t n) {
    gmd::System system;
    system.resize(n, n);
    gmd::Box box;
    box.set_lengths({30.0, 30.0, 30.0});
    system.set_box(box);
    for (auto& mass : system.mutable_masses()) {
        mass = 1.0;
    }
    return system;
}

void test_two_atom_projection(int& failures) {
    auto system = make_system(2);
    system.mutable_coordinates()[0] = {10.0, 10.0, 10.0};
    system.mutable_coordinates()[1] = {11.2, 10.0, 10.0};
    system.mutable_velocities()[0] = {0.2, 0.1, 0.0};
    system.mutable_velocities()[1] = {-0.1, -0.1, 0.0};

    gmd::ConstraintSolver solver({gmd::BondConstraint{0, 1, 1.0}},
                                 gmd::ConstraintSettings{1.0e-10, 100, true});
    const auto shake = solver.apply_shake(system);
    const auto rattle = solver.apply_rattle(system);
    check(shake.converged && std::abs(distance(system, 0, 1) - 1.0) < 1.0e-9,
          "two-atom SHAKE did not enforce target length", failures);
    check(rattle.converged && std::abs(velocity_dot_bond(system, 0, 1)) < 1.0e-9,
          "two-atom RATTLE did not remove radial velocity", failures);
}

void test_water_geometry(int& failures) {
    auto system = make_system(3);
    system.mutable_masses()[0] = 15.999;
    system.mutable_masses()[1] = 1.008;
    system.mutable_masses()[2] = 1.008;
    system.mutable_coordinates()[0] = {10.0, 10.0, 10.0};
    system.mutable_coordinates()[1] = {10.98, 10.02, 10.0};
    system.mutable_coordinates()[2] = {9.74, 10.91, 10.0};

    const double oh = 0.9572;
    const double angle = 104.52 * std::numbers::pi / 180.0;
    const double hh = std::sqrt(2.0 * oh * oh * (1.0 - std::cos(angle)));
    gmd::ConstraintSolver solver({
        gmd::BondConstraint{0, 1, oh},
        gmd::BondConstraint{0, 2, oh},
        gmd::BondConstraint{1, 2, hh},
    }, gmd::ConstraintSettings{1.0e-9, 200, true});
    const auto shake = solver.apply_shake(system);
    check(shake.converged, "water SHAKE did not converge", failures);
    check(std::abs(distance(system, 0, 1) - oh) < 5.0e-8,
          "water O-H(1) length error too large", failures);
    check(std::abs(distance(system, 0, 2) - oh) < 5.0e-8,
          "water O-H(2) length error too large", failures);
    check(std::abs(distance(system, 1, 2) - hh) < 5.0e-8,
          "water H-H length error too large", failures);
}

void test_nve_stability(int& failures) {
    auto system = make_system(2);
    system.mutable_coordinates()[0] = {10.0, 10.0, 10.0};
    system.mutable_coordinates()[1] = {11.0, 10.0, 10.0};
    system.mutable_velocities()[0] = {0.0, 0.1, 0.0};
    system.mutable_velocities()[1] = {0.0, -0.1, 0.0};

    auto solver = std::make_shared<gmd::ConstraintSolver>(
        std::vector<gmd::BondConstraint>{gmd::BondConstraint{0, 1, 1.0}},
        gmd::ConstraintSettings{1.0e-9, 100, true});
    gmd::VelocityVerletIntegrator integrator(0.05);
    integrator.set_constraint_solver(solver);
    ZeroForceProvider force;
    gmd::RuntimeContext runtime;
    integrator.initialize(system, runtime);
    force.initialize(runtime);

    const double initial_ke = 0.5 * gmd::compute_twice_ke(system);
    for (std::uint64_t step = 0; step < 1000; ++step) {
        integrator.step(system, force, gmd::IntegratorStepContext{step, 0.05}, runtime);
    }
    const double final_ke = 0.5 * gmd::compute_twice_ke(system);
    check(std::abs(distance(system, 0, 1) - 1.0) < 1.0e-7,
          "constrained NVE bond length drifted", failures);
    check(std::abs(final_ke - initial_ke) < 1.0e-2,
          "constrained NVE kinetic energy drift too large: initial=" +
              std::to_string(initial_ke) + " final=" + std::to_string(final_ke),
          failures);
}

void test_small_unconstrained_vs_large_constrained(int& failures) {
    auto unconstrained = make_system(2);
    unconstrained.mutable_coordinates()[0] = {10.0, 10.0, 10.0};
    unconstrained.mutable_coordinates()[1] = {11.1, 10.0, 10.0};

    HarmonicBondForceProvider force;
    gmd::RuntimeContext runtime;
    gmd::VelocityVerletIntegrator small_dt(0.001);
    small_dt.initialize(unconstrained, runtime);
    force.initialize(runtime);
    for (std::uint64_t step = 0; step < 1000; ++step) {
        small_dt.step(unconstrained, force, gmd::IntegratorStepContext{step, 0.001}, runtime);
    }
    check(std::isfinite(distance(unconstrained, 0, 1)) &&
              std::abs(distance(unconstrained, 0, 1) - 1.0) < 0.2,
          "unconstrained small timestep harmonic bond became unstable", failures);

    auto constrained = make_system(2);
    constrained.mutable_coordinates()[0] = {10.0, 10.0, 10.0};
    constrained.mutable_coordinates()[1] = {11.1, 10.0, 10.0};
    auto solver = std::make_shared<gmd::ConstraintSolver>(
        std::vector<gmd::BondConstraint>{gmd::BondConstraint{0, 1, 1.0}},
        gmd::ConstraintSettings{1.0e-9, 100, true});
    gmd::VelocityVerletIntegrator large_dt(0.05);
    large_dt.set_constraint_solver(solver);
    large_dt.initialize(constrained, runtime);
    for (std::uint64_t step = 0; step < 1000; ++step) {
        large_dt.step(constrained, force, gmd::IntegratorStepContext{step, 0.05}, runtime);
    }
    check(std::abs(distance(constrained, 0, 1) - 1.0) < 1.0e-7,
          "constrained large timestep did not hold the harmonic bond length", failures);
}

}  // namespace

int main() {
    int failures = 0;
    test_two_atom_projection(failures);
    test_water_geometry(failures);
    test_nve_stability(failures);
    test_small_unconstrained_vs_large_constrained(failures);
    return failures == 0 ? 0 : 1;
}
