// Regression tests for the state left behind after a barostat volume change.
//
// VelocityVerletIntegrator::step() evaluates forces, stores them, and only then
// applies the barostat. A barostat that rescales the box and the coordinates
// therefore used to leave the system holding forces computed for the *previous*
// geometry, which the next step's first half-kick would then integrate.
//
// step() now detects an accepted volume change and re-establishes forces, virial
// and neighbor-list state for the final geometry before returning. Conversely, a
// step that did not change the cell -- no barostat, a rejected Monte Carlo trial,
// or a Monte Carlo step that was not a move attempt -- must not pay for a second
// force evaluation.

#include <array>
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
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/berendsen_barostat.hpp"
#include "gmd/integrator/mc_barostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[barostat refresh] " << message << '\n';
        ++failures;
    }
}

constexpr double kBoxLength = 14.0;
constexpr double kCutoff = 4.0;

// Wraps a real LJ provider and counts how many times compute() is called, so the
// tests can assert that no redundant force evaluation is performed.
class CountingForceProvider final : public gmd::ForceProvider {
public:
    CountingForceProvider()
        : inner_(0.25, 2.5, kCutoff) {}

    std::string_view name() const noexcept override { return "counting_lj"; }
    void initialize(gmd::RuntimeContext& runtime) override { inner_.initialize(runtime); }
    void finalize(gmd::RuntimeContext& runtime) override { inner_.finalize(runtime); }

    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext& runtime) override {
        ++compute_count;
        inner_.compute(request, result, runtime);
    }

    int compute_count = 0;

private:
    gmd::ClassicalForceProvider inner_;
};

gmd::System make_system() {
    // A slightly irregular cluster: a perfect lattice can produce forces that
    // cancel to zero, which would make the comparison below vacuous.
    gmd::System system;
    const std::size_t n = 8;
    system.resize(n, n);
    gmd::Box box;
    box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    system.set_box(box);

    const double offsets[8][3] = {
        {3.0, 3.0, 3.0}, {6.1, 3.2, 3.0}, {3.0, 6.3, 3.1}, {6.2, 6.0, 3.0},
        {3.1, 3.0, 6.2}, {6.0, 3.1, 6.0}, {3.0, 6.1, 6.3}, {6.3, 6.2, 6.1},
    };
    for (std::size_t i = 0; i < n; ++i) {
        system.mutable_masses()[i] = 12.0;
        system.mutable_coordinates()[i] = {offsets[i][0], offsets[i][1], offsets[i][2]};
        system.mutable_velocities()[i] = {0.001 * static_cast<double>(i + 1), -0.0005, 0.0007};
    }
    return system;
}

// Evaluates forces from scratch for whatever geometry `system` currently holds.
gmd::ForceResult fresh_force_evaluation(gmd::System& system,
                                        gmd::ForceProvider& provider,
                                        gmd::RuntimeContext& runtime) {
    gmd::ForceResult result;
    const auto coordinates = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .step = 1,
        .time = 1.0,
        .coordinates = std::span<const gmd::Coordinate3D>(coordinates.data(),
                                                          coordinates.size()),
        .neighbor_list = nullptr,
    };
    provider.compute(request, result, runtime);
    return result;
}

double max_force_difference(const gmd::System& system, const gmd::ForceResult& fresh) {
    double worst = 0.0;
    const auto stored = system.forces();
    for (std::size_t i = 0; i < system.num_local_atoms(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            worst = std::max(worst, std::abs(stored[i][d] - fresh.forces[i][d]));
        }
    }
    return worst;
}

// --- Berendsen: the box always changes, so forces must always be refreshed ---

void test_berendsen_refreshes_forces() {
    gmd::System system = make_system();
    CountingForceProvider provider;
    gmd::RuntimeContext runtime;

    // Coupling chosen so that mu^3 = 1 - beta*(dt/tau)*(P_target - P) lands a
    // little below 1: enough to rescale the cell measurably, but comfortably
    // clear of the mu^3 <= 0 guard that would make the barostat a no-op.
    //   1 - 4.5e-5 * (0.5/10) * 4000 ~= 0.991
    auto barostat = std::make_shared<gmd::BerendsenBarostat>(10.0, 4.5e-5);

    gmd::VelocityVerletIntegrator integrator(0.5);
    integrator.set_barostat(barostat);
    integrator.set_target_pressure(4000.0);
    integrator.set_target_temperature(300.0);
    integrator.initialize(system, runtime);

    const gmd::Box box_before = system.box();

    const gmd::IntegratorStepContext ctx{.step = 0, .dt = 0.5};
    integrator.step(system, provider, ctx, runtime);

    const gmd::Box box_after = system.box();
    check(!gmd::box_lengths_equal(box_before, box_after),
          "Berendsen barostat should have rescaled the box; the test needs it to");

    // The heart of the matter: the forces the system now holds must be the
    // forces of the geometry it now holds.
    const int count_before_check = provider.compute_count;
    gmd::ForceResult fresh = fresh_force_evaluation(system, provider, runtime);
    provider.compute_count = count_before_check;  // don't count the probe

    const double worst = max_force_difference(system, fresh);
    check(worst < 1.0e-12,
          "stored forces do not match a fresh evaluation at the post-barostat "
          "coordinates and box (max component difference " + std::to_string(worst) + ")");

    // The cached virial trace used by the next barostat call must be refreshed
    // from the same evaluation, not left over from the pre-scaling geometry.
    const double fresh_trace = fresh.virial[0] + fresh.virial[4] + fresh.virial[8];
    const double stored_trace = system.last_virial()[0] + system.last_virial()[4] +
                                system.last_virial()[8];
    check(std::abs(stored_trace - fresh_trace) < 1.0e-12,
          "stored virial does not match a fresh evaluation at the post-barostat "
          "geometry");

    // One evaluation for the step itself, one to refresh after the volume
    // change. Anything more means a redundant evaluation crept in.
    check(provider.compute_count == 2,
          "expected exactly 2 force evaluations for a volume-changing step, got " +
              std::to_string(provider.compute_count));
}

// --- No barostat: nothing to refresh, and no extra evaluation --------------

void test_no_barostat_is_untouched() {
    gmd::System system = make_system();
    CountingForceProvider provider;
    gmd::RuntimeContext runtime;

    gmd::VelocityVerletIntegrator integrator(0.5);
    integrator.set_target_temperature(300.0);
    integrator.initialize(system, runtime);

    const gmd::Box box_before = system.box();
    const gmd::IntegratorStepContext ctx{.step = 0, .dt = 0.5};
    integrator.step(system, provider, ctx, runtime);

    check(gmd::box_lengths_equal(box_before, system.box()),
          "a step without a barostat must not change the box");
    check(provider.compute_count == 1,
          "a step without a barostat must evaluate forces exactly once, got " +
              std::to_string(provider.compute_count));

    const int count_before_check = provider.compute_count;
    gmd::ForceResult fresh = fresh_force_evaluation(system, provider, runtime);
    provider.compute_count = count_before_check;
    check(max_force_difference(system, fresh) < 1.0e-12,
          "stored forces must match the current geometry even without a barostat");
}

// --- Monte Carlo: rejected trials must cost nothing extra ------------------

void test_mc_rejected_trial_costs_nothing_extra() {
    gmd::System system = make_system();
    CountingForceProvider provider;
    gmd::RuntimeContext runtime;

    // frequency = 100 means no volume move is attempted on step 1, so this
    // stands in for every step on which the Monte Carlo barostat does nothing.
    auto barostat = std::make_shared<gmd::MCBarostat>(100, 0.01);

    gmd::VelocityVerletIntegrator integrator(0.5);
    integrator.set_barostat(barostat);
    integrator.set_target_pressure(1.0);
    integrator.set_target_temperature(300.0);
    integrator.initialize(system, runtime);

    const gmd::Box box_before = system.box();
    const gmd::IntegratorStepContext ctx{.step = 1, .dt = 0.5};
    integrator.step(system, provider, ctx, runtime);

    check(gmd::box_lengths_equal(box_before, system.box()),
          "no volume move should have been attempted on this step");
    check(provider.compute_count == 1,
          "a step with no accepted volume move must evaluate forces exactly once, "
          "got " + std::to_string(provider.compute_count));

    const int count_before_check = provider.compute_count;
    gmd::ForceResult fresh = fresh_force_evaluation(system, provider, runtime);
    provider.compute_count = count_before_check;
    check(max_force_difference(system, fresh) < 1.0e-12,
          "stored forces must match the current geometry after a no-op Monte "
          "Carlo step");
}

// --- Monte Carlo: an accepted move must refresh forces ---------------------

void test_mc_accepted_move_refreshes_forces() {
    gmd::System system = make_system();
    CountingForceProvider provider;
    gmd::RuntimeContext runtime;

    // frequency = 1 attempts a move on every step. Acceptance is stochastic, so
    // step until the box actually changes, then assert the invariant.
    auto barostat = std::make_shared<gmd::MCBarostat>(1, 0.02);

    gmd::VelocityVerletIntegrator integrator(0.5);
    integrator.set_barostat(barostat);
    integrator.set_target_pressure(1.0);
    integrator.set_target_temperature(300.0);
    integrator.initialize(system, runtime);

    bool observed_accepted_move = false;
    for (std::uint64_t step = 0; step < 40 && !observed_accepted_move; ++step) {
        const gmd::Box box_before = system.box();
        const gmd::IntegratorStepContext ctx{.step = step, .dt = 0.5};
        integrator.step(system, provider, ctx, runtime);

        if (gmd::box_lengths_equal(box_before, system.box())) {
            continue;  // trial rejected, nothing to check
        }
        observed_accepted_move = true;

        gmd::ForceResult fresh = fresh_force_evaluation(system, provider, runtime);
        const double worst = max_force_difference(system, fresh);
        check(worst < 1.0e-12,
              "after an accepted Monte Carlo volume move the stored forces do not "
              "match a fresh evaluation (max component difference " +
                  std::to_string(worst) + ")");
    }

    check(observed_accepted_move,
          "expected at least one accepted Monte Carlo volume move in 40 steps; "
          "the test did not exercise the accepted path");
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    // An MPI-enabled build still runs this as a single-process test, but the
    // force providers issue collectives that require an initialised MPI.
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif
    test_berendsen_refreshes_forces();
    test_no_barostat_is_untouched();
    test_mc_rejected_trial_costs_nothing_extra();
    test_mc_accepted_move_refreshes_forces();

    if (failures != 0) {
        std::cerr << "[barostat refresh] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[barostat refresh] all checks passed\n";
    return 0;
}
