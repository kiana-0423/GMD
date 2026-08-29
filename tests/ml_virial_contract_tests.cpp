// MLForceProvider's virial contract.
//
// No ML backend in this engine reports a stress tensor, so the provider has no
// virial to offer and must say so. It previously said nothing at all: it left
// `ForceResult::virial` and `virial_valid` exactly as it found them.
//
// Every call site in the engine happens to hand it a freshly constructed
// ForceResult, whose `virial_valid` defaults to false, so the omission was not
// reachable today. That is a property of the callers, not of the provider, and
// it is the kind of property that stops holding quietly: a caller that reuses
// one ForceResult across providers -- an obvious optimisation -- would have
// handed the model a previous provider's tensor and had the pressure computed
// from it. The provider now clears both fields explicitly.
//
// WHAT IS DELIBERATELY NOT DONE. `sum_i r_i (x) F_i` is not synthesised as a
// stand-in. For a periodic, cell-dependent model that expression is not the
// virial -- the same reason it is wrong for the Ewald reciprocal term, where
// the energy depends on the cell explicitly through 1/V and through
// k = 2 pi n / L. Manufacturing one would let a pressure-coupled barostat
// consume a number that does not mean what it claims. The correct behaviour is
// to report nothing and let the machinery downstream refuse to run.

#include <array>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/composite_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/force/ml_force_provider.hpp"
#include "gmd/force/model_runtime_adapter.hpp"
#include "gmd/integrator/berendsen_barostat.hpp"
#include "gmd/integrator/mc_barostat.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[ml virial contract] " << message << '\n';
        ++failures;
    }
}

// A stand-in ML runtime: an energy and forces, and no stress tensor. That is
// the shape of every backend this provider currently supports.
class StubModelAdapter final : public gmd::ModelRuntimeAdapter {
public:
    std::string_view name() const noexcept override { return "stub_model"; }
    void load_model(const std::filesystem::path&, gmd::RuntimeContext&) override {}
    void unload_model(gmd::RuntimeContext&) override {}
    void evaluate(const gmd::ModelEvaluationRequest& request,
                  gmd::ModelEvaluationResult& result,
                  gmd::RuntimeContext&) override {
        result.success = true;
        result.total_energy = -1.25;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.1, -0.2, 0.3});
    }
};

// A provider with a known, valid virial, so the composite behaviour below is
// about propagation rather than about the numbers.
class StubVirialProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "stub_virial"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.success = true;
        result.potential_energy = 2.5;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.virial = {1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0};
        result.virial_valid = true;
    }
};

gmd::System make_system() {
    gmd::System system;
    system.resize(2, 2);
    gmd::Box box;
    box.set_lengths({20.0, 20.0, 20.0});
    system.set_box(box);
    system.mutable_coordinates()[0] = {8.0, 8.0, 8.0};
    system.mutable_coordinates()[1] = {10.13, 9.47, 10.71};
    system.mutable_masses()[0] = 12.0;
    system.mutable_masses()[1] = 12.0;
    return system;
}

gmd::ForceRequest make_request(gmd::System& system,
                               std::span<const gmd::Coordinate3D> coordinates) {
    return gmd::ForceRequest{
        .system = &system,
        .box = &system.box(),
        .step = 0,
        .time = 0.0,
        .coordinates = coordinates,
        .neighbor_list = nullptr,
    };
}

double tensor_scale(const std::array<double, 9>& tensor) {
    double scale = 0.0;
    for (const double value : tensor) scale = std::max(scale, std::fabs(value));
    return scale;
}

// ---------------------------------------------------------------------------

void test_fresh_result_reports_no_virial() {
    gmd::System system = make_system();
    gmd::MLForceProvider provider("stub.pt", std::make_shared<StubModelAdapter>());
    gmd::RuntimeContext runtime;

    const auto coordinates = system.coordinates();
    gmd::ForceResult result;
    provider.compute(make_request(system, {coordinates.data(), coordinates.size()}),
                     result, runtime);

    check(result.success, "the stub model evaluation should succeed");
    check(result.potential_energy == -1.25,
          "the provider should pass the model's energy through");
    check(!result.virial_valid,
          "MLForceProvider must report virial_valid == false: it has no stress "
          "tensor to report");
    check(tensor_scale(result.virial) == 0.0,
          "MLForceProvider must leave the virial zeroed rather than carrying a "
          "value a caller might read");
}

// The regression that motivated the fix.
void test_stale_result_is_cleared() {
    gmd::System system = make_system();
    gmd::MLForceProvider provider("stub.pt", std::make_shared<StubModelAdapter>());
    gmd::RuntimeContext runtime;

    gmd::ForceResult reused;
    reused.virial = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    reused.virial_valid = true;
    reused.potential_energy = 123.0;

    const auto coordinates = system.coordinates();
    provider.compute(make_request(system, {coordinates.data(), coordinates.size()}),
                     reused, runtime);

    check(!reused.virial_valid,
          "a ForceResult that already carried virial_valid == true must come back "
          "false: the model did not produce that tensor and must not inherit it");
    check(tensor_scale(reused.virial) == 0.0,
          "a stale virial must be cleared, not left in place (largest remaining "
          "component " + std::to_string(tensor_scale(reused.virial)) + ")");
    check(reused.potential_energy == -1.25,
          "the stale energy must also be replaced, not accumulated onto");
}

// Nothing may quietly turn the absent tensor into a present one.
void test_absent_virial_is_not_synthesised_from_the_force_moment() {
    gmd::System system = make_system();
    gmd::MLForceProvider provider("stub.pt", std::make_shared<StubModelAdapter>());
    gmd::RuntimeContext runtime;

    const auto coordinates = system.coordinates();
    gmd::ForceResult result;
    provider.compute(make_request(system, {coordinates.data(), coordinates.size()}),
                     result, runtime);

    // The stub returns non-zero forces, so sum_i r_i (x) F_i is non-zero. If
    // the provider ever started reporting it, this fixture would notice.
    double moment = 0.0;
    for (std::size_t i = 0; i < result.forces.size(); ++i) {
        for (std::size_t a = 0; a < 3; ++a) {
            for (std::size_t b = 0; b < 3; ++b) {
                moment += std::fabs(coordinates[i][a] * result.forces[i][b]);
            }
        }
    }
    check(moment > 1.0,
          "the fixture's force moment is negligible, so it could not detect the "
          "provider synthesising one");
    check(tensor_scale(result.virial) == 0.0,
          "the reported virial is non-zero, which for a provider with no stress "
          "tensor can only mean a force moment was synthesised");
}

// A composite containing an ML child must be invalid as a whole, whatever else
// it contains and in whatever order.
void test_composite_with_an_ml_child_is_invalid() {
    for (const bool ml_first : {true, false}) {
        auto composite = std::make_shared<gmd::CompositeForceProvider>();
        if (ml_first) {
            composite->add(std::make_shared<gmd::MLForceProvider>(
                "stub.pt", std::make_shared<StubModelAdapter>()));
            composite->add(std::make_shared<StubVirialProvider>());
        } else {
            composite->add(std::make_shared<StubVirialProvider>());
            composite->add(std::make_shared<gmd::MLForceProvider>(
                "stub.pt", std::make_shared<StubModelAdapter>()));
        }

        gmd::System system = make_system();
        gmd::RuntimeContext runtime;
        const auto coordinates = system.coordinates();
        gmd::ForceResult result;
        composite->compute(make_request(system, {coordinates.data(), coordinates.size()}),
                           result, runtime);

        check(result.success,
              std::string("an absent virial is not a failed evaluation (ML ") +
                  (ml_first ? "first" : "second") + ")");
        check(!result.virial_valid,
              std::string("a composite containing a provider with no virial must "
                          "report virial_valid == false (ML ") +
                  (ml_first ? "first" : "second") + ")");
        check(result.potential_energy == 1.25,
              std::string("energies must still add when one child has no virial (ML ") +
                  (ml_first ? "first" : "second") + ")");
    }
}

// The gate that stops the absent value being consumed.
void test_pressure_barostats_cannot_consume_an_absent_virial() {
    gmd::BerendsenBarostat berendsen(2000.0, 4.5e-5);
    check(berendsen.requires_virial(),
          "the Berendsen barostat must declare that it needs a virial, otherwise "
          "the integrator has nothing to gate on");

    gmd::MCBarostat monte_carlo(25, 0.01, 100, 12345);
    check(!monte_carlo.requires_virial(),
          "the Monte Carlo barostat re-evaluates the energy and does not consume a "
          "virial, so it is not gated");

    // And the value a pressure-coupled barostat would have consumed is
    // genuinely absent once an ML result is installed on the System.
    gmd::System system = make_system();
    gmd::MLForceProvider provider("stub.pt", std::make_shared<StubModelAdapter>());
    gmd::RuntimeContext runtime;
    const auto coordinates = system.coordinates();
    gmd::ForceResult result;
    provider.compute(make_request(system, {coordinates.data(), coordinates.size()}),
                     result, runtime);

    system.set_provider_virial(result.virial, result.virial_valid);
    check(!system.last_virial_valid(),
          "a System fed an ML result must report no usable virial, which is what "
          "stops VelocityVerletIntegrator from running a barostat whose "
          "requires_virial() is true");

    // The same System must accept a real virial, or the check above would pass
    // for a System that never reports one at all.
    StubVirialProvider real_provider;
    gmd::ForceResult real_result;
    real_provider.compute(make_request(system, {coordinates.data(), coordinates.size()}),
                          real_result, runtime);
    system.set_provider_virial(real_result.virial, real_result.virial_valid);
    check(system.last_virial_valid(),
          "the System refuses a valid provider virial, so the negative result above "
          "says nothing about the ML path");
}

}  // namespace

int main() {
    test_fresh_result_reports_no_virial();
    test_stale_result_is_cleared();
    test_absent_virial_is_not_synthesised_from_the_force_moment();
    test_composite_with_an_ml_child_is_invalid();
    test_pressure_barostats_cannot_consume_an_absent_virial();

    if (failures != 0) {
        std::cerr << "ML virial contract tests failed: " << failures << '\n';
        return 1;
    }
    std::cout << "ML virial contract tests passed\n";
    return 0;
}
