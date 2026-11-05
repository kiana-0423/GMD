// MLForceProvider refuses atoms of unknown species.
//
// GMD sets Z = 0 when it cannot identify an atom's element: an element symbol
// outside its table, an integer type with no element mapping, or (before
// checkpoint version 3) every atom of a restarted run. ML models use 0 as the
// padding species, so evaluating such an atom does not fail -- it silently
// returns a different, wrong energy. The provider must therefore refuse to
// call the model at all, with an error that names the offending atoms.

#include <filesystem>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/force/ml_force_provider.hpp"
#include "gmd/force/model_runtime_adapter.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[ml species guard] " << message << '\n';
        ++failures;
    }
}

class CountingAdapter final : public gmd::ModelRuntimeAdapter {
public:
    int evaluations = 0;
    std::string_view name() const noexcept override { return "counting"; }
    void load_model(const std::filesystem::path&, gmd::RuntimeContext&) override {}
    void unload_model(gmd::RuntimeContext&) override {}
    void evaluate(const gmd::ModelEvaluationRequest& request,
                  gmd::ModelEvaluationResult& result,
                  gmd::RuntimeContext&) override {
        ++evaluations;
        result.success = true;
        result.total_energy = -1.0;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
    }
};

gmd::System make_system(const std::vector<int>& z) {
    gmd::System system;
    system.resize(z.size(), z.size());
    gmd::Box box;
    box.set_lengths({20.0, 20.0, 20.0});
    system.set_box(box);
    for (std::size_t i = 0; i < z.size(); ++i) {
        system.mutable_coordinates()[i] = {1.0 + 2.0 * static_cast<double>(i), 5.0, 5.0};
        system.mutable_masses()[i] = 1.0;
        system.mutable_atomic_numbers()[i] = z[i];
        system.mutable_atom_tags()[i] = 100 + static_cast<int>(i);
    }
    return system;
}

std::string evaluate(const std::vector<int>& z, int& evaluations) {
    gmd::System system = make_system(z);
    auto adapter = std::make_shared<CountingAdapter>();
    gmd::MLForceProvider provider("stub.pt", adapter);
    gmd::RuntimeContext runtime;
    const auto coordinates = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .step = 0,
        .time = 0.0,
        .coordinates = {coordinates.data(), coordinates.size()},
        .neighbor_list = nullptr,
    };
    gmd::ForceResult result;
    std::string error;
    try {
        provider.compute(request, result, runtime);
    } catch (const std::runtime_error& e) {
        error = e.what();
    }
    evaluations = adapter->evaluations;
    return error;
}

}  // namespace

int main() {
    int evaluations = 0;

    std::string error = evaluate({8, 1, 1}, evaluations);
    check(error.empty(), "known species were refused: " + error);
    check(evaluations == 1, "known species were not evaluated");

    error = evaluate({8, 0, 1, 0}, evaluations);
    check(!error.empty(), "atoms with Z = 0 were accepted");
    check(evaluations == 0, "the model was called despite unknown species");
    check(error.find("unknown species") != std::string::npos,
          "the error does not say 'unknown species': " + error);
    check(error.find("2 of 4") != std::string::npos,
          "the error does not count the offending atoms: " + error);
    check(error.find("tag 101") != std::string::npos &&
              error.find("tag 103") != std::string::npos,
          "the error does not name the offending atoms by tag: " + error);

    error = evaluate({-3}, evaluations);
    check(!error.empty() && evaluations == 0, "a negative atomic number was accepted");

    // The static check on its own, including a short species array.
    bool threw = false;
    try {
        const std::vector<int> z = {1, 1};
        gmd::MLForceProvider::validate_species(z, 3);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    check(threw, "a species array shorter than the atom count was accepted");

    if (failures != 0) {
        std::cerr << "[ml species guard] " << failures << " failure(s)\n";
        return 1;
    }
    std::cout << "[ml species guard] all checks passed\n";
    return 0;
}
