#include "gmd/force/ml_force_provider.hpp"

#include <stdexcept>
#include <string>
#include <utility>

#include "gmd/force/model_runtime_adapter.hpp"
#include "gmd/core/runtime_context.hpp"

namespace gmd {

MLForceProvider::MLForceProvider(std::filesystem::path model_path,
                                 std::shared_ptr<ModelRuntimeAdapter> adapter) noexcept
    : model_path_(std::move(model_path)),
      adapter_(adapter ? std::move(adapter) : CreateUnavailableModelRuntimeAdapter()) {}

std::string_view MLForceProvider::name() const noexcept {
    return "ml_force_provider";
}

float MLForceProvider::cutoff() const noexcept {
    if (adapter_) return adapter_->cutoff();
    return 0.0f;
}

void MLForceProvider::initialize(RuntimeContext& runtime) {
    if (runtime.size() > 1) {
        throw std::runtime_error(
            "MLForceProvider does not support MPI domain decomposition: "
            "local-plus-ghost model energy ownership and message-passing halo depth "
            "are not defined");
    }

    if (!adapter_) {
        adapter_ = CreateUnavailableModelRuntimeAdapter();
    }
    if (!model_loaded_ && !model_path_.empty()) {
        adapter_->load_model(model_path_, runtime);
        model_loaded_ = true;
    }
}

void MLForceProvider::compute(const ForceRequest& request,
                              ForceResult& result,
                              RuntimeContext& runtime) {
    if (runtime.size() > 1) {
        throw std::runtime_error(
            "MLForceProvider does not support MPI domain decomposition: "
            "local-plus-ghost model energy ownership and message-passing halo depth "
            "are not defined");
    }

    result.success = false;
    result.potential_energy = 0.0;
    result.forces.clear();

    // No ML backend in this engine reports a stress tensor, so this provider
    // has no virial to offer. Say so explicitly rather than leaving the caller's
    // fields untouched: a reused ForceResult would otherwise carry a previous
    // provider's tensor and have it attributed to the model.
    //
    // sum_i r_i (x) F_i is deliberately NOT synthesised here. For a periodic,
    // cell-dependent model that expression is not the virial -- the same reason
    // it is wrong for the Ewald reciprocal term -- and manufacturing one would
    // let a pressure-coupled barostat consume a number that does not mean what
    // it claims. CompositeForceProvider propagates virial_valid == false, and
    // VelocityVerletIntegrator refuses to run a barostat whose
    // requires_virial() is true while none is available.
    result.virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    result.virial_valid = false;

    if (!adapter_ || model_path_.empty()) {
        return;
    }

    if (!model_loaded_) {
        initialize(runtime);
    }

    // Refuse unknown species BEFORE the model sees them. The first force
    // evaluation happens in Simulation::initialize(), so this fires before
    // any MD step is taken.
    if (request.system != nullptr) {
        validate_species(request.system->atomic_numbers(),
                         request.coordinates.size(),
                         request.system->atom_tags());
    } else {
        validate_species({}, request.coordinates.size());
    }

    ModelEvaluationRequest model_request{
        .model_path = model_path_,
        .model_format = model_format_,
        .coordinates = request.coordinates,
        .box = request.box,
        .atomic_numbers = request.system ? request.system->atomic_numbers()
                                         : std::span<const int>{},
        .neighbor_list = request.neighbor_list,
    };
    ModelEvaluationResult model_result;
    adapter_->evaluate(model_request, model_result, runtime);

    result.success = model_result.success;
    result.potential_energy = model_result.total_energy;
    result.forces = std::move(model_result.forces);
}

void MLForceProvider::validate_species(std::span<const int> atomic_numbers,
                                       std::size_t atom_count,
                                       std::span<const int> tags) {
    if (atomic_numbers.size() < atom_count) {
        throw std::runtime_error(
            "ML force field: atomic numbers are unavailable for " +
            std::to_string(atom_count - atomic_numbers.size()) + " of " +
            std::to_string(atom_count) +
            " atom(s); the model needs the species of every atom");
    }
    std::size_t bad = 0;
    std::string examples;
    for (std::size_t i = 0; i < atom_count; ++i) {
        if (atomic_numbers[i] > 0) continue;
        if (bad < 5) {
            const int id = (i < tags.size()) ? tags[i] : static_cast<int>(i);
            examples += (bad == 0 ? "" : ", ") + std::string("tag ") + std::to_string(id) +
                        " (Z=" + std::to_string(atomic_numbers[i]) + ")";
        }
        ++bad;
    }
    if (bad == 0) return;
    throw std::runtime_error(
        "ML force field: " + std::to_string(bad) + " of " + std::to_string(atom_count) +
        " atom(s) have unknown species (atomic number 0): " + examples +
        (bad > 5 ? ", ..." : "") +
        ". Give every atom an element symbol GMD recognises in the xyz file (or a force "
        "field type mapping to one), and restart ML runs only from version 3+ checkpoints.");
}

void MLForceProvider::finalize(RuntimeContext& runtime) {
    if (adapter_ && model_loaded_) {
        adapter_->unload_model(runtime);
        model_loaded_ = false;
    }
}

}  // namespace gmd
