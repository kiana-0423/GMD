#include "gmd/force/ml_force_provider.hpp"

#include <stdexcept>
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

void MLForceProvider::finalize(RuntimeContext& runtime) {
    if (adapter_ && model_loaded_) {
        adapter_->unload_model(runtime);
        model_loaded_ = false;
    }
}

}  // namespace gmd
