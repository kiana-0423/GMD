#include "gmd/ml/ml_force_provider.hpp"

#include <utility>

#include "gmd/ml/model_runtime_adapter.hpp"

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
    result.success = false;
    result.potential_energy = 0.0;
    result.forces.clear();

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