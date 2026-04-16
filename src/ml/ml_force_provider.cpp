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

void MLForceProvider::set_model_path(std::filesystem::path model_path) noexcept {
    model_path_ = std::move(model_path);
    model_loaded_ = false;
}

const std::filesystem::path& MLForceProvider::model_path() const noexcept {
    return model_path_;
}

void MLForceProvider::set_model_format(std::string model_format) noexcept {
    model_format_ = std::move(model_format);
}

std::string_view MLForceProvider::model_format() const noexcept {
    return model_format_;
}

void MLForceProvider::set_runtime_adapter(std::shared_ptr<ModelRuntimeAdapter> adapter) noexcept {
    adapter_ = adapter ? std::move(adapter) : CreateUnavailableModelRuntimeAdapter();
    model_loaded_ = false;
}

const std::shared_ptr<ModelRuntimeAdapter>& MLForceProvider::runtime_adapter() const noexcept {
    return adapter_;
}

}  // namespace gmd