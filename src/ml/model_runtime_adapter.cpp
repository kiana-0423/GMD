#include "gmd/ml/model_runtime_adapter.hpp"

#include <memory>

namespace gmd {
namespace {

class UnavailableModelRuntimeAdapter final : public ModelRuntimeAdapter {
public:
    std::string_view name() const noexcept override {
        return "unavailable_model_runtime_adapter";
    }

    void load_model(const std::filesystem::path& model_path, RuntimeContext& runtime) override {
        (void)runtime;
        model_path_ = model_path;
        loaded_ = !model_path_.empty();
    }

    void evaluate(const ModelEvaluationRequest& request,
                  ModelEvaluationResult& result,
                  RuntimeContext& runtime) override {
        (void)runtime;
        result.success = loaded_ && request.model_path == model_path_;
        result.total_energy = 0.0;
        result.forces.assign(request.coordinates.size(), Force3D{0.0, 0.0, 0.0});
    }

    void unload_model(RuntimeContext& runtime) override {
        (void)runtime;
        loaded_ = false;
        model_path_.clear();
    }

private:
    std::filesystem::path model_path_;
    bool loaded_ = false;
};

}  // namespace

std::shared_ptr<ModelRuntimeAdapter> CreateUnavailableModelRuntimeAdapter() {
    return std::make_shared<UnavailableModelRuntimeAdapter>();
}

}  // namespace gmd