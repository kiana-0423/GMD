#pragma once

#include <filesystem>
#include <memory>
#include <span>
#include <string_view>
#include <vector>

#include "../force/force_provider.hpp"

namespace gmd {

class RuntimeContext;

// Request passed to an ML runtime after the model artifact has been selected.
struct ModelEvaluationRequest {
	std::filesystem::path model_path;
	std::string_view model_format;
	std::span<const Coordinate3D> coordinates{};
	const Box* box = nullptr;
};

// Raw outputs returned by an ML model runtime.
struct ModelEvaluationResult {
	bool success = false;
	double total_energy = 0.0;
	std::vector<Force3D> forces;
};

// Backend adapter for TorchScript, ONNX, custom CUDA runtimes, and similar model containers.
class ModelRuntimeAdapter {
public:
	virtual ~ModelRuntimeAdapter() = default;

	virtual std::string_view name() const noexcept = 0;

	virtual void load_model(const std::filesystem::path& model_path, RuntimeContext& runtime) = 0;
	virtual void evaluate(const ModelEvaluationRequest& request,
						  ModelEvaluationResult& result,
						  RuntimeContext& runtime) = 0;
	virtual void unload_model(RuntimeContext& runtime) = 0;
};

// Returns a placeholder runtime that accepts a model path and produces zeroed outputs.
std::shared_ptr<ModelRuntimeAdapter> CreateUnavailableModelRuntimeAdapter();

}  // namespace gmd
