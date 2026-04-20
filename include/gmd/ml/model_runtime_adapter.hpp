#pragma once

#include <filesystem>
#include <memory>
#include <span>
#include <string_view>
#include <vector>

#include "../force/force_provider.hpp"
#include "../system/system.hpp"

namespace gmd {

class RuntimeContext;

// Request passed to an ML runtime after the model artifact has been selected.
struct ModelEvaluationRequest {
	std::filesystem::path model_path;
	std::string_view model_format;
	std::span<const Coordinate3D> coordinates{};
	const Box* box = nullptr;
	// Atomic numbers (Z) per atom, required by SE3-GNN and similar models.
	std::span<const int> atomic_numbers{};
	// Pre-built neighbor list; non-null when a NeighborBuilder is active.
	const NeighborList* neighbor_list = nullptr;
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

	// Returns the force cutoff radius [Å] read from the model artifact (0 = unknown).
	virtual float cutoff() const noexcept { return 0.0f; }

	virtual void load_model(const std::filesystem::path& model_path, RuntimeContext& runtime) = 0;
	virtual void evaluate(const ModelEvaluationRequest& request,
						  ModelEvaluationResult& result,
						  RuntimeContext& runtime) = 0;
	virtual void unload_model(RuntimeContext& runtime) = 0;
};

// Returns a placeholder runtime that accepts a model path and produces zeroed outputs.
std::shared_ptr<ModelRuntimeAdapter> CreateUnavailableModelRuntimeAdapter();

}  // namespace gmd
