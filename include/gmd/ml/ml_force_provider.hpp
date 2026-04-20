#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>

#include "../force/force_provider.hpp"

namespace gmd {

class ModelRuntimeAdapter;

// Force provider that evaluates a user-supplied ML force field artifact.
// The wrapped model consumes atomic coordinates and returns per-atom forces
// together with the total potential energy of the system.
class MLForceProvider final : public ForceProvider {
public:
	explicit MLForceProvider(std::filesystem::path model_path,
							 std::shared_ptr<ModelRuntimeAdapter> adapter = nullptr) noexcept;
	~MLForceProvider() override = default;

	std::string_view name() const noexcept override;

	// Returns the force cutoff radius [Å] as read from the model artifact.
	// Returns 0 if the backend does not expose this information.
	float cutoff() const noexcept;

	void initialize(RuntimeContext& runtime) override;
	void compute(const ForceRequest& request, ForceResult& result, RuntimeContext& runtime) override;
	void finalize(RuntimeContext& runtime) override;

private:
	std::filesystem::path model_path_;
	std::string model_format_;
	std::shared_ptr<ModelRuntimeAdapter> adapter_;
	bool model_loaded_ = false;
};

}  // namespace gmd
