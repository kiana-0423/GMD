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

	void initialize(RuntimeContext& runtime) override;
	void compute(const ForceRequest& request, ForceResult& result, RuntimeContext& runtime) override;
	void finalize(RuntimeContext& runtime) override;

	void set_model_path(std::filesystem::path model_path) noexcept;
	const std::filesystem::path& model_path() const noexcept;

	void set_model_format(std::string model_format) noexcept;
	std::string_view model_format() const noexcept;

	void set_runtime_adapter(std::shared_ptr<ModelRuntimeAdapter> adapter) noexcept;
	const std::shared_ptr<ModelRuntimeAdapter>& runtime_adapter() const noexcept;

private:
	std::filesystem::path model_path_;
	std::string model_format_;
	std::shared_ptr<ModelRuntimeAdapter> adapter_;
	bool model_loaded_ = false;
};

}  // namespace gmd
