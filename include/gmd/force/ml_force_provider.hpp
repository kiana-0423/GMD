#pragma once

#include <filesystem>
#include <cstddef>
#include <memory>
#include <span>
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

	// Throws std::runtime_error naming the offending atoms if any of the first
	// `atom_count` atomic numbers is <= 0, or if fewer than `atom_count` are
	// supplied. Z = 0 means GMD could not identify the species (an element
	// symbol outside its table, or an integer type with no element mapping, or
	// a pre-v3 checkpoint); ML models treat 0 as padding, so evaluating it
	// would silently produce wrong energies and forces. `tags` (optional, same
	// length) is used only to name atoms in the message.
	static void validate_species(std::span<const int> atomic_numbers,
	                             std::size_t atom_count,
	                             std::span<const int> tags = {});

private:
	std::filesystem::path model_path_;
	std::string model_format_;
	std::shared_ptr<ModelRuntimeAdapter> adapter_;
	bool model_loaded_ = false;
};

}  // namespace gmd
