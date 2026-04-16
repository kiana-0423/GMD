#pragma once

#include <string_view>

#include "gmd/force/force_provider.hpp"

namespace gmd {

// Minimal placeholder for traditional analytical force fields.
class ClassicalForceProvider final : public ForceProvider {
public:
	ClassicalForceProvider() = default;
	~ClassicalForceProvider() override = default;

	std::string_view name() const noexcept override;

	void initialize(RuntimeContext& runtime) override;
	void compute(const ForceRequest& request, ForceResult& result, RuntimeContext& runtime) override;
	void finalize(RuntimeContext& runtime) override;
};

}  // namespace gmd
