#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include "gmd/force/force_provider.hpp"

namespace gmd {

// Holds an ordered list of ForceProvider instances and accumulates their
// contributions into a single ForceResult.  Energies are summed; force vectors
// are summed element-wise.  Used to combine, e.g., LJ + Ewald/PME.
class CompositeForceProvider final : public ForceProvider {
public:
    CompositeForceProvider() = default;

    void add(std::shared_ptr<ForceProvider> provider) {
        providers_.push_back(std::move(provider));
    }

    std::string_view name() const noexcept override {
        return "composite_force_provider";
    }

    void initialize(RuntimeContext& runtime) override {
        for (auto& p : providers_) p->initialize(runtime);
    }

    void compute(const ForceRequest& request,
                 ForceResult& result,
                 RuntimeContext& runtime) override {
        const std::size_t n = request.coordinates.size();
        result.success = true;
        result.potential_energy = 0.0;
        result.forces.assign(n, Force3D{0.0, 0.0, 0.0});
        result.virial_valid = false;

        for (auto& p : providers_) {
            ForceResult sub;
            p->compute(request, sub, runtime);
            if (!sub.success) {
                result.success = false;
                return;
            }
            result.potential_energy += sub.potential_energy;
            for (std::size_t i = 0; i < n; ++i) {
                result.forces[i][0] += sub.forces[i][0];
                result.forces[i][1] += sub.forces[i][1];
                result.forces[i][2] += sub.forces[i][2];
            }
        }
    }

    void finalize(RuntimeContext& runtime) override {
        for (auto& p : providers_) p->finalize(runtime);
    }

private:
    std::vector<std::shared_ptr<ForceProvider>> providers_;
};

}  // namespace gmd
