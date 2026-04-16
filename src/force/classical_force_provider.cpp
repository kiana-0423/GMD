#include "gmd/force/classical_force_provider.hpp"

namespace gmd {

std::string_view ClassicalForceProvider::name() const noexcept {
    return "classical_force_provider";
}

void ClassicalForceProvider::initialize(RuntimeContext& runtime) {
    (void)runtime;
}

void ClassicalForceProvider::compute(const ForceRequest& request,
                                     ForceResult& result,
                                     RuntimeContext& runtime) {
    (void)runtime;
    result.success = true;
    result.potential_energy = 0.0;
    result.forces.assign(request.coordinates.size(), Force3D{0.0, 0.0, 0.0});
    result.virial_valid = false;
}

void ClassicalForceProvider::finalize(RuntimeContext& runtime) {
    (void)runtime;
}

}  // namespace gmd