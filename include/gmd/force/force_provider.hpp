#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace gmd {

struct Box;
class RuntimeContext;
class System;

using Coordinate3D = std::array<double, 3>;
using Force3D = std::array<double, 3>;

// Immutable inputs supplied to a force evaluation at a given simulation step.
struct ForceRequest {
    const System* system = nullptr;
    const Box* box = nullptr;
    std::uint64_t step = 0;
    double time = 0.0;
    std::span<const Coordinate3D> coordinates;
};

// Aggregate outputs produced by a force evaluation.
struct ForceResult {
    bool success = false;
    double potential_energy = 0.0;
    std::vector<Force3D> forces;
    std::array<double, 9> virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    bool virial_valid = false;
};

// Abstract interface for classical, ML, or hybrid force backends.
class ForceProvider {
public:
    virtual ~ForceProvider() = default;

    virtual std::string_view name() const noexcept = 0;

    // Prepares backend-owned state, buffers, or model handles.
    virtual void initialize(RuntimeContext& runtime) = 0;
    // Computes forces and energy for the provided system state.
    virtual void compute(const ForceRequest& request, ForceResult& result, RuntimeContext& runtime) = 0;
    // Releases backend resources tied to the runtime.
    virtual void finalize(RuntimeContext& runtime) = 0;
};

}  // namespace gmd
