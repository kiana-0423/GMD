#pragma once

#include <cstdint>
#include <string_view>

namespace gmd {

class RuntimeContext;
class System;

// Optional diagnostics returned after rebuilding the neighbor data structure.
struct NeighborBuildStats {
    bool rebuilt = false;
    std::uint64_t pair_count = 0;
};

// Abstract interface for Verlet list or cell-list construction policies.
class NeighborBuilder {
public:
    virtual ~NeighborBuilder() = default;

    virtual std::string_view name() const noexcept = 0;

    // Performs one-time setup against the initial system state.
    virtual void initialize(System& system, RuntimeContext& runtime) = 0;
    // Checks whether particle motion invalidated the current neighbor representation.
    virtual bool needs_rebuild(const System& system, std::uint64_t step) const = 0;
    // Rebuilds neighbor data and optionally reports lightweight statistics.
    virtual void rebuild(System& system, RuntimeContext& runtime, NeighborBuildStats* stats) = 0;
};

}  // namespace gmd
