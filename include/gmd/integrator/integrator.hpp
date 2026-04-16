#pragma once

#include <cstdint>
#include <string_view>

namespace gmd {

struct ForceResult;
class RuntimeContext;
class System;

// Per-step metadata passed into an integration scheme.
struct IntegratorStepContext {
    std::uint64_t step = 0;
    double dt = 0.0;
};

// Abstract interface for time integration algorithms.
class Integrator {
public:
    virtual ~Integrator() = default;

    virtual std::string_view name() const noexcept = 0;

    // Initializes scheme-specific state before the first time step.
    virtual void initialize(System& system, RuntimeContext& runtime) = 0;
    // Updates positions and velocities using the latest force evaluation.
    virtual void step(System& system,
                      const ForceResult& force,
                      const IntegratorStepContext& ctx,
                      RuntimeContext& runtime) = 0;
};

}  // namespace gmd
