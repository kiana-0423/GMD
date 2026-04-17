#pragma once

#include "gmd/integrator/integrator.hpp"

namespace gmd {

// Default second-order integrator for deterministic molecular dynamics.
class VelocityVerletIntegrator final : public Integrator {
public:
    explicit VelocityVerletIntegrator(double dt) noexcept;
    ~VelocityVerletIntegrator() override = default;

    std::string_view name() const noexcept override;

    void initialize(System& system, RuntimeContext& runtime) override;
    void step(System& system,
              ForceProvider& force_provider,
              const IntegratorStepContext& ctx,
              RuntimeContext& runtime) override;

    // Returns the configured default time step for this integrator instance.
    double dt() const noexcept;
    void set_dt(double dt) noexcept;

private:
    double dt_ = 0.0;
};

}  // namespace gmd
