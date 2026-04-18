#pragma once

#include <cstdint>

namespace gmd {

class ForceProvider;
class RuntimeContext;
class System;

// Abstract interface for pressure control (barostats).
//
// apply() is called once per MD step with full access to the force provider
// and runtime so that implementations that require energy re-evaluation
// (e.g. Monte Carlo barostat) can do so without coupling to the integrator.
class Barostat {
public:
    virtual ~Barostat() = default;

    // Apply pressure control after a completed velocity-Verlet step.
    //
    // Parameters:
    //   system          — the simulation system (may be modified in-place)
    //   provider        — force backend (used by MC barostat for trial moves)
    //   runtime         — runtime context forwarded to force provider
    //   step            — current MD step index
    //   dt              — integration time step [fs]
    //   temperature     — target temperature [K]  (used as kT in MC criterion)
    //   target_pressure — target pressure [bar]
    //   virial_trace    — W_xx + W_yy + W_zz from the last force evaluation
    //                     (valid only when ForceResult::virial_valid is true)
    virtual void apply(System& system,
                       ForceProvider& provider,
                       RuntimeContext& runtime,
                       std::uint64_t step,
                       double dt,
                       double temperature,
                       double target_pressure,
                       double virial_trace) = 0;

    // Returns true when this barostat requires a valid virial tensor to function
    // correctly (e.g. Berendsen).  Returns false when it works without virial
    // (e.g. Monte Carlo).  The integrator uses this to decide whether to skip
    // the barostat call when the active force provider has virial_valid == false.
    virtual bool requires_virial() const noexcept { return true; }

    // Optional: reset any internal state (e.g. between equilibration and production).
    virtual void reset() {}
};

}  // namespace gmd
