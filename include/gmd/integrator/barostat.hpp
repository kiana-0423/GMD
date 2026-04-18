#pragma once

namespace gmd {

class System;

// Abstract interface for pressure control (barostats).
class Barostat {
public:
    virtual ~Barostat() = default;

    // Apply barostat scaling after a completed MD step.
    // virial_trace: sum of the diagonal components of the virial tensor (W_xx + W_yy + W_zz).
    virtual void apply(System& system, double dt, double target_pressure,
                       double virial_trace) = 0;

    // Optional: reset any internal state (e.g. between equilibration and production).
    virtual void reset() {}
};

}  // namespace gmd
