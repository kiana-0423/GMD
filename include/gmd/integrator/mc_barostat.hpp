#pragma once

#include <cstdint>
#include <random>

#include "gmd/integrator/barostat.hpp"

namespace gmd {

// Monte Carlo (isotropic) barostat for the NPT ensemble.
//
// Every `frequency` MD steps, a trial isotropic volume move is proposed and
// accepted or rejected via the Metropolis criterion:
//
//   w = dU + P_ext * dV - N * kT * ln(V_new / V_old)
//   P(accept) = min(1, exp(-w / kT))
//
// where kT = k_B * temperature, P_ext is the target pressure, and the ln term
// accounts for the change in translational phase-space volume.
//
// This formulation does NOT require the virial and therefore works correctly
// with any combination of force providers (LJ, Ewald, PME, bonded, ML).
//
// Step-size adaptation: every `adjust_interval` attempts, the maximum
// ln-volume displacement is scaled to drive the acceptance rate toward 50 %.
//
// Reference: Frenkel & Smit, "Understanding Molecular Simulation", 2nd ed.,
//   Algorithm 5.1 (isotropic NPT MC move).
class MCBarostat final : public Barostat {
public:
    // frequency       — attempt a volume move every this many MD steps
    // max_delta_ln_V  — initial maximum |Δ ln V| for trial moves (dimensionless)
    // adjust_interval — number of attempts between step-size adaptation events
    // seed            — RNG seed (deterministic by default)
    explicit MCBarostat(std::uint32_t frequency      = 25,
                        double        max_delta_ln_V  = 0.01,
                        std::uint32_t adjust_interval = 100,
                        std::uint32_t seed            = 12345) noexcept;

    ~MCBarostat() override = default;

    // --- Barostat interface ---
    void apply(System& system,
               ForceProvider& provider,
               RuntimeContext& runtime,
               std::uint64_t step,
               double dt,
               double temperature,
               double target_pressure,
               double virial_trace) override;

    // MC barostat does not require the virial tensor.
    bool requires_virial() const noexcept override { return false; }

    void reset() override;

    // --- Diagnostics ---
    std::uint64_t attempts()  const noexcept { return n_attempts_; }
    std::uint64_t accepted()  const noexcept { return n_accepted_; }
    double acceptance_rate()  const noexcept {
        return n_attempts_ > 0
            ? static_cast<double>(n_accepted_) / static_cast<double>(n_attempts_)
            : 0.0;
    }
    double max_delta_ln_V()   const noexcept { return max_delta_ln_V_; }

    void set_frequency(std::uint32_t f)       noexcept { frequency_       = f; }
    void set_adjust_interval(std::uint32_t n) noexcept { adjust_interval_ = n; }

private:
    // Physical constants (internal unit system: eV, Å, amu, fs).
    static constexpr double kB_eV           = 8.617333262e-5;  // Boltzmann [eV/K]
    static constexpr double bar_to_eV_per_A3 = 6.2415091e-7;   // 1 bar → eV/Å³

    // Adaptive step-size bounds [dimensionless ln-volume displacement].
    static constexpr double min_delta_ = 1.0e-5;
    static constexpr double max_delta_hard_ = 0.3;
    static constexpr double adapt_factor_   = 1.1;   // multiply/divide by this

    std::uint32_t frequency_       = 25;
    double        max_delta_ln_V_  = 0.01;
    std::uint32_t adjust_interval_ = 100;

    // Counters (since last reset).
    std::uint64_t n_attempts_  = 0;
    std::uint64_t n_accepted_  = 0;
    // Counters within the current adapt window.
    std::uint32_t window_attempts_ = 0;
    std::uint32_t window_accepted_ = 0;

    std::mt19937                      rng_;
    std::uniform_real_distribution<double> uniform_{-1.0, 1.0};
    std::uniform_real_distribution<double> uniform01_{0.0, 1.0};
};

}  // namespace gmd
