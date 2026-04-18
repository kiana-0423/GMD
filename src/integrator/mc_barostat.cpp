#include "gmd/integrator/mc_barostat.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <vector>

#include "gmd/force/force_provider.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/system.hpp"

namespace gmd {

MCBarostat::MCBarostat(std::uint32_t frequency,
                       double        max_delta_ln_V,
                       std::uint32_t adjust_interval,
                       std::uint32_t seed) noexcept
    : frequency_(frequency),
      max_delta_ln_V_(max_delta_ln_V),
      adjust_interval_(adjust_interval),
      rng_(seed)
{}

void MCBarostat::reset() {
    n_attempts_     = 0;
    n_accepted_     = 0;
    window_attempts_ = 0;
    window_accepted_ = 0;
}

void MCBarostat::apply(System& system,
                       ForceProvider& provider,
                       RuntimeContext& runtime,
                       std::uint64_t step,
                       double /*dt*/,
                       double temperature,
                       double target_pressure,
                       double /*virial_trace*/)
{
    // Only attempt a volume move every `frequency_` steps.
    if (frequency_ == 0 || (step % static_cast<std::uint64_t>(frequency_) != 0))
        return;

    const std::size_t n = system.atom_count();
    if (n == 0) return;

    // --- Physical constants for this move ---
    const double kT       = kB_eV * temperature;
    if (kT <= 0.0) return;
    const double P_ext    = target_pressure * bar_to_eV_per_A3;   // [eV/Å³]

    // --- Save current state ---
    const Box old_box          = system.box();
    const double V_old         = old_box.lengths[0] * old_box.lengths[1] * old_box.lengths[2];
    if (V_old < 1.0e-30) return;

    // Copy positions.
    const std::span<const std::array<double,3>> coords_view = system.coordinates();
    const std::vector<std::array<double,3>> old_coords(coords_view.begin(), coords_view.end());

    const double U_old = system.potential_energy();
    const bool   nl_valid_old = system.neighbor_list().valid;

    // --- Propose isotropic volume change ---
    const double delta     = max_delta_ln_V_ * uniform_(rng_);   // Δ ln V ~ U(-Δ_max, +Δ_max)
    const double ln_ratio  = delta;                               // ln(V_new / V_old)
    const double mu        = std::exp(ln_ratio / 3.0);            // linear scaling factor

    // Scale box.
    const std::array<double,3> new_lengths = {
        old_box.lengths[0] * mu,
        old_box.lengths[1] * mu,
        old_box.lengths[2] * mu
    };
    system.mutable_box().set_lengths(new_lengths);

    // Scale all coordinates.
    {
        auto coords = system.mutable_coordinates();
        for (auto& c : coords) {
            c[0] *= mu;
            c[1] *= mu;
            c[2] *= mu;
        }
    }

    // Invalidate neighbor list so the force provider does not use stale data.
    system.mutable_neighbor_list().valid = false;

    // --- Evaluate energy at trial volume ---
    const std::span<const std::array<double,3>> new_coords_view = system.coordinates();
    ForceRequest req;
    req.system        = &system;
    req.box           = &system.box();
    req.step          = step;
    req.time          = 0.0;
    req.coordinates   = std::span<const std::array<double,3>>(
                            new_coords_view.data(), new_coords_view.size());
    req.neighbor_list = nullptr;   // force O(N²) or provider-internal rebuild

    ForceResult res;
    res.success = false;
    provider.compute(req, res, runtime);

    const double U_new = res.success ? res.potential_energy : U_old + 1.0e30;

    // --- Metropolis acceptance criterion ---
    // w = ΔU + P_ext * ΔV - N * kT * ln(V_new/V_old)
    const double V_new = new_lengths[0] * new_lengths[1] * new_lengths[2];
    const double dU    = U_new - U_old;
    const double dV    = V_new - V_old;
    const double w     = dU + P_ext * dV - static_cast<double>(n) * kT * ln_ratio;

    const bool accept = (w <= 0.0) || (uniform01_(rng_) < std::exp(-w / kT));

    ++n_attempts_;
    ++window_attempts_;

    if (accept) {
        // Keep new state.  Forces were not evaluated with a proper neighbor list,
        // so set potential_energy from the trial result and leave nl invalid so
        // the next neighbor build is triggered at the start of the next step.
        system.set_potential_energy(U_new);
        ++n_accepted_;
        ++window_accepted_;
    } else {
        // Restore old box and coordinates.
        system.mutable_box().set_lengths(old_box.lengths);
        {
            auto coords = system.mutable_coordinates();
            for (std::size_t i = 0; i < n; ++i)
                coords[i] = old_coords[i];
        }
        system.set_potential_energy(U_old);
        // Restore neighbor list validity — positions are back to the reference.
        system.mutable_neighbor_list().valid = nl_valid_old;
    }

    // --- Adaptive step-size control (target 50 % acceptance) ---
    if (adjust_interval_ > 0 && window_attempts_ >= adjust_interval_) {
        const double frac = static_cast<double>(window_accepted_)
                          / static_cast<double>(window_attempts_);
        if (frac > 0.5) {
            max_delta_ln_V_ = std::min(max_delta_ln_V_ * adapt_factor_, max_delta_hard_);
        } else {
            max_delta_ln_V_ = std::max(max_delta_ln_V_ / adapt_factor_, min_delta_);
        }
        window_attempts_ = 0;
        window_accepted_ = 0;
    }
}

}  // namespace gmd
