#pragma once

#include <cstdint>
#include <string_view>

namespace gmd_next::model {

// Truncation conventions, following docs/operator_contracts.md section 3:
//   none             u = phi(r) everywhere; only defined without periodicity
//   potential_shift  u = phi(r) - phi(rc), the force is unchanged below rc
//   force_shift      u = phi(r) - phi(rc) - (r - rc) phi'(rc), radial derivative
//                    reduced by phi'(rc), so the force vanishes at rc
// These are different models. Nothing may switch between them to fit a backend;
// comparisons with the current GMD Lennard-Jones use potential_shift.
enum class CutoffMode : std::uint8_t { none, potential_shift, force_shift };

std::string_view name(CutoffMode mode);

// Single-type Lennard-Jones: one epsilon and one sigma for every pair.
// Under UnitSystem::metal epsilon is eV and sigma, cutoff and skin are Angstrom.
struct LjModel {
    double epsilon = 0.0;
    double sigma = 0.0;
    CutoffMode cutoff_mode = CutoffMode::none;
    // Exactly zero when cutoff_mode is none, otherwise finite and positive.
    double cutoff = 0.0;
    // Neighbour-list margin. Positive in the production scope; exactly zero for
    // static reference evaluation, which builds no list and cannot consume it.
    double skin = 0.0;
};

// The list radius rc + skin that a neighbour build must cover.
double list_radius(const LjModel& model);

}  // namespace gmd_next::model
