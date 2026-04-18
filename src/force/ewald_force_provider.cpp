#include "gmd/force/ewald_force_provider.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "gmd/boundary/minimum_image.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace gmd {

// Coulomb constant k_e = e² / (4π ε₀) in units of [eV · Å / e²].
static constexpr double kEwaldCoulomb = 14.3996;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

EwaldForceProvider::EwaldForceProvider(double alpha, int kmax,
                                       double real_cutoff) noexcept
    : alpha_(alpha),
      alpha_sq_(alpha * alpha),
      real_cutoff_(real_cutoff),
      real_cutoff_sq_(real_cutoff * real_cutoff),
      kmax_request_(kmax),
      kmax_resolved_(kmax) {}

std::string_view EwaldForceProvider::name() const noexcept {
    return "ewald_force_provider";
}
void EwaldForceProvider::initialize(RuntimeContext&) {}
void EwaldForceProvider::finalize(RuntimeContext&) {}

// ---------------------------------------------------------------------------
// Parameter auto-selection
// ---------------------------------------------------------------------------

void EwaldForceProvider::resolve_params(const Box& box) noexcept {
    // Auto real-space cutoff: 45 % of the shortest box dimension.
    if (real_cutoff_ <= 0.0) {
        real_cutoff_ = std::min({box.lengths[0], box.lengths[1], box.lengths[2]}) * 0.45;
        real_cutoff_sq_ = real_cutoff_ * real_cutoff_;
    }
    // Auto alpha: ~3.2 / r_cut gives good accuracy with modest kmax.
    if (alpha_ <= 0.0) {
        alpha_    = 3.2 / real_cutoff_;
        alpha_sq_ = alpha_ * alpha_;
    }
    // Auto kmax: ceil(alpha * L_max * 3.5 / pi), at least 3.
    if (kmax_request_ <= 0) {
        const double L_max = std::max({box.lengths[0], box.lengths[1], box.lengths[2]});
        kmax_resolved_ = std::max(3,
            static_cast<int>(std::ceil(alpha_ * L_max * 3.5 / std::numbers::pi)));
    } else {
        kmax_resolved_ = kmax_request_;
    }
}

// ---------------------------------------------------------------------------
// Main compute entry-point
// ---------------------------------------------------------------------------

void EwaldForceProvider::compute(const ForceRequest& req,
                                  ForceResult& res,
                                  RuntimeContext&) {
    const std::size_t n = req.coordinates.size();
    res.success = true;
    res.potential_energy = 0.0;
    res.forces.assign(n, Force3D{0.0, 0.0, 0.0});
    res.virial_valid = false;

    if (n == 0 || req.box == nullptr || req.system == nullptr) return;

    const auto charges = req.system->charges();
    if (charges.size() != n) return;

    // Early exit when all charges are zero.
    bool has_charges = false;
    for (std::size_t i = 0; i < n; ++i) {
        if (charges[i] != 0.0) { has_charges = true; break; }
    }
    if (!has_charges) return;

    resolve_params(*req.box);

    compute_real_space(req, res);
    compute_reciprocal(req, res);
    compute_self_correction(req, res);
}

// ---------------------------------------------------------------------------
// Real-space contribution
// ---------------------------------------------------------------------------
//
//   U_real = k_e * Σ_{i<j} q_i q_j * erfc(α r_ij) / r_ij
//
//   F_i from pair (i,j):
//     f_factor = k_e q_i q_j / r² * [erfc(αr)/r + (2α/√π) exp(-α²r²)]
//     F_i += f_factor * dr_ij    (dr_ij = r_i - r_j, minimum-image)

void EwaldForceProvider::compute_real_space(const ForceRequest& req,
                                             ForceResult& res) const {
    const auto& coords  = req.coordinates;
    const auto  charges = req.system->charges();
    const Box&  box     = *req.box;
    const std::size_t n = coords.size();

    const double two_alpha_over_sqrt_pi =
        2.0 * alpha_ / std::sqrt(std::numbers::pi);

    auto eval_pair = [&](std::size_t i, std::size_t j) {
        const double qi = charges[i], qj = charges[j];
        if (qi == 0.0 && qj == 0.0) return;

        Force3D dr = {coords[i][0] - coords[j][0],
                      coords[i][1] - coords[j][1],
                      coords[i][2] - coords[j][2]};
        apply_minimum_image(dr, box);

        const double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        if (r2 >= real_cutoff_sq_ || r2 < 1e-12) return;

        const double r       = std::sqrt(r2);
        const double ar      = alpha_ * r;
        const double erfc_ar = std::erfc(ar);
        const double exp_ar2 = std::exp(-ar * ar);

        // Energy contribution (Newton III: count once).
        res.potential_energy += kEwaldCoulomb * qi * qj * erfc_ar / r;

        // Force factor (∂/∂r of k_e*qi*qj*erfc(αr)/r) · 1/r
        const double ff = kEwaldCoulomb * qi * qj / r2
                          * (erfc_ar / r + two_alpha_over_sqrt_pi * exp_ar2);
        res.forces[i][0] += ff * dr[0];
        res.forces[i][1] += ff * dr[1];
        res.forces[i][2] += ff * dr[2];
        res.forces[j][0] -= ff * dr[0];
        res.forces[j][1] -= ff * dr[1];
        res.forces[j][2] -= ff * dr[2];
    };

    if (req.neighbor_list != nullptr && req.neighbor_list->valid) {
        const NeighborList& nl = *req.neighbor_list;
        for (std::size_t i = 0; i < n; ++i) {
            const int start = nl.offsets[i];
            const int count = nl.counts[i];
            for (int k = 0; k < count; ++k)
                eval_pair(i, static_cast<std::size_t>(nl.neighbors[start + k]));
        }
    } else {
        for (std::size_t i = 0; i < n - 1; ++i)
            for (std::size_t j = i + 1; j < n; ++j)
                eval_pair(i, j);
    }
}

// ---------------------------------------------------------------------------
// Reciprocal-space contribution
// ---------------------------------------------------------------------------
//
//   S(k) = Σ_j q_j exp(i k·r_j)   (structure factor)
//
//   U_recip = (k_e / 2V) Σ_{k≠0} (4π/k²) exp(-k²/4α²) |S(k)|²
//
//   F_i = (k_e / V) Σ_{k≠0} (4π k / k²) exp(-k²/4α²)
//            × q_i [S_re sin(k·r_i) - S_im cos(k·r_i)]
//
// The sum runs over all integer triples (nx,ny,nz) in [-kmax,kmax]³ \ {0}.
// The k±k symmetry means both (+k) and (-k) are included, so the 1/2 in the
// energy formula is already accounted for by iterating all k-vectors.

void EwaldForceProvider::compute_reciprocal(const ForceRequest& req,
                                             ForceResult& res) const {
    const auto& coords  = req.coordinates;
    const auto  charges = req.system->charges();
    const Box&  box     = *req.box;
    const std::size_t n = coords.size();

    const double Lx = box.lengths[0], Ly = box.lengths[1], Lz = box.lengths[2];
    const double V  = Lx * Ly * Lz;
    const double tpLx = 2.0 * std::numbers::pi / Lx;
    const double tpLy = 2.0 * std::numbers::pi / Ly;
    const double tpLz = 2.0 * std::numbers::pi / Lz;
    const double inv_4a2 = 1.0 / (4.0 * alpha_sq_);
    const int km = kmax_resolved_;

    for (int nx = -km; nx <= km; ++nx) {
        const double kx = nx * tpLx;
        for (int ny = -km; ny <= km; ++ny) {
            const double ky = ny * tpLy;
            for (int nz = -km; nz <= km; ++nz) {
                if (nx == 0 && ny == 0 && nz == 0) continue;
                const double kz  = nz * tpLz;
                const double k2  = kx*kx + ky*ky + kz*kz;

                // Damping and prefactor.
                const double gfactor = kEwaldCoulomb
                                     * (4.0 * std::numbers::pi)
                                     / (V * k2)
                                     * std::exp(-k2 * inv_4a2);

                // Structure factor S(k).
                double S_re = 0.0, S_im = 0.0;
                for (std::size_t j = 0; j < n; ++j) {
                    const double phi = kx*coords[j][0]
                                     + ky*coords[j][1]
                                     + kz*coords[j][2];
                    S_re += charges[j] * std::cos(phi);
                    S_im += charges[j] * std::sin(phi);
                }

                // Energy: ½ · gfactor · |S|²
                res.potential_energy += 0.5 * gfactor * (S_re*S_re + S_im*S_im);

                // Force on each atom.
                for (std::size_t i = 0; i < n; ++i) {
                    if (charges[i] == 0.0) continue;
                    const double phi = kx*coords[i][0]
                                     + ky*coords[i][1]
                                     + kz*coords[i][2];
                    const double c = std::cos(phi), s = std::sin(phi);
                    const double fi = gfactor * charges[i]
                                    * (S_re * s - S_im * c);
                    res.forces[i][0] += fi * kx;
                    res.forces[i][1] += fi * ky;
                    res.forces[i][2] += fi * kz;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Self-energy and net-charge corrections (energy only, no forces)
// ---------------------------------------------------------------------------
//
//   U_self = -k_e · (α/√π) · Σ_i q_i²
//   U_net  = -k_e · π / (2 V α²) · (Σ_i q_i)²   [non-neutral systems only]

void EwaldForceProvider::compute_self_correction(const ForceRequest& req,
                                                  ForceResult& res) const {
    const auto charges = req.system->charges();
    const std::size_t n = charges.size();

    double q2_sum = 0.0, Q_net = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        q2_sum += charges[i] * charges[i];
        Q_net  += charges[i];
    }

    res.potential_energy -=
        kEwaldCoulomb * (alpha_ / std::sqrt(std::numbers::pi)) * q2_sum;

    if (Q_net != 0.0) {
        const double V = req.box->lengths[0]
                       * req.box->lengths[1]
                       * req.box->lengths[2];
        res.potential_energy -=
            kEwaldCoulomb * std::numbers::pi / (2.0 * V * alpha_sq_) * Q_net * Q_net;
    }
}

}  // namespace gmd
