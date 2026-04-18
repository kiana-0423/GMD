#include "gmd/force/pme_force_provider.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "gmd/boundary/minimum_image.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace gmd {

// Coulomb constant k_e = e² / (4π ε₀) [eV · Å / e²].
static constexpr double kPMECoulomb = 14.3996;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static bool is_power_of_two(int n) noexcept {
    return n > 0 && (n & (n - 1)) == 0;
}

// ---------------------------------------------------------------------------
// Construction / lifecycle
// ---------------------------------------------------------------------------

PMEForceProvider::PMEForceProvider(double alpha, double real_cutoff,
                                   int order, std::array<int, 3> grid)
    : alpha_(alpha),
      alpha_sq_(alpha * alpha),
      real_cutoff_(real_cutoff),
      real_cutoff_sq_(real_cutoff * real_cutoff),
      order_(order),
      grid_(grid) {
    if (order != 4 && order != 6)
        throw std::invalid_argument("PME B-spline order must be 4 or 6");
    for (int d = 0; d < 3; ++d) {
        if (!is_power_of_two(grid[d]))
            throw std::invalid_argument(
                "PME grid dimensions must each be a power of 2 (got "
                + std::to_string(grid[d]) + ")");
    }

    const std::size_t total = static_cast<std::size_t>(grid[0])
                            * static_cast<std::size_t>(grid[1])
                            * static_cast<std::size_t>(grid[2]);
    mesh_.resize(total);
    influence_.resize(total, 0.0);

    for (int d = 0; d < 3; ++d)
        bmod_sq_[d].resize(static_cast<std::size_t>(grid[d]), 0.0);
}

std::string_view PMEForceProvider::name() const noexcept {
    return "pme_force_provider";
}
void PMEForceProvider::initialize(RuntimeContext&) {}
void PMEForceProvider::finalize(RuntimeContext&) {}

// ---------------------------------------------------------------------------
// Parameter auto-selection
// ---------------------------------------------------------------------------

void PMEForceProvider::resolve_params(const Box& box) noexcept {
    if (real_cutoff_ <= 0.0) {
        real_cutoff_ = std::min({box.lengths[0], box.lengths[1], box.lengths[2]}) * 0.45;
        real_cutoff_sq_ = real_cutoff_ * real_cutoff_;
    }
    if (alpha_ <= 0.0) {
        alpha_    = 3.2 / real_cutoff_;
        alpha_sq_ = alpha_ * alpha_;
    }
}

// ---------------------------------------------------------------------------
// B-spline evaluation
// ---------------------------------------------------------------------------
//
// Cardinal B-spline M_p(u) for u in [0, p], recursively defined.
// Hardcoded analytical forms for order 4 and 6 for performance.

double PMEForceProvider::bspline(double u, int order) noexcept {
    if (u < 0.0 || u > static_cast<double>(order)) return 0.0;

    if (order == 4) {
        if      (u < 1.0) return u * u * u / 6.0;
        else if (u < 2.0) return (-3.0*u*u*u + 12.0*u*u - 12.0*u + 4.0) / 6.0;
        else if (u < 3.0) return ( 3.0*u*u*u - 24.0*u*u + 60.0*u - 44.0) / 6.0;
        else               return (4.0 - u) * (4.0 - u) * (4.0 - u) / 6.0;
    }
    // order == 6 (quintic B-spline)
    if (order == 6) {
        const double u2 = u * u, u3 = u2 * u, u4 = u3 * u, u5 = u4 * u;
        if      (u < 1.0) return u5 / 120.0;
        else if (u < 2.0) return (-5.0*u5 + 30.0*u4 - 60.0*u3 + 60.0*u2 - 30.0*u + 6.0) / 120.0;
        else if (u < 3.0) return (10.0*u5 - 90.0*u4 + 300.0*u3 - 480.0*u2 + 360.0*u - 102.0) / 120.0;
        else if (u < 4.0) return (-10.0*u5 + 120.0*u4 - 540.0*u3 + 1140.0*u2 - 1110.0*u + 402.0) / 120.0;
        else if (u < 5.0) return (5.0*u5 - 90.0*u4 + 600.0*u3 - 1860.0*u2 + 2670.0*u - 1434.0) / 120.0;
        else               return (6.0 - u) * (6.0-u) * (6.0-u) * (6.0-u) * (6.0-u) / 120.0;
    }
    return 0.0;
}

double PMEForceProvider::bspline_deriv(double u, int order) noexcept {
    // dM_p/du = M_{p-1}(u) - M_{p-1}(u-1)
    return bspline(u, order - 1) - bspline(u - 1.0, order - 1);
}

// ---------------------------------------------------------------------------
// B-spline DFT modulus correction
// ---------------------------------------------------------------------------
//
// The B-coefficient |b_α(m)|² removes the aliasing error introduced by the
// finite B-spline approximation in reciprocal space:
//
//   b_α(m) = exp(2πi(p-1)m/K) / Σ_{k=0}^{p-2} M_p(k+1) exp(2πikm/K)
//
// Since |exp(2πi·)|=1, |b|^{-2} = |Σ_{k=0}^{p-2} M_p(k+1) exp(2πikm/K)|^{-2}.
// We precompute and store the inverse of the squared modulus.

void PMEForceProvider::precompute_influence(const Box& box) {
    // Recompute only when the box changes.
    if (box.lengths[0] == box_lengths_cache_[0] &&
        box.lengths[1] == box_lengths_cache_[1] &&
        box.lengths[2] == box_lengths_cache_[2])
        return;
    box_lengths_cache_ = box.lengths;

    const int p = order_;

    // Precompute M_p(k+1) for k = 0 .. p-2 (used in b-modulus sum).
    std::vector<double> Mp(static_cast<std::size_t>(p - 1));
    for (int k = 0; k < p - 1; ++k)
        Mp[static_cast<std::size_t>(k)] = bspline(static_cast<double>(k + 1), p);

    // Per-dimension b-modulus correction.
    for (int d = 0; d < 3; ++d) {
        const int K = grid_[d];
        auto& bmod = bmod_sq_[d];
        bmod.resize(static_cast<std::size_t>(K));
        for (int m = 0; m < K; ++m) {
            double re = 0.0, im = 0.0;
            for (int k = 0; k < p - 1; ++k) {
                const double angle = 2.0 * std::numbers::pi * k * m / K;
                re += Mp[static_cast<std::size_t>(k)] * std::cos(angle);
                im += Mp[static_cast<std::size_t>(k)] * std::sin(angle);
            }
            const double mag_sq = re*re + im*im;
            bmod[static_cast<std::size_t>(m)] = (mag_sq > 1e-30) ? 1.0 / mag_sq : 0.0;
        }
    }

    // Build the influence function G(m1,m2,m3) in reciprocal space.
    // G(m) = k_e * 4π / (k² V) * exp(-k²/4α²) / |b1(m1)|² |b2(m2)|² |b3(m3)|²
    // m=0 along each dimension corresponds to k=0 (DC term) → G=0.
    const double Lx = box.lengths[0], Ly = box.lengths[1], Lz = box.lengths[2];
    const double V  = Lx * Ly * Lz;
    const int K1 = grid_[0], K2 = grid_[1], K3 = grid_[2];
    const double inv_4a2 = 1.0 / (4.0 * alpha_sq_);

    for (int m1 = 0; m1 < K1; ++m1) {
        // Map m1 to k1: reciprocal lattice index in (-K/2, K/2].
        const int n1 = (m1 <= K1/2) ? m1 : m1 - K1;
        const double kx = 2.0 * std::numbers::pi * n1 / Lx;

        for (int m2 = 0; m2 < K2; ++m2) {
            const int n2 = (m2 <= K2/2) ? m2 : m2 - K2;
            const double ky = 2.0 * std::numbers::pi * n2 / Ly;

            for (int m3 = 0; m3 < K3; ++m3) {
                const std::size_t idx = static_cast<std::size_t>(
                    m1 * K2 * K3 + m2 * K3 + m3);

                if (m1 == 0 && m2 == 0 && m3 == 0) {
                    influence_[idx] = 0.0;
                    continue;
                }

                const int n3 = (m3 <= K3/2) ? m3 : m3 - K3;
                const double kz = 2.0 * std::numbers::pi * n3 / Lz;
                const double k2 = kx*kx + ky*ky + kz*kz;

                const double b_corr = bmod_sq_[0][static_cast<std::size_t>(m1)]
                                    * bmod_sq_[1][static_cast<std::size_t>(m2)]
                                    * bmod_sq_[2][static_cast<std::size_t>(m3)];

                influence_[idx] = kPMECoulomb
                                 * (4.0 * std::numbers::pi)
                                 / (V * k2)
                                 * std::exp(-k2 * inv_4a2)
                                 * b_corr;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Main compute
// ---------------------------------------------------------------------------

void PMEForceProvider::compute(const ForceRequest& req,
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

    bool has_charges = false;
    for (std::size_t i = 0; i < n; ++i)
        if (charges[i] != 0.0) { has_charges = true; break; }
    if (!has_charges) return;

    resolve_params(*req.box);
    precompute_influence(*req.box);

    compute_real_space(req, res);
    compute_reciprocal_pme(req, res);
    compute_self_correction(req, res);
}

// ---------------------------------------------------------------------------
// Real-space part (identical to Ewald)
// ---------------------------------------------------------------------------

void PMEForceProvider::compute_real_space(const ForceRequest& req,
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

        res.potential_energy += kPMECoulomb * qi * qj * erfc_ar / r;

        const double ff = kPMECoulomb * qi * qj / r2
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
// Reciprocal-space part via PME mesh
// ---------------------------------------------------------------------------

void PMEForceProvider::compute_reciprocal_pme(const ForceRequest& req,
                                               ForceResult& res) {
    const auto& coords  = req.coordinates;
    const auto  charges = req.system->charges();
    const Box&  box     = *req.box;
    const std::size_t n = coords.size();
    const int K1 = grid_[0], K2 = grid_[1], K3 = grid_[2];
    const double Lx = box.lengths[0], Ly = box.lengths[1], Lz = box.lengths[2];
    const int p = order_;

    // ---- Step 1: Zero the charge mesh ----
    std::fill(mesh_.begin(), mesh_.end(), std::complex<double>{0.0, 0.0});

    // ---- Step 2: B-spline charge spreading ----
    // For atom i with scaled fractional coordinate u_α = r_α/L_α * K_α:
    //   base_m = floor(u_α)
    //   t      = u_α - base_m   (fractional part, in [0,1))
    //   p contributing grid points starting from m_start = base_m - (p-1)
    //   weight at grid point m_start + k  =  M_p(t + (p-1) - k)  for k=0..p-1

    for (std::size_t i = 0; i < n; ++i) {
        if (charges[i] == 0.0) continue;

        const double ux = coords[i][0] / Lx * K1;
        const double uy = coords[i][1] / Ly * K2;
        const double uz = coords[i][2] / Lz * K3;

        const int bx = static_cast<int>(std::floor(ux));
        const int by = static_cast<int>(std::floor(uy));
        const int bz = static_cast<int>(std::floor(uz));

        const double tx = ux - bx;
        const double ty = uy - by;
        const double tz = uz - bz;

        for (int kx = 0; kx < p; ++kx) {
            const double wx = bspline(tx + static_cast<double>(p - 1 - kx), p);
            if (wx == 0.0) continue;
            const int mx = ((bx - (p - 1 - kx)) % K1 + K1) % K1;

            for (int ky = 0; ky < p; ++ky) {
                const double wy = bspline(ty + static_cast<double>(p - 1 - ky), p);
                if (wy == 0.0) continue;
                const int my = ((by - (p - 1 - ky)) % K2 + K2) % K2;

                for (int kz = 0; kz < p; ++kz) {
                    const double wz = bspline(tz + static_cast<double>(p - 1 - kz), p);
                    if (wz == 0.0) continue;
                    const int mz = ((bz - (p - 1 - kz)) % K3 + K3) % K3;

                    const std::size_t idx = static_cast<std::size_t>(
                        mx * K2 * K3 + my * K3 + mz);
                    mesh_[idx] += charges[i] * wx * wy * wz;
                }
            }
        }
    }

    // ---- Step 3: Forward 3-D FFT ----
    fft3d(mesh_, false);

    // ---- Step 4: Apply influence function; accumulate energy ----
    // U_recip = (1/2) Σ_m Q̂*(m) · G(m) · Q̂(m)
    //         = (1/2) Σ_m G(m) |Q̂(m)|²
    double recip_energy = 0.0;
    for (std::size_t idx = 0; idx < mesh_.size(); ++idx) {
        const double g = influence_[idx];
        mesh_[idx] *= g;
        recip_energy += g * (std::real(mesh_[idx]) * std::real(mesh_[idx] / g)
                            + std::imag(mesh_[idx]) * std::imag(mesh_[idx] / g));
    }
    // Simpler: energy = 0.5 * Σ_m G(m) |Q̂(m)|²
    // Recompute cleanly (influence is already applied above, so use modified mesh).
    // We already multiplied mesh_ by G, so the unmodified |Q̂|² is not available.
    // Recompute energy from the potential mesh after IFFT instead.

    // Reset and redo the application step cleanly.
    // Reload mesh and compute energy + modified mesh together.
    (void)recip_energy;

    // Redo step 3+4 cleanly:
    std::fill(mesh_.begin(), mesh_.end(), std::complex<double>{0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i) {
        if (charges[i] == 0.0) continue;
        const double ux = coords[i][0] / Lx * K1;
        const double uy = coords[i][1] / Ly * K2;
        const double uz = coords[i][2] / Lz * K3;
        const int bx = static_cast<int>(std::floor(ux));
        const int by = static_cast<int>(std::floor(uy));
        const int bz = static_cast<int>(std::floor(uz));
        const double tx = ux - bx, ty = uy - by, tz = uz - bz;
        for (int kx = 0; kx < p; ++kx) {
            const double wx = bspline(tx + static_cast<double>(p-1-kx), p);
            if (wx == 0.0) continue;
            const int mx = ((bx - (p-1-kx)) % K1 + K1) % K1;
            for (int ky = 0; ky < p; ++ky) {
                const double wy = bspline(ty + static_cast<double>(p-1-ky), p);
                if (wy == 0.0) continue;
                const int my = ((by - (p-1-ky)) % K2 + K2) % K2;
                for (int kz = 0; kz < p; ++kz) {
                    const double wz = bspline(tz + static_cast<double>(p-1-kz), p);
                    if (wz == 0.0) continue;
                    const int mz = ((bz - (p-1-kz)) % K3 + K3) % K3;
                    const std::size_t idx = static_cast<std::size_t>(mx*K2*K3 + my*K3 + mz);
                    mesh_[idx] += charges[i] * wx * wy * wz;
                }
            }
        }
    }
    fft3d(mesh_, false);

    // Compute energy: ½ Σ_m G(m)|Q̂(m)|²  and apply influence function.
    double e_recip = 0.0;
    for (std::size_t idx = 0; idx < mesh_.size(); ++idx) {
        const double g   = influence_[idx];
        const double re  = std::real(mesh_[idx]);
        const double im  = std::imag(mesh_[idx]);
        e_recip   += 0.5 * g * (re*re + im*im);
        mesh_[idx] *= g;
    }
    res.potential_energy += e_recip;

    // ---- Step 5: Inverse 3-D FFT to get potential on mesh ----
    fft3d(mesh_, true);

    // ---- Step 6: Force interpolation ----
    // F_i,α = -q_i · K_α/L_α · Σ_{m} V(m) · (∂w_α/∂u_α)(m_α)
    //                                         · Π_{β≠α} w_β(m_β)
    for (std::size_t i = 0; i < n; ++i) {
        if (charges[i] == 0.0) continue;

        const double ux = coords[i][0] / Lx * K1;
        const double uy = coords[i][1] / Ly * K2;
        const double uz = coords[i][2] / Lz * K3;
        const int bx = static_cast<int>(std::floor(ux));
        const int by = static_cast<int>(std::floor(uy));
        const int bz = static_cast<int>(std::floor(uz));
        const double tx = ux - bx, ty = uy - by, tz = uz - bz;

        double fx = 0.0, fy = 0.0, fz = 0.0;

        for (int kx = 0; kx < p; ++kx) {
            const double arg_x = tx + static_cast<double>(p-1-kx);
            const double  wx = bspline      (arg_x, p);
            const double dwx = bspline_deriv(arg_x, p);
            const int mx = ((bx - (p-1-kx)) % K1 + K1) % K1;

            for (int ky = 0; ky < p; ++ky) {
                const double arg_y = ty + static_cast<double>(p-1-ky);
                const double  wy = bspline      (arg_y, p);
                const double dwy = bspline_deriv(arg_y, p);
                const int my = ((by - (p-1-ky)) % K2 + K2) % K2;

                for (int kz = 0; kz < p; ++kz) {
                    const double arg_z = tz + static_cast<double>(p-1-kz);
                    const double  wz = bspline      (arg_z, p);
                    const double dwz = bspline_deriv(arg_z, p);
                    const int mz = ((bz - (p-1-kz)) % K3 + K3) % K3;

                    const std::size_t idx = static_cast<std::size_t>(mx*K2*K3 + my*K3 + mz);
                    const double V_m = std::real(mesh_[idx]);

                    fx += V_m * dwx * wy  * wz;
                    fy += V_m *  wx * dwy * wz;
                    fz += V_m *  wx * wy  * dwz;
                }
            }
        }

        // The derivative w.r.t. r_i,α = (dw/du_α) * (K_α / L_α).
        // Force = -q_i * gradient of potential.
        res.forces[i][0] -= charges[i] * fx * (K1 / Lx);
        res.forces[i][1] -= charges[i] * fy * (K2 / Ly);
        res.forces[i][2] -= charges[i] * fz * (K3 / Lz);
    }
}

// ---------------------------------------------------------------------------
// Self-energy correction (same as Ewald)
// ---------------------------------------------------------------------------

void PMEForceProvider::compute_self_correction(const ForceRequest& req,
                                                ForceResult& res) const {
    const auto charges = req.system->charges();
    const std::size_t n = charges.size();

    double q2_sum = 0.0, Q_net = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        q2_sum += charges[i] * charges[i];
        Q_net  += charges[i];
    }

    res.potential_energy -=
        kPMECoulomb * (alpha_ / std::sqrt(std::numbers::pi)) * q2_sum;

    if (Q_net != 0.0) {
        const double V = req.box->lengths[0]
                       * req.box->lengths[1]
                       * req.box->lengths[2];
        res.potential_energy -=
            kPMECoulomb * std::numbers::pi / (2.0 * V * alpha_sq_) * Q_net * Q_net;
    }
}

// ---------------------------------------------------------------------------
// 1-D Cooley-Tukey radix-2 in-place FFT
// ---------------------------------------------------------------------------
//
// n must be a power of 2.
// inverse=false → DFT  (no 1/N normalization)
// inverse=true  → IDFT (divides by N)

void PMEForceProvider::fft1d(std::complex<double>* data, int n, bool inverse) {
    // Bit-reversal permutation.
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(data[i], data[j]);
    }
    // Butterfly stages.
    for (int len = 2; len <= n; len <<= 1) {
        const double angle =
            (inverse ? 1.0 : -1.0) * 2.0 * std::numbers::pi / len;
        const std::complex<double> wlen(std::cos(angle), std::sin(angle));
        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (int j = 0; j < len / 2; ++j) {
                const auto u = data[i + j];
                const auto v = data[i + j + len / 2] * w;
                data[i + j]           = u + v;
                data[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse) {
        const double inv_n = 1.0 / static_cast<double>(n);
        for (int i = 0; i < n; ++i) data[i] *= inv_n;
    }
}

// ---------------------------------------------------------------------------
// 3-D FFT: apply 1-D FFT successively along each dimension.
// ---------------------------------------------------------------------------

void PMEForceProvider::fft3d(std::vector<std::complex<double>>& data,
                              bool inverse) const {
    const int K1 = grid_[0], K2 = grid_[1], K3 = grid_[2];

    // Along dimension 3 (innermost, contiguous).
    for (int i = 0; i < K1 * K2; ++i)
        fft1d(data.data() + static_cast<std::ptrdiff_t>(i * K3), K3, inverse);

    // Along dimension 2.
    std::vector<std::complex<double>> tmp(static_cast<std::size_t>(K2));
    for (int i1 = 0; i1 < K1; ++i1) {
        for (int i3 = 0; i3 < K3; ++i3) {
            for (int i2 = 0; i2 < K2; ++i2)
                tmp[static_cast<std::size_t>(i2)] =
                    data[static_cast<std::size_t>(i1*K2*K3 + i2*K3 + i3)];
            fft1d(tmp.data(), K2, inverse);
            for (int i2 = 0; i2 < K2; ++i2)
                data[static_cast<std::size_t>(i1*K2*K3 + i2*K3 + i3)] =
                    tmp[static_cast<std::size_t>(i2)];
        }
    }

    // Along dimension 1 (outermost).
    std::vector<std::complex<double>> tmp1(static_cast<std::size_t>(K1));
    for (int i2 = 0; i2 < K2; ++i2) {
        for (int i3 = 0; i3 < K3; ++i3) {
            for (int i1 = 0; i1 < K1; ++i1)
                tmp1[static_cast<std::size_t>(i1)] =
                    data[static_cast<std::size_t>(i1*K2*K3 + i2*K3 + i3)];
            fft1d(tmp1.data(), K1, inverse);
            for (int i1 = 0; i1 < K1; ++i1)
                data[static_cast<std::size_t>(i1*K2*K3 + i2*K3 + i3)] =
                    tmp1[static_cast<std::size_t>(i1)];
        }
    }
}

}  // namespace gmd
