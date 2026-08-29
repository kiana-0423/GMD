// Test-only reference implementation of the PME reciprocal machinery.
//
// Nothing here is used by the engine. Every piece is written from its own
// definition so that comparing it against the production path is a second
// opinion rather than a restatement:
//
//   * the cardinal B-spline comes from the order recursion, not from a
//     polynomial table, so a transcription error in the engine's hard-coded
//     polynomials shows up as a disagreement;
//   * the transforms are direct DFTs, not Cooley-Tukey, and both directions
//     are UNNORMALISED. The normalisation convention is the whole point of
//     pinning the FFT path down, so it is stated explicitly here rather than
//     inherited from an FFT routine:
//
//         forward:  Qhat(m)  = sum_j Q(j)    exp(-2 pi i j.m / K)
//         inverse:  theta(j) = sum_m Ghat(m) exp(+2 pi i j.m / K)
//
//     Note that the inverse carries NO 1/N. The engine's fft3d(..., true) does
//     divide by N, which is exactly why its force interpolation has to
//     multiply by K1*K2*K3 to recover the gradient. A reference that hid the
//     same convention inside the same helper could not detect a missing or
//     doubled factor.
//
// Two independent expressions for the reciprocal energy are provided, related
// by Parseval:
//
//     E = 1/2 sum_m G(m) |Qhat(m)|^2        (reciprocal space)
//     E = 1/2 sum_j Q(j) theta(j)           (real space, with the convention above)
//
// They share no arithmetic beyond Qhat, so agreement between them and the
// production energy is a genuine cross-check of the transform convention.

#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <vector>

namespace pme_ref {

using Real = long double;

// CODATA 2022 k_e = E_h * a0, matching gmd::kCoulombConstant.
inline constexpr Real kCoulomb = 14.3996454686836L;
// Declared here rather than imported from
// gmd/core/physical_constants.hpp on purpose: this header is an
// independent reference, and sharing the production symbol would make
// every comparison built on it self-confirming. It must track the
// production value, and tests/electrostatic_constant_tests.cpp is what
// fails if it stops doing so.

// --- cardinal B-spline, from the recursion --------------------------------
//
//   M_2(u) = 1 - |u - 1|                       on (0, 2)
//   M_p(u) = [ u M_(p-1)(u) + (p - u) M_(p-1)(u - 1) ] / (p - 1)
//
// The engine evaluates hard-coded piecewise polynomials instead; this is the
// definition they are supposed to expand to.
inline Real cardinal_bspline(Real u, int order) {
    if (u <= 0.0L || u >= static_cast<Real>(order)) return 0.0L;
    if (order == 2) return 1.0L - std::fabs(u - 1.0L);
    return (u * cardinal_bspline(u, order - 1)
            + (static_cast<Real>(order) - u) * cardinal_bspline(u - 1.0L, order - 1))
           / static_cast<Real>(order - 1);
}

inline Real cardinal_bspline_derivative(Real u, int order) {
    return cardinal_bspline(u, order - 1) - cardinal_bspline(u - 1.0L, order - 1);
}

// |b_d(m)|^2 = |sum_{k=0}^{p-2} M_p(k+1) exp(2 pi i k m / K)|^-2 (Essmann 1995).
inline std::vector<Real> bspline_modulus_squared(int grid, int order) {
    std::vector<Real> moduli(static_cast<std::size_t>(grid), 0.0L);
    for (int m = 0; m < grid; ++m) {
        Real re = 0.0L;
        Real im = 0.0L;
        for (int k = 0; k <= order - 2; ++k) {
            const Real weight = cardinal_bspline(static_cast<Real>(k + 1), order);
            const Real angle = 2.0L * std::numbers::pi_v<Real> * k * m / grid;
            re += weight * std::cos(angle);
            im += weight * std::sin(angle);
        }
        const Real magnitude = re * re + im * im;
        moduli[static_cast<std::size_t>(m)] =
            magnitude > 1.0e-30L ? 1.0L / magnitude : 0.0L;
    }
    return moduli;
}

// --- direct transforms ----------------------------------------------------

// Separable direct DFT. `sign` is -1 for the forward transform and +1 for the
// inverse; neither direction is normalised (see the header comment).
inline std::vector<std::complex<Real>> direct_dft3(
        const std::vector<std::complex<Real>>& input,
        const std::array<int, 3>& grid,
        int sign) {
    const int k1 = grid[0], k2 = grid[1], k3 = grid[2];
    std::vector<std::complex<Real>> data = input;

    auto transform_axis = [&](int length, int outer1, int outer2, auto index_of) {
        std::vector<std::complex<Real>> in(static_cast<std::size_t>(length));
        std::vector<std::complex<Real>> out(static_cast<std::size_t>(length));
        for (int p = 0; p < outer1; ++p) {
            for (int q = 0; q < outer2; ++q) {
                for (int j = 0; j < length; ++j) {
                    in[static_cast<std::size_t>(j)] = data[index_of(p, q, j)];
                }
                for (int m = 0; m < length; ++m) {
                    std::complex<Real> sum(0.0L, 0.0L);
                    for (int j = 0; j < length; ++j) {
                        const Real angle = static_cast<Real>(sign) * 2.0L
                                         * std::numbers::pi_v<Real> * j * m / length;
                        sum += in[static_cast<std::size_t>(j)]
                             * std::complex<Real>(std::cos(angle), std::sin(angle));
                    }
                    out[static_cast<std::size_t>(m)] = sum;
                }
                for (int m = 0; m < length; ++m) {
                    data[index_of(p, q, m)] = out[static_cast<std::size_t>(m)];
                }
            }
        }
    };

    transform_axis(k3, k1, k2, [&](int p, int q, int j) {
        return static_cast<std::size_t>(p * k2 * k3 + q * k3 + j);
    });
    transform_axis(k2, k1, k3, [&](int p, int q, int j) {
        return static_cast<std::size_t>(p * k2 * k3 + j * k3 + q);
    });
    transform_axis(k1, k2, k3, [&](int p, int q, int j) {
        return static_cast<std::size_t>(j * k2 * k3 + p * k3 + q);
    });
    return data;
}

// --- charge mesh ----------------------------------------------------------

// The B-spline charge mesh Q(m1,m2,m3) from fractional coordinates. The index
// arithmetic mirrors the engine's (base = floor(u), weights running from
// M_p(t + p - 1) down), because the mesh is the shared definition, not the
// thing under test here.
inline std::vector<Real> charge_mesh(const std::vector<std::array<Real, 3>>& fractional,
                                     const std::vector<Real>& charges,
                                     const std::array<int, 3>& grid,
                                     int order) {
    const std::size_t total = static_cast<std::size_t>(grid[0])
                            * static_cast<std::size_t>(grid[1])
                            * static_cast<std::size_t>(grid[2]);
    std::vector<Real> mesh(total, 0.0L);

    for (std::size_t i = 0; i < fractional.size(); ++i) {
        const Real q = charges[i];
        if (q == 0.0L) continue;

        std::array<int, 3> base{};
        std::array<Real, 3> frac{};
        for (std::size_t d = 0; d < 3; ++d) {
            const Real u = fractional[i][d] * static_cast<Real>(grid[d]);
            base[d] = static_cast<int>(std::floor(static_cast<double>(u)));
            frac[d] = u - static_cast<Real>(base[d]);
        }

        for (int a = 0; a < order; ++a) {
            const Real wx = cardinal_bspline(frac[0] + static_cast<Real>(order - 1 - a), order);
            if (wx == 0.0L) continue;
            const int mx = ((base[0] - (order - 1 - a)) % grid[0] + grid[0]) % grid[0];
            for (int b = 0; b < order; ++b) {
                const Real wy = cardinal_bspline(frac[1] + static_cast<Real>(order - 1 - b), order);
                if (wy == 0.0L) continue;
                const int my = ((base[1] - (order - 1 - b)) % grid[1] + grid[1]) % grid[1];
                for (int c = 0; c < order; ++c) {
                    const Real wz =
                        cardinal_bspline(frac[2] + static_cast<Real>(order - 1 - c), order);
                    if (wz == 0.0L) continue;
                    const int mz = ((base[2] - (order - 1 - c)) % grid[2] + grid[2]) % grid[2];
                    const auto index = static_cast<std::size_t>(
                        mx * grid[1] * grid[2] + my * grid[2] + mz);
                    mesh[index] += q * wx * wy * wz;
                }
            }
        }
    }
    return mesh;
}

// --- influence function ---------------------------------------------------

// G(m) = k_e 4 pi / (V k^2) exp(-k^2 / 4 alpha^2) |b1|^2 |b2|^2 |b3|^2,
// for an orthorhombic cell of the given edge lengths. G(0) = 0.
inline std::vector<Real> influence(const std::array<Real, 3>& lengths,
                                   Real alpha,
                                   const std::array<int, 3>& grid,
                                   int order) {
    const std::array<std::vector<Real>, 3> moduli = {
        bspline_modulus_squared(grid[0], order),
        bspline_modulus_squared(grid[1], order),
        bspline_modulus_squared(grid[2], order)};

    const Real volume = lengths[0] * lengths[1] * lengths[2];
    const Real inv_four_alpha_sq = 1.0L / (4.0L * alpha * alpha);

    const auto total = static_cast<std::size_t>(grid[0]) * static_cast<std::size_t>(grid[1])
                     * static_cast<std::size_t>(grid[2]);
    std::vector<Real> g(total, 0.0L);

    for (int m1 = 0; m1 < grid[0]; ++m1) {
        const int n1 = (m1 <= grid[0] / 2) ? m1 : m1 - grid[0];
        const Real kx = 2.0L * std::numbers::pi_v<Real> * n1 / lengths[0];
        for (int m2 = 0; m2 < grid[1]; ++m2) {
            const int n2 = (m2 <= grid[1] / 2) ? m2 : m2 - grid[1];
            const Real ky = 2.0L * std::numbers::pi_v<Real> * n2 / lengths[1];
            for (int m3 = 0; m3 < grid[2]; ++m3) {
                const auto index = static_cast<std::size_t>(
                    m1 * grid[1] * grid[2] + m2 * grid[2] + m3);
                if (m1 == 0 && m2 == 0 && m3 == 0) continue;
                const int n3 = (m3 <= grid[2] / 2) ? m3 : m3 - grid[2];
                const Real kz = 2.0L * std::numbers::pi_v<Real> * n3 / lengths[2];
                const Real k_sq = kx * kx + ky * ky + kz * kz;
                g[index] = kCoulomb * 4.0L * std::numbers::pi_v<Real> / (volume * k_sq)
                         * std::exp(-k_sq * inv_four_alpha_sq)
                         * moduli[0][static_cast<std::size_t>(m1)]
                         * moduli[1][static_cast<std::size_t>(m2)]
                         * moduli[2][static_cast<std::size_t>(m3)];
            }
        }
    }
    return g;
}

// --- assembled reference --------------------------------------------------

struct ReciprocalReference {
    Real energy_reciprocal_space = 0.0L;  // 1/2 sum_m G |Qhat|^2
    Real energy_real_space = 0.0L;        // 1/2 sum_j Q(j) theta(j)
    std::vector<std::array<Real, 3>> forces;
};

// Evaluates the PME reciprocal energy both ways and the reciprocal force, all
// from the unnormalised transforms above.
//
// The force is the analytic gradient of E with respect to the atomic position,
// taken through the spline weights:
//
//   F_i,a = -q_i (K_a / L_a) sum_m theta(m) (dw_a/du_a)(m_a) prod_{b != a} w_b(m_b)
//
// with theta the UNNORMALISED inverse transform of G*Qhat. Written this way,
// the mesh-point count never appears as a separate factor: it is absorbed into
// the convention, which is precisely what makes this able to detect the
// engine either omitting or double-applying it.
inline ReciprocalReference reciprocal_reference(
        const std::vector<std::array<Real, 3>>& positions,
        const std::vector<Real>& charges,
        const std::array<Real, 3>& lengths,
        Real alpha,
        const std::array<int, 3>& grid,
        int order) {
    std::vector<std::array<Real, 3>> fractional(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            fractional[i][d] = positions[i][d] / lengths[d];
        }
    }

    const std::vector<Real> mesh = charge_mesh(fractional, charges, grid, order);
    const std::vector<Real> g = influence(lengths, alpha, grid, order);

    std::vector<std::complex<Real>> spectral(mesh.size());
    for (std::size_t i = 0; i < mesh.size(); ++i) {
        spectral[i] = std::complex<Real>(mesh[i], 0.0L);
    }
    const auto transformed = direct_dft3(spectral, grid, -1);

    ReciprocalReference reference;

    std::vector<std::complex<Real>> filtered(transformed.size());
    for (std::size_t i = 0; i < transformed.size(); ++i) {
        filtered[i] = transformed[i] * g[i];
        const Real re = transformed[i].real();
        const Real im = transformed[i].imag();
        reference.energy_reciprocal_space += 0.5L * g[i] * (re * re + im * im);
    }

    // theta = unnormalised inverse transform of G * Qhat.
    const auto potential = direct_dft3(filtered, grid, +1);
    for (std::size_t i = 0; i < mesh.size(); ++i) {
        reference.energy_real_space += 0.5L * mesh[i] * potential[i].real();
    }

    reference.forces.assign(positions.size(), std::array<Real, 3>{0.0L, 0.0L, 0.0L});
    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (charges[i] == 0.0L) continue;

        std::array<int, 3> base{};
        std::array<Real, 3> frac{};
        for (std::size_t d = 0; d < 3; ++d) {
            const Real u = fractional[i][d] * static_cast<Real>(grid[d]);
            base[d] = static_cast<int>(std::floor(static_cast<double>(u)));
            frac[d] = u - static_cast<Real>(base[d]);
        }

        std::array<Real, 3> gradient = {0.0L, 0.0L, 0.0L};
        for (int a = 0; a < order; ++a) {
            const Real ax = frac[0] + static_cast<Real>(order - 1 - a);
            const Real wx = cardinal_bspline(ax, order);
            const Real dx = cardinal_bspline_derivative(ax, order);
            const int mx = ((base[0] - (order - 1 - a)) % grid[0] + grid[0]) % grid[0];
            for (int b = 0; b < order; ++b) {
                const Real ay = frac[1] + static_cast<Real>(order - 1 - b);
                const Real wy = cardinal_bspline(ay, order);
                const Real dy = cardinal_bspline_derivative(ay, order);
                const int my = ((base[1] - (order - 1 - b)) % grid[1] + grid[1]) % grid[1];
                for (int c = 0; c < order; ++c) {
                    const Real az = frac[2] + static_cast<Real>(order - 1 - c);
                    const Real wz = cardinal_bspline(az, order);
                    const Real dz = cardinal_bspline_derivative(az, order);
                    const int mz = ((base[2] - (order - 1 - c)) % grid[2] + grid[2]) % grid[2];

                    const auto index = static_cast<std::size_t>(
                        mx * grid[1] * grid[2] + my * grid[2] + mz);
                    const Real value = potential[index].real();
                    gradient[0] += value * dx * wy * wz;
                    gradient[1] += value * wx * dy * wz;
                    gradient[2] += value * wx * wy * dz;
                }
            }
        }
        for (std::size_t d = 0; d < 3; ++d) {
            reference.forces[i][d] =
                -charges[i] * gradient[d] * static_cast<Real>(grid[d]) / lengths[d];
        }
    }

    return reference;
}

// U_self = -k_e (alpha / sqrt(pi)) sum_i q_i^2. No cell or position
// dependence, so it shifts the energy and contributes no force.
inline Real self_energy(const std::vector<Real>& charges, Real alpha) {
    Real sum = 0.0L;
    for (const Real q : charges) sum += q * q;
    return -kCoulomb * (alpha / std::sqrt(std::numbers::pi_v<Real>)) * sum;
}

// U_net = -k_e pi / (2 V alpha^2) Q^2 for a non-neutral cell.
inline Real net_charge_energy(const std::vector<Real>& charges, Real volume, Real alpha) {
    Real total = 0.0L;
    for (const Real q : charges) total += q;
    if (total == 0.0L) return 0.0L;
    return -kCoulomb * std::numbers::pi_v<Real> / (2.0L * volume * alpha * alpha)
           * total * total;
}

}  // namespace pme_ref
