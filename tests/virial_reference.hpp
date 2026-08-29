// Test-only independent references for virial validation.
//
// Nothing in this header is used by the engine. Every quantity is recomputed
// from its own definition so that a comparison against a production tensor is
// a genuine second opinion rather than a restatement of the same code.
//
// WHY A TRICLINIC REFERENCE EXISTS HERE
//
// `gmd::Box` stores three edge lengths, so the engine can represent
// orthorhombic cells only and cannot apply a shear strain. That is a limit of
// the *engine*, not of a reference implementation. The Ewald and PME
// reciprocal energies are written here against a general 3x3 cell matrix h,
// which can be sheared freely. Differentiating that reference energy with
// respect to a general strain therefore yields all nine virial components,
// including the off-diagonal ones the engine's own finite-difference test
// cannot reach. Comparing the engine's analytic tensor -- evaluated at an
// orthorhombic h, which is all it supports -- against that derivative is a
// real validation of the off-diagonal components.
//
// The strain convention throughout is
//
//     h -> (I + e) h,   r_i -> (I + e) r_i   (fractional coordinates fixed)
//     W_ab = -dE/de_ab |_(e=0)
//
// which is the convention that makes W_ab = sum_i r_ia F_ib for a pair
// potential, and P = (2K + tr W) / (3V).
//
// `long double` is used wherever it is free. Note that it is only 53-bit on
// Apple arm64 and 64-bit on x86-64 Linux, so every tolerance in the tests is
// calibrated for the 53-bit case.

#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <numbers>
#include <vector>

#include "pme_reference.hpp"

namespace virial_ref {

using Real = long double;
using Vec3 = std::array<Real, 3>;
// Row-major 3x3: M[a * 3 + b] is M_ab.
using Mat3 = std::array<Real, 9>;

inline constexpr Real kCoulomb = 14.3996L;  // matches kEwaldCoulomb / kPMECoulomb

// --- small linear algebra -------------------------------------------------

inline Mat3 identity3() noexcept {
    return {1.0L, 0.0L, 0.0L, 0.0L, 1.0L, 0.0L, 0.0L, 0.0L, 1.0L};
}

inline Mat3 mat_mul(const Mat3& a, const Mat3& b) noexcept {
    Mat3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            Real sum = 0.0L;
            for (std::size_t k = 0; k < 3; ++k) sum += a[i * 3 + k] * b[k * 3 + j];
            out[i * 3 + j] = sum;
        }
    }
    return out;
}

inline Vec3 mat_vec(const Mat3& m, const Vec3& v) noexcept {
    return {m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
            m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
            m[6] * v[0] + m[7] * v[1] + m[8] * v[2]};
}

inline Mat3 transpose3(const Mat3& m) noexcept {
    return {m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]};
}

inline Real det3(const Mat3& m) noexcept {
    return m[0] * (m[4] * m[8] - m[5] * m[7])
         - m[1] * (m[3] * m[8] - m[5] * m[6])
         + m[2] * (m[3] * m[7] - m[4] * m[6]);
}

inline Mat3 inverse3(const Mat3& m) noexcept {
    const Real d = det3(m);
    Mat3 inv{};
    inv[0] = (m[4] * m[8] - m[5] * m[7]) / d;
    inv[1] = (m[2] * m[7] - m[1] * m[8]) / d;
    inv[2] = (m[1] * m[5] - m[2] * m[4]) / d;
    inv[3] = (m[5] * m[6] - m[3] * m[8]) / d;
    inv[4] = (m[0] * m[8] - m[2] * m[6]) / d;
    inv[5] = (m[2] * m[3] - m[0] * m[5]) / d;
    inv[6] = (m[3] * m[7] - m[4] * m[6]) / d;
    inv[7] = (m[1] * m[6] - m[0] * m[7]) / d;
    inv[8] = (m[0] * m[4] - m[1] * m[3]) / d;
    return inv;
}

inline Vec3 sub(const Vec3& a, const Vec3& b) noexcept {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline Real dot(const Vec3& a, const Vec3& b) noexcept {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline Real norm(const Vec3& v) noexcept { return std::sqrt(dot(v, v)); }

// Outer product accumulated as W_ab += u_a v_b.
inline void accumulate_outer(Mat3& w, const Vec3& u, const Vec3& v) noexcept {
    for (std::size_t a = 0; a < 3; ++a) {
        for (std::size_t b = 0; b < 3; ++b) {
            w[a * 3 + b] += u[a] * v[b];
        }
    }
}

// R W R^T, the covariance a rank-2 Cartesian tensor must obey under rotation.
inline Mat3 rotate_tensor(const Mat3& w, const Mat3& r) noexcept {
    return mat_mul(mat_mul(r, w), transpose3(r));
}

// Rotation by `angle` about a unit axis (Rodrigues).
inline Mat3 rotation_matrix(const Vec3& axis, Real angle) noexcept {
    const Real len = norm(axis);
    const Vec3 u = {axis[0] / len, axis[1] / len, axis[2] / len};
    const Real c = std::cos(angle);
    const Real s = std::sin(angle);
    const Real t = 1.0L - c;
    return {t * u[0] * u[0] + c,        t * u[0] * u[1] - s * u[2], t * u[0] * u[2] + s * u[1],
            t * u[0] * u[1] + s * u[2], t * u[1] * u[1] + c,        t * u[1] * u[2] - s * u[0],
            t * u[0] * u[2] - s * u[1], t * u[1] * u[2] + s * u[0], t * u[2] * u[2] + c};
}

// --- generic strain derivative -------------------------------------------

// W_ab = -dE/de_ab, central difference, Richardson-extrapolated from steps
// h and h/2 to suppress the leading O(h^2) truncation term.
//
// `energy` receives the deformation matrix (I + e) and must apply it to both
// the cell and the Cartesian coordinates, holding fractional coordinates
// fixed. Returning the energy of the *deformed* configuration is the whole
// contract; the caller decides what "the configuration" means.
inline Mat3 strain_derivative_virial(
        const std::function<Real(const Mat3& deformation)>& energy,
        Real step = 4.0e-5L) {
    Mat3 w{};
    for (std::size_t a = 0; a < 3; ++a) {
        for (std::size_t b = 0; b < 3; ++b) {
            auto derivative_at = [&](Real h) {
                Mat3 plus = identity3();
                Mat3 minus = identity3();
                plus[a * 3 + b] += h;
                minus[a * 3 + b] -= h;
                return (energy(plus) - energy(minus)) / (2.0L * h);
            };
            const Real coarse = derivative_at(step);
            const Real fine = derivative_at(step * 0.5L);
            // Richardson: (4 * fine - coarse) / 3 cancels the h^2 term.
            w[a * 3 + b] = -((4.0L * fine - coarse) / 3.0L);
        }
    }
    return w;
}

// --- Ewald reciprocal space, general cell ---------------------------------

// A configuration expressed the way a strain deformation wants it: a cell
// matrix whose COLUMNS are the three cell vectors, and fractional coordinates.
struct CellConfiguration {
    Mat3 cell = identity3();          // h; column c is cell vector c
    std::vector<Vec3> fractional;     // s_i, with r_i = h s_i
    std::vector<Real> charges;

    Vec3 cartesian(std::size_t i) const noexcept {
        return mat_vec(cell, fractional[i]);
    }
    Real volume() const noexcept { return det3(cell); }
    // Reciprocal lattice: k_n = 2 pi h^{-T} n.
    Mat3 reciprocal_basis() const noexcept {
        Mat3 inv_t = transpose3(inverse3(cell));
        for (auto& value : inv_t) value *= 2.0L * std::numbers::pi_v<Real>;
        return inv_t;
    }
};

// Builds a CellConfiguration for an orthorhombic cell from Cartesian
// coordinates -- the only cell shape the engine can produce.
inline CellConfiguration make_orthorhombic(const std::array<Real, 3>& lengths,
                                           const std::vector<Vec3>& coordinates,
                                           const std::vector<Real>& charges) {
    CellConfiguration configuration;
    configuration.cell = {lengths[0], 0.0L, 0.0L,
                          0.0L, lengths[1], 0.0L,
                          0.0L, 0.0L, lengths[2]};
    configuration.charges = charges;
    configuration.fractional.reserve(coordinates.size());
    for (const auto& r : coordinates) {
        configuration.fractional.push_back(
            {r[0] / lengths[0], r[1] / lengths[1], r[2] / lengths[2]});
    }
    return configuration;
}

// Applies (I + e) to the cell. Fractional coordinates are unchanged by
// construction, which is exactly the strain protocol.
inline CellConfiguration deform(const CellConfiguration& configuration,
                                const Mat3& deformation) {
    CellConfiguration out = configuration;
    out.cell = mat_mul(deformation, configuration.cell);
    return out;
}

// U_recip = (k_e / 2V) sum_{n != 0} (4 pi / k^2) exp(-k^2 / 4 alpha^2) |S(k)|^2
//
// The index range is the cube [-kmax, kmax]^3 minus the origin, matching the
// truncation EwaldForceProvider uses. A spherical cutoff would be a different
// sum and would not be comparable.
inline Real ewald_reciprocal_energy(const CellConfiguration& configuration,
                                    Real alpha, int kmax) {
    const Real volume = configuration.volume();
    const Mat3 basis = configuration.reciprocal_basis();
    const Real inv_four_alpha_sq = 1.0L / (4.0L * alpha * alpha);
    const std::size_t n = configuration.fractional.size();

    std::vector<Vec3> positions(n);
    for (std::size_t i = 0; i < n; ++i) positions[i] = configuration.cartesian(i);

    Real energy = 0.0L;
    for (int nx = -kmax; nx <= kmax; ++nx) {
        for (int ny = -kmax; ny <= kmax; ++ny) {
            for (int nz = -kmax; nz <= kmax; ++nz) {
                if (nx == 0 && ny == 0 && nz == 0) continue;
                const Vec3 index = {static_cast<Real>(nx), static_cast<Real>(ny),
                                    static_cast<Real>(nz)};
                const Vec3 k = mat_vec(basis, index);
                const Real k_sq = dot(k, k);

                Real structure_re = 0.0L;
                Real structure_im = 0.0L;
                for (std::size_t i = 0; i < n; ++i) {
                    const Real phase = dot(k, positions[i]);
                    structure_re += configuration.charges[i] * std::cos(phase);
                    structure_im += configuration.charges[i] * std::sin(phase);
                }

                const Real prefactor = kCoulomb * 4.0L * std::numbers::pi_v<Real>
                                     / (volume * k_sq)
                                     * std::exp(-k_sq * inv_four_alpha_sq);
                energy += 0.5L * prefactor
                        * (structure_re * structure_re + structure_im * structure_im);
            }
        }
    }
    return energy;
}

// F_i = -dU_recip/dr_i, derived independently of the provider:
//   d|S|^2/dr_i = 2 q_i k (S_im cos(k.r_i) - S_re sin(k.r_i))
inline std::vector<Vec3> ewald_reciprocal_forces(const CellConfiguration& configuration,
                                                 Real alpha, int kmax) {
    const Real volume = configuration.volume();
    const Mat3 basis = configuration.reciprocal_basis();
    const Real inv_four_alpha_sq = 1.0L / (4.0L * alpha * alpha);
    const std::size_t n = configuration.fractional.size();

    std::vector<Vec3> positions(n);
    for (std::size_t i = 0; i < n; ++i) positions[i] = configuration.cartesian(i);
    std::vector<Vec3> forces(n, Vec3{0.0L, 0.0L, 0.0L});

    for (int nx = -kmax; nx <= kmax; ++nx) {
        for (int ny = -kmax; ny <= kmax; ++ny) {
            for (int nz = -kmax; nz <= kmax; ++nz) {
                if (nx == 0 && ny == 0 && nz == 0) continue;
                const Vec3 index = {static_cast<Real>(nx), static_cast<Real>(ny),
                                    static_cast<Real>(nz)};
                const Vec3 k = mat_vec(basis, index);
                const Real k_sq = dot(k, k);

                Real structure_re = 0.0L;
                Real structure_im = 0.0L;
                for (std::size_t i = 0; i < n; ++i) {
                    const Real phase = dot(k, positions[i]);
                    structure_re += configuration.charges[i] * std::cos(phase);
                    structure_im += configuration.charges[i] * std::sin(phase);
                }

                const Real prefactor = kCoulomb * 4.0L * std::numbers::pi_v<Real>
                                     / (volume * k_sq)
                                     * std::exp(-k_sq * inv_four_alpha_sq);

                for (std::size_t i = 0; i < n; ++i) {
                    if (configuration.charges[i] == 0.0L) continue;
                    const Real phase = dot(k, positions[i]);
                    const Real scale = prefactor * configuration.charges[i]
                                     * (structure_re * std::sin(phase)
                                        - structure_im * std::cos(phase));
                    for (std::size_t d = 0; d < 3; ++d) forces[i][d] += scale * k[d];
                }
            }
        }
    }
    return forces;
}

// The analytic reciprocal virial, derived here from scratch:
//
//   E_k depends on the cell through 1/V and through k = 2 pi h^{-T} n.
//   dV/de_ab = V d_ab and dk_c/de_ab = -d_ca k_b, so dk^2/de_ab = -2 k_a k_b.
//   With E_k proportional to (1/V)(1/k^2) exp(-k^2/4a^2),
//
//     dE_k/de_ab = -d_ab E_k + 2 E_k (1/k^2 + 1/(4a^2)) k_a k_b
//     W_ab = -dE_k/de_ab = E_k [ d_ab - 2 (1/k^2 + 1/(4a^2)) k_a k_b ]
//
// This is an independent derivation of the same identity the provider uses;
// the finite-difference reference above is the check that does not share the
// derivation at all.
inline Mat3 ewald_reciprocal_virial_analytic(const CellConfiguration& configuration,
                                             Real alpha, int kmax) {
    const Real volume = configuration.volume();
    const Mat3 basis = configuration.reciprocal_basis();
    const Real inv_four_alpha_sq = 1.0L / (4.0L * alpha * alpha);
    const std::size_t n = configuration.fractional.size();

    std::vector<Vec3> positions(n);
    for (std::size_t i = 0; i < n; ++i) positions[i] = configuration.cartesian(i);

    Mat3 w{};
    for (int nx = -kmax; nx <= kmax; ++nx) {
        for (int ny = -kmax; ny <= kmax; ++ny) {
            for (int nz = -kmax; nz <= kmax; ++nz) {
                if (nx == 0 && ny == 0 && nz == 0) continue;
                const Vec3 index = {static_cast<Real>(nx), static_cast<Real>(ny),
                                    static_cast<Real>(nz)};
                const Vec3 k = mat_vec(basis, index);
                const Real k_sq = dot(k, k);

                Real structure_re = 0.0L;
                Real structure_im = 0.0L;
                for (std::size_t i = 0; i < n; ++i) {
                    const Real phase = dot(k, positions[i]);
                    structure_re += configuration.charges[i] * std::cos(phase);
                    structure_im += configuration.charges[i] * std::sin(phase);
                }

                const Real prefactor = kCoulomb * 4.0L * std::numbers::pi_v<Real>
                                     / (volume * k_sq)
                                     * std::exp(-k_sq * inv_four_alpha_sq);
                const Real e_k = 0.5L * prefactor
                               * (structure_re * structure_re
                                  + structure_im * structure_im);
                const Real coefficient = 2.0L * (inv_four_alpha_sq + 1.0L / k_sq);

                for (std::size_t a = 0; a < 3; ++a) {
                    for (std::size_t b = 0; b < 3; ++b) {
                        const Real delta = (a == b) ? 1.0L : 0.0L;
                        w[a * 3 + b] += e_k * (delta - coefficient * k[a] * k[b]);
                    }
                }
            }
        }
    }
    return w;
}

// U_self = -k_e (alpha / sqrt(pi)) sum_i q_i^2. No cell dependence at all, so
// its virial is identically zero -- that is the claim the tests check.
inline Real ewald_self_energy(const std::vector<Real>& charges, Real alpha) {
    Real sum = 0.0L;
    for (Real q : charges) sum += q * q;
    return -kCoulomb * (alpha / std::sqrt(std::numbers::pi_v<Real>)) * sum;
}

// U_net = -k_e pi / (2 V alpha^2) Q^2, the uniform neutralising background.
inline Real ewald_net_charge_energy(const std::vector<Real>& charges,
                                    Real volume, Real alpha) {
    Real total = 0.0L;
    for (Real q : charges) total += q;
    if (total == 0.0L) return 0.0L;
    return -kCoulomb * std::numbers::pi_v<Real> / (2.0L * volume * alpha * alpha)
           * total * total;
}

// U_net depends on the cell only through V, and V^-1 has strain derivative
// -d_ab / V, so W_ab = U_net d_ab: isotropic, with the energy on the diagonal.
inline Mat3 ewald_net_charge_virial(const std::vector<Real>& charges,
                                    Real volume, Real alpha) {
    const Real energy = ewald_net_charge_energy(charges, volume, alpha);
    Mat3 w{};
    w[0] = w[4] = w[8] = energy;
    return w;
}

// --- PME reciprocal mesh, general cell ------------------------------------

// The B-spline kernel, the b-moduli and the direct DFT live in
// tests/pme_reference.hpp, alongside the direct-DFT tests that pin the
// engine's transform convention. They are shared rather than duplicated so
// that there is exactly one test-side definition of "what PME is supposed to
// compute" for both files to disagree with the engine about.
// E_rec = 1/2 sum_m G(m) |Qhat(m)|^2, with
//   G(m) = k_e 4 pi / (V k_m^2) exp(-k_m^2 / 4 alpha^2) |b1|^2 |b2|^2 |b3|^2
// and k_m built from the GENERAL reciprocal basis, so this works for a sheared
// cell that `gmd::Box` could never represent. That is the whole point: the
// strain derivative of this function reaches the off-diagonal components.
//
// The charge mesh is rebuilt at every deformed configuration rather than
// carried over. It happens to be invariant -- it depends only on fractional
// coordinates -- but assuming that here would be assuming the very argument
// the production comment makes, and this reference exists to check it
// independently.
inline Real pme_reciprocal_energy(const CellConfiguration& configuration,
                                  Real alpha,
                                  const std::array<int, 3>& grid,
                                  int order) {
    std::vector<std::array<Real, 3>> fractional;
    fractional.reserve(configuration.fractional.size());
    for (const auto& s : configuration.fractional) {
        fractional.push_back({s[0], s[1], s[2]});
    }

    const std::vector<Real> mesh =
        pme_ref::charge_mesh(fractional, configuration.charges, grid, order);

    std::vector<std::complex<Real>> spectral(mesh.size());
    for (std::size_t i = 0; i < mesh.size(); ++i) {
        spectral[i] = std::complex<Real>(mesh[i], 0.0L);
    }
    const auto transformed = pme_ref::direct_dft3(spectral, grid, -1);

    const std::array<std::vector<Real>, 3> moduli = {
        pme_ref::bspline_modulus_squared(grid[0], order),
        pme_ref::bspline_modulus_squared(grid[1], order),
        pme_ref::bspline_modulus_squared(grid[2], order)};

    const Real volume = configuration.volume();
    const Mat3 basis = configuration.reciprocal_basis();
    const Real inv_four_alpha_sq = 1.0L / (4.0L * alpha * alpha);

    Real energy = 0.0L;
    for (int m1 = 0; m1 < grid[0]; ++m1) {
        const int n1 = (m1 <= grid[0] / 2) ? m1 : m1 - grid[0];
        for (int m2 = 0; m2 < grid[1]; ++m2) {
            const int n2 = (m2 <= grid[1] / 2) ? m2 : m2 - grid[1];
            for (int m3 = 0; m3 < grid[2]; ++m3) {
                if (m1 == 0 && m2 == 0 && m3 == 0) continue;
                const int n3 = (m3 <= grid[2] / 2) ? m3 : m3 - grid[2];

                const Vec3 index = {static_cast<Real>(n1), static_cast<Real>(n2),
                                    static_cast<Real>(n3)};
                const Vec3 k = mat_vec(basis, index);
                const Real k_sq = dot(k, k);

                const Real correction = moduli[0][static_cast<std::size_t>(m1)]
                                      * moduli[1][static_cast<std::size_t>(m2)]
                                      * moduli[2][static_cast<std::size_t>(m3)];
                const Real g = kCoulomb * 4.0L * std::numbers::pi_v<Real>
                             / (volume * k_sq)
                             * std::exp(-k_sq * inv_four_alpha_sq) * correction;

                const auto position = static_cast<std::size_t>(
                    m1 * grid[1] * grid[2] + m2 * grid[2] + m3);
                const Real re = transformed[position].real();
                const Real im = transformed[position].imag();
                energy += 0.5L * g * (re * re + im * im);
            }
        }
    }
    return energy;
}

}  // namespace virial_ref
