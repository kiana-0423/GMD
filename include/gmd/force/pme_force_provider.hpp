#pragma once

#include <array>
#include <complex>
#include <string_view>
#include <vector>

#include "gmd/force/force_provider.hpp"

namespace gmd {

// Particle-Mesh Ewald (PME) for long-range Coulomb interactions.
//
// The real-space part is identical to Ewald summation with erfc damping.
// The reciprocal-space part is evaluated efficiently on a 3-D mesh:
//
//   1. Fractional coordinates → B-spline charge assignment onto a K₁×K₂×K₃ mesh.
//   2. 3-D FFT of the charge mesh.
//   3. Multiply each reciprocal-space element by the influence function G(m),
//      which encodes the Coulomb kernel and the B-spline DFT correction.
//   4. Inverse 3-D FFT to obtain the electrostatic potential on the mesh.
//   5. Back-interpolation of forces and energy using B-spline derivative weights.
//
// The internal FFT implementation uses radix-2 Cooley-Tukey; each mesh
// dimension must therefore be a power of 2.
//
// Physical units: Distances [Å], Energies [eV], Charges [e],
//   k_e = 14.3996 eV·Å/e².
//
// Supported B-spline orders: 4 (default, cubic) or 6 (quintic).
class PMEForceProvider final : public ForceProvider {
public:
    // alpha       – Ewald splitting [1/Å]; 0 = auto-select
    // real_cutoff – real-space cutoff [Å]; 0 = auto-select
    // order       – B-spline order (4 or 6)
    // grid        – mesh dimensions; each must be a power of 2
    PMEForceProvider(double alpha, double real_cutoff,
                     int order, std::array<int, 3> grid);

    std::string_view name() const noexcept override;
    void initialize(RuntimeContext& runtime) override;
    void compute(const ForceRequest& request, ForceResult& result,
                 RuntimeContext& runtime) override;
    void finalize(RuntimeContext& runtime) override;

    double alpha()       const noexcept { return alpha_; }
    double real_cutoff() const noexcept { return real_cutoff_; }

private:
    double alpha_;
    double alpha_sq_;
    double real_cutoff_;
    double real_cutoff_sq_;
    int    order_;
    std::array<int, 3> grid_;   // K1, K2, K3

    // Precomputed reciprocal-space influence function (size K1*K2*K3, real).
    std::vector<double> influence_;

    // Working charge/potential mesh (size K1*K2*K3).
    std::vector<std::complex<double>> mesh_;

    // B-spline DFT modulus corrections b_α²(m) for each dimension.
    std::array<std::vector<double>, 3> bmod_sq_;

    // Box lengths stored at last precompute call.
    std::array<double, 3> box_lengths_cache_ = {-1.0, -1.0, -1.0};

    void resolve_params(const Box& box) noexcept;

    // Precompute influence function and b-moduli for the current box.
    void precompute_influence(const Box& box);

    // B-spline evaluation: M_p(u) for u in [0, p].
    static double bspline(double u, int order) noexcept;
    // First derivative: dM_p/du.
    static double bspline_deriv(double u, int order) noexcept;

    // Real-space part (identical to Ewald).
    void compute_real_space(const ForceRequest& req, ForceResult& res) const;

    // Reciprocal-space part via PME mesh.
    void compute_reciprocal_pme(const ForceRequest& req, ForceResult& res);

    // Self-energy correction (energy only).
    void compute_self_correction(const ForceRequest& req, ForceResult& res) const;

    // 3-D in-place FFT / IFFT.
    void fft3d(std::vector<std::complex<double>>& data, bool inverse) const;
    static void fft1d(std::complex<double>* data, int n, bool inverse);
};

}  // namespace gmd
