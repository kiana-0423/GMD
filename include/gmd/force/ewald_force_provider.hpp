#pragma once

#include <string_view>

#include "gmd/force/force_provider.hpp"

namespace gmd {

// Ewald summation for long-range Coulomb interactions in a 3-D periodic box.
//
// The Coulomb energy is split into three parts:
//   U = U_real + U_recip + U_self
//
// U_real  – short-range; computed in real space up to `real_cutoff` using a
//           complementary error function (erfc) damping.  Uses the Verlet
//           neighbor list when available.
//
// U_recip – long-range; computed as a sum over k-vectors in reciprocal space
//           up to |k| = kmax along each box axis.  Scales as O(N · K³).
//
// U_self  – analytic self-interaction correction (no forces).
//
// Physical units:
//   Distances [Å],  Energies [eV],  Charges [e],
//   Coulomb constant k_e = 14.3996 eV·Å/e².
//
// Parameter selection:
//   alpha        Controls the real/reciprocal split.  Larger alpha puts more
//                work in k-space.  A good default is alpha ≈ 3.2 / r_cut.
//   kmax         Max k-vector index per dimension.  Auto-selected when ≤ 0 as
//                ceil(alpha · L_max · 3.5 / π).
//   real_cutoff  Real-space cutoff [Å].  When ≤ 0, defaults to
//                0.45 × min(L_x, L_y, L_z) on first use.
class EwaldForceProvider final : public ForceProvider {
public:
    // alpha [1/Å]  (0 = auto-select from real_cutoff)
    // kmax         (0 = auto-select from alpha and box size)
    // real_cutoff  (0 = auto-select from box)
    EwaldForceProvider(double alpha, int kmax, double real_cutoff) noexcept;

    std::string_view name() const noexcept override;
    void initialize(RuntimeContext& runtime) override;
    void compute(const ForceRequest& request, ForceResult& result,
                 RuntimeContext& runtime) override;
    void finalize(RuntimeContext& runtime) override;

    double alpha()       const noexcept { return alpha_; }
    int    kmax()        const noexcept { return kmax_resolved_; }
    double real_cutoff() const noexcept { return real_cutoff_; }

private:
    double alpha_;
    double alpha_sq_;
    double real_cutoff_;
    double real_cutoff_sq_;
    int    kmax_request_;      // user-provided (0 = auto)
    int    kmax_resolved_;     // resolved value used in the last compute() call

    // Resolves auto-parameters from the current box dimensions.
    void resolve_params(const Box& box) noexcept;

    void compute_real_space   (const ForceRequest& req, ForceResult& res) const;
    void compute_reciprocal   (const ForceRequest& req, ForceResult& res) const;
    void compute_self_correction(const ForceRequest& req, ForceResult& res) const;
};

}  // namespace gmd
