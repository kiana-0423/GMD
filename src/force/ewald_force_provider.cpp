#include "gmd/force/ewald_force_provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

#include "gmd/boundary/minimum_image.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

// Coulomb constant k_e = e² / (4π ε₀) in units of [eV · Å / e²].
static constexpr double kEwaldCoulomb = 14.3996;

namespace {

void accumulate_coordinate_virial(const std::span<const Coordinate3D> coordinates,
                                  const std::vector<Force3D>& forces,
                                  std::array<double, 9>& virial) noexcept {
    const std::size_t n = std::min(coordinates.size(), forces.size());
    for (std::size_t i = 0; i < n; ++i) {
        const auto& r = coordinates[i];
        const auto& f = forces[i];
        virial[0] += r[0] * f[0];
        virial[1] += r[0] * f[1];
        virial[2] += r[0] * f[2];
        virial[3] += r[1] * f[0];
        virial[4] += r[1] * f[1];
        virial[5] += r[1] * f[2];
        virial[6] += r[2] * f[0];
        virial[7] += r[2] * f[1];
        virial[8] += r[2] * f[2];
    }
}

#ifdef GMD_ENABLE_MPI
bool mpi_is_available() noexcept {
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    return is_initialized != 0 && is_finalized == 0;
}

int mpi_size() noexcept {
    if (!mpi_is_available()) return 1;
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

void allreduce_complex_sum(double& re, double& im) noexcept {
    if (!mpi_is_available()) return;
    double buf[2] = {re, im};
    double global[2] = {0.0, 0.0};
    MPI_Allreduce(buf, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    re = global[0];
    im = global[1];
}

double allreduce_scalar_ewald(double local) noexcept {
    if (!mpi_is_available()) return local;
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

void allreduce_virial(std::array<double, 9>& virial) noexcept {
    if (!mpi_is_available()) return;
    auto local = virial;
    MPI_Allreduce(local.data(),
                  virial.data(),
                  static_cast<int>(virial.size()),
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD);
}
#else
int mpi_size() noexcept {
    return 1;
}
#endif

}  // namespace

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
    res.virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    res.virial_valid = true;

    // Determine whether this rank can meaningfully contribute.
    // (n may be zero on ranks with no atoms in their domain.)
    bool can_compute = (n > 0 && req.box != nullptr && req.system != nullptr);
    if (can_compute) {
        const auto charges = req.system->charges();
        if (charges.size() != n) {
            res.virial_valid = false;
            can_compute = false;
        }
    }

    // Coordinated charge-existence check — MPI_Allreduce ensures every rank
    // reaches the same decision, even ranks with zero atoms.
    bool has_charges = false;
    if (can_compute) {
        const std::size_t num_local = req.system->num_local_atoms();
        const auto charges = req.system->charges();
        for (std::size_t i = 0; i < num_local; ++i) {
            if (charges[i] != 0.0) { has_charges = true; break; }
        }
    }
#ifdef GMD_ENABLE_MPI
    {
        int local_has = has_charges ? 1 : 0;
        int global_has = 0;
        MPI_Allreduce(&local_has, &global_has, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        has_charges = global_has != 0;
    }
#endif

    if (can_compute && has_charges) {
        resolve_params(*req.box);

        compute_real_space(req, res);
        compute_reciprocal(req, res);
        compute_self_correction(req, res);
        accumulate_coordinate_virial(req.coordinates, res.forces, res.virial);
    }

#ifdef GMD_ENABLE_MPI
    allreduce_virial(res.virial);
#endif
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

        // In MPI mode, a local/ghost boundary pair can exist on both ranks.
        // Use global tags to evaluate exactly once.
        if (req.system != nullptr) {
            const std::size_t num_local = req.system->num_local_atoms();
            const bool i_is_local = i < num_local;
            const bool j_is_local = j < num_local;
            if (i_is_local && !j_is_local) {
                if (req.system->atom_tag(i) > req.system->atom_tag(j)) return;
            } else if (!i_is_local && j_is_local) {
                if (req.system->atom_tag(j) < req.system->atom_tag(i)) return;
            }
        }

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
        const std::size_t num_local = req.system != nullptr
            ? req.system->num_local_atoms() : n;
        for (std::size_t i = 0; i < num_local; ++i) {
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
    const std::size_t num_local = req.system != nullptr
        ? req.system->num_local_atoms() : n;

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

                // Local structure factor S_local(k) — only local atoms.
                double S_re = 0.0, S_im = 0.0;
                for (std::size_t j = 0; j < num_local; ++j) {
                    const double phi = kx*coords[j][0]
                                     + ky*coords[j][1]
                                     + kz*coords[j][2];
                    S_re += charges[j] * std::cos(phi);
                    S_im += charges[j] * std::sin(phi);
                }

#ifdef GMD_ENABLE_MPI
                allreduce_complex_sum(S_re, S_im);
#endif

                // Energy: ½ · gfactor · |S|²
                res.potential_energy +=
                    0.5 * gfactor * (S_re*S_re + S_im*S_im)
                    / static_cast<double>(mpi_size());

                // Force on each local atom.
                for (std::size_t i = 0; i < num_local; ++i) {
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

    // Reciprocal-space virial via coordinate-force outer product.
    for (std::size_t i = 0; i < num_local; ++i) {
        const auto& r = coords[i];
        const auto& f = res.forces[i];
        res.virial[0] += r[0] * f[0];
        res.virial[1] += r[0] * f[1];
        res.virial[2] += r[0] * f[2];
        res.virial[3] += r[1] * f[0];
        res.virial[4] += r[1] * f[1];
        res.virial[5] += r[1] * f[2];
        res.virial[6] += r[2] * f[0];
        res.virial[7] += r[2] * f[1];
        res.virial[8] += r[2] * f[2];
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
    const std::size_t num_local = req.system != nullptr
        ? req.system->num_local_atoms() : charges.size();

    double q2_sum = 0.0, Q_net = 0.0;
    for (std::size_t i = 0; i < num_local; ++i) {
        q2_sum += charges[i] * charges[i];
        Q_net  += charges[i];
    }

#ifdef GMD_ENABLE_MPI
    q2_sum = allreduce_scalar_ewald(q2_sum);
    Q_net  = allreduce_scalar_ewald(Q_net);
#endif

    res.potential_energy -=
        kEwaldCoulomb * (alpha_ / std::sqrt(std::numbers::pi)) * q2_sum
        / static_cast<double>(mpi_size());

    if (Q_net != 0.0) {
        const double V = req.box->lengths[0]
                       * req.box->lengths[1]
                       * req.box->lengths[2];
        res.potential_energy -=
            kEwaldCoulomb * std::numbers::pi / (2.0 * V * alpha_sq_) * Q_net * Q_net
            / static_cast<double>(mpi_size());
    }
}

}  // namespace gmd
