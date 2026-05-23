#include "gmd/force/classical_force_provider.hpp"

#include <cmath>
#include <cstddef>

#include "gmd/system/minimum_image.hpp"
#include "gmd/io/config_loader.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

namespace {
struct LjEval {
    double energy;
    double force_factor;
};

LjEval lj_eval(double r2, double eps4, double sig2) noexcept {
    const double s2 = sig2 / r2;
    const double s6 = s2 * s2 * s2;
    const double s12 = s6 * s6;
    const double energy = eps4 * (s12 - s6);
    const double force_factor = eps4 * (12.0 * s12 - 6.0 * s6) / r2;
    return {energy, force_factor};
}

#ifdef GMD_ENABLE_MPI
bool mpi_is_available() noexcept {
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    return is_initialized != 0 && is_finalized == 0;
}
#endif
}  // namespace

ClassicalForceProvider::ClassicalForceProvider(double epsilon,
                                               double sigma,
                                               double cutoff) noexcept
    : cutoff_(cutoff), cutoff_sq_(cutoff * cutoff) {
    const double eps4 = 4.0 * epsilon;
    const double sig2 = sigma * sigma;
    const double eshift = lj_eval(cutoff_sq_, eps4, sig2).energy;
    pair_table_ = {{PairCache{eps4, sig2, eshift}}};
}

ClassicalForceProvider::ClassicalForceProvider(const LJForceFieldConfig& config) noexcept
    : cutoff_(config.cutoff), cutoff_sq_(config.cutoff * config.cutoff) {
    build_pair_table(config);
}

void ClassicalForceProvider::build_pair_table(const LJForceFieldConfig& config) noexcept {
    const std::size_t n = config.elements.size();
    pair_table_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        pair_table_[i].resize(n);
        for (std::size_t j = 0; j < n; ++j) {
            double eps_ij, sig_ij;
            config.pair_params(static_cast<int>(i), static_cast<int>(j), eps_ij, sig_ij);
            const double eps4 = 4.0 * eps_ij;
            const double sig2 = sig_ij * sig_ij;
            pair_table_[i][j] = PairCache{eps4, sig2, lj_eval(cutoff_sq_, eps4, sig2).energy};
        }
    }
}

void ClassicalForceProvider::set_params(const LJForceFieldConfig& config) noexcept {
    cutoff_    = config.cutoff;
    cutoff_sq_ = cutoff_ * cutoff_;
    build_pair_table(config);
}

std::string_view ClassicalForceProvider::name() const noexcept {
    return "classical_force_provider";
}

void ClassicalForceProvider::initialize(RuntimeContext& runtime) {
    (void)runtime;
}

void ClassicalForceProvider::compute(const ForceRequest& request,
                                     ForceResult& result,
                                     RuntimeContext& runtime) {
    (void)runtime;

    const std::size_t n = request.coordinates.size();
    const std::size_t local_atom_count = request.system != nullptr
        ? request.system->num_local_atoms()
        : n;
    result.success = true;
    result.potential_energy = 0.0;
    result.forces.assign(n, Force3D{0.0, 0.0, 0.0});
    result.virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    result.virial_valid = true;

    // Force computation (skipped on ranks with no atoms/box).
    if (n != 0 && local_atom_count != 0 && request.box != nullptr) {
        const Box& box = *request.box;
        double total_pe = 0.0;

        // Determine if we have per-atom type information.
        const bool multi_element = pair_table_.size() > 1
                                   && request.system != nullptr
                                   && request.system->atom_types().size() == n;

        // Evaluate one half-pair (i, j) and accumulate.
        auto eval_pair = [&](std::size_t i, std::size_t j) {
            const int atom_i = static_cast<int>(i);
            const int atom_j = static_cast<int>(j);
            const bool evaluate_pair = request.system != nullptr
                ? should_evaluate_pair(*request.system, atom_i, atom_j)
                : should_evaluate_pair(atom_i, atom_j);
            if (!evaluate_pair) {
                return;
            }

            Force3D dr = {
                request.coordinates[i][0] - request.coordinates[j][0],
                request.coordinates[i][1] - request.coordinates[j][1],
                request.coordinates[i][2] - request.coordinates[j][2]
            };
            apply_minimum_image(dr, box);

            const double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if (r2 >= cutoff_sq_ || r2 < 1e-12) return;

            const double lj_scale = request.system != nullptr
                ? request.system->nonbonded_scale(i, j).lj
                : 1.0;
            if (lj_scale == 0.0) return;

            // Select pair cache by atom types.
            const PairCache& pc = multi_element
                ? pair_table_[static_cast<std::size_t>(request.system->atom_types()[i])]
                             [static_cast<std::size_t>(request.system->atom_types()[j])]
                : pair_table_[0][0];

            const auto [energy, ff] = lj_eval(r2, pc.eps4, pc.sig2);
            const double scaled_ff = lj_scale * ff;
            total_pe += lj_scale * (energy - pc.energy_shift);

            result.forces[i][0] += scaled_ff * dr[0];
            result.forces[i][1] += scaled_ff * dr[1];
            result.forces[i][2] += scaled_ff * dr[2];
            result.forces[j][0] -= scaled_ff * dr[0];
            result.forces[j][1] -= scaled_ff * dr[1];
            result.forces[j][2] -= scaled_ff * dr[2];

            // Pair virial tensor contribution W_ab = r_a * F_b.
            result.virial[0] += dr[0] * (scaled_ff * dr[0]);
            result.virial[1] += dr[0] * (scaled_ff * dr[1]);
            result.virial[2] += dr[0] * (scaled_ff * dr[2]);
            result.virial[3] += dr[1] * (scaled_ff * dr[0]);
            result.virial[4] += dr[1] * (scaled_ff * dr[1]);
            result.virial[5] += dr[1] * (scaled_ff * dr[2]);
            result.virial[6] += dr[2] * (scaled_ff * dr[0]);
            result.virial[7] += dr[2] * (scaled_ff * dr[1]);
            result.virial[8] += dr[2] * (scaled_ff * dr[2]);
        };

        if (request.neighbor_list != nullptr && request.neighbor_list->valid) {
            const NeighborList& nl = *request.neighbor_list;
            for (std::size_t i = 0; i < local_atom_count; ++i) {
                const int start = nl.offsets[i];
                const int count = nl.counts[i];
                for (int k = 0; k < count; ++k) {
                    eval_pair(i, static_cast<std::size_t>(nl.neighbors[start + k]));
                }
            }
        } else {
            for (std::size_t i = 0; i < local_atom_count; ++i) {
                for (std::size_t j = i + 1; j < n; ++j) {
                    eval_pair(i, j);
                }
            }
        }

        result.potential_energy = total_pe;
    }

    // MPI: virial allreduce — always called on all ranks because MPI_Allreduce
    // is a collective operation that must be matched by every rank in the
    // communicator, even those with zero atoms.
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        auto local_virial = result.virial;
        MPI_Allreduce(local_virial.data(),
                      result.virial.data(),
                      static_cast<int>(result.virial.size()),
                      MPI_DOUBLE,
                      MPI_SUM,
                      MPI_COMM_WORLD);
    }
#endif
}

bool ClassicalForceProvider::should_evaluate_pair(int i, int j) const noexcept {
    return i >= 0 && j > i;
}

bool ClassicalForceProvider::should_evaluate_pair(const System& system,
                                                 int i,
                                                 int j) const noexcept {
    if (!should_evaluate_pair(i, j)) {
        return false;
    }

    const auto atom_i = static_cast<std::size_t>(i);
    const auto atom_j = static_cast<std::size_t>(j);
    if (system.is_local_atom(atom_j)) {
        return true;
    }

    // A local/ghost boundary pair can exist on both ranks. Global tags choose
    // the rank that owns the lower-tagged local side for energy and virial.
    return system.atom_tag(atom_i) < system.atom_tag(atom_j);
}

void ClassicalForceProvider::finalize(RuntimeContext& runtime) {
    (void)runtime;
}

}  // namespace gmd
