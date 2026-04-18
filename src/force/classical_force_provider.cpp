#include "gmd/force/classical_force_provider.hpp"

#include <cmath>
#include <cstddef>

#include "gmd/boundary/minimum_image.hpp"
#include "gmd/io/config_loader.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

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

double ClassicalForceProvider::epsilon() const noexcept {
    return pair_table_.empty() ? 0.0 : pair_table_[0][0].eps4 / 4.0;
}

double ClassicalForceProvider::sigma() const noexcept {
    return pair_table_.empty() ? 0.0 : std::sqrt(pair_table_[0][0].sig2);
}

double ClassicalForceProvider::energy_shift() const noexcept {
    return pair_table_.empty() ? 0.0 : pair_table_[0][0].energy_shift;
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
    result.success = true;
    result.potential_energy = 0.0;
    result.forces.assign(n, Force3D{0.0, 0.0, 0.0});
    result.virial_valid = false;

    if (n == 0 || request.box == nullptr) {
        return;
    }

    const Box& box = *request.box;
    double total_pe = 0.0;

    // Determine if we have per-atom type information.
    const bool multi_element = pair_table_.size() > 1
                               && request.system != nullptr
                               && request.system->atom_types().size() == n;

    // Evaluate one half-pair (i, j) and accumulate.
    auto eval_pair = [&](std::size_t i, std::size_t j) {
        Force3D dr = {
            request.coordinates[i][0] - request.coordinates[j][0],
            request.coordinates[i][1] - request.coordinates[j][1],
            request.coordinates[i][2] - request.coordinates[j][2]
        };
        apply_minimum_image(dr, box);

        const double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        if (r2 >= cutoff_sq_ || r2 < 1e-12) return;

        // Select pair cache by atom types.
        const PairCache& pc = multi_element
            ? pair_table_[static_cast<std::size_t>(request.system->atom_types()[i])]
                         [static_cast<std::size_t>(request.system->atom_types()[j])]
            : pair_table_[0][0];

        const auto [energy, ff] = lj_eval(r2, pc.eps4, pc.sig2);
        total_pe += energy - pc.energy_shift;

        result.forces[i][0] += ff * dr[0];
        result.forces[i][1] += ff * dr[1];
        result.forces[i][2] += ff * dr[2];
        result.forces[j][0] -= ff * dr[0];
        result.forces[j][1] -= ff * dr[1];
        result.forces[j][2] -= ff * dr[2];
    };

    if (request.neighbor_list != nullptr && request.neighbor_list->valid) {
        const NeighborList& nl = *request.neighbor_list;
        for (std::size_t i = 0; i < n; ++i) {
            const int start = nl.offsets[i];
            const int count = nl.counts[i];
            for (int k = 0; k < count; ++k) {
                eval_pair(i, static_cast<std::size_t>(nl.neighbors[start + k]));
            }
        }
    } else {
        for (std::size_t i = 0; i < n - 1; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                eval_pair(i, j);
            }
        }
    }

    result.potential_energy = total_pe;
}

void ClassicalForceProvider::finalize(RuntimeContext& runtime) {
    (void)runtime;
}

}  // namespace gmd
