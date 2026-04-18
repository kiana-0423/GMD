#pragma once

#include <string_view>
#include <vector>

#include "gmd/force/force_provider.hpp"

namespace gmd {

// Forward declaration — full definition in gmd/io/config_loader.hpp
struct LJForceFieldConfig;

// Lennard-Jones (12-6) force field with a potential shift at the cutoff.
// V(r) = 4*epsilon * [(sigma/r)^12 - (sigma/r)^6] - V(r_cut)
// Supports multiple atom types via Lorentz-Berthelot combining rules.
class ClassicalForceProvider final : public ForceProvider {
public:
    // Per-pair precomputed constants.
    struct PairCache {
        double eps4;          // 4 * epsilon_ij
        double sig2;          // sigma_ij^2
        double energy_shift;  // V(r_cut) for the ij pair
    };

    // Construct with explicit single-element parameters (epsilon [eV], sigma [Ang], cutoff [Ang]).
    explicit ClassicalForceProvider(double epsilon = 0.01032,
                                    double sigma   = 3.405,
                                    double cutoff  = 8.5) noexcept;

    // Construct directly from a loaded config struct (supports multi-element).
    explicit ClassicalForceProvider(const LJForceFieldConfig& config) noexcept;

    ~ClassicalForceProvider() override = default;

    std::string_view name() const noexcept override;

    void initialize(RuntimeContext& runtime) override;
    void compute(const ForceRequest& request, ForceResult& result, RuntimeContext& runtime) override;
    void finalize(RuntimeContext& runtime) override;

    // --- Parameter accessors (single-element convenience) ---
    double epsilon()      const noexcept;
    double sigma()        const noexcept;
    double cutoff()       const noexcept { return cutoff_; }
    double energy_shift() const noexcept;

    // Replace all parameters at once and rebuild the pair table.
    void set_params(const LJForceFieldConfig& config) noexcept;

    // Access the full pair table (num_types × num_types).
    const std::vector<std::vector<PairCache>>& pair_table() const noexcept { return pair_table_; }

private:
    double cutoff_;        // cutoff radius [Ang]
    double cutoff_sq_;     // cutoff^2 (cached)

    // Pair parameters for all type × type combinations.
    // pair_table_[type_i][type_j] is valid for type_i <= type_j (symmetric).
    std::vector<std::vector<PairCache>> pair_table_;

    // Build pair_table_ from an LJForceFieldConfig.
    void build_pair_table(const LJForceFieldConfig& config) noexcept;
};

}  // namespace gmd
