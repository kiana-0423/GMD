#pragma once

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "gmd/system/topology.hpp"

namespace gmd {

// Scale factors applied to non-bonded interactions between topology-related
// atoms. A value of zero is an exclusion; one is an unmodified interaction.
struct NonbondedScale {
    double lj = 1.0;
    double coulomb = 1.0;
};

struct SpecialPairScaleConfig {
    NonbondedScale pair_12{0.0, 0.0};
    NonbondedScale pair_13{0.0, 0.0};
    NonbondedScale pair_14{0.5, 0.833333333333};
};

struct SpecialPair {
    int atom_tag_a = 0;
    int atom_tag_b = 0;
    NonbondedScale scale;
};

class SpecialPairMap {
public:
    SpecialPairMap() = default;

    explicit SpecialPairMap(const Topology& topology,
                            const SpecialPairScaleConfig& scales = {}) {
        // Insert shortest graph separations first: an atom pair that also
        // appears in a longer topology term remains governed by its closest
        // bonded relationship.
        for (const auto& bond : topology.bonds) {
            insert_if_absent(bond.i, bond.j, scales.pair_12);
        }
        for (const auto& angle : topology.angles) {
            insert_if_absent(angle.i, angle.k, scales.pair_13);
        }
        for (const auto& dihedral : topology.dihedrals) {
            insert_if_absent(dihedral.i, dihedral.l, scales.pair_14);
        }
    }

    NonbondedScale scale_for(int atom_tag_a, int atom_tag_b) const noexcept {
        const auto found = pairs_.find(make_key(atom_tag_a, atom_tag_b));
        return found == pairs_.end() ? NonbondedScale{} : found->second;
    }

    bool empty() const noexcept {
        return pairs_.empty();
    }

    std::size_t size() const noexcept {
        return pairs_.size();
    }

    const std::vector<SpecialPair>& entries() const noexcept {
        return entries_;
    }

private:
    static std::uint64_t make_key(int atom_tag_a, int atom_tag_b) noexcept {
        const auto lo = static_cast<std::uint32_t>(std::min(atom_tag_a, atom_tag_b));
        const auto hi = static_cast<std::uint32_t>(std::max(atom_tag_a, atom_tag_b));
        return (static_cast<std::uint64_t>(lo) << 32U) |
               static_cast<std::uint64_t>(hi);
    }

    void insert_if_absent(int atom_tag_a,
                          int atom_tag_b,
                          const NonbondedScale& scale) {
        const int lo = std::min(atom_tag_a, atom_tag_b);
        const int hi = std::max(atom_tag_a, atom_tag_b);
        const auto [position, inserted] = pairs_.try_emplace(make_key(lo, hi), scale);
        (void)position;
        if (inserted) {
            entries_.push_back(SpecialPair{lo, hi, scale});
        }
    }

    std::unordered_map<std::uint64_t, NonbondedScale> pairs_;
    std::vector<SpecialPair> entries_;
};

}  // namespace gmd
