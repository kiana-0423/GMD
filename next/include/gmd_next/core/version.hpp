#pragma once

#include "gmd_next/core/status.hpp"

#include <array>
#include <bit>
#include <cstdint>
#include <initializer_list>
#include <string_view>

namespace gmd_next::core {

// What can change underneath a cached result.
//   coordinates  every position update, including each drift
//   cell         box edits; the first production scope keeps the box fixed
//   model        potential parameters, cutoff mode, cutoff or skin
//   permutation  device reordering, so slot-indexed data stops lining up
enum class VersionDomain : std::uint8_t { coordinates, cell, model, permutation };

inline constexpr std::array<VersionDomain, 4> kVersionDomains{
    VersionDomain::coordinates, VersionDomain::cell, VersionDomain::model,
    VersionDomain::permutation};

std::string_view name(VersionDomain domain);

// A counter for one domain. The domain is part of the type so that a cell
// version cannot be compared with, or assigned from, a coordinate version.
// Counter 0 means "unset" and matches nothing, including another unset value.
template<VersionDomain D>
class Version {
public:
    static constexpr VersionDomain domain = D;

    constexpr Version() = default;
    constexpr explicit Version(std::uint64_t counter) : counter_(counter) {}

    constexpr bool is_set() const { return counter_ != 0; }
    constexpr std::uint64_t counter() const { return counter_; }
    constexpr Version next() const { return Version{counter_ + 1}; }

    friend constexpr bool operator==(const Version&, const Version&) = default;

private:
    std::uint64_t counter_ = 0;
};

using CoordinateVersion = Version<VersionDomain::coordinates>;
using CellVersion = Version<VersionDomain::cell>;
using ModelVersion = Version<VersionDomain::model>;
using PermutationVersion = Version<VersionDomain::permutation>;

class VersionDomainSet {
public:
    constexpr VersionDomainSet() = default;
    constexpr VersionDomainSet(std::initializer_list<VersionDomain> domains) {
        for (const auto domain : domains) insert(domain);
    }

    constexpr VersionDomainSet& insert(VersionDomain domain) {
        bits_ = static_cast<std::uint8_t>(bits_ | bit(domain));
        return *this;
    }
    constexpr bool contains(VersionDomain domain) const { return (bits_ & bit(domain)) != 0; }
    constexpr bool empty() const { return bits_ == 0; }
    constexpr std::size_t size() const { return static_cast<std::size_t>(std::popcount(bits_)); }
    constexpr VersionDomainSet operator|(VersionDomainSet other) const {
        VersionDomainSet result;
        result.bits_ = static_cast<std::uint8_t>(bits_ | other.bits_);
        return result;
    }
    constexpr VersionDomainSet operator&(VersionDomainSet other) const {
        VersionDomainSet result;
        result.bits_ = static_cast<std::uint8_t>(bits_ & other.bits_);
        return result;
    }
    friend constexpr bool operator==(VersionDomainSet, VersionDomainSet) = default;

private:
    static constexpr std::uint8_t bit(VersionDomain domain) {
        return static_cast<std::uint8_t>(1u << static_cast<unsigned>(domain));
    }
    std::uint8_t bits_ = 0;
};

// A slot-indexed force array is stale if any of these moved.
inline constexpr VersionDomainSet kCachedForceDomains{
    VersionDomain::coordinates, VersionDomain::cell, VersionDomain::model,
    VersionDomain::permutation};
// A neighbour list is not invalidated by coordinates alone: that is what the
// skin and the tracked displacement bound are for. Any of these does end it.
inline constexpr VersionDomainSet kNeighborListDomains{
    VersionDomain::cell, VersionDomain::model, VersionDomain::permutation};

// The four counters as recorded at one moment, stored next to a cached result.
struct VersionStamp {
    CoordinateVersion coordinates{};
    CellVersion cell{};
    ModelVersion model{};
    PermutationVersion permutation{};

    bool is_complete() const;
    // Domains that differ, or that are unset on either side.
    VersionDomainSet mismatches(const VersionStamp& other) const;
    // Complete on both sides and equal in every domain.
    bool matches(const VersionStamp& other) const;

    friend bool operator==(const VersionStamp&, const VersionStamp&) = default;
};

// Owner-side counters. The state owner advances a domain whenever it changes
// that state; consumers only ever copy stamps.
class StateVersions {
public:
    StateVersions() = default;

    CoordinateVersion coordinates() const { return CoordinateVersion{counters_[0]}; }
    CellVersion cell() const { return CellVersion{counters_[1]}; }
    ModelVersion model() const { return ModelVersion{counters_[2]}; }
    PermutationVersion permutation() const { return PermutationVersion{counters_[3]}; }

    void advance(VersionDomain domain);
    VersionStamp stamp() const;

private:
    // The first valid version is 1, so a default VersionStamp never matches.
    std::array<std::uint64_t, 4> counters_{1, 1, 1, 1};
};

// Consumer-side guard: report every listed domain in which the recorded stamp
// no longer describes the current state.
ValidationReport check_versions(const VersionStamp& recorded, const VersionStamp& current,
                                VersionDomainSet required);

}  // namespace gmd_next::core
