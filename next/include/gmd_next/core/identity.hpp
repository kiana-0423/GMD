#pragma once

#include "gmd_next/core/status.hpp"

#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace gmd_next::core {

// Stable global identity. It survives sorting, device permutation, migration
// and restart; it is the only key an output record may order rows by.
class AtomId {
public:
    using value_type = std::int64_t;

    constexpr AtomId() = default;
    constexpr explicit AtomId(value_type value) : value_(value) {}

    constexpr value_type value() const { return value_; }
    // Negative values mark "no identity"; inputs must carry nonnegative ids.
    constexpr bool is_valid() const { return value_ >= 0; }

    friend constexpr auto operator<=>(const AtomId&, const AtomId&) = default;

private:
    value_type value_ = -1;
};

// Slot in the current layout. Sorting, migration or capacity growth invalidates
// it, so it must never be stored as an identity or written to a file.
class LocalIndex {
public:
    using value_type = std::uint32_t;
    static constexpr value_type kInvalid = std::numeric_limits<value_type>::max();
    // The highest addressable slot; one value stays reserved for kInvalid.
    static constexpr std::size_t kCapacityLimit = static_cast<std::size_t>(kInvalid);

    constexpr LocalIndex() = default;
    constexpr explicit LocalIndex(value_type value) : value_(value) {}

    // Explicit capacity check, because the device index space is narrower than
    // the host size_t that produced the slot.
    static LocalIndex from_size(std::size_t slot) {
        if (slot >= kCapacityLimit) {
            throw ContractError({ErrorCode::capacity_overflow, "local_index",
                                 "slot " + std::to_string(slot) +
                                     " exceeds the 32-bit device index space"});
        }
        return LocalIndex{static_cast<value_type>(slot)};
    }

    constexpr bool is_valid() const { return value_ != kInvalid; }
    constexpr value_type value() const { return value_; }

    friend constexpr auto operator<=>(const LocalIndex&, const LocalIndex&) = default;

private:
    value_type value_ = kInvalid;
};

// Whether an atom count can be addressed with LocalIndex at all.
inline bool fits_local_index_space(std::size_t count) {
    return count < LocalIndex::kCapacityLimit;
}

}  // namespace gmd_next::core
