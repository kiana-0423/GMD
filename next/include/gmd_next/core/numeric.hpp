#pragma once

#include <array>
#include <cmath>
#include <cstddef>

namespace gmd_next::core {

// Host-side value types for validation, records and contracts. Device-resident
// SoA layouts are introduced by the storage milestone and are not these types.
using Vec3 = std::array<double, 3>;
using Tensor3 = std::array<double, 9>;  // row-major, component (a,b) at 3*a+b

inline bool all_finite(double value) {
    return std::isfinite(value);
}

template<std::size_t N>
bool all_finite(const std::array<double, N>& values) {
    for (const double value : values) {
        if (!std::isfinite(value)) return false;
    }
    return true;
}

}  // namespace gmd_next::core
