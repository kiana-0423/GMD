#pragma once

#include "gmd_next/core/numeric.hpp"

#include <array>
#include <optional>

namespace gmd_next::model {

// Fixed orthorhombic cell. The production scope keeps it periodic on all three
// axes and constant for the whole run; tilted cells are a separate model.
struct CellSpec {
    core::Vec3 lengths{0.0, 0.0, 0.0};  // Angstrom under UnitSystem::metal
    std::array<bool, 3> periodic{true, true, true};

    friend bool operator==(const CellSpec&, const CellSpec&) = default;
};

bool is_fully_periodic(const CellSpec& cell);
// Shortest periodic edge, absent when no axis is periodic.
std::optional<double> shortest_periodic_length(const CellSpec& cell);
double volume(const CellSpec& cell);

}  // namespace gmd_next::model
