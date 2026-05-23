#pragma once

#include <array>

#include "gmd/system/box.hpp"

namespace gmd {

class System;

double wrap_coordinate(double coordinate, double length) noexcept;
void wrap_position(std::array<double, 3>& coordinate, const Box& box) noexcept;
void wrap_positions(System& system) noexcept;

}  // namespace gmd