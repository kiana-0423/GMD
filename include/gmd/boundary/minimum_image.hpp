#pragma once

#include <array>

#include "gmd/system/box.hpp"

namespace gmd {

double apply_minimum_image_component(double displacement,
                                     double length,
                                     double half_length) noexcept;
void apply_minimum_image(std::array<double, 3>& displacement, const Box& box) noexcept;

}  // namespace gmd