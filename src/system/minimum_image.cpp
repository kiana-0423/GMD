#include "gmd/system/minimum_image.hpp"

#include <cmath>

namespace gmd {

double apply_minimum_image_component(double displacement,
                                     double length,
                                     double half_length) noexcept {
    if (length <= 0.0) {
        return displacement;
    }

    // Coordinates normally remain in the primary cell because the integrator
    // wraps them after every drift. Initial configurations and direct API users
    // are not required to satisfy that precondition, however, so a displacement
    // can span more than one box. Reduce it first, then retain the historical
    // tie convention at exactly +/- half a box.
    if (!std::isfinite(displacement)) {
        return displacement;
    }
    displacement = std::fmod(displacement, length);
    if (displacement < -half_length) {
        displacement += length;
    } else if (displacement > half_length) {
        displacement -= length;
    }
    return displacement;
}

void apply_minimum_image(std::array<double, 3>& displacement, const Box& box) noexcept {
    for (std::size_t dim = 0; dim < 3; ++dim) {
        displacement[dim] = apply_minimum_image_component(
            displacement[dim], box.lengths[dim], box.half_lengths[dim]);
    }
}

}  // namespace gmd
