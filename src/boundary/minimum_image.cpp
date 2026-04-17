#include "gmd/boundary/minimum_image.hpp"

namespace gmd {

double apply_minimum_image_component(double displacement,
                                     double length,
                                     double half_length) noexcept {
    if (length <= 0.0) {
        return displacement;
    }
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