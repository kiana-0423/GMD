#include "gmd/boundary/periodic_boundary.hpp"

#include "gmd/system/system.hpp"

namespace gmd {

double wrap_coordinate(double coordinate, double length) noexcept {
    if (length <= 0.0) {
        return coordinate;
    }

    while (coordinate < 0.0) {
        coordinate += length;
    }
    while (coordinate >= length) {
        coordinate -= length;
    }
    return coordinate;
}

void wrap_position(std::array<double, 3>& coordinate, const Box& box) noexcept {
    for (std::size_t dim = 0; dim < 3; ++dim) {
        coordinate[dim] = wrap_coordinate(coordinate[dim], box.lengths[dim]);
    }
}

void wrap_positions(System& system) noexcept {
    auto coordinates = system.mutable_coordinates();
    const auto& box = system.box();
    for (auto& coordinate : coordinates) {
        wrap_position(coordinate, box);
    }
}

}  // namespace gmd