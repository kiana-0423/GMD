#pragma once

#include <array>

namespace gmd {

struct Box {
	std::array<double, 3> lengths = {0.0, 0.0, 0.0};
	std::array<double, 3> half_lengths = {0.0, 0.0, 0.0};  //Recording hlaf_lengths for later use in minimum image convention.

	void set_lengths(const std::array<double, 3>& values) noexcept {
		lengths = values;
		half_lengths = {values[0] * 0.5, values[1] * 0.5, values[2] * 0.5};
	}
};

}  // namespace gmd
