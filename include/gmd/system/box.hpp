#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace gmd {

struct Box {
	std::array<double, 3> lengths = {0.0, 0.0, 0.0};
	std::array<double, 3> half_lengths = {0.0, 0.0, 0.0};  //Recording hlaf_lengths for later use in minimum image convention.

	void set_lengths(const std::array<double, 3>& values) {
		validate_lengths(values);
		lengths = values;
		half_lengths = {values[0] * 0.5, values[1] * 0.5, values[2] * 0.5};
	}

	// Every box length must be finite and strictly positive. Zero or negative
	// lengths would make periodic wrapping divide by zero and would turn the
	// ghost-shift loops in the neighbor builder into infinite loops, so they are
	// rejected here at the single point where box dimensions are installed.
	static void validate_lengths(const std::array<double, 3>& values) {
		static constexpr const char* axis_name[3] = {"x", "y", "z"};
		for (std::size_t dim = 0; dim < values.size(); ++dim) {
			if (!std::isfinite(values[dim]) || values[dim] <= 0.0) {
				throw std::invalid_argument(
					std::string("Box length along ") + axis_name[dim] +
					" must be finite and strictly positive, got " +
					std::to_string(values[dim]));
			}
		}
	}
};

// True when two boxes have identical edge lengths. Used to detect whether a
// barostat actually changed the cell (an accepted move) as opposed to leaving
// it untouched (a rejected Monte Carlo trial).
inline bool box_lengths_equal(const Box& a, const Box& b) noexcept {
	return a.lengths[0] == b.lengths[0]
	    && a.lengths[1] == b.lengths[1]
	    && a.lengths[2] == b.lengths[2];
}

}  // namespace gmd
