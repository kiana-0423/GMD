#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "../system/box.hpp"

namespace gmd {

// Owns per-atom simulation state. In the simpleMD reference program this role
// is carried by Atom: storing coordinates loaded from input together with the
// arrays updated during dynamics, such as velocities and forces. Velocity
// values are typically populated by the system initializer module after input
// loading completes.
class System {
public:
	using Vec3 = std::array<double, 3>;

	System() = default;
	~System() = default;

	void resize(std::size_t atom_count) {
		masses_.assign(atom_count, 0.0);
		coordinates_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
		velocities_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
		forces_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
	}

	std::size_t atom_count() const noexcept {
		return coordinates_.size();
	}

	void set_box(const Box& box) noexcept {
		box_ = box;
	}

	const Box& box() const noexcept {
		return box_;
	}

	Box& mutable_box() noexcept {
		return box_;
	}

	std::span<const double> masses() const noexcept {
		return masses_;
	}

	std::span<double> mutable_masses() noexcept {
		return masses_;
	}

	std::span<const Vec3> coordinates() const noexcept {
		return coordinates_;
	}

	std::span<Vec3> mutable_coordinates() noexcept {
		return coordinates_;
	}

	std::span<const Vec3> velocities() const noexcept {
		return velocities_;
	}

	std::span<Vec3> mutable_velocities() noexcept {
		return velocities_;
	}

	std::span<const Vec3> forces() const noexcept {
		return forces_;
	}

	std::span<Vec3> mutable_forces() noexcept {
		return forces_;
	}

	double potential_energy() const noexcept {
		return potential_energy_;
	}

	void set_potential_energy(double value) noexcept {
		potential_energy_ = value;
	}

private:
	Box box_;
	std::vector<double> masses_;
	std::vector<Vec3> coordinates_;
	std::vector<Vec3> velocities_;
	std::vector<Vec3> forces_;
	double potential_energy_ = 0.0;
};

}  // namespace gmd
