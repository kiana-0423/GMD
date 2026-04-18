#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "../system/box.hpp"

namespace gmd {

// Compressed-sparse-row neighbor list stored inside System.
// After a NeighborBuilder::rebuild() call:
//   - atom i's neighbors are neighbors[offsets[i] .. offsets[i]+counts[i])
//   - each entry is a half-pair (j > i for Verlet lists built with Newton III)
//   - ref_coordinates holds the positions at last rebuild, used for skin check
struct NeighborList {
	std::vector<int>   counts;           // counts[i]  = number of neighbors of atom i
	std::vector<int>   offsets;          // offsets[i] = start index in `neighbors`
	std::vector<int>   neighbors;        // flat neighbor index storage
	std::vector<std::array<double, 3>> ref_coordinates;  // positions at last rebuild
	bool valid = false;                  // false until first rebuild

	void clear() noexcept {
		counts.clear();
		offsets.clear();
		neighbors.clear();
		ref_coordinates.clear();
		valid = false;
	}
};

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
		charges_.assign(atom_count, 0.0);
		coordinates_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
		velocities_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
		forces_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
		atom_types_.assign(atom_count, 0);
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

	std::span<const double> charges() const noexcept {
		return charges_;
	}

	std::span<double> mutable_charges() noexcept {
		return charges_;
	}

	std::span<const int> atom_types() const noexcept {
		return atom_types_;
	}

	std::span<int> mutable_atom_types() noexcept {
		return atom_types_;
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

	const NeighborList& neighbor_list() const noexcept {
		return neighbor_list_;
	}

	NeighborList& mutable_neighbor_list() noexcept {
		return neighbor_list_;
	}

private:
	Box box_;
	std::vector<double> masses_;
	std::vector<double> charges_;
	std::vector<int>    atom_types_;
	std::vector<Vec3> coordinates_;
	std::vector<Vec3> velocities_;
	std::vector<Vec3> forces_;
	double potential_energy_ = 0.0;
	NeighborList neighbor_list_;
};

}  // namespace gmd
