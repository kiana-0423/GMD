#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include "../system/box.hpp"
#include "../integrator/constraint_solver.hpp"
#include "../system/special_pair_map.hpp"

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
	// Per-pair integer image shift vectors S such that edge_shift = S * box.lengths.
	// image_flags[k] corresponds to neighbors[k].
	std::vector<std::array<int, 3>> image_flags;
	std::vector<std::array<double, 3>> ref_coordinates;  // positions at last rebuild
	bool valid = false;                  // false until first rebuild

	void clear() noexcept {
		counts.clear();
		offsets.clear();
		neighbors.clear();
		image_flags.clear();
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

	void resize(std::size_t atom_count,
	            std::optional<std::size_t> local_count = std::nullopt) {
		const std::size_t new_local_count = local_count.value_or(atom_count);
		if (new_local_count > atom_count) {
			throw std::invalid_argument("Local atom count cannot exceed total atom count");
		}

		masses_.assign(atom_count, 0.0);
		charges_.assign(atom_count, 0.0);
		coordinates_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
		velocities_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
		forces_.assign(atom_count, Vec3{0.0, 0.0, 0.0});
		atom_types_.assign(atom_count, 0);
		molecule_ids_.assign(atom_count, 0);
		atomic_numbers_.assign(atom_count, 0);
		atom_tags_.resize(atom_count);
		for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
			atom_tags_[atom_index] = static_cast<int>(atom_index);
		}
		atom_owners_.assign(atom_count, 0);
		num_local_atoms_ = new_local_count;
		clear_reverse_ghost_virial();
		neighbor_list_.clear();
	}

	std::size_t atom_count() const noexcept {
		return coordinates_.size();
	}

	std::size_t num_local_atoms() const noexcept {
		return num_local_atoms_;
	}

	std::size_t num_ghost_atoms() const noexcept {
		return atom_count() - num_local_atoms_;
	}

	bool is_local_atom(std::size_t atom_index) const noexcept {
		return atom_index < num_local_atoms_;
	}

	int atom_tag(std::size_t atom_index) const noexcept {
		return atom_tags_[atom_index];
	}

	int atom_owner(std::size_t atom_index) const noexcept {
		return atom_owners_[atom_index];
	}

	std::span<int> mutable_atom_tags() noexcept {
		return atom_tags_;
	}

	std::span<int> mutable_atom_owners() noexcept {
		return atom_owners_;
	}

	void mark_atoms_local(int count) {
		if (count < 0 || static_cast<std::size_t>(count) > atom_count()) {
			throw std::invalid_argument("Local atom count must fit the atom storage");
		}

		num_local_atoms_ = static_cast<std::size_t>(count);
	}

	void add_ghost_atom(double mass,
	                    double charge,
	                    const Vec3& position,
	                    int tag,
	                    int owner) {
		masses_.push_back(mass);
		charges_.push_back(charge);
		coordinates_.push_back(position);
		velocities_.push_back(Vec3{0.0, 0.0, 0.0});
		forces_.push_back(Vec3{0.0, 0.0, 0.0});
		atom_types_.push_back(0);
		molecule_ids_.push_back(0);
		atomic_numbers_.push_back(0);
		atom_tags_.push_back(tag);
		atom_owners_.push_back(owner);
		neighbor_list_.clear();
	}

	void clear_ghost_atoms() {
		masses_.resize(num_local_atoms_);
		charges_.resize(num_local_atoms_);
		coordinates_.resize(num_local_atoms_);
		velocities_.resize(num_local_atoms_);
		forces_.resize(num_local_atoms_);
		atom_types_.resize(num_local_atoms_);
		molecule_ids_.resize(num_local_atoms_);
		atomic_numbers_.resize(num_local_atoms_);
		atom_tags_.resize(num_local_atoms_);
		atom_owners_.resize(num_local_atoms_);
		neighbor_list_.clear();
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

	std::span<const int> molecule_ids() const noexcept {
		return molecule_ids_;
	}

	std::span<int> mutable_molecule_ids() noexcept {
		return molecule_ids_;
	}

	std::span<const int> atomic_numbers() const noexcept {
		return atomic_numbers_;
	}

	std::span<int> mutable_atomic_numbers() noexcept {
		return atomic_numbers_;
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

	const std::array<double, 9>& reverse_ghost_virial() const noexcept {
		return reverse_ghost_virial_;
	}

	void clear_reverse_ghost_virial() noexcept {
		reverse_ghost_virial_ = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	}

	void accumulate_reverse_ghost_virial(const std::array<double, 9>& virial) noexcept {
		for (std::size_t index = 0; index < reverse_ghost_virial_.size(); ++index) {
			reverse_ghost_virial_[index] += virial[index];
		}
	}

	double potential_energy() const noexcept {
		return potential_energy_;
	}

	void set_potential_energy(double value) noexcept {
		potential_energy_ = value;
	}

	const std::array<double, 9>& last_virial() const noexcept {
		return last_virial_;
	}

	void set_last_virial(const std::array<double, 9>& virial, bool valid) noexcept {
		last_virial_ = virial;
		last_virial_valid_ = valid;
	}

	bool last_virial_valid() const noexcept {
		return last_virial_valid_;
	}

	void set_last_shake_stats(const ConstraintProjectionStats& stats) noexcept {
		last_shake_stats_ = stats;
	}

	const ConstraintProjectionStats& last_shake_stats() const noexcept {
		return last_shake_stats_;
	}

	void set_last_rattle_stats(const ConstraintProjectionStats& stats) noexcept {
		last_rattle_stats_ = stats;
	}

	const ConstraintProjectionStats& last_rattle_stats() const noexcept {
		return last_rattle_stats_;
	}

	void set_special_pair_map(std::shared_ptr<const SpecialPairMap> special_pairs) noexcept {
		special_pairs_ = std::move(special_pairs);
	}

	const SpecialPairMap* special_pair_map() const noexcept {
		return special_pairs_.get();
	}

	NonbondedScale nonbonded_scale(std::size_t atom_i,
	                               std::size_t atom_j) const noexcept {
		return special_pairs_ != nullptr
			? special_pairs_->scale_for(atom_tag(atom_i), atom_tag(atom_j))
			: NonbondedScale{};
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
	std::vector<int>    molecule_ids_;
	std::vector<int>    atomic_numbers_;
	std::vector<int>    atom_tags_;
	std::vector<int>    atom_owners_;
	std::vector<Vec3> coordinates_;
	std::vector<Vec3> velocities_;
	std::vector<Vec3> forces_;
	std::size_t num_local_atoms_ = 0;
	std::array<double, 9> reverse_ghost_virial_ = {
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0
	};
	std::array<double, 9> last_virial_ = {
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0
	};
	bool last_virial_valid_ = false;
	ConstraintProjectionStats last_shake_stats_;
	ConstraintProjectionStats last_rattle_stats_;
	double potential_energy_ = 0.0;
	std::shared_ptr<const SpecialPairMap> special_pairs_;
	NeighborList neighbor_list_;
};

}  // namespace gmd
