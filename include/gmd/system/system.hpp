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

	void set_box(const Box& box) {
		// set_box() copies a whole Box, bypassing Box::set_lengths(), so the
		// dimensions are validated here as well. This is the other entry point
		// through which a degenerate box could reach wrapping and ghost shifts.
		Box::validate_lengths(box.lengths);
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

	// --- Provider virial + constraint virial ---------------------------------
	//
	// last_virial() is the COMBINED tensor everything that reports pressure
	// reads. It is assembled here, in one place, from two parts kept separately
	// so that either can be replaced without losing or double-counting the
	// other:
	//
	//   provider_virial_    the force providers' tensor, evaluated at the
	//                       current coordinates;
	//   constraint_virial_  the step's constraint contribution, recovered from
	//                       the converged RATTLE multipliers.
	//
	// TIME LEVELS. The two must belong to the same state. Within a step the
	// providers run at r(t+dt) and RATTLE runs afterwards on those same
	// coordinates, so they match. Anything that installs a provider virial for a
	// *different* geometry -- a barostat rescale, or the initial evaluation of a
	// fresh run -- goes through set_provider_virial(), which drops the constraint
	// part to Unavailable. While it is Unavailable the combined tensor is
	// reported invalid rather than being passed off as a complete pressure
	// virial, because a provider-only virial is not one when constraints act.
	//
	// With no constraints in the run the state is NotApplicable and the combined
	// value is the provider virial unchanged, so unconstrained runs are
	// bit-for-bit unaffected by any of this.
	const std::array<double, 9>& provider_virial() const noexcept {
		return provider_virial_;
	}

	bool provider_virial_valid() const noexcept { return provider_virial_valid_; }

	// Installs the force-provider virial for the current coordinates. Any
	// constraint contribution that was attached to a previous evaluation is
	// dropped: it does not necessarily belong to this geometry, and deciding
	// that here is not possible. RATTLE re-attaches the contemporaneous one.
	void set_provider_virial(const std::array<double, 9>& virial, bool valid) noexcept {
		provider_virial_ = virial;
		provider_virial_valid_ = valid;
		if (constraint_virial_state_ == ConstraintVirialState::Valid) {
			constraint_virial_state_ = ConstraintVirialState::Unavailable;
		}
		constraint_virial_ = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		refresh_combined_virial();
	}

	// Declares whether constraints act in this run. Switching them off makes the
	// provider virial complete again; switching them on only records that a
	// constraint term is now required, without inventing one.
	void set_constraints_active(bool active) noexcept {
		if (!active) {
			constraint_virial_ = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			constraint_virial_state_ = ConstraintVirialState::NotApplicable;
			refresh_combined_virial();
		} else if (constraint_virial_state_ == ConstraintVirialState::NotApplicable) {
			constraint_virial_state_ = ConstraintVirialState::Unavailable;
			refresh_combined_virial();
		}
	}

	const std::array<double, 9>& constraint_virial() const noexcept {
		return constraint_virial_;
	}

	ConstraintVirialState constraint_virial_state() const noexcept {
		return constraint_virial_state_;
	}

	bool constraint_virial_valid() const noexcept {
		return constraint_virial_state_ == ConstraintVirialState::Valid;
	}

	// Attaches the constraint contribution belonging to the currently installed
	// provider virial. Called by the integrator right after RATTLE, and by a
	// restart restoring the value the checkpointed step ended with.
	void set_constraint_virial(const std::array<double, 9>& virial) noexcept {
		constraint_virial_ = virial;
		constraint_virial_state_ = ConstraintVirialState::Valid;
		refresh_combined_virial();
	}

	// Constraints act, but no multiplier belongs to the installed provider
	// virial. The combined tensor becomes invalid.
	void mark_constraint_virial_unavailable() noexcept {
		constraint_virial_ = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		if (constraint_virial_state_ == ConstraintVirialState::Valid) {
			constraint_virial_state_ = ConstraintVirialState::Unavailable;
		}
		refresh_combined_virial();
	}

	// --- Completed-step thermodynamics ---------------------------------------
	//
	// The thermodynamic state OF THE STEP THAT JUST FINISHED, and the force
	// state attached to the CURRENT GEOMETRY, are different things, and a
	// barostat separates them: it rescales the cell after the step is complete,
	// and the forces are then re-evaluated so the next step starts consistent.
	// That re-evaluation produces a virial and a potential energy for a geometry
	// no dynamics ever integrated, at a volume the completed step never had.
	//
	// So the completed step's numbers are captured together, once, at the end of
	// finish_step() and before any barostat runs, and are NOT touched by a later
	// force evaluation. Reported as a set they are mutually consistent:
	// P = (2K + tr W) / 3V holds among exactly these values, and E_total = PE + K
	// is the energy of the state that produced them. The pressure here is also
	// bit-for-bit the value the barostat itself consumed.
	//
	// last_virial() remains the current-geometry tensor and is not this.
	struct StepThermodynamics {
		bool valid = false;
		double pressure = 0.0;               // (2K + tr W) / 3V  [eV/A^3]
		std::array<double, 9> virial{};      // provider + constraint, both at t+dt
		double twice_kinetic_energy = 0.0;
		double volume = 0.0;
		double potential_energy = 0.0;
	};

	const StepThermodynamics& step_thermodynamics() const noexcept {
		return step_thermodynamics_;
	}

	void set_step_thermodynamics(const StepThermodynamics& value) noexcept {
		step_thermodynamics_ = value;
	}

	void clear_step_thermodynamics() noexcept {
		step_thermodynamics_ = StepThermodynamics{};
	}

	// Everything a trajectory frame needs that is NOT per-atom.
	//
	// This exists for the MPI output path, which writes a separate System holding
	// gathered global coordinates. That System takes part in no dynamics, so
	// without this it would report the state it was copied from -- the initial
	// box, no constraint diagnostics, no completed-step thermodynamics -- while
	// the real System moved on. Keeping the list here, next to the members it
	// copies, is what stops a newly added frame field from being forgotten there.
	//
	// Nothing here is reduced across ranks: every field is either replicated (the
	// box, the constraint virial) or already global (the provider virial, which
	// the providers allreduce, and the completed-step record, which is built from
	// globally reduced quantities). Reducing again would multiply them.
	void copy_frame_state_from(const System& other) {
		set_box(other.box());
		set_potential_energy(other.potential_energy());
		set_step_thermodynamics(other.step_thermodynamics());
		set_last_shake_stats(other.last_shake_stats());
		set_last_rattle_stats(other.last_rattle_stats());
		// Provider first: installing one drops any constraint term, so the
		// constraint state has to be re-established afterwards.
		set_provider_virial(other.provider_virial_, other.provider_virial_valid_);
		switch (other.constraint_virial_state_) {
			case ConstraintVirialState::NotApplicable:
				set_constraints_active(false);
				break;
			case ConstraintVirialState::Unavailable:
				set_constraints_active(true);
				mark_constraint_virial_unavailable();
				break;
			case ConstraintVirialState::Valid:
				set_constraints_active(true);
				set_constraint_virial(other.constraint_virial_);
				break;
		}
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
	// The single place the reported virial is assembled. Invalid whenever a
	// constraint term is required but not available for this geometry.
	void refresh_combined_virial() noexcept {
		std::array<double, 9> combined = provider_virial_;
		if (constraint_virial_state_ == ConstraintVirialState::Valid) {
			for (std::size_t index = 0; index < combined.size(); ++index) {
				combined[index] += constraint_virial_[index];
			}
		}
		set_last_virial(combined,
		                provider_virial_valid_ &&
		                constraint_virial_state_ != ConstraintVirialState::Unavailable);
	}

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
	std::array<double, 9> constraint_virial_ = {
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0
	};
	std::array<double, 9> provider_virial_ = {
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0
	};
	bool provider_virial_valid_ = false;
	ConstraintVirialState constraint_virial_state_ = ConstraintVirialState::NotApplicable;
	StepThermodynamics step_thermodynamics_;
	ConstraintProjectionStats last_shake_stats_;
	ConstraintProjectionStats last_rattle_stats_;
	double potential_energy_ = 0.0;
	std::shared_ptr<const SpecialPairMap> special_pairs_;
	NeighborList neighbor_list_;
};

}  // namespace gmd
