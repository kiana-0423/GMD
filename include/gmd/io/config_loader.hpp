#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace gmd {

class System;

struct RunConfig {
	std::uint64_t num_steps = 0;
	double time_step = 0.0;
	double temperature = 0.0;
	std::string velocity_init_mode = "random";
	std::uint32_t velocity_seed = 5489u;
	bool remove_center_of_mass_velocity = true;
};

// Per-element LJ parameters.
struct LJElementParams {
	double epsilon = 0.0;  // well depth [eV]
	double sigma   = 0.0;  // zero-crossing distance [Angstrom]
};

// Parameters for a Lennard-Jones (12-6) force field.
//
// Single-element file format:
//   force_field  lj
//   epsilon      0.01032
//   sigma        3.405
//   cutoff       8.5
//
// Multi-element file format (Lorentz-Berthelot mixing rules applied automatically):
//   force_field  lj
//   cutoff       8.5
//   element  Ar   epsilon 0.01032  sigma 3.405
//   element  Ne   epsilon 0.00312  sigma 2.749
//
// In the single-element format `epsilon`/`sigma` are stored under element type 0.
struct LJForceFieldConfig {
	// Cutoff radius shared by all pairs [Angstrom].
	double cutoff = 8.5;

	// Per-element parameters indexed by element type integer (0-based).
	// For single-element configs this always has exactly one entry (index 0).
	std::vector<LJElementParams> elements;

	// Maps element name (e.g. "Ar") → type index into `elements`.
	std::unordered_map<std::string, int> element_name_to_type;

	// Convenience: true when only one element class was specified.
	bool is_single_element() const noexcept { return elements.size() == 1; }

	// Returns epsilon[eV] and sigma[Ang] for the pair (type_i, type_j)
	// using Lorentz-Berthelot mixing rules.
	void pair_params(int type_i, int type_j,
	                 double& eps_out, double& sig_out) const noexcept;

	// Legacy single-element accessors for backward compatibility.
	double epsilon() const noexcept { return elements.empty() ? 0.0 : elements[0].epsilon; }
	double sigma()   const noexcept { return elements.empty() ? 0.0 : elements[0].sigma; }
};

class ConfigLoader {
public:
	ConfigLoader() = default;

	void load_xyz(const std::filesystem::path& xyz_path, System& system) const;
	RunConfig load_run(const std::filesystem::path& run_path) const;

	// Reads a force field parameter file and returns an LJForceFieldConfig.
	// Throws std::runtime_error on parse errors or unsupported force_field type.
	LJForceFieldConfig load_force_field(const std::filesystem::path& ff_path) const;
};

}  // namespace gmd
