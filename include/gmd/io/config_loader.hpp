#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gmd/force/bonded_params.hpp"
#include "gmd/system/topology.hpp"

namespace gmd {

class System;

// Per-type LJ parameters (one entry per type number, stored 0-based internally).
struct LJElementParams {
	std::string element;   // element symbol (e.g. "Ar"), used for mass lookup
	double mass    = 0.0;  // atomic mass [amu]
	double epsilon = 0.0;  // LJ well depth [eV]
	double sigma   = 0.0;  // LJ zero-crossing distance [Angstrom]
	double charge  = 0.0;  // partial charge [e] (elementary charge units)
};

// Parameters for a Lennard-Jones (12-6) force field.
//
// Single-element file format (legacy, no explicit type number):
//   force_field  lj
//   epsilon      0.01032
//   sigma        3.405
//   cutoff       8.5
//
// Multi-type format (type numbers are 1-based, like LAMMPS):
//   force_field  lj
//   cutoff       8.5
//   type 1  Ar  epsilon 0.01032  sigma 3.405
//   type 2  Ne  epsilon 0.00312  sigma 2.749
//   # Same element, different parameters:
//   type 3  Ar  epsilon 0.00800  sigma 3.200
//
// Lorentz-Berthelot mixing rules are applied automatically for cross-type pairs.
struct LJForceFieldConfig {
	// Cutoff radius shared by all pairs [Angstrom].
	double cutoff = 8.5;

	// Per-type parameters, indexed 0-based (type N in the file → index N-1).
	std::vector<LJElementParams> elements;

	// Maps element symbol → last assigned type index (0-based).
	// Useful when xyz input uses element symbols and each symbol is unambiguous.
	std::unordered_map<std::string, int> element_name_to_type;

	// Returns epsilon[eV] and sigma[Ang] for the pair (type_i, type_j)
	// using Lorentz-Berthelot mixing rules. Indices are 0-based.
	void pair_params(int type_i, int type_j,
	                 double& eps_out, double& sig_out) const noexcept;
};

// Configuration for long-range electrostatics.
//
// Ewald summation:
//   coulomb        ewald
//   ewald_alpha    0.3        # splitting parameter [1/Å]; 0 = auto
//   ewald_kmax     7          # max k-vector index; 0 = auto
//   ewald_cutoff   10.0       # real-space cutoff [Å]; 0 = use LJ cutoff
//
// Particle-Mesh Ewald (PME):
//   coulomb        pme
//   pme_alpha      0.3        # splitting parameter [1/Å]; 0 = auto
//   pme_cutoff     10.0       # real-space cutoff [Å]; 0 = use LJ cutoff
//   pme_order      4          # B-spline order (4 or 6)
//   pme_grid       32 32 32   # mesh dimensions (each must be a power of 2)
struct CoulombConfig {
	std::string method = "ewald";         // "ewald" or "pme"
	double alpha        = 0.0;            // Ewald splitting [1/Å]; 0 = auto
	double real_cutoff  = 0.0;            // real-space cutoff [Å]; 0 = use LJ cutoff
	// Ewald-specific
	int kmax = 0;                         // max k-vector index; 0 = auto
	// PME-specific
	int pme_order = 4;                    // B-spline order (4 or 6)
	std::array<int, 3> pme_grid = {32, 32, 32};
};

struct RunConfig {
	std::uint64_t num_steps = 0;
	double time_step = 0.0;            // internal integration time unit
	double time_step_fs = 0.0;         // user-facing time step in femtoseconds
	double temperature = 0.0;
	std::string velocity_init_mode = "random";
	std::uint32_t velocity_seed = 5489u;
	bool remove_center_of_mass_velocity = true;

	// Inline force field parameters parsed directly from the run input file.
	// When present, no separate .ff file is required.
	std::optional<LJForceFieldConfig> force_field;

	// Long-range Coulomb method.  Absent means no electrostatics.
	std::optional<CoulombConfig> coulomb;

	// Thermostat.  Empty string means NVE (no temperature coupling).
	// Supported values: "velocity_rescaling", "nose_hoover".
	std::string thermostat_type = "";
	double      thermostat_tau  = 100.0;   // Nosé-Hoover coupling time [fs]

	// Barostat.  Empty string means no pressure coupling.
	// Supported values: "berendsen", "monte_carlo".
	std::string   barostat_type    = "";
	double        target_pressure  = 1.0;    // target pressure [bar]
	double        barostat_tau     = 2000.0; // pressure bath relaxation time [fs] (Berendsen)
	double        compressibility  = 4.5e-5; // isothermal compressibility [1/bar]  (Berendsen)
	std::uint32_t mc_frequency     = 25;     // attempt volume move every N steps   (MC)
	double        mc_volume_step   = 0.01;   // initial max |Δ ln V| for trial moves (MC)
};

// ---------------------------------------------------------------------------
// Explicit LJ cross-pair override (overrides Lorentz-Berthelot mixing).
// ---------------------------------------------------------------------------
struct ExplicitPairOverride {
	double epsilon = 0.0;  // well depth [eV]
	double sigma   = 0.0;  // zero-crossing distance [Å]
};

// ---------------------------------------------------------------------------
// Full molecular force field:
//   - LJ nonbonded section  (atom types, cutoff, pair overrides)
//   - bonded parameter tables  (bonds, angles, dihedrals, impropers)
//
// Designed to be loaded from a .ff file by ConfigLoader::load_molecular_ff().
// All internal energy/length units follow the project convention
// (eV / Å / amu).  The parser converts kcal/mol → eV and deg → rad
// automatically, so .ff files may use the conventional force-field units.
// ---------------------------------------------------------------------------
struct MolecularForceFieldConfig {
	// ---- nonbonded ----
	LJForceFieldConfig lj;   // atom types (mass, ε, σ, q) + cutoff

	// Explicit cross-pair overrides keyed by (min_type, max_type), 0-based.
	// When present for a pair (i,j) these values replace Lorentz-Berthelot.
	std::map<std::pair<int,int>, ExplicitPairOverride> pair_overrides;

	// ---- bonded parameter tables ----
	// Indices in these vectors correspond to the type_idx stored in Topology terms.
	std::vector<BondParams>     bond_types;
	std::vector<AngleParams>    angle_types;
	std::vector<DihedralParams> dihedral_types;
	std::vector<ImproperParams> improper_types;
};

class ConfigLoader {
public:
	ConfigLoader() = default;

	// Loads atom positions from an XYZ-format file into `system`.
	//
	// Atom lines accept two first-column formats:
	//
	//   Integer type (1-based, LAMMPS style):
	//     1  x  y  z              (4 or 7 columns — mass from ff->elements[0])
	//     1  x  y  z  mass        (5 or 8 columns — explicit mass)
	//
	//   Element symbol (when each symbol maps to a unique type):
	//     Ar  x  y  z             (4 or 7 columns — mass from built-in table)
	//     Ar  x  y  z  mass       (5 or 8 columns — explicit mass overrides table)
	//
	// When `ff` is non-null and integer type IDs are used, mass is looked up from
	// ff->elements[type-1].mass (which was populated from the built-in element table
	// at FF parse time).  Element symbols in xyz are resolved via
	// ff->element_name_to_type.  If `ff` is null both formats still work but
	// integer-type atoms must carry an explicit mass column.
	void load_xyz(const std::filesystem::path& xyz_path, System& system,
	              const LJForceFieldConfig* ff = nullptr) const;

	RunConfig load_run(const std::filesystem::path& run_path) const;

	// Reads a force field parameter file and returns an LJForceFieldConfig.
	// Throws std::runtime_error on parse errors or unsupported force_field type.
	LJForceFieldConfig load_force_field(const std::filesystem::path& ff_path) const;

	// Reads a molecular force field file (.ff) and returns a
	// MolecularForceFieldConfig.  The file supports:
	//   - Per-type mass, LJ epsilon/sigma/charge via  type  directives
	//   - LJ cutoff via  lj_cutoff
	//   - Explicit cross-pair overrides via  pair  directives
	//   - Bond, angle, dihedral and improper parameter tables
	// Throws std::runtime_error on parse errors.
	MolecularForceFieldConfig load_molecular_ff(const std::filesystem::path& ff_path) const;

	// Reads a topology file (.top) and returns a Topology.
	// Each section begins with its keyword followed by an entry count, then one
	// entry per line.  Atom indices in the file are 1-based and stored 0-based in
	// the returned Topology.  Example:
	//
	//   bonds 2
	//   1 2  bond_type 1
	//   2 3  bond_type 1
	//
	//   angles 1
	//   1 2 3  angle_type 1
	//
	//   dihedrals 0
	//
	//   impropers 0
	//
	// Throws std::runtime_error on parse errors.
	std::shared_ptr<Topology> load_topology(const std::filesystem::path& top_path) const;
};

}  // namespace gmd
