#include "gmd/io/config_loader.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gmd/system/topology.hpp"

#include "gmd/system/system.hpp"

namespace gmd {

namespace {

// The integrator uses a reduced internal time unit derived from the code's
// eV / Angstrom / amu convention. Keep the conversion explicit so input/output
// can continue to use femtoseconds.
constexpr double kInternalTimeUnitsPerFs = 1.018051e+1;

// Standard atomic masses [amu] for common elements.
// Used when xyz input provides element symbols instead of explicit masses.
const std::unordered_map<std::string, double> kElementMasses = {
    {"H",  1.008},  {"He", 4.003},  {"Li", 6.941},  {"Be", 9.012},
    {"B",  10.811}, {"C",  12.011}, {"N",  14.007},  {"O",  15.999},
    {"F",  18.998}, {"Ne", 20.180}, {"Na", 22.990},  {"Mg", 24.305},
    {"Al", 26.982}, {"Si", 28.086}, {"P",  30.974},  {"S",  32.065},
    {"Cl", 35.453}, {"Ar", 39.948}, {"K",  39.098},  {"Ca", 40.078},
    {"Ti", 47.867}, {"Fe", 55.845}, {"Ni", 58.693},  {"Cu", 63.546},
    {"Zn", 65.38},  {"Kr", 83.798}, {"Ag", 107.868}, {"Sn", 118.710},
    {"Xe", 131.29}, {"Au", 196.967},{"Pb", 207.2},
};

// Atomic numbers (Z) for common elements.
const std::unordered_map<std::string, int> kElementAtomicNumbers = {
    {"H",  1},  {"He", 2},  {"Li", 3},  {"Be", 4},  {"B",  5},
    {"C",  6},  {"N",  7},  {"O",  8},  {"F",  9},  {"Ne", 10},
    {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P",  15},
    {"S",  16}, {"Cl", 17}, {"Ar", 18}, {"K",  19}, {"Ca", 20},
    {"Ti", 22}, {"Fe", 26}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
    {"Kr", 36}, {"Ag", 47}, {"Sn", 50}, {"Xe", 54}, {"Au", 79},
    {"Pb", 82},
};

// Returns true if the string looks like an element symbol (starts with a letter).
bool is_element_symbol(const std::string& s) {
    return !s.empty() && std::isalpha(static_cast<unsigned char>(s[0]));
}

std::string to_lower_ascii(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

std::string parse_mixing_rule(const std::string& token) {
    const std::string rule = to_lower_ascii(token);
    if (rule == "lb" || rule == "lorentz_berthelot" || rule == "lorentz-berthelot") {
        return "lorentz_berthelot";
    }
    if (rule == "geometric" || rule == "geom") {
        return "geometric";
    }
    if (rule == "waldman_hagler" || rule == "waldman-hagler" || rule == "wh") {
        return "waldman_hagler";
    }
    throw std::runtime_error(
        "Unsupported mixing_rule: " + token
        + " (supported: lorentz_berthelot|lb, geometric|geom, waldman_hagler|wh)");
}

std::vector<std::string> tokenize_line(std::istream& input) {
    std::string line;
    std::getline(input, line);

    std::istringstream line_stream(line);
    std::vector<std::string> tokens;
    std::string token;
    while (line_stream >> token) {
        // Stop at inline comment.
        if (token.starts_with('#')) break;
        tokens.push_back(token);
    }
    return tokens;
}

int parse_int(const std::string& token, const char* field_name) {
    try {
        return std::stoi(token);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Failed to parse integer field: ") + field_name);
    }
}

std::uint32_t parse_uint32(const std::string& token, const char* field_name) {
    const auto value = parse_int(token, field_name);
    if (value < 0) {
        throw std::runtime_error(std::string("Expected non-negative integer for field: ") + field_name);
    }
    return static_cast<std::uint32_t>(value);
}

bool parse_bool(const std::string& token, const char* field_name) {
    if (token == "1" || token == "true" || token == "on" || token == "yes") {
        return true;
    }
    if (token == "0" || token == "false" || token == "off" || token == "no") {
        return false;
    }
    throw std::runtime_error(std::string("Failed to parse boolean field: ") + field_name);
}

double parse_double(const std::string& token, const char* field_name) {
    try {
        return std::stod(token);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Failed to parse floating-point field: ") + field_name);
    }
}

// Parses a "type <N> <elem> epsilon <val> sigma <val>" directive into `ff`.
// `tokens` must be the full tokenized line starting with "type".
// N is 1-based; the entry is stored at index N-1 (vector is resized if needed).
void parse_type_into_ff(const std::vector<std::string>& tokens, LJForceFieldConfig& ff) {
    // type <N> <elem> epsilon <val> sigma <val>
    if (tokens.size() < 7) {
        throw std::runtime_error(
            "type directive expects: type <N> <element> epsilon <val> sigma <val>");
    }
    const int type_num = parse_int(tokens[1], "type_number");
    if (type_num < 1) {
        throw std::runtime_error("type number must be >= 1, got " + tokens[1]);
    }
    const int idx = type_num - 1;
    const std::string& elem_name = tokens[2];

    LJElementParams ep;
    ep.element = elem_name;
    // Look up mass from the built-in table.
    auto mit = kElementMasses.find(elem_name);
    if (mit != kElementMasses.end()) {
        ep.mass = mit->second;
    }
    // else mass stays 0.0; load_xyz will require explicit mass column in that case.

    bool got_eps = false, got_sig = false;
    for (std::size_t k = 3; k + 1 < tokens.size(); k += 2) {
        if (tokens[k] == "epsilon") {
            ep.epsilon = parse_double(tokens[k + 1], "epsilon");
            if (ep.epsilon <= 0.0) throw std::runtime_error("type epsilon must be positive");
            got_eps = true;
        } else if (tokens[k] == "sigma") {
            ep.sigma = parse_double(tokens[k + 1], "sigma");
            if (ep.sigma <= 0.0) throw std::runtime_error("type sigma must be positive");
            got_sig = true;
        } else if (tokens[k] == "charge") {
            ep.charge = parse_double(tokens[k + 1], "charge");
        }
    }
    if (!got_eps || !got_sig) {
        throw std::runtime_error(
            "type " + std::to_string(type_num) + " must specify both epsilon and sigma");
    }

    // Resize to accommodate out-of-order type numbers.
    if (static_cast<int>(ff.elements.size()) <= idx) {
        ff.elements.resize(static_cast<std::size_t>(idx + 1));
    }
    ff.elements[static_cast<std::size_t>(idx)] = ep;
    ff.element_name_to_type[elem_name] = idx;  // last assignment wins
}

}  // namespace

void ConfigLoader::load_xyz(const std::filesystem::path& xyz_path, System& system,
                            const LJForceFieldConfig* ff) const {
    std::ifstream input(xyz_path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open xyz input file: " + xyz_path.string());
    }

    auto tokens = tokenize_line(input);
    if (tokens.size() != 1) {
        throw std::runtime_error("The first line of xyz input must contain exactly one atom count value");
    }

    const auto parsed_atom_count = parse_int(tokens[0], "atom_count");
    if (parsed_atom_count < 0) {
        throw std::runtime_error("Atom count must be non-negative");
    }
    const auto atom_count = static_cast<std::size_t>(parsed_atom_count);
    system.resize(atom_count);

    tokens = tokenize_line(input);
    if (tokens.size() != 3) {
        throw std::runtime_error("The second line of xyz input must contain exactly three box lengths");
    }

    Box box;
    const auto box_x = parse_double(tokens[0], "box_x");
    const auto box_y = parse_double(tokens[1], "box_y");
    const auto box_z = parse_double(tokens[2], "box_z");
    if (box_x <= 0.0 || box_y <= 0.0 || box_z <= 0.0) {
        throw std::runtime_error("Box lengths must be positive");
    }
    box.set_lengths({box_x, box_y, box_z});
    system.set_box(box);

    auto coordinates = system.mutable_coordinates();
    auto masses      = system.mutable_masses();
    auto charges     = system.mutable_charges();
    auto atom_types  = system.mutable_atom_types();
    // Tracks which atoms received an explicit charge column in the XYZ file.
    // These will NOT be overwritten by the FF-type default charge pass below.
    std::vector<bool> has_explicit_charge(atom_count, false);

    // Local element→type map built on first-appearance order when ff is null.
    std::unordered_map<std::string, int> local_element_map;
    int next_type_id = 0;

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        tokens = tokenize_line(input);

        if (tokens.empty()) {
            throw std::runtime_error(
                "Unexpected end of xyz input at atom line " + std::to_string(atom_index + 1));
        }

        if (is_element_symbol(tokens[0])) {
            // ---- Element-symbol mode ----
            // Accepted column counts:
            //   4: elem x y z                      (mass from ff or built-in table)
            //   5: elem x y z mass                 (explicit mass)
            //   6: elem x y z mass charge          (explicit mass + charge)
            //   7: elem x y z vx vy vz             (mass from ff or built-in table)
            //   8: elem x y z vx vy vz mass
            //   9: elem x y z vx vy vz mass charge (explicit mass + charge)
            if (tokens.size() != 4 && tokens.size() != 5 && tokens.size() != 6 &&
                tokens.size() != 7 && tokens.size() != 8 && tokens.size() != 9) {
                throw std::runtime_error(
                    "Atom line " + std::to_string(atom_index + 1)
                    + " (element-symbol format) must have 4, 5, 6, 7, 8, or 9 columns, got "
                    + std::to_string(tokens.size()));
            }

            const std::string& elem = tokens[0];

            // Resolve 0-based type index.
            if (ff != nullptr) {
                auto it = ff->element_name_to_type.find(elem);
                if (it == ff->element_name_to_type.end()) {
                    throw std::runtime_error(
                        "Element '" + elem + "' in xyz input is not defined in the force field");
                }
                atom_types[atom_index] = it->second;
            } else {
                auto [it, inserted] = local_element_map.emplace(elem, next_type_id);
                if (inserted) { ++next_type_id; }
                atom_types[atom_index] = it->second;
            }

            coordinates[atom_index] = {parse_double(tokens[1], "x"),
                                       parse_double(tokens[2], "y"),
                                       parse_double(tokens[3], "z")};

            auto velocities = system.mutable_velocities();
            if (tokens.size() == 7 || tokens.size() == 8 || tokens.size() == 9) {
                velocities[atom_index] = {parse_double(tokens[4], "vx"),
                                          parse_double(tokens[5], "vy"),
                                          parse_double(tokens[6], "vz")};
            }

            // Mass: explicit column → ff type entry → built-in table.
            // Layouts with explicit mass: 5,6 (no vel) and 8,9 (with vel).
            if (tokens.size() == 5 || tokens.size() == 6) {
                masses[atom_index] = parse_double(tokens[4], "mass");
            } else if (tokens.size() == 8 || tokens.size() == 9) {
                masses[atom_index] = parse_double(tokens[7], "mass");
            } else if (ff != nullptr) {
                const int tidx = atom_types[atom_index];
                const double m = ff->elements[static_cast<std::size_t>(tidx)].mass;
                if (m <= 0.0) {
                    throw std::runtime_error(
                        "No mass defined for element '" + elem
                        + "' in the force field. Provide mass explicitly as the last column.");
                }
                masses[atom_index] = m;
            } else {
                auto mit = kElementMasses.find(elem);
                if (mit == kElementMasses.end()) {
                    throw std::runtime_error(
                        "No built-in mass for element '" + elem
                        + "'. Provide mass explicitly as the last column.");
                }
                masses[atom_index] = mit->second;
            }

            // Charge: explicit last column overrides FF-type default.
            // Layouts with explicit charge: 6 (no vel) and 9 (with vel).
            if (tokens.size() == 6) {
                charges[atom_index] = parse_double(tokens[5], "charge");
                has_explicit_charge[atom_index] = true;
            } else if (tokens.size() == 9) {
                charges[atom_index] = parse_double(tokens[8], "charge");
                has_explicit_charge[atom_index] = true;
            }
        } else {
            // ---- Integer type mode (1-based, LAMMPS style) ----
            // When ff is provided mass can be omitted (4 or 7 columns).
            // Without ff, explicit mass is required (5, 6, 8, or 9 columns).
            // Column layouts:
            //   4: type x y z                      (mass from ff)
            //   5: type x y z mass
            //   6: type x y z mass charge
            //   7: type x y z vx vy vz             (mass from ff)
            //   8: type x y z vx vy vz mass
            //   9: type x y z vx vy vz mass charge
            const bool has_ff = (ff != nullptr);
            if (has_ff) {
                if (tokens.size() != 4 && tokens.size() != 5 && tokens.size() != 6 &&
                    tokens.size() != 7 && tokens.size() != 8 && tokens.size() != 9) {
                    throw std::runtime_error(
                        "Atom line " + std::to_string(atom_index + 1)
                        + " must have 4, 5, 6, 7, 8, or 9 columns, got "
                        + std::to_string(tokens.size()));
                }
            } else {
                if (tokens.size() != 5 && tokens.size() != 6 &&
                    tokens.size() != 8 && tokens.size() != 9) {
                    throw std::runtime_error(
                        "Atom line " + std::to_string(atom_index + 1)
                        + " must have 5/6 columns (type x y z mass [charge]) or 8/9 columns"
                          " (type x y z vx vy vz mass [charge]), got "
                        + std::to_string(tokens.size()));
                }
            }

            const int type_num = parse_int(tokens[0], "atom_type");
            if (type_num < 1) {
                throw std::runtime_error(
                    "Atom type on line " + std::to_string(atom_index + 1)
                    + " must be >= 1 (got " + std::to_string(type_num) + ")");
            }
            const int type_idx = type_num - 1;  // convert to 0-based
            atom_types[atom_index] = type_idx;

            coordinates[atom_index] = {parse_double(tokens[1], "x"),
                                       parse_double(tokens[2], "y"),
                                       parse_double(tokens[3], "z")};
            auto velocities = system.mutable_velocities();
            if (tokens.size() == 7 || tokens.size() == 8 || tokens.size() == 9) {
                velocities[atom_index] = {parse_double(tokens[4], "vx"),
                                          parse_double(tokens[5], "vy"),
                                          parse_double(tokens[6], "vz")};
            }

            // Mass: explicit column or from ff type definition.
            if (tokens.size() == 5 || tokens.size() == 6) {
                masses[atom_index] = parse_double(tokens[4], "mass");
            } else if (tokens.size() == 8 || tokens.size() == 9) {
                masses[atom_index] = parse_double(tokens[7], "mass");
            } else {
                // 4 or 7 columns: mass comes from ff.
                if (type_idx >= static_cast<int>(ff->elements.size())) {
                    throw std::runtime_error(
                        "Type " + std::to_string(type_num)
                        + " on atom line " + std::to_string(atom_index + 1)
                        + " is not defined in the force field");
                }
                const double m = ff->elements[static_cast<std::size_t>(type_idx)].mass;
                if (m <= 0.0) {
                    throw std::runtime_error(
                        "No mass for type " + std::to_string(type_num)
                        + " in the force field. Provide mass explicitly as the last column.");
                }
                masses[atom_index] = m;
            }

            // Charge: explicit last column overrides FF-type default.
            if (tokens.size() == 6) {
                charges[atom_index] = parse_double(tokens[5], "charge");
                has_explicit_charge[atom_index] = true;
            } else if (tokens.size() == 9) {
                charges[atom_index] = parse_double(tokens[8], "charge");
                has_explicit_charge[atom_index] = true;
            }
        }

        if (masses[atom_index] <= 0.0) {
            throw std::runtime_error(
                "Mass of atom " + std::to_string(atom_index + 1) + " must be positive");
        }
    }

    // Populate per-atom charges from FF type definitions for atoms that did NOT
    // receive an explicit charge column in the XYZ file.
    if (ff != nullptr) {
        const auto types = system.atom_types();
        for (std::size_t i = 0; i < atom_count; ++i) {
            if (has_explicit_charge[i]) continue;  // explicit column takes precedence
            const int tidx = types[i];
            if (tidx >= 0 && static_cast<std::size_t>(tidx) < ff->elements.size()) {
                charges[i] = ff->elements[static_cast<std::size_t>(tidx)].charge;
            }
        }
    }

    // Populate per-atom atomic numbers from element symbols.
    // For element-symbol mode the element comes from the xyz line (stored in
    // local_element_map or ff->element_name_to_type).  For integer-type mode
    // the element symbol is taken from ff->elements[type_idx].element.
    {
        auto atomic_numbers = system.mutable_atomic_numbers();
        const auto types    = system.atom_types();

        // Build a per-type atomic number lookup from the ff (if available).
        std::unordered_map<int, int> type_to_z;
        if (ff != nullptr) {
            for (auto& [sym, tidx] : ff->element_name_to_type) {
                auto zit = kElementAtomicNumbers.find(sym);
                if (zit != kElementAtomicNumbers.end()) {
                    type_to_z[tidx] = zit->second;
                }
            }
        }
        // Also include entries built from the local element map (no-ff element-symbol mode).
        for (auto& [sym, tidx] : local_element_map) {
            auto zit = kElementAtomicNumbers.find(sym);
            if (zit != kElementAtomicNumbers.end()) {
                type_to_z[tidx] = zit->second;
            }
        }

        for (std::size_t i = 0; i < atom_count; ++i) {
            const int tidx = types[i];
            auto it = type_to_z.find(tidx);
            atomic_numbers[i] = (it != type_to_z.end()) ? it->second : 0;
        }
    }
}

RunConfig ConfigLoader::load_run(const std::filesystem::path& run_path) const {
    std::ifstream input(run_path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open run input file: " + run_path.string());
    }

    RunConfig config;

    // Inline force field state (populated when FF directives appear in the run file).
    bool has_inline_ff    = false;
    double single_epsilon = 0.0;
    double single_sigma   = 0.0;
    bool   has_single_eps = false;
    bool   has_single_sig = false;
    LJForceFieldConfig inline_ff;

    while (input.peek() != EOF) {
        const auto tokens = tokenize_line(input);
        if (tokens.empty() || tokens[0].starts_with('#')) {
            continue;
        }

        if (tokens.size() < 2) {
            throw std::runtime_error("Run input directive is missing a value: " + tokens[0]);
        }

        if (tokens[0] == "time_step") {
            config.time_step_fs = parse_double(tokens[1], "time_step");
            config.time_step = config.time_step_fs / kInternalTimeUnitsPerFs;
            if (config.time_step_fs < 0.0) {
                throw std::runtime_error("time_step must be non-negative");
            }
        } else if (tokens[0] == "run") {
            config.num_steps = static_cast<std::uint64_t>(parse_int(tokens[1], "run"));
            if (config.num_steps == 0) {
                throw std::runtime_error("run must be greater than zero");
            }
        } else if (tokens[0] == "velocity") {
            config.temperature = parse_double(tokens[1], "velocity");
            if (config.temperature < 0.0) {
                throw std::runtime_error("velocity must be non-negative");
            }
        } else if (tokens[0] == "velocity_init") {
            if (tokens[1] != "random" && tokens[1] != "input") {
                throw std::runtime_error("velocity_init must be either 'random' or 'input'");
            }
            config.velocity_init_mode = tokens[1];
        } else if (tokens[0] == "velocity_seed") {
            config.velocity_seed = parse_uint32(tokens[1], "velocity_seed");
        } else if (tokens[0] == "remove_com_velocity") {
            config.remove_center_of_mass_velocity = parse_bool(tokens[1], "remove_com_velocity");
        } else if (tokens[0] == "molecular_nonbonded") {
            if (tokens[1] != "none" && tokens[1] != "lj_unsafe") {
                throw std::runtime_error(
                    "molecular_nonbonded must be 'none' or 'lj_unsafe'");
            }
            config.molecular_nonbonded_mode = tokens[1];
        // ---- Inline force field directives ----
        } else if (tokens[0] == "force_field") {
            if (tokens[1] == "lj") {
                has_inline_ff = true;
                config.force_field_type = "lj";
            } else if (tokens[1] == "ml") {
                config.force_field_type = "ml";
            } else {
                throw std::runtime_error("Unsupported force_field type: " + tokens[1]
                                         + " (supported: lj, ml)");
            }
        } else if (tokens[0] == "model_path") {
            config.ml_model_path = tokens[1];
        } else if (tokens[0] == "cutoff") {
            inline_ff.cutoff = parse_double(tokens[1], "cutoff");
            if (inline_ff.cutoff <= 0.0) throw std::runtime_error("cutoff must be positive");
        } else if (tokens[0] == "mixing_rule" || tokens[0] == "mixing") {
            inline_ff.mixing_rule = parse_mixing_rule(tokens[1]);
        } else if (tokens[0] == "epsilon") {
            single_epsilon = parse_double(tokens[1], "epsilon");
            if (single_epsilon <= 0.0) throw std::runtime_error("epsilon must be positive");
            has_single_eps = true;
        } else if (tokens[0] == "sigma") {
            single_sigma = parse_double(tokens[1], "sigma");
            if (single_sigma <= 0.0) throw std::runtime_error("sigma must be positive");
            has_single_sig = true;
        } else if (tokens[0] == "element") {
            parse_type_into_ff(tokens, inline_ff);
        } else if (tokens[0] == "type") {
            parse_type_into_ff(tokens, inline_ff);
        // ---- Long-range Coulomb directives ----
        } else if (tokens[0] == "coulomb") {
            if (tokens[1] != "ewald" && tokens[1] != "pme") {
                throw std::runtime_error("coulomb method must be 'ewald' or 'pme'");
            }
            if (!config.coulomb.has_value()) config.coulomb = CoulombConfig{};
            config.coulomb->method = tokens[1];
        } else if (tokens[0] == "ewald_alpha") {
            if (!config.coulomb.has_value()) config.coulomb = CoulombConfig{};
            config.coulomb->alpha = parse_double(tokens[1], "ewald_alpha");
        } else if (tokens[0] == "ewald_kmax") {
            if (!config.coulomb.has_value()) config.coulomb = CoulombConfig{};
            config.coulomb->kmax = parse_int(tokens[1], "ewald_kmax");
        } else if (tokens[0] == "ewald_cutoff") {
            if (!config.coulomb.has_value()) config.coulomb = CoulombConfig{};
            config.coulomb->real_cutoff = parse_double(tokens[1], "ewald_cutoff");
        } else if (tokens[0] == "pme_alpha") {
            if (!config.coulomb.has_value()) config.coulomb = CoulombConfig{};
            config.coulomb->alpha = parse_double(tokens[1], "pme_alpha");
        } else if (tokens[0] == "pme_cutoff") {
            if (!config.coulomb.has_value()) config.coulomb = CoulombConfig{};
            config.coulomb->real_cutoff = parse_double(tokens[1], "pme_cutoff");
        } else if (tokens[0] == "pme_order") {
            if (!config.coulomb.has_value()) config.coulomb = CoulombConfig{};
            const int ord = parse_int(tokens[1], "pme_order");
            if (ord != 4 && ord != 6) throw std::runtime_error("pme_order must be 4 or 6");
            config.coulomb->pme_order = ord;
        } else if (tokens[0] == "pme_grid") {
            if (tokens.size() < 4) {
                throw std::runtime_error("pme_grid expects three integers: pme_grid Nx Ny Nz");
            }
            if (!config.coulomb.has_value()) config.coulomb = CoulombConfig{};
            config.coulomb->pme_grid[0] = parse_int(tokens[1], "pme_grid_x");
            config.coulomb->pme_grid[1] = parse_int(tokens[2], "pme_grid_y");
            config.coulomb->pme_grid[2] = parse_int(tokens[3], "pme_grid_z");
        // ---- Thermostat directives ----
        } else if (tokens[0] == "thermostat") {
            if (tokens[1] != "velocity_rescaling" && tokens[1] != "nose_hoover") {
                throw std::runtime_error(
                    "thermostat must be 'velocity_rescaling' or 'nose_hoover'");
            }
            config.thermostat_type = tokens[1];
        } else if (tokens[0] == "thermostat_tau") {
            config.thermostat_tau = parse_double(tokens[1], "thermostat_tau");
            if (config.thermostat_tau <= 0.0)
                throw std::runtime_error("thermostat_tau must be positive");
        // ---- Barostat directives ----
        } else if (tokens[0] == "barostat") {
            if (tokens[1] != "berendsen" && tokens[1] != "monte_carlo") {
                throw std::runtime_error("barostat must be 'berendsen' or 'monte_carlo'");
            }
            config.barostat_type = tokens[1];
        } else if (tokens[0] == "pressure") {
            config.target_pressure = parse_double(tokens[1], "pressure");
        } else if (tokens[0] == "barostat_tau") {
            config.barostat_tau = parse_double(tokens[1], "barostat_tau");
            if (config.barostat_tau <= 0.0)
                throw std::runtime_error("barostat_tau must be positive");
        } else if (tokens[0] == "compressibility") {
            config.compressibility = parse_double(tokens[1], "compressibility");
            if (config.compressibility <= 0.0)
                throw std::runtime_error("compressibility must be positive");
        } else if (tokens[0] == "mc_frequency") {
            const double v = parse_double(tokens[1], "mc_frequency");
            if (v < 1.0) throw std::runtime_error("mc_frequency must be >= 1");
            config.mc_frequency = static_cast<std::uint32_t>(v);
        } else if (tokens[0] == "mc_volume_step") {
            config.mc_volume_step = parse_double(tokens[1], "mc_volume_step");
            if (config.mc_volume_step <= 0.0)
                throw std::runtime_error("mc_volume_step must be positive");
        } else {
            throw std::runtime_error("Unsupported run input directive: " + tokens[0]);
        }
    }

    if (has_inline_ff) {
        // Promote legacy single-element epsilon/sigma if no type directives were given.
        if (inline_ff.elements.empty()) {
            if (!has_single_eps || !has_single_sig) {
                throw std::runtime_error(
                    "Inline LJ force field must specify 'type' directives "
                    "or both 'epsilon' and 'sigma'");
            }
            LJElementParams ep;
            ep.epsilon = single_epsilon;
            ep.sigma   = single_sigma;
            inline_ff.elements.push_back(ep);
            inline_ff.element_name_to_type["default"] = 0;
        }
        config.force_field = std::move(inline_ff);
    }

    return config;
}

void LJForceFieldConfig::pair_params(int type_i, int type_j,
                                     double& eps_out, double& sig_out) const noexcept {
    // Check explicit cross-pair overrides first.
    const auto key = std::make_pair(std::min(type_i, type_j), std::max(type_i, type_j));
    const auto it  = pair_overrides.find(key);
    if (it != pair_overrides.end()) {
        eps_out = it->second.epsilon;
        sig_out = it->second.sigma;
        return;
    }
    // Fall back to configurable mixing rules.
    const auto& ei = elements[static_cast<std::size_t>(type_i)];
    const auto& ej = elements[static_cast<std::size_t>(type_j)];

    if (mixing_rule == "geometric") {
        eps_out = std::sqrt(ei.epsilon * ej.epsilon);
        sig_out = std::sqrt(ei.sigma * ej.sigma);
        return;
    }
    if (mixing_rule == "waldman_hagler") {
        const double si6 = std::pow(ei.sigma, 6.0);
        const double sj6 = std::pow(ej.sigma, 6.0);
        const double denom = si6 + sj6;
        sig_out = std::pow(0.5 * denom, 1.0 / 6.0);
        eps_out = (2.0 * std::sqrt(ei.epsilon * ej.epsilon) * std::pow(ei.sigma, 3.0)
                   * std::pow(ej.sigma, 3.0)) / denom;
        return;
    }

    // Default: Lorentz-Berthelot.
    eps_out = std::sqrt(ei.epsilon * ej.epsilon);  // geometric mean
    sig_out = 0.5 * (ei.sigma + ej.sigma);         // arithmetic mean
}

LJForceFieldConfig ConfigLoader::load_force_field(const std::filesystem::path& ff_path) const {
    std::ifstream input(ff_path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open force field file: " + ff_path.string());
    }

    LJForceFieldConfig config;
    bool has_force_field_type = false;
    // Temporary single-element epsilon/sigma (legacy format).
    double single_epsilon = 0.0;
    double single_sigma   = 0.0;
    bool   has_single_eps = false;
    bool   has_single_sig = false;

    while (input.peek() != EOF) {
        const auto tokens = tokenize_line(input);
        if (tokens.empty() || tokens[0].starts_with('#')) {
            continue;
        }
        if (tokens.size() < 2) {
            throw std::runtime_error("Force field directive missing value: " + tokens[0]);
        }

        if (tokens[0] == "force_field") {
            if (tokens[1] != "lj") {
                throw std::runtime_error("Unsupported force_field type: " + tokens[1]
                                         + " (only 'lj' is supported)");
            }
            has_force_field_type = true;
        } else if (tokens[0] == "cutoff") {
            config.cutoff = parse_double(tokens[1], "cutoff");
            if (config.cutoff <= 0.0) {
                throw std::runtime_error("cutoff must be positive");
            }
        } else if (tokens[0] == "mixing_rule" || tokens[0] == "mixing") {
            config.mixing_rule = parse_mixing_rule(tokens[1]);
        } else if (tokens[0] == "epsilon") {
            // Legacy single-element format.
            single_epsilon = parse_double(tokens[1], "epsilon");
            if (single_epsilon <= 0.0) throw std::runtime_error("epsilon must be positive");
            has_single_eps = true;
        } else if (tokens[0] == "sigma") {
            single_sigma = parse_double(tokens[1], "sigma");
            if (single_sigma <= 0.0) throw std::runtime_error("sigma must be positive");
            has_single_sig = true;
        } else if (tokens[0] == "element") {
            // Backward compat: treat old "element Ar epsilon ... sigma ..." as type 1, 2, ...
            // Rewrite tokens to "type <N+1> ..." and delegate.
            std::vector<std::string> rewritten;
            rewritten.push_back("type");
            rewritten.push_back(std::to_string(static_cast<int>(config.elements.size()) + 1));
            rewritten.insert(rewritten.end(), tokens.begin() + 1, tokens.end());
            parse_type_into_ff(rewritten, config);
        } else if (tokens[0] == "type") {
            parse_type_into_ff(tokens, config);
        } else {
            throw std::runtime_error("Unsupported force field directive: " + tokens[0]);
        }
    }

    if (!has_force_field_type) {
        throw std::runtime_error("Force field file must specify 'force_field' type");
    }

    // If legacy single-element format was used, promote to elements vector.
    if (config.elements.empty()) {
        if (!has_single_eps || !has_single_sig) {
            throw std::runtime_error(
                "Single-element LJ force field must specify both 'epsilon' and 'sigma'");
        }
        LJElementParams ep;
        ep.epsilon = single_epsilon;
        ep.sigma   = single_sigma;
        config.elements.push_back(ep);
        config.element_name_to_type["default"] = 0;
    }

    return config;
}

// ============================================================================
// ConfigLoader::load_molecular_ff
// ============================================================================
//
// File format (.ff):
//
//   force_field  molecular
//   mixing_rule  lorentz_berthelot   # optional
//
//   # --- Atom types ---------------------------------------------------------
//   # type <N>  <elem>  mass <m>  epsilon <ε>  sigma <σ>  [charge <q>]
//   # Units: mass [amu], epsilon [kcal/mol], sigma [Å], charge [e]
//   type 1  CT  mass 12.011  epsilon 0.0860  sigma 3.400  charge -0.180
//   type 2  HC  mass  1.008  epsilon 0.0157  sigma 2.650  charge  0.060
//
//   lj_cutoff  12.0          # LJ cutoff [Å]
//
//   # --- Explicit pair overrides (optional) ---------------------------------
//   # pair <type_i> <type_j>  epsilon <ε>  sigma <σ>   (kcal/mol, Å)
//   pair 1 2  epsilon 0.0367  sigma 3.025
//
//   # --- Bond types ----------------------------------------------------------
//   # bond_type <N>  k <k>  r0 <r0>   (k: kcal/mol/Å², r0: Å)
//   bond_type 1  k 310.0  r0 1.526
//
//   # --- Angle types ---------------------------------------------------------
//   # angle_type <N>  k <k>  theta0 <θ_deg>  (k: kcal/mol/rad²)
//   angle_type 1  k 40.0  theta0 112.7
//
//   # --- Dihedral types ------------------------------------------------------
//   # dihedral_type <N>  k <k>  n <n>  delta <δ_deg>  (k: kcal/mol)
//   dihedral_type 1  k 1.400  n 3  delta 0.0
//
//   # --- Improper types ------------------------------------------------------
//   # improper_type <N>  k <k>  phi0 <φ_deg>  (k: kcal/mol/rad²)
//   improper_type 1  k 10.5  phi0 0.0
//
// ============================================================================

MolecularForceFieldConfig ConfigLoader::load_molecular_ff(
        const std::filesystem::path& ff_path) const {

    std::ifstream input(ff_path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open molecular FF file: " + ff_path.string());
    }

    MolecularForceFieldConfig cfg;
    bool has_ff_type = false;

    // kcal/mol → eV and deg → rad conversion constants (mirrored locally)
    constexpr double kcal_to_eV = 4.336410e-2;
    constexpr double d2r        = 3.14159265358979323846 / 180.0;

    while (input.peek() != EOF) {
        const auto tokens = tokenize_line(input);
        if (tokens.empty() || tokens[0].starts_with('#')) continue;

        if (tokens.size() < 2) {
            throw std::runtime_error("Molecular FF directive missing value: " + tokens[0]);
        }

        // ---- Top-level keyword ----
        if (tokens[0] == "force_field") {
            if (tokens[1] != "molecular") {
                throw std::runtime_error(
                    "load_molecular_ff: expected 'force_field molecular', got: " + tokens[1]);
            }
            has_ff_type = true;

        } else if (tokens[0] == "lj_cutoff") {
            cfg.lj.cutoff = parse_double(tokens[1], "lj_cutoff");
            if (cfg.lj.cutoff <= 0.0)
                throw std::runtime_error("lj_cutoff must be positive");
        } else if (tokens[0] == "mixing_rule" || tokens[0] == "mixing") {
            cfg.lj.mixing_rule = parse_mixing_rule(tokens[1]);

        // ---- Atom type ----
        } else if (tokens[0] == "type") {
            // type <N>  <elem>  mass <m>  epsilon <ε>  sigma <σ>  [charge <q>]
            if (tokens.size() < 8) {
                throw std::runtime_error(
                    "type directive: type <N> <elem> mass <m> epsilon <ε> sigma <σ> "
                    "[charge <q>]");
            }
            const int type_num = parse_int(tokens[1], "type_number");
            if (type_num < 1)
                throw std::runtime_error("type number must be >= 1");

            const int idx = type_num - 1;
            if (static_cast<int>(cfg.lj.elements.size()) <= idx)
                cfg.lj.elements.resize(static_cast<std::size_t>(idx + 1));

            LJElementParams& ep = cfg.lj.elements[static_cast<std::size_t>(idx)];
            ep.element = tokens[2];

            // Look up built-in mass default
            {
                auto mit = kElementMasses.find(ep.element);
                if (mit != kElementMasses.end()) ep.mass = mit->second;
            }

            // Parse key-value pairs after element symbol
            bool got_eps = false, got_sig = false;
            for (std::size_t k = 3; k + 1 < tokens.size(); k += 2) {
                const auto& key = tokens[k];
                const auto& val = tokens[k + 1];
                if      (key == "mass")    { ep.mass    = parse_double(val, "mass");
                                              if (ep.mass <= 0.0) throw std::runtime_error("mass must be positive"); }
                else if (key == "epsilon") { ep.epsilon = parse_double(val, "epsilon") * kcal_to_eV;
                                              if (ep.epsilon <= 0.0) throw std::runtime_error("epsilon must be positive");
                                              got_eps = true; }
                else if (key == "sigma")  { ep.sigma   = parse_double(val, "sigma");
                                              if (ep.sigma <= 0.0) throw std::runtime_error("sigma must be positive");
                                              got_sig = true; }
                else if (key == "charge") { ep.charge  = parse_double(val, "charge"); }
                else {
                    throw std::runtime_error("Unknown key in type directive: " + key);
                }
            }
            if (!got_eps || !got_sig)
                throw std::runtime_error(
                    "type " + std::to_string(type_num) + " must specify both epsilon and sigma");

            cfg.lj.element_name_to_type[ep.element] = idx;

        // ---- Explicit pair override ----
        } else if (tokens[0] == "pair") {
            // pair <type_i> <type_j>  epsilon <ε>  sigma <σ>
            if (tokens.size() < 6) {
                throw std::runtime_error(
                    "pair directive: pair <type_i> <type_j> epsilon <ε> sigma <σ>");
            }
            const int ti = parse_int(tokens[1], "pair_type_i") - 1;
            const int tj = parse_int(tokens[2], "pair_type_j") - 1;
            if (ti < 0 || tj < 0)
                throw std::runtime_error("pair type numbers must be >= 1");

            ExplicitPairOverride ov;
            bool got_eps = false, got_sig = false;
            for (std::size_t k = 3; k + 1 < tokens.size(); k += 2) {
                const auto& key = tokens[k];
                if      (key == "epsilon") { ov.epsilon = parse_double(tokens[k+1], "pair_epsilon") * kcal_to_eV; got_eps = true; }
                else if (key == "sigma")   { ov.sigma   = parse_double(tokens[k+1], "pair_sigma");                got_sig = true; }
                else throw std::runtime_error("Unknown key in pair directive: " + key);
            }
            if (!got_eps || !got_sig)
                throw std::runtime_error("pair directive must specify both epsilon and sigma");

            const auto key = std::make_pair(std::min(ti, tj), std::max(ti, tj));
            cfg.lj.pair_overrides[key] = ov;

        // ---- Bond types ----
        } else if (tokens[0] == "bond_type") {
            // bond_type <N>  k <k>  r0 <r0>
            if (tokens.size() < 6) {
                throw std::runtime_error("bond_type: bond_type <N> k <k> r0 <r0>");
            }
            const int num = parse_int(tokens[1], "bond_type_number");
            if (num < 1) throw std::runtime_error("bond_type number must be >= 1");
            const int idx = num - 1;
            if (static_cast<int>(cfg.bond_types.size()) <= idx)
                cfg.bond_types.resize(static_cast<std::size_t>(idx + 1));

            BondParams& bp = cfg.bond_types[static_cast<std::size_t>(idx)];
            bool got_k = false, got_r0 = false;
            for (std::size_t k = 2; k + 1 < tokens.size(); k += 2) {
                if      (tokens[k] == "k")  { bp.k  = parse_double(tokens[k+1], "bond_k") * kcal_to_eV; got_k  = true; }
                else if (tokens[k] == "r0") { bp.r0 = parse_double(tokens[k+1], "bond_r0");              got_r0 = true; }
                else throw std::runtime_error("Unknown key in bond_type: " + tokens[k]);
            }
            if (!got_k || !got_r0)
                throw std::runtime_error("bond_type " + std::to_string(num) + " must specify k and r0");

        // ---- Angle types ----
        } else if (tokens[0] == "angle_type") {
            // angle_type <N>  k <k>  theta0 <θ_deg>
            if (tokens.size() < 6) {
                throw std::runtime_error("angle_type: angle_type <N> k <k> theta0 <θ_deg>");
            }
            const int num = parse_int(tokens[1], "angle_type_number");
            if (num < 1) throw std::runtime_error("angle_type number must be >= 1");
            const int idx = num - 1;
            if (static_cast<int>(cfg.angle_types.size()) <= idx)
                cfg.angle_types.resize(static_cast<std::size_t>(idx + 1));

            AngleParams& ap = cfg.angle_types[static_cast<std::size_t>(idx)];
            bool got_k = false, got_t0 = false;
            for (std::size_t k = 2; k + 1 < tokens.size(); k += 2) {
                if      (tokens[k] == "k")      { ap.k      = parse_double(tokens[k+1], "angle_k") * kcal_to_eV; got_k  = true; }
                else if (tokens[k] == "theta0") { ap.theta0 = parse_double(tokens[k+1], "angle_theta0") * d2r;   got_t0 = true; }
                else throw std::runtime_error("Unknown key in angle_type: " + tokens[k]);
            }
            if (!got_k || !got_t0)
                throw std::runtime_error("angle_type " + std::to_string(num) + " must specify k and theta0");

        // ---- Dihedral types ----
        } else if (tokens[0] == "dihedral_type") {
            // dihedral_type <N>  k <k>  n <n>  delta <δ_deg>
            if (tokens.size() < 8) {
                throw std::runtime_error(
                    "dihedral_type: dihedral_type <N> k <k> n <n> delta <δ_deg>");
            }
            const int num = parse_int(tokens[1], "dihedral_type_number");
            if (num < 1) throw std::runtime_error("dihedral_type number must be >= 1");
            const int idx = num - 1;
            if (static_cast<int>(cfg.dihedral_types.size()) <= idx)
                cfg.dihedral_types.resize(static_cast<std::size_t>(idx + 1));

            DihedralParams& dp = cfg.dihedral_types[static_cast<std::size_t>(idx)];
            bool got_k = false, got_n = false, got_d = false;
            for (std::size_t k = 2; k + 1 < tokens.size(); k += 2) {
                if      (tokens[k] == "k")     { dp.k     = parse_double(tokens[k+1], "dihedral_k") * kcal_to_eV; got_k = true; }
                else if (tokens[k] == "n")     { dp.n     = parse_int(tokens[k+1], "dihedral_n");                 got_n = true; }
                else if (tokens[k] == "delta") { dp.delta = parse_double(tokens[k+1], "dihedral_delta") * d2r;    got_d = true; }
                else throw std::runtime_error("Unknown key in dihedral_type: " + tokens[k]);
            }
            if (!got_k || !got_n || !got_d)
                throw std::runtime_error(
                    "dihedral_type " + std::to_string(num) + " must specify k, n and delta");

        // ---- Improper types ----
        } else if (tokens[0] == "improper_type") {
            // improper_type <N>  k <k>  phi0 <φ_deg>
            if (tokens.size() < 6) {
                throw std::runtime_error(
                    "improper_type: improper_type <N> k <k> phi0 <φ_deg>");
            }
            const int num = parse_int(tokens[1], "improper_type_number");
            if (num < 1) throw std::runtime_error("improper_type number must be >= 1");
            const int idx = num - 1;
            if (static_cast<int>(cfg.improper_types.size()) <= idx)
                cfg.improper_types.resize(static_cast<std::size_t>(idx + 1));

            ImproperParams& ip = cfg.improper_types[static_cast<std::size_t>(idx)];
            bool got_k = false, got_p0 = false;
            for (std::size_t k = 2; k + 1 < tokens.size(); k += 2) {
                if      (tokens[k] == "k")    { ip.k    = parse_double(tokens[k+1], "improper_k") * kcal_to_eV; got_k  = true; }
                else if (tokens[k] == "phi0") { ip.phi0 = parse_double(tokens[k+1], "improper_phi0") * d2r;     got_p0 = true; }
                else throw std::runtime_error("Unknown key in improper_type: " + tokens[k]);
            }
            if (!got_k || !got_p0)
                throw std::runtime_error("improper_type " + std::to_string(num) + " must specify k and phi0");

        } else {
            throw std::runtime_error("Unknown directive in molecular FF file: " + tokens[0]);
        }
    }

    if (!has_ff_type)
        throw std::runtime_error("Molecular FF file must start with 'force_field molecular'");
    if (cfg.lj.elements.empty())
        throw std::runtime_error("Molecular FF file must define at least one 'type' entry");

    return cfg;
}

// ============================================================================
// ConfigLoader::load_topology
// ============================================================================
//
// File format (.top):
//
//   # Comment lines start with #.
//
//   bonds  <count>
//   <i> <j>  bond_type <N>
//   ...
//
//   angles  <count>
//   <i> <j> <k>  angle_type <N>
//   ...
//
//   dihedrals  <count>
//   <i> <j> <k> <l>  dihedral_type <N>
//   ...
//
//   impropers  <count>
//   <i> <j> <k> <l>  improper_type <N>
//   ...
//
// All atom indices are 1-based in the file; stored 0-based in Topology.
// Type numbers are 1-based in the file; stored 0-based in Topology terms.
// Sections may appear in any order and may be repeated (entries are appended).
// A count of 0 is allowed (section header with no following entries).
// ============================================================================

std::shared_ptr<Topology> ConfigLoader::load_topology(
        const std::filesystem::path& top_path) const {

    std::ifstream input(top_path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open topology file: " + top_path.string());
    }

    auto topo = std::make_shared<Topology>();

    while (input.peek() != EOF) {
        const auto tokens = tokenize_line(input);
        if (tokens.empty() || tokens[0].starts_with('#')) continue;

        // Every non-comment line is a section header: <keyword> <count>
        if (tokens.size() < 2) {
            throw std::runtime_error("Topology: expected '<section> <count>', got: " + tokens[0]);
        }

        const std::string& section = tokens[0];
        const int count = parse_int(tokens[1], "section_count");
        if (count < 0) throw std::runtime_error("Topology section count must be >= 0");

        if (section == "bonds") {
            for (int entry = 0; entry < count; ++entry) {
                const auto row = tokenize_line(input);
                // <i> <j>  bond_type <N>
                if (row.size() < 4 || row[2] != "bond_type")
                    throw std::runtime_error(
                        "bonds entry must be: <i> <j> bond_type <N>");
                BondTerm bt;
                bt.i        = parse_int(row[0], "bond_i") - 1;
                bt.j        = parse_int(row[1], "bond_j") - 1;
                bt.type_idx = parse_int(row[3], "bond_type") - 1;
                if (bt.i < 0 || bt.j < 0 || bt.type_idx < 0)
                    throw std::runtime_error("Bond atom indices and type must be >= 1");
                topo->bonds.push_back(bt);
            }
        } else if (section == "angles") {
            for (int entry = 0; entry < count; ++entry) {
                const auto row = tokenize_line(input);
                // <i> <j> <k>  angle_type <N>
                if (row.size() < 5 || row[3] != "angle_type")
                    throw std::runtime_error(
                        "angles entry must be: <i> <j> <k> angle_type <N>");
                AngleTerm at;
                at.i        = parse_int(row[0], "angle_i") - 1;
                at.j        = parse_int(row[1], "angle_j") - 1;
                at.k        = parse_int(row[2], "angle_k") - 1;
                at.type_idx = parse_int(row[4], "angle_type") - 1;
                if (at.i < 0 || at.j < 0 || at.k < 0 || at.type_idx < 0)
                    throw std::runtime_error("Angle atom indices and type must be >= 1");
                topo->angles.push_back(at);
            }
        } else if (section == "dihedrals") {
            for (int entry = 0; entry < count; ++entry) {
                const auto row = tokenize_line(input);
                // <i> <j> <k> <l>  dihedral_type <N>
                if (row.size() < 6 || row[4] != "dihedral_type")
                    throw std::runtime_error(
                        "dihedrals entry must be: <i> <j> <k> <l> dihedral_type <N>");
                DihedralTerm dt;
                dt.i        = parse_int(row[0], "dihedral_i") - 1;
                dt.j        = parse_int(row[1], "dihedral_j") - 1;
                dt.k        = parse_int(row[2], "dihedral_k") - 1;
                dt.l        = parse_int(row[3], "dihedral_l") - 1;
                dt.type_idx = parse_int(row[5], "dihedral_type") - 1;
                if (dt.i < 0 || dt.j < 0 || dt.k < 0 || dt.l < 0 || dt.type_idx < 0)
                    throw std::runtime_error("Dihedral atom indices and type must be >= 1");
                topo->dihedrals.push_back(dt);
            }
        } else if (section == "impropers") {
            for (int entry = 0; entry < count; ++entry) {
                const auto row = tokenize_line(input);
                // <i> <j> <k> <l>  improper_type <N>
                if (row.size() < 6 || row[4] != "improper_type")
                    throw std::runtime_error(
                        "impropers entry must be: <i> <j> <k> <l> improper_type <N>");
                ImproperTerm it;
                it.i        = parse_int(row[0], "improper_i") - 1;
                it.j        = parse_int(row[1], "improper_j") - 1;
                it.k        = parse_int(row[2], "improper_k") - 1;
                it.l        = parse_int(row[3], "improper_l") - 1;
                it.type_idx = parse_int(row[5], "improper_type") - 1;
                if (it.i < 0 || it.j < 0 || it.k < 0 || it.l < 0 || it.type_idx < 0)
                    throw std::runtime_error("Improper atom indices and type must be >= 1");
                topo->impropers.push_back(it);
            }
        } else {
            throw std::runtime_error("Unknown topology section: " + section);
        }
    }

    return topo;
}

}  // namespace gmd
