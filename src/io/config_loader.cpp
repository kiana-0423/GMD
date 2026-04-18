#include "gmd/io/config_loader.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "gmd/system/system.hpp"

namespace gmd {

namespace {

constexpr double kTimeUnitConversion = 1.018051e+1;

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

// Returns true if the string looks like an element symbol (starts with a letter).
bool is_element_symbol(const std::string& s) {
    return !s.empty() && std::isalpha(static_cast<unsigned char>(s[0]));
}

std::vector<std::string> tokenize_line(std::istream& input) {
    std::string line;
    std::getline(input, line);

    std::istringstream line_stream(line);
    std::vector<std::string> tokens;
    std::string token;
    while (line_stream >> token) {
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
    auto atom_types  = system.mutable_atom_types();

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
            //   4: elem x y z            (mass from ff or built-in table)
            //   5: elem x y z mass       (explicit mass)
            //   7: elem x y z vx vy vz   (mass from ff or built-in table)
            //   8: elem x y z vx vy vz mass
            if (tokens.size() != 4 && tokens.size() != 5 &&
                tokens.size() != 7 && tokens.size() != 8) {
                throw std::runtime_error(
                    "Atom line " + std::to_string(atom_index + 1)
                    + " (element-symbol format) must have 4, 5, 7, or 8 columns, got "
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
            if (tokens.size() == 7 || tokens.size() == 8) {
                velocities[atom_index] = {parse_double(tokens[4], "vx"),
                                          parse_double(tokens[5], "vy"),
                                          parse_double(tokens[6], "vz")};
            }

            // Mass: explicit last column → ff type entry → built-in table.
            if (tokens.size() == 5) {
                masses[atom_index] = parse_double(tokens[4], "mass");
            } else if (tokens.size() == 8) {
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
        } else {
            // ---- Integer type mode (1-based, LAMMPS style) ----
            // When ff is provided mass can be omitted (4 or 7 columns).
            // Without ff, explicit mass is required (5 or 8 columns).
            const bool has_ff = (ff != nullptr);
            if (has_ff) {
                if (tokens.size() != 4 && tokens.size() != 5 &&
                    tokens.size() != 7 && tokens.size() != 8) {
                    throw std::runtime_error(
                        "Atom line " + std::to_string(atom_index + 1)
                        + " must have 4, 5, 7, or 8 columns, got "
                        + std::to_string(tokens.size()));
                }
            } else {
                if (tokens.size() != 5 && tokens.size() != 8) {
                    throw std::runtime_error(
                        "Atom line " + std::to_string(atom_index + 1)
                        + " must have 5 columns (type x y z mass) or 8 columns"
                          " (type x y z vx vy vz mass), got "
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
            if (tokens.size() == 7 || tokens.size() == 8) {
                velocities[atom_index] = {parse_double(tokens[4], "vx"),
                                          parse_double(tokens[5], "vy"),
                                          parse_double(tokens[6], "vz")};
            }

            // Mass: explicit column or from ff type definition.
            if (tokens.size() == 5) {
                masses[atom_index] = parse_double(tokens[4], "mass");
            } else if (tokens.size() == 8) {
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
        }

        if (masses[atom_index] <= 0.0) {
            throw std::runtime_error(
                "Mass of atom " + std::to_string(atom_index + 1) + " must be positive");
        }
    }

    // Populate per-atom charges from the force field type definitions.
    // Charges default to 0 when no FF is provided or a type has no charge defined.
    auto charges = system.mutable_charges();
    if (ff != nullptr) {
        const auto types = system.atom_types();
        for (std::size_t i = 0; i < atom_count; ++i) {
            const int tidx = types[i];
            if (tidx >= 0 && static_cast<std::size_t>(tidx) < ff->elements.size()) {
                charges[i] = ff->elements[static_cast<std::size_t>(tidx)].charge;
            }
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
            config.time_step = parse_double(tokens[1], "time_step") / kTimeUnitConversion;
            if (config.time_step < 0.0) {
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
        // ---- Inline force field directives ----
        } else if (tokens[0] == "force_field") {
            if (tokens[1] != "lj") {
                throw std::runtime_error("Unsupported force_field type: " + tokens[1]
                                         + " (only 'lj' is supported)");
            }
            has_inline_ff = true;
        } else if (tokens[0] == "cutoff") {
            inline_ff.cutoff = parse_double(tokens[1], "cutoff");
            if (inline_ff.cutoff <= 0.0) throw std::runtime_error("cutoff must be positive");
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
            if (tokens[1] != "berendsen") {
                throw std::runtime_error("barostat must be 'berendsen'");
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
    const auto& ei = elements[static_cast<std::size_t>(type_i)];
    const auto& ej = elements[static_cast<std::size_t>(type_j)];
    eps_out = std::sqrt(ei.epsilon * ej.epsilon);   // Lorentz-Berthelot (geometric)
    sig_out = 0.5 * (ei.sigma + ej.sigma);           // Lorentz-Berthelot (arithmetic)
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

}  // namespace gmd