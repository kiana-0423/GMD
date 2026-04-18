#include "gmd/io/config_loader.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gmd/system/system.hpp"

namespace gmd {

namespace {

constexpr double kTimeUnitConversion = 1.018051e+1;

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

}  // namespace

void ConfigLoader::load_xyz(const std::filesystem::path& xyz_path, System& system) const {
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
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        tokens = tokenize_line(input);
        if (tokens.size() != 5 && tokens.size() != 8) {
            throw std::runtime_error(
                "Atom line " + std::to_string(atom_index + 1)
                + " must have 5 columns (type x y z mass) or 8 columns (type x y z vx vy vz mass),"
                  " got " + std::to_string(tokens.size()));
        }

        // Column 0: integer atom type (stored in system so force providers can use it).
        atom_types[atom_index] = parse_int(tokens[0], "atom_type");

        coordinates[atom_index] = {parse_double(tokens[1], "x"),
                                   parse_double(tokens[2], "y"),
                                   parse_double(tokens[3], "z")};
        auto velocities = system.mutable_velocities();
        if (tokens.size() == 8) {
            velocities[atom_index] = {parse_double(tokens[4], "vx"),
                                      parse_double(tokens[5], "vy"),
                                      parse_double(tokens[6], "vz")};
            masses[atom_index] = parse_double(tokens[7], "mass");
        } else {
            masses[atom_index] = parse_double(tokens[4], "mass");
        }
        if (masses[atom_index] <= 0.0) {
            throw std::runtime_error(
                "Mass of atom " + std::to_string(atom_index + 1) + " must be positive");
        }
    }
}

RunConfig ConfigLoader::load_run(const std::filesystem::path& run_path) const {
    std::ifstream input(run_path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open run input file: " + run_path.string());
    }

    RunConfig config;
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
        } else {
            throw std::runtime_error("Unsupported run input directive: " + tokens[0]);
        }
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
            // Multi-element format:
            //   element  Ar  epsilon  0.01032  sigma  3.405
            if (tokens.size() < 6) {
                throw std::runtime_error(
                    "element directive expects: element <name> epsilon <val> sigma <val>");
            }
            const std::string& elem_name = tokens[1];
            // Parse key-value pairs after the name (order-independent).
            LJElementParams ep;
            bool got_eps = false, got_sig = false;
            for (std::size_t k = 2; k + 1 < tokens.size(); k += 2) {
                if (tokens[k] == "epsilon") {
                    ep.epsilon = parse_double(tokens[k + 1], "epsilon");
                    if (ep.epsilon <= 0.0) throw std::runtime_error("element epsilon must be positive");
                    got_eps = true;
                } else if (tokens[k] == "sigma") {
                    ep.sigma = parse_double(tokens[k + 1], "sigma");
                    if (ep.sigma <= 0.0) throw std::runtime_error("element sigma must be positive");
                    got_sig = true;
                }
            }
            if (!got_eps || !got_sig) {
                throw std::runtime_error("element '" + elem_name + "' must specify both epsilon and sigma");
            }
            const int type_idx = static_cast<int>(config.elements.size());
            config.elements.push_back(ep);
            config.element_name_to_type[elem_name] = type_idx;
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
        config.elements.push_back({single_epsilon, single_sigma});
        config.element_name_to_type["default"] = 0;
    }

    return config;
}

}  // namespace gmd