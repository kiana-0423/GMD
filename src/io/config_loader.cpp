#include "gmd/io/config_loader.hpp"

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
    auto masses = system.mutable_masses();
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        tokens = tokenize_line(input);
        if (tokens.size() != 5 && tokens.size() != 8) {
            throw std::runtime_error("Each atom line in xyz input must contain either 5 items (type x y z mass) or 8 items (type x y z vx vy vz mass)");
        }

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
            throw std::runtime_error("Atom masses must be positive");
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

}  // namespace gmd