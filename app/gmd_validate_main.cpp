#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/bonded_force_provider.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/io/config_loader.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/system/special_pair_map.hpp"
#include "gmd/system/system.hpp"

namespace {

enum class ForceFieldFileKind {
    LJ,
    Molecular,
};

struct CommandLine {
    std::string xyz_path;
    std::string run_path;
    std::optional<std::string> ff_path;
    std::optional<std::string> top_path;
    std::string json_path;
};

CommandLine parse_command_line(int argc, char** argv) {
    std::vector<std::string> positionals;
    CommandLine cli;

    for (int arg_index = 1; arg_index < argc; ++arg_index) {
        const std::string arg = argv[arg_index];
        if (arg == "--json") {
            if (arg_index + 1 >= argc) {
                throw std::runtime_error("--json requires an output path");
            }
            cli.json_path = argv[++arg_index];
            continue;
        }
        if (arg.starts_with("--")) {
            throw std::runtime_error("Unknown command-line option: " + arg);
        }
        positionals.push_back(arg);
    }

    if (cli.json_path.empty()) {
        throw std::runtime_error("Usage: gmd_validate input.xyz run.in [ff.ff] [top.top] --json result.json");
    }
    if (positionals.size() < 2 || positionals.size() > 4) {
        throw std::runtime_error("Usage: gmd_validate input.xyz run.in [ff.ff] [top.top] --json result.json");
    }

    cli.xyz_path = positionals[0];
    cli.run_path = positionals[1];
    if (positionals.size() > 2) cli.ff_path = positionals[2];
    if (positionals.size() > 3) cli.top_path = positionals[3];
    return cli;
}

ForceFieldFileKind detect_force_field_file_kind(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open force field file: " + path.string());
    }

    std::string line;
    while (std::getline(input, line)) {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;
        if (!key.empty() && key[0] == '#') continue;
        if (key != "force_field") continue;

        std::string value;
        if (!(iss >> value)) {
            throw std::runtime_error(
                "force_field directive in " + path.string() + " is missing its value");
        }
        if (value == "lj") return ForceFieldFileKind::LJ;
        if (value == "molecular") return ForceFieldFileKind::Molecular;
        throw std::runtime_error("Unsupported force_field type in " + path.string() + ": " + value);
    }

    throw std::runtime_error("Could not determine force_field type from file: " + path.string());
}

gmd::ForceRequest make_request(const gmd::System& system) {
    const auto coords = system.coordinates();
    return {
        .system = &system,
        .box = &system.box(),
        .coordinates = std::span<const gmd::Coordinate3D>(coords.data(), coords.size()),
        .neighbor_list = nullptr,
    };
}

gmd::ForceResult compute_component(gmd::ForceProvider& provider,
                                   const gmd::System& system,
                                   gmd::RuntimeContext& runtime) {
    gmd::ForceResult result;
    provider.compute(make_request(system), result, runtime);
    if (!result.success) {
        throw std::runtime_error("Validation force evaluation failed");
    }
    return result;
}

double force_rms(const std::vector<gmd::Force3D>& forces) {
    if (forces.empty()) return 0.0;
    double sum_sq = 0.0;
    for (const auto& force : forces) {
        sum_sq += force[0] * force[0] + force[1] * force[1] + force[2] * force[2];
    }
    return std::sqrt(sum_sq / static_cast<double>(forces.size()));
}

double max_force_norm(const std::vector<gmd::Force3D>& forces) {
    double max_norm = 0.0;
    for (const auto& force : forces) {
        const double norm = std::sqrt(force[0] * force[0] +
                                      force[1] * force[1] +
                                      force[2] * force[2]);
        max_norm = std::max(max_norm, norm);
    }
    return max_norm;
}

void write_json_array(std::ostream& output, const std::vector<gmd::Force3D>& forces) {
    output << "[\n";
    for (std::size_t index = 0; index < forces.size(); ++index) {
        const auto& force = forces[index];
        output << "    [" << force[0] << ", " << force[1] << ", " << force[2] << "]";
        if (index + 1 != forces.size()) {
            output << ",";
        }
        output << "\n";
    }
    output << "  ]";
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment mpi_env(argc, argv);
#endif

    try {
        const CommandLine cli = parse_command_line(argc, argv);

        gmd::RuntimeContext runtime;
        gmd::ConfigLoader loader;
        gmd::System system;
        const gmd::RunConfig run_config = loader.load_run(cli.run_path);

        std::optional<gmd::LJForceFieldConfig> external_lj_ff;
        std::optional<gmd::MolecularForceFieldConfig> molecular_ff;
        std::shared_ptr<gmd::Topology> topology;
        const gmd::LJForceFieldConfig* xyz_ff = nullptr;

        if (run_config.force_field.has_value()) {
            xyz_ff = &run_config.force_field.value();
        } else if (cli.ff_path.has_value()) {
            const auto ff_kind = detect_force_field_file_kind(*cli.ff_path);
            if (ff_kind == ForceFieldFileKind::Molecular) {
                if (!cli.top_path.has_value()) {
                    throw std::runtime_error(
                        "Molecular force fields require a topology file as the fourth CLI argument");
                }
                molecular_ff = loader.load_molecular_ff(*cli.ff_path);
                topology = loader.load_topology(*cli.top_path);
                xyz_ff = &molecular_ff->lj;
            } else {
                external_lj_ff = loader.load_force_field(*cli.ff_path);
                xyz_ff = &external_lj_ff.value();
            }
        }

        loader.load_xyz(cli.xyz_path, system, xyz_ff);
        if (topology != nullptr) {
            system.set_special_pair_map(std::make_shared<gmd::SpecialPairMap>(
                *topology, run_config.special_pair_scales));
        }

        double lj_energy = 0.0;
        double bonded_energy = 0.0;
        double coulomb_energy = 0.0;
        std::string coulomb_method = "none";
        std::vector<gmd::Force3D> total_forces(system.atom_count(), gmd::Force3D{0.0, 0.0, 0.0});

        auto accumulate_component = [&](const gmd::ForceResult& result) {
            if (result.forces.size() != total_forces.size()) {
                throw std::runtime_error("Validation force buffer size mismatch");
            }
            for (std::size_t atom_index = 0; atom_index < total_forces.size(); ++atom_index) {
                total_forces[atom_index][0] += result.forces[atom_index][0];
                total_forces[atom_index][1] += result.forces[atom_index][1];
                total_forces[atom_index][2] += result.forces[atom_index][2];
            }
        };

        if (run_config.force_field.has_value()) {
            gmd::ClassicalForceProvider provider(run_config.force_field.value());
            const auto result = compute_component(provider, system, runtime);
            lj_energy = result.potential_energy;
            accumulate_component(result);
        } else if (external_lj_ff.has_value()) {
            gmd::ClassicalForceProvider provider(external_lj_ff.value());
            const auto result = compute_component(provider, system, runtime);
            lj_energy = result.potential_energy;
            accumulate_component(result);
        } else if (molecular_ff.has_value()) {
            gmd::BondedForceProvider bonded_provider(topology);
            for (const auto& params : molecular_ff->bond_types) bonded_provider.add_bond_type(params);
            for (const auto& params : molecular_ff->angle_types) bonded_provider.add_angle_type(params);
            for (const auto& params : molecular_ff->dihedral_types) bonded_provider.add_dihedral_type(params);
            for (const auto& params : molecular_ff->improper_types) bonded_provider.add_improper_type(params);
            const auto bonded_result = compute_component(bonded_provider, system, runtime);
            bonded_energy = bonded_result.potential_energy;
            accumulate_component(bonded_result);

            if (run_config.molecular_nonbonded_mode != "none") {
                gmd::ClassicalForceProvider lj_provider(molecular_ff->lj);
                const auto lj_result = compute_component(lj_provider, system, runtime);
                lj_energy = lj_result.potential_energy;
                accumulate_component(lj_result);
            }
        }

        if (run_config.coulomb.has_value()) {
            coulomb_method = run_config.coulomb->method;
            if (coulomb_method == "pme") {
                gmd::PMEForceProvider provider(run_config.coulomb->alpha,
                                               run_config.coulomb->real_cutoff,
                                               run_config.coulomb->pme_order,
                                               run_config.coulomb->pme_grid);
                provider.initialize(runtime);
                const auto result = compute_component(provider, system, runtime);
                provider.finalize(runtime);
                coulomb_energy = result.potential_energy;
                accumulate_component(result);
            } else {
                gmd::EwaldForceProvider provider(run_config.coulomb->alpha,
                                                 run_config.coulomb->kmax,
                                                 run_config.coulomb->real_cutoff);
                const auto result = compute_component(provider, system, runtime);
                coulomb_energy = result.potential_energy;
                accumulate_component(result);
            }
        }

        const auto json_parent = std::filesystem::path(cli.json_path).parent_path();
        if (!json_parent.empty()) {
            std::filesystem::create_directories(json_parent);
        }
        std::ofstream output(cli.json_path);
        if (!output.is_open()) {
            throw std::runtime_error("Failed to open validation JSON output: " + cli.json_path);
        }
        output << std::setprecision(15);

        output << "{\n"
               << "  \"metadata\": {\n"
               << "    \"xyz\": \"" << cli.xyz_path << "\",\n"
               << "    \"run\": \"" << cli.run_path << "\",\n"
               << "    \"atom_count\": " << system.atom_count() << ",\n"
               << "    \"coulomb_method\": \"" << coulomb_method << "\"\n"
               << "  },\n"
               << "  \"energy\": {\n"
               << "    \"total\": " << (lj_energy + bonded_energy + coulomb_energy) << ",\n"
               << "    \"components\": {\n"
               << "      \"lj\": " << lj_energy << ",\n"
               << "      \"bonded\": " << bonded_energy << ",\n"
               << "      \"coulomb\": " << coulomb_energy << "\n"
               << "    }\n"
               << "  },\n"
               << "  \"force\": {\n"
               << "    \"rms_norm\": " << force_rms(total_forces) << ",\n"
               << "    \"max_norm\": " << max_force_norm(total_forces) << ",\n"
               << "    \"atoms\": ";
        write_json_array(output, total_forces);
        output << "\n"
               << "  }\n"
               << "}\n";
    } catch (const std::exception& error) {
        std::cerr << "gmd_validate error: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
