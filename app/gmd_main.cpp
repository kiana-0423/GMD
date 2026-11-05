#include <algorithm>
#include <array>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#include "gmd/force/bonded_force_provider.hpp"
#include "gmd/core/simulation.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/composite_force_provider.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/integrator/berendsen_barostat.hpp"
#include "gmd/integrator/constraint_solver.hpp"
#include "gmd/integrator/mc_barostat.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/velocity_rescaling_thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/io/checkpoint.hpp"
#include "gmd/io/config_loader.hpp"
#include "gmd/io/file_utils.hpp"
#include "gmd/io/frame_writer.hpp"
#include "gmd/io/trajectory_writer.hpp"
#include "gmd/force/ml_force_provider.hpp"
#include "gmd/system/verlet_neighbor_builder.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/core/runtime_context.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/special_pair_map.hpp"
#include "gmd/system/system.hpp"
#ifdef GMD_ENABLE_TORCH
#include "gmd/force/torchscript_adapter.hpp"
#endif

namespace {

// --- Cooperative signal handling -------------------------------------------
//
// SIGTERM / SIGINT only record the signal; the MD loop acts on it at its next
// control poll (every control_poll_interval steps) exactly as on a stop file:
// final frame, final checkpoint, status file, exit 0. The handler then
// restores the default action, so a SECOND signal terminates immediately.
volatile std::sig_atomic_t g_stop_signal = 0;

extern "C" void gmd_record_stop_signal(int signal_number) {
    g_stop_signal = signal_number;
    std::signal(signal_number, SIG_DFL);
}

std::string signal_name(int signal_number) {
    if (signal_number == SIGTERM) return "SIGTERM";
    if (signal_number == SIGINT) return "SIGINT";
    return "signal " + std::to_string(signal_number);
}

// --- JSON status file --------------------------------------------------------
std::string json_string(const std::string& value) {
    std::string out = "\"";
    for (const char c : value) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buffer[8];
                    std::snprintf(buffer, sizeof(buffer), "\\u%04x", c);
                    out += buffer;
                } else {
                    out.push_back(c);
                }
        }
    }
    out.push_back('"');
    return out;
}

std::string json_string_or_null(const std::string& value) {
    return value.empty() ? std::string("null") : json_string(value);
}

std::string json_number(double value) {
    if (!std::isfinite(value)) return "null";
    char buffer[40];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    return buffer;
}

// Everything the status file reports. Filled in as the run progresses so the
// failure path can still write what is known.
struct RunStatus {
    std::filesystem::path path;      // empty -> no status file
    std::uint64_t start_step = 0;
    std::uint64_t step = 0;
    double time_step_fs = 0.0;
    std::string checkpoint;          // last checkpoint written (absolute)
    std::string frame_file;          // absolute, empty when off
    std::uint64_t frames_written = 0;
    std::string force_field_kind;
    std::string force_field_sha256;
    std::string model_sha256;
    std::string model_changed_from;
    std::string run_id;
    std::string stop_reason;
};

inline constexpr const char* kStatusFormat = "gmd-status-v1";

void write_status_file(const RunStatus& status, const std::string& state,
                       const std::string& message = {}) {
    if (status.path.empty()) return;
    std::ostringstream out;
    out << "{\n"
        << "  \"format_version\": " << json_string(kStatusFormat) << ",\n"
        << "  \"state\": " << json_string(state) << ",\n"
        << "  \"step\": " << status.step << ",\n"
        << "  \"start_step\": " << status.start_step << ",\n"
        << "  \"time_fs\": "
        << json_number(static_cast<double>(status.step) * status.time_step_fs) << ",\n"
        << "  \"checkpoint\": " << json_string_or_null(status.checkpoint) << ",\n"
        << "  \"checkpoint_version\": " << gmd::kCheckpointVersion << ",\n"
        << "  \"frame_file\": " << json_string_or_null(status.frame_file) << ",\n"
        << "  \"frame_format\": "
        << (status.frame_file.empty() ? std::string("null") : json_string(gmd::kFrameFormat))
        << ",\n"
        << "  \"frames_written\": " << status.frames_written << ",\n"
        << "  \"force_field_kind\": " << json_string_or_null(status.force_field_kind) << ",\n"
        << "  \"force_field_sha256\": " << json_string_or_null(status.force_field_sha256)
        << ",\n"
        << "  \"model_sha256\": " << json_string_or_null(status.model_sha256) << ",\n"
        << "  \"model_changed_from\": " << json_string_or_null(status.model_changed_from)
        << ",\n"
        << "  \"run_id\": " << json_string_or_null(status.run_id) << ",\n"
        << "  \"stop_reason\": " << json_string_or_null(status.stop_reason) << ",\n"
        << "  \"message\": " << json_string_or_null(message) << "\n"
        << "}\n";
    gmd::atomic_write_file(status.path, out.str());
}

std::string absolute_string(const std::filesystem::path& path) {
    if (path.empty()) return {};
    std::error_code ec;
    const auto absolute = std::filesystem::absolute(path, ec);
    return (ec ? path : absolute).lexically_normal().string();
}

enum class ForceFieldFileKind {
    LJ,
    Molecular,
};

struct CommandLine {
    std::string xyz_path = "xyz.in";
    std::string run_path = "run.in";
    std::optional<std::string> ff_path;
    std::optional<std::string> top_path;
    std::optional<int> expected_process_count;
    std::optional<std::array<int, 3>> process_grid;
};

int parse_positive_int_arg(const std::string& value, const char* name) {
    std::size_t parsed_chars = 0;
    const int parsed = std::stoi(value, &parsed_chars);
    if (parsed_chars != value.size() || parsed <= 0) {
        throw std::runtime_error(std::string(name) + " requires a positive integer");
    }
    return parsed;
}

CommandLine parse_command_line(int argc, char** argv) {
    std::vector<std::string> positionals;
    CommandLine cli;

    for (int arg_index = 1; arg_index < argc; ++arg_index) {
        const std::string arg = argv[arg_index];
        if (arg == "--np") {
            if (arg_index + 1 >= argc) {
                throw std::runtime_error("--np requires a positive process count");
            }

            cli.expected_process_count = parse_positive_int_arg(argv[++arg_index], "--np");
            continue;
        }
        if (arg == "--proc-grid") {
            if (arg_index + 3 >= argc) {
                throw std::runtime_error("--proc-grid requires Px Py Pz extents");
            }
            cli.process_grid = {
                parse_positive_int_arg(argv[++arg_index], "--proc-grid Px"),
                parse_positive_int_arg(argv[++arg_index], "--proc-grid Py"),
                parse_positive_int_arg(argv[++arg_index], "--proc-grid Pz"),
            };
            continue;
        }
        if (arg.starts_with("--")) {
            throw std::runtime_error("Unknown command-line option: " + arg);
        }

        positionals.push_back(arg);
    }

    if (positionals.size() > 4) {
        throw std::runtime_error(
            "Usage: gmd input.xyz run.in [ff.ff] [top.top] [--np N]"
            " [--proc-grid Px Py Pz]");
    }
    if (!positionals.empty()) cli.xyz_path = positionals[0];
    if (positionals.size() > 1) cli.run_path = positionals[1];
    if (positionals.size() > 2) cli.ff_path = positionals[2];
    if (positionals.size() > 3) cli.top_path = positionals[3];
    return cli;
}

std::string rank_prefix(int rank) {
    return "[gmd rank " + std::to_string(rank) + "] ";
}

gmd::PmeExecutionMode parse_pme_mode(const std::string& mode) {
    if (mode == "replicated") return gmd::PmeExecutionMode::Replicated;
    if (mode == "distributed") return gmd::PmeExecutionMode::Distributed;
    if (mode == "auto") return gmd::PmeExecutionMode::Auto;
    throw std::runtime_error("Unsupported pme_mode: " + mode);
}

int process_grid_size(const std::array<int, 3>& grid) {
    return grid[0] * grid[1] * grid[2];
}

void keep_rank_local_atoms(gmd::System& system,
                           const gmd::DomainDecomposition& decomposition,
                           int my_rank) {
    const auto box = system.box();
    const auto masses = system.masses();
    const auto charges = system.charges();
    const auto atom_types = system.atom_types();
    const auto molecule_ids = system.molecule_ids();
    const auto atomic_numbers = system.atomic_numbers();
    const auto coordinates = system.coordinates();
    const auto velocities = system.velocities();

    struct LocalAtom {
        double mass;
        double charge;
        int atom_type;
        int molecule_id;
        int atomic_number;
        int tag;
        gmd::System::Vec3 position;
        gmd::System::Vec3 velocity;
    };

    std::vector<LocalAtom> local_atoms;
    local_atoms.reserve(system.atom_count());
    for (std::size_t atom_index = 0; atom_index < system.atom_count(); ++atom_index) {
        if (decomposition.owner_rank(box, coordinates[atom_index]) != my_rank) {
            continue;
        }

        local_atoms.push_back(LocalAtom{
            .mass = masses[atom_index],
            .charge = charges[atom_index],
            .atom_type = atom_types[atom_index],
            .molecule_id = molecule_ids[atom_index],
            .atomic_number = atomic_numbers[atom_index],
            .tag = system.atom_tag(atom_index),
            .position = coordinates[atom_index],
            .velocity = velocities[atom_index],
        });
    }

    system.resize(local_atoms.size(), local_atoms.size());
    system.set_box(box);
    auto local_masses = system.mutable_masses();
    auto local_charges = system.mutable_charges();
    auto local_atom_types = system.mutable_atom_types();
    auto local_molecule_ids = system.mutable_molecule_ids();
    auto local_atomic_numbers = system.mutable_atomic_numbers();
    auto local_coordinates = system.mutable_coordinates();
    auto local_velocities = system.mutable_velocities();
    auto local_tags = system.mutable_atom_tags();
    auto local_owners = system.mutable_atom_owners();

    for (std::size_t atom_index = 0; atom_index < local_atoms.size(); ++atom_index) {
        const LocalAtom& atom = local_atoms[atom_index];
        local_masses[atom_index] = atom.mass;
        local_charges[atom_index] = atom.charge;
        local_atom_types[atom_index] = atom.atom_type;
        local_molecule_ids[atom_index] = atom.molecule_id;
        local_atomic_numbers[atom_index] = atom.atomic_number;
        local_coordinates[atom_index] = atom.position;
        local_velocities[atom_index] = atom.velocity;
        local_tags[atom_index] = atom.tag;
        local_owners[atom_index] = my_rank;
    }
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
        throw std::runtime_error(
            "Unsupported force_field type in " + path.string() + ": " + value);
    }

    throw std::runtime_error(
        "Could not determine force_field type from file: " + path.string());
}

}  // namespace

int main(int argc, char** argv)
{
#ifdef GMD_ENABLE_MPI
    gmd::MpiEnvironment mpi_env(argc, argv);
#endif

    gmd::RuntimeContext runtime;
    const int my_rank = runtime.rank();
    const int nprocs = runtime.size();
    auto mpi_comm = std::make_shared<gmd::MpiCommunicator>();
    const bool is_root_rank = my_rank == 0;
    const std::string log_prefix = rank_prefix(my_rank);

    CommandLine cli;
    try {
        cli = parse_command_line(argc, argv);
        if (cli.expected_process_count.has_value() &&
            *cli.expected_process_count != nprocs) {
            throw std::runtime_error(
                "--np " + std::to_string(*cli.expected_process_count) +
                " does not match the MPI world size " + std::to_string(nprocs));
        }
    } catch (const std::exception& error) {
        std::cerr << log_prefix << "Error: " << error.what() << "\n";
        return 1;
    }

    std::filesystem::path xyz_path = cli.xyz_path;
    std::filesystem::path run_path = cli.run_path;
    std::optional<std::filesystem::path> ff_path =
        cli.ff_path.has_value() ? std::optional<std::filesystem::path>(*cli.ff_path)
                                : std::nullopt;
    std::optional<std::filesystem::path> top_path =
        cli.top_path.has_value() ? std::optional<std::filesystem::path>(*cli.top_path)
                                 : std::nullopt;

    // Output stem: same directory as xyz input, base name "output".
    const std::filesystem::path output_stem =
        xyz_path.parent_path() / "output";

    // Cooperative SIGTERM / SIGINT: acted on at the next control poll.
    std::signal(SIGTERM, gmd_record_stop_signal);
    std::signal(SIGINT, gmd_record_stop_signal);

    RunStatus status;

    try {
        gmd::ConfigLoader loader;
        gmd::System system;

        // Load run config first so inline force field is available before the xyz file is parsed.
        const gmd::RunConfig run_config = loader.load_run(run_path);
        if (is_root_rank) {
            status.path = run_config.status_file;
        }
        status.run_id = run_config.run_id;
        status.time_step_fs = run_config.time_step_fs;
        const std::uint64_t output_interval = run_config.output_interval;
        const bool is_restart = !run_config.restart_from.empty();
        std::optional<gmd::CheckpointMetadata> restart_metadata;

        std::optional<gmd::LJForceFieldConfig> external_lj_ff;
        std::optional<gmd::MolecularForceFieldConfig> molecular_ff;
        std::shared_ptr<gmd::Topology> topology;

        // Resolve external FF before xyz loading so type→mass mapping works.
        const gmd::LJForceFieldConfig* xyz_ff = nullptr;
        if (is_restart) {
            auto checkpoint_topology = std::make_shared<gmd::Topology>();
            restart_metadata = gmd::read_checkpoint(run_config.restart_from,
                                                    system,
                                                    checkpoint_topology.get());
            if (run_config.force_field_type == "ml" && !restart_metadata->has_atomic_numbers) {
                throw std::runtime_error(
                    "Checkpoint " + run_config.restart_from.string() + " is version " +
                    std::to_string(restart_metadata->version) +
                    ", which does not store atomic numbers; an ML force field needs the "
                    "species of every atom. Restart ML runs only from version " +
                    std::to_string(gmd::kFirstCheckpointVersionWithAtomicNumbers) +
                    "+ checkpoints (or start again from the xyz file).");
            }
            if (!restart_metadata->force_field_file.empty() && !ff_path.has_value()) {
                ff_path = restart_metadata->force_field_file;
            }
            if (!restart_metadata->topology_file.empty() && !top_path.has_value()) {
                top_path = restart_metadata->topology_file;
            }
            if (!checkpoint_topology->bonds.empty() ||
                !checkpoint_topology->angles.empty() ||
                !checkpoint_topology->dihedrals.empty() ||
                !checkpoint_topology->impropers.empty() ||
                !checkpoint_topology->constraints.empty()) {
                topology = checkpoint_topology;
            }
            if (is_root_rank) {
                std::cout << log_prefix << "Restarting from checkpoint "
                          << run_config.restart_from.string()
                          << " at step " << restart_metadata->step << "\n";
            }
        }

        if (run_config.force_field.has_value()) {
            xyz_ff = &run_config.force_field.value();
        } else if (ff_path.has_value() && std::filesystem::exists(*ff_path)) {
            const auto ff_kind = detect_force_field_file_kind(*ff_path);
            if (ff_kind == ForceFieldFileKind::Molecular) {
                if (!top_path.has_value() && topology == nullptr) {
                    throw std::runtime_error(
                        "Molecular force fields require a topology file as the fourth CLI argument");
                }
                molecular_ff = loader.load_molecular_ff(*ff_path);
                if (topology == nullptr) {
                    topology = loader.load_topology(*top_path);
                }
                xyz_ff = &molecular_ff->lj;
            } else {
                external_lj_ff = loader.load_force_field(*ff_path);
                xyz_ff = &external_lj_ff.value();
            }
        }

        if (!is_restart) {
            loader.load_xyz(xyz_path, system, xyz_ff);
        }
        if (topology != nullptr) {
            system.set_special_pair_map(std::make_shared<gmd::SpecialPairMap>(
                *topology, run_config.special_pair_scales));
        }

        // --- Force provider ---
        std::shared_ptr<gmd::ClassicalForceProvider> lj_provider;
        std::shared_ptr<gmd::ForceProvider> active_provider;
        double short_range_cutoff = 0.0;
        bool need_neighbor_builder = false;
        if (run_config.force_field.has_value()) {
            const auto& ff_config = run_config.force_field.value();
            lj_provider = std::make_shared<gmd::ClassicalForceProvider>(ff_config);
            active_provider = lj_provider;
            short_range_cutoff = lj_provider->cutoff();
            need_neighbor_builder = true;
            if (is_root_rank) {
                std::cout << log_prefix << "Loaded inline force field from " << run_path
                          << " (" << ff_config.elements.size() << " element type(s))\n";
            }
        } else if (molecular_ff.has_value()) {
            auto bonded_provider = std::make_shared<gmd::BondedForceProvider>(topology);
            for (const auto& bp : molecular_ff->bond_types) {
                bonded_provider->add_bond_type(bp);
            }
            for (const auto& ap : molecular_ff->angle_types) {
                bonded_provider->add_angle_type(ap);
            }
            for (const auto& dp : molecular_ff->dihedral_types) {
                bonded_provider->add_dihedral_type(dp);
            }
            for (const auto& ip : molecular_ff->improper_types) {
                bonded_provider->add_improper_type(ip);
            }

            active_provider = bonded_provider;
            short_range_cutoff = molecular_ff->lj.cutoff;

            if (is_root_rank) {
                std::cout << log_prefix << "Loaded force field from " << ff_path->string()
                          << " (" << molecular_ff->lj.elements.size() << " atom type(s), "
                          << topology->bonds.size() << " bond(s), "
                          << topology->angles.size() << " angle(s), "
                          << topology->dihedrals.size() << " dihedral(s), "
                          << topology->impropers.size() << " improper(s))\n";
                if (top_path.has_value()) {
                    std::cout << log_prefix << "Loaded topology from " << top_path->string() << "\n";
                } else {
                    std::cout << log_prefix << "Loaded topology from checkpoint\n";
                }
            }
            if (run_config.molecular_nonbonded_mode != "none") {
                lj_provider = std::make_shared<gmd::ClassicalForceProvider>(molecular_ff->lj);
                auto composite = std::make_shared<gmd::CompositeForceProvider>();
                composite->add(lj_provider);
                composite->add(bonded_provider);
                active_provider = composite;
                short_range_cutoff = lj_provider->cutoff();
                need_neighbor_builder = true;
                if (is_root_rank) {
                    std::cout << log_prefix
                              << "Molecular non-bonded mode: topology special pairs "
                              << "(1-2/1-3 exclusions, 1-4 scaling).\n";
                }
            } else {
                if (is_root_rank) {
                    std::cout << log_prefix << "Molecular non-bonded mode: bonded-only.\n";
                }
            }
        } else if (external_lj_ff.has_value()) {
            const auto& ff_config = external_lj_ff.value();
            lj_provider = std::make_shared<gmd::ClassicalForceProvider>(ff_config);
            active_provider = lj_provider;
            short_range_cutoff = lj_provider->cutoff();
            need_neighbor_builder = true;
            if (is_root_rank) {
                std::cout << log_prefix << "Loaded force field from " << ff_path->string()
                          << " (" << ff_config.elements.size() << " element type(s))\n";
            }
        } else if (run_config.force_field_type == "ml") {
            if (nprocs > 1) {
                throw std::runtime_error(
                    "ML force provider cannot run with MPI domain decomposition: "
                    "local-plus-ghost model energy ownership and message-passing halo depth "
                    "are not defined");
            }
#ifdef GMD_ENABLE_TORCH
            if (run_config.ml_model_path.empty()) {
                throw std::runtime_error(
                    "force_field ml requires a 'model_path' directive in the run file");
            }
            auto ts_adapter = std::make_shared<gmd::TorchScriptModelRuntimeAdapter>();
            auto ml_provider = std::make_shared<gmd::MLForceProvider>(
                run_config.ml_model_path, ts_adapter);
            // Load the model now so cutoff() is available before NeighborBuilder is set up.
            gmd::RuntimeContext tmp_runtime;
            ml_provider->initialize(tmp_runtime);
            short_range_cutoff = static_cast<double>(ml_provider->cutoff());
            if (short_range_cutoff <= 0.0) {
                throw std::runtime_error(
                    "ML model 'local_cutoff' attribute is missing or <= 0");
            }
            active_provider = ml_provider;
            need_neighbor_builder = true;
            if (is_root_rank) {
                std::cout << log_prefix << "Loaded ML model from " << run_config.ml_model_path.string()
                          << "  cutoff=" << short_range_cutoff << " \u00c5\n";
            }
#else
            throw std::runtime_error(
                "force_field ml requires GMD to be built with -DGMD_ENABLE_TORCH=ON");
#endif
        } else {
            lj_provider = std::make_shared<gmd::ClassicalForceProvider>();
            active_provider = lj_provider;
            short_range_cutoff = lj_provider->cutoff();
            need_neighbor_builder = true;
            if (is_root_rank) {
                std::cout << log_prefix << "No force field supplied; using default Ar LJ parameters.\n";
            }
        }

        // --- Long-range Coulomb (Ewald or PME) ---
        // If a Coulomb section is present in run.in, add it to the active provider.
        if (run_config.coulomb.has_value()) {
            const auto& cc = *run_config.coulomb;
            const double coulomb_cutoff =
                cc.real_cutoff > 0.0 ? cc.real_cutoff
                                     : (short_range_cutoff > 0.0 ? short_range_cutoff : 8.5);
            short_range_cutoff = std::max(short_range_cutoff, coulomb_cutoff);
            need_neighbor_builder = true;
            auto composite = std::dynamic_pointer_cast<gmd::CompositeForceProvider>(active_provider);
            if (!composite) {
                composite = std::make_shared<gmd::CompositeForceProvider>();
                composite->add(active_provider);
                active_provider = composite;
            }

            if (cc.method == "pme") {
                auto pme = std::make_shared<gmd::PMEForceProvider>(
                    cc.alpha,
                    cc.real_cutoff,
                    cc.pme_order,
                    cc.pme_grid,
                    parse_pme_mode(cc.pme_mode),
                    cc.pme_benchmark);
                composite->add(pme);
                if (is_root_rank) {
                    std::cout << log_prefix << "Coulomb: PME  order=" << cc.pme_order
                              << "  grid=" << cc.pme_grid[0] << "x"
                              << cc.pme_grid[1] << "x" << cc.pme_grid[2]
                              << "  mode=" << cc.pme_mode
                              << (cc.pme_benchmark ? "  benchmark=on" : "")
                              << "\n";
                }
            } else {
                // Default to Ewald.
                auto ewald = std::make_shared<gmd::EwaldForceProvider>(
                    cc.alpha, cc.kmax, cc.real_cutoff);
                composite->add(ewald);
                if (is_root_rank) {
                    std::cout << log_prefix << "Coulomb: Ewald  alpha=" << cc.alpha
                              << "  kmax=" << cc.kmax
                              << "  r_cut=" << cc.real_cutoff << "\n";
                }
            }
        }

        // --- Force-field provenance -------------------------------------------
        //
        // What produced the forces, recorded in every checkpoint, frame and the
        // status file, and compared on restart. The fingerprint is the SHA-256
        // of the model file (ML), of the force-field file (file-based fields),
        // or of the canonical parameter summary (inline / default LJ), with
        // any Coulomb settings folded in.
        std::string force_field_summary;
        std::string force_field_kind;
        std::string force_field_fingerprint;
        std::string model_sha256;
        if (run_config.force_field.has_value()) {
            std::ostringstream ff_summary;
            ff_summary.precision(17);
            ff_summary << "inline_lj cutoff " << run_config.force_field->cutoff
                       << " mixing_rule " << run_config.force_field->mixing_rule
                       << " types " << run_config.force_field->elements.size();
            for (std::size_t type_index = 0;
                 type_index < run_config.force_field->elements.size();
                 ++type_index) {
                const auto& type = run_config.force_field->elements[type_index];
                ff_summary << " type " << type_index
                           << ' ' << type.element
                           << " mass " << type.mass
                           << " epsilon " << type.epsilon
                           << " sigma " << type.sigma
                           << " charge " << type.charge;
            }
            for (const auto& [pair, override_params] :
                 run_config.force_field->pair_overrides) {
                ff_summary << " pair " << pair.first << ' ' << pair.second
                           << " epsilon " << override_params.epsilon
                           << " sigma " << override_params.sigma;
            }
            force_field_summary = ff_summary.str();
            force_field_kind = "inline_lj";
            force_field_fingerprint = gmd::sha256_hex(force_field_summary);
        } else if (molecular_ff.has_value() || external_lj_ff.has_value()) {
            force_field_summary = "file " + ff_path->string();
            force_field_kind = molecular_ff.has_value() ? "molecular" : "lj_file";
            force_field_fingerprint = gmd::sha256_file_hex(*ff_path);
        } else if (run_config.force_field_type == "ml") {
            model_sha256 = gmd::sha256_file_hex(run_config.ml_model_path);
            force_field_summary = "ml model " + run_config.ml_model_path.string();
            force_field_kind = "ml";
            force_field_fingerprint = model_sha256;
        } else {
            force_field_summary = "default_lj";
            force_field_kind = "default_lj";
            force_field_fingerprint = gmd::sha256_hex(force_field_summary);
        }
        if (run_config.coulomb.has_value()) {
            const auto& cc = *run_config.coulomb;
            std::ostringstream coulomb;
            coulomb.precision(17);
            coulomb << cc.method << " alpha " << cc.alpha << " real_cutoff " << cc.real_cutoff
                    << " kmax " << cc.kmax << " order " << cc.pme_order << " grid "
                    << cc.pme_grid[0] << ' ' << cc.pme_grid[1] << ' ' << cc.pme_grid[2];
            force_field_kind += "+" + cc.method;
            force_field_fingerprint =
                gmd::sha256_hex(force_field_fingerprint + " coulomb " + coulomb.str());
        }
        status.force_field_kind = force_field_kind;
        status.force_field_sha256 = force_field_fingerprint;
        status.model_sha256 = model_sha256;
        if (is_root_rank) {
            std::cout << log_prefix << "Force field provenance: kind=" << force_field_kind
                      << " sha256=" << force_field_fingerprint << "\n";
        }

        // --- Restart provenance check ------------------------------------------
        std::string model_changed_from;
        if (is_restart) {
            const auto& md = *restart_metadata;
            if (!md.has_atomic_numbers) {
                if (is_root_rank) {
                    std::cout << log_prefix << "Checkpoint version " << md.version
                              << " records no force-field provenance; the force field "
                                 "is not validated against it.\n";
                }
            } else if (md.force_field_kind != force_field_kind ||
                       md.force_field_sha256 != force_field_fingerprint) {
                const std::string old_id =
                    md.force_field_kind + " sha256=" + md.force_field_sha256;
                const std::string new_id =
                    force_field_kind + " sha256=" + force_field_fingerprint;
                if (!run_config.allow_model_change) {
                    throw std::runtime_error(
                        "Force field differs from the checkpoint " +
                        run_config.restart_from.string() + ": checkpoint " + old_id +
                        ", this run " + new_id +
                        ". Add 'model_change allow' to the run file to restart with a "
                        "different force field or model.");
                }
                model_changed_from = md.model_sha256.empty() ? md.force_field_sha256
                                                             : md.model_sha256;
                status.model_changed_from = model_changed_from;
                if (is_root_rank) {
                    std::cout << log_prefix << "Model change on restart at step " << md.step
                              << " (allowed by 'model_change allow'): " << old_id << " -> "
                              << new_id << "\n";
                }
            }
        }

        // --- Neighbor builder (uses the largest active short-range cutoff, r_skin = 2.0 Å) ---
        constexpr double r_skin = 2.0;
        std::shared_ptr<gmd::VerletNeighborBuilder> neighbor_builder;
        if (need_neighbor_builder) {
            neighbor_builder = std::make_shared<gmd::VerletNeighborBuilder>(
                short_range_cutoff, r_skin);
        }

        const std::size_t global_atom_count = system.atom_count();
        gmd::System output_system;
        if (is_root_rank) {
            output_system = system;
        }
        std::shared_ptr<gmd::DomainDecomposition> domain_decomposition;
        if (nprocs > 1) {
            domain_decomposition = std::make_shared<gmd::DomainDecomposition>();
            const auto proc_grid = cli.process_grid.has_value()
                ? cli.process_grid
                : run_config.mpi_grid;
            if (proc_grid.has_value()) {
                if (process_grid_size(*proc_grid) != nprocs) {
                    throw std::runtime_error(
                        "MPI processor grid product does not match the MPI world size");
                }
                domain_decomposition->create_decomposition(system.box(),
                                                           *proc_grid,
                                                           my_rank,
                                                           short_range_cutoff,
                                                           r_skin,
                                                           {true, true, true});
            } else {
                domain_decomposition->create_decomposition(system.box(),
                                                           nprocs,
                                                           my_rank,
                                                           short_range_cutoff,
                                                           r_skin,
                                                           {true, true, true});
            }
            keep_rank_local_atoms(system, *domain_decomposition, my_rank);
            const auto grid = domain_decomposition->info().proc_grid;
            std::cout << log_prefix << "Domain decomposition owns "
                      << system.num_local_atoms() << " of " << global_atom_count
                      << " atoms before ghost exchange on grid "
                      << grid[0] << "x" << grid[1] << "x" << grid[2] << "\n";
        }

        // --- Integrator ---
        auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(run_config.time_step);
        integrator->set_target_temperature(run_config.temperature);
        if (run_config.constraints_enabled) {
            if (topology == nullptr || !molecular_ff.has_value()) {
                throw std::runtime_error(
                    "constraints require a molecular topology and force-field bond parameters");
            }
            std::vector<double> bond_type_distances;
            bond_type_distances.reserve(molecular_ff->bond_types.size());
            for (const auto& params : molecular_ff->bond_types) {
                bond_type_distances.push_back(params.r0);
            }
            auto constraints = gmd::constraints_from_bond_types(
                *topology, run_config.constrained_bond_types, bond_type_distances);
            if (constraints.empty()) {
                throw std::runtime_error(
                    "constraints are enabled, but no explicit constraints or constrained bond types were found");
            }
            auto constraint_solver = std::make_shared<gmd::ConstraintSolver>(
                std::move(constraints), run_config.constraint_settings);
            integrator->set_constraint_solver(constraint_solver);
            if (is_root_rank) {
                // ConstraintSolver does not log; it exposes what it normalised
                // and the application reports it here, once.
                const auto& diagnostics = constraint_solver->normalization_diagnostics();
                if (!diagnostics.empty()) {
                    std::cout << log_prefix << "Constraints: collapsed "
                              << diagnostics.exact_duplicates << " exact duplicate(s) and "
                              << diagnostics.tolerance_equivalent_duplicates
                              << " tolerance-equivalent duplicate(s), of which "
                              << diagnostics.reversed_duplicates
                              << " were the same pair listed in reverse order\n";
                    for (const auto& entry : diagnostics.discarded_targets) {
                        std::cout << log_prefix << "  discarded constraint target: "
                                  << entry << "\n";
                    }
                }
                std::cout << log_prefix << "Constraints: SHAKE/RATTLE enabled  tolerance="
                          << run_config.constraint_settings.tolerance
                          << "  max_iterations="
                          << run_config.constraint_settings.max_iterations
                          << "  rattle="
                          << (run_config.constraint_settings.enable_rattle ? "on" : "off")
                          << "\n";
            }
        }

        // --- Thermostat ---
        std::shared_ptr<gmd::Thermostat> thermostat;
        if (run_config.thermostat_type == "nose_hoover") {
            auto tstat = std::make_shared<gmd::NoseHooverThermostat>(run_config.thermostat_tau);
            thermostat = tstat;
            integrator->set_thermostat(tstat);
            if (is_root_rank) {
                std::cout << log_prefix << "Thermostat: Nose-Hoover  tau=" << run_config.thermostat_tau_fs << " fs\n";
            }
        } else if (run_config.thermostat_type == "velocity_rescaling") {
            auto tstat = std::make_shared<gmd::VelocityRescalingThermostat>();
            thermostat = tstat;
            integrator->set_thermostat(tstat);
            if (is_root_rank) {
                std::cout << log_prefix << "Thermostat: velocity rescaling\n";
            }
        }

        // --- Barostat ---
        std::shared_ptr<gmd::Barostat> barostat;
        if (run_config.barostat_type == "berendsen") {
            auto bstat = std::make_shared<gmd::BerendsenBarostat>(
                run_config.barostat_tau, run_config.compressibility);
            barostat = bstat;
            integrator->set_barostat(bstat);
            integrator->set_target_pressure(run_config.target_pressure);
            if (is_root_rank) {
                std::cout << log_prefix << "Barostat: Berendsen  P=" << run_config.target_pressure
                          << " bar  tau=" << run_config.barostat_tau_fs << " fs\n";
            }
        } else if (run_config.barostat_type == "monte_carlo") {
            if (nprocs > 1) {
                throw std::runtime_error(
                    "Monte Carlo barostat cannot run with MPI domain decomposition: "
                    "trial volume moves need coordinated ghost refresh, global trial energy, "
                    "and one accept/reject decision");
            }
            auto bstat = std::make_shared<gmd::MCBarostat>(
                run_config.mc_frequency,
                run_config.mc_volume_step);
            barostat = bstat;
            integrator->set_barostat(bstat);
            integrator->set_target_pressure(run_config.target_pressure);
            if (is_root_rank) {
                std::cout << log_prefix << "Barostat: Monte Carlo NPT  P=" << run_config.target_pressure
                          << " bar  freq=" << run_config.mc_frequency
                          << "  max_delta_ln_V=" << run_config.mc_volume_step << "\n";
            }
        }

        // --- Velocity initializer ---
        const auto velocity_mode =
            run_config.velocity_init_mode == "input" ? gmd::VelocityInitMode::FromInput
                                                     : gmd::VelocityInitMode::Random;
        auto velocity_initializer = std::make_shared<gmd::VelocityInitializer>(run_config.velocity_seed);

        // --- Assemble simulation ---
        gmd::Simulation simulation(&system);
        if (!is_restart) {
            simulation.set_velocity_initializer(velocity_initializer);
            simulation.set_velocity_init_mode(velocity_mode);
        }
        simulation.set_remove_center_of_mass_velocity(run_config.remove_center_of_mass_velocity);
        simulation.set_initial_temperature(run_config.temperature);
        simulation.set_force_provider(active_provider);
        simulation.set_mpi_communicator(mpi_comm);
        if (domain_decomposition) {
            simulation.set_domain_decomposition(domain_decomposition);
        }
        if (neighbor_builder) {
            simulation.set_neighbor_builder(neighbor_builder);
        }
        simulation.set_integrator(integrator);
        simulation.set_time_step(run_config.time_step);

        // --- Trajectory writer ---
        gmd::TrajectoryWriter writer;
        if (is_root_rank) {
            writer.open(output_stem);
            std::cout << log_prefix << "Writing trajectory to " << output_stem.string() << ".xyz"
                      << " and energy log to " << output_stem.string() << ".log\n";
            if (!model_changed_from.empty()) {
                writer.write_log_comment(
                    "model_change at step " + std::to_string(restart_metadata->step) + ": " +
                    restart_metadata->force_field_kind + " " + model_changed_from + " -> " +
                    force_field_kind + " " +
                    (model_sha256.empty() ? force_field_fingerprint : model_sha256));
            }
        }

        // --- ASE-readable frame file (optional) ---
        const bool frame_output_enabled = !run_config.frame_file.empty();
        const std::uint64_t frame_interval =
            run_config.frame_interval.value_or(output_interval);
        gmd::FrameWriter frame_writer;
        if (frame_output_enabled) {
            status.frame_file = absolute_string(run_config.frame_file);
            if (is_root_rank) {
                frame_writer.open(run_config.frame_file, run_config.frame_velocities,
                                  gmd::FrameProvenance{
                                      .run_id = run_config.run_id,
                                      .force_field_kind = force_field_kind,
                                      .model_sha256 = model_sha256,
                                      .model_changed_from = model_changed_from,
                                  });
                std::cout << log_prefix << "Writing " << gmd::kFrameFormat << " frames to "
                          << run_config.frame_file.string() << " every "
                          << frame_interval << " step(s)"
                          << (run_config.frame_velocities ? " with velocities" : "") << "\n";
            }
        }

        simulation.initialize(runtime);

        if (is_restart) {
            // Must FOLLOW initialize(): that call evaluates the forces at the
            // checkpointed coordinates and installs the provider virial, which
            // drops any constraint term attached to an earlier geometry and
            // clears the completed-step pressure.
            //
            // The constraint virial is the ENDPOINT RATTLE value of the step
            // that produced this state, so it belongs to exactly these
            // coordinates and is attached to that freshly recomputed provider
            // virial. The completed-step pressure is installed as recorded, so
            // the restarted run reports for this frame the identical number the
            // uninterrupted run reported, rather than reconstructing it.
            //
            // metadata.provider_virial is deliberately NOT installed: the run
            // has just recomputed it from these coordinates, and a value summed
            // under a different rank decomposition would only differ in
            // round-off. It is persisted for diagnosis, not for restoration.
            if (restart_metadata->constraint_virial_state == "valid") {
                if (restart_metadata->constraint_virial_time_level !=
                    "endpoint_rattle_t_plus_dt") {
                    throw std::runtime_error(
                        "Checkpoint carries a constraint virial at time level '" +
                        restart_metadata->constraint_virial_time_level +
                        "', which this build cannot pair with the endpoint provider "
                        "virial it evaluates at the checkpointed coordinates");
                }
                system.set_constraint_virial(restart_metadata->constraint_virial);
            }
            if (restart_metadata->step_pressure_valid) {
                gmd::System::StepThermodynamics record;
                record.valid = true;
                record.pressure = restart_metadata->step_pressure;
                record.virial = restart_metadata->step_pressure_virial;
                record.twice_kinetic_energy = restart_metadata->step_pressure_twice_ke;
                record.volume = restart_metadata->step_pressure_volume;
                record.potential_energy = restart_metadata->step_pressure_potential_energy;
                system.set_step_thermodynamics(record);
            }
        }

        // Degrees of freedom for every temperature the run reports. Read back
        // from the integrator rather than recomputed here, so the trajectory
        // log and the thermostat are guaranteed to use the same count: 3N,
        // less 3 when the COM velocity is removed, less one per constraint.
        // Independence warnings from the check initialize() performed. A set that
        // is dependent has already thrown; these are the ones that are
        // independent but close enough to degenerate to be worth saying so.
        if (is_root_rank) {
            for (const auto& warning : integrator->constraint_rank_report().warnings) {
                std::cout << log_prefix << "Constraints: WARNING: " << warning << "\n";
            }
        }

        const std::size_t dof = integrator->degrees_of_freedom(system);
        if (is_root_rank) {
            std::cout << log_prefix << "Degrees of freedom: " << dof
                      << "  (3N=" << 3 * global_atom_count
                      << (run_config.remove_center_of_mass_velocity ? ", -3 COM" : ", COM kept")
                      << ", -" << integrator->constraint_count() << " constraints)\n";
        }
        if (dof == 0) {
            throw std::runtime_error(
                "System has zero degrees of freedom after removing centre-of-mass "
                "motion and constraints; temperature is undefined");
        }
        if (is_restart) {
            if (restart_metadata->thermostat_type != run_config.thermostat_type) {
                throw std::runtime_error(
                    "Checkpoint thermostat type does not match run input");
            }
            if (restart_metadata->barostat_type != run_config.barostat_type) {
                throw std::runtime_error(
                    "Checkpoint barostat type does not match run input");
            }
            if (thermostat != nullptr) {
                thermostat->load_checkpoint_state(restart_metadata->thermostat_state);
            }
            if (barostat != nullptr) {
                barostat->load_checkpoint_state(restart_metadata->barostat_state);
            }
            simulation.set_current_step(restart_metadata->step);
        }

        // Writes the legacy trajectory frame (output.xyz/.log) and/or a frame-file
        // frame for `step`. Collective under MPI: every rank must call it with
        // the same flags, which are derived from the step and the run config.
        std::optional<std::uint64_t> last_legacy_frame_step;
        std::optional<std::uint64_t> last_frame_file_step;
        auto write_global_frame = [&](std::uint64_t step, double time, bool legacy,
                                      bool frame) {
            if (!legacy && !frame) {
                return;
            }
            // compute_twice_ke already performs MPI_Allreduce internally to
            // return the global kinetic energy; do NOT double-wrap here.
            const double twice_ke = gmd::compute_twice_ke(system);
            if (legacy) last_legacy_frame_step = step;
            if (frame) last_frame_file_step = step;
            if (nprocs == 1) {
                if (is_root_rank) {
                    if (legacy) writer.write_frame(system, step, time, twice_ke, dof);
                    if (frame) {
                        frame_writer.write_frame(system, step, time, twice_ke, dof);
                        status.frames_written = frame_writer.frame_count();
                    }
                }
                return;
            }

            // Coordinates always; forces and velocities only for a frame-file
            // frame, which reports them.
            const std::size_t fields_per_atom = frame ? 9 : 3;
            std::vector<double> local_fields(global_atom_count * fields_per_atom, 0.0);
            const auto coordinates = system.coordinates();
            const auto forces = system.forces();
            const auto velocities = system.velocities();
            for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
                const int tag = system.atom_tag(atom_index);
                if (tag < 0 || static_cast<std::size_t>(tag) >= global_atom_count) {
                    throw std::runtime_error("Local atom tag is outside the output coordinate map");
                }

                const std::size_t offset = static_cast<std::size_t>(tag) * fields_per_atom;
                for (int d = 0; d < 3; ++d) {
                    local_fields[offset + d] = coordinates[atom_index][d];
                    if (frame) {
                        local_fields[offset + 3 + d] = forces[atom_index][d];
                        local_fields[offset + 6 + d] = velocities[atom_index][d];
                    }
                }
            }

            std::vector<double> global_fields;
            mpi_comm->allreduce_vector(local_fields, global_fields);
            if (is_root_rank) {
                auto frame_coordinates = output_system.mutable_coordinates();
                auto frame_forces = output_system.mutable_forces();
                auto frame_velocities = output_system.mutable_velocities();
                for (std::size_t atom_index = 0; atom_index < global_atom_count; ++atom_index) {
                    const std::size_t offset = atom_index * fields_per_atom;
                    frame_coordinates[atom_index] = {
                        global_fields[offset],
                        global_fields[offset + 1],
                        global_fields[offset + 2]
                    };
                    if (frame) {
                        frame_forces[atom_index] = {
                            global_fields[offset + 3],
                            global_fields[offset + 4],
                            global_fields[offset + 5]
                        };
                        frame_velocities[atom_index] = {
                            global_fields[offset + 6],
                            global_fields[offset + 7],
                            global_fields[offset + 8]
                        };
                    }
                }
                // output_system holds only the gathered global coordinates; it
                // takes part in no dynamics, so every non-atomic field the
                // writer reads has to be brought across from the real System.
                // Without this the frame would report the state output_system
                // was copied from at set-up: the initial box, no constraint
                // diagnostics and no completed-step thermodynamics.
                //
                // Nothing here is reduced. Every field is either replicated (the
                // box, the constraint virial) or already global (the provider
                // virial, allreduced by the providers, and the completed-step
                // record, built from globally reduced quantities), so each rank
                // already holds the identical value. Reducing would multiply
                // them; tests/mpi_constraint_virial.cpp asserts the identity
                // across ranks rather than leaving it assumed.
                output_system.copy_frame_state_from(system);
                if (legacy) writer.write_frame(output_system, step, time, twice_ke, dof);
                if (frame) {
                    frame_writer.write_frame(output_system, step, time, twice_ke, dof);
                    status.frames_written = frame_writer.frame_count();
                }
            }
        };

        auto gather_global_system = [&]() {
            if (nprocs <= 1) {
                return system;
            }

            // Atomic numbers travel with the rest of the per-atom state: a
            // version 3 checkpoint stores them, and an ML restart needs them.
            constexpr std::size_t fields_per_atom = 12;
            std::vector<double> local(global_atom_count * fields_per_atom, 0.0);
            const auto atom_types = system.atom_types();
            const auto molecule_ids = system.molecule_ids();
            const auto masses = system.masses();
            const auto charges = system.charges();
            const auto coordinates = system.coordinates();
            const auto velocities = system.velocities();
            const auto atomic_numbers = system.atomic_numbers();
            for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
                const int tag = system.atom_tag(atom_index);
                if (tag < 0 || static_cast<std::size_t>(tag) >= global_atom_count) {
                    throw std::runtime_error("Local atom tag is outside the checkpoint atom map");
                }
                const std::size_t offset = static_cast<std::size_t>(tag) * fields_per_atom;
                local[offset] = static_cast<double>(tag);
                local[offset + 1] = static_cast<double>(atom_types[atom_index]);
                local[offset + 2] = static_cast<double>(molecule_ids[atom_index]);
                local[offset + 3] = masses[atom_index];
                local[offset + 4] = charges[atom_index];
                local[offset + 5] = coordinates[atom_index][0];
                local[offset + 6] = coordinates[atom_index][1];
                local[offset + 7] = coordinates[atom_index][2];
                local[offset + 8] = velocities[atom_index][0];
                local[offset + 9] = velocities[atom_index][1];
                local[offset + 10] = velocities[atom_index][2];
                local[offset + 11] = static_cast<double>(atomic_numbers[atom_index]);
            }

            std::vector<double> global;
            mpi_comm->allreduce_vector(local, global);
            gmd::System checkpoint_system;
            checkpoint_system.resize(global_atom_count);
            checkpoint_system.set_box(system.box());
            auto out_atom_types = checkpoint_system.mutable_atom_types();
            auto out_molecule_ids = checkpoint_system.mutable_molecule_ids();
            auto out_masses = checkpoint_system.mutable_masses();
            auto out_charges = checkpoint_system.mutable_charges();
            auto out_coordinates = checkpoint_system.mutable_coordinates();
            auto out_velocities = checkpoint_system.mutable_velocities();
            auto out_tags = checkpoint_system.mutable_atom_tags();
            auto out_atomic_numbers = checkpoint_system.mutable_atomic_numbers();
            for (std::size_t atom_index = 0; atom_index < global_atom_count; ++atom_index) {
                const std::size_t offset = atom_index * fields_per_atom;
                out_tags[atom_index] = static_cast<int>(global[offset]);
                out_atom_types[atom_index] = static_cast<int>(global[offset + 1]);
                out_molecule_ids[atom_index] = static_cast<int>(global[offset + 2]);
                out_masses[atom_index] = global[offset + 3];
                out_charges[atom_index] = global[offset + 4];
                out_coordinates[atom_index] = {
                    global[offset + 5],
                    global[offset + 6],
                    global[offset + 7]
                };
                out_velocities[atom_index] = {
                    global[offset + 8],
                    global[offset + 9],
                    global[offset + 10]
                };
                out_atomic_numbers[atom_index] = static_cast<int>(global[offset + 11]);
            }
            checkpoint_system.set_potential_energy(system.potential_energy());
            return checkpoint_system;
        };

        auto checkpoint_path_for = [&](std::uint64_t step) {
            return run_config.checkpoint_step_suffix
                       ? gmd::step_suffixed_checkpoint_path(run_config.checkpoint_file, step)
                       : run_config.checkpoint_file;
        };

        // Collective (the MPI gather). Returns true when a checkpoint was written.
        std::optional<std::uint64_t> last_checkpoint_step;
        auto write_checkpoint_if_needed = [&](std::uint64_t step, bool force) {
            if (run_config.checkpoint_file.empty()) {
                return false;
            }
            if (!force &&
                (run_config.write_checkpoint_every == 0 ||
                 step % run_config.write_checkpoint_every != 0)) {
                return false;
            }
            if (last_checkpoint_step == step) {
                return true;  // this exact state is already on disk
            }
            last_checkpoint_step = step;
            gmd::System checkpoint_system = gather_global_system();
            const std::filesystem::path checkpoint_path = checkpoint_path_for(step);
            status.checkpoint = absolute_string(checkpoint_path);
            if (!is_root_rank) {
                return true;
            }

            gmd::CheckpointMetadata metadata;
            metadata.step = step;
            metadata.time_fs = static_cast<double>(step) * run_config.time_step_fs;
            metadata.xyz_file = xyz_path.string();
            metadata.run_file = run_path.string();
            metadata.force_field_file = ff_path.has_value() ? ff_path->string() : "";
            metadata.topology_file = top_path.has_value() ? top_path->string() : "";
            metadata.velocity_seed = run_config.velocity_seed;
            metadata.force_field_summary = force_field_summary;
            metadata.config_summary =
                "dt_fs=" + std::to_string(run_config.time_step_fs) +
                " output_interval=" + std::to_string(run_config.output_interval) +
                " mpi_size=" + std::to_string(nprocs);
            metadata.thermostat_type = run_config.thermostat_type;
            metadata.thermostat_state =
                thermostat != nullptr ? thermostat->checkpoint_state() : "stateless";
            metadata.barostat_type = run_config.barostat_type;
            metadata.barostat_state =
                barostat != nullptr ? barostat->checkpoint_state() : "stateless";
            // Virial and pressure state, each named explicitly; see
            // CheckpointMetadata for what a restart does with each one.
            metadata.constraint_virial = system.constraint_virial();
            switch (system.constraint_virial_state()) {
                case gmd::ConstraintVirialState::Valid:
                    metadata.constraint_virial_state = "valid";
                    metadata.constraint_virial_time_level = "endpoint_rattle_t_plus_dt";
                    break;
                case gmd::ConstraintVirialState::Unavailable:
                    metadata.constraint_virial_state = "unavailable";
                    metadata.constraint_virial_time_level = "none";
                    break;
                case gmd::ConstraintVirialState::NotApplicable:
                    metadata.constraint_virial_state = "not_applicable";
                    metadata.constraint_virial_time_level = "none";
                    break;
            }
            metadata.provider_virial = system.provider_virial();
            metadata.provider_virial_valid = system.provider_virial_valid();
            const auto& completed = system.step_thermodynamics();
            metadata.step_pressure_valid = completed.valid;
            metadata.step_pressure = completed.pressure;
            metadata.step_pressure_virial = completed.virial;
            metadata.step_pressure_twice_ke = completed.twice_kinetic_energy;
            metadata.step_pressure_volume = completed.volume;
            metadata.step_pressure_potential_energy = completed.potential_energy;
            // Provenance (version 3).
            metadata.force_field_kind = force_field_kind;
            metadata.force_field_sha256 = force_field_fingerprint;
            metadata.model_path = run_config.force_field_type == "ml"
                                      ? run_config.ml_model_path.string()
                                      : std::string();
            metadata.model_sha256 = model_sha256;
            metadata.run_id = run_config.run_id;

            gmd::CheckpointData checkpoint{
                .metadata = metadata,
                .system = &checkpoint_system,
                .topology = topology.get(),
            };
            gmd::write_checkpoint(checkpoint_path, checkpoint);
            return true;
        };

        auto output_due = [&](std::uint64_t step) {
            return output_interval == 0 || step % output_interval == 0;
        };
        auto frame_due = [&](std::uint64_t step) {
            return frame_output_enabled && (frame_interval == 0 || step % frame_interval == 0);
        };

        // Write the first frame (t=0, or the restart state).
        const std::uint64_t start_step = simulation.current_step();
        status.start_step = start_step;
        status.step = start_step;
        write_global_frame(start_step, static_cast<double>(start_step) * run_config.time_step_fs,
                           true, frame_output_enabled);
        if (is_root_rank) {
            write_status_file(status, "running");
        }

        if (is_root_rank) {
            std::cout << log_prefix << "Running " << run_config.num_steps << " more steps with "
                      << global_atom_count << " atoms";
            if (!run_config.checkpoint_file.empty()) {
                std::cout << "  checkpoint=" << run_config.checkpoint_file.string()
                          << " every=" << run_config.write_checkpoint_every
                          << (run_config.checkpoint_step_suffix ? " (step-suffixed)" : "");
            }
            std::cout << "...\n";
            if (!run_config.stop_file.empty() || !run_config.checkpoint_request_file.empty()) {
                std::cout << log_prefix << "Control: polling every "
                          << run_config.control_poll_interval << " step(s)";
                if (!run_config.stop_file.empty()) {
                    std::cout << "  stop_file=" << run_config.stop_file.string();
                }
                if (!run_config.checkpoint_request_file.empty()) {
                    std::cout << "  checkpoint_request_file="
                              << run_config.checkpoint_request_file.string();
                }
                std::cout << "\n";
            }
        }

        // --- Cooperative control -------------------------------------------
        //
        // Polled after every control_poll_interval-th global step, once that
        // step's frame and periodic checkpoint are written. Rank 0 checks the
        // files; every rank checks its own signal flag; the decisions are
        // OR-reduced so all ranks act identically.
        struct ControlRequests {
            bool stop_file = false;
            bool signal = false;
            bool checkpoint = false;
        };
        auto poll_control = [&]() {
            ControlRequests requests;
            if (is_root_rank) {
                std::error_code ec;
                requests.stop_file = !run_config.stop_file.empty() &&
                                     std::filesystem::exists(run_config.stop_file, ec);
                requests.checkpoint =
                    !run_config.checkpoint_request_file.empty() &&
                    std::filesystem::exists(run_config.checkpoint_request_file, ec);
            }
            requests.signal = g_stop_signal != 0;
            if (nprocs > 1) {
                requests.stop_file = mpi_comm->allreduce_logical_or(requests.stop_file);
                requests.signal = mpi_comm->allreduce_logical_or(requests.signal);
                requests.checkpoint = mpi_comm->allreduce_logical_or(requests.checkpoint);
            }
            return requests;
        };
        auto acknowledge_checkpoint_request = [&]() {
            if (!is_root_rank || run_config.checkpoint_request_file.empty()) return;
            std::error_code ec;
            std::filesystem::remove(run_config.checkpoint_request_file, ec);
        };

        bool stopped = false;
        for (std::uint64_t s = 1; s <= run_config.num_steps; ++s) {
            simulation.step(runtime);
            const std::uint64_t global_step = simulation.current_step();
            status.step = global_step;
            write_global_frame(global_step, global_step * run_config.time_step_fs,
                               output_due(global_step), frame_due(global_step));
            write_checkpoint_if_needed(global_step, false);

            if (global_step % run_config.control_poll_interval != 0) {
                continue;
            }
            const ControlRequests requests = poll_control();
            if (requests.stop_file || requests.signal) {
                if (requests.stop_file) {
                    status.stop_reason = "stop_file";
                } else {
                    // The flag may be set on another rank only; name it there.
                    status.stop_reason = g_stop_signal != 0 ? signal_name(g_stop_signal)
                                                            : std::string("signal");
                }
                stopped = true;
                break;
            }
            if (requests.checkpoint) {
                write_checkpoint_if_needed(global_step, true);
                if (is_root_rank) {
                    std::cout << log_prefix << "Checkpoint written by request at step "
                              << global_step << ": " << status.checkpoint << "\n";
                    // Status first, then remove the request: once the request
                    // file is gone, the status file already names the checkpoint.
                    write_status_file(status, "running");
                }
                acknowledge_checkpoint_request();
            }
        }

        // Always write the final frame (unless this step's is already out).
        const std::uint64_t final_step = simulation.current_step();
        status.step = final_step;
        write_global_frame(final_step, final_step * run_config.time_step_fs,
                           last_legacy_frame_step != final_step,
                           frame_output_enabled && last_frame_file_step != final_step);
        write_checkpoint_if_needed(final_step, true);
        if (stopped) {
            // A checkpoint request that arrived with the stop is satisfied by
            // the final checkpoint just written.
            acknowledge_checkpoint_request();
        }

        if (is_root_rank) {
            if (stopped) {
                const std::string line = "stopped by request at step " +
                                         std::to_string(final_step) + " (" +
                                         status.stop_reason + ")";
                writer.write_log_comment(line);
                std::cout << log_prefix << "Run " << line << "\n";
            }
            writer.close();
            frame_writer.close();
            write_status_file(status, stopped ? "stopped" : "completed");

            std::cout << log_prefix << "Done. " << writer.frame_count() << " frames written.\n"
                      << log_prefix << "Final PE = " << system.potential_energy() << " eV\n";
        }

    } catch (const std::exception& error) {
        std::cerr << log_prefix << "Error: " << error.what() << "\n";
        if (!status.path.empty()) {
            try {
                write_status_file(status, "failed", error.what());
            } catch (const std::exception& status_error) {
                std::cerr << log_prefix << "Error: could not write status file: "
                          << status_error.what() << "\n";
            }
        }
        return 1;
    }

    return 0;
}
