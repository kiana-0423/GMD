#include <algorithm>
#include <array>
#include <cmath>
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

    try {
        gmd::ConfigLoader loader;
        gmd::System system;

        // Load run config first so inline force field is available before the xyz file is parsed.
        const gmd::RunConfig run_config = loader.load_run(run_path);
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
                              << " tolerance-equivalent duplicate(s)\n";
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
                std::cout << log_prefix << "Thermostat: Nose-Hoover  tau=" << run_config.thermostat_tau << " fs\n";
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
                          << " bar  tau=" << run_config.barostat_tau << " fs\n";
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
        }

        simulation.initialize(runtime);

        if (is_restart) {
            // Must FOLLOW initialize(): that call evaluates the forces at the
            // checkpointed coordinates and installs the provider virial, which
            // drops any constraint term attached to an earlier geometry.
            //
            // The constraint virial is the ENDPOINT RATTLE value of the step
            // that produced this state, so it belongs to exactly these
            // coordinates and is attached to that freshly recomputed provider
            // virial, so the restarted run reports for this frame the same
            // pressure the uninterrupted run reported.
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
        }


        // Degrees of freedom for every temperature the run reports. Read back
        // from the integrator rather than recomputed here, so the trajectory
        // log and the thermostat are guaranteed to use the same count: 3N,
        // less 3 when the COM velocity is removed, less one per constraint.
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

        auto write_global_frame = [&](std::uint64_t step, double time) {
            // compute_twice_ke already performs MPI_Allreduce internally to
            // return the global kinetic energy; do NOT double-wrap here.
            const double twice_ke = gmd::compute_twice_ke(system);
            if (is_root_rank) {
                if (nprocs == 1) {
                    writer.write_frame(system, step, time, twice_ke, dof);
                }
            }

            if (nprocs <= 1) {
                return;
            }

            std::vector<double> local_coordinates(global_atom_count * 3, 0.0);
            const auto coordinates = system.coordinates();
            for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
                const int tag = system.atom_tag(atom_index);
                if (tag < 0 || static_cast<std::size_t>(tag) >= global_atom_count) {
                    throw std::runtime_error("Local atom tag is outside the output coordinate map");
                }

                const std::size_t offset = static_cast<std::size_t>(tag) * 3;
                local_coordinates[offset] = coordinates[atom_index][0];
                local_coordinates[offset + 1] = coordinates[atom_index][1];
                local_coordinates[offset + 2] = coordinates[atom_index][2];
            }

            std::vector<double> global_coordinates;
            mpi_comm->allreduce_vector(local_coordinates, global_coordinates);
            if (is_root_rank) {
                auto frame_coordinates = output_system.mutable_coordinates();
                for (std::size_t atom_index = 0; atom_index < global_atom_count; ++atom_index) {
                    const std::size_t offset = atom_index * 3;
                    frame_coordinates[atom_index] = {
                        global_coordinates[offset],
                        global_coordinates[offset + 1],
                        global_coordinates[offset + 2]
                    };
                }
                output_system.set_potential_energy(system.potential_energy());
                writer.write_frame(output_system, step, time, twice_ke, dof);
            }
        };

        auto gather_global_system = [&]() {
            if (nprocs <= 1) {
                return system;
            }

            constexpr std::size_t fields_per_atom = 11;
            std::vector<double> local(global_atom_count * fields_per_atom, 0.0);
            const auto atom_types = system.atom_types();
            const auto molecule_ids = system.molecule_ids();
            const auto masses = system.masses();
            const auto charges = system.charges();
            const auto coordinates = system.coordinates();
            const auto velocities = system.velocities();
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
            }
            checkpoint_system.set_potential_energy(system.potential_energy());
            return checkpoint_system;
        };

        auto write_checkpoint_if_needed = [&](std::uint64_t step, bool force) {
            if (run_config.checkpoint_file.empty()) {
                return;
            }
            if (!force &&
                (run_config.write_checkpoint_every == 0 ||
                 step % run_config.write_checkpoint_every != 0)) {
                return;
            }
            gmd::System checkpoint_system = gather_global_system();
            if (!is_root_rank) {
                return;
            }

            gmd::CheckpointMetadata metadata;
            metadata.step = step;
            metadata.time_fs = static_cast<double>(step) * run_config.time_step_fs;
            metadata.xyz_file = xyz_path.string();
            metadata.run_file = run_path.string();
            metadata.force_field_file = ff_path.has_value() ? ff_path->string() : "";
            metadata.topology_file = top_path.has_value() ? top_path->string() : "";
            metadata.velocity_seed = run_config.velocity_seed;
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
                metadata.force_field_summary = ff_summary.str();
            } else if (ff_path.has_value()) {
                metadata.force_field_summary = "file " + ff_path->string();
            } else {
                metadata.force_field_summary = "default_lj";
            }
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

            // The constraint contribution and which multiplier it is; see
            // CheckpointMetadata for what a restart does with it.
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
            gmd::CheckpointData checkpoint{
                .metadata = metadata,
                .system = &checkpoint_system,
                .topology = topology.get(),
            };
            gmd::write_checkpoint(run_config.checkpoint_file, checkpoint);
        };

        // Write t=0 frame.
        const std::uint64_t start_step = simulation.current_step();
        write_global_frame(start_step, static_cast<double>(start_step) * run_config.time_step_fs);

        if (is_root_rank) {
            std::cout << log_prefix << "Running " << run_config.num_steps << " more steps with "
                      << global_atom_count << " atoms";
            if (!run_config.checkpoint_file.empty()) {
                std::cout << "  checkpoint=" << run_config.checkpoint_file.string()
                          << " every=" << run_config.write_checkpoint_every;
            }
            std::cout << "...\n";
        }

        for (std::uint64_t s = 1; s <= run_config.num_steps; ++s) {
            simulation.step(runtime);
            const std::uint64_t global_step = simulation.current_step();
            if (output_interval == 0 || global_step % output_interval == 0) {
                write_global_frame(global_step, global_step * run_config.time_step_fs);
            }
            write_checkpoint_if_needed(global_step, false);
        }

        // Always write final frame.
        const std::uint64_t final_step = simulation.current_step();
        if (output_interval != 0 && final_step % output_interval != 0) {
            write_global_frame(final_step,
                               final_step * run_config.time_step_fs);
        }
        write_checkpoint_if_needed(final_step, true);

        if (is_root_rank) {
            writer.close();

            std::cout << log_prefix << "Done. " << writer.frame_count() << " frames written.\n"
                      << log_prefix << "Final PE = " << system.potential_energy() << " eV\n";
        }

    } catch (const std::exception& error) {
        std::cerr << log_prefix << "Error: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
