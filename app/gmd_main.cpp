#include <algorithm>
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
#include "gmd/integrator/mc_barostat.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/velocity_rescaling_thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/integrator/thermostat.hpp"
#include "gmd/io/config_loader.hpp"
#include "gmd/io/trajectory_writer.hpp"
#include "gmd/ml/ml_force_provider.hpp"
#include "gmd/neighbor/verlet_neighbor_builder.hpp"
#include "gmd/parallel/domain_decomposition.hpp"
#include "gmd/parallel/mpi_communicator.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"
#ifdef GMD_ENABLE_TORCH
#include "gmd/ml/torchscript_adapter.hpp"
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
};

CommandLine parse_command_line(int argc, char** argv) {
    std::vector<std::string> positionals;
    CommandLine cli;

    for (int arg_index = 1; arg_index < argc; ++arg_index) {
        const std::string arg = argv[arg_index];
        if (arg == "--np") {
            if (arg_index + 1 >= argc) {
                throw std::runtime_error("--np requires a positive process count");
            }

            const std::string count_arg = argv[++arg_index];
            std::size_t parsed_chars = 0;
            const int count = std::stoi(count_arg, &parsed_chars);
            if (parsed_chars != count_arg.size() || count <= 0) {
                throw std::runtime_error("--np requires a positive integer process count");
            }
            cli.expected_process_count = count;
            continue;
        }
        if (arg.starts_with("--")) {
            throw std::runtime_error("Unknown command-line option: " + arg);
        }

        positionals.push_back(arg);
    }

    if (positionals.size() > 4) {
        throw std::runtime_error(
            "Usage: gmd input.xyz run.in [ff.ff] [top.top] [--np N]");
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

int owner_rank_for_x(const gmd::Box& box, const gmd::System::Vec3& position, int nprocs) {
    double x = std::fmod(position[0], box.lengths[0]);
    if (x < 0.0) {
        x += box.lengths[0];
    }

    const double domain_width = box.lengths[0] / static_cast<double>(nprocs);
    const int owner = static_cast<int>(x / domain_width);
    return std::min(owner, nprocs - 1);
}

void keep_rank_local_atoms(gmd::System& system, int my_rank, int nprocs) {
    const auto box = system.box();
    const auto masses = system.masses();
    const auto charges = system.charges();
    const auto atom_types = system.atom_types();
    const auto atomic_numbers = system.atomic_numbers();
    const auto coordinates = system.coordinates();
    const auto velocities = system.velocities();

    struct LocalAtom {
        double mass;
        double charge;
        int atom_type;
        int atomic_number;
        int tag;
        gmd::System::Vec3 position;
        gmd::System::Vec3 velocity;
    };

    std::vector<LocalAtom> local_atoms;
    local_atoms.reserve(system.atom_count());
    for (std::size_t atom_index = 0; atom_index < system.atom_count(); ++atom_index) {
        if (owner_rank_for_x(box, coordinates[atom_index], nprocs) != my_rank) {
            continue;
        }

        local_atoms.push_back(LocalAtom{
            .mass = masses[atom_index],
            .charge = charges[atom_index],
            .atom_type = atom_types[atom_index],
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

    const char* xyz_path = cli.xyz_path.c_str();
    const char* run_path = cli.run_path.c_str();
    const char* ff_path = cli.ff_path.has_value() ? cli.ff_path->c_str() : nullptr;
    const char* top_path = cli.top_path.has_value() ? cli.top_path->c_str() : nullptr;

    // Output stem: same directory as xyz input, base name "output".
    const std::filesystem::path output_stem =
        std::filesystem::path(xyz_path).parent_path() / "output";

    // Write a frame every this many steps (0 = only first and last).
    constexpr std::uint64_t output_interval = 100;

    try {
        gmd::ConfigLoader loader;
        gmd::System system;

        // Load run config first so inline force field is available before the xyz file is parsed.
        const gmd::RunConfig run_config = loader.load_run(run_path);

        std::optional<gmd::LJForceFieldConfig> external_lj_ff;
        std::optional<gmd::MolecularForceFieldConfig> molecular_ff;
        std::shared_ptr<gmd::Topology> topology;

        // Resolve external FF before xyz loading so type→mass mapping works.
        const gmd::LJForceFieldConfig* xyz_ff = nullptr;
        if (run_config.force_field.has_value()) {
            xyz_ff = &run_config.force_field.value();
        } else if (ff_path != nullptr && std::filesystem::exists(ff_path)) {
            const auto ff_kind = detect_force_field_file_kind(ff_path);
            if (ff_kind == ForceFieldFileKind::Molecular) {
                if (top_path == nullptr) {
                    throw std::runtime_error(
                        "Molecular force fields require a topology file as the fourth CLI argument");
                }
                molecular_ff = loader.load_molecular_ff(ff_path);
                topology = loader.load_topology(top_path);
                xyz_ff = &molecular_ff->lj;
            } else {
                external_lj_ff = loader.load_force_field(ff_path);
                xyz_ff = &external_lj_ff.value();
            }
        }

        loader.load_xyz(xyz_path, system, xyz_ff);

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
                std::cout << log_prefix << "Loaded force field from " << ff_path
                          << " (" << molecular_ff->lj.elements.size() << " atom type(s), "
                          << topology->bonds.size() << " bond(s), "
                          << topology->angles.size() << " angle(s), "
                          << topology->dihedrals.size() << " dihedral(s), "
                          << topology->impropers.size() << " improper(s))\n";
                std::cout << log_prefix << "Loaded topology from " << top_path << "\n";
            }
            if (run_config.molecular_nonbonded_mode == "lj_unsafe") {
                lj_provider = std::make_shared<gmd::ClassicalForceProvider>(molecular_ff->lj);
                auto composite = std::make_shared<gmd::CompositeForceProvider>();
                composite->add(lj_provider);
                composite->add(bonded_provider);
                active_provider = composite;
                short_range_cutoff = lj_provider->cutoff();
                need_neighbor_builder = true;
                if (is_root_rank) {
                    std::cerr << log_prefix << "WARNING: molecular_nonbonded=lj_unsafe enables LJ without 1-2/1-3 exclusions.\n"
                              << log_prefix << "         This is unphysical for most molecular force fields and is intended\n"
                              << log_prefix << "         only for diagnostics until exclusion lists are implemented.\n";
                }
            } else {
                if (is_root_rank) {
                    std::cout << log_prefix << "Molecular non-bonded mode: bonded-only (default).\n"
                              << log_prefix << "Set 'molecular_nonbonded lj_unsafe' in run.in to explicitly enable LJ.\n";
                }
            }
        } else if (external_lj_ff.has_value()) {
            const auto& ff_config = external_lj_ff.value();
            lj_provider = std::make_shared<gmd::ClassicalForceProvider>(ff_config);
            active_provider = lj_provider;
            short_range_cutoff = lj_provider->cutoff();
            need_neighbor_builder = true;
            if (is_root_rank) {
                std::cout << log_prefix << "Loaded force field from " << ff_path
                          << " (" << ff_config.elements.size() << " element type(s))\n";
            }
        } else if (run_config.force_field_type == "ml") {
            if (nprocs > 1) {
                throw std::runtime_error(
                    "ML force provider is not yet supported with MPI domain decomposition");
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
                    cc.alpha, cc.real_cutoff, cc.pme_order, cc.pme_grid);
                composite->add(pme);
                if (is_root_rank) {
                    std::cout << log_prefix << "Coulomb: PME  order=" << cc.pme_order
                              << "  grid=" << cc.pme_grid[0] << "x"
                              << cc.pme_grid[1] << "x" << cc.pme_grid[2] << "\n";
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
            domain_decomposition->create_1d_decomposition(
                system.box(), nprocs, my_rank, short_range_cutoff, r_skin);
            keep_rank_local_atoms(system, my_rank, nprocs);
            std::cout << log_prefix << "Domain decomposition owns "
                      << system.num_local_atoms() << " of " << global_atom_count
                      << " atoms before ghost exchange\n";
        }

        // --- Integrator ---
        auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(run_config.time_step);
        integrator->set_target_temperature(run_config.temperature);

        // --- Thermostat ---
        if (run_config.thermostat_type == "nose_hoover") {
            auto tstat = std::make_shared<gmd::NoseHooverThermostat>(run_config.thermostat_tau);
            integrator->set_thermostat(tstat);
            if (is_root_rank) {
                std::cout << log_prefix << "Thermostat: Nose-Hoover  tau=" << run_config.thermostat_tau << " fs\n";
            }
        } else if (run_config.thermostat_type == "velocity_rescaling") {
            auto tstat = std::make_shared<gmd::VelocityRescalingThermostat>();
            integrator->set_thermostat(tstat);
            if (is_root_rank) {
                std::cout << log_prefix << "Thermostat: velocity rescaling\n";
            }
        }

        // --- Barostat ---
        if (run_config.barostat_type == "berendsen") {
            auto bstat = std::make_shared<gmd::BerendsenBarostat>(
                run_config.barostat_tau, run_config.compressibility);
            integrator->set_barostat(bstat);
            integrator->set_target_pressure(run_config.target_pressure);
            if (is_root_rank) {
                std::cout << log_prefix << "Barostat: Berendsen  P=" << run_config.target_pressure
                          << " bar  tau=" << run_config.barostat_tau << " fs\n";
            }
        } else if (run_config.barostat_type == "monte_carlo") {
            if (nprocs > 1) {
                throw std::runtime_error(
                    "Monte Carlo barostat is not yet supported with MPI domain decomposition");
            }
            auto bstat = std::make_shared<gmd::MCBarostat>(
                run_config.mc_frequency,
                run_config.mc_volume_step);
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
        simulation.set_velocity_initializer(velocity_initializer);
        simulation.set_velocity_init_mode(velocity_mode);
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

        // Degrees of freedom = 3N - 3 (after COM velocity removal).
        const std::size_t dof = global_atom_count > 1 ? 3 * global_atom_count - 3 : 3;

        simulation.initialize(runtime);

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

        // Write t=0 frame.
        write_global_frame(0, 0.0);

        if (is_root_rank) {
            std::cout << log_prefix << "Running " << run_config.num_steps << " steps with "
                      << global_atom_count << " atoms...\n";
        }

        for (std::uint64_t s = 1; s <= run_config.num_steps; ++s) {
            simulation.step(runtime);
            if (output_interval == 0 || s % output_interval == 0) {
                write_global_frame(s, s * run_config.time_step_fs);
            }
        }

        // Always write final frame.
        if (output_interval != 0 && run_config.num_steps % output_interval != 0) {
            write_global_frame(run_config.num_steps,
                               run_config.num_steps * run_config.time_step_fs);
        }

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
