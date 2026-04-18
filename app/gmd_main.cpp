#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <stdexcept>

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
#include "gmd/io/config_loader.hpp"
#include "gmd/io/trajectory_writer.hpp"
#include "gmd/neighbor/verlet_neighbor_builder.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

namespace {

enum class ForceFieldFileKind {
    LJ,
    Molecular,
};

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
    const char* xyz_path = argc > 1 ? argv[1] : "xyz.in";
    const char* run_path = argc > 2 ? argv[2] : "run.in";
    const char* ff_path  = argc > 3 ? argv[3] : nullptr;
    const char* top_path = argc > 4 ? argv[4] : nullptr;

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
            std::cout << "[gmd] Loaded inline force field from " << run_path
                      << " (" << ff_config.elements.size() << " element type(s))\n";
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

            std::cout << "[gmd] Loaded force field from " << ff_path
                      << " (" << molecular_ff->lj.elements.size() << " atom type(s), "
                      << topology->bonds.size() << " bond(s), "
                      << topology->angles.size() << " angle(s), "
                      << topology->dihedrals.size() << " dihedral(s), "
                      << topology->impropers.size() << " improper(s))\n";
            std::cout << "[gmd] Loaded topology from " << top_path << "\n";
            if (run_config.molecular_nonbonded_mode == "lj_unsafe") {
                lj_provider = std::make_shared<gmd::ClassicalForceProvider>(molecular_ff->lj);
                auto composite = std::make_shared<gmd::CompositeForceProvider>();
                composite->add(lj_provider);
                composite->add(bonded_provider);
                active_provider = composite;
                short_range_cutoff = lj_provider->cutoff();
                need_neighbor_builder = true;
                std::cerr << "[gmd] WARNING: molecular_nonbonded=lj_unsafe enables LJ without 1-2/1-3 exclusions.\n"
                          << "[gmd]          This is unphysical for most molecular force fields and is intended\n"
                          << "[gmd]          only for diagnostics until exclusion lists are implemented.\n";
            } else {
                std::cout << "[gmd] Molecular non-bonded mode: bonded-only (default).\n"
                          << "[gmd] Set 'molecular_nonbonded lj_unsafe' in run.in to explicitly enable LJ.\n";
            }
        } else if (external_lj_ff.has_value()) {
            const auto& ff_config = external_lj_ff.value();
            lj_provider = std::make_shared<gmd::ClassicalForceProvider>(ff_config);
            active_provider = lj_provider;
            short_range_cutoff = lj_provider->cutoff();
            need_neighbor_builder = true;
            std::cout << "[gmd] Loaded force field from " << ff_path
                      << " (" << ff_config.elements.size() << " element type(s))\n";
        } else {
            lj_provider = std::make_shared<gmd::ClassicalForceProvider>();
            active_provider = lj_provider;
            short_range_cutoff = lj_provider->cutoff();
            need_neighbor_builder = true;
            std::cout << "[gmd] No force field supplied; using default Ar LJ parameters.\n";
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
                std::cout << "[gmd] Coulomb: PME  order=" << cc.pme_order
                          << "  grid=" << cc.pme_grid[0] << "x"
                          << cc.pme_grid[1] << "x" << cc.pme_grid[2] << "\n";
            } else {
                // Default to Ewald.
                auto ewald = std::make_shared<gmd::EwaldForceProvider>(
                    cc.alpha, cc.kmax, cc.real_cutoff);
                composite->add(ewald);
                std::cout << "[gmd] Coulomb: Ewald  alpha=" << cc.alpha
                          << "  kmax=" << cc.kmax
                          << "  r_cut=" << cc.real_cutoff << "\n";
            }
        }

        // --- Neighbor builder (uses the largest active short-range cutoff, r_skin = 2.0 Å) ---
        std::shared_ptr<gmd::VerletNeighborBuilder> neighbor_builder;
        if (need_neighbor_builder) {
            constexpr double r_skin = 2.0;
            neighbor_builder = std::make_shared<gmd::VerletNeighborBuilder>(
                short_range_cutoff, r_skin);
        }

        // --- Integrator ---
        auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(run_config.time_step);
        integrator->set_target_temperature(run_config.temperature);

        // --- Thermostat ---
        if (run_config.thermostat_type == "nose_hoover") {
            auto tstat = std::make_shared<gmd::NoseHooverThermostat>(run_config.thermostat_tau);
            integrator->set_thermostat(tstat);
            std::cout << "[gmd] Thermostat: Nose-Hoover  tau=" << run_config.thermostat_tau << " fs\n";
        } else if (run_config.thermostat_type == "velocity_rescaling") {
            auto tstat = std::make_shared<gmd::VelocityRescalingThermostat>();
            integrator->set_thermostat(tstat);
            std::cout << "[gmd] Thermostat: velocity rescaling\n";
        }

        // --- Barostat ---
        if (run_config.barostat_type == "berendsen") {
            auto bstat = std::make_shared<gmd::BerendsenBarostat>(
                run_config.barostat_tau, run_config.compressibility);
            integrator->set_barostat(bstat);
            integrator->set_target_pressure(run_config.target_pressure);
            std::cout << "[gmd] Barostat: Berendsen  P=" << run_config.target_pressure
                      << " bar  tau=" << run_config.barostat_tau << " fs\n";
        } else if (run_config.barostat_type == "monte_carlo") {
            auto bstat = std::make_shared<gmd::MCBarostat>(
                run_config.mc_frequency,
                run_config.mc_volume_step);
            integrator->set_barostat(bstat);
            integrator->set_target_pressure(run_config.target_pressure);
            std::cout << "[gmd] Barostat: Monte Carlo NPT  P=" << run_config.target_pressure
                      << " bar  freq=" << run_config.mc_frequency
                      << "  max_delta_ln_V=" << run_config.mc_volume_step << "\n";
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
        if (neighbor_builder) {
            simulation.set_neighbor_builder(neighbor_builder);
        }
        simulation.set_integrator(integrator);
        simulation.set_time_step(run_config.time_step);

        // --- Trajectory writer ---
        gmd::TrajectoryWriter writer;
        writer.open(output_stem);
        std::cout << "[gmd] Writing trajectory to " << output_stem.string() << ".xyz"
                  << " and energy log to " << output_stem.string() << ".log\n";

        // Degrees of freedom = 3N - 3 (after COM velocity removal).
        const std::size_t dof = system.atom_count() > 1 ? 3 * system.atom_count() - 3 : 3;

        gmd::RuntimeContext runtime;
        simulation.initialize(runtime);

        // Write t=0 frame.
        writer.write_frame(system, 0, 0.0, dof);

        std::cout << "[gmd] Running " << run_config.num_steps << " steps with "
                  << system.atom_count() << " atoms...\n";

        for (std::uint64_t s = 1; s <= run_config.num_steps; ++s) {
            simulation.step(runtime);
            writer.write_frame_if(system, s, s * run_config.time_step_fs, dof, output_interval);
        }

        // Always write final frame.
        if (run_config.num_steps % output_interval != 0) {
            writer.write_frame(system, run_config.num_steps,
                               run_config.num_steps * run_config.time_step_fs, dof);
        }

        writer.close();

        std::cout << "[gmd] Done. " << writer.frame_count() << " frames written.\n"
                  << "      Final PE = " << system.potential_energy() << " eV\n";

    } catch (const std::exception& error) {
        std::cerr << "[gmd] Error: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
