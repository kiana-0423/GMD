#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "gmd/core/simulation.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/io/config_loader.hpp"
#include "gmd/io/trajectory_writer.hpp"
#include "gmd/neighbor/verlet_neighbor_builder.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

int main(int argc, char** argv)
{
    const char* xyz_path = argc > 1 ? argv[1] : "xyz.in";
    const char* run_path = argc > 2 ? argv[2] : "run.in";
    const char* ff_path  = argc > 3 ? argv[3] : nullptr;

    // Output stem: same directory as xyz input, base name "output".
    const std::filesystem::path output_stem =
        std::filesystem::path(xyz_path).parent_path() / "output";

    // Write a frame every this many steps (0 = only first and last).
    constexpr std::uint64_t output_interval = 100;

    try {
        gmd::ConfigLoader loader;
        gmd::System system;
        loader.load_xyz(xyz_path, system);
        const gmd::RunConfig run_config = loader.load_run(run_path);

        // --- Force provider ---
        std::shared_ptr<gmd::ClassicalForceProvider> force_provider;
        if (ff_path != nullptr && std::filesystem::exists(ff_path)) {
            const auto ff_config = loader.load_force_field(ff_path);
            force_provider = std::make_shared<gmd::ClassicalForceProvider>(ff_config);
            std::cout << "[gmd] Loaded force field from " << ff_path
                      << " (" << ff_config.elements.size() << " element type(s))\n";
        } else {
            force_provider = std::make_shared<gmd::ClassicalForceProvider>();
            std::cout << "[gmd] No force field file supplied; using default Ar LJ parameters.\n";
        }

        // --- Neighbor builder (r_cut from force provider, r_skin = 2.0 Å) ---
        constexpr double r_skin = 2.0;
        auto neighbor_builder = std::make_shared<gmd::VerletNeighborBuilder>(
            force_provider->cutoff(), r_skin);

        // --- Integrator ---
        auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(run_config.time_step);

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
        simulation.set_force_provider(force_provider);
        simulation.set_neighbor_builder(neighbor_builder);
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
            writer.write_frame_if(system, s, s * run_config.time_step, dof, output_interval);
        }

        // Always write final frame.
        if (run_config.num_steps % output_interval != 0) {
            writer.write_frame(system, run_config.num_steps,
                               run_config.num_steps * run_config.time_step, dof);
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
