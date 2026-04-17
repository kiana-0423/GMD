#include <iostream>
#include <memory>
#include <stdexcept>

#include "gmd/core/simulation.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/io/config_loader.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

int main(int argc, char** argv)
{
    const char* xyz_path = argc > 1 ? argv[1] : "xyz.in";
    const char* run_path = argc > 2 ? argv[2] : "run.in";

    try {
        gmd::ConfigLoader loader;
        gmd::System system;
        loader.load_xyz(xyz_path, system);
        const gmd::RunConfig run_config = loader.load_run(run_path);

        const auto velocity_mode =
            run_config.velocity_init_mode == "input" ? gmd::VelocityInitMode::FromInput
                                                      : gmd::VelocityInitMode::Random;

        auto velocity_initializer = std::make_shared<gmd::VelocityInitializer>(run_config.velocity_seed);
        auto force_provider = std::make_shared<gmd::ClassicalForceProvider>();
        auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(run_config.time_step);

        gmd::Simulation simulation(&system);
        simulation.set_velocity_initializer(velocity_initializer);
        simulation.set_velocity_init_mode(velocity_mode);
        simulation.set_remove_center_of_mass_velocity(run_config.remove_center_of_mass_velocity);
        simulation.set_initial_temperature(run_config.temperature);
        simulation.set_force_provider(force_provider);
        simulation.set_integrator(integrator);
        simulation.set_time_step(run_config.time_step);

        gmd::RuntimeContext runtime;
        simulation.initialize(runtime);
        simulation.run(runtime, run_config.num_steps);

        std::cout << "Loaded " << system.atom_count() << " atoms from " << xyz_path
                  << ", initialized velocities at " << run_config.temperature << " K"
                  << " using mode " << run_config.velocity_init_mode
                  << ", and completed " << run_config.num_steps << " integration steps" << std::endl;
    } catch (const std::exception& error) {
        std::cerr << "Initialization failed: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}
