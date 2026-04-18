// =============================================================================
// Ethane demo: bonded molecular force field run
//
// Demonstrates:
//  1. load_molecular_ff()  — reads ethane.ff (mass, LJ, bond/angle/dihedral params)
//  2. load_topology()      — reads ethane.top (bond/angle/dihedral connectivity)
//  3. BondedForceProvider  — computes intramolecular bonded energy & forces
//  4. ClassicalForceProvider (LJ non-bonded)
//  5. CompositeForceProvider — sums both; total energy printed every step
//
// Build (from project root):
//   see run_ethane_demo.sh
// =============================================================================

#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "gmd/core/simulation.hpp"
#include "gmd/force/bonded_force_provider.hpp"
#include "gmd/force/classical_force_provider.hpp"
#include "gmd/force/composite_force_provider.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/io/config_loader.hpp"
#include "gmd/io/trajectory_writer.hpp"
#include "gmd/neighbor/verlet_neighbor_builder.hpp"
#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

namespace {

// Print a banner line to stdout.
void banner(const std::string& msg) {
    std::cout << "\n=== " << msg << " ===\n";
}

// Convert eV to kcal/mol for readable output.
constexpr double eV_to_kcal = 23.0609;

}  // namespace

int main(int argc, char** argv)
{
    // Default paths relative to the directory containing this binary.
    // Override via command-line arguments:
    //   ethane_demo  [xyz]  [run]  [ff]  [top]  [out_stem]
    const std::string xyz_path = argc > 1 ? argv[1] : "ethane.xyz";
    const std::string run_path = argc > 2 ? argv[2] : "ethane.run";
    const std::string ff_path  = argc > 3 ? argv[3] : "ethane.ff";
    const std::string top_path = argc > 4 ? argv[4] : "ethane.top";
    const std::string out_stem = argc > 5 ? argv[5] : "ethane_out";

    // Print step/energy every this many MD steps.
    constexpr std::uint64_t print_interval = 10;

    try {
        banner("Loading inputs");

        gmd::ConfigLoader loader;

        // 1. Load molecular force field -------------------------------------
        const gmd::MolecularForceFieldConfig mff =
            loader.load_molecular_ff(ff_path);

        std::cout << "  FF file       : " << ff_path << "\n"
                  << "  Atom types    : " << mff.lj.elements.size() << "\n"
                  << "  Bond types    : " << mff.bond_types.size() << "\n"
                  << "  Angle types   : " << mff.angle_types.size() << "\n"
                  << "  Dihedral types: " << mff.dihedral_types.size() << "\n"
                  << "  Improper types: " << mff.improper_types.size() << "\n"
                  << "  LJ cutoff     : " << mff.lj.cutoff << " Å\n";

        // 2. Load topology --------------------------------------------------
        const std::shared_ptr<gmd::Topology> topo =
            loader.load_topology(top_path);

        std::cout << "  Topology file : " << top_path << "\n"
                  << "  Bonds         : " << topo->bonds.size() << "\n"
                  << "  Angles        : " << topo->angles.size() << "\n"
                  << "  Dihedrals     : " << topo->dihedrals.size() << "\n"
                  << "  Impropers     : " << topo->impropers.size() << "\n";

        // 3. Load run parameters + XYZ structure ----------------------------
        const gmd::RunConfig run_config = loader.load_run(run_path);

        gmd::System system;
        loader.load_xyz(xyz_path, system, &mff.lj);

        std::cout << "  Atoms         : " << system.atom_count() << "\n"
                  << "  Box           : "
                  << system.box().lengths[0] << " x "
                  << system.box().lengths[1] << " x "
                  << system.box().lengths[2] << " Å\n"
                  << "  Steps         : " << run_config.num_steps << "\n"
                  << "  dt            : " << run_config.time_step_fs << " fs\n"
                  << "  T_target      : " << run_config.temperature << " K\n";

        // Print per-type masses
        std::cout << "\n  Atom type masses:\n";
        for (std::size_t t = 0; t < mff.lj.elements.size(); ++t) {
            const auto& ep = mff.lj.elements[t];
            std::cout << "    type " << t + 1 << " (" << ep.element << ")"
                      << "  mass=" << ep.mass << " amu"
                      << "  ε=" << ep.epsilon / gmd::kcal_per_mol_to_eV << " kcal/mol"
                      << "  σ=" << ep.sigma << " Å\n";
        }

        // Print bond type parameters
        std::cout << "\n  Bond type parameters:\n";
        for (std::size_t i = 0; i < mff.bond_types.size(); ++i) {
            const auto& bp = mff.bond_types[i];
            std::cout << "    bond_type " << i + 1
                      << "  k=" << bp.k / gmd::kcal_per_mol_to_eV << " kcal/mol/Å²"
                      << "  r0=" << bp.r0 << " Å\n";
        }

        // Print angle type parameters
        std::cout << "\n  Angle type parameters:\n";
        for (std::size_t i = 0; i < mff.angle_types.size(); ++i) {
            const auto& ap = mff.angle_types[i];
            std::cout << "    angle_type " << i + 1
                      << "  k=" << ap.k / gmd::kcal_per_mol_to_eV << " kcal/mol/rad²"
                      << "  θ0=" << ap.theta0 / gmd::deg_to_rad << "°\n";
        }

        // Print dihedral type parameters
        std::cout << "\n  Dihedral type parameters:\n";
        for (std::size_t i = 0; i < mff.dihedral_types.size(); ++i) {
            const auto& dp = mff.dihedral_types[i];
            std::cout << "    dihedral_type " << i + 1
                      << "  k=" << dp.k / gmd::kcal_per_mol_to_eV << " kcal/mol"
                      << "  n=" << dp.n
                      << "  δ=" << dp.delta / gmd::deg_to_rad << "°\n";
        }

        // 4. Build BondedForceProvider --------------------------------------
        banner("Building force providers");

        auto bonded = std::make_shared<gmd::BondedForceProvider>(topo);
        for (const auto& bp : mff.bond_types)     bonded->add_bond_type(bp);
        for (const auto& ap : mff.angle_types)    bonded->add_angle_type(ap);
        for (const auto& dp : mff.dihedral_types) bonded->add_dihedral_type(dp);
        for (const auto& ip : mff.improper_types) bonded->add_improper_type(ip);

        std::cout << "  BondedForceProvider ready\n";

        // 5. Build LJ (non-bonded) provider ---------------------------------
        auto lj = std::make_shared<gmd::ClassicalForceProvider>(mff.lj);
        std::cout << "  ClassicalForceProvider (LJ) ready  cutoff=" << mff.lj.cutoff << " Å\n";

        // 6. Composite = LJ + bonded ----------------------------------------
        auto composite = std::make_shared<gmd::CompositeForceProvider>();
        composite->add(lj);
        composite->add(bonded);
        std::cout << "  CompositeForceProvider: LJ + bonded\n";

        // 7. Neighbor builder -----------------------------------------------
        constexpr double r_skin = 2.0;
        auto neighbor_builder =
            std::make_shared<gmd::VerletNeighborBuilder>(mff.lj.cutoff, r_skin);

        // 8. Integrator + Nosé-Hoover thermostat ----------------------------
        auto integrator =
            std::make_shared<gmd::VelocityVerletIntegrator>(run_config.time_step);
        integrator->set_target_temperature(run_config.temperature);

        auto thermostat =
            std::make_shared<gmd::NoseHooverThermostat>(run_config.thermostat_tau);
        integrator->set_thermostat(thermostat);

        // 9. Assemble simulation --------------------------------------------
        gmd::Simulation simulation(&system);
        simulation.set_force_provider(composite);
        simulation.set_neighbor_builder(neighbor_builder);
        simulation.set_integrator(integrator);
        simulation.set_time_step(run_config.time_step);
        simulation.set_initial_temperature(run_config.temperature);
        auto vel_init = std::make_shared<gmd::VelocityInitializer>(run_config.velocity_seed);
        simulation.set_velocity_initializer(vel_init);
        simulation.set_remove_center_of_mass_velocity(
            run_config.remove_center_of_mass_velocity);

        // 10. Trajectory writer --------------------------------------------
        gmd::TrajectoryWriter writer;
        writer.open(out_stem);
        std::cout << "\n  Output trajectory : " << out_stem << ".xyz\n"
                  << "  Output energy log : " << out_stem << ".log\n";

        // 11. Initialize (t=0 force evaluation) ----------------------------
        gmd::RuntimeContext runtime;
        simulation.initialize(runtime);

        const std::size_t dof =
            system.atom_count() > 1 ? 3 * system.atom_count() - 3 : 3;

        writer.write_frame(system, 0, 0.0, dof);

        // 12. Print t=0 state ----------------------------------------------
        banner("Initial state");

        // Compute initial bonded energies separately for diagnostics
        {
            gmd::ForceResult bonded_result;
            gmd::ForceRequest req{&system, &system.box(), 0, 0.0,
                                  system.coordinates(), nullptr};
            bonded->compute(req, bonded_result, runtime);

            gmd::ForceResult lj_result;
            lj->compute(req, lj_result, runtime);

            std::cout << std::fixed << std::setprecision(6)
                      << "  Bonded PE  = " << bonded_result.potential_energy
                      << " eV  ("
                      << bonded_result.potential_energy * eV_to_kcal << " kcal/mol)\n"
                      << "  LJ PE      = " << lj_result.potential_energy
                      << " eV  ("
                      << lj_result.potential_energy * eV_to_kcal << " kcal/mol)\n"
                      << "  Total PE   = " << system.potential_energy()
                      << " eV  ("
                      << system.potential_energy() * eV_to_kcal << " kcal/mol)\n";
        }

        // 13. Run MD --------------------------------------------------------
        banner("Running MD  (" + std::to_string(run_config.num_steps) + " steps)");

        std::cout << std::left
                  << std::setw(8)  << "Step"
                  << std::setw(12) << "Time(fs)"
                  << std::setw(16) << "PE(kcal/mol)"
                  << std::setw(16) << "KE(kcal/mol)"
                  << std::setw(16) << "E_tot(kcal/mol)"
                  << std::setw(10) << "T(K)"
                  << "\n"
                  << std::string(78, '-') << "\n";

        auto print_row = [&](std::uint64_t s) {
            const double pe   = system.potential_energy();
            // Kinetic energy via VelocityInitializer helper
            const double ke   = vel_init->kinetic_energy(system);
            const double etot = pe + ke;
            const double temp = dof > 0
                ? ke * 2.0 / (static_cast<double>(dof) * 8.617333e-5) : 0.0;

            std::cout << std::fixed << std::setprecision(3)
                      << std::left
                      << std::setw(8)  << s
                      << std::setw(12) << s * run_config.time_step_fs
                      << std::setw(16) << pe   * eV_to_kcal
                      << std::setw(16) << ke   * eV_to_kcal
                      << std::setw(16) << etot * eV_to_kcal
                      << std::setw(10) << temp
                      << "\n";
        };

        print_row(0);

        for (std::uint64_t s = 1; s <= run_config.num_steps; ++s) {
            simulation.step(runtime);
            writer.write_frame_if(system, s, s * run_config.time_step_fs, dof, 50);
            if (s % print_interval == 0) print_row(s);
        }

        // Always print final state
        if (run_config.num_steps % print_interval != 0) {
            print_row(run_config.num_steps);
        }

        if (run_config.num_steps % 50 != 0) {
            writer.write_frame(system, run_config.num_steps,
                               run_config.num_steps * run_config.time_step_fs, dof);
        }
        writer.close();

        banner("Done");
        std::cout << "  Frames written : " << writer.frame_count() << "\n"
                  << "  Final PE       : "
                  << std::fixed << std::setprecision(6)
                  << system.potential_energy() << " eV  ("
                  << system.potential_energy() * eV_to_kcal << " kcal/mol)\n";

    } catch (const std::exception& e) {
        std::cerr << "[ethane_demo] Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
