#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "gmd/core/runtime_context.hpp"
#include "gmd/core/simulation.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/integrator/mc_barostat.hpp"
#include "gmd/integrator/nose_hoover_thermostat.hpp"
#include "gmd/integrator/velocity_verlet_integrator.hpp"
#include "gmd/io/checkpoint.hpp"
#include "gmd/system/system.hpp"

namespace {

class ZeroForceProvider final : public gmd::ForceProvider {
public:
    std::string_view name() const noexcept override { return "zero_force"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request,
                 gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.success = true;
        result.potential_energy = 0.0;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
        result.virial_valid = true;
    }
};

void check(bool condition, const std::string& message, int& failures) {
    if (!condition) {
        std::cerr << "[checkpoint] " << message << "\n";
        ++failures;
    }
}

gmd::System make_system() {
    gmd::System system;
    system.resize(4, 4);
    gmd::Box box;
    box.set_lengths({20.0, 21.0, 22.0});
    system.set_box(box);
    auto masses = system.mutable_masses();
    auto charges = system.mutable_charges();
    auto types = system.mutable_atom_types();
    auto molecules = system.mutable_molecule_ids();
    auto coords = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();
    for (std::size_t i = 0; i < system.atom_count(); ++i) {
        masses[i] = 10.0 + static_cast<double>(i);
        charges[i] = 0.1 * static_cast<double>(i);
        types[i] = static_cast<int>(i % 2);
        molecules[i] = static_cast<int>(i / 2);
        coords[i] = {1.0 + static_cast<double>(i),
                     2.0 + 0.5 * static_cast<double>(i),
                     3.0 + 0.25 * static_cast<double>(i)};
        velocities[i] = {0.01 * static_cast<double>(i + 1),
                         -0.02 * static_cast<double>(i + 1),
                         0.03 * static_cast<double>(i + 1)};
    }
    return system;
}

void run_steps(gmd::System& system, std::uint64_t start_step, std::uint64_t steps) {
    auto integrator = std::make_shared<gmd::VelocityVerletIntegrator>(0.5);
    auto force = std::make_shared<ZeroForceProvider>();
    gmd::RuntimeContext runtime;
    gmd::Simulation simulation(&system);
    simulation.set_force_provider(force);
    simulation.set_integrator(integrator);
    simulation.set_time_step(0.5);
    simulation.initialize(runtime);
    simulation.set_current_step(start_step);
    simulation.run(runtime, steps);
}

void compare_systems(const gmd::System& a,
                     const gmd::System& b,
                     double tolerance,
                     int& failures) {
    check(a.atom_count() == b.atom_count(), "atom counts differ", failures);
    const auto ca = a.coordinates();
    const auto cb = b.coordinates();
    const auto va = a.velocities();
    const auto vb = b.velocities();
    for (std::size_t i = 0; i < a.atom_count(); ++i) {
        for (std::size_t dim = 0; dim < 3; ++dim) {
            check(std::abs(ca[i][dim] - cb[i][dim]) < tolerance,
                  "positions differ after restart", failures);
            check(std::abs(va[i][dim] - vb[i][dim]) < tolerance,
                  "velocities differ after restart", failures);
        }
    }
    check(std::abs(a.potential_energy() - b.potential_energy()) < tolerance,
          "potential energies differ after restart", failures);
}

void test_serial_restart_equivalence(int& failures) {
    const std::filesystem::path checkpoint_path =
        std::filesystem::current_path() / "checkpoint_restart_test.gmdchk";

    gmd::System continuous = make_system();
    run_steps(continuous, 0, 100);

    gmd::System split = make_system();
    run_steps(split, 0, 50);
    gmd::Topology topology;
    topology.bonds.push_back(gmd::BondTerm{0, 1, 0});
    gmd::CheckpointMetadata metadata;
    metadata.step = 50;
    metadata.time_fs = 25.0;
    metadata.xyz_file = "input.xyz";
    metadata.run_file = "run.in";
    metadata.force_field_file = "ff.ff";
    metadata.topology_file = "top.top";
    metadata.velocity_seed = 1234;
    metadata.config_summary = "test";
    gmd::CheckpointData checkpoint{metadata, &split, &topology};
    gmd::write_checkpoint(checkpoint_path, checkpoint);

    gmd::System restarted;
    gmd::Topology restored_topology;
    const auto restored_metadata =
        gmd::read_checkpoint(checkpoint_path, restarted, &restored_topology);
    check(restored_metadata.step == 50, "checkpoint step was not restored", failures);
    check(restored_metadata.time_fs == 25.0, "checkpoint time was not restored", failures);
    check(restored_topology.bonds.size() == 1, "checkpoint topology was not restored", failures);
    run_steps(restarted, restored_metadata.step, 50);

    compare_systems(continuous, restarted, 1.0e-12, failures);
}

void test_state_roundtrip(int& failures) {
    gmd::System system = make_system();
    gmd::NoseHooverThermostat nh(50.0);
    nh.initialize(system);
    const std::string nh_state = nh.checkpoint_state();
    gmd::NoseHooverThermostat nh_restored;
    nh_restored.load_checkpoint_state(nh_state);
    check(nh_restored.checkpoint_state() == nh_state,
          "Nose-Hoover checkpoint state did not round-trip", failures);

    gmd::MCBarostat mc(5, 0.02, 10, 123);
    const std::string mc_state = mc.checkpoint_state();
    gmd::MCBarostat mc_restored;
    mc_restored.load_checkpoint_state(mc_state);
    check(mc_restored.checkpoint_state() == mc_state,
          "MC barostat RNG/counter state did not round-trip", failures);
}

void test_bad_checkpoint_errors(int& failures) {
    const auto dir = std::filesystem::current_path();
    const auto corrupt_path = dir / "checkpoint_corrupt.gmdchk";
    {
        std::ofstream out(corrupt_path);
        out << "not a checkpoint\n";
    }
    gmd::System system;
    bool threw = false;
    try {
        (void)gmd::read_checkpoint(corrupt_path, system, nullptr);
    } catch (const std::exception& error) {
        threw = std::string(error.what()).find("Invalid checkpoint header") != std::string::npos;
    }
    check(threw, "corrupt checkpoint did not produce a clear error", failures);

    const auto version_path = dir / "checkpoint_bad_version.gmdchk";
    {
        std::ofstream out(version_path);
        out << "GMD_CHECKPOINT 999\n";
    }
    threw = false;
    try {
        (void)gmd::read_checkpoint(version_path, system, nullptr);
    } catch (const std::exception& error) {
        threw = std::string(error.what()).find("Unsupported checkpoint version") != std::string::npos;
    }
    check(threw, "bad checkpoint version did not produce a clear error", failures);
}

}  // namespace

int main() {
    int failures = 0;
    test_serial_restart_equivalence(failures);
    test_state_roundtrip(failures);
    test_bad_checkpoint_errors(failures);
    if (failures != 0) {
        std::cerr << failures << " checkpoint test(s) failed\n";
        return 1;
    }
    std::cout << "checkpoint/restart tests passed\n";
    return 0;
}
