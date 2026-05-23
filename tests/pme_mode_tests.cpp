#include <array>
#include <cmath>
#include <iostream>
#include <string>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/system/system.hpp"

namespace {

void check(bool condition, const std::string& message, int& failures) {
    if (!condition) {
        std::cerr << "[pme-mode] " << message << "\n";
        ++failures;
    }
}

gmd::System make_charged_system() {
    gmd::System system;
    system.resize(4, 4);
    gmd::Box box;
    box.set_lengths({18.0, 21.0, 25.0});
    system.set_box(box);
    auto masses = system.mutable_masses();
    auto charges = system.mutable_charges();
    auto coords = system.mutable_coordinates();
    for (std::size_t i = 0; i < system.atom_count(); ++i) {
        masses[i] = 1.0;
    }
    charges[0] = 1.0;
    charges[1] = -1.0;
    charges[2] = 0.5;
    charges[3] = -0.5;
    coords[0] = {2.0, 3.0, 4.0};
    coords[1] = {7.0, 5.5, 9.0};
    coords[2] = {11.0, 14.0, 13.5};
    coords[3] = {15.0, 18.0, 20.0};
    return system;
}

gmd::ForceResult compute_pme(gmd::System& system, gmd::PmeExecutionMode mode) {
    gmd::RuntimeContext runtime;
    gmd::PMEForceProvider pme(0.28, 8.0, 4, {16, 16, 16}, mode, false);
    pme.initialize(runtime);
    gmd::ForceResult result;
    const auto coords = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .step = 0,
        .time = 0.0,
        .coordinates = std::span<const gmd::Coordinate3D>(coords.data(), coords.size()),
        .neighbor_list = nullptr,
    };
    pme.compute(request, result, runtime);
    return result;
}

}  // namespace

int main() {
    int failures = 0;
    auto replicated_system = make_charged_system();
    auto distributed_system = make_charged_system();
    const auto replicated = compute_pme(replicated_system, gmd::PmeExecutionMode::Replicated);
    const auto distributed = compute_pme(distributed_system, gmd::PmeExecutionMode::Distributed);

    check(replicated.success && distributed.success, "PME compute failed", failures);
    check(std::abs(replicated.potential_energy - distributed.potential_energy) < 1.0e-12,
          "distributed PME energy differs from replicated PME", failures);
    for (std::size_t i = 0; i < replicated.forces.size(); ++i) {
        for (int dim = 0; dim < 3; ++dim) {
            check(std::abs(replicated.forces[i][dim] - distributed.forces[i][dim]) < 1.0e-12,
                  "distributed PME force differs from replicated PME", failures);
        }
    }

    if (failures != 0) {
        std::cerr << failures << " PME mode test(s) failed\n";
        return 1;
    }
    std::cout << "PME mode tests passed\n";
    return 0;
}
