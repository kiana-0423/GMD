#pragma once

#include "gmd_next/core/observables.hpp"
#include "gmd_next/core/units.hpp"
#include "gmd_next/model/cell.hpp"
#include "gmd_next/model/lj_model.hpp"

#include <cstdint>
#include <string_view>

namespace gmd_next::model {

// Only nve is implemented. The other two are listed so that a run asking for
// them is rejected with a named reason instead of being read as nve.
enum class Ensemble : std::uint8_t { nve, nvt, npt };

std::string_view name(Ensemble ensemble);

// A complete run request, before validation. Every field is consumed: nothing
// here may be ignored because the selected configuration cannot honour it.
struct RunSpec {
    core::UnitSystem units = core::UnitSystem::metal;
    Ensemble ensemble = Ensemble::nve;
    double timestep = 0.0;  // expressed in timestep_unit, finite and positive
    core::TimeUnit timestep_unit = core::TimeUnit::femtoseconds;
    std::int64_t steps = 0;  // nonnegative; zero means evaluate and stop
    CellSpec cell{};
    LjModel model{};
    // Quantities to report. Forces are computed every step regardless; energy
    // and virial are computed only when dynamics or output asks for them.
    core::ObservableSet output_observables{};
    bool remove_initial_center_of_mass_motion = false;
};

}  // namespace gmd_next::model
