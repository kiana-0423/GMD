#include "gmd_next/model/cell.hpp"
#include "gmd_next/model/lj_model.hpp"
#include "gmd_next/model/run_spec.hpp"

#include <algorithm>

namespace gmd_next::model {

std::string_view name(CutoffMode mode) {
    switch (mode) {
    case CutoffMode::none: return "none";
    case CutoffMode::potential_shift: return "potential_shift";
    case CutoffMode::force_shift: return "force_shift";
    }
    return "unknown";
}

std::string_view name(Ensemble ensemble) {
    switch (ensemble) {
    case Ensemble::nve: return "nve";
    case Ensemble::nvt: return "nvt";
    case Ensemble::npt: return "npt";
    }
    return "unknown";
}

double list_radius(const LjModel& model) {
    return model.cutoff + model.skin;
}

bool is_fully_periodic(const CellSpec& cell) {
    return std::all_of(cell.periodic.begin(), cell.periodic.end(), [](bool flag) { return flag; });
}

std::optional<double> shortest_periodic_length(const CellSpec& cell) {
    std::optional<double> shortest;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (!cell.periodic[axis]) continue;
        if (!shortest || cell.lengths[axis] < *shortest) shortest = cell.lengths[axis];
    }
    return shortest;
}

double volume(const CellSpec& cell) {
    return cell.lengths[0] * cell.lengths[1] * cell.lengths[2];
}

}  // namespace gmd_next::model
