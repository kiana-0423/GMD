#include "gmd/io/trajectory_writer.hpp"

#include <iomanip>
#include <limits>
#include <stdexcept>

#include "gmd/core/physical_constants.hpp"
#include "gmd/integrator/thermostat.hpp"  // compute_twice_ke, temperature_from_twice_ke
#include "gmd/system/system.hpp"

namespace gmd {

namespace {

// Single definition in gmd/core/physical_constants.hpp; not a second literal.
// The eV/A^3 -> bar direction is the one this file needs, and it is a
// multiplication rather than a division by the forward constant: the two are
// exact reciprocals, so the result is the same bits either way, but the
// multiplication says which way the conversion runs.
constexpr double kEvPerA3ToBar = kEVPerAngstromCubedToBar;

double volume_from_box(const Box& box) noexcept {
    return box.lengths[0] * box.lengths[1] * box.lengths[2];
}

// The thermodynamic state a frame reports, and whether its pressure means
// anything.
//
// PREFERENCE ORDER. The completed-step record is used when there is one. It was
// captured at the end of the step, before any barostat rescaled the cell, and
// every field in it was taken at that same instant, so the reported P, V, KE, PE
// and T describe ONE state: P = (2K + tr W) / 3V holds among exactly these
// numbers. Its pressure is also the value the barostat consumed. Using it is
// what stops a rescale from destroying the pressure the step actually measured.
//
// The current-geometry state is the fallback for a frame that no step produced
// -- the initial frame of a run -- and its pressure is used only when the
// current virial is complete.
//
// A numerical zero must never stand for "unavailable": it is a perfectly
// ordinary pressure. When no complete pressure exists the pressure is NaN and
// `valid` is false, and the log carries an explicit validity column beside it.
struct FrameThermodynamics {
    double pressure_bar = 0.0;
    bool pressure_valid = false;
    double volume = 0.0;
    double twice_ke = 0.0;
    double potential_energy = 0.0;
};

FrameThermodynamics frame_thermodynamics(const System& system,
                                         double current_twice_ke) noexcept {
    FrameThermodynamics frame;

    const auto& completed = system.step_thermodynamics();
    if (completed.valid) {
        frame.pressure_bar = completed.pressure * kEvPerA3ToBar;
        frame.pressure_valid = true;
        frame.volume = completed.volume;
        frame.twice_ke = completed.twice_kinetic_energy;
        frame.potential_energy = completed.potential_energy;
        return frame;
    }

    frame.volume = volume_from_box(system.box());
    frame.twice_ke = current_twice_ke;
    frame.potential_energy = system.potential_energy();
    if (system.last_virial_valid() && frame.volume > 0.0) {
        const auto& virial = system.last_virial();
        const double virial_trace = virial[0] + virial[4] + virial[8];
        frame.pressure_bar =
            ((current_twice_ke + virial_trace) / (3.0 * frame.volume)) * kEvPerA3ToBar;
        frame.pressure_valid = true;
        return frame;
    }

    frame.pressure_bar = std::numeric_limits<double>::quiet_NaN();
    frame.pressure_valid = false;
    return frame;
}

}  // namespace

TrajectoryWriter::~TrajectoryWriter() {
    close();
}

void TrajectoryWriter::open(const std::filesystem::path& stem) {
    close();

    const auto xyz_path = std::filesystem::path(stem).replace_extension(".xyz");
    const auto log_path = std::filesystem::path(stem).replace_extension(".log");

    xyz_.open(xyz_path, std::ios::out | std::ios::trunc);
    if (!xyz_.is_open()) {
        throw std::runtime_error("TrajectoryWriter: failed to open " + xyz_path.string());
    }

    log_.open(log_path, std::ios::out | std::ios::trunc);
    if (!log_.is_open()) {
        xyz_.close();
        throw std::runtime_error("TrajectoryWriter: failed to open " + log_path.string());
    }

    frame_count_ = 0;
    write_log_header();
}

void TrajectoryWriter::close() {
    if (xyz_.is_open()) xyz_.close();
    if (log_.is_open()) log_.close();
}

void TrajectoryWriter::write_log_header() {
    log_ << "# step  time[fs]  PE[eV]  KE[eV]  E_total[eV]  T[K]  P[bar]  V[A^3]"
         << "  shake_iter  shake_error[A]  rattle_iter  rattle_error[A/fs]"
         << "  P_valid\n";
    log_ << "# P_valid is 0 when no complete pressure exists for the frame, in "
            "which case P[bar] is nan. Zero is a pressure, not a sentinel.\n";
    log_ << "# PE, KE, E_total, T, P and V describe the completed step, taken "
            "together before any barostat rescale, so they are mutually "
            "consistent. Under a barostat the coordinates in the .xyz are the "
            "rescaled ones the next step starts from.\n";
}

void TrajectoryWriter::write_frame(const System& system, std::uint64_t step, double time,
                                    double twice_ke, std::size_t dof) {
    const std::size_t n = system.atom_count();
    // One coherent set: see frame_thermodynamics(). PE, KE, T, P and V all come
    // from the same state, which for a barostat run is the completed step rather
    // than the rescaled geometry the coordinates below describe.
    const FrameThermodynamics frame = frame_thermodynamics(system, twice_ke);
    const double pe     = frame.potential_energy;
    const double ke     = 0.5 * frame.twice_ke;
    const double temp   = (dof > 0) ? temperature_from_twice_ke(frame.twice_ke, dof) : 0.0;
    const bool pressure_valid = frame.pressure_valid;
    const double pressure_bar = frame.pressure_bar;
    const double volume = frame.volume;

    // --- XYZ frame ---
    xyz_ << n << '\n';
    xyz_ << std::fixed << std::setprecision(6)
         << "step=" << step
         << " time=" << time
         << " PE=" << pe
         << " KE=" << ke
         << " T=" << temp
         << " P=" << pressure_bar
         << " P_valid=" << (pressure_valid ? 1 : 0)
         << " V=" << volume
         << " SHAKE_iter=" << system.last_shake_stats().iterations
         << " SHAKE_error=" << system.last_shake_stats().max_error
         << " RATTLE_iter=" << system.last_rattle_stats().iterations
         << " RATTLE_error=" << system.last_rattle_stats().max_error
         << '\n';

    const auto coords     = system.coordinates();
    const auto atom_types = system.atom_types();
    xyz_ << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < n; ++i) {
        const int type = (atom_types.size() == n) ? atom_types[i] : 0;
        xyz_ << type
             << "  " << coords[i][0]
             << "  " << coords[i][1]
             << "  " << coords[i][2]
             << '\n';
    }

    // --- log line ---
    log_ << std::fixed << std::setprecision(6)
         << step
         << "  " << time
         << "  " << pe
         << "  " << ke
         << "  " << (pe + ke)
         << "  " << temp
         << "  " << pressure_bar
         << "  " << volume
         << "  " << system.last_shake_stats().iterations
         << "  " << system.last_shake_stats().max_error
         << "  " << system.last_rattle_stats().iterations
         << "  " << system.last_rattle_stats().max_error
         << "  " << (pressure_valid ? 1 : 0)
         << '\n';

    ++frame_count_;
}

void TrajectoryWriter::write_frame(const System& system, std::uint64_t step, double time,
                                    std::size_t dof) {
    const double twice_ke = compute_twice_ke(system);
    write_frame(system, step, time, twice_ke, dof);
}

void TrajectoryWriter::write_frame_if(const System& system, std::uint64_t step, double time,
                                       std::size_t dof, std::uint64_t interval) {
    if (interval == 0 || step % interval == 0) {
        write_frame(system, step, time, dof);
    }
}

}  // namespace gmd
