#include "gmd/io/trajectory_writer.hpp"

#include <iomanip>
#include <stdexcept>

#include "gmd/integrator/thermostat.hpp"  // compute_twice_ke, temperature_from_twice_ke
#include "gmd/system/system.hpp"

namespace gmd {

namespace {

constexpr double kBarToEvPerA3 = 6.2415091e-7;

double volume_from_box(const Box& box) noexcept {
    return box.lengths[0] * box.lengths[1] * box.lengths[2];
}

double pressure_bar_from_system(const System& system, double twice_ke) noexcept {
    if (!system.last_virial_valid()) {
        return 0.0;
    }

    const double volume = volume_from_box(system.box());
    if (volume <= 0.0) {
        return 0.0;
    }

    const auto& virial = system.last_virial();
    const double virial_trace = virial[0] + virial[4] + virial[8];
    const double pressure_ev_per_a3 = (twice_ke + virial_trace) / (3.0 * volume);
    return pressure_ev_per_a3 / kBarToEvPerA3;
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
         << "  shake_iter  shake_error[A]  rattle_iter  rattle_error[A/fs]\n";
}

void TrajectoryWriter::write_frame(const System& system, std::uint64_t step, double time,
                                    double twice_ke, std::size_t dof) {
    const std::size_t n = system.atom_count();
    const double pe     = system.potential_energy();
    const double ke     = 0.5 * twice_ke;
    const double temp   = (dof > 0) ? temperature_from_twice_ke(twice_ke, dof) : 0.0;
    const double pressure_bar = pressure_bar_from_system(system, twice_ke);
    const double volume = volume_from_box(system.box());

    // --- XYZ frame ---
    xyz_ << n << '\n';
    xyz_ << std::fixed << std::setprecision(6)
         << "step=" << step
         << " time=" << time
         << " PE=" << pe
         << " KE=" << ke
         << " T=" << temp
         << " P=" << pressure_bar
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
