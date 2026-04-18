#include "gmd/io/trajectory_writer.hpp"

#include <iomanip>
#include <stdexcept>

#include "gmd/integrator/thermostat.hpp"  // compute_twice_ke, temperature_from_twice_ke
#include "gmd/system/system.hpp"

namespace gmd {

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
    log_ << "# step  time[fs]  PE[eV]  KE[eV]  E_total[eV]  T[K]\n";
}

void TrajectoryWriter::write_frame(const System& system, std::uint64_t step, double time,
                                    double twice_ke, std::size_t dof) {
    const std::size_t n = system.atom_count();
    const double pe     = system.potential_energy();
    const double ke     = 0.5 * twice_ke;
    const double temp   = (dof > 0) ? temperature_from_twice_ke(twice_ke, dof) : 0.0;

    // --- XYZ frame ---
    xyz_ << n << '\n';
    xyz_ << std::fixed << std::setprecision(6)
         << "step=" << step
         << " time=" << time
         << " PE=" << pe
         << " KE=" << ke
         << " T=" << temp
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
