#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>

namespace gmd {

class System;

// Writes simulation output to disk.
//
// Two files are produced:
//   - XYZ trajectory file (*.xyz)  — one frame per write_frame() call,
//     extended XYZ format: atom-count line, comment line with step/time/PE/KE/T,
//     then one atom line per atom: "type  x  y  z"
//   - Energy log file (*.log)      — space-delimited table with columns:
//     step  time  PE  KE  E_total  temperature
//
// Usage:
//   TrajectoryWriter writer;
//   writer.open("output");          // creates output.xyz and output.log
//   writer.write_frame(system, step, time);
//   writer.close();
//
// If write_interval > 1, callers are responsible for calling write_frame only
// on the appropriate steps, or use write_frame_if() for automatic filtering.
class TrajectoryWriter {
public:
    TrajectoryWriter() = default;
    ~TrajectoryWriter();

    // Opens (or truncates) output files "<stem>.xyz" and "<stem>.log".
    // Throws std::runtime_error on failure.
    void open(const std::filesystem::path& stem);

    // Closes output files (flushed automatically on destruction if not closed).
    void close();

    bool is_open() const noexcept { return xyz_.is_open(); }

    // Writes one trajectory frame.
    // twice_ke: sum of m*v^2 (caller provides; avoids re-scanning velocities).
    // dof: degrees of freedom for temperature (typically 3N-3).
    void write_frame(const System& system, std::uint64_t step, double time,
                     double twice_ke, std::size_t dof);

    // Convenience overload that computes twice_ke internally.
    void write_frame(const System& system, std::uint64_t step, double time,
                     std::size_t dof);

    // Writes a frame only when (step % interval == 0).
    void write_frame_if(const System& system, std::uint64_t step, double time,
                        std::size_t dof, std::uint64_t interval);

    // Number of frames written so far.
    std::uint64_t frame_count() const noexcept { return frame_count_; }

private:
    std::ofstream xyz_;
    std::ofstream log_;
    std::uint64_t frame_count_ = 0;

    void write_log_header();
};

}  // namespace gmd
