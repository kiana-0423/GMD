#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>

namespace gmd {

class System;

inline constexpr const char* kFrameFormat = "gmd-frame-v1";

// Provenance written into every frame's comment line.
struct FrameProvenance {
    std::string run_id;              // empty -> written as run_id=none
    std::string force_field_kind;    // e.g. "ml", "inline_lj"
    std::string model_sha256;        // empty for classical -> model_sha256=none
    // Set when this run restarted from a checkpoint produced by a different
    // force field and `model_change allow` accepted it: the checkpoint's model
    // (or force-field) SHA-256. Omitted from the frame otherwise.
    std::string model_changed_from;
};

// ASE-readable extended-XYZ frame output (format "gmd-frame-v1").
//
// One frame per write_frame() call:
//
//   <N>
//   Lattice="Lx 0 0 0 Ly 0 0 0 Lz" Properties=id:I:1:species:S:1:pos:R:3:forces:R:3[:vel:R:3]
//     pbc="T T T" step=<int> time_fs=<f> energy=<PE> potential_energy=<PE>
//     kinetic_energy=<KE> temperature=<K> volume=<A^3> format=gmd-frame-v1
//     units="..." run_id=<id|none> force_field=<kind> model_sha256=<hex|none>
//     [model_changed_from=<hex>]
//   <id> <symbol> x y z fx fy fz [vx vy vz]      (one line per atom)
//
// (The comment is a single line; it is wrapped above for readability.)
//
//  * id is the GMD atom tag: 0-based, stable for the whole run and across
//    restarts, equal to the atom's 0-based line index in the input xyz. Atom
//    lines are always written in tag order.
//  * species is the element symbol from the atomic number; Z = 0 (unknown) is
//    written as "X", ASE's dummy element.
//  * Positions [A] are wrapped into the box, forces [eV/A], velocities [A/fs]
//    (only with frame_velocities on). Every real number is written with 17
//    significant digits (%.17g), which round-trips an IEEE double exactly.
//  * energy and potential_energy are the same number [eV]: `energy` is the key
//    ASE attaches to its SinglePointCalculator, potential_energy is explicit.
//    All per-frame values describe the written geometry: positions, forces,
//    potential energy, and kinetic energy / temperature from the written
//    velocities.
//
// FLUSHING. Each frame is formatted completely in memory, written with one
// write and flushed before write_frame() returns, so whenever GMD is between
// frames the file ends on a complete frame. A reader polling concurrently with
// a write may briefly see an incomplete trailing frame and should use the
// atom-count line to ignore it.
class FrameWriter {
public:
    FrameWriter() = default;
    ~FrameWriter();

    // Opens (truncates) `path`. Throws std::runtime_error on failure.
    void open(const std::filesystem::path& path, bool write_velocities,
              FrameProvenance provenance);
    void close();
    bool is_open() const noexcept { return out_.is_open(); }

    // `twice_ke` is sum m v^2 [eV] for the written velocities; `dof` the
    // degrees of freedom used for the temperature.
    void write_frame(const System& system, std::uint64_t step, double time_fs,
                     double twice_ke, std::size_t dof);

    std::uint64_t frame_count() const noexcept { return frame_count_; }

    static std::string format_frame(const System& system, std::uint64_t step, double time_fs,
                                    double twice_ke, std::size_t dof, bool write_velocities,
                                    const FrameProvenance& provenance);

private:
    std::ofstream out_;
    bool write_velocities_ = false;
    FrameProvenance provenance_;
    std::uint64_t frame_count_ = 0;
};

// Element symbol for an atomic number; "X" for 0 or anything out of range.
const char* element_symbol(int atomic_number) noexcept;

}  // namespace gmd
