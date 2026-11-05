#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>

#include "gmd/system/topology.hpp"

namespace gmd {

class System;

// Version 2 adds the constraint virial. Version 1 files are still read: the
// field is simply absent and the restarted run reports its first frame's
// pressure without the constraint contribution, exactly as before.
//
// Version 3 adds the atomic number Z of every atom (an extra `atomic_number`
// column) and a provenance block naming the force field that produced the
// state: its kind, a SHA-256 fingerprint, the ML model path and the SHA-256 of
// the model file, and an optional run id. Version 1 and 2 files are still read;
// they carry no Z, so `has_atomic_numbers` is false and every Z is 0. That is
// harmless for a classical force field, which never reads Z, and is refused by
// the application for an ML force field, which would otherwise be fed species 0.
inline constexpr int kCheckpointVersion = 3;
inline constexpr int kMinReadableCheckpointVersion = 1;
inline constexpr int kFirstCheckpointVersionWithAtomicNumbers = 3;

struct CheckpointMetadata {
    // Version the file was READ from. write_checkpoint always writes
    // kCheckpointVersion and ignores this field.
    int version = kCheckpointVersion;
    // True when the file carried a per-atom atomic number column (version 3+).
    bool has_atomic_numbers = true;

    std::uint64_t step = 0;
    double time_fs = 0.0;
    std::string boundary = "periodic periodic periodic";
    std::string xyz_file;
    std::string run_file;
    std::string force_field_file;
    std::string topology_file;
    std::uint32_t velocity_seed = 5489u;
    std::string force_field_summary;
    std::string config_summary;
    std::string thermostat_type;
    std::string thermostat_state = "stateless";
    std::string barostat_type;
    std::string barostat_state = "stateless";

    // --- Virial and pressure state (version 2 and later) ---------------------
    //
    // Three explicitly named quantities, because they are three different things
    // and a restart has to reproduce each of them.
    //
    // 1. CONSTRAINT CONTRIBUTION, and its validity and time level. `state`
    //    is one of "not_applicable" (the run has no constraints),
    //    "unavailable" (constraints act but no multiplier belongs to the stored
    //    provider virial) or "valid". `time_level` names which multiplier it is;
    //    the only value this build writes is "endpoint_rattle_t_plus_dt", the
    //    RATTLE multiplier of the step that produced the checkpointed state,
    //    which belongs to exactly the coordinates stored here. Nothing in the
    //    state itself determines it: it is a property of the step that reached
    //    the state, which is why it has to be persisted. A restart re-evaluates
    //    forces at these coordinates and attaches this to that provider virial.
    std::array<double, 9> constraint_virial{};
    std::string constraint_virial_state = "not_applicable";
    std::string constraint_virial_time_level = "none";

    // 2. CURRENT-GEOMETRY PROVIDER VIRIAL, as it stood when the checkpoint was
    //    written. This is RECORDED BUT NOT INSTALLED on restart: the restarted
    //    run recomputes it from the checkpointed coordinates, which is where it
    //    came from in the first place, and reinstalling a value summed under a
    //    different rank decomposition would only introduce a discrepancy. It is
    //    persisted so the two can be compared when diagnosing a restart.
    std::array<double, 9> provider_virial{};
    bool provider_virial_valid = false;

    // 3. COMPLETED-STEP THERMODYNAMICS: the pressure of the step that produced
    //    this state, captured before any barostat rescaled the cell, together
    //    with the virial, 2K, volume and potential energy it was taken with --
    //    one mutually consistent set. This is what the frame reported, and it is
    //    installed on restart so the restarted run reports the identical numbers
    //    for the frame rather than reconstructing them.
    bool step_pressure_valid = false;
    double step_pressure = 0.0;                  // [eV/A^3]
    std::array<double, 9> step_pressure_virial{};
    double step_pressure_twice_ke = 0.0;
    double step_pressure_volume = 0.0;
    double step_pressure_potential_energy = 0.0;

    // --- Provenance (version 3 and later; empty when read from older files) --
    //
    // force_field_kind is one of "inline_lj", "lj_file", "molecular",
    // "default_lj" or "ml", optionally suffixed "+ewald" / "+pme".
    // force_field_sha256 fingerprints the force field: for "ml" it is the
    // SHA-256 of the model file, for a file-based field the SHA-256 of that
    // file, and otherwise the SHA-256 of the canonical parameter summary. A
    // restart compares it against the current run. Values never contain
    // whitespace except run_id and model_path, which are stored to end of line.
    std::string force_field_kind;
    std::string force_field_sha256;
    std::string model_path;
    std::string model_sha256;
    std::string run_id;
};

struct CheckpointData {
    CheckpointMetadata metadata;
    System* system = nullptr;
    const Topology* topology = nullptr;
};

// Serialises the checkpoint to text (always kCheckpointVersion).
std::string format_checkpoint(const CheckpointData& checkpoint);

// Writes the checkpoint ATOMICALLY: the text goes to a temporary file in the
// same directory, is fsync'ed and then renamed over `path` (see
// atomic_write_file). A crash mid-write leaves the previous checkpoint intact.
void write_checkpoint(const std::filesystem::path& path,
                      const CheckpointData& checkpoint);

// `state.gmdchk` -> `state_step000000001000.gmdchk`: the step, zero-padded to
// 12 digits so names sort in step order, inserted before the extension.
std::filesystem::path step_suffixed_checkpoint_path(const std::filesystem::path& path,
                                                    std::uint64_t step);

CheckpointMetadata read_checkpoint(const std::filesystem::path& path,
                                   System& system,
                                   Topology* topology = nullptr);

}  // namespace gmd
