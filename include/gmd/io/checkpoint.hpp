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
inline constexpr int kCheckpointVersion = 2;
inline constexpr int kMinReadableCheckpointVersion = 1;

struct CheckpointMetadata {
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
};

struct CheckpointData {
    CheckpointMetadata metadata;
    System* system = nullptr;
    const Topology* topology = nullptr;
};

void write_checkpoint(const std::filesystem::path& path,
                      const CheckpointData& checkpoint);

CheckpointMetadata read_checkpoint(const std::filesystem::path& path,
                                   System& system,
                                   Topology* topology = nullptr);

}  // namespace gmd
