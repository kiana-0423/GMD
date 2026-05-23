#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

#include "gmd/system/topology.hpp"

namespace gmd {

class System;

inline constexpr int kCheckpointVersion = 1;

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
