#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

namespace gmd {

class System;

struct RunConfig {
	std::uint64_t num_steps = 0;
	double time_step = 0.0;
	double temperature = 0.0;
	std::string velocity_init_mode = "random";
	std::uint32_t velocity_seed = 5489u;
	bool remove_center_of_mass_velocity = true;
};

class ConfigLoader {
public:
	ConfigLoader() = default;

	void load_xyz(const std::filesystem::path& xyz_path, System& system) const;
	RunConfig load_run(const std::filesystem::path& run_path) const;
};

}  // namespace gmd
