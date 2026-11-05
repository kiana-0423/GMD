// run.in directives for external-controller integration, and run-file-relative
// path resolution.
//
// Every relative path in a run file resolves against the directory holding
// the run file; absolute paths are untouched; a run file given without a
// directory leaves relative paths relative to the working directory (which is
// then that directory). New directives default to "off" so a legacy run file
// parses to exactly the configuration it always did.

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "gmd/io/config_loader.hpp"
#include "scoped_temp_dir.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[run directives] " << message << "\n";
        ++failures;
    }
}

gmd::RunConfig load(const std::filesystem::path& path, const std::string& text) {
    std::ofstream(path) << text;
    return gmd::ConfigLoader{}.load_run(path);
}

std::string load_error(const std::filesystem::path& path, const std::string& text) {
    try {
        load(path, text);
    } catch (const std::runtime_error& error) {
        return error.what();
    }
    return {};
}

void test_defaults(const std::filesystem::path& dir) {
    const auto config = load(dir / "legacy.in", "run 10\ntime_step 1.0\n");
    check(config.run_id.empty(), "run_id should default to empty");
    check(!config.allow_model_change, "model_change should default to forbid");
    check(config.frame_file.empty(), "frame output should default to off");
    check(!config.frame_interval.has_value(), "frame_interval should default to unset");
    check(!config.frame_velocities, "frame_velocities should default to off");
    check(config.stop_file.empty() && config.checkpoint_request_file.empty() &&
              config.status_file.empty(),
          "control files should default to off");
    check(config.control_poll_interval == 1, "control_poll_interval should default to 1");
    check(!config.checkpoint_step_suffix, "checkpoint_step_suffix should default to off");
    check(config.checkpoint_file.empty(), "checkpoint_file should default to empty");
    check(config.ml_jit_optimize, "ml_jit_optimize should default to on");
}

void test_all_directives(const std::filesystem::path& dir) {
    const auto config = load(dir / "full.in",
                             "run 10\n"
                             "time_step 0.5\n"
                             "run_id al-r3.seg_7\n"
                             "model_change allow   # inline comment\n"
                             "frame_file frames/out.extxyz\n"
                             "frame_interval 5\n"
                             "frame_velocities on\n"
                             "stop_file STOP\n"
                             "checkpoint_request_file CHECKPOINT\n"
                             "control_poll_interval 3\n"
                             "status_file status.json\n"
                             "checkpoint_file state.gmdchk\n"
                             "checkpoint_step_suffix on\n"
                             "restart_from ../previous/state.gmdchk\n"
                             "force_field ml\n"
                             "model_path model.pt\n"
                             "ml_jit_optimize off\n");
    check(config.run_id == "al-r3.seg_7", "run_id not parsed");
    check(config.allow_model_change, "model_change allow not parsed");
    check(config.frame_file == dir / "frames/out.extxyz", "frame_file not resolved: " +
                                                            config.frame_file.string());
    check(config.frame_interval == 5u, "frame_interval not parsed");
    check(config.frame_velocities, "frame_velocities not parsed");
    check(config.stop_file == dir / "STOP", "stop_file not resolved");
    check(config.checkpoint_request_file == dir / "CHECKPOINT",
          "checkpoint_request_file not resolved");
    check(config.control_poll_interval == 3, "control_poll_interval not parsed");
    check(config.status_file == dir / "status.json", "status_file not resolved");
    check(config.checkpoint_file == dir / "state.gmdchk", "checkpoint_file not resolved");
    check(config.checkpoint_step_suffix, "checkpoint_step_suffix not parsed");
    check(config.restart_from == dir / "../previous/state.gmdchk", "restart_from not resolved");
    check(!config.ml_jit_optimize, "ml_jit_optimize off not parsed");
    check(config.ml_model_path == dir / "model.pt", "model_path not resolved against run dir: " +
                                                    config.ml_model_path.string());
}

void test_absolute_paths_unchanged(const std::filesystem::path& dir) {
    const std::filesystem::path absolute = dir / "elsewhere" / "model.pt";
    const auto config = load(dir / "abs.in",
                             "run 1\nforce_field ml\nmodel_path " + absolute.string() +
                                 "\nstatus_file " + (dir / "s.json").string() + "\n");
    check(config.ml_model_path == absolute, "absolute model_path was rewritten");
    check(config.status_file == dir / "s.json", "absolute status_file was rewritten");
}

void test_bare_run_file_name(const std::filesystem::path& dir) {
    // A run file named without a directory: paths stay relative, i.e. relative
    // to the working directory, which is where the run file is.
    const auto previous = std::filesystem::current_path();
    std::filesystem::current_path(dir);
    std::ofstream("bare.in") << "run 1\nframe_file f.xyz\n";
    const auto config = gmd::ConfigLoader{}.load_run("bare.in");
    std::filesystem::current_path(previous);
    check(config.frame_file == "f.xyz", "bare run file name changed a relative path");
}

void test_errors(const std::filesystem::path& dir) {
    const auto path = dir / "bad.in";
    const auto expect = [&](const std::string& text, const std::string& fragment) {
        const std::string error = load_error(path, text);
        check(!error.empty(), "accepted invalid run file:\n" + text);
        check(error.find(fragment) != std::string::npos,
              "error for:\n" + text + "does not mention '" + fragment + "': " + error);
    };
    expect("run 1\nmodel_change maybe\n", "model_change");
    expect("run 1\nframe_interval -1\nframe_file f\n", "frame_interval");
    expect("run 1\nframe_interval 5\n", "require frame_file");
    expect("run 1\nframe_velocities on\n", "require frame_file");
    expect("run 1\nframe_file f\nframe_velocities maybe\n", "frame_velocities");
    expect("run 1\ncontrol_poll_interval 0\n", "control_poll_interval");
    expect("run 1\ncheckpoint_request_file REQ\n", "requires checkpoint_file");
    expect("run 1\ncheckpoint_step_suffix on\n", "requires checkpoint_file");
    expect("run 1\nrun_id has\"quote\n", "run_id");
    expect("run 1\nstop_file\n", "missing a value");
    expect("run 1\nframe_fiel x\n", "Unsupported run input directive");
}

}  // namespace

int main() {
    gmd_test::ScopedTempDir temp("gmd_run_directives");
    // Canonical path: the temporary directory may sit behind a symlink
    // (/tmp -> /private/tmp), and resolution is purely lexical.
    const auto dir = std::filesystem::canonical(temp.path());
    test_defaults(dir);
    test_all_directives(dir);
    test_absolute_paths_unchanged(dir);
    test_bare_run_file_name(dir);
    test_errors(dir);

    if (failures != 0) {
        std::cerr << "[run directives] " << failures << " failure(s)\n";
        return 1;
    }
    std::cout << "[run directives] all checks passed\n";
    return 0;
}
