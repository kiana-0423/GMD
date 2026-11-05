// Checkpoint format version 3 and the file utilities it is written with.
//
//  * SHA-256 against the FIPS 180-4 / NIST example vectors (provenance hashes
//    must be the standard digest, so `sha256sum model.pt` agrees with GMD).
//  * atomic_write_file: replaces in place, leaves no temporary behind, and on
//    failure leaves the previous file untouched.
//  * Version 3 round trip: atomic numbers and the provenance block survive
//    bit for bit, alongside everything version 2 already stored.
//  * Version 1 and 2 compatibility: still readable, reported as carrying no
//    atomic numbers (every Z is 0), so the application can refuse them for ML.
//  * step-suffixed checkpoint names.

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "gmd/io/checkpoint.hpp"
#include "gmd/io/file_utils.hpp"
#include "gmd/system/system.hpp"
#include "scoped_temp_dir.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[checkpoint v3] " << message << "\n";
        ++failures;
    }
}

std::string read_file(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

void test_sha256_vectors() {
    check(gmd::sha256_hex("") ==
              "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
          "sha256('') is wrong");
    check(gmd::sha256_hex("abc") ==
              "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
          "sha256('abc') is wrong");
    check(gmd::sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
              "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
          "sha256(448-bit message) is wrong");
    const std::string million(1000000, 'a');
    check(gmd::sha256_hex(million) ==
              "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
          "sha256(10^6 x 'a') is wrong");
}

void test_atomic_write(const std::filesystem::path& dir) {
    const auto path = dir / "atomic.txt";
    gmd::atomic_write_file(path, "first\n");
    check(read_file(path) == "first\n", "atomic write did not create the file");
    gmd::atomic_write_file(path, "second version\n");
    check(read_file(path) == "second version\n", "atomic write did not replace the file");
    check(gmd::sha256_file_hex(path) == gmd::sha256_hex("second version\n"),
          "sha256_file_hex disagrees with sha256_hex");

    std::size_t entries = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        (void)entry;
        ++entries;
    }
    check(entries == 1, "atomic write left a temporary file behind (" +
                            std::to_string(entries) + " entries)");

    // A failing write (directory does not exist) throws and creates nothing.
    bool threw = false;
    try {
        gmd::atomic_write_file(dir / "missing_dir" / "x.txt", "never\n");
    } catch (const std::runtime_error&) {
        threw = true;
    }
    check(threw, "atomic write into a missing directory did not throw");
    check(!std::filesystem::exists(dir / "missing_dir"), "failed atomic write created a path");

    // A failing rename (target is a non-empty directory) leaves the target
    // untouched and removes the temporary.
    const auto blocker = dir / "blocker";
    std::filesystem::create_directory(blocker);
    gmd::atomic_write_file(blocker / "inside.txt", "keep\n");
    threw = false;
    try {
        gmd::atomic_write_file(blocker, "cannot replace a directory\n");
    } catch (const std::runtime_error&) {
        threw = true;
    }
    check(threw, "atomic write over a non-empty directory did not throw");
    check(read_file(blocker / "inside.txt") == "keep\n", "failed rename damaged the target");
    entries = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        (void)entry;
        ++entries;
    }
    check(entries == 2, "failed atomic write left a temporary behind");
}

gmd::System make_system() {
    gmd::System system;
    system.resize(3, 3);
    gmd::Box box;
    box.set_lengths({10.25, 11.5, 12.125});
    system.set_box(box);
    const int z[3] = {8, 1, 1};
    for (std::size_t i = 0; i < 3; ++i) {
        system.mutable_masses()[i] = (i == 0) ? 15.999 : 1.008;
        system.mutable_charges()[i] = (i == 0) ? -0.8476 : 0.4238;
        system.mutable_atom_types()[i] = (i == 0) ? 0 : 1;
        system.mutable_molecule_ids()[i] = 7;
        system.mutable_atomic_numbers()[i] = z[i];
        // Values with no short decimal representation: the round trip must be
        // exact, not merely close.
        system.mutable_coordinates()[i] = {1.0 / 3.0 + static_cast<double>(i),
                                           2.0 / 7.0, 0.1 + 0.2};
        system.mutable_velocities()[i] = {1e-3 / 3.0, -2.0 / 11.0 * 1e-2,
                                          static_cast<double>(i) * 1e-17};
    }
    return system;
}

void test_v3_round_trip(const std::filesystem::path& dir) {
    gmd::System system = make_system();
    gmd::CheckpointMetadata md;
    md.step = 1234;
    md.time_fs = 617.0;
    md.thermostat_type = "nose_hoover";
    md.thermostat_state = "nose_hoover 0.125 3.5";
    md.force_field_kind = "ml";
    md.force_field_sha256 = std::string(64, 'a');
    md.model_path = "/models/with space/model.pt";
    md.model_sha256 = std::string(64, 'a');
    md.run_id = "al round 3 segment 7";

    const auto path = dir / "v3.gmdchk";
    gmd::write_checkpoint(path, {.metadata = md, .system = &system, .topology = nullptr});
    const std::string text = read_file(path);
    check(text.rfind("GMD_CHECKPOINT 3\n", 0) == 0, "checkpoint does not start with version 3");
    check(text.find("atoms tag type molecule mass charge x y z vx vy vz atomic_number\n") !=
              std::string::npos,
          "version 3 atoms header is missing");

    gmd::System restored;
    const auto read = gmd::read_checkpoint(path, restored);
    check(read.version == 3, "read version is not 3");
    check(read.has_atomic_numbers, "version 3 checkpoint reports no atomic numbers");
    check(read.step == 1234 && read.time_fs == 617.0, "step/time not restored");
    check(read.force_field_kind == "ml", "force_field_kind not restored");
    check(read.force_field_sha256 == md.force_field_sha256, "force_field_sha256 not restored");
    check(read.model_path == md.model_path, "model_path (with a space) not restored");
    check(read.model_sha256 == md.model_sha256, "model_sha256 not restored");
    check(read.run_id == md.run_id, "run_id (with spaces) not restored");
    check(read.thermostat_state == md.thermostat_state, "thermostat state not restored");
    check(restored.atom_count() == 3, "atom count not restored");
    for (std::size_t i = 0; i < 3; ++i) {
        check(restored.atomic_numbers()[i] == system.atomic_numbers()[i],
              "atomic number of atom " + std::to_string(i) + " not restored");
        for (int d = 0; d < 3; ++d) {
            check(restored.coordinates()[i][d] == system.coordinates()[i][d],
                  "coordinate not restored bitwise");
            check(restored.velocities()[i][d] == system.velocities()[i][d],
                  "velocity not restored bitwise");
        }
        check(restored.masses()[i] == system.masses()[i], "mass not restored");
        check(restored.charges()[i] == system.charges()[i], "charge not restored");
        check(restored.atom_types()[i] == system.atom_types()[i], "type not restored");
        check(restored.atom_tag(i) == system.atom_tag(i), "tag not restored");
    }
    check(restored.box().lengths == system.box().lengths, "box not restored");

    // Empty provenance (a classical run without a run id) round-trips as empty.
    gmd::CheckpointMetadata classical;
    classical.force_field_kind = "default_lj";
    gmd::write_checkpoint(path, {.metadata = classical, .system = &system, .topology = nullptr});
    gmd::System restored_classical;
    const auto read_classical = gmd::read_checkpoint(path, restored_classical);
    check(read_classical.force_field_kind == "default_lj", "classical kind not restored");
    check(read_classical.model_path.empty() && read_classical.model_sha256.empty() &&
              read_classical.run_id.empty(),
          "empty provenance fields did not round-trip as empty");
}

std::string legacy_checkpoint(int version) {
    std::ostringstream out;
    out << "GMD_CHECKPOINT " << version << "\n"
        << "step 50\n"
        << "time_fs 25\n"
        << "boundary periodic periodic periodic\n"
        << "box 20 20 20\n"
        << "xyz_file \n"
        << "run_file \n"
        << "force_field_file \n"
        << "topology_file \n"
        << "velocity_seed 5489\n"
        << "force_field_summary default_lj\n"
        << "config_summary legacy\n"
        << "thermostat_type \n"
        << "thermostat_state stateless\n"
        << "barostat_type \n"
        << "barostat_state stateless\n";
    if (version >= 2) {
        out << "constraint_virial not_applicable none 0 0 0 0 0 0 0 0 0\n"
            << "provider_virial 0 0 0 0 0 0 0 0 0 0\n"
            << "step_pressure 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
    }
    out << "atom_count 2\n"
        << "atoms tag type molecule mass charge x y z vx vy vz\n"
        << "0 0 0 39.948 0 1 2 3 0.1 0.2 0.3\n"
        << "1 0 0 39.948 0 4 5 6 -0.1 -0.2 -0.3\n"
        << "topology_counts 0 0 0 0 0\n"
        << "end\n";
    return out.str();
}

void test_legacy_versions(const std::filesystem::path& dir) {
    for (int version : {1, 2}) {
        const auto path = dir / ("v" + std::to_string(version) + ".gmdchk");
        std::ofstream(path) << legacy_checkpoint(version);
        gmd::System system;
        bool threw = false;
        gmd::CheckpointMetadata md;
        try {
            md = gmd::read_checkpoint(path, system);
        } catch (const std::exception& error) {
            threw = true;
            std::cerr << "[checkpoint v3] v" << version << " read failed: " << error.what() << "\n";
        }
        check(!threw, "version " + std::to_string(version) + " checkpoint is no longer readable");
        check(md.version == version, "legacy version not reported");
        check(!md.has_atomic_numbers,
              "version " + std::to_string(version) + " claims to carry atomic numbers");
        check(system.atom_count() == 2 && system.atomic_numbers()[0] == 0 &&
                  system.atomic_numbers()[1] == 0,
              "legacy checkpoint atoms should read with Z = 0");
        check(md.force_field_kind.empty() && md.model_sha256.empty(),
              "legacy checkpoint should carry no provenance");
        check(system.coordinates()[1][2] == 6.0, "legacy coordinates not restored");
    }

    // A version 3 file with a version 2 atoms header is corrupt, not legacy.
    std::string bad = legacy_checkpoint(2);
    bad.replace(bad.find("GMD_CHECKPOINT 2"), 16, "GMD_CHECKPOINT 3");
    bad.replace(bad.find("atom_count"), 0,
                "force_field_kind default_lj\nforce_field_sha256 \nmodel_path \n"
                "model_sha256 \nrun_id \n");
    const auto bad_path = dir / "bad_v3.gmdchk";
    std::ofstream(bad_path) << bad;
    gmd::System system;
    bool threw = false;
    try {
        gmd::read_checkpoint(bad_path, system);
    } catch (const std::runtime_error& error) {
        threw = std::string(error.what()).find("atoms header") != std::string::npos;
    }
    check(threw, "version 3 checkpoint without the atomic_number column was accepted");

    // Version 4 does not exist yet.
    std::string future = legacy_checkpoint(2);
    future.replace(future.find("GMD_CHECKPOINT 2"), 16, "GMD_CHECKPOINT 4");
    std::ofstream(dir / "v4.gmdchk") << future;
    threw = false;
    try {
        gmd::read_checkpoint(dir / "v4.gmdchk", system);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    check(threw, "unsupported checkpoint version 4 was accepted");
}

void test_step_suffix() {
    check(gmd::step_suffixed_checkpoint_path("run/state.gmdchk", 1000) ==
              std::filesystem::path("run/state_step000000001000.gmdchk"),
          "step suffix with extension is wrong");
    check(gmd::step_suffixed_checkpoint_path("state", 7) ==
              std::filesystem::path("state_step000000000007"),
          "step suffix without extension is wrong");
}

}  // namespace

int main() {
    gmd_test::ScopedTempDir temp("gmd_checkpoint_v3");
    test_sha256_vectors();
    std::filesystem::create_directories(temp.path() / "atomic");
    test_atomic_write(temp.path() / "atomic");
    test_v3_round_trip(temp.path());
    test_legacy_versions(temp.path());
    test_step_suffix();
    // ScopedTempDir removes flat files only; the atomic-write fixture made a
    // subdirectory tree of known shape inside the private temp directory.
    std::error_code ignored;
    std::filesystem::remove_all(temp.path() / "atomic", ignored);

    if (failures != 0) {
        std::cerr << "[checkpoint v3] " << failures << " failure(s)\n";
        return 1;
    }
    std::cout << "[checkpoint v3] all checks passed\n";
    return 0;
}
