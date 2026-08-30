#pragma once

// A private, uniquely-named temporary directory that cleans itself up.
//
// WHY THIS EXISTS. tests/pressure_unit_tests.cpp used to write
// <system temp>/gmd_pressure_unit_probe.log with that exact name, every time.
// One process at a time per build directory is fine, but the same test binary
// exists in every configured build tree -- Release and Debug, serial and
// MPI-enabled -- and running two of them at once has each truncating and then
// deleting the other's file underneath it. The failure would be intermittent
// and would look like a defect in the code under test rather than in the
// harness.
//
// mkdtemp() is the fix rather than a longer fixed name or a PID suffix: it
// creates the directory atomically with mode 0700 and fails if the name is
// taken, so there is no window between choosing a name and owning it, and no
// predictable path another process could have created first.
//
// CLEANUP is deliberately not recursive. The destructor removes the regular
// files it finds in its own directory and then the directory itself, and gives
// up quietly if anything else is in there. Tests create flat files here and
// nothing else, so this is sufficient, and it cannot become a recursive delete
// of a path that turned out not to be what was expected.

#include <cerrno>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace gmd_test {

class ScopedTempDir {
public:
    // `prefix` names the directory for a human reading /tmp; uniqueness comes
    // from mkdtemp, not from the prefix.
    explicit ScopedTempDir(const std::string& prefix) {
        const std::filesystem::path base = std::filesystem::temp_directory_path();
        std::string pattern = (base / (prefix + ".XXXXXX")).string();
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        if (::mkdtemp(buffer.data()) == nullptr) {
            throw std::runtime_error("could not create a temporary directory from " +
                                     pattern + ": " + std::strerror(errno));
        }
        path_ = std::filesystem::path(buffer.data());
    }

    ScopedTempDir(const ScopedTempDir&) = delete;
    ScopedTempDir& operator=(const ScopedTempDir&) = delete;

    ~ScopedTempDir() {
        std::error_code ec;
        for (const auto& entry : std::filesystem::directory_iterator(path_, ec)) {
            if (ec) break;
            std::error_code remove_ec;
            if (entry.is_regular_file(remove_ec) && !remove_ec) {
                std::filesystem::remove(entry.path(), remove_ec);
            }
        }
        // Non-recursive: fails harmlessly if anything unexpected remains.
        std::filesystem::remove(path_, ec);
    }

    const std::filesystem::path& path() const noexcept { return path_; }

    // Convenience for the common "give me a file in here" case.
    std::filesystem::path file(const std::string& name) const {
        return path_ / name;
    }

private:
    std::filesystem::path path_;
};

}  // namespace gmd_test
