#pragma once

#include <filesystem>
#include <string>
#include <string_view>

namespace gmd {

// Lower-case hexadecimal SHA-256 (FIPS 180-4) of a byte string.
std::string sha256_hex(std::string_view data);

// Lower-case hexadecimal SHA-256 of a file's bytes. Throws std::runtime_error
// if the file cannot be read.
std::string sha256_file_hex(const std::filesystem::path& path);

// Replaces `path` with `content` atomically.
//
// The bytes go to a uniquely named temporary file in the SAME directory, are
// flushed and fsync'ed, and the temporary is then rename()d over `path`. On a
// POSIX file system rename() within one directory is atomic, so a reader (or a
// restart after a crash) sees either the complete previous file or the complete
// new one, never a truncated mixture. The directory is fsync'ed afterwards on a
// best-effort basis so the rename itself survives a power loss.
//
// On any failure the temporary is removed, `path` is left untouched, and
// std::runtime_error is thrown.
void atomic_write_file(const std::filesystem::path& path, std::string_view content);

}  // namespace gmd
