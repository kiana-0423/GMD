#include "gmd/io/file_utils.hpp"

#include <array>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

#if defined(_WIN32)
#include <process.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

namespace gmd {

namespace {

// ---------------------------------------------------------------------------
// SHA-256, FIPS 180-4. Small, dependency-free, used for model provenance.
// ---------------------------------------------------------------------------
constexpr std::array<std::uint32_t, 64> kRoundConstants = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
    0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
    0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
    0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
    0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
    0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u,
};

constexpr std::uint32_t rotr(std::uint32_t value, int bits) noexcept {
    return (value >> bits) | (value << (32 - bits));
}

class Sha256 {
public:
    void update(const unsigned char* data, std::size_t size) {
        for (std::size_t i = 0; i < size; ++i) {
            block_[block_size_++] = data[i];
            if (block_size_ == 64) {
                compress();
                block_size_ = 0;
            }
        }
        total_bytes_ += size;
    }

    std::string hex_digest() {
        const std::uint64_t total_bits = total_bytes_ * 8u;
        const unsigned char pad_start = 0x80;
        const unsigned char zero = 0x00;
        update(&pad_start, 1);
        while (block_size_ != 56) {
            update(&zero, 1);
        }
        unsigned char length[8];
        for (int i = 0; i < 8; ++i) {
            length[i] = static_cast<unsigned char>(total_bits >> (56 - 8 * i));
        }
        update(length, 8);

        static constexpr char kHex[] = "0123456789abcdef";
        std::string out;
        out.reserve(64);
        for (std::uint32_t word : state_) {
            for (int shift = 28; shift >= 0; shift -= 4) {
                out.push_back(kHex[(word >> shift) & 0xFu]);
            }
        }
        return out;
    }

private:
    void compress() {
        std::array<std::uint32_t, 64> w{};
        for (int i = 0; i < 16; ++i) {
            w[i] = (static_cast<std::uint32_t>(block_[4 * i]) << 24) |
                   (static_cast<std::uint32_t>(block_[4 * i + 1]) << 16) |
                   (static_cast<std::uint32_t>(block_[4 * i + 2]) << 8) |
                   static_cast<std::uint32_t>(block_[4 * i + 3]);
        }
        for (int i = 16; i < 64; ++i) {
            const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        std::uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        std::uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (int i = 0; i < 64; ++i) {
            const std::uint32_t s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const std::uint32_t ch = (e & f) ^ (~e & g);
            const std::uint32_t t1 = h + s1 + ch + kRoundConstants[i] + w[i];
            const std::uint32_t s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const std::uint32_t t2 = s0 + maj;
            h = g; g = f; f = e; e = d + t1;
            d = c; c = b; b = a; a = t1 + t2;
        }
        state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
        state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
    }

    std::array<std::uint32_t, 8> state_ = {
        0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
        0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u,
    };
    std::array<unsigned char, 64> block_{};
    std::size_t block_size_ = 0;
    std::uint64_t total_bytes_ = 0;
};

int current_process_id() noexcept {
#if defined(_WIN32)
    return _getpid();
#else
    return static_cast<int>(::getpid());
#endif
}

}  // namespace

std::string sha256_hex(std::string_view data) {
    Sha256 hash;
    hash.update(reinterpret_cast<const unsigned char*>(data.data()), data.size());
    return hash.hex_digest();
}

std::string sha256_file_hex(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open file for hashing: " + path.string());
    }
    Sha256 hash;
    std::vector<char> buffer(1 << 16);
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize got = input.gcount();
        if (got > 0) {
            hash.update(reinterpret_cast<const unsigned char*>(buffer.data()),
                        static_cast<std::size_t>(got));
        }
    }
    if (input.bad()) {
        throw std::runtime_error("Failed while reading file for hashing: " + path.string());
    }
    return hash.hex_digest();
}

void atomic_write_file(const std::filesystem::path& path, std::string_view content) {
    const std::filesystem::path directory =
        path.has_parent_path() ? path.parent_path() : std::filesystem::path(".");
    static unsigned long counter = 0;
    const std::filesystem::path temp_path =
        directory / ("." + path.filename().string() + ".tmp." +
                     std::to_string(current_process_id()) + "." + std::to_string(counter++));

#if defined(_WIN32)
    {
        std::ofstream out(temp_path, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error("Failed to open temporary file for writing: " +
                                     temp_path.string());
        }
        out.write(content.data(), static_cast<std::streamsize>(content.size()));
        out.flush();
        if (!out) {
            out.close();
            std::error_code ignored;
            std::filesystem::remove(temp_path, ignored);
            throw std::runtime_error("Failed to write temporary file: " + temp_path.string());
        }
    }
    std::error_code rename_error;
    std::filesystem::rename(temp_path, path, rename_error);
    if (rename_error) {
        std::error_code ignored;
        std::filesystem::remove(temp_path, ignored);
        throw std::runtime_error("Failed to rename " + temp_path.string() + " to " +
                                 path.string() + ": " + rename_error.message());
    }
#else
    const int fd = ::open(temp_path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_TRUNC, 0644);
    if (fd < 0) {
        throw std::runtime_error("Failed to open temporary file for writing: " +
                                 temp_path.string() + ": " + std::strerror(errno));
    }
    auto fail = [&](const std::string& what) {
        const std::string reason = std::strerror(errno);
        ::close(fd);
        ::unlink(temp_path.c_str());
        throw std::runtime_error(what + " " + temp_path.string() + ": " + reason);
    };
    const char* cursor = content.data();
    std::size_t remaining = content.size();
    while (remaining > 0) {
        const ssize_t written = ::write(fd, cursor, remaining);
        if (written < 0) {
            if (errno == EINTR) continue;
            fail("Failed to write temporary file");
        }
        cursor += written;
        remaining -= static_cast<std::size_t>(written);
    }
    if (::fsync(fd) != 0) {
        fail("Failed to fsync temporary file");
    }
    if (::close(fd) != 0) {
        const std::string reason = std::strerror(errno);
        ::unlink(temp_path.c_str());
        throw std::runtime_error("Failed to close temporary file " + temp_path.string() +
                                 ": " + reason);
    }
    if (::rename(temp_path.c_str(), path.c_str()) != 0) {
        const std::string reason = std::strerror(errno);
        ::unlink(temp_path.c_str());
        throw std::runtime_error("Failed to rename " + temp_path.string() + " to " +
                                 path.string() + ": " + reason);
    }
    // Persist the directory entry. Best effort: some file systems refuse to
    // fsync a directory, and the data itself is already durable.
    const int dir_fd = ::open(directory.c_str(), O_RDONLY);
    if (dir_fd >= 0) {
        ::fsync(dir_fd);
        ::close(dir_fd);
    }
#endif
}

}  // namespace gmd
