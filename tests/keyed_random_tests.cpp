// The keyed random generator, pinned against independently derived vectors.
//
// gmd/core/keyed_random.hpp exists so that a physical atom's random draw is a
// function of its identity rather than of the order something happened to walk
// an array. That property is only worth anything if the mapping itself is
// stable, so this file pins it exactly.
//
// The reference values below were produced by an independent reimplementation
// of the documented algorithm in Python -- a different language, arbitrary-
// precision integers, and no shared code with the header -- and are reproduced
// in tests/keyed_random_reference.py so the derivation can be re-run rather
// than taken on trust. They are not captured output: had the C++ and the Python
// disagreed, the vectors would have exposed it.
//
// The integer path is exact everywhere: key derivation, the SplitMix64 mixer
// and the uniform mapping are integer arithmetic and a power-of-two division.
// Those are checked for bitwise equality. The normal transform goes through
// sqrt, log and cos, which are not guaranteed bit-identical across libm
// implementations, so those vectors are checked to a few ulp and the file says
// so rather than pretending otherwise.

#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "gmd/core/keyed_random.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[keyed random] " << message << '\n';
        ++failures;
    }
}

std::string hex(std::uint64_t value) {
    std::ostringstream out;
    out << "0x" << std::hex << std::uppercase << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

std::string number(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

// --- exact key vectors -----------------------------------------------------

struct KeyVector {
    std::uint64_t seed;
    std::uint64_t stream;
    std::uint64_t identity;
    std::uint64_t component;
    std::uint64_t expected;
};

void test_key_reference_vectors() {
    static const std::array<KeyVector, 8> kVectors = {{
        {0u, 0u, 0u, 0u, 0x4534183A9817A9A3ULL},
        {1u, 1u, 0u, 0u, 0x6C07BF2F695B8072ULL},
        {20260830u, 1u, 0u, 0u, 0x452F533D779C3A30ULL},
        {20260830u, 1u, 0u, 1u, 0xAEBC99A1CA4C3561ULL},
        {20260830u, 1u, 0u, 2u, 0x716154D692A37F53ULL},
        {20260830u, 1u, 1u, 0u, 0x40966F4BABBCEB12ULL},
        {20260830u, 1u, 4294967295u, 2u, 0xEBE6F8B20AD1813FULL},
        {4294967295u, 1u, 123456789u, 1u, 0x7ECFA976750DE3D6ULL},
    }};
    for (const auto& v : kVectors) {
        const std::uint64_t produced =
            gmd::random_key(v.seed, v.stream, v.identity, v.component);
        check(produced == v.expected,
              "random_key(" + std::to_string(v.seed) + ", " + std::to_string(v.stream) +
                  ", " + std::to_string(v.identity) + ", " + std::to_string(v.component) +
                  ") = " + hex(produced) + ", expected " + hex(v.expected));
    }

    // An all-zero input must not produce an all-zero key: that is what the
    // domain constant in random_key is for.
    check(gmd::random_key(0, 0, 0, 0) != 0ULL,
          "an all-zero input produced an all-zero key");
}

void test_stream_reference_vectors() {
    const std::uint64_t key = gmd::random_key(20260830u, 1u, 0u, 0u);
    check(gmd::random_bits(key, 0) == 0xCF13FFC908274D4CULL,
          "random_bits(key, 0) = " + hex(gmd::random_bits(key, 0)));
    check(gmd::random_bits(key, 1) == 0x4F2A87459ACD50F4ULL,
          "random_bits(key, 1) = " + hex(gmd::random_bits(key, 1)));
    check(gmd::random_bits(key, 0) != gmd::random_bits(key, 1),
          "successive stream outputs are identical");
}

// --- the uniform mapping ---------------------------------------------------

void test_uniform_is_strictly_inside_the_unit_interval() {
    // The endpoints are where this matters: a zero would make log(u) infinite
    // and a one would make it exactly zero, and both are reachable if the
    // mapping is half-open.
    const double lowest = gmd::uniform_open_unit(0ULL);
    const double highest = gmd::uniform_open_unit(~0ULL);
    check(lowest > 0.0, "uniform_open_unit(0) is " + number(lowest) + ", not > 0");
    check(highest < 1.0, "uniform_open_unit(~0) is " + number(highest) + ", not < 1");
    check(std::isfinite(std::log(lowest)),
          "log of the smallest uniform is not finite");
    check(std::log(highest) < 0.0,
          "log of the largest uniform is not strictly negative, so the variate "
          "it produces would be exactly zero");

    // Exact, because both are computed by exact arithmetic.
    check(lowest == 1.1102230246251565e-16,
          "the smallest uniform is " + number(lowest));
    check(highest == 0.9999999999999999,
          "the largest uniform is " + number(highest));

    // Discarding eleven bits instead of twelve would put the top of the range
    // at exactly 1.0. This states why the header discards twelve.
    const double naive = (static_cast<double>((~0ULL) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
    check(naive == 1.0,
          "the 53-bit mapping was expected to round up to exactly 1.0, which is "
          "the reason the header uses 52; it produced " + number(naive));
}

void test_uniform_has_no_modulo_bias() {
    // The mapping is a bit shift and a power-of-two division, so the only way
    // to bias it is to get the arithmetic wrong. A coarse histogram over a
    // large sample is enough to catch that.
    constexpr int kBuckets = 16;
    constexpr std::uint64_t kSamples = 200000;
    std::array<std::uint64_t, kBuckets> histogram{};
    const std::uint64_t key = gmd::random_key(4242u, 1u, 7u, 1u);
    for (std::uint64_t i = 0; i < kSamples; ++i) {
        const double u = gmd::uniform_open_unit(gmd::random_bits(key, i));
        check(u > 0.0 && u < 1.0, "a uniform escaped (0,1): " + number(u));
        auto bucket = static_cast<std::size_t>(u * kBuckets);
        if (bucket >= kBuckets) bucket = kBuckets - 1;
        ++histogram[bucket];
    }
    const double expected = static_cast<double>(kSamples) / kBuckets;
    // 5 sigma on a multinomial bucket of this size is about 3.3 percent.
    for (int b = 0; b < kBuckets; ++b) {
        const double deviation =
            std::abs(static_cast<double>(histogram[static_cast<std::size_t>(b)]) - expected) /
            expected;
        check(deviation < 0.05,
              "uniform bucket " + std::to_string(b) + " holds " +
                  std::to_string(histogram[static_cast<std::size_t>(b)]) +
                  " of " + std::to_string(kSamples) + " samples, " +
                  number(deviation * 100.0) + " percent from uniform");
    }
}

// --- normals ---------------------------------------------------------------

struct NormalVector {
    std::uint64_t identity;
    std::uint64_t component;
    double expected;
};

void test_normal_reference_vectors() {
    static const std::array<NormalVector, 6> kVectors = {{
        {0u, 0u, -0.23686721361008206},
        {0u, 1u, -2.0443383893440688},
        {0u, 2u, 1.3107738684555905},
        {1u, 0u, -0.765644214806505},
        {7u, 2u, -0.648866582295988},
        {23u, 1u, 2.967502265772108},
    }};
    for (const auto& v : kVectors) {
        const double produced = gmd::standard_normal(
            20260830u, gmd::RandomStream::VelocityInitialization, v.identity, v.component);
        // A few ulp, not equality: sqrt is correctly rounded by IEEE-754 but
        // log and cos are library-dependent. On this platform they agree
        // exactly; the bound is what the algorithm can promise elsewhere.
        const double tolerance = 8.0 * 2.220446049250313e-16 * std::abs(v.expected);
        check(std::abs(produced - v.expected) <= tolerance,
              "standard_normal(identity " + std::to_string(v.identity) + ", component " +
                  std::to_string(v.component) + ") = " + number(produced) +
                  ", expected " + number(v.expected));
    }
}

void test_components_are_independent_draws() {
    // The failure this rules out is a single hash output reused for x, y and z,
    // or a shared uniform pair with cos and sin. Either would make components
    // equal or exactly related.
    for (std::uint64_t identity = 0; identity < 64; ++identity) {
        const double x = gmd::standard_normal(
            20260830u, gmd::RandomStream::VelocityInitialization, identity, 0);
        const double y = gmd::standard_normal(
            20260830u, gmd::RandomStream::VelocityInitialization, identity, 1);
        const double z = gmd::standard_normal(
            20260830u, gmd::RandomStream::VelocityInitialization, identity, 2);
        check(x != y && y != z && x != z,
              "atom " + std::to_string(identity) +
                  " received equal components: " + number(x) + ", " + number(y) +
                  ", " + number(z));
        // sin^2 + cos^2 = 1 would show up as x^2 + y^2 being suspiciously
        // related for a shared pair; check they are not negatives either.
        check(x != -y && y != -z && x != -z,
              "atom " + std::to_string(identity) +
                  " received exactly opposite components");
    }

    // Sample correlation between components, over enough atoms to mean
    // something. Independent draws give |r| ~ 1/sqrt(n).
    constexpr std::uint64_t kAtoms = 20000;
    double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (std::uint64_t identity = 0; identity < kAtoms; ++identity) {
        const double x = gmd::standard_normal(
            20260830u, gmd::RandomStream::VelocityInitialization, identity, 0);
        const double y = gmd::standard_normal(
            20260830u, gmd::RandomStream::VelocityInitialization, identity, 1);
        sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
    }
    const auto n = static_cast<double>(kAtoms);
    const double covariance = sxy / n - (sx / n) * (sy / n);
    const double sigma_x = std::sqrt(sxx / n - (sx / n) * (sx / n));
    const double sigma_y = std::sqrt(syy / n - (sy / n) * (sy / n));
    const double correlation = covariance / (sigma_x * sigma_y);
    std::cout << "  x/y correlation over " << kAtoms << " atoms  " << correlation << '\n';
    // 5/sqrt(n) is about 0.035.
    check(std::abs(correlation) < 0.035,
          "components x and y are correlated: r = " + number(correlation));
}

void test_normals_have_the_right_moments() {
    constexpr std::uint64_t kAtoms = 200000;
    double sum = 0.0, square_sum = 0.0;
    double lowest = 0.0, highest = 0.0;
    for (std::uint64_t identity = 0; identity < kAtoms; ++identity) {
        const double z = gmd::standard_normal(
            20260830u, gmd::RandomStream::VelocityInitialization, identity, 0);
        check(std::isfinite(z), "a non-finite variate appeared at identity " +
                                    std::to_string(identity));
        sum += z;
        square_sum += z * z;
        lowest = std::min(lowest, z);
        highest = std::max(highest, z);
    }
    const auto n = static_cast<double>(kAtoms);
    const double mean = sum / n;
    const double variance = square_sum / n - mean * mean;
    std::cout << "  normal mean " << mean << "  variance " << variance
              << "  range [" << lowest << ", " << highest << "]\n";
    // Standard error of the mean is 1/sqrt(n) ~ 0.0022; of the variance,
    // sqrt(2/n) ~ 0.0032. Five sigma each.
    check(std::abs(mean) < 0.012,
          "the sample mean is " + number(mean) + ", not zero");
    check(std::abs(variance - 1.0) < 0.017,
          "the sample variance is " + number(variance) + ", not one");
    // A Box-Muller pass that lost its logarithm would be bounded by 1.
    check(highest > 3.0 && lowest < -3.0,
          "the sample never reached +/-3 sigma over " + std::to_string(kAtoms) +
              " draws, so the tails are missing");
}

void test_streams_are_separated() {
    // A second stream must not reproduce the first's numbers for the same
    // identity: that is what stops a future random feature from silently
    // shifting velocity initialization.
    for (std::uint64_t identity = 0; identity < 32; ++identity) {
        const std::uint64_t a = gmd::random_key(20260830u, 1u, identity, 0u);
        const std::uint64_t b = gmd::random_key(20260830u, 2u, identity, 0u);
        check(a != b, "streams 1 and 2 produced the same key for identity " +
                          std::to_string(identity));
    }
}

void test_neighbouring_inputs_decorrelate() {
    // Adjacent seeds, tags and components must not produce adjacent keys: a
    // weak mixer would leave visible structure and neighbouring atoms would get
    // similar velocities.
    std::vector<std::uint64_t> keys;
    for (std::uint64_t identity = 0; identity < 256; ++identity) {
        keys.push_back(gmd::random_key(1u, 1u, identity, 0u));
    }
    std::size_t low_hamming = 0;
    for (std::size_t i = 1; i < keys.size(); ++i) {
        const int distance = __builtin_popcountll(keys[i] ^ keys[i - 1]);
        // Two independent 64-bit values differ in 32 bits on average; fewer
        // than 12 would be a strong sign of structure.
        if (distance < 12) ++low_hamming;
    }
    check(low_hamming == 0,
          std::to_string(low_hamming) +
              " adjacent identities produced keys differing in fewer than 12 bits");
}

}  // namespace

int main() {
    std::cout << "[keyed random] pinning the keyed generator\n";
    test_key_reference_vectors();
    test_stream_reference_vectors();
    test_uniform_is_strictly_inside_the_unit_interval();
    test_uniform_has_no_modulo_bias();
    test_normal_reference_vectors();
    test_components_are_independent_draws();
    test_normals_have_the_right_moments();
    test_streams_are_separated();
    test_neighbouring_inputs_decorrelate();

    if (failures == 0) {
        std::cout << "[keyed random] all checks passed\n";
        return 0;
    }
    std::cerr << "[keyed random] " << failures << " check(s) failed\n";
    return 1;
}
