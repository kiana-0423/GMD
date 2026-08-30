#pragma once

// Counter-based random numbers keyed on physical identity.
//
// WHY NOT A CONVENTIONAL GENERATOR. A std::mt19937 advanced in a loop assigns
// its draws by traversal order. That is fine in serial and wrong under domain
// decomposition: each rank starts from the same seeded state and gives its own
// first local atom the generator's first draws, so which numbers a physical
// atom receives depends on the rank count, the decomposition, the local
// storage order, and any migration that happened first. Nothing about the atom
// itself enters.
//
// Here the draw is a pure function of a key instead:
//
//     value = f(seed, stream, identity, component, counter)
//
// with no state carried between atoms and no dependence on the order they are
// visited. The same atom gets the same number whether it is alone on one rank
// or one of thousands, first in the array or last.
//
// THE MIXER is SplitMix64 (Steele, Lea and Flood, "Fast splittable
// pseudorandom number generators", OOPSLA 2014), the finalizer used by
// java.util.SplittableRandom and by Vigna's reference implementation. It is
// two multiply-xorshift rounds with documented constants and full 64-bit
// avalanche, and it is exact integer arithmetic: identical on every platform,
// in every optimisation mode, at every rank count.
//
//   random_bits(key, i) reproduces the SplitMix64 stream seeded at `key`: the
//   generator's i-th output is mix(key + (i+1)*gamma), which is what the
//   reference implementation computes by adding gamma to its state before each
//   mix.
//
// KEY DERIVATION absorbs each input through a full mixing round, so every input
// avalanches over the whole key and no two distinct inputs collide except by
// 64-bit accident. The inputs are absorbed in a fixed order, and the initial
// value is an arbitrary domain constant so that a key of all zeros is not the
// identity.
//
// UNIFORMS land in the OPEN interval (0,1):
//
//     u = ((bits >> 12) + 0.5) / 2^52    in  [2^-53, 1 - 2^-53]
//
// so log(u) is always finite. A closed or half-open mapping would eventually
// hand Box-Muller a zero and produce infinity. The divisor is a power of two
// and the numerator is a plain bit extraction, so there is no modulo bias and
// the division is exact.
//
// TWELVE bits are discarded, not eleven, and the difference matters. With
// 53 bits the largest value is (2^53 - 1) + 0.5, which is NOT representable --
// doubles are spaced 1.0 at that magnitude -- so it rounds half-to-even up to
// 2^53 and the quotient comes out as exactly 1.0, reintroducing the endpoint
// the offset was there to remove. With 52 bits the spacing at the top of the
// range is 0.5, the sum is exact, and the largest quotient is 1 - 2^-53. The
// cost is one bit of resolution out of 52, which is nothing next to a variate
// that then goes through a logarithm and a cosine.
//
// NORMALS use Box-Muller. Each component draws its OWN pair of uniforms from
// its own key, rather than taking the cosine and sine of a shared pair: a
// shared pair would make a component's two outputs exactly dependent, and
// reusing one uniform across components is the correlation this is meant to
// avoid. Half the entropy per pair is discarded, which costs nothing that
// matters here.
//
// REPRODUCIBILITY. The integer path -- key derivation, the mixer, the uniform
// mapping -- is bit-exact everywhere. sqrt, log and cos are not guaranteed
// bit-identical across libm implementations, so the guarantee is: bitwise
// reproducible for a given platform and standard library across rank counts,
// decompositions, storage orders, build types and optimisation levels; and
// numerically reproducible to within a few ulp across platforms. sqrt is
// correctly rounded by IEEE-754, and log and cos are typically within an ulp,
// but that is a property of the library rather than a promise made here.
// tests/keyed_random_tests.cpp pins exact reference vectors so any change in
// the integer path is caught immediately.

#include <cmath>
#include <cstdint>

namespace gmd {

// Independent draw streams. A new random feature takes a new value here rather
// than sharing an existing one, so that adding it cannot shift the numbers an
// existing feature produces.
enum class RandomStream : std::uint64_t {
    VelocityInitialization = 1,
};

// The odd increment SplitMix64 adds to its state: 2^64 / phi, rounded to odd.
inline constexpr std::uint64_t kSplitMixGamma = 0x9E3779B97F4A7C15ULL;

// SplitMix64's finalizer. Two multiply-xorshift rounds.
inline constexpr std::uint64_t splitmix64_mix(std::uint64_t z) noexcept {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// Folds one more input into a key.
inline constexpr std::uint64_t key_absorb(std::uint64_t key,
                                          std::uint64_t value) noexcept {
    return splitmix64_mix(key + kSplitMixGamma + value);
}

// The key for one (seed, stream, identity, component) tuple. `identity` is the
// stable global atom tag, never an array position or a rank.
inline constexpr std::uint64_t random_key(std::uint64_t seed,
                                          std::uint64_t stream,
                                          std::uint64_t identity,
                                          std::uint64_t component) noexcept {
    // Fractional bits of pi: an arbitrary constant, present so that an
    // all-zero input does not produce an all-zero key.
    std::uint64_t key = 0x243F6A8885A308D3ULL;
    key = key_absorb(key, seed);
    key = key_absorb(key, stream);
    key = key_absorb(key, identity);
    key = key_absorb(key, component);
    return key;
}

// The i-th output of the SplitMix64 stream seeded at `key`.
inline constexpr std::uint64_t random_bits(std::uint64_t key,
                                           std::uint64_t index) noexcept {
    return splitmix64_mix(key + (index + 1ULL) * kSplitMixGamma);
}

// 2^52. See the note above on why 52 bits rather than 53.
inline constexpr double kUniformDenominator = 4503599627370496.0;
inline constexpr int kUniformDiscardedBits = 12;

// Uniform on the OPEN interval (0,1). See the note above on why open.
inline constexpr double uniform_open_unit(std::uint64_t bits) noexcept {
    return (static_cast<double>(bits >> kUniformDiscardedBits) + 0.5) *
           (1.0 / kUniformDenominator);
}

inline constexpr double kTwoPi = 6.283185307179586476925286766559;

// One standard normal variate for this (seed, stream, identity, component).
inline double standard_normal(std::uint64_t seed,
                              RandomStream stream,
                              std::uint64_t identity,
                              std::uint64_t component) noexcept {
    const std::uint64_t key =
        random_key(seed, static_cast<std::uint64_t>(stream), identity, component);
    const double u1 = uniform_open_unit(random_bits(key, 0));
    const double u2 = uniform_open_unit(random_bits(key, 1));
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(kTwoPi * u2);
}

}  // namespace gmd
