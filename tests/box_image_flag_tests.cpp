// Regression tests for box-dimension validation and for the neighbor list's
// per-pair image flags.
//
// Box validation: a zero or negative edge length used to flow straight into
// periodic wrapping (division by zero) and into the neighbor builder's
// ghost-shift loops (which advance by one box length per iteration and so would
// never terminate). Box::set_lengths()/System::set_box() now reject it.
//
// Image flags: image_flags[k] is the integer shift S for the pair stored at
// neighbors[k]. The convention every consumer relies on -- including the
// TorchScript adapter, which feeds edge_shift = S * L to the model -- is
//
//     r_j + S * L - r_i  ==  minimum_image(r_j - r_i)
//
// These tests pin that convention down across each periodic boundary.

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#ifdef GMD_ENABLE_MPI
#include "gmd/parallel/mpi_environment.hpp"
#endif
#include "gmd/system/box.hpp"
#include "gmd/system/minimum_image.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/verlet_neighbor_builder.hpp"

namespace {

int failures = 0;

void check(bool value, const std::string& message) {
    if (!value) {
        std::cerr << "[box/image] " << message << '\n';
        ++failures;
    }
}

// Runs `body` and reports whether it threw std::invalid_argument.
template <typename Body>
bool throws_invalid_argument(Body body) {
    try {
        body();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

// --- Box validation -------------------------------------------------------

void test_box_validation() {
    for (int dim = 0; dim < 3; ++dim) {
        {   // Zero length on each axis in turn.
            std::array<double, 3> lengths = {10.0, 10.0, 10.0};
            lengths[static_cast<std::size_t>(dim)] = 0.0;
            check(throws_invalid_argument([&] { gmd::Box b; b.set_lengths(lengths); }),
                  "zero box length on axis " + std::to_string(dim) + " must throw");
        }
        {   // Negative length on each axis in turn.
            std::array<double, 3> lengths = {10.0, 10.0, 10.0};
            lengths[static_cast<std::size_t>(dim)] = -4.0;
            check(throws_invalid_argument([&] { gmd::Box b; b.set_lengths(lengths); }),
                  "negative box length on axis " + std::to_string(dim) + " must throw");
        }
    }

    // Non-finite lengths are rejected too: NaN would silently poison every
    // wrapped coordinate rather than looping forever.
    check(throws_invalid_argument([] {
              gmd::Box b;
              b.set_lengths({std::nan(""), 10.0, 10.0});
          }),
          "NaN box length must throw");
    check(throws_invalid_argument([] {
              gmd::Box b;
              b.set_lengths({std::numeric_limits<double>::infinity(), 10.0, 10.0});
          }),
          "infinite box length must throw");

    // A valid box still works and still derives its half-lengths.
    gmd::Box good;
    good.set_lengths({10.0, 12.0, 14.0});
    check(good.half_lengths[0] == 5.0 && good.half_lengths[1] == 6.0 &&
              good.half_lengths[2] == 7.0,
          "valid box must still set half_lengths");

    // System::set_box() copies a whole Box and so bypasses set_lengths(); it
    // must validate as well, otherwise a default-constructed Box (all zeros)
    // would reach the neighbor builder.
    check(throws_invalid_argument([] {
              gmd::System system;
              system.resize(2, 2);
              gmd::Box degenerate;  // default-constructed: all lengths zero
              system.set_box(degenerate);
          }),
          "System::set_box must reject a degenerate default-constructed box");
}

// --- Image flags ----------------------------------------------------------

constexpr double kBoxLength = 12.0;
constexpr double kCutoff = 3.0;
constexpr double kSkin = 1.0;

// Builds a two-atom system, runs the neighbor builder, and verifies the stored
// image flag against the minimum-image displacement.
void check_pair_image_flag(const std::array<double, 3>& r_i,
                           const std::array<double, 3>& r_j,
                           const std::string& label) {
    gmd::System system;
    system.resize(2, 2);
    gmd::Box box;
    box.set_lengths({kBoxLength, kBoxLength, kBoxLength});
    system.set_box(box);
    system.mutable_masses()[0] = 1.0;
    system.mutable_masses()[1] = 1.0;
    system.mutable_coordinates()[0] = r_i;
    system.mutable_coordinates()[1] = r_j;

    gmd::VerletNeighborBuilder builder(kCutoff, kSkin);
    gmd::RuntimeContext runtime;
    builder.rebuild(system, runtime, nullptr);

    const auto& nl = system.neighbor_list();
    check(nl.valid, label + ": neighbor list must be valid");
    check(nl.image_flags.size() == nl.neighbors.size(),
          label + ": image_flags must be parallel to neighbors");

    bool found = false;
    for (std::size_t i = 0; i < 2; ++i) {
        const int start = nl.offsets[i];
        for (int k = 0; k < nl.counts[i]; ++k) {
            const std::size_t idx = static_cast<std::size_t>(start + k);
            const std::size_t j = static_cast<std::size_t>(nl.neighbors[idx]);
            const std::array<int, 3>& S = nl.image_flags[idx];

            // Expected: the minimum-image displacement from i to j.
            std::array<double, 3> expected = {
                system.coordinates()[j][0] - system.coordinates()[i][0],
                system.coordinates()[j][1] - system.coordinates()[i][1],
                system.coordinates()[j][2] - system.coordinates()[i][2],
            };
            gmd::apply_minimum_image(expected, box);

            // Convention under test: r_j + S * L - r_i.
            for (std::size_t d = 0; d < 3; ++d) {
                const double actual = system.coordinates()[j][d]
                                    + static_cast<double>(S[d]) * box.lengths[d]
                                    - system.coordinates()[i][d];
                check(std::abs(actual - expected[d]) < 1.0e-9,
                      label + ": image flag mismatch on axis " + std::to_string(d) +
                          " (r_j + S*L - r_i = " + std::to_string(actual) +
                          ", minimum image = " + std::to_string(expected[d]) + ")");
            }
            found = true;
        }
    }
    check(found, label + ": expected the pair to be in the neighbor list");
}

void test_image_flags() {
    // A pair well inside the box: no wrapping, S must be zero.
    check_pair_image_flag({6.0, 6.0, 6.0}, {7.5, 6.0, 6.0}, "interior pair");

    // One pair crossing each periodic boundary in turn. The atoms sit either
    // side of the x/y/z face, so the true separation goes through the boundary
    // and the stored shift has to account for it.
    check_pair_image_flag({0.5, 6.0, 6.0}, {11.5, 6.0, 6.0}, "pair crossing x");
    check_pair_image_flag({6.0, 0.5, 6.0}, {6.0, 11.5, 6.0}, "pair crossing y");
    check_pair_image_flag({6.0, 6.0, 0.5}, {6.0, 6.0, 11.5}, "pair crossing z");

    // The opposite crossing direction, to catch a sign error that happens to
    // cancel in one direction only.
    check_pair_image_flag({11.5, 6.0, 6.0}, {0.5, 6.0, 6.0}, "pair crossing x (reversed)");
    check_pair_image_flag({6.0, 11.5, 6.0}, {6.0, 0.5, 6.0}, "pair crossing y (reversed)");
    check_pair_image_flag({6.0, 6.0, 11.5}, {6.0, 6.0, 0.5}, "pair crossing z (reversed)");

    // A corner pair crossing all three boundaries at once.
    check_pair_image_flag({0.5, 0.5, 0.5}, {11.5, 11.5, 11.5}, "pair crossing xyz corner");
}

}  // namespace

int main(int argc, char** argv) {
#ifdef GMD_ENABLE_MPI
    // An MPI-enabled build still runs this as a single-process test, but the
    // force providers issue collectives that require an initialised MPI.
    gmd::MpiEnvironment environment(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif
    test_box_validation();
    test_image_flags();

    if (failures != 0) {
        std::cerr << "[box/image] " << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "[box/image] all checks passed\n";
    return 0;
}
