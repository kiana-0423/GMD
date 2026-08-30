// Collective validation of the atom tags the random draw is keyed on.
//
// The tag became load-bearing when the draw stopped depending on storage
// order, so a duplicate now means two physical atoms receive the same velocity
// and a negative one has no defined mapping. Both are rejected -- but under
// MPI, "rejected" has a second requirement beyond throwing: EVERY rank must
// throw. A rank that noticed nothing locally and carried on would be left
// waiting in the next collective while the others unwound, and a bad input
// would present as a hang.
//
// That is not hypothetical. The first version of this validation threw the
// moment a rank saw a negative tag, before the gather that finds duplicates,
// which is exactly the deadlock described above. The negative-tag scan is now
// a separate MPI_Allreduce ahead of the gather, so all ranks agree before any
// of them unwinds.
//
// Duplicates are the more interesting case, because a duplicate can span
// ranks: two ranks each holding one atom tagged 5 are locally consistent and
// only a global check can see it. The np=2 and np=4 cases below are built that
// way deliberately.

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "gmd/system/box.hpp"
#include "gmd/system/initializer.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;
int global_rank = 0;
int global_size = 1;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[mpi velocity tags][rank " << global_rank << "] " << message << '\n';
        ++failures;
    }
}

constexpr std::size_t kTotalAtoms = 8;
constexpr double kTargetTemperature = 300.0;

double mass_for(std::size_t tag) {
    static const std::array<double, 8> kMasses = {
        1.008, 12.011, 15.999, 4.0026, 6.941, 10.811, 14.007, 32.065};
    return kMasses[tag % kMasses.size()];
}

// At np>=4 the last rank owns nothing.
std::vector<int> owned_tags() {
    std::vector<int> owned;
    const int contributing = (global_size >= 4) ? global_size - 1 : global_size;
    for (int tag = 0; tag < static_cast<int>(kTotalAtoms); ++tag) {
        if (tag % contributing == global_rank) owned.push_back(tag);
    }
    return owned;
}

gmd::System build(const std::vector<int>& tags) {
    gmd::System system;
    system.resize(tags.size(), tags.size());
    gmd::Box box;
    box.set_lengths({50.0, 50.0, 50.0});
    system.set_box(box);
    for (std::size_t slot = 0; slot < tags.size(); ++slot) {
        const auto identity = static_cast<std::size_t>(std::abs(tags[slot]));
        system.mutable_masses()[slot] = mass_for(identity % kTotalAtoms);
        system.mutable_coordinates()[slot] = {2.0 + 3.0 * static_cast<double>(slot),
                                              3.0, 4.0};
        system.mutable_atom_tags()[slot] = tags[slot];
        system.mutable_atom_owners()[slot] = global_rank;
    }
    return system;
}

// Runs initialization and reports, collectively, whether THIS rank threw.
// Returns the number of ranks that threw.
int ranks_that_threw(std::vector<int> tags) {
    gmd::System system = build(tags);
    gmd::VelocityInitializer initializer(20260830u);
    int threw = 0;
    try {
        initializer.initialize(system, kTargetTemperature,
                               gmd::VelocityInitMode::Random, true);
    } catch (const std::exception&) {
        threw = 1;
    }
    int total = 0;
    MPI_Allreduce(&threw, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return total;
}

void test_valid_tags_are_accepted() {
    const int total = ranks_that_threw(owned_tags());
    check(total == 0,
          "a correctly tagged system was rejected on " + std::to_string(total) +
              " of " + std::to_string(global_size) + " rank(s)");
}

void test_empty_rank_participates() {
    std::vector<int> tags = owned_tags();
    if (global_size >= 4 && global_rank == global_size - 1) {
        check(tags.empty(),
              "this fixture expects the last rank to own no atoms at np>=4");
    }
    // Reaching here at all is the test: a rank owning nothing must still enter
    // the tag gather and the reductions inside initialize().
    const int total = ranks_that_threw(tags);
    check(total == 0, "an empty rank made a valid system fail");
}

void test_cross_rank_duplicate_is_rejected_everywhere() {
    // Every rank claims tag 0. At np=1 that is a local duplicate only if the
    // rank holds more than one atom, so the fixture gives every rank two atoms:
    // one genuine and one colliding.
    std::vector<int> tags = owned_tags();
    tags.push_back(0);
    const int total = ranks_that_threw(tags);
    check(total == global_size,
          "a duplicated atom tag was rejected on only " + std::to_string(total) +
              " of " + std::to_string(global_size) +
              " rank(s). Every rank must throw, or the ones that did not would "
              "wait alone in the next collective");
}

void test_duplicate_spanning_two_ranks_is_rejected() {
    // The case a local check cannot see: each rank's own tags are unique, but
    // two ranks share one. At np=1 there is no second rank, so the fixture
    // degenerates to a valid system and must be accepted -- which is worth
    // asserting, because it says the check is not simply rejecting everything.
    std::vector<int> tags = owned_tags();
    if (global_size > 1 && global_rank == global_size - 1) {
        tags.push_back(0);          // tag 0 is owned by rank 0
    }
    const int total = ranks_that_threw(tags);
    if (global_size == 1) {
        check(total == 0, "the np=1 degenerate case should be a valid system");
    } else {
        check(total == global_size,
              "a duplicate spanning two ranks was rejected on only " +
                  std::to_string(total) + " of " + std::to_string(global_size) +
                  " rank(s); each rank's own tags are unique, so only a global "
                  "check can see this one");
    }
}

void test_negative_tag_is_rejected_everywhere() {
    // On ONE rank only. The negative scan is local, so this is the case that
    // deadlocks if the rank that sees it throws before the collective.
    std::vector<int> tags = owned_tags();
    if (global_rank == 0 && !tags.empty()) {
        tags[0] = -3;
    }
    const int total = ranks_that_threw(tags);
    check(total == global_size,
          "a negative atom tag on rank 0 was rejected on only " +
              std::to_string(total) + " of " + std::to_string(global_size) +
              " rank(s). Reaching this assertion at all means no rank hung, "
              "which is the other half of the requirement");
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    test_valid_tags_are_accepted();
    test_empty_rank_participates();
    test_cross_rank_duplicate_is_rejected_everywhere();
    test_duplicate_spanning_two_ranks_is_rejected();
    test_negative_tag_is_rejected_everywhere();

    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_rank == 0) {
        if (total == 0) {
            std::cout << "[mpi velocity tags] all checks passed on " << global_size
                      << " rank(s)\n";
        } else {
            std::cerr << "[mpi velocity tags] " << total << " check(s) failed on "
                      << global_size << " rank(s)\n";
        }
    }
    MPI_Finalize();
    return total == 0 ? 0 : 1;
}
