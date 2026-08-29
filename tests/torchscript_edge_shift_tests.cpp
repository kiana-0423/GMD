// Direct validation of the periodic `edge_shift` tensor the TorchScript
// interface hands to an ML model.
//
// WHAT MAKES THIS DIFFERENT FROM THE EXISTING COVERAGE
//
// tests/box_image_flag_tests.cpp already pins down NeighborList::image_flags,
// the integer shift S that VerletNeighborBuilder stores per half-pair. That is
// a helper. It says nothing about what TorchScriptModelRuntimeAdapter does with
// S afterwards: the half-to-full-graph expansion, the reverse edge, the
// multiplication by the box lengths, the float32 conversion, the src/dst
// ordering inside edge_index, or whether the tensor that reaches forward() is
// the one that was built. Every test here observes edge_shift from INSIDE a
// real scripted model, reached through the production MLForceProvider ->
// TorchScriptModelRuntimeAdapter -> torch::jit path.
//
// THE CONTRACT (established by reading the production path, then confirmed here)
//
//   edge_index   int64   [2, 2E]   row 0 = src, row 1 = dst
//   edge_shift   float32 [2E, 3]   CARTESIAN translation in Angstrom, not an
//                                  integer image offset
//   device                         whatever the loaded module is on (CPU here)
//
//   For a directed edge (src, dst) with shift s, the periodic displacement is
//
//       d(src -> dst) = r_dst + s - r_src
//
//   and that equals the minimum-image separation. The builder stores
//   S = round((r_i - r_j) / L) for the half-pair (i, j), and the adapter emits
//   s = S * L for the edge i -> j and s = -S * L for the edge j -> i. Every
//   half-pair therefore produces exactly two directed edges, and the two
//   shifts are exact negatives of each other.
//
//   Sign convention, stated the way it is easy to get backwards: a NEGATIVE
//   x-shift on the edge i -> j means j's relevant periodic image sits one box
//   length in the -x direction from where j is stored.
//
//   The shift depends on the STORED coordinates, not only on the physical
//   configuration. The same dimer written with one atom wrapped into the cell
//   and written with it unwrapped produces different shifts -- and the same
//   displacement. That is the contract behaving correctly, not a defect, and
//   test_wrapped_and_unwrapped_agree_on_displacement() pins both halves of it.
//
//   All three directions are always periodic. gmd::Box carries three edge
//   lengths and nothing else; there is no per-axis periodicity flag anywhere in
//   the minimum-image or neighbour-construction path, so a non-periodic
//   direction cannot be requested. See test_all_three_axes_are_periodic().
//
// HOW THE MODEL IS OBSERVED
//
// The adapter returns only {"energy", "forces"}, so the model has to encode
// what it saw into those. Two scripted models are built here, from C++ via
// torch::jit::Module::define, so no Python and no checked-in .pt file is
// needed and nothing is written into the source tree:
//
//   graph observer  forces[i] = (number of edges with src == i,
//                                sum of species[dst],
//                                sum of species[dst]^2)
//   shift observer  forces[i] = sum over edges with src == i of
//                                 edge_shift * species[dst]
//                   energy    = sum over edges of (shift . w) * key,
//                               w = (1, 1e3, 1e6), key = 31*species[src] + species[dst]
//
// `species` carries a UNIQUE identity per atom, so every assertion is written
// against a stable physical identity rather than an index into whatever order
// the neighbour list happened to produce. Fixtures are built from well
// separated dimers, so each atom is the source of exactly one edge and
// forces[i] / species[j] recovers that edge's shift exactly -- a value, not a
// pass/fail. The graph observer establishes that one-edge-per-atom structure
// from the model's own point of view, so the recovery does not assume it.
//
// The expected shifts are written out from the fixture geometry and the box,
// never obtained from image_flags or apply_minimum_image. Reusing the helper
// under test to predict its own output would prove nothing.
//
// TOLERANCE. Every comparison of a shift component is EXACT (==). The box
// lengths are 10, 12 and 14, the species are small integers, and the shifts are
// integer multiples of the box lengths, so every value in the chain is exactly
// representable in float32 and in float64 and no rounding occurs anywhere. A
// tolerance would only hide a defect. The one place a tolerance appears is the
// energy checksum, and it is justified where it is used.

#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <span>
#include <string>
#include <vector>

#include <torch/script.h>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/composite_force_provider.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/force/ml_force_provider.hpp"
#include "gmd/force/torchscript_adapter.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/verlet_neighbor_builder.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[edge shift] " << message << '\n';
        ++failures;
    }
}

// --- the box -------------------------------------------------------------
//
// Non-cubic, every edge a different length, and every length exactly
// representable in float32 so the shift arithmetic is exact end to end.
constexpr double kLx = 10.0;
constexpr double kLy = 12.0;
constexpr double kLz = 14.0;
constexpr double kCutoff = 3.5;
constexpr double kSkin = 0.5;

// ===========================================================================
// Scripted observer models, built from C++ and saved next to the test binary
// ===========================================================================

// The shift observer, optionally mutated. The mutation is applied to the
// tensor the model received, which from the test's point of view is
// indistinguishable from the production path having delivered a mutated
// tensor -- which is exactly what a negative control needs.
enum class Mutation {
    None,
    ForceZero,      // all shifts zeroed
    ReverseSign,    // every shift negated
    PermuteAxes,    // (x, y, z) <- (z, x, y)
};

std::string mutation_source(Mutation mutation) {
    switch (mutation) {
        case Mutation::None:        return "";
        case Mutation::ForceZero:   return "    shift = shift * 0.0\n";
        case Mutation::ReverseSign: return "    shift = -shift\n";
        case Mutation::PermuteAxes:
            return "    shift = shift.index_select(1, torch.tensor([2, 0, 1]))\n";
    }
    return "";
}

std::string mutation_name(Mutation mutation) {
    switch (mutation) {
        case Mutation::None:        return "unmutated";
        case Mutation::ForceZero:   return "shifts forced to zero";
        case Mutation::ReverseSign: return "shift sign reversed";
        case Mutation::PermuteAxes: return "x/y/z components permuted";
    }
    return "";
}

// Writes a scripted module to `path` and returns it. Built with
// torch::jit::Module::define rather than exported from Python: the test then
// needs only LibTorch, and the artifact lives in the test's working directory,
// never in the source tree.
void write_shift_observer(const std::filesystem::path& path, Mutation mutation) {
    torch::jit::Module module("EdgeShiftObserver");
    module.register_attribute("local_cutoff", c10::FloatType::get(),
                              c10::IValue(static_cast<double>(kCutoff)));
    const std::string source = std::string(R"JIT(
def forward(self,
            species: Tensor,
            positions: Tensor,
            edge_index: Tensor,
            edge_shift: Tensor) -> Dict[str, Tensor]:
    src = edge_index[0]
    dst = edge_index[1]
    sp = species.double()
    shift = edge_shift.double()
)JIT") + mutation_source(mutation) + R"JIT(
    w = torch.tensor([1.0, 1000.0, 1000000.0]).double()
    key = sp.index_select(0, src) * 31.0 + sp.index_select(0, dst)
    energy = torch.sum(torch.matmul(shift, w) * key)
    tag = sp.index_select(0, dst).unsqueeze(1)
    forces = torch.zeros_like(positions).double()
    forces = forces.index_add(0, src, shift * tag)
    return {"energy": energy, "forces": forces.float()}
)JIT";
    module.define(source);
    module.eval();
    module.save(path.string());
}

// Reports the directed graph as the model sees it, so the structural
// assumptions below are established from inside forward() rather than by
// consulting the production neighbour list a second time.
void write_graph_observer(const std::filesystem::path& path) {
    torch::jit::Module module("EdgeGraphObserver");
    module.register_attribute("local_cutoff", c10::FloatType::get(),
                              c10::IValue(static_cast<double>(kCutoff)));
    module.define(R"JIT(
def forward(self,
            species: Tensor,
            positions: Tensor,
            edge_index: Tensor,
            edge_shift: Tensor) -> Dict[str, Tensor]:
    src = edge_index[0]
    dst = edge_index[1]
    sp = species.double()
    ones = torch.ones_like(src).double()
    neighbour = sp.index_select(0, dst)
    degree = torch.zeros_like(sp).index_add(0, src, ones)
    total = torch.zeros_like(sp).index_add(0, src, neighbour)
    squares = torch.zeros_like(sp).index_add(0, src, neighbour * neighbour)
    forces = torch.stack([degree, total, squares], 1)
    energy = torch.sum(degree)
    return {"energy": energy, "forces": forces.float()}
)JIT");
    module.eval();
    module.save(path.string());
}

// ===========================================================================
// Fixtures
// ===========================================================================

struct Atom {
    std::array<double, 3> position;
    int identity;  // unique species value, the stable handle for assertions
};

// One dimer and the shift its i -> j edge must carry, written out from the
// geometry rather than derived from the code under test.
struct Dimer {
    Atom first;
    Atom second;
    std::array<double, 3> expected_shift_first_to_second;
    std::string description;
};

gmd::System build_system(const std::vector<Atom>& atoms) {
    gmd::System system;
    system.resize(atoms.size(), atoms.size());
    gmd::Box box;
    box.set_lengths({kLx, kLy, kLz});
    system.set_box(box);
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        system.mutable_coordinates()[i] = atoms[i].position;
        system.mutable_masses()[i] = 1.0;
        system.mutable_atomic_numbers()[i] = atoms[i].identity;
    }
    return system;
}

struct Observation {
    double energy = 0.0;
    std::vector<std::array<double, 3>> forces;
};

// Runs the real provider over the real neighbour list and the real scripted
// module. Nothing here shortcuts the production path.
Observation observe(const std::vector<Atom>& atoms,
                    const std::filesystem::path& model,
                    bool through_composite = false) {
    gmd::System system = build_system(atoms);
    gmd::RuntimeContext runtime;

    gmd::VerletNeighborBuilder builder(kCutoff, kSkin);
    builder.initialize(system, runtime);

    auto ml = std::make_shared<gmd::MLForceProvider>(
        model, std::make_shared<gmd::TorchScriptModelRuntimeAdapter>());

    std::shared_ptr<gmd::ForceProvider> provider = ml;
    if (through_composite) {
        auto composite = std::make_shared<gmd::CompositeForceProvider>();
        composite->add(ml);
        provider = composite;
    }
    provider->initialize(runtime);

    const auto coordinates = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .step = 0,
        .time = 0.0,
        .coordinates = std::span<const gmd::Coordinate3D>(coordinates.data(),
                                                          coordinates.size()),
        .neighbor_list = &system.neighbor_list(),
    };
    gmd::ForceResult result;
    provider->compute(request, result, runtime);

    Observation observation;
    observation.energy = result.potential_energy;
    observation.forces.assign(result.forces.size(), {0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < result.forces.size(); ++i) {
        observation.forces[i] = {result.forces[i][0], result.forces[i][1],
                                 result.forces[i][2]};
    }
    check(result.success, "the ML evaluation did not succeed");
    return observation;
}

// Recovers each edge's shift from the shift observer's per-atom sums, keyed by
// the SPECIES of the two endpoints. Valid because every fixture is a set of
// well separated dimers, which the graph observer confirms independently.
std::map<std::pair<int, int>, std::array<double, 3>> recover_shifts(
        const std::vector<Atom>& atoms, const Observation& observation) {
    std::map<std::pair<int, int>, std::array<double, 3>> shifts;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        // The partner is the only other atom within the cutoff; identified from
        // the fixture's construction, which pairs atoms 2k and 2k+1.
        const std::size_t partner = (i % 2 == 0) ? i + 1 : i - 1;
        const double tag = static_cast<double>(atoms[partner].identity);
        std::array<double, 3> shift{};
        for (std::size_t d = 0; d < 3; ++d) {
            shift[d] = observation.forces[i][d] / tag;
        }
        shifts[{atoms[i].identity, atoms[partner].identity}] = shift;
    }
    return shifts;
}

bool shift_equals(const std::array<double, 3>& observed,
                  const std::array<double, 3>& expected) {
    // Exact. See the tolerance note in the file header.
    return observed[0] == expected[0] && observed[1] == expected[1]
        && observed[2] == expected[2];
}

std::string to_string(const std::array<double, 3>& v) {
    return "(" + std::to_string(v[0]) + ", " + std::to_string(v[1]) + ", "
         + std::to_string(v[2]) + ")";
}

// The dimer set. Each pair sits far from every other pair, so the graph is a
// perfect matching and no atom has two neighbours.
std::vector<Dimer> boundary_dimers() {
    return {
        // 1. No boundary crossing: both shifts must be exactly zero.
        {{{4.0, 5.0, 6.0}, 1}, {{5.0, 5.0, 6.0}, 2},
         {0.0, 0.0, 0.0}, "interior, no crossing"},

        // 2. Each face in turn. The first atom sits just inside the low face
        //    and the second just inside the high face, so the pair is 1 A apart
        //    through the boundary and S = -1 on that axis.
        {{{0.5, 5.0, 6.0}, 3}, {{9.5, 5.0, 6.0}, 4},
         {-kLx, 0.0, 0.0}, "across the x face, negative shift"},
        {{{4.0, 0.5, 6.0}, 5}, {{4.0, 11.5, 6.0}, 6},
         {0.0, -kLy, 0.0}, "across the y face, negative shift"},
        {{{4.0, 5.0, 0.5}, 7}, {{4.0, 5.0, 13.5}, 8},
         {0.0, 0.0, -kLz}, "across the z face, negative shift"},

        // 3. The same three faces with the atoms swapped, so the shift on the
        //    i -> j edge is POSITIVE. Without these, every non-zero expectation
        //    in the suite would have the same sign and a global sign error
        //    would survive.
        {{{9.5, 5.0, 6.0}, 9}, {{0.5, 5.0, 6.0}, 10},
         {kLx, 0.0, 0.0}, "across the x face, positive shift"},
        {{{4.0, 11.5, 6.0}, 11}, {{4.0, 0.5, 6.0}, 12},
         {0.0, kLy, 0.0}, "across the y face, positive shift"},
        {{{4.0, 5.0, 13.5}, 13}, {{4.0, 5.0, 0.5}, 14},
         {0.0, 0.0, kLz}, "across the z face, positive shift"},

        // 4. A corner, with all three components non-zero AND of mixed sign,
        //    so an axis permutation or a partial sign error cannot survive.
        {{{0.5, 11.5, 0.5}, 15}, {{9.5, 0.5, 13.5}, 16},
         {-kLx, kLy, -kLz}, "across the corner, mixed signs"},
    };
}

// A set of dimers that can all coexist in one cell: every atom of one dimer is
// further than r_list = cutoff + skin = 4.0 from every atom of every other,
// measured under the minimum image, while each dimer still crosses the face it
// is meant to cross. The closest inter-dimer approach is about 5.2 A.
//
// The dimers in boundary_dimers() deliberately overlap each other -- several of
// them occupy the same region of the cell -- so they are only ever evaluated one
// at a time. Combining them would produce atoms with two neighbours and the
// per-atom readout would stop being a per-edge readout. test_many_edges_at_once()
// checks the degree explicitly rather than trusting either fixture.
std::vector<Dimer> separated_dimers() {
    return {
        {{{0.5, 1.0, 1.0}, 21}, {{9.5, 1.0, 1.0}, 22},
         {-kLx, 0.0, 0.0}, "x face, negative"},
        {{{5.0, 0.5, 7.0}, 23}, {{5.0, 11.5, 7.0}, 24},
         {0.0, -kLy, 0.0}, "y face, negative"},
        {{{7.0, 6.0, 3.0}, 25}, {{7.0, 6.0, 4.0}, 26},
         {0.0, 0.0, 0.0}, "interior, no crossing"},
        {{{2.0, 8.0, 0.5}, 27}, {{2.0, 8.0, 13.5}, 28},
         {0.0, 0.0, -kLz}, "z face, negative"},
        {{{9.5, 4.0, 10.0}, 29}, {{0.5, 4.0, 10.0}, 30},
         {kLx, 0.0, 0.0}, "x face, positive"},
    };
}

std::vector<Atom> flatten(const std::vector<Dimer>& dimers) {
    std::vector<Atom> atoms;
    for (const auto& dimer : dimers) {
        atoms.push_back(dimer.first);
        atoms.push_back(dimer.second);
    }
    return atoms;
}

// ===========================================================================
// Tests
// ===========================================================================

// Each dimer is exercised on its own, which keeps the graph trivially a single
// edge and makes the per-atom readout an exact per-edge readout.
void test_each_boundary_crossing(const std::filesystem::path& shift_model,
                                 const std::filesystem::path& graph_model) {
    for (const auto& dimer : boundary_dimers()) {
        const std::vector<Atom> atoms = {dimer.first, dimer.second};

        // Structure first, from the model's own view of edge_index.
        const Observation graph = observe(atoms, graph_model);
        check(graph.forces.size() == 2,
              dimer.description + ": expected two atoms in the model output");
        for (std::size_t i = 0; i < 2; ++i) {
            const std::size_t partner = 1 - i;
            check(graph.forces[i][0] == 1.0,
                  dimer.description + ": atom " + std::to_string(i) +
                      " is the source of " + std::to_string(graph.forces[i][0]) +
                      " edges, expected exactly 1 -- every half-pair must produce "
                      "one edge in each direction");
            check(graph.forces[i][1] == static_cast<double>(atoms[partner].identity),
                  dimer.description + ": atom " + std::to_string(i) +
                      "'s edge points at species " + std::to_string(graph.forces[i][1]) +
                      ", expected " + std::to_string(atoms[partner].identity) +
                      " -- src/dst are not paired as edge_index claims");
        }

        // Then the values.
        const Observation shift = observe(atoms, shift_model);
        const auto recovered = recover_shifts(atoms, shift);

        const auto forward_key = std::make_pair(dimer.first.identity, dimer.second.identity);
        const auto reverse_key = std::make_pair(dimer.second.identity, dimer.first.identity);

        const std::array<double, 3>& forward = recovered.at(forward_key);
        const std::array<double, 3>& reverse = recovered.at(reverse_key);

        check(shift_equals(forward, dimer.expected_shift_first_to_second),
              dimer.description + ": edge " + std::to_string(dimer.first.identity) +
                  " -> " + std::to_string(dimer.second.identity) + " carries shift " +
                  to_string(forward) + ", expected " +
                  to_string(dimer.expected_shift_first_to_second));

        const std::array<double, 3> expected_reverse = {
            -dimer.expected_shift_first_to_second[0],
            -dimer.expected_shift_first_to_second[1],
            -dimer.expected_shift_first_to_second[2]};
        check(shift_equals(reverse, expected_reverse),
              dimer.description + ": the reverse edge carries shift " +
                  to_string(reverse) + ", expected the exact negative " +
                  to_string(expected_reverse));

        // The displacement the contract promises, reconstructed from the shift
        // the model actually received and the stored coordinates.
        for (std::size_t d = 0; d < 3; ++d) {
            const double displacement =
                dimer.second.position[d] + forward[d] - dimer.first.position[d];
            const double length = (d == 0) ? kLx : (d == 1) ? kLy : kLz;
            check(std::fabs(displacement) <= 0.5 * length,
                  dimer.description + ": r_dst + shift - r_src = " +
                      std::to_string(displacement) + " on axis " + std::to_string(d) +
                      " is not the minimum image (half box = " +
                      std::to_string(0.5 * length) + ")");
        }
    }
}

// All eight dimers at once: several edges, several distinct shifts, and the
// per-atom readout still exact because the graph is a perfect matching.
void test_many_edges_at_once(const std::filesystem::path& shift_model,
                             const std::filesystem::path& graph_model) {
    const auto dimers = separated_dimers();
    const auto atoms = flatten(dimers);

    const Observation graph = observe(atoms, graph_model);
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        check(graph.forces[i][0] == 1.0,
              "multi-edge fixture: atom with species " +
                  std::to_string(atoms[i].identity) + " is the source of " +
                  std::to_string(graph.forces[i][0]) +
                  " edges, expected 1; the dimers are not isolated and the "
                  "per-atom readout would not be a per-edge readout");
    }

    const Observation shift = observe(atoms, shift_model);
    const auto recovered = recover_shifts(atoms, shift);
    for (const auto& dimer : dimers) {
        const auto& forward = recovered.at({dimer.first.identity, dimer.second.identity});
        check(shift_equals(forward, dimer.expected_shift_first_to_second),
              "multi-edge fixture, " + dimer.description + ": shift " +
                  to_string(forward) + ", expected " +
                  to_string(dimer.expected_shift_first_to_second));
    }
}

// Atom ordering must not matter. The assertions are keyed by species, so a
// permutation that changed which shift belongs to which physical pair would be
// caught; one that merely renumbered indices must not disturb anything.
void test_atom_order_permutation(const std::filesystem::path& shift_model) {
    const auto dimers = separated_dimers();
    const auto atoms = flatten(dimers);

    const Observation baseline = observe(atoms, shift_model);
    const auto expected = recover_shifts(atoms, baseline);

    // A reversal, and a rotation by one dimer. Both keep the 2k / 2k+1 pairing
    // that recover_shifts() relies on while changing every atom index.
    std::vector<std::vector<Atom>> permutations;
    {
        std::vector<Atom> reversed;
        for (std::size_t pair = dimers.size(); pair-- > 0;) {
            reversed.push_back(atoms[2 * pair]);
            reversed.push_back(atoms[2 * pair + 1]);
        }
        permutations.push_back(reversed);
    }
    {
        std::vector<Atom> rotated;
        for (std::size_t pair = 0; pair < dimers.size(); ++pair) {
            const std::size_t source = (pair + 2) % dimers.size();
            rotated.push_back(atoms[2 * source]);
            rotated.push_back(atoms[2 * source + 1]);
        }
        permutations.push_back(rotated);
    }
    // Swapping the two atoms within every dimer, which also flips the sign of
    // every expected shift -- a stronger permutation than a reordering.
    {
        std::vector<Atom> swapped;
        for (std::size_t pair = 0; pair < dimers.size(); ++pair) {
            swapped.push_back(atoms[2 * pair + 1]);
            swapped.push_back(atoms[2 * pair]);
        }
        permutations.push_back(swapped);
    }

    for (std::size_t index = 0; index < permutations.size(); ++index) {
        const Observation permuted = observe(permutations[index], shift_model);
        const auto recovered = recover_shifts(permutations[index], permuted);
        check(recovered.size() == expected.size(),
              "permutation " + std::to_string(index) +
                  ": recovered a different number of edges");
        for (const auto& [key, value] : expected) {
            const auto found = recovered.find(key);
            check(found != recovered.end(),
                  "permutation " + std::to_string(index) + ": edge " +
                      std::to_string(key.first) + " -> " + std::to_string(key.second) +
                      " disappeared; edges are being identified by position, not by "
                      "physical identity");
            if (found == recovered.end()) continue;
            check(shift_equals(found->second, value),
                  "permutation " + std::to_string(index) + ": edge " +
                      std::to_string(key.first) + " -> " + std::to_string(key.second) +
                      " changed from " + to_string(value) + " to " +
                      to_string(found->second));
        }
        check(permuted.energy == baseline.energy,
              "permutation " + std::to_string(index) +
                  ": the whole-graph energy checksum changed, so the edge set is "
                  "not permutation invariant");
    }
}

// The shift is a function of the STORED coordinates. Writing the same physical
// dimer with one atom unwrapped changes the shift and must not change the
// displacement -- both halves are asserted, because checking only the
// displacement would pass for an implementation that always reported zero.
void test_wrapped_and_unwrapped_agree_on_displacement(
        const std::filesystem::path& shift_model) {
    // Wrapped: atom B at x = 9.5, one box length to the right of its image.
    const std::vector<Atom> wrapped = {{{0.5, 5.0, 6.0}, 1}, {{9.5, 5.0, 6.0}, 2}};
    // Unwrapped: the same image written directly, at x = -0.5.
    const std::vector<Atom> unwrapped = {{{0.5, 5.0, 6.0}, 1}, {{-0.5, 5.0, 6.0}, 2}};

    const auto wrapped_shift = recover_shifts(wrapped, observe(wrapped, shift_model));
    const auto unwrapped_shift = recover_shifts(unwrapped, observe(unwrapped, shift_model));

    const std::array<double, 3> from_wrapped = wrapped_shift.at({1, 2});
    const std::array<double, 3> from_unwrapped = unwrapped_shift.at({1, 2});

    check(shift_equals(from_wrapped, {-kLx, 0.0, 0.0}),
          "wrapped placement: expected shift (-10, 0, 0), got " + to_string(from_wrapped));
    check(shift_equals(from_unwrapped, {0.0, 0.0, 0.0}),
          "unwrapped placement: expected shift (0, 0, 0) because the stored "
          "coordinate is already the correct image, got " + to_string(from_unwrapped));

    // Different shifts, same physics.
    check(!shift_equals(from_wrapped, from_unwrapped),
          "the two placements produced the same shift, so this fixture is not "
          "testing the dependence on the stored coordinates");

    for (std::size_t d = 0; d < 3; ++d) {
        const double displacement_wrapped =
            wrapped[1].position[d] + from_wrapped[d] - wrapped[0].position[d];
        const double displacement_unwrapped =
            unwrapped[1].position[d] + from_unwrapped[d] - unwrapped[0].position[d];
        check(displacement_wrapped == displacement_unwrapped,
              "axis " + std::to_string(d) + ": the wrapped placement gives "
              "displacement " + std::to_string(displacement_wrapped) +
              " and the unwrapped placement " + std::to_string(displacement_unwrapped) +
              "; physically equivalent coordinates must give the same displacement");
    }
}

// There is no non-periodic direction to test. gmd::Box holds three edge lengths
// and nothing else, apply_minimum_image() wraps all three unconditionally, and
// VerletNeighborBuilder applies periodic cell wrapping on all three axes. The
// per-axis `periodic` flags on DomainDecomposition govern MPI halo exchange,
// not neighbour construction, and never reach the ML path (which rejects MPI).
//
// So the honest test is the documented behaviour: every axis is periodic, and a
// crossing on any one of them produces a shift. That is asserted positively
// here so the limitation cannot quietly change.
void test_all_three_axes_are_periodic(const std::filesystem::path& shift_model) {
    const std::array<double, 3> lengths = {kLx, kLy, kLz};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        std::array<double, 3> low = {4.0, 5.0, 6.0};
        std::array<double, 3> high = {4.0, 5.0, 6.0};
        low[axis] = 0.5;
        high[axis] = lengths[axis] - 0.5;

        const std::vector<Atom> atoms = {{low, 1}, {high, 2}};
        const auto recovered = recover_shifts(atoms, observe(atoms, shift_model));
        const std::array<double, 3>& shift = recovered.at({1, 2});

        std::array<double, 3> expected = {0.0, 0.0, 0.0};
        expected[axis] = -lengths[axis];
        check(shift_equals(shift, expected),
              "axis " + std::to_string(axis) +
                  " is expected to be periodic (no axis can be disabled), so a "
                  "crossing must produce shift " + to_string(expected) + ", got " +
                  to_string(shift));
    }
}

// A provider with a fixed, known answer, so the composite arithmetic below is
// about propagation rather than about the numbers.
class StubProvider final : public gmd::ForceProvider {
public:
    explicit StubProvider(double value) : value_(value) {}
    std::string_view name() const noexcept override { return "stub"; }
    void initialize(gmd::RuntimeContext&) override {}
    void finalize(gmd::RuntimeContext&) override {}
    void compute(const gmd::ForceRequest& request, gmd::ForceResult& result,
                 gmd::RuntimeContext&) override {
        result.success = true;
        result.potential_energy = value_;
        result.forces.assign(request.coordinates.size(),
                             gmd::Force3D{value_, value_, value_});
        result.virial = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        result.virial_valid = true;
    }
private:
    double value_;
};

// The same contract must reach an ML child inside a composite. The composite
// hands every child the same ForceRequest and accumulates into one result, so
// this is where a stale or overwritten tensor would show up.
void test_composite_with_ml_child(const std::filesystem::path& shift_model) {
    const auto dimers = separated_dimers();
    const auto atoms = flatten(dimers);

    const Observation standalone = observe(atoms, shift_model);
    const auto expected = recover_shifts(atoms, standalone);

    // 1. ML alone inside a composite: identical to running it directly.
    const Observation wrapped = observe(atoms, shift_model, /*through_composite=*/true);
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        for (std::size_t d = 0; d < 3; ++d) {
            check(wrapped.forces[i][d] == standalone.forces[i][d],
                  "composite: atom " + std::to_string(i) + " component " +
                      std::to_string(d) + " changed from " +
                      std::to_string(standalone.forces[i][d]) + " to " +
                      std::to_string(wrapped.forces[i][d]) +
                      " merely by being wrapped in a CompositeForceProvider");
        }
    }
    check(wrapped.energy == standalone.energy,
          "composite: the energy changed when the ML provider was wrapped");

    // 2. ML between two other children, run twice. The ML contribution is
    //    recovered by subtracting the known stub forces, and must still decode
    //    to exactly the same shifts. A tensor left over from a previous
    //    evaluation, or a result buffer a sibling overwrote, breaks this.
    gmd::System system = build_system(atoms);
    gmd::RuntimeContext runtime;
    gmd::VerletNeighborBuilder builder(kCutoff, kSkin);
    builder.initialize(system, runtime);

    auto composite = std::make_shared<gmd::CompositeForceProvider>();
    composite->add(std::make_shared<StubProvider>(2.0));
    composite->add(std::make_shared<gmd::MLForceProvider>(
        shift_model, std::make_shared<gmd::TorchScriptModelRuntimeAdapter>()));
    composite->add(std::make_shared<StubProvider>(-0.5));
    composite->initialize(runtime);

    const auto coordinates = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .step = 0,
        .time = 0.0,
        .coordinates = std::span<const gmd::Coordinate3D>(coordinates.data(),
                                                          coordinates.size()),
        .neighbor_list = &system.neighbor_list(),
    };

    for (int evaluation = 0; evaluation < 3; ++evaluation) {
        gmd::ForceResult result;
        composite->compute(request, result, runtime);
        check(result.success, "composite evaluation " + std::to_string(evaluation) +
                                  " did not succeed");

        Observation isolated;
        isolated.energy = result.potential_energy - 1.5;  // 2.0 + (-0.5)
        isolated.forces.assign(result.forces.size(), {0.0, 0.0, 0.0});
        for (std::size_t i = 0; i < result.forces.size(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) {
                isolated.forces[i][d] = result.forces[i][d] - 1.5;
            }
        }

        const auto recovered = recover_shifts(atoms, isolated);
        for (const auto& [key, value] : expected) {
            check(shift_equals(recovered.at(key), value),
                  "composite evaluation " + std::to_string(evaluation) + ": edge " +
                      std::to_string(key.first) + " -> " + std::to_string(key.second) +
                      " decoded to " + to_string(recovered.at(key)) + ", expected " +
                      to_string(value));
        }
        check(std::fabs(isolated.energy - standalone.energy)
                  < 1.0e-9 * std::max(1.0, std::fabs(standalone.energy)),
              "composite evaluation " + std::to_string(evaluation) +
                  ": the isolated ML energy is " + std::to_string(isolated.energy) +
                  ", expected " + std::to_string(standalone.energy));
    }
}

// ===========================================================================
// Negative controls
// ===========================================================================

// Each mutated model reports what a mutated production path would have
// delivered. The suite must reject every one of them; a control that passes
// means the corresponding assertion is not actually load-bearing.
void test_negative_controls(const std::filesystem::path& directory) {
    const auto dimers = separated_dimers();
    const auto atoms = flatten(dimers);

    for (const Mutation mutation : {Mutation::ForceZero, Mutation::ReverseSign,
                                    Mutation::PermuteAxes}) {
        const std::filesystem::path model =
            directory / ("edge_shift_mutant_" +
                         std::to_string(static_cast<int>(mutation)) + ".pt");
        write_shift_observer(model, mutation);

        const auto recovered = recover_shifts(atoms, observe(atoms, model));

        // A zero shift is a fixed point of all three mutations -- negating it,
        // zeroing it or permuting its components all leave it unchanged -- so
        // the interior dimer is not detectable by construction and is excluded
        // from the count rather than weakening the assertion for the others.
        // Every edge that CAN be detected must be.
        int detectable = 0;
        int mismatches = 0;
        for (const auto& dimer : dimers) {
            const auto& expected = dimer.expected_shift_first_to_second;
            const bool is_zero = shift_equals(expected, {0.0, 0.0, 0.0});
            if (!is_zero) ++detectable;
            const auto& observed = recovered.at({dimer.first.identity,
                                                 dimer.second.identity});
            if (!shift_equals(observed, expected)) ++mismatches;
        }
        check(detectable > 0,
              "negative control (" + mutation_name(mutation) +
                  "): no fixture carries a non-zero shift, so nothing could be "
                  "detected and this control is vacuous");
        check(mismatches == detectable,
              "negative control (" + mutation_name(mutation) + ") rejected " +
                  std::to_string(mismatches) + " of the " +
                  std::to_string(detectable) +
                  " edges it should have rejected; the assertions are not catching "
                  "that defect on every affected edge");
        std::cout << "    negative control, " << mutation_name(mutation) << ": "
                  << mismatches << " of " << detectable
                  << " detectable edges rejected (" << (dimers.size() - detectable)
                  << " zero-shift edge is a fixed point of this mutation)\n";
    }
}

}  // namespace

int main() {
    // Artifacts go in the current working directory, which CTest sets to the
    // build tree. Nothing is written into the source tree and no absolute path
    // is baked in.
    const std::filesystem::path directory = std::filesystem::current_path();
    const std::filesystem::path shift_model = directory / "edge_shift_observer.pt";
    const std::filesystem::path graph_model = directory / "edge_graph_observer.pt";

    try {
        write_shift_observer(shift_model, Mutation::None);
        write_graph_observer(graph_model);
    } catch (const std::exception& error) {
        std::cerr << "[edge shift] could not build the scripted observer models: "
                  << error.what() << '\n';
        return 1;
    }

    std::cout << "TorchScript edge_shift, observed inside a scripted model\n";
    test_each_boundary_crossing(shift_model, graph_model);
    test_many_edges_at_once(shift_model, graph_model);
    test_atom_order_permutation(shift_model);
    test_wrapped_and_unwrapped_agree_on_displacement(shift_model);
    test_all_three_axes_are_periodic(shift_model);
    test_composite_with_ml_child(shift_model);
    std::cout << "  negative controls\n";
    test_negative_controls(directory);

    if (failures != 0) {
        std::cerr << "TorchScript edge shift tests failed: " << failures << '\n';
        return 1;
    }
    std::cout << "TorchScript edge shift tests passed\n";
    return 0;
}
