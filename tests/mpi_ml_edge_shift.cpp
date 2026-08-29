// The ML path under MPI: proving the limitation rather than pretending to test
// around it.
//
// `edge_shift` never reaches a model under MPI domain decomposition, because
// MLForceProvider refuses to run there at all. Both initialize() and compute()
// throw when RuntimeContext::size() > 1. So there is no multi-rank edge_shift
// behaviour to validate, and a "np=2 edge_shift test" would either be testing a
// path that cannot execute or quietly running everything on one rank and
// calling it parallel.
//
// What is worth testing is the refusal itself, because it is what stands
// between a user and a silently wrong answer. The reason it exists:
// TorchScriptModelRuntimeAdapter builds the graph from the LOCAL neighbour list
// and hands the model every atom in `request.coordinates`, which under domain
// decomposition is local atoms plus ghosts. Nothing in that path decides which
// rank owns the energy of an edge that straddles a boundary, and nothing
// defines how deep the halo must be for a model whose receptive field is
// several message-passing layers wide rather than one cutoff. Summing each
// rank's model output would double-count every cross-boundary contribution.
//
// This test therefore asserts:
//   np = 1  the provider runs normally -- so the refusal is about rank count,
//           not about the provider being broken;
//   np > 1  initialize() and compute() both throw, on EVERY rank, with a
//           message that names the limitation.
//
// It deliberately uses a stub adapter rather than TorchScript, so the proof
// holds in every MPI build whether or not LibTorch is available: the refusal
// happens before any adapter is touched.

#include <cstddef>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include <mpi.h>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/composite_force_provider.hpp"
#include "gmd/force/ml_force_provider.hpp"
#include "gmd/force/model_runtime_adapter.hpp"
#include "gmd/parallel/mpi_environment.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/system.hpp"

namespace {

int failures = 0;
int global_rank = 0;
int global_size = 1;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[mpi ml edge shift][rank " << global_rank << "] " << message
                  << '\n';
        ++failures;
    }
}

// Stands in for a real model runtime. Reaching this at all under MPI would
// itself be the defect, so it records whether it was ever asked to do anything.
class RecordingAdapter final : public gmd::ModelRuntimeAdapter {
public:
    bool loaded = false;
    bool evaluated = false;

    std::string_view name() const noexcept override { return "recording_stub"; }
    float cutoff() const noexcept override { return 3.5f; }
    void load_model(const std::filesystem::path&, gmd::RuntimeContext&) override {
        loaded = true;
    }
    void unload_model(gmd::RuntimeContext&) override {}
    void evaluate(const gmd::ModelEvaluationRequest& request,
                  gmd::ModelEvaluationResult& result,
                  gmd::RuntimeContext&) override {
        evaluated = true;
        result.success = true;
        result.total_energy = 1.0;
        result.forces.assign(request.coordinates.size(), gmd::Force3D{0.0, 0.0, 0.0});
    }
};

// A dimer straddling the x face, split across ranks when there is more than
// one: a genuine cross-rank periodic edge, which is exactly the configuration
// whose ownership the ML path does not define.
gmd::System make_cross_boundary_system() {
    const std::array<std::array<double, 3>, 2> positions = {{
        {0.5, 5.0, 6.0}, {9.5, 5.0, 6.0}}};

    // Rank 0 owns the low-x atom, rank 1 (when present) the high-x one.
    std::vector<int> owned;
    for (int tag = 0; tag < 2; ++tag) {
        const int owner = (global_size > 1) ? tag % global_size : 0;
        if (owner == global_rank) owned.push_back(tag);
    }

    gmd::System system;
    system.resize(owned.size(), owned.size());
    gmd::Box box;
    box.set_lengths({10.0, 12.0, 14.0});
    system.set_box(box);
    for (std::size_t i = 0; i < owned.size(); ++i) {
        const auto tag = static_cast<std::size_t>(owned[i]);
        system.mutable_coordinates()[i] = positions[tag];
        system.mutable_masses()[i] = 1.0;
        system.mutable_atomic_numbers()[i] = owned[i] + 1;
        system.mutable_atom_tags()[i] = owned[i];
        system.mutable_atom_owners()[i] = global_rank;
    }
    return system;
}

bool mentions_mpi_limitation(const std::string& message) {
    return message.find("MPI") != std::string::npos;
}

void test_single_rank_runs_and_multi_rank_refuses() {
    gmd::System system = make_cross_boundary_system();
    gmd::RuntimeContext runtime;

    auto adapter = std::make_shared<RecordingAdapter>();
    gmd::MLForceProvider provider("unused_by_the_stub.pt", adapter);

    const auto coordinates = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .step = 0,
        .time = 0.0,
        .coordinates = std::span<const gmd::Coordinate3D>(coordinates.data(),
                                                          coordinates.size()),
        .neighbor_list = nullptr,
    };

    check(runtime.size() == global_size,
          "RuntimeContext does not agree with MPI_Comm_size, so the gate under "
          "test is not seeing the real rank count");

    if (global_size == 1) {
        bool threw = false;
        std::string message;
        try {
            provider.initialize(runtime);
            gmd::ForceResult result;
            provider.compute(request, result, runtime);
            check(result.success, "np=1: the ML provider should evaluate normally");
        } catch (const std::exception& error) {
            threw = true;
            message = error.what();
        }
        check(!threw, "np=1: the ML provider must run; it threw \"" + message + "\"");
        check(adapter->evaluated,
              "np=1: the adapter was never asked to evaluate, so this run proves "
              "nothing about the provider working at one rank");
        return;
    }

    // np > 1: both entry points must refuse, on this rank.
    std::string initialize_message;
    bool initialize_threw = false;
    try {
        provider.initialize(runtime);
    } catch (const std::exception& error) {
        initialize_threw = true;
        initialize_message = error.what();
    }
    check(initialize_threw,
          "np=" + std::to_string(global_size) +
              ": MLForceProvider::initialize() must refuse to run under MPI domain "
              "decomposition; it returned normally");
    check(mentions_mpi_limitation(initialize_message),
          "np=" + std::to_string(global_size) +
              ": the refusal message does not mention MPI: \"" +
              initialize_message + "\"");

    std::string compute_message;
    bool compute_threw = false;
    try {
        gmd::ForceResult result;
        provider.compute(request, result, runtime);
    } catch (const std::exception& error) {
        compute_threw = true;
        compute_message = error.what();
    }
    check(compute_threw,
          "np=" + std::to_string(global_size) +
              ": MLForceProvider::compute() must refuse as well; initialize() alone "
              "is not a gate, because a provider can be constructed and computed "
              "without an explicit initialize()");
    check(mentions_mpi_limitation(compute_message),
          "np=" + std::to_string(global_size) +
              ": the compute refusal message does not mention MPI: \"" +
              compute_message + "\"");

    // The adapter must never have been reached, so no partial or duplicated
    // model contribution can exist on any rank.
    check(!adapter->evaluated,
          "np=" + std::to_string(global_size) +
              ": the model adapter WAS evaluated under MPI. Every rank would then "
              "contribute its own model energy for the same cross-boundary edge, "
              "which is precisely the double count the refusal exists to prevent");
    check(!adapter->loaded,
          "np=" + std::to_string(global_size) +
              ": the model was loaded despite the refusal");
}

// The refusal must not be escapable by burying the provider in a composite.
void test_composite_with_ml_child_also_refuses() {
    if (global_size == 1) return;

    gmd::System system = make_cross_boundary_system();
    gmd::RuntimeContext runtime;
    auto adapter = std::make_shared<RecordingAdapter>();

    auto composite = std::make_shared<gmd::CompositeForceProvider>();
    composite->add(std::make_shared<gmd::MLForceProvider>("unused.pt", adapter));

    bool threw = false;
    try {
        composite->initialize(runtime);
    } catch (const std::exception&) {
        threw = true;
    }
    check(threw,
          "np=" + std::to_string(global_size) +
              ": a CompositeForceProvider containing an ML child must propagate the "
              "MPI refusal rather than swallowing it");
    check(!adapter->evaluated,
          "np=" + std::to_string(global_size) +
              ": the composite reached the model adapter under MPI");
}

// Every rank must reach the same verdict. A refusal that fired on some ranks
// and not others would deadlock a real run at the next collective.
void test_every_rank_agrees() {
    int local = failures;
    int minimum = 0;
    int maximum = 0;
    MPI_Allreduce(&local, &minimum, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local, &maximum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    check(minimum == maximum,
          "ranks disagree about the ML/MPI verdict (min " + std::to_string(minimum) +
              ", max " + std::to_string(maximum) +
              "); a refusal that fires on some ranks only would deadlock at the "
              "next collective");
}

}  // namespace

int main(int argc, char** argv) {
    gmd::MpiEnvironment environment(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    test_single_rank_runs_and_multi_rank_refuses();
    test_composite_with_ml_child_also_refuses();
    test_every_rank_agrees();

    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (global_rank == 0) {
        if (total == 0) {
            std::cout << "[mpi ml edge shift] all checks passed on " << global_size
                      << " rank(s)"
                      << (global_size > 1
                              ? "; the ML path refuses MPI, so edge_shift never "
                                "reaches a model there"
                              : "; the ML path runs normally at one rank")
                      << '\n';
        } else {
            std::cerr << "[mpi ml edge shift] " << total << " check(s) failed on "
                      << global_size << " rank(s)\n";
        }
    }
    return total == 0 ? 0 : 1;
}
