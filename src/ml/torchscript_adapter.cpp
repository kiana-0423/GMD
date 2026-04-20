#ifdef GMD_ENABLE_TORCH

#include "gmd/ml/torchscript_adapter.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <torch/script.h>

#include "gmd/runtime/runtime_context.hpp"
#include "gmd/system/system.hpp"

namespace gmd {

// ---------------------------------------------------------------------------
// Private implementation
// ---------------------------------------------------------------------------
struct TorchScriptModelRuntimeAdapter::Impl {
    torch::jit::script::Module module;
    torch::Device device{torch::kCPU};
    float cutoff{0.0f};
    bool loaded{false};
};

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------
TorchScriptModelRuntimeAdapter::TorchScriptModelRuntimeAdapter() noexcept
    : impl_(std::make_unique<Impl>()) {}

TorchScriptModelRuntimeAdapter::~TorchScriptModelRuntimeAdapter() = default;

std::string_view TorchScriptModelRuntimeAdapter::name() const noexcept {
    return "torchscript_model_runtime_adapter";
}

float TorchScriptModelRuntimeAdapter::cutoff() const noexcept {
    return impl_->cutoff;
}

void TorchScriptModelRuntimeAdapter::load_model(
        const std::filesystem::path& model_path,
        RuntimeContext& /*runtime*/) {
    try {
        impl_->module = torch::jit::load(model_path.string(), impl_->device);
    } catch (const c10::Error& e) {
        throw std::runtime_error(
            "TorchScript: failed to load model from " + model_path.string()
            + ": " + e.what());
    }
    impl_->module.eval();

    // Read cutoff from model attribute.  The attribute may be stored as
    // float or double by the export script.
    try {
        auto attr = impl_->module.attr("local_cutoff");
        if (attr.isDouble()) {
            impl_->cutoff = static_cast<float>(attr.toDouble());
        } else if (attr.isInt()) {
            impl_->cutoff = static_cast<float>(attr.toInt());
        } else {
            impl_->cutoff = static_cast<float>(attr.toTensor().item<double>());
        }
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("TorchScript: model is missing 'local_cutoff' attribute: ")
            + e.what());
    }

    impl_->loaded = true;
}

void TorchScriptModelRuntimeAdapter::evaluate(
        const ModelEvaluationRequest& request,
        ModelEvaluationResult& result,
        RuntimeContext& /*runtime*/) {
    result.success = false;

    if (!impl_->loaded) {
        return;
    }

    const std::size_t N = request.coordinates.size();
    if (N == 0) {
        result.total_energy = 0.0;
        result.forces.clear();
        result.success = true;
        return;
    }

    // -----------------------------------------------------------------------
    // 1. Build input tensors
    // -----------------------------------------------------------------------

    // species: int64 [N]
    auto species_t = torch::zeros({static_cast<long>(N)}, torch::kInt64);
    {
        auto acc = species_t.accessor<int64_t, 1>();
        for (std::size_t i = 0; i < N; ++i) {
            acc[static_cast<long>(i)] = (request.atomic_numbers.size() == N)
                                            ? static_cast<int64_t>(request.atomic_numbers[i])
                                            : 0LL;
        }
    }

    // positions: float32 [N, 3]
    auto pos_t = torch::zeros({static_cast<long>(N), 3}, torch::kFloat32);
    {
        auto acc = pos_t.accessor<float, 2>();
        for (std::size_t i = 0; i < N; ++i) {
            acc[static_cast<long>(i)][0] = static_cast<float>(request.coordinates[i][0]);
            acc[static_cast<long>(i)][1] = static_cast<float>(request.coordinates[i][1]);
            acc[static_cast<long>(i)][2] = static_cast<float>(request.coordinates[i][2]);
        }
    }

    // -----------------------------------------------------------------------
    // 2. Build edge_index [2, 2E] and edge_shift [2E, 3] from NeighborList
    //    Half-pairs (j > i) are expanded to a directed full graph.
    // -----------------------------------------------------------------------
    const NeighborList* nl = request.neighbor_list;

    torch::Tensor edge_index_t;
    torch::Tensor edge_shift_t;

    if (nl != nullptr && nl->valid && !nl->neighbors.empty()) {
        const std::size_t E_half = nl->neighbors.size();
        const std::size_t E_full = 2 * E_half;

        const double* Lx = request.box ? &request.box->lengths[0] : nullptr;

        // Collect all directed edges in plain vectors first.
        std::vector<int64_t> ei_src, ei_dst;
        std::vector<float> es_data;
        ei_src.reserve(E_full);
        ei_dst.reserve(E_full);
        es_data.reserve(E_full * 3);

        for (std::size_t i = 0; i < N; ++i) {
            const int offset = nl->offsets[i];
            for (int k = 0; k < nl->counts[i]; ++k) {
                const int j = nl->neighbors[static_cast<std::size_t>(offset + k)];
                const std::array<int, 3>& S =
                    nl->image_flags[static_cast<std::size_t>(offset + k)];

                // Edge i → j
                ei_src.push_back(static_cast<int64_t>(i));
                ei_dst.push_back(static_cast<int64_t>(j));
                if (Lx) {
                    es_data.push_back(static_cast<float>(S[0]) * static_cast<float>(Lx[0]));
                    es_data.push_back(static_cast<float>(S[1]) * static_cast<float>(Lx[1]));
                    es_data.push_back(static_cast<float>(S[2]) * static_cast<float>(Lx[2]));
                } else {
                    es_data.push_back(0.0f);
                    es_data.push_back(0.0f);
                    es_data.push_back(0.0f);
                }

                // Edge j → i (negated shift)
                ei_src.push_back(static_cast<int64_t>(j));
                ei_dst.push_back(static_cast<int64_t>(i));
                if (Lx) {
                    es_data.push_back(-static_cast<float>(S[0]) * static_cast<float>(Lx[0]));
                    es_data.push_back(-static_cast<float>(S[1]) * static_cast<float>(Lx[1]));
                    es_data.push_back(-static_cast<float>(S[2]) * static_cast<float>(Lx[2]));
                } else {
                    es_data.push_back(0.0f);
                    es_data.push_back(0.0f);
                    es_data.push_back(0.0f);
                }
            }
        }

        const long E = static_cast<long>(E_full);

        // Stack src/dst into [2, E]
        auto src_t = torch::from_blob(ei_src.data(), {E},
                                       torch::TensorOptions().dtype(torch::kInt64))
                         .clone();
        auto dst_t = torch::from_blob(ei_dst.data(), {E},
                                       torch::TensorOptions().dtype(torch::kInt64))
                         .clone();
        edge_index_t = torch::stack({src_t, dst_t}, /*dim=*/0);

        edge_shift_t = torch::from_blob(es_data.data(), {E, 3},
                                         torch::TensorOptions().dtype(torch::kFloat32))
                           .clone();
    } else {
        // No neighbor list available: empty edges (non-periodic or single-atom).
        edge_index_t = torch::zeros({2, 0}, torch::kInt64);
        edge_shift_t = torch::zeros({0, 3}, torch::kFloat32);
    }

    // Move tensors to the model device if needed.
    species_t     = species_t.to(impl_->device);
    pos_t         = pos_t.to(impl_->device);
    edge_index_t  = edge_index_t.to(impl_->device);
    edge_shift_t  = edge_shift_t.to(impl_->device);

    // -----------------------------------------------------------------------
    // 3. Call model forward()
    // -----------------------------------------------------------------------
    c10::IValue output;
    try {
        output = impl_->module.forward({species_t, pos_t, edge_index_t, edge_shift_t});
    } catch (const c10::Error& e) {
        throw std::runtime_error(
            std::string("TorchScript: model.forward() failed: ") + e.what());
    }

    // -----------------------------------------------------------------------
    // 4. Extract energy and forces from the returned dict
    // -----------------------------------------------------------------------
    auto dict = output.toGenericDict();

    const double energy = dict.at("energy").toTensor().item<double>();

    auto forces_t = dict.at("forces").toTensor().cpu().contiguous();
    const auto forces_acc = forces_t.accessor<float, 2>();

    result.total_energy = energy;
    result.forces.resize(N);
    for (std::size_t i = 0; i < N; ++i) {
        result.forces[i][0] = static_cast<double>(forces_acc[static_cast<long>(i)][0]);
        result.forces[i][1] = static_cast<double>(forces_acc[static_cast<long>(i)][1]);
        result.forces[i][2] = static_cast<double>(forces_acc[static_cast<long>(i)][2]);
    }
    result.success = true;
}

void TorchScriptModelRuntimeAdapter::unload_model(RuntimeContext& /*runtime*/) {
    impl_->loaded = false;
    impl_->cutoff = 0.0f;
}

}  // namespace gmd

#endif  // GMD_ENABLE_TORCH
