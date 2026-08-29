#pragma once

#ifdef GMD_ENABLE_TORCH

#include <filesystem>
#include <memory>
#include <string_view>

#include "gmd/force/model_runtime_adapter.hpp"

// No forward declaration of the TorchScript module type is needed, and none is
// possible: `torch::jit::script` is a namespace ALIAS for `torch::jit`, so
// `struct script::Module;` is an elaborated-type-specifier carrying a
// nested-name-specifier, which is ill-formed. It never compiled. The handle
// lives behind the opaque `Impl` below, which is the mechanism that actually
// keeps torch/script.h out of every translation unit including this header, so
// the declaration was not doing the job its comment claimed either.

namespace gmd {

// ModelRuntimeAdapter backed by a TorchScript (*.pt) model exported by
// gmd_se3gnn's export.py.
//
// Expected model interface (set by the shared contract):
//   forward(species:     Tensor int64  [N]     -- atomic numbers
//           positions:   Tensor float32 [N,3]  -- Å
//           edge_index:  Tensor int64  [2,2E]  -- directed full graph
//           edge_shift:  Tensor float32 [2E,3] -- Cartesian shift (Å))
//   -> Dict{"energy": Tensor float64 scalar,
//            "forces": Tensor float32 [N,3]}
//
// The model must also expose an attribute `local_cutoff` (float / double)
// that ConfigLoader uses to set the VerletNeighborBuilder cutoff radius.
class TorchScriptModelRuntimeAdapter final : public ModelRuntimeAdapter {
public:
    TorchScriptModelRuntimeAdapter() noexcept;
    ~TorchScriptModelRuntimeAdapter() override;

    std::string_view name() const noexcept override;

    // Returns the local_cutoff attribute read from the model at load time.
    // Returns 0 if the model has not been loaded yet or the attribute is absent.
    float cutoff() const noexcept override;

    void load_model(const std::filesystem::path& model_path,
                    RuntimeContext& runtime) override;

    void evaluate(const ModelEvaluationRequest& request,
                  ModelEvaluationResult& result,
                  RuntimeContext& runtime) override;

    void unload_model(RuntimeContext& runtime) override;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace gmd

#endif  // GMD_ENABLE_TORCH
