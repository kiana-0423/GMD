#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include "gmd/force/bonded_params.hpp"
#include "gmd/force/force_provider.hpp"
#include "gmd/system/topology.hpp"

namespace gmd {

// ---------------------------------------------------------------------------
// Unit-conversion helpers
// ---------------------------------------------------------------------------
inline constexpr double kcal_per_mol_to_eV = 4.336410e-2;  // 1 kcal/mol in eV
inline constexpr double deg_to_rad         = 3.14159265358979323846 / 180.0;

// ---------------------------------------------------------------------------
// BondedForceProvider
//
// Implements the ForceProvider interface for intramolecular bonded terms:
//   - harmonic bonds
//   - harmonic angles
//   - periodic proper dihedrals
//   - harmonic improper dihedrals
//
// Usage:
//   auto bonded = std::make_shared<BondedForceProvider>();
//   bonded->set_topology(topo);
//   int bt = bonded->add_bond_type({k, r0});
//   // ... populate topology BondTerm entries with bt ...
//   composite.add(bonded);
// ---------------------------------------------------------------------------
class BondedForceProvider final : public ForceProvider {
public:
    BondedForceProvider() = default;
    explicit BondedForceProvider(std::shared_ptr<Topology> topology)
        : topology_(std::move(topology)) {}

    // Register a parameter set; returns 0-based type index.
    int add_bond_type(BondParams p);
    int add_angle_type(AngleParams p);
    int add_dihedral_type(DihedralParams p);
    int add_improper_type(ImproperParams p);

    // Replace the full topology (can be called before or after initialize).
    void set_topology(std::shared_ptr<Topology> topology);

    const Topology* topology() const noexcept { return topology_.get(); }

    // ForceProvider interface ------------------------------------------------
    std::string_view name() const noexcept override {
        return "bonded_force_provider";
    }

    void initialize(RuntimeContext& runtime) override;
    void compute(const ForceRequest& request,
                 ForceResult& result,
                 RuntimeContext& runtime) override;
    void finalize(RuntimeContext& runtime) override;

private:
    std::shared_ptr<Topology> topology_;

    std::vector<BondParams>     bond_types_;
    std::vector<AngleParams>    angle_types_;
    std::vector<DihedralParams> dihedral_types_;
    std::vector<ImproperParams> improper_types_;

    // Internal compute kernels (each accumulates into result).
    void compute_bonds    (const ForceRequest& req, ForceResult& result) const;
    void compute_angles   (const ForceRequest& req, ForceResult& result) const;
    void compute_dihedrals(const ForceRequest& req, ForceResult& result) const;
    void compute_impropers(const ForceRequest& req, ForceResult& result) const;
};

}  // namespace gmd
