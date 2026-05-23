#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "gmd/force/force_provider.hpp"
#include "gmd/system/minimum_image.hpp"
#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

// Applies the direct-space correction that converts an unmodified
// Ewald/PME reciprocal calculation into scaled topology special pairs.
// MPI gathers positions for only locally owned atoms and applies force only
// to owned endpoints, so correctness does not depend on halo width.
inline void apply_special_pair_coulomb_corrections(const ForceRequest& req,
                                                    ForceResult& res,
                                                    double coulomb_constant) {
    if (req.system == nullptr || req.box == nullptr ||
        req.system->special_pair_map() == nullptr ||
        req.system->special_pair_map()->empty()) {
        return;
    }

    const auto& entries = req.system->special_pair_map()->entries();
    int maximum_tag = -1;
    for (const auto& pair : entries) {
        maximum_tag = std::max(maximum_tag, pair.atom_tag_b);
    }
    if (maximum_tag < 0) return;

    std::vector<Coordinate3D> positions(static_cast<std::size_t>(maximum_tag + 1));
    std::vector<double> tag_charges(static_cast<std::size_t>(maximum_tag + 1), 0.0);
    std::vector<bool> position_present(static_cast<std::size_t>(maximum_tag + 1), false);
    std::vector<int> local_index(static_cast<std::size_t>(maximum_tag + 1), -1);
    std::vector<double> owned_records;
    const auto coordinates = req.coordinates;
    const auto charges = req.system->charges();
    for (std::size_t i = 0; i < req.system->num_local_atoms(); ++i) {
        const int tag = req.system->atom_tag(i);
        if (tag < 0 || tag > maximum_tag) continue;
        local_index[static_cast<std::size_t>(tag)] = static_cast<int>(i);
        owned_records.push_back(static_cast<double>(tag));
        owned_records.push_back(coordinates[i][0]);
        owned_records.push_back(coordinates[i][1]);
        owned_records.push_back(coordinates[i][2]);
        owned_records.push_back(charges[i]);
    }

    std::vector<double> all_records = owned_records;
#ifdef GMD_ENABLE_MPI
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized != 0 && finalized == 0) {
        const int send_count = static_cast<int>(owned_records.size());
        int size = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        std::vector<int> counts(static_cast<std::size_t>(size), 0);
        MPI_Allgather(&send_count, 1, MPI_INT,
                      counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        std::vector<int> displacements(static_cast<std::size_t>(size), 0);
        int total = 0;
        for (int rank = 0; rank < size; ++rank) {
            displacements[static_cast<std::size_t>(rank)] = total;
            total += counts[static_cast<std::size_t>(rank)];
        }
        all_records.resize(static_cast<std::size_t>(total));
        MPI_Allgatherv(owned_records.data(), send_count, MPI_DOUBLE,
                       all_records.data(), counts.data(), displacements.data(),
                       MPI_DOUBLE, MPI_COMM_WORLD);
    }
#endif
    for (std::size_t offset = 0; offset + 4 < all_records.size(); offset += 5) {
        const int tag = static_cast<int>(all_records[offset]);
        if (tag < 0 || tag > maximum_tag) continue;
        positions[static_cast<std::size_t>(tag)] = {
            all_records[offset + 1], all_records[offset + 2], all_records[offset + 3]};
        tag_charges[static_cast<std::size_t>(tag)] = all_records[offset + 4];
        position_present[static_cast<std::size_t>(tag)] = true;
    }

    for (const auto& pair : entries) {
        const double delta = pair.scale.coulomb - 1.0;
        if (delta == 0.0 ||
            !position_present[static_cast<std::size_t>(pair.atom_tag_a)] ||
            !position_present[static_cast<std::size_t>(pair.atom_tag_b)]) {
            continue;
        }
        const int index_a = local_index[static_cast<std::size_t>(pair.atom_tag_a)];
        const int index_b = local_index[static_cast<std::size_t>(pair.atom_tag_b)];
        if (index_a < 0 && index_b < 0) continue;

        const double qa = tag_charges[static_cast<std::size_t>(pair.atom_tag_a)];
        const double qb = tag_charges[static_cast<std::size_t>(pair.atom_tag_b)];
        if (qa == 0.0 || qb == 0.0) continue;

        Force3D dr = {
            positions[static_cast<std::size_t>(pair.atom_tag_a)][0] -
                positions[static_cast<std::size_t>(pair.atom_tag_b)][0],
            positions[static_cast<std::size_t>(pair.atom_tag_a)][1] -
                positions[static_cast<std::size_t>(pair.atom_tag_b)][1],
            positions[static_cast<std::size_t>(pair.atom_tag_a)][2] -
                positions[static_cast<std::size_t>(pair.atom_tag_b)][2],
        };
        apply_minimum_image(dr, *req.box);
        const double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        if (r2 < 1.0e-12) continue;
        const double r = std::sqrt(r2);
        const double correction = delta * coulomb_constant * qa * qb;
        if (index_a >= 0) {
            res.potential_energy += correction / r;
            const double ff = correction / (r2 * r);
            for (std::size_t d = 0; d < 3; ++d) {
                res.forces[static_cast<std::size_t>(index_a)][d] += ff * dr[d];
            }
        }
        if (index_b >= 0) {
            const double ff = correction / (r2 * r);
            for (std::size_t d = 0; d < 3; ++d) {
                res.forces[static_cast<std::size_t>(index_b)][d] -= ff * dr[d];
            }
        }
    }
}

}  // namespace gmd
