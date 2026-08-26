#include "gmd/integrator/constraint_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "gmd/system/minimum_image.hpp"
#include "gmd/system/periodic_boundary.hpp"
#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {
namespace {

struct AtomRecord {
    int tag = 0;
    double mass = 0.0;
    System::Vec3 coordinate{0.0, 0.0, 0.0};
    System::Vec3 velocity{0.0, 0.0, 0.0};
};

#ifdef GMD_ENABLE_MPI
bool mpi_is_available() noexcept {
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    return initialized != 0 && finalized == 0;
}

int mpi_size() noexcept {
    if (!mpi_is_available()) return 1;
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

void allgather_owned_atoms(const System& system,
                           std::vector<AtomRecord>& records) {
    constexpr int width = 8;
    std::vector<double> local;
    local.reserve(system.num_local_atoms() * width);
    const auto masses = system.masses();
    const auto coordinates = system.coordinates();
    const auto velocities = system.velocities();
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        local.push_back(static_cast<double>(system.atom_tag(atom_index)));
        local.push_back(masses[atom_index]);
        local.push_back(coordinates[atom_index][0]);
        local.push_back(coordinates[atom_index][1]);
        local.push_back(coordinates[atom_index][2]);
        local.push_back(velocities[atom_index][0]);
        local.push_back(velocities[atom_index][1]);
        local.push_back(velocities[atom_index][2]);
    }

    const int size = mpi_size();
    const int send_count = static_cast<int>(local.size());
    std::vector<int> recv_counts(static_cast<std::size_t>(size), 0);
    MPI_Allgather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> displs(static_cast<std::size_t>(size), 0);
    int total = 0;
    for (int rank = 0; rank < size; ++rank) {
        displs[static_cast<std::size_t>(rank)] = total;
        total += recv_counts[static_cast<std::size_t>(rank)];
    }
    if (total < 0 || total % width != 0) {
        throw std::runtime_error("Constraint MPI atom gather received malformed data");
    }

    std::vector<double> global(static_cast<std::size_t>(total), 0.0);
    MPI_Allgatherv(local.data(),
                   send_count,
                   MPI_DOUBLE,
                   global.data(),
                   recv_counts.data(),
                   displs.data(),
                   MPI_DOUBLE,
                   MPI_COMM_WORLD);

    records.clear();
    records.reserve(global.size() / width);
    for (std::size_t offset = 0; offset < global.size(); offset += width) {
        records.push_back(AtomRecord{
            .tag = static_cast<int>(global[offset]),
            .mass = global[offset + 1],
            .coordinate = {global[offset + 2], global[offset + 3], global[offset + 4]},
            .velocity = {global[offset + 5], global[offset + 6], global[offset + 7]},
        });
    }
}
#endif

void load_atoms(const System& system, std::vector<AtomRecord>& records) {
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available() && mpi_size() > 1) {
        allgather_owned_atoms(system, records);
        return;
    }
#endif
    records.clear();
    records.reserve(system.num_local_atoms());
    const auto masses = system.masses();
    const auto coordinates = system.coordinates();
    const auto velocities = system.velocities();
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        records.push_back(AtomRecord{
            .tag = system.atom_tag(atom_index),
            .mass = masses[atom_index],
            .coordinate = coordinates[atom_index],
            .velocity = velocities[atom_index],
        });
    }
}

std::unordered_map<int, std::size_t> make_tag_index(const std::vector<AtomRecord>& records) {
    std::unordered_map<int, std::size_t> index;
    index.reserve(records.size());
    for (std::size_t atom_index = 0; atom_index < records.size(); ++atom_index) {
        index.emplace(records[atom_index].tag, atom_index);
    }
    return index;
}

void write_local_atoms(System& system, const std::vector<AtomRecord>& records) {
    auto coordinates = system.mutable_coordinates();
    auto velocities = system.mutable_velocities();
    const auto index = make_tag_index(records);
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        const auto found = index.find(system.atom_tag(atom_index));
        if (found == index.end()) {
            throw std::runtime_error("Constraint projection lost a local atom tag");
        }
        coordinates[atom_index] = records[found->second].coordinate;
        wrap_position(coordinates[atom_index], system.box());
        velocities[atom_index] = records[found->second].velocity;
    }
    system.mutable_neighbor_list().valid = false;
}

ConstraintProjectionStats make_stats(const char* stage) {
    ConstraintProjectionStats stats;
    stats.enabled = true;
    stats.stage = stage;
    return stats;
}

}  // namespace

ConstraintSolver::ConstraintSolver(std::vector<BondConstraint> constraints,
                                   ConstraintSettings settings)
    : settings_(settings) {
    if (settings_.tolerance <= 0.0) {
        throw std::invalid_argument("Constraint tolerance must be positive");
    }
    if (settings_.max_iterations <= 0) {
        throw std::invalid_argument("Constraint max_iterations must be positive");
    }

    // Normalise the incoming list into distinct constraints. See the class
    // comment for what this does and does not guarantee: pairs are made
    // order-independent and duplicates removed, but independence of the
    // resulting set is assumed, not proven.
    constraints_.reserve(constraints.size());
    std::unordered_map<std::uint64_t, std::size_t> seen;
    seen.reserve(constraints.size());

    for (const auto& incoming : constraints) {
        if (incoming.i < 0 || incoming.j < 0) {
            throw std::invalid_argument(
                "BondConstraint requires non-negative atom tags, got (" +
                std::to_string(incoming.i) + ", " + std::to_string(incoming.j) + ")");
        }
        if (incoming.i == incoming.j) {
            throw std::invalid_argument(
                "BondConstraint cannot constrain atom " + std::to_string(incoming.i) +
                " to itself");
        }
        if (!std::isfinite(incoming.target_distance) || incoming.target_distance <= 0.0) {
            throw std::invalid_argument(
                "BondConstraint requires a finite, strictly positive target distance, got " +
                std::to_string(incoming.target_distance));
        }

        // (i, j) and (j, i) denote the same constraint.
        const int lo = std::min(incoming.i, incoming.j);
        const int hi = std::max(incoming.i, incoming.j);
        const auto key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32U) |
                          static_cast<std::uint64_t>(static_cast<std::uint32_t>(hi));

        const auto found = seen.find(key);
        if (found == seen.end()) {
            seen.emplace(key, constraints_.size());
            constraints_.push_back(BondConstraint{lo, hi, incoming.target_distance});
            continue;
        }

        // The pair is already present. Classify the repeat by how far its
        // target sits from the one already kept.
        const double existing = constraints_[found->second].target_distance;
        const double difference = std::abs(existing - incoming.target_distance);

        if (difference > settings_.tolerance) {
            // No geometry satisfies both; silently keeping one would hide a
            // bad input.
            throw std::invalid_argument(
                "Conflicting constraints for atom pair (" + std::to_string(lo) + ", " +
                std::to_string(hi) + "): target distances " + std::to_string(existing) +
                " and " + std::to_string(incoming.target_distance) + " differ by " +
                std::to_string(difference) +
                ", which exceeds the constraint tolerance " +
                std::to_string(settings_.tolerance));
        }

        if (difference == 0.0) {
            ++diagnostics_.exact_duplicates;
            continue;
        }

        // Tolerance-equivalent: SHAKE converges to within `tolerance`, so the
        // two targets are not distinguishable by the solver. Keep the first
        // value -- which makes the outcome deterministic and independent of the
        // orientation each pair was supplied in -- and record what was dropped
        // so the caller can report it.
        ++diagnostics_.tolerance_equivalent_duplicates;
        diagnostics_.discarded_targets.push_back(
            "atom pair (" + std::to_string(lo) + ", " + std::to_string(hi) +
            "): kept target distance " + std::to_string(existing) +
            ", discarded " + std::to_string(incoming.target_distance) +
            " (differ by " + std::to_string(difference) +
            ", within constraint tolerance " + std::to_string(settings_.tolerance) + ")");
    }
}

ConstraintProjectionStats ConstraintSolver::apply_shake(System& system) const {
    auto stats = make_stats("SHAKE");
    if (constraints_.empty()) {
        stats.enabled = false;
        return stats;
    }

    std::vector<AtomRecord> atoms;
    load_atoms(system, atoms);
    auto index = make_tag_index(atoms);
    const Box& box = system.box();

    for (int iteration = 1; iteration <= settings_.max_iterations; ++iteration) {
        stats.iterations = iteration;
        stats.max_error = 0.0;
        for (const auto& constraint : constraints_) {
            const auto found_i = index.find(constraint.i);
            const auto found_j = index.find(constraint.j);
            if (found_i == index.end() || found_j == index.end()) {
                throw std::runtime_error("SHAKE constraint references an unknown atom tag");
            }

            auto& atom_i = atoms[found_i->second];
            auto& atom_j = atoms[found_j->second];
            System::Vec3 dr{
                atom_i.coordinate[0] - atom_j.coordinate[0],
                atom_i.coordinate[1] - atom_j.coordinate[1],
                atom_i.coordinate[2] - atom_j.coordinate[2],
            };
            apply_minimum_image(dr, box);

            const double r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if (r2 < 1.0e-20) {
                throw std::runtime_error("SHAKE constraint encountered coincident atoms");
            }

            const double target2 = constraint.target_distance * constraint.target_distance;
            const double violation = r2 - target2;
            stats.max_error = std::max(
                stats.max_error,
                std::abs(std::sqrt(r2) - constraint.target_distance));

            if (std::abs(std::sqrt(r2) - constraint.target_distance) <= settings_.tolerance) {
                continue;
            }

            const double wi = 1.0 / atom_i.mass;
            const double wj = 1.0 / atom_j.mass;
            const double lambda = -violation / (2.0 * (wi + wj) * r2);
            for (std::size_t dim = 0; dim < 3; ++dim) {
                const double correction = lambda * dr[dim];
                atom_i.coordinate[dim] += wi * correction;
                atom_j.coordinate[dim] -= wj * correction;
            }
        }

        if (stats.max_error <= settings_.tolerance) {
            stats.converged = true;
            write_local_atoms(system, atoms);
            return stats;
        }
    }

    stats.converged = false;
    write_local_atoms(system, atoms);
    throw std::runtime_error(
        "SHAKE failed to converge within " + std::to_string(settings_.max_iterations) +
        " iterations; max bond error = " + std::to_string(stats.max_error));
}

ConstraintProjectionStats ConstraintSolver::apply_rattle(System& system) const {
    auto stats = make_stats("RATTLE");
    if (constraints_.empty() || !settings_.enable_rattle) {
        stats.enabled = !constraints_.empty();
        return stats;
    }

    std::vector<AtomRecord> atoms;
    load_atoms(system, atoms);
    auto index = make_tag_index(atoms);
    const Box& box = system.box();

    for (int iteration = 1; iteration <= settings_.max_iterations; ++iteration) {
        stats.iterations = iteration;
        stats.max_error = 0.0;
        for (const auto& constraint : constraints_) {
            const auto found_i = index.find(constraint.i);
            const auto found_j = index.find(constraint.j);
            if (found_i == index.end() || found_j == index.end()) {
                throw std::runtime_error("RATTLE constraint references an unknown atom tag");
            }

            auto& atom_i = atoms[found_i->second];
            auto& atom_j = atoms[found_j->second];
            System::Vec3 dr{
                atom_i.coordinate[0] - atom_j.coordinate[0],
                atom_i.coordinate[1] - atom_j.coordinate[1],
                atom_i.coordinate[2] - atom_j.coordinate[2],
            };
            apply_minimum_image(dr, box);
            const double r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if (r2 < 1.0e-20) {
                throw std::runtime_error("RATTLE constraint encountered coincident atoms");
            }

            const System::Vec3 dv{
                atom_i.velocity[0] - atom_j.velocity[0],
                atom_i.velocity[1] - atom_j.velocity[1],
                atom_i.velocity[2] - atom_j.velocity[2],
            };
            const double dot = dr[0] * dv[0] + dr[1] * dv[1] + dr[2] * dv[2];
            const double velocity_error = std::abs(dot) / std::sqrt(r2);
            stats.max_error = std::max(stats.max_error, velocity_error);
            if (velocity_error <= settings_.tolerance) {
                continue;
            }

            const double wi = 1.0 / atom_i.mass;
            const double wj = 1.0 / atom_j.mass;
            const double lambda = -dot / ((wi + wj) * r2);
            for (std::size_t dim = 0; dim < 3; ++dim) {
                const double correction = lambda * dr[dim];
                atom_i.velocity[dim] += wi * correction;
                atom_j.velocity[dim] -= wj * correction;
            }
        }

        if (stats.max_error <= settings_.tolerance) {
            stats.converged = true;
            write_local_atoms(system, atoms);
            return stats;
        }
    }

    stats.converged = false;
    write_local_atoms(system, atoms);
    throw std::runtime_error(
        "RATTLE failed to converge within " + std::to_string(settings_.max_iterations) +
        " iterations; max velocity constraint error = " + std::to_string(stats.max_error));
}

std::vector<BondConstraint> constraints_from_bond_types(
    const Topology& topology,
    const std::vector<int>& constrained_bond_types,
    const std::vector<double>& bond_type_distances) {
    std::vector<BondConstraint> constraints = topology.constraints;
    for (const auto& bond : topology.bonds) {
        if (std::find(constrained_bond_types.begin(),
                      constrained_bond_types.end(),
                      bond.type_idx) == constrained_bond_types.end()) {
            continue;
        }
        if (bond.type_idx < 0 ||
            static_cast<std::size_t>(bond.type_idx) >= bond_type_distances.size()) {
            throw std::runtime_error("Constraint bond type has no target distance");
        }
        constraints.push_back(BondConstraint{
            .i = bond.i,
            .j = bond.j,
            .target_distance = bond_type_distances[static_cast<std::size_t>(bond.type_idx)],
        });
    }
    return constraints;
}

}  // namespace gmd
