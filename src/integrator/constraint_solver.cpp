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
    // Total SHAKE displacement applied to this atom, accumulated across the
    // iterations. The matching half-step velocity impulse is displacement/dt.
    System::Vec3 displacement{0.0, 0.0, 0.0};
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

// Builds the constraint virial from the converged RATTLE multipliers.
//
// See ConstraintVirialResult for the full derivation. In short: the RATTLE
// correction v_i += w_i Lambda_c r_c is the constraint share of the SECOND
// half-kick, the SHAKE impulse having already carried the first. Matching it
// against (dt/2m_i) G_i gives the ENDPOINT constraint force at t+dt,
// G_i = (2 Lambda_c / dt) r_c, and
//
//     W(t+dt) = sum_c (2 Lambda_c / dt) (r_c (x) r_c).
//
// This is an endpoint quantity, not a step average: it pairs with the provider
// virial at t+dt with no time-level mixing.
//
// r_c is fixed for the whole RATTLE solve, so summing the scalar multipliers
// gives an exactly central, exactly antisymmetric pair force and a symmetric
// tensor, with nothing projected away.
//
// MPI: this runs identically on every rank. load_atoms() allgathers all owned
// atoms when the communicator has more than one rank, and the constraint list is
// replicated against global atom tags, so every rank iterates the same
// constraints over the same coordinates and velocities and reaches the same
// multipliers. The tensor produced here is therefore ALREADY the global one and
// must not be reduced again -- an allreduce would multiply it by the rank count.
// There is no per-rank ownership split to make and no extra collective to enter,
// so ranks that own no atoms cannot diverge.
void accumulate_constraint_virial(const std::vector<AtomRecord>& atoms,
                                  const std::unordered_map<int, std::size_t>& index,
                                  const Box& box,
                                  const std::vector<double>& lambda_sum,
                                  double time_step,
                                  const std::vector<BondConstraint>& constraints,
                                  ConstraintVirialResult& out) {
    out.virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    out.contributions.clear();
    out.contributions.reserve(constraints.size());

    // Endpoint conversion: matching the RATTLE impulse against the second
    // half-kick gives (dt/2m_i) G_i = w_i Lambda_c r_c, hence G = (2/dt) Lambda r.
    const double force_scale = 2.0 / time_step;

    for (std::size_t c = 0; c < constraints.size(); ++c) {
        const auto& constraint = constraints[c];
        const auto found_i = index.find(constraint.i);
        const auto found_j = index.find(constraint.j);
        if (found_i == index.end() || found_j == index.end()) {
            throw std::runtime_error(
                "Constraint virial references an unknown atom tag");
        }
        const AtomRecord& atom_i = atoms[found_i->second];
        const AtomRecord& atom_j = atoms[found_j->second];

        // Bond vector at the geometry RATTLE ran on -- r_c(t+dt), which is the
        // geometry the force providers were evaluated at. Minimum image, so the result does
        // not depend on the origin or on how either atom happens to be wrapped.
        System::Vec3 bond{
            atom_i.coordinate[0] - atom_j.coordinate[0],
            atom_i.coordinate[1] - atom_j.coordinate[1],
            atom_i.coordinate[2] - atom_j.coordinate[2],
        };
        apply_minimum_image(bond, box);

        const double r2 = bond[0]*bond[0] + bond[1]*bond[1] + bond[2]*bond[2];
        if (!(r2 > 1.0e-20) || !std::isfinite(r2)) {
            throw std::runtime_error(
                "Constraint virial encountered a degenerate bond vector for atom pair (" +
                std::to_string(constraint.i) + ", " + std::to_string(constraint.j) + ")");
        }

        ConstraintContribution contribution;
        contribution.atom_tag_i = constraint.i;
        contribution.atom_tag_j = constraint.j;
        contribution.bond = bond;
        contribution.multiplier = lambda_sum[c];

        const double magnitude = force_scale * lambda_sum[c];
        for (std::size_t d = 0; d < 3; ++d) {
            contribution.force[d] = magnitude * bond[d];
            if (!std::isfinite(contribution.force[d])) {
                throw std::runtime_error(
                    "Constraint virial produced a non-finite force for atom pair (" +
                    std::to_string(constraint.i) + ", " + std::to_string(constraint.j) + ")");
            }
        }

        // W += r_c (x) G_c, exactly one contribution per physical constraint.
        // Symmetric by construction: G_c is parallel to r_c.
        for (std::size_t a = 0; a < 3; ++a) {
            for (std::size_t b = 0; b < 3; ++b) {
                out.virial[a * 3 + b] += bond[a] * contribution.force[b];
            }
        }
        out.contributions.push_back(contribution);
    }

    out.valid = true;
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

ConstraintReference ConstraintSolver::capture_reference(const System& system) const {
    ConstraintReference reference;
    if (constraints_.empty()) {
        return reference;
    }

    // Collective under MPI: load_atoms() allgathers, so every rank computes the
    // same reference bonds from the same coordinates.
    std::vector<AtomRecord> atoms;
    load_atoms(system, atoms);
    const auto index = make_tag_index(atoms);
    const Box& box = system.box();

    reference.bonds_.reserve(constraints_.size());
    for (const auto& constraint : constraints_) {
        const auto found_i = index.find(constraint.i);
        const auto found_j = index.find(constraint.j);
        if (found_i == index.end() || found_j == index.end()) {
            throw std::runtime_error(
                "SHAKE reference capture references an unknown atom tag");
        }
        const AtomRecord& atom_i = atoms[found_i->second];
        const AtomRecord& atom_j = atoms[found_j->second];
        System::Vec3 bond{
            atom_i.coordinate[0] - atom_j.coordinate[0],
            atom_i.coordinate[1] - atom_j.coordinate[1],
            atom_i.coordinate[2] - atom_j.coordinate[2],
        };
        apply_minimum_image(bond, box);
        const double r2 = bond[0]*bond[0] + bond[1]*bond[1] + bond[2]*bond[2];
        if (!(r2 > 1.0e-20) || !std::isfinite(r2)) {
            throw std::runtime_error(
                "SHAKE reference capture found a degenerate bond for atom pair (" +
                std::to_string(constraint.i) + ", " + std::to_string(constraint.j) + ")");
        }
        reference.bonds_.push_back(bond);
    }
    return reference;
}

ConstraintProjectionStats ConstraintSolver::apply_shake(System& system) const {
    return shake_impl(system, nullptr, 0.0);
}

ConstraintProjectionStats ConstraintSolver::apply_shake(
        System& system,
        const ConstraintReference& reference,
        double time_step) const {
    if (!(time_step > 0.0) || !std::isfinite(time_step)) {
        throw std::invalid_argument(
            "A dynamical SHAKE projection needs the strictly positive time step the "
            "drift used, since the matching velocity impulse is dr/dt; got " +
            std::to_string(time_step));
    }
    if (reference.size() != constraints_.size()) {
        throw std::invalid_argument(
            "SHAKE reference holds " + std::to_string(reference.size()) +
            " bonds but the solver has " + std::to_string(constraints_.size()) +
            " constraints; capture_reference() must be called on the same solver "
            "immediately before the position update");
    }
    return shake_impl(system, &reference, time_step);
}

// SHAKE.
//
// With a reference this is the textbook projection: the constraint gradient is
// evaluated at the configuration the step started from, so atoms are displaced
// along r_c(t) and the iteration solves for the amplitude. Writing
// s = r_ij(current) and r0 = r_c(t), one Gauss-Seidel sweep takes
//
//     g   = (|s|^2 - d^2) / (2 (w_i + w_j) (s . r0))
//     r_i -= w_i g r0        r_j += w_j g r0
//
// which linearises to |s|^2 - (|s|^2 - d^2) = d^2 as required. The correction is
// along r0 for BOTH atoms and mass-weighted oppositely, so it is a central pair
// impulse and conserves total momentum exactly. The accumulated displacement is
// then converted into the matching half-step velocity impulse, v += dr/dt.
//
// Without a reference this degenerates to displacing along the current bond,
// which is a plain geometric projection onto the manifold. That is correct for
// putting an initial or a rescaled state onto the manifold, and is NOT correct
// as a step of the dynamics -- see the ConstraintReference comment.
ConstraintProjectionStats ConstraintSolver::shake_impl(
        System& system,
        const ConstraintReference* reference,
        double time_step) const {
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
        for (std::size_t constraint_index = 0; constraint_index < constraints_.size();
             ++constraint_index) {
            const auto& constraint = constraints_[constraint_index];
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

            // Direction of the correction: the reference gradient when this is a
            // step of the dynamics, the current bond when it is a geometric
            // projection.
            const System::Vec3& direction =
                reference != nullptr ? reference->bond(constraint_index) : dr;
            const double projection =
                dr[0] * direction[0] + dr[1] * direction[1] + dr[2] * direction[2];
            if (!(std::abs(projection) > 1.0e-12 * std::sqrt(r2))) {
                // s . r0 -> 0 means the bond turned through ~90 degrees in one
                // step. SHAKE's linearisation has no solution there and silently
                // continuing would produce a wild correction.
                throw std::runtime_error(
                    "SHAKE cannot project atom pair (" + std::to_string(constraint.i) +
                    ", " + std::to_string(constraint.j) +
                    "): the bond rotated too far in one step for the reference "
                    "gradient to be usable (r(t+dt) . r(t) = " +
                    std::to_string(projection) + "). Reduce the time step.");
            }

            const double wi = 1.0 / atom_i.mass;
            const double wj = 1.0 / atom_j.mass;
            const double g = violation / (2.0 * (wi + wj) * projection);
            for (std::size_t dim = 0; dim < 3; ++dim) {
                const double correction = g * direction[dim];
                atom_i.coordinate[dim] -= wi * correction;
                atom_j.coordinate[dim] += wj * correction;
                atom_i.displacement[dim] -= wi * correction;
                atom_j.displacement[dim] += wj * correction;
            }
        }

        if (stats.max_error <= settings_.tolerance) {
            stats.converged = true;
            if (reference != nullptr) {
                // The half-step velocity impulse that belongs to the position
                // correction. Omitting it is what made the previous scheme
                // non-symplectic; see the ConstraintReference comment.
                const double inverse_dt = 1.0 / time_step;
                for (auto& atom : atoms) {
                    for (std::size_t dim = 0; dim < 3; ++dim) {
                        atom.velocity[dim] += atom.displacement[dim] * inverse_dt;
                    }
                }
            }
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
    return rattle_impl(system, 0.0, nullptr);
}

ConstraintProjectionStats ConstraintSolver::apply_rattle(
        System& system,
        double time_step,
        ConstraintVirialResult& virial_out) const {
    return rattle_impl(system, time_step, &virial_out);
}

ConstraintProjectionStats ConstraintSolver::rattle_impl(
        System& system,
        double time_step,
        ConstraintVirialResult* virial_out) const {
    auto stats = make_stats("RATTLE");
    if (virial_out != nullptr) {
        *virial_out = ConstraintVirialResult{};
        if (!(time_step > 0.0) || !std::isfinite(time_step)) {
            throw std::invalid_argument(
                "Recovering the constraint virial needs the strictly positive time step "
                "the step used; got " + std::to_string(time_step));
        }
    }
    if (constraints_.empty() || !settings_.enable_rattle) {
        // With RATTLE disabled there is no contemporaneous multiplier to report,
        // so the caller is left with an invalid result rather than a zero one.
        stats.enabled = !constraints_.empty();
        return stats;
    }

    std::vector<AtomRecord> atoms;
    load_atoms(system, atoms);
    auto index = make_tag_index(atoms);
    const Box& box = system.box();

    // Sum of each constraint's multipliers over the iterations. RATTLE does not
    // move atoms, so the bond vector each multiplier acts along is the same for
    // every iteration and a SCALAR sum is exact -- which is what makes the
    // recovered pair force central and antisymmetric by construction.
    std::vector<double> lambda_sum;
    if (virial_out != nullptr) {
        lambda_sum.assign(constraints_.size(), 0.0);
    }

    for (int iteration = 1; iteration <= settings_.max_iterations; ++iteration) {
        stats.iterations = iteration;
        stats.max_error = 0.0;
        for (std::size_t constraint_index = 0; constraint_index < constraints_.size();
             ++constraint_index) {
            const auto& constraint = constraints_[constraint_index];
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
            if (virial_out != nullptr) {
                lambda_sum[constraint_index] += lambda;
            }
        }

        if (stats.max_error <= settings_.tolerance) {
            stats.converged = true;
            if (virial_out != nullptr) {
                accumulate_constraint_virial(atoms, index, box, lambda_sum, time_step,
                                             constraints_, *virial_out);
            }
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
