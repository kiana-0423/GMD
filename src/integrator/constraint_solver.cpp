#include "gmd/integrator/constraint_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
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

// --- Independence analysis -------------------------------------------------

// Greedy rank-revealing selection with column pivoting.
//
// Modified Gram-Schmidt: repeatedly take the column whose residual after
// projecting out the basis so far is largest, and accept it if that residual
// clears `tolerance`. Returns the accepted column indices; everything else is
// redundant given them.
//
// TIE-BREAKING IS CANONICAL, NOT POSITIONAL. Symmetric molecules produce
// genuinely tied pivot candidates -- four equally spaced collinear atoms of
// equal mass, a regular ring -- and breaking those ties by column index would
// make the answer depend on the order the constraints happened to be listed in,
// which is exactly the property this is documented to have. Ties are settled by
// the canonical atom pair (min tag, max tag) instead, which is a property of the
// physical constraint.
//
// A tie is decided on a RELATIVE criterion. Two columns that are physically
// equivalent reach their residual norms by different arithmetic paths, so they
// agree only to within round-off accumulated over the projections; an absolute
// comparison would call them distinct and hand the choice back to input order.
// sqrt(eps) is the usual threshold for "these two computations of the same
// quantity are indistinguishable", and being generous here is safe: it only ever
// affects WHICH of two near-equivalent columns is named, and both are equally
// valid answers.
//
// This decides WHICH rows are redundant. HOW MANY is decided by the singular
// values, which is the more reliable question; the caller cross-checks the two
// and falls back to naming the whole component if they disagree.
std::vector<std::size_t> select_independent_columns(
        const std::vector<std::vector<double>>& columns,
        const std::vector<std::array<int, 2>>& canonical_keys,
        double tolerance) {
    const std::size_t n = columns.size();
    std::vector<std::size_t> chosen;
    if (n == 0) return chosen;

    std::vector<std::vector<double>> residual = columns;   // working copies
    std::vector<std::vector<double>> basis;                // orthonormal
    std::vector<bool> taken(n, false);

    const double tie_relative_tolerance =
        std::sqrt(std::numeric_limits<double>::epsilon());

    for (std::size_t step = 0; step < n; ++step) {
        std::vector<double> norms(n, -1.0);
        double best_norm = -1.0;
        for (std::size_t c = 0; c < n; ++c) {
            if (taken[c]) continue;
            double norm2 = 0.0;
            for (double v : residual[c]) norm2 += v * v;
            norms[c] = std::sqrt(norm2);
            best_norm = std::max(best_norm, norms[c]);
        }
        if (!(best_norm > tolerance)) break;

        // Everything within a relative tie of the leader is an equally good
        // pivot; take the one with the smallest canonical atom pair.
        const double cutoff = best_norm * (1.0 - tie_relative_tolerance);
        std::size_t best = n;
        for (std::size_t c = 0; c < n; ++c) {
            if (taken[c] || norms[c] < cutoff) continue;
            if (best == n || canonical_keys[c] < canonical_keys[best]) best = c;
        }
        if (best == n) break;
        best_norm = norms[best];

        taken[best] = true;
        chosen.push_back(best);
        std::vector<double> unit = residual[best];
        for (double& v : unit) v /= best_norm;
        // Re-orthogonalise once: MGS loses orthogonality on ill-conditioned
        // input, and this analysis is run precisely on such input.
        for (int pass = 0; pass < 2; ++pass) {
            for (const auto& b : basis) {
                double projection = 0.0;
                for (std::size_t k = 0; k < unit.size(); ++k) projection += unit[k] * b[k];
                for (std::size_t k = 0; k < unit.size(); ++k) unit[k] -= projection * b[k];
            }
            double norm2 = 0.0;
            for (double v : unit) norm2 += v * v;
            const double norm = std::sqrt(norm2);
            if (!(norm > 0.0)) break;
            for (double& v : unit) v /= norm;
        }
        basis.push_back(unit);

        for (std::size_t c = 0; c < n; ++c) {
            if (taken[c]) continue;
            double projection = 0.0;
            for (std::size_t k = 0; k < residual[c].size(); ++k) {
                projection += residual[c][k] * basis.back()[k];
            }
            for (std::size_t k = 0; k < residual[c].size(); ++k) {
                residual[c][k] -= projection * basis.back()[k];
            }
        }
    }
    std::sort(chosen.begin(), chosen.end());
    return chosen;
}

// Connected components of the constraint graph, as index lists.
std::vector<std::vector<std::size_t>> constraint_components(
        const std::vector<BondConstraint>& constraints) {
    // Union-find over atom tags.
    std::unordered_map<int, int> parent;
    std::function<int(int)> find = [&](int x) {
        auto it = parent.find(x);
        if (it == parent.end()) { parent.emplace(x, x); return x; }
        if (it->second == x) return x;
        const int root = find(it->second);
        it->second = root;
        return root;
    };
    auto unite = [&](int a, int b) {
        const int ra = find(a);
        const int rb = find(b);
        if (ra != rb) parent[ra] = rb;
    };
    for (const auto& constraint : constraints) {
        find(constraint.i);
        find(constraint.j);
        unite(constraint.i, constraint.j);
    }

    std::map<int, std::vector<std::size_t>> grouped;
    for (std::size_t index = 0; index < constraints.size(); ++index) {
        grouped[find(constraints[index].i)].push_back(index);
    }
    std::vector<std::vector<std::size_t>> components;
    components.reserve(grouped.size());
    for (auto& entry : grouped) components.push_back(std::move(entry.second));
    return components;
}

// std::to_string gives six fixed decimals, which prints every singular value and
// tolerance in this analysis as "0.000000". These numbers span many orders of
// magnitude and their exponent is the whole point, so they are formatted in
// scientific notation.
std::string format_number(double value) {
    std::ostringstream out;
    out << std::setprecision(6) << std::scientific << value;
    return out.str();
}

std::string format_tags(const std::vector<int>& tags) {
    std::string text = "{";
    for (std::size_t k = 0; k < tags.size(); ++k) {
        if (k != 0) text += ", ";
        if (k == 12) { text += "..."; break; }
        text += std::to_string(tags[k]);
    }
    return text + "}";
}

}  // namespace

ConstraintRankReport ConstraintSolver::analyze_independence(const System& system) const {
    ConstraintRankReport report;
    report.constraint_count = constraints_.size();
    report.independent = true;
    if (constraints_.empty()) return report;

    std::vector<AtomRecord> atoms;
    load_atoms(system, atoms);
    const auto index = make_tag_index(atoms);
    const Box& box = system.box();

    for (const auto& member_indices : constraint_components(constraints_)) {
        ConstraintComponent component;
        component.constraint_indices = member_indices;

        // Local atom numbering for this component.
        std::vector<int> tags;
        for (std::size_t c : member_indices) {
            tags.push_back(constraints_[c].i);
            tags.push_back(constraints_[c].j);
        }
        std::sort(tags.begin(), tags.end());
        tags.erase(std::unique(tags.begin(), tags.end()), tags.end());
        component.atom_tags = tags;
        std::unordered_map<int, std::size_t> local;
        for (std::size_t k = 0; k < tags.size(); ++k) local.emplace(tags[k], k);

        // Columns of J_M^T: one mass-weighted gradient per constraint.
        std::vector<std::vector<double>> columns;
        columns.reserve(member_indices.size());
        for (std::size_t c : member_indices) {
            const auto& constraint = constraints_[c];
            const auto found_i = index.find(constraint.i);
            const auto found_j = index.find(constraint.j);
            if (found_i == index.end() || found_j == index.end()) {
                throw std::runtime_error(
                    "Constraint independence analysis references an unknown atom tag");
            }
            const AtomRecord& atom_i = atoms[found_i->second];
            const AtomRecord& atom_j = atoms[found_j->second];
            if (!(atom_i.mass > 0.0) || !(atom_j.mass > 0.0)) {
                throw std::runtime_error(
                    "Constraint independence analysis requires strictly positive masses "
                    "for atom pair (" + std::to_string(constraint.i) + ", " +
                    std::to_string(constraint.j) + ")");
            }

            System::Vec3 bond{
                atom_i.coordinate[0] - atom_j.coordinate[0],
                atom_i.coordinate[1] - atom_j.coordinate[1],
                atom_i.coordinate[2] - atom_j.coordinate[2],
            };
            apply_minimum_image(bond, box);
            const double length = std::sqrt(bond[0]*bond[0] + bond[1]*bond[1] +
                                            bond[2]*bond[2]);

            // A gradient of 2 r_c vanishes with the bond. Scale the test against
            // the constraint's own target distance, which is the only length in
            // the problem, so it means the same thing in any unit system.
            if (!std::isfinite(length) ||
                length <= 1.0e-8 * constraint.target_distance) {
                report.problems.push_back(
                    "constraint " + std::to_string(c) + " on atom pair (" +
                    std::to_string(constraint.i) + ", " + std::to_string(constraint.j) +
                    ") in component " + format_tags(component.atom_tags) +
                    " has a degenerate gradient: the atoms are separated by " +
                    format_number(length) + " A against a target of " +
                    format_number(constraint.target_distance) +
                    ", so grad sigma = 2 r_ij is numerically zero and the constraint "
                    "cannot remove a degree of freedom");
                report.independent = false;
            }

            const double wi = 1.0 / std::sqrt(atom_i.mass);
            const double wj = 1.0 / std::sqrt(atom_j.mass);
            std::vector<double> column(3 * tags.size(), 0.0);
            const std::size_t li = local.at(constraint.i);
            const std::size_t lj = local.at(constraint.j);
            for (std::size_t d = 0; d < 3; ++d) {
                column[3 * li + d] = 2.0 * bond[d] * wi;
                column[3 * lj + d] = -2.0 * bond[d] * wj;
            }
            columns.push_back(std::move(column));
        }

        const std::size_t rows = 3 * tags.size();
        const std::size_t cols = member_indices.size();
        const auto values = jacobi_singular_values(
            columns, "component " + format_tags(component.atom_tags));
        component.largest_singular_value = values.empty() ? 0.0 : values.front();
        component.smallest_singular_value = values.empty() ? 0.0 : values.back();

        // SCALE-AWARE RANK TOLERANCE.
        //
        //     tol = max(rows, cols) * eps * sigma_max
        //
        // Every part of that is needed. Multiplying by sigma_max makes the test
        // RELATIVE: the entries of J_M are 2 r / sqrt(m), so they carry units of
        // length over root mass and their magnitude depends entirely on the unit
        // system and on how long the bonds are. An absolute threshold would mean
        // something different for a 1 A bond than for a 10 A one, and something
        // different again in another unit system. eps is the unit round-off, the
        // floor below which a computed singular value carries no information.
        // max(rows, cols) accounts for error accumulating over the matrix: each
        // singular value is the result of O(max(rows, cols)) arithmetic
        // operations, so the round-off floor grows with the larger dimension.
        // This is the standard LAPACK rank convention.
        component.rank_tolerance = static_cast<double>(std::max(rows, cols)) *
                                   std::numeric_limits<double>::epsilon() *
                                   component.largest_singular_value;
        component.rank = 0;
        for (double value : values) {
            if (value > component.rank_tolerance) ++component.rank;
        }
        component.condition_number =
            component.smallest_singular_value > component.rank_tolerance
                ? component.largest_singular_value / component.smallest_singular_value
                : std::numeric_limits<double>::infinity();

        report.rank += component.rank;
        if (!component.independent()) {
            report.independent = false;

            // Name the actual redundant rows when a pivoted rank-revealing
            // selection agrees with the singular values about how many there
            // are, and say so honestly when it does not. This is recorded PER
            // COMPONENT: with more than one deficient component, one may be
            // identified and another not, and a single report-level flag could
            // not describe that.
            std::vector<std::array<int, 2>> canonical_keys;
            canonical_keys.reserve(member_indices.size());
            for (std::size_t c : member_indices) {
                canonical_keys.push_back({std::min(constraints_[c].i, constraints_[c].j),
                                          std::max(constraints_[c].i, constraints_[c].j)});
            }
            const auto kept = select_independent_columns(columns, canonical_keys,
                                                         component.rank_tolerance);
            std::vector<std::size_t> redundant;
            for (std::size_t k = 0; k < member_indices.size(); ++k) {
                if (std::find(kept.begin(), kept.end(), k) == kept.end()) {
                    redundant.push_back(member_indices[k]);
                }
            }
            const bool identified = kept.size() == component.rank;
            component.redundant_constraints_identified = identified;
            component.redundant_constraint_indices = identified ? redundant : member_indices;

            // The theoretical ceiling, for a message that says what is wrong
            // rather than only that something is. A rigid body of K >= 3
            // non-collinear atoms has 6 rigid-body motions, so at most 3K - 6 of
            // its internal coordinates can be frozen; a pair has 5, giving 1.
            const std::size_t atom_count = tags.size();
            const std::size_t ceiling = atom_count >= 3
                ? (3 * atom_count >= 6 ? 3 * atom_count - 6 : 0)
                : (atom_count == 2 ? 1 : 0);
            std::string message =
                "component " + format_tags(component.atom_tags) + " has " +
                std::to_string(cols) + " constraint(s) over " +
                std::to_string(atom_count) + " atom(s) but the mass-weighted "
                "constraint Jacobian has rank " + std::to_string(component.rank) +
                ": " + std::to_string(cols - component.rank) +
                " of them are linearly dependent on the others and remove no degree "
                "of freedom. Singular values run from " +
                format_number(component.largest_singular_value) + " down to " +
                format_number(component.smallest_singular_value) +
                ", against a rank tolerance of " +
                format_number(component.rank_tolerance) + ".";
            if (cols > ceiling) {
                message += " A rigid body of " + std::to_string(atom_count) +
                           " atoms has at most " + std::to_string(ceiling) +
                           " independent distance constraints, so this component is "
                           "over-constrained by construction.";
            } else {
                message += " The count is within the rigid-body limit of " +
                           std::to_string(ceiling) +
                           ", so the dependence is geometric: at this configuration "
                           "the constraint gradients are coplanar or collinear.";
            }
            if (identified) {
                message += " Redundant given the others:";
                for (std::size_t c : redundant) {
                    message += " constraint " + std::to_string(c) + " on (" +
                               std::to_string(constraints_[c].i) + "," +
                               std::to_string(constraints_[c].j) + ")";
                }
                message += " (which member of a dependent group is named is not unique; "
                           "these are the ones a column-pivoted rank-revealing "
                           "factorisation leaves over, which is deterministic).";
            }
            message += " Constrained atom pairs:";
            for (std::size_t c : member_indices) {
                message += " (" + std::to_string(constraints_[c].i) + "," +
                           std::to_string(constraints_[c].j) + ")";
            }
            report.problems.push_back(message);
        } else if (component.condition_number >
                   1.0 / std::sqrt(std::numeric_limits<double>::epsilon())) {
            // Independent, but only just. Past this point roughly half the
            // significant digits of the multipliers are lost, so say so rather
            // than let the solver converge slowly to a poorly determined answer.
            report.warnings.push_back(
                "component " + format_tags(component.atom_tags) +
                " is independent but ill-conditioned: the mass-weighted constraint "
                "Jacobian has condition number " +
                format_number(component.condition_number) +
                ", so the constraint multipliers are poorly determined and SHAKE and "
                "RATTLE will converge slowly. The configuration is close to one where "
                "these constraints become dependent.");
        }

        report.components.push_back(std::move(component));
    }

    return report;
}

ConstraintRankReport ConstraintSolver::require_independent(
        const System& system, const std::string& context) const {
    ConstraintRankReport report = analyze_independence(system);
    if (report.independent) return report;

    std::string message =
        "Constraint set is not independent: " + std::to_string(report.constraint_count) +
        " constraint(s) supplied, but the mass-weighted constraint Jacobian has rank " +
        std::to_string(report.rank) + " at " + context + ", so " +
        std::to_string(report.redundant_count()) +
        " of them remove no degree of freedom. Counting them against the degrees of "
        "freedom would report a temperature that is too high, and the redundant rows "
        "leave the SHAKE and RATTLE multipliers undetermined. Remove the dependent "
        "constraints.";
    for (const auto& problem : report.problems) {
        message += "\n  - " + problem;
    }
    throw std::invalid_argument(message);
}

namespace {
// Rejects target distances that no non-degenerate configuration can satisfy.
//
// This is a property of the TARGETS ALONE, so it needs no coordinates and runs
// at construction, before any geometry exists. For three atoms all constrained
// to each other with targets a, b, c the triangle inequality decides everything:
//
//   c > a + b            no configuration satisfies all three at all;
//   c = a + b  (or the
//   reversed form
//   c = |a - b|)         the only configurations are COLLINEAR, where the three
//                        constraint gradients span two dimensions instead of
//                        three and the set is rank deficient.
//
// The second case is why this exists. The supplied coordinates can be a perfectly
// well conditioned triangle while the targets describe a collinear body, so a
// rank check on the input geometry passes and the projection then converges
// towards a configuration it cannot solve -- stalling without ever reaching a
// geometry degenerate enough for the rank check to catch. Catching it here says
// what is actually wrong, and says it before 500 wasted iterations.
void reject_degenerate_target_triangles(const std::vector<BondConstraint>& constraints,
                                        double tolerance) {
    std::unordered_map<std::uint64_t, double> target_of;
    std::unordered_map<int, std::vector<int>> neighbours;
    auto key_of = [](int a, int b) {
        const int lo = std::min(a, b);
        const int hi = std::max(a, b);
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32U) |
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(hi));
    };
    for (const auto& constraint : constraints) {
        target_of.emplace(key_of(constraint.i, constraint.j), constraint.target_distance);
        neighbours[constraint.i].push_back(constraint.j);
        neighbours[constraint.j].push_back(constraint.i);
    }

    for (auto& entry : neighbours) {
        const int j = entry.first;
        auto& adjacent = entry.second;
        std::sort(adjacent.begin(), adjacent.end());
        adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
        for (std::size_t p = 0; p + 1 < adjacent.size(); ++p) {
            for (std::size_t q = p + 1; q < adjacent.size(); ++q) {
                const int i = adjacent[p];
                const int k = adjacent[q];
                const auto closing = target_of.find(key_of(i, k));
                if (closing == target_of.end()) continue;   // not a constrained triple
                // Visit each triple once.
                if (!(j < i)) continue;

                const double a = target_of.at(key_of(i, j));
                const double b = target_of.at(key_of(j, k));
                const double c = closing->second;
                // Order the three so `longest` is the side the other two must
                // reach; only that arrangement can fail.
                double sides[3] = {a, b, c};
                std::sort(sides, sides + 3);
                const double shorter_sum = sides[0] + sides[1];
                const double longest = sides[2];
                const double slack = shorter_sum - longest;
                const std::string triple =
                    "(" + std::to_string(i) + ", " + std::to_string(j) + ", " +
                    std::to_string(k) + ") with target distances " +
                    format_number(a) + ", " + format_number(b) + " and " +
                    format_number(c);

                if (slack < -tolerance) {
                    throw std::invalid_argument(
                        "Infeasible constraint targets for the constrained triple " +
                        triple + ": " + format_number(longest) +
                        " exceeds the sum of the other two, " +
                        format_number(shorter_sum) +
                        ", so no configuration satisfies all three");
                }

                // How degenerate a nearly-collinear triple is must be measured as
                // a LENGTH, not as the triangle-inequality slack. The slack is a
                // second-order quantity: a triangle whose apex stands off the
                // base by a height h has slack ~ h^2 (a+b)/(2ab), so comparing it
                // against a distance tolerance would reject triangles that are
                // thin but perfectly resolvable. Inverting that relation gives
                // the height the targets imply,
                //
                //     h = sqrt(2 a b slack / (a + b))
                //
                // which is compared against the solver's own distance tolerance:
                // below it, SHAKE cannot tell the triangle from a straight line,
                // and the three gradients span two dimensions rather than three.
                const double product = sides[0] * sides[1];
                const double implied_height =
                    shorter_sum > 0.0
                        ? std::sqrt(2.0 * product * std::max(0.0, slack) / shorter_sum)
                        : 0.0;
                if (implied_height <= tolerance) {
                    throw std::invalid_argument(
                        "Degenerate constraint targets for the constrained triple " +
                        triple + ": " + format_number(longest) +
                        " reaches the sum of the other two, " +
                        format_number(shorter_sum) +
                        ", to within a triangle height of " +
                        format_number(implied_height) +
                        " A, at or below the constraint tolerance " +
                        format_number(tolerance) +
                        ". The only configurations satisfying these targets are "
                        "COLLINEAR to the solver's own resolution, and three collinear "
                        "distances have only two independent gradients, so the set is "
                        "rank deficient at every geometry it can reach -- whatever "
                        "geometry it starts from. Remove one of the three constraints");
                }
            }
        }
    }
}

ConstraintProjectionStats make_stats(const char* stage) {
    ConstraintProjectionStats stats;
    stats.enabled = true;
    stats.stage = stage;
    return stats;
}

}  // namespace

// Singular values of a small matrix given by its COLUMNS, by one-sided Jacobi.
//
// The columns here are the mass-weighted constraint gradients of one connected
// component, so the matrix is (3 * atoms) x (constraints) and the constraint
// count is the small dimension. One-sided Jacobi rotates pairs of columns until
// they are mutually orthogonal; the column norms are then the singular values.
//
// This is used rather than forming the Gram matrix G = J_M J_M^T and taking its
// eigenvalues. Squaring the matrix squares the condition number, which halves
// the number of correct digits in exactly the small singular values that decide
// the rank. One-sided Jacobi works on the matrix itself and is accurate to
// relative precision in each singular value, which is what a rank decision
// needs.
//
// CONVERGENCE IS CHECKED, NOT ASSUMED. After every sweep the scale-free
// off-orthogonality residual
//
//     max_{p<q} |a_p . a_q| / (|a_p| |a_q|)
//
// is measured, and the routine returns only once it falls below
//
//     cols * eps
//
// The threshold is dimension aware for a reason: a sweep rotates each column
// against every other, and orthogonalising a later pair perturbs an earlier one
// by O(eps) each time, so the residual an entire converged sweep can leave
// behind grows with the number of columns. eps alone would be unreachable for
// anything but a 2-column matrix; anything larger than O(cols * eps) would be
// slack. A pair involving a numerically zero column is orthogonal by definition
// and contributes nothing to the residual.
//
// Reaching the sweep limit without meeting it throws rather than returning a
// rank derived from an unconverged decomposition.
std::vector<double> jacobi_singular_values(std::vector<std::vector<double>> columns,
                                           const std::string& label) {
    const std::size_t n = columns.size();
    std::vector<double> values(n, 0.0);
    if (n == 0) return values;

    auto column_norms = [&columns](std::size_t count) {
        std::vector<double> norms(count, 0.0);
        for (std::size_t p = 0; p < count; ++p) {
            double norm2 = 0.0;
            for (double v : columns[p]) norm2 += v * v;
            norms[p] = std::sqrt(norm2);
        }
        return norms;
    };

    if (n == 1) {
        values[0] = column_norms(1).front();
        return values;
    }

    const std::size_t rows = columns.front().size();
    for (std::size_t p = 0; p < n; ++p) {
        if (columns[p].size() != rows) {
            throw std::runtime_error("Constraint Jacobian of " + label +
                                     " has columns of unequal length");
        }
        for (double v : columns[p]) {
            if (!std::isfinite(v)) {
                throw std::runtime_error(
                    "Constraint Jacobian of " + label + " contains a non-finite entry, so "
                    "its rank cannot be decided; check the coordinates and masses of the "
                    "atoms involved");
            }
        }
    }
    const double required_residual =
        static_cast<double>(n) * std::numeric_limits<double>::epsilon();
    // Jacobi converges quadratically once the columns are nearly orthogonal, so
    // this limit is far above what any well-posed matrix of this size needs. It
    // exists to bound the work, not to be reached.
    constexpr int max_sweeps = 60;

    double residual = std::numeric_limits<double>::infinity();
    int sweep = 0;
    for (; sweep < max_sweeps; ++sweep) {
        for (std::size_t p = 0; p + 1 < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                double alpha = 0.0;
                double beta = 0.0;
                double gamma = 0.0;
                for (std::size_t k = 0; k < rows; ++k) {
                    alpha += columns[p][k] * columns[p][k];
                    beta += columns[q][k] * columns[q][k];
                    gamma += columns[p][k] * columns[q][k];
                }
                if (alpha <= 0.0 || beta <= 0.0) continue;
                if (std::abs(gamma) <=
                    std::numeric_limits<double>::epsilon() * std::sqrt(alpha * beta)) {
                    continue;
                }

                const double zeta = (beta - alpha) / (2.0 * gamma);
                const double t = (zeta >= 0.0 ? 1.0 : -1.0) /
                                 (std::abs(zeta) + std::sqrt(1.0 + zeta * zeta));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double sn = c * t;
                for (std::size_t k = 0; k < rows; ++k) {
                    const double vp = columns[p][k];
                    const double vq = columns[q][k];
                    columns[p][k] = c * vp - sn * vq;
                    columns[q][k] = sn * vp + c * vq;
                }
            }
        }

        // Scale-free off-orthogonality of the whole matrix after this sweep.
        const auto norms = column_norms(n);
        residual = 0.0;
        for (std::size_t p = 0; p + 1 < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                if (!(norms[p] > 0.0) || !(norms[q] > 0.0)) continue;
                double gamma = 0.0;
                for (std::size_t k = 0; k < rows; ++k) gamma += columns[p][k] * columns[q][k];
                residual = std::max(residual, std::abs(gamma) / (norms[p] * norms[q]));
            }
        }
        if (residual <= required_residual) break;
    }

    if (!(residual <= required_residual)) {
        std::ostringstream message;
        message << std::setprecision(6) << std::scientific
                << "One-sided Jacobi failed to orthogonalise the constraint Jacobian of "
                << label << ": " << n << " column(s) of length " << rows
                << " after " << sweep << " sweep(s) still have off-orthogonality residual "
                << residual << ", above the required " << required_residual
                << ". The rank cannot be decided from an unconverged decomposition.";
        throw std::runtime_error(message.str());
    }

    values = column_norms(n);
    std::sort(values.begin(), values.end(), std::greater<double>());
    return values;
}

ConstraintSolver::ConstraintSolver(std::vector<BondConstraint> constraints,
                                   ConstraintSettings settings)
    : settings_(settings) {
    if (settings_.tolerance <= 0.0) {
        throw std::invalid_argument("Constraint tolerance must be positive");
    }
    if (settings_.max_iterations <= 0) {
        throw std::invalid_argument("Constraint max_iterations must be positive");
    }

    // Normalise the incoming list into distinct constraints: pairs are made
    // order-independent and duplicates removed. Independence of the resulting
    // set is not decided here -- it depends on a geometry this constructor does
    // not have -- and is verified by analyze_independence() before any dynamics
    // run. What IS decided here is target consistency, which needs no geometry:
    // conflicting targets for one pair, and target triples that no
    // non-degenerate configuration can satisfy.
    constraints_.reserve(constraints.size());
    std::unordered_map<std::uint64_t, std::size_t> seen;
    seen.reserve(constraints.size());
    // The orientation each pair was FIRST supplied in. "Reversed" is a property
    // of a repeat relative to that first occurrence, not of either entry alone:
    // (1, 0) after (0, 1) is reversed, and so is (0, 1) after (1, 0), while the
    // same orientation twice is not.
    std::unordered_map<std::uint64_t, int> first_orientation_i;

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
            first_orientation_i.emplace(key, incoming.i);
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
                std::to_string(hi) + "): target distances " + format_number(existing) +
                " and " + format_number(incoming.target_distance) + " differ by " +
                format_number(difference) +
                ", which exceeds the constraint tolerance " +
                format_number(settings_.tolerance));
        }

        if (incoming.i != first_orientation_i.at(key)) {
            ++diagnostics_.reversed_duplicates;
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

    // Targets that no non-degenerate configuration can satisfy are rejected here,
    // before any geometry exists to be projected.
    reject_degenerate_target_triangles(constraints_, settings_.tolerance);
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
