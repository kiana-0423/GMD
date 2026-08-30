#include "gmd/system/initializer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "gmd/core/keyed_random.hpp"
#include "gmd/core/physical_constants.hpp"
#include "gmd/system/system.hpp"

#ifdef GMD_ENABLE_MPI
#include <mpi.h>
#endif

namespace gmd {

namespace {

// One authoritative definition, with its derivation and rounding policy,
// lives in gmd/core/physical_constants.hpp. This file previously carried
// its own literal, 1.13e-06 above the one the thermostats and the barostat
// used, so a system initialised to 300 K reported 300.000339 K.
constexpr double kBoltzmannConstant = kBoltzmannConstantEVPerKelvin;

#ifdef GMD_ENABLE_MPI
bool mpi_is_available() noexcept {
    int is_initialized = 0;
    int is_finalized = 0;
    MPI_Initialized(&is_initialized);
    MPI_Finalized(&is_finalized);
    return is_initialized != 0 && is_finalized == 0;
}
#endif

}  // namespace

VelocityInitializer::VelocityInitializer(std::uint32_t seed) noexcept
    : seed_(seed) {}

double VelocityInitializer::kinetic_energy(const System& system) const noexcept {
    const auto masses = system.masses();
    const auto velocities = system.velocities();

    double kinetic_energy = 0.0;
    // Owned atoms only. A ghost is another rank's atom borrowed for force
    // evaluation; adding its energy here would count that atom twice.
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        const auto& velocity = velocities[atom_index];
        const double velocity_squared = velocity[0] * velocity[0] +
                                        velocity[1] * velocity[1] +
                                        velocity[2] * velocity[2];
        kinetic_energy += masses[atom_index] * velocity_squared;
    }
    return kinetic_energy * 0.5;
}

void VelocityInitializer::sample_random_velocities(System& system, double target_temperature) const {
    auto masses = system.masses();
    auto velocities = system.mutable_velocities();

    // Owned atoms only. Ghosts are copies of atoms another rank owns and
    // receive their velocities through the normal communication path; sampling
    // them here would be harmless for the value -- the key is the tag, so a
    // ghost would draw its owner's velocity -- but it would double that atom's
    // contribution to the reductions below.
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        if (masses[atom_index] <= 0.0) {
            throw std::runtime_error("All particle masses must be positive before velocity initialization");
        }

        if (target_temperature == 0.0) {
            velocities[atom_index] = {0.0, 0.0, 0.0};
            continue;
        }

        // The atom's stable global tag is the random identity -- never its
        // index in this array, and never the rank. validate_atom_tags() has
        // already established that the tag is non-negative and globally
        // unique.
        const auto identity = static_cast<std::uint64_t>(system.atom_tag(atom_index));
        const double sigma =
            std::sqrt(kBoltzmannConstant * target_temperature / masses[atom_index]);
        velocities[atom_index] = {
            sigma * standard_normal(seed_, RandomStream::VelocityInitialization, identity, 0),
            sigma * standard_normal(seed_, RandomStream::VelocityInitialization, identity, 1),
            sigma * standard_normal(seed_, RandomStream::VelocityInitialization, identity, 2),
        };
    }
}

void VelocityInitializer::validate_atom_tags(const System& system) const {
    // The tag is the random draw's identity, so a duplicate would hand two
    // physical atoms the same velocity and a negative one has no defined
    // mapping. Both are rejected loudly rather than silently falling back to
    // an array index, which is the behaviour this replaced.
    // Negative tags are found first and reported COLLECTIVELY. Throwing here
    // the moment one is seen would leave the other ranks waiting in the gather
    // below, turning a bad input into a hang; every rank has to learn the
    // verdict and every rank has to throw.
    constexpr long long kNoNegative = std::numeric_limits<long long>::max();
    long long negative_tag = kNoNegative;
    std::vector<long long> local_tags;
    local_tags.reserve(system.num_local_atoms());
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        const int tag = system.atom_tag(atom_index);
        if (tag < 0) {
            negative_tag = std::min(negative_tag, static_cast<long long>(tag));
            continue;
        }
        local_tags.push_back(static_cast<long long>(tag));
    }

#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        long long global_negative = kNoNegative;
        MPI_Allreduce(&negative_tag, &global_negative, 1, MPI_LONG_LONG,
                      MPI_MIN, MPI_COMM_WORLD);
        negative_tag = global_negative;
    }
#endif

    if (negative_tag != kNoNegative) {
        throw std::runtime_error(
            "Atom tag " + std::to_string(negative_tag) + " is negative. Random velocity "
            "initialization keys each atom's draw on its global tag, which must "
            "be a non-negative identifier");
    }

    long long duplicate = -1;

#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        int rank = 0;
        int size = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        const int local_count = static_cast<int>(local_tags.size());
        std::vector<int> counts(static_cast<std::size_t>(size), 0);
        MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<int> displacements(static_cast<std::size_t>(size), 0);
        int total = 0;
        if (rank == 0) {
            for (int r = 0; r < size; ++r) {
                displacements[static_cast<std::size_t>(r)] = total;
                total += counts[static_cast<std::size_t>(r)];
            }
        }
        // Gathered on rank 0 rather than allgathered, so no rank has to hold a
        // buffer the size of the whole system just to check the tags.
        std::vector<long long> all_tags(rank == 0 ? static_cast<std::size_t>(total) : 0);
        MPI_Gatherv(local_tags.data(), local_count, MPI_LONG_LONG,
                    all_tags.data(), counts.data(), displacements.data(),
                    MPI_LONG_LONG, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::sort(all_tags.begin(), all_tags.end());
            for (std::size_t i = 1; i < all_tags.size(); ++i) {
                if (all_tags[i] == all_tags[i - 1]) {
                    duplicate = all_tags[i];
                    break;
                }
            }
        }
        // Every rank must learn the verdict and every rank must throw, or the
        // ranks that did not would run on into the next collective alone.
        MPI_Bcast(&duplicate, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    } else
#endif
    {
        std::vector<long long> sorted = local_tags;
        std::sort(sorted.begin(), sorted.end());
        for (std::size_t i = 1; i < sorted.size(); ++i) {
            if (sorted[i] == sorted[i - 1]) {
                duplicate = sorted[i];
                break;
            }
        }
    }

    if (duplicate >= 0) {
        throw std::runtime_error(
            "Global atom tag " + std::to_string(duplicate) + " appears more than once. "
            "Random velocity initialization keys each atom's draw on its global tag, so "
            "duplicated tags would give two physical atoms the same velocity");
    }
}

void VelocityInitializer::remove_center_of_mass_velocity(System& system) const {
    auto masses = system.masses();
    auto velocities = system.mutable_velocities();
    double total_mass = 0.0;
    System::Vec3 center_of_mass_velocity = {0.0, 0.0, 0.0};
    // Owned atoms only, for the same reason as in kinetic_energy(): a ghost
    // would contribute its owner's momentum a second time.
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        const auto mass = masses[atom_index];
        if (mass <= 0.0) {
            throw std::runtime_error("All particle masses must be positive before velocity initialization");
        }
        total_mass += mass;
        const auto& velocity = velocities[atom_index];
        center_of_mass_velocity[0] += masses[atom_index] * velocity[0];
        center_of_mass_velocity[1] += masses[atom_index] * velocity[1];
        center_of_mass_velocity[2] += masses[atom_index] * velocity[2];
    }

#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        double global_mass = 0.0;
        MPI_Allreduce(&total_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        total_mass = global_mass;

        double momentum[3] = {
            center_of_mass_velocity[0],
            center_of_mass_velocity[1],
            center_of_mass_velocity[2]
        };
        double global_momentum[3] = {0.0, 0.0, 0.0};
        MPI_Allreduce(momentum, global_momentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        center_of_mass_velocity[0] = global_momentum[0];
        center_of_mass_velocity[1] = global_momentum[1];
        center_of_mass_velocity[2] = global_momentum[2];
    }
#endif

    if (total_mass <= 0.0) {
        throw std::runtime_error("Total system mass must be positive before velocity initialization");
    }

    center_of_mass_velocity[0] /= total_mass;
    center_of_mass_velocity[1] /= total_mass;
    center_of_mass_velocity[2] /= total_mass;

    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        velocities[atom_index][0] -= center_of_mass_velocity[0];
        velocities[atom_index][1] -= center_of_mass_velocity[1];
        velocities[atom_index][2] -= center_of_mass_velocity[2];
    }
}

void VelocityInitializer::rescale_temperature(System& system,
                                              double target_temperature,
                                              bool center_of_mass_removed,
                                              std::size_t override_dof) const {
    if (target_temperature == 0.0) {
        auto velocities = system.mutable_velocities();
        for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
            velocities[atom_index] = {0.0, 0.0, 0.0};
        }
        return;
    }

    double current_kinetic_energy = kinetic_energy(system);
    double atom_count = static_cast<double>(system.num_local_atoms());

#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        double global_ke = 0.0;
        MPI_Allreduce(&current_kinetic_energy, &global_ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        current_kinetic_energy = global_ke;

        long long local_n = static_cast<long long>(system.num_local_atoms());
        long long global_n = 0;
        MPI_Allreduce(&local_n, &global_n, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        atom_count = static_cast<double>(global_n);
    }
#endif

    if (current_kinetic_energy <= 0.0) {
        throw std::runtime_error("Velocity initialization produced zero kinetic energy");
    }

    // The authoritative count when one was supplied. It is 3N - rank - 3, the
    // same number temperature reporting, the thermostats and the pressure use;
    // computing 3N-3 here instead is what made a constrained run start at the
    // wrong temperature. The fallback is only for callers with no constraints.
    const auto dof = override_dof != 0
                         ? static_cast<double>(override_dof)
                         : (center_of_mass_removed ? 3.0 * atom_count - 3.0
                                                   : 3.0 * atom_count);
    if (dof <= 0.0) {
        throw std::runtime_error("Not enough degrees of freedom to define a temperature");
    }

    const double current_temperature = current_kinetic_energy * 2.0 /
                                       (dof * kBoltzmannConstant);
    const double scale_factor = std::sqrt(target_temperature / current_temperature);

    auto velocities = system.mutable_velocities();
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        velocities[atom_index][0] *= scale_factor;
        velocities[atom_index][1] *= scale_factor;
        velocities[atom_index][2] *= scale_factor;
    }
}


double VelocityInitializer::global_twice_kinetic_energy(const System& system) const {
    double local = 2.0 * kinetic_energy(system);
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        double global = 0.0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return global;
    }
#endif
    return local;
}

namespace {

double global_total_mass(const System& system) {
    const auto masses = system.masses();
    double total = 0.0;
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        total += masses[atom_index];
    }
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        double global = 0.0;
        MPI_Allreduce(&total, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return global;
    }
#endif
    return total;
}

}  // namespace

double VelocityInitializer::global_momentum_magnitude(const System& system) const {
    const auto masses = system.masses();
    const auto velocities = system.velocities();
    double momentum[3] = {0.0, 0.0, 0.0};
    for (std::size_t atom_index = 0; atom_index < system.num_local_atoms(); ++atom_index) {
        for (std::size_t dim = 0; dim < 3; ++dim) {
            momentum[dim] += masses[atom_index] * velocities[atom_index][dim];
        }
    }
#ifdef GMD_ENABLE_MPI
    if (mpi_is_available()) {
        double global[3] = {0.0, 0.0, 0.0};
        MPI_Allreduce(momentum, global, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        momentum[0] = global[0];
        momentum[1] = global[1];
        momentum[2] = global[2];
    }
#endif
    return std::sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1] +
                     momentum[2] * momentum[2]);
}

void VelocityInitializer::initialize(System& system,
                                     double target_temperature,
                                     VelocityInitMode mode,
                                     bool remove_center_of_mass_velocity_flag) const {
    initialize(system, target_temperature, mode, remove_center_of_mass_velocity_flag,
               VelocityConstraintContext{});
}

// THE SEQUENCE, and why it is this one.
//
// Sampling, centre-of-mass removal, constraint projection and rescaling all
// touch the same velocities, so the order matters and not every order works.
// This one does, and it does so in a single pass:
//
//   1  sample tag-keyed Gaussians
//   2  remove the global centre-of-mass velocity, if requested
//   3  project into the constraint tangent space  (RATTLE)
//   4  rescale so that 2K = dof * k_B * T, with the AUTHORITATIVE dof
//
// The reason no iteration is needed is that each step preserves what the
// previous ones established:
//
//   * Centre-of-mass removal preserves tangency. A distance constraint's
//     Jacobian row is (+r_ij, -r_ij), so a uniform velocity shift contributes
//     r_ij.c - r_ij.c = 0: rigid translation lies in the null space of J.
//
//   * Projection preserves zero momentum. RATTLE's correction is
//     dv_i = +lambda r_ij / m_i, dv_j = -lambda r_ij / m_j, whose momentum
//     change is lambda r_ij - lambda r_ij = 0.
//
//   * Rescaling by a scalar preserves both: J(alpha v) = alpha J v = 0, and
//     alpha * 0 = 0.
//
// So after step 4 all three properties hold simultaneously and exactly, which
// is why the loop below normally runs once. It is still a loop because the
// projection is iterative and converges to a tolerance rather than exactly, and
// because a future constraint type might not have the null-space property that
// makes step 2 harmless. It verifies rather than assumes, and it reports a
// failure instead of accepting an inconsistent field.
void VelocityInitializer::initialize(System& system,
                                     double target_temperature,
                                     VelocityInitMode mode,
                                     bool remove_center_of_mass_velocity_flag,
                                     const VelocityConstraintContext& constraints) const {
    if (target_temperature < 0.0) {
        throw std::runtime_error("Target temperature must be non-negative");
    }

    // In MPI mode we must NOT return early when atom_count() == 0, because
    // other ranks may participate in collective MPI_Allreduce calls inside
    // remove_center_of_mass_velocity() and rescale_temperature(). Skipping
    // those calls on this rank would cause MPI_ERR_TRUNCATE.
    if (system.num_local_atoms() == 0) {
#ifdef GMD_ENABLE_MPI
        if (!mpi_is_available()) return;
        // Fall through — participate in allreduces with zero contributions.
#else
        return;
#endif
    }

    if (mode == VelocityInitMode::Random) {
        // Collective, and before any sampling: the tags are the draw's
        // identity, so they must be shown to be usable before they are used.
        validate_atom_tags(system);
        sample_random_velocities(system, target_temperature);
    }

    // Nothing below can make a field at rest anything other than at rest, and
    // the residual checks would divide by a zero target.
    if (target_temperature == 0.0) {
        rescale_temperature(system, target_temperature,
                            remove_center_of_mass_velocity_flag,
                            constraints.degrees_of_freedom);
        return;
    }

    // Relative tolerances on the three properties the field must satisfy at
    // once. The kinetic-energy target is set exactly by the final rescale, so
    // its bound only has to absorb the reduction round-off in a global sum; the
    // momentum bound is scaled by the momentum the field would have had before
    // removal, so it is dimensionless in the same way.
    constexpr double kKineticTolerance = 1.0e-12;
    constexpr double kMomentumTolerance = 1.0e-10;
    constexpr int kMaxPasses = 8;

    // The momentum residual is measured against the system's natural momentum
    // scale, sqrt(2K * M), and NOT against the net momentum the field happened
    // to start with. That distinction matters: velocities read from an input
    // file are often already very nearly momentum-free, so the starting net
    // momentum can be nine orders below the per-atom scale, and dividing a
    // round-off residual by it would report an enormous relative error for a
    // field that is in fact perfectly balanced.
    const double momentum_scale =
        remove_center_of_mass_velocity_flag
            ? std::sqrt(global_twice_kinetic_energy(system) * global_total_mass(system))
            : 0.0;

    int pass = 0;
    for (; pass < kMaxPasses; ++pass) {
        if (remove_center_of_mass_velocity_flag) {
            remove_center_of_mass_velocity(system);
        }
        if (constraints.has_constraints()) {
            constraints.project_velocities(system);
        }
        rescale_temperature(system, target_temperature,
                            remove_center_of_mass_velocity_flag,
                            constraints.degrees_of_freedom);

        // Verify rather than assume. Tangency is not re-measured here -- that
        // would mean reimplementing the constraint Jacobian, which is the
        // solver's job and is asserted by the tests instead; what is checked is
        // that the two GLOBAL scalar properties survived the projection.
        const double twice_ke = global_twice_kinetic_energy(system);
        const std::size_t dof = constraints.degrees_of_freedom;
        double kinetic_error = 0.0;
        if (dof != 0) {
            const double target =
                static_cast<double>(dof) * kBoltzmannConstant * target_temperature;
            kinetic_error = std::abs(twice_ke - target) / target;
        }
        double momentum_error = 0.0;
        if (remove_center_of_mass_velocity_flag && momentum_scale > 0.0) {
            momentum_error = global_momentum_magnitude(system) / momentum_scale;
        }
        if (kinetic_error <= kKineticTolerance && momentum_error <= kMomentumTolerance) {
            return;
        }
    }

    throw std::runtime_error(
        "Velocity initialization did not converge in " + std::to_string(kMaxPasses) +
        " passes: the centre-of-mass removal, the constraint projection and the "
        "temperature rescale could not be satisfied simultaneously. Rigid "
        "translation lies in the null space of a distance-constraint Jacobian and "
        "the projection conserves momentum, so one pass normally suffices; a "
        "failure here means a constraint type that breaks one of those "
        "properties, or a projection that is not converging");
}

}  // namespace gmd