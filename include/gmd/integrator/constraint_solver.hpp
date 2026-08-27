#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "gmd/system/topology.hpp"

namespace gmd {

class System;

struct ConstraintSettings {
    double tolerance = 1.0e-6;
    int max_iterations = 100;
    bool enable_rattle = true;
};

struct ConstraintProjectionStats {
    bool enabled = false;
    bool converged = true;
    int iterations = 0;
    double max_error = 0.0;
    std::string stage;
};

// Whether the reported virial has a constraint contribution attached to it.
//
//   NotApplicable  the run has no active constraints, so the force-provider
//                  virial is complete as it stands.
//   Unavailable    constraints are active but no constraint multiplier belongs
//                  to the currently stored provider virial. The reported virial
//                  is then INCOMPLETE, and is marked invalid rather than being
//                  passed off as a complete pressure virial.
//   Valid          a contemporaneous constraint virial is stored and included.
enum class ConstraintVirialState { NotApplicable, Unavailable, Valid };

// One constrained pair's contribution to the constraint virial, kept so the
// total can be audited per constraint.
struct ConstraintContribution {
    int atom_tag_i = 0;
    int atom_tag_j = 0;
    // r_c(t+dt) = r_i - r_j at the geometry RATTLE ran on, minimum image [A].
    std::array<double, 3> bond{};
    // ENDPOINT constraint force on atom i at t+dt: G_i = (2 Lambda / dt) r_c.
    // Not a step average. G_j = -G_i.
    std::array<double, 3> force{};
    double multiplier = 0.0;             // Lambda: RATTLE multipliers, summed
};

// The constraint contribution to the configurational virial at time level
// t+dt, recovered from the converged RATTLE multipliers.
//
// ---------------------------------------------------------------------------
// 1. THE ALGORITHM THIS IS DERIVED FROM
// ---------------------------------------------------------------------------
// VelocityVerletIntegrator implements the standard constrained velocity-Verlet
// (SHAKE/RATTLE) splitting. Over one step, with w = 1/m, r_c(t) the bond vector
// of constraint c BEFORE the position update, r_c(t+dt) the one after it, and
// s_ic = +1 when atom i is the first atom of c and -1 when it is the second:
//
//   (1)  v_i += (dt / 2 m_i) F_i(t)                          first half-kick
//   (2)  v_i += w_i sum_c s_ic Gamma_c r_c(t)                SHAKE impulse
//   (3)  r_i += dt v_i                                       drift
//   (4)  v_i += (dt / 2 m_i) F_i(t+dt)                       second half-kick
//   (5)  v_i += w_i sum_c s_ic Lambda_c r_c(t+dt)            RATTLE impulse
//
// Steps (2) and (3) are what apply_shake() performs together: SHAKE solves for
// the multipliers Gamma_c that put r(t+dt) on the constraint manifold, and the
// same impulse is applied to the velocity as well, which is the `v += dr/dt`
// correction of the textbook scheme. RATTLE then solves for the Lambda_c that
// make the velocities tangent at r(t+dt).
//
// TWO DISTINCT CONSTRAINT HALF-IMPULSES. (2) carries the constraint force at
// time level t, paired with F(t); (5) carries it at t+dt, paired with F(t+dt).
// They are separate quantities and are NOT summed here.
//
// ---------------------------------------------------------------------------
// 2. WHICH ONE THE REPORTED PRESSURE NEEDS, AND THE CONVERSION FACTOR
// ---------------------------------------------------------------------------
// The force providers are evaluated at r(t+dt), so the provider virial belongs
// to time level t+dt. Its constraint partner is therefore the RATTLE impulse of
// step (5), not the SHAKE impulse of step (2). Matching (5) against the second
// half-kick (4) term by term,
//
//     (dt / 2 m_i) G_i(t+dt) = w_i Lambda_c r_c(t+dt)
//     =>  G_i(t+dt) = (2 / dt) sum_c s_ic Lambda_c r_c(t+dt)
//
// THE FACTOR IS 2/dt, and it is exact: it is not a step average and involves no
// approximation. Both G and the provider force are ENDPOINT quantities at t+dt,
// so the combined virial is a single-time-level estimator with no O(dt) mixing.
//
//     W_constraint(t+dt) = sum_c (2 Lambda_c / dt) (r_c(t+dt) (x) r_c(t+dt))
//
// The factor was confirmed against an independent physical result rather than
// assumed: for a freely rotating rigid dimer the endpoint constraint force must
// equal the centripetal force mu omega^2 d, and it does, to O(dt^2). Halving it
// to 1/dt -- the value that would be right if the SHAKE velocity impulse (2)
// were omitted, so that RATTLE had to carry the whole step -- undershoots it by
// exactly two. See the rigid-rotor cases in tests/constraint_virial_tests.cpp.
//
// Units need no conversion. G is an impulse divided by a time, both read
// directly off the integrator's own update, so it carries the same force units
// the providers use and W is directly additive to their virial.
//
// ---------------------------------------------------------------------------
// 3. CENTRAL AND SYMMETRIC BY CONSTRUCTION
// ---------------------------------------------------------------------------
// RATTLE does not move atoms, so r_c(t+dt) is fixed for the entire solve.
// Summing the SCALAR multipliers therefore gives a pair force exactly along
// r_c and exactly antisymmetric between the two atoms -- for coupled
// constraints such as a rigid triangle just as much as for isolated ones.
// Nothing is projected onto anything and no residual is discarded.
//
// In the project's convention W = sum_i r_i (x) F_i, using G_j = -G_i, the
// tensor above follows, and it is manifestly symmetric. Being built from bond
// vectors alone it is independent of the coordinate origin and of how the atoms
// happen to be wrapped.
//
// ---------------------------------------------------------------------------
// 4. SIGN (worked, not fitted)
// ---------------------------------------------------------------------------
// Each RATTLE iteration takes lambda = -(r_c . v_c) / ((w_i + w_j) |r_c|^2). A
// separating pair has r_c . v_c > 0, hence lambda < 0, hence G_i pointing from i
// towards j: an attractive restoring force and a negative contribution to tr W.
// That matches the negative virial an attractive Lennard-Jones pair produces
// under the same convention, and it matches the rotor, whose constraint force
// must be centripetal.
//
// ---------------------------------------------------------------------------
// 5. WHY THE POST-RATTLE KINETIC TERM DOES NOT ALREADY COVER THIS
// ---------------------------------------------------------------------------
// The estimator is P = (2K + tr W) / 3V with 2K = sum m v^2 taken from the
// post-RATTLE velocities. The kinetic term and the configurational constraint
// virial are different quantities; projecting the velocities does not remove the
// need for the second. The rigid rotor makes it concrete. A dimer of reduced
// mass mu rotating at omega with bond length d and no centre-of-mass motion has
// relative speed omega d, so
//
//     2K = mu omega^2 d^2,        tr W_constraint = -mu omega^2 d^2,
//
// and 2K + tr W = 0. That is the physically required answer: the bond is frozen,
// so the rigid pair contributes nothing to the pressure. Omitting the constraint
// virial leaves 2K uncancelled and reports a spurious pressure of 2K / 3V.
// tests/constraint_virial_tests.cpp asserts this cancellation directly, and
// shows the residual falling as dt^2.
struct ConstraintVirialResult {
    bool valid = false;
    std::array<double, 9> virial{};      // row-major, W_ab = virial[a*3 + b]
    std::vector<ConstraintContribution> contributions;

    double trace() const noexcept { return virial[0] + virial[4] + virial[8]; }
};

// The constraint geometry as it stood BEFORE a position update, which is what a
// standard SHAKE projection corrects along.
//
// WHY THIS EXISTS. sigma_c = |r_c|^2 - d_c^2 has gradient 2 r_c, evaluated at
// the configuration the step started from. Textbook SHAKE therefore displaces
// atoms along r_c(t) and solves for the amplitude; it does NOT displace them
// along the drifted bond r_c'(t+dt). Both land on the manifold, but at different
// points, and only the first one reproduces the constrained equations of motion.
// Correcting along the drifted bond is a pure rescaling of the bond, which is
// not symplectic and drains energy secularly: a free rigid rotor loses ~86% of
// its kinetic energy over 4000 steps at omega*dt = 0.04, where the reference
// form conserves it to 3e-13.
//
// Populated by ConstraintSolver::capture_reference() before the drift, and it
// is a globally replicated quantity: the gather happens inside it, so every rank
// holds the same reference bonds.
class ConstraintReference {
public:
    bool empty() const noexcept { return bonds_.empty(); }
    std::size_t size() const noexcept { return bonds_.size(); }
    const std::array<double, 3>& bond(std::size_t constraint_index) const {
        return bonds_.at(constraint_index);
    }

private:
    friend class ConstraintSolver;
    std::vector<std::array<double, 3>> bonds_;
};

// Diagnostics produced while normalising a constraint list.
//
// A repeated atom pair falls into exactly one of three categories, decided by
// how far apart the two target distances are:
//
//   |d1 - d2| == 0                    exact duplicate
//   0 < |d1 - d2| <= tolerance        tolerance-equivalent duplicate
//   |d1 - d2| >  tolerance            conflict (rejected outright)
//
// The middle category is a judgement call, so it is counted separately and each
// occurrence is recorded in `discarded_targets` with both distances and the
// tolerance that was applied. The first value encountered wins.
struct ConstraintNormalizationDiagnostics {
    std::size_t exact_duplicates = 0;
    std::size_t tolerance_equivalent_duplicates = 0;
    // Repeats supplied with the atoms the other way round FROM THE FIRST
    // OCCURRENCE of that pair: a later (j, i) against an earlier (i, j). This is
    // a property of the two entries relative to each other, not of either one on
    // its own -- (1, 0) followed by (0, 1) is just as reversed as (0, 1)
    // followed by (1, 0), while (1, 0) twice is not reversed at all. The
    // orientation of each pair's first occurrence is therefore remembered and
    // compared against. Counted separately because the pair normalisation makes
    // the two indistinguishable afterwards, and a topology that lists a bond
    // twice in opposite orders is usually a generation bug worth surfacing.
    // Every reversed repeat is also counted in one of the two categories above.
    std::size_t reversed_duplicates = 0;

    // One entry per tolerance-equivalent duplicate whose target differed from
    // the value that was kept. Human-readable; the application decides how to
    // surface these (gmd prints them once at start-up).
    std::vector<std::string> discarded_targets;

    bool empty() const noexcept {
        return exact_duplicates == 0 && tolerance_equivalent_duplicates == 0;
    }
};

// One connected component of the constraint graph, and the rank of its own
// constraint rows.
//
// Independence is a per-component property. Two constraints that share no atom
// have Jacobian rows with disjoint support, so they are trivially independent
// and the total rank is the sum over components. Working component by component
// keeps the linear algebra on matrices the size of a molecule rather than of the
// whole system.
struct ConstraintComponent {
    std::vector<int> atom_tags;                   // sorted, globally tagged
    std::vector<std::size_t> constraint_indices;  // into ConstraintSolver::constraints()

    std::size_t rank = 0;              // numerical rank of this component's rows
    double largest_singular_value = 0.0;
    double smallest_singular_value = 0.0;
    double rank_tolerance = 0.0;       // the threshold rank was decided against
    // largest / smallest. Infinite when the component is rank deficient.
    double condition_number = 0.0;

    // Which of THIS component's constraints are redundant, as indices into
    // ConstraintSolver::constraints().
    //
    // WHICH ones are redundant is not unique: within a dependent group any one
    // member could be called the extra one. These are the columns left over by a
    // rank-revealing factorisation with column pivoting -- at each step it takes
    // the constraint whose gradient is most nearly orthogonal to those already
    // chosen, with numerically tied candidates settled by the canonical atom pair
    // (min tag, max tag) rather than by position in the input. The answer is
    // therefore deterministic and independent of the order the constraints were
    // supplied in.
    //
    // `redundant_constraints_identified` is false when the pivoted factorisation
    // did not agree with the singular values about HOW MANY rows are redundant.
    // In that case `redundant_constraint_indices` lists EVERY constraint in this
    // component instead, since nothing better can be said about which.
    bool redundant_constraints_identified = false;
    std::vector<std::size_t> redundant_constraint_indices;

    bool independent() const noexcept { return rank == constraint_indices.size(); }
};

// Result of testing a constraint set for independence at a given geometry.
//
// THE CRITERION. A holonomic constraint sigma_c(r) = |r_c|^2 - d_c^2 removes one
// degree of freedom only if its gradient is linearly independent of the others.
// The gradients must be compared in the metric the dynamics actually uses, which
// is mass-weighted: the constrained equations of motion involve J M^-1 J^T, so
// the relevant matrix is the mass-weighted Jacobian
//
//     J_M = J M^(-1/2),      J_M[c, 3i+a] = +2 r_c[a] / sqrt(m_i)
//                            J_M[c, 3j+a] = -2 r_c[a] / sqrt(m_j)
//
// with i and j the two atoms of constraint c. The number of degrees of freedom
// the set actually removes is rank(J_M), never the number of constraints
// supplied. J M^-1 J^T is singular exactly when J_M is rank deficient, which is
// also when the SHAKE/RATTLE iteration has no unique multiplier to converge to.
//
// The rank is CONFIGURATION DEPENDENT: a set can be independent at one geometry
// and degenerate at another (three collinear atoms, a flattened ring). This
// report describes the geometry it was computed at, which is the initial one.
struct ConstraintRankReport {
    bool independent = false;          // rank == constraint count everywhere
    std::size_t constraint_count = 0;
    std::size_t rank = 0;              // summed over components
    std::vector<ConstraintComponent> components;

    // Human-readable, each naming the atom tags and the component involved.
    // `problems` are hard failures; `warnings` are sets that are independent but
    // close enough to degenerate that the multipliers will be poorly determined.
    std::vector<std::string> problems;
    std::vector<std::string> warnings;

    std::size_t redundant_count() const noexcept { return constraint_count - rank; }

    // Whether EVERY rank-deficient component had its redundant rows identified.
    //
    // Derived from the components rather than stored, so a report describing a
    // mixed result -- one component identified, another fallen back to naming
    // its whole membership -- cannot claim to be fully identified. Vacuously
    // true when the set is independent, where there is nothing to identify.
    bool redundant_identified() const noexcept {
        for (const auto& component : components) {
            if (!component.independent() && !component.redundant_constraints_identified) {
                return false;
            }
        }
        return true;
    }

    // Every redundant constraint across all components, as indices into
    // ConstraintSolver::constraints(). Empty when the set is independent. Also
    // derived, so it cannot disagree with the per-component lists it summarises.
    // Where a component fell back, its entire membership appears here.
    std::vector<std::size_t> dependent_constraints() const {
        std::vector<std::size_t> all;
        for (const auto& component : components) {
            all.insert(all.end(), component.redundant_constraint_indices.begin(),
                       component.redundant_constraint_indices.end());
        }
        return all;
    }
};

// Singular values, in descending order, of a small dense matrix given by its
// COLUMNS, computed by one-sided Jacobi.
//
// This is the kernel the rank analysis is built on. It is exposed so that it can
// be validated against matrices whose singular values are known independently --
// a decomposition that is only ever compared against itself is not validated.
//
// Throws std::runtime_error if the sweeps do not reach mutual orthogonality; a
// rank must never be decided from a silently unconverged decomposition. `label`
// is folded into that message so the caller can say which constraints were
// involved.
std::vector<double> jacobi_singular_values(std::vector<std::vector<double>> columns,
                                           const std::string& label = "matrix");

// Holonomic bond-length constraints solved with SHAKE/RATTLE.
//
// Constraint list normalisation (performed by the constructor):
//   - atom pairs are stored as (min, max), so (i, j) and (j, i) are the same
//     constraint and are normalised identically regardless of the orientation
//     they were supplied in;
//   - repeated pairs are collapsed to one entry, keeping the FIRST target
//     distance encountered, which makes the result deterministic and
//     independent of pair orientation;
//   - a repeat whose target differs by more than the constraint tolerance is
//     rejected as a conflict, since no geometry satisfies both;
//   - negative atom tags, self-constraints (i == j) and non-positive or
//     non-finite target distances are rejected.
//
// A repeat whose target differs by more than zero but at most the tolerance is
// a "tolerance-equivalent duplicate": SHAKE converges to within `tolerance`, so
// the two targets are not distinguishable by the solver. It is accepted, the
// first value is kept, and it is reported through
// normalization_diagnostics().
//
// INDEPENDENCE IS VERIFIED, NOT ASSUMED, and a dependent set is REJECTED.
// Distinct is not the same as independent: a closed topology such as every pair
// among four or more atoms, or any three collinear atoms, contains more distance
// constraints than the rigid body has removable degrees of freedom. Counting
// those against the degrees of freedom would over-subtract and report a
// temperature that is too high, and the redundant rows make J M^-1 J^T singular,
// so the multipliers SHAKE and RATTLE converge to are not unique either.
//
// analyze_independence() computes the rank of the mass-weighted constraint
// Jacobian at a given geometry, per connected component; require_independent()
// throws unless the rank equals the constraint count. The integrator calls the
// latter once, before any dynamics run. THE POLICY IS TO REJECT: after it
// returns, active_constraint_count() IS the number of degrees of freedom the set
// removes, so the DOF subtraction is exact rather than assumed.
//
// Rejecting rather than silently using the rank is the deliberate choice. A
// dependent constraint set is nearly always a topology or input error, and
// running it with a quietly corrected DOF count would hide that while leaving
// the non-unique multipliers in place. See ConstraintRankReport.
class ConstraintSolver {
public:
    ConstraintSolver() = default;
    ConstraintSolver(std::vector<BondConstraint> constraints,
                     ConstraintSettings settings);

    bool enabled() const noexcept { return !constraints_.empty(); }
    const ConstraintSettings& settings() const noexcept { return settings_; }

    // The normalised list: distinct pairs, each stored as (min, max).
    const std::vector<BondConstraint>& constraints() const noexcept { return constraints_; }

    // Number of distinct active constraints, i.e. the number of degrees of
    // freedom removed *assuming the configured constraints are independent*.
    std::size_t active_constraint_count() const noexcept { return constraints_.size(); }

    // Total repeats dropped, of either kind.
    std::size_t dropped_duplicate_count() const noexcept {
        return diagnostics_.exact_duplicates +
               diagnostics_.tolerance_equivalent_duplicates;
    }

    // Per-category counts and the human-readable record of any non-identical
    // target that was discarded. ConstraintSolver deliberately does not print
    // anything itself: it is a low-level library class and the project has no
    // logging abstraction, so diagnostics are exposed to the caller instead.
    // `gmd` reports them on stdout once, right after the constraint summary.
    const ConstraintNormalizationDiagnostics& normalization_diagnostics() const noexcept {
        return diagnostics_;
    }

    // Rank of the mass-weighted constraint Jacobian at `system`'s geometry,
    // component by component, with diagnostics naming the atoms involved. See
    // ConstraintRankReport for the criterion. Collective under MPI: the atom
    // gather happens inside, and every rank computes the identical report.
    ConstraintRankReport analyze_independence(const System& system) const;

    // analyze_independence(), throwing a diagnostic listing every problem unless
    // the set is independent. Warnings (ill-conditioned but independent sets) are
    // returned in the report rather than thrown, for the caller to surface.
    //
    // `context` names the geometry in the message. The rank is a property of the
    // configuration, so saying which one was analysed matters: the integrator
    // calls this on the PROJECTED geometry, not the supplied one.
    ConstraintRankReport require_independent(
        const System& system,
        const std::string& context = "the analysed geometry") const;

    // The constraint bond vectors at the current geometry, to be handed to the
    // SHAKE projection that follows a position update. Call this BEFORE the
    // drift. Collective under MPI; the result is replicated on every rank.
    ConstraintReference capture_reference(const System& system) const;

    // Standard SHAKE: solves for the multipliers Gamma_c that put the drifted
    // positions back on the constraint manifold, displacing along the REFERENCE
    // bonds, and applies the matching half-step velocity impulse dr/dt. This is
    // steps (2) and (3) of the splitting documented on ConstraintVirialResult.
    // `time_step` is the dt the drift used and must be strictly positive.
    //
    // Contributes no virial: its multiplier is the constraint force at time
    // level t, while the reported virial is evaluated at t+dt.
    ConstraintProjectionStats apply_shake(System& system,
                                          const ConstraintReference& reference,
                                          double time_step) const;

    // Geometric projection of positions onto the constraint manifold, with no
    // reference geometry and no velocity change. This is NOT a step of the
    // dynamics: it is for putting an initial or a barostat-rescaled state onto
    // the manifold, where no preceding configuration exists to correct along.
    ConstraintProjectionStats apply_shake(System& system) const;

    // Velocity projection (RATTLE). The three-argument form additionally
    // recovers the step's ENDPOINT constraint virial from the converged
    // multipliers; `time_step` is the dt the step used and must be strictly
    // positive, since the constraint force is an impulse divided by it.
    ConstraintProjectionStats apply_rattle(System& system) const;
    ConstraintProjectionStats apply_rattle(System& system,
                                           double time_step,
                                           ConstraintVirialResult& virial_out) const;

private:
    ConstraintProjectionStats shake_impl(System& system,
                                         const ConstraintReference* reference,
                                         double time_step) const;
    ConstraintProjectionStats rattle_impl(System& system,
                                          double time_step,
                                          ConstraintVirialResult* virial_out) const;

    std::vector<BondConstraint> constraints_;
    ConstraintSettings settings_;
    ConstraintNormalizationDiagnostics diagnostics_;
};

std::vector<BondConstraint> constraints_from_bond_types(
    const Topology& topology,
    const std::vector<int>& constrained_bond_types,
    const std::vector<double>& bond_type_distances);

}  // namespace gmd
