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

    // One entry per tolerance-equivalent duplicate whose target differed from
    // the value that was kept. Human-readable; the application decides how to
    // surface these (gmd prints them once at start-up).
    std::vector<std::string> discarded_targets;

    bool empty() const noexcept {
        return exact_duplicates == 0 && tolerance_equivalent_duplicates == 0;
    }
};

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
// INDEPENDENCE IS ASSUMED, NOT VERIFIED. After normalisation the list contains
// distinct constraints, but distinct is not the same as independent: a closed
// topology such as an all-pairs cage over five or more atoms contains more
// distance constraints than the rigid body has removable degrees of freedom.
// Deciding that in general means computing the rank of the 3N x M constraint
// Jacobian, which is configuration dependent (the rank can drop at particular
// geometries) and is not attempted here. A redundant set is accepted and every
// distinct constraint is counted, so degrees_of_freedom() over-subtracts and
// the reported temperature comes out high. Configure independent constraints.
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
