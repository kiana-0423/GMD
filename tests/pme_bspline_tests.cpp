// Direct regression tests for the PRODUCTION cardinal B-spline kernel:
//
//     PMEForceProvider::bspline(u, order)
//     PMEForceProvider::bspline_deriv(u, order)
//
// These are the functions the engine actually evaluates, not a reference
// standing in for them. That distinction is the whole reason this file exists:
// the reference implementation in tests/pme_reference.hpp was correct all
// along, and the production polynomials were not.
//
// WHAT WENT WRONG, AND WHAT WOULD HAVE CAUGHT IT
//
//   1. bspline() implemented orders 4 and 6 only and returned 0.0 for anything
//      else. bspline_deriv(u, p) evaluates M_(p-1), so order 4 asked for order
//      3 and order 6 asked for order 5 -- both got zero. Every derivative
//      weight in the force interpolation was zero, so the PME reciprocal force
//      vanished entirely, at every supported order. Caught here by
//      test_orders_3_and_5_are_supported() and by the derivative identity.
//
//   2. The order-6 polynomials on [2,3), [3,4) and [4,5) were not M_6. The
//      spline took negative values, partition of unity summed to 0.9 instead
//      of 1, and M_6(u) = M_6(6-u) failed. Caught here by non-negativity,
//      partition of unity, symmetry, the recursion comparison, and by the
//      explicit regression values in test_order_six_regression_values().
//
// Orders 3, 4, 5 and 6 are all covered. 3 and 5 are not spare capacity: they
// are what the derivatives of the supported production orders are made of, and
// leaving them untested is what allowed defect 1 to survive.

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "gmd/force/pme_force_provider.hpp"

#include "pme_reference.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[pme bspline] " << message << '\n';
        ++failures;
    }
}

std::string format(double value, int precision = 3) {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(precision) << value;
    return stream.str();
}

// The production entry points, named once so every test below is unambiguously
// exercising the shipped code.
double production_bspline(double u, int order) {
    return gmd::PMEForceProvider::bspline(u, order);
}

double production_bspline_deriv(double u, int order) {
    return gmd::PMEForceProvider::bspline_deriv(u, order);
}

constexpr std::array<int, 5> kOrders = {2, 3, 4, 5, 6};

// Evaluation tolerance, per order.
//
// These polynomials are expanded in the monomial basis with coefficients that
// grow fast: order 6 on [4,5) carries 12270 and 10974, and M_6 there is about
// 0.2. Evaluating it therefore cancels five significant figures away, and the
// result carries an absolute round-off of order 1e-15 * 12270 ~ 1e-11. A
// tolerance tighter than that would be measuring the monomial basis, not the
// spline. Observed worst deviations are 2e-13 (order 6) and 4e-14 (order 5).
double evaluation_tolerance(int order) {
    return order >= 5 ? 5.0e-12 : 1.0e-14;
}

// The d-th derivative of M_p, from the difference formula
//   M_p^(d)(u) = sum_{j=0}^{d} (-1)^j C(d,j) M_(p-d)(u - j)
// evaluated through the PRODUCTION spline of order p-d. Defined for
// d <= p - 2, where M_(p-d) is still continuous.
double production_derivative(double u, int order, int d) {
    double sum = 0.0;
    double binomial = 1.0;
    for (int j = 0; j <= d; ++j) {
        sum += ((j % 2 == 0) ? 1.0 : -1.0) * binomial
             * gmd::PMEForceProvider::bspline(u - j, order - d);
        binomial = binomial * (d - j) / (j + 1);
    }
    return sum;
}

// ===========================================================================
// Support, boundaries, and the intervals themselves
// ===========================================================================

void test_zero_outside_support() {
    for (const int order : kOrders) {
        const double p = order;
        for (const double u : {-5.0, -1.0, -1e-9, 0.0, p, p + 1e-9, p + 1.0, p + 5.0}) {
            check(production_bspline(u, order) == 0.0,
                  "order " + std::to_string(order) + ": M_p(" + format(u) +
                      ") must be exactly zero outside (0, p), got " +
                      format(production_bspline(u, order)));
        }
    }
}

// Every interval [k, k+1) is a different polynomial branch. Sampling each one
// densely is what makes "the polynomials are right" a statement about all of
// them rather than about whichever the fixtures happened to hit.
void test_every_interval_against_the_recursion() {
    constexpr int kSamplesPerInterval = 401;
    for (const int order : kOrders) {
        for (int interval = 0; interval < order; ++interval) {
            double worst = 0.0;
            double worst_at = 0.0;
            for (int sample = 0; sample < kSamplesPerInterval; ++sample) {
                // Strictly inside the interval, so this test is about the
                // branch and the boundary test below is about the seams.
                const double t = (sample + 0.5) / kSamplesPerInterval;
                const double u = interval + t;
                const double produced = production_bspline(u, order);
                const auto expected = static_cast<double>(
                    pme_ref::cardinal_bspline(static_cast<pme_ref::Real>(u), order));
                const double error = std::fabs(produced - expected);
                if (error > worst) {
                    worst = error;
                    worst_at = u;
                }
            }
            check(worst < evaluation_tolerance(order),
                  "order " + std::to_string(order) + ", interval [" +
                      std::to_string(interval) + "," + std::to_string(interval + 1) +
                      "): production polynomial disagrees with the recursion by " +
                      format(worst) + " at u = " + format(worst_at));
        }
    }
}

// The seams. A branch whose domain test is off by one, or a polynomial that is
// right in the interior and wrong at the join, shows up here and nowhere else.
void test_interval_boundaries_and_both_sides() {
    for (const int order : kOrders) {
        for (int knot = 0; knot <= order; ++knot) {
            const double u = knot;
            const auto expected = static_cast<double>(
                pme_ref::cardinal_bspline(static_cast<pme_ref::Real>(u), order));
            check(std::fabs(production_bspline(u, order) - expected)
                      < evaluation_tolerance(order),
                  "order " + std::to_string(order) + ": M_p at the exact knot u = " +
                      std::to_string(knot) + " is " +
                      format(production_bspline(u, order)) + ", expected " +
                      format(expected));

            // Immediately either side: the two branches must agree with each
            // other in the limit, which is continuity across the knot.
            for (const double epsilon : {1.0e-11, 1.0e-9, 1.0e-7}) {
                const double below = production_bspline(u - epsilon, order);
                const double above = production_bspline(u + epsilon, order);
                check(std::fabs(below - above) < 1.0e-5,
                      "order " + std::to_string(order) + ": M_p is discontinuous at knot " +
                          std::to_string(knot) + " (below " + format(below) +
                          ", above " + format(above) + ", epsilon " + format(epsilon) + ")");
            }
        }
    }
}

// ===========================================================================
// Structural invariants
// ===========================================================================

void test_non_negativity() {
    for (const int order : kOrders) {
        double most_negative = 0.0;
        double at = 0.0;
        for (int sample = 0; sample <= 20000; ++sample) {
            const double u = order * sample / 20000.0;
            const double value = production_bspline(u, order);
            if (value < most_negative) {
                most_negative = value;
                at = u;
            }
        }
        check(most_negative >= -1.0e-15,
              "order " + std::to_string(order) +
                  ": a cardinal B-spline is non-negative on its support, but M_p(" +
                  format(at) + ") = " + format(most_negative));
    }
}

void test_symmetry() {
    for (const int order : kOrders) {
        double worst = 0.0;
        double at = 0.0;
        for (int sample = 0; sample <= 10000; ++sample) {
            const double u = order * sample / 10000.0;
            const double error =
                std::fabs(production_bspline(u, order) - production_bspline(order - u, order));
            if (error > worst) {
                worst = error;
                at = u;
            }
        }
        check(worst < evaluation_tolerance(order),
              "order " + std::to_string(order) +
                  ": M_p(u) = M_p(p-u) fails by " + format(worst) + " at u = " + format(at));
    }
}

// sum_k M_p(t + k) == 1 for every fractional offset t. This is the property
// the charge-spreading step depends on -- a spline that violates it does not
// conserve the charge it deposits on the mesh.
void test_partition_of_unity() {
    constexpr int kOffsets = 1997;  // prime, so the offsets do not align with knots
    for (const int order : kOrders) {
        double worst = 0.0;
        double at = 0.0;
        for (int index = 0; index < kOffsets; ++index) {
            const double t = static_cast<double>(index) / kOffsets;
            double sum = 0.0;
            for (int k = 0; k < order; ++k) sum += production_bspline(t + k, order);
            const double error = std::fabs(sum - 1.0);
            if (error > worst) {
                worst = error;
                at = t;
            }
        }
        check(worst < evaluation_tolerance(order),
              "order " + std::to_string(order) +
                  ": sum_k M_p(t+k) departs from 1 by " + format(worst) +
                  " at t = " + format(at));
    }
}

// M_p is a piecewise polynomial of degree p-1 joined with C^(p-2) continuity.
// Checking derivatives 0 .. p-2 across every interior knot pins that down; the
// (p-1)-th derivative is genuinely discontinuous and is not checked.
void test_continuity_through_expected_derivative_order() {
    for (const int order : kOrders) {
        for (int derivative = 0; derivative <= order - 2; ++derivative) {
            for (int knot = 1; knot < order; ++knot) {
                // One-sided limits of the d-th derivative, taken from the
                // closed-form difference formula rather than by repeated
                // numerical differentiation, so the test measures continuity
                // and not finite-difference noise.
                // Continuity means the gap closes with epsilon, so the bound
                // has to close with it too. A fixed tolerance would either
                // admit a genuine jump at small epsilon or reject the ordinary
                // O(epsilon) variation of a continuous function at large
                // epsilon. The floor is the round-off of the monomial
                // evaluation, below which the gap stops shrinking.
                for (const double epsilon : {1.0e-10, 1.0e-8, 1.0e-6}) {
                    const double below = production_derivative(knot - epsilon, order, derivative);
                    const double above = production_derivative(knot + epsilon, order, derivative);
                    const double gap = std::fabs(below - above);
                    const double bound = std::max(100.0 * epsilon, 1.0e-9);
                    check(gap < bound,
                          "order " + std::to_string(order) + ": derivative " +
                              std::to_string(derivative) + " is discontinuous at knot " +
                              std::to_string(knot) + " (one-sided limits differ by " +
                              format(gap) + " at epsilon " + format(epsilon) +
                              ", bound " + format(bound) + ")");
                }
            }
        }

        // The continuity claim above is only worth anything if the highest
        // derivative it checks is not identically zero. M_1 is not implemented
        // by this kernel, so the first genuinely discontinuous derivative
        // (order p-1) cannot be evaluated through the difference formula; what
        // can be established is that derivative p-2 carries real signal.
        if (order >= 3) {
            const int highest = order - 2;
            double largest = 0.0;
            for (int sample = 0; sample <= 4000; ++sample) {
                largest = std::max(largest,
                                   std::fabs(production_derivative(
                                       order * sample / 4000.0, order, highest)));
            }
            check(largest > 0.1,
                  "order " + std::to_string(order) + ": derivative " +
                      std::to_string(highest) +
                      " is identically zero over the whole support (largest " +
                      format(largest) + "), so the continuity check above passes "
                      "vacuously");
        }
    }
}

void test_derivative_identity() {
    for (const int order : kOrders) {
        // M_1 is the indicator function, not implemented here, so order 2's
        // derivative is out of scope. Every order the provider can reach, and
        // every order those reach through their derivatives, is covered.
        if (order < 3) continue;
        double worst = 0.0;
        double at = 0.0;
        for (int sample = 0; sample <= 20000; ++sample) {
            const double u = order * sample / 20000.0;
            const double produced = production_bspline_deriv(u, order);
            const double identity =
                production_bspline(u, order - 1) - production_bspline(u - 1.0, order - 1);
            const double error = std::fabs(produced - identity);
            if (error > worst) {
                worst = error;
                at = u;
            }
        }
        check(worst < 1.0e-15,
              "order " + std::to_string(order) +
                  ": bspline_deriv does not equal M_(p-1)(u) - M_(p-1)(u-1), worst " +
                  format(worst) + " at u = " + format(at));
    }
}

// The derivative must also be the derivative of the production spline itself,
// established by a coordinate central difference rather than by the identity
// above. Both checks were needed: the identity held trivially when both sides
// were zero.
void test_derivative_against_finite_difference() {
    // h = 1e-4, not 1e-6. The polynomial evaluation itself carries ~1e-12 of
    // absolute round-off at order 6 (see evaluation_tolerance), and a central
    // difference divides that by 2h: at h = 1e-6 the noise floor alone is
    // ~5e-7, which is larger than anything this test could hope to resolve.
    // At h = 1e-4 the noise floor is ~5e-9 and truncation is ~1e-8|M'''|.
    const double step = 1.0e-4;
    for (const int order : kOrders) {
        if (order < 3) continue;
        double worst_absolute = 0.0;
        double scale = 0.0;
        double at = 0.0;
        for (int sample = 1; sample < 20000; ++sample) {
            const double u = order * sample / 20000.0;
            // Skip a neighbourhood of the knots: the derivative of a spline of
            // order 3 is only C^0 there, so a central difference straddling a
            // knot measures the kink, not the derivative.
            bool near_knot = false;
            for (int knot = 0; knot <= order; ++knot) {
                if (std::fabs(u - knot) < 10.0 * step) near_knot = true;
            }
            if (near_knot) continue;

            const double numerical =
                (production_bspline(u + step, order) - production_bspline(u - step, order))
                / (2.0 * step);
            const double analytic = production_bspline_deriv(u, order);
            scale = std::max(scale, std::fabs(analytic));
            if (std::fabs(numerical - analytic) > worst_absolute) {
                worst_absolute = std::fabs(numerical - analytic);
                at = u;
            }
        }
        check(scale > 1.0e-3,
              "order " + std::to_string(order) +
                  ": the derivative is identically zero over the whole support, so this "
                  "comparison would pass vacuously");
        check(worst_absolute < 1.0e-6,
              "order " + std::to_string(order) +
                  ": bspline_deriv disagrees with a central difference of bspline by " +
                  format(worst_absolute) + " at u = " + format(at) + " (scale " +
                  format(scale) + ")");
        std::cout << "    order " << order << ": max |M'_p - central difference| = "
                  << format(worst_absolute) << " (|M'_p| up to " << format(scale) << ")\n";
    }
}

// The specific defect: orders 3 and 5 returning zero made every derivative
// weight vanish. State it directly, so a regression cannot hide behind an
// aggregate.
void test_orders_3_and_5_are_supported() {
    for (const int order : {2, 3, 5}) {
        double largest = 0.0;
        for (int sample = 0; sample <= 1000; ++sample) {
            largest = std::max(largest,
                               std::fabs(production_bspline(order * sample / 1000.0, order)));
        }
        check(largest > 0.1,
              "order " + std::to_string(order) +
                  " is not implemented (largest value over its support is " +
                  format(largest) + "). Orders 4 and 6 need it: bspline_deriv(u, p) "
                  "evaluates M_(p-1), and an unimplemented M_(p-1) zeroes every PME "
                  "reciprocal force.");
    }

    // And the consequence, stated at the level that actually broke: the
    // derivative weights of the supported production orders must not be zero.
    for (const int order : {4, 6}) {
        double largest = 0.0;
        for (int sample = 0; sample <= 1000; ++sample) {
            largest = std::max(largest,
                               std::fabs(production_bspline_deriv(order * sample / 1000.0, order)));
        }
        check(largest > 0.1,
              "order " + std::to_string(order) +
                  ": every derivative weight is zero (largest " + format(largest) +
                  "), which is exactly the state in which PME applies no reciprocal force");
    }
}

// ===========================================================================
// Explicit regression values
// ===========================================================================

// The three order-6 intervals that were wrong, with values computed from
// M_6(u) = 1/120 sum_j (-1)^j C(6,j) (u-j)_+^5 by hand. Written out as literals
// rather than derived at run time, so that this test still fails if both the
// production code and the reference recursion were changed together.
//
// For reference, the superseded polynomials produced -0.183333 at u = 2,
// 0.350000 at u = 3 and -0.283333 at u = 4 -- negative values for a
// non-negative function.
void test_order_six_regression_values() {
    struct Expectation {
        double u;
        double value;
    };
    const std::vector<Expectation> expectations = {
        // [2,3): (10u^5 - 120u^4 + 540u^3 - 1140u^2 + 1170u - 474) / 120
        {2.00, 26.0 / 120.0},                 // 0.21666666666666667
        {2.50, 0.43802083333333333},
        {2.75, 0.51964518229166667},
        // [3,4): (-10u^5 + 180u^4 - 1260u^3 + 4260u^2 - 6930u + 4386) / 120
        {3.00, 66.0 / 120.0},                 // 0.55
        {3.25, 0.51964518229166667},
        {3.50, 0.43802083333333333},
        // [4,5): (5u^5 - 120u^4 + 1140u^3 - 5340u^2 + 12270u - 10974) / 120
        {4.00, 26.0 / 120.0},
        {4.25, 0.12491048177083333},
        {4.75, 0.02538248697916667},
    };

    for (const auto& expectation : expectations) {
        const double produced = production_bspline(expectation.u, 6);
        const double error = std::fabs(produced - expectation.value);
        check(error < 1.0e-14,
              "order 6 regression: M_6(" + format(expectation.u, 2) + ") = " +
                  format(produced, 12) + ", expected " + format(expectation.value, 12) +
                  " (error " + format(error) + ")");
    }

    // The middle three intervals must not be negative anywhere -- the single
    // most visible symptom of the superseded polynomials.
    double most_negative = 0.0;
    for (int sample = 0; sample <= 30000; ++sample) {
        const double u = 2.0 + 3.0 * sample / 30000.0;
        most_negative = std::min(most_negative, production_bspline(u, 6));
    }
    check(most_negative >= -1.0e-15,
          "order 6: M_6 is negative on [2,5], minimum " + format(most_negative) +
              " -- the signature of the superseded quintic polynomials");
}

// Order 4 was correct throughout and must stay that way; pinning a few values
// keeps a future "cleanup" of the polynomial table honest.
void test_order_four_regression_values() {
    struct Expectation {
        double u;
        double value;
    };
    const std::vector<Expectation> expectations = {
        {0.50, 1.0 / 48.0},        // u^3/6
        {1.00, 1.0 / 6.0},
        {1.50, 23.0 / 48.0},
        {2.00, 4.0 / 6.0},
        {2.50, 23.0 / 48.0},
        {3.00, 1.0 / 6.0},
        {3.50, 1.0 / 48.0},
    };
    for (const auto& expectation : expectations) {
        const double produced = production_bspline(expectation.u, 4);
        check(std::fabs(produced - expectation.value) < 1.0e-15,
              "order 4 regression: M_4(" + format(expectation.u, 2) + ") = " +
                  format(produced, 12) + ", expected " + format(expectation.value, 12));
    }
}

}  // namespace

int main() {
    std::cout << "PME production B-spline kernel\n";
    test_zero_outside_support();
    test_every_interval_against_the_recursion();
    test_interval_boundaries_and_both_sides();
    test_non_negativity();
    test_symmetry();
    test_partition_of_unity();
    test_continuity_through_expected_derivative_order();
    test_derivative_identity();
    test_derivative_against_finite_difference();
    test_orders_3_and_5_are_supported();
    test_order_six_regression_values();
    test_order_four_regression_values();

    if (failures != 0) {
        std::cerr << "PME B-spline tests failed: " << failures << '\n';
        return 1;
    }
    std::cout << "PME B-spline tests passed\n";
    return 0;
}
