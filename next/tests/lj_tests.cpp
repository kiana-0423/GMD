#include "gmd_next/reference.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <vector>

namespace {
using namespace gmd_next::reference;

void compare(const PairResult& a, const PairResult& b) {
    near(a.energy, b.energy);
    expect(a.forces.size() == b.forces.size(), "Force size mismatch");
    for (std::size_t i = 0; i < a.forces.size(); ++i) {
        for (std::size_t axis = 0; axis < 3; ++axis) near(a.forces[i][axis], b.forces[i][axis]);
    }
    for (std::size_t c = 0; c < 9; ++c) near(a.virial[c], b.virial[c]);
}

std::vector<Pair> half_pairs(std::size_t size) {
    std::vector<Pair> pairs;
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = i + 1; j < size; ++j) pairs.push_back({i, j});
    }
    return pairs;
}

void analytic_pair() {
    const std::array<Vec3, 2> x{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}};
    const std::array<Pair, 1> pairs{{{0, 1}}};
    const auto result = evaluate_lj(x, pairs, {1.0, 1.0});
    near(result.energy, 0.0);
    near(result.forces[0][0], -24.0);
    near(result.forces[1][0], 24.0);
    near(result.virial[0], 24.0);
    for (std::size_t c = 1; c < 9; ++c) near(result.virial[c], 0.0);
    auto minimum = x;
    minimum[1][0] = std::pow(2.0, 1.0 / 6.0);
    const auto at_minimum = evaluate_lj(minimum, pairs, {1.0, 1.0});
    near(at_minimum.energy, -1.0);
    near(at_minimum.forces[0][0], 0.0);
}

void gradients_and_representation() {
    const std::vector<Vec3> x{{0.1, 0.2, 0.0}, {1.4, 0.1, 0.2}, {0.6, 1.45, -0.1},
                             {1.8, 1.4, 0.65}, {-0.8, 0.7, 0.5}};
    const auto pairs = half_pairs(x.size());
    const LennardJones plain{0.8, 0.9};
    const auto baseline = evaluate_lj(x, pairs, plain);

    // Force finite differences exercise all coordinates and both cutoff models.
    // Configurations avoid the nonsmooth potential-shift cutoff surface.
    for (const auto model : {plain, LennardJones{0.8, 0.9, CutoffMode::potential_shift, 2.5},
                             LennardJones{0.8, 0.9, CutoffMode::force_shift, 2.5}}) {
        const auto result = evaluate_lj(x, pairs, model);
        constexpr double h = 1e-6;
        for (std::size_t i = 0; i < x.size(); ++i) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                auto plus = x;
                auto minus = x;
                plus[i][axis] += h;
                minus[i][axis] -= h;
                const double fd = -(evaluate_lj(plus, pairs, model).energy -
                                    evaluate_lj(minus, pairs, model).energy) / (2.0 * h);
                near(result.forces[i][axis], fd, 2e-8, 2e-8);
            }
        }
        // Uniform deformation X -> X(I+epsilon), including shear components.
        for (std::size_t a = 0; a < 3; ++a) {
            for (std::size_t b = 0; b < 3; ++b) {
                auto plus = x;
                auto minus = x;
                for (std::size_t i = 0; i < x.size(); ++i) {
                    plus[i][b] += h * x[i][a];
                    minus[i][b] -= h * x[i][a];
                }
                const double fd = -(evaluate_lj(plus, pairs, model).energy -
                                    evaluate_lj(minus, pairs, model).energy) / (2.0 * h);
                near(result.virial[3 * a + b], fd, 2e-8, 2e-8);
            }
        }
        for (std::size_t axis = 0; axis < 3; ++axis) {
            double total = 0.0;
            for (const auto& force : result.forces) total += force[axis];
            near(total, 0.0);
        }
    }

    auto reversed = pairs;
    for (auto& pair : reversed) std::swap(pair.source, pair.target);
    std::reverse(reversed.begin(), reversed.end());
    compare(baseline, evaluate_lj(x, reversed, plain));
    auto scaled = pairs;
    for (auto& pair : scaled) pair.scale = 0.25;
    auto quarter = baseline;
    quarter.energy *= 0.25;
    for (auto& force : quarter.forces) for (auto& v : force) v *= 0.25;
    for (auto& v : quarter.virial) v *= 0.25;
    compare(quarter, evaluate_lj(x, scaled, plain));

    // Independent source-local full-list formula: half energy/virial, full force.
    PairResult full;
    full.forces.resize(x.size(), Vec3{});
    for (std::size_t i = 0; i < x.size(); ++i) {
        for (std::size_t j = 0; j < x.size(); ++j) {
            if (i == j) continue;
            Vec3 d{};
            double r2 = 0.0;
            for (std::size_t a = 0; a < 3; ++a) {
                d[a] = x[j][a] - x[i][a];
                r2 += d[a] * d[a];
            }
            const double s6 = std::pow(plain.sigma * plain.sigma / r2, 3);
            full.energy += 2.0 * plain.epsilon * (s6 * s6 - s6);
            const double factor = 24.0 * plain.epsilon * (s6 - 2.0 * s6 * s6) / r2;
            for (std::size_t a = 0; a < 3; ++a) {
                full.forces[i][a] += factor * d[a];
                for (std::size_t b = 0; b < 3; ++b) {
                    full.virial[3 * a + b] -= 0.5 * d[a] * factor * d[b];
                }
            }
        }
    }
    compare(baseline, full);

    // Permute coordinates and index metadata together.
    const std::array<std::size_t, 5> old_to_new{3, 0, 4, 1, 2};
    auto permuted_x = x;
    auto permuted_pairs = pairs;
    auto expected = baseline;
    for (std::size_t i = 0; i < x.size(); ++i) {
        permuted_x[old_to_new[i]] = x[i];
        expected.forces[old_to_new[i]] = baseline.forces[i];
    }
    for (auto& pair : permuted_pairs) {
        pair.source = old_to_new[pair.source];
        pair.target = old_to_new[pair.target];
    }
    compare(expected, evaluate_lj(permuted_x, permuted_pairs, plain));
}

void cutoff_periodic_and_errors() {
    const std::array<Pair, 1> pairs{{{0, 1}}};
    const LennardJones shifted{1.0, 1.0, CutoffMode::potential_shift, 2.5};
    const LennardJones force_shifted{1.0, 1.0, CutoffMode::force_shift, 2.5};
    std::array<Vec3, 2> x{{{0.0, 0.0, 0.0}, {2.5, 0.0, 0.0}}};
    for (const auto model : {shifted, force_shifted}) {
        const auto result = evaluate_lj(x, pairs, model);
        near(result.energy, 0.0);
        near(result.forces[0][0], 0.0);
        x[1][0] = 3.0;
        near(evaluate_lj(x, pairs, model).energy, 0.0);
        x[1][0] = 2.5;
    }
    x[1][0] = 2.5 - 1e-7;
    expect(std::abs(evaluate_lj(x, pairs, shifted).forces[0][0]) > 1e-3,
           "Potential shift must retain the original force below cutoff");
    near(evaluate_lj(x, pairs, force_shifted).forces[0][0], 0.0, 1e-7, 0.0);

    const OrthorhombicCell cell{{10.0, 10.0, 10.0}};
    x = {{{0.2, 0.3, 0.4}, {9.0, 0.6, 0.8}}};
    auto unwrapped = x;
    unwrapped[1][0] -= 10.0;
    const auto periodic = evaluate_lj(x, pairs, shifted, cell);
    compare(periodic, evaluate_lj(unwrapped, pairs, shifted));
    x[1][0] += 20.0;
    x[0][1] -= 10.0;
    compare(periodic, evaluate_lj(x, pairs, shifted, cell));

    const auto empty = evaluate_lj({}, {}, {1.0, 1.0});
    near(empty.energy, 0.0);
    expect(empty.forces.empty(), "Empty system forces must be empty");
    x = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    const std::array<Pair, 1> excluded{{{0, 1, 0.0}}};
    near(evaluate_lj(x, excluded, shifted).energy, 0.0);
    throws<std::domain_error>([&] { evaluate_lj(x, pairs, shifted); });
    x[1][0] = 1.0;
    const std::array<Pair, 2> duplicate{{{0, 1}, {1, 0}}};
    throws<std::invalid_argument>([&] { evaluate_lj(x, duplicate, shifted); });
    const std::array<Pair, 1> self{{{0, 0}}};
    throws<std::invalid_argument>([&] { evaluate_lj(x, self, shifted); });
    const std::array<Pair, 1> invalid{{{0, 2}}};
    throws<std::out_of_range>([&] { evaluate_lj(x, invalid, shifted); });
    throws<std::invalid_argument>([&] { evaluate_lj(x, pairs, {-1.0, 1.0}); });
    throws<std::invalid_argument>([&] { evaluate_lj(x, pairs, {1.0, 1.0, CutoffMode::none, 2.5}); });
    throws<std::invalid_argument>([&] { evaluate_lj(x, pairs, {1.0, 1.0}, cell); });
    throws<std::invalid_argument>([&] {
        evaluate_lj(x, pairs, shifted, OrthorhombicCell{{5.0, 10.0, 10.0}});
    });
    x[0][0] = std::numeric_limits<double>::quiet_NaN();
    throws<std::invalid_argument>([&] { evaluate_lj(x, pairs, shifted); });
    std::cout << "Analytic LJ, gradients, virial, counting, permutation, PBC and errors passed\n";
}
}  // namespace

int main() {
    try {
        analytic_pair();
        gradients_and_representation();
        cutoff_periodic_and_errors();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
