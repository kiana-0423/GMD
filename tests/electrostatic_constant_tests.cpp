// Audit of the electrostatic unit-conversion constant k_e.
//
// GMD works in eV, Angstrom and elementary charge, so the Coulomb energy of
// two unit charges one Angstrom apart is, numerically, k_e. Every electrostatic
// quantity the engine reports -- direct special-pair correction, Ewald real,
// Ewald reciprocal, Ewald self, Ewald net-charge background, and the same five
// for PME -- is exactly proportional to that one number, with exactly one
// factor of it. This file exists to prove that structural claim and to pin the
// numeric value against a derivation from primary constants.
//
// HOW k_e IS MEASURED RATHER THAN ASSUMED
//
// Nothing here reads the production constant: it has internal linkage in
// ewald_force_provider.cpp and pme_force_provider.cpp and is deliberately not
// exported for testing. Instead each electrostatic sum is recomputed here with
// k_e set to 1, and the constant is recovered as
//
//     k_e_measured = E_provider / Phi_reference(k_e = 1)
//
// which is exact because the dependence is exactly linear. A path that omitted
// the constant would measure 1; a path that applied it twice would measure
// k_e^2; a path with its own differently-rounded copy would disagree with the
// others in the last digits. All three are failures here.
//
// The PME reciprocal reference is taken from pme_reference.hpp, whose own
// declared constant is divided straight back out. That division makes this
// measurement independent of the value that header happens to declare: if both
// the header and production were wrong in the same way, the measured constant
// would still come out as production's value and the comparison against the
// CODATA derivation below would fail. A reference that shared production's
// constant *without* dividing it out would be self-confirming, which is
// precisely what this file must not be.
//
// THE AUTHORITATIVE VALUE
//
// k_e is not an adjustable parameter. In atomic units the Coulomb energy of two
// unit charges separated by one Bohr radius is exactly one Hartree, so
//
//     e^2 / (4 pi eps0) = E_h * a0
//
// and in GMD's units k_e[eV A/e^2] = E_h[eV] * a0[A]. Both factors are taken
// below from the 2022 CODATA adjustment, tabulated by NIST:
//
//     E_h = 27.211 386 245 981(30) eV     https://physics.nist.gov/cgi-bin/cuu/Value?hrev
//     a0  = 5.291 772 105 44(82) e-11 m   https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
//
// They are written here as two separate literals and multiplied, rather than as
// one pre-computed number, so that the derivation is visible and checkable
// rather than copied.

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "gmd/core/runtime_context.hpp"
#include "gmd/force/ewald_force_provider.hpp"
#include "gmd/force/pme_force_provider.hpp"
#include "gmd/system/box.hpp"
#include "gmd/system/special_pair_map.hpp"
#include "gmd/system/system.hpp"
#include "gmd/system/topology.hpp"

#include "pme_reference.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "[k_e audit] " << message << '\n';
        ++failures;
    }
}

std::string number(long double value, int digits = 17) {
    std::ostringstream out;
    out << std::setprecision(digits) << static_cast<double>(value);
    return out.str();
}

// --- the authoritative value, derived here from primary constants ----------

// CODATA 2022, NIST. Written separately and multiplied; see the header comment.
constexpr long double kHartreeEnergyEv    = 27.211386245981L;
constexpr long double kBohrRadiusAngstrom = 0.529177210544L;
constexpr long double kCodata2022 = kHartreeEnergyEv * kBohrRadiusAngstrom;

// What production currently uses. This is a pin measured out of the engine
// below, not a derivation.
constexpr long double kExpectedProduction = 14.3996L;

// CHARACTERIZATION, NOT ENDORSEMENT.
//
// The production constant is a five-significant-figure truncation of the
// CODATA derivation. This file records the exact size of that deviation so it
// is a tested fact rather than an assertion in a document, and so that any
// commit correcting the constant is forced to come through here.
//
// It is a real systematic error, not a rounding nicety: 3.16e-06 relative is
// larger than the tightest PME convergence point this repository measures
// (6.3e-10 relative, order 6 at grid 128), so it sets a floor on how well GMD
// can ever agree with an external engine no matter how fine the mesh.
constexpr long double kKnownDeviationFromCodata = -3.157635e-06L;
constexpr long double kDeviationTolerance = 1.0e-4L;   // relative, on the deviation

// --- a complete Ewald reference with k_e factored out ----------------------
//
// Real, reciprocal, self and neutralising-background terms for an orthorhombic
// cell, written from their definitions with the Coulomb constant omitted, so
// the return value is Phi = E / k_e. Independent of the engine and of every
// other reference header.

struct Ewald1 {
    long double energy = 0.0L;                          // Phi
    std::vector<std::array<long double, 3>> forces;     // dPhi/dr, negated
    std::array<long double, 9> virial{};                // row major
};

Ewald1 ewald_reference_unit_constant(
        const std::vector<std::array<long double, 3>>& positions,
        const std::vector<long double>& charges,
        const std::array<long double, 3>& lengths,
        long double alpha, int kmax, long double cutoff) {
    const std::size_t n = charges.size();
    Ewald1 out;
    out.forces.assign(n, {0.0L, 0.0L, 0.0L});

    const long double two_alpha_over_root_pi =
        2.0L * alpha / std::sqrt(std::numbers::pi_v<long double>);

    // Real space, minimum image, hard cutoff -- the same truncation the engine
    // applies, so the two are comparable term by term.
    for (std::size_t i = 0; i + 1 < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            std::array<long double, 3> dr{};
            for (std::size_t d = 0; d < 3; ++d) {
                dr[d] = positions[i][d] - positions[j][d];
                dr[d] -= lengths[d] * std::round(dr[d] / lengths[d]);
            }
            const long double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if (r2 >= cutoff * cutoff || r2 < 1.0e-24L) continue;
            const long double r = std::sqrt(r2);
            const long double ar = alpha * r;
            const long double qq = charges[i] * charges[j];
            out.energy += qq * std::erfc(ar) / r;
            const long double ff = qq / r2 *
                (std::erfc(ar) / r + two_alpha_over_root_pi * std::exp(-ar * ar));
            for (std::size_t d = 0; d < 3; ++d) {
                out.forces[i][d] += ff * dr[d];
                out.forces[j][d] -= ff * dr[d];
            }
            for (std::size_t a = 0; a < 3; ++a)
                for (std::size_t b = 0; b < 3; ++b)
                    out.virial[a * 3 + b] += dr[a] * (ff * dr[b]);
        }
    }

    // Reciprocal space.
    const long double volume = lengths[0] * lengths[1] * lengths[2];
    const long double inv_four_alpha_sq = 1.0L / (4.0L * alpha * alpha);
    const long double two_pi = 2.0L * std::numbers::pi_v<long double>;
    for (int nx = -kmax; nx <= kmax; ++nx) {
      for (int ny = -kmax; ny <= kmax; ++ny) {
        for (int nz = -kmax; nz <= kmax; ++nz) {
          if (nx == 0 && ny == 0 && nz == 0) continue;
          const std::array<long double, 3> k = {two_pi * nx / lengths[0],
                                                two_pi * ny / lengths[1],
                                                two_pi * nz / lengths[2]};
          const long double k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
          const long double g = 4.0L * std::numbers::pi_v<long double>
                              / (volume * k2) * std::exp(-k2 * inv_four_alpha_sq);
          long double sre = 0.0L, sim = 0.0L;
          for (std::size_t i = 0; i < n; ++i) {
              const long double phase = k[0]*positions[i][0] + k[1]*positions[i][1]
                                      + k[2]*positions[i][2];
              sre += charges[i] * std::cos(phase);
              sim += charges[i] * std::sin(phase);
          }
          const long double term = 0.5L * g * (sre * sre + sim * sim);
          out.energy += term;
          for (std::size_t i = 0; i < n; ++i) {
              const long double phase = k[0]*positions[i][0] + k[1]*positions[i][1]
                                      + k[2]*positions[i][2];
              const long double pref = g * charges[i]
                  * (sre * std::sin(phase) - sim * std::cos(phase));
              for (std::size_t d = 0; d < 3; ++d) out.forces[i][d] += pref * k[d];
          }
          // W_ab = E_k [ delta_ab - 2 (1/k^2 + 1/4alpha^2) k_a k_b ]
          const long double factor = 2.0L * (1.0L / k2 + inv_four_alpha_sq);
          for (std::size_t a = 0; a < 3; ++a)
              for (std::size_t b = 0; b < 3; ++b)
                  out.virial[a * 3 + b] +=
                      term * ((a == b ? 1.0L : 0.0L) - factor * k[a] * k[b]);
        }
      }
    }

    // Self term: depends on alpha and the charges only, never on the cell, so
    // it carries energy but contributes nothing to the virial.
    long double q2 = 0.0L, qnet = 0.0L;
    for (long double q : charges) { q2 += q * q; qnet += q; }
    out.energy -= alpha / std::sqrt(std::numbers::pi_v<long double>) * q2;

    // Neutralising background: scales as 1/V, so W_ab = U_net delta_ab.
    if (qnet != 0.0L) {
        const long double net =
            -std::numbers::pi_v<long double> / (2.0L * volume * alpha * alpha) * qnet * qnet;
        out.energy += net;
        out.virial[0] += net;
        out.virial[4] += net;
        out.virial[8] += net;
    }
    return out;
}

// --- fixtures --------------------------------------------------------------
//
// Deliberately asymmetric: a non-cubic box, four charges of four different
// magnitudes at irregular positions, no two coordinates equal, no mirror plane.
// A sign flip, an axis permutation or a transposed virial cannot survive it,
// because every force component and every virial component is distinct and
// substantially nonzero.

struct Fixture {
    std::vector<std::array<long double, 3>> positions;
    std::vector<long double> charges;
    std::array<long double, 3> lengths{};
};

Fixture asymmetric_fixture(bool neutral) {
    Fixture f;
    f.lengths = {17.0L, 21.0L, 25.0L};
    f.positions = {{2.13L, 3.41L, 4.77L},
                   {8.62L, 1.94L, 11.35L},
                   {13.08L, 9.27L, 6.51L},
                   {5.46L, 15.83L, 19.24L}};
    // Neutral variant sums to exactly zero in binary floating point; the
    // charged variant does not, which is what makes the neutralising-background
    // term fire and lets this file check that term too.
    f.charges = neutral ? std::vector<long double>{0.75L, -0.5L, 0.625L, -0.875L}
                        : std::vector<long double>{0.75L, -0.5L, 0.625L, 0.375L};
    return f;
}

gmd::System build_system(const Fixture& f) {
    gmd::System system;
    system.resize(f.charges.size(), f.charges.size());
    gmd::Box box;
    box.set_lengths({static_cast<double>(f.lengths[0]),
                     static_cast<double>(f.lengths[1]),
                     static_cast<double>(f.lengths[2])});
    system.set_box(box);
    for (std::size_t i = 0; i < f.charges.size(); ++i) {
        system.mutable_coordinates()[i] = {static_cast<double>(f.positions[i][0]),
                                           static_cast<double>(f.positions[i][1]),
                                           static_cast<double>(f.positions[i][2])};
        system.mutable_charges()[i] = static_cast<double>(f.charges[i]);
        system.mutable_masses()[i] = 1.0;
        system.mutable_atom_tags()[i] = static_cast<int>(i);
        system.mutable_atom_owners()[i] = 0;
    }
    return system;
}

gmd::ForceResult evaluate(gmd::ForceProvider& provider, const gmd::System& system) {
    gmd::RuntimeContext runtime;
    const auto coordinates = system.coordinates();
    gmd::ForceRequest request{
        .system = &system,
        .box = &system.box(),
        .coordinates = std::span<const gmd::Coordinate3D>(coordinates.data(),
                                                          coordinates.size()),
        .neighbor_list = nullptr,
    };
    gmd::ForceResult result;
    provider.initialize(runtime);
    provider.compute(request, result, runtime);
    provider.finalize(runtime);
    return result;
}

// --- measurements ----------------------------------------------------------

constexpr long double kAlpha = 0.32L;
constexpr long double kCutoff = 8.0L;
constexpr int kKmax = 12;

struct Measured {
    long double from_energy = 0.0L;
    long double from_force = 0.0L;
    long double from_virial = 0.0L;
};

// Recovers k_e from a provider result and a k_e = 1 reference. The force and
// virial channels use the largest-magnitude component so the division is never
// near a zero crossing.
Measured recover(const gmd::ForceResult& result, const Ewald1& reference) {
    Measured m;
    m.from_energy = static_cast<long double>(result.potential_energy) / reference.energy;

    std::size_t best_atom = 0, best_axis = 0;
    long double best = 0.0L;
    for (std::size_t i = 0; i < reference.forces.size(); ++i)
        for (std::size_t d = 0; d < 3; ++d)
            if (std::fabs(reference.forces[i][d]) > best) {
                best = std::fabs(reference.forces[i][d]);
                best_atom = i; best_axis = d;
            }
    m.from_force = static_cast<long double>(result.forces[best_atom][best_axis])
                 / reference.forces[best_atom][best_axis];

    std::size_t best_component = 0;
    best = 0.0L;
    for (std::size_t c = 0; c < 9; ++c)
        if (std::fabs(reference.virial[c]) > best) {
            best = std::fabs(reference.virial[c]); best_component = c;
        }
    m.from_virial = static_cast<long double>(result.virial[best_component])
                  / reference.virial[best_component];
    return m;
}

// Round-off budget. The reference sums thousands of k-vectors in long double
// (53-bit on Apple arm64) and the provider sums the same terms in double in a
// different order, so agreement is limited by cancellation in the sum, not by
// the constant. 1e-12 relative is roughly four orders above the observed
// spread and still six orders below the 3.16e-06 change under audit.
constexpr long double kRelativeTolerance = 1.0e-12L;

bool close(long double a, long double b, long double tolerance = kRelativeTolerance) {
    return std::fabs(a / b - 1.0L) <= tolerance;
}

void report(const std::string& label, const Measured& m) {
    std::cout << "    " << std::left << std::setw(34) << label
              << "  E " << number(m.from_energy)
              << "   F " << number(m.from_force)
              << "   W " << number(m.from_virial) << '\n';
}

Measured measure_ewald(bool neutral) {
    const Fixture f = asymmetric_fixture(neutral);
    gmd::System system = build_system(f);
    gmd::EwaldForceProvider provider(static_cast<double>(kAlpha), kKmax,
                                     static_cast<double>(kCutoff));
    const auto result = evaluate(provider, system);
    const Ewald1 reference =
        ewald_reference_unit_constant(f.positions, f.charges, f.lengths,
                                      kAlpha, kKmax, kCutoff);
    return recover(result, reference);
}

// PME differs from Ewald only in the reciprocal term, so its measurement needs
// a PME reciprocal reference. pme_reference.hpp supplies one; its own constant
// is divided straight back out, which is what keeps this independent of the
// value that header declares.
Measured measure_pme(bool neutral, int order, int grid) {
    const Fixture f = asymmetric_fixture(neutral);
    gmd::System system = build_system(f);
    gmd::PMEForceProvider provider(static_cast<double>(kAlpha),
                                   static_cast<double>(kCutoff), order,
                                   {grid, grid, grid});
    const auto result = evaluate(provider, system);

    // Real, self and background from the k_e = 1 reference; reciprocal from the
    // PME reference with its constant removed.
    Ewald1 reference =
        ewald_reference_unit_constant(f.positions, f.charges, f.lengths,
                                      kAlpha, /*kmax=*/0, kCutoff);
    std::vector<std::array<pme_ref::Real, 3>> pme_positions;
    std::vector<pme_ref::Real> pme_charges;
    for (std::size_t i = 0; i < f.charges.size(); ++i) {
        pme_positions.push_back({f.positions[i][0], f.positions[i][1], f.positions[i][2]});
        pme_charges.push_back(f.charges[i]);
    }
    const auto mesh = pme_ref::reciprocal_reference(
        pme_positions, pme_charges,
        {f.lengths[0], f.lengths[1], f.lengths[2]},
        kAlpha, std::array<int, 3>{grid, grid, grid}, order);
    reference.energy += mesh.energy_reciprocal_space / pme_ref::kCoulomb;
    for (std::size_t i = 0; i < f.charges.size(); ++i)
        for (std::size_t d = 0; d < 3; ++d)
            reference.forces[i][d] += mesh.forces[i][d] / pme_ref::kCoulomb;

    Measured m;
    m.from_energy = static_cast<long double>(result.potential_energy) / reference.energy;
    std::size_t best_atom = 0, best_axis = 0;
    long double best = 0.0L;
    for (std::size_t i = 0; i < reference.forces.size(); ++i)
        for (std::size_t d = 0; d < 3; ++d)
            if (std::fabs(reference.forces[i][d]) > best) {
                best = std::fabs(reference.forces[i][d]); best_atom = i; best_axis = d;
            }
    m.from_force = static_cast<long double>(result.forces[best_atom][best_axis])
                 / reference.forces[best_atom][best_axis];
    // The PME reciprocal virial is the analytic mesh tensor, not a sum over the
    // reference's k-vectors, so it is not recovered from this reference. The
    // Ewald channel above covers the virial, and the PME virial has its own
    // component-wise coverage in tests/pme_reciprocal_virial_tests.cpp.
    m.from_virial = m.from_energy;
    return m;
}

// The special-pair correction is a bare k_e q_a q_b (scale - 1) / r term with
// no lattice sum at all, so it measures the constant on the shortest possible
// path. It is isolated by differencing two runs that differ only in whether a
// special-pair map is attached.
Measured measure_special_pair() {
    const Fixture f = asymmetric_fixture(true);
    gmd::EwaldForceProvider provider(static_cast<double>(kAlpha), kKmax,
                                     static_cast<double>(kCutoff));

    gmd::System plain = build_system(f);
    const auto without = evaluate(provider, plain);

    gmd::System scaled = build_system(f);
    auto topology = std::make_shared<gmd::Topology>();
    topology->bonds.push_back({0, 2, 0});   // 1-2 pair: atoms 0 and 2
    gmd::SpecialPairScaleConfig scales;
    scales.pair_12 = {0.0, 0.0};            // full exclusion
    scaled.set_special_pair_map(
        std::make_shared<gmd::SpecialPairMap>(*topology, scales));
    const auto with = evaluate(provider, scaled);

    // Difference is exactly (0 - 1) * k_e * q0 * q2 / r02.
    std::array<long double, 3> dr{};
    for (std::size_t d = 0; d < 3; ++d) {
        dr[d] = f.positions[0][d] - f.positions[2][d];
        dr[d] -= f.lengths[d] * std::round(dr[d] / f.lengths[d]);
    }
    const long double r = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
    const long double phi = -f.charges[0] * f.charges[2] / r;
    const long double ff = -f.charges[0] * f.charges[2] / (r * r * r);

    Measured m;
    m.from_energy =
        (static_cast<long double>(with.potential_energy) -
         static_cast<long double>(without.potential_energy)) / phi;
    std::size_t axis = 0;
    for (std::size_t d = 1; d < 3; ++d)
        if (std::fabs(dr[d]) > std::fabs(dr[axis])) axis = d;
    m.from_force =
        (static_cast<long double>(with.forces[0][axis]) -
         static_cast<long double>(without.forces[0][axis])) / (ff * dr[axis]);
    std::size_t a = 0, b = 0;
    long double best = 0.0L;
    for (std::size_t p = 0; p < 3; ++p)
        for (std::size_t q = 0; q < 3; ++q)
            if (std::fabs(dr[p] * dr[q]) > best) { best = std::fabs(dr[p]*dr[q]); a = p; b = q; }
    m.from_virial =
        (static_cast<long double>(with.virial[a * 3 + b]) -
         static_cast<long double>(without.virial[a * 3 + b])) / (ff * dr[a] * dr[b]);
    return m;
}

// --- tests -----------------------------------------------------------------

std::vector<std::pair<std::string, Measured>> all_paths;

void test_every_path_measures_the_same_constant() {
    std::cout << "\n  Constant recovered from each electrostatic path\n";
    all_paths = {
        {"Ewald, neutral", measure_ewald(true)},
        {"Ewald, net-charged", measure_ewald(false)},
        {"PME order 4, neutral", measure_pme(true, 4, 32)},
        {"PME order 6, net-charged", measure_pme(false, 6, 32)},
        {"special-pair direct correction", measure_special_pair()},
    };
    for (const auto& [label, m] : all_paths) report(label, m);

    for (const auto& [label, m] : all_paths) {
        check(close(m.from_energy, kExpectedProduction),
              label + ": energy path measures k_e = " + number(m.from_energy) +
                  ", expected " + number(kExpectedProduction) +
                  " (relative deviation " +
                  number(m.from_energy / kExpectedProduction - 1.0L, 4) + ")");
        check(close(m.from_force, kExpectedProduction),
              label + ": FORCE path measures k_e = " + number(m.from_force) +
                  " but the energy path measures " + number(m.from_energy) +
                  "; the force is not the gradient of the energy the engine reports");
        check(close(m.from_virial, kExpectedProduction),
              label + ": VIRIAL path measures k_e = " + number(m.from_virial) +
                  " against energy " + number(m.from_energy));
    }
}

void test_no_squared_or_missing_factor() {
    // A path that dropped the constant would measure 1; one that applied it
    // twice would measure k_e^2. Both are far outside any round-off budget, so
    // they are called out by name rather than left to a generic tolerance.
    for (const auto& [label, m] : all_paths) {
        for (const auto& [channel, value] :
             {std::pair{"energy", m.from_energy}, std::pair{"force", m.from_force},
              std::pair{"virial", m.from_virial}}) {
            check(std::fabs(value - 1.0L) > 1.0L,
                  std::string(label) + " " + channel +
                      " measures 1, i.e. the Coulomb constant is missing from that path");
            check(std::fabs(value / (kExpectedProduction * kExpectedProduction) - 1.0L) > 1e-6L,
                  std::string(label) + " " + channel + " measures k_e squared (" +
                      number(value) + "), i.e. the constant is applied twice");
        }
    }
}

void test_paths_agree_with_each_other() {
    // Stronger than each matching the pin: no two production paths may carry
    // differently rounded copies of the constant.
    const long double base = all_paths.front().second.from_energy;
    for (const auto& [label, m] : all_paths) {
        check(close(m.from_energy, base),
              label + " measures k_e = " + number(m.from_energy) + " but " +
                  all_paths.front().first + " measures " + number(base) +
                  "; two production paths are using different constants");
    }
}

void test_linear_in_the_constant() {
    // k_e enters linearly, so scaling every charge by lambda must scale the
    // energy by lambda^2 exactly and leave the recovered constant unchanged.
    const Fixture f = asymmetric_fixture(true);
    const long double lambda = 3.0L;

    gmd::System plain = build_system(f);
    gmd::EwaldForceProvider provider(static_cast<double>(kAlpha), kKmax,
                                     static_cast<double>(kCutoff));
    const auto base = evaluate(provider, plain);

    gmd::System scaled = build_system(f);
    for (std::size_t i = 0; i < f.charges.size(); ++i)
        scaled.mutable_charges()[i] *= static_cast<double>(lambda);
    const auto boosted = evaluate(provider, scaled);

    const long double ratio =
        static_cast<long double>(boosted.potential_energy) /
        static_cast<long double>(base.potential_energy);
    check(close(ratio, lambda * lambda, 1.0e-13L),
          "scaling every charge by " + number(lambda, 3) + " scaled the energy by " +
              number(ratio) + ", expected exactly " + number(lambda * lambda, 3) +
              "; the constant is not entering linearly");
}

void test_sign_conventions() {
    // Two like charges must repel and two unlike charges must attract, with the
    // force on the first atom pointing along +dr for repulsion. A global sign
    // error in k_e would invert both and is not detectable from magnitudes.
    const std::array<long double, 3> lengths = {17.0L, 21.0L, 25.0L};
    for (const bool like : {true, false}) {
        Fixture f;
        f.lengths = lengths;
        f.positions = {{4.0L, 5.0L, 6.0L}, {7.5L, 5.0L, 6.0L}};
        f.charges = {0.8L, like ? 0.8L : -0.8L};
        gmd::System system = build_system(f);
        gmd::EwaldForceProvider provider(static_cast<double>(kAlpha), kKmax,
                                         static_cast<double>(kCutoff));
        const auto result = evaluate(provider, system);
        // dr = r0 - r1 points in -x, so a repulsive force on atom 0 is -x.
        const double fx = result.forces[0][0];
        check(like ? (fx < 0.0) : (fx > 0.0),
              std::string(like ? "like" : "unlike") +
                  " charges: force on atom 0 is x = " + number(fx) +
                  ", which has the wrong sign for " +
                  (like ? "repulsion" : "attraction"));
    }
}

void test_two_charge_energy_at_known_separation() {
    // The most direct statement of what the constant means: two unit charges
    // one Angstrom apart have Coulomb energy k_e. Ewald cannot be evaluated at
    // 1 A in a small box without the periodic images mattering, so this uses
    // the special-pair correction, which is the bare 1/r term.
    Fixture f;
    f.lengths = {40.0L, 44.0L, 48.0L};
    f.positions = {{10.0L, 12.0L, 14.0L}, {11.0L, 12.0L, 14.0L}};  // exactly 1 A apart
    f.charges = {1.0L, 1.0L};

    gmd::EwaldForceProvider provider(static_cast<double>(kAlpha), kKmax,
                                     static_cast<double>(kCutoff));
    gmd::System plain = build_system(f);
    const auto without = evaluate(provider, plain);

    gmd::System excluded = build_system(f);
    auto topology = std::make_shared<gmd::Topology>();
    topology->bonds.push_back({0, 1, 0});
    gmd::SpecialPairScaleConfig scales;
    scales.pair_12 = {0.0, 0.0};   // full exclusion
    excluded.set_special_pair_map(
        std::make_shared<gmd::SpecialPairMap>(*topology, scales));
    const auto with = evaluate(provider, excluded);

    const long double delta =
        static_cast<long double>(without.potential_energy) -
        static_cast<long double>(with.potential_energy);
    std::cout << "\n  Two unit charges 1 A apart: k_e = " << number(delta) << '\n';
    check(close(delta, kExpectedProduction),
          "two unit charges 1 Angstrom apart give " + number(delta) +
              " eV; k_e is " + number(kExpectedProduction));

    // And the force magnitude, which must be k_e / r^2 = k_e at r = 1 A.
    const long double force =
        static_cast<long double>(without.forces[0][0]) -
        static_cast<long double>(with.forces[0][0]);
    check(close(std::fabs(force), kExpectedProduction),
          "force magnitude between two unit charges 1 Angstrom apart is " +
              number(std::fabs(force)) + " eV/A; k_e / r^2 is " +
              number(kExpectedProduction));
    check(force < 0.0L,
          "force on atom 0 from a like charge at +x must point in -x, got " +
              number(force));

    // The virial of a single 1/r pair: W_ab = r_a F_b, so W_xx = -k_e here
    // (dr = r0 - r1 = -1 A in x, F_x on atom 0 = -k_e).
    const long double wxx =
        static_cast<long double>(without.virial[0]) -
        static_cast<long double>(with.virial[0]);
    check(close(wxx, kExpectedProduction),
          "virial W_xx of two unit charges 1 Angstrom apart is " + number(wxx) +
              ", expected k_e = " + number(kExpectedProduction));
    for (const std::size_t off : {1u, 2u, 3u, 5u, 6u, 7u}) {
        const long double value =
            static_cast<long double>(without.virial[off]) -
            static_cast<long double>(with.virial[off]);
        check(std::fabs(value) < 1.0e-12L,
              "virial component " + std::to_string(off) +
                  " of an x-aligned pair must vanish, got " + number(value));
    }
}

void test_deviation_from_codata_is_as_documented() {
    const long double relative = kExpectedProduction / kCodata2022 - 1.0L;
    std::cout << "\n  CODATA 2022  E_h * a0 = " << number(kCodata2022, 17) << '\n';
    std::cout << "  production            = " << number(kExpectedProduction, 17)
              << "   relative deviation " << number(relative, 4) << '\n';
    check(std::fabs(relative / kKnownDeviationFromCodata - 1.0L) <= kDeviationTolerance,
          "the production constant's deviation from the CODATA 2022 derivation is " +
              number(relative, 6) + ", not the documented " +
              number(kKnownDeviationFromCodata, 6) +
              ". If the constant was corrected, update kExpectedProduction and "
              "replace this characterization with a rounding-policy assertion.");
}

}  // namespace

int main() {
    std::cout << "Electrostatic unit-conversion constant audit\n";
    test_deviation_from_codata_is_as_documented();
    test_two_charge_energy_at_known_separation();
    test_every_path_measures_the_same_constant();
    test_no_squared_or_missing_factor();
    test_paths_agree_with_each_other();
    test_linear_in_the_constant();
    test_sign_conventions();

    if (failures != 0) {
        std::cerr << "\nElectrostatic constant audit failed: " << failures << '\n';
        return 1;
    }
    std::cout << "\nElectrostatic constant audit passed\n";
    return 0;
}
