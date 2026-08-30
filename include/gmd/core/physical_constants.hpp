#pragma once

// Physical constants in GMD's internal unit system: eV for energy, Angstrom
// for distance, elementary charge for charge, and femtosecond-derived time.
//
// One definition per constant. Before this header the Coulomb constant existed
// as two independent literals, one in ewald_force_provider.cpp and one in
// pme_force_provider.cpp, with nothing keeping them equal.

namespace gmd {

// Coulomb constant k_e = e^2 / (4 pi eps0), in eV * Angstrom / e^2.
//
// Numerically it is the Coulomb energy, in eV, of two unit charges one
// Angstrom apart -- which is exactly how tests/electrostatic_constant_tests.cpp
// measures it back out of the engine.
//
// DERIVATION. In atomic units the Coulomb energy of two unit charges separated
// by one Bohr radius is exactly one Hartree, so
//
//     e^2 / (4 pi eps0) = E_h * a0
//
// and therefore, in GMD's units,
//
//     k_e [eV A / e^2] = E_h [eV] * a0 [A]
//                      = 27.211386245981 * 0.529177210544
//                      = 14.399645468683593...
//
// Both factors are 2022 CODATA recommended values, tabulated by NIST:
//
//     E_h = 27.211 386 245 981(30) eV
//           https://physics.nist.gov/cgi-bin/cuu/Value?hrev
//     a0  = 5.291 772 105 44(82) x 10^-11 m
//           https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
//
// The same value follows from the elementary charge and the vacuum permittivity
// directly, k_e = e / (4 pi eps0 * 1e-10) with e = 1.602176634e-19 C (exact by
// the SI definition) and eps0 = 8.8541878188(14)e-12 F/m; the two routes agree
// to 1.1e-12 relative, which is the uncertainty in eps0 rather than an error in
// either.
//
// ROUNDING POLICY. Enough digits that the literal reproduces the derivation to
// double-precision representation (4.6e-16 relative). This is not excess
// precision: CODATA's own relative uncertainty here is 1.5e-10, and the
// tightest PME convergence point this repository measures is 6.3e-10 relative,
// so a constant rounded to, say, LAMMPS' six decimals (3.3e-08 relative) would
// be the coarsest quantity in that comparison and would appear as a spurious
// convergence floor. Digits beyond the tenth decimal are below the CODATA
// uncertainty and are carried only so the literal is exact as a double.
//
// COMPATIBILITY. This is the CODATA-derived value, not a match to any external
// engine. Nothing in this repository documents a requirement to reproduce
// another code's constant, and the engines do not agree with each other in any
// case: LAMMPS `units metal` uses 14.399645 (its own six-decimal rounding,
// 3.3e-08 low) and OpenMM's ONE_4PI_EPS0 works out to 14.399645478 (6.8e-10
// high). validation/pme_external converts each engine's result by the exactly
// linear ratio k_e^GMD / k_e^engine rather than adopting either value.
inline constexpr double kCoulombConstant = 14.3996454686836;

// Boltzmann constant, in eV / K.
//
// Numerically it is the kinetic energy in eV that corresponds to one kelvin per
// degree of freedom -- which is exactly how tests/boltzmann_constant_tests.cpp
// measures it back out of the engine, from 2K = dof * k_B * T.
//
// DERIVATION. Unlike most physical constants this one is EXACT and carries no
// uncertainty, because both of its ingredients have been exact by definition
// since the 2019 SI redefinition:
//
//     k_B = 1.380649e-23 J/K          https://physics.nist.gov/cgi-bin/cuu/Value?k
//     e   = 1.602176634e-19 C         https://physics.nist.gov/cgi-bin/cuu/Value?e
//
// One electronvolt is e joules exactly, so
//
//     k_B [eV/K] = 1.380649e-23 / 1.602176634e-19
//                = 1380649 / 16021766340
//                = 8.6173332621451774336636593340806392...e-5
//
// The denominator has factors other than 2 and 5, so this exact rational has a
// non-terminating decimal expansion: every literal is a truncation of it. NIST
// tabulates the value as "8.617 333 262... x 10^-5 eV K^-1, exact", where the
// ellipsis is their notation for precisely that.
// (https://physics.nist.gov/cgi-bin/cuu/Value?kev, CODATA 2022.)
//
// ROUNDING POLICY. Since there is no physical uncertainty to hide behind, the
// only sensible cut-off is double precision. The literal below is the shortest
// decimal that parses to the nearest double to the exact rational; it sits
// 1.1e-18 from the true value, which is below one ulp. Writing fewer digits
// would be an arbitrary truncation with nothing to justify it: the previous
// 8.617333262e-5 was 1.7e-11 low and 8.617343e-5 was 1.13e-06 low.
inline constexpr double kBoltzmannConstantEVPerKelvin = 8.617333262145177e-5;

// Pressure conversion between bar and GMD's internal unit.
//
// GMD computes pressure as P = (2K + tr W) / 3V. With energies in eV and
// lengths in Angstrom that is an energy density in eV/A^3, and it is the unit
// every internal pressure is stored in: System::StepThermodynamics::pressure,
// the checkpoint's step_pressure, and the value each barostat compares against.
// bar is an interface unit only -- it is what a run input asks for and what the
// P[bar] log column reports -- so exactly one conversion stands between them.
//
// DERIVATION. Unlike k_e and k_B this is not a measured quantity at all. It is a
// pure unit identity, and all four of its ingredients are exact by definition:
//
//     1 bar = 100000 Pa            (definition of the bar)
//     1 Pa  = 1 J/m^3              (definition of the pascal)
//     1 A   = 1e-10 m              (definition of the angstrom)
//     1 eV  = 1.602176634e-19 J    (exact since the 2019 SI redefinition;
//                                   https://physics.nist.gov/cgi-bin/cuu/Value?e
//                                   and https://physics.nist.gov/cgi-bin/cuu/Value?evj)
//
// so
//
//     1 bar = 1e5 J/m^3
//           = 1e5 * 1e-30 J/A^3
//           = 1e-25 / 1.602176634e-19  eV/A^3
//           = 1e-6 / 1.602176634       eV/A^3
//           = 1000 / 1602176634  =  500 / 801088317  eV/A^3
//           = 6.24150907446076260777624098...e-7 eV/A^3
//
// WHICH DIRECTION IS PRIMARY. The reverse conversion is the reciprocal,
// 801088317 / 500, and unlike the forward one it TERMINATES:
//
//     1 eV/A^3 = 1602176.634 bar, exactly.
//
// So that direction is written as the literal and the forward one is derived
// from it. This is not cosmetic. Writing the forward direction as a literal and
// the reverse as its reciprocal, or computing either through the SI constants,
// rounds more than once and lands one ulp off; taken from the terminating
// decimal instead, each direction is the nearest double to its exact rational
// AND the two are exact reciprocals of each other in double arithmetic. The
// static_asserts below hold the pair to both properties, so a conversion applied
// forwards and then backwards returns the original bits.
//
// PRECISION. The old value, 6.2415091e-7, was a rounding of the forward
// direction to eight significant figures and sat 4.091837e-09 relative above the
// exact value. Since nothing here is uncertain there is no measurement precision
// to round to, and the only defensible cut-off is the double itself.
inline constexpr double kEVPerAngstromCubedToBar = 1602176.634;
inline constexpr double kBarToEVPerAngstromCubed = 1.0 / kEVPerAngstromCubedToBar;

static_assert(kEVPerAngstromCubedToBar * kBarToEVPerAngstromCubed == 1.0,
              "the two pressure conversion directions must be exact reciprocals");
static_assert(kBarToEVPerAngstromCubed == 500.0 / 801088317.0,
              "bar -> eV/A^3 must be the nearest double to the exact rational "
              "500/801088317");
static_assert(kEVPerAngstromCubedToBar == 801088317.0 / 500.0,
              "eV/A^3 -> bar must be the nearest double to the exact rational "
              "801088317/500");

// Internal time unit, in femtoseconds.
//
// GMD integrates v += (F/m)*dt and r += v*dt with F in eV/A and m in amu.
// Neither line mentions seconds, so the time unit is not a free choice: it is
// whatever makes F/m an acceleration in this unit system. F/m carries units of
// eV/(A*amu), and requiring that to equal A/T^2 gives
//
//     T = A * sqrt(amu / eV)
//
// One internal velocity unit is correspondingly sqrt(eV/amu), which is
// 0.0982269474... A/fs -- the number a free particle actually travels per
// femtosecond, and how tests/time_unit_tests.cpp measures this back out of the
// integrator.
//
// DERIVATION. Two of the ingredients are exact and one is not:
//
//     e   = 1.602176634e-19 J per eV     exact, SI 2019
//           https://physics.nist.gov/cgi-bin/cuu/Value?evj
//     m_u = 1.66053906892(52)e-27 kg     CODATA 2022, relative 3.1e-10
//           https://physics.nist.gov/cgi-bin/cuu/Value?ukg
//     1 A  = 1e-10 m                     exact, by definition
//     1 fs = 1e-15 s                     exact, by definition
//
// so
//
//     T [fs] = 1e-10 * sqrt(m_u / e) / 1e-15 = sqrt(m_u / e) * 1e5
//            = 10.1805057178711931077510010336...
//
// The square root halves the mass constant's uncertainty, so T carries 1.57e-10
// relative -- the only constant in this header that has a real uncertainty at
// all.
//
// ROUNDING POLICY. Enough digits to reproduce the derivation as a double. That
// is well below the 1.57e-10 physical uncertainty and is not a claim of
// physical precision; it is so that the literal does not ADD error to a
// quantity that already has some. The superseded value, 1.018051e+1, was this
// rounded to seven significant figures and sat 4.206204e-07 relative high --
// about 2700 times the CODATA uncertainty, so it was not a defensible
// truncation of it.
//
// NAMING. The superseded constant was called kInternalTimeUnitsPerFs but was
// used as a divisor of a femtosecond timestep, which makes it femtoseconds per
// internal time unit -- the reciprocal of what its name said. Both directions
// are given here so that neither call site has to divide by a constant whose
// name reads the wrong way round.
inline constexpr double kFemtosecondsPerInternalTime = 10.180505717871194;
inline constexpr double kInternalTimePerFemtosecond = 1.0 / kFemtosecondsPerInternalTime;

static_assert(kFemtosecondsPerInternalTime * kInternalTimePerFemtosecond == 1.0,
              "the two internal-time conversion directions must be exact "
              "reciprocals");

}  // namespace gmd
