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

}  // namespace gmd
