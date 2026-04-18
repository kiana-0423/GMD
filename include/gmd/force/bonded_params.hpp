#pragma once

namespace gmd {

// All energy units follow the project convention: eV / Å / amu.
// Use the helpers in bonded_force_provider.hpp (kcal_per_mol_to_eV, deg_to_rad)
// to convert from typical AMBER/CHARMM/OPLS input units.

// Harmonic bond:  V = k (r - r0)^2
struct BondParams {
    double k  = 0.0;   // force constant [eV/Å²]
    double r0 = 0.0;   // equilibrium bond length [Å]
};

// Harmonic angle:  V = k (θ - θ0)^2    (θ in radians internally)
struct AngleParams {
    double k      = 0.0;  // force constant [eV/rad²]
    double theta0 = 0.0;  // equilibrium angle [rad]
};

// Periodic proper dihedral:  V = k [1 + cos(n φ - δ)]
// Multiple entries with the same atom quad are summed (Fourier series).
struct DihedralParams {
    double k     = 0.0;  // barrier height [eV]
    int    n     = 1;    // periodicity (multiplicity ≥ 1)
    double delta = 0.0;  // phase shift [rad]
};

// Harmonic improper dihedral:  V = k (φ - φ0)^2
struct ImproperParams {
    double k    = 0.0;   // force constant [eV/rad²]
    double phi0 = 0.0;   // equilibrium angle [rad]
};

}  // namespace gmd
