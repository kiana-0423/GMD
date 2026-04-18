#pragma once

#include <vector>

namespace gmd {

// ---------------------------------------------------------------------------
// Per-interaction topology entries
// Each entry stores 1-based atom indices (converted to 0-based internally) and
// an index into the corresponding parameter-type table held by
// BondedForceProvider.
// ---------------------------------------------------------------------------

// Harmonic bond:  V = k (r - r0)^2
struct BondTerm {
    int i, j;       // atom indices (0-based)
    int type_idx;   // index into BondedForceProvider::bond_types_
};

// Harmonic angle:  V = k (θ - θ0)^2   (j is the vertex atom)
struct AngleTerm {
    int i, j, k;    // atom indices (0-based), j = vertex
    int type_idx;
};

// Periodic dihedral:  V = k [1 + cos(n φ - δ)]
// i-j-k-l torsion angle
struct DihedralTerm {
    int i, j, k, l; // atom indices (0-based)
    int type_idx;
};

// Harmonic improper:  V = k (φ - φ0)^2
// Uses the same dihedral-angle geometry (i-j-k-l) but a harmonic well.
struct ImproperTerm {
    int i, j, k, l; // atom indices (0-based)
    int type_idx;
};

// Holds the full bonded connectivity for a molecular system.
struct Topology {
    std::vector<BondTerm>     bonds;
    std::vector<AngleTerm>    angles;
    std::vector<DihedralTerm> dihedrals;
    std::vector<ImproperTerm> impropers;
};

}  // namespace gmd
