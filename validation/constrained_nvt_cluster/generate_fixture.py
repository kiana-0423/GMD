#!/usr/bin/env python3
"""Regenerate the constrained NVT cluster fixture deterministically.

Run from the repository root:

    python3 validation/constrained_nvt_cluster/generate_fixture.py

Rewrites cluster.xyz and cluster.top in place. No randomness is involved: the
molecular placements and velocities below are fixed literals, so the fixture is
reproducible from this script alone. See README.md for the provenance record.

PHYSICAL INTENT
---------------
FOUR rigid-water-like molecules that DO interact. Unlike the single molecule of
the constrained_nve_water case, the intermolecular Lennard-Jones terms here are
real, so the force providers produce a non-zero virial and the constraint virial
has a genuine partner to be added to. Each molecule carries three distance
constraints (two O-H and the H-H that fixes the angle), giving twelve
constraints over four connected components of the constraint graph.

DELIBERATE PROPERTIES
---------------------
  * IRREGULAR placement and orientation. Every molecule has its own generic
    rotation and its own centre, none of them related by symmetry, so no error
    can cancel between molecules. The velocities are likewise unequal.
  * MPI RANK BOUNDARIES. GMD splits the cell over a balanced process grid:
    (2,1,1) at np=2 and (2,2,1) at np=4. Molecule C is placed astride
    x = Lx/2 and molecule D astride y = Ly/2, so at both rank counts at least
    one constraint has its two atoms owned by different ranks.
  * AN EMPTY RANK AT np=4. The (2,2,1) grid divides the cell into four
    quadrants in x-y. All four molecules are placed so that the
    (x > Lx/2, y > Ly/2) quadrant contains no atom at all, so one of the four
    ranks owns nothing and the collective paths are exercised with an empty
    local domain.
  * Molecules are >= 3 A apart, well inside the 8 A cutoff, so they interact
    without the initial configuration being a hard clash that would blow up.
"""
from __future__ import annotations

import math
import pathlib

R_OH = 0.9572
ANGLE_HOH_DEG = 104.52
D_HH = 2.0 * R_OH * math.sin(math.radians(ANGLE_HOH_DEG) / 2.0)

BOX = (20.0, 22.0, 24.0)
MASSES = (15.9994, 1.0080, 1.0080)
TYPES = (1, 2, 2)

# Each molecule starts this far off its target distances, so the initial SHAKE
# projection does real work on every component. The values differ per molecule
# so the four components are not in the same state.
SLACK = (0.0018, -0.0022, 0.0031, -0.0014)

# centre, rotation (rad), centre-of-mass velocity [A/fs], angular velocity [rad/fs]
#
# C straddles x = 10.0 and D straddles y = 11.0; the quadrant x>10, y>11 is left
# empty so that np=4 has a rank owning no atoms.
MOLECULES = [
    # A: lower-left quadrant, well inside its own domain.
    ((4.7310, 5.2170, 7.3140), (0.6109, 0.4014, 0.9599),
     (0.0031, -0.0024, 0.0017), (0.0022, 0.0036, -0.0029)),
    # B: lower-left quadrant, a different corner of it.
    ((7.1930, 3.8410, 15.6270), (1.2217, 0.8727, 0.3491),
     (-0.0027, 0.0019, -0.0033), (-0.0031, 0.0017, 0.0024)),
    # C: astride x = Lx/2 = 10.0  -> cross-rank at np=2 and np=4.
    ((10.0000, 4.6120, 11.2830), (0.3491, 1.0472, 1.3963),
     (0.0022, 0.0029, -0.0015), (0.0018, -0.0026, 0.0033)),
    # D: astride y = Ly/2 = 11.0, kept at x < Lx/2 so the x>10,y>11 quadrant
    #    stays empty -> cross-rank at np=4.
    ((5.8470, 11.0000, 19.4160), (0.9599, 0.1745, 0.7854),
     (-0.0019, -0.0026, 0.0028), (0.0027, 0.0021, -0.0018)),
]


def rotate(point, angles):
    a, b, c = angles
    x, y, z = point
    x, y = x * math.cos(a) - y * math.sin(a), x * math.sin(a) + y * math.cos(a)
    y, z = y * math.cos(b) - z * math.sin(b), y * math.sin(b) + z * math.cos(b)
    z, x = z * math.cos(c) - x * math.sin(c), z * math.sin(c) + x * math.cos(c)
    return (x, y, z)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def molecule(center, angles, v_com, omega, slack):
    half = math.radians(ANGLE_HOH_DEG) / 2.0
    local = [
        (0.0, 0.0, 0.0),
        (R_OH * math.cos(half), R_OH * math.sin(half), 0.0),
        (R_OH * math.cos(half), -R_OH * math.sin(half), 0.0),
    ]
    local = [rotate(p, angles) for p in local]
    local = [tuple((1.0 + slack) * v for v in p) for p in local]
    total = sum(MASSES)
    com = tuple(sum(MASSES[i] * local[i][d] for i in range(3)) / total for d in range(3))
    offsets = [tuple(local[i][d] - com[d] for d in range(3)) for i in range(3)]

    rows = []
    for i in range(3):
        position = tuple(center[d] + offsets[i][d] for d in range(3))
        spin = cross(omega, offsets[i])
        velocity = tuple(v_com[d] + spin[d] for d in range(3))
        rows.append((TYPES[i], position, velocity))
    return rows


def main() -> int:
    here = pathlib.Path(__file__).resolve().parent

    rows = []
    for (center, angles, v_com, omega), slack in zip(MOLECULES, SLACK):
        rows.extend(molecule(center, angles, v_com, omega, slack))

    lines = [str(len(rows)), f"{BOX[0]} {BOX[1]} {BOX[2]}"]
    for atom_type, position, velocity in rows:
        lines.append(
            f"{atom_type}  "
            f"{position[0]:.14f} {position[1]:.14f} {position[2]:.14f}  "
            f"{velocity[0]:.14f} {velocity[1]:.14f} {velocity[2]:.14f}"
        )
    (here / "cluster.xyz").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Topology: per molecule, two exclusion-only bonds, one exclusion-only angle,
    # and three distance constraints.
    bonds, angles_out, constraints = [], [], []
    for m in range(len(MOLECULES)):
        o, h1, h2 = 3 * m + 1, 3 * m + 2, 3 * m + 3      # 1-based
        bonds += [f"  {o} {h1}  bond_type 1", f"  {o} {h2}  bond_type 1"]
        angles_out.append(f"  {h1} {o} {h2}  angle_type 1")
        constraints += [
            f"  {o} {h1}  {R_OH:.14f}",
            f"  {o} {h2}  {R_OH:.14f}",
            f"  {h1} {h2}  {D_HH:.14f}",
        ]

    top = [
        "# Four rigid-water-like molecules; see generate_fixture.py.",
        "#",
        "# The bonds and angles carry NO force (their types have k = 0 in",
        "# cluster.ff). They exist to exclude the intramolecular non-bonded pairs,",
        "# so the only real forces are INTERMOLECULAR. The rigid geometry comes",
        "# from the constraints section.",
        "",
        f"bonds {len(bonds)}", *bonds, "",
        f"angles {len(angles_out)}", *angles_out, "",
        "dihedrals 0",
        "impropers 0",
        "",
        f"constraints {len(constraints)}", *constraints, "",
    ]
    (here / "cluster.top").write_text("\n".join(top) + "\n", encoding="utf-8")

    # --- report, so the deliberate properties are visible when regenerating ---
    print(f"atoms={len(rows)}  constraints={len(constraints)}  components={len(MOLECULES)}")
    print(f"d_HH = {D_HH:.14f} A")
    for m, ((center, _a, _v, _o), slack) in enumerate(zip(MOLECULES, SLACK)):
        tag = "ABCD"[m]
        sub = rows[3 * m:3 * m + 3]
        d01 = math.dist(sub[0][1], sub[1][1])
        print(f"  molecule {tag}: center={center} slack={slack*100:+.2f}%  |O-H1|={d01:.10f}")
    xs = [r[1][0] for r in rows]
    ys = [r[1][1] for r in rows]
    print(f"  x range {min(xs):.3f}..{max(xs):.3f}   split at {BOX[0]/2}")
    print(f"  y range {min(ys):.3f}..{max(ys):.3f}   split at {BOX[1]/2}")
    for m in range(len(MOLECULES)):
        sub = rows[3 * m:3 * m + 3]
        sx = [1 if r[1][0] > BOX[0] / 2 else 0 for r in sub]
        sy = [1 if r[1][1] > BOX[1] / 2 else 0 for r in sub]
        note = []
        if len(set(sx)) > 1: note.append("crosses x=Lx/2")
        if len(set(sy)) > 1: note.append("crosses y=Ly/2")
        print(f"  molecule {'ABCD'[m]}: {', '.join(note) if note else 'inside one domain'}")
    occupied = {(1 if r[1][0] > BOX[0]/2 else 0, 1 if r[1][1] > BOX[1]/2 else 0) for r in rows}
    print(f"  occupied (x,y) quadrants at np=4: {sorted(occupied)}")
    print(f"  empty quadrants: {sorted({(0,0),(0,1),(1,0),(1,1)} - occupied)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
