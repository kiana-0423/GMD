#!/usr/bin/env python3
"""Regenerate the constrained NVE water fixtures deterministically.

Run from the repository root:

    python3 validation/constrained_nve_water/generate_fixture.py

This rewrites water.xyz and water_wrapped.xyz in place. It takes no random
input: every number below is either a documented physical constant of the rigid
water geometry or a fixed literal, so the fixture is reproducible from this
script alone. See README.md for the provenance record.

PHYSICAL INTENT
---------------
One rigid-water-like molecule (O, H, H) with THREE distance constraints: the two
O-H bonds and the H-H distance that is equivalent to fixing the H-O-H angle.
Because the molecule is alone in the cell and every intramolecular pair is
excluded (1-2 by the bonds, 1-3 by the angle) with zero force constants, the
molecule feels NO force at all. It is therefore a free rigid body, whose motion
is analytically known: the centre of mass travels in a straight line at constant
speed, and the total energy, the linear momentum and the angular momentum about
the centre of mass are all exactly conserved. That is what makes this fixture
worth checking against rather than only against a previous GMD run.

DELIBERATE PROPERTIES
---------------------
  * unequal masses (O is ~16x H), so no error cancels by symmetry;
  * a generic rotation, so no bond is parallel to a box axis and every one of
    the nine virial components is non-trivial;
  * a non-lattice translation, so no coordinate is a round number;
  * INITIAL SLACK: every distance starts 0.25% away from its target, so the
    initial SHAKE projection has real work to do and "distances match their
    targets after projection" is a meaningful assertion rather than a tautology;
  * a NON-TANGENT initial velocity component, so the initial RATTLE projection
    also has real work to do. The projected state is a pure rigid-body motion,
    which is what the analytic conservation laws above then apply to.
"""
from __future__ import annotations

import math
import pathlib

# --- Rigid water geometry (TIP3P-like) -------------------------------------
R_OH = 0.9572                                      # [A]
ANGLE_HOH_DEG = 104.52                             # [deg]
D_HH = 2.0 * R_OH * math.sin(math.radians(ANGLE_HOH_DEG) / 2.0)

BOX = (24.0, 26.0, 28.0)

# Every distance starts this far ABOVE its target, so SHAKE must contract the
# molecule onto the manifold before the first step.
INITIAL_SLACK = 0.0025

# Rigid-body motion of the projected state.
V_COM = (0.0042, -0.0031, 0.0019)                  # [A/fs]
OMEGA = (0.0021, 0.0035, -0.0028)                  # [rad/fs]
# A radial (non-tangent) perturbation the initial RATTLE has to remove. It is
# applied along each atom's offset from the centre of mass, which is exactly the
# direction a rigid-body velocity field cannot contain.
RADIAL_PERTURBATION = 0.0013                       # [1/fs]

# Generic Euler-like rotation, chosen so no molecular axis lines up with a box
# axis. Radians, ~35 / 23 / 55 degrees.
ROTATION = (0.6109, 0.4014, 0.9599)

MASSES = (15.9994, 1.0080, 1.0080)
TYPES = (1, 2, 2)

# Placement of the centre of mass.
#
# The wrapped variant has to straddle a periodic face AT THE STEP THE CHECKS
# LOOK AT, which is the end of the run, not merely at the start. The molecule
# carries a centre-of-mass velocity, so a variant placed on the face at t=0
# simply drifts off it and the minimum-image path stops being exercised at all
# -- the fixture would still pass while testing nothing. Its x is therefore set
# so that the molecule arrives ON the x = 0 face after the 400 steps of nve.run.
# analyze.py ASSERTS that it is genuinely straddling there, so this cannot
# silently rot if the run length or the velocities change.
COM_IN_BOX = (9.3137, 11.6411, 13.2290)
COM_ON_FACE = (22.6566, 11.6411, 13.2290)


def rotate(point: tuple[float, float, float]) -> tuple[float, float, float]:
    a, b, c = ROTATION
    x, y, z = point
    x, y = x * math.cos(a) - y * math.sin(a), x * math.sin(a) + y * math.cos(a)
    y, z = y * math.cos(b) - z * math.sin(b), y * math.sin(b) + z * math.cos(b)
    z, x = z * math.cos(c) - x * math.sin(c), z * math.sin(c) + x * math.cos(c)
    return (x, y, z)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def molecule_frame() -> list[tuple[float, float, float]]:
    """Atom offsets from the centre of mass, rotated and with the slack applied."""
    half = math.radians(ANGLE_HOH_DEG) / 2.0
    local = [
        (0.0, 0.0, 0.0),
        (R_OH * math.cos(half), R_OH * math.sin(half), 0.0),
        (R_OH * math.cos(half), -R_OH * math.sin(half), 0.0),
    ]
    local = [rotate(p) for p in local]
    # Scaling the whole molecule scales every distance by the same factor, so a
    # single number puts all three constraints equally far off target.
    scale = 1.0 + INITIAL_SLACK
    local = [tuple(scale * v for v in p) for p in local]
    total = sum(MASSES)
    com = tuple(sum(MASSES[i] * local[i][d] for i in range(3)) / total for d in range(3))
    return [tuple(local[i][d] - com[d] for d in range(3)) for i in range(3)]


def build(com: tuple[float, float, float]) -> list[tuple]:
    offsets = molecule_frame()
    rows = []
    for i in range(3):
        position = tuple(com[d] + offsets[i][d] for d in range(3))
        spin = cross(OMEGA, offsets[i])
        velocity = tuple(
            V_COM[d] + spin[d] + RADIAL_PERTURBATION * offsets[i][d] for d in range(3)
        )
        rows.append((TYPES[i], position, velocity))
    return rows


def write_xyz(path: pathlib.Path, rows: list[tuple], comment: str) -> None:
    lines = [str(len(rows)), comment]
    for atom_type, position, velocity in rows:
        lines.append(
            f"{atom_type}  "
            f"{position[0]:.14f} {position[1]:.14f} {position[2]:.14f}  "
            f"{velocity[0]:.14f} {velocity[1]:.14f} {velocity[2]:.14f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    here = pathlib.Path(__file__).resolve().parent
    box = f"{BOX[0]} {BOX[1]} {BOX[2]}"
    write_xyz(here / "water.xyz", build(COM_IN_BOX), box)
    write_xyz(here / "water_wrapped.xyz", build(COM_ON_FACE), box)
    print(f"d_HH target = {D_HH:.14f} A")
    print(f"initial slack = {INITIAL_SLACK * 100:.4f}% on every constraint")
    for name, com in (("water.xyz", COM_IN_BOX), ("water_wrapped.xyz", COM_ON_FACE)):
        rows = build(com)
        print(f"{name}: com={com}")
        for (i, j, target) in ((0, 1, R_OH), (0, 2, R_OH), (1, 2, D_HH)):
            d = math.dist(rows[i][1], rows[j][1])
            print(f"    |{i}-{j}| = {d:.14f}  target {target:.14f}  "
                  f"slack {(d - target) / target * 100:+.4f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
