"""Evaluate the pme_external fixture with OpenMM PME and print JSON on stdout.

Run with an interpreter that has OpenMM installed. It is deliberately a
separate script rather than part of generate_external_reference.py so that the
OpenMM installation can live in its own environment, outside this repository,
and so that this file records exactly and only what was asked of OpenMM.

Everything is reported in OpenMM's own native units (kJ/mol, nm). No conversion
to GMD units happens here: the conversion is a separate, auditable step in
generate_external_reference.py, and mixing the two would make it impossible to
tell an algorithm difference from a units mistake.
"""

from __future__ import annotations

import argparse
import json
import sys

import openmm as mm
import openmm.unit as u
from openmm import Vec3


def read_fixture(path: str) -> tuple[list[float], list[list[float]], list[float]]:
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line for line in handle.read().splitlines() if line.strip()]
    count = int(lines[0].split()[0])
    box = [float(value) for value in lines[1].split()[:3]]
    positions: list[list[float]] = []
    charges: list[float] = []
    for line in lines[2 : 2 + count]:
        fields = line.split()
        positions.append([float(fields[1]), float(fields[2]), float(fields[3])])
        charges.append(float(fields[5]))
    return box, positions, charges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xyz", required=True)
    parser.add_argument("--alpha", type=float, required=True,
                        help="Ewald splitting parameter in 1/Angstrom")
    parser.add_argument("--cutoff", type=float, required=True,
                        help="real-space cutoff in Angstrom")
    parser.add_argument("--grid", type=int, nargs=3, required=True)
    parser.add_argument("--platform", default="Reference")
    args = parser.parse_args()

    box, positions, charges = read_fixture(args.xyz)

    system = mm.System()
    # Box vectors are the fixture box converted to nm. The fixture is
    # orthorhombic, which is all GMD's Box can express, so the off-diagonal
    # entries are exactly zero rather than approximately zero.
    system.setDefaultPeriodicBoxVectors(
        Vec3(box[0] / 10.0, 0.0, 0.0) * u.nanometer,
        Vec3(0.0, box[1] / 10.0, 0.0) * u.nanometer,
        Vec3(0.0, 0.0, box[2] / 10.0) * u.nanometer,
    )

    force = mm.NonbondedForce()
    force.setNonbondedMethod(mm.NonbondedForce.PME)
    force.setCutoffDistance(args.cutoff / 10.0 * u.nanometer)
    # sigma and epsilon are zero for every particle: this fixture is Coulomb
    # only, so there is no Lennard-Jones shift, switch, mixing rule or
    # long-range dispersion convention to reconcile between the two codes.
    force.setUseDispersionCorrection(False)
    force.setUseSwitchingFunction(False)
    # setPMEParameters takes alpha in 1/nm. Passing a nonzero alpha and an
    # explicit grid disables OpenMM's own error-estimate-driven selection, so
    # the splitting and the mesh are the ones asked for and not ones OpenMM
    # chose. The values actually in force are read back below and reported.
    force.setPMEParameters(args.alpha * 10.0, args.grid[0], args.grid[1], args.grid[2])
    for charge in charges:
        force.addParticle(charge, 0.0, 0.0)
    # No exclusions and no exceptions: every pair interacts at full strength,
    # matching a GMD system with no topology and therefore no special pairs.
    system.addForce(force)

    for _ in charges:
        system.addParticle(1.0 * u.dalton)

    platform = mm.Platform.getPlatformByName(args.platform)
    context = mm.Context(system, mm.VerletIntegrator(1.0 * u.femtosecond), platform)
    context.setPositions([Vec3(*[c / 10.0 for c in p]) for p in positions] * u.nanometer)

    state = context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    forces = [
        [component for component in vector.value_in_unit(u.kilojoule_per_mole / u.nanometer)]
        for vector in state.getForces()
    ]

    # What OpenMM actually used, which is not necessarily what was requested.
    used_alpha, used_nx, used_ny, used_nz = force.getPMEParametersInContext(context)

    json.dump(
        {
            "engine": "OpenMM",
            "version": mm.version.version,
            "git_revision": mm.version.git_revision,
            "platform": platform.getName(),
            "platform_precision": "double (Reference platform is double precision)",
            "energy_kj_per_mol": energy,
            "forces_kj_per_mol_per_nm": forces,
            "requested": {
                "alpha_inv_angstrom": args.alpha,
                "cutoff_angstrom": args.cutoff,
                "grid": list(args.grid),
            },
            "used": {
                "alpha_inv_nm": used_alpha.value_in_unit(u.nanometer ** -1)
                if u.is_quantity(used_alpha) else used_alpha,
                "grid": [used_nx, used_ny, used_nz],
                "bspline_order": 5,
                "bspline_order_note": (
                    "OpenMM's PME interpolation order is fixed at 5 in the source "
                    "(PME_ORDER in ReferencePME.cpp and the CPU/CUDA kernels). It is "
                    "not exposed through any API, so it cannot be matched to GMD's 4 "
                    "or 6. This is a real, unremovable misalignment, not an omission."
                ),
            },
            "reports_virial": False,
            "reports_virial_note": (
                "OpenMM exposes no virial tensor. State objects carry energy, forces "
                "and positions only; pressure is available solely through Monte Carlo "
                "barostat volume moves, which do not evaluate a virial."
            ),
        },
        sys.stdout,
        indent=2,
        sort_keys=True,
    )
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
