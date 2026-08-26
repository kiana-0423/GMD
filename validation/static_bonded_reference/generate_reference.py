#!/usr/bin/env python3
"""Regenerate reference.json for the bonded validation case from LAMMPS.

This is the provenance script for the checked-in reference values. It is NOT
run as part of the validation suite (validation must not depend on LAMMPS being
installed); it is run by hand when the fixture or the reference engine changes.

    python3 generate_reference.py --lammps $(which lmp_serial)

Everything about the conventions, units and the improper sign mapping is
documented in README.md next to this file.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys


# The same constant GMD uses internally (src/io/config_loader.cpp) to convert
# the kcal/mol coefficients in a .ff file to eV. LAMMPS runs in `units real`
# with numerically identical coefficients, so this converts its output to the
# units GMD reports. Because the identical constant is applied on both sides,
# its finite precision does not degrade the comparison.
KCAL_TO_EV = 4.336410e-2

CASE_DIR = pathlib.Path(__file__).resolve().parent


def run_lammps(lammps: str, input_file: str, log: str) -> None:
    completed = subprocess.run(
        [lammps, "-in", input_file, "-log", log],
        cwd=CASE_DIR, text=True, capture_output=True, check=False,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise RuntimeError(f"LAMMPS failed on {input_file}")


def lammps_version(log: str) -> str:
    text = (CASE_DIR / log).read_text()
    match = re.search(r"LAMMPS \(([^)]+)\)", text)
    if not match:
        raise RuntimeError(f"could not read the LAMMPS version from {log}")
    return match.group(1)


def read_energies(log: str) -> dict:
    """Pull the single thermo row produced by `run 0`."""
    lines = (CASE_DIR / log).read_text().splitlines()
    for index, line in enumerate(lines):
        if line.split()[:1] == ["Step"]:
            header = line.split()
            values = [float(v) for v in lines[index + 1].split()]
            row = dict(zip(header, values))
            return {
                "bond": row["E_bond"] * KCAL_TO_EV,
                "angle": row["E_angle"] * KCAL_TO_EV,
                "dihedral": row["E_dihed"] * KCAL_TO_EV,
                "improper": row["E_impro"] * KCAL_TO_EV,
                "pair": row["E_pair"] * KCAL_TO_EV,
                "total": row["PotEng"] * KCAL_TO_EV,
            }
    raise RuntimeError(f"no thermo row found in {log}")


def read_forces(dump: str) -> list[list[float]]:
    lines = (CASE_DIR / dump).read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS")) + 1
    forces = []
    for line in lines[start:]:
        fields = line.split()
        if len(fields) != 4:
            break
        forces.append([float(fields[1]) * KCAL_TO_EV,
                       float(fields[2]) * KCAL_TO_EV,
                       float(fields[3]) * KCAL_TO_EV])
    return forces


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lammps", required=True,
                        help="path to the LAMMPS executable, e.g. lmp_serial")
    args = parser.parse_args()

    run_lammps(args.lammps, "lammps.in", "lammps.log")
    run_lammps(args.lammps, "lammps_wrapped.in", "lammps_wrapped.log")

    intact = read_energies("lammps.log")
    wrapped = read_energies("lammps_wrapped.log")

    reference = {
        "source": "lammps",
        "engine": {
            "name": "LAMMPS",
            "version": lammps_version("lammps.log"),
            "units": "real (kcal/mol, Angstrom); converted to eV with "
                     f"kcal_to_eV = {KCAL_TO_EV!r}",
            "styles": {
                "pair": "zero 2.0 (no non-bonded contribution)",
                "bond": "harmonic       E = K (r - r0)^2",
                "angle": "harmonic      E = K (theta - theta0)^2",
                "dihedral": "charmm    E = K [1 + cos(n phi - d)]",
                "improper": "harmonic  E = K (chi - chi0)^2",
            },
            "command": "lmp_serial -in lammps.in -log lammps.log",
        },
        "atom_order": "GMD .xyz line order == LAMMPS atom-ID order (1..5)",
        "energy_ev": {
            "intact": intact,
            "wrapped": wrapped,
        },
        "force_ev_per_angstrom": {
            "intact": read_forces("forces_intact.txt"),
            "wrapped": read_forces("forces_wrapped.txt"),
        },
        "externally_validated_terms": ["bond", "angle", "dihedral", "improper"],
        "notes": (
            "The improper term is compared under a derived sign mapping: GMD "
            "phi0 = -chi0. See README.md."
        ),
    }

    path = CASE_DIR / "reference.json"
    with path.open("w", encoding="utf-8") as handle:
        json.dump(reference, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
