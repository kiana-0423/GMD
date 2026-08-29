from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


# CODATA 2022 k_e = E_h * a0 = 27.211386245981 * 0.529177210544, matching
# gmd::kCoulombConstant in include/gmd/core/physical_constants.hpp. Declared
# here rather than read from the C++ source on purpose: this module is an
# independent analytic reference, and importing the production value would
# make the comparisons it feeds self-confirming. It must track production,
# and tests/electrostatic_constant_tests.cpp is what fails if it stops.
COULOMB = 14.3996454686836
KCAL_PER_MOL_TO_EV = 4.336410e-2


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def minimum_image(dr: list[float], box: list[float]) -> list[float]:
    out = dr[:]
    for dim, length in enumerate(box):
        while out[dim] > 0.5 * length:
            out[dim] -= length
        while out[dim] < -0.5 * length:
            out[dim] += length
    return out


def add_force(forces: list[list[float]], i: int, j: int, factor: float, dr: list[float]) -> None:
    for dim in range(3):
        value = factor * dr[dim]
        forces[i][dim] += value
        forces[j][dim] -= value


def shifted_lj(
    coords: list[list[float]],
    box: list[float],
    epsilon: float,
    sigma: float,
    cutoff: float,
    pair_scales: dict[tuple[int, int], float] | None = None,
) -> tuple[float, list[list[float]]]:
    n = len(coords)
    forces = [[0.0, 0.0, 0.0] for _ in range(n)]
    energy = 0.0
    cutoff_sq = cutoff * cutoff
    eps4 = 4.0 * epsilon
    sig2 = sigma * sigma
    s2_cut = sig2 / cutoff_sq
    s6_cut = s2_cut * s2_cut * s2_cut
    shift = eps4 * (s6_cut * s6_cut - s6_cut)
    pair_scales = pair_scales or {}

    for i in range(n - 1):
        for j in range(i + 1, n):
            scale = pair_scales.get((i, j), 1.0)
            if scale == 0.0:
                continue
            dr = minimum_image([coords[i][d] - coords[j][d] for d in range(3)], box)
            r2 = sum(value * value for value in dr)
            if r2 >= cutoff_sq or r2 < 1.0e-12:
                continue
            s2 = sig2 / r2
            s6 = s2 * s2 * s2
            s12 = s6 * s6
            pair_energy = eps4 * (s12 - s6)
            factor = eps4 * (12.0 * s12 - 6.0 * s6) / r2
            energy += scale * (pair_energy - shift)
            add_force(forces, i, j, scale * factor, dr)
    return energy, forces


def ewald(
    coords: list[list[float]],
    charges: list[float],
    box: list[float],
    alpha: float,
    kmax: int,
    cutoff: float,
    special_coulomb_scales: dict[tuple[int, int], float] | None = None,
) -> tuple[float, list[list[float]]]:
    n = len(coords)
    forces = [[0.0, 0.0, 0.0] for _ in range(n)]
    energy = 0.0
    cutoff_sq = cutoff * cutoff
    two_alpha_over_sqrt_pi = 2.0 * alpha / math.sqrt(math.pi)
    special_coulomb_scales = special_coulomb_scales or {}

    # Real-space Ewald sum.
    for i in range(n - 1):
        for j in range(i + 1, n):
            qi = charges[i]
            qj = charges[j]
            if qi == 0.0 and qj == 0.0:
                continue
            dr = minimum_image([coords[i][d] - coords[j][d] for d in range(3)], box)
            r2 = sum(value * value for value in dr)
            if r2 >= cutoff_sq or r2 < 1.0e-12:
                continue
            r = math.sqrt(r2)
            ar = alpha * r
            erfc_ar = math.erfc(ar)
            exp_ar2 = math.exp(-ar * ar)
            energy += COULOMB * qi * qj * erfc_ar / r
            factor = (
                COULOMB
                * qi
                * qj
                / r2
                * (erfc_ar / r + two_alpha_over_sqrt_pi * exp_ar2)
            )
            add_force(forces, i, j, factor, dr)

    # Reciprocal-space Ewald sum.
    lx, ly, lz = box
    volume = lx * ly * lz
    twopi_over_l = [2.0 * math.pi / lx, 2.0 * math.pi / ly, 2.0 * math.pi / lz]
    inv_4a2 = 1.0 / (4.0 * alpha * alpha)
    for nx in range(-kmax, kmax + 1):
        kx = nx * twopi_over_l[0]
        for ny in range(-kmax, kmax + 1):
            ky = ny * twopi_over_l[1]
            for nz in range(-kmax, kmax + 1):
                if nx == 0 and ny == 0 and nz == 0:
                    continue
                kz = nz * twopi_over_l[2]
                k2 = kx * kx + ky * ky + kz * kz
                gfactor = (
                    COULOMB
                    * 4.0
                    * math.pi
                    / (volume * k2)
                    * math.exp(-k2 * inv_4a2)
                )
                s_re = 0.0
                s_im = 0.0
                for atom, q in enumerate(charges):
                    phi = kx * coords[atom][0] + ky * coords[atom][1] + kz * coords[atom][2]
                    s_re += q * math.cos(phi)
                    s_im += q * math.sin(phi)
                energy += 0.5 * gfactor * (s_re * s_re + s_im * s_im)
                for atom, q in enumerate(charges):
                    if q == 0.0:
                        continue
                    phi = kx * coords[atom][0] + ky * coords[atom][1] + kz * coords[atom][2]
                    c = math.cos(phi)
                    s = math.sin(phi)
                    pref = gfactor * q * (s_re * s - s_im * c)
                    forces[atom][0] += pref * kx
                    forces[atom][1] += pref * ky
                    forces[atom][2] += pref * kz

    # Self and neutralizing-background corrections.
    q2_sum = sum(q * q for q in charges)
    q_net = sum(charges)
    energy -= COULOMB * (alpha / math.sqrt(math.pi)) * q2_sum
    if q_net != 0.0:
        energy -= COULOMB * math.pi / (2.0 * volume * alpha * alpha) * q_net * q_net

    # Special-pair correction: convert the full unscaled periodic Coulomb sum
    # into topology-scaled direct interactions for 1-2/1-3/1-4 pairs.
    for (i, j), scale in special_coulomb_scales.items():
        delta = scale - 1.0
        if delta == 0.0:
            continue
        qi = charges[i]
        qj = charges[j]
        if qi == 0.0 or qj == 0.0:
            continue
        dr = minimum_image([coords[i][d] - coords[j][d] for d in range(3)], box)
        r2 = sum(value * value for value in dr)
        if r2 < 1.0e-12:
            continue
        r = math.sqrt(r2)
        correction = delta * COULOMB * qi * qj
        energy += correction / r
        add_force(forces, i, j, correction / (r2 * r), dr)

    return energy, forces


def sum_forces(*force_terms: list[list[float]]) -> list[list[float]]:
    n = len(force_terms[0])
    out = [[0.0, 0.0, 0.0] for _ in range(n)]
    for term in force_terms:
        for atom in range(n):
            for dim in range(3):
                out[atom][dim] += term[atom][dim]
    return out


def special_pair_reference(lj14: float, coul14: float) -> dict:
    coords = [
        [5.0, 5.0, 5.0],
        [6.0, 5.0, 5.0],
        [7.0, 5.0, 5.0],
        [8.0, 5.0, 5.0],
    ]
    charges = [1.0, -1.0, 1.0, -1.0]
    box = [30.0, 30.0, 30.0]
    lj_scales = {
        (0, 1): 0.0,
        (1, 2): 0.0,
        (2, 3): 0.0,
        (0, 2): 0.0,
        (1, 3): 0.0,
        (0, 3): lj14,
    }
    coul_scales = {
        (0, 1): 0.0,
        (1, 2): 0.0,
        (2, 3): 0.0,
        (0, 2): 0.0,
        (1, 3): 0.0,
        (0, 3): coul14,
    }
    lj_energy, lj_forces = shifted_lj(
        coords,
        box,
        epsilon=0.2 * KCAL_PER_MOL_TO_EV,
        sigma=1.0,
        cutoff=9.0,
        pair_scales=lj_scales,
    )
    coulomb_energy, coulomb_forces = ewald(
        coords,
        charges,
        box,
        alpha=0.3,
        kmax=7,
        cutoff=9.0,
        special_coulomb_scales=coul_scales,
    )
    total_forces = sum_forces(lj_forces, coulomb_forces)
    return {
        "energy": {
            "components": {
                "bonded": 0.0,
                "coulomb": coulomb_energy,
                "lj": lj_energy,
            },
            "total": lj_energy + coulomb_energy,
        },
        "force": {"atoms": total_forces},
    }


def static_coulomb_ewald_reference() -> dict:
    coords = []
    charges = []
    rows = [
        (2.0, 2.0, 2.0, 0.4),
        (8.0, 2.0, 2.0, -0.4),
        (2.0, 8.0, 2.0, 0.4),
        (8.0, 8.0, 2.0, -0.4),
        (2.0, 2.0, 8.0, 0.4),
        (8.0, 2.0, 8.0, -0.4),
        (2.0, 8.0, 8.0, 0.4),
        (8.0, 8.0, 8.0, -0.4),
    ]
    for x, y, z, q in rows:
        coords.append([x, y, z])
        charges.append(q)
    box = [20.0, 20.0, 20.0]
    lj_energy, lj_forces = shifted_lj(coords, box, epsilon=0.00774, sigma=2.5, cutoff=8.0)
    coulomb_energy, coulomb_forces = ewald(coords, charges, box, alpha=0.3, kmax=3, cutoff=8.0)
    total_forces = sum_forces(lj_forces, coulomb_forces)
    return {
        "energy": {
            "components": {
                "bonded": 0.0,
                "coulomb": coulomb_energy,
                "lj": lj_energy,
            },
            "total": lj_energy + coulomb_energy,
        },
        "force": {"atoms": total_forces},
    }


def with_metadata(payload: dict, source: str, description: str, formulas: list[str]) -> dict:
    out = dict(payload)
    out["source"] = source
    out["reference_details"] = {
        "description": description,
        "units": "GMD internal units: eV, Angstrom, elementary charge",
        "formulas": formulas,
        # Emitted rather than hand-added, so it survives the next --write. Every
        # Coulomb quantity in this file is exactly proportional to this constant;
        # the Lennard-Jones ones do not depend on it at all, which is why the
        # totals do not scale when it changes.
        "coulomb_constant_ev_angstrom": COULOMB,
        "coulomb_constant_source": (
            "CODATA 2022 via E_h * a0 = 27.211386245981 eV * 0.529177210544 A, "
            "matching gmd::kCoulombConstant in "
            "include/gmd/core/physical_constants.hpp"
        ),
        "coulomb_constant_superseded": 14.3996,
        "coulomb_constant_change": (
            "Regenerated 2026-08-30. The constant was corrected from 14.3996, a "
            "five-significant-figure truncation, to the CODATA 2022 value; every "
            "Coulomb energy and force here changed by +3.157635e-06 relative and "
            "every Lennard-Jones one is bit-identical."
        ),
        "regeneration_command": (
            "python3 validation/analytic_references.py --root validation --write, "
            "run from <repo>"
        ),
    }
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(Path(__file__).resolve().parent))
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    root = Path(args.root)

    special = with_metadata(
        special_pair_reference(lj14=0.5, coul14=0.833333333333),
        "analytic_reference_python_ewald_shifted_lj_special_pairs",
        "Four-atom 1-2-3-4 chain with 1-2/1-3 exclusions and scaled 1-4 direct corrections.",
        [
            "shifted LJ: 4 epsilon [(sigma/r)^12 - (sigma/r)^6] - U(rc)",
            "periodic Ewald real + reciprocal + self/background corrections",
            "special Coulomb correction: (scale - 1) * k_e q_i q_j / r",
        ],
    )
    special["scale_variants"] = {
        "lj14_0.25_coul14_0.5": with_metadata(
            special_pair_reference(lj14=0.25, coul14=0.5),
            "analytic_reference_python_ewald_shifted_lj_special_pairs",
            "Same four-atom chain with modified 1-4 LJ/Coulomb scale factors.",
            [
                "shifted LJ 1-4 contribution scales linearly with lj_scale_14",
                "direct special-pair Coulomb correction changes linearly with coul_scale_14",
            ],
        )
    }

    ewald_ref = with_metadata(
        static_coulomb_ewald_reference(),
        "analytic_reference_python_periodic_ewald",
        "Eight-charge periodic system evaluated with explicit Ewald real, reciprocal, self, and background terms.",
        [
            "U_real = k_e sum_{i<j} q_i q_j erfc(alpha r_ij) / r_ij",
            "U_recip = 1/2 sum_{k!=0} k_e 4 pi exp(-k^2/4 alpha^2) |S(k)|^2 / (V k^2)",
            "U_self = -k_e alpha / sqrt(pi) sum_i q_i^2",
            "U_background = -k_e pi Q^2 / (2 V alpha^2)",
            "shifted LJ: 4 epsilon [(sigma/r)^12 - (sigma/r)^6] - U(rc)",
        ],
    )

    if args.write:
        write_json(root / "static_special_pairs" / "reference.json", special)
        write_json(root / "static_coulomb" / "reference_ewald.json", ewald_ref)
    else:
        print(json.dumps({"static_special_pairs": special, "static_coulomb_ewald": ewald_ref}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
