"""Compare GMD PME against the checked-in external reference.

This runs in normal CI. It needs Python and a built gmd_validate, and nothing
else: OpenMM and LAMMPS are required only to regenerate the reference, which is
a separate manual step (see generate_external_reference.py).

What makes this an external validation rather than a regression baseline is the
reference file, not this script. So the first thing checked is that the file is
still what it claims to be: an external engine's output, for this geometry, at
these settings. Only then are the numbers compared.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import shutil
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import load_json, write_json  # noqa: E402

CASE = pathlib.Path(__file__).resolve().parent
REFERENCE_RUN = "pme_g64_o6.run"
EWALD_RUN = "ewald_converged.run"
FIXTURES = ["charges.xyz", "charges_wrapped.xyz", "charges_translated.xyz"]
VIRIAL_ORDER = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


class Checks:
    """Collects pass/fail rows so a failing run reports every problem at once."""

    def __init__(self) -> None:
        self.rows: list[dict] = []

    def record(self, name: str, ok: bool, detail: str, **extra) -> bool:
        self.rows.append({"check": name, "passed": bool(ok), "detail": detail, **extra})
        return bool(ok)

    def close(self, name: str, actual: float, limit: float, detail: str, **extra) -> bool:
        return self.record(
            name, actual <= limit,
            f"{detail}: {actual:.6e} (limit {limit:.6e})",
            actual=actual, limit=limit, **extra)

    @property
    def passed(self) -> bool:
        return all(row["passed"] for row in self.rows)


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_gmd(binary: str, work: pathlib.Path, xyz: str, run_file: str, tag: str) -> dict:
    output = work / f"gmd_{tag}.json"
    completed = subprocess.run([binary, xyz, run_file, "--json", str(output)],
                               cwd=work, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(f"gmd_validate failed for {tag}:\n{completed.stderr}")
    return json.loads(output.read_text())


def worst_force_difference(actual: list[list[float]],
                           expected: list[list[float]]) -> tuple[float, int, int]:
    """Largest single force-component difference, and where it is.

    Identifying the atom and the axis is the whole point: "max force error
    3e-06" says nothing about whether one atom is wrong or all of them are.
    Atom 0 axis 0 is returned when every difference is exactly zero, which is
    reported as a failure by the independence check rather than here.
    """
    worst, atom, axis = 0.0, 0, 0
    for index, (row_a, row_b) in enumerate(zip(actual, expected)):
        for component in range(3):
            delta = abs(row_a[component] - row_b[component])
            if delta > worst:
                worst, atom, axis = delta, index, component
    return worst, atom, axis


def force_rms(actual: list[list[float]], expected: list[list[float]]) -> float:
    total = sum((a - b) ** 2 for row_a, row_b in zip(actual, expected)
                for a, b in zip(row_a, row_b))
    return math.sqrt(total / (3 * len(expected)))


def net_force(forces: list[list[float]]) -> float:
    return max(abs(sum(row[component] for row in forces)) for component in range(3))


def read_run_parameters(path: pathlib.Path) -> dict:
    values: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        stripped = line.split("#", 1)[0].split()
        if stripped:
            values[stripped[0]] = stripped[1:]
    return values


def check_reference_is_external(checks: Checks, reference: dict) -> None:
    """The reference must still be an external engine's output for this fixture."""
    source = reference.get("source", "")
    checks.record(
        "reference_source_is_external", source.startswith("external_"),
        f"reference source is {source!r}; an external reference must declare "
        f"a source beginning with 'external_'")

    block = reference.get("reference", {})
    engine, version = block.get("engine", ""), block.get("version", "")
    checks.record(
        "reference_names_an_external_engine",
        engine in {"OpenMM", "LAMMPS"} and bool(version),
        f"reference engine {engine!r} version {version!r}")

    provenance = reference.get("provenance", {})
    recorded = provenance.get("fixture_sha256", {})
    for name in FIXTURES:
        actual = sha256(CASE / name)
        checks.record(
            f"fixture_hash_{name}", recorded.get(name) == actual,
            f"{name} sha256 {actual} vs recorded {recorded.get(name)}; the "
            f"reference was computed for a different geometry if these differ")

    recorded_inputs = provenance.get("input_sha256", {})
    for name in [REFERENCE_RUN, EWALD_RUN]:
        actual = sha256(CASE / name)
        checks.record(
            f"input_hash_{name}", recorded_inputs.get(name) == actual,
            f"{name} sha256 {actual} vs recorded {recorded_inputs.get(name)}")

    # The settings GMD is about to be run at must be the settings the external
    # engine was run at. A hash cannot catch this on its own, because the two
    # sides read different files.
    parameters = read_run_parameters(CASE / REFERENCE_RUN)
    checks.record(
        "run_file_matches_reference_alpha",
        float(parameters["pme_alpha"][0]) == block.get("alpha_inv_angstrom"),
        f"{REFERENCE_RUN} pme_alpha {parameters['pme_alpha'][0]} vs reference "
        f"{block.get('alpha_inv_angstrom')}")
    checks.record(
        "run_file_matches_reference_cutoff",
        float(parameters["pme_cutoff"][0]) == block.get("cutoff_angstrom"),
        f"{REFERENCE_RUN} pme_cutoff {parameters['pme_cutoff'][0]} vs reference "
        f"{block.get('cutoff_angstrom')}")
    checks.record(
        "run_file_matches_reference_grid",
        [int(value) for value in parameters["pme_grid"]] == block.get("grid"),
        f"{REFERENCE_RUN} pme_grid {parameters['pme_grid']} vs reference "
        f"{block.get('grid')}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()

    work = pathlib.Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    for name in FIXTURES:
        shutil.copy2(CASE / name, work / name)
    for run_file in CASE.glob("*.run"):
        shutil.copy2(run_file, work / run_file.name)

    reference = load_json(CASE / "reference_external.json")
    tolerance = load_json(CASE / "tolerance.json")
    checks = Checks()

    check_reference_is_external(checks, reference)

    external = reference["reference"]
    external_energy = external["energy_ev"]
    external_forces = external["forces_ev_per_angstrom"]

    pme = run_gmd(args.gmd_validate, work, "charges.xyz", REFERENCE_RUN, "pme")
    ewald = run_gmd(args.gmd_validate, work, "charges.xyz", EWALD_RUN, "ewald")
    wrapped = run_gmd(args.gmd_validate, work, "charges_wrapped.xyz", REFERENCE_RUN, "wrapped")
    translated = run_gmd(args.gmd_validate, work, "charges_translated.xyz",
                         REFERENCE_RUN, "translated")

    gmd_energy = pme["energy"]["total"]
    gmd_forces = pme["force"]["atoms"]

    # --- the reference must not be GMD's own output -------------------------
    energy_gap = abs(gmd_energy - external_energy)
    force_gap, _, _ = worst_force_difference(gmd_forces, external_forces)
    checks.record(
        "reference_is_independent_of_gmd",
        energy_gap > tolerance["independence_floor_energy"]
        and force_gap > tolerance["independence_floor_force"],
        f"|GMD - reference| energy {energy_gap:.6e}, force {force_gap:.6e}; both "
        f"must exceed {tolerance['independence_floor_energy']:.0e}. Two different "
        f"PME implementations at different interpolation orders cannot agree this "
        f"closely, so a smaller gap means reference_external.json holds "
        f"GMD-generated numbers and is not an external reference at all")

    # --- energy and forces against the external reference -------------------
    checks.close("energy_vs_external", energy_gap,
                 tolerance["energy_abs_vs_external"],
                 f"GMD {gmd_energy:.12f} eV vs {external['engine']} "
                 f"{external_energy:.12f} eV, absolute error")
    checks.close("energy_rel_vs_external", energy_gap / abs(external_energy),
                 tolerance["energy_rel_vs_external"], "relative energy error")

    worst, atom, axis = worst_force_difference(gmd_forces, external_forces)
    checks.close("force_max_vs_external", worst,
                 tolerance["force_max_abs_vs_external"],
                 f"worst force component is atom {atom} {'xyz'[axis]} "
                 f"(GMD {gmd_forces[atom][axis]:.9f}, reference "
                 f"{external_forces[atom][axis]:.9f}); error",
                 atom=atom, component="xyz"[axis])
    checks.close("force_rms_vs_external", force_rms(gmd_forces, external_forces),
                 tolerance["force_rms_abs_vs_external"], "force RMS error")

    # A mesh breaks translational invariance, so PME's net force is nonzero at
    # the mesh-error level in both codes. Exact Ewald has no mesh, so its net
    # force is held to round-off instead -- that is the sharp check, and it is
    # what would catch a sign or index error in the pair loop.
    checks.close("pme_net_force_gmd", net_force(gmd_forces),
                 tolerance["pme_net_force_max_abs"], "GMD PME net force")
    checks.close("pme_net_force_external", net_force(external_forces),
                 tolerance["pme_net_force_max_abs"],
                 f"{external['engine']} PME net force")
    checks.close("ewald_net_force_gmd", net_force(ewald["force"]["atoms"]),
                 tolerance["ewald_net_force_max_abs"],
                 "GMD exact-Ewald net force, which has no mesh and so must "
                 "vanish to round-off")

    # --- second external engine ---------------------------------------------
    lammps = reference["cross_reference"]["lammps_pppm"]
    checks.close("energy_vs_lammps_pppm", abs(gmd_energy - lammps["energy_ev"]),
                 tolerance["energy_abs_vs_lammps_pppm"],
                 f"GMD vs LAMMPS PPPM ({lammps['version']}) energy")
    worst, atom, axis = worst_force_difference(
        gmd_forces, lammps["forces_ev_per_angstrom"])
    checks.close("force_max_vs_lammps_pppm", worst,
                 tolerance["force_max_abs_vs_lammps_pppm"],
                 f"GMD vs LAMMPS PPPM worst force component, atom {atom} {'xyz'[axis]}",
                 atom=atom, component="xyz"[axis])

    # --- the exact limit both mesh methods approximate ----------------------
    analytic = reference["analytic_ewald"]
    checks.close("pme_energy_vs_analytic_ewald",
                 abs(gmd_energy - analytic["energy_ev"]),
                 tolerance["energy_abs_vs_analytic_ewald"],
                 "GMD PME vs exact Ewald energy (mesh error)")
    worst, atom, axis = worst_force_difference(
        gmd_forces, analytic["forces_ev_per_angstrom"])
    checks.close("pme_force_vs_analytic_ewald", worst,
                 tolerance["force_max_abs_vs_analytic_ewald"],
                 f"GMD PME vs exact Ewald worst force component, atom {atom} "
                 f"{'xyz'[axis]} (mesh error)", atom=atom, component="xyz"[axis])
    checks.close("ewald_energy_vs_analytic_ewald",
                 abs(ewald["energy"]["total"] - analytic["energy_ev"]),
                 tolerance["gmd_ewald_energy_abs_vs_analytic_ewald"],
                 "GMD Ewald vs the same sum in Python energy")
    worst, atom, axis = worst_force_difference(
        ewald["force"]["atoms"], analytic["forces_ev_per_angstrom"])
    checks.close("ewald_force_vs_analytic_ewald", worst,
                 tolerance["gmd_ewald_force_max_abs_vs_analytic_ewald"],
                 f"GMD Ewald vs the same sum in Python, worst component atom "
                 f"{atom} {'xyz'[axis]}", atom=atom, component="xyz"[axis])

    # --- virial, against LAMMPS only ----------------------------------------
    # OpenMM exposes no virial, so this is the one quantity the primary
    # reference cannot cover.
    for label, result, lammps_block, limit_key in [
        ("ewald", ewald, reference["cross_reference"]["lammps_ewald"],
         "virial_max_abs_ewald_vs_lammps"),
        ("pme", pme, reference["cross_reference"]["lammps_pppm"],
         "virial_max_abs_pme_vs_lammps"),
    ]:
        tensor = result["virial"]["tensor"]
        checks.record(f"virial_valid_{label}", result["virial"]["valid"],
                      f"{label} reported virial_valid={result['virial']['valid']}")
        asymmetry = max(abs(tensor[a * 3 + b] - tensor[b * 3 + a])
                        for a in range(3) for b in range(3))
        # LAMMPS reports six components because the tensor is symmetric; the
        # remaining three GMD components are covered by checking that GMD's own
        # tensor is symmetric rather than by assuming it.
        checks.close(f"virial_symmetry_{label}", asymmetry,
                     tolerance["virial_max_asymmetry"],
                     f"{label} GMD virial asymmetry, which is what makes the six "
                     f"LAMMPS components cover all nine GMD components")
        for (a, b), lammps_value in zip(
                VIRIAL_ORDER, lammps_block["virial_ev_xx_yy_zz_xy_xz_yz"]):
            name = f"W_{'xyz'[a]}{'xyz'[b]}"
            checks.close(f"virial_{label}_{name}",
                         abs(tensor[a * 3 + b] - lammps_value),
                         tolerance[limit_key],
                         f"{label} {name}: GMD {tensor[a * 3 + b]:.9f} vs LAMMPS "
                         f"{lammps_value:.9f}", component=name)

    # --- same configuration, different periodic images ----------------------
    checks.close("wrapped_energy", abs(wrapped["energy"]["total"] - gmd_energy),
                 tolerance["wrapped_energy_abs"],
                 "energy change when every atom is written in a different "
                 "periodic image")
    worst, atom, axis = worst_force_difference(wrapped["force"]["atoms"], gmd_forces)
    checks.close("wrapped_force", worst, tolerance["wrapped_force_max_abs"],
                 f"worst force change under re-imaging, atom {atom} {'xyz'[axis]}",
                 atom=atom, component="xyz"[axis])

    # --- rigid translation by a non-lattice vector --------------------------
    # Exact for an exact method; for a mesh method it is only invariant to the
    # mesh error, because the atoms move relative to a grid fixed to the box.
    checks.close("translated_energy", abs(translated["energy"]["total"] - gmd_energy),
                 tolerance["translated_energy_abs"],
                 "energy change under rigid translation by a non-lattice vector")
    worst, atom, axis = worst_force_difference(translated["force"]["atoms"], gmd_forces)
    checks.close("translated_force", worst, tolerance["translated_force_max_abs"],
                 f"worst force change under rigid translation, atom {atom} {'xyz'[axis]}",
                 atom=atom, component="xyz"[axis])

    summary = {
        "case": "pme_external",
        "reference_source": reference.get("source"),
        "reference_engine": f"{external['engine']} {external['version']}",
        "cross_reference_engine": lammps["version"],
        "settings": {
            "alpha_inv_angstrom": external["alpha_inv_angstrom"],
            "cutoff_angstrom": external["cutoff_angstrom"],
            "grid": external["grid"],
            "gmd_bspline_order": 6,
            "external_bspline_order": external["bspline_order"],
        },
        "checks": checks.rows,
        "passed": checks.passed,
    }
    write_json(work / "summary.json", summary)

    if not checks.passed:
        for row in checks.rows:
            if not row["passed"]:
                sys.stderr.write(f"FAIL {row['check']}: {row['detail']}\n")
        return 1
    print(f"pme_external: {len(checks.rows)} checks passed against "
          f"{external['engine']} {external['version']} and {lammps['version']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
