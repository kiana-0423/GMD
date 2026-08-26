from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import compare_scalar, force_error_metrics, load_json, run_command, write_json


CASE = "static_bonded_reference"

# GMD reports only a single aggregated "bonded" energy component, so each term
# is isolated with its own topology file and compared against the matching
# per-term energy LAMMPS reports from one run of the full fixture.
PER_TERM_TOPOLOGY = {
    "bond": "bonded_bond_only.top",
    "angle": "bonded_angle_only.top",
    "dihedral": "bonded_dihedral_only.top",
    "improper": "bonded_improper_only.top",
}

INPUTS = [
    "bonded.xyz", "bonded_wrapped.xyz", "bonded.ff", "bonded.run", "bonded.top",
    *PER_TERM_TOPOLOGY.values(),
]


def evaluate(gmd_validate: str, work_dir: pathlib.Path, xyz: str, topology: str,
             tag: str) -> dict:
    output = work_dir / f"result_{tag}.json"
    run_command(
        [gmd_validate, xyz, "bonded.run", "bonded.ff", topology,
         "--json", str(output)],
        work_dir,
    )
    return load_json(output)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()

    case_dir = pathlib.Path(__file__).resolve().parent
    work_dir = pathlib.Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    for name in INPUTS:
        shutil.copy2(case_dir / name, work_dir / name)

    reference = load_json(case_dir / "reference.json")
    tolerance = load_json(case_dir / "tolerance.json")

    summary: dict = {
        "case": CASE,
        "reference_engine": reference["engine"],
        "externally_validated_terms": reference["externally_validated_terms"],
    }
    passed = True

    # --- per-term energies, each against its own LAMMPS term ---------------
    per_term = {}
    for term, topology in PER_TERM_TOPOLOGY.items():
        result = evaluate(args.gmd_validate, work_dir, "bonded.xyz", topology, term)
        # Only the bonded component is under test; the fixture is built so the
        # other two are identically zero, which is asserted below.
        per_term[term] = compare_scalar(
            result["energy"]["components"]["bonded"],
            reference["energy_ev"]["intact"][term],
            tolerance["energy_abs"],
        )
        passed = passed and per_term[term]["passed"]
    summary["per_term_energy"] = per_term

    # --- full fixture: total energy and every force component --------------
    variants = {}
    for tag, xyz in (("intact", "bonded.xyz"), ("wrapped", "bonded_wrapped.xyz")):
        result = evaluate(args.gmd_validate, work_dir, xyz, "bonded.top", f"full_{tag}")

        components = result["energy"]["components"]
        # The fixture isolates the bonded terms; if this ever stops holding the
        # comparison against a bonded-only reference is meaningless.
        non_bonded_clean = components["lj"] == 0.0 and components["coulomb"] == 0.0

        total = compare_scalar(
            components["bonded"],
            reference["energy_ev"][tag]["total"],
            tolerance["energy_abs"],
        )
        forces = force_error_metrics(
            result["force"]["atoms"],
            reference["force_ev_per_angstrom"][tag],
        )
        forces["passed"] = (
            forces["rms_error"] <= tolerance["force_rms_abs"]
            and forces["max_component_error"] <= tolerance["force_max_abs"]
        )

        variants[tag] = {
            "total_energy": total,
            "force_error": forces,
            "non_bonded_is_zero": non_bonded_clean,
            "passed": total["passed"] and forces["passed"] and non_bonded_clean,
        }
        passed = passed and variants[tag]["passed"]
    summary["variants"] = variants

    # --- wrapping must not change anything ---------------------------------
    wrap_shift = abs(
        variants["intact"]["total_energy"]["actual"]
        - variants["wrapped"]["total_energy"]["actual"]
    )
    summary["wrapped_matches_intact"] = {
        "abs_error": wrap_shift,
        "tolerance": tolerance["energy_abs"],
        "passed": wrap_shift <= tolerance["energy_abs"],
    }
    passed = passed and summary["wrapped_matches_intact"]["passed"]

    summary["passed"] = passed
    write_json(work_dir / "summary.json", summary)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
