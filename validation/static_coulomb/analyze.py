from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import compare_scalar, force_error_metrics, load_json, run_command, write_json


def evaluate(work_dir: pathlib.Path, gmd_validate: str, run_name: str, reference_name: str, tolerance: dict) -> dict:
    actual_path = work_dir / f"{run_name}.json"
    run_command(
        [
            gmd_validate,
            "charges.xyz",
            f"{run_name}.run",
            "--json",
            str(actual_path),
        ],
        work_dir,
    )
    actual = load_json(actual_path)
    reference = load_json(pathlib.Path(__file__).resolve().parent / reference_name)

    components = {}
    for key in ["lj", "bonded", "coulomb"]:
        components[key] = compare_scalar(
            actual["energy"]["components"][key],
            reference["energy"]["components"][key],
            tolerance["component_abs"],
        )
    total = compare_scalar(
        actual["energy"]["total"],
        reference["energy"]["total"],
        tolerance["total_abs"],
    )
    force_metrics = force_error_metrics(
        actual["force"]["atoms"],
        reference["force"]["atoms"],
    )
    force_metrics["passed"] = (
        force_metrics["rms_error"] <= tolerance["force_rms_abs"]
        and force_metrics["max_component_error"] <= tolerance["force_max_abs"]
    )

    return {
        "reference_source": reference.get("source", "unknown"),
        "components": components,
        "total_energy": total,
        "force_error": force_metrics,
        "passed": total["passed"]
        and all(metric["passed"] for metric in components.values())
        and force_metrics["passed"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()

    case_dir = pathlib.Path(__file__).resolve().parent
    work_dir = pathlib.Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    for name in ["charges.xyz", "ewald.run", "pme.run"]:
        shutil.copy2(case_dir / name, work_dir / name)

    tolerance = load_json(case_dir / "tolerance.json")
    summary = {
        "case": "static_coulomb",
        "ewald": evaluate(work_dir, args.gmd_validate, "ewald", "reference_ewald.json", tolerance),
        "pme": evaluate(work_dir, args.gmd_validate, "pme", "reference_pme.json", tolerance),
    }
    summary["passed"] = summary["ewald"]["passed"] and summary["pme"]["passed"]
    write_json(work_dir / "summary.json", summary)
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
