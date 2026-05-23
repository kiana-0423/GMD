from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import compare_scalar, force_error_metrics, load_json, run_command, write_json


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()

    case_dir = pathlib.Path(__file__).resolve().parent
    work_dir = pathlib.Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    for name in ["cluster.xyz", "cluster.run"]:
        shutil.copy2(case_dir / name, work_dir / name)

    actual_path = work_dir / "result.json"
    run_command(
        [
            args.gmd_validate,
            "cluster.xyz",
            "cluster.run",
            "--json",
            str(actual_path),
        ],
        work_dir,
    )

    actual = load_json(actual_path)
    reference = load_json(case_dir / "reference.json")
    tolerance = load_json(case_dir / "tolerance.json")

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

    summary = {
        "case": "static_lj_cluster",
        "components": components,
        "total_energy": total,
        "force_error": force_metrics,
        "passed": total["passed"]
        and all(metric["passed"] for metric in components.values())
        and force_metrics["passed"],
    }
    write_json(work_dir / "summary.json", summary)
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
