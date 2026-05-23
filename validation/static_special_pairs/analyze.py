from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import compare_scalar, force_error_metrics, load_json, run_command, write_json


def evaluate(work_dir: pathlib.Path,
             gmd_validate: str,
             run_name: str,
             reference: dict,
             tolerance: dict) -> dict:
    actual_path = work_dir / f"{run_name}.json"
    run_command(
        [
            gmd_validate,
            "chain.xyz",
            f"{run_name}.run",
            "chain.ff",
            "chain.top",
            "--json",
            str(actual_path),
        ],
        work_dir,
    )

    actual = load_json(actual_path)
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

    for name in ["chain.xyz", "chain.run", "chain.ff", "chain.top"]:
        shutil.copy2(case_dir / name, work_dir / name)

    variant_run = (case_dir / "chain.run").read_text(encoding="utf-8")
    variant_run = variant_run.replace("lj_scale_14 0.5", "lj_scale_14 0.25")
    variant_run = variant_run.replace("coul_scale_14 0.833333333333", "coul_scale_14 0.5")
    (work_dir / "chain_variant.run").write_text(variant_run, encoding="utf-8")

    reference = load_json(case_dir / "reference.json")
    tolerance = load_json(case_dir / "tolerance.json")

    default = evaluate(work_dir, args.gmd_validate, "chain", reference, tolerance)
    variant_reference = reference["scale_variants"]["lj14_0.25_coul14_0.5"]
    variant = evaluate(work_dir, args.gmd_validate, "chain_variant", variant_reference, tolerance)

    summary = {
        "case": "static_special_pairs",
        "reference_source": reference["source"],
        "default_scales": default,
        "scale_variant_lj14_0.25_coul14_0.5": variant,
        "passed": default["passed"] and variant["passed"],
    }
    write_json(work_dir / "summary.json", summary)
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
