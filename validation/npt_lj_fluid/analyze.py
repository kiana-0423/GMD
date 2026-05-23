from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import compare_scalar, load_json, mean, parse_energy_log, run_command, stddev, write_json, write_series_csv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()

    case_dir = pathlib.Path(__file__).resolve().parent
    work_dir = pathlib.Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    for name in ["lj32.xyz", "npt.run"]:
        shutil.copy2(case_dir / name, work_dir / name)

    run_command([args.gmd, "lj32.xyz", "npt.run"], work_dir)
    rows = parse_energy_log(work_dir / "output.log")
    write_series_csv(
        work_dir / "thermo_series.csv",
        rows,
        ["step", "time_fs", "pe", "ke", "total_energy", "temperature", "pressure_bar", "volume_a3"],
    )
    temperatures = [row["temperature"] for row in rows]
    pressures = [row["pressure_bar"] for row in rows]

    reference = load_json(case_dir / "reference.json")
    tolerance = load_json(case_dir / "tolerance.json")
    metrics = {
        "temperature_mean": compare_scalar(
            mean(temperatures),
            reference["metrics"]["temperature_mean"],
            tolerance["temperature_mean_abs"],
        ),
        "temperature_stddev": compare_scalar(
            stddev(temperatures),
            reference["metrics"]["temperature_stddev"],
            tolerance["temperature_stddev_abs"],
        ),
        "pressure_mean": compare_scalar(
            mean(pressures),
            reference["metrics"]["pressure_mean"],
            tolerance["pressure_mean_abs"],
        ),
        "pressure_stddev": compare_scalar(
            stddev(pressures),
            reference["metrics"]["pressure_stddev"],
            tolerance["pressure_stddev_abs"],
        ),
    }
    summary = {
        "case": "npt_lj_fluid",
        "metrics": metrics,
        "passed": all(metric["passed"] for metric in metrics.values()),
    }
    write_json(work_dir / "summary.json", summary)
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
