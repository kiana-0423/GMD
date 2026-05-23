from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import compare_scalar, load_json, parse_energy_log, run_command, write_json, write_series_csv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()

    case_dir = pathlib.Path(__file__).resolve().parent
    work_dir = pathlib.Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    for name in ["lj32.xyz", "nve.run"]:
        shutil.copy2(case_dir / name, work_dir / name)

    run_command([args.gmd, "lj32.xyz", "nve.run"], work_dir)
    rows = parse_energy_log(work_dir / "output.log")
    write_series_csv(
        work_dir / "thermo_series.csv",
        rows,
        ["step", "time_fs", "pe", "ke", "total_energy", "temperature", "pressure_bar", "volume_a3"],
    )
    atom_count = 32.0
    total_time_ps = (rows[-1]["time_fs"] - rows[0]["time_fs"]) * 1.0e-3
    drift = (rows[-1]["total_energy"] - rows[0]["total_energy"]) / (atom_count * total_time_ps)

    reference = load_json(case_dir / "reference.json")
    tolerance = load_json(case_dir / "tolerance.json")
    metric = compare_scalar(
        drift,
        reference["metrics"]["energy_drift_per_atom_ps"],
        tolerance["energy_drift_per_atom_ps_abs"],
    )
    summary = {
        "case": "nve_lj_fluid",
        "metrics": {"energy_drift_per_atom_ps": metric},
        "passed": metric["passed"],
    }
    write_json(work_dir / "summary.json", summary)
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
