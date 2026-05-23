from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import compare_scalar, compute_msd_series, fit_diffusion_coefficient, load_json, parse_energy_log, parse_xyz_frames, run_command, write_json, write_series_csv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()

    case_dir = pathlib.Path(__file__).resolve().parent
    work_dir = pathlib.Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    for name in ["lj32.xyz", "diffusion.run"]:
        shutil.copy2(case_dir / name, work_dir / name)

    run_command([args.gmd, "lj32.xyz", "diffusion.run"], work_dir)
    rows = parse_energy_log(work_dir / "output.log")
    frames = parse_xyz_frames(work_dir / "output.xyz")
    times = [row["time_fs"] for row in rows]
    msd_series = compute_msd_series(frames, [18.0, 18.0, 18.0], times)
    diffusion = fit_diffusion_coefficient(msd_series)

    reference = load_json(case_dir / "reference.json")
    tolerance = load_json(case_dir / "tolerance.json")
    metric = compare_scalar(
        diffusion,
        reference["metrics"]["diffusion_coefficient_a2_per_ps"],
        tolerance["diffusion_coefficient_a2_per_ps_abs"],
    )
    write_series_csv(work_dir / "msd_series.csv", msd_series, ["time_fs", "msd_a2"])
    summary = {
        "case": "diffusion_lj_fluid",
        "metrics": {"diffusion_coefficient_a2_per_ps": metric},
        "passed": metric["passed"],
    }
    write_json(work_dir / "summary.json", summary)
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
