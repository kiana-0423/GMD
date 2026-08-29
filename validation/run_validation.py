from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys


SHORT_CASES = [
    "static_lj_cluster",
    "static_special_pairs",
    "static_coulomb",
    "static_bonded_reference",
    "pme_external",
]

LONG_CASES = [
    "nve_lj_fluid",
    "nvt_lj_fluid",
    "npt_lj_fluid",
    "diffusion_lj_fluid",
]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GMD validation cases.")
    parser.add_argument("--case", choices=SHORT_CASES + LONG_CASES + ["short", "long", "all"], required=True)
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-root", required=True)
    args = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parent
    work_root = pathlib.Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    if args.case == "short":
        cases = SHORT_CASES
    elif args.case == "long":
        cases = LONG_CASES
    elif args.case == "all":
        cases = SHORT_CASES + LONG_CASES
    else:
        cases = [args.case]

    for case in cases:
        case_dir = root / case
        script = case_dir / "analyze.py"
        completed = subprocess.run(
            [
                sys.executable,
                str(script),
                "--gmd",
                args.gmd,
                "--gmd-validate",
                args.gmd_validate,
                "--work-dir",
                str(work_root / case),
            ],
            cwd=root.parent,
            text=True,
        )
        if completed.returncode != 0:
            return completed.returncode

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
