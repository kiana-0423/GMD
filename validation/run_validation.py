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
    "berendsen_npt_lj",
    "diffusion_lj_fluid",
]

# Constrained dynamics through the real CLI. These are grouped separately
# because they accept the extra MPI arguments below: the constraint solver
# replicates its state across ranks, so serial/MPI agreement is part of what
# they validate rather than a separate concern.
CONSTRAINED_CASES = [
    "constrained_nve_water",
    "constrained_nvt_cluster",
]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GMD validation cases.")
    parser.add_argument(
        "--case",
        choices=SHORT_CASES + LONG_CASES + CONSTRAINED_CASES
        + ["short", "long", "constrained", "all"],
        required=True,
    )
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-root", required=True)
    # Only the constrained cases understand these. They are optional: without
    # them those cases run serially and record their MPI check as skipped.
    parser.add_argument("--mpiexec", default="",
                        help="mpiexec/mpirun to use for the constrained cases")
    parser.add_argument("--mpi-ranks", default="",
                        help="comma-separated rank counts, e.g. 1,2,4")
    args = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parent
    work_root = pathlib.Path(args.work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    if args.case == "short":
        cases = SHORT_CASES
    elif args.case == "long":
        cases = LONG_CASES
    elif args.case == "constrained":
        cases = CONSTRAINED_CASES
    elif args.case == "all":
        cases = SHORT_CASES + LONG_CASES + CONSTRAINED_CASES
    else:
        cases = [args.case]

    for case in cases:
        case_dir = root / case
        script = case_dir / "analyze.py"
        command = [
            sys.executable,
            str(script),
            "--gmd",
            args.gmd,
            "--gmd-validate",
            args.gmd_validate,
            "--work-dir",
            str(work_root / case),
        ]
        if case in CONSTRAINED_CASES and args.mpiexec and args.mpi_ranks:
            command += ["--mpiexec", args.mpiexec, "--mpi-ranks", args.mpi_ranks]
        completed = subprocess.run(command, cwd=root.parent, text=True)
        if completed.returncode != 0:
            return completed.returncode

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
