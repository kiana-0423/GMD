"""Berendsen NPT trajectory validation.

WHAT THIS IS. A dynamics regression case, not an external scientific reference.
No other engine was run for it and none of its numbers is derived from theory.
Its value is that it exercises the Berendsen barostat over a real trajectory,
which nothing else in the repository did: npt_lj_fluid uses the Monte Carlo
barostat, and the Berendsen path had only direct unit tests behind it.

WHAT IT IS SENSITIVE TO. The case is a sign-sensitive pair. Two runs differ
only in their target pressure, both on a compressed 32-atom Lennard-Jones
fixture whose mean pressure over the run is roughly 1000-1500 bar:

    expand.run    target  200 bar, BELOW the fixture's mean pressure -> the
                  cell must grow
    compress.run  target 4000 bar, ABOVE it -> the cell must shrink

Two production defects are caught by construction rather than by tolerance:

  * A barostat comparing a bar target against an eV/A^3 pressure -- the defect
    corrected in "fix: convert the Berendsen barostat's pressure comparison to
    bar" -- reads 200 bar as 200 eV/A^3, which is 3.2e8 bar. Both targets then
    sit enormously above any instantaneous pressure and BOTH runs compress. The
    expansion assertion fails outright, whatever the tolerances say.

  * A relaxation time consumed in the wrong unit -- the defect corrected in
    "fix: convert thermostat and barostat relaxation times to internal units" --
    leaves the coupling about ten times too weak, so the cell barely moves. The
    response-magnitude assertions below are what catch that; they are set from
    the measured response of the correct code and of the defective code, not
    from a round number.

The reference values are regression values. The sign and magnitude assertions
are not: they hold for any correct barostat and are the reason this case exists.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import shutil
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import compare_scalar, load_json, mean, parse_energy_log, run_command, write_json, write_series_csv

CASE = "berendsen_npt_lj"
FIXTURE = "lj32_compressed.xyz"


def finite(values) -> bool:
    return all(math.isfinite(v) for v in values)


def series(rows, key):
    return [row[key] for row in rows]


def run_case(gmd: str, work_dir: pathlib.Path, run_file: str, log_name: str,
             mpiexec: str | None = None, np: int = 1) -> list[dict]:
    command: list[str] = []
    if np > 1:
        if not mpiexec:
            raise RuntimeError("MPI comparison requires --mpiexec")
        command += [mpiexec, "-n", str(np)]
    command += [gmd, FIXTURE, run_file]
    if np > 1:
        command += ["--np", str(np)]
    run_command(command, work_dir)
    produced = work_dir / "output.log"
    target = work_dir / log_name
    if target.exists():
        target.unlink()
    produced.rename(target)
    return parse_energy_log(target)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=False, default="")
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--mpiexec", default="",
                        help="if given, the expand run is repeated under MPI and "
                             "the two trajectories are compared")
    parser.add_argument("--np", type=int, default=2)
    args = parser.parse_args()

    case_dir = pathlib.Path(__file__).resolve().parent
    work_dir = pathlib.Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    for name in [FIXTURE, "expand.run", "compress.run",
                 "restart_first.run", "restart_second.run"]:
        shutil.copy2(case_dir / name, work_dir / name)

    expand = run_case(args.gmd, work_dir, "expand.run", "expand.log")
    compress = run_case(args.gmd, work_dir, "compress.run", "compress.log")

    write_series_csv(work_dir / "expand_series.csv", expand,
                     ["step", "time_fs", "pe", "ke", "total_energy",
                      "temperature", "pressure_bar", "volume_a3"])
    write_series_csv(work_dir / "compress_series.csv", compress,
                     ["step", "time_fs", "pe", "ke", "total_energy",
                      "temperature", "pressure_bar", "volume_a3"])

    reference = load_json(case_dir / "reference.json")
    tolerance = load_json(case_dir / "tolerance.json")

    initial_volume = expand[0]["volume_a3"]
    expand_ratio = expand[-1]["volume_a3"] / initial_volume
    compress_ratio = compress[-1]["volume_a3"] / initial_volume

    # --- structural assertions: true of any correct barostat -----------------
    structural: dict[str, dict] = {}

    def assert_true(name: str, condition: bool, detail: str) -> None:
        structural[name] = {"passed": bool(condition), "detail": detail}

    assert_true(
        "expand_direction", expand_ratio > 1.0,
        f"target 200 bar is below the fixture's mean pressure "
        f"({mean(series(expand, 'pressure_bar')):.1f} bar), so the cell must grow; "
        f"volume ratio was {expand_ratio!r}")
    assert_true(
        "compress_direction", compress_ratio < 1.0,
        f"target 4000 bar is above the fixture's mean pressure "
        f"({mean(series(compress, 'pressure_bar')):.1f} bar), so the cell must "
        f"shrink; volume ratio was {compress_ratio!r}")
    assert_true(
        "directions_differ", expand_ratio > 1.0 > compress_ratio,
        f"the two targets must move the cell in opposite directions; got "
        f"{expand_ratio!r} and {compress_ratio!r}")

    # Response magnitude. A coupling ten times too weak leaves the expansion at
    # about 1.01 instead of 1.11; the bound sits between the two.
    min_response = tolerance["min_volume_response"]
    assert_true(
        "expand_magnitude", abs(expand_ratio - 1.0) > min_response,
        f"the cell expanded by only {abs(expand_ratio - 1.0):.4f}, below the "
        f"{min_response} a correctly coupled barostat produces; a relaxation "
        f"time consumed in the wrong unit looks like this")
    assert_true(
        "compress_magnitude", abs(compress_ratio - 1.0) > min_response,
        f"the cell compressed by only {abs(compress_ratio - 1.0):.4f}, below "
        f"{min_response}")

    # --- sanity: nothing diverged --------------------------------------------
    for label, rows in (("expand", expand), ("compress", compress)):
        volumes = series(rows, "volume_a3")
        assert_true(f"{label}_volume_positive", all(v > 0.0 for v in volumes),
                    f"{label}: a non-positive volume appeared")
        assert_true(f"{label}_volume_finite", finite(volumes),
                    f"{label}: a non-finite volume appeared")
        assert_true(
            f"{label}_volume_bounded",
            all(0.25 * initial_volume < v < 4.0 * initial_volume for v in volumes),
            f"{label}: the cell left the range [0.25, 4.0] x initial volume, "
            f"which is a runaway rather than coupling")
        for key in ("pe", "ke", "total_energy", "temperature", "pressure_bar"):
            assert_true(f"{label}_{key}_finite", finite(series(rows, key)),
                        f"{label}: a non-finite {key} appeared")
        assert_true(f"{label}_temperature_positive",
                    all(t >= 0.0 for t in series(rows, "temperature")),
                    f"{label}: a negative temperature appeared")

    # --- restart continuity ---------------------------------------------------
    run_case(args.gmd, work_dir, "restart_first.run", "restart_first.log")
    split = run_case(args.gmd, work_dir, "restart_second.run", "restart_second.log")
    continuous_final = expand[-1]["volume_a3"]
    split_final = split[-1]["volume_a3"]
    assert_true(
        "restart_continuity", split_final == continuous_final,
        f"a run split across a checkpoint ended at volume {split_final!r} where "
        f"the continuous run ended at {continuous_final!r}; these must be "
        f"bit-identical, since the split run resumes the same trajectory")

    # --- serial/MPI: reported, and deliberately not asserted ------------------
    #
    # A trajectory started from `velocity_init random` does not reproduce across
    # rank counts, and the reason is not the barostat. VelocityInitializer draws
    # from one generator sequentially by LOCAL atom index, so each rank gives its
    # own first atom the generator's first three draws; under decomposition that
    # is a different physical atom than in a serial run. The total kinetic energy
    # is then rescaled to the requested temperature, which is why step 0 reports
    # an identical temperature and an identical potential energy while the
    # underlying velocity field differs.
    #
    # It is a pre-existing limitation, not something this case introduced, and it
    # is demonstrable without any barostat at all: the same fixture run as plain
    # NVE diverges between np=1 and np=2 in exactly the same way, while the same
    # fixture run with `velocity 0.0` stays bit-identical for 50 steps.
    #
    # So this reports the difference rather than asserting on it. What IS
    # asserted about the barostat under MPI lives in tests/mpi_berendsen_barostat.cpp,
    # which gives every rank the same velocity field by construction and requires
    # the coupling factor to be bit-identical at np=1, 2 and 4.
    mpi_detail = "not run (no --mpiexec given)"
    if args.mpiexec:
        mpi_rows = run_case(args.gmd, work_dir, "expand.run", "expand_mpi.log",
                            mpiexec=args.mpiexec, np=args.np)
        mpi_final = mpi_rows[-1]["volume_a3"]
        difference = abs(mpi_final - continuous_final)
        mpi_detail = (f"np={args.np} ended at {mpi_final!r} against the serial "
                      f"{continuous_final!r} (difference {difference:.3e}); reported "
                      f"only -- see the note in analyze.py on rank-dependent "
                      f"random velocity initialisation")
        # The one thing that must hold: both runs stayed physical.
        assert_true("mpi_volume_bounded",
                    all(0.25 * initial_volume < r["volume_a3"] < 4.0 * initial_volume
                        for r in mpi_rows),
                    f"the MPI run left the physical volume range: {mpi_detail}")
    structural["serial_mpi"] = {"passed": True, "detail": mpi_detail}

    # --- regression metrics ---------------------------------------------------
    metrics = {
        "expand_volume_ratio": compare_scalar(
            expand_ratio, reference["metrics"]["expand_volume_ratio"],
            tolerance["volume_ratio_abs"]),
        "compress_volume_ratio": compare_scalar(
            compress_ratio, reference["metrics"]["compress_volume_ratio"],
            tolerance["volume_ratio_abs"]),
        # Means exclude frame 0. That frame is the untouched input lattice
        # before any dynamics, and its pressure -- about 22600 bar against a
        # run mean near 900 -- dominates both the mean and its spread without
        # describing anything the barostat did.
        "expand_temperature_mean": compare_scalar(
            mean(series(expand[1:], "temperature")),
            reference["metrics"]["expand_temperature_mean"],
            tolerance["temperature_mean_abs"]),
        "expand_pressure_mean": compare_scalar(
            mean(series(expand[1:], "pressure_bar")),
            reference["metrics"]["expand_pressure_mean"],
            tolerance["pressure_mean_abs"]),
    }

    passed = (all(m["passed"] for m in metrics.values())
              and all(s["passed"] for s in structural.values()))
    summary = {
        "case": CASE,
        "metrics": metrics,
        "structural": structural,
        "passed": passed,
    }
    write_json(work_dir / "summary.json", summary)
    if not passed:
        for name, entry in structural.items():
            if not entry["passed"]:
                print(f"[{CASE}] {name}: {entry['detail']}", file=sys.stderr)
        for name, entry in metrics.items():
            if not entry["passed"]:
                print(f"[{CASE}] {name}: actual={entry['actual']!r} "
                      f"expected={entry['expected']!r} tolerance={entry['tolerance']!r}",
                      file=sys.stderr)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
