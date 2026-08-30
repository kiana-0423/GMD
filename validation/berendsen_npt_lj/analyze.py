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



def read_raw_log(path: pathlib.Path) -> list[list[str]]:
    """Every data row of the log as printed text, so comparison is exact."""
    return [line.split() for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]


LOG_COLUMNS = ["step", "time_fs", "pe", "ke", "total_energy", "temperature",
               "pressure_bar", "volume_a3", "shake_iter", "shake_error",
               "rattle_iter", "rattle_error", "p_valid"]


def first_log_difference(reference, other):
    """The first differing field, named by column and step, or None."""
    for row_index, (a, b) in enumerate(zip(reference, other)):
        for column_index, (x, y) in enumerate(zip(a, b)):
            if x != y:
                column = (LOG_COLUMNS[column_index]
                          if column_index < len(LOG_COLUMNS)
                          else f"column {column_index}")
                return (f"log column '{column}' differs at step {a[0]} "
                        f"(row {row_index}): serial {x}, MPI {y}")
    return None


def read_checkpoint_atoms(path: pathlib.Path) -> dict:
    """tag -> (x, y, z, vx, vy, vz) as printed text."""
    atoms: dict = {}
    reading = False
    if not path.exists():
        return atoms
    for line in path.read_text().splitlines():
        if line.startswith("atoms tag type molecule mass charge x y z vx vy vz"):
            reading = True
            continue
        if not reading:
            continue
        fields = line.split()
        if len(fields) < 11:
            break
        atoms[int(fields[0])] = tuple(fields[5:11])
    return atoms


ATOM_FIELDS = ["x", "y", "z", "vx", "vy", "vz"]


def first_atom_difference(reference, other, tolerance):
    """The worst per-atom difference above `tolerance`, named by tag and
    component, or None. `tolerance` of zero demands bitwise equality."""
    if set(reference) != set(other):
        return (f"checkpoint atom tags differ "
                f"(missing {sorted(set(reference) - set(other))}, "
                f"extra {sorted(set(other) - set(reference))})")
    worst, detail = 0.0, None
    for tag in sorted(reference):
        for index, (x, y) in enumerate(zip(reference[tag], other[tag])):
            if tolerance == 0.0:
                if x != y:
                    return (f"atom tag {tag} component '{ATOM_FIELDS[index]}' "
                            f"differs: serial {x}, MPI {y}")
                continue
            difference = abs(float(x) - float(y))
            if difference > worst:
                worst = difference
                detail = (f"atom tag {tag} component '{ATOM_FIELDS[index]}' differs "
                          f"by {difference:.6e}: serial {x}, MPI {y}")
    if tolerance > 0.0 and worst > tolerance:
        return f"{detail} (bound {tolerance:.1e})"
    return None

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
    # expand.run writes expand_final.chk. Every arrangement runs that same input
    # file, so the checkpoint has to be moved aside with its log or the next run
    # overwrites the state this one is about to be compared against. Only
    # expand.run writes it -- the other run files have their own checkpoint
    # names -- so the rename is scoped to it, or a compress run would carry the
    # serial expand state away with it.
    checkpoint = work_dir / "expand_final.chk"
    if run_file == "expand.run" and checkpoint.exists():
        stem = log_name[:-len(".log")] if log_name.endswith(".log") else log_name
        renamed = work_dir / f"{stem}_final.chk"
        # The serial run's log is expand.log, so its checkpoint already has the
        # name it wants; renaming it onto itself would delete it.
        if renamed != checkpoint:
            if renamed.exists():
                renamed.unlink()
            checkpoint.rename(renamed)
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

    # --- serial/MPI: enforced --------------------------------------------------
    #
    # This was report-only. It was written when a run started from
    # `velocity_init random` genuinely did not reproduce across rank counts --
    # VelocityInitializer drew from one generator in local storage order, so
    # every rank handed its own first local atom the seeded generator's first
    # draws. That is fixed: draws are now keyed on the global atom tag, and the
    # trajectories agree BIT FOR BIT.
    #
    # So the comparison is exact, on the printed text of every column of every
    # frame, and on the final per-atom state keyed by TAG. Exact is the right
    # bar because it is what the code achieves; a tolerance here would only hide
    # a future regression. If exact agreement ever becomes unavailable that is a
    # finding to explain, not a number to widen.
    mpi_detail = "not run (no --mpiexec given)"
    if args.mpiexec:
        reference_rows = read_raw_log(work_dir / "expand.log")
        reference_atoms = read_checkpoint_atoms(work_dir / "expand_final.chk")
        assert_true("serial_checkpoint_present", bool(reference_atoms),
                    "the serial expand run wrote no comparable checkpoint")
        compared = []
        for np_count in (1, 2, 4):
            label = f"np{np_count}"
            rows = run_case(args.gmd, work_dir, "expand.run", f"expand_{label}.log",
                            mpiexec=args.mpiexec, np=np_count)
            # Structural first: a NaN or a missing frame must be reported as
            # itself rather than folded into a difference metric.
            assert_true(f"{label}_frames",
                        len(rows) == len(expand),
                        f"{label}: {len(rows)} frames against the serial run's "
                        f"{len(expand)}")
            finite_columns = all(
                finite(series(rows, key))
                for key in ("pe", "ke", "total_energy", "temperature",
                            "pressure_bar", "volume_a3"))
            assert_true(f"{label}_finite", finite_columns,
                        f"{label}: a non-finite value appeared in the log")
            if len(rows) != len(expand) or not finite_columns:
                continue

            raw = read_raw_log(work_dir / f"expand_{label}.log")
            difference = first_log_difference(reference_rows, raw)
            assert_true(f"{label}_log_identical", difference is None,
                        f"{label}: {difference}" if difference else "")

            # np=1 has no cross-rank reduction to reorder, so the final state
            # must be BITWISE identical. Beyond that the Berendsen coupling
            # factor is built from globally reduced quantities whose summation
            # order changes with the rank count, and 500 steps of chaotic
            # dynamics amplify that: measured worst case 1.7e-14 at np=2 and
            # 2.1e-14 at np=4, on coordinates of order 1 to 15 Angstrom. The
            # bound is roughly fifty times that. It is not slack for the
            # trajectory -- the LOG is compared exactly at every rank count, and
            # every failure mode this case is built to catch shows up there.
            atoms = read_checkpoint_atoms(work_dir / f"expand_{label}_final.chk")
            state_tolerance = 0.0 if np_count == 1 else tolerance["serial_mpi_state_abs"]
            atom_difference = first_atom_difference(reference_atoms, atoms,
                                                    state_tolerance)
            assert_true(f"{label}_state_identical", atom_difference is None,
                        f"{label}: {atom_difference}" if atom_difference else "")

            # The property the case exists for must survive decomposition too.
            ratio = rows[-1]["volume_a3"] / rows[0]["volume_a3"]
            assert_true(f"{label}_expand_direction", ratio > 1.0,
                        f"{label}: the cell did not expand (volume ratio {ratio!r})")
            compared.append(label)
        mpi_detail = ("identical logs and final states at " + ", ".join(compared)
                      if compared else "no rank count completed")
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
