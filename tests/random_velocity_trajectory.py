#!/usr/bin/env python3
"""Serial/MPI trajectory equivalence for runs that start from random velocities.

tests/velocity_init_equivalence.py compares the initial velocity field. This
compares what the field is FOR: a real trajectory, driven through the gmd
executable, at np=1, 2 and 4 plus a plain serial run with no mpiexec at all.

Why this is a regression test and not a smoke test. Before the initializer was
keyed on the atom tag, every rank handed its own first local atom the seeded
generator's first draws, so the runs below started from different velocity
fields and diverged immediately -- the potential energy differed in the sixth
decimal by step 10. The comparisons here are exact, so that failure is
unmissable.

What is compared:

  * every column of the energy log, as printed text, for every frame;
  * the final checkpoint's per-atom coordinates AND velocities, keyed by TAG
    rather than by line, since the gathered ordering is not the input ordering;
  * restart continuity: a run split across a checkpoint must reach the same
    final state as the continuous one, at np=1 and under MPI.

Exact comparison is the right bar here rather than a tolerance. The log is
six-decimal text and the checkpoint is written at full precision, and the
preceding work established that this fixture's force path is bit-identical
across rank counts: with `velocity 0.0` it reproduces exactly for fifty steps.
The one thing that was not identical was the initial field, which is what this
change fixes. If a future change makes exact agreement impossible, that is a
finding to explain rather than a tolerance to widen.
"""

from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import sys

problems: list[str] = []

STEPS = 50
HALF = 25
TIME_STEP = 2.0
SEED = 20260523


def check(condition: bool, message: str) -> None:
    if not condition:
        problems.append(message)


def write_fixture(directory: pathlib.Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    # 32 atoms on a slightly irregular lattice in an 18 A box: dense enough to
    # interact through the 8.5 A cutoff, loose enough not to blow up.
    lines = ["32", "18.0 18.0 18.0"]
    index = 0
    for i in range(4):
        for j in range(4):
            for k in range(2):
                x = 2.0 + 4.0 * i + 0.13 * ((index * 7) % 5)
                y = 2.0 + 4.0 * j + 0.11 * ((index * 3) % 7)
                z = 3.0 + 8.0 * k + 0.17 * ((index * 5) % 3)
                lines.append(f"1  {x:.6f}  {y:.6f}  {z:.6f}")
                index += 1
    (directory / "lj32.xyz").write_text("\n".join(lines) + "\n")


def run_file(steps: int, *, checkpoint: str, restart_from: str | None,
             velocity_init: str) -> str:
    lines = [
        "velocity 120.0",
        f"time_step {TIME_STEP}",
        f"run {steps}",
        "output_interval 10",
        f"velocity_init {velocity_init}",
        f"velocity_seed {SEED}",
        "remove_com_velocity true",
        "",
        "force_field lj",
        "cutoff 8.5",
        "type 1 Ar epsilon 0.01032 sigma 3.405 charge 0.0",
        "",
        f"checkpoint_file {checkpoint}",
        "write_checkpoint_every 0",
    ]
    if restart_from:
        lines.append(f"restart_from {restart_from}")
    return "\n".join(lines) + "\n"


def run_gmd(args, directory: pathlib.Path, run_name: str, np_count: int | None) -> None:
    command: list[str] = []
    if np_count is not None:
        command += [args.mpiexec, args.mpiexec_np_flag, str(np_count)]
    command += [args.gmd, "lj32.xyz", run_name]
    if np_count is not None and np_count > 1:
        command += ["--np", str(np_count)]
    result = subprocess.run(command, cwd=directory, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    (directory / "stdout.txt").write_text(result.stdout)
    (directory / "stderr.txt").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"gmd failed ({result.returncode}) in {directory}\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}")


def read_log_lines(path: pathlib.Path) -> list[list[str]]:
    return [line.split() for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]


def read_checkpoint_atoms(path: pathlib.Path) -> dict[int, tuple[float, ...]]:
    """tag -> (x, y, z, vx, vy, vz)."""
    atoms: dict[int, tuple[float, ...]] = {}
    reading = False
    for line in path.read_text().splitlines():
        if line.startswith("atoms tag type molecule mass charge x y z vx vy vz"):
            reading = True
            continue
        if not reading:
            continue
        fields = line.split()
        # tag type molecule mass charge x y z vx vy vz
        if len(fields) < 11:
            break
        tag = int(fields[0])
        if tag in atoms:
            problems.append(f"{path.name}: atom tag {tag} appears more than once")
        atoms[tag] = tuple(float(value) for value in fields[5:11])
    return atoms


def compare_logs(label: str, reference: list[list[str]], other: list[list[str]]) -> None:
    check(len(reference) == len(other),
          f"{label}: log has {len(other)} frames, the serial run has {len(reference)}")
    for row_index, (a, b) in enumerate(zip(reference, other)):
        check(a == b,
              f"{label}: log row {row_index} differs from serial.\n"
              f"      serial: {' '.join(a)}\n"
              f"      {label:<6s}: {' '.join(b)}")
        if a != b:
            return          # one diff is enough; the rest are consequences


def compare_atoms(label: str, reference: dict, other: dict, tolerance: float) -> float:
    """Compares by TAG. Returns the worst absolute difference seen."""
    check(set(reference) == set(other),
          f"{label}: checkpoint atom tags differ from serial "
          f"(missing {sorted(set(reference) - set(other))}, "
          f"extra {sorted(set(other) - set(reference))})")
    worst, worst_tag = 0.0, -1
    for tag in sorted(set(reference) & set(other)):
        for a, b in zip(reference[tag], other[tag]):
            if abs(a - b) > worst:
                worst, worst_tag = abs(a - b), tag
    if worst > tolerance:
        check(False,
              f"{label}: the final state differs from serial by {worst:.6e} at atom "
              f"tag {worst_tag}, above the {tolerance:.1e} bound.\n"
              f"      serial: {reference[worst_tag]}\n"
              f"      {label}: {other[worst_tag]}")
    return worst


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--mpiexec", default="")
    parser.add_argument("--mpiexec-np-flag", default="-n")
    parser.add_argument("--max-np", type=int, default=4)
    # Graded, and every bound is measured rather than chosen.
    #
    # The LOG is compared exactly at every rank count, as printed text.
    #
    # The CHECKPOINT is full precision, so it exposes what the log's six
    # decimals hide. A plain serial run and an np=1 run agree BITWISE, and so
    # does a serial run split across a checkpoint -- there is nothing for a
    # reduction to reorder in either. At np>=2 the centre-of-mass and
    # kinetic-energy sums are combined by MPI_Allreduce over a different tree,
    # which moves the initial field by ~2.8e-17; fifty steps of chaotic
    # dynamics amplify that to a measured worst case of 8.88e-16, about 28 ulp.
    # 1e-13 is ~113x that. It is not slack for the random draws, which are
    # identical by tag; before this fix the runs differed in the sixth printed
    # decimal of the potential energy by step 10, which the exact log
    # comparison catches outright.
    parser.add_argument("--mpi-state-tolerance", type=float, default=1.0e-13)
    args = parser.parse_args()

    work = pathlib.Path(args.work_root).resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    # --- the reference: a plain serial run, no mpiexec involved ------------
    serial_dir = work / "serial"
    write_fixture(serial_dir)
    (serial_dir / "run.in").write_text(
        run_file(STEPS, checkpoint="final.chk", restart_from=None, velocity_init="random"))
    run_gmd(args, serial_dir, "run.in", None)
    reference_log = read_log_lines(serial_dir / "output.log")
    reference_atoms = read_checkpoint_atoms(serial_dir / "final.chk")
    check(len(reference_atoms) == 32,
          f"the serial checkpoint holds {len(reference_atoms)} atoms, expected 32")
    check(len(reference_log) > 1,
          "the serial run produced no trajectory")

    arrangements = []
    if args.mpiexec:
        arrangements = [n for n in (1, 2, 4) if n <= args.max_np]

    for np_count in arrangements:
        label = f"np{np_count}"
        directory = work / label
        write_fixture(directory)
        (directory / "run.in").write_text(
            run_file(STEPS, checkpoint="final.chk", restart_from=None,
                     velocity_init="random"))
        run_gmd(args, directory, "run.in", np_count)
        compare_logs(label, reference_log, read_log_lines(directory / "output.log"))
        # np=1 has no cross-rank reduction to reorder, so it must be bitwise.
        tolerance = 0.0 if np_count == 1 else args.mpi_state_tolerance
        worst = compare_atoms(label, reference_atoms,
                              read_checkpoint_atoms(directory / "final.chk"), tolerance)
        print(f"  {label:<8s} log identical, final state within {worst:.3e}"
              + ("  (bitwise)" if worst == 0.0 else ""))

    # --- restart continuity ------------------------------------------------
    # A run split across a checkpoint must land where the continuous one did.
    # `run` counts additional steps, so the second half asks for HALF more.
    restart_counts = [None] + ([2] if args.mpiexec else [])
    for np_count in restart_counts:
        label = "serial" if np_count is None else f"np{np_count}"
        directory = work / f"restart_{label}"
        write_fixture(directory)
        (directory / "first.in").write_text(
            run_file(HALF, checkpoint="half.chk", restart_from=None,
                     velocity_init="random"))
        run_gmd(args, directory, "first.in", np_count)
        (directory / "second.in").write_text(
            run_file(STEPS - HALF, checkpoint="split.chk", restart_from="half.chk",
                     velocity_init="input"))
        run_gmd(args, directory, "second.in", np_count)
        tolerance = 0.0 if np_count is None else args.mpi_state_tolerance
        worst = compare_atoms(f"restart({label})", reference_atoms,
                              read_checkpoint_atoms(directory / "split.chk"), tolerance)
        print(f"  restart({label}) final state within {worst:.3e}"
              + ("  (bitwise)" if worst == 0.0 else ""))

    if problems:
        for message in problems:
            print(f"[random velocity trajectory] {message}", file=sys.stderr)
        print(f"[random velocity trajectory] {len(problems)} problem(s)", file=sys.stderr)
        return 1

    covered = ", ".join(["serial"] + [f"np={n}" for n in arrangements]) or "serial"
    print(f"[random velocity trajectory] identical trajectories and final states "
          f"across {covered}, and across a checkpoint split")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
