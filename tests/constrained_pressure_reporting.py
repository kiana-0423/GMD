#!/usr/bin/env python3
"""Constrained pressure reporting, through the real `gmd` CLI.

Two things are checked, both of which used to fail.

NVE
    A constrained run with no thermostat and no barostat must show only ordinary
    velocity-Verlet discretisation error in its total energy -- bounded, and
    second order in dt -- and no secular drain from the constraint projection.
    The projection that puts the drifted positions back on the constraint
    manifold has to displace along the constraint gradient at the START of the
    step; displacing along the drifted bond also lands on the manifold but is not
    symplectic and bleeds energy away steadily.

    Two things separate the one from the other, and an absolute bound on the
    drift is neither of them, so neither is used here. First, the drift must fall
    as dt^2 over a fixed span of physical time; a projection drain is first order
    and would fall only as dt. Second, it must be no worse than the SAME system
    run unconstrained, which is the control for how much drift the integrator and
    force field produce on their own.

NPT
    A barostat rescales the cell after the step is complete, and the forces are
    then re-evaluated so the next step starts consistent. That re-evaluation
    produces a virial for a geometry no dynamics integrated, and with constraints
    active it has no constraint partner. It must not be allowed to destroy the
    pressure the completed step measured.

    The specific failure this guards against is a reported pressure that
    silently becomes a column of zeros. A numerical zero is a perfectly ordinary
    pressure, so it must never double as "unavailable": the log carries an
    explicit P_valid column, and an unavailable pressure is written as nan.

MPI
    The MPI output path writes a separate System that holds only the gathered
    global coordinates; every non-atomic field the writer reads has to be brought
    across from the distributed System. Run with --np N this test drives the real
    executable under mpiexec and compares the whole log against a serial run of
    the same fixture, column by column, which is the only thing that covers that
    path end to end.

The sharpest assertion here is that the NPT run and an otherwise identical run
with no barostat report the EXACT same pressure for step 1. The barostat first
acts at the end of step 1, after the completed-step pressure has been captured,
so the two runs are identical up to that point and the reported numbers must
agree bit for bit. They cannot if the reported value is taken from the
post-rescale re-evaluation, or if the constraint term is dropped on the way.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

# A four-atom chain in a non-cubic box, bonds constrained, angles not. Velocities
# are given explicitly and `velocity_init input` is used, so every run here is
# deterministic with no RNG involved.
BOX = (18.0, 21.0, 24.0)
ATOMS = [
    (1, 8.00, 9.00, 10.00, 0.0030, -0.0011, 0.0007),
    (2, 9.42, 9.18, 10.31, -0.0016, 0.0024, -0.0009),
    (2, 10.11, 10.42, 10.94, 0.0021, 0.0008, 0.0013),
    (1, 11.50, 10.27, 11.44, -0.0025, -0.0019, -0.0006),
]
BOND_R0 = 1.500
STEPS = 40
TIME_STEP_FS = 0.5
# The NVE cases run over a fixed span of physical time at two timesteps, so the
# drift can be compared between them rather than against an arbitrary bound.
NVE_TIME_STEP = 0.5
NVE_STEPS = 200

LOG_COLUMNS = [
    "step", "time", "pe", "ke", "etot", "temperature", "pressure", "volume",
    "shake_iter", "shake_error", "rattle_iter", "rattle_error", "p_valid",
]
INT_COLUMNS = {"step", "shake_iter", "rattle_iter", "p_valid"}


def write_inputs(directory: Path, *, barostat: bool, time_step: float = TIME_STEP_FS,
                 steps: int = STEPS, constrain: bool = True) -> None:
    directory.mkdir(parents=True, exist_ok=True)

    lines = [str(len(ATOMS)), " ".join(str(value) for value in BOX)]
    for atom_type, x, y, z, vx, vy, vz in ATOMS:
        lines.append(f"{atom_type}  {x:.10f} {y:.10f} {z:.10f}  {vx:.10f} {vy:.10f} {vz:.10f}")
    (directory / "input.xyz").write_text("\n".join(lines) + "\n")

    (directory / "ff.ff").write_text("\n".join([
        "force_field molecular",
        "lj_cutoff 8.0",
        "type 1 C mass 12.011 epsilon 0.0700 sigma 3.400 charge 0.0",
        "type 2 N mass 14.007 epsilon 0.0715 sigma 3.310 charge 0.0",
        f"bond_type 1 k 300.0 r0 {BOND_R0}",
        "angle_type 1 k 60.0 theta0 112.0",
        "",
    ]))

    (directory / "top.top").write_text("\n".join([
        "bonds 3",
        "  1 2  bond_type 1",
        "  2 3  bond_type 1",
        "  3 4  bond_type 1",
        "",
        "angles 2",
        "  1 2 3  angle_type 1",
        "  2 3 4  angle_type 1",
        "",
        "dihedrals 0",
        "impropers 0",
        "",
    ]))

    run = [
        "velocity 300.0",
        f"time_step {time_step}",
        f"run {steps}",
        "output_interval 1",
        "velocity_init input",
        "remove_com_velocity true",
        "molecular_nonbonded special",
        "",
        # No thermostat: this is NVE, so the energy check means something.
        "constraint_tolerance 1e-12",
        "constraint_max_iterations 500",
    ]
    if constrain:
        run.append("constrain_bond_type 1")
    if barostat:
        run += ["", "barostat berendsen", "pressure 1.0", "barostat_tau 100.0"]
    (directory / "run.in").write_text("\n".join(run) + "\n")


def run_gmd(args: argparse.Namespace, directory: Path, np: int) -> None:
    command: list[str] = []
    if np > 1:
        command += [args.mpiexec, args.mpiexec_np_flag, str(np)]
    command += [args.gmd, "input.xyz", "run.in", "ff.ff", "top.top"]
    if np > 1:
        command += ["--np", str(np)]
    result = subprocess.run(
        command, cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False,
    )
    (directory / "stdout.txt").write_text(result.stdout)
    (directory / "stderr.txt").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"gmd failed ({result.returncode}, np={np}) in {directory}\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )


def read_log(path: Path) -> tuple[list[dict], list[str]]:
    rows: list[dict] = []
    headers: list[str] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            headers.append(line)
            continue
        fields = line.split()
        if len(fields) < len(LOG_COLUMNS):
            raise RuntimeError(
                f"{path}: expected at least {len(LOG_COLUMNS)} columns "
                f"({' '.join(LOG_COLUMNS)}), got {len(fields)}: {line}"
            )
        row = {}
        for name, raw in zip(LOG_COLUMNS, fields):
            row[name] = int(raw) if name in INT_COLUMNS else float(raw)
        rows.append(row)
    return rows, headers


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--mpiexec-np-flag", default="-n")
    # Serial-vs-MPI agreement is limited by the log format, not by the physics:
    # every column is written with six decimals, so two runs that agree to full
    # double precision can still differ by one unit in the last printed place.
    # That is the bound, and it is a statement about %.6f rather than a tolerance
    # on the numbers.
    parser.add_argument("--log-print-resolution", type=float, default=1.0e-6)
    # Absolute, in bar. The two runs are identical through step 1, so the only
    # difference permitted here is decimal rounding in the log itself.
    parser.add_argument("--pressure-match-tolerance", type=float, default=1.0e-9)
    # How much more the constrained run may drift than the same system run
    # unconstrained. Anything above 1 is slack; the constrained run is normally
    # the quieter of the two.
    parser.add_argument("--control-drift-factor", type=float, default=2.0)
    parser.add_argument("--constraint-tolerance", type=float, default=1.0e-9)
    args = parser.parse_args()

    work = Path(args.work_root)
    problems: list[str] = []

    # NVE at two timesteps over the same span of physical time, plus the same
    # system with the constraints switched off as a control.
    nve_dir = work / "constrained_nve"
    nve_half_dir = work / "constrained_nve_half_dt"
    control_dir = work / "unconstrained_nve"
    npt_dir = work / "constrained_npt"

    write_inputs(nve_dir, barostat=False, time_step=NVE_TIME_STEP, steps=NVE_STEPS)
    write_inputs(nve_half_dir, barostat=False,
                 time_step=0.5 * NVE_TIME_STEP, steps=2 * NVE_STEPS)
    write_inputs(control_dir, barostat=False, time_step=NVE_TIME_STEP, steps=NVE_STEPS,
                 constrain=False)
    write_inputs(npt_dir, barostat=True, time_step=TIME_STEP_FS, steps=STEPS)
    for directory in (nve_dir, nve_half_dir, control_dir, npt_dir):
        run_gmd(args, directory, args.np)

    nve, nve_headers = read_log(nve_dir / "output.log")
    nve_half, _ = read_log(nve_half_dir / "output.log")
    control, _ = read_log(control_dir / "output.log")
    npt, _ = read_log(npt_dir / "output.log")

    # --- the log must declare the validity column -------------------------
    if not any("P_valid" in line for line in nve_headers):
        problems.append("the log header does not mention the P_valid column")

    # --- constrained NVE: bounded, second-order, no worse than the control -
    # Every value is checked for finiteness before it reaches a metric. max() over
    # a nan keeps its running value -- `max(0.0, nan)` is 0.0 -- so a column that
    # was never reported would aggregate to "no drift" and pass by vacuity.
    def relative_drift(label: str, rows: list[dict]) -> float:
        for row in rows:
            if not math.isfinite(row["etot"]):
                problems.append(
                    f"{label}: total energy is {row['etot']} at step {row['step']}; a "
                    f"non-finite value cannot be aggregated into a drift metric")
                return float("inf")
        initial = rows[0]["etot"]
        return max(abs(row["etot"] - initial) for row in rows) / max(abs(initial), 1e-30)

    coarse = relative_drift("constrained NVE", nve)
    fine = relative_drift("constrained NVE at half dt", nve_half)
    unconstrained = relative_drift("unconstrained control", control)
    ratio = coarse / fine if fine > 0.0 else float("inf")
    span = NVE_TIME_STEP * NVE_STEPS
    print(f"constrained NVE over {span:.1f} fs: max relative energy drift "
          f"{coarse:.3e} at dt={NVE_TIME_STEP}, {fine:.3e} at dt={0.5 * NVE_TIME_STEP} "
          f"(ratio {ratio:.2f}); unconstrained control {unconstrained:.3e}")

    if not (3.0 <= ratio <= 5.0):
        problems.append(
            f"halving dt changed the constrained NVE energy drift by a factor of "
            f"{ratio:.3f}, not the ~4 of a second-order scheme. A factor near 2 means a "
            f"first-order secular drain, which is what a SHAKE projection that does not "
            f"use the reference gradient produces")
    if coarse > args.control_drift_factor * unconstrained:
        problems.append(
            f"constrained NVE drifts {coarse:.3e} against {unconstrained:.3e} for the same "
            f"system unconstrained, more than the permitted factor of "
            f"{args.control_drift_factor}; the constraints are adding a drain of their own")

    worst_shake = max(row["shake_error"] for row in nve + npt)
    worst_rattle = max(row["rattle_error"] for row in nve + npt)
    if worst_shake > args.constraint_tolerance:
        problems.append(f"SHAKE residual {worst_shake} exceeds {args.constraint_tolerance}")
    if worst_rattle > args.constraint_tolerance:
        problems.append(f"RATTLE residual {worst_rattle} exceeds {args.constraint_tolerance}")

    # --- the premise: the barostat really does rescale, repeatedly --------
    volumes = [row["volume"] for row in npt]
    consecutive = 0
    best_run = 0
    for previous, current in zip(volumes, volumes[1:]):
        if current != previous:
            consecutive += 1
            best_run = max(best_run, consecutive)
        else:
            consecutive = 0
    print(f"constrained NPT: volume {volumes[0]:.6f} -> {volumes[-1]:.6f} A^3, "
          f"longest run of consecutive rescaling steps = {best_run}")
    if best_run < 10:
        problems.append(
            f"the NPT fixture rescaled on at most {best_run} consecutive steps; this test "
            f"has to cover many consecutive rescales to mean anything")

    # --- the pressure column must not silently become zeros ---------------
    stepped = [row for row in npt if row["step"] >= 1]
    zeros = [row["step"] for row in stepped if row["pressure"] == 0.0]
    invalid = [row["step"] for row in stepped if row["p_valid"] != 1]
    nans = [row["step"] for row in stepped if math.isnan(row["pressure"])]
    distinct = {row["pressure"] for row in stepped}
    print(f"constrained NPT: {len(stepped)} stepped frames, "
          f"{len(distinct)} distinct pressures, "
          f"range [{min(r['pressure'] for r in stepped):.3f}, "
          f"{max(r['pressure'] for r in stepped):.3f}] bar, "
          f"{len(zeros)} exact zeros, {len(invalid)} marked invalid")
    if zeros:
        problems.append(
            f"constrained NPT reported an exactly-zero pressure on steps {zeros[:8]}; "
            f"a rescale must not blank the completed step's pressure")
    if invalid:
        problems.append(f"constrained NPT marked the pressure invalid on steps {invalid[:8]}")
    if nans:
        problems.append(f"constrained NPT reported nan on steps {nans[:8]}")
    if len(distinct) < len(stepped) // 2:
        problems.append(
            f"constrained NPT pressure took only {len(distinct)} distinct values over "
            f"{len(stepped)} frames, which looks frozen rather than measured")

    # --- an unavailable pressure is nan and flagged, never zero -----------
    # A constrained run's initial frame has no completed step behind it, and its
    # current-geometry virial has no constraint term, so there is no complete
    # pressure for it. That must show as nan + P_valid 0.
    for name, rows in (("NVE", nve), ("NPT", npt)):
        first = rows[0]
        if first["step"] != 0:
            problems.append(f"{name} log does not start at step 0")
        elif first["p_valid"] != 0:
            problems.append(
                f"{name} initial frame claims a valid pressure, but a constrained run has "
                f"no constraint multiplier before its first step")
        elif not math.isnan(first["pressure"]):
            problems.append(
                f"{name} initial frame reports {first['pressure']} for an unavailable "
                f"pressure; it must be nan, never a number that could be mistaken for one")

    # The control run has no constraints, so its provider virial is complete on
    # its own and even the initial frame reports a pressure. This pins that none
    # of the above changed anything for an unconstrained run.
    if control[0]["p_valid"] != 1 or math.isnan(control[0]["pressure"]):
        problems.append(
            "the unconstrained control's initial frame lost its pressure; nothing here "
            "should affect a run without constraints")

    # --- the decisive one: identical runs must report identical step-1 P ---
    npt_first = next(row for row in npt if row["step"] == 1)
    nve_first = next(row for row in nve if row["step"] == 1)
    delta = abs(npt_first["pressure"] - nve_first["pressure"])
    print(f"step 1 pressure: NPT {npt_first['pressure']:.6f} bar, "
          f"no-barostat {nve_first['pressure']:.6f} bar, difference {delta:.3e}")
    if delta > args.pressure_match_tolerance:
        problems.append(
            f"step-1 pressure differs by {delta} bar between the NPT run and an otherwise "
            f"identical run with no barostat. The barostat first acts at the END of step 1, "
            f"after the completed-step pressure is captured, so these must agree exactly; a "
            f"difference means the reported value came from the post-rescale re-evaluation "
            f"or lost its constraint term")

    # --- serial vs MPI, over the whole log --------------------------------
    if args.np > 1:
        for label, directory in (("NVE", nve_dir), ("NPT", npt_dir)):
            serial_dir = work / f"serial_reference_{label.lower()}"
            write_inputs(serial_dir, barostat=(label == "NPT"),
                         time_step=NVE_TIME_STEP if label == "NVE" else TIME_STEP_FS,
                         steps=NVE_STEPS if label == "NVE" else STEPS)
            run_gmd(args, serial_dir, 1)
            serial_rows, _ = read_log(serial_dir / "output.log")
            mpi_rows, _ = read_log(directory / "output.log")

            if len(serial_rows) != len(mpi_rows):
                problems.append(
                    f"{label}: serial produced {len(serial_rows)} frames, np={args.np} "
                    f"produced {len(mpi_rows)}")
                continue

            worst = {name: 0.0 for name in LOG_COLUMNS}
            scale = {name: 0.0 for name in LOG_COLUMNS}
            for serial_row, mpi_row in zip(serial_rows, mpi_rows):
                for name in LOG_COLUMNS:
                    a, b = serial_row[name], mpi_row[name]
                    if name in INT_COLUMNS:
                        if a != b:
                            problems.append(
                                f"{label} step {serial_row['step']}: {name} is {a} in serial "
                                f"and {b} at np={args.np}; this must match exactly")
                        continue
                    # Finiteness first: a difference involving nan is nan, and
                    # max() would keep its running value and read as agreement.
                    if not math.isfinite(a) or not math.isfinite(b):
                        if math.isnan(a) and math.isnan(b) and name == "pressure" \
                                and serial_row["p_valid"] == 0 and mpi_row["p_valid"] == 0:
                            # Both sides agree there is no pressure for this frame,
                            # which the P_valid column states explicitly.
                            continue
                        problems.append(
                            f"{label} step {serial_row['step']}: {name} is {a} in serial and "
                            f"{b} at np={args.np}; a non-finite value here is either a "
                            f"disagreement or a quantity that was never reported")
                        continue
                    worst[name] = max(worst[name], abs(a - b))
                    # Track the magnitude too: the bound below is expressed in
                    # units of the last printed decimal, and converting those
                    # decimals back to binary is only exact up to a few ULP of
                    # the value itself.
                    scale[name] = max(scale[name], abs(a), abs(b))
            reported = "  ".join(
                f"{name}={worst[name]:.1e}"
                for name in ("pe", "ke", "etot", "temperature", "pressure", "volume",
                             "shake_error", "rattle_error"))
            print(f"{label} serial vs np={args.np} over {len(mpi_rows)} frames, "
                  f"max |difference|: {reported}")
            for name, value in worst.items():
                if name in INT_COLUMNS:
                    continue
                # The bound is "at most one unit in the last printed place",
                # which is a statement about %.6f and not about the numbers. Two
                # decimals one unit apart do not parse back to a difference of
                # exactly 1e-6: neither is representable in binary, so the
                # parsed gap can land a few ULP either side of it. Comparing
                # against the bare resolution therefore rejects the very case
                # the bound is meant to permit -- as it did for a pressure pair
                # printed -1.794485 and -1.794484, whose parsed difference is
                # 1.000000000139778e-06, over by 35% of a single ULP. The slack
                # is scaled by the column's own magnitude, so it stays a
                # statement about decimal conversion rather than a tolerance
                # that grows to hide real disagreement.
                allowed = args.log_print_resolution + 8.0 * sys.float_info.epsilon * scale[name]
                if value > allowed:
                    problems.append(
                        f"{label}: {name} differs by {value} between serial and np={args.np}, "
                        f"more than one unit in the last printed place "
                        f"({args.log_print_resolution}, allowing {allowed - args.log_print_resolution:.1e} "
                        f"for decimal-to-binary conversion at magnitude {scale[name]:.3g}). "
                        f"The MPI output path writes a separate System holding only gathered "
                        f"coordinates; a field that was not carried across from the "
                        f"distributed System shows up here")

    if problems:
        raise RuntimeError("constrained pressure reporting failed:\n  - " +
                           "\n  - ".join(problems))
    print(f"constrained pressure reporting (np={args.np}): all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
