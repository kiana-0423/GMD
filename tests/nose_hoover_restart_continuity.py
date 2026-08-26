#!/usr/bin/env python3
"""End-to-end CLI checkpoint/restart continuity for a CONSTRAINED Nose-Hoover run.

This drives the real `gmd` executable through the whole application path --

    run-file / force-field / topology parsing
      -> Simulation::initialize()
      -> checkpoint loading
      -> Nose-Hoover thermostat-state validation
      -> continued integration

-- rather than poking the thermostat's checkpoint methods directly.

  continuous : xyz + run -> N_TOTAL steps
  split      : xyz + run -> N_FIRST steps -> checkpoint -> restart -> rest

The two must agree on coordinates, velocities, forces, energies, temperature,
pressure, the full Nose-Hoover extended-system state (tau, xi, Q, dof, cached
temperature), step number, simulation time and constraint residuals.

Two negative cases are also driven through the CLI: restarting with a changed
centre-of-mass removal setting, and with a changed constraint configuration.
Both must be rejected before a single step is taken.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path


# --- Fixture ---------------------------------------------------------------
#
# A four-atom chain in a non-cubic box. Bonds are constrained (SHAKE/RATTLE),
# angles are not, so the run exercises constraint-aware degrees of freedom:
#
#   DOF = 3N - 3 (COM removed) - 3 (constraints) = 12 - 3 - 3 = 6
#
# Velocities are given explicitly and `velocity_init input` is used, so the run
# is fully deterministic with no RNG involved.

BOX = (18.0, 21.0, 24.0)

ATOMS = [
    # type, x, y, z, vx, vy, vz
    (1, 8.00, 9.00, 10.00, 0.0030, -0.0011, 0.0007),
    (2, 9.42, 9.18, 10.31, -0.0016, 0.0024, -0.0009),
    (2, 10.11, 10.42, 10.94, 0.0021, 0.0008, 0.0013),
    (1, 11.50, 10.27, 11.44, -0.0025, -0.0019, -0.0006),
]

BOND_R0 = 1.500          # constrained target distance [A]
EXPECTED_DOF = 6
N_TOTAL = 60
N_FIRST = 30
TIME_STEP_FS = 0.5
TARGET_TEMPERATURE = 300.0   # [K]; must be > 0 or the thermostat is inert


def write_xyz(path: Path) -> None:
    lines = [str(len(ATOMS)), " ".join(str(value) for value in BOX)]
    for atom_type, x, y, z, vx, vy, vz in ATOMS:
        lines.append(f"{atom_type}  {x:.10f} {y:.10f} {z:.10f}  {vx:.10f} {vy:.10f} {vz:.10f}")
    path.write_text("\n".join(lines) + "\n")


def write_ff(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "force_field molecular",
                "lj_cutoff 8.0",
                "type 1 C mass 12.011 epsilon 0.0700 sigma 3.400 charge 0.0",
                "type 2 N mass 14.007 epsilon 0.0715 sigma 3.310 charge 0.0",
                f"bond_type 1 k 300.0 r0 {BOND_R0}",
                "angle_type 1 k 60.0 theta0 112.0",
                "",
            ]
        )
    )


def write_top(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
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
            ]
        )
    )


def write_run(
    path: Path,
    *,
    steps: int,
    checkpoint_file: Path,
    restart_from: Path | None = None,
    remove_com: bool = True,
    constrain: bool = True,
) -> None:
    lines = [
        # A positive target temperature is what makes the thermostat live:
        # apply_half_kick() returns immediately when it is zero, which would
        # leave xi and Q at zero and make this whole comparison vacuous.
        f"velocity {TARGET_TEMPERATURE}",
        f"time_step {TIME_STEP_FS}",
        f"run {steps}",
        "output_interval 1",
        "velocity_init input",
        f"remove_com_velocity {'true' if remove_com else 'false'}",
        "molecular_nonbonded special",
        "",
        "thermostat nose_hoover",
        "thermostat_tau 100.0",
        "",
        "constraint_tolerance 1e-10",
        "constraint_max_iterations 200",
    ]
    if constrain:
        lines.append("constrain_bond_type 1")
    lines.append(f"checkpoint_file {checkpoint_file}")
    lines.append("write_checkpoint_every 0")
    if restart_from is not None:
        lines.append(f"restart_from {restart_from}")
    path.write_text("\n".join(lines) + "\n")


# --- Running ---------------------------------------------------------------

def run_gmd(args: argparse.Namespace, cwd: Path, check: bool = True):
    command = []
    if args.np > 1:
        command.extend([args.mpiexec, args.mpiexec_np_flag, str(args.np)])
    command.extend([args.gmd, "input.xyz", "run.in", "ff.ff", "top.top"])
    if args.np > 1:
        command.extend(["--np", str(args.np)])

    result = subprocess.run(
        command, cwd=cwd, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    (cwd / "stdout.txt").write_text(result.stdout)
    (cwd / "stderr.txt").write_text(result.stderr)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"gmd failed ({result.returncode}) in {cwd}\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )
    return result


def prepare(directory: Path, **run_kwargs) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    write_xyz(directory / "input.xyz")
    write_ff(directory / "ff.ff")
    write_top(directory / "top.top")
    write_run(directory / "run.in", **run_kwargs)


# --- Parsing ---------------------------------------------------------------

def read_checkpoint(path: Path) -> dict:
    lines = path.read_text().splitlines()
    if not lines or not lines[0].startswith("GMD_CHECKPOINT"):
        raise RuntimeError(f"{path} is not a GMD checkpoint")

    data: dict = {"atoms": {}}
    index = 1
    while index < len(lines):
        line = lines[index]
        if line.startswith("atoms tag"):
            count = int(data["atom_count"])
            for atom_line in lines[index + 1 : index + 1 + count]:
                f = atom_line.split()
                data["atoms"][int(f[0])] = {
                    "pos": [float(f[5]), float(f[6]), float(f[7])],
                    "vel": [float(f[8]), float(f[9]), float(f[10])],
                    "mass": float(f[3]),
                }
            index += count + 1
            continue
        if line == "end":
            break
        parts = line.split()
        if parts:
            key = parts[0]
            if key in {"step", "atom_count"}:
                data[key] = int(parts[1])
            elif key == "time_fs":
                data[key] = float(parts[1])
            elif key == "box":
                data[key] = [float(parts[1]), float(parts[2]), float(parts[3])]
            elif key == "thermostat_state":
                data[key] = " ".join(parts[1:])
            elif key == "thermostat_type":
                data[key] = " ".join(parts[1:])
        index += 1
    return data


def parse_thermostat_state(state: str) -> dict:
    """`tau <v> xi <v> Q <v> dof <v> current_temperature <v>` -> dict."""
    tokens = state.split()
    if len(tokens) != 10:
        raise RuntimeError(f"unexpected Nose-Hoover state: {state!r}")
    out = {}
    for i in range(0, len(tokens), 2):
        key, value = tokens[i], tokens[i + 1]
        out[key] = int(value) if key == "dof" else float(value)
    return out


LOG_COLUMNS = [
    "step", "time", "pe", "ke", "etot", "temperature", "pressure", "volume",
    "shake_iter", "shake_error", "rattle_iter", "rattle_error",
]


def read_log(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        row = {}
        for name, raw in zip(LOG_COLUMNS, fields):
            row[name] = int(raw) if name in {"step", "shake_iter", "rattle_iter"} else float(raw)
        rows.append(row)
    return rows


# --- Force comparison ------------------------------------------------------

def forces_from_checkpoint(args, work: Path, name: str, state: dict) -> list[list[float]]:
    """Re-evaluate forces at a checkpoint's full-precision coordinates.

    The trajectory .xyz is written at six decimal places, which is far too
    coarse for the tolerances used here, so the geometry is rebuilt from the
    checkpoint and handed to `gmd_validate`, which reports per-atom forces.
    """
    directory = work / f"forces_{name}"
    directory.mkdir(parents=True, exist_ok=True)

    tags = sorted(state["atoms"])
    lines = [str(len(tags)), " ".join(str(v) for v in state["box"])]
    for position, tag in enumerate(tags):
        atom = state["atoms"][tag]
        atom_type = ATOMS[position][0]
        p = atom["pos"]
        v = atom["vel"]
        lines.append(
            f"{atom_type}  {p[0]!r} {p[1]!r} {p[2]!r}  {v[0]!r} {v[1]!r} {v[2]!r}"
        )
    (directory / "input.xyz").write_text("\n".join(lines) + "\n")
    write_ff(directory / "ff.ff")
    write_top(directory / "top.top")
    write_run(directory / "run.in", steps=1, checkpoint_file=directory / "unused.gmdchk")

    result = subprocess.run(
        [args.gmd_validate, "input.xyz", "run.in", "ff.ff", "top.top",
         "--json", "result.json"],
        cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"gmd_validate failed for {name}:\n{result.stdout}")
    return json.loads((directory / "result.json").read_text())["force"]["atoms"]


# --- Checks ----------------------------------------------------------------

def max_abs(lhs, rhs) -> float:
    return max(abs(a - b) for a, b in zip(lhs, rhs))


def check_continuity(args, work: Path, continuous_dir: Path, split2_dir: Path) -> dict:
    continuous = read_checkpoint(continuous_dir / "final.gmdchk")
    split = read_checkpoint(split2_dir / "final.gmdchk")

    problems: list[str] = []
    report: dict = {}

    # --- step and simulation time ---
    expected_time = N_TOTAL * TIME_STEP_FS
    for name, state in (("continuous", continuous), ("split", split)):
        if state["step"] != N_TOTAL:
            problems.append(f"{name} final step is {state['step']}, expected {N_TOTAL}")
        if not math.isclose(state["time_fs"], expected_time, abs_tol=1e-12):
            problems.append(f"{name} final time is {state['time_fs']}, expected {expected_time}")

    # --- coordinates and velocities ---
    if set(continuous["atoms"]) != set(split["atoms"]):
        problems.append("atom tag sets differ between runs")
    else:
        pos_error = vel_error = 0.0
        for tag in sorted(continuous["atoms"]):
            pos_error = max(pos_error,
                            max_abs(continuous["atoms"][tag]["pos"], split["atoms"][tag]["pos"]))
            vel_error = max(vel_error,
                            max_abs(continuous["atoms"][tag]["vel"], split["atoms"][tag]["vel"]))
        report["max_position_error"] = pos_error
        report["max_velocity_error"] = vel_error

        # The run must have gone somewhere, or comparing endpoints proves nothing.
        travelled = 0.0
        tags = sorted(continuous["atoms"])
        for position, tag in enumerate(tags):
            start = [ATOMS[position][3], ATOMS[position][4], ATOMS[position][5]]
            travelled = max(travelled, max_abs(continuous["atoms"][tag]["pos"], start))
        report["max_displacement_from_start"] = travelled
        if travelled < 1.0e-3:
            problems.append(f"atoms barely moved ({travelled} A); the comparison is vacuous")
        if pos_error > args.position_tolerance:
            problems.append(f"position error {pos_error} exceeds {args.position_tolerance}")
        if vel_error > args.velocity_tolerance:
            problems.append(f"velocity error {vel_error} exceeds {args.velocity_tolerance}")

    # --- Nose-Hoover extended-system state ---
    a = parse_thermostat_state(continuous["thermostat_state"])
    b = parse_thermostat_state(split["thermostat_state"])
    report["thermostat_continuous"] = a
    report["thermostat_split"] = b

    # Guard against a vacuous pass. If the thermostat never engaged, every
    # comparison below would trivially hold at zero.
    if a["Q"] <= 0.0:
        problems.append("Nose-Hoover mass Q is not positive; the thermostat never engaged")
    if a["xi"] == 0.0:
        problems.append("Nose-Hoover friction xi is exactly zero; the thermostat never engaged")
    if a["current_temperature"] <= 0.0:
        problems.append("Nose-Hoover cached temperature is not positive; nothing was integrated")

    if a["dof"] != EXPECTED_DOF:
        problems.append(f"continuous DOF is {a['dof']}, expected {EXPECTED_DOF} "
                        f"(3N - 3 COM - 3 constraints)")
    if b["dof"] != a["dof"]:
        problems.append(f"restart DOF {b['dof']} differs from continuous {a['dof']}")
    for key, tolerance in (("tau", 0.0), ("xi", args.thermostat_tolerance),
                           ("Q", args.thermostat_tolerance),
                           ("current_temperature", args.thermostat_tolerance)):
        delta = abs(a[key] - b[key])
        report[f"thermostat_{key}_error"] = delta
        if delta > tolerance:
            problems.append(f"Nose-Hoover {key} differs by {delta} (tolerance {tolerance}): "
                            f"continuous {a[key]}, split {b[key]}")

    # --- logged thermodynamics, over the whole post-restart overlap ---
    continuous_log = {row["step"]: row for row in read_log(continuous_dir / "output.log")}
    split_log = {row["step"]: row for row in read_log(split2_dir / "output.log")}

    overlap = sorted(set(continuous_log) & set(split_log))
    if not overlap:
        problems.append("continuous and restarted logs share no steps")
    report["compared_steps"] = len(overlap)
    if len(overlap) < (N_TOTAL - N_FIRST) // 2:
        problems.append(f"only {len(overlap)} shared log steps; expected the whole "
                        f"post-restart region")

    worst = {key: 0.0 for key in ("pe", "ke", "etot", "temperature", "pressure")}
    for step in overlap:
        for key in worst:
            worst[key] = max(worst[key], abs(continuous_log[step][key] - split_log[step][key]))
    report["log_errors"] = worst

    for key, tolerance in (("pe", args.energy_tolerance),
                           ("ke", args.energy_tolerance),
                           ("etot", args.energy_tolerance),
                           ("temperature", args.temperature_tolerance),
                           ("pressure", args.pressure_tolerance)):
        if worst[key] > tolerance:
            problems.append(f"log {key} differs by {worst[key]} over {len(overlap)} shared "
                            f"steps (tolerance {tolerance})")

    # The restarted run must actually pick up where it left off.
    if split_log and min(split_log) <= N_FIRST - 1:
        problems.append(f"restarted log starts at step {min(split_log)}, expected > {N_FIRST - 1}")

    # --- constraint residuals ---
    for name, directory in (("continuous", continuous_dir), ("split", split2_dir)):
        rows = read_log(directory / "output.log")
        residual = max(row["shake_error"] for row in rows)
        rattle_residual = max(row["rattle_error"] for row in rows)
        report[f"{name}_max_shake_error"] = residual
        report[f"{name}_max_rattle_error"] = rattle_residual
        if residual > args.constraint_tolerance:
            problems.append(f"{name} SHAKE residual {residual} exceeds {args.constraint_tolerance}")
        if rattle_residual > args.constraint_tolerance:
            problems.append(f"{name} RATTLE residual {rattle_residual} exceeds "
                            f"{args.constraint_tolerance}")

    # Constrained bond lengths must still hold at their target after restart.
    worst_bond = 0.0
    for state_name, state in (("continuous", continuous), ("split", split)):
        tags = sorted(state["atoms"])
        for first, second in ((0, 1), (1, 2), (2, 3)):
            p = state["atoms"][tags[first]]["pos"]
            q = state["atoms"][tags[second]]["pos"]
            length = math.dist(p, q)
            worst_bond = max(worst_bond, abs(length - BOND_R0))
    report["max_constrained_bond_deviation"] = worst_bond
    if worst_bond > args.constraint_tolerance:
        problems.append(f"constrained bond length deviates by {worst_bond} from {BOND_R0} "
                        f"(tolerance {args.constraint_tolerance})")

    # --- forces ---
    continuous_forces = forces_from_checkpoint(args, work, "continuous", continuous)
    split_forces = forces_from_checkpoint(args, work, "split", split)
    if len(continuous_forces) != len(split_forces):
        problems.append("force vector counts differ")
    else:
        force_error = 0.0
        for lhs, rhs in zip(continuous_forces, split_forces):
            force_error = max(force_error, max_abs(lhs, rhs))
        report["max_force_error"] = force_error
        if force_error > args.force_tolerance:
            problems.append(f"force error {force_error} exceeds {args.force_tolerance}")

    if problems:
        raise RuntimeError("constrained Nose-Hoover restart continuity failed:\n  - "
                           + "\n  - ".join(problems))
    return report


def check_rejected_restart(args, directory: Path, description: str,
                           checkpoint: Path, expect_in_error: list[str]) -> dict:
    """A restart that must be refused before any step is taken."""
    result = run_gmd(args, directory, check=False)
    combined = result.stdout + result.stderr

    problems = []
    if result.returncode == 0:
        problems.append(f"{description}: gmd exited 0, expected a rejection")
    for needle in expect_in_error:
        if needle not in combined:
            problems.append(f"{description}: error text does not mention {needle!r}\n{combined}")

    # Nothing may have been advanced: no checkpoint, and no data rows logged.
    if checkpoint.exists():
        problems.append(f"{description}: a checkpoint was written despite rejection")
    log_path = directory / "output.log"
    if log_path.exists():
        rows = read_log(log_path)
        if rows:
            problems.append(f"{description}: {len(rows)} log rows written despite rejection")
    xyz_path = directory / "output.xyz"
    if xyz_path.exists() and xyz_path.stat().st_size > 0:
        problems.append(f"{description}: trajectory frames written despite rejection")

    if problems:
        raise RuntimeError("\n  - ".join(["negative restart check failed:"] + problems))
    return {"description": description, "returncode": result.returncode}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--mpiexec")
    parser.add_argument("--mpiexec-np-flag", default="-n")
    parser.add_argument("--position-tolerance", type=float, default=1e-10)
    parser.add_argument("--velocity-tolerance", type=float, default=1e-10)
    parser.add_argument("--force-tolerance", type=float, default=1e-9)
    parser.add_argument("--energy-tolerance", type=float, default=1e-6)
    parser.add_argument("--temperature-tolerance", type=float, default=1e-6)
    parser.add_argument("--pressure-tolerance", type=float, default=1e-4)
    parser.add_argument("--thermostat-tolerance", type=float, default=1e-10)
    parser.add_argument("--constraint-tolerance", type=float, default=1e-8)
    args = parser.parse_args()

    if args.np > 1 and not args.mpiexec:
        parser.error("--mpiexec is required when --np > 1")

    work = Path(args.work_root)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    # --- continuous reference run ---
    continuous_dir = work / "continuous"
    prepare(continuous_dir, steps=N_TOTAL,
            checkpoint_file=continuous_dir / "final.gmdchk")
    run_gmd(args, continuous_dir)

    # --- split run: first leg, then restart ---
    split1_dir = work / "split1"
    prepare(split1_dir, steps=N_FIRST,
            checkpoint_file=split1_dir / "restart.gmdchk")
    run_gmd(args, split1_dir)

    restart_checkpoint = split1_dir / "restart.gmdchk"
    if not restart_checkpoint.exists():
        raise RuntimeError("first leg did not write a checkpoint")

    split2_dir = work / "split2"
    prepare(split2_dir, steps=N_TOTAL - N_FIRST,
            checkpoint_file=split2_dir / "final.gmdchk",
            restart_from=restart_checkpoint)
    run_gmd(args, split2_dir)

    report = check_continuity(args, work, continuous_dir, split2_dir)

    # --- negative case 1: COM-removal setting changed across the restart ---
    com_dir = work / "reject_com"
    prepare(com_dir, steps=N_TOTAL - N_FIRST,
            checkpoint_file=com_dir / "final.gmdchk",
            restart_from=restart_checkpoint,
            remove_com=False)
    report["reject_com"] = check_rejected_restart(
        args, com_dir, "changed COM-removal setting",
        com_dir / "final.gmdchk",
        # DOF goes 6 -> 9 when the three COM modes are kept.
        ["degrees of freedom", "centre-of-mass", str(EXPECTED_DOF), "9"],
    )

    # --- negative case 2: constraint configuration changed ---
    constraint_dir = work / "reject_constraints"
    prepare(constraint_dir, steps=N_TOTAL - N_FIRST,
            checkpoint_file=constraint_dir / "final.gmdchk",
            restart_from=restart_checkpoint,
            constrain=False)
    report["reject_constraints"] = check_rejected_restart(
        args, constraint_dir, "changed constraint configuration",
        constraint_dir / "final.gmdchk",
        # DOF goes 6 -> 9 when the three bond constraints are dropped.
        ["degrees of freedom", "constraint", str(EXPECTED_DOF), "9"],
    )

    (work / "summary.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
