#!/usr/bin/env python3
"""Run a real CLI checkpoint/restart continuity test for GMD.

The test drives the `gmd` executable, not the checkpoint unit-test helpers:

  continuous: initial checkpoint -> 100 MD steps
  split:      initial checkpoint -> 50 MD steps -> checkpoint -> 50 MD steps

Final checkpoints and energy logs are compared by global atom tag.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path


ATOMS = [
    # tag, type, molecule, mass, charge, x, y, z, vx, vy, vz
    (0, 0, 101, 39.948, 0.0, 3.0, 3.0, 3.0, 0.003, 0.001, 0.000),
    (1, 0, 101, 39.948, 0.0, 7.2, 3.1, 3.0, -0.002, 0.000, 0.001),
    (2, 0, 202, 39.948, 0.0, 12.8, 12.0, 3.2, 0.000, -0.002, 0.001),
    (3, 0, 202, 39.948, 0.0, 17.0, 12.1, 3.1, -0.001, 0.002, -0.001),
]


def write_initial_checkpoint(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "GMD_CHECKPOINT 1",
                "step 0",
                "time_fs 0",
                "boundary periodic periodic periodic",
                "box 20 20 20",
                "xyz_file ",
                "run_file ",
                "force_field_file ",
                "topology_file ",
                "velocity_seed 5489",
                "force_field_summary inline_lj restart_continuity",
                "config_summary restart_continuity_initial",
                "thermostat_type ",
                "thermostat_state stateless",
                "barostat_type ",
                "barostat_state stateless",
                f"atom_count {len(ATOMS)}",
                "atoms tag type molecule mass charge x y z vx vy vz",
                *(" ".join(str(value) for value in atom) for atom in ATOMS),
                "topology_counts 0 0 0 0 0",
                "end",
                "",
            ]
        )
    )


def write_dummy_xyz(path: Path) -> None:
    lines = [
        str(len(ATOMS)),
        "20.0 20.0 20.0",
    ]
    for _, atom_type, _, _, _, x, y, z, vx, vy, vz in ATOMS:
        # XYZ type ids are 1-based; checkpoint stores internal 0-based types.
        lines.append(f"{atom_type + 1} {x} {y} {z} {vx} {vy} {vz}")
    path.write_text("\n".join(lines) + "\n")


def write_run(
    path: Path,
    *,
    steps: int,
    restart_from: Path,
    checkpoint_file: Path,
    np: int,
) -> None:
    lines = [
        "velocity 0.0",
        "time_step 1.0",
        f"run {steps}",
        "output_interval 100",
        "velocity_init input",
        "remove_com_velocity false",
        "",
        "force_field lj",
        "cutoff 8.5",
        "type 1 Ar epsilon 0.01032 sigma 3.405 charge 0.0",
        "",
        f"restart_from {restart_from}",
        f"checkpoint_file {checkpoint_file}",
        "write_checkpoint_every 0",
    ]
    if np == 2:
        lines.append("mpi_grid 2 1 1")
    elif np == 4:
        lines.append("mpi_grid 4 1 1")
    path.write_text("\n".join(lines) + "\n")


def run_gmd(args: argparse.Namespace, xyz: Path, run_file: Path, cwd: Path) -> None:
    command = []
    if args.np > 1:
        if not args.mpiexec:
            raise RuntimeError("--mpiexec is required for MPI restart continuity")
        command.extend([args.mpiexec, args.mpiexec_np_flag, str(args.np)])
    command.extend([args.gmd, str(xyz), str(run_file)])
    if args.np > 1:
        command.extend(["--np", str(args.np)])

    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (cwd / "command.out").write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {result.returncode}: {' '.join(command)}\n"
            f"{result.stdout}"
        )


def read_checkpoint(path: Path) -> dict:
    lines = path.read_text().splitlines()
    if not lines or lines[0] != "GMD_CHECKPOINT 1":
        raise RuntimeError(f"{path} is not a version-1 GMD checkpoint")

    data: dict = {"atoms": {}}
    index = 1
    while index < len(lines):
        line = lines[index]
        if line == "atoms tag type molecule mass charge x y z vx vy vz":
            count = int(data["atom_count"])
            for atom_line in lines[index + 1 : index + 1 + count]:
                fields = atom_line.split()
                tag = int(fields[0])
                data["atoms"][tag] = {
                    "type": int(fields[1]),
                    "molecule": int(fields[2]),
                    "mass": float(fields[3]),
                    "charge": float(fields[4]),
                    "pos": [float(fields[5]), float(fields[6]), float(fields[7])],
                    "vel": [float(fields[8]), float(fields[9]), float(fields[10])],
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
        index += 1
    return data


def read_final_energy(log_path: Path) -> tuple[int, float, float]:
    final = None
    for line in log_path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        final = (int(fields[0]), float(fields[1]), float(fields[4]))
    if final is None:
        raise RuntimeError(f"{log_path} did not contain an energy row")
    return final


def max_abs_delta(lhs: list[float], rhs: list[float]) -> float:
    return max(abs(a - b) for a, b in zip(lhs, rhs))


def compare(args: argparse.Namespace, continuous_dir: Path, split2_dir: Path) -> dict:
    continuous = read_checkpoint(continuous_dir / "final.gmdchk")
    split = read_checkpoint(split2_dir / "final.gmdchk")
    continuous_energy = read_final_energy(continuous_dir / "output.log")
    split_energy = read_final_energy(split2_dir / "output.log")

    errors = {
        "max_position_error": 0.0,
        "max_velocity_error": 0.0,
        "energy_error": abs(continuous_energy[2] - split_energy[2]),
    }

    expected_step = 100
    expected_time = 100.0
    for name, state, energy in (
        ("continuous", continuous, continuous_energy),
        ("split", split, split_energy),
    ):
        if state["step"] != expected_step or energy[0] != expected_step:
            raise RuntimeError(f"{name} final step is not {expected_step}")
        if not math.isclose(state["time_fs"], expected_time, abs_tol=1.0e-12):
            raise RuntimeError(f"{name} checkpoint time is not {expected_time} fs")
        if not math.isclose(energy[1], expected_time, abs_tol=1.0e-12):
            raise RuntimeError(f"{name} log time is not {expected_time} fs")
        if state["box"] != [20.0, 20.0, 20.0]:
            raise RuntimeError(f"{name} box changed unexpectedly: {state['box']}")

    continuous_tags = set(continuous["atoms"].keys())
    split_tags = set(split["atoms"].keys())
    if continuous_tags != split_tags:
        raise RuntimeError(f"atom tag mismatch: {continuous_tags} vs {split_tags}")
    if continuous_tags != {atom[0] for atom in ATOMS}:
        raise RuntimeError(f"final atom tags do not match initial tags: {continuous_tags}")

    expected_molecules = {atom[0]: atom[2] for atom in ATOMS}
    for tag in sorted(continuous_tags):
        lhs = continuous["atoms"][tag]
        rhs = split["atoms"][tag]
        if lhs["molecule"] != expected_molecules[tag]:
            raise RuntimeError(f"continuous molecule id changed for tag {tag}")
        if rhs["molecule"] != expected_molecules[tag]:
            raise RuntimeError(f"split molecule id changed for tag {tag}")
        if lhs["type"] != rhs["type"] or lhs["molecule"] != rhs["molecule"]:
            raise RuntimeError(f"atom metadata mismatch for tag {tag}")
        errors["max_position_error"] = max(
            errors["max_position_error"], max_abs_delta(lhs["pos"], rhs["pos"])
        )
        errors["max_velocity_error"] = max(
            errors["max_velocity_error"], max_abs_delta(lhs["vel"], rhs["vel"])
        )

    if errors["max_position_error"] > args.position_tolerance:
        raise RuntimeError(
            f"position error {errors['max_position_error']} exceeds "
            f"{args.position_tolerance}"
        )
    if errors["max_velocity_error"] > args.velocity_tolerance:
        raise RuntimeError(
            f"velocity error {errors['max_velocity_error']} exceeds "
            f"{args.velocity_tolerance}"
        )
    if errors["energy_error"] > args.energy_tolerance:
        raise RuntimeError(
            f"energy error {errors['energy_error']} exceeds {args.energy_tolerance}"
        )
    return {
        "np": args.np,
        "position_tolerance": args.position_tolerance,
        "velocity_tolerance": args.velocity_tolerance,
        "energy_tolerance": args.energy_tolerance,
        **errors,
        "continuous_energy": continuous_energy[2],
        "split_energy": split_energy[2],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--mpiexec")
    parser.add_argument("--mpiexec-np-flag", default="-n")
    parser.add_argument("--position-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--velocity-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--energy-tolerance", type=float, default=1.0e-6)
    args = parser.parse_args()

    work_root = Path(args.work_root)
    if work_root.exists():
        shutil.rmtree(work_root)
    work_root.mkdir(parents=True)

    initial_checkpoint = work_root / "initial.gmdchk"
    write_initial_checkpoint(initial_checkpoint)

    continuous_dir = work_root / "continuous"
    split1_dir = work_root / "split1"
    split2_dir = work_root / "split2"
    for directory in (continuous_dir, split1_dir, split2_dir):
        directory.mkdir()
        write_dummy_xyz(directory / "input.xyz")

    write_run(
        continuous_dir / "run.in",
        steps=100,
        restart_from=initial_checkpoint,
        checkpoint_file=continuous_dir / "final.gmdchk",
        np=args.np,
    )
    write_run(
        split1_dir / "run.in",
        steps=50,
        restart_from=initial_checkpoint,
        checkpoint_file=split1_dir / "restart_50.gmdchk",
        np=args.np,
    )
    write_run(
        split2_dir / "run.in",
        steps=50,
        restart_from=split1_dir / "restart_50.gmdchk",
        checkpoint_file=split2_dir / "final.gmdchk",
        np=args.np,
    )

    try:
        run_gmd(args, continuous_dir / "input.xyz", continuous_dir / "run.in", continuous_dir)
        run_gmd(args, split1_dir / "input.xyz", split1_dir / "run.in", split1_dir)
        run_gmd(args, split2_dir / "input.xyz", split2_dir / "run.in", split2_dir)
        metrics = compare(args, continuous_dir, split2_dir)
    except Exception as exc:
        print(f"restart continuity failed: {exc}", file=sys.stderr)
        return 1

    metrics_path = work_root / "restart_continuity_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metrics, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
