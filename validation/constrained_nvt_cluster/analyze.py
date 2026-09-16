#!/usr/bin/env python3
"""Constrained NVT validation: four interacting rigid molecules, real CLI.

WHAT THIS CASE ADDS OVER constrained_nve_water
----------------------------------------------
The NVE case is a single force-free molecule, so its provider virial is exactly
zero and the constraint virial is the whole of it. Here the four molecules
interact through intermolecular Lennard-Jones terms, so the provider virial is
non-zero and the constraint term has a genuine partner: the pressure identity
below tests the SUM of two independently produced tensors at the same time level.

The constraint graph has four components and twelve constraints, two of the
molecules straddle an MPI rank boundary, and at np=4 one rank owns no atoms at
all.

WHAT IS NOT CLAIMED
-------------------
This is a 400-step deterministic run. It is DYNAMICS AND REGRESSION COVERAGE --
that the thermostat holds the temperature near its target with bounded
excursions, that the constraints stay satisfied, that the virial and pressure
stay consistent, and that a restart reproduces the trajectory. It is NOT a
statistical validation of the canonical ensemble: 400 steps of one deterministic
trajectory cannot establish an ensemble distribution, and nothing here claims it
does.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import shutil
import statistics
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import load_json, write_json, write_series_csv                  # noqa: E402
from constrained_common import (                                            # noqa: E402
    ValidationFailure, check_tensor_finite, compare_tensors,
    constraint_residuals, describe_constraint, max_asymmetry, parse_checkpoint,
    parse_constrained_log, trace, TENSOR_COMPONENT_NAMES,
)

INPUTS = ["cluster.xyz", "cluster.ff", "cluster.top", "nvt.run"]

R_OH = 0.95720000000000
D_HH = 1.51390065452732
MOLECULE_COUNT = 4
ATOM_COUNT = 12
# Global atom tags are 0-based: molecule m owns tags 3m, 3m+1, 3m+2 (O, H, H).
CONSTRAINTS: list[tuple[int, int, float]] = []
for _m in range(MOLECULE_COUNT):
    _o, _h1, _h2 = 3 * _m, 3 * _m + 1, 3 * _m + 2
    CONSTRAINTS += [(_o, _h1, R_OH), (_o, _h2, R_OH), (_h1, _h2, D_HH)]

EXPECTED_DOF = 21          # 3N - constraints - 3 = 36 - 12 - 3
UNCONSTRAINED_DOF = 33     # 3N - 3, the value a DOF regression would produce
TARGET_TEMPERATURE = 300.0
TOTAL_STEPS = 400
RESIDUAL_SAMPLE_STEPS = [100, 200, 300, 400]


def run_gmd(gmd: str, work: pathlib.Path, run_file: str, np: int = 1,
            mpiexec: str = "mpiexec") -> str:
    command = []
    if np > 1:
        command += [mpiexec, "-np", str(np)]
    command += [gmd, "cluster.xyz", run_file, "cluster.ff", "cluster.top"]
    completed = subprocess.run(command, cwd=work, text=True, capture_output=True)
    if completed.returncode != 0:
        raise ValidationFailure(
            f"gmd failed (exit {completed.returncode}, np={np}) on {run_file} "
            f"in {work}\n--- stdout ---\n{completed.stdout}\n"
            f"--- stderr ---\n{completed.stderr}")
    return completed.stdout


def write_run(work: pathlib.Path, name: str, *, steps: int, checkpoint: str,
              restart_from: str | None = None) -> None:
    lines = []
    for line in (work / "nvt.run").read_text(encoding="utf-8").splitlines():
        head = line.strip().split()
        key = head[0] if head else ""
        if key in ("run", "write_checkpoint_every", "checkpoint_file", "restart_from"):
            continue
        lines.append(line)
    lines += [f"run {steps}",
              f"write_checkpoint_every {steps}",
              f"checkpoint_file {checkpoint}"]
    if restart_from is not None:
        lines.append(f"restart_from {restart_from}")
    (work / name).write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_dof(stdout: str) -> int:
    for line in stdout.splitlines():
        if "Degrees of freedom:" in line:
            return int(line.split("Degrees of freedom:")[1].split()[0])
    raise ValidationFailure(
        "gmd did not report a degree-of-freedom count on stdout:\n" + stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--mpiexec", default="")
    parser.add_argument("--mpi-ranks", default="")
    args = parser.parse_args()

    case_dir = pathlib.Path(__file__).resolve().parent
    work = pathlib.Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    for name in INPUTS:
        shutil.copy2(case_dir / name, work / name)

    reference = load_json(case_dir / "reference.json")
    tolerance = load_json(case_dir / "tolerance.json")
    checks: dict = {}

    def record(name: str, passed: bool, **detail) -> None:
        checks[name] = {"passed": bool(passed), **detail}

    # ---------------------------------------------------------------- baseline
    write_run(work, "base.run", steps=TOTAL_STEPS, checkpoint="base.gmdchk")
    stdout = run_gmd(args.gmd, work, "base.run")

    # ---- authoritative constrained DOF -----------------------------------
    dof = parse_dof(stdout)
    record("authoritative_constrained_dof", dof == EXPECTED_DOF,
           actual=dof, expected=EXPECTED_DOF,
           unconstrained_would_be=UNCONSTRAINED_DOF,
           note="3N - constraints - 3 = 36 - 12 - 3 = 21. A run whose DOF "
                "reverted to the unconstrained 3N-3 = 33 would report a "
                "temperature low by 33/21 and thermostat to the wrong energy, "
                "so this value is load bearing rather than cosmetic")

    rows = parse_constrained_log(work / "output.log")
    write_series_csv(work / "thermo_series.csv", rows, list(rows[0].keys()))

    # ---- temperature after constraint-aware initialization ---------------
    # The initializer projects into the constraint tangent space and only then
    # rescales, using the constrained DOF. Both must be right for frame 0 to sit
    # on the target: projecting after rescaling removes kinetic energy and lands
    # low, and using 3N-3 lands high by 33/21.
    initial_temperature = rows[0]["temperature"]
    record("initial_temperature_after_constrained_initialization",
           abs(initial_temperature - TARGET_TEMPERATURE)
           <= tolerance["initial_temperature_abs_k"],
           actual=initial_temperature, target=TARGET_TEMPERATURE,
           tolerance=tolerance["initial_temperature_abs_k"],
           dof_used=dof)

    # ---- temperature mean and bounded fluctuation ------------------------
    temperatures = [row["temperature"] for row in rows]
    mean_temperature = statistics.fmean(temperatures)
    record("temperature_mean",
           abs(mean_temperature - TARGET_TEMPERATURE)
           <= tolerance["temperature_mean_abs_k"],
           mean=mean_temperature, target=TARGET_TEMPERATURE,
           tolerance=tolerance["temperature_mean_abs_k"],
           note="dynamics coverage over a short deterministic run, NOT a "
                "statistical test of the canonical ensemble")
    excursion = max(abs(value - TARGET_TEMPERATURE) for value in temperatures)
    record("temperature_bounded",
           excursion <= tolerance["temperature_max_excursion_k"],
           max_excursion_k=excursion, minimum=min(temperatures),
           maximum=max(temperatures),
           tolerance=tolerance["temperature_max_excursion_k"])

    # ---- per-frame log bounds --------------------------------------------
    log_bound = tolerance["log_residual_column_bound_a"]
    frame_failures = []
    for row in rows:
        if row["shake_error_a"] > log_bound:
            frame_failures.append(
                f"step {row['step']}: log column shake_error[A] = "
                f"{row['shake_error_a']} exceeds {log_bound}")
        if row["rattle_error_a_per_fs"] > log_bound:
            frame_failures.append(
                f"step {row['step']}: log column rattle_error = "
                f"{row['rattle_error_a_per_fs']} exceeds {log_bound}")
        if not (row["temperature"] > 0.0):
            frame_failures.append(
                f"step {row['step']}: T[K] = {row['temperature']} is not positive")
    record("per_frame_log_bounds", not frame_failures, frames=len(rows),
           failures=frame_failures[:8])

    record("initial_frame_pressure_invalid",
           rows[0]["pressure_valid"] == 0 and rows[0]["pressure_bar"] is None,
           note="no RATTLE multiplier belongs to the initial provider virial")
    later = [row for row in rows if row["step"] > 0]
    record("completed_steps_pressure_valid",
           all(row["pressure_valid"] == 1 for row in later),
           invalid=[row["step"] for row in later if row["pressure_valid"] != 1])

    # ---- tight residuals sampled along the trajectory --------------------
    residual_series = []
    worst_position = {"value": 0.0, "step": None, "constraint": None}
    worst_velocity = {"value": 0.0, "step": None, "constraint": None}
    for steps in RESIDUAL_SAMPLE_STEPS:
        write_run(work, f"sample{steps}.run", steps=steps,
                  checkpoint=f"sample{steps}.gmdchk")
        run_gmd(args.gmd, work, f"sample{steps}.run")
        state = parse_checkpoint(work / f"sample{steps}.gmdchk")
        residuals = constraint_residuals(state, CONSTRAINTS)
        position = residuals["max_position_residual_a"]
        velocity = residuals["max_velocity_residual_internal"]
        residual_series.append({
            "step": steps,
            "max_position_residual_a": position["value"],
            "max_velocity_residual_internal": velocity["value"],
        })
        if position["value"] > worst_position["value"]:
            worst_position = {"value": position["value"], "step": steps,
                              "constraint": position["constraint"]}
        if velocity["value"] > worst_velocity["value"]:
            worst_velocity = {"value": velocity["value"], "step": steps,
                              "constraint": velocity["constraint"]}
    write_series_csv(work / "residual_series.csv", residual_series,
                     ["step", "max_position_residual_a",
                      "max_velocity_residual_internal"])

    record("max_position_residual",
           worst_position["value"] <= tolerance["position_residual_abs_a"],
           value=worst_position["value"],
           tolerance=tolerance["position_residual_abs_a"],
           step=worst_position["step"],
           where=describe_constraint(worst_position["constraint"]))
    record("max_velocity_tangency_residual",
           worst_velocity["value"] <= tolerance["velocity_residual_abs_internal"],
           value=worst_velocity["value"],
           tolerance=tolerance["velocity_residual_abs_internal"],
           step=worst_velocity["step"],
           where=describe_constraint(worst_velocity["constraint"]))

    # ---- thermostat state -------------------------------------------------
    final = parse_checkpoint(work / "base.gmdchk")
    thermostat = final["thermostat"]
    friction_values = []
    for token in thermostat["state"].split():
        try:
            friction_values.append(float(token))
        except ValueError:
            continue
    record("thermostat_state_finite",
           thermostat["type"] == "nose_hoover"
           and bool(friction_values)
           and all(math.isfinite(value) for value in friction_values),
           type=thermostat["type"], state=thermostat["state"],
           note="the Nose-Hoover friction variable is the extended-system "
                "quantity this build exposes; it must be finite for the "
                "extended dynamics to mean anything")

    # ---- constraint virial and pressure ----------------------------------
    constraint_virial = final["constraint_virial"]
    record("constraint_virial_state_valid",
           constraint_virial["state"] == "valid"
           and constraint_virial["time_level"] == "endpoint_rattle_t_plus_dt",
           state=constraint_virial["state"],
           time_level=constraint_virial["time_level"])

    tensor = constraint_virial["tensor"]
    check_tensor_finite(tensor, "constraint virial")
    provider = final["provider_virial"]
    check_tensor_finite(provider["tensor"], "provider virial")
    record("provider_virial_valid_and_nonzero",
           provider["valid"]
           and max(abs(v) for v in provider["tensor"])
           >= tolerance["provider_virial_min_abs"],
           valid=provider["valid"],
           largest_component=max(abs(v) for v in provider["tensor"]),
           threshold=tolerance["provider_virial_min_abs"],
           note="unlike the NVE water case the molecules interact here, so the "
                "constraint virial is added to a real provider tensor")

    record("constraint_virial_symmetric",
           max_asymmetry(tensor) <= tolerance["virial_symmetry_abs"],
           max_asymmetry=max_asymmetry(tensor),
           tolerance=tolerance["virial_symmetry_abs"])
    smallest_offdiagonal = min(abs(tensor[i]) for i in (1, 2, 5))
    record("constraint_virial_offdiagonals_nontrivial",
           smallest_offdiagonal >= tolerance["virial_offdiagonal_min_abs"],
           smallest_offdiagonal=smallest_offdiagonal,
           threshold=tolerance["virial_offdiagonal_min_abs"],
           components={name: value
                       for name, value in zip(TENSOR_COMPONENT_NAMES, tensor)})

    step_pressure = final["step_pressure"]
    if not step_pressure["valid"]:
        raise ValidationFailure(
            "the completed-step pressure must be valid at the end of a "
            "constrained NVT run, but the checkpoint reports it invalid")
    # The stored completed-step virial must be the SUM of the two tensors, at the
    # same time level. Checking the sum rather than the total alone is what makes
    # a one-step shift of the constraint term detectable.
    expected_total = [provider["tensor"][i] + tensor[i] for i in range(9)]
    total_match = compare_tensors(step_pressure["tensor"], expected_total,
                                  tolerance["virial_sum_abs"],
                                  "completed-step virial vs provider + constraint")
    record("total_virial_is_provider_plus_constraint", total_match["passed"],
           **{k: v for k, v in total_match.items() if k != "passed"},
           note="both tensors belong to t+dt; this is the check that a "
                "constraint term taken from the wrong step would fail")

    identity = (step_pressure["twice_kinetic_energy"] + trace(step_pressure["tensor"])) \
        / (3.0 * step_pressure["volume"])
    record("pressure_kinetic_plus_virial_identity",
           abs(identity - step_pressure["pressure"])
           <= tolerance["pressure_identity_rel"] * abs(step_pressure["pressure"]),
           reported=step_pressure["pressure"], recomputed=identity,
           relative_difference=abs(identity - step_pressure["pressure"])
           / abs(step_pressure["pressure"]),
           tolerance=tolerance["pressure_identity_rel"])

    # ---- checkpoint / restart continuity ---------------------------------
    half = TOTAL_STEPS // 2
    write_run(work, "restart_first.run", steps=half, checkpoint="restart_half.gmdchk")
    run_gmd(args.gmd, work, "restart_first.run")
    shutil.copy2(work / "output.log", work / "restart_first.log")
    write_run(work, "restart_second.run", steps=half,
              checkpoint="restart_final.gmdchk", restart_from="restart_half.gmdchk")
    run_gmd(args.gmd, work, "restart_second.run")
    restarted = parse_checkpoint(work / "restart_final.gmdchk")
    restarted_rows = parse_constrained_log(work / "output.log")

    position_difference, velocity_difference, worst_tag = 0.0, 0.0, None
    by_tag = {a["tag"]: a for a in restarted["atoms"]}
    for atom in final["atoms"]:
        other = by_tag[atom["tag"]]
        for d in range(3):
            dp = abs(atom["position"][d] - other["position"][d])
            if dp > position_difference:
                position_difference, worst_tag = dp, atom["tag"]
            velocity_difference = max(
                velocity_difference, abs(atom["velocity"][d] - other["velocity"][d]))
    restart_virial = compare_tensors(
        restarted["constraint_virial"]["tensor"], tensor,
        tolerance["restart_virial_abs"], "restart vs continuous constraint virial")
    record("restart_continuity",
           position_difference <= tolerance["restart_position_abs_a"]
           and velocity_difference <= tolerance["restart_velocity_abs_internal"]
           and restart_virial["passed"]
           and restarted["constraint_virial"]["state"] == "valid",
           max_position_difference_a=position_difference,
           max_velocity_difference_internal=velocity_difference,
           worst_atom_tag=worst_tag, constraint_virial=restart_virial,
           restored_state=restarted["constraint_virial"]["state"])

    # THE RESTORED COMPLETED-STEP PRESSURE. A fresh constrained run reports its
    # initial frame as having NO valid pressure, because no step has completed
    # and no RATTLE multiplier belongs to the initial provider virial. A
    # RESTARTED run is different: the checkpoint carries the completed-step
    # pressure of the step it was written after, so the first frame must report
    # that value and must report it as valid. This is what catches a restart
    # that drops the recorded pressure state and silently falls back to the
    # current geometry.
    half_state = parse_checkpoint(work / "restart_half.gmdchk")
    first_restarted = restarted_rows[0]
    expected_bar = half_state["step_pressure"]["pressure"] * 1.602176634e6
    record("restart_restores_completed_step_pressure",
           first_restarted["pressure_valid"] == 1
           and first_restarted["pressure_bar"] is not None
           and abs(first_restarted["pressure_bar"] - expected_bar)
           <= tolerance["restart_pressure_abs_bar"],
           step=first_restarted["step"],
           reported_bar=first_restarted["pressure_bar"],
           checkpoint_bar=expected_bar,
           tolerance=tolerance["restart_pressure_abs_bar"],
           note="a restarted run's first frame carries the checkpoint's "
                "completed-step pressure, unlike a fresh run's invalid one")

    # NOTE ON WHAT THIS CANNOT SEE. The checkpoint also stores the constraint
    # virial itself, with its time level. Dropping THAT alone is not observable
    # through any CLI output: the reported pressure of the restored frame comes
    # from the completed-step record checked just above, and the next step's
    # RATTLE recomputes the constraint virial from scratch. The stored tensor is
    # a reporting quantity for the restored frame, and the CLI does not expose a
    # current-geometry virial. tests/mpi_constraint_virial.cpp covers the tensor
    # itself at the library level.

    # No systematic jump: the first frames after the restart must continue the
    # continuous run's temperature and energy rather than step away from them.
    continuous_by_step = {row["step"]: row for row in rows}
    jump_failures = []
    for row in restarted_rows:
        absolute_step = row["step"]
        if absolute_step not in continuous_by_step:
            continue
        reference_row = continuous_by_step[absolute_step]
        temperature_jump = abs(row["temperature"] - reference_row["temperature"])
        energy_jump = abs(row["total_energy"] - reference_row["total_energy"])
        if temperature_jump > tolerance["restart_temperature_abs_k"]:
            jump_failures.append(
                f"step {absolute_step}: temperature {row['temperature']} after "
                f"restart against {reference_row['temperature']} continuous "
                f"(jump {temperature_jump:.3e} K)")
        if energy_jump > tolerance["restart_energy_abs_ev"]:
            jump_failures.append(
                f"step {absolute_step}: total energy {row['total_energy']} after "
                f"restart against {reference_row['total_energy']} continuous "
                f"(jump {energy_jump:.3e} eV)")
    record("no_restart_energy_or_temperature_jump", not jump_failures,
           compared_frames=len([r for r in restarted_rows
                                if r["step"] in continuous_by_step]),
           failures=jump_failures[:8])

    # ---- MPI -------------------------------------------------------------
    ranks = [int(r) for r in args.mpi_ranks.split(",") if r.strip()] if args.mpi_ranks else []
    if args.mpiexec and ranks:
        differences = []
        for np in ranks:
            directory = work / f"np{np}"
            directory.mkdir(exist_ok=True)
            for name in INPUTS:
                shutil.copy2(case_dir / name, directory / name)
            write_run(directory, "base.run", steps=TOTAL_STEPS,
                      checkpoint="base.gmdchk")
            rank_stdout = run_gmd(args.gmd, directory, "base.run", np=np,
                                  mpiexec=args.mpiexec)
            if parse_dof(rank_stdout) != EXPECTED_DOF:
                differences.append(
                    f"np={np}: reported DOF {parse_dof(rank_stdout)} against "
                    f"{EXPECTED_DOF} serial")
            state = parse_checkpoint(directory / "base.gmdchk")
            comparison = compare_tensors(
                state["constraint_virial"]["tensor"], tensor,
                tolerance["mpi_virial_abs"], f"np={np} constraint virial")
            if not comparison["passed"]:
                differences.append(
                    f"np={np}: constraint virial component "
                    f"{comparison['worst_component']} differs from serial by "
                    f"{comparison['max_component_difference']:.3e} "
                    f"(tolerance {tolerance['mpi_virial_abs']:.3e}); a tensor "
                    f"multiplied by the rank count would differ by ~{np - 1}x "
                    f"its own magnitude")
            others = {a["tag"]: a for a in state["atoms"]}
            worst = 0.0
            worst_tag_mpi = None
            for atom in final["atoms"]:
                for d in range(3):
                    difference = abs(atom["position"][d]
                                     - others[atom["tag"]]["position"][d])
                    if difference > worst:
                        worst, worst_tag_mpi = difference, atom["tag"]
            if worst > tolerance["mpi_position_abs_a"]:
                differences.append(
                    f"np={np}: atom tag {worst_tag_mpi} position differs from "
                    f"serial by {worst:.3e} A")
            residuals = constraint_residuals(state, CONSTRAINTS)
            if residuals["max_position_residual_a"]["value"] > \
                    tolerance["position_residual_abs_a"]:
                differences.append(
                    f"np={np}: "
                    + describe_constraint(
                        residuals["max_position_residual_a"]["constraint"])
                    + f" residual {residuals['max_position_residual_a']['value']:.3e}")
        record("mpi_agreement", not differences, ranks=ranks, differences=differences,
               note="molecule C straddles x=Lx/2 and molecule D straddles "
                    "y=Ly/2, so constraints cross rank boundaries at np=2 and "
                    "np=4; at np=4 one rank owns no atoms")
    else:
        checks["mpi_agreement"] = {"passed": True, "skipped": True,
                                   "note": "no mpiexec supplied"}

    summary = {
        "case": "constrained_nvt_cluster",
        "ensemble": "NVT (Nose-Hoover)",
        "claim": "dynamics and regression coverage, not statistical ensemble "
                 "validation",
        "baseline_kind": reference["provenance"]["baseline_kind"],
        "checks": checks,
        "passed": all(check["passed"] for check in checks.values()),
    }
    write_json(work / "summary.json", summary)

    failed = [name for name, check in checks.items() if not check["passed"]]
    if failed:
        sys.stderr.write("FAILED constrained_nvt_cluster checks: "
                         + ", ".join(failed) + "\n")
        for name in failed:
            sys.stderr.write(f"  {name}: "
                             + json.dumps(checks[name], indent=2, sort_keys=True) + "\n")
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationFailure as failure:
        sys.stderr.write(f"constrained_nvt_cluster: {failure}\n")
        raise SystemExit(1)
