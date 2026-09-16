#!/usr/bin/env python3
"""Constrained NVE validation: a free rigid water molecule through the real CLI.

Everything here is driven through the production `gmd` executable. Nothing calls
the constraint solver directly, and no expected value is taken from the solver's
own convergence flag: the residuals are RECOMPUTED from checkpointed state at 17
significant digits (see validation/constrained_common.py).

WHAT MAKES THIS CASE WORTH HAVING
---------------------------------
The molecule is alone in the cell and every intramolecular pair is excluded with
zero force constants, so it feels no force whatsoever. It is a FREE RIGID BODY,
and a free rigid body has analytically known behaviour that does not depend on
any previous GMD run:

  * the total energy is exactly constant;
  * the linear momentum is exactly constant, and the centre of mass travels in a
    straight line at constant velocity;
  * the angular momentum about the centre of mass is exactly constant;
  * the constraint virial cancels the ROTATIONAL kinetic energy, so the reported
    pressure collapses to the ideal-gas pressure of the centre of mass,
    M |v_com|^2 / 3V.

That last one is the sharpest check in this file. It pins the constraint virial's
magnitude (the 2/dt endpoint conversion), its sign and its time level all at
once, against a number derived from mechanics rather than from GMD.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import shutil
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from common import load_json, write_json, write_series_csv                  # noqa: E402
from constrained_common import (                                            # noqa: E402
    INTERNAL_TIME_PER_FEMTOSECOND, ValidationFailure, angular_momentum,
    minimum_image,
    center_of_mass_velocity, check_tensor_finite, compare_tensors,
    constraint_residuals, describe_constraint, group_center_of_mass,
    max_asymmetry, parse_checkpoint, parse_constrained_log,
    separation_modulo_box, total_momentum, trace, twice_kinetic_energy,
    TENSOR_COMPONENT_NAMES,
)

INPUTS = ["water.xyz", "water_wrapped.xyz", "water.ff", "water.top", "nve.run"]

R_OH = 0.95720000000000
D_HH = 1.51390065452732
# Global atom tags are 0-based; the topology file lists them 1-based.
CONSTRAINTS = [(0, 1, R_OH), (0, 2, R_OH), (1, 2, D_HH)]

EXPECTED_DOF = 6          # 3N - constraints = 9 - 3, with the COM velocity KEPT
ATOM_COUNT = 3
TIME_STEP_FS = 0.5
TOTAL_STEPS = 400
# Steps at which the tight, full-precision residuals are sampled. `gmd` keeps one
# checkpoint file and overwrites it, so each sample is its own run of that length.
RESIDUAL_SAMPLE_STEPS = [100, 200, 300, 400]

# eV/A^3 -> bar, matching src/core/physical_constants.hpp.
EV_PER_A3_TO_BAR = 1.602176634e6


def run_gmd(gmd: str, work: pathlib.Path, run_file: str, np: int = 1,
            mpiexec: str = "mpiexec") -> str:
    command = []
    if np > 1:
        command += [mpiexec, "-np", str(np)]
    command += [gmd, "water.xyz", run_file, "water.ff", "water.top"]
    completed = subprocess.run(command, cwd=work, text=True, capture_output=True)
    if completed.returncode != 0:
        raise ValidationFailure(
            f"gmd failed (exit {completed.returncode}, np={np}) on {run_file} "
            f"in {work}\n--- stdout ---\n{completed.stdout}\n"
            f"--- stderr ---\n{completed.stderr}")
    return completed.stdout


def write_run(work: pathlib.Path, name: str, *, steps: int, checkpoint: str,
              xyz: str = "water.xyz", restart_from: str | None = None,
              base: str = "nve.run") -> None:
    """Derive a run file from the checked-in nve.run, changing only what varies."""
    lines = []
    for line in (work / base).read_text(encoding="utf-8").splitlines():
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
    if xyz != "water.xyz":
        raise ValidationFailure("the run file always names water.xyz on the CLI")


def parse_dof(stdout: str) -> int:
    for line in stdout.splitlines():
        if "Degrees of freedom:" in line:
            return int(line.split("Degrees of freedom:")[1].split()[0])
    raise ValidationFailure(
        "gmd did not report a degree-of-freedom count on stdout:\n" + stdout)


def parse_xyz_frame_zero(path: pathlib.Path) -> list[list[float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[0].strip())
    return [[float(v) for v in lines[2 + i].split()[1:4]] for i in range(count)]


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

    dof = parse_dof(stdout)
    record("constrained_dof", dof == EXPECTED_DOF, actual=dof, expected=EXPECTED_DOF,
           note="3N - constraints = 9 - 3, with the COM velocity kept")

    rows = parse_constrained_log(work / "output.log")
    write_series_csv(work / "thermo_series.csv", rows, list(rows[0].keys()))

    # ---- per-frame log checks -------------------------------------------
    # The log's residual columns carry six fixed decimals, so they can only
    # BOUND the residual, not measure it. That bound is still worth asserting on
    # every frame; the tight measurement happens from checkpoints below.
    log_bound = tolerance["log_residual_column_bound_a"]
    frame_failures = []
    for row in rows:
        if row["shake_error_a"] > log_bound:
            frame_failures.append(
                f"step {row['step']}: log column shake_error[A] = "
                f"{row['shake_error_a']} exceeds {log_bound}")
        if row["rattle_error_a_per_fs"] > log_bound:
            frame_failures.append(
                f"step {row['step']}: log column rattle_error[A/fs] = "
                f"{row['rattle_error_a_per_fs']} exceeds {log_bound}")
        if row["temperature"] <= 0.0:
            frame_failures.append(
                f"step {row['step']}: log column T[K] = {row['temperature']} "
                f"is not a positive finite temperature")
        if row["shake_iterations"] > reference["solver"]["max_iterations"]:
            frame_failures.append(
                f"step {row['step']}: shake_iter = {row['shake_iterations']} "
                f"exceeds the configured maximum")
    record("per_frame_log_bounds", not frame_failures, frames=len(rows),
           failures=frame_failures[:8],
           note="bounds every frame from the log; see checkpoint residuals for "
                "the tight measurement")

    # ---- pressure validity semantics ------------------------------------
    initial = rows[0]
    later = [row for row in rows if row["step"] > 0]
    record("initial_frame_pressure_invalid",
           initial["pressure_valid"] == 0 and initial["pressure_bar"] is None,
           step=initial["step"], pressure_valid=initial["pressure_valid"],
           note="no step has completed at frame 0, so no RATTLE multiplier "
                "belongs to the stored provider virial and the pressure must be "
                "reported invalid rather than as a number")
    record("completed_steps_pressure_valid",
           all(row["pressure_valid"] == 1 for row in later),
           frames=len(later),
           invalid=[row["step"] for row in later if row["pressure_valid"] != 1])

    # ---- initial projection, from the step-0 trajectory frame ------------
    # Six decimals in the .xyz, so this resolves 1e-6 A. The fixture starts
    # 0.25% (about 2.4e-3 A) off target, so a projection that did nothing would
    # miss by 2400x this bound. The tight check is the checkpoint one below.
    frame_zero = parse_xyz_frame_zero(work / "output.xyz")
    projection_failures = []
    for index, (tag_i, tag_j, target) in enumerate(CONSTRAINTS):
        distance = math.dist(frame_zero[tag_i], frame_zero[tag_j])
        error = abs(distance - target)
        if error > tolerance["initial_projection_xyz_abs_a"]:
            projection_failures.append(
                f"constraint {index} on atom tags ({tag_i}, {tag_j}): "
                f"measured {distance:.6f} A against target {target:.6f} A, "
                f"error {error:.3e} A")
    record("initial_target_distances_after_projection", not projection_failures,
           failures=projection_failures,
           tolerance=tolerance["initial_projection_xyz_abs_a"],
           note="measured from the step-0 .xyz frame, which carries six decimals")

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

    # ---- energy drift ----------------------------------------------------
    span_ps = (rows[-1]["time_fs"] - rows[0]["time_fs"]) * 1.0e-3
    drift = (rows[-1]["total_energy"] - rows[0]["total_energy"]) / (ATOM_COUNT * span_ps)
    record("energy_drift_per_atom_ps",
           abs(drift) <= tolerance["energy_drift_per_atom_ps_abs"],
           value=drift, tolerance=tolerance["energy_drift_per_atom_ps_abs"],
           note="a free rigid body has no potential energy, so this is the "
                "integrator's own drift and nothing else")

    # ---- analytic rigid-body results ------------------------------------
    final = parse_checkpoint(work / "base.gmdchk")
    first = parse_checkpoint(work / f"sample{RESIDUAL_SAMPLE_STEPS[0]}.gmdchk")

    momentum_initial = total_momentum(first)
    momentum_final = total_momentum(final)
    momentum_drift = max(abs(a - b) for a, b in zip(momentum_initial, momentum_final))
    momentum_scale = math.sqrt(sum(p * p for p in momentum_final)) or 1.0
    record("linear_momentum_conserved",
           momentum_drift / momentum_scale <= tolerance["momentum_rel"],
           relative_change=momentum_drift / momentum_scale,
           tolerance=tolerance["momentum_rel"],
           initial=momentum_initial, final=momentum_final)

    angular_initial = angular_momentum(first)
    angular_final = angular_momentum(final)
    angular_drift = max(abs(a - b) for a, b in zip(angular_initial, angular_final))
    angular_scale = math.sqrt(sum(l * l for l in angular_final)) or 1.0
    record("angular_momentum_conserved",
           angular_drift / angular_scale <= tolerance["angular_momentum_rel"],
           relative_change=angular_drift / angular_scale,
           tolerance=tolerance["angular_momentum_rel"],
           initial=angular_initial, final=angular_final,
           note="a free rigid body conserves L about its centre of mass exactly")

    # COM straight-line motion: r_com(t2) - r_com(t1) = v_com * (t2 - t1).
    #
    # Two wrinkles, both handled rather than assumed away. The molecule's atoms
    # are wrapped INDEPENDENTLY, so the centre of mass has to be formed from
    # positions unwrapped into a single image; and the centre of mass itself
    # drifts out of the cell over 150 fs, so the predicted and observed values
    # are compared under the minimum image convention.
    tags = [tag for tag, _, _ in [(0, None, None), (1, None, None), (2, None, None)]]
    velocity_com = center_of_mass_velocity(final)
    elapsed_fs = final["time_fs"] - first["time_fs"]
    # Velocities are in A per INTERNAL time unit, so the elapsed time has to be
    # converted before it can multiply one. Getting this wrong is a factor of
    # 10.18, which is why it is spelled out rather than inlined.
    elapsed = elapsed_fs * INTERNAL_TIME_PER_FEMTOSECOND
    predicted = [group_center_of_mass(first, tags)[d] + velocity_com[d] * elapsed
                 for d in range(3)]
    observed = group_center_of_mass(final, tags)
    com_error = separation_modulo_box(predicted, observed, final["box"])
    record("com_straight_line_motion",
           com_error <= tolerance["com_trajectory_abs_a"],
           max_component_error_a=com_error,
           tolerance=tolerance["com_trajectory_abs_a"],
           elapsed_fs=elapsed_fs, predicted=predicted, observed=observed,
           note="analytic: a force-free body's COM moves at constant velocity")

    record("finite_temperature",
           all(math.isfinite(row["temperature"]) and row["temperature"] > 0.0
               for row in rows),
           minimum=min(row["temperature"] for row in rows),
           maximum=max(row["temperature"] for row in rows))

    # ---- constraint virial and pressure ----------------------------------
    constraint_virial = final["constraint_virial"]
    record("constraint_virial_state_valid",
           constraint_virial["state"] == "valid"
           and constraint_virial["time_level"] == "endpoint_rattle_t_plus_dt",
           state=constraint_virial["state"],
           time_level=constraint_virial["time_level"],
           note="the stored multiplier must be the endpoint RATTLE one, at the "
                "same time level as the provider virial it is added to")

    tensor = constraint_virial["tensor"]
    check_tensor_finite(tensor, "constraint virial")
    record("constraint_virial_components_finite", True,
           tensor={name: value for name, value in zip(TENSOR_COMPONENT_NAMES, tensor)})
    record("constraint_virial_symmetric",
           max_asymmetry(tensor) <= tolerance["virial_symmetry_abs"],
           max_asymmetry=max_asymmetry(tensor),
           tolerance=tolerance["virial_symmetry_abs"])
    smallest_offdiagonal = min(abs(tensor[i]) for i in (1, 2, 5))
    record("constraint_virial_offdiagonals_nontrivial",
           smallest_offdiagonal >= tolerance["virial_offdiagonal_min_abs"],
           smallest_offdiagonal=smallest_offdiagonal,
           threshold=tolerance["virial_offdiagonal_min_abs"],
           note="the molecule is generically rotated, so no off-diagonal may "
                "vanish; a fixture aligned with the box axes would hide errors")

    step_pressure = final["step_pressure"]
    if not step_pressure["valid"]:
        raise ValidationFailure(
            "the completed-step pressure must be valid at the end of a "
            "constrained NVE run, but the checkpoint reports it invalid")
    identity = (step_pressure["twice_kinetic_energy"] + trace(step_pressure["tensor"])) \
        / (3.0 * step_pressure["volume"])
    record("pressure_kinetic_plus_virial_identity",
           abs(identity - step_pressure["pressure"])
           <= tolerance["pressure_identity_rel"] * abs(step_pressure["pressure"]),
           reported=step_pressure["pressure"], recomputed=identity,
           relative_difference=abs(identity - step_pressure["pressure"])
           / abs(step_pressure["pressure"]),
           tolerance=tolerance["pressure_identity_rel"],
           note="P = (2K + tr W) / 3V, recomputed from the checkpoint's own "
                "2K, volume and total virial")

    # THE ANALYTIC PRESSURE. For a free rigid body the constraint virial cancels
    # the rotational kinetic energy exactly, leaving the ideal-gas pressure of a
    # single particle of the molecule's total mass moving at v_com.
    total_mass = sum(a["mass"] for a in final["atoms"])
    twice_ke_com = total_mass * sum(v * v for v in velocity_com)
    twice_ke_total = twice_kinetic_energy(final)
    twice_ke_rotational = twice_ke_total - twice_ke_com
    cancellation = abs(twice_ke_rotational + trace(tensor))
    record("constraint_virial_cancels_rotational_ke",
           cancellation <= tolerance["rotational_cancellation_abs"],
           twice_ke_rotational=twice_ke_rotational,
           trace_constraint_virial=trace(tensor),
           residual=cancellation,
           tolerance=tolerance["rotational_cancellation_abs"],
           note="analytic: a rigid body's frozen internal motion must not "
                "contribute to the pressure, so tr W = -2K_rot to O(dt^2)")

    analytic_pressure = twice_ke_com / (3.0 * step_pressure["volume"])
    pressure_relative_error = abs(step_pressure["pressure"] - analytic_pressure) \
        / abs(analytic_pressure)
    record("pressure_matches_analytic_ideal_gas_com",
           pressure_relative_error <= tolerance["analytic_pressure_rel"],
           reported_ev_per_a3=step_pressure["pressure"],
           analytic_ev_per_a3=analytic_pressure,
           reported_bar=step_pressure["pressure"] * EV_PER_A3_TO_BAR,
           analytic_bar=analytic_pressure * EV_PER_A3_TO_BAR,
           relative_error=pressure_relative_error,
           tolerance=tolerance["analytic_pressure_rel"],
           note="analytic reference: M |v_com|^2 / 3V, independent of GMD")

    # ---- wrapped / unwrapped equivalence ---------------------------------
    wrapped = work / "wrapped"
    wrapped.mkdir(exist_ok=True)
    for name in INPUTS:
        shutil.copy2(case_dir / name, wrapped / name)
    shutil.copy2(case_dir / "water_wrapped.xyz", wrapped / "water.xyz")
    write_run(wrapped, "base.run", steps=TOTAL_STEPS, checkpoint="base.gmdchk")
    run_gmd(args.gmd, wrapped, "base.run")
    wrapped_state = parse_checkpoint(wrapped / "base.gmdchk")
    wrapped_rows = parse_constrained_log(wrapped / "output.log")

    # The wrapped fixture is the SAME molecule translated so that it straddles a
    # periodic face. Coordinates therefore differ, but every quantity that does
    # not depend on the choice of origin must agree.
    wrapped_residuals = constraint_residuals(wrapped_state, CONSTRAINTS)

    # THE FIXTURE MUST ACTUALLY STRADDLE A FACE, and at THIS step -- the one the
    # comparison below looks at -- not merely at t=0. The molecule carries a
    # centre-of-mass velocity, so a variant that starts on the face drifts off
    # it, after which every raw distance equals its minimum-image distance and
    # the wrapping path is not exercised at all. The case would keep passing
    # while testing nothing. This asserts the precondition instead of assuming
    # it: at least one constrained pair must have its atoms in different images.
    straddling = []
    for index, (tag_i, tag_j, _target) in enumerate(CONSTRAINTS):
        by_tag = {atom["tag"]: atom for atom in wrapped_state["atoms"]}
        delta = [by_tag[tag_i]["position"][d] - by_tag[tag_j]["position"][d]
                 for d in range(3)]
        raw = math.sqrt(sum(v * v for v in delta))
        imaged = minimum_image(delta, wrapped_state["box"])
        folded = math.sqrt(sum(v * v for v in imaged))
        if abs(raw - folded) > 1.0e-6:
            straddling.append(index)
    record("wrapped_fixture_actually_straddles_a_face", bool(straddling),
           straddling_constraints=straddling,
           note="precondition for the equivalence check below; without it the "
                "minimum-image path is never taken and the comparison is vacuous")

    scalar_differences = []
    for key in ("total_energy", "ke", "temperature"):
        difference = abs(wrapped_rows[-1][key] - rows[-1][key])
        if difference > tolerance["wrapped_scalar_abs"]:
            scalar_differences.append(
                f"log column {key}: wrapped {wrapped_rows[-1][key]} vs "
                f"unwrapped {rows[-1][key]} (difference {difference:.3e})")
    virial_match = compare_tensors(
        wrapped_state["constraint_virial"]["tensor"], tensor,
        tolerance["wrapped_virial_abs"], "wrapped vs unwrapped constraint virial")
    record("wrapped_unwrapped_equivalence",
           not scalar_differences and virial_match["passed"],
           scalar_differences=scalar_differences,
           virial=virial_match,
           wrapped_max_position_residual_a=
               wrapped_residuals["max_position_residual_a"]["value"],
           note="the molecule straddles a periodic face; energies, temperature "
                "and the constraint virial are origin independent and must match")

    # ---- checkpoint / restart continuity ---------------------------------
    half = TOTAL_STEPS // 2
    write_run(work, "restart_first.run", steps=half, checkpoint="restart_half.gmdchk")
    run_gmd(args.gmd, work, "restart_first.run")
    write_run(work, "restart_second.run", steps=half,
              checkpoint="restart_final.gmdchk", restart_from="restart_half.gmdchk")
    run_gmd(args.gmd, work, "restart_second.run")
    restarted = parse_checkpoint(work / "restart_final.gmdchk")

    if restarted["step"] != final["step"]:
        raise ValidationFailure(
            f"restarted run ended at step {restarted['step']} but the continuous "
            f"run ended at step {final['step']}")
    position_difference, velocity_difference, worst_tag = 0.0, 0.0, None
    by_tag = {a["tag"]: a for a in restarted["atoms"]}
    for atom in final["atoms"]:
        other = by_tag[atom["tag"]]
        for d in range(3):
            dp = abs(atom["position"][d] - other["position"][d])
            dv = abs(atom["velocity"][d] - other["velocity"][d])
            if dp > position_difference:
                position_difference, worst_tag = dp, atom["tag"]
            velocity_difference = max(velocity_difference, dv)
    restart_virial = compare_tensors(
        restarted["constraint_virial"]["tensor"], tensor,
        tolerance["restart_virial_abs"], "restart vs continuous constraint virial")
    record("restart_continuity",
           position_difference <= tolerance["restart_position_abs_a"]
           and velocity_difference <= tolerance["restart_velocity_abs_a_per_fs"]
           and restart_virial["passed"]
           and restarted["constraint_virial"]["state"] == "valid",
           max_position_difference_a=position_difference,
           max_velocity_difference_a_per_fs=velocity_difference,
           worst_atom_tag=worst_tag,
           constraint_virial=restart_virial,
           restored_state=restarted["constraint_virial"]["state"],
           note=f"{half} + {half} steps against {TOTAL_STEPS} continuous")

    # ---- MPI -------------------------------------------------------------
    ranks = [int(r) for r in args.mpi_ranks.split(",") if r.strip()] if args.mpi_ranks else []
    if args.mpiexec and ranks:
        mpi_results = {}
        differences = []
        for np in ranks:
            directory = work / f"np{np}"
            directory.mkdir(exist_ok=True)
            for name in INPUTS:
                shutil.copy2(case_dir / name, directory / name)
            write_run(directory, "base.run", steps=TOTAL_STEPS,
                      checkpoint="base.gmdchk")
            run_gmd(args.gmd, directory, "base.run", np=np, mpiexec=args.mpiexec)
            state = parse_checkpoint(directory / "base.gmdchk")
            mpi_results[np] = state
            comparison = compare_tensors(
                state["constraint_virial"]["tensor"], tensor,
                tolerance["mpi_virial_abs"],
                f"np={np} vs serial constraint virial")
            if not comparison["passed"]:
                differences.append(
                    f"np={np}: constraint virial component "
                    f"{comparison['worst_component']} differs by "
                    f"{comparison['max_component_difference']:.3e} "
                    f"(tolerance {tolerance['mpi_virial_abs']:.3e})")
            worst = 0.0
            others = {a["tag"]: a for a in state["atoms"]}
            for atom in final["atoms"]:
                for d in range(3):
                    worst = max(worst,
                                abs(atom["position"][d] - others[atom["tag"]]["position"][d]))
            if worst > tolerance["mpi_position_abs_a"]:
                differences.append(
                    f"np={np}: positions differ from serial by {worst:.3e} A")
        record("mpi_agreement", not differences, ranks=ranks,
               differences=differences,
               note="the constraint virial is already global; agreement without "
                    "a factor of the rank count is the point of this check")
    else:
        checks["mpi_agreement"] = {"passed": True, "skipped": True,
                                   "note": "no mpiexec supplied"}

    # ------------------------------------------------------------------ done
    summary = {
        "case": "constrained_nve_water",
        "ensemble": "NVE",
        "baseline_kind": reference["provenance"]["baseline_kind"],
        "checks": checks,
        "passed": all(check["passed"] for check in checks.values()),
    }
    write_json(work / "summary.json", summary)

    failed = [name for name, check in checks.items() if not check["passed"]]
    if failed:
        sys.stderr.write("FAILED constrained_nve_water checks: "
                         + ", ".join(failed) + "\n")
        for name in failed:
            sys.stderr.write(f"  {name}: "
                             + json.dumps(checks[name], indent=2, sort_keys=True) + "\n")
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationFailure as failure:
        sys.stderr.write(f"constrained_nve_water: {failure}\n")
        raise SystemExit(1)
