#!/usr/bin/env python3
"""External-controller integration test for the `gmd` CLI.

Drives the real executable the way an active-learning controller does: in
segments, reading frames and status files, requesting checkpoints and stops.

Classical part (always; standard library only, ASE checks when importable):
  * continuous 2N steps vs N + N segmented restart (v3 checkpoint): final
    checkpoint state, energy log and frame-file frames must agree BITWISE;
  * restart with a changed force field is refused without `model_change
    allow`, and accepted and recorded with it;
  * the gmd-frame-v1 frame file: structure, 0-based ids, and forces equal to
    an independent Lennard-Jones evaluation; parsed with ase.io.read(":")
    when ASE is importable (required with --require-ase);
  * run-file-relative paths: gmd started from another working directory;
  * stop_file, checkpoint_request_file and SIGTERM on a long run, and the
    status file contents in every case.

ML part (with --ml; needs a Torch-enabled gmd, torch, ase and GMD_SGT):
  * exports two small GMD_SGT models (TorchScript);
  * continuous 2N vs N + N segmented restart with the same model;
  * restart with model B refused without the directive, accepted and
    recorded (old -> new sha) with it;
  * a pre-v3 (no atomic numbers) checkpoint is refused for ML;
  * atoms of unknown species (Z = 0) are refused before any step;
  * frame forces equal a Python re-evaluation of the exported model.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True

problems: list[str] = []


def check(condition: bool, message: str) -> bool:
    if not condition:
        problems.append(message)
        print(f"  FAIL: {message}")
    return condition


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------
# Running gmd
# --------------------------------------------------------------------------

def run_gmd(gmd: Path, xyz: Path, run: Path, *, cwd: Path | None = None,
            expect_ok: bool = True, env: dict | None = None) -> subprocess.CompletedProcess:
    result = subprocess.run([str(gmd), str(xyz), str(run)], cwd=cwd or run.parent,
                            text=True, capture_output=True, env=env, timeout=600)
    if expect_ok and result.returncode != 0:
        raise RuntimeError(f"gmd failed ({result.returncode}) for {run}:\n"
                           f"{result.stdout}\n{result.stderr}")
    return result


def start_gmd(gmd: Path, xyz: Path, run: Path, env: dict | None = None) -> subprocess.Popen:
    return subprocess.Popen([str(gmd), str(xyz), str(run)], cwd=run.parent, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)


def wait_for(predicate, timeout: float, what: str) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.02)
    check(False, f"timed out waiting for {what}")
    return False


def read_status(path: Path) -> dict:
    return json.loads(path.read_text())


# --------------------------------------------------------------------------
# Parsing GMD outputs (standard library only)
# --------------------------------------------------------------------------

def read_checkpoint(path: Path) -> dict:
    lines = path.read_text().splitlines()
    data: dict = {"version": int(lines[0].split()[1]), "atoms": {}}
    index = 1
    while index < len(lines):
        line = lines[index]
        if line.startswith("atoms tag "):
            columns = line.split()[1:]
            count = int(data["atom_count"])
            for atom_line in lines[index + 1:index + 1 + count]:
                fields = atom_line.split()
                data["atoms"][int(fields[0])] = dict(zip(columns, fields))
            index += count + 1
            continue
        key, _, value = line.partition(" ")
        data[key] = value
        index += 1
    return data


def parse_frames(path: Path) -> list[dict]:
    """Minimal gmd-frame-v1 reader, independent of ASE."""
    lines = path.read_text().splitlines()
    frames = []
    index = 0
    while index < len(lines):
        n = int(lines[index])
        comment = lines[index + 1]
        info = {}
        token = ""
        in_quotes = False
        for char in comment + " ":
            if char == '"':
                in_quotes = not in_quotes
            elif char == " " and not in_quotes:
                if token:
                    key, _, value = token.partition("=")
                    info[key] = value
                token = ""
                continue
            token += char if char != '"' else ""
        atoms = [line.split() for line in lines[index + 2:index + 2 + n]]
        check(len(atoms) == n, f"{path.name}: truncated frame at line {index + 1}")
        frames.append({"info": info, "atoms": atoms, "text": "\n".join(lines[index:index + 2 + n])})
        index += 2 + n
    return frames


def read_log(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if not line.startswith("#")]


def frame_by_step(frames: list[dict]) -> dict[int, dict]:
    return {int(frame["info"]["step"]): frame for frame in frames}


# --------------------------------------------------------------------------
# Independent Lennard-Jones reference (minimum image, plain cutoff)
# --------------------------------------------------------------------------

def lj_forces(positions, box, epsilon, sigma, cutoff):
    n = len(positions)
    forces = [[0.0, 0.0, 0.0] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = [positions[i][k] - positions[j][k] for k in range(3)]
            d = [d[k] - box[k] * round(d[k] / box[k]) for k in range(3)]
            r2 = sum(x * x for x in d)
            if r2 >= cutoff * cutoff:
                continue
            s6 = (sigma * sigma / r2) ** 3
            f_over_r = 24.0 * epsilon * (2.0 * s6 * s6 - s6) / r2
            for k in range(3):
                forces[i][k] += f_over_r * d[k]
                forces[j][k] -= f_over_r * d[k]
    return forces


# --------------------------------------------------------------------------
# Classical (LJ) fixtures
# --------------------------------------------------------------------------

LJ_EPS = 0.01032
LJ_SIGMA = 3.405
LJ_CUTOFF = 8.5
LJ_BOX = 18.0


def lj_xyz(path: Path) -> None:
    # 3 x 3 x 3 argon lattice, slightly perturbed so forces are non-trivial.
    lines = ["27", f"{LJ_BOX} {LJ_BOX} {LJ_BOX}"]
    for a in range(3):
        for b in range(3):
            for c in range(3):
                jitter = 0.05 * math.sin(7.0 * a + 3.0 * b + 1.3 * c)
                lines.append(f"Ar {1.0 + 6.0 * a + jitter:.6f} {1.0 + 6.0 * b - jitter:.6f} "
                             f"{1.0 + 6.0 * c + 0.5 * jitter:.6f}")
    path.write_text("\n".join(lines) + "\n")


def lj_run(path: Path, steps: int, *, epsilon: float = LJ_EPS, extra: list[str] = ()) -> None:
    path.write_text("\n".join([
        f"run {steps}",
        "time_step 2.0",
        "velocity 120.0",
        "velocity_init random",
        "velocity_seed 77",
        "output_interval 5",
        "force_field lj",
        f"cutoff {LJ_CUTOFF}",
        f"type 1 Ar epsilon {epsilon} sigma {LJ_SIGMA} charge 0.0",
        *extra,
    ]) + "\n")


def compare_states(label: str, a: dict, b: dict, tolerance: float) -> float:
    worst = 0.0
    check(set(a["atoms"]) == set(b["atoms"]), f"{label}: atom tags differ")
    for tag in a["atoms"]:
        for key in ("x", "y", "z", "vx", "vy", "vz"):
            diff = abs(float(a["atoms"][tag][key]) - float(b["atoms"][tag][key]))
            worst = max(worst, diff)
        if "atomic_number" in a["atoms"][tag] and "atomic_number" in b["atoms"][tag]:
            check(a["atoms"][tag]["atomic_number"] == b["atoms"][tag]["atomic_number"],
                  f"{label}: atomic number of tag {tag} differs")
    check(worst <= tolerance, f"{label}: state differs by {worst:.3e} > {tolerance:.1e}")
    return worst


def segment_test(gmd: Path, root: Path, make_inputs, steps: int, label: str,
                 tolerance: float, env: dict | None = None) -> dict:
    """Continuous 2N vs N + N. make_inputs(dir, steps, extra) writes in.xyz/run.in."""
    cont = root / f"{label}_continuous"
    seg1 = root / f"{label}_seg1"
    seg2 = root / f"{label}_seg2"
    for d in (cont, seg1, seg2):
        d.mkdir(parents=True)
    common = ["frame_file frames.extxyz", "frame_interval 5", "frame_velocities on",
              f"run_id {label}", "status_file status.json"]
    make_inputs(cont, 2 * steps, common + ["checkpoint_file final.gmdchk"])
    run_gmd(gmd, cont / "in.xyz", cont / "run.in", env=env)
    make_inputs(seg1, steps, common + ["checkpoint_file half.gmdchk"])
    run_gmd(gmd, seg1 / "in.xyz", seg1 / "run.in", env=env)
    make_inputs(seg2, steps, common + ["checkpoint_file final.gmdchk",
                                       f"restart_from {seg1 / 'half.gmdchk'}",
                                       "velocity_init input"])
    run_gmd(gmd, seg2 / "in.xyz", seg2 / "run.in", env=env)

    half = read_checkpoint(seg1 / "half.gmdchk")
    check(half["version"] == 3, f"{label}: segment checkpoint is not version 3")
    worst = compare_states(f"{label} continuous vs segmented", read_checkpoint(cont / "final.gmdchk"),
                           read_checkpoint(seg2 / "final.gmdchk"), tolerance)

    # Energies: the log rows of the second segment must equal the continuous
    # run's rows for the same steps (6 printed decimals, compared as text).
    # With a non-zero tolerance, rows are compared numerically to within the
    # printed resolution (1e-6) plus the tolerance.
    cont_log = {row.split()[0]: row for row in read_log(cont / "output.log")}
    seg_rows = read_log(seg1 / "output.log") + read_log(seg2 / "output.log")

    def rows_agree(a: str | None, b: str) -> bool:
        if a is None:
            return False
        if tolerance == 0.0:
            return a == b
        return all(x == y or abs(float(x) - float(y)) <= 1.5e-6 + tolerance
                   for x, y in zip(a.split(), b.split()))

    mismatched = [row.split()[0] for row in seg_rows
                  if not rows_agree(cont_log.get(row.split()[0]), row)]
    check(not mismatched, f"{label}: energy log rows differ at steps {mismatched[:5]}")

    # Frames at full precision: positions, forces, velocities and energies.
    cont_frames = frame_by_step(parse_frames(cont / "frames.extxyz"))
    seg_frames = frame_by_step(parse_frames(seg1 / "frames.extxyz"))
    seg_frames.update(frame_by_step(parse_frames(seg2 / "frames.extxyz")))
    check(sorted(cont_frames) == sorted(seg_frames),
          f"{label}: frame steps differ: {sorted(cont_frames)} vs {sorted(seg_frames)}")
    worst_frame = 0.0
    worst_energy = 0.0
    for step, frame in cont_frames.items():
        other = seg_frames.get(step)
        if other is None:
            continue
        for row_a, row_b in zip(frame["atoms"], other["atoms"]):
            check(row_a[:2] == row_b[:2], f"{label}: id/species differ at step {step}")
            for x, y in zip(row_a[2:], row_b[2:]):
                worst_frame = max(worst_frame, abs(float(x) - float(y)))
        for key in ("potential_energy", "kinetic_energy"):
            worst_energy = max(worst_energy, abs(float(frame["info"][key]) -
                                                 float(other["info"][key])))
    check(worst_frame <= tolerance, f"{label}: frame values differ by {worst_frame:.3e}")
    check(worst_energy <= tolerance, f"{label}: frame energies differ by {worst_energy:.3e}")
    status = read_status(seg2 / "status.json")
    check(status["state"] == "completed" and status["step"] == 2 * steps and
          status["start_step"] == steps,
          f"{label}: unexpected segment-2 status {status}")
    print(f"  {label}: continuous {2 * steps} vs {steps}+{steps}: state max |diff| = {worst:.3e}, "
          f"frames max |diff| = {worst_frame:.3e}, energies max |diff| = {worst_energy:.3e}"
          + ("  (bitwise)" if worst == worst_frame == worst_energy == 0.0 else ""))
    return {"cont": cont, "seg1": seg1, "seg2": seg2}


# --------------------------------------------------------------------------
# Classical tests
# --------------------------------------------------------------------------

def lj_inputs(directory: Path, steps: int, extra: list[str], epsilon: float = LJ_EPS) -> None:
    lj_xyz(directory / "in.xyz")
    lj_run(directory / "run.in", steps, epsilon=epsilon, extra=extra)


def lj_multicell_xyz(path: Path) -> None:
    # 5 x 5 x 5 argon in a 33 A box: a 3 x 3 x 3 neighbour-cell grid.
    lines = ["125", "33.0 33.0 33.0"]
    for a in range(5):
        for b in range(5):
            for c in range(5):
                jitter = 0.05 * math.sin(7.0 * a + 3.0 * b + 1.3 * c)
                lines.append(f"Ar {0.5 + 6.6 * a + jitter:.6f} {0.5 + 6.6 * b - jitter:.6f} "
                             f"{0.5 + 6.6 * c + 0.5 * jitter:.6f}")
    path.write_text("\n".join(lines) + "\n")


def lj_multicell_inputs(directory: Path, steps: int, extra: list[str]) -> None:
    lj_multicell_xyz(directory / "in.xyz")
    lj_run(directory / "run.in", steps, extra=extra)


def test_lj_segments(gmd: Path, root: Path, steps: int) -> dict:
    print("[LJ] continuous vs segmented restart")
    # A single neighbour cell (18 A box, 10.5 A list radius): the Verlet list
    # order cannot depend on history, so the split run is bitwise identical.
    dirs = segment_test(gmd, root, lj_inputs, steps, "lj", 0.0)
    # Several cells: the classical Verlet list visits neighbours in an order
    # that depends on the cells atoms occupied when it was last built, and a
    # restart rebuilds it at step N. The float64 force sums can then differ in
    # the last bit (~1e-16 relative) and chaotic dynamics amplify that; measured
    # 3.6e-14 A after 1500 + 1500 steps at 600 K. 1e-9 bounds that with margin
    # and is far below anything physical.
    segment_test(gmd, root, lj_multicell_inputs, steps, "lj_multicell", 1e-9)
    return dirs


def test_lj_model_change(gmd: Path, root: Path, dirs: dict, steps: int) -> None:
    print("[LJ] force-field change on restart")
    checkpoint = dirs["seg1"] / "half.gmdchk"
    old_sha = read_checkpoint(checkpoint)["force_field_sha256"]
    refused = root / "lj_change_refused"
    refused.mkdir()
    lj_inputs(refused, steps, [f"restart_from {checkpoint}", "status_file status.json"],
              epsilon=2 * LJ_EPS)
    result = run_gmd(gmd, refused / "in.xyz", refused / "run.in", expect_ok=False)
    check(result.returncode != 0, "LJ restart with a changed epsilon was accepted")
    check("model_change allow" in result.stderr,
          f"refusal does not name the directive: {result.stderr.strip()}")
    status = read_status(refused / "status.json")
    check(status["state"] == "failed" and "differs from the checkpoint" in (status["message"] or ""),
          f"refused restart status is not 'failed' with the reason: {status}")

    allowed = root / "lj_change_allowed"
    allowed.mkdir()
    lj_inputs(allowed, steps, [f"restart_from {checkpoint}", "model_change allow",
                               "frame_file frames.extxyz", "status_file status.json",
                               "checkpoint_file final.gmdchk"], epsilon=2 * LJ_EPS)
    result = run_gmd(gmd, allowed / "in.xyz", allowed / "run.in")
    check("Model change on restart" in result.stdout and old_sha in result.stdout,
          "allowed model change is not reported on stdout")
    log_comments = [l for l in (allowed / "output.log").read_text().splitlines()
                    if l.startswith("# model_change")]
    check(len(log_comments) == 1 and old_sha in log_comments[0],
          f"output.log does not record the model change: {log_comments}")
    frames = parse_frames(allowed / "frames.extxyz")
    check(all(f["info"].get("model_changed_from") == old_sha for f in frames),
          "frames do not record model_changed_from=<old sha>")
    status = read_status(allowed / "status.json")
    check(status["model_changed_from"] == old_sha, f"status lacks model_changed_from: {status}")
    new_sha = read_checkpoint(allowed / "final.gmdchk")["force_field_sha256"]
    check(new_sha != old_sha and status["force_field_sha256"] == new_sha,
          "the new checkpoint does not carry the new force-field fingerprint")
    print(f"  refused without directive; accepted with it: {old_sha[:12]} -> {new_sha[:12]}")


def test_lj_frames(gmd: Path, dirs: dict, require_ase: bool) -> None:
    print("[LJ] frame file contents")
    path = dirs["cont"] / "frames.extxyz"
    frames = parse_frames(path)
    worst = 0.0
    for frame in frames:
        info = frame["info"]
        check(info["Lattice"] == "18.0 0 0 0 18.0 0 0 0 18.0", f"bad Lattice {info['Lattice']}")
        check(info["pbc"] == "T T T", "bad pbc")
        check(info["Properties"] ==
              "id:I:1:species:S:1:pos:R:3:forces:R:3:vel:R:3", "bad Properties")
        check(info["format"] == "gmd-frame-v1", "bad format")
        check(info["run_id"] == "lj" and info["model_sha256"] == "none" and
              info["force_field"] == "inline_lj", f"bad provenance {info}")
        for key in ("step", "time_fs", "energy", "potential_energy", "kinetic_energy",
                    "temperature", "volume", "units"):
            check(key in info, f"frame lacks info key {key}")
        check([int(a[0]) for a in frame["atoms"]] == list(range(27)), "ids are not 0..N-1")
        check(all(a[1] == "Ar" for a in frame["atoms"]), "species are not Ar")
        positions = [[float(v) for v in a[2:5]] for a in frame["atoms"]]
        reference = lj_forces(positions, [LJ_BOX] * 3, LJ_EPS, LJ_SIGMA, LJ_CUTOFF)
        for row, ref in zip(frame["atoms"], reference):
            for k in range(3):
                worst = max(worst, abs(float(row[5 + k]) - ref[k]))
    check(worst < 1e-12, f"frame forces differ from the LJ reference by {worst:.3e}")
    print(f"  {len(frames)} frames; forces vs independent LJ: max |diff| = {worst:.3e} eV/A")
    check_frames_with_ase(path, frames, require_ase)


def check_frames_with_ase(path: Path, frames: list[dict], require_ase: bool) -> None:
    try:
        import ase.io  # noqa: PLC0415
        import numpy as np  # noqa: PLC0415
    except ImportError:
        check(not require_ase, "ASE is required (--require-ase) but not importable")
        print("  ASE not importable: ase.io.read check SKIPPED")
        return
    images = ase.io.read(path, index=":")
    check(len(images) == len(frames), f"ASE read {len(images)} frames, expected {len(frames)}")
    worst = 0.0
    for atoms, frame in zip(images, frames):
        text_forces = np.array([[float(v) for v in a[5:8]] for a in frame["atoms"]])
        text_pos = np.array([[float(v) for v in a[2:5]] for a in frame["atoms"]])
        worst = max(worst, float(np.abs(atoms.get_forces() - text_forces).max()))
        check(np.array_equal(atoms.get_positions(), text_pos), "ASE positions differ from text")
        check(atoms.get_potential_energy() == float(frame["info"]["energy"]),
              "ASE potential energy differs from the frame")
        check(list(atoms.arrays["id"]) == [int(a[0]) for a in frame["atoms"]], "ASE ids differ")
        check(atoms.get_chemical_symbols() == [a[1] for a in frame["atoms"]], "ASE species differ")
        check(bool(atoms.pbc.all()), "ASE pbc is not all True")
        check(atoms.info["step"] == int(frame["info"]["step"]), "ASE step differs")
    check(worst == 0.0, f"ASE forces differ from GMD's written forces by {worst:.3e}")
    print(f"  ase.io.read(index=':') -> {len(images)} frames; forces identical to GMD's "
          f"(max |diff| = {worst:.1e})")


def test_relative_paths(gmd: Path, root: Path) -> None:
    print("[LJ] run-file-relative paths")
    run_dir = root / "relative" / "case"
    run_dir.mkdir(parents=True)
    lj_inputs(run_dir, 4, ["frame_file out/frames.extxyz", "status_file status.json",
                           "checkpoint_file state.gmdchk"])
    (run_dir / "out").mkdir()
    # Start from the work root, naming the inputs by relative paths.
    run_gmd(gmd, Path("relative/case/in.xyz"), Path("relative/case/run.in"), cwd=root)
    check((run_dir / "out" / "frames.extxyz").exists(), "frame_file not created next to run.in")
    check((run_dir / "status.json").exists(), "status_file not created next to run.in")
    check((run_dir / "state.gmdchk").exists(), "checkpoint_file not created next to run.in")
    check(not (root / "status.json").exists() and not (root / "state.gmdchk").exists(),
          "a relative path was resolved against the working directory")
    status = read_status(run_dir / "status.json")
    check(Path(status["checkpoint"]) == (run_dir / "state.gmdchk").resolve()
          or Path(status["checkpoint"]).resolve() == (run_dir / "state.gmdchk").resolve(),
          f"status checkpoint path is not the resolved absolute path: {status['checkpoint']}")
    print("  frame/status/checkpoint files land next to run.in when started from elsewhere")


def control_case(gmd: Path, root: Path, name: str, poll: int, extra: list[str]) -> tuple[Path, subprocess.Popen]:
    directory = root / name
    directory.mkdir()
    lj_inputs(directory, 100_000_000, ["frame_file frames.extxyz", "frame_interval 10",
                                       "status_file status.json", "checkpoint_file state.gmdchk",
                                       f"control_poll_interval {poll}", "run_id control",
                                       *extra])
    process = start_gmd(gmd, directory / "in.xyz", directory / "run.in")
    status_path = directory / "status.json"
    wait_for(lambda: status_path.exists() and read_status(status_path)["state"] == "running",
             60, f"{name}: running status")
    return directory, process


def finish(process: subprocess.Popen, timeout: float = 60) -> tuple[int, str, str]:
    stdout, stderr = process.communicate(timeout=timeout)
    return process.returncode, stdout, stderr


def check_stopped_run(directory: Path, rc: int, stdout: str, reason: str, poll: int) -> None:
    check(rc == 0, f"{directory.name}: exit code {rc}, expected 0")
    status = read_status(directory / "status.json")
    step = status["step"]
    check(status["state"] == "stopped", f"{directory.name}: state {status['state']}")
    check(status["stop_reason"] == reason, f"{directory.name}: stop_reason {status['stop_reason']}")
    check(step % poll == 0, f"{directory.name}: stopped at step {step}, not a poll boundary")
    check(f"stopped by request at step {step}" in stdout,
          f"{directory.name}: no 'stopped by request at step {step}' line on stdout")
    log_text = (directory / "output.log").read_text()
    check(f"# stopped by request at step {step}" in log_text,
          f"{directory.name}: output.log lacks the stop line")
    checkpoint = read_checkpoint(Path(status["checkpoint"]))
    check(int(checkpoint["step"]) == step, f"{directory.name}: checkpoint step != status step")
    frames = parse_frames(directory / "frames.extxyz")
    check(int(frames[-1]["info"]["step"]) == step, f"{directory.name}: last frame is not step {step}")
    check(status["frames_written"] == len(frames), f"{directory.name}: frames_written mismatch")
    check(status["frame_file"] and Path(status["frame_file"]).resolve() ==
          (directory / "frames.extxyz").resolve(), f"{directory.name}: frame_file in status")
    check(abs(status["time_fs"] - 2.0 * step) < 1e-9, f"{directory.name}: time_fs wrong")
    check(status["format_version"] == "gmd-status-v1" and status["checkpoint_version"] == 3,
          f"{directory.name}: status format fields wrong")
    # Legacy output ends on a complete frame too.
    xyz_lines = (directory / "output.xyz").read_text().splitlines()
    check(len(xyz_lines) % 29 == 0, f"{directory.name}: output.xyz ends mid-frame")
    print(f"  {directory.name}: stopped at step {step} ({reason}); status, checkpoint, "
          f"final frame and logs consistent")


def test_control(gmd: Path, root: Path) -> None:
    print("[LJ] cooperative control")
    poll = 7
    directory, process = control_case(gmd, root, "control_stop", poll, ["stop_file STOP"])
    # Tail the frame file while the run is live: every snapshot must parse up
    # to its last complete frame, and complete frames must be well formed.
    frames_path = directory / "frames.extxyz"
    for _ in range(20):
        text = frames_path.read_text() if frames_path.exists() else ""
        lines = text.split("\n")
        index = 0
        while index < len(lines) and lines[index]:
            n = int(lines[index])
            if index + 2 + n > len(lines) - (0 if text.endswith("\n") else 1):
                break  # an in-progress trailing frame; ignored as documented
            index += 2 + n
        time.sleep(0.01)
    (directory / "STOP").touch()
    rc, stdout, _ = finish(process)
    check_stopped_run(directory, rc, stdout, "stop_file", poll)

    directory, process = control_case(gmd, root, "control_request", 5,
                                      ["checkpoint_request_file CHECKPOINT",
                                       "checkpoint_step_suffix on"])
    request = directory / "CHECKPOINT"
    request.touch()
    wait_for(lambda: not request.exists(), 60, "checkpoint request acknowledgement")
    status = read_status(directory / "status.json")
    check(status["state"] == "running", f"status after request is {status['state']}")
    requested = Path(status["checkpoint"]) if status["checkpoint"] else None
    if check(requested is not None and requested.exists(),
             f"requested checkpoint not reported/written: {status['checkpoint']}"):
        chk = read_checkpoint(requested)
        check(requested.name == f"state_step{int(chk['step']):012d}.gmdchk",
              f"step-suffixed name wrong: {requested.name}")
        check(int(chk["step"]) % 5 == 0, "requested checkpoint not at a poll boundary")
        print(f"  checkpoint request served at step {chk['step']} -> {requested.name}; "
              "request file removed")
    process.send_signal(signal.SIGTERM)
    rc, stdout, _ = finish(process)
    check_stopped_run(directory, rc, stdout, "SIGTERM", 5)
    final = Path(read_status(directory / "status.json")["checkpoint"])
    check(requested is not None and final != requested and final.exists() and requested.exists(),
          "SIGTERM stop did not write a separate step-suffixed checkpoint")

    directory, process = control_case(gmd, root, "control_sigint", 3, [])
    process.send_signal(signal.SIGINT)
    rc, stdout, _ = finish(process)
    check_stopped_run(directory, rc, stdout, "SIGINT", 3)


# --------------------------------------------------------------------------
# ML tests
# --------------------------------------------------------------------------

ML_BOX = 12.0
# Two water molecules. The second straddles the x face of the box, with its O
# 0.001 A inside it, so the O crosses the face early in the run. Atoms are wrapped into the box every step while the
# Verlet list's image flags date from the last rebuild, so this exercises the
# shift correction for atoms that wrapped since the list was built (without it
# the frame forces stop matching the Python evaluation below).
ML_ATOMS = [
    ("O", 3.000, 3.000, 3.000), ("H", 3.957, 3.000, 3.000), ("H", 2.760, 3.927, 3.000),
    ("O", 11.999, 4.400, 3.400), ("H", 0.650, 5.000, 3.100), ("H", 11.450, 3.600, 2.950),
]
Z = {"H": 1, "O": 8}


def export_models(sgt_root: Path, out_dir: Path) -> list[Path]:
    sys.path.insert(0, str(sgt_root))
    import torch  # noqa: PLC0415
    from gmd_sgt.inference.export import script_model  # noqa: PLC0415
    from gmd_sgt.models import AllegroStyleBackbone  # noqa: PLC0415

    paths = []
    for name, seed in (("modelA.pt", 1), ("modelB.pt", 2)):
        torch.manual_seed(seed)
        model = AllegroStyleBackbone(n_species=10, hidden_channels=16, num_layers=2, n_basis=4,
                                     cutoff=3.0, l_max=2, avg_neighbors=4.0)
        for parameter in model.readout.net.parameters():
            torch.nn.init.normal_(parameter, std=0.3)
        model.eval()
        path = out_dir / name
        script_model(model).save(str(path))
        paths.append(path)
    return paths


def ml_inputs(directory: Path, steps: int, extra: list[str], model: Path,
              atoms=ML_ATOMS) -> None:
    lines = [str(len(atoms)), f"{ML_BOX} {ML_BOX} {ML_BOX}"]
    for symbol, x, y, z, *mass in atoms:
        lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}" + "".join(f" {m}" for m in mass))
    (directory / "in.xyz").write_text("\n".join(lines) + "\n")
    (directory / "run.in").write_text("\n".join([
        f"run {steps}", "time_step 0.25", "velocity 200.0", "velocity_init random",
        "velocity_seed 11", "output_interval 5", "force_field ml", f"model_path {model}",
        *extra]) + "\n")


def eval_model(model_path: Path, symbols: list[str], positions) -> tuple[float, list]:
    import torch  # noqa: PLC0415
    model = torch.jit.load(str(model_path))
    cutoff = float(model.local_cutoff)
    n = len(symbols)
    src, dst, shifts = [], [], []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = [positions[j][k] - positions[i][k] for k in range(3)]
            s = [-ML_BOX * round(d[k] / ML_BOX) for k in range(3)]
            if sum((d[k] + s[k]) ** 2 for k in range(3)) < cutoff * cutoff:
                src.append(i)
                dst.append(j)
                shifts.append(s)
    out = model(torch.tensor([Z[s] for s in symbols], dtype=torch.long),
                torch.tensor(positions, dtype=torch.float32),
                torch.tensor([src, dst], dtype=torch.long),
                torch.tensor(shifts, dtype=torch.float32).reshape(-1, 3))
    return float(out["energy"]), out["forces"].tolist()


def test_ml(gmd: Path, root: Path, sgt_root: Path, steps: int, require_ase: bool) -> None:
    print("[ML] TorchScript (GMD_SGT export)")
    models_dir = root / "models"
    models_dir.mkdir()
    model_a, model_b = export_models(sgt_root, models_dir)
    sha_a, sha_b = sha256_file(model_a), sha256_file(model_b)

    # ml_jit_optimize off: every model evaluation is a pure function of its
    # inputs, so the segmented run must reproduce the continuous one BITWISE.
    def inputs(directory: Path, n: int, extra: list[str]) -> None:
        ml_inputs(directory, n, extra + ["ml_jit_optimize off"], model_a)

    dirs = segment_test(gmd, root, inputs, steps, "ml", 0.0)

    # Default (ml_jit_optimize on): TorchScript's profiling executor switches
    # to a fused graph after its first calls, so the restarted run's first
    # evaluation (call 1) and the continuous run's (call N+1) differ by float32
    # round-off: one ulp of a float32 force is ~6e-8 relative, i.e. <= ~1e-7
    # eV/A for forces of order 1 eV/A. Over N short steps that perturbs
    # positions far less. 1e-6 bounds the force/energy round-off with margin
    # and is still ~5 orders below any physical effect.
    def inputs_default(directory: Path, n: int, extra: list[str]) -> None:
        ml_inputs(directory, n, extra, model_a)

    segment_test(gmd, root, inputs_default, steps, "ml_jit_default", 1e-6)
    half = read_checkpoint(dirs["seg1"] / "half.gmdchk")
    check(half["force_field_kind"] == "ml" and half["model_sha256"] == sha_a and
          half["force_field_sha256"] == sha_a and Path(half["model_path"]) == model_a,
          f"ML checkpoint provenance wrong: {half.get('model_sha256')}")
    check([half["atoms"][t]["atomic_number"] for t in sorted(half["atoms"])] ==
          [str(Z[a[0]]) for a in ML_ATOMS], "ML checkpoint atomic numbers wrong")

    # Frame forces vs a Python re-evaluation of the exported model.
    frames = parse_frames(dirs["cont"] / "frames.extxyz")
    worst_f, worst_e = 0.0, 0.0
    for frame in frames:
        check(frame["info"]["model_sha256"] == sha_a, "frame model_sha256 is not model A's")
        symbols = [a[1] for a in frame["atoms"]]
        positions = [[float(v) for v in a[2:5]] for a in frame["atoms"]]
        energy, forces = eval_model(model_a, symbols, positions)
        worst_e = max(worst_e, abs(energy - float(frame["info"]["potential_energy"])))
        for row, ref in zip(frame["atoms"], forces):
            for k in range(3):
                worst_f = max(worst_f, abs(float(row[5 + k]) - ref[k]))
    # GMD passes float32 positions and edge shifts and the model returns float32
    # forces; GMD's own shift vectors (S*L from its neighbour list) and Python's
    # are equal in exact arithmetic, so the residual is float32 round-off.
    check(worst_f < 1e-5 and worst_e < 1e-5,
          f"ML frame forces/energy differ from Python evaluation: {worst_f:.3e}, {worst_e:.3e}")
    print(f"  frames vs Python model evaluation: max |dF| = {worst_f:.3e} eV/A, "
          f"max |dE| = {worst_e:.3e} eV")
    check_frames_with_ase(dirs["cont"] / "frames.extxyz", frames, require_ase)

    checkpoint = dirs["seg1"] / "half.gmdchk"
    refused = root / "ml_swap_refused"
    refused.mkdir()
    ml_inputs(refused, steps, [f"restart_from {checkpoint}", "status_file status.json"], model_b)
    result = run_gmd(gmd, refused / "in.xyz", refused / "run.in", expect_ok=False)
    check(result.returncode != 0 and "model_change allow" in result.stderr,
          f"ML restart with model B was not refused: rc={result.returncode} {result.stderr}")
    check(read_status(refused / "status.json")["state"] == "failed", "refused swap status")

    allowed = root / "ml_swap_allowed"
    allowed.mkdir()
    ml_inputs(allowed, steps, [f"restart_from {checkpoint}", "model_change allow",
                               "frame_file frames.extxyz", "status_file status.json",
                               "checkpoint_file final.gmdchk"], model_b)
    result = run_gmd(gmd, allowed / "in.xyz", allowed / "run.in")
    check(sha_a in result.stdout and sha_b in result.stdout,
          "model swap not reported with old and new sha on stdout")
    comment = [l for l in (allowed / "output.log").read_text().splitlines()
               if l.startswith("# model_change")]
    check(len(comment) == 1 and sha_a in comment[0] and sha_b in comment[0],
          f"output.log model_change line wrong: {comment}")
    frames_b = parse_frames(allowed / "frames.extxyz")
    check(all(f["info"]["model_sha256"] == sha_b and f["info"].get("model_changed_from") == sha_a
              for f in frames_b), "swap frames do not record model B and the change from A")
    status = read_status(allowed / "status.json")
    check(status["model_sha256"] == sha_b and status["model_changed_from"] == sha_a,
          f"swap status wrong: {status}")
    final_b = read_checkpoint(allowed / "final.gmdchk")
    check(final_b["model_sha256"] == sha_b, "post-swap checkpoint does not record model B")
    # The restarted state is evaluated with the correct species: the first
    # frame's energy is model B at the checkpoint geometry with real Z.
    first = frames_b[0]
    energy_b, _ = eval_model(model_b, [a[1] for a in first["atoms"]],
                             [[float(v) for v in a[2:5]] for a in first["atoms"]])
    check(abs(energy_b - float(first["info"]["potential_energy"])) < 1e-5,
          "restarted ML energy is not model B with the checkpointed species")
    print(f"  swap A->B refused without directive; accepted and recorded with it "
          f"({sha_a[:12]} -> {sha_b[:12]})")

    # A pre-v3 checkpoint (no atomic numbers) must be refused for ML.
    legacy = root / "ml_legacy"
    legacy.mkdir()
    lines = checkpoint.read_text().splitlines()
    out = []
    for line in lines:
        if line.startswith("GMD_CHECKPOINT"):
            out.append("GMD_CHECKPOINT 2")
        elif line.split(" ")[0] in ("force_field_kind", "force_field_sha256", "model_path",
                                    "model_sha256", "run_id"):
            continue
        elif line.startswith("atoms tag"):
            out.append("atoms tag type molecule mass charge x y z vx vy vz")
        elif len(line.split()) == 12 and line.split()[0].isdigit():
            out.append(" ".join(line.split()[:11]))
        else:
            out.append(line)
    (legacy / "v2.gmdchk").write_text("\n".join(out) + "\n")
    ml_inputs(legacy, steps, [f"restart_from {legacy / 'v2.gmdchk'}", "model_change allow"],
              model_a)
    result = run_gmd(gmd, legacy / "in.xyz", legacy / "run.in", expect_ok=False)
    check(result.returncode != 0 and "does not store atomic numbers" in result.stderr,
          f"ML restart from a v2 checkpoint was not refused: {result.stderr.strip()}")
    print("  v2 checkpoint refused for ML (no atomic numbers)")

    # Unknown species: refused before any step.
    unknown = root / "ml_unknown_species"
    unknown.mkdir()
    atoms = [("Xx", 3.0, 3.0, 3.0, 20.0)] + [(*a, 1.008 if a[0] == "H" else 15.999)
                                             for a in ML_ATOMS[1:]]
    ml_inputs(unknown, steps, ["status_file status.json"], model_a, atoms=atoms)
    result = run_gmd(gmd, unknown / "in.xyz", unknown / "run.in", expect_ok=False)
    check(result.returncode != 0 and "unknown species" in result.stderr,
          f"atom with unknown species was not refused: {result.stderr.strip()}")
    status = read_status(unknown / "status.json")
    check(status["state"] == "failed" and "unknown species" in status["message"],
          f"unknown-species status wrong: {status}")
    check(not (unknown / "output.log").exists() or len(read_log(unknown / "output.log")) == 0,
          "a frame was written despite unknown species")
    print("  atom with Z = 0 refused before the first step")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gmd", required=True, type=Path)
    parser.add_argument("--work-root", required=True, type=Path)
    parser.add_argument("--steps", type=int, default=40, help="N for 2N vs N+N")
    parser.add_argument("--require-ase", action="store_true")
    parser.add_argument("--ml", action="store_true", help="run the ML part")
    parser.add_argument("--sgt-root", type=Path, help="GMD_SGT checkout (for --ml)")
    parser.add_argument("--skip-classical", action="store_true")
    args = parser.parse_args()

    root = args.work_root.resolve()
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    gmd = args.gmd.resolve()

    if not args.skip_classical:
        dirs = test_lj_segments(gmd, root, args.steps)
        test_lj_model_change(gmd, root, dirs, args.steps)
        test_lj_frames(gmd, dirs, args.require_ase)
        test_relative_paths(gmd, root)
        test_control(gmd, root)
    if args.ml:
        if args.sgt_root is None:
            parser.error("--ml requires --sgt-root")
        test_ml(gmd, root, args.sgt_root.resolve(), args.steps, args.require_ase)

    if problems:
        print(f"[controller integration] {len(problems)} problem(s)", file=sys.stderr)
        for message in problems:
            print(f"  {message}", file=sys.stderr)
        return 1
    print("[controller integration] all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
