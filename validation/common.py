from __future__ import annotations

import csv
import json
import math
import pathlib
import subprocess
import sys
from typing import Iterable


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def load_json(path: pathlib.Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: pathlib.Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def run_command(args: list[str], cwd: pathlib.Path) -> None:
    completed = subprocess.run(args, cwd=cwd, text=True, capture_output=True)
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise RuntimeError(f"Command failed: {' '.join(args)}")


def compare_scalar(actual: float, expected: float, abs_tol: float, rel_tol: float = 0.0) -> dict:
    diff = actual - expected
    limit = max(abs_tol, rel_tol * abs(expected))
    return {
        "actual": actual,
        "expected": expected,
        "abs_error": abs(diff),
        "passed": abs(diff) <= limit,
        "tolerance": limit,
    }


def force_error_metrics(actual: list[list[float]], expected: list[list[float]]) -> dict:
    if len(actual) != len(expected):
        raise RuntimeError("Force array length mismatch")

    sum_sq = 0.0
    max_abs = 0.0
    max_norm = 0.0
    count = 0
    for actual_force, expected_force in zip(actual, expected):
        for actual_component, expected_component in zip(actual_force, expected_force):
            delta = actual_component - expected_component
            sum_sq += delta * delta
            max_abs = max(max_abs, abs(delta))
            count += 1
        dx = actual_force[0] - expected_force[0]
        dy = actual_force[1] - expected_force[1]
        dz = actual_force[2] - expected_force[2]
        max_norm = max(max_norm, math.sqrt(dx * dx + dy * dy + dz * dz))

    rms = math.sqrt(sum_sq / count) if count > 0 else 0.0
    return {
        "rms_error": rms,
        "max_component_error": max_abs,
        "max_force_norm_error": max_norm,
    }


def parse_energy_log(path: pathlib.Path) -> list[dict]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            rows.append(
                {
                    "step": int(tokens[0]),
                    "time_fs": float(tokens[1]),
                    "pe": float(tokens[2]),
                    "ke": float(tokens[3]),
                    "total_energy": float(tokens[4]),
                    "temperature": float(tokens[5]),
                    "pressure_bar": float(tokens[6]) if len(tokens) > 6 else 0.0,
                    "volume_a3": float(tokens[7]) if len(tokens) > 7 else 0.0,
                }
            )
    if not rows:
        raise RuntimeError(f"No data rows found in {path}")
    return rows


def parse_xyz_frames(path: pathlib.Path) -> list[dict]:
    frames = []
    with path.open("r", encoding="utf-8") as handle:
        while True:
            count_line = handle.readline()
            if not count_line:
                break
            atom_count = int(count_line.strip())
            comment = handle.readline().strip()
            atoms = []
            for _ in range(atom_count):
                tokens = handle.readline().split()
                atoms.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
            frames.append({"comment": comment, "coordinates": atoms})
    if not frames:
        raise RuntimeError(f"No XYZ frames found in {path}")
    return frames


def compute_msd_series(frames: list[dict], box_lengths: list[float], time_fs: list[float]) -> list[dict]:
    if len(frames) != len(time_fs):
        raise RuntimeError("Frame/time length mismatch")
    if not frames:
        return []

    reference = frames[0]["coordinates"]
    previous = [coords[:] for coords in reference]
    unwrapped = [coords[:] for coords in reference]
    series = []

    for frame_index, frame in enumerate(frames):
        current = frame["coordinates"]
        if frame_index > 0:
            for atom_index in range(len(current)):
                for dim in range(3):
                    delta = current[atom_index][dim] - previous[atom_index][dim]
                    length = box_lengths[dim]
                    if delta > 0.5 * length:
                        delta -= length
                    elif delta < -0.5 * length:
                        delta += length
                    unwrapped[atom_index][dim] += delta
            previous = [coords[:] for coords in current]

        msd = 0.0
        for atom_index in range(len(unwrapped)):
            dx = unwrapped[atom_index][0] - reference[atom_index][0]
            dy = unwrapped[atom_index][1] - reference[atom_index][1]
            dz = unwrapped[atom_index][2] - reference[atom_index][2]
            msd += dx * dx + dy * dy + dz * dz
        msd /= len(unwrapped)
        series.append({"time_fs": time_fs[frame_index], "msd_a2": msd})
    return series


def fit_diffusion_coefficient(msd_series: list[dict], fit_fraction: float = 0.5) -> float:
    if len(msd_series) < 2:
        return 0.0
    start_index = max(0, int(len(msd_series) * (1.0 - fit_fraction)))
    subset = msd_series[start_index:]
    xs = [row["time_fs"] * 1.0e-3 for row in subset]
    ys = [row["msd_a2"] for row in subset]
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    denominator = sum((x - x_mean) * (x - x_mean) for x in xs)
    slope = numerator / denominator if denominator > 0.0 else 0.0
    return slope / 6.0


def write_series_csv(path: pathlib.Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def stddev(values: list[float]) -> float:
    if not values:
        return 0.0
    avg = mean(values)
    return math.sqrt(sum((value - avg) * (value - avg) for value in values) / len(values))
