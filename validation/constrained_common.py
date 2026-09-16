"""Shared helpers for the constrained-dynamics validation cases.

WHY THIS EXISTS SEPARATELY FROM common.py
-----------------------------------------
The unconstrained cases read `output.log`, whose columns are written with six
fixed decimals. That is fine for energies of order 1e-2 eV, but it is useless for
a constraint residual: a bond error of 1e-13 A prints as "0.000000", and so does
a bond error of 1e-7 A. Validating that constraints are satisfied to 1e-10 CANNOT
be done from the log.

The CHECKPOINT can do it. `gmd` writes checkpoints with 17 significant digits and
they carry exactly what is needed:

  * every atom's global tag, mass, position and velocity, so the position
    residual |r_ij| - d_ij and the velocity tangency residual |r_ij . v_ij| /
    |r_ij| can be recomputed here, independently of whatever the solver believed;
  * the constraint virial, its validity state and its TIME LEVEL;
  * the provider virial and its validity;
  * the completed-step pressure together with the 2K, volume and virial it was
    formed from, so the kinetic-plus-virial identity can be checked rather than
    taken on trust.

Recomputing the residuals here from the checkpointed state is the point: it is an
independent measurement of what the production path actually produced, not a
read-back of the solver's own convergence flag.

NaN POLICY
----------
Every numeric accessor below refuses non-finite values at the point of parsing,
BEFORE any error metric is formed. This is deliberate: `abs(nan - x) <= tol` is
False in Python, so a NaN would often fail a comparison anyway -- but `nan` fed
into a mean, a max or a drift denominator can silently produce a passing number.
The one place a NaN is legitimate is the pressure column of a frame whose
P_valid flag is 0, which is how `gmd` spells "no complete pressure exists"; that
case is handled explicitly and never reaches an error metric.
"""
from __future__ import annotations

import math
import pathlib
from typing import Iterable

# GMD stores velocities in its INTERNAL time unit, not in femtoseconds. With
# energies in eV, masses in amu and distances in Angstrom the time unit is fixed
# at A sqrt(amu/eV); see include/gmd/core/physical_constants.hpp. Every velocity
# read out of a checkpoint is therefore in A per internal time unit, and so is
# the velocity-tangency residual |r.v|/|r| that the solver converges.
#
# NOTE: the `output.log` header labels that column "rattle_error[A/fs]", which
# is wrong by this factor. The number itself is the internal-unit one.
FEMTOSECONDS_PER_INTERNAL_TIME = 10.180505717871194
INTERNAL_TIME_PER_FEMTOSECOND = 1.0 / FEMTOSECONDS_PER_INTERNAL_TIME

# The constrained log has five more columns than the unconstrained one.
LOG_COLUMNS = [
    "step", "time_fs", "pe", "ke", "total_energy", "temperature",
    "pressure_bar", "volume_a3",
    "shake_iterations", "shake_error_a", "rattle_iterations",
    "rattle_error_a_per_fs", "pressure_valid",
]
INTEGER_COLUMNS = {"step", "shake_iterations", "rattle_iterations", "pressure_valid"}


class ValidationFailure(Exception):
    """Raised with a message that names exactly what failed and where."""


def _finite(value: float, what: str) -> float:
    if not math.isfinite(value):
        raise ValidationFailure(f"{what} is not finite (got {value!r})")
    return value


# ---------------------------------------------------------------------------
# Constrained energy log
# ---------------------------------------------------------------------------
def parse_constrained_log(path: pathlib.Path) -> list[dict]:
    """Parse the 13-column constrained `output.log`.

    A frame whose pressure_valid flag is 0 keeps `pressure_bar = None` rather
    than NaN, so that a caller cannot accidentally average it.
    """
    if not path.is_file():
        raise ValidationFailure(f"missing log file {path}")

    rows: list[dict] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        if len(tokens) != len(LOG_COLUMNS):
            raise ValidationFailure(
                f"{path}:{line_number}: expected {len(LOG_COLUMNS)} columns "
                f"({', '.join(LOG_COLUMNS)}), got {len(tokens)}: {line!r}"
            )
        row: dict = {}
        for column, token in zip(LOG_COLUMNS, tokens):
            if column in INTEGER_COLUMNS:
                row[column] = int(token)
                continue
            value = float(token)
            if column == "pressure_bar":
                row[column] = value          # validated below against the flag
                continue
            row[column] = _finite(value, f"{path}:{line_number}: column {column}")

        valid = row["pressure_valid"]
        if valid not in (0, 1):
            raise ValidationFailure(
                f"{path}:{line_number}: pressure_valid must be 0 or 1, got {valid}")
        if valid == 1:
            row["pressure_bar"] = _finite(
                row["pressure_bar"], f"{path}:{line_number}: column pressure_bar "
                                     f"(step {row['step']}, pressure_valid=1)")
        else:
            # `gmd` documents an unavailable pressure as nan. Anything else --
            # a zero in particular -- would be a real pressure masquerading as a
            # sentinel, which is the failure the P_valid column exists to stop.
            if not math.isnan(row["pressure_bar"]):
                raise ValidationFailure(
                    f"{path}:{line_number}: step {row['step']} has "
                    f"pressure_valid=0 but pressure_bar={row['pressure_bar']!r}; "
                    f"an unavailable pressure must be written as nan, never as a "
                    f"number (a zero is a pressure, not a sentinel)")
            row["pressure_bar"] = None
        rows.append(row)

    if not rows:
        raise ValidationFailure(f"no data rows in {path}")
    return rows


# ---------------------------------------------------------------------------
# Checkpoint
# ---------------------------------------------------------------------------
def parse_checkpoint(path: pathlib.Path) -> dict:
    """Parse a `.gmdchk` checkpoint into a dict.

    Returns keys: step, time_fs, box, atoms (list of dicts with tag/mass/
    position/velocity), constraint_virial{state,time_level,tensor},
    provider_virial{valid,tensor}, step_pressure{...}, thermostat{type,state}.
    """
    if not path.is_file():
        raise ValidationFailure(f"missing checkpoint {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("GMD_CHECKPOINT"):
        raise ValidationFailure(f"{path}: not a GMD checkpoint")

    def find(key: str) -> list[str]:
        for line in lines:
            if line.startswith(key + " ") or line == key:
                return line.split()
        raise ValidationFailure(f"{path}: missing required key {key!r}")

    result: dict = {}
    result["step"] = int(find("step")[1])
    result["time_fs"] = _finite(float(find("time_fs")[1]), f"{path}: time_fs")
    result["box"] = [_finite(float(v), f"{path}: box") for v in find("box")[1:4]]

    tokens = find("constraint_virial")
    state, time_level = tokens[1], tokens[2]
    if state not in ("not_applicable", "unavailable", "valid"):
        raise ValidationFailure(f"{path}: unknown constraint_virial state {state!r}")
    result["constraint_virial"] = {
        "state": state,
        "time_level": time_level,
        "tensor": [_finite(float(v), f"{path}: constraint_virial component {i}")
                   for i, v in enumerate(tokens[3:12])],
    }

    tokens = find("provider_virial")
    result["provider_virial"] = {
        "valid": tokens[1] == "1",
        "tensor": [_finite(float(v), f"{path}: provider_virial component {i}")
                   for i, v in enumerate(tokens[2:11])],
    }

    tokens = find("step_pressure")
    valid = tokens[1] == "1"
    step_pressure = {"valid": valid}
    if valid:
        names = ["pressure", "twice_kinetic_energy", "volume", "potential_energy"]
        for offset, name in enumerate(names):
            step_pressure[name] = _finite(
                float(tokens[2 + offset]), f"{path}: step_pressure {name}")
        step_pressure["tensor"] = [
            _finite(float(v), f"{path}: step_pressure virial component {i}")
            for i, v in enumerate(tokens[6:15])]
    result["step_pressure"] = step_pressure

    thermostat_type = find("thermostat_type")
    result["thermostat"] = {
        "type": thermostat_type[1] if len(thermostat_type) > 1 else "",
        "state": " ".join(find("thermostat_state")[1:]),
    }

    count = int(find("atom_count")[1])
    start = next(i for i, line in enumerate(lines) if line.startswith("atoms ")) + 1
    atoms = []
    for offset in range(count):
        tokens = lines[start + offset].split()
        if len(tokens) < 11:
            raise ValidationFailure(
                f"{path}: atom row {offset} has {len(tokens)} fields, expected >= 11")
        tag = int(tokens[0])
        atoms.append({
            "tag": tag,
            "mass": _finite(float(tokens[3]), f"{path}: atom tag {tag} mass"),
            "position": [_finite(float(tokens[5 + d]), f"{path}: atom tag {tag} position")
                         for d in range(3)],
            "velocity": [_finite(float(tokens[8 + d]), f"{path}: atom tag {tag} velocity")
                         for d in range(3)],
        })
    if len(atoms) != count:
        raise ValidationFailure(f"{path}: expected {count} atoms, parsed {len(atoms)}")
    result["atoms"] = atoms
    return result


# ---------------------------------------------------------------------------
# Constraint residuals, recomputed from checkpointed state
# ---------------------------------------------------------------------------
def minimum_image(delta: list[float], box: list[float]) -> list[float]:
    out = []
    for value, length in zip(delta, box):
        while value > 0.5 * length:
            value -= length
        while value < -0.5 * length:
            value += length
        out.append(value)
    return out


def constraint_residuals(state: dict, constraints: list[tuple[int, int, float]]) -> dict:
    """Recompute both constraint residuals from a checkpointed state.

    `constraints` are (tag_i, tag_j, target_distance) with GLOBAL atom tags, so
    the check does not depend on atom ordering or on MPI decomposition.

    Returns the worst of each residual together with the constraint and the atom
    tags that produced it, so a failure can name them.

      position residual  | |r_ij| - d_ij |            [A]
      velocity residual  | r_ij . v_ij | / |r_ij|     [A / internal time unit]

    The velocity residual is the per-constraint row of |J v| normalised by the
    bond length, which is the same quantity the solver converges and the same one
    the log calls rattle_error[A/fs].
    """
    by_tag = {atom["tag"]: atom for atom in state["atoms"]}
    box = state["box"]

    worst_position = {"value": 0.0, "constraint": None}
    worst_velocity = {"value": 0.0, "constraint": None}

    for index, (tag_i, tag_j, target) in enumerate(constraints):
        for tag in (tag_i, tag_j):
            if tag not in by_tag:
                raise ValidationFailure(
                    f"constraint {index} references atom tag {tag}, which is not "
                    f"present in the checkpoint (tags present: "
                    f"{sorted(by_tag)[:16]}...)")
        atom_i, atom_j = by_tag[tag_i], by_tag[tag_j]

        dr = minimum_image(
            [atom_i["position"][d] - atom_j["position"][d] for d in range(3)], box)
        distance = math.sqrt(sum(value * value for value in dr))
        if not math.isfinite(distance) or distance <= 0.0:
            raise ValidationFailure(
                f"constraint {index} on tags ({tag_i}, {tag_j}) has a degenerate "
                f"bond vector {dr!r}")
        position_error = abs(distance - target)

        dv = [atom_i["velocity"][d] - atom_j["velocity"][d] for d in range(3)]
        velocity_error = abs(sum(dr[d] * dv[d] for d in range(3))) / distance

        detail = {"index": index, "tags": (tag_i, tag_j), "target": target,
                  "distance": distance}
        if position_error > worst_position["value"]:
            worst_position = {"value": position_error, "constraint": detail}
        if velocity_error > worst_velocity["value"]:
            worst_velocity = {"value": velocity_error, "constraint": detail}

    return {"max_position_residual_a": worst_position,
            "max_velocity_residual_internal": worst_velocity}


def describe_constraint(detail: dict | None) -> str:
    if detail is None:
        return "(no constraint)"
    tag_i, tag_j = detail["tags"]
    return (f"constraint {detail['index']} on atom tags ({tag_i}, {tag_j}), "
            f"target {detail['target']:.14f} A, measured {detail['distance']:.14f} A")


# ---------------------------------------------------------------------------
# Momentum, angular momentum, and the pressure identity
# ---------------------------------------------------------------------------
def unwrapped_positions(state: dict, tags: list[int]) -> dict[int, list[float]]:
    """Positions of `tags`, unwrapped into one continuous image.

    Atoms are wrapped into the cell INDEPENDENTLY, so a molecule straddling a
    periodic face has its atoms in different images and any centre of mass
    computed from the raw coordinates is meaningless. Each atom is therefore
    placed in the image nearest the first tag, which reconstructs the molecule
    as a single connected object. This is well defined for a group whose extent
    is under half the box, which every molecule here is.
    """
    by_tag = {atom["tag"]: atom for atom in state["atoms"]}
    box = state["box"]
    anchor = by_tag[tags[0]]["position"]
    out = {}
    for tag in tags:
        delta = minimum_image(
            [by_tag[tag]["position"][d] - anchor[d] for d in range(3)], box)
        out[tag] = [anchor[d] + delta[d] for d in range(3)]
    return out


def group_center_of_mass(state: dict, tags: list[int]) -> list[float]:
    """Centre of mass of a molecule, computed from unwrapped positions."""
    by_tag = {atom["tag"]: atom for atom in state["atoms"]}
    positions = unwrapped_positions(state, tags)
    total_mass = sum(by_tag[tag]["mass"] for tag in tags)
    return [sum(by_tag[tag]["mass"] * positions[tag][d] for tag in tags) / total_mass
            for d in range(3)]


def separation_modulo_box(a: list[float], b: list[float],
                          box: list[float]) -> float:
    """Largest component of a - b, taken under the minimum image convention.

    Used to compare a predicted position against an observed one when the
    observed one may have been wrapped an arbitrary number of times.
    """
    return max(abs(v) for v in minimum_image([a[d] - b[d] for d in range(3)], box))


def total_momentum(state: dict) -> list[float]:
    return [sum(a["mass"] * a["velocity"][d] for a in state["atoms"]) for d in range(3)]


def center_of_mass_velocity(state: dict) -> list[float]:
    total_mass = sum(a["mass"] for a in state["atoms"])
    return [p / total_mass for p in total_momentum(state)]


def angular_momentum(state: dict) -> list[float]:
    """L = sum_i m_i (r_i - r_com) x (v_i - v_com), in the COM frame."""
    total_mass = sum(a["mass"] for a in state["atoms"])
    com = [sum(a["mass"] * a["position"][d] for a in state["atoms"]) / total_mass
           for d in range(3)]
    vcom = center_of_mass_velocity(state)
    out = [0.0, 0.0, 0.0]
    for atom in state["atoms"]:
        r = [atom["position"][d] - com[d] for d in range(3)]
        v = [atom["velocity"][d] - vcom[d] for d in range(3)]
        out[0] += atom["mass"] * (r[1] * v[2] - r[2] * v[1])
        out[1] += atom["mass"] * (r[2] * v[0] - r[0] * v[2])
        out[2] += atom["mass"] * (r[0] * v[1] - r[1] * v[0])
    return out


def twice_kinetic_energy(state: dict) -> float:
    return sum(a["mass"] * sum(v * v for v in a["velocity"]) for a in state["atoms"])


def trace(tensor: list[float]) -> float:
    return tensor[0] + tensor[4] + tensor[8]


def max_asymmetry(tensor: list[float]) -> float:
    return max(abs(tensor[1] - tensor[3]),
               abs(tensor[2] - tensor[6]),
               abs(tensor[5] - tensor[7]))


TENSOR_COMPONENT_NAMES = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]


def check_tensor_finite(tensor: list[float], label: str) -> None:
    if len(tensor) != 9:
        raise ValidationFailure(f"{label}: expected 9 components, got {len(tensor)}")
    for index, value in enumerate(tensor):
        if not math.isfinite(value):
            raise ValidationFailure(
                f"{label}: component {TENSOR_COMPONENT_NAMES[index]} "
                f"(index {index}) is not finite (got {value!r})")


def compare_tensors(actual: list[float], expected: list[float], tolerance: float,
                    label: str) -> dict:
    """Component-by-component comparison that names the component that failed."""
    check_tensor_finite(actual, f"{label} (actual)")
    check_tensor_finite(expected, f"{label} (expected)")
    worst_index, worst = 0, 0.0
    for index in range(9):
        difference = abs(actual[index] - expected[index])
        if difference > worst:
            worst, worst_index = difference, index
    return {
        "max_component_difference": worst,
        "worst_component": TENSOR_COMPONENT_NAMES[worst_index],
        "tolerance": tolerance,
        "passed": worst <= tolerance,
    }


def mean(values: Iterable[float]) -> float:
    values = list(values)
    return sum(values) / len(values) if values else 0.0
