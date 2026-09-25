#!/usr/bin/env python3
"""Small FP64 checks for the revised graph--tensor manuscript.

This is a new reference implementation, not the missing script behind the old
draft. It imports no GMD engine code. NumPy operator expressions are compared
with separately written scalar/indexed formulas and finite differences of the
scalar energies. Shared model definitions are not an external-code validation.
Run from any directory; default inputs and outputs are next to this file.
"""

import argparse
import hashlib
import itertools
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def incidence(indices, n):
    """Tiny explicit matrices used only for mathematical verification."""
    indices = np.asarray(indices, dtype=np.int64)
    source = np.zeros((len(indices), n))
    target = np.zeros_like(source)
    source[np.arange(len(indices)), indices[:, 0]] = 1.0
    target[np.arange(len(indices)), indices[:, 1]] = 1.0
    return target - source, source, target


def assemble(indices, local, n):
    out = np.zeros((n, 3))
    np.add.at(out, np.asarray(indices).ravel(), local.reshape(-1, 3))
    return out


def lj_values(r, epsilon, sigma, cutoff=None):
    """Array energy and radial derivative, optionally force shifted."""
    a = (sigma / r) ** 6
    u = 4.0 * epsilon * (a * a - a)
    du = 24.0 * epsilon / r * (a - 2.0 * a * a)
    if cutoff is not None:
        ac = (sigma / cutoff) ** 6
        uc = 4.0 * epsilon * (ac * ac - ac)
        duc = 24.0 * epsilon / cutoff * (ac - 2.0 * ac * ac)
        u = np.where(r < cutoff, u - uc - (r - cutoff) * duc, 0.0)
        du = np.where(r < cutoff, du - duc, 0.0)
    return u, du


def lj_operator(x, pairs, parameters, box=None, scales=None, cutoff=None):
    pairs = np.asarray(pairs, dtype=np.int64)
    weights = np.ones(len(pairs)) if scales is None else np.asarray(scales)
    # Drop excluded pairs BEFORE evaluating their distance or inverse powers.
    active = weights != 0.0
    pairs, weights = pairs[active], weights[active]
    b, _, _ = incidence(pairs, len(x))
    d = b @ x
    if box is not None:
        d = d - np.rint(d / box) * box
    r = np.linalg.norm(d, axis=1)
    if np.any(r <= 0.0):
        raise ValueError("Coincident particles in an active LJ pair")
    u, du = lj_values(r, parameters["epsilon"], parameters["sigma"], cutoff)
    g = (weights * du / r)[:, None] * d
    return float(np.sum(weights * u)), -b.T @ g, -d.T @ g


def bonded_operator(x, parameters, kind):
    indices = np.asarray(parameters["indices"], dtype=np.int64)
    q = x[indices]
    if kind == "bonds":
        d = q[:, 1] - q[:, 0]
        r = np.linalg.norm(d, axis=1)
        k = np.asarray(parameters["k"])
        dr = r - parameters["r0"]
        g = (k * dr / r)[:, None] * d
        return float(np.sum(0.5 * k * dr ** 2)), assemble(indices, np.stack((g, -g), axis=1), len(x))
    if kind == "angles":
        a, b = q[:, 0] - q[:, 1], q[:, 2] - q[:, 1]
        alpha, beta = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
        c = np.einsum("ij,ij->i", a, b) / (alpha * beta)
        theta = np.arccos(c)
        k = np.asarray(parameters["k"])
        delta = theta - parameters["theta0"]
        prefactor = -k * delta / np.sqrt(1.0 - c * c)
        ga = prefactor[:, None] * (b / (alpha * beta)[:, None] - (c / alpha ** 2)[:, None] * a)
        gb = prefactor[:, None] * (a / (alpha * beta)[:, None] - (c / beta ** 2)[:, None] * b)
        return float(np.sum(0.5 * k * delta ** 2)), assemble(indices, np.stack((-ga, ga + gb, -gb), axis=1), len(x))
    if kind != "dihedrals":
        raise ValueError(kind)
    b1, b2, b3 = q[:, 1] - q[:, 0], q[:, 2] - q[:, 1], q[:, 3] - q[:, 2]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    ell = np.linalg.norm(b2, axis=1)
    tangent = b2 / ell[:, None]
    cosine = np.einsum("ij,ij->i", n1, n2)
    sine = np.einsum("ij,ij->i", tangent, np.cross(n1, n2))
    chi = np.arctan2(sine, cosine)
    kappa, mult, phase = (np.asarray(parameters[key]) for key in ("kappa", "multiplicity", "phase"))
    u = kappa * (1.0 + np.cos(mult * chi - phase))
    duchi = -kappa * mult * np.sin(mult * chi - phase)
    xb = -duchi * sine / (cosine ** 2 + sine ** 2)
    yb = duchi * cosine / (cosine ** 2 + sine ** 2)
    n1b = xb[:, None] * n2 + yb[:, None] * np.cross(n2, tangent)
    n2b = xb[:, None] * n1 + yb[:, None] * np.cross(tangent, n1)
    tb = yb[:, None] * np.cross(n1, n2)
    b1b = np.cross(b2, n1b)
    b2b = (tb - tangent * np.einsum("ij,ij->i", tangent, tb)[:, None]) / ell[:, None]
    b2b += np.cross(n1b, b1) + np.cross(b3, n2b)
    b3b = np.cross(n2b, b2)
    local_force = np.stack((b1b, -b1b + b2b, -b2b + b3b, -b3b), axis=1)
    return float(np.sum(u)), assemble(indices, local_force, len(x))


def eam_operator(x, pairs, parameters):
    b, source, target = incidence(pairs, len(x))
    d = b @ x
    r = np.linalg.norm(d, axis=1)
    z = np.asarray(parameters["types"])
    amplitude = np.asarray(parameters["density_amplitude"])
    embedding = np.asarray(parameters["embedding_coefficient"])[z]
    pairs = np.asarray(pairs)
    decay = parameters["density_decay"]
    p = amplitude[z[pairs[:, 0]], z[pairs[:, 1]]] * np.exp(-decay * r)
    q = amplitude[z[pairs[:, 1]], z[pairs[:, 0]]] * np.exp(-decay * r)
    rho = source.T @ p + target.T @ q
    node_derivative = 2.0 * embedding * rho
    phi = parameters["pair_amplitude"] * np.exp(-parameters["pair_decay"] * r)
    du = -parameters["pair_decay"] * phi - decay * (p * (source @ node_derivative) + q * (target @ node_derivative))
    g = (du / r)[:, None] * d
    return float(np.sum(embedding * rho ** 2) + np.sum(phi)), -b.T @ g, -d.T @ g


# The loop path uses Python math, indexed sums, and direct atom updates. It does
# not call the operator energy, local-gradient, incidence, or assembly functions.
def sub(a, b):
    return [float(a[k]) - float(b[k]) for k in range(3)]


def dot(a, b):
    return sum(a[k] * b[k] for k in range(3))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def norm(a):
    return math.sqrt(dot(a, a))


def lj_loop(x, pairs, parameters, box=None, scales=None, cutoff=None):
    energy = 0.0
    force = [[0.0] * 3 for _ in x]
    epsilon, sigma = parameters["epsilon"], parameters["sigma"]
    for e, (i, j) in enumerate(pairs):
        scale = 1.0 if scales is None else scales[e]
        if scale == 0.0:
            continue
        d = sub(x[j], x[i])
        if box is not None:
            d = [d[k] - float(box[k]) * round(d[k] / float(box[k])) for k in range(3)]
        r = norm(d)
        if r == 0.0:
            raise ValueError("Coincident particles in an active LJ pair")
        if cutoff is not None and r >= cutoff:
            continue
        u = 4.0 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
        derivative = 24.0 * epsilon / r * ((sigma / r) ** 6 - 2.0 * (sigma / r) ** 12)
        if cutoff is not None:
            uc = 4.0 * epsilon * ((sigma / cutoff) ** 12 - (sigma / cutoff) ** 6)
            duc = 24.0 * epsilon / cutoff * ((sigma / cutoff) ** 6 - 2.0 * (sigma / cutoff) ** 12)
            u -= uc + (r - cutoff) * duc
            derivative -= duc
        energy += scale * u
        for k in range(3):
            value = scale * derivative * d[k] / r
            force[i][k] += value
            force[j][k] -= value
    return energy, np.array(force)


def bonded_loop(x, parameters, kind):
    energy = 0.0
    force = [[0.0] * 3 for _ in x]
    for term, atoms in enumerate(parameters["indices"]):
        if kind == "bonds":
            i, j = atoms
            d = sub(x[j], x[i])
            r = norm(d)
            k, r0 = parameters["k"][term], parameters["r0"][term]
            energy += 0.5 * k * (r - r0) ** 2
            for axis in range(3):
                value = k * (r - r0) * d[axis] / r
                force[i][axis] += value
                force[j][axis] -= value
        elif kind == "angles":
            i, j, k = atoms
            a, b = sub(x[i], x[j]), sub(x[k], x[j])
            alpha, beta = norm(a), norm(b)
            cosine = dot(a, b) / (alpha * beta)
            theta = math.acos(cosine)
            spring, theta0 = parameters["k"][term], parameters["theta0"][term]
            energy += 0.5 * spring * (theta - theta0) ** 2
            prefactor = spring * (theta - theta0) / math.sin(theta)
            for axis in range(3):
                fi = prefactor / alpha * (b[axis] / beta - cosine * a[axis] / alpha)
                fk = prefactor / beta * (a[axis] / alpha - cosine * b[axis] / beta)
                force[i][axis] += fi
                force[k][axis] += fk
                force[j][axis] -= fi + fk
        else:
            i, j, k, l = atoms
            b1, b2, b3 = sub(x[j], x[i]), sub(x[k], x[j]), sub(x[l], x[k])
            n1, n2 = cross(b1, b2), cross(b2, b3)
            ell = norm(b2)
            # Triple-product atan2 and closed-form Cartesian torsion gradient;
            # this path does not reuse the reverse chain above.
            chi = math.atan2(ell * dot(b1, n2), dot(n1, n2))
            kap, mult, phase = (parameters[key][term] for key in ("kappa", "multiplicity", "phase"))
            energy += kap * (1.0 + math.cos(mult * chi - phase))
            du = -kap * mult * math.sin(mult * chi - phase)
            c1, c3 = dot(b1, b2) / ell ** 2, dot(b3, b2) / ell ** 2
            for axis in range(3):
                gi = -ell * n1[axis] / dot(n1, n1)
                gl = ell * n2[axis] / dot(n2, n2)
                gj = -(1.0 + c1) * gi + c3 * gl
                gk = c1 * gi - (1.0 + c3) * gl
                for atom, gradient in zip(atoms, (gi, gj, gk, gl)):
                    force[atom][axis] -= du * gradient
    return energy, np.array(force)


def eam_loop(x, pairs, parameters):
    neighbors = [[] for _ in x]
    for i, j in pairs:
        neighbors[i].append(j)
        neighbors[j].append(i)
    z = parameters["types"]
    amp = parameters["density_amplitude"]
    a = parameters["embedding_coefficient"]
    decay, pdecay = parameters["density_decay"], parameters["pair_decay"]
    rho = [sum(amp[z[i]][z[j]] * math.exp(-decay * norm(sub(x[j], x[i]))) for j in neighbors[i]) for i in range(len(x))]
    energy = sum(a[z[i]] * rho[i] ** 2 for i in range(len(x)))
    energy += sum(parameters["pair_amplitude"] * math.exp(-pdecay * norm(sub(x[j], x[i]))) for i, j in pairs)
    force = [[0.0] * 3 for _ in x]
    for i in range(len(x)):
        for j in neighbors[i]:
            d = sub(x[j], x[i])
            r = norm(d)
            derivative = -pdecay * parameters["pair_amplitude"] * math.exp(-pdecay * r)
            derivative -= decay * math.exp(-decay * r) * (2.0 * a[z[i]] * rho[i] * amp[z[i]][z[j]] + 2.0 * a[z[j]] * rho[j] * amp[z[j]][z[i]])
            for axis in range(3):
                force[i][axis] += derivative * d[axis] / r
    return energy, np.array(force)


def finite_difference(energy, x, step):
    force = np.zeros_like(x)
    for index in np.ndindex(x.shape):
        plus, minus = x.copy(), x.copy()
        plus[index] += step
        minus[index] -= step
        force[index] = -(energy(plus) - energy(minus)) / (2.0 * step)
    return force


def max_abs(value):
    return float(np.max(np.abs(value)))


def serialize(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=serialize, allow_nan=False) + "\n", encoding="utf-8")


def verify(config):
    x = np.asarray(config["positions"], dtype=np.float64)
    pairs = np.asarray(config["pairs"], dtype=np.int64)
    fd_config = config["finite_difference"]
    checks, models = {}, {}

    def check(name, residual, tolerance):
        residual = float(residual)
        checks[name] = {"residual": residual, "tolerance": tolerance, "passed": math.isfinite(residual) and abs(residual) <= tolerance}

    def model(name, coordinates, operator, loop):
        energy, force = operator(coordinates)[:2]
        reference_energy, reference_force = loop(coordinates)[:2]
        scale = max(fd_config["normalization_force_scale"], max_abs(force))
        sweep = []
        for step in fd_config["steps"]:
            fd_force = finite_difference(lambda y: loop(y)[0], coordinates, step)
            absolute = max_abs(force - fd_force)
            sweep.append({"step": step, "max_abs_force_error": absolute, "normalized_force_error": absolute / scale, "finite_difference_force": fd_force})
        table_row = next(row for row in sweep if row["step"] == fd_config["table_step"])
        check(name + "/energy_loop_operator", abs(energy - reference_energy), 1e-12 * max(1.0, abs(energy)))
        check(name + "/force_loop_operator", max_abs(force - reference_force), 1e-12 * scale)
        check(name + "/finite_difference", table_row["normalized_force_error"], fd_config["normalized_error_tolerance_at_table_step"])
        check(name + "/net_force", max_abs(np.sum(force, axis=0)), 1e-12 * scale)
        models[name] = {"operator_energy": energy, "loop_energy": reference_energy, "operator_force": force, "loop_force": reference_force, "normalization_denominator": scale, "finite_difference_sweep": sweep, "table_step_result": table_row}

    lj_op = lambda y: lj_operator(y, pairs, config["lj"])
    lj_ref = lambda y: lj_loop(y, pairs, config["lj"])
    model("lj", x, lj_op, lj_ref)
    for kind in ("bonds", "angles", "dihedrals"):
        parameters = config[kind]
        model(kind, x, lambda y, p=parameters, k=kind: bonded_operator(y, p, k), lambda y, p=parameters, k=kind: bonded_loop(y, p, k))
    for name in ("eam_single", "eam_multi"):
        parameters = config[name]
        model(name, x, lambda y, p=parameters: eam_operator(y, pairs, p), lambda y, p=parameters: eam_loop(y, pairs, p))
    periodic = config["periodic"]
    xp = np.asarray(periodic["positions"], dtype=np.float64)
    box = np.asarray(periodic["box_lengths"], dtype=np.float64)
    ppairs = np.asarray(periodic["pairs"], dtype=np.int64)
    model("lj_periodic", xp, lambda y: lj_operator(y, ppairs, config["lj"], box), lambda y: lj_loop(y, ppairs, config["lj"], box))
    bp, _, _ = incidence(ppairs, len(xp))
    image = np.rint((bp @ xp) / box).astype(np.int64)
    dp = bp @ xp - image * box
    branch_margin = float(np.min(box / 2.0 - np.abs(dp)))
    # All +/- FD perturbations remain on the same image branch.
    check("periodic/stencil_branch_margin", max(0.0, 2.0 * max(fd_config["steps"]) - branch_margin), 0.0)

    scales = config["supplemental"]["pair_scales"]
    rc = config["supplemental"]["force_shift_cutoff"]
    model("lj_scaled_force_shifted", x,
          lambda y: lj_operator(y, pairs, config["lj"], scales=scales, cutoff=rc),
          lambda y: lj_loop(y, pairs, config["lj"], scales=scales, cutoff=rc))
    coincident = np.zeros((2, 3))
    excluded_u, excluded_f, _ = lj_operator(coincident, [[0, 1]], config["lj"], scales=[0.0])
    check("excluded_coincident_pair", abs(excluded_u) + max_abs(excluded_f), 0.0)
    uc, duc = lj_values(np.array([rc]), **config["lj"], cutoff=rc)
    check("force_shift/at_cutoff", max_abs(uc) + max_abs(duc), 0.0)
    for offset in (1e-3, 1e-4, 1e-5):
        left_u, left_du = lj_values(np.array([rc - offset]), **config["lj"], cutoff=rc)
        # Record one-sided approach without asserting smooth second derivatives.
        checks["force_shift/left_" + str(offset)] = {"offset": offset, "energy": float(left_u[0]), "radial_derivative": float(left_du[0])}

    u, f, w = lj_op(x)
    full = np.vstack((pairs, pairs[:, ::-1]))
    bf, sf, _ = incidence(full, len(x))
    df = bf @ x
    rf = np.linalg.norm(df, axis=1)
    uf, duf = lj_values(rf, **config["lj"])
    gf = (duf / rf)[:, None] * df
    check("half_full/energy", abs(u - 0.5 * np.sum(uf)), 1e-12)
    check("half_full/adjoint_force", max_abs(f + 0.5 * bf.T @ gf), 1e-12)
    check("half_full/source_force", max_abs(f - sf.T @ gf), 1e-12)
    adjoint = config["adjoint"]
    indices = np.asarray(adjoint["indices"], dtype=np.int64)
    z = np.random.default_rng(adjoint["seed"]).standard_normal(tuple(adjoint["shape"]))
    assembled = assemble(indices, z, len(x))
    lhs, rhs = float(np.sum(x[indices] * z)), float(np.sum(x * assembled))
    check("gather/adjoint", abs(lhs - rhs), 1e-12)

    # Reversal and permutation change only representation, not the pair model.
    for name, edges in (("reversal", pairs[:, ::-1]), ("permutation", pairs[::-1])):
        ur, fr, _ = lj_operator(x, edges, config["lj"])
        check("invariance/edge_" + name, abs(ur - u) + max_abs(fr - f), 1e-12)
    permutation = np.array([2, 4, 0, 3, 1])  # old label -> new label
    relabeled = np.empty_like(x)
    relabeled[permutation] = x
    ur, fr, _ = lj_operator(relabeled, permutation[pairs], config["lj"])
    check("invariance/atom_relabeling", abs(ur - u) + max_abs(fr[permutation] - f), 1e-12)
    for name in ("lj", "bonds", "angles", "dihedrals", "eam_single", "eam_multi"):
        fm = models[name]["operator_force"]
        check(name + "/nonperiodic_torque", max_abs(np.sum(np.cross(x, fm), axis=0)), 1e-12 * max(1.0, max_abs(fm)))
    check("virial/nonperiodic_XTF", max_abs(w - x.T @ f), 1e-12)

    # Deform both particle coordinates and the periodic cell, at fixed images.
    # The scalar reference energy uses the deformed displacements explicitly.
    strain_results = {}
    hstrain = config["supplemental"]["strain_step"]
    for name, coordinates, edges, lengths in (("lj", x, pairs, None), ("lj_periodic", xp, ppairs, box)):
        b, _, _ = incidence(edges, len(coordinates))
        displacement = b @ coordinates
        if lengths is not None:
            displacement -= np.rint(displacement / lengths) * lengths
        analytic_w = lj_operator(coordinates, edges, config["lj"], lengths)[2]
        def strain_energy(strain):
            deformed = displacement @ (np.eye(3) + strain)
            return sum(4.0 * config["lj"]["epsilon"] * ((config["lj"]["sigma"] / norm(row)) ** 12 - (config["lj"]["sigma"] / norm(row)) ** 6) for row in deformed)
        fd_w = finite_difference(strain_energy, np.zeros((3, 3)), hstrain)
        residual = max_abs(analytic_w - fd_w)
        check("virial/strain_" + name, residual / max(1.0, max_abs(analytic_w)), 1e-7)
        strain_results[name] = {"analytic_virial": analytic_w, "finite_difference_virial": fd_w, "max_abs_error": residual, "step": hstrain}

    vv = config["verlet"]
    v = np.asarray(vv["velocities"], dtype=np.float64)
    mass = np.asarray(vv["masses"])
    dt = vv["time_step"]
    vh = v + 0.5 * dt * f / mass[:, None]
    x1 = x + dt * vh
    f1 = lj_op(x1)[1]
    v1 = vh + 0.5 * dt * f1 / mass[:, None]
    fr0 = lj_ref(x)[1]
    vrh, xr1 = np.empty_like(v), np.empty_like(x)
    for i in range(len(x)):
        for axis in range(3):
            vrh[i, axis] = v[i, axis] + dt / (2.0 * mass[i]) * fr0[i, axis]
            xr1[i, axis] = x[i, axis] + dt * vrh[i, axis]
    fr1 = lj_ref(xr1)[1]
    vr1 = np.empty_like(v)
    for i in range(len(x)):
        for axis in range(3):
            vr1[i, axis] = vrh[i, axis] + dt / (2.0 * mass[i]) * fr1[i, axis]
    check("verlet/position", max_abs(x1 - xr1), 1e-12)
    check("verlet/velocity", max_abs(v1 - vr1), 1e-12)
    check("verlet/final_force", max_abs(f1 - fr1), 1e-12)

    # Explicit counterexample to an energy jump at a nearest-image tie.
    # Length 4, epsilon=sigma=1: continuous energy, opposite one-sided slopes.
    tie_offsets = [1e-2, 1e-4, 1e-6]
    def tie_energy(displacement):
        radius = abs(displacement - 4.0 * round(displacement / 4.0))
        return 4.0 * (radius ** -12 - radius ** -6)
    tie_samples = [{"offset": e, "left_energy": tie_energy(2.0 - e), "right_energy": tie_energy(2.0 + e)} for e in tie_offsets]
    check("periodic/tie_equal_energy", max(abs(row["left_energy"] - row["right_energy"]) for row in tie_samples), 1e-12)
    # A skew-cell counterexample to unrestricted fractional rounding.
    cell = np.array([[1.0, 0.0, 0.0], [0.9, 0.4, 0.0], [0.0, 0.0, 5.0]])
    displacement = np.array([0.5, 0.18, 0.0])
    rounded_image = np.rint(displacement @ np.linalg.inv(cell)).astype(int)
    candidates = np.array(list(itertools.product(range(-2, 3), repeat=3)))
    distances = np.linalg.norm(displacement - candidates @ cell, axis=1)
    closest = int(np.argmin(distances))
    # Enumeration is a bounded illustrative check, not a general image solver.
    check("periodic/skew_counterexample", 0.0 if distances[closest] < np.linalg.norm(displacement - rounded_image @ cell) else 1.0, 0.0)

    return {
        "schema_version": 1,
        "provenance": config["provenance"],
        "inputs": config,
        "models": models,
        "checks": checks,
        "adjoint": {"random_tensor": z, "assembled_tensor": assembled, "lhs": lhs, "rhs": rhs},
        "representation_checks": {"atom_permutation_old_to_new": permutation, "edge_permutation": "reverse the configured pair order", "edge_reversal": "swap both endpoints of every configured pair", "full_list": full},
        "supplemental_fixtures": {"excluded_pair_positions": coincident, "excluded_pair_indices": [[0, 1]], "excluded_pair_scales": [0.0], "cutoff_approach_offsets": [1e-3, 1e-4, 1e-5], "strain_reference": np.zeros((3, 3))},
        "periodic_geometry": {"image_vectors": image, "displacements": dp, "minimum_distance_to_image_tie": branch_margin},
        "virial": strain_results,
        "verlet": {"operator_positions": x1, "loop_positions": xr1, "operator_velocities": v1, "loop_velocities": vr1, "operator_final_force": f1, "loop_final_force": fr1},
        "image_tie_example": {"box_length": 4.0, "epsilon": 1.0, "sigma": 1.0, "energy_at_tie": tie_energy(2.0), "left_derivative_limit": 0.181640625, "right_derivative_limit": -0.181640625, "samples": tie_samples},
        "skew_example": {"cell": cell, "displacement": displacement, "rounded_image": rounded_image, "rounded_distance": float(np.linalg.norm(displacement - rounded_image @ cell)), "enumeration_range": [-2, 2], "shorter_image": candidates[closest], "shorter_distance": float(distances[closest])},
        "all_checks_passed": all(item["passed"] for item in checks.values() if "passed" in item),
    }


def latex_number(value):
    if value == 0.0:
        return "0"
    mantissa, exponent = f"{value:.4e}".split("e")
    return mantissa + r"\times10^{" + str(int(exponent)) + "}"


def table_rows(results):
    names = [("lj", "Lennard--Jones"), ("bonds", "Harmonic bonds"), ("angles", "Harmonic angles"), ("dihedrals", "Periodic dihedral potential"), ("eam_single", "Single-element EAM form"), ("eam_multi", "Multi-element EAM form"), ("lj_periodic", "LJ in a fixed orthorhombic box")]
    lines = []
    for key, label in names:
        row = results["models"][key]["table_step_result"]
        lines.append(label + " & $" + latex_number(row["max_abs_force_error"]) + "$ & $" + latex_number(row["normalized_force_error"]) + r"$\\")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.json")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--check-manuscript", type=Path, help="Check the marked numerical table against this run; never edit the manuscript")
    args = parser.parse_args()
    config_bytes = args.config.read_bytes()
    config = json.loads(config_bytes)
    np.seterr(divide="raise", invalid="raise", over="raise")
    results = verify(config)
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    config_hash = hashlib.sha256(config_bytes).hexdigest()
    results["source_sha256"] = {"verify_equivalence.py": script_hash, "config.json": config_hash}
    environment = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "python_implementation": platform.python_implementation(),
        "numpy_version": np.__version__,
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "float64_epsilon": float(np.finfo(np.float64).eps),
        "numpy_build_configuration": np.show_config(mode="dicts"),
        "source_sha256": results["source_sha256"],
        "execution": "CPU only; no GMD engine, CUDA execution, LaTeX compilation, or external MD package.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = table_rows(results)
    write_json(args.output_dir / "verification_results.json", results)
    environment["verification_results_sha256"] = hashlib.sha256((args.output_dir / "verification_results.json").read_bytes()).hexdigest()
    write_json(args.output_dir / "environment.json", environment)
    (args.output_dir / "table_rows.tex").write_text(rows, encoding="utf-8")
    if args.check_manuscript is not None:
        manuscript = args.check_manuscript.read_text(encoding="utf-8")
        recorded_rows = manuscript.split("% BEGIN GENERATED VERIFICATION ROWS\n", 1)[1].split("% END GENERATED VERIFICATION ROWS", 1)[0]
        if recorded_rows != rows:
            raise SystemExit("Manuscript table differs from this run. Inspect results and environment; do not assume portable bitwise agreement.")
    counted = [item for item in results["checks"].values() if "passed" in item]
    print(f"{sum(item['passed'] for item in counted)}/{len(counted)} checks passed")
    for key, value in results["models"].items():
        row = value["table_step_result"]
        print(f"{key:25s} absolute={row['max_abs_force_error']:.4e} normalized={row['normalized_force_error']:.4e}")
    if not results["all_checks_passed"]:
        for name, item in results["checks"].items():
            if item.get("passed") is False:
                print("FAILED:", name, item)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
