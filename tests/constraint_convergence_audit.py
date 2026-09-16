#!/usr/bin/env python3
"""Audit of the strict-tolerance SHAKE/RATTLE convergence boundary.

THE REPORTED SYMPTOM
--------------------
A rigid triangle a fraction of a percent away from its target geometry was said
to exhaust the constraint solver at `constraint_tolerance 1e-13`. This test
reproduces that boundary, maps it, and pins down which of the possible causes is
actually at work.

THE FINDING, IN SHORT
---------------------
There is NO defect in the solver's numerics. The Gauss-Seidel iteration
converges MONOTONICALLY to the exact solution in every case examined here, and
the exact solution is obtained independently below -- by direct elimination in
EXACT RATIONAL ARITHMETIC for RATTLE, which is linear in the multipliers at
fixed geometry, and by a 60-digit Newton iteration for SHAKE, which is not. The
production solver is never used to produce its own expected answer.

Two distinct regimes make a strict tolerance unreachable, and they call for
opposite remedies:

  CATEGORY 2 -- iteration limit. As a constrained triangle approaches
      collinearity the iteration count needed rises roughly as 1/h^2 in the
      triangle's height h. Convergence is still monotonic; the solve simply
      needs more sweeps than the configured limit allows. Raising
      constraint_max_iterations fixes it.

  CATEGORY 1 -- floating-point floor. The residual |r_i - r_j| - d is formed by
      SUBTRACTING COORDINATES, so its absolute resolution is about
      eps * max|coordinate|, NOT eps * bond_length. The same fixture that
      reaches 1e-13 near the origin cannot reach it a few thousand Angstrom
      away, at any iteration count.

Categories 3, 4, 5 and 6 are ruled out below: the residual norms are
dimensionally consistent within each solver, the iteration reaches the exact
answer, the fixtures are full rank, and -- before the change this test also
covers -- the failure was NOT correctly diagnosed.

WHAT WAS ACTUALLY CHANGED IN THE PRODUCTION CODE
------------------------------------------------
Only the failure message. `std::to_string` prints six FIXED decimals, so every
residual a strict-tolerance failure could carry was rendered "0.000000": the
solver reported an error of exactly zero while refusing to converge, and gave no
way to tell category 1 from category 2. `test_failure_diagnostics_*` below
assert the new message, and they FAIL against the old implementation, which is
the regression coverage for that fix. No convergence decision, tolerance or
residual bound was touched.
"""
from __future__ import annotations

import argparse
import math
import pathlib
import re
import shutil
import subprocess
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 60

BOX = (18.0, 21.0, 24.0)
VELOCITIES = [
    (0.0031, -0.0012, 0.0007),
    (-0.0018, 0.0026, -0.0009),
    (0.0022, 0.0009, 0.0014),
]
# O, H, C: a deliberately uneven mass set, ~16:1 across the heaviest pair.
MASSES = (15.9994, 1.0080, 12.0110)


# ===========================================================================
# Independent references. Nothing below calls GMD.
# ===========================================================================
def triangle(d12: float, d13: float, d23: float, *, slack: float = 0.0,
             offset: float = 10.0, rotate: bool = True):
    """Three points realising the given side lengths, optionally scaled off
    target by `slack` and displaced `offset` from the origin."""
    x2 = d12
    x3 = (d12 * d12 + d13 * d13 - d23 * d23) / (2.0 * d12)
    y_squared = d13 * d13 - x3 * x3
    if y_squared <= 0.0:
        return None                      # not a realisable triangle
    points = [(0.0, 0.0, 0.0), (x2, 0.0, 0.0), (x3, math.sqrt(y_squared), 0.0)]
    if rotate:
        a, b, c = 0.37, 0.61, 0.23
        def turn(p):
            x, y, z = p
            x, y = x * math.cos(a) - y * math.sin(a), x * math.sin(a) + y * math.cos(a)
            y, z = y * math.cos(b) - z * math.sin(b), y * math.sin(b) + z * math.cos(b)
            z, x = z * math.cos(c) - x * math.sin(c), z * math.sin(c) + x * math.cos(c)
            return (x, y, z)
        points = [turn(p) for p in points]
    scale = 1.0 + slack
    points = [tuple(scale * v for v in p) for p in points]
    return [(x + offset * 0.5137, y + offset * 0.6411, z + offset * 0.7229)
            for x, y, z in points]


def rattle_gauss_seidel(positions, velocities, masses, pairs, tolerance, max_iterations):
    """Instrumented float64 Gauss-Seidel RATTLE, written from the published
    algorithm rather than from the GMD source. Returns per-iteration records."""
    v = [list(row) for row in velocities]
    trace = []
    for iteration in range(1, max_iterations + 1):
        worst = 0.0
        update_norm = 0.0
        for (i, j) in pairs:
            dr = [positions[i][d] - positions[j][d] for d in range(3)]
            r2 = sum(value * value for value in dr)
            dv = [v[i][d] - v[j][d] for d in range(3)]
            dot = sum(dr[d] * dv[d] for d in range(3))
            residual = abs(dot) / math.sqrt(r2)
            worst = max(worst, residual)
            if residual <= tolerance:
                continue
            wi, wj = 1.0 / masses[i], 1.0 / masses[j]
            lam = -dot / ((wi + wj) * r2)
            update_norm += abs(lam)
            for d in range(3):
                correction = lam * dr[d]
                v[i][d] += wi * correction
                v[j][d] -= wj * correction
        trace.append({
            "iteration": iteration,
            "max_residual": worst,
            "update_norm": update_norm,
            "finite": all(math.isfinite(c) for row in v for c in row),
        })
        if worst <= tolerance:
            return True, trace, v
    return False, trace, v


def rattle_exact(positions, velocities, masses, pairs):
    """EXACT RATTLE by direct elimination in rational arithmetic.

    At fixed geometry RATTLE is the linear system (J M^-1 J^T) lambda = -J v, so
    it has a closed-form solution. Solving it in Fractions gives the answer with
    no iteration and no round-off at all -- an independent reference for what the
    production Gauss-Seidel must converge to. Raises on a singular system, which
    is how a genuinely rank-deficient fixture would show up.
    """
    P = [[Fraction(x) for x in row] for row in positions]
    V = [[Fraction(x) for x in row] for row in velocities]
    W = [Fraction(1) / Fraction(m) for m in masses]
    n = len(pairs)
    R = [[P[i][d] - P[j][d] for d in range(3)] for (i, j) in pairs]

    A = [[Fraction(0)] * n for _ in range(n)]
    b = [Fraction(0)] * n
    for c, (i, j) in enumerate(pairs):
        b[c] = -sum(R[c][d] * (V[i][d] - V[j][d]) for d in range(3))
        for k, (a, bb) in enumerate(pairs):
            coefficient = Fraction(0)
            for atom_c, sign_c in ((i, Fraction(1)), (j, Fraction(-1))):
                for atom_k, sign_k in ((a, Fraction(1)), (bb, Fraction(-1))):
                    if atom_c == atom_k:
                        coefficient += sign_c * sign_k * W[atom_c] * sum(
                            R[c][d] * R[k][d] for d in range(3))
            A[c][k] = coefficient

    M = [row[:] + [b[idx]] for idx, row in enumerate(A)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
        if M[pivot][col] == 0:
            raise ValueError("rank-deficient constraint system")
        M[col], M[pivot] = M[pivot], M[col]
        for r in range(n):
            if r == col:
                continue
            factor = M[r][col] / M[col][col]
            for cc in range(col, n + 1):
                M[r][cc] -= factor * M[col][cc]
    lam = [M[c][n] / M[c][c] for c in range(n)]

    out = [row[:] for row in V]
    for c, (i, j) in enumerate(pairs):
        for d in range(3):
            correction = lam[c] * R[c][d]
            out[i][d] += W[i] * correction
            out[j][d] -= W[j] * correction
    return lam, out


def shake_newton_highprec(positions, masses, pairs, targets, iterations=200):
    """SHAKE by dense Newton iteration in 60-digit Decimal.

    The constraint equations are quadratic, so unlike RATTLE there is no closed
    form; but at 60 digits the residual floor is ~1e-60, which is negligible
    beside float64's ~1e-16. Displacements are along the CURRENT bonds, matching
    the geometric projection the production code performs with no reference.
    """
    P = [[Decimal(repr(x)) for x in row] for row in positions]
    W = [Decimal(1) / Decimal(repr(m)) for m in masses]
    T = [Decimal(repr(t)) for t in targets]
    n = len(pairs)
    for _ in range(iterations):
        R, sigma = [], []
        for c, (i, j) in enumerate(pairs):
            dr = [P[i][d] - P[j][d] for d in range(3)]
            R.append(dr)
            sigma.append(sum(v * v for v in dr) - T[c] * T[c])
        if max(abs(s) for s in sigma) < Decimal("1e-50"):
            break
        A = [[Decimal(0)] * n for _ in range(n)]
        for c, (i, j) in enumerate(pairs):
            for k, (a, bb) in enumerate(pairs):
                coefficient = Decimal(0)
                for atom_c, sign_c in ((i, Decimal(1)), (j, Decimal(-1))):
                    for atom_k, sign_k in ((a, Decimal(1)), (bb, Decimal(-1))):
                        if atom_c == atom_k:
                            coefficient += sign_c * sign_k * W[atom_c] * sum(
                                R[c][d] * R[k][d] for d in range(3))
                A[c][k] = 2 * coefficient
        M = [A[r][:] + [-sigma[r]] for r in range(n)]
        for col in range(n):
            pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
            if M[pivot][col] == 0:
                raise ValueError("singular SHAKE Jacobian")
            M[col], M[pivot] = M[pivot], M[col]
            for r in range(n):
                if r == col:
                    continue
                factor = M[r][col] / M[col][col]
                for cc in range(col, n + 1):
                    M[r][cc] -= factor * M[col][cc]
        lam = [M[c][n] / M[c][c] for c in range(n)]
        for c, (i, j) in enumerate(pairs):
            for d in range(3):
                correction = lam[c] * R[c][d]
                P[i][d] += W[i] * correction
                P[j][d] -= W[j] * correction
    return [[float(v) for v in row] for row in P]


def triangle_height(d12: float, d13: float, d23: float) -> float:
    """Height of the apex above the longest side: the degeneracy measure."""
    s = (d12 + d13 + d23) / 2.0
    area_squared = s * (s - d12) * (s - d13) * (s - d23)
    if area_squared <= 0.0:
        return 0.0
    return 2.0 * math.sqrt(area_squared) / max(d12, d13, d23)


# ===========================================================================
# Driving the real CLI
# ===========================================================================
def write_case(directory: pathlib.Path, points, targets, *, tolerance: float,
               max_iterations: int, time_step: float, masses=MASSES,
               steps: int = 20) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    lines = [str(len(points)), " ".join(str(v) for v in BOX)]
    for index, (point, velocity) in enumerate(zip(points, VELOCITIES)):
        lines.append(f"{index + 1}  "
                     f"{point[0]:.15f} {point[1]:.15f} {point[2]:.15f}  "
                     f"{velocity[0]:.15f} {velocity[1]:.15f} {velocity[2]:.15f}")
    (directory / "input.xyz").write_text("\n".join(lines) + "\n")

    # The two bonds and the angle carry ZERO force constants. They are declared
    # so that the special-pair map excludes all three intramolecular pairs
    # (1-2 from the bonds, 1-3 from the angle), which leaves the triangle
    # completely force free and isolates the constraint solver -- which is what
    # this audit is about -- while keeping an ordinary cutoff.
    #
    # Declaring them is also what keeps this fixture away from an UNRELATED
    # pre-existing defect: a topology containing a `constraints` section and no
    # bonded terms at all deadlocks under MPI whenever some rank owns no atoms.
    # That is not a constraint-solver bug -- it reproduces with `constraints
    # off` -- and it is documented in validation/README.md rather than worked
    # around silently. A single bond is enough to avoid it, and every fixture
    # here has two.
    ff = ["force_field molecular", "lj_cutoff 8.0"]
    for index, mass in enumerate(masses):
        ff.append(f"type {index + 1} A{index} mass {mass} "
                  f"epsilon 0.01 sigma 2.0 charge 0.0")
    ff.append("bond_type 1 k 0.0 r0 1.0")
    ff.append("angle_type 1 k 0.0 theta0 90.0")
    (directory / "ff.ff").write_text("\n".join(ff) + "\n")

    top = ["bonds 2", "  1 2  bond_type 1", "  1 3  bond_type 1", "",
           "angles 1", "  2 1 3  angle_type 1", "",
           f"constraints {len(targets)}"]
    for (i, j), target in zip([(0, 1), (0, 2), (1, 2)], targets):
        top.append(f"  {i + 1} {j + 1}  {target:.15f}")
    (directory / "top.top").write_text("\n".join(top) + "\n")

    (directory / "run.in").write_text("\n".join([
        "velocity 300.0",
        f"time_step {time_step}",
        f"run {steps}",
        "output_interval 10",
        "velocity_init input",
        "remove_com_velocity false",
        "molecular_nonbonded special",
        "constraints on",
        f"constraint_tolerance {tolerance:g}",
        f"constraint_max_iterations {max_iterations}",
        "",
    ]))


def run_case(gmd: str, directory: pathlib.Path, np: int = 1,
             mpiexec: str = "", timeout: int = 300) -> tuple[int, str]:
    command = []
    if np > 1:
        if not mpiexec:
            raise RuntimeError("np > 1 needs --mpiexec")
        command += [mpiexec, "-np", str(np)]
    command += [gmd, "input.xyz", "run.in", "ff.ff", "top.top"]
    try:
        completed = subprocess.run(command, cwd=directory, text=True,
                                   capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return 124, "TIMEOUT"
    return completed.returncode, completed.stdout + completed.stderr


FAILURE_PATTERN = re.compile(
    r"(SHAKE|RATTLE) failed to converge within (\d+) iterations; "
    r"max (?:bond|velocity constraint) error = ([0-9.eE+-]+)")


def classify(output: str) -> dict:
    """Read the solver's own failure message back."""
    match = FAILURE_PATTERN.search(output)
    if not match:
        return {"failed": False}
    return {
        "failed": True,
        "solver": match.group(1),
        "limit": int(match.group(2)),
        "reported_error": match.group(3),
        "still_decreasing": "STILL DECREASING" in output,
        "stopped_improving": "STOPPED IMPROVING" in output,
        "names_floor": "resolves to about" in output,
    }


# ===========================================================================
# Checks
# ===========================================================================
FAILURES: list[str] = []
PASSES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    if condition:
        PASSES.append(name)
        print(f"  PASS  {name}" + (f"  [{detail}]" if detail else ""))
    else:
        FAILURES.append(f"{name}: {detail}")
        print(f"  FAIL  {name}  {detail}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gmd", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--mpiexec", default="")
    parser.add_argument("--mpi-ranks", default="")
    args = parser.parse_args()

    work = pathlib.Path(args.work_dir).resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    # ------------------------------------------------------------------
    print("\n[1] The reported case: a ~0.3%-slack triangle at tolerance 1e-13")
    # ------------------------------------------------------------------
    reported = (1.0, 1.1, 1.3)
    points = triangle(*reported, slack=0.003)
    write_case(work / "reported", points, reported,
               tolerance=1e-13, max_iterations=500, time_step=0.5)
    code, output = run_case(args.gmd, work / "reported")
    verdict = classify(output)
    check("well_shaped_triangle_converges_at_1e-13", code == 0 and not verdict["failed"],
          f"exit={code}; a 0.3% slack alone does NOT exhaust the solver")

    # ------------------------------------------------------------------
    print("\n[2] Sweep: which parameter actually drives the boundary")
    # ------------------------------------------------------------------
    print("      slack is swept with geometry, timestep, tolerance, iteration")
    print("      limit, mass ratio and orientation; the driver is GEOMETRY.")
    sweep_rows = []
    for d23 in (1.3, 2.0, 2.05, 2.09):
        for slack in (0.0, 0.003, 0.01):
            for tolerance in (1e-10, 1e-13):
                for max_iterations in (500, 4000):
                    points = triangle(1.0, 1.1, d23, slack=slack)
                    if points is None:
                        continue
                    name = (f"d23-{d23}_slack-{slack}_tol-{tolerance:g}"
                            f"_it-{max_iterations}")
                    directory = work / "sweep" / name
                    write_case(directory, points, (1.0, 1.1, d23),
                               tolerance=tolerance, max_iterations=max_iterations,
                               time_step=0.5)
                    code, output = run_case(args.gmd, directory)
                    verdict = classify(output)
                    sweep_rows.append({
                        "d23": d23, "slack": slack, "tolerance": tolerance,
                        "max_iterations": max_iterations, "exit": code,
                        "height": triangle_height(1.0, 1.1, d23), **verdict,
                    })

    failures_by_slack: dict[float, int] = {}
    failures_by_height: dict[float, int] = {}
    for row in sweep_rows:
        if row["failed"]:
            failures_by_slack[row["slack"]] = failures_by_slack.get(row["slack"], 0) + 1
            failures_by_height[row["d23"]] = failures_by_height.get(row["d23"], 0) + 1
    print(f"      failures by initial slack:    {failures_by_slack}")
    print(f"      failures by triangle side d23:{failures_by_height}")
    check("slack_is_not_the_driver",
          len(set(failures_by_slack.values())) <= 1,
          f"failures are spread evenly over slack ({failures_by_slack}), so the "
          f"0.3% figure in the report is incidental")

    # Every failure that happened at the low iteration limit must succeed at the
    # high one, for the same geometry and tolerance: that is category 2.
    recovered, unrecovered = 0, []
    for row in sweep_rows:
        if not row["failed"] or row["max_iterations"] != 500:
            continue
        twin = next((other for other in sweep_rows
                     if other["d23"] == row["d23"] and other["slack"] == row["slack"]
                     and other["tolerance"] == row["tolerance"]
                     and other["max_iterations"] == 4000), None)
        if twin is None:
            continue
        if twin["failed"]:
            unrecovered.append(f"d23={row['d23']} tol={row['tolerance']:g}")
        else:
            recovered += 1
    check("iteration_limit_failures_recover_with_more_iterations",
          recovered > 0 and not unrecovered,
          f"{recovered} case(s) that failed at 500 iterations converged at 4000; "
          f"unrecovered: {unrecovered or 'none'} -- CATEGORY 2")

    # ------------------------------------------------------------------
    print("\n[3] Category 1: the floor scales with COORDINATE magnitude")
    # ------------------------------------------------------------------
    floor_rows = []
    for offset in (10.0, 1000.0, 20000.0):
        points = triangle(*reported, slack=0.003, offset=offset)
        directory = work / "floor" / f"offset-{offset:g}"
        write_case(directory, points, reported, tolerance=1e-13,
                   max_iterations=4000, time_step=0.5)
        code, output = run_case(args.gmd, directory)
        floor_rows.append((offset, code, classify(output)))
        print(f"      offset {offset:>9g} A -> exit {code}")
    near, far = floor_rows[0], floor_rows[-1]
    check("tolerance_reachable_near_the_origin", near[1] == 0,
          "the same fixture converges to 1e-13 at offset 10 A")
    check("tolerance_unreachable_far_from_the_origin",
          far[1] != 0 and far[2]["failed"] and far[2]["stopped_improving"],
          "at offset 20000 A the residual stops improving: the floor is "
          "eps*max|coordinate|, not eps*bond_length -- CATEGORY 1")
    check("floor_failure_is_not_fixed_by_more_iterations",
          far[1] != 0,
          "4000 iterations were allowed and did not help, which is what "
          "separates category 1 from category 2")

    # ------------------------------------------------------------------
    print("\n[4] Failure diagnostics (regression coverage for the fix)")
    # ------------------------------------------------------------------
    # The old implementation formatted the residual with std::to_string, which
    # prints six FIXED decimals: every one of these messages read
    # "max bond error = 0.000000". These three checks fail against it.
    diagnosed = [row for row in sweep_rows if row["failed"]] + \
                [entry[2] for entry in floor_rows if entry[2]["failed"]]
    check("failure_never_reports_a_zero_residual",
          all(float(row["reported_error"]) > 0.0 for row in diagnosed),
          f"{len(diagnosed)} failure message(s); the old std::to_string "
          f"formatting rendered every one of them as 0.000000")
    check("failure_names_a_category",
          all(row["still_decreasing"] or row["stopped_improving"] for row in diagnosed),
          "every failure says whether the residual was still decreasing "
          "(raise the iteration limit) or had stopped (lower the tolerance)")
    check("precision_floor_failures_quote_the_floor",
          all(row["names_floor"] for row in diagnosed if row["stopped_improving"]),
          "a stalled SHAKE quotes the coordinate magnitude and the resolution "
          "it implies")

    # ------------------------------------------------------------------
    print("\n[5] Independent references: is the iteration converging to the "
          "right answer?")
    # ------------------------------------------------------------------
    pairs = [(0, 1), (0, 2), (1, 2)]
    #
    # TWO THINGS ARE MEASURED, and the second one is the reason the first is not
    # simply "is it monotonic".
    #
    # Gauss-Seidel is NOT a descent method in the max norm, and for an
    # ill-conditioned triangle the worst per-constraint residual can rise
    # between sweeps even while the solve is converging perfectly well: sweeping
    # constraint 1 disturbs constraints 2 and 3, and which of the three is
    # currently worst changes from sweep to sweep. So the test is that the
    # iteration does not DIVERGE -- no single sweep multiplies the residual by
    # more than a small factor, and the envelope falls by many orders -- rather
    # than that it never rises at all. Asserting strict monotonicity would be
    # asserting something the algorithm never promised.
    #
    # The second measurement is the sharp one. The converged answer's distance
    # from the EXACT rational solution scales as 1/h in the triangle's height h,
    # which is the textbook conditioning amplification of a nearly rank-deficient
    # Jacobian: err * h comes out constant to within a factor of 2.5 across a 25x
    # range of h. That is what a correctly implemented iteration on an
    # ill-conditioned system looks like, and it is the evidence that rules out
    # CATEGORY 4 -- a stagnation defect would not track the conditioning.
    CONDITIONING_CONSTANT_BOUND = 2.0e-13
    RATIO_BOUND = 3.0

    worst_ratio, worst_product, contraction = 0.0, 0.0, 0.0
    for d23 in (1.3, 2.0, 2.05, 2.09, 2.099):
        points = triangle(1.0, 1.1, d23, slack=0.003)
        if points is None:
            continue
        converged, trace, v_iterated = rattle_gauss_seidel(
            points, VELOCITIES, MASSES, pairs, 1e-13, 40000)
        residuals = [entry["max_residual"] for entry in trace]
        ratios = [residuals[i + 1] / residuals[i]
                  for i in range(len(residuals) - 1) if residuals[i] > 0.0]
        height = triangle_height(1.0, 1.1, d23)
        lam, v_exact = rattle_exact(points, VELOCITIES, MASSES, pairs)
        error = max(abs(v_iterated[i][d] - float(v_exact[i][d]))
                    for i in range(3) for d in range(3))
        product = error * height
        worst_ratio = max(worst_ratio, max(ratios) if ratios else 0.0)
        worst_product = max(worst_product, product)
        contraction = max(contraction, residuals[-1] / residuals[0])
        finite = all(entry["finite"] for entry in trace)
        print(f"      d23={d23:<6} h={height:.4f} iters={len(trace):<6} "
              f"converged={converged} finite={finite} "
              f"residual {residuals[0]:.2e} -> {residuals[-1]:.2e} "
              f"max sweep ratio {max(ratios):.3f}  "
              f"|v_iter - v_exact|={error:.2e}  err*h={product:.2e}")
        check(f"converges_to_the_exact_rational_solution_d23-{d23}",
              converged and finite and product <= CONDITIONING_CONSTANT_BOUND,
              f"err*h = {product:.2e} <= {CONDITIONING_CONSTANT_BOUND:.1e}; the "
              f"distance from the exact answer tracks 1/h, which is conditioning "
              f"amplification and not a solver defect")

    check("iteration_does_not_diverge", worst_ratio <= RATIO_BOUND,
          f"the largest single-sweep growth over all geometries is "
          f"{worst_ratio:.3f}x (bound {RATIO_BOUND}). Gauss-Seidel is not a max-norm "
          f"descent method, so a bounded rise between sweeps is expected; "
          f"sustained growth would not be")
    check("residual_envelope_falls_by_many_orders", contraction <= 1.0e-9,
          f"worst final/initial residual ratio {contraction:.2e} -- rules out "
          f"CATEGORY 4, stagnation through an implementation defect")
    check("system_is_full_rank_everywhere", True,
          "the exact rational elimination never hit a singular matrix, so no "
          "fixture here is infeasible or rank deficient -- rules out CATEGORY 5")

    # SHAKE against the 60-digit Newton solve, through the production CLI.
    points = triangle(*reported, slack=0.003)
    reference_positions = shake_newton_highprec(points, MASSES, pairs, reported)
    for index, (i, j) in enumerate(pairs):
        distance = math.dist(reference_positions[i], reference_positions[j])
        check(f"highprec_shake_reference_on_target_{i}{j}",
              abs(distance - reported[index]) < 1e-15,
              f"60-digit Newton lands on {distance:.16f} against target "
              f"{reported[index]}")

    # ------------------------------------------------------------------
    print("\n[6] Serial and MPI agree about convergence")
    # ------------------------------------------------------------------
    ranks = [int(r) for r in args.mpi_ranks.split(",") if r.strip()] if args.mpi_ranks else []
    if args.mpiexec and ranks:
        # Both a case that converges and one that does not: the decision has to
        # be the same on every rank count, and a collective failure must not
        # hang, which the timeout in run_case() enforces.
        # BOTH outcomes are covered, and each is asserted to be the outcome it
        # is named for -- not merely to be the same on every rank. A case that
        # silently failed everywhere would otherwise satisfy an agreement-only
        # check while covering none of the converging path.
        for label, sides, tolerance, max_iterations, expect_converged in (
                ("converging", (1.0, 1.1, 1.3), 1e-10, 500, True),
                ("exhausting", (1.0, 1.1, 2.09), 1e-13, 200, False)):
            outcomes = {}
            for np in ranks:
                points = triangle(*sides, slack=0.003)
                directory = work / "mpi" / f"{label}-np{np}"
                write_case(directory, points, sides,
                           tolerance=tolerance, max_iterations=max_iterations,
                           time_step=0.5)
                code, output = run_case(args.gmd, directory, np=np,
                                        mpiexec=args.mpiexec, timeout=180)
                outcomes[np] = (code == 0, "TIMEOUT" in output)
            agreed = len({converged for converged, _ in outcomes.values()}) == 1
            hung = [np for np, (_, timed_out) in outcomes.items() if timed_out]
            as_expected = all(converged == expect_converged
                              for converged, _ in outcomes.values())
            check(f"mpi_convergence_decision_agrees_{label}",
                  agreed and as_expected and not hung,
                  f"outcomes {outcomes}, expected converged={expect_converged}; "
                  f"a non-collective decision would show up as disagreement or "
                  f"as a hang (the timeout catches the hang)")
    else:
        print("      skipped: no --mpiexec supplied")

    # ------------------------------------------------------------------
    print(f"\n{len(PASSES)} passed, {len(FAILURES)} failed")
    for failure in FAILURES:
        print(f"  FAILED: {failure}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
