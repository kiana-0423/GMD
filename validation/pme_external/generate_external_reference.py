#!/usr/bin/env python3
"""Generate the independent external PME reference for validation/pme_external.

This script runs OpenMM and LAMMPS. It does NOT run at test time: the normal CI
comparison reads the JSON this produces. Regenerating is a deliberate, manual
step.

It refuses to run if an engine is missing. It never falls back to GMD output to
stand in for an external result -- that failure mode is exactly what this whole
directory exists to rule out, so a missing engine is a hard error.

Nothing is written into the repository. All scratch files go to --work-dir, and
only --output (which the caller points at a path of its choosing) is a result.

Reproduce with:

    python3 validation/pme_external/generate_external_reference.py \
        --work-dir <work> \
        --gmd-validate <build>/gmd_validate \
        --openmm-python <python-with-openmm> \
        --lammps <lammps-binary> \
        --output <work>/reference_external.json

run from <repo>. Copy <work>/reference_external.json over
validation/pme_external/reference_external.json to update the checked-in
reference.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import pathlib
import platform
import shutil
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
CASE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(REPO / "validation"))

from analytic_references import ewald as analytic_ewald  # noqa: E402

# --- the settings every engine is pinned to ---------------------------------
ALPHA = 0.35          # Ewald splitting, 1/Angstrom
CUTOFF = 8.0          # real-space cutoff, Angstrom
GRIDS = [16, 32, 64, 128]
# LAMMPS PPPM is swept only up to this grid. `kspace_modify mesh N N N` with any
# N > 64 aborts this LAMMPS build with SIGSEGV before it prints anything, on this
# fixture, with or without gewald and order pinned, and with LAMMPS' own
# auto-selected splitting as well. It is a limit of the external engine, not of
# the comparison: GMD and OpenMM are still swept to 128, and the reference grid
# is 64, so LAMMPS covers the settings the reference is actually taken at.
# Reproduce with `kspace_style pppm 1.0e-10` + `kspace_modify mesh 72 72 72`.
LAMMPS_MAX_GRID = 64
GMD_ORDERS = [4, 6]
REFERENCE_GRID = 64   # the grid the checked-in CI reference is taken at
REFERENCE_GMD_ORDER = 6
ANALYTIC_KMAX = 24    # exp(-k^2/4 alpha^2) < 1e-20 on every axis of this box

# GMD's Coulomb constant, from gmd::kCoulombConstant in
# include/gmd/core/physical_constants.hpp. Every engine's constant is MEASURED
# below rather than taken from its documentation, because the whole comparison
# is scaled by it.
KE_GMD = 14.3996454686836


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fixture(path: pathlib.Path) -> dict:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    count = int(lines[0].split()[0])
    box = [float(value) for value in lines[1].split()[:3]]
    elements, positions, masses, charges = [], [], [], []
    for line in lines[2 : 2 + count]:
        fields = line.split()
        elements.append(fields[0])
        positions.append([float(fields[1]), float(fields[2]), float(fields[3])])
        masses.append(float(fields[4]))
        charges.append(float(fields[5]))
    return {"box": box, "elements": elements, "positions": positions,
            "masses": masses, "charges": charges}


# ---------------------------------------------------------------------------
# LAMMPS
# ---------------------------------------------------------------------------

def write_lammps_data(fixture: dict, path: pathlib.Path) -> None:
    box = fixture["box"]
    elements = fixture["elements"]
    kinds = sorted(set(elements))
    type_of = {element: index + 1 for index, element in enumerate(kinds)}
    mass_of = {element: fixture["masses"][elements.index(element)] for element in kinds}

    out = ["pme_external fixture", "",
           f"{len(elements)} atoms", f"{len(kinds)} atom types", ""]
    for axis, low_high in zip("xyz", box):
        out.append(f"0.0 {low_high!r} {axis}lo {axis}hi")
    out += ["", "Masses", ""]
    for element in kinds:
        out.append(f"{type_of[element]} {mass_of[element]!r}")
    out += ["", "Atoms # charge", ""]
    for index, element in enumerate(elements):
        position = fixture["positions"][index]
        # Image flags rather than raw out-of-box coordinates, so that the
        # wrapped fixture reaches LAMMPS as the same periodic images the file
        # names instead of images LAMMPS picks for itself.
        images = [int(math.floor(position[d] / box[d])) for d in range(3)]
        folded = [position[d] - images[d] * box[d] for d in range(3)]
        out.append(
            f"{index + 1} {type_of[element]} {fixture['charges'][index]!r} "
            f"{folded[0]!r} {folded[1]!r} {folded[2]!r} "
            f"{images[0]} {images[1]} {images[2]}")
    path.write_text("\n".join(out) + "\n")


def parse_lammps_thermo(log_text: str) -> dict:
    lines = log_text.splitlines()
    for index, line in enumerate(lines):
        fields = line.split()
        if fields[:1] == ["Step"]:
            values = lines[index + 1].split()
            if len(values) != len(fields):
                continue
            return {name: float(value) for name, value in zip(fields, values)}
    raise RuntimeError("no thermo row found in LAMMPS output")


def parse_lammps_dump(path: pathlib.Path, count: int) -> list[list[float]]:
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("ITEM: ATOMS")) + 1
    forces = [None] * count
    for line in lines[start : start + count]:
        fields = line.split()
        forces[int(fields[0]) - 1] = [float(fields[1]), float(fields[2]), float(fields[3])]
    if any(force is None for force in forces):
        raise RuntimeError("LAMMPS dump did not cover every atom id")
    return forces


def run_lammps(binary: str, script: pathlib.Path, work: pathlib.Path,
               variables: dict, tag: str) -> tuple[dict, list[list[float]], str]:
    force_file = work / f"forces_{tag}.dump"
    if force_file.exists():
        force_file.unlink()
    command = [binary, "-in", str(script), "-log", "none", "-screen", str(work / f"{tag}.log")]
    for name, value in variables.items():
        command += ["-var", name, str(value)]
    command += ["-var", "forcefile", str(force_file)]
    completed = subprocess.run(command, cwd=work, text=True, capture_output=True)
    log_text = (work / f"{tag}.log").read_text() if (work / f"{tag}.log").exists() else ""
    if completed.returncode != 0:
        raise RuntimeError(
            f"LAMMPS failed for {tag}:\n{completed.stdout}\n{completed.stderr}\n{log_text}")
    return parse_lammps_thermo(log_text), force_file, log_text


def measure_lammps_coulomb_constant(binary: str, work: pathlib.Path) -> tuple[float, float]:
    """Measure k_e and the pressure unit factor from LAMMPS itself.

    Two opposite unit charges 2 Angstrom apart under plain coul/cut in a large
    box. Analytically E = k_e q1 q2 / r and the configurational virial of a 1/r
    pair potential is W_xx = E with W_yy = W_zz = 0, so both the energy constant
    and LAMMPS' pressure-to-internal factor fall out of one run. Measuring beats
    quoting the documentation: a wrong constant here would rescale the entire
    comparison by a factor that looks exactly like a physics disagreement.
    """
    separation, half_box = 2.0, 100.0
    data = work / "ke_probe.data"
    data.write_text(
        "ke probe\n\n2 atoms\n1 atom types\n\n"
        f"0.0 {half_box} xlo xhi\n0.0 {half_box} ylo yhi\n0.0 {half_box} zlo zhi\n\n"
        "Masses\n\n1 1.0\n\nAtoms # charge\n\n"
        f"1 1 1.0 10.0 10.0 10.0 0 0 0\n"
        f"2 1 -1.0 {10.0 + separation} 10.0 10.0 0 0 0\n")
    script = work / "in.ke_probe"
    script.write_text(
        "units metal\natom_style charge\nboundary p p p\n"
        f"read_data {data}\n"
        "pair_style coul/cut 20.0\npair_coeff * *\n"
        "compute virial all pressure NULL virial\n"
        "thermo_style custom step pe c_virial[1]\n"
        "thermo_modify format float %.17g\nrun 0\n")
    completed = subprocess.run(
        [binary, "-in", str(script), "-log", "none", "-screen", str(work / "ke_probe.log")],
        cwd=work, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(f"LAMMPS k_e probe failed:\n{completed.stderr}")
    row = parse_lammps_thermo((work / "ke_probe.log").read_text())
    energy = row["PotEng"]
    coulomb_constant = -energy * separation           # E = -k_e/r for q = +1,-1
    volume = half_box ** 3
    # P_xx = W_xx / V * nktv2p and W_xx = E for this pair, so nktv2p follows.
    pressure_factor = row["c_virial[1]"] * volume / energy
    return coulomb_constant, pressure_factor


# ---------------------------------------------------------------------------
# OpenMM
# ---------------------------------------------------------------------------

def run_openmm(python: str, xyz: pathlib.Path, grid: int) -> dict:
    completed = subprocess.run(
        [python, str(CASE / "openmm" / "run_openmm.py"),
         "--xyz", str(xyz), "--alpha", str(ALPHA), "--cutoff", str(CUTOFF),
         "--grid", str(grid), str(grid), str(grid)],
        text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(f"OpenMM failed:\n{completed.stdout}\n{completed.stderr}")
    return json.loads(completed.stdout)


def measure_openmm_coulomb_constant(python: str, work: pathlib.Path) -> float:
    """Measure OpenMM's ONE_4PI_EPS0 the same way, in kJ/mol nm /e^2."""
    probe = work / "ke_probe_openmm.py"
    probe.write_text(
        "import openmm as mm, openmm.unit as u\n"
        "from openmm import Vec3\n"
        "s = mm.System(); s.addParticle(1.0); s.addParticle(1.0)\n"
        "f = mm.NonbondedForce(); f.setNonbondedMethod(mm.NonbondedForce.NoCutoff)\n"
        "f.addParticle(1.0, 1.0, 0.0); f.addParticle(-1.0, 1.0, 0.0); s.addForce(f)\n"
        "c = mm.Context(s, mm.VerletIntegrator(1.0),\n"
        "               mm.Platform.getPlatformByName('Reference'))\n"
        "c.setPositions([Vec3(0,0,0), Vec3(0.2,0,0)]*u.nanometer)\n"
        "e = c.getState(getEnergy=True).getPotentialEnergy()\n"
        "print(repr(-e.value_in_unit(u.kilojoule_per_mole) * 0.2))\n")
    completed = subprocess.run([python, str(probe)], text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(f"OpenMM k_e probe failed:\n{completed.stderr}")
    return float(completed.stdout.strip().splitlines()[-1])


# ---------------------------------------------------------------------------
# GMD
# ---------------------------------------------------------------------------

def run_gmd(binary: str, work: pathlib.Path, xyz: str, run_file: str, tag: str) -> dict:
    output = work / f"gmd_{tag}.json"
    completed = subprocess.run(
        [binary, xyz, run_file, "--json", str(output)],
        cwd=work, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(f"gmd_validate failed for {tag}:\n{completed.stderr}")
    return json.loads(output.read_text())


# ---------------------------------------------------------------------------
# comparison helpers
# ---------------------------------------------------------------------------

def force_errors(actual: list[list[float]], expected: list[list[float]]) -> dict:
    worst_abs, worst_atom, worst_axis = 0.0, -1, -1
    sum_sq, reference_sq = 0.0, 0.0
    for atom, (row_a, row_b) in enumerate(zip(actual, expected)):
        for axis in range(3):
            delta = abs(row_a[axis] - row_b[axis])
            sum_sq += delta * delta
            reference_sq += row_b[axis] * row_b[axis]
            if delta > worst_abs:
                worst_abs, worst_atom, worst_axis = delta, atom, axis
    count = 3 * len(expected)
    largest = max(abs(value) for row in expected for value in row)
    return {
        "max_abs": worst_abs,
        "max_abs_atom": worst_atom,
        "max_abs_component": "xyz"[worst_axis] if worst_axis >= 0 else None,
        "max_abs_reference_value": expected[worst_atom][worst_axis] if worst_atom >= 0 else None,
        # Relative to the largest force component present, not to the local
        # component: dividing by a component that happens to sit near a zero
        # crossing manufactures a huge "relative error" that means nothing.
        "max_rel_to_largest_component": worst_abs / largest if largest > 0 else 0.0,
        "rms_abs": math.sqrt(sum_sq / count),
        "rms_rel": math.sqrt(sum_sq / reference_sq) if reference_sq > 0 else 0.0,
        "largest_reference_component": largest,
    }


def scale_forces(forces: list[list[float]], factor: float) -> list[list[float]]:
    return [[value * factor for value in row] for row in forces]


# LAMMPS reports the virial as six independent components in this order; GMD
# reports the full row-major 3x3, of which these are the same six entries.
VIRIAL_ORDER = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def virial_comparison(gmd_tensor: list[float], lammps_six: list[float]) -> dict:
    """Compare all nine GMD components against LAMMPS' six independent ones.

    LAMMPS reports six because the tensor is symmetric. That is not assumed
    here: GMD's own asymmetry is measured and reported, so a GMD tensor that
    was quietly non-symmetric could not pass by being read through a symmetric
    six-component view.
    """
    asymmetry = max(abs(gmd_tensor[a * 3 + b] - gmd_tensor[b * 3 + a])
                    for a in range(3) for b in range(3))
    rows, worst, worst_name = [], 0.0, None
    for (a, b), lammps_value in zip(VIRIAL_ORDER, lammps_six):
        name = f"W_{'xyz'[a]}{'xyz'[b]}"
        delta = gmd_tensor[a * 3 + b] - lammps_value
        rows.append({"component": name, "gmd": gmd_tensor[a * 3 + b],
                     "lammps": lammps_value, "abs_error": abs(delta)})
        if abs(delta) > worst:
            worst, worst_name = abs(delta), name
    largest = max(abs(value) for value in lammps_six)
    return {
        "components": rows,
        "max_abs_error": worst,
        "max_abs_error_component": worst_name,
        "max_rel_to_largest_component": worst / largest if largest > 0 else 0.0,
        "gmd_max_asymmetry": asymmetry,
        "note": "all nine GMD components are covered: the six independent ones "
                "are compared against LAMMPS and the remaining three are the "
                "transposes, checked through gmd_max_asymmetry",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--gmd-validate", required=True)
    parser.add_argument("--openmm-python", required=True,
                        help="interpreter that can import openmm")
    parser.add_argument("--lammps", default="lmp_serial")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    work = pathlib.Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    if REPO in work.parents or work == REPO:
        raise SystemExit(
            f"--work-dir must be outside the repository; got {work} inside {REPO}")

    lammps = shutil.which(args.lammps) or args.lammps
    if not pathlib.Path(lammps).exists():
        raise SystemExit(
            f"LAMMPS binary not found: {args.lammps}. This script will not substitute "
            "GMD output for an external reference; install LAMMPS or pass --lammps.")
    if not pathlib.Path(args.gmd_validate).exists():
        raise SystemExit(f"gmd_validate not found: {args.gmd_validate}")
    probe = subprocess.run([args.openmm_python, "-c", "import openmm"],
                           text=True, capture_output=True)
    if probe.returncode != 0:
        raise SystemExit(
            f"{args.openmm_python} cannot import openmm:\n{probe.stderr}\n"
            "This script will not substitute GMD output for an external reference.")

    fixtures = {name: read_fixture(CASE / f"{name}.xyz")
                for name in ["charges", "charges_wrapped", "charges_translated"]}
    base = fixtures["charges"]
    count = len(base["charges"])

    for name in fixtures:
        shutil.copy2(CASE / f"{name}.xyz", work / f"{name}.xyz")
    for run_file in CASE.glob("*.run"):
        shutil.copy2(run_file, work / run_file.name)

    # --- engine constants, measured -----------------------------------------
    ke_lammps, lammps_pressure_factor = measure_lammps_coulomb_constant(lammps, work)
    ke_openmm_native = measure_openmm_coulomb_constant(args.openmm_python, work)
    # Coulomb energy is exactly proportional to k_e in every term (real,
    # reciprocal, self and background each carry one factor), so converting an
    # engine's result into GMD's constant is an exact rescaling and not an
    # approximation. Note that for OpenMM the eV <-> kJ/mol convention cancels:
    # the factor is k_e^GMD / (10 * ONE_4PI_EPS0), which contains no energy-unit
    # constant at all.
    #
    # This rescaling used to be compensating GMD's own truncated constant, which
    # was 3.16e-06 below CODATA and therefore the largest single discrepancy in
    # the whole comparison. GMD now uses the CODATA 2022 value, so what remains
    # is the engines' own rounding of the same physical constant -- LAMMPS to six
    # decimals, OpenMM to whatever CODATA release it was built against -- which is
    # genuinely engine-specific and has to stay. It is now a 3e-08 correction
    # rather than a 3e-06 one.
    openmm_energy_to_ev = KE_GMD / (10.0 * ke_openmm_native)
    openmm_force_to_ev_per_angstrom = KE_GMD / (100.0 * ke_openmm_native)
    lammps_scale = KE_GMD / ke_lammps

    # --- exact analytic Ewald, in GMD's own constant -------------------------
    analytic_energy, analytic_forces = analytic_ewald(
        base["positions"], base["charges"], base["box"],
        alpha=ALPHA, kmax=ANALYTIC_KMAX, cutoff=CUTOFF)

    # --- GMD -----------------------------------------------------------------
    gmd_runs = {}
    for grid in GRIDS:
        for order in GMD_ORDERS:
            tag = f"g{grid}_o{order}"
            gmd_runs[tag] = run_gmd(args.gmd_validate, work, "charges.xyz",
                                    f"pme_{tag}.run", tag)
    gmd_runs["ewald"] = run_gmd(args.gmd_validate, work, "charges.xyz",
                                "ewald_converged.run", "ewald")
    reference_tag = f"g{REFERENCE_GRID}_o{REFERENCE_GMD_ORDER}"
    for variant in ["charges_wrapped", "charges_translated"]:
        gmd_runs[f"{reference_tag}_{variant}"] = run_gmd(
            args.gmd_validate, work, f"{variant}.xyz",
            f"pme_{reference_tag}.run", f"{reference_tag}_{variant}")

    # --- OpenMM ---------------------------------------------------------------
    openmm_runs = {}
    for grid in GRIDS:
        openmm_runs[f"g{grid}"] = run_openmm(args.openmm_python, work / "charges.xyz", grid)
    for variant in ["charges_wrapped", "charges_translated"]:
        openmm_runs[f"g{REFERENCE_GRID}_{variant}"] = run_openmm(
            args.openmm_python, work / f"{variant}.xyz", REFERENCE_GRID)

    # --- LAMMPS ---------------------------------------------------------------
    lammps_runs = {}
    lammps_version = None
    for name, fixture in fixtures.items():
        write_lammps_data(fixture, work / f"{name}.data")
    for grid in [value for value in GRIDS if value <= LAMMPS_MAX_GRID]:
        for order in GMD_ORDERS:
            tag = f"g{grid}_o{order}"
            row, dump, log = run_lammps(
                lammps, CASE / "lammps" / "in.pppm", work,
                {"datafile": work / "charges.data", "rcut": CUTOFF, "alpha": ALPHA,
                 "nx": grid, "ny": grid, "nz": grid, "order": order},
                f"pppm_{tag}")
            lammps_runs[f"pppm_{tag}"] = (row, parse_lammps_dump(dump, count))
            lammps_version = lammps_version or log.splitlines()[0].strip()
    row, dump, log = run_lammps(
        lammps, CASE / "lammps" / "in.ewald", work,
        {"datafile": work / "charges.data", "rcut": CUTOFF, "alpha": ALPHA}, "ewald")
    lammps_runs["ewald"] = (row, parse_lammps_dump(dump, count))
    for variant in ["charges_wrapped", "charges_translated"]:
        tag = f"pppm_g{REFERENCE_GRID}_o{REFERENCE_GMD_ORDER}_{variant}"
        row, dump, _ = run_lammps(
            lammps, CASE / "lammps" / "in.pppm", work,
            {"datafile": work / f"{variant}.data", "rcut": CUTOFF, "alpha": ALPHA,
             "nx": REFERENCE_GRID, "ny": REFERENCE_GRID, "nz": REFERENCE_GRID,
             "order": REFERENCE_GMD_ORDER}, tag)
        lammps_runs[tag] = (row, parse_lammps_dump(dump, count))

    def lammps_energy(tag: str) -> float:
        return lammps_runs[tag][0]["PotEng"] * lammps_scale

    def lammps_forces(tag: str) -> list[list[float]]:
        return scale_forces(lammps_runs[tag][1], lammps_scale)

    def lammps_virial(tag: str) -> list[float]:
        row = lammps_runs[tag][0]
        volume = base["box"][0] * base["box"][1] * base["box"][2]
        # P_ab = W_ab / V * nktv2p, inverted with the factor measured above.
        return [row[f"c_virial[{i}]"] * volume / lammps_pressure_factor * lammps_scale
                for i in range(1, 7)]

    def openmm_energy(tag: str) -> float:
        return openmm_runs[tag]["energy_kj_per_mol"] * openmm_energy_to_ev

    def openmm_forces(tag: str) -> list[list[float]]:
        return scale_forces(openmm_runs[tag]["forces_kj_per_mol_per_nm"],
                            openmm_force_to_ev_per_angstrom)

    # --- convergence toward the exact Ewald limit -----------------------------
    convergence = []
    for grid in GRIDS:
        row = {"grid": [grid, grid, grid]}
        for order in GMD_ORDERS:
            result = gmd_runs[f"g{grid}_o{order}"]
            row[f"gmd_order{order}"] = {
                "energy_ev": result["energy"]["total"],
                "energy_abs_error_vs_analytic_ewald":
                    abs(result["energy"]["total"] - analytic_energy),
                "force": force_errors(result["force"]["atoms"], analytic_forces),
            }
            if grid <= LAMMPS_MAX_GRID:
                row[f"lammps_pppm_order{order}"] = {
                    "energy_ev": lammps_energy(f"pppm_g{grid}_o{order}"),
                    "energy_abs_error_vs_analytic_ewald":
                        abs(lammps_energy(f"pppm_g{grid}_o{order}") - analytic_energy),
                    "force": force_errors(lammps_forces(f"pppm_g{grid}_o{order}"),
                                          analytic_forces),
                }
            else:
                row[f"lammps_pppm_order{order}"] = {
                    "skipped": f"LAMMPS PPPM aborts for mesh > {LAMMPS_MAX_GRID}"
                }
        row["openmm_order5"] = {
            "energy_ev": openmm_energy(f"g{grid}"),
            "energy_abs_error_vs_analytic_ewald":
                abs(openmm_energy(f"g{grid}") - analytic_energy),
            "force": force_errors(openmm_forces(f"g{grid}"), analytic_forces),
        }
        convergence.append(row)

    reference_forces = openmm_forces(f"g{REFERENCE_GRID}")
    reference_energy = openmm_energy(f"g{REFERENCE_GRID}")
    gmd_reference = gmd_runs[reference_tag]

    payload = {
        "source": "external_reference_openmm_pme",
        "schema": "gmd/pme_external/1",
        "reference": {
            "engine": "OpenMM",
            "version": openmm_runs[f"g{REFERENCE_GRID}"]["version"],
            "git_revision": openmm_runs[f"g{REFERENCE_GRID}"]["git_revision"],
            "platform": openmm_runs[f"g{REFERENCE_GRID}"]["platform"],
            "grid": [REFERENCE_GRID] * 3,
            "bspline_order": 5,
            "alpha_inv_angstrom": ALPHA,
            "cutoff_angstrom": CUTOFF,
            "energy_ev": reference_energy,
            "forces_ev_per_angstrom": reference_forces,
            "energy_native_kj_per_mol":
                openmm_runs[f"g{REFERENCE_GRID}"]["energy_kj_per_mol"],
            "wrapped": {
                "energy_ev": openmm_energy(f"g{REFERENCE_GRID}_charges_wrapped"),
                "forces_ev_per_angstrom":
                    openmm_forces(f"g{REFERENCE_GRID}_charges_wrapped"),
            },
            "translated": {
                "energy_ev": openmm_energy(f"g{REFERENCE_GRID}_charges_translated"),
                "forces_ev_per_angstrom":
                    openmm_forces(f"g{REFERENCE_GRID}_charges_translated"),
            },
        },
        "cross_reference": {
            "lammps_pppm": {
                "engine": "LAMMPS",
                "version": lammps_version,
                "grid": [REFERENCE_GRID] * 3,
                "order": REFERENCE_GMD_ORDER,
                "energy_ev": lammps_energy(f"pppm_{reference_tag}"),
                "forces_ev_per_angstrom": lammps_forces(f"pppm_{reference_tag}"),
                "virial_ev_xx_yy_zz_xy_xz_yz": lammps_virial(f"pppm_{reference_tag}"),
            },
            "lammps_ewald": {
                "engine": "LAMMPS",
                "version": lammps_version,
                "energy_ev": lammps_energy("ewald"),
                "forces_ev_per_angstrom": lammps_forces("ewald"),
                "virial_ev_xx_yy_zz_xy_xz_yz": lammps_virial("ewald"),
            },
        },
        "analytic_ewald": {
            "source": "validation/analytic_references.py ewald(), kmax "
                      f"{ANALYTIC_KMAX}, same alpha and same real-space cutoff",
            "energy_ev": analytic_energy,
            "forces_ev_per_angstrom": analytic_forces,
        },
        "gmd": {
            tag: {"energy_ev": result["energy"]["total"],
                  "forces_ev_per_angstrom": result["force"]["atoms"],
                  "virial_ev_row_major": result["virial"]["tensor"],
                  "virial_valid": result["virial"]["valid"]}
            for tag, result in gmd_runs.items()
        },
        "virial": {
            "note":
                "OpenMM exposes no virial at all, so the virial is validated "
                "against LAMMPS only. GMD's convention is W_ab = sum_i r_ia F_ib "
                "row major; LAMMPS' compute pressure with a NULL temperature is "
                "the same configurational sum, reported as pressure and converted "
                "back with the volume and the measured nktv2p.",
            "exact_ewald_vs_lammps_ewald": virial_comparison(
                gmd_runs["ewald"]["virial"]["tensor"],
                lammps_virial("ewald")),
            "pme_vs_lammps_pppm_at_reference_settings": virial_comparison(
                gmd_reference["virial"]["tensor"],
                lammps_virial(f"pppm_{reference_tag}")),
            "trace_versus_energy": {
                "gmd_ewald_trace": sum(gmd_runs["ewald"]["virial"]["tensor"][i]
                                       for i in (0, 4, 8)),
                "gmd_ewald_energy": gmd_runs["ewald"]["energy"]["total"],
                "note":
                    "For an untruncated 1/r lattice sum the trace equals the "
                    "energy exactly. These differ at the 1e-4 level because the "
                    "real-space sum is truncated at the 8 Angstrom cutoff, which "
                    "breaks the homogeneity the identity relies on. The identity "
                    "is therefore NOT used as a check; it is recorded so the "
                    "discrepancy is not later mistaken for a virial defect.",
            },
        },
        "convergence": convergence,
        "agreement_at_reference_settings": {
            "gmd_run": reference_tag,
            "gmd_energy_ev": gmd_reference["energy"]["total"],
            "openmm_energy_ev": reference_energy,
            "energy_abs_error": abs(gmd_reference["energy"]["total"] - reference_energy),
            "energy_rel_error":
                abs(gmd_reference["energy"]["total"] - reference_energy) / abs(reference_energy),
            "force": force_errors(gmd_reference["force"]["atoms"], reference_forces),
            "gmd_energy_error_vs_analytic_ewald":
                abs(gmd_reference["energy"]["total"] - analytic_energy),
            "openmm_energy_error_vs_analytic_ewald":
                abs(reference_energy - analytic_energy),
        },
        "provenance": {
            "generated_utc": datetime.datetime.now(datetime.timezone.utc)
                .replace(microsecond=0).isoformat(),
            "source_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
                capture_output=True).stdout.strip(),
            "platform": f"{platform.system()} {platform.release()} {platform.machine()}",
            "engines": {
                "openmm": {
                    "version": openmm_runs[f"g{REFERENCE_GRID}"]["version"],
                    "git_revision": openmm_runs[f"g{REFERENCE_GRID}"]["git_revision"],
                    "platform": openmm_runs[f"g{REFERENCE_GRID}"]["platform"],
                    "precision": "double (Reference platform)",
                    "installation": "pip install openmm (PyPI wheel) into a virtual "
                                    "environment outside this repository",
                    "measured_one_4pi_eps0_kj_per_mol_nm": ke_openmm_native,
                    "reports_virial": False,
                },
                "lammps": {
                    "version": lammps_version,
                    "installation": "Homebrew formula 'lammps'",
                    "precision": "double",
                    "units": "metal (eV, Angstrom, e)",
                    "measured_coulomb_constant_ev_angstrom": ke_lammps,
                    "measured_pressure_unit_factor_nktv2p": lammps_pressure_factor,
                    "reports_virial": True,
                    "reports_virial_note":
                        "six independent components (xx yy zz xy xz yz) via "
                        "compute pressure with NULL temperature; the tensor is "
                        "symmetric so this is the full tensor",
                    "max_usable_pppm_grid": LAMMPS_MAX_GRID,
                    "max_usable_pppm_grid_note":
                        "kspace_modify mesh N N N aborts this build with SIGSEGV "
                        "for every N > 64 on this fixture, with pinned or "
                        "auto-selected gewald and order alike, so the LAMMPS PPPM "
                        "convergence sweep stops at 64 while GMD and OpenMM "
                        "continue to 128.",
                },
            },
            "unit_conversion": {
                "gmd_coulomb_constant_ev_angstrom": KE_GMD,
                "openmm_energy_kj_per_mol_to_ev": openmm_energy_to_ev,
                "openmm_force_kj_per_mol_nm_to_ev_per_angstrom":
                    openmm_force_to_ev_per_angstrom,
                "lammps_scale": lammps_scale,
                "note":
                    "Each engine's result is rescaled by k_e^GMD / k_e^engine. "
                    "Coulomb energy is exactly proportional to k_e in every term, "
                    "so this is exact. GMD uses the CODATA 2022 value "
                    f"{KE_GMD}; the residual against LAMMPS is "
                    f"{abs(KE_GMD / ke_lammps - 1.0):.3e} relative and against "
                    f"OpenMM {abs(KE_GMD / (ke_openmm_native * 10.0 / (6.02214076e23 * 1.602176634e-19 / 1000.0)) - 1.0):.3e}, "
                    "which is those engines' own rounding of the same physical "
                    "constant rather than a disagreement about physics. Before "
                    "GMD's constant was corrected this factor was 3.16e-06 and "
                    "was compensating GMD's truncation, which would otherwise "
                    "have exceeded the grid-128 mesh error and been misread as a "
                    "convergence floor.",
            },
            "fixture_sha256": {
                path.name: sha256(path)
                for path in sorted(CASE.glob("*.xyz"))
            },
            "input_sha256": {
                f"{path.parent.name}/{path.name}" if path.parent != CASE else path.name:
                    sha256(path)
                for path in sorted(
                    list(CASE.glob("*.run")) + list((CASE / "lammps").glob("in.*"))
                    + [CASE / "openmm" / "run_openmm.py",
                       CASE / "generate_external_reference.py"])
            },
            "command": (
                "python3 validation/pme_external/generate_external_reference.py "
                "--work-dir <work> --gmd-validate <build>/gmd_validate "
                "--openmm-python <python-with-openmm> --lammps <lammps-binary> "
                "--output <work>/reference_external.json"),
            "working_directory": "<repo>",
            "generated_file": "<work>/reference_external.json",
            "gmd_build": "Release (-DCMAKE_BUILD_TYPE=Release)",
            "not_gmd_generated":
                "The reference and cross_reference blocks come from OpenMM and "
                "LAMMPS. The gmd block is GMD's own output, recorded alongside so "
                "the comparison is auditable; it is never the reference.",
        },
    }

    output = pathlib.Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    agreement = payload["agreement_at_reference_settings"]
    print(f"wrote {output}")
    print(f"  GMD {reference_tag} vs OpenMM grid {REFERENCE_GRID}:")
    print(f"    energy  abs {agreement['energy_abs_error']:.4e} eV"
          f"  rel {agreement['energy_rel_error']:.4e}")
    print(f"    force   max abs {agreement['force']['max_abs']:.4e} eV/A"
          f"  at atom {agreement['force']['max_abs_atom']}"
          f" component {agreement['force']['max_abs_component']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
