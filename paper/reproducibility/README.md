# Reference checks for the graph–tensor manuscript

This directory accompanies `../md_graph_tensor_preprint.tex`. It contains a **new reference implementation and recorded run prepared for the September 17, 2026 revision**. The original draft's verification script, JSON, and environment metadata were unavailable. These files are not recovered original-run artifacts. The manuscript's numerical table and residuals have been replaced by the results actually produced here.

The checks use CPU double precision. The script does not import the existing GMD engine or `next/`, require CUDA, compile LaTeX, or generate a PDF. AI assistance was used in preparing this implementation and the manuscript revision; author review remains pending.

## Files

| File | Purpose |
| --- | --- |
| `verify_equivalence.py` | NumPy operator implementation, separately written scalar-loop formulas, finite differences, and consistency checks |
| `config.json` | Complete main-test coordinates, ordered indices, potential parameters, reference units, velocities, random-input specification, and tolerances |
| `requirements.txt` | Recorded NumPy version |
| `verification_results.json` | All inputs, energies, analytic and finite-difference force arrays, step sweeps, residuals, pass/fail results, and supplemental geometry examples |
| `environment.json` | Actual interpreter, NumPy, operating system, architecture, NumPy build/BLAS information, run timestamp, and source hashes |
| `table_rows.tex` | Generated rows for the numerical table; the main TeX keeps its own marked copy and has no new `\input` dependency |

## Run

The recorded environment is CPython 3.9.6, NumPy 2.0.2, macOS arm64, with Accelerate BLAS. `environment.json` gives the details. Using an existing environment with the recorded dependency is sufficient. Otherwise, from this directory, create an isolated environment using an interpreter compatible with the pinned NumPy version:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python verify_equivalence.py --output-dir /tmp/gmd-paper-reproduction
```

To run from the repository root using an environment that already provides NumPy:

```sh
python3 paper/reproducibility/verify_equivalence.py --output-dir /tmp/gmd-paper-reproduction
```

The output directory option keeps the checked-in reference artifacts unchanged. Without it, the three generated files are overwritten **in this directory**, regardless of the working directory. There is no network access during verification.

The command exits nonzero if a numerical check fails. The recorded run passes 56 checks. Values approaching the force-shift cutoff are recorded as diagnostic samples rather than counted as extra pass/fail tests.

The following optional command also verifies that the manuscript's marked table rows exactly match the newly computed, rounded rows:

```sh
python3 paper/reproducibility/verify_equivalence.py --output-dir /tmp/gmd-paper-reproduction --check-manuscript paper/md_graph_tensor_preprint.tex
```

Exact table matching is useful when synchronizing an edited manuscript on the same environment. A different platform may pass all numerical tolerances yet produce different last digits and fail this editorial check. Inspect the results and environment before updating a table; do not tune fixtures to reproduce old rounding errors. The script never edits the manuscript.

## What is checked

- Uncut LJ, harmonic bonds, harmonic angles, periodic angular dihedrals, and single- and multi-element constructed EAM functions on the five-particle fixture.
- Uncut nearest-image LJ on a separately specified four-particle orthorhombic fixture, strictly away from image ties. This is a local check of a minimum-image model, not an infinite periodic LJ lattice sum.
- Scalar-loop versus matrix/Gather energy and force evaluation for every model. The dihedral loop uses closed-form Cartesian angle derivatives; the operator implementation uses the appendix's reverse chain.
- Central finite differences of scalar-loop energies at steps `1e-4`, `1e-5`, `1e-6`, and `1e-7`. All arrays and errors are saved, including steps with larger errors. The manuscript uses `1e-6` for every model.
- Half/full-list energies and force assembly, including source-only assembly, Gather–assembly adjointness with a repeated slot, edge reversal/permutation, and atom relabeling.
- Net forces, nonperiodic torques, LJ virial versus `X.T @ F`, and nonperiodic/periodic virial strain derivatives. The periodic strain check deforms coordinates and cell together at fixed image labels.
- One velocity Verlet step with separate loop and operator force evaluations and fixed masses and initial velocities.
- A force-shifted LJ fixture with fixed exclusions and scaling, an excluded coincident pair, and cutoff values/one-sided samples.
- Equal radial energy on the two sides of an image tie and a concrete skew-cell rounding counterexample. The bounded image enumeration in that example is not a general closest-lattice-vector algorithm.

No neighbor-list builder, capacity overflow mechanism, long trajectory, thermostat, variable-cell dynamics, PME solver, constraint solver, GPU performance, or third-party MD package is tested. References in the manuscript do not imply an executed comparison with those packages.

## Error and independence conventions

`max_abs_force_error` is the infinity norm of the difference between the operator force and the central finite difference of the loop energy. `normalized_force_error` divides it by `max(1, max(abs(operator_force)))` in the configured reference force unit. This is a scale-normalized error, not an elementwise relative error.

At step `1e-6`, the acceptance tolerance is `1e-7` for this normalized error. Loop/operator forces use `1e-12 * max(1, max(abs(force)))`; energies use the analogous energy scale. Other tests store their individual tolerances alongside their residuals. The table is not a selection of the lowest error at each step size. Numerical tolerances are intended for these fixed, nondegenerate fixtures, not universal guarantees.

The two evaluation paths share physical definitions and input data, but do not call each other's energy, force, or assembly routines. They can still share a mistaken convention. Agreement is evidence of internal mathematical/implementation consistency, not validation against an external simulator or of material properties. Floating-point summation, transcendental functions, and BLAS choices can change final digits.

Values use fixed reference length, energy, and mass units. Time is measured in `length * sqrt(mass/energy)` and force in `energy/length`. Angles are in radians. The constructed EAM densities use a chosen density unit; exponential arguments use the numerical distance in the reference length unit. No mapping to a particular material or experimental unit system is claimed.

The sole random draw is a `(3, 3, 3)` standard-normal tensor from `numpy.random.default_rng(17)` with no previous draws. It is saved in the results along with its index array. Initial velocities are explicit, not random. Source and configuration SHA-256 hashes in both generated JSON files identify the run inputs; changing either requires a new recorded run.

## Manuscript revisions supported by this package

The revision clarifies the difference between local differentiation at a fixed image branch and the cutoff condition for equivalence to a periodic image sum. Smooth radial minimum-image energy is continuous at image ties, although force may jump. A cutoff is not necessary for local differentiation away from ties.

It also makes pair scaling and exclusions explicit, documents all formerly missing fixtures, states the limits of the numerical evidence, and replaces unsupported availability claims with the files in this directory. Primary references and JAX-MD's structure-mapping documentation have been added to distinguish the paper's contract-based organization from established physical and computational ingredients.
