# Constrained NVE: a free rigid water molecule

The first checked-in validation of constrained dynamics on the production path.
Everything is driven through the real `gmd` executable; nothing calls the
constraint solver directly.

## Fixture and physical intent

One rigid-water-like molecule in a 24 x 26 x 28 A cell:

| | |
|---|---|
| atoms | O (15.9994 amu), H, H (1.008 amu) — **unequal masses** |
| constraints | 2 x O-H at 0.9572 A, 1 x H-H at 1.51390065452732 A |
| H-H target | `2 r_OH sin(104.52 deg / 2)`, i.e. the H-O-H angle expressed as a distance |
| orientation | a generic rotation; no bond is parallel to a box axis |
| initial slack | every distance starts **0.25% above** its target |
| initial velocities | a rigid-body field plus a deliberate **non-tangent radial** component |

`water.ff` gives the two bonds and the angle **zero force constants**. They are
declared only so the special-pair map excludes the intramolecular non-bonded
pairs. The molecule is alone in the cell, so it feels **no force at all**.

That is the point of the fixture. A force-free rigid body has analytically known
motion, so most of this case is checked against mechanics rather than against a
previous GMD run:

* total energy exactly constant;
* linear momentum constant, and the centre of mass on a straight line;
* angular momentum about the centre of mass constant;
* the constraint virial cancels the **rotational** kinetic energy, so the
  reported pressure collapses to the ideal-gas value `M |v_com|^2 / 3V`.

The slack and the non-tangent velocity component exist so that the initial SHAKE
and RATTLE projections each have real work to do. "The distances match their
targets after projection" would otherwise be a tautology.

`water_wrapped.xyz` is the same molecule translated onto a periodic face, so its
atoms wrap to opposite sides of the cell. Every origin-independent quantity must
match the unwrapped run.

## What is checked

Residuals are **recomputed here** from 17-digit checkpoint state, not read back
from the solver's convergence flag:

* `|r_ij| - d_ij` — position residual, in A;
* `|r_ij . v_ij| / |r_ij|` — velocity tangency residual, the per-constraint row
  of `|Jv|`, in **A per internal time unit** (see "Units" below).

Also: the constrained degree-of-freedom count, per-frame log bounds, the
pressure-validity semantics of the initial frame, all nine constraint-virial
components, tensor symmetry, non-trivial off-diagonals, the kinetic-plus-virial
identity, wrapped/unwrapped equivalence, checkpoint/restart continuity, and
serial/MPI agreement at np = 1, 2, 4.

## Units

GMD stores velocities in its **internal time unit**, `A sqrt(amu/eV)`, which is
`10.180505717871194` fs — not in femtoseconds. The velocity tangency residual is
therefore in A per internal time unit.

Note that the `output.log` header labels that column `rattle_error[A/fs]`, which
is wrong by this factor of 10.18. The number itself is the internal-unit one.

## Tolerance derivation

Every bound in `tolerance.json` comes from a measurement. The pattern is: at
least ~10x above what the build actually produces, so it is not brittle, and
orders of magnitude below the failure it exists to catch.

| bound | value | observed | derivation |
|---|---|---|---|
| `position_residual_abs_a` | 2e-10 | 9.7e-11 | the configured solver tolerance is 1e-10; 2x absorbs the independent recomputation's own round-off |
| `velocity_residual_abs_internal` | 2e-10 | 9.9e-11 | as above |
| `initial_projection_xyz_abs_a` | 5e-6 | — | the `.xyz` carries six decimals, so ~1e-6 is the floor of the measurement. The fixture starts 2.4e-3 A off target, so a projection that did nothing would miss by 480x this bound |
| `log_residual_column_bound_a` | 1e-6 | 0.0 | the log's six fixed decimals can only *bound* the residual; the tight measurement is the checkpoint one |
| `energy_drift_per_atom_ps_abs` | 1e-8 | 0.0 | a force-free body has no potential energy, so this is integrator drift alone |
| `momentum_rel` | 1e-11 | 4.0e-15 | round-off only; exactly conserved by the pair-wise constraint impulse |
| `angular_momentum_rel` | 1e-10 | 6.5e-14 | as above |
| `com_trajectory_abs_a` | 1e-9 | 1.1e-13 | analytic straight-line prediction |
| `virial_symmetry_abs` | 1e-15 | 2.2e-19 | symmetric by construction; this catches a transposition |
| `virial_offdiagonal_min_abs` | 1e-5 | 4.0e-4 | a *floor*, not a ceiling: it fails if the fixture ever becomes axis-aligned enough to hide errors |
| `pressure_identity_rel` | 1e-12 | 0.0 | `P = (2K + tr W)/3V` recomputed from the checkpoint's own terms |
| `rotational_cancellation_abs` | 5e-7 | 3.6e-8 | second order in dt (1.2e-7 at dt=1.0, 3.6e-8 at dt=0.5, 1.8e-8 at dt=0.25, then a ~2e-8 floor). Removing the constraint virial makes this 7.7e-3 — five orders above the bound |
| `analytic_pressure_rel` | 5e-6 | 2.4e-7 | against `M |v_com|^2 / 3V`. Removing the constraint virial gives ~5e-2 |
| `restart_*` | 1e-12 | 0.0 | restart reproduces the continuous run bit for bit here |
| `mpi_*` | 1e-14 / 1e-12 | 0.0 | this fixture is force free, so there is no summation order to differ; the interacting NVT case has the looser, derived MPI bounds |

## Reproducing

```sh
python3 validation/constrained_nve_water/generate_fixture.py     # rewrite the fixtures
python3 validation/run_validation.py --case constrained_nve_water \
    --gmd <build>/gmd --gmd-validate <build>/gmd_validate --work-root <work> \
    [--mpiexec mpirun --mpi-ranks 1,2,4]
```

`reference.json` records the provenance: source commit, SHA-256 of every input,
build type, rank counts, and which numbers are analytic and which are GMD
regression results. The values under `recorded_observations` are **not**
compared against; they are there so a reviewer can see what this build produced.

## What is not claimed

No independent external engine was used. LAMMPS, GROMACS and OpenMM do not
expose a constraint virial in a form that could be compared component by
component, so no cross-engine reference exists for the quantity this case is
mainly about. The analytic checks listed above are the independent part; the
rest is regression coverage.
