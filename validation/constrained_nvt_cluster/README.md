# Constrained NVT: four interacting rigid molecules under Nose-Hoover

Constrained dynamics with a thermostat, real intermolecular forces, and MPI rank
boundaries that cut through molecules. Driven entirely through the real `gmd`
executable.

## Fixture and physical intent

Four rigid-water-like molecules in a 20 x 22 x 24 A cell — 12 atoms, 12
constraints, 4 connected components of the constraint graph.

| | |
|---|---|
| per molecule | 2 x O-H at 0.9572 A, 1 x H-H at 1.51390065452732 A |
| masses | O 15.9994 amu, H 1.008 amu |
| initial slack | +0.18%, -0.22%, +0.31%, -0.14% — **different per molecule**, and signed both ways |
| orientation | each molecule has its own generic rotation; none are related by symmetry |
| forces | **intermolecular Lennard-Jones only** (bonded terms have k = 0 and exist solely for exclusions) |

### Why this case exists alongside the NVE one

In `constrained_nve_water` the molecule is alone and force free, so the provider
virial is exactly zero and the constraint virial is the whole of it. Here the
molecules interact, so the provider virial is genuinely non-zero and the
completed-step virial is a **sum of two independently produced tensors at the
same time level**. The check that the total equals provider + constraint is what
a constraint term taken from the wrong step would fail.

### MPI layout, chosen deliberately

GMD splits the cell over a balanced process grid: `(2,1,1)` at np=2 and
`(2,2,1)` at np=4.

* molecule **C** straddles `x = Lx/2 = 10.0` — cross-rank at np=2 and np=4;
* molecule **D** straddles `y = Ly/2 = 11.0` — cross-rank at np=4;
* the `x > Lx/2, y > Ly/2` quadrant is **empty**, so at np=4 one rank owns no
  atoms and the collective paths run with an empty local domain.

`generate_fixture.py` prints all of this when run, so the properties can be
re-checked rather than taken on trust.

## Degrees of freedom

`3N - constraints - 3 = 36 - 12 - 3 = **21**`, with COM removal on.

This is load bearing. The unconstrained count would be `3N - 3 = 33`; a run that
used it would report a temperature low by `33/21` and, because Nose-Hoover's
thermostat mass is `Q = dof kB T tau^2`, would thermostat to the wrong kinetic
energy. The case asserts the reported value directly and asserts that the
initial temperature lands exactly on 300 K, which only happens when the
projection-then-rescale order and the constrained DOF are both right.

## What is checked

Residuals are recomputed here from 17-digit checkpoint state. Beyond the NVE
case, this one adds: the authoritative DOF, the initial temperature after
constraint-aware initialization, the temperature mean and bounded excursion,
Nose-Hoover state finiteness, a non-zero provider virial, the
provider-plus-constraint sum, the absence of any energy or temperature jump
across a restart, and serial/MPI agreement at np = 1, 2, 4.

## Tolerance derivation

| bound | value | observed | derivation |
|---|---|---|---|
| `initial_temperature_abs_k` | 1e-9 | 0.0 | the initializer rescales to the target exactly, using the constrained DOF. Using 3N-3 would land at ~471 K |
| `temperature_mean_abs_k` | 1.0 | 0.031 | dynamics coverage over 400 steps |
| `temperature_max_excursion_k` | 5.0 | 0.20 | bounded fluctuation, not an ensemble width |
| `position_residual_abs_a` | 2e-10 | 1.0e-10 | solver tolerance 1e-10, 2x for recomputation round-off |
| `velocity_residual_abs_internal` | 2e-10 | 9.7e-11 | as above |
| `provider_virial_min_abs` | 1e-6 | 5.1e-3 | a *floor*: it fails if the molecules ever stop interacting, which would silently turn this into the NVE case |
| `virial_sum_abs` | 1e-14 | 0.0 | total = provider + constraint, exactly |
| `pressure_identity_rel` | 1e-12 | 0.0 | `P = (2K + tr W)/3V` |
| `restart_*` | 1e-12 / 1e-6 K / 1e-9 eV | 0.0 | restart reproduces the continuous run bit for bit |
| `mpi_virial_abs` | 1e-9 | 1.03e-11 | **derived, not guessed.** Domain decomposition changes the order of the force summation, so MPI results differ from serial at round-off and that difference grows over 400 steps. 1.03e-11 absolute is 5.8e-10 relative to the largest component. The failure this bound exists to catch — a tensor multiplied by the rank count — would shift `zz` by 9.6e-3, nine orders larger. np=2 and np=4 give the *identical* difference, which is itself evidence that it is reordering and not rank scaling |
| `mpi_position_abs_a` | 1e-10 | < 1e-10 | as above |

## What is NOT claimed

This is a **400-step deterministic run**. It is dynamics and regression
coverage: that the thermostat holds the temperature near its target with bounded
excursions, that the constraints stay satisfied, that the virial and pressure
stay mutually consistent, and that a restart reproduces the trajectory.

It is **not** a statistical validation of the canonical ensemble. 400 steps of
one deterministic trajectory cannot establish an ensemble distribution, and
nothing here claims it does. The temperature bounds are dynamics bounds, not
fluctuation-theorem predictions.

Every expected value is a **GMD regression result**, with two exceptions that
rest on identities rather than recorded numbers: the kinetic-plus-virial
pressure identity, and total virial = provider + constraint. No external engine
was used; see `reference.json` for the full provenance statement.

## Reproducing

```sh
python3 validation/constrained_nvt_cluster/generate_fixture.py
python3 validation/run_validation.py --case constrained_nvt_cluster \
    --gmd <build>/gmd --gmd-validate <build>/gmd_validate --work-root <work> \
    [--mpiexec mpirun --mpi-ranks 1,2,4]
```
