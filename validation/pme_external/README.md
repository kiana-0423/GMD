# PME External Validation Plan

Status: **TODO / design only**.

This directory records the planned independent validation for replicated PME.
It must not be used to claim that PME is independently validated until LAMMPS
PPPM, OpenMM PME, or another external reference result is committed with exact
version, units, settings, outputs, and tolerances.

## Purpose

Upgrade the current `validation/static_coulomb/reference_pme.json` from
`provisional_gmd_baseline_2026-05-23` regression evidence to an independent
scientific reference.

## Candidate Reference Engines

- LAMMPS PPPM with an explicitly recorded version, units, boundary conditions,
  cutoff, mesh accuracy, PPPM order, and force/energy output.
- OpenMM PME with recorded version, platform, precision mode, cutoff, PME
  tolerance, box vectors, and force/energy output.

Use one reference engine first; adding a second cross-check is useful but not
required for `v2.4`.

## Suggested Cases

1. `small_charged_particles`

   Fixed-position charged particles in an orthorhombic periodic box. No bonded
   terms. Compare Coulomb total energy and per-atom forces.

2. `neutral_ion_pair`

   A neutral periodic ion-pair system with nonzero LJ parameters. Compare total
   energy, Coulomb component where the reference engine exposes it, and forces.

3. `molecular_14_scaling`

   Four-atom molecular chain with 1-2/1-3 exclusions and 1-4 Coulomb scaling.
   This case should reuse the existing special-pair topology inputs where
   possible and document how the reference engine expresses special bonds.

## Required Comparisons

- Total Coulomb energy.
- Total potential energy when LJ is present.
- Per-atom force RMS error.
- Maximum force component error.
- `np=1`, `np=2`, and `np=4` GMD consistency for the same PME settings.
- Grid/order/cutoff convergence, at minimum:
  - fixed cutoff with increasing PME grid;
  - fixed grid with PME order 4 vs 6 where both are supported;
  - one tighter reference setting to show convergence direction.

## Proposed Directory Layout

```text
validation/pme_external/
├── README.md
├── small_charged_particles/
│   ├── gmd/
│   ├── lammps/          # or openmm/
│   ├── reference.json
│   ├── tolerance.json
│   └── analyze.py
├── neutral_ion_pair/
└── molecular_14_scaling/
```

## Acceptance Target

Set tolerances only after the external reference is generated. Expected first
targets for a small fixed-coordinate PME case are:

- energy absolute error: documented per case, likely `1e-5` to `1e-4` eV until
  mesh convergence is characterized;
- force RMS error: documented per case, likely `1e-5` to `1e-4` eV/Angstrom;
- max force error: documented per case and justified by grid/order/cutoff.

These numbers are placeholders for planning, not current validation claims.

## Reuse Existing Framework

Prefer extending `validation/run_validation.py` and `validation/common.py` so
the output remains JSON/CSV-compatible with the existing static validation
cases. Avoid adding heavy Python dependencies beyond the current numpy-style
analysis stack unless the selected external engine requires its own driver.
