# pme_external

An independent external reference for replicated PME: GMD's PME energy, forces
and electrostatic virial compared against **OpenMM PME** and **LAMMPS PPPM**.

This directory previously held a design plan and no result. It now holds a
generated reference, the inputs that produced it, and a comparison that runs in
normal CI.

## Which reference is which

The three kinds of reference in this repository are not interchangeable, and
this case uses all three at once:

| Kind | File | What it proves |
|---|---|---|
| **External** | `reference_external.json` → `reference` | OpenMM 8.6.0 PME, an unrelated implementation by unrelated authors |
| **External, second engine** | → `cross_reference` | LAMMPS PPPM and LAMMPS exact Ewald; the only source of a virial tensor |
| **Analytic** | → `analytic_ewald` | `validation/analytic_references.py`, an explicit k-vector Ewald sum in Python |
| **GMD's own output** | → `gmd` | recorded so the comparison is auditable; **never** the reference |

`analytic_ewald` uses the same formulas and the same Coulomb constant as GMD, so
it is an independent *implementation* but not an independent *derivation*. It is
cross-checked here against LAMMPS' exact Ewald, which agrees with it to 8.6e-08
eV in energy and 4.0e-08 eV/Å in force.

## The fixture

`charges.xyz` — 12 atoms, `Na`/`Cl` labels used only to make atom identity
legible in each engine's output.

- **Exactly charge neutral.** Charges are multiples of 1/16, so the sum is
  exactly `0.0` in binary floating point and no neutralising-background term is
  triggered in any engine. The background correction is therefore *not* under
  test here.
- **Non-cubic**: 18 × 22 × 26 Å. A cubic box hides axis-permutation mistakes.
- **No symmetry.** No mirror plane, no inversion centre, no two atoms sharing a
  coordinate. Every one of the 36 force components is substantially nonzero, so
  nothing passes by being zero on both sides. The off-diagonal virial components
  are likewise all nonzero.
- Minimum image separation 4.68 Å; cutoff 8.0 Å is strictly below `L_min/2 = 9`.
- Small enough that an exact Ewald sum at `kmax = 24` is instant.

`charges_wrapped.xyz` is the same physical configuration with each atom written
in a different periodic image — mixed signs on every axis, several atoms moved on
more than one axis at once. Energy and forces must be identical.

`charges_translated.xyz` is a rigid translation by `(3.7, −5.3, 8.9)` Å, which is
not a lattice vector. This is *not* an exact symmetry of a mesh method and is
checked only to the mesh error.

## Settings, and how aligned each one is

α = 0.35 Å⁻¹, real-space cutoff 8.0 Å, grids 16³/32³/64³/128³.

**Exactly aligned:** splitting α, reciprocal grid, real-space cutoff, box
vectors, periodicity on all three axes, absence of exclusions, charge neutrality,
double precision.

**Not alignable:**

- **B-spline order.** OpenMM's PME interpolation order is fixed at 5 in its
  sources and exposed by no API. GMD supports 4 and 6. LAMMPS' `kspace_modify
  order` does match, and both 4 and 6 are run there.
- **Reciprocal influence function.** LAMMPS PPPM uses an optimised Green's
  function chosen to minimise RMS force error; GMD and OpenMM both use the
  Essmann smooth-PME kernel. LAMMPS is a different mesh approximation of the same
  exact sum, not a different implementation of the same approximation.
- **Coulomb constant.** GMD uses `14.3996` eV·Å/e²; both engines use the CODATA
  value (LAMMPS `14.399645`, OpenMM `138.935457` kJ/mol·nm/e², measured). GMD is
  low by **3.16e-06 relative**. Because every Ewald term carries exactly one
  factor of `k_e`, each engine's result is rescaled by `k_e^GMD / k_e^engine`,
  which is exact. Both constants are measured from the engines at generation
  time rather than quoted. Uncorrected, this offset would exceed the grid-128
  mesh error and look like a convergence floor.

**Cannot be compared at all:** the decomposition into real / reciprocal / self /
background. Neither engine exposes a split with a proven-equivalent definition,
so only the totals are compared.

## Reproducing the reference

Requires OpenMM and LAMMPS. Neither is needed to *run* the validation.

```bash
python3 validation/pme_external/generate_external_reference.py \
    --work-dir <work> \
    --gmd-validate <build>/gmd_validate \
    --openmm-python <python-with-openmm> \
    --lammps <lammps-binary> \
    --output <work>/reference_external.json
```

Run from `<repo>`. Writes only into `<work>`; nothing is written into the source
tree. To adopt the result, copy `<work>/reference_external.json` over
`validation/pme_external/reference_external.json`.

The script **fails** if either engine is missing or if `--work-dir` is inside the
repository. It never falls back to GMD output to stand in for an external
result: that failure mode is exactly what this directory exists to rule out.

OpenMM was installed with `pip install openmm` into a virtual environment outside
the repository; LAMMPS from the Homebrew `lammps` formula. Exact versions, git
revision, platform, measured constants and input hashes are recorded in the
`provenance` block of the reference file.

## Running the validation

```bash
python3 validation/run_validation.py --case pme_external \
    --gmd <build>/gmd --gmd-validate <build>/gmd_validate --work-root <work>
```

or `ctest -R gmd_validation_pme_external`. 44 checks, roughly a tenth of a
second, no external engine required.

## What is checked

- The reference is still an external engine's output: `source` prefix, engine
  name and version present, fixture and run-file SHA-256 hashes match the ones
  recorded at generation, and the run file's α, cutoff and grid still match the
  reference block.
- **The reference is not GMD's own output.** Two PME implementations at
  different interpolation orders cannot agree below 1e-08; a smaller gap means
  the file was replaced with GMD-generated numbers. Observed gap: 1.13e-06 eV.
- Total energy vs OpenMM, absolute and relative.
- Every force component of every atom vs OpenMM; failures name the atom and axis.
- Force RMS, and net force for both codes.
- The same against LAMMPS PPPM.
- GMD PME vs the exact Ewald sum (mesh error), and GMD's own Ewald against the
  same sum in Python (6.2e-15 eV).
- All nine virial components against LAMMPS, for both GMD Ewald and GMD PME.
  LAMMPS reports six because the tensor is symmetric; GMD's symmetry is measured
  rather than assumed, which is what makes six cover nine.
- Wrapped-image equivalence, exact to 1.7e-13 eV.
- Translation invariance, to the mesh error.

## Tolerances

Derived in `tolerance.json`, which carries the measurement each one comes from.
None is an observed error rounded up. The reference-settings tolerance is the
triangle-inequality bound `|GMD − exact| + |OpenMM − exact|`, doubled; that bound
is tight here rather than loose, because the two codes deviate from exact Ewald
in opposite directions.

## Limits

- **No virial from OpenMM.** OpenMM exposes no virial tensor in any form, so the
  tensor rests on LAMMPS alone.
- **No decomposition.** Only total electrostatic energy, forces and virial.
- **LAMMPS PPPM stops at grid 64.** `kspace_modify mesh N N N` aborts that build
  with SIGSEGV for every `N > 64` on this fixture, with pinned or auto-selected
  gewald and order alike. GMD and OpenMM are swept to 128.
- **Coulomb-only, static, serial.** No LJ, no bonded terms, no exclusions or 1-4
  scaling, no dynamics, no MPI.
- **Orthorhombic only**, because `gmd::Box` cannot represent anything else.
- **GMD's `k_e` is 3.16e-06 low.** Corrected exactly in this comparison, but it
  is a real systematic offset in every Coulomb energy and force GMD reports, and
  it is *not* fixed here — changing it would move every checked-in baseline in
  the repository.
