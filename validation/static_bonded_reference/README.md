# static_bonded_reference

Independent external validation of GMD's bonded terms — bonds, angles, proper
dihedrals and improper dihedrals — against LAMMPS.

This case exists because the internal finite-difference tests can only show that
GMD's forces are consistent with GMD's *own* energy. They cannot show that the
energy function itself is the intended one. A proper dihedral force bug survived
in this codebase precisely because it was self-consistent under the checks that
existed at the time.

## Reference engine

| | |
|---|---|
| Engine | LAMMPS |
| Version | `22 Jul 2025 - Update 5` (Homebrew bottle, `lmp_serial`) |
| Packages | `MOLECULE` (bond/angle/dihedral/improper styles) |
| Units | `real` — energy kcal/mol, distance Å, force kcal/mol/Å |
| Command | `lmp_serial -in lammps.in -log lammps.log` |

Regenerate with:

```bash
python3 validation/static_bonded_reference/generate_reference.py --lammps $(which lmp_serial)
```

`generate_reference.py` is provenance tooling and is **not** invoked by the
validation suite — the suite reads the checked-in `reference.json`, so running
validation never requires LAMMPS to be installed.

## Units and conversion

A GMD molecular `.ff` file takes energies in **kcal/mol**, distances in
**Å** and angles in **degrees** (see the format documentation in
`src/io/config_loader.cpp`), converting to eV internally with

```
kcal_to_eV = 4.336410e-2
```

LAMMPS runs in `units real`, which takes exactly the same units, so **every
coefficient in `lammps.in` is numerically identical to the one in `bonded.ff`**.
`generate_reference.py` multiplies LAMMPS energies by `kcal_to_eV` and forces by
the same factor (Å is common to both).

Because the identical constant is applied to both sides, its seven-digit
precision does not limit the comparison — the agreement achieved is ~1e-15 eV.

**This is also a limit on what the comparison proves.** The reference is
converted into eV using GMD's *own* declared constant, so that constant appears
on both sides of every comparison and cancels. If `kcal_to_eV` were itself
wrong, LAMMPS's converted values would be wrong by exactly the same factor and
this case would still pass. **The numerical value of the unit-conversion
constant is therefore not independently verified here.** What is verified is
everything downstream of it: the functional forms, the angle and sign
conventions, periodic wrapping behaviour, and the energies and forces they
produce. Verifying the constant itself would require an independent physical
determination of kcal/mol in eV, which is outside the scope of this case.

## Functional forms and convention mapping

| Term | GMD | LAMMPS style | Equivalence |
|---|---|---|---|
| Bond | `V = k (r - r0)²` | `bond_style harmonic`, `E = K (r - r0)²` | identical |
| Angle | `V = k (θ - θ0)²` | `angle_style harmonic`, `E = K (θ - θ0)²` | identical |
| Proper dihedral | `V = k [1 + cos(n φ - δ)]` | `dihedral_style charmm`, `E = K [1 + cos(n φ - d)]` | identical, **including the sign of φ** (proven below) |
| Improper | `V = k (φ - φ0)²`, φ = torsion(i,j,k,l) | `improper_style harmonic`, `E = K (χ - χ0)²` | same form, **opposite angle sign**: `φ = -χ`, so `φ0 = -χ0` |

`dihedral_style charmm` takes a fourth coefficient, the 1-4 pair weighting
factor. It is set to `0.0` so the style contributes no non-bonded term, matching
GMD, where 1-4 scaling is handled separately by the special-pair map.

### Proper dihedral: sign convention *is* equivalent, and this was proven

`E = K[1 + cos(nφ - δ)]` is an even function of φ when `δ = 0`, so a δ=0
comparison agrees whichever sign convention each engine uses and proves nothing
about the sign. The fixture therefore uses **δ = 60°**, which makes the energy
sign-sensitive. Both engines agree to 0.00e+00 relative error at δ = 60°, which
establishes that GMD's φ and the CHARMM φ have the same sign.

Verification performed while deriving this mapping:

| Experiment | GMD [eV] | LAMMPS·C [eV] | Relative error |
|---|---|---|---|
| proper dihedral, δ = 0 (sign-insensitive) | 0.0324842894903 | 0.0324842894903 | 6.4e-16 |
| proper dihedral, δ = 60° (**sign-sensitive**) | 0.0975304240662 | 0.0975304240662 | 0.0 |

### Improper sign convention: derived, not assumed

For this fixture GMD reports the improper torsion angle of atoms (1,2,3,5) as
**φ = −54.7915°**, while LAMMPS's `χ` for the same quadruple is **+54.79°**.
The two definitions differ by a sign.

That was established by measurement, not by flipping a sign until the numbers
agreed:

| Experiment | GMD [eV] | LAMMPS·C [eV] | Relative error |
|---|---|---|---|
| φ0 = 0 (even function, sign-insensitive) | 1.58625051305 | 1.58625051305 | 9.8e-16 |
| φ0 = +10° vs χ0 = +10° (naive, same sign) | 2.21810148833 | 1.06007521978 | **5.2e-01** |
| φ0 = −10° vs χ0 = +10° (**mapped**) | 1.06007521978 | 1.06007521978 | 6.3e-16 |

Since `E = k(φ - φ0)²` and `φ = -χ`, we have
`(φ - φ0) = (-χ + φ0) = -(χ - (-φ0))`, so squaring gives an identical energy
provided `φ0 = -χ0`. The fixture uses GMD `phi0 = -10.0` against LAMMPS
`chi0 = 10.0`, and agrees to 6.3e-16 in energy and 4.4e-15 in force.

The φ0 = 0 row is the control: it agrees under *either* convention, so it alone
would not have exposed the difference.

## Fixture

Five atoms in a **non-cubic orthorhombic** box, 18 × 21 × 24 Å, in a
deliberately non-symmetric, non-degenerate geometry (no collinear triples, no
planar torsions, all three box lengths distinct so an axis mix-up cannot hide):

```
bonds      1-2, 2-3, 3-4, 3-5
angles     1-2-3, 2-3-4, 2-3-5
dihedral   1-2-3-4   (proper,  δ = 60°)
improper   1-2-3-5   (harmonic, φ0 = -10° <-> χ0 = +10°)
```

`bonded_wrapped.xyz` is the same molecule translated by −7 Å in x and wrapped,
so bond 1-2 crosses the periodic boundary. Both engines must return the same
energy and forces for it as for the intact molecule; LAMMPS agrees with itself
across the two variants to 2.0e-15 eV.

**Atom ordering**: line order in the GMD `.xyz` is the same as LAMMPS atom-ID
order, 1 through 5. `write_dump ... modify sort id` guarantees the dump is in
atom-ID order.

Non-bonded terms are eliminated rather than matched: all charges are zero, and
`lj_cutoff` (1.2 Å) is shorter than the closest interatomic distance, so GMD
reports exactly `lj = 0, coulomb = 0`. LAMMPS uses `pair_style zero`. `analyze.py`
asserts the GMD non-bonded components really are zero, so the comparison against
a bonded-only reference cannot silently become meaningless.

## What is compared

- **Per-term energies.** GMD aggregates all bonded terms into one component, so
  each term is isolated with its own topology file (`bonded_*_only.top`) and
  compared against the matching per-term energy from a single LAMMPS run.
- **Total bonded energy**, intact and wrapped.
- **Every force component on every atom**, intact and wrapped.
- **Wrapping invariance**: intact and wrapped totals must agree.

## Tolerances

| Quantity | Tolerance | Why |
|---|---|---|
| Energy | 1e-12 eV | Observed agreement is ~2e-15 eV, three orders tighter. The limit is double-precision accumulation over a handful of terms, not any physical approximation — both engines evaluate the same closed-form expressions. |
| Force (max component, RMS) | 1e-11 eV/Å | Observed agreement is ~6e-15 eV/Å. Forces involve more arithmetic than energies (cross products, normalisations), so the bound is set an order looser than the energy one while still being ~4 orders tighter than observed. |

These are absolute tolerances on quantities of order 1 eV and 1–10 eV/Å, so they
are effectively relative tolerances of 1e-12.

## Sensitivity

This case fails with the pre-fix dihedral force formula:

```
variant intact   energy_ok=True force_ok=False max_force_err=5.385e+00
variant wrapped  energy_ok=True force_ok=False max_force_err=5.385e+00
```

Note that the **energies still pass** — the defect was in the forces only. An
energy-only reference would not have caught it, which is why per-atom forces are
compared here.

## Scope of the claim

External validation is claimed, for **bond, angle, proper dihedral and
improper**, of:

- the **functional forms** — that GMD evaluates the intended expressions;
- the **angle and sign conventions**, including the proper-dihedral phase sign
  (proven with a sign-sensitive δ = 60° case) and the improper sign mapping
  (derived, with a disagreeing control);
- **periodic wrapping behaviour** — the wrapped and intact variants agree, in
  both engines;
- the resulting **energies** and **per-atom forces**.

Explicitly **not** claimed:

- The **numerical value of the kcal/mol → eV conversion constant.** The
  reference is converted with GMD's own declared constant, so it cancels out of
  the comparison and an error in it would not be detected. See *Units and
  conversion* above.

Also not covered here:

- Non-bonded (LJ / Coulomb) terms — covered by the other static cases.
- The virial. LAMMPS's per-term virial decomposition and GMD's do not have
  proven-equivalent semantics for this fixture, so no virial comparison is made.
  GMD's virial is validated internally by finite difference instead
  (`tests/virial_finite_difference_tests.cpp`), and only its diagonal.
- Any dihedral or improper style other than the two mapped above.

The internal finite-difference gradient test
(`tests/bonded_force_gradient_tests.cpp`) is kept as a separate, independent
check: it verifies self-consistency of forces with energy, while this case
verifies that the energy function is the intended one.
