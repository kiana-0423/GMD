# GMD — Good Molecular Dynamics  `v2.4`

Release history is tracked in [CHANGELOG.md](CHANGELOG.md). The current documented release is `v2.4` from `2026-05-24`.

GMD is a C++20 molecular dynamics engine built around a small set of composable runtime abstractions:

- `ForceProvider` — force and energy evaluation (LJ, bonded, Ewald/PME, ML)
- `CompositeForceProvider` — sums multiple providers; total energy is automatically accumulated
- `NeighborBuilder` — Verlet-list construction with skin-distance rebuild check
- `Integrator` — velocity-Verlet time stepping with optional thermostat / barostat hooks
- `Simulation` — orchestrator that wires everything together

---

## What's New in v2.4

| Feature | Files |
|---|---|
| **Validation-progress release status** — README and validation docs now distinguish analytic validation, regression baselines, tested workflows, and prototype interfaces | `README.md` `validation/README.md` |
| **Special-pair / 1-4 scaling validation** — topology-derived 1-2/1-3 exclusions and 1-4 LJ/Coulomb scaling are documented as analytic-validation covered | `include/gmd/system/special_pair_map.hpp` `validation/static_special_pairs/` |
| **SHAKE/RATTLE and checkpoint/restart status corrected** — both are documented as tested features, with MPI constraint and checkpoint scalability limits called out | `README.md` `tests/restart_continuity.py` `tests/constraint_solver_tests.cpp` |
| **Independent external PME validation** — replicated PME is now compared against OpenMM PME and LAMMPS PPPM on a Coulomb-only fixture, with a convergence study showing all three codes approaching the same exact Ewald limit | `validation/pme_external/` |
| **Release evidence tooling** — collect serial/MPI configure, build, CTest, validation summaries, and environment info without touching existing build dirs | `scripts/collect_release_logs.sh` `docs/release_logs/README.md` |
| **Manual/nightly full MPI CI** — full MPI CTest workflow separate from ordinary PR CI | `.github/workflows/full-mpi.yml` |

### Correctness fixes (post-v2.4)

These change simulation results for the configurations they affect.

| Fix | Files |
|---|---|
| **Degrees of freedom now account for constraints and the COM setting** — one shared `compute_degrees_of_freedom()` replaces the hardcoded `3N - 3` in both thermostats and in trajectory temperature output. A constrained run previously reported a temperature that was too low, so the thermostat drove the system hotter than its target | `include/gmd/integrator/thermostat.hpp` `src/integrator/thermostat.cpp` `src/integrator/{nose_hoover,velocity_rescaling}_thermostat.cpp` `src/integrator/velocity_verlet_integrator.cpp` `app/gmd_main.cpp` |
| **Neighbor-list rebuild decided collectively under MPI** *(correctness safeguard, not a result-changing fix)* — local rebuild flags are combined with a logical OR so every rank rebuilds during the same force evaluation, or none does. Under the current domain-decomposition path ghost cleanup already clears the list every step, so this is expected to be a no-op there and may not alter existing decomposed trajectories; it prevents divergent control flow wherever a valid list survives between steps, and protects future implementations | `src/core/simulation.cpp` `src/parallel/mpi_communicator.cpp` |
| **Forces refreshed after a barostat volume change** — `VelocityVerletIntegrator::step()` no longer returns holding forces computed for the pre-scaling geometry; `Simulation::step()` skips the refresh when the cell did not change | `src/integrator/velocity_verlet_integrator.cpp` `src/core/simulation.cpp` |
| **Physically derived virial for Ewald, PME and bonded terms** — see *Virial and pressure* below | `src/force/{ewald,pme,bonded}_force_provider.cpp` `include/gmd/force/special_pair_coulomb.hpp` |
| **Box dimensions validated** — zero, negative and non-finite edge lengths are rejected where a box is installed, instead of reaching periodic wrapping or the neighbor builder | `include/gmd/system/box.hpp` `include/gmd/system/system.hpp` `src/system/verlet_neighbor_builder.cpp` |
| **Proper dihedral and improper forces corrected** — the terminal-atom force had a flipped sign and the middle-atom projection used the coefficients for the opposite `b1` convention. Any run with proper dihedrals or impropers was integrating incorrect torsional forces. Found by component-wise virial validation: a torsion angle is invariant under isotropic scaling, so the error is invisible in `tr(W)` | `src/force/bonded_force_provider.cpp` |
| **Nose-Hoover restart validates, never overwrites, the DOF** — a checkpoint's `dof` is checked against the run's authoritative count and an incompatible restart is rejected rather than silently continued with a mismatched thermostat mass | `src/integrator/nose_hoover_thermostat.cpp` `include/gmd/integrator/thermostat.hpp` |
| **Constraint lists normalised** — `(i,j)` and `(j,i)` collapse to one constraint, exact duplicates are dropped, conflicting target distances and self-constraints are rejected, so the DOF subtraction counts distinct constraints | `src/integrator/constraint_solver.cpp` `include/gmd/integrator/constraint_solver.hpp` |
| **PME reciprocal forces were identically zero** — `bspline_deriv(u, p)` evaluates `M_(p-1)`, but `bspline()` implemented only orders 4 and 6 and returned `0.0` for anything else. Every derivative weight was therefore zero and **`coulomb pme` applied no reciprocal electrostatic force at all**, at any order. Orders 2, 3 and 5 are now implemented — every supported order needs its predecessor, all the way down. Energies were unaffected, which is why the existing PME regression baseline never noticed | `src/force/pme_force_provider.cpp` |
| **PME force interpolation was missing the mesh-point-count factor** — the energy uses an unnormalised forward transform while `fft3d(..., true)` divides by `K1·K2·K3`, so `dE/dQ(j) = N·IFFT[G·Q̂](j)`. The factor `N` was absent, leaving every PME reciprocal force `N` times too small. Masked by the defect above, which zeroed the term outright | `src/force/pme_force_provider.cpp` |
| **PME B-spline order 6 was wrong on three of its six intervals** — the hard-coded quintic polynomials on `[2,3)`, `[3,4)` and `[4,5)` did not match `M_6`; the spline went negative, summed to 0.9 instead of 1 under partition of unity, and broke the symmetry `M(u) = M(6-u)`. `pme_order 6` produced nonsense energies (13401 eV where the correct value is -3.05 eV) | `src/force/pme_force_provider.cpp` |
| **`velocity_init random` was rank-dependent** — one `std::mt19937` advanced in local storage order, so an atom's draw was decided by its array position and an MPI run never reproduced the serial trajectory for the same seed. The global temperature rescale hid it. Draws are now keyed on `(seed, stream, global atom tag, component)` | `include/gmd/core/keyed_random.hpp` `src/system/initializer.cpp` `include/gmd/system/initializer.hpp` |
| **Nosé–Hoover and Berendsen relaxation times were consumed in internal time units while being written, printed and defaulted in femtoseconds** — `tau = 100 fs` relaxed on 1018.05 fs and `tau_P = 2000 fs` coupled on 20361 fs, both exactly `T` too long. Neither was visible from inside its own file: every quantity was self-consistent | `src/io/config_loader.cpp` `include/gmd/io/config_loader.hpp` `app/gmd_main.cpp` |
| **The internal time constant was rounded to seven figures and inversely named** — `kInternalTimeUnitsPerFs = 1.018051e+1` was 4.206204e-07 above `Å·√(amu/eV)`, ~2700× the CODATA uncertainty, and was *divided into* a femtosecond timestep despite its name | `include/gmd/core/physical_constants.hpp` `src/io/config_loader.cpp` |
| **The Berendsen barostat compared a target in bar against a pressure in eV/Å³** — it subtracted the two directly with no conversion on either side, so a run asking for 1 bar was asking for 1 eV/Å³, or 1602176.634 bar. The comparison sets the *sign* of the coupling, so for any ordinary target the box was pushed the same direction regardless of the true pressure. The instantaneous pressure is now converted to bar, where `beta` already lived | `src/integrator/berendsen_barostat.cpp` `include/gmd/integrator/berendsen_barostat.hpp` |
| **One authoritative bar ⇄ eV/Å³ conversion** — the factor existed as two independent literals, both `6.2415091e-7` and both 4.091837e-09 relative high. It is a pure unit identity with four exact ingredients, so there was no precision to round to; the terminating reverse direction `1 eV/Å³ = 1602176.634 bar` is now the single definition | `include/gmd/core/physical_constants.hpp` `src/io/trajectory_writer.cpp` `include/gmd/integrator/mc_barostat.hpp` |
| **`MLForceProvider` states its virial contract** — it left `virial`/`virial_valid` untouched, so a reused `ForceResult` would carry a previous provider's tensor and have it attributed to the model. It now clears both explicitly. `Σ r ⊗ F` is deliberately *not* synthesised: for a periodic, cell-dependent model that expression is not the virial | `src/force/ml_force_provider.cpp` |

---

## Virial and pressure

The virial tensor `W` feeds the barostats through `P = (2*KE + tr(W)) / 3V`, and
each force term now contributes it in the form that term actually requires:

| Term | Virial |
|---|---|
| LJ, Ewald/PME real space, special-pair Coulomb correction | pair virial `r_ij ⊗ F_ij` from the minimum-image separation |
| Bonds, angles, dihedrals, impropers | per-interaction `Σ_a (r_a - r_ref) ⊗ F_a`, with positions taken relative to one atom of the interaction |
| Ewald/PME reciprocal space | `Σ_k E_k [δ_ab - 2(1/4α² + 1/k²) k_a k_b]` |
| Ewald/PME net-charge correction | isotropic `U_net · δ_ab` (the term scales as `1/V`) |
| Ewald/PME self-energy | none — it does not depend on the cell |
| SHAKE/RATTLE constraints | `Σ_c (2Λ_c/dt)·(r_c ⊗ r_c)`, the endpoint force recovered from the converged RATTLE multipliers |
| MLForceProvider | none — the provider reports `virial_valid == false` |

### Virial validation coverage

Every source above is validated component by component against a reference that
does not share the engine's arithmetic. The evidence is not uniform across
sources, so it is spelled out rather than summarised as "validated":

| Source | Analytical reference | Independent numerical reference | Nine components | PBC | MPI np=1/2/4 | External engine |
|---|---|---|---|---|---|---|
| LJ pair | `W_ab = r_a F_b`, force cross-checked against `-dV/dr` | rotation covariance `R W Rᵀ` | yes | x/y/z faces, corner, translate-and-wrap | yes | **no** |
| LJ mixed types, explicit pair override | yes | — | yes | — | — | **no** |
| 1-4 LJ scaling (0, partial, 1) | yes, plus exact linearity in the scale | — | yes | — | — | **no** |
| Bond, angle, proper dihedral, improper | force moment `Σ_a (r_a - r_ref) ⊗ F_a`, built without `accum_interaction_virial()` | forces first verified against a central difference of an independently written energy; rotation covariance | yes | intact vs wrapped molecule | yes, chain spanning rank boundaries | **no** (LAMMPS covers the *forces*, not the tensor) |
| Ewald real space | pair identity, independent regrouping of `-dV/dr` | — | yes | minimum image | yes | **yes** — LAMMPS `kspace_style ewald`, all six independent components to 1.6e-7 eV |
| Ewald reciprocal | independent re-derivation of the cell derivative | **shear strain derivative** of a test-only triclinic reciprocal energy | yes | k-space, not applicable | yes | **yes** — same LAMMPS comparison; the reported tensor is the sum of all Ewald terms |
| Ewald net-charge correction | isotropic `U_net·δ_ab`, incl. off-diagonals required to vanish | — | yes | not applicable | yes | **no** |
| Ewald self term | required to be identically zero | — | yes | not applicable | yes | **no** |
| Special-pair Coulomb correction | `(scale-1)` times the bare pair tensor | — | yes | minimum image | yes, corrections spanning ranks counted once | **no** |
| PME real space | compared directly against the Ewald pair sum | — | yes | minimum image | yes | **yes** — included in the LAMMPS PPPM tensor comparison |
| PME reciprocal mesh | — | PME → Ewald convergence under mesh refinement, **and** the reciprocal metric derivative of an independently written PME energy | yes | not applicable | yes | **yes** — LAMMPS PPPM, all six independent components to 8.6e-7 eV; energy and forces additionally against OpenMM PME |
| SHAKE/RATTLE constraints | endpoint force from the RATTLE multipliers; rigid-rotor closed form | rigid-body second moment `-ω² Σ m a ⊗ a` | yes | yes | yes | **no** |
| CompositeForceProvider | exact component-wise addition | — | yes | not applicable | yes | not applicable |
| MLForceProvider | reports no virial; `virial_valid == false` asserted | — | not applicable | not applicable | rejects MPI outright | not applicable |

**External-engine coverage of the virial is real but partial.** The
electrostatic virial — every Ewald and PME term, summed — is compared against
LAMMPS in `validation/pme_external`: all six independent components (and so, via
an asserted symmetry, all nine) agree to 1.6e-7 eV against `kspace_style ewald`
and 8.6e-7 eV against PPPM, on a fixture whose off-diagonals are all
substantially nonzero. That is a genuine cross-code check of a tensor.

Every other row is still marked **no**, and means it. LJ, bonded, special-pair
and constraint virials have no external reference: LAMMPS reports one combined
pressure tensor, so isolating a single contribution would require reconstructing
it by subtraction under an assumption of equivalence that has not been proven
here. OpenMM exposes no virial at all, in any decomposition. The constraint
virial is additionally a dynamical quantity that cannot be compared from a
static configuration.

**Why the off-diagonal claim is not a finite-difference claim.** `Box` stores
three edge lengths, so the engine represents orthorhombic cells only and no
shear strain is expressible. `tests/virial_finite_difference_tests.cpp`
therefore validates the trace and the three diagonal components and explicitly
does not validate the rest. The off-diagonal components are covered here
instead by references that do not need the engine to shear: exact pair
identities, bonded force moments, and — for the reciprocal-space terms, where
no pair identity exists — a test-only reciprocal energy written against a
general 3×3 cell matrix, which *can* be sheared. Differentiating that reference
under a general strain yields all nine components, and the engine's analytic
tensor is compared against it at the orthorhombic configuration both agree on.

**Rotation covariance is not a substitute for a shear derivative.**
`W → R W Rᵀ` is a necessary condition on any Cartesian rank-2 tensor and it does
constrain the off-diagonal components, but it tests how the tensor *transforms*,
not that it is the derivative of anything. It also only applies where the
periodic lattice is irrelevant, since rotating a configuration inside a fixed
orthorhombic box does not rotate the lattice — which is why the reciprocal-space
terms use a 90° axis permutation (exact for an orthorhombic cell) plus the
strain-derivative reference, and not a general rotation.

Tests: `tests/virial_source_inventory_tests.cpp`,
`tests/pme_reciprocal_virial_tests.cpp`, `tests/mpi_virial_sources.cpp`, with
the shared references in `tests/virial_reference.hpp`.

### Reproducible random velocity initialization

`velocity_init random` used to be **rank-dependent**: the same seed gave a
different velocity field at every rank count, so an MPI run never reproduced
the serial trajectory.

**Root cause.** `VelocityInitializer` advanced one `std::mt19937` in a loop
over local storage, so the draw a physical atom received was decided by its
position in the array. Under decomposition every rank starts from the same
seeded state and hands it to its own first local atom — a different physical
atom on every rank and at every rank count. The global temperature rescale then
hid it: total kinetic energy is forced to the target either way, so step 0
reported an identical temperature *and* an identical potential energy while the
underlying field differed.

**The fix.** The draw is now a pure function of identity, with no state carried
between atoms and no dependence on visitation order:

```
value = f(seed, stream, global atom tag, component)
```

`include/gmd/core/keyed_random.hpp` implements it with SplitMix64 (Steele, Lea
and Flood, OOPSLA 2014 — the finalizer behind `java.util.SplittableRandom`):
two multiply-xorshift rounds over documented constants, exact integer
arithmetic. Each input is absorbed through a full mixing round so it avalanches
over the whole key.

| Detail | Choice, and why |
|---|---|
| Uniform | `((bits >> 12) + 0.5) / 2^52`, strictly inside `(0,1)`. **Twelve** bits, not eleven: with 53 the largest value is `(2^53 − 1) + 0.5`, which is not representable — doubles are spaced 1.0 there — so it rounds half-to-even up to `2^53` and the quotient is exactly `1.0`, putting `log(u)` at exactly `0`. A test asserts that the 53-bit variant *does* round to 1.0, so the reason is recorded rather than remembered. |
| Bias | a bit shift and a power-of-two division; there is no modulus anywhere |
| Components | x, y and z each draw their **own** pair of uniforms from their own key — never `cos` and `sin` of a shared pair, which makes two outputs exactly dependent. Half the entropy per pair is discarded, which costs nothing here. |
| Streams | `RandomStream` separates draws, so adding a random feature later cannot shift the numbers this one produces |

**Tag validation.** The tag is now load-bearing, so it is checked before use.
Negative tags are rejected, and duplicates are found by a collective
gather-sort whose verdict is broadcast. A duplicate can **span ranks** — each
rank's own tags unique, two ranks sharing one — and only a global check sees
that. Every rank throws the same diagnostic naming the same tag: a rank that
noticed nothing would wait alone in the next collective, and a bad input would
present as a hang rather than an error. There is no silent fallback to the
array index; that fallback was the defect.

**Ghost atoms are not initialized.** A ghost carries its owner's tag, so keyed
on the tag it would draw the owner's velocity — the right value. What would be
wrong is the reductions: that atom's momentum and energy would enter the
centre-of-mass sum and the rescale twice. Sampling and all three reductions
iterate owned atoms only, and the tag validation likewise, so a halo's repeated
tags are not mistaken for duplicates.

#### What is guaranteed

| Property | Guarantee |
|---|---|
| Same seed, same tags/masses/configuration | same field |
| np=1 vs np=2 vs np=4, and reversed local storage | identical **to 2.8e-17** |
| The draws themselves | **bitwise** identical by tag |
| Trajectory: serial vs np=1/2/4, 50 steps | **identical logs**; final state bitwise identical at np=1, within 8.9e-16 (~28 ulp) beyond |
| Serial run split across a checkpoint | **bitwise** identical to the continuous run |
| Changing the seed, or any atom's tag | changes that atom's draw |
| Cross-platform | a few ulp, not bitwise — see below |

The residual is **reduction round-off**, not the draws: the centre-of-mass and
kinetic-energy sums are accumulated in storage order within a rank and combined
by `MPI_Allreduce` across ranks, and neither order is the same at every rank
count. That reaches every atom through one shared shift and one shared factor.
Before the fix all 24 atoms of the equivalence fixture differed with a worst
component of **4.7e-01** — the improvement is sixteen orders of magnitude.

Key derivation, the mixer and the uniform mapping are exact integer arithmetic
and identical on every platform and in every optimisation mode. `sqrt`, `log`
and `cos` are not guaranteed bit-identical across libm implementations, so
cross-platform agreement is to a few ulp. The integer path is pinned by
`tests/keyed_random_tests.cpp` against vectors derived independently in
`tests/keyed_random_reference.py` — a third implementation, in another
language, that the vectors came *from* rather than were captured from.

#### Serial compatibility and baselines

**Serial results change.** This is a different random stream, and the old one
is deliberately not preserved: keeping it would mean keeping traversal order as
the random identity. The **distribution** does not change — target temperature,
the `sqrt(k_B T/m)` width, centre-of-mass removal and the 3N−3 rescale are
untouched, and `tests/velocity_init_tests.cpp` predicts the entire field from
an independent reimplementation of the specification and matches it to 2.8e-17.

| Case | Change |
|---|---|
| `nve_lj_fluid` | drift +4.687e-08 → −1.563e-08 eV/atom/ps |
| `nvt_lj_fluid` | T mean 120.41 → 119.74 K, stddev 22.00 → 18.72 K |
| `npt_lj_fluid` | T mean 119.92 → 120.14 K, pressure 35.56 → 47.32 bar |
| `berendsen_npt_lj` | expand ratio 1.1102 → 1.1046 |
| `diffusion_lj_fluid` | 1.796 → 3.273 Å²/ps — **the only metric put outside its tolerance** |

That tolerance, 0.05, is a *fixed-seed reproducibility* bound and not a
statement about how well D is determined: running the unchanged fixture under
five seeds gives 3.273, 2.664, 2.260, 2.224 and 2.205 Å²/ps, a spread of ~48%.
The new value is an ordinary draw from that spread; the old one sat at its low
end. **Static energy, force and virial references are untouched** — none of
them assigns velocities.

**Restart** is unaffected: a restart does not install the velocity initializer
at all, so checkpoint velocities are resumed rather than resampled, and no RNG
state is serialized. `velocity_init input` and a checkpoint restart are
distinct: the former reads velocities from the `.xyz`, the latter from the
checkpoint.

### The internal time unit

GMD integrates

```
v += (F/m)·dt        r += v·dt
```

with `F` in eV/Å and `m` in amu. Neither line mentions seconds, so the time
unit is **not a free choice**: requiring `F/m` to be an acceleration in this
unit system fixes it.

| | |
|---|---|
| **Definition** | `T = Å·√(amu/eV) = √(m_u/e)·1e5 fs` |
| **Value** | `10.1805057178711931...` fs |
| **Velocity unit** | `√(eV/amu) = 0.0982269474...` Å/fs |
| **Uncertainty** | 1.57e-10 relative |

| Ingredient | Value | Source |
|---|---|---|
| `e` | 1.602176634e-19 J, exact | [physics.nist.gov/cgi-bin/cuu/Value?evj](https://physics.nist.gov/cgi-bin/cuu/Value?evj) |
| `m_u` | 1.66053906892(52)e-27 kg | [physics.nist.gov/cgi-bin/cuu/Value?ukg](https://physics.nist.gov/cgi-bin/cuu/Value?ukg) |
| Å | 1e-10 m, by definition | — |
| fs | 1e-15 s, by definition | — |

Unlike `k_B` and the pressure conversion, this one **is** uncertain: `m_u`
carries 3.1e-10 relative, halved by the square root. Both directions live in
`include/gmd/core/physical_constants.hpp` as
`gmd::kFemtosecondsPerInternalTime` and `gmd::kInternalTimePerFemtosecond`,
held to being exact reciprocals by `static_assert`.

**Units, stated.** Coordinates Å · velocities Å/`T` · masses amu · forces eV/Å ·
accelerations Å/`T²` · internal time `T` · **displayed time fs** · diffusion
coefficients **Å²/ps**.

**Three defects, one audit.**

| # | Was | Effect |
|---|---|---|
| 1 | `kInternalTimeUnitsPerFs = 1.018051e+1` | **+4.206204e-07** on `dt`; ~2700× the CODATA uncertainty. Its *name* was also the reciprocal of its use — it was **divided into** a femtosecond timestep. |
| 2 | Nosé–Hoover `tau` never converted | `Q = dof·k_B·T·tau²` built from femtoseconds, `ξ` integrated against an internal `dt`. **`tau = 100 fs` relaxed on 1018.05 fs**; `Q` was `T² = 103.6427`× too large. |
| 3 | Berendsen `tau_P` never converted | Same in `mu³ = 1 − beta·(dt/tau)·ΔP`. **`tau = 2000 fs` coupled on 20361 fs.** |

Neither tau defect was visible from inside its own file — every quantity was
self-consistent, and both objects are constructed straight from `RunConfig`.
`RunConfig` now carries `thermostat_tau_fs`/`thermostat_tau` and
`barostat_tau_fs`/`barostat_tau`, exactly as it already carried
`time_step_fs`/`time_step`: the `_fs` field is what the user wrote and what the
CLI prints, the unsuffixed one is what the object consumes. **No thermostat or
barostat equation changed.**

**What was already correct.** The trajectory `time[fs]` column
(`step × time_step_fs`, computed by the caller), the checkpoint's `time_fs`,
and the femtosecond-to-picosecond conversion behind the reported diffusion
coefficient. `Simulation` passes `force_time` to force providers in *internal*
units — a different quantity from the reported column, and not what any log or
checkpoint carries.

**Baseline impact.** Results-changing for every dynamics run, and substantially
for thermostatted ones — the thermostat is now about ten times more strongly
coupled.

| Case | Metric | Change | Tolerance |
|---|---|---|---|
| `nvt_lj_fluid` | T mean / stddev | 123.81 → 120.41 K / 15.66 → **22.00** K | 5.0 K |
| `npt_lj_fluid` | T mean / stddev | 126.44 → 119.92 K / 17.04 → 21.17 K | 10.0 K |
| `npt_lj_fluid` | pressure mean | 45.99 → 35.56 bar | 2000 bar |
| `diffusion_lj_fluid` | D | 1.79598 → 1.79603 Å²/ps | 0.05 |
| `nve_lj_fluid` | drift | **0** — bit-identical | 1e-06 |

The nvt standard deviation is the only metric in the repository this puts
*outside* its tolerance. Tighter coupling raising the fluctuation while pulling
the mean toward the run's 120 K target is the expected signature, not drift.
`diffusion_lj_fluid` moves only 2.8e-05 because that run configures no
thermostat at all. `nve_lj_fluid` has neither thermostat nor barostat, and its
drift metric — a difference of two six-decimal energies — does not resolve
4.2e-07. **No static energy, force or virial baseline involves time; none is
regenerated and none changes.**

**What is asserted.** `tests/time_unit_tests.cpp` reads no production constant
for its reference — it derives `T` from the SI definitions and cross-checks the
grouping, which matters: applying the metre and second conversions separately
rounds twice and lands one ulp low. It then measures the unit back out of the
integrator three ways: free flight (the distance a particle carrying one
velocity unit covers per femtosecond), a harmonic oscillator's period against
the exact `2π√(m/k)`, and exact `t`/`t²` scaling under a constant force in all
three components and both signs. Second-order convergence is checked separately,
so "the unit is right" is distinguished from "the integrator is accurate at one
step size". The relaxation times are recovered from the thermostat mass and
from the observable Berendsen coupling factor and required to match what was
asked for, at three values each.

### Berendsen NPT trajectory validation

`validation/berendsen_npt_lj` is the only trajectory-level coverage of the
Berendsen barostat — `npt_lj_fluid` uses the Monte Carlo barostat, so the
Berendsen path previously had only unit tests. It is a **dynamics regression
case, not an external reference**: every number is GMD's own.

It is a sign-sensitive pair on a compressed 32-atom LJ fixture whose mean
pressure is ~2300 bar. `expand.run` targets 200 bar and the cell must grow;
`compress.run` targets 4000 bar and it must shrink. The direction and magnitude
assertions do not read `reference.json`.

Both previously-corrected defects are caught by construction, with thresholds
**measured from the defective code**:

| Restored defect | expand | compress | Caught by |
|---|---|---|---|
| *(correct)* | 1.1102 | 0.9510 | — |
| bar compared against eV/Å³ | 0.9911 | 0.8355 | both runs compress → `expand_direction` |
| `tau` unconverted | 1.0129 | 0.9962 | response 0.0129 < `min_volume_response` 0.03 |

It also checks that volume, energies, temperature and pressure stay finite,
that the cell never leaves `[0.25, 4.0]×` its initial volume, and that a run
split across a checkpoint ends at a **bit-identical** volume. Unlike the other
long-dynamics cases it is registered as a CTest test: four 500-step runs of 32
atoms take well under a second, so the path is genuinely gated.

**Serial/MPI is reported, not asserted**, and the reason is a pre-existing
defect this case surfaced: see *Known limitations*.

### The pressure unit conversion

GMD computes pressure internally as `P = (2K + tr W) / 3V`, which in eV and
Ångström is an energy density in **eV/Å³**. Users never see that: they set a
target in **bar** and read a `P[bar]` log column. Exactly one conversion stands
between the two, and two things were wrong with it.

**One authoritative conversion, and it is exact.** The factor existed as two
independent literals — `6.2415091e-7` in `src/io/trajectory_writer.cpp` and
`6.2415091e-7` in `include/gmd/integrator/mc_barostat.hpp` — with nothing
keeping them equal, both **4.091837e-09 relative high**. They are now

```
gmd::kEVPerAngstromCubedToBar = 1602176.634          // include/gmd/core/physical_constants.hpp
gmd::kBarToEVPerAngstromCubed = 1.0 / kEVPerAngstromCubedToBar
```

Unlike `k_e` and `k_B` this is not a measured quantity — it is a pure unit
identity whose four ingredients are all exact by definition:

| Ingredient | Value | Source |
|---|---|---|
| bar | 100000 Pa, by definition | — |
| pascal | 1 J/m³, by definition | — |
| ångström | 1e-10 m, by definition | — |
| eV | 1.602176634e-19 J, exact | [physics.nist.gov/cgi-bin/cuu/Value?evj](https://physics.nist.gov/cgi-bin/cuu/Value?evj) |

```
1 bar = 1e5 J/m³ = 1e-25 J/Å³ = 1e-25 / 1.602176634e-19 eV/Å³
      = 500/801088317 = 6.24150907446076260777624098...e-7 eV/Å³
```

There is no uncertainty to round to, so eight significant figures had nothing
behind it. The **reverse** direction is the one written as a literal, because
unlike the forward direction it **terminates**: `1 eV/Å³ = 1602176.634 bar`,
exactly. Taken from that decimal, each direction is the nearest `double` to its
exact rational *and* the two are exact reciprocals in `double` arithmetic, so a
pressure converted out and back returns the original bits. Writing the forward
direction as the literal instead, or routing either through the SI constants,
rounds more than once and lands one ulp off. Three `static_assert`s hold the
pair to all of it.

**The Berendsen barostat was comparing bar against eV/Å³.** It computed
`(2K + tr W) / 3V` and subtracted the target pressure from it with no
conversion on either side. The target arrives in bar — it is the run input's
`pressure` field, set through the same `set_target_pressure()` the Monte Carlo
barostat reads and converts. So **a run asking for 1 bar was asking for 1 eV/Å³,
which is 1602176.634 bar.** That is not a scale error: the comparison sets the
*sign* of the coupling, so for any ordinary target the barostat pushed the box
the same direction regardless of the true pressure. `beta` is a compressibility
in bar⁻¹ — its default 4.5e-5 is liquid water's value in those units — so the
comparison belongs in bar, and it is the instantaneous pressure that is now
converted. No algorithm, sign convention or virial handling changed.

**Affected production paths.** The `P[bar]` column of the log and the `.xyz`
comment line; the Monte Carlo barostat's `P_ext·ΔV` term; the Berendsen
barostat's entire coupling. **Results-changing for reported pressures and for
pressure-controlled runs; Berendsen runs change qualitatively.**

**What does not change.** Forces, energies and the virial tensor. Pressure here
is reported and controlled, never an input to a force, so every static energy,
force and virial baseline is untouched by construction.

**Baseline impact.** One reference moved, and the trajectory behind it did not.

| Case | Metric | Relative change | Tolerance |
|---|---|---|---|
| `npt_lj_fluid` | pressure mean / stddev | +4.327e-09 / +4.752e-09 | 2000 bar |
| `npt_lj_fluid` | temperature mean / stddev | **0** — bit-identical | 10.0 K |
| `nve` / `nvt` / `diffusion` | all metrics | **0** — bit-identical | — |

Every log column except pressure — step, time, PE, KE, total energy,
temperature and volume — is bit-identical across all 201 frames, so the
barostat accepted exactly the same sequence of volume moves. The pressure
metrics do **not** move by the constant's own 4.0918e-09 ratio, and the reason
is the log rather than the physics: at ~46 bar printed with six decimals one
printed unit is 1e-6 bar, while the conversion moves each value by ~1.9e-7. The
metric is a mean of quantized values. Frame by frame the effect is bounded and
fully accounted for: 62 of 201 frames changed, every one by exactly one unit in
the last printed place, 51 up (all with positive pressure) and 11 down (all with
negative pressure), no exceptions. The net shift those counts predict,
`(51−11)·1e-6/201 = 1.9900497512e-07` bar, matches the observed shift in the
mean, `1.9900497250e-07` bar, to 2.6e-15.

**Restart implications.** Checkpoints store pressure in eV/Å³, the internal
unit, so the conversion is **not** serialized and an old checkpoint resumes
under the corrected one. The stored number does not change; the bar value
reported from it moves by 4.09e-09. Existing logs and trajectories carry
`P[bar]` columns written with the superseded factor and are not comparable with
new ones below that level.

**What is asserted.** `tests/pressure_unit_tests.cpp` and
`tests/mpi_pressure_units.cpp` read no production constant — a test that
imported one would agree with a wrong one. The reference is derived from the SI
definitions as the rational `500/801088317` and cross-checked against the naive
floating-point route. Each path is then measured back out of the engine:
reporting by writing a frame with a known internal pressure and reading the bar
column; the Monte Carlo barostat by bisecting for the target pressure at which a
fixed-seed trial move flips from accepted to rejected, run twice at different
atom counts so that the subtraction cancels the unobservable uniform draw *and*
the Boltzmann factor exactly; Berendsen by inverting its observable coupling
factor to recover the pressure it believed it had. Omitting the conversion
(measures 1), applying it twice (`c²`), inverting it (`1/c`), reversing the
pressure-work sign and restoring the superseded literal are each named failures.
Under MPI the measurement is repeated at 1, 2 and 4 ranks with one empty rank,
since the pressure it converts is itself a reduced quantity.

### The Boltzmann constant

GMD carried **four** different values for one constant. They are now one:

```
gmd::kBoltzmannConstantEVPerKelvin = 8.617333262145177e-5   // eV/K
```

| Where | Was | Relative to the correct value |
|---|---|---|
| `src/system/initializer.cpp` | `8.617343e-5` | **+1.130031e-06** |
| `include/gmd/integrator/thermostat.hpp` | `8.617333262e-5` | −1.685e-11 |
| `include/gmd/integrator/mc_barostat.hpp` | `8.617333262e-5` | −1.685e-11 |
| `examples/ethane_demo/ethane_demo.cpp` | `8.617333e-5` | −3.042e-08 |

**This is a results-changing correction.** Velocity initialization moves by
−1.130030e-06, which changes every trajectory that starts from randomly sampled
velocities.

**The observable bug.** Velocity initialization and temperature reporting used
different constants, so the engine did not agree with itself: a system
initialised to 300 K reported **300.000339014 K**. It now reports 300 K to
2.2e-16. The Nosé–Hoover thermostat mass and the MC barostat's Metropolis
exponent were on the reporting side of the split; the initializer was the odd
one out.

**Derivation.** Unusually, `k_B` in eV/K is **exact** and carries no
uncertainty, because both ingredients have been exact by definition since the
2019 SI redefinition:

| | | |
|---|---|---|
| `k_B` | 1.380649e-23 J/K, exact | [physics.nist.gov/cgi-bin/cuu/Value?k](https://physics.nist.gov/cgi-bin/cuu/Value?k) |
| `e` | 1.602176634e-19 C, exact | [physics.nist.gov/cgi-bin/cuu/Value?e](https://physics.nist.gov/cgi-bin/cuu/Value?e) |

One electronvolt is `e` joules exactly, so

```
k_B [eV/K] = 1.380649e-23 / 1.602176634e-19
           = 1380649 / 16021766340
           = 8.6173332621451774336636593340806392...e-5
```

The denominator has factors other than 2 and 5, so this exact rational has a
**non-terminating** decimal expansion — every literal is a truncation of it.
NIST tabulates it as *"8.617 333 262… × 10⁻⁵ eV K⁻¹, exact"*
([CODATA 2022](https://physics.nist.gov/cgi-bin/cuu/Value?kev)); the ellipsis is
their notation for precisely that.

**Rounding policy.** With no physical uncertainty to hide behind, the only
defensible cut-off is double precision. The literal is the shortest decimal that
parses to the nearest `double` to the exact rational, 1.1e-18 away — below one
ulp. Fewer digits would be an arbitrary truncation.

**Affected production paths.** Maxwell–Boltzmann velocity sampling
(`σ = √(k_B T/m)`), the rescale-to-target that follows it, instantaneous
temperature reporting, the Nosé–Hoover thermostat mass `Q = dof·k_B·T·τ²` and
its friction force, the velocity-rescaling thermostat's target, and the MC
barostat's `β = 1/(k_B T)`.

**What is *not* affected.** Static energy, force and virial results are
untouched: no Coulomb, Lennard-Jones or bonded quantity depends on `k_B`. The
`static_lj_cluster`, `static_coulomb`, `static_special_pairs`,
`static_bonded_reference` and `pme_external` baselines are unchanged by
construction.

**Baseline impact.** Three of the four long dynamics baselines moved and were
regenerated; one provably did not.

| Case | Metric | Relative change | Tolerance |
|---|---|---|---|
| `nvt_lj_fluid` | temperature mean / stddev | +1.079e-07 / +1.912e-06 | 5.0 K |
| `npt_lj_fluid` | temperature mean / stddev | +3.306e-07 / +2.192e-06 | 10.0 K |
| `npt_lj_fluid` | pressure mean / stddev | −1.828e-05 / −1.174e-06 | 2000 bar |
| `diffusion_lj_fluid` | diffusion coefficient | +2.593e-04 | 0.05 |
| `nve_lj_fluid` | energy drift | **0** — bit-identical | 1e-06 |

Every one stayed inside its existing tolerance, so all would have passed without
regeneration; they were regenerated because they had been reproducing bit-for-bit
and no longer did. The diffusion coefficient moves three orders more than the
5.65e-07 velocity perturbation — Lyapunov amplification over the run, not an
error. The NVE metric is a difference of two energies printed at six decimals,
and the correction shifts the total energy by ~3e-08 eV on ~0.29 eV, two orders
below the last printed digit, so neither printed value moves.

**Restart implications.** The constant is not serialized by name anywhere. But
the Nosé–Hoover thermostat mass **is**, and `Q = dof·k_B·T·τ²` carries the
constant inside it. `load_checkpoint_state()` restores `Q` verbatim by design, so
that a continued run is deterministic — which means **a checkpoint written before
this change keeps the old constant's thermostat mass after a resume**, while
every other path uses the new one. The discrepancy is 1.7e-11 for a checkpoint
from a Nosé–Hoover run (the thermostat was already on the correct side of the
split). Restart from the input rather than the checkpoint if that matters.

**No bitwise backward compatibility.** Any run that samples initial velocities
diverges from a pre-correction run of the same input and seed.

**What is asserted.** `tests/boltzmann_constant_tests.cpp` and
`tests/mpi_boltzmann_constant.cpp` do not read a production constant — the
initializer's is file-local and the barostat's is private. Each path is driven
with known inputs and the constant recovered from its output: `k_B = 2K/(dof·T)`
after initialization, the inverse relation for reporting, `Q/(dof·T·τ²)` for
Nosé–Hoover, and for the MC barostat a bisection for the temperature at which
the first trial move's accept/reject decision flips, which is exactly inversely
proportional to its constant. Initialization is asserted to agree with reporting
directly, so two paths drifting together could not pass. Omitting the constant
(measures 1), applying it twice (`k_B²`) and inverting the conversion (`1/k_B`)
are each named failures. Under MPI the measurement is repeated at 1, 2 and 4
ranks, where one rank deliberately owns no atoms.

### The electrostatic constant

GMD's Coulomb constant was `14.3996` eV·Å/e² through v2.4. It is now

```
gmd::kCoulombConstant = 14.3996454686836   // include/gmd/core/physical_constants.hpp
```

a **results-changing correction of +3.157635e-06 relative**. Every Coulomb,
Ewald and PME energy, force and virial the engine produces changes by that
factor. Nothing is bitwise compatible with a previous release.

**Derivation.** In atomic units the Coulomb energy of two unit charges one Bohr
radius apart is exactly one Hartree, so `e²/(4πε₀) = E_h·a₀` and, in GMD's
units, `k_e[eV·Å/e²] = E_h[eV]·a₀[Å]`. Both factors are 2022 CODATA
recommended values from NIST:

| | | |
|---|---|---|
| `E_h` | 27.211386245981(30) eV | [physics.nist.gov/cgi-bin/cuu/Value?hrev](https://physics.nist.gov/cgi-bin/cuu/Value?hrev) |
| `a₀` | 5.29177210544(82) × 10⁻¹¹ m | [physics.nist.gov/cgi-bin/cuu/Value?bohrrada0](https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0) |
| product | **14.399645468683593…** eV·Å/e² | |

Deriving the same number straight from the elementary charge and the vacuum
permittivity — `k_e = e / (4πε₀ × 10⁻¹⁰)`, with `e = 1.602176634e-19` C exact
by the SI definition and `ε₀ = 8.8541878188(14)e-12` F/m — agrees to 1.1e-12
relative, which is the uncertainty in `ε₀` rather than an error in either route.

**Why it was worth changing.** 3.16e-06 is not a cosmetic digit here. It is
larger than the tightest PME convergence point this repository measures —
6.3e-10 relative, order 6 at grid 128 — so it was a floor on agreement with any
external code no matter how fine the mesh, and `validation/pme_external` had to
correct for it explicitly to keep its convergence study meaningful.

**Rounding policy.** Enough digits to reproduce the derivation to
double-precision representation (4.6e-16 relative). Not excess precision:
rounding to LAMMPS' six decimals would be 3.3e-08, still coarser than that
6.3e-10 convergence point, and would reintroduce a smaller version of the same
floor. Digits past the tenth decimal are below CODATA's own 1.5e-10
uncertainty and are carried only so the literal is exact as a `double`.

**Compatibility policy: CODATA, not an engine.** Nothing in this repository
documents a requirement to reproduce another code's constant, and the engines
do not agree with each other in any case — LAMMPS `units metal` uses
`14.399645`, OpenMM's `ONE_4PI_EPS0` works out to `14.399645478`, both
*measured* rather than quoted. `validation/pme_external` converts each engine's
result by the exactly linear ratio `k_e^GMD / k_e^engine` instead of adopting
either.

**One definition.** The constant previously existed as two independent literals,
in `ewald_force_provider.cpp` and `pme_force_provider.cpp`, with nothing keeping
them equal. Both now derive from `gmd::kCoulombConstant`. The test-side
references (`tests/virial_reference.hpp`, `tests/pme_reference.hpp`,
`tests/special_pair_tests.cpp`, `tests/mpi_special_pair.cpp`,
`validation/analytic_references.py`) each keep their **own** declaration on
purpose: importing the production symbol would turn every comparison built on
them into a restatement of production rather than a second opinion.

**What is asserted.** `tests/electrostatic_constant_tests.cpp` and
`tests/mpi_electrostatic_constant.cpp` do not read the production constant —
it has internal linkage. They recompute each electrostatic sum with `k_e = 1`
and recover `k_e_measured = E_provider / Φ_reference`, which is exact because
the dependence is exactly linear. Covered: Ewald real, reciprocal, self and
background; the same four for PME at orders 4 and 6; and the bare special-pair
correction. Energy, force and virial are recovered *separately*, so a missing
factor in the force path cannot hide behind a correct energy. A path that
dropped the constant would measure 1, one that applied it twice would measure
`k_e²`, and two paths with differently rounded copies would disagree — each is
a named failure. Plus two unit charges exactly 1 Å apart, where the energy, the
force magnitude and `W_xx` must each equal `k_e` outright; exact quadratic
scaling under a charge rescale; repulsion and attraction signs; and the same
measurement at 1, 2 and 4 MPI ranks, where the providers return rank-local
partials and a term applied once per rank instead of once per pair would stay
self-consistent within any single rank.

**Restart and baseline implications.** Checkpoints do not store the constant, so
a checkpoint written before this change and resumed after it will continue with
the corrected value; the trajectory diverges from the original run at the
3.16e-06 level. Every stored Coulomb baseline was regenerated — see
`validation/static_coulomb`, `validation/static_special_pairs` and
`validation/pme_external` — with the Coulomb components verified to scale by
exactly the constant ratio and the Lennard-Jones components verified
bit-identical. Totals and stored per-atom forces do **not** scale exactly,
because they are sums of a quantity that depends on `k_e` and one that does not.

### Independent external PME validation

Everything above compares GMD against references written for GMD. That is
sufficient to catch an implementation mistake and insufficient to catch a shared
misunderstanding. `validation/pme_external` closes that gap by running two
unrelated engines on the same configuration.

| | |
|---|---|
| Primary reference | **OpenMM 8.6.0**, `NonbondedForce` with `PME`, Reference platform (double precision), installed from the PyPI wheel |
| Second reference | **LAMMPS 22 Jul 2025 – Update 5**, `kspace_style pppm` and `kspace_style ewald`, `units metal`, Homebrew build |
| Fixture | `validation/pme_external/charges.xyz` — 12 atoms, exactly charge neutral, irregular positions, **non-cubic** 18 × 22 × 26 Å box, no symmetry that makes any force component vanish |
| Coupling | Coulomb only. No Lennard-Jones, no bonded terms, no exclusions, no 1-4 scaling — so there is no shift, switch, mixing rule or dispersion-correction convention to reconcile |
| Settings | α = 0.35 Å⁻¹, real-space cutoff 8.0 Å (< L_min/2), grids 16³ … 128³ |

**What is exactly aligned, and what is not.** This distinction is the whole
substance of a cross-code comparison, so it is stated per parameter rather than
asserted in aggregate:

| Parameter | GMD | OpenMM | LAMMPS | Aligned? |
|---|---|---|---|---|
| Splitting α | set directly | `setPMEParameters`, read back from the context as 3.5 nm⁻¹ | `kspace_modify gewald` | **exactly** |
| Reciprocal grid | set directly | `setPMEParameters`, read back | `kspace_modify mesh` | **exactly** |
| Real-space cutoff | 8.0 Å | 0.8 nm | 8.0 Å | **exactly** |
| Periodicity, box | orthorhombic, all three axes | same, off-diagonals exactly zero | same | **exactly** |
| Exclusions, net charge | none; system neutral | none; system neutral | none; system neutral | **exactly** |
| B-spline order | 4 or 6 | **fixed at 5**, not exposed by any API | `kspace_modify order`, 4 and 6 both run | **impossible with OpenMM** |
| Reciprocal influence function | Essmann smooth PME | Essmann smooth PME | PPPM *optimised* Green's function — a different mesh approximation by construction | **not aligned with LAMMPS** |
| Coulomb constant `k_e` | `14.3996454686836` eV·Å/e² (CODATA 2022) | 138.935457 kJ/mol·nm/e² → 14.399645478 | 14.399645 eV·Å/e² | **not aligned** — each engine rounds it differently; see below |

**The Coulomb constant is corrected exactly, not tolerated.** All three codes
use a *different* number for the same physical constant, because each rounds it
somewhere different: GMD `14.3996454686836` (CODATA 2022, see *The
electrostatic constant* below), LAMMPS `14.399645` (its own six decimals,
3.3e-08 low) and OpenMM `14.399645478` (6.8e-10 high). Every Ewald term carries
exactly one factor of `k_e`, so each engine's result is rescaled by
`k_e^GMD / k_e^engine`, which is an exact correction rather than an
approximation. Both engine constants are *measured* at generation time with a
two-charge probe rather than quoted from documentation. For OpenMM the eV ↔
kJ/mol convention cancels out of the combined factor entirely, so no
energy-unit constant enters the comparison.

Until GMD's constant was corrected this rescaling was doing much more work: the
old `14.3996` was 3.16e-06 low, which exceeded the grid-128 mesh error and
would have been misread as a convergence floor. What it compensates now is only
the engines' own rounding, two orders of magnitude smaller.

**Convergence, not a single tolerance.** All three codes are swept over the same
grids and measured against the exact Ewald sum at the same α and the same
cutoff. Absolute energy error in eV:

| Grid | GMD order 4 | GMD order 6 | LAMMPS PPPM order 4 | LAMMPS PPPM order 6 | OpenMM order 5 |
|---|---|---|---|---|---|
| 16³ | 1.54e-02 | 1.87e-03 | 4.51e-03 | 4.22e-04 | 4.68e-03 |
| 32³ | 8.04e-04 | 1.19e-05 | 1.17e-04 | 2.18e-06 | 4.76e-05 |
| 64³ | 5.96e-05 | 1.97e-07 | 6.92e-06 | 1.25e-07 | 9.34e-07 |
| 128³ | 3.54e-06 | 2.65e-09 | — | — | 1.17e-08 |

Every column falls monotonically, at the rate its interpolation order predicts
(GMD order 4 ≈ h⁴, order 6 ≈ h⁶). The three codes converge toward the *same*
limit from different mesh approximations, which is the actual claim being made.
LAMMPS PPPM stops at 64³ because `kspace_modify mesh N N N` aborts that build
with SIGSEGV for every N > 64 on this fixture, with pinned or auto-selected
parameters alike — an external-engine limit, recorded in the reference file.

**Result at the reference settings** (grid 64³, GMD order 6 vs OpenMM order 5):

| Quantity | Error | Where |
|---|---|---|
| Energy | 1.13e-06 eV absolute, 2.70e-07 relative | — |
| Force, worst component | 1.53e-06 eV/Å | atom 0, z |
| Force RMS | 4.46e-07 eV/Å | — |
| Against LAMMPS PPPM | 7.21e-08 eV, 2.69e-07 eV/Å | — |
| Virial, all six independent components vs LAMMPS PPPM | ≤ 8.58e-07 eV | `W_zz` |
| Virial vs LAMMPS exact Ewald | ≤ 1.65e-07 eV | `W_zz` |

The tolerances are derived, not observed-and-rounded: at these settings GMD's
own mesh error is 1.97e-07 eV and OpenMM's is 9.34e-07 eV, so two mesh methods
can differ by at most their sum, 1.13e-06 eV. The observed difference *is*
1.1313e-06 — the triangle-inequality bound is tight, because the two codes
deviate from exact Ewald in opposite directions. The committed tolerance is that
bound doubled.

**Also checked, at test time:** that the reference is not GMD's own output (two
different PME implementations at different interpolation orders cannot agree
below 1e-8, so a smaller gap is a replaced file, not a passing test); that the
fixture and run-file hashes match the ones the reference was generated for; that
the run file's α, cutoff and grid still match the reference block; that the same
configuration written in different periodic images gives identical energy
(1.7e-13 eV) and forces (1.9e-14 eV/Å); and that a rigid translation by a
non-lattice vector changes the answer only at the mesh-error level, since a mesh
fixed to the box makes exact translation invariance false by construction.

**What this does *not* establish.** OpenMM exposes no virial in any form, so the
tensor comparison rests on LAMMPS alone. Neither engine exposes a proven
equivalent of GMD's *decomposition* into real, reciprocal, self and background
parts, so only the total electrostatic energy, the forces and the total
electrostatic virial are cross-validated — not the internal split. Nothing here
covers LJ, bonded, special-pair or constraint contributions, and nothing here
runs under MPI.

Regenerating the reference requires both engines and is a deliberate manual
step; the CI comparison reads the checked-in JSON and needs neither. See
`validation/pme_external/README.md`.

### PME: the production kernel, tested directly

The virial coverage above validates a tensor. It does not, on its own, say
anything about whether the force the integrator receives is the gradient of the
energy the engine reports — and that is precisely where PME was broken. Three
further test files cover the production path itself:

| File | What it holds against what |
|---|---|
| `tests/pme_bspline_tests.cpp` | `PMEForceProvider::bspline()` and `bspline_deriv()` **directly** (they are public static members for this reason), at orders 2–6: every polynomial interval against the order recursion, exact knots and both sides of them, zero outside `[0,p]`, non-negativity, `M_p(u) = M_p(p-u)`, partition of unity over 1997 offsets, continuity through derivative `p-2`, the identity `M'_p = M_(p-1)(u) - M_(p-1)(u-1)`, a coordinate finite difference, and explicit regression values on `[2,3)`, `[3,4)`, `[4,5)` |
| `tests/pme_energy_force_tests.cpp` | the transform convention, against an **explicit direct DFT** with both directions unnormalised, on an 8³ mesh — which pins the `K1·K2·K3` factor absolutely rather than through convergence; and `F = -dU/dr` by central differences over orders 4/6, meshes 16/32, three alphas, neutral and net-charged, an asymmetric configuration, atoms on and across the faces, and lattice-translation invariance |
| `tests/ml_virial_contract_tests.cpp` | `MLForceProvider` reports no virial, clears a stale `ForceResult`, does not synthesise a force moment, invalidates any composite it sits in, and cannot feed a pressure-coupled barostat |
| `tests/electrostatic_constant_tests.cpp`, `tests/mpi_electrostatic_constant.cpp` | the Coulomb constant itself, measured out of the engine rather than read: each electrostatic sum is recomputed with `k_e = 1` and the constant recovered as `E_provider / Φ_reference`. Ewald real/reciprocal/self/background, the same four for PME at orders 4 and 6, and the bare special-pair correction; energy, force and virial recovered separately; two unit charges 1 Å apart; quadratic charge scaling; repulsion and attraction signs; and the same measurement at 1, 2 and 4 MPI ranks |

**Why convergence was not enough.** `tests/pme_reciprocal_virial_tests.cpp`
shows PME converging to Ewald as the mesh is refined, and that test passed
throughout the period when the reciprocal force was zero — because it compared
tensors and energies, and the tensor is computed before the inverse transform.
Convergence to another method also cannot detect an internally inconsistent
energy and force: both defects in the force path left the energy exactly right.
The direct tests above are absolute comparisons at a fixed mesh, so they fail
immediately rather than in the limit.

**Lattice translation is exact; arbitrary translation is not.** Moving every
atom by a whole number of box vectors describes the identical periodic system,
and the energy and forces come back unchanged to `1e-11`. An arbitrary
translation is *not* an exact symmetry of PME — it moves the charges relative
to the mesh and changes the B-spline aliasing error — so what the tests require
there is that the deviation shrinks with mesh refinement, which it does:
`1.1e-2 → 2.2e-3 → 3.9e-4` at order 4 and `1.2e-3 → 2.9e-5 → 1.2e-6` at order 6
over meshes 16, 32, 64.

### TorchScript `edge_shift`

The tensor an ML model receives for periodic edges. The contract, established by
reading the production path and then confirmed from inside a scripted model:

| | |
|---|---|
| `edge_index` | `int64 [2, 2E]`, row 0 = `src`, row 1 = `dst` |
| `edge_shift` | `float32 [2E, 3]`, a **Cartesian translation in Å** — not an integer image offset |
| device | whatever the loaded module is on |
| displacement | `d(src → dst) = r_dst + edge_shift − r_src`, equal to the minimum-image separation |
| reverse edge | every half-pair yields both directions; the two shifts are exact negatives |
| sign | a negative x-shift on `i → j` means j's relevant image sits one box length in −x from where j is stored |
| periodicity | all three axes, always; there is no per-axis flag to disable one |

`VerletNeighborBuilder` stores `S = round((r_i − r_j)/L)` per half-pair and the
adapter emits `S·L` for `i → j` and `−S·L` for `j → i`. The shift is a function
of the **stored** coordinates, so the same physical dimer written wrapped and
unwrapped produces different shifts and the same displacement — both halves are
asserted.

**What is covered, and where.** `tests/box_image_flag_tests.cpp` covers
`image_flags`, the helper. `tests/torchscript_edge_shift_tests.cpp` covers what
reaches `forward()`: the half-to-full-graph expansion, the reverse edge, the
box-length multiplication, the `src`/`dst` ordering, and the float32 conversion.
It builds its observer models from C++ at run time, so no Python and no
checked-in `.pt` are needed and no artifact lands in the source tree.

Serial fixtures: no crossing; each face in both sign directions; a corner with
three non-zero mixed-sign components; the reverse edge of every fixture;
wrapped versus unwrapped equivalence; every axis periodic; atom-order
permutation including swapping each pair. Composite: an ML child alone and
between two other providers, evaluated repeatedly, its contribution recovered by
subtraction. Shift components are compared **exactly** — the box lengths and
species are chosen so nothing in the chain rounds.

Negative controls run mutated observer models (shifts zeroed, sign reversed,
x/y/z permuted) and require every detectable edge to be rejected.

**MPI.** `edge_shift` never reaches a model under MPI: `MLForceProvider` throws
from both `initialize()` and `compute()` when the rank count exceeds one,
because nothing in the adapter decides which rank owns a cross-boundary edge's
energy or how deep the halo must be for a multi-layer model.
`tests/mpi_ml_edge_shift.cpp` proves that at np=1/2/4 rather than adding a
single-rank test dressed as parallel: at one rank the provider must work, at
more it must refuse on every rank, a composite must propagate the refusal, and
the adapter must never be reached.

**This is an interface-contract test, not a scientific validation.** It shows
the model receives the shifts the contract specifies. It says nothing about
whether any particular model's energies or forces are correct.

**Availability.** The test exists only when `GMD_ENABLE_TORCH=ON`. A build
without LibTorch prints a configure-time `STATUS` line saying the test is not
registered and the contract is unverified in that build.

### Constraint (SHAKE/RATTLE) virial

**The splitting.** Constrained dynamics uses the standard velocity-Verlet
SHAKE/RATTLE splitting, in which the constraint force enters **twice** per step:

```
(1)  v_i += (dt/2m_i)·F_i(t)                    first half-kick
(2)  v_i += w_i·Σ_c s_ic·Γ_c·r_c(t)             SHAKE impulse      } solved
(3)  r_i += dt·v_i                              drift              } together
(4)  v_i += (dt/2m_i)·F_i(t+dt)                 second half-kick
(5)  v_i += w_i·Σ_c s_ic·Λ_c·r_c(t+dt)          RATTLE impulse
```

(2) is the constraint force at time level `t`, paired with `F(t)`; (5) is the
same force at `t+dt`, paired with `F(t+dt)`. They are different quantities.

**Which one the pressure needs, and the factor.** The providers are evaluated at
`r(t+dt)`, so the provider virial belongs to `t+dt` and its partner is the
RATTLE impulse. Matching (5) against (4) term by term,

```
(dt/2m_i)·G_i(t+dt) = w_i·Λ_c·r_c(t+dt)   ⇒   G_i(t+dt) = (2/dt)·Σ_c s_ic·Λ_c·r_c
W_constraint(t+dt) = Σ_c (2Λ_c/dt)·(r_c ⊗ r_c)
```

The factor is **`2/dt`**, and it is exact — this is an endpoint quantity, not a
step average, so **both halves of the reported virial live at the same time
level and there is no `O(dt)` mixing**. Pinned from the outside rather than
assumed: a freely rotating rigid dimer's endpoint constraint force must equal
the centripetal force `μω²d`, and it does, to `O(dt²)`. `1/dt` — the value that
would be right if SHAKE omitted its velocity impulse (2), leaving RATTLE to
carry the whole step — undershoots it by exactly two.

**Central and symmetric by construction.** RATTLE does not move atoms, so
`r_c(t+dt)` is fixed for the whole solve and summing the *scalar* multipliers
gives a pair force exactly along `r_c` and exactly antisymmetric — for coupled
constraints such as a rigid triangle as much as for isolated ones. Nothing is
projected and no residual is discarded.

**Sign.** `λ = -(r_c·v_c)/((w_i+w_j)|r_c|²)`, so a separating pair gives `λ < 0`
and an attractive restoring force with a negative trace — matching the negative
virial an attractive Lennard-Jones pair produces under the same convention.

**Why the kinetic term does not already cover it.** `2K` and the configurational
constraint virial are different quantities, and projecting the velocities does
not remove the need for the second. For a rigid rotor of reduced mass `μ`,

```
2K = μω²d²        tr W_constraint = -μω²d²        2K + tr W = 0
```

which is the physically required answer — a rigid body's frozen internal
coordinate contributes nothing to the pressure. Drop the constraint virial and
`2K` is left uncancelled, reporting a spurious `2K/3V`.

### Constraint independence

A constraint set is rejected before any dynamics run unless it removes as many
degrees of freedom as it has constraints. **Distinct is not independent**: every
pair among five atoms is ten distinct constraints over a body with at most nine
internal degrees of freedom, and any three collinear atoms have three distinct
pair distances of which only two are free. Neither is visible to a duplicate
check.

**The criterion.** A constraint `σ_c = |r_c|² − d_c²` removes a degree of freedom
only if its gradient is independent of the others *in the metric the dynamics
uses*. The constrained equations of motion involve `J M⁻¹ Jᵀ`, so the matrix that
matters is the mass-weighted Jacobian

```
J_M = J M^(−1/2)      J_M[c, 3i+a] = +2 r_c[a] / √m_i
                      J_M[c, 3j+a] = −2 r_c[a] / √m_j
```

and the number of degrees of freedom removed is `rank(J_M)`, never the number of
constraints supplied. `J M⁻¹ Jᵀ` is singular exactly when `J_M` is rank
deficient, which is also when SHAKE and RATTLE have no unique multiplier to
converge to — so a dependent set is not merely a bookkeeping problem.

**Per component.** Constraints sharing no atom have Jacobian rows with disjoint
support, so rank is additive over connected components of the constraint graph.
The analysis works component by component, which keeps the linear algebra on
matrices the size of a molecule rather than of the system.

**Numerics.** Singular values come from a one-sided Jacobi SVD of the gradient
columns, not from the eigenvalues of `J_M J_Mᵀ`: squaring the matrix squares the
condition number and halves the digits in exactly the small singular values that
decide the rank. **Convergence is checked, not assumed** — after every sweep the
scale-free off-orthogonality residual `max_{p<q} |a_p·a_q| / (|a_p||a_q|)` must
fall below `cols · ε`, and reaching the sweep limit throws rather than returning a
rank derived from an unconverged decomposition. The threshold is dimension aware
because orthogonalising a later pair perturbs an earlier one by `O(ε)` each time,
so what a converged sweep can leave behind grows with the column count. The rank
tolerance is

```
tol = max(rows, cols) · ε · σ_max
```

Each factor earns its place. `σ_max` makes the test **relative**: the entries of
`J_M` are `2r/√m`, whose magnitude depends entirely on the unit system and on how
long the bonds are, so an absolute threshold would mean one thing for a 1 Å bond
and another for a 10 Å one. `ε` is the unit round-off, the floor below which a
computed singular value carries no information. `max(rows, cols)` accounts for
error accumulating over the `O(max(rows, cols))` operations behind each singular
value. This is the standard LAPACK rank convention.

**The geometry that decides is the projected one.** Rank is a property of the
configuration, and the configuration the dynamics start from is the one on the
constraint manifold, not the one supplied — which can differ by a lot. The
supplied geometry is analysed only for early diagnostics; the authoritative check
runs on the converged projection, and velocities are projected only once it has
passed. Degrees of freedom are computed from that accepted state.

**Degenerate targets are caught earlier still, from the targets alone.** For three
atoms all constrained to each other, the triangle inequality decides everything:
`c > a + b` admits no configuration at all, and `c = a + b` admits only collinear
ones, where three distances have two independent gradients. That is checked at
construction, before any geometry exists, because a set like it is rank deficient
at *every* geometry it can reach and the projection would simply stall trying.
How degenerate a thin triple is gets measured as a **length** — the triangle
height `h = √(2ab·slack/(a+b))` the targets imply, compared against the solver's
own distance tolerance — because the triangle-inequality slack itself is a
second-order quantity, `slack ≈ h²(a+b)/2ab`, and comparing it to a distance
tolerance would reject thin but perfectly resolvable triangles.

**The policy is to reject.** `require_independent()` throws, naming the connected
component, its atom tags, the rank against the count, the singular values and
whether the dependence is structural (more constraints than the rigid-body limit
`3K − 6`) or geometric (within the limit, but degenerate at this configuration).
It also names **which** rows are redundant, when a column-pivoted rank-revealing
factorisation agrees with the singular values about how many there are — which
member of a dependent group gets named is not unique, so the pivoted choice is
used because it is deterministic and independent of the order the constraints
were supplied in. When the two disagree the report says so and falls back to
listing the whole component. Silently substituting the rank would keep the run
going with non-unique multipliers and hide what is almost always a topology
error. After `require_independent()` returns, the constraint count **is** the
rank, so the DOF subtraction is exact.

A set that is independent but close to degenerate is **accepted with a warning**
— rejecting a valid geometry for being awkward would be worse than the problem.
The threshold is `σ_max/σ_min > 1/√ε ≈ 6.7×10⁷`, the point at which roughly half
the significant digits of the multipliers are lost.

Detected and reported: exact duplicates, the same pair listed in reverse order
*relative to its first occurrence*, tolerance-equivalent duplicates, conflicting
targets, infeasible and degenerate target triples, linearly dependent
constraints, over-constrained components, degenerate (zero-length) gradients,
non-finite Jacobian entries and ill-conditioned sets.

The SVD kernel is exposed as `jacobi_singular_values()` so it can be validated
against matrices with independently known singular values — built as `A = U S Vᵀ`
from orthonormal factors, and cross-checked against a long-double Gram-matrix
eigensolver — rather than only against itself.

### Reported pressure vs current-geometry virial

These are two different things and the engine stores them separately.

- **`System::last_virial()`** is the *current geometry's* tensor: provider plus,
  when one belongs to that geometry, the constraint term. After a barostat
  rescale the forces are re-evaluated for a geometry no dynamics integrated, and
  that tensor has no constraint partner — it is then marked **invalid** rather
  than passed off as a complete pressure virial.
- **`System::step_thermodynamics()`** is the state *of the step that just
  finished* — pressure, virial, `2K`, volume and potential energy — captured
  together at the end of `finish_step()`, before any barostat runs. Its pressure
  is bit-for-bit the number the barostat itself consumed, and a later force
  evaluation deliberately does not touch the record.

The log row is that record: **`PE`, `KE`, `E_total`, `T`, `P` and `V` all come
from one state**, so `P = (2K + tr W)/3V` holds among exactly those numbers and
`E_total` is the energy of the configuration that produced them. Keeping the
completed step and the current geometry in one slot would force a rescale to
destroy the pressure the step actually measured, which is why constrained NPT
previously had no pressure to report.

The one thing that is *not* from that state is the **`.xyz` configuration**:
those are the current coordinates, which under a barostat are the rescaled ones
the next step starts from. Snapshotting coordinates every step to match would
cost an `N x 3` copy per step for a difference no thermodynamic quantity in the
row depends on. Under a barostat, therefore, the row's `V` lags the `.xyz`
coordinates by one rescale, by design.

**Validity is explicit, never a sentinel.** A numerical zero is an ordinary
pressure. When no complete pressure exists — the initial frame of a constrained
run, before any step has produced a RATTLE multiplier — the log writes `nan` and
sets the `P_valid` column to `0`.

**MPI output.** The MPI output path writes a separate `System` holding the
gathered global coordinates. It takes part in no dynamics, so every non-atomic
field the writer reads is brought across from the distributed `System` by
`System::copy_frame_state_from()` — box, potential energy, completed-step
record, SHAKE/RATTLE diagnostics and the virial state used as the initial-frame
fallback. Nothing there is reduced: each field is either replicated (box,
constraint virial) or already global (provider virial, and the completed-step
record built from globally reduced quantities), so every rank already holds the
identical value and reducing would multiply it. That identity is asserted
directly in `tests/mpi_constraint_virial.cpp` rather than assumed, and
`tests/constrained_pressure_reporting.py` compares whole logs from serial,
`np=2` and `np=4` runs of the real executable column by column.

**MPI constraint solving.** `ConstraintSolver` allgathers every owned atom and works from a
constraint list replicated against global atom tags, so every rank converges to
the same multipliers and the tensor each rank computes is *already global*. It
is therefore **not** reduced — an allreduce would multiply it by the rank count.
`tests/mpi_constraint_virial.cpp` runs the same fixture at 1, 2 and 4 ranks,
which catches both a dropped and a rank-multiplied contribution.

**Validated scope, all nine components.** Against the centripetal force of a
rigid rotor — a reference independent of the implementation — for a dimer and for
a rigid triangle of three coupled constraints, both in a **general 3D
orientation** so that no tensor entry is trivially zero. For the tilted dimer the
reference is `W_ab = −μω²d² û_a û_b` with `û` the endpoint bond direction; for the
tilted triangle it is `W = −ω² Σ_i m_i a_i ⊗ a_i`, the body's second-moment
tensor, both written down from rigid-body mechanics without touching the
accumulation under test. Every component, diagonal and off-diagonal, converges to
its reference as `dt²`. The tensor is also checked to transform covariantly,
`W → R W Rᵀ`, under a rotation with no special relationship to the axes. Plus the
magnitude of the recovered endpoint force, the sign, the trace, the identity
`2K + tr W = 0`, and second-order convergence of both errors in `dt`. Plus exact tensor symmetry, exact centrality
of every pair force, the two constraint half-impulses being measurably distinct,
energy conservation of a free rotor, per-constraint additivity, translation
invariance, invariance across a periodic boundary, rank invariance, and restart
continuity of the reported pressure. Constrained pressure is complete and
analytically validated for its **trace and diagonal** on orthorhombic cells. It
has **no independent external reference**: LAMMPS does not expose a separately
extractable constraint virial with a proven-equivalent definition, and the
quantity is dynamical, so it cannot be compared from a static configuration.

None of these depend on the coordinate origin or on how atoms happen to be
wrapped into the cell. In particular the reciprocal sum is **not** computed as
`Σ_i r_i ⊗ F_i`: the reciprocal energy depends on the cell explicitly, through
both the `1/V` prefactor and `k = 2πn/L`, so that expression is not the virial
at all there.

All of this is checked against the definition rather than against a closed form,
in two ways:

- **Isotropic.** `tr(W)` against `-dU/ds` under `r → s·r, L → s·L`.
- **Per-axis.** Each diagonal component `W_aa` against `-dU/dε` under a normal
  strain on that axis alone (`L_a → (1+ε)L_a`, fractional coordinates fixed),
  on deliberately non-cubic cells so no error hides behind cubic symmetry.

`tests/virial_finite_difference_tests.cpp` runs both for LJ, bonded (intact and
wrapped across a boundary), Ewald and PME, neutral and net-charged.
`tests/mpi_dof_virial.cpp` repeats both under 4-way domain decomposition, where a
virial reduced the wrong number of times fails by a factor of the rank count.

**Scope of the claim.** `Box` stores three edge lengths, so the engine
represents orthorhombic cells only and no shear strain can be applied.
**Only the three diagonal components are validated. The off-diagonal components
are not.** They are checked for the symmetry `W_ab == W_ba`, which is necessary
but not sufficient. Validating them would need a triclinic box representation
and a shear deformation, neither of which exists here.

The per-axis check earns its keep: it is what exposed the dihedral force bug
listed above, which the isotropic trace could not see because a torsion angle is
unchanged by isotropic scaling and so contributes zero to `tr(W)` either way.

---

## What's New in v2.3

| Feature | Files |
|---|---|
| **Directory restructuring** — merged 12 subdirectories into 6 for better cohesion; no code logic changed | `include/gmd/` `src/` `CMakeLists.txt` |
| **Include path cleanup** — old `gmd/boundary`, `gmd/neighbor`, `gmd/runtime`, and `gmd/ml` paths moved under `system`, `core`, and `force` | `include/gmd/` `src/` |

---

## What's New in v2.2

| Feature | Files |
|---|---|
| **3D MPI domain decomposition** — full 3D process grids with face/edge/corner ghost exchange, 3D periodic wraparound, `MPI_Alltoallv` reverse force accumulation, `MPI_Allgatherv` atom redistribution, and per-dimension periodicity control | `include/gmd/parallel/domain_decomposition.hpp` `include/gmd/parallel/mpi_communicator.hpp` `src/parallel/` |
| **`--proc-grid Px Py Pz` CLI flag** and **`mpi_grid` run.in directive** — user-selectable 3D process grids validated against MPI world size | `app/gmd_main.cpp` `include/gmd/io/config_loader.hpp` |
| **Comprehensive MPI test suite** — 17+ CTest targets covering 1D periodic ghost exchange/force/migration, 3D face/edge/corner ghost exchange, 3D reverse force, 8-rank LJ/Ewald/PME consistency, and 2-rank smoke tests | `tests/mpi_periodic_1d.cpp` `tests/mpi_ghost_exchange_3d.cpp` `tests/mpi_domain_decomposition_3d.cpp` `CMakeLists.txt` |
| **3D ghost image flags in Verlet lists** — `S[dim]` computed for all three dimensions, fixing neighbor images in 2D/3D process grids | `src/system/verlet_neighbor_builder.cpp` |
| **Periodic boundary migration test data** — atoms placed to exercise wraparound ghost exchange and redistribution | `tests/smoke_mpi_periodic_lj.{xyz,run}` `tests/smoke_mpi_boundary_migration.{xyz,run}` |

---

## What's New in v2.1

| Feature | Files |
|---|---|
| **MPI spatial domain decomposition** — balanced or user-selected 3D process grids, ghost-atom coordinate exchange, reverse force accumulation, migration, and allreduce helpers | `include/gmd/parallel/domain_decomposition.hpp` `include/gmd/parallel/mpi_communicator.hpp` `include/gmd/parallel/mpi_environment.hpp` `src/parallel/` |
| **PME execution-mode interface** — replicated PME remains the numerical backend; `pme_mode distributed` is a distributed-PME interface/prototype path for workflow testing and future FFTW-MPI integration, not a completed distributed mesh/FFT implementation | `include/gmd/force/pme_force_provider.hpp` `src/force/pme_force_provider.cpp` `CMakeLists.txt` |
| **`GMD_ENABLE_MPI` CMake option** — opt-in MPI build with `find_package(MPI REQUIRED)`, `MPI::MPI_CXX` linkage, and `GMD_ENABLE_MPI` compile definition | `CMakeLists.txt` `cmake/MPIOptions.cmake` |
| **MPI smoke tests** — multi-rank LJ/Ewald/PME runs and serial-vs-parallel logged energy consistency checks | `tests/smoke_mpi_lj.xyz` `tests/compare_energy_logs.cpp` `cmake/RunMpiLJConsistency.cmake` `cmake/RunMpiConsistency.cmake` |
| **CLI `--np N` flag** — validates expected process count against MPI world size | `app/gmd_main.cpp` |
| **Rank-aware I/O** — only rank 0 writes trajectory/log files; global coordinate gather for output frames | `app/gmd_main.cpp` |

---

## What's New in v2.0

| Feature | Files |
|---|---|
| **ML force provider integration** via TorchScript (LibTorch) — loads `.pt` models, reads `local_cutoff`, calls `forward(species, positions, edge_index, edge_shift)` | `include/gmd/force/torchscript_adapter.hpp` `include/gmd/force/ml_force_provider.hpp` `src/force/` |
| **`edge_index` and `edge_shift` tensors** — `NeighborList::image_flags` stores per-pair integer shift vectors; `VerletNeighborBuilder` computes them | `include/gmd/system/system.hpp` `src/system/verlet_neighbor_builder.cpp` |
| **Atomic numbers (`Z`)** — populated per atom during XYZ loading via built-in element table; forwarded to ML models | `include/gmd/system/system.hpp` `src/io/config_loader.cpp` |
| **`GMD_ENABLE_TORCH` CMake option** — opt-in LibTorch linkage; `force_field ml` and `model_path` directives in run files | `CMakeLists.txt` `app/gmd_main.cpp` |
| **PME smoke test** — validates Particle-Mesh Ewald path | `tests/smoke_pme.run` |

---

## What's New in v1.1

| Feature | Files |
|---|---|
| **Selectable LJ mixing rules** — `mixing_rule` / `mixing` now works for inline LJ in `run.in`, standalone LJ `.ff`, and molecular `.ff` files | `include/gmd/io/config_loader.hpp` `src/io/config_loader.cpp` |
| **Three supported combining rules** — `lorentz_berthelot` (default), `geometric`, and `waldman_hagler` for cross-type LJ pairs without explicit overrides | same files above |
| **Rule aliases and flexible spelling** — accepts `lb`, `geom`, `wh`, plus hyphenated names such as `lorentz-berthelot` and `waldman-hagler` | same files above |
| **Clearer non-bonded documentation** — examples and input format docs now describe `mixing_rule` and note that explicit `pair` entries override the configured fallback rule | `README.md` `examples/ethane_demo/ethane.ff` |

---

## What's New in v1.0

| Feature | Files |
|---|---|
| **Bonded force provider** — harmonic bonds, harmonic angles, periodic dihedrals, harmonic impropers | `include/gmd/force/bonded_params.hpp` `include/gmd/force/bonded_force_provider.hpp` `src/force/bonded_force_provider.cpp` |
| **Topology** — per-molecule connectivity (`BondTerm`, `AngleTerm`, `DihedralTerm`, `ImproperTerm`) | `include/gmd/system/topology.hpp` |
| **Molecular force field loader** (`load_molecular_ff`) — atom types with mass/ε/σ/q, explicit pair overrides, bonded parameter tables, inline comment support | `include/gmd/io/config_loader.hpp` `src/io/config_loader.cpp` |
| **Topology loader** (`load_topology`) — reads `.top` files with bond/angle/dihedral/improper sections | same files above |
| **Ethane demo** — complete working example with `.ff`, `.top`, `.xyz`, `.run` and a standalone build script | `examples/ethane_demo/` `run_ethane_demo.sh` |

---

## Current Scope

Implemented:

- periodic boundary conditions and minimum-image convention
- cell-list Verlet neighbor lists with skin-distance rebuild check and per-pair image shift vectors
- shifted Lennard-Jones pair interactions (multi-element, selectable mixing rules, explicit pair overrides that take priority over mixing); the static LJ cluster has an analytic validation reference
- **harmonic bonds** (`V = k(r − r₀)²`)
- **harmonic angles** (`V = k(θ − θ₀)²`)
- **periodic proper dihedrals** (`V = k[1 + cos(nφ − δ)]`)
- **harmonic improper dihedrals** (`V = k(φ − φ₀)²`)
- topology-derived **special pairs**: 1-2/1-3 LJ and Coulomb exclusions plus configurable 1-4 scaling; the four-atom chain validation uses an analytic reference
- Velocity Verlet integration
- velocity-rescaling thermostat
- Nosé-Hoover thermostat (VVNH splitting)
- Berendsen barostat (weak-coupling, requires virial)
- **Monte Carlo barostat** (isotropic NPT, Metropolis criterion, no virial required, adaptive step-size)
- long-range Coulomb via Ewald and replicated PME (self-contained 3D FFT, B-spline orders 4/6); Ewald has an analytic static validation reference, and replicated PME has an **independent external reference** — OpenMM PME with LAMMPS PPPM as a second engine, see *Independent external PME validation* above
- **ML force provider** — TorchScript backend (`GMD_ENABLE_TORCH=ON`), SE3-GNN compatible, reads `local_cutoff` from model
- **MPI spatial domain decomposition** — balanced or user-selected 3D process grids, face/edge/corner ghost-atom exchange, reverse force accumulation via `MPI_Alltoallv`, atom redistribution via `MPI_Allgatherv`, and allreduce collectives (`GMD_ENABLE_MPI=ON`); covered by 2/4/8-rank tests
- replicated-mesh PME under MPI; the `pme_mode distributed` path currently exercises an interface/prototype while still using the replicated numerical backend
- SHAKE / RATTLE bond constraints; serial constraints are tested, and MPI constraints currently use a correctness-first global-gather projection rather than a scalable owner-based solver
- checkpoint / restart through the CLI; restart continuity is tested by comparing continuous 100-step runs against 50-step + checkpoint + restart + 50-step runs in serial and MPI
- extended XYZ trajectory output and energy logging (rank-0 I/O with global coordinate gather in MPI mode)
- inline `run.in` LJ force-field definitions
- external `.ff` files (LJ and molecular; molecular runs apply topology special pairs by default)
- external `.top` topology files
- atomic numbers (`Z`) populated from element symbols during XYZ loading

Not yet implemented:

- CUDA compute backend (CMake option present but CPU-only)
- Python bindings
- full unit / regression test coverage
- true distributed PME/FFT: distributed charge assignment, grid decomposition, MPI FFT communication, and force interpolation from a distributed reciprocal grid

---

## Requirements

- C++ compiler with C++20 support (GCC 11+, Clang 14+, MSVC 2019+)
- CMake 3.24+

---

## Build

```bash
cd /path/to/GMD
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The main executable is produced at `build/gmd`.

Available CMake options:

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` / `Debug` |
| `GMD_ENABLE_CUDA` | `OFF` | Enable CUDA language (runtime still CPU-only) |
| `GMD_BUILD_PYTHON` | `OFF` | Python bindings stub |
| `GMD_ENABLE_TORCH` | `OFF` | Enable TorchScript ML force provider (requires LibTorch) |
| `GMD_ENABLE_MPI` | `OFF` | Enable MPI spatial domain decomposition (requires MPI) |

MPI build example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGMD_ENABLE_MPI=ON
cmake --build build --parallel
```

Torch + MPI build example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DGMD_ENABLE_TORCH=ON -DCMAKE_PREFIX_PATH=/path/to/libtorch \
      -DGMD_ENABLE_MPI=ON
cmake --build build --parallel
```

---

## Quick Start — Lennard-Jones Smoke Case

```bash
./build/gmd tests/smoke_lj.xyz tests/smoke_lj.run
```

Output files:

- `output.xyz` — extended XYZ trajectory
- `output.log` — tabular thermodynamic log with `PE`, `KE`, `E_total`, `T`, `P`, `V`

The write cadence defaults to every 100 steps and can be overridden with:

```ini
output_interval 10
```

---

## MPI Parallel Run

Build with MPI enabled:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGMD_ENABLE_MPI=ON
cmake --build build --parallel
```

Run with multiple processes:

```bash
mpirun -np 4 ./build/gmd tests/smoke_mpi_lj.xyz tests/smoke_mpi_lj.run
```

The `--np N` CLI flag validates that the MPI world size matches expectations:

```bash
mpirun -np 4 ./build/gmd --np 4 input.xyz run.in
```

In MPI mode:
- Atoms are distributed across ranks via a balanced 3D process grid (computed by `MPI_Dims_create` or a balanced factorisation); `--proc-grid Px Py Pz` or `mpi_grid Px Py Pz` can request a compatible grid
- Atoms that cross domain boundaries are redistributed to their owning rank after each drift step
- Ghost atoms are exchanged across face, edge, and corner neighbors each step based on the neighbor-list cutoff + skin
- Forces on ghost atoms are reverse-accumulated back to their home ranks via `MPI_Alltoallv`
- Only rank 0 writes trajectory and energy log files
- Global coordinates are gathered for output frames
- Ewald and PME run with the 3D domain decomposition. PME currently sums rank-local charge assignment onto a replicated full mesh and performs a full FFT on every rank; it is not distributed PME.
- TorchScript ML force fields and the Monte Carlo barostat are rejected in multi-rank runs because their distributed ownership/trial-move contracts are not implemented.

---

## Ethane Demo — Molecular Force Field

A fully self-contained example in `examples/ethane_demo/` demonstrating the bonded force field:

```bash
cd /path/to/GMD
bash run_ethane_demo.sh
```

The script compiles all sources with `c++` (no CMake required), copies the input files into `build_demo/`, runs the simulation, and prints a per-step energy table.

The main CLI also supports molecular runs directly. External molecular FF
runs enable LJ with topology-derived special pairs by default: bonded 1-2
and angled 1-3 pairs are excluded, while dihedral endpoint 1-4 pairs use
configurable scaling.

```bash
./build/gmd \
  examples/ethane_demo/ethane.xyz \
  examples/ethane_demo/ethane.run \
  examples/ethane_demo/ethane.ff \
  examples/ethane_demo/ethane.top
```

To make that behavior explicit, or to customize OPLS-style 1-4 factors, add:

```ini
molecular_nonbonded special
lj_scale_14 0.5
coul_scale_14 0.833333333333
```

`molecular_nonbonded none` remains available for bonded-only diagnostics.

### SHAKE / RATTLE Constraints

Molecular runs can constrain selected bond lengths with SHAKE after the
Velocity-Verlet position update and RATTLE after the final velocity update:

```ini
constraints on
constraint_tolerance 1e-6
constraint_max_iterations 100
rattle on
constrain_bond_type 1
```

`constrain_bond_type N` uses the `r0` value from `bond_type N` as the target
distance. Multiple constrained bond types can be listed with repeated
`constrain_bond_type` lines or one `constrain_bond_types 1 2 ...` line.
Topology files may also include explicit constraints:

```text
constraints 2
  1 2 0.9572
  1 3 0.9572
```

The thermodynamic log appends SHAKE/RATTLE iteration counts and maximum
constraint errors for every written frame. In MPI runs, constraints are applied
with a correctness-first global gather over owned atoms keyed by stable atom
tags, so constraints may cross domain boundaries. This is intentionally simple;
it is not yet the scalable local constraint-communication algorithm needed for
very large constrained biomolecular systems.

### Checkpoint / Restart

Long simulations can write a restartable text checkpoint that is separate from
the normal `.xyz` trajectory and `.log` thermodynamic output:

```ini
write_checkpoint_every 1000
checkpoint_file state.gmdchk
```

The checkpoint starts with `GMD_CHECKPOINT 1` and stores the current step,
simulation time, box, atom tags/types/molecule ids/masses/charges, positions,
velocities, topology entries, force-field/topology file references, boundary
summary, run/config summary, thermostat state, barostat state, and RNG state
where the active module owns one. Neighbor lists, ghost atoms, PME grids,
forces, and domain decomposition are intentionally rebuilt after restart.

To continue from a checkpoint, use a run input containing:

```ini
restart_from state.gmdchk
run 50000
```

`run` means "additional steps" after the checkpoint step. Restart writes a
fresh `output.xyz` and `output.log` with continuous step/time values beginning
at the checkpoint step; it does not append to an existing trajectory. MPI
restart reads the same global checkpoint on every rank, then rebuilds domain
decomposition for the current MPI size. This supports changing the process
count for the currently replicated text checkpoint format, assuming the same
force-field/topology inputs are available.

Expected output (truncated):

```
=== Loading inputs ===
  FF file       : ethane.ff
  Atom types    : 2
  Bond types    : 2
  Angle types   : 3
  Dihedral types: 1
  Bonds         : 7
  Angles        : 12
  Dihedrals     : 9

=== Building force providers ===
  BondedForceProvider ready
  LJ non-bonded: enabled (1-2/1-3 excluded, 1-4 scaled)

=== Initial state ===
  Bond PE     = 0.011972 eV  (0.276080 kcal/mol)

=== Running MD  (500 steps) ===
Step    Time(fs)    PE(kcal/mol)    KE(kcal/mol)    E_tot(kcal/mol)  T(K)
------------------------------------------------------------------------------
0       0.000       0.276           6.260           6.536           300.000
10      10.000      1.603           5.057           6.660           242.346
...
500     500.000     4.119           3.478           7.597           166.691

=== Done ===
  Frames written : 11
  Final PE       : 0.178612 eV  (4.118959 kcal/mol)
```

> **Note:** The standalone ethane demo driver may choose its provider set
> independently. The main CLI molecular workflow now enables non-bonded LJ
> safely through its topology-derived special-pair table.

---

## Input Files

GMD accepts up to four plain-text input files. Lines starting with `#` and inline comments (everything after `#` on a line) are ignored.

### `xyz.in` — initial structure

```text
<N>
<Lx> <Ly> <Lz>
<atom-line> ...
```

Atom-line formats (`type` is 1-based; element symbol also accepted):

| Columns | Content |
|---|---|
| `type x y z` | mass from FF or built-in element table |
| `type x y z mass` | explicit mass |
| `type x y z vx vy vz` | velocities from file |
| `type x y z vx vy vz mass` | both |
| replace `type` with `Ar`, `C`, `H`, … | same variants |

### `run.in` — simulation parameters

Core directives:

| Directive | Meaning |
|---|---|
| `run <N>` | number of MD steps |
| `time_step <dt>` | time step [fs] |
| `velocity <T>` | target / initial temperature [K] |
| `velocity_init random\|input` | randomize or read from xyz |
| `velocity_seed <s>` | RNG seed |
| `remove_com_velocity true\|false` | remove center-of-mass drift |
| `molecular_nonbonded special\|none` | molecular LJ mode; default `special` (`lj_unsafe` accepted as a legacy alias) |
| `lj_scale_12`, `coul_scale_12` | 1-2 pair factors; default `0.0`, `0.0` |
| `lj_scale_13`, `coul_scale_13` | 1-3 pair factors; default `0.0`, `0.0` |
| `lj_scale_14`, `coul_scale_14` | 1-4 pair factors; default `0.5`, `0.833333333333` |

Inline LJ force field:

```ini
force_field lj
cutoff 8.5
mixing_rule lorentz_berthelot   # or geometric / waldman_hagler
type 1 Ar epsilon 0.01032 sigma 3.405 charge 0.0
type 2 Ne epsilon 0.00312 sigma 2.749
```

Supported mixing-rule spellings:

- `lorentz_berthelot`, `lorentz-berthelot`, `lb`
- `geometric`, `geom`
- `waldman_hagler`, `waldman-hagler`, `wh`

Thermostat:

```ini
thermostat nose_hoover
thermostat_tau 100.0    # coupling time [fs]
```

or `thermostat velocity_rescaling`.

Barostat:

```ini
barostat berendsen
pressure 1.0
barostat_tau 2000.0
compressibility 4.5e-5
```

or `barostat monte_carlo`:

```ini
barostat monte_carlo
pressure 1.0
mc_frequency 25
mc_volume_step 0.01
```

Long-range Coulomb (Ewald or PME):

```ini
coulomb ewald
ewald_alpha 0.3
ewald_kmax 7
ewald_cutoff 10.0
```

```ini
coulomb pme
pme_alpha 0.3
pme_cutoff 10.0
pme_order 4
pme_grid 32 32 32
pme_mode replicated
pme_benchmark off
```

PME modes:

```ini
pme_mode replicated    # default, full PME mesh and FFT on every MPI rank
pme_mode distributed   # interface/prototype path; still uses replicated PME numerics
pme_mode auto          # reserved for future automatic backend choice
pme_benchmark on       # print charge/FFT/Green/interpolation timings
```

The distributed-PME prototype interface is guarded by:

```bash
cmake -S . -B build-mpi-pme \
  -DGMD_ENABLE_MPI=ON \
  -DGMD_ENABLE_DISTRIBUTED_PME=ON
```

`GMD_ENABLE_DISTRIBUTED_PME=ON` performs FFTW-MPI dependency checks
(`fftw3-mpi.h`, `libfftw3_mpi`, `libfftw3`) and fails at configure time if
they are missing. The default build keeps replicated PME and does not require
FFTW. In the current implementation, requesting `pme_mode distributed` selects
the distributed-PME interface/prototype path and timing hooks, but the
numerical work is still performed by the replicated PME backend and a warning
is printed. It is useful for exercising input/configuration, CMake dependency
checks, tests, and workflow compatibility; it does not reduce PME mesh memory
per rank and should not be used as evidence of scalable distributed-PME
performance.

Future true distributed PME still needs distributed charge assignment,
distributed grid ownership, an MPI FFT backend using slab or pencil
decomposition, reciprocal-space Green's-function application on distributed
mesh data, and force interpolation back from the distributed reciprocal grid
to domain-decomposed atoms.

ML force field (requires `GMD_ENABLE_TORCH=ON`):

```ini
force_field ml
model_path /path/to/model.pt
```

### `molecule.ff` — molecular force field

Used with `ConfigLoader::load_molecular_ff()` for bonded systems.  
All `kcal/mol` and degree values are **converted automatically** to internal `eV / rad` units.

```ini
force_field  molecular

mixing_rule  lorentz_berthelot   # optional; default is lorentz_berthelot
lj_cutoff  10.0

# Atom types: mass [amu], epsilon [kcal/mol], sigma [Å], charge [e]
type 1  C  mass 12.011  epsilon 0.0660  sigma 3.500  charge -0.180
type 2  H  mass  1.008  epsilon 0.0300  sigma 2.500  charge  0.060

# Explicit cross-pair override (overrides mixing_rule fallback)
pair 1 2  epsilon 0.0447  sigma 3.000

# Bond types: k [kcal/mol/Å²], r0 [Å]
bond_type 1  k 310.0  r0 1.526   # C-C
bond_type 2  k 340.0  r0 1.090   # C-H

# Angle types: k [kcal/mol/rad²], theta0 [deg]
angle_type 1  k 40.0  theta0 112.7   # C-C-C
angle_type 2  k 50.0  theta0 110.7   # H-C-H
angle_type 3  k 50.0  theta0 110.7   # H-C-C

# Dihedral types: k [kcal/mol], n (periodicity), delta [deg]
dihedral_type 1  k 0.150  n 3  delta 0.0

# Improper types: k [kcal/mol/rad²], phi0 [deg]
improper_type 1  k 10.5  phi0 0.0
```

### `molecule.top` — topology

Used with `ConfigLoader::load_topology()`. All indices are 1-based.

```text
bonds  7
1 2  bond_type 2   # C1-H1a
1 5  bond_type 1   # C1-C2

angles  3
2 1 5  angle_type 3   # j = vertex atom

dihedrals  9
2 1 5 6  dihedral_type 1

impropers  0
```

---

## Combining Force Providers

`CompositeForceProvider` accumulates energy and forces from an ordered list of providers:

```cpp
gmd::ConfigLoader loader;

auto mff  = loader.load_molecular_ff("molecule.ff");
auto topo = loader.load_topology("molecule.top");

auto bonded = std::make_shared<gmd::BondedForceProvider>(topo);
for (auto& bp : mff.bond_types)     bonded->add_bond_type(bp);
for (auto& ap : mff.angle_types)    bonded->add_angle_type(ap);
for (auto& dp : mff.dihedral_types) bonded->add_dihedral_type(dp);
for (auto& ip : mff.improper_types) bonded->add_improper_type(ip);

system.set_special_pair_map(std::make_shared<gmd::SpecialPairMap>(*topo));
auto lj = std::make_shared<gmd::ClassicalForceProvider>(mff.lj);
auto composite = std::make_shared<gmd::CompositeForceProvider>();
composite->add(lj); composite->add(bonded);
simulation.set_force_provider(composite);
// Total PE = scaled/excluded LJ + bond + angle + dihedral + improper
```

---

## Bonded Potential Functions

| Term | Potential | Parameters |
|---|---|---|
| Bond | $V = k(r - r_0)^2$ | `k` [eV/Å²], `r0` [Å] |
| Angle | $V = k(\theta - \theta_0)^2$ | `k` [eV/rad²], `theta0` [rad] |
| Dihedral | $V = k[1 + \cos(n\varphi - \delta)]$ | `k` [eV], `n` (int), `delta` [rad] |
| Improper | $V = k(\varphi - \varphi_0)^2$ | `k` [eV/rad²], `phi0` [rad] |

Forces are derived analytically. Dihedral and improper forces use the Blondel–Karplus gradient distribution to all four atoms.

Unit conversion helpers available in `bonded_force_provider.hpp`:

```cpp
gmd::kcal_per_mol_to_eV   // = 4.336410e-2
gmd::deg_to_rad           // = π / 180
```

---

## Testing

```bash
# Non-MPI tests
cd build
ctest --output-on-failure

# MPI tests (requires GMD_ENABLE_MPI=ON)
ctest --output-on-failure -L mpi

# All tests
ctest --output-on-failure
```

When `GMD_ENABLE_MPI=ON`, CTest additionally registers MPI smoke tests for LJ,
Ewald, PME, and the minimal bonded molecular path, plus serial-vs-MPI LJ,
Ewald, and replicated-PME consistency checks. The consistency checks validate
the full logged trajectory row-by-row as well as the final energy drift.


| Test name | Exercises |
|---|---|
| `gmd_smoke_inline_lj` | inline LJ force field, NVT |
| `gmd_smoke_ewald` | Ewald electrostatics, NVT |
| `gmd_smoke_pme` | PME electrostatics, NVT |
| `gmd_smoke_mc_barostat` | MC barostat NPT, Nosé-Hoover |
| `gmd_smoke_molecular` | molecular FF (`.ff` + `.top`), bonded forces and special-pair LJ path |
| `gmd_special_pair_rules` | four-atom chain: hand-checked LJ exclusion/scaling plus Ewald/PME Coulomb corrections |
| `gmd_constraint_solver` | SHAKE/RATTLE two-atom, water geometry, and constrained NVE stability checks |
| `gmd_checkpoint_restart` | serial 100-step continuity vs 50+checkpoint+restart+50, state round-trip, bad-file errors |
| `gmd_pme_modes` | fixed charged system: replicated PME vs distributed-PME interface/prototype selection path energy/force equivalence, both using replicated numerics |
| `gmd_validation_static_lj_cluster` | validation suite: fixed-coordinate LJ cluster vs stored baseline |
| `gmd_validation_static_special_pairs` | validation suite: analytic 1-2/1-3 exclusion and 1-4 scaling reference, including a modified-scale variant |
| `gmd_validation_static_coulomb` | validation suite: charged Ewald analytic single-point reference plus provisional PME regression baseline |
| `gmd_mpi_periodic_1d_ghost_exchange_2proc` | 1D periodic wraparound ghost exchange (rank 0 ↔ rank 1) |
| `gmd_mpi_periodic_1d_force_consistency_2proc` | distributed vs serial LJ force/energy across a periodic x boundary |
| `gmd_mpi_periodic_1d_migration_2proc` | atom migration across periodic x high→low and low→high boundaries |
| `gmd_mpi_special_pair_consistency_2proc` | cross-rank four-atom chain: tag-based LJ/Ewald special-pair behavior |
| `gmd_mpi_constraint_solver_2proc` | cross-rank SHAKE/RATTLE constraint projection by global atom tag |
| `gmd_mpi_checkpoint_restart_2proc` | MPI ranks read the same global checkpoint and recover atom state consistently |
| `gmd_mpi_domain_decomposition_3d_grid_4proc` | 3D process grid creation and coordinate mapping (4 ranks) |
| `gmd_mpi_domain_decomposition_3d_migration_8proc` | 3D multi-axis atom migration with periodic wrapping (8 ranks) |
| `gmd_mpi_ghost_exchange_3d_face_8proc` | 3D face ghost exchange (8 ranks) |
| `gmd_mpi_ghost_exchange_3d_edge_8proc` | 3D edge ghost exchange (8 ranks) |
| `gmd_mpi_ghost_exchange_3d_corner_8proc` | 3D corner ghost exchange (8 ranks) |
| `gmd_mpi_ghost_exchange_3d_periodic_8proc` | 3D periodic corner wraparound ghost exchange (8 ranks) |
| `gmd_mpi_reverse_force_3d_face_8proc` | 3D face reverse force accumulation (8 ranks) |
| `gmd_mpi_reverse_force_3d_edge_8proc` | 3D edge reverse force accumulation (8 ranks) |
| `gmd_mpi_reverse_force_3d_corner_8proc` | 3D corner reverse force accumulation (8 ranks) |
| `gmd_mpi_lj_3d_periodic_consistency_8proc` | 8-rank periodic LJ energy/force vs serial consistency |
| `gmd_smoke_mpi_lj_2proc` | 2-process MPI LJ NVT run (verifies execution does not crash) |
| `gmd_smoke_mpi_ewald_2proc` | 2-process MPI Ewald smoke run on a true cross-rank charged system |
| `gmd_smoke_mpi_pme_2proc` | 2-process MPI PME smoke run on a true cross-rank charged system |
| `gmd_smoke_mpi_molecular_2proc` | 2-process MPI bonded molecular smoke run across a rank boundary |
| `gmd_smoke_mpi_lj_consistency` | serial vs 2-process MPI energy consistency check (verifies numerical equivalence) |
| `gmd_smoke_mpi_ewald_consistency_8proc` | serial vs 8-process Ewald consistency on a 3D process grid |
| `gmd_smoke_mpi_pme_consistency_8proc` | serial vs 8-process replicated-PME consistency on a 3D process grid, including empty local domains |
| `gmd_mpi_ml_edge_shift_1proc` | ML provider runs normally at one rank (proves the MPI refusal below is about rank count) |
| `gmd_mpi_ml_edge_shift_2proc` | ML provider refuses MPI on every rank; adapter never reached |
| `gmd_mpi_ml_edge_shift_4proc` | same at 4 ranks, on a dimer straddling a periodic face |
| `gmd_torchscript_edge_shift` | **`GMD_ENABLE_TORCH=ON` only** — `edge_shift` observed inside a real scripted model |

---

## Validation

GMD now ships a lightweight validation suite under `validation/`:

- `static_*` cases are short, deterministic, and registered in CTest.
- `nve_*`, `nvt_*`, `npt_*`, and `diffusion_*` are optional longer workflows.
- Static single-point checks use `gmd_validate`, which writes JSON containing component energies and per-atom forces.
- Dynamics cases post-process `output.log` and `output.xyz` into `summary.json`, `thermo_series.csv`, and, where relevant, `msd_series.csv`.

Examples:

```bash
# short validation cases
python3 validation/run_validation.py \
  --case short \
  --gmd build/gmd \
  --gmd-validate build/gmd_validate \
  --work-root build/validation

# optional long workflow validation
cmake --build build --target gmd_validation_long
```

Each case directory contains:

- GMD input files
- a reference result JSON
- where available, a LAMMPS/OpenMM reproduction note or input stub
- a Python analysis script
- an explicit tolerance file
- a case README documenting intent, settings, and limitations

Reference provenance is stated per case. Entries marked `analytic_reference_*`
are independent analytic checks. Entries marked `provisional_gmd_baseline` are
repeatable regression baselines only; they still need an external software
result to become scientific cross-code validation.

Current validation status:

- **LJ static cluster**: analytic reference, suitable for short CI.
- **Special-pair / 1-4 scaling static chain**: analytic shifted-LJ + Ewald special-pair reference, including a modified-scale variant.
- **Ewald static Coulomb**: analytic periodic Ewald reference.
- **Replicated PME static Coulomb**: `validation/static_coulomb` remains a GMD regression baseline. The independent external validation is `validation/pme_external` — OpenMM 8.6.0 PME as the primary reference and LAMMPS 22 Jul 2025 PPPM plus exact Ewald as a second engine, covering total energy, every per-atom force component and all nine virial components, with a four-grid convergence study toward the exact Ewald limit.
- **Distributed PME mode**: interface/prototype selection path only; it uses replicated PME numerics and is not evidence of scalable distributed FFT performance.
- **SHAKE/RATTLE**: tested for serial projection, water geometry, constrained NVE stability, and a cross-rank MPI projection; the MPI algorithm is correctness-first global gather.
- **Checkpoint/restart**: tested through the real `gmd` CLI restart-continuity workflow in serial and MPI.

---

## Continuous Integration

GitHub Actions runs a minimal CI workflow in `.github/workflows/ci.yml`:

- serial `RelWithDebInfo` configure/build with Ninja
- full serial CTest suite, including short static validation cases
- MPI `RelWithDebInfo` configure/build with OpenMPI
- MPI test registration check via `ctest -N`
- selected 2-rank MPI smoke/unit tests, including ghost exchange, migration,
  special pairs, constraints, checkpoint/restart, restart continuity, LJ,
  Ewald, PME, and molecular smoke tests

Long validation workflows (`nve_*`, `nvt_*`, `npt_*`, `diffusion_*`) are not
part of every CI run. They remain optional/manual validation because they are
slower and many currently use `provisional_gmd_baseline` regression references:

```bash
cmake --build build --target gmd_validation_long
```

A separate manual/nightly workflow can run the full MPI CTest suite. It is kept
out of ordinary pull-request CI because the 4/8-rank tests are slower and more
sensitive to runner MPI availability.

Local CI reproduction:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --parallel
ctest --test-dir build --output-on-failure --verbose

cmake -S . -B build-mpi -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DGMD_ENABLE_MPI=ON
cmake --build build-mpi --parallel
ctest --test-dir build-mpi -N
ctest --test-dir build-mpi --output-on-failure --verbose \
  -R "gmd_mpi_periodic_1d_ghost_exchange_2proc|gmd_mpi_periodic_1d_force_consistency_2proc|gmd_mpi_periodic_1d_migration_2proc|gmd_mpi_special_pair_consistency_2proc|gmd_mpi_constraint_solver_2proc|gmd_mpi_checkpoint_restart_2proc|gmd_mpi_restart_continuity_cli_2proc|gmd_smoke_mpi_lj_2proc|gmd_smoke_mpi_ewald_2proc|gmd_smoke_mpi_pme_2proc|gmd_smoke_mpi_molecular_2proc"
```

Release evidence can be collected without touching an existing `build/`
directory:

```bash
scripts/collect_release_logs.sh /tmp/gmd-v2.4-release-logs
```

---

## Project Layout

```
GMD/
├── app/gmd_main.cpp               CLI entry point (LJ / Ewald / PME / ML / MPI workflow)
├── cmake/                         Build helper scripts (CompilerOptions, CUDAOptions, MPIOptions, …)
├── examples/
│   └── ethane_demo/               Complete molecular FF example
│       ├── ethane_demo.cpp        Standalone demo driver
│       ├── ethane.ff              Molecular force field (OPLS-AA subset)
│       ├── ethane.top             Bond/angle/dihedral topology
│       ├── ethane.xyz             Initial coordinates (C2H6, 8 atoms)
│       └── ethane.run             MD parameters (500 steps, 300 K, NHC)
├── include/gmd/
│   ├── boundary/                  PBC & minimum-image
│   ├── core/                      Simulation orchestrator
│   ├── cuda/                      CUDA backend (reserved)
│   ├── force/
│   │   ├── force_provider.hpp         abstract interface
│   │   ├── composite_force_provider.hpp
│   │   ├── classical_force_provider.hpp   (LJ)
│   │   ├── bonded_params.hpp             (BondParams, AngleParams, …)
│   │   ├── bonded_force_provider.hpp     (bonds / angles / dihedrals)
│   │   ├── ewald_force_provider.hpp
│   │   └── pme_force_provider.hpp
│   ├── integrator/                VV + thermostat + barostat
│   ├── io/                        ConfigLoader (xyz/run/ff/top), TrajectoryWriter
│   ├── ml/                        ML force provider + TorchScript adapter
│   ├── neighbor/                  VerletNeighborBuilder
│   ├── parallel/                  MPI domain decomposition + communication
│   │   ├── domain_decomposition.hpp
│   │   ├── mpi_communicator.hpp
│   │   ├── mpi_environment.hpp
│   │   └── pme_parallel.hpp
│   ├── runtime/                   RuntimeContext
│   ├── system/                    System, Box, Topology
│   └── utils/                     Utilities (reserved)
├── src/                           Implementation (mirrors include/gmd/)
├── tests/                         Smoke / integration / MPI unit / consistency tests
│   ├── mpi_periodic_1d.cpp            1D periodic MPI unit tests
│   ├── mpi_ghost_exchange_3d.cpp      3D ghost exchange + reverse force unit tests
│   ├── mpi_domain_decomposition_3d.cpp 3D grid + migration unit tests
│   ├── compare_energy_logs.cpp        Serial-vs-MPI energy log comparator
│   └── smoke_*                        Smoke test inputs and run files
├── run_ethane_demo.sh             One-command build + run for ethane demo
└── CMakeLists.txt
```

---

## Architecture Overview

```
ConfigLoader ─────────────────────────────────────────────────┐
  load_xyz()          → System (coords, masses, types, Z)      │
  load_run()          → RunConfig                              │
  load_molecular_ff() → MolecularForceFieldConfig              ├──► Simulation
  load_topology()     → Topology                               │       │
                                                               │  DomainDecomposition (MPI)
ForceProvider (interface)                                      │  MpiCommunicator (MPI)
  ClassicalForceProvider   (LJ non-bonded)                     │  NeighborBuilder::rebuild()
  BondedForceProvider      (bonds/angles/dihedrals/impropers)  │  ForceProvider::compute()
  EwaldForceProvider       (k-space Coulomb)                   │  Integrator::step()
  PMEForceProvider         (mesh Ewald)                        │       │
  MLForceProvider          (TorchScript)                       │  TrajectoryWriter::write_frame()
  CompositeForceProvider   (sums all above) ───────────────────┘
```

---

## Known Limitations

- SHAKE / RATTLE constraints use a correctness-first MPI global-gather projection; scalable local constraint communication and SETTLE/LINCS are not implemented
- Constraint forces contribute to the virial as an endpoint quantity at `t+dt` (see *Virial and pressure*), so constrained pressure is complete and single-time-level. It has **no independent external reference**, only an analytical one
- Constraint **independence is verified at the initial geometry only**. The rank of the mass-weighted constraint Jacobian is configuration dependent, and a set that is independent at the start can become degenerate later in a trajectory (three constrained atoms drifting collinear, a ring flattening). That is not re-checked during the run, because the analysis is collective and per-step cost would be significant; a set that is heading that way is flagged as ill-conditioned at start-up, which is the warning signal
- The **initial frame of a constrained run** has no completed step behind it and no RATTLE multiplier for its force evaluation, so it has no complete pressure. It is reported as `nan` with `P_valid 0`, never as a number
- The SHAKE reference-gradient linearisation has no solution when a bond turns through ~90° in a single step (`r(t+dt)·r(t) → 0`). That is diagnosed as a hard error naming the pair and asking for a smaller time step, rather than being allowed to produce a wild correction
- Off-diagonal virial components are validated for the **constraint** term only (see *Constraint independence*, below, and *Virial and pressure*). The force providers' own off-diagonal virial is still unvalidated: `Box` stores three edge lengths, so the engine is orthorhombic-only and no shear strain can be applied to finite-difference it
- A Nose-Hoover checkpoint can only be restarted into a run with the same degrees of freedom. Changing the constraint set, the centre-of-mass removal setting or the atom count is rejected, as is a checkpoint predating constraint-aware DOF accounting (the old `3N-3` rule) whenever the two disagree. There is no migration path: the thermostat mass `Q` and friction variable `xi` belong to the DOF they were generated under. Restart from the input instead
- **The initial velocity field is reproducible across rank counts to reduction round-off, not bitwise.** The draws are bitwise identical by tag, but the centre-of-mass and kinetic-energy sums are accumulated in storage order within a rank and combined by `MPI_Allreduce` across ranks, and neither order is the same at every rank count. That reaches every atom through one shared shift and one shared factor: 2.8e-17 on the initial field, amplifying to 8.9e-16 (~28 ulp) over 50 steps. Making it bitwise would need a fixed-order global reduction — an O(N) gather on every initialization — for a difference far below anything physical
- **The keyed generator's normal transform is not bit-portable across standard libraries.** Key derivation, the SplitMix64 mixer and the uniform mapping are exact integer arithmetic and identical everywhere; `sqrt`, `log` and `cos` are not guaranteed bit-identical across libm implementations, so cross-platform agreement is to a few ulp
- **The velocity initializer's degrees-of-freedom convention is its own.** It rescales against 3N−3 (or 3N) directly, while thermostats and temperature reporting use the constraint-aware `compute_degrees_of_freedom()`. For an unconstrained system the two agree exactly; for a constrained one a system initialized to T will not report exactly T. Pre-existing and not addressed here
- **`validation/berendsen_npt_lj` still reports rather than asserts its serial/MPI difference.** With the initializer corrected the two now agree to reduction round-off, but the case's comparison was written when they did not and has not been tightened; the rank-independence asserted about the barostat itself lives in `tests/mpi_berendsen_barostat.cpp`
- The internal time unit is now audited, but the **Monte Carlo barostat's `mc_frequency` is a step count, not a time**, so it does not scale with the timestep: halving `time_step` halves the physical interval between volume-move attempts. That is the documented behaviour rather than a defect, but it is the one scheduling parameter in the code that is not expressed in femtoseconds
- No GPU execution (CUDA option present but CPU-only)
- Checkpoint/restart uses a readable replicated text file; large-scale binary/parallel checkpoint I/O is not implemented
- Berendsen barostat requires virial from every active force term; barostat pressure is computed with MPI-allreduced kinetic energy and virial. Every provider now reports a physically derived virial rather than a coordinate approximation (see *Virial and pressure* below), and each contribution is reduced exactly once across the communicator
- `pme_mode distributed` currently uses the replicated PME numerical backend. It is an interface/prototype for workflow compatibility, dependency checks, and timing hooks; it does not reduce PME grid memory per rank, does not perform distributed FFT communication, and should not be used as evidence of scalable distributed-PME performance
- Ewald/PME special-pair Coulomb scaling is applied with analytical `(scale - 1) q_i q_j/r` corrections; PME retains its normal mesh discretization error and still needs broader accuracy regression coverage
- ML force provider requires `GMD_ENABLE_TORCH=ON` and a compatible TorchScript model; MPI domain decomposition is rejected because local-plus-ghost model energy ownership and message-passing halo depth are not defined. The rejection is proven at np=1/2/4 by `tests/mpi_ml_edge_shift.cpp`, and the periodic `edge_shift` tensor the model receives is directly covered by `tests/torchscript_edge_shift_tests.cpp` (see *TorchScript `edge_shift`*) — but only in a build configured with LibTorch; without it neither the adapter nor its contract is exercised at all
- MPI NVE/NVT is supported for implemented classical short-range, Ewald, replicated-PME, and current bonded paths; Berendsen NPT is also supported, but the Monte Carlo barostat is serial-only because MPI trial-volume coordination is not implemented
- The MPI LJ NVE regression tracks serial logged energies and drift within test tolerances (`gmd_smoke_mpi_lj_consistency`)
- Limited automated test coverage

---

## Contributing

Contributions are welcome. If you change behavior:

- Keep public headers in `include/gmd/` and implementations in `src/`
- Add or expand tests in `tests/` where possible
- Update this README when input formats or runtime behavior change
- Keep `ForceProvider`, `Integrator`, and `NeighborBuilder` as the primary extension points

---

## License

This project is licensed under the terms in [LICENSE](LICENSE).

## Citation

If you use GMD in research, you can cite it as software:

```bibtex
@software{gmd2026,
  title={GMD: Good Molecular Dynamics},
  author={Contributors},
  year={2026},
  url={https://github.com/kiana-0423/GMD}
}
```
