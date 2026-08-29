# Changelog

All notable user-facing changes in GMD are documented here.

## [Unreleased]

### ⚠️ Simulation-results-changing corrections

- **`GMD_ENABLE_TORCH=ON` did not compile.** `include/gmd/force/torchscript_adapter.hpp`
  carried `namespace torch::jit { struct script::Module; }`. That is ill-formed:
  an elaborated-type-specifier cannot carry a nested-name-specifier, and
  `torch::jit::script` is in any case a namespace *alias* for `torch::jit`, so
  the name being declared did not exist to declare. Every Torch-enabled build
  failed in the first translation unit that included the header, which is
  `src/force/torchscript_adapter.cpp` itself -- so the entire TorchScript path
  had never been compiled, let alone run.

  The declaration was also unnecessary. The header names no Torch type: the
  module handle lives behind an opaque `Impl`, which is what actually keeps
  `torch/script.h` out of including translation units. It is removed rather
  than corrected. No default build changes, since the option is off by default.

- **PME reciprocal-space forces were absent, and are now correct.** Three
  independent defects, found while building component-wise virial validation
  for the reciprocal mesh. **Every run using `coulomb pme` changes.** Ewald runs
  are unaffected.

  *No reciprocal force at all, at any B-spline order.* `bspline_deriv(u, p)`
  evaluates `M_(p-1)(u) - M_(p-1)(u-1)`, but `bspline()` implemented only orders
  4 and 6 and fell through to `return 0.0` for every other order. Order 4 asked
  for order 3 and order 6 asked for order 5; both got zero. Every derivative
  weight in the force interpolation was therefore zero, so PME contributed the
  real-space erfc force and nothing else. Orders 2, 3 and 5 are now implemented:
  every supported order needs its predecessor, and order 3 needs order 2, so the
  chain is closed rather than extended by one link.

  *A missing factor of the mesh point count.* The reciprocal energy is
  `E = ½ Σ_m G(m)|Q̂(m)|²` with an unnormalised forward transform, while
  `fft3d(..., true)` divides by `K1·K2·K3`. Differentiating gives
  `dE/dQ(j) = N·IFFT[G·Q̂](j)`, so the array left in the mesh is `1/N` of the
  potential the gradient needs. Without the factor every PME reciprocal force
  was `N` times too small — invisible while the weights above were zero.

  *The order-6 B-spline was wrong.* The hard-coded quintic polynomials on
  `[2,3)`, `[3,4)` and `[4,5)` were not `M_6`: the spline took negative values,
  `Σ_k M_6(t+k)` came to 0.9 instead of 1, and `M_6(u) = M_6(6-u)` failed.
  `pme_order 6` returned 13401 eV where the correct energy is -3.05 eV. The
  order-4 polynomials were correct and are unchanged.

  `PMEForceProvider::bspline()` and `bspline_deriv()` are now public static
  members so they can be tested directly rather than through the provider's
  behaviour. They are pure functions of `(u, order)`; leaving them private is
  what let a total loss of the reciprocal force sit behind a public surface that
  looked healthy.

  The corrected PME now converges to a high-accuracy Ewald reference in energy,
  force and all nine virial components as the mesh is refined, with order 6
  converging faster than order 4. `validation/static_coulomb/reference_pme.json`
  has been regenerated: the superseded baseline had recorded a reciprocal force
  of zero, so its stored forces were ~94% too small. The regenerated forces
  agree with that case's analytic Ewald reference to 8.2e-4 — the expected mesh
  error at grid 16³, order 4, alpha 0.3 — where the old ones disagreed by
  8.9e-2. The stored energies are unchanged, which is why an energy-only
  regression never caught any of this.

- **Constrained dynamics now uses the standard velocity-Verlet SHAKE/RATTLE
  splitting.** Two defects are corrected, and **every constrained trajectory
  changes** as a result.

  *SHAKE projected along the wrong gradient.* `sigma_c = |r_c|^2 - d_c^2` has
  gradient `2 r_c` at the configuration the step starts from, so the projection
  must displace atoms along `r_c(t)`. It displaced them along the drifted bond
  `r_c(t+dt)` instead. Both land on the constraint manifold, but at different
  points, and only the first reproduces the constrained equations of motion. The
  second is a plain rescaling of the bond: not symplectic, and it bleeds energy
  away steadily. A free rigid rotor lost **86% of its kinetic energy over 4000
  steps** at `omega*dt = 0.04`, where the corrected form conserves it to 1e-13.

  *The SHAKE half-step velocity impulse was missing.* The position correction
  carries a matching velocity impulse `v += dr/dt`; without it the velocities
  were inconsistent with the constrained positions and RATTLE had to absorb the
  whole step's constraint impulse rather than the second half-kick's share.

  A bond that turns through ~90 degrees in one step now raises a hard error
  naming the pair, instead of producing a wild correction: the reference-gradient
  linearisation has no solution there.

- **Holonomic constraint forces now contribute to the virial, and hence to the
  pressure.** Constrained runs previously reported a pressure that was simply
  missing the constraint term, so **every reported pressure for a system with
  SHAKE/RATTLE constraints changes**. Unconstrained runs are unaffected.

  The contribution is the **endpoint** constraint force at `t+dt`, recovered
  from the converged RATTLE multipliers. In the splitting above the constraint
  force enters twice per step: through the SHAKE impulse, paired with `F(t)`,
  and through the RATTLE impulse, paired with `F(t+dt)`. The providers are
  evaluated at `r(t+dt)`, so the virial takes the second. Matching it against the
  second half-kick term by term gives

  ```
  (dt/2m_i) G_i(t+dt) = w_i Lambda_c r_c(t+dt)
  W_constraint(t+dt)  = sum_c (2 Lambda_c / dt) (r_c (x) r_c)
  ```

  in the project's `W = sum r (x) F` convention. The factor `2/dt` is exact, and
  **both halves of the reported virial live at the same time level** — there is
  no step-averaged quantity paired with an endpoint one anywhere. The factor was
  established against an independent physical result, not assumed: a freely
  rotating rigid dimer's endpoint constraint force must equal the centripetal
  force `mu omega^2 d`, and `1/dt` undershoots it by exactly two.

  Because RATTLE does not move atoms, the bond vector is fixed for the whole
  solve and summing the *scalar* multipliers gives a pair force that is exactly
  central and exactly antisymmetric — for coupled constraints as much as for
  isolated ones. Nothing is projected and no residual is discarded, and the
  tensor is symmetric by construction. It is built from minimum-image bond
  vectors, so it is independent of the origin and of how atoms are wrapped.

  **RATTLE is not redundant with the kinetic term.** `2K` and the configurational
  constraint virial are different quantities; projecting the velocities does not
  remove the need for the second. For a rigid rotor `2K = mu omega^2 d^2` and
  `tr W_constraint = -mu omega^2 d^2`, so `2K + tr W = 0` — a rigid body's frozen
  internal coordinate contributes nothing to the pressure. Omitting the
  constraint virial leaves `2K` uncancelled and reports a spurious `2K/3V`.

  Validated against the centripetal force of a rigid rotor, a reference
  independent of the implementation, for a dimer and for a rigid triangle of
  three coupled constraints: endpoint force magnitude, sign, trace, the identity
  `2K + tr W = 0`, and second-order convergence of both errors in `dt`. Plus
  exact symmetry and centrality, the two constraint half-impulses being measured
  separately and shown to differ, free-rotor energy conservation, additivity,
  translation and periodic-wrapping invariance, rank invariance at 1, 2 and 4
  ranks, and restart continuity of the reported pressure. No independent
  external reference is claimed; see *Remaining limitations*.

- **The reported pressure is now the completed step's, not the post-rescale
  geometry's.** The pressure of the step that just finished and the virial
  attached to the geometry a barostat has since rescaled are different things,
  and they are now stored separately. The completed-step pressure is captured at
  the end of the step, before any barostat runs, and is bit-for-bit the value the
  barostat itself consumed; a later force evaluation cannot destroy it.

  This is what makes constrained NPT report a pressure at all. The post-rescale
  re-evaluation has no constraint partner, so under the previous arrangement a
  constrained run under a cell-rescaling barostat had nothing complete to report.

  For unconstrained runs the two definitions differ only on a step that actually
  rescales; all eight validation cases are unchanged, and in the NPT case this
  was checked rather than assumed — with `mc_frequency 20` against
  `output_interval 10`, no logged frame is itself a rescale step.

Every item in this section changes numerical results for the configurations it
affects. **Trajectories, energies, temperatures and pressures produced by
earlier versions are not reproducible after upgrading** for systems using
constraints, proper dihedrals, impropers, charged NPT, or MPI barostat runs.
Re-run anything you intend to compare against; do not mix old and new output in
one analysis.

- **Proper dihedral and improper forces were wrong.** *This is the most serious
  item here.* **Earlier versions could produce incorrect trajectories for any
  system using proper dihedral or improper terms.** The terminal-atom force had
  a flipped overall sign, and the middle-atom projection used the coefficients
  belonging to the opposite `b1 = r_i - r_j` convention while `b1` is built as
  `r_j - r_i`. Energies were unaffected — only forces — so the error was
  invisible to energy-only checks and drove torsions toward the wrong states.
  It is also invisible in the virial trace, because a torsion angle is
  unchanged by isotropic scaling. Now validated externally against LAMMPS
  (`validation/static_bonded_reference`) and internally against the energy
  gradient (`tests/bonded_force_gradient_tests.cpp`).
- **Temperature and thermostat degrees of freedom now account for
  constraints.** All temperature consumers previously hardcoded `3N - 3`,
  ignoring both SHAKE/RATTLE constraints and the centre-of-mass removal setting.
  A constrained run therefore reported a temperature that was too low, and the
  thermostat compensated by driving the system hotter than its target. The DOF
  count is now `3N`, less 3 when the COM velocity is removed, less one per
  distinct constraint, computed in one shared place and used by both thermostats
  and by trajectory output.
- **Periodic bonded and reciprocal-space virial accounting corrected.** The
  bonded virial was formed from absolute wrapped coordinates and silently
  changed value when a molecule crossed a periodic boundary; it is now built
  per interaction from minimum-image separations. The Ewald/PME reciprocal
  virial was computed as `Σ r_i ⊗ F_i`, which is not the virial for an energy
  that depends on the cell explicitly; it now uses the analytic k-space form.
  Ewald was additionally double-counting its reciprocal contribution. **Pressure
  and NPT behaviour change for any charged system**, and for molecular systems
  whose molecules straddle a boundary.
- **Forces are refreshed after a barostat volume change.**
  `VelocityVerletIntegrator::step()` previously returned holding forces computed
  for the pre-scaling geometry, which the next step's half-kick then integrated.
  NPT trajectories change. `Simulation::step()` conversely no longer re-evaluates
  forces when the cell did not change (a rejected Monte Carlo trial, or a step
  with no volume move), removing a redundant force evaluation.
- **Incompatible Nosé-Hoover checkpoints are now rejected.** See *Checkpoint
  compatibility* below. Runs that previously restarted and silently continued
  with a mismatched thermostat mass now stop with an error.

### Dependent constraint sets are now rejected

**Results-changing for any run whose constraints were not independent** — such a
run now stops at start-up instead of reporting a temperature that was too high.

Distinct was not independent. The solver collapsed duplicates and rejected
conflicts, but a redundant closed topology (every pair among five atoms: ten
constraints over a body with nine internal degrees of freedom) or any three
collinear atoms was accepted, and every distinct constraint was subtracted from
the degrees of freedom. The redundant rows also make `J M⁻¹ Jᵀ` singular, so the
multipliers SHAKE and RATTLE converged to were not unique.

`ConstraintSolver::analyze_independence()` computes the rank of the mass-weighted
constraint Jacobian `J_M = J M^(−1/2)`, per connected component of the constraint
graph, using a one-sided Jacobi SVD with the standard `max(rows, cols)·ε·σ_max`
relative rank tolerance. The decomposition's convergence is verified against a
scale-free off-orthogonality residual and throws rather than returning a rank
from an unconverged factorisation.

`VelocityVerletIntegrator::initialize()` analyses the supplied geometry for early
diagnostics, projects positions onto the constraint manifold, and then calls
`require_independent()` on the **converged projected geometry**, which is the
authoritative report; velocities are projected only after it passes, and degrees
of freedom are computed from that accepted state. Rank is a property of the
configuration, and a set that is full rank as supplied can project onto a
degenerate one. Targets that admit no non-degenerate configuration at all — a
constrained triple whose distances violate or exactly meet the triangle
inequality — are rejected at construction, before any geometry exists.

The diagnostic names the component, its atom tags, the singular values, whether
the dependence is structural or geometric, and which rows are redundant when a
column-pivoted rank-revealing factorisation agrees with the singular values about
how many there are. After `require_independent()` returns the constraint count
**is** the rank, so the degrees-of-freedom subtraction is exact rather than
assumed.

A set that is independent but ill-conditioned (`σ_max/σ_min > 1/√ε`) is accepted
with a warning rather than rejected. Reversed duplicate pairs — the same bond
listed as `(i, j)` and `(j, i)` — are now counted and reported separately, since
a topology that does that is usually a generation bug.

### MPI trajectory output reported stale frame state

**Bug fix, MPI runs only.** The MPI output path writes a separate `System`
holding the gathered global coordinates. That object was copied from the real
system at start-up and thereafter received only coordinates and the potential
energy, so every other frame field it reported was frozen at its start-up value.
In practice an MPI log reported the **initial box** — so the volume column of an
MPI NPT run never moved — and **zeroed SHAKE/RATTLE iteration counts and
residuals**. Both predate this release.

All non-atomic frame state is now carried across explicitly, and serial, `np=2`
and `np=4` runs of the real executable are compared column by column over whole
logs.

### Trajectory-log format

The energy log gains a trailing **`P_valid`** column, and an unavailable pressure
is now written as **`nan`** rather than `0`.

**`PE`, `KE`, `E_total`, `T`, `P` and `V` in a row now describe one state**: the
completed step, captured together before any barostat rescale. They are mutually
consistent, and `E_total` is the energy of the configuration that produced the
reported pressure. Previously `P` came from one state and `PE`/`V` from the
post-rescale geometry. For a run without a barostat nothing changes — there is
no rescale to be on either side of. Under a barostat the row's `V` now lags the
`.xyz` coordinates by one rescale: the coordinates are the current ones, since
snapshotting them per step would cost an `N x 3` copy for a difference no
thermodynamic quantity in the row depends on.

A numerical zero is a perfectly ordinary pressure, so it must never double as a
sentinel for "no pressure available". The only case that produces one today is
the initial frame of a *constrained* run: no step has completed, so there is no
completed-step pressure, and the initial force evaluation has no RATTLE
multiplier to pair with. Unconstrained runs report a valid pressure on every
frame, initial one included.

The column is appended, so parsers that read the log by position are unaffected.

### Checkpoint compatibility

The checkpoint format moves to **version 2**, which appends three explicitly
named blocks. Reading stays backward-compatible — version 1 files are still
accepted, and every existing field keeps its name, order and meaning.

- `constraint_virial <state> <time_level> <9 components>` — the constraint
  contribution, its validity (`not_applicable` / `unavailable` / `valid`) and
  **which multiplier it is**. The only time level this build writes is
  `endpoint_rattle_t_plus_dt`, and a restart refuses any other rather than
  pairing an unknown time level with the endpoint provider virial it computes.
  Nothing in the checkpointed state determines this value — it is a property of
  the step that reached the state, which is why it has to be persisted.
- `provider_virial <valid> <9 components>` — the current-geometry provider
  virial. Recorded for diagnosis but **not** installed on restart: the restarted
  run recomputes it from the checkpointed coordinates, which is where it came
  from, and reinstalling a value summed under a different rank decomposition
  would only introduce a discrepancy.
- `step_pressure <valid> <P> <2K> <V> <9 components>` — the completed step's
  pressure and the pieces it was built from. Installed on restart, so the
  restarted run reports for that frame the identical number the uninterrupted
  run reported rather than reconstructing it.

A version-1 checkpoint carries none of these. A constrained run restarted from
one reports its first frame's pressure as unavailable (`nan`, `P_valid 0`); every
subsequent step is unaffected, because each closes with its own RATTLE.

Older builds cannot read a version-2 checkpoint.

What changed is validation of the Nosé-Hoover `dof` field. It is now checked
against the degrees of freedom the current run computes, and **never installed**.

| Old checkpoint | Outcome |
|---|---|
| Nosé-Hoover, unconstrained, COM velocity removed | **Accepted.** The old `3N - 3` equals the new count, so these restart unchanged and deterministically. |
| Nosé-Hoover, with constraints active | **Rejected.** The old `3N - 3` disagrees with the constraint-aware count. |
| Nosé-Hoover, COM velocity kept (`remove_com_velocity false`) | **Rejected.** The old rule subtracted 3 regardless. |
| Restart into a different constraint set, COM setting, or atom count | **Rejected.** |
| Velocity-rescaling thermostat, or no thermostat | **Accepted.** These carry no extended-system state to invalidate. |
| Any non-thermostat checkpoint content (positions, velocities, box, topology) | **Accepted**, unchanged. |

A rejection names both the checkpoint's DOF and the run's DOF, and lists the
likely causes.

**Why mismatched Nosé-Hoover state cannot be migrated.** The thermostat mass
`Q = dof · k_B · T_target · τ²` and the friction variable `ξ` were both
generated under the checkpoint's DOF, and `ξ` carries the accumulated history of
an extended system defined by that value. Changing the DOF defines a *different*
extended system. Rescaling `Q` would not reconstruct the `ξ` trajectory that the
new system would have produced, so no rescaling makes the continued run a
continuation of the original one. Silently rescaling would produce a plausible
but unfaithful trajectory, which is worse than stopping. Restart from the input
instead.

### Added

- **External bonded reference case** `validation/static_bonded_reference`,
  independently validating, for bond, angle, proper dihedral and improper terms:
  the **functional forms**, the **angle and sign conventions** (including the
  proper-dihedral phase sign, proven with a sign-sensitive δ = 60° case, and the
  derived improper sign mapping GMD `φ0 = -χ0`), **periodic wrapping
  behaviour**, and the resulting **energies and per-atom forces** — against
  LAMMPS `22 Jul 2025 - Update 5`, intact and wrapped across a periodic
  boundary, in a non-cubic box.

  The reference is converted to eV using GMD's own declared
  `kcal_to_eV = 4.336410e-2`, so that constant cancels out of the comparison.
  **The numerical value of the unit-conversion constant is therefore not
  independently verified by this case**; everything downstream of it is. Units,
  conventions, atom ordering and tolerances are documented in the case README.

  The suite reads the checked-in `reference.json`, so running validation does
  not invoke or require LAMMPS.
- **Component-wise virial validation.** Each diagonal component is checked
  independently against a per-axis finite difference on non-cubic cells, serial
  and under 4-rank decomposition.
- **Full nine-component validation of every virial source.** `tests/virial_reference.hpp`
  holds test-only references written from each term's own definition:
  analytical pair identities with the force cross-checked against `-dV/dr`;
  bonded force moments built without `accum_interaction_virial()`, from forces
  first verified against a central difference of an independently written
  energy; and — for the reciprocal-space terms, where no pair identity exists —
  an Ewald and a PME reciprocal energy written against a **general 3×3 cell
  matrix**, complete with its own B-spline recursion and direct DFT.

  That last piece is what makes the off-diagonal components reachable. `Box`
  stores three edge lengths, so the engine cannot apply a shear strain and the
  existing finite-difference test validates only the trace and the diagonal. A
  *reference* implementation is under no such restriction: differentiating it
  under a general strain, shear included, yields all nine components, and the
  engine's analytic tensor is compared against that derivative at the
  orthorhombic configuration both agree on.

  Covered: LJ (attractive and repulsive branches, mixed types, explicit pair
  overrides, 1-4 scaling at 0 / partial / 1 with exact linearity); bond, angle,
  sign-sensitive proper dihedral with non-zero phase, and harmonic improper,
  each intact and wrapped; Ewald real space, reciprocal space, self term
  (required to contribute exactly zero) and net-charge correction (required to
  be exactly isotropic); the special-pair Coulomb correction; PME real space
  and reciprocal mesh at B-spline orders 4 and 6 over several meshes, alphas and
  a non-cubic mesh, reported as component-wise convergence tables rather than a
  tolerance copied from one run; `CompositeForceProvider` addition and validity
  propagation; and `MLForceProvider`'s absent-virial contract. Every check
  carries a non-zero guard so a blank tensor cannot pass, and periodic coverage
  includes each face, the box corner, and translate-and-wrap invariance.

  `tests/mpi_virial_sources.cpp` runs the same fixtures at 1, 2 and 4 ranks
  against the same rank-count-independent references, with explicit guards
  against a tensor multiplied by the rank count and against special-pair or
  bonded terms spanning a rank boundary being counted more than once.
- **Direct regression tests for the production PME kernel.**
  `tests/pme_bspline_tests.cpp` evaluates `PMEForceProvider::bspline()` and
  `bspline_deriv()` themselves at orders 2 through 6 — every polynomial
  interval against the order recursion, the exact knots and both sides of them,
  zero outside the support, non-negativity, `M_p(u) = M_p(p-u)`, partition of
  unity, continuity through derivative `p-2`, the derivative identity, a
  coordinate finite difference, and explicit regression values on the three
  intervals that were wrong.

  `tests/pme_energy_force_tests.cpp` pins the transform convention against an
  explicit direct DFT with both directions unnormalised, on a small mesh, which
  makes a missing or doubled `K1·K2·K3` an immediate absolute failure rather
  than a convergence one; and checks `F = -dU/dr` by central differences of the
  provider's own energy across orders 4 and 6, two mesh sizes, three alphas,
  neutral and net-charged systems, an asymmetric configuration, and atoms on
  and across the periodic faces. Alpha and the cutoff are always passed
  strictly positive so `resolve_params()` cannot change anything between the
  two sides of a difference, and repeated evaluation is asserted to be bitwise
  reproducible.

  `tests/ml_virial_contract_tests.cpp` covers the `MLForceProvider` contract:
  no virial reported, a stale `ForceResult` cleared, no force moment
  synthesised, any composite containing it invalidated in either order, and a
  pressure-coupled barostat unable to consume the absent value.
- **Direct coverage of the TorchScript `edge_shift` tensor.**
  `tests/torchscript_edge_shift_tests.cpp` observes `edge_shift` from *inside* a
  real scripted model, reached through the production
  `MLForceProvider` -> `TorchScriptModelRuntimeAdapter` -> `torch::jit` path.
  The existing `tests/box_image_flag_tests.cpp` covers
  `NeighborList::image_flags`, which is the helper; it says nothing about the
  half-to-full-graph expansion, the reverse edge, the multiplication by the box
  lengths, the `src`/`dst` ordering inside `edge_index`, or whether the tensor
  that was built is the one that reaches `forward()`.

  Two observer models are built from C++ via `torch::jit::Module::define` at run
  time, so the test needs no Python and no checked-in `.pt` file, and the
  artifacts are written to the CTest working directory rather than the source
  tree. Atoms carry unique `species` values, so every assertion is keyed to a
  stable physical identity instead of a position in whatever order the
  neighbour list produced. Expected shifts are written out from the fixture
  geometry, never obtained from `image_flags` or `apply_minimum_image`.

  Serial cases: no boundary crossing (all shifts exactly zero); each of the x,
  y and z faces in both sign directions; a corner crossing with all three
  components non-zero and of mixed sign; the reverse directed edge of every
  fixture, required to be the exact negative and to correspond to the right
  `edge_index` entry; wrapped versus unwrapped placement of the same physical
  dimer, which must produce *different* shifts and the *same* displacement;
  every axis confirmed periodic, since no per-axis periodicity flag exists to
  disable one; and atom-order permutation, including swapping the two atoms of
  every pair. Composite: an ML child inside `CompositeForceProvider`, alone and
  between two other providers, evaluated repeatedly, with the ML contribution
  recovered by subtraction and required to decode to the same shifts. Every
  shift component is compared exactly; the box lengths and species are chosen
  so no rounding occurs anywhere in the chain.

  Negative controls run mutated observer models -- shifts zeroed, sign
  reversed, x/y/z permuted -- and require every detectable edge to be rejected.
  A zero shift is a fixed point of all three mutations, so the interior fixture
  is excluded from the count rather than weakening the assertion.

- **The ML/MPI limitation is proven rather than asserted.**
  `tests/mpi_ml_edge_shift.cpp` runs at 1, 2 and 4 ranks on a dimer straddling a
  periodic face. At one rank the provider must evaluate normally; at more, both
  `initialize()` and `compute()` must throw on *every* rank with a message
  naming MPI, a composite containing an ML child must propagate the refusal
  rather than swallow it, and the model adapter must never be reached -- which
  is what rules out a duplicated per-rank contribution for a cross-boundary
  edge. It uses a stub adapter, so the proof holds in MPI builds without
  LibTorch.
- **Bonded force-gradient tests** verifying every bonded force against
  `-dU/dx`.
- **End-to-end constrained Nosé-Hoover restart tests** driving the real CLI
  through configuration parsing, `Simulation::initialize()`, checkpoint loading,
  thermostat validation and continued integration — serial and MPI — plus
  negative cases for a changed COM setting and a changed constraint set.
- **Constraint-list normalisation diagnostics.** Repeated atom pairs are
  classified as exact duplicates, tolerance-equivalent duplicates, or conflicts,
  and counted separately. `ConstraintSolver` performs no logging itself; it
  exposes `normalization_diagnostics()` and `gmd` reports them at start-up.
- **Box-dimension validation.** Zero, negative and non-finite edge lengths are
  rejected where a box is installed.
- **`.gitignore`**, and 295 generated build artifacts removed from Git tracking.

### Changed

- **Neighbor-list rebuild decisions are collective under MPI** — a correctness
  safeguard, not a result-changing fix. `VerletNeighborBuilder::needs_rebuild()`
  inspects only locally owned atoms, so in principle one rank could rebuild
  while another kept a stale list. Local flags are now combined with a logical
  OR, so every rank rebuilds during the same force evaluation or none does.

  In the current domain-decomposition path this is expected to be a no-op:
  `reverse_accumulate_ghost_forces()` ends in `clear_ghost_atoms()`, which
  clears the neighbor list, so the `!valid` term already forces a rebuild on
  every rank every step and `needs_rebuild()` is never reached. **This change
  therefore may not alter existing domain-decomposed trajectories.** Its value
  is that it makes the all-or-nothing invariant explicit and enforced: it
  prevents divergent control flow on any path where a valid neighbor list
  survives between steps (as it does today when a communicator is used without
  a domain decomposition), and it protects future implementations — for example
  one that stops discarding the list during ghost cleanup — from silently
  reintroducing the divergence.
- `(i, j)` and `(j, i)` are now the same constraint; repeated pairs collapse to
  one entry with deterministic first-value-wins semantics, and a repeat whose
  target distance differs by more than the constraint tolerance is rejected.
- Constraint terminology no longer claims independence the implementation cannot
  prove: counts are of *distinct* constraints.

### Remaining limitations

- **The constraint virial has no independent external reference.** It is
  validated against the analytical centripetal-force result for a rigid rotor,
  plus invariance and convergence properties, but not against LAMMPS or OpenMM:
  neither exposes a separately extractable constraint virial whose definition
  could be *proven* equivalent to this one, and the quantity is dynamical, so it
  cannot be compared from a static configuration. Only the analytical claim is
  made.
- **The initial frame of a constrained run has no complete pressure.** No step
  has completed, and the initial force evaluation has no RATTLE multiplier to
  pair with. It is reported as `nan` with `P_valid 0`, never as a number.
- **A bond may not turn through ~90 degrees in one step.** The SHAKE
  reference-gradient linearisation has no solution there; this is diagnosed as a
  hard error naming the pair and asking for a smaller time step.
- `Box` stores three edge lengths, so the engine remains orthorhombic-only and
  cannot itself apply shear strain. All nine virial components are nevertheless
  validated against independent general-cell references, including shear
  derivatives, as described under *Added* above.
- Constraint independence is validated on the converged projected geometry.
  The defensive rank-selection fallback remains API-tested because no physical
  fixture has been found that makes the primary and fallback rank counts disagree.
- The virial is not compared against LAMMPS; its per-term decomposition and
  GMD's do not have proven-equivalent semantics for the reference fixture.
- TorchScript `edge_shift` is directly covered (see *Added*), but only in a build configured with `GMD_ENABLE_TORCH=ON`. A build without LibTorch does not register the test and reports that at configure time; in such a build the contract is unverified.

## [v2.4] - 2026-05-24

### Added

- **Validation-progress release documentation** for the current GMD research-platform status.
  - README now states which paths are analytic validation, regression baseline, prototype interface, or tested workflow.
  - `validation/README.md` records release-facing validation maturity for LJ, Ewald, replicated PME, distributed PME mode, SHAKE/RATTLE, and checkpoint/restart.
- **Release evidence tooling**.
  - Added `scripts/collect_release_logs.sh` to collect serial/MPI configure, build, CTest, validation summaries, and environment information without deleting or reusing the repository `build/` directories.
  - Added `docs/release_logs/README.md` describing expected release artifacts.
- **PME external validation plan** under `validation/pme_external/`.
  - Documents planned LAMMPS PPPM or OpenMM PME reference cases and required comparisons.
  - Explicitly keeps replicated PME marked as tested/regression-only until an external reference is committed.
- **Full MPI GitHub Actions workflow**.
  - Added a manual/nightly workflow that runs the full MPI CTest suite and uploads CTest logs.

### Changed

- Bumped the project version to `2.4.0`.
- Corrected release-facing documentation that still listed SHAKE/RATTLE and checkpoint/restart as not implemented.
- Clarified that `pme_mode distributed` is an interface/prototype using the replicated PME numerical backend, not true distributed mesh ownership or distributed FFT.
- Updated project inventory and architecture documents to reflect v2.4 status and current file/test counts.

### Validation Status

- LJ static cluster: analytic validation.
- Special-pair / 1-4 scaling static chain: analytic validation.
- Static Ewald Coulomb: analytic validation.
- Replicated PME: tested, but independent external validation is still pending.
- Restart continuity: tested through the real CLI workflow in serial and MPI.

## [v2.3] - 2026-05-23

### Changed

- **Directory restructuring** — merged 12 subdirectories into 6 for better cohesion:
  - `boundary/` → `system/` (PBC and minimum-image are system geometry properties)
  - `neighbor/` → `system/` (neighbor lists are part of the system data structure)
  - `runtime/` → `core/` (runtime context is a simulation-core utility)
  - `ml/` → `force/` (ML force providers are force-field implementations)
  - Removed empty `cuda/` and `utils/` placeholder directories.
- All `#include` paths updated accordingly (e.g. `gmd/boundary/` → `gmd/system/`).
- `CMakeLists.txt` and `run_ethane_demo.sh` source file lists updated.
- No code logic or architecture changed — only file locations and include paths.

## [v2.2] - 2026-05-23

### Added

- **3D MPI domain decomposition** — the 1D x-axis decomposition introduced in v2.1 is now a full 3D process grid.
  - `DomainDecomposition` supports arbitrary `{Px, Py, Pz}` process grids with face, edge, and corner neighbor communication.
  - `choose_processor_grid()` uses `MPI_Dims_create` when MPI is available, falling back to a balanced factorisation for non-MPI builds.
  - `create_1d_decomposition()` remains available for backward compatibility and simple x-only splitting.
  - `neighbor_rank(offset)` correctly wraps periodic neighbors and returns `no_rank` (-1) for non-periodic boundaries.
  - `owner_rank(pos)` wraps periodic coordinates before determining the owning process.
- **`--proc-grid Px Py Pz` CLI flag** — explicitly specify a custom 3D process grid at runtime.
- **`mpi_grid Px Py Pz` run.in directive** — configure the process grid from the run file.
- **Comprehensive MPI unit tests** (17+ CTest targets when `GMD_ENABLE_MPI=ON`):
  - `gmd_mpi_periodic_1d_*` — 1D periodic ghost exchange, force consistency, and atom migration (2 ranks).
  - `gmd_mpi_domain_decomposition_3d_*` — 3D grid creation, coordinate mapping, and multi-axis atom migration (4/8 ranks).
  - `gmd_mpi_ghost_exchange_3d_*` — 3D face/edge/corner ghost exchange, periodic corner wrapping, and reverse force accumulation (8 ranks).
  - `gmd_mpi_lj_3d_periodic_consistency` — 8-rank periodic LJ energy/force consistency vs serial (8 ranks).
  - `gmd_smoke_mpi_ewald_2proc` / `gmd_smoke_mpi_pme_2proc` — 2-rank Ewald/PME smoke tests on cross-rank charged systems.
  - `gmd_smoke_mpi_molecular_2proc` — 2-rank bonded molecular smoke test across a rank boundary.
  - `gmd_smoke_mpi_ewald_consistency_8proc` / `gmd_smoke_mpi_pme_consistency_8proc` — serial vs 8-rank Ewald/PME consistency on 3D process grids.
- **Periodic 1D MPI boundary migration test data** — `tests/smoke_mpi_periodic_lj.{xyz,run}` and `tests/smoke_mpi_boundary_migration.{xyz,run}` with atoms placed to exercise wraparound ghost exchange and redistribution.

### Changed

- **`MpiCommunicator` ghost exchange now uses 3D neighbor offsets** — `neighbor_offsets()` iterates up to 26 face/edge/corner neighbors, automatically skipping offsets along axes where `proc_grid[dim] == 1`. `periodic_shift()` wraps coordinates for periodic boundaries in all three dimensions.
- **`MpiCommunicator` reverse force accumulation uses `MPI_Alltoallv`** — forces on ghost atoms are packed per home rank and exchanged in a single all-to-all collective, replacing the previous pair-wise `MPI_Sendrecv` scheme.
- **`MpiCommunicator` atom redistribution uses `MPI_Allgatherv`** — all atom states are gathered globally, periodic coordinates are normalized, and `owner_rank()` determines the new owning process.
- **`VerletNeighborBuilder` ghost image flags extended to 3D** — `S[dim]` is now computed for all three dimensions (previously only x), correctly handling ghost atoms that are periodically shifted along y or z in 3D decompositions.
- **`DomainDecomposition` periodic state stored as `std::array<bool, 3>`** — replaces the previous single `periodic_x` flag, enabling per-dimension periodicity control.
- **`validate_rank_grid()` added** — verifies that the MPI rank matches the domain process grid coordinate at the start of every collective operation.
- **`MpiCommunicator` exposes `pack_send_buffer()` and `unpack_recv_buffer()`** as private helpers with clear 3D offset semantics.

### Fixed

- Ghost image shift vectors (`S[dim]`) in `VerletNeighborBuilder` now correctly computed for all three spatial dimensions, fixing incorrect neighbor images in 2D and 3D process grids.
- `MPI_Allreduce` collectives in force providers (`ClassicalForceProvider`, `EwaldForceProvider`, `PMEForceProvider`, `BondedForceProvider`) are now placed outside early-return paths so ranks with zero atoms still participate.
- `compute_twice_ke()` and thermostat initializers use globally allreduced atom counts for consistent `dof_` across all ranks.
- `VelocityInitializer` COM velocity removal and temperature rescaling use `MPI_Allreduce` for consistent scaling factors.
- Double-allreduce of kinetic energy removed from `write_global_frame` in `gmd_main.cpp`.

---

## [v2.1] - 2026-05-22

### Added

- **MPI spatial domain decomposition** for distributed parallel simulations (initial 1D x-axis release, upgraded to full 3D in v2.2).
  - `DomainDecomposition` (`include/gmd/parallel/domain_decomposition.hpp`, `src/parallel/domain_decomposition.cpp`) provides 1D x-axis box splitting with configurable ghost width derived from cutoff + skin.
  - `MpiCommunicator` (`include/gmd/parallel/mpi_communicator.hpp`, `src/parallel/mpi_communicator.cpp`) encapsulates allreduce (scalar/vector), broadcast, barrier, ghost-atom coordinate exchange, and reverse force accumulation.
  - `MpiEnvironment` (`include/gmd/parallel/mpi_environment.hpp`, `src/parallel/mpi_environment.cpp`) provides RAII-style `MPI_Init`/`MPI_Finalize` with non-MPI fallback.
- **PME pencil decomposition** infrastructure.
  - `PmeParallelDecomposition` (`include/gmd/parallel/pme_parallel.hpp`, `src/parallel/pme_parallel.cpp`) defines a 2D process grid (Py × Pz) with transpose collectives (`x→y`, `y→z`, `z→y`, `y→x`). Currently the PME mesh is replicated; distributed FFT will be activated in a future release.
- **`GMD_ENABLE_MPI` CMake option** — opt-in MPI build with `find_package(MPI REQUIRED)`, `MPI::MPI_CXX` linkage, and `GMD_ENABLE_MPI` compile definition. See `cmake/MPIOptions.cmake` for the `gmd_configure_mpi_target` helper.
- **MPI smoke tests** — two new CTest targets (only registered when `GMD_ENABLE_MPI=ON`):
  - `gmd_smoke_mpi_lj_2proc`: runs a 2-process LJ simulation.
  - `gmd_smoke_mpi_lj_consistency`: compares serial and 2-process energy logs to verify numerical consistency.
  - New test collateral: `tests/smoke_mpi_lj.xyz`, `tests/smoke_mpi_lj.run`, `tests/compare_energy_logs.cpp`, `cmake/RunMpiLJConsistency.cmake`.
- **CLI `--np N` flag** — validates that the MPI world size matches the expected process count.
- **Rank-aware I/O** — only rank 0 writes trajectory and energy log files; global coordinates are gathered via `allreduce_vector` for output frames.
- **System MPI extensions** — `System` gains `atom_tags_`, `atom_owners_`, `num_local_atoms_` with corresponding accessors (`atom_tag()`, `atom_owner()`, `num_local_atoms()`, `mutable_atom_tags()`, `mutable_atom_owners()`). `resize()` accepts an optional `local_count` parameter.

### Changed (MPI correctness fixes)

- **`compute_twice_ke()`** (`src/integrator/thermostat.cpp`) — added `MPI_Allreduce` so all ranks return the same global kinetic energy. This ensures that temperature-dependent operations (thermostat scaling, barostat pressure) are consistent across ranks.
- **Nose-Hoover thermostat `initialize()`** (`src/integrator/nose_hoover_thermostat.cpp`) — `dof_` is now computed from the globally allreduced atom count so that the friction coefficient xi evolves identically on every rank.
- **Velocity rescaling thermostat `initialize()`** (`src/integrator/velocity_rescaling_thermostat.cpp`) — same global `dof_` fix.
- **`VelocityInitializer`** (`src/system/initializer.cpp`):
  - `remove_center_of_mass_velocity()` — added `MPI_Allreduce` for total mass and COM momentum so that all ranks subtract the same COM velocity.
  - `rescale_temperature()` — added `MPI_Allreduce` for kinetic energy and atom count so that the scaling factor is consistent.
  - `initialize()` — in MPI mode no longer early-returns when `atom_count() == 0`, because doing so would skip collective MPI calls and cause `MPI_ERR_TRUNCATE`.
- **`ClassicalForceProvider::compute()`** (`src/force/classical_force_provider.cpp`) — virial `MPI_Allreduce` moved outside the early-return path so that ranks with zero atoms still participate in the collective call.
- **`BondedForceProvider::compute()`** (`src/force/bonded_force_provider.cpp`) — added virial `MPI_Allreduce` for global bonded virial, also placed unconditionally after the computation block.
- **`EwaldForceProvider`** (`src/force/ewald_force_provider.cpp`):
  - `compute_reciprocal()` — structure factor S(k) built from local atoms only, then allreduced to obtain the global value; forces computed on local atoms only.
  - `compute_self_correction()` — q²_sum and Q_net computed from local atoms only, then allreduced.
  - `compute_real_space()` — added pair deduplication via atom_tag comparison for ghost-atom boundary pairs (same pattern as `ClassicalForceProvider`); outer i-loop limited to local atoms.
  - `compute()` — refactored `has_charges` early-exit into a coordinated check using `MPI_LOR` allreduce, ensuring all ranks reach the same decision.
- **`write_global_frame`** (`app/gmd_main.cpp`) — removed redundant `allreduce_scalar` wrapper around `compute_twice_ke`, which already performs its own allreduce. Previously this double-allreduced the kinetic energy in MPI builds.

### Changed

- `Simulation` gains `set_mpi_communicator()` and `set_domain_decomposition()`; `step()` now performs ghost exchange before force evaluation and reverse force accumulation afterward, with allreduce on potential energy.
- `ClassicalForceProvider::compute()` respects `num_local_atoms()` — only local atoms are iterated as primary index `i`; ghost atoms `j` still participate in pair evaluation.
- `VerletNeighborBuilder::rebuild()` builds neighbor lists for local atoms only; ghost atoms are included as secondary partners.
- `BondedForceProvider::compute()` uses owner-rank resolution for bonded terms spanning process boundaries.
- `app/gmd_main.cpp` initializes `MpiEnvironment`, distributes atoms via `keep_rank_local_atoms()`, sets up `DomainDecomposition`, and uses rank-gated logging throughout.

### Compatibility Notes

- All existing serial workflows are unaffected. Omitting `-DGMD_ENABLE_MPI=ON` produces an identical binary to previous releases.
- MPI and non-MPI builds share the same input file format; no input changes are required to run in parallel.
- The `--np N` flag is optional and only validated when MPI is active.

---

## [v2.0] - 2026-04-20

### Added

- **ML force provider integration** via TorchScript (LibTorch) backend.
  - New `TorchScriptModelRuntimeAdapter` (`include/gmd/force/torchscript_adapter.hpp`, `src/force/torchscript_adapter.cpp`) loads a `.pt` model exported by `gmd_se3gnn`, reads the `local_cutoff` attribute, and calls `forward(species, positions, edge_index, edge_shift)`.
  - `MLForceProvider` now exposes `float cutoff() const noexcept` (delegated to the adapter) so the correct `VerletNeighborBuilder` cutoff is set automatically.
  - Run files accept two new directives: `force_field ml` and `model_path /path/to/model.pt`.
- **`edge_index` and `edge_shift` tensors** for periodic-boundary ML models.
  - `NeighborList` stores per-pair integer image shift vectors (`image_flags [E,3]`) alongside the existing CSR neighbor array.
  - `VerletNeighborBuilder::rebuild()` computes and stores these shift vectors; `TorchScriptModelRuntimeAdapter` expands CSR half-pairs into a directed full graph and computes Cartesian `edge_shift = S * box.lengths` (Å).
- **Atomic numbers** (`Z`) populated per atom during XYZ loading.
  - `System` gains an `atomic_numbers_` array with `atomic_numbers()` / `mutable_atomic_numbers()` accessors.
  - `ConfigLoader::load_xyz` resolves element symbols to atomic numbers via a built-in `kElementAtomicNumbers` table and writes them into `System`.
  - `ModelEvaluationRequest` now carries an `atomic_numbers` span and a `neighbor_list` pointer, both forwarded by `MLForceProvider::compute()`.
- **`GMD_ENABLE_TORCH` CMake option**: when `ON`, `find_package(Torch REQUIRED)` is called, `src/force/torchscript_adapter.cpp` is compiled into `gmd_core`, and `${TORCH_LIBRARIES}` is linked. Build with `-DGMD_ENABLE_TORCH=ON -DCMAKE_PREFIX_PATH=/path/to/libtorch`.

### Changed

- `ModelRuntimeAdapter` interface extended with a default-impl `virtual float cutoff() const noexcept` and `ModelEvaluationRequest` fields for `atomic_numbers` and `neighbor_list`.
- `RunConfig` gains `force_field_type` and `ml_model_path` fields; the `force_field` parser in `ConfigLoader::load_run` now accepts `ml` in addition to `lj`.

### Compatibility Notes

- All existing `lj`, `molecular`, and inline force-field inputs are unaffected.
- The Torch backend is fully opt-in: omitting `-DGMD_ENABLE_TORCH=ON` produces an identical binary to previous releases, and a clear runtime error is emitted if a run file requests `force_field ml` on a non-Torch build.
- `image_flags` is a new field on `NeighborList`; code that constructs or clears `NeighborList` directly must call `clear()` (which now also clears `image_flags`).

## [v1.1] - 2026-04-20

### Added

- Configurable LJ cross-type mixing via `mixing_rule` / `mixing` for inline `run.in` force fields, standalone LJ `.ff` files, and molecular `.ff` files.
- Support for three built-in combining rules: `lorentz_berthelot` (default), `geometric`, and `waldman_hagler`.
- Flexible rule aliases and spellings, including `lb`, `geom`, `wh`, `lorentz-berthelot`, and `waldman-hagler`.

### Changed

- Explicit `pair` entries remain the highest-priority source for cross-type LJ parameters; the configured mixing rule is only used as a fallback.
- Non-bonded examples and documentation now describe `mixing_rule` directly, including the override behavior for explicit cross-pair definitions.

### Compatibility Notes

- Existing inputs remain compatible because the default mixing rule is still `lorentz_berthelot`.
- This release does not change the current molecular non-bonded safety model: external molecular force fields still default to bonded-only mode unless `molecular_nonbonded lj_unsafe` is requested explicitly.

## [v1.0] - 2026-04-19

### Added

- Bonded molecular mechanics support with harmonic bonds, harmonic angles, periodic dihedrals, and harmonic impropers.
- Topology parsing for `.top` files with bond, angle, dihedral, and improper sections.
- Molecular force-field loading for `.ff` files with per-type mass, epsilon, sigma, charge, bonded parameter tables, and explicit pair overrides.
- A complete ethane demo with `.xyz`, `.run`, `.ff`, and `.top` inputs plus a standalone build-and-run script.

### Notes

- `v1.0` established the first end-to-end molecular workflow in GMD on top of the existing LJ, neighbor-list, and integration infrastructure.
