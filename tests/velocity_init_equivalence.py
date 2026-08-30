#!/usr/bin/env python3
"""Compare the random velocity field across rank counts and storage orders.

Runs tests/mpi_velocity_init.cpp at np=1, np=2, np=4 and again at np=2 with
each rank's local storage reversed, then compares every field against the np=1
one BY ATOM TAG.

Comparing by tag is the whole point. Under domain decomposition, array index i
is a different physical atom at every rank count, so an index-keyed comparison
would be meaningless -- and an index-keyed *generator* is precisely the defect
this exists to measure.

--require-match turns the comparison into an assertion. Without it the
differences are reported and only the structural invariants are enforced: the
same set of tags is present, every value is finite, the global momentum is
zero, and the kinetic energy matches the target temperature. Those hold however
the draws are assigned, which is why they can be checked before the fix lands.
"""

import argparse
import pathlib
import subprocess
import sys

problems: list[str] = []


def check(condition: bool, message: str) -> None:
    if not condition:
        problems.append(message)


def read_field(path: pathlib.Path) -> tuple[dict[int, tuple[float, float, float]], dict]:
    velocities: dict[int, tuple[float, float, float]] = {}
    meta: dict = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            continue
        fields = line.split()
        if fields[0] == "momentum":
            meta["momentum"] = tuple(float(v) for v in fields[1:4])
        elif fields[0] == "twice_kinetic_energy":
            meta["twice_ke"] = float(fields[1])
        elif fields[0] == "degrees_of_freedom":
            # Constrained fixtures emit the authoritative count. Unconstrained
            # ones omit it and the 3N-3 default below applies, which for them is
            # the same number.
            meta["dof"] = int(fields[1])
        else:
            tag = int(fields[0])
            if tag in velocities:
                problems.append(f"{path.name}: atom tag {tag} appears more than once")
            velocities[tag] = (float(fields[1]), float(fields[2]), float(fields[3]))
    return velocities, meta


def run(args, np_count: int, out: pathlib.Path, reverse: bool) -> None:
    command = [args.mpiexec, args.mpiexec_np_flag, str(np_count), args.executable,
               "--out", str(out)]
    command += args.extra_arg
    if reverse:
        command.append("--reverse-storage")
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"np={np_count} exited {result.returncode}\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}")


def worst_difference(a, b) -> tuple[float, int, int]:
    worst, worst_tag, differing = 0.0, -1, 0
    for tag, va in a.items():
        vb = b[tag]
        atom_worst = max(abs(x - y) for x, y in zip(va, vb))
        if atom_worst > 0.0:
            differing += 1
        if atom_worst > worst:
            worst, worst_tag = atom_worst, tag
    return worst, worst_tag, differing


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--mpiexec-np-flag", default="-n")
    parser.add_argument("--require-match", action="store_true",
                        help="fail when any field differs from the np=1 field by "
                             "more than --tolerance")
    # Not a physics tolerance and not slack for the random draws, which are
    # bitwise identical by tag. The centre-of-mass and kinetic-energy sums are
    # accumulated in storage order within a rank and combined by MPI_Allreduce
    # across ranks, and neither order is the same at every rank count, so the
    # one shared scale factor and the one shared centre-of-mass shift differ in
    # their last bits. Measured worst case across np=1/2/4 and a reversed
    # storage order is 2.8e-17 on velocities averaging 5.1e-02 -- well under an
    # ulp of a typical component. 1e-15 is ~36x that, and fourteen orders below
    # the 4.7e-01 the storage-ordered generator produced.
    parser.add_argument("--tolerance", type=float, default=1.0e-15)
    # k_B from the exact SI definitions, as elsewhere in this repository.
    parser.add_argument("--target-temperature", type=float, default=300.0)
    parser.add_argument("--extra-arg", action="append", default=[],
                        help="passed through to the executable on every run")
    args = parser.parse_args()

    boltzmann = 1.380649e-23 / 1.602176634e-19

    work = pathlib.Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)

    arrangements = [(1, False), (2, False), (4, False), (2, True)]
    fields = {}
    for np_count, reverse in arrangements:
        label = f"np{np_count}" + ("_reversed" if reverse else "")
        out = work / f"velocities_{label}.txt"
        run(args, np_count, out, reverse)
        fields[label] = read_field(out)

    reference, reference_meta = fields["np1"]
    check(len(reference) > 0, "the np=1 run produced no atoms")

    for label, (field, meta) in fields.items():
        # Structural: the same physical system, whatever the decomposition.
        check(set(field.keys()) == set(reference.keys()),
              f"{label}: the set of atom tags differs from np=1 "
              f"(missing {sorted(set(reference) - set(field))}, "
              f"extra {sorted(set(field) - set(reference))})")
        for tag, v in field.items():
            check(all(x == x and abs(x) != float("inf") for x in v),
                  f"{label}: atom {tag} has a non-finite velocity {v}")

        # The centre-of-mass constraint and the target temperature are global
        # properties and must hold at every rank count.
        momentum = meta["momentum"]
        scale = max(abs(x) for v in field.values() for x in v) or 1.0
        check(all(abs(p) < 1.0e-12 * scale * len(field) for p in momentum),
              f"{label}: residual centre-of-mass momentum {momentum}")
        dof = meta.get("dof", 3 * len(field) - 3)
        temperature = meta["twice_ke"] / (dof * boltzmann)
        check(abs(temperature / args.target_temperature - 1.0) < 1.0e-12,
              f"{label}: field carries {temperature!r} K, not "
              f"{args.target_temperature}")

    print(f"{'arrangement':16s} {'differing':>10s} {'worst component':>18s}")
    all_match = True
    for label, (field, _) in fields.items():
        if label == "np1":
            print(f"{label:16s} {'reference':>10s} {'-':>18s}")
            continue
        worst, worst_tag, differing = worst_difference(reference, field)
        print(f"{label:16s} {differing:>10d} {worst:>18.6e}"
              + (f"   (atom {worst_tag})" if worst_tag >= 0 else ""))
        if worst != 0.0:
            all_match = False
        if args.require_match and worst > args.tolerance:
            problems.append(
                f"{label}: the velocity field differs from np=1 in {differing} "
                f"of {len(reference)} atoms, worst component {worst:.6e} at "
                f"atom {worst_tag}, above the {args.tolerance:.1e} reduction "
                f"round-off bound. The random draw a physical atom receives must "
                f"not depend on rank count or local storage order")

    if problems:
        for message in problems:
            print(f"[velocity equivalence] {message}", file=sys.stderr)
        print(f"[velocity equivalence] {len(problems)} problem(s)", file=sys.stderr)
        return 1

    if all_match:
        print("[velocity equivalence] every arrangement reproduces the np=1 field "
              "bit for bit")
    elif args.require_match:
        print(f"[velocity equivalence] every arrangement reproduces the np=1 field "
              f"to within the {args.tolerance:.1e} reduction round-off bound")
    else:
        print("[velocity equivalence] structural invariants hold; field differences "
              "reported above are not asserted (no --require-match)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
