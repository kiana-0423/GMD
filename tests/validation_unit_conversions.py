#!/usr/bin/env python3
"""Unit tests for the time conversions inside validation/common.py.

The validation analyzers read a log whose time column is in femtoseconds and
report a diffusion coefficient in A^2/ps. Exactly one factor of 1000 stands
between those two, inside fit_diffusion_coefficient(), and nothing else in the
repository exercises it: the diffusion case compares one fitted number against a
stored one, so a conversion that was wrong by 1000 in both the reference and the
run would agree with itself forever.

These cases feed synthetic series with an analytically known slope, so the
expected value is arithmetic rather than a previous output.
"""

import argparse
import pathlib
import sys

failures = []


def check(condition: bool, message: str) -> None:
    if not condition:
        failures.append(message)


def relative(actual: float, expected: float) -> float:
    if expected == 0.0:
        return abs(actual)
    return abs(actual - expected) / abs(expected)


def test_diffusion_converts_fs_to_ps(common) -> None:
    # MSD = 6 D t with t in ps. Build a series in FEMTOSECONDS whose slope
    # corresponds to a known D in A^2/ps, and require the fit to recover it.
    #
    # D = 2.5 A^2/ps  =>  MSD grows by 6*2.5 = 15 A^2 per ps = 0.015 A^2 per fs.
    expected_d = 2.5
    per_fs = 6.0 * expected_d * 1.0e-3
    series = [{"time_fs": float(i) * 10.0, "msd_a2": per_fs * float(i) * 10.0}
              for i in range(200)]
    fitted = common.fit_diffusion_coefficient(series)
    check(relative(fitted, expected_d) < 1.0e-12,
          f"fit_diffusion_coefficient returned {fitted!r} A^2/ps for a series "
          f"built with D = {expected_d} A^2/ps")

    # The named failure modes, each a specific wrong constant rather than a
    # generic tolerance: treating fs as ps is 1000x low, the reciprocal is
    # 1000x high, and omitting the 1/6 is 6x high.
    check(relative(fitted, expected_d * 1.0e-3) > 0.1,
          "the diffusion fit is 1000x low: femtoseconds are being used as if "
          "they were picoseconds")
    check(relative(fitted, expected_d * 1.0e3) > 0.1,
          "the diffusion fit is 1000x high: the fs-to-ps conversion is inverted")
    check(relative(fitted, expected_d * 6.0) > 0.1,
          "the diffusion fit is 6x high: the Einstein relation's 1/6 is missing")


def test_diffusion_is_linear_in_the_slope(common) -> None:
    # Doubling the slope must double D exactly; an additive offset in the time
    # axis must not change it at all.
    def fit(scale: float, offset_fs: float) -> float:
        series = [{"time_fs": offset_fs + float(i) * 10.0,
                   "msd_a2": scale * 0.015 * float(i) * 10.0}
                  for i in range(200)]
        return common.fit_diffusion_coefficient(series)

    base = fit(1.0, 0.0)
    check(relative(fit(2.0, 0.0), 2.0 * base) < 1.0e-12,
          "doubling the MSD slope did not double the diffusion coefficient")
    check(relative(fit(1.0, 5000.0), base) < 1.0e-12,
          "shifting the time origin changed the diffusion coefficient")


def test_zero_and_degenerate_series(common) -> None:
    check(common.fit_diffusion_coefficient([]) == 0.0,
          "an empty series must fit to exactly zero")
    check(common.fit_diffusion_coefficient([{"time_fs": 0.0, "msd_a2": 0.0}]) == 0.0,
          "a single-point series must fit to exactly zero")
    flat = [{"time_fs": float(i), "msd_a2": 3.0} for i in range(50)]
    check(common.fit_diffusion_coefficient(flat) == 0.0,
          "a flat MSD must fit to exactly zero diffusion")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation-dir", required=True,
                        help="directory containing common.py")
    args = parser.parse_args()

    sys.path.insert(0, str(pathlib.Path(args.validation_dir).resolve()))
    import common  # noqa: E402

    test_diffusion_converts_fs_to_ps(common)
    test_diffusion_is_linear_in_the_slope(common)
    test_zero_and_degenerate_series(common)

    if failures:
        for message in failures:
            print(f"[validation units] {message}", file=sys.stderr)
        print(f"[validation units] {len(failures)} check(s) failed", file=sys.stderr)
        return 1
    print("[validation units] diffusion time-unit conversion verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
