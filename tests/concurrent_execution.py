#!/usr/bin/env python3
"""Run several copies of a test executable at once and require all to succeed.

A test that writes a fixed filename into the shared system temporary directory
passes when run alone and fails intermittently when it does not. That is the
worst kind of failure: it appears in CI, blames the code under test, and does
not reproduce. The same test binary exists in every configured build tree --
Release and Debug, serial and MPI-enabled -- so "run alone" is not something
the harness can promise.

This driver makes the concurrent case explicit rather than incidental.
"""

import argparse
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--copies", type=int, default=4)
    args = parser.parse_args()

    processes = [
        subprocess.Popen([args.executable], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True)
        for _ in range(args.copies)
    ]
    failed = []
    for index, process in enumerate(processes):
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            failed.append((index, process.returncode, stdout, stderr))

    if failed:
        for index, code, stdout, stderr in failed:
            print(f"[concurrent] copy {index} exited {code}", file=sys.stderr)
            if stdout:
                print(f"--- stdout ---\n{stdout}", file=sys.stderr)
            if stderr:
                print(f"--- stderr ---\n{stderr}", file=sys.stderr)
        print(f"[concurrent] {len(failed)} of {args.copies} concurrent copies of "
              f"{args.executable} failed; the test is not safe to run alongside "
              f"itself", file=sys.stderr)
        return 1

    print(f"[concurrent] {args.copies} concurrent copies of {args.executable} "
          f"all passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
