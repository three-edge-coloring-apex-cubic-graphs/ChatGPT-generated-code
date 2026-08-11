#!/usr/bin/env python3
"""Run the five center-degree cases of Lemma B.2 in parallel processes."""

from __future__ import annotations

import argparse
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

EXPECTED = {7: 4438, 8: 4939, 9: 2409, 10: 567, 11: 38}
METRIC = re.compile(r"^possible_bad_wheels_degree_(\d+)=(\d+)$")


def run_degree(degree: int, executable: Path, data_root: Path, output_dir: Path) -> tuple[int, int]:
    log = output_dir / f"b2-degree-{degree}.log"
    command = [
        str(executable),
        "b2",
        "--data-root",
        str(data_root),
        "--first-degree",
        str(degree),
        "--last-degree",
        str(degree),
    ]
    with log.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    return degree, completed.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, default=Path("build/apex_verify"))
    parser.add_argument("--data-root", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("b2-parallel-output"))
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--check-metrics", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    failed = False
    degrees = sorted(EXPECTED)
    with ThreadPoolExecutor(max_workers=min(args.jobs, len(degrees))) as executor:
        futures = {
            executor.submit(run_degree, d, args.executable, args.data_root, args.output_dir): d
            for d in degrees
        }
        for future in as_completed(futures):
            degree, returncode = future.result()
            print(f"degree {degree} exited with status {returncode}")
            failed |= returncode != 0
    if failed:
        return 1

    actual: dict[int, int] = {}
    for degree in degrees:
        for line in (args.output_dir / f"b2-degree-{degree}.log").read_text().splitlines():
            match = METRIC.match(line)
            if match:
                actual[int(match.group(1))] = int(match.group(2))
    for degree in degrees:
        if degree not in actual:
            print(f"missing metric for degree {degree}")
            failed = True
            continue
        print(f"possible_bad_wheels_degree_{degree}={actual[degree]}")
        if args.check_metrics and actual[degree] != EXPECTED[degree]:
            print(f"mismatch for degree {degree}: expected {EXPECTED[degree]}")
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
