#!/usr/bin/env python3
"""Compare apex_verify output with the expected reproducibility metrics."""

from __future__ import annotations

import pathlib
import sys

EXPECTED = {
    "combined_rule_count": 747,
    "maximum_combined_charge": 5,
    "possible_bad_wheels_degree_7": 4438,
    "possible_bad_wheels_degree_8": 4939,
    "possible_bad_wheels_degree_9": 2409,
    "possible_bad_wheels_degree_10": 567,
    "possible_bad_wheels_degree_11": 38,
    "generated_island_occurrences_with_0_rings": 254,
    "generated_island_occurrences_with_1_rings": 88393,
    "generated_island_occurrences_with_2_rings": 20836,
    "generated_island_occurrences_with_3_rings": 18,
    "generated_island_occurrences_with_at_least_4_rings": 0,
}


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} METRICS.txt", file=sys.stderr)
        return 2
    values: dict[str, int] = {}
    for raw in pathlib.Path(sys.argv[1]).read_text(encoding="utf-8").splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        try:
            values[key.strip()] = int(value.strip())
        except ValueError:
            pass

    failed = False
    for key, expected in EXPECTED.items():
        actual = values.get(key)
        if actual != expected:
            print(f"MISMATCH {key}: expected {expected}, got {actual}")
            failed = True
        else:
            print(f"OK       {key}: {actual}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
