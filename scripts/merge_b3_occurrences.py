#!/usr/bin/env python3
"""Merge raw island-occurrence records from sharded Lemma B.3 runs.

Each input line contains only the number of boundary rings of one generated
island occurrence. Duplicate and isomorphic islands are deliberately counted
again. The default search may memoize an isomorphic subproblem, but it replays
the complete result
for every state and dart-pair occurrence, so these raw occurrence records retain
the Appendix B.4 multiplicity.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()

    counts: Counter[int] = Counter()
    total = 0
    for path in args.files:
        with path.open(encoding="utf-8") as stream:
            for line_number, raw in enumerate(stream, 1):
                raw = raw.strip()
                if not raw:
                    continue
                try:
                    rings = int(raw)
                except ValueError as error:
                    raise SystemExit(
                        f"{path}:{line_number}: malformed island-occurrence record"
                    ) from error
                if rings < 0:
                    raise SystemExit(
                        f"{path}:{line_number}: ring count must be nonnegative"
                    )
                counts[rings] += 1
                total += 1

    for rings in range(4):
        print(f"generated_island_occurrences_with_{rings}_rings={counts[rings]}")
    at_least_four = sum(count for rings, count in counts.items() if rings >= 4)
    print(
        "generated_island_occurrences_with_at_least_4_rings="
        f"{at_least_four}"
    )
    print(f"total_generated_island_occurrences={total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
