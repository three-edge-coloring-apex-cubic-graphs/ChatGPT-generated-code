#!/usr/bin/env python3
"""Run Lemma B.3 in independent filename-order shards and sum raw outputs."""

from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def run_one(
    shard: int,
    first: int,
    last: int,
    executable: Path,
    data_root: Path,
    output_dir: Path,
    check_reducibility: bool,
    literal_search: bool,
    use_reducibility_cache: bool,
) -> tuple[int, int]:
    log = output_dir / f"b3-shard-{shard:02d}.log"
    occurrences = output_dir / f"b3-shard-{shard:02d}.occurrences"
    command = [
        str(executable),
        "b3",
        "--data-root",
        str(data_root),
        "--first-configuration",
        str(first),
        "--last-configuration",
        str(last),
        "--island-occurrences",
        str(occurrences),
    ]
    if not check_reducibility:
        command.append("--skip-reducibility")
    if literal_search:
        command.append("--literal-search")
    if not use_reducibility_cache:
        command.append("--no-reducibility-cache")
    with log.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    return shard, completed.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, default=Path("build/apex_verify"))
    parser.add_argument("--data-root", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("b3-parallel-output"))
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--skip-reducibility", action="store_true")
    parser.add_argument(
        "--literal-search",
        action="store_true",
        help="disable safe B.4 pruning/memoization for a literal regression run",
    )
    parser.add_argument(
        "--no-reducibility-cache",
        action="store_true",
        help="recompute the semi-reducibility verdict for every occurrence",
    )
    args = parser.parse_args()

    if args.jobs < 1:
        parser.error("--jobs must be positive")
    configuration_dir = args.data_root / "configurations" / "K"
    count = len(sorted(configuration_dir.glob("K*.conf")))
    if count == 0:
        parser.error(f"no configurations found in {configuration_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    shard_count = min(args.jobs, count)
    ranges: list[tuple[int, int, int]] = []
    for shard in range(shard_count):
        first = (shard * count) // shard_count + 1
        last = ((shard + 1) * count) // shard_count
        ranges.append((shard + 1, first, last))

    failed = False
    with ThreadPoolExecutor(max_workers=shard_count) as executor:
        futures = {
            executor.submit(
                run_one,
                shard,
                first,
                last,
                args.executable,
                args.data_root,
                args.output_dir,
                not args.skip_reducibility,
                args.literal_search,
                not args.no_reducibility_cache,
            ): (shard, first, last)
            for shard, first, last in ranges
        }
        for future in as_completed(futures):
            shard, first, last = futures[future]
            _, returncode = future.result()
            print(f"shard {shard} ({first}-{last}) exited with status {returncode}")
            failed |= returncode != 0
    if failed:
        return 1

    occurrence_files = [
        args.output_dir / f"b3-shard-{shard:02d}.occurrences"
        for shard, _, _ in ranges
    ]
    merge = Path(__file__).with_name("merge_b3_occurrences.py")
    command = [sys.executable, str(merge)]
    command.extend(str(path) for path in occurrence_files)
    return subprocess.run(command).returncode


if __name__ == "__main__":
    raise SystemExit(main())
