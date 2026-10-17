#!/usr/bin/env python3
"""Execute one or more indexed tasks from a polarity1d campaign manifest."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import load_manifest, run_manifest_task  # noqa: E402


def _indices(specification: str, count: int) -> list[int]:
    if specification == "all":
        return list(range(count))
    if ":" in specification:
        pieces = specification.split(":")
        if len(pieces) not in {2, 3}:
            raise ValueError("range must be START:STOP or START:STOP:STEP")
        start = int(pieces[0]) if pieces[0] else 0
        stop = int(pieces[1]) if pieces[1] else count
        step = int(pieces[2]) if len(pieces) == 3 and pieces[2] else 1
        return list(range(start, stop, step))
    return [int(specification)]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument(
        "--task",
        help="integer index, START:STOP[:STEP], or 'all'; defaults to SLURM_ARRAY_TASK_ID",
    )
    parser.add_argument("--force", action="store_true")
    arguments = parser.parse_args(argv)
    manifest = load_manifest(arguments.manifest)
    specification = arguments.task
    if specification is None:
        specification = os.environ.get("SLURM_ARRAY_TASK_ID")
    if specification is None:
        raise SystemExit("provide --task or run inside a SLURM array")
    for index in _indices(specification, manifest["task_count"]):
        destination = run_manifest_task(
            manifest,
            index,
            arguments.output_root,
            force=arguments.force,
        )
        print(destination)


if __name__ == "__main__":
    main()
