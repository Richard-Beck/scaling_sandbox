#!/usr/bin/env python3
"""Collect a polarity1d production campaign into task and group tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import collect_campaign, load_manifest  # noqa: E402


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_root", type=Path)
    arguments = parser.parse_args(argv)
    report = collect_campaign(load_manifest(arguments.manifest), arguments.output_root)
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["failure_count"] or report["invalid_count"]:
        raise SystemExit(1)
    if report["missing_count"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
