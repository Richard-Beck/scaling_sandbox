#!/usr/bin/env python3
"""Run one frozen-boundary, controlled-wave, or deformable campaign case."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from taniguchi.geometry_campaign import run_geometry_campaign  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    path = args.config if args.config.is_absolute() else ROOT / args.config
    config = json.loads(path.read_text(encoding="utf-8"))
    output = Path(config["output_path"])
    if not output.is_absolute():
        output = ROOT / output
    started = time.perf_counter()
    summary = run_geometry_campaign(config, output)
    summary["elapsed_seconds"] = time.perf_counter() - started
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(output),
        "arm": summary["arm"],
        "volume_ratio": summary["volume_ratio"],
        "seed": summary["seed"],
        "event_count": summary["event_count"],
        "wave_tracks": summary["wave_track_count"],
        "elapsed_seconds": summary["elapsed_seconds"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
