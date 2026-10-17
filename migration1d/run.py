#!/usr/bin/env python3
"""Single command-line entry point for all migration1d campaigns."""

from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT / "src"))

from migration1d.runner import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
