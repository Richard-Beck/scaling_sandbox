#!/usr/bin/env python3
"""Build an aggregate HTML report from completed single-species analyses."""

from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT / "src"))

from migration1d.single_species.report import main  # noqa: E402


if __name__ == "__main__":
    main()
