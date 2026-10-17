#!/usr/bin/env python3
"""Source-tree entry point for the polarity1d config runner."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.runner import main  # noqa: E402


if __name__ == "__main__":
    main()

