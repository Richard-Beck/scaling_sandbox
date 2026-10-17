"""Shared output-path handling for config-driven campaigns."""

from __future__ import annotations

import os
from pathlib import Path


def output_dir(default_name: str) -> Path:
    """Return the runner-selected output directory or a local default."""
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    if configured:
        return Path(configured).expanduser().resolve()
    return Path.cwd() / "results" / default_name
