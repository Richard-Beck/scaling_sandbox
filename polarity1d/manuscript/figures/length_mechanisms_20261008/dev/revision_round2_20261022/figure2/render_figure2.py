"""Compatibility entry point; canonical runtime source is under manuscript/code.

The pending-data renderer is preserved as supporting metadata provenance.
"""
from pathlib import Path
import runpy
root = next(p for p in Path(__file__).resolve().parents if (p / "polarity1d/manuscript").is_dir())
runpy.run_path(str(root / "polarity1d/manuscript/code/mechanisms/revision_round2_figure2_20261022.py"), run_name="__main__")
