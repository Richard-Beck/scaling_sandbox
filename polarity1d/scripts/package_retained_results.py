#!/usr/bin/env python3
"""Assemble retained PNGs and compact derived data after figure generation."""

from __future__ import annotations

import json
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
DATA = RESULTS / "data"
MANUSCRIPT = RESULTS / "manuscript_wp1_v3"

FIGURE_SOURCES = {
    "Figure1_WP1_cross_model.png": MANUSCRIPT / "Figure1_WP1_cross_model.png",
    "FigureS1_Otsuji_reversal_sweep.png": MANUSCRIPT / "FigureS1_Otsuji_reversal_sweep.png",
    "FigureS2_tracking_examples.png": MANUSCRIPT / "FigureS2_tracking_examples.png",
    "Figure2_standalone_v2_julia_activation_defaults_linear50.png": (
        MANUSCRIPT / "standalone_alternatives/v2_julia_activation_defaults_linear50/"
        "Figure2_standalone_v2_julia_activation_defaults_linear50.png"
    ),
    "Figure2_supplement_full_model_kymographs.png": (
        MANUSCRIPT / "standalone_alternatives/v2_julia_activation_defaults_linear50/supplement/"
        "Figure2_supplement_full_model_kymographs.png"
    ),
    "debelly_qualification.png": RESULTS / "qualification/debelly_qualification.png",
    "jilkine_qualification.png": RESULTS / "qualification/jilkine_qualification.png",
}


def copy_files(source: Path, destination: Path, names: tuple[str, ...] | None = None) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    paths = sorted(source.iterdir()) if names is None else [source / name for name in names]
    for path in paths:
        if path.is_file():
            shutil.copy2(path, destination / path.name)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for name, source in FIGURE_SOURCES.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, FIGURES / name)

    copy_files(MANUSCRIPT / "source_data", DATA / "manuscript")
    linear = MANUSCRIPT / "standalone_alternatives/v2_julia_activation_defaults_linear50"
    copy_files(linear, DATA / "figure2", (
        "README.md", "figure2_v2_countercue_repolarization_audit.csv",
        "figure2_v2_minimum_countercue_duration.csv",
        "figure2_v2_polarity_trajectories.csv.gz",
        "figure2_v2_steady_state_metrics.csv",
    ))
    copy_files(linear / "supplement", DATA / "figure2_supplement", (
        "README.md", "figure2_supplement_kymograph_data.npz",
        "figure2_supplement_orientation_trajectories.csv.gz",
    ))
    copy_files(RESULTS / "qualification", DATA / "qualification", (
        "qualification_summary.json",
    ))
    report = {
        "checkpoint": "accepted",
        "comparison_basis": "decoded RGBA pixel arrays",
        "pixel_identical": [
            "Figure1_WP1_cross_model.png",
            "FigureS1_Otsuji_reversal_sweep.png",
            "FigureS2_tracking_examples.png",
            "Figure2_standalone_v2_julia_activation_defaults_linear50.png",
            "Figure2_supplement_full_model_kymographs.png",
            "jilkine_qualification.png",
        ],
        "accepted_difference": {
            "debelly_qualification.png": {"differing_pixels": 787769},
        },
    }
    (DATA / "checkpoint_pixel_comparison.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(FIGURES)


if __name__ == "__main__":
    main()
