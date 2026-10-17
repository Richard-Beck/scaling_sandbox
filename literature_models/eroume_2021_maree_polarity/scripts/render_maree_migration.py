#!/usr/bin/env python3
"""Render a saved moving Eroume/Maree trajectory without rerunning it."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from maree_migration import (ChemistryParameters, MotilityParameters,  # noqa: E402
                             SimulationConfig, SimulationResult, StimulusProtocol)
from maree_migration.render import render_outputs  # noqa: E402


def load_result(directory: Path) -> SimulationResult:
    summary = json.loads((directory / "summary.json").read_text())
    raw = summary["configuration"]
    chemistry = ChemistryParameters(**raw.pop("chemistry"))
    motility = MotilityParameters(**raw.pop("motility"))
    stimulus = StimulusProtocol(**raw.pop("stimulus"))
    config = SimulationConfig(**raw, chemistry=chemistry, motility=motility,
                              stimulus=stimulus)
    with np.load(directory / "trajectory.npz", allow_pickle=False) as data:
        material = (data["material"] if "material" in data.files else
                    (data["phi"][:, None] * data["concentrations"]).sum((2, 3))
                    * config.dx**2)
        exterior = (data["exterior_material"] if "exterior_material" in data.files else
                    np.asarray([(phi[None] * fields)[:, phi < 0.05].sum(axis=1)
                                * config.dx**2
                                for phi, fields in zip(data["phi"], data["concentrations"])]))
        return SimulationResult(config, data["x"], data["y"], data["times"],
                                data["phi"], data["concentrations"], material, exterior)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fps", type=float, default=15.0)
    parser.add_argument("--pixels", type=int, default=640)
    parser.add_argument("--trajectory-overlay", action="store_true")
    args = parser.parse_args()
    result = load_result(args.input)
    output = args.output or args.input / "rendered"
    for path in render_outputs(result, output, fps=args.fps, pixels=args.pixels,
                               trajectory=args.trajectory_overlay):
        print(path)


if __name__ == "__main__":
    main()
