"""Config-driven runner for individual qualification runs (not production screens)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .metrics import phenotype_summary
from .protocols import (
    AlternatingGradient,
    GaussianPulse,
    FinitePulseSequence,
    LinearGradient,
    LocalizedPulse,
    NoStimulus,
    ReversingGradient,
    ReversingFiniteGradient,
)
from .registry import make_model
from .simulation import simulate


PROTOCOLS = {
    "none": NoStimulus,
    "localized_pulse": LocalizedPulse,
    "finite_pulse_sequence": FinitePulseSequence,
    "gaussian_pulse": GaussianPulse,
    "linear_gradient": LinearGradient,
    "reversing_gradient": ReversingGradient,
    "reversing_finite_gradient": ReversingFiniteGradient,
    "alternating_gradient": AlternatingGradient,
}


def run_config(config: dict[str, Any], output: Path) -> tuple[Path, Path]:
    model = make_model(config["model"], config.get("model_parameters"))
    protocol_config = dict(config.get("protocol", {"name": "none"}))
    protocol_name = protocol_config.pop("name")
    try:
        protocol = PROTOCOLS[protocol_name](**protocol_config)
    except KeyError as error:
        raise ValueError(f"unknown protocol {protocol_name!r}") from error
    simulation_config = dict(config.get("simulation", {}))
    trajectory = simulate(model, protocol=protocol, **simulation_config)
    metric_config = dict(config.get("metrics", {}))
    summary = phenotype_summary(trajectory, **metric_config)
    summary["trajectory_metadata"] = trajectory.metadata
    output.mkdir(parents=True, exist_ok=True)
    trajectory_path = trajectory.save(output / "trajectory.npz")
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return trajectory_path, summary_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args(argv)
    config = json.loads(arguments.config.read_text())
    trajectory, summary = run_config(config, arguments.output)
    print(trajectory)
    print(summary)


if __name__ == "__main__":
    main()
