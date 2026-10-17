#!/usr/bin/env python3
"""Create focused configs for the representative profiles in WP1 Figure 1."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "campaigns" / "wp1_figure1_examples_v1"


def _source_task(campaign: str, *, assay: str, model: str, length: float, period: float | None = None) -> dict:
    """Select a task from its immutable campaign manifest, not prior results."""
    manifest = json.loads((ROOT / "campaigns" / f"{campaign}.json").read_text())
    for task in manifest["tasks"]:
        if task["assay"] != assay or task["model"] != model:
            continue
        if float(task["simulation"]["length"]) != length or int(task.get("replicate", 0)) != 0:
            continue
        if period is not None and float(task["protocol"]["period"]) != period:
            continue
        return task
    raise RuntimeError(f"no source task for {campaign}: {assay}, {model}, L={length}, period={period}")


def main() -> None:
    models = ("wave_pinning", "otsuji", "goryachev", "legi")
    configs = []
    for model in models:
        competition = _source_task(
            "wp1_revision_audit_v1", assay="competing_fronts_extended_4x",
            model=model, length=15.0,
        )
        competition["simulation"]["end_time"] = 20_000.0
        competition["simulation"]["output_interval"] = 20.0
        configs.append(("competition", model, competition))

        reversal = _source_task(
            "jilkine_reversal_long_conditioning_v1",
            assay="gradient_reversal_long_conditioning", model=model, length=15.0,
        )
        reversal_time = float(reversal["protocol"]["reversal_time"])
        reversal["simulation"]["end_time"] = reversal_time + 20_000.0
        reversal["simulation"]["output_interval"] = 20.0
        configs.append(("reversal", model, reversal))

        tracking = _source_task(
            "phase1_production_v1", assay="dynamic_tracking",
            model=model, length=15.0, period=200.0,
        )
        configs.append(("tracking", model, tracking))

    index = []
    for task_index, (assay, model, task) in enumerate(configs):
        config = {
            "model": task["model"],
            "model_parameters": task.get("model_parameters"),
            "protocol": task["protocol"],
            "simulation": task["simulation"],
            "metrics": task["metrics"],
        }
        destination = OUTPUT / "configs" / f"{task_index:02d}_{assay}_{model}.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
        index.append({"task_index": task_index, "assay": assay, "model": model,
                      "config": str(destination.relative_to(ROOT))})
    (OUTPUT / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
