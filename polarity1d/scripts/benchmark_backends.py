#!/usr/bin/env python3
"""Benchmark candidate WP1 numerical strategies on representative single runs.

This deliberately does not launch a screen.  It compares adaptive sparse BDF,
an exact-diffusion Strang discretization, explicit RK45 on a small stiff case,
and the optimized de Belly finite-volume port.  The retained Julia/Gridap
smoke runner is optional because each fresh Julia process has a long startup.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
from time import perf_counter

import numpy as np
import scipy


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import (  # noqa: E402
    DeBelly,
    GaussianPulse,
    Goryachev,
    LEGI,
    LocalizedPulse,
    Otsuji,
    WavePinning,
    simulate,
)


def _timed(run, repeats: int) -> tuple[object, dict[str, float]]:
    elapsed = []
    result = None
    for _ in range(repeats):
        start = perf_counter()
        result = run()
        elapsed.append(perf_counter() - start)
    assert result is not None
    return result, {
        "median_wall_seconds": float(np.median(elapsed)),
        "minimum_wall_seconds": float(np.min(elapsed)),
        "maximum_wall_seconds": float(np.max(elapsed)),
    }


def _nrmse(candidate: np.ndarray, reference: np.ndarray) -> float:
    return float(
        np.sqrt(np.mean((candidate - reference) ** 2))
        / max(float(np.ptp(reference)), 1e-12)
    )


def benchmark_jilkine(repeats: int) -> list[dict[str, object]]:
    protocol = LocalizedPulse(
        amplitude=0.2,
        stop=40.0,
        taper=20.0,
        sides=("left", "right"),
        asymmetry=0.02,
    )
    records: list[dict[str, object]] = []
    for model in (WavePinning(), Otsuji(), Goryachev(), LEGI()):
        reference = simulate(
            model,
            protocol=protocol,
            cells=101,
            end_time=300.0,
            output_interval=5.0,
            backend="bdf",
            rtol=2e-8,
            atol=1e-11,
        )
        for backend in ("bdf", "split"):
            trajectory, timing = _timed(
                lambda backend=backend, model=model: simulate(
                    model,
                    protocol=protocol,
                    cells=101,
                    end_time=300.0,
                    output_interval=5.0,
                    backend=backend,
                ),
                repeats,
            )
            records.append(
                {
                    "case": "jilkine_dual_300s",
                    "model": model.name,
                    "strategy": backend,
                    "cells": 101,
                    "fixed_step": model.preferred_step if backend == "split" else None,
                    "simulated_seconds": 300.0,
                    "terminal_nrmse_vs_tight_bdf": _nrmse(
                        trajectory.observable[-1], reference.observable[-1]
                    ),
                    **timing,
                }
            )

    small_protocol = LocalizedPulse(amplitude=0.2, stop=10.0, taper=5.0)
    for model in (WavePinning(), Otsuji()):
        for backend in ("bdf", "rk45"):
            trajectory, timing = _timed(
                lambda backend=backend, model=model: simulate(
                    model,
                    protocol=small_protocol,
                    cells=41,
                    end_time=20.0,
                    output_interval=2.0,
                    backend=backend,
                ),
                1 if backend == "rk45" else repeats,
            )
            records.append(
                {
                    "case": "small_stiff_20s",
                    "model": model.name,
                    "strategy": backend,
                    "cells": 41,
                    "simulated_seconds": 20.0,
                    "rhs_evaluations": trajectory.metadata.get("rhs_evaluations"),
                    **timing,
                }
            )
    return records


def benchmark_debelly(repeats: int) -> list[dict[str, object]]:
    protocol = GaussianPulse()
    records: list[dict[str, object]] = []
    trajectories = {}
    for step in (0.5, 0.25, 0.125):
        trajectory, timing = _timed(
            lambda step=step: simulate(
                DeBelly(),
                protocol=protocol,
                cells=100,
                end_time=350.0,
                output_interval=5.0,
                step=step,
            ),
            repeats,
        )
        trajectories[step] = trajectory
        records.append(
            {
                "case": "debelly_single_350s",
                "model": trajectory.model,
                "strategy": "finite_volume_banded_staggered_be",
                "cells": 100,
                "fixed_step": step,
                "simulated_seconds": 350.0,
                "mca_conservation_error": trajectory.metadata["mca_conservation_error"],
                **timing,
            }
        )
    for step in (0.5, 0.25):
        record = next(
            item
            for item in records
            if item.get("fixed_step") == step and item["case"] == "debelly_single_350s"
        )
        record["terminal_nrmse_vs_next_refinement"] = _nrmse(
            trajectories[step].observable[-1], trajectories[step / 2].observable[-1]
        )

    spatial = {}
    for cells in (50, 100, 200):
        trajectory, timing = _timed(
            lambda cells=cells: simulate(
                DeBelly(),
                protocol=protocol,
                cells=cells,
                end_time=350.0,
                output_interval=5.0,
                step=0.25,
            ),
            repeats,
        )
        spatial[cells] = trajectory
        records.append(
            {
                "case": "debelly_spatial_refinement_350s",
                "model": trajectory.model,
                "strategy": "finite_volume_banded_staggered_be",
                "cells": cells,
                "fixed_step": 0.25,
                "simulated_seconds": 350.0,
                **timing,
            }
        )
    for cells in (50, 100):
        coarse, fine = spatial[cells], spatial[2 * cells]
        interpolated = np.interp(fine.x, coarse.x, coarse.observable[-1])
        record = next(
            item
            for item in records
            if item["case"] == "debelly_spatial_refinement_350s" and item["cells"] == cells
        )
        record["terminal_nrmse_vs_next_refinement"] = _nrmse(
            interpolated, fine.observable[-1]
        )
    return records


def benchmark_julia() -> dict[str, object]:
    runner = (
        REPOSITORY
        / "literature_models"
        / "debelly_2026_mechanochemical"
        / "scripts"
        / "run_debelly.sh"
    )
    environment = os.environ.copy()
    environment.update(
        {
            "JULIA": "/app/eb/software/Julia/1.8.2-linux-x86_64/bin/julia",
            "DEBELLY_QUIET": "1",
        }
    )
    start = perf_counter()
    completed = subprocess.run(
        [str(runner), "mechanochemical", "smoke"],
        cwd=REPOSITORY,
        env=environment,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    elapsed = perf_counter() - start
    runner_time = re.search(r"wall_seconds=(\d+)", completed.stdout)
    return {
        "case": "debelly_smoke_6s",
        "model": "debelly_1d_mechanochemical",
        "strategy": "julia_gridap_finite_element_fresh_process",
        "cells": 20,
        "fixed_step": 1.0,
        "simulated_seconds": 6.0,
        "median_wall_seconds": elapsed,
        "reported_runner_seconds": None if runner_time is None else int(runner_time.group(1)),
    }


def write_markdown(payload: dict[str, object], path: Path) -> None:
    records = payload["records"]
    lines = [
        "# WP1 backend benchmark",
        "",
        "Representative single trajectories only; no parameter screen was run.",
        "",
        "| Case | Model | Strategy | Cells | dt | Wall (s) | Error / next refinement |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for record in records:
        error = record.get(
            "terminal_nrmse_vs_tight_bdf",
            record.get("terminal_nrmse_vs_next_refinement", ""),
        )
        lines.append(
            "| {case} | {model} | {strategy} | {cells} | {step} | {wall:.4g} | {error} |".format(
                case=record["case"],
                model=record["model"],
                strategy=record["strategy"],
                cells=record["cells"],
                step=record.get("fixed_step") or "adaptive",
                wall=record["median_wall_seconds"],
                error="" if error == "" else f"{error:.4g}",
            )
        )
    lines.extend(
        [
            "",
            "Selected defaults: sparse adaptive BDF for Jilkine models and the",
            "conservative banded finite-volume solver at dt=0.5 s for de Belly.",
            "The exact-diffusion split solver remains as an independent regression oracle.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--include-julia", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "benchmarks" / "backend_benchmark.json",
    )
    arguments = parser.parse_args()
    if arguments.repeats < 1:
        raise SystemExit("--repeats must be positive")
    records = benchmark_jilkine(arguments.repeats) + benchmark_debelly(arguments.repeats)
    if arguments.include_julia:
        records.append(benchmark_julia())
    payload = {
        "environment": {
            "host": platform.node(),
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "conda_environment": os.environ.get("CONDA_DEFAULT_ENV"),
        },
        "repeats": arguments.repeats,
        "records": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    markdown = arguments.output.with_suffix(".md")
    write_markdown(payload, markdown)
    print(arguments.output)
    print(markdown)


if __name__ == "__main__":
    main()
