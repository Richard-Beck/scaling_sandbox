#!/usr/bin/env python3
"""Run and render the full spatial Buttenschön mechanochemical model."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from buttenschoen import SimulationParameters, figure9_regimes, simulate  # noqa: E402
from buttenschoen.render import plot_four_regimes, plot_kymograph  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--regime", default="all",
        choices=("all", *figure9_regimes().keys()),
        help="Figure 9 reconstruction regime to run (default: all)",
    )
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/buttenschoen_spatial/figure9")
    parser.add_argument("--n-cells", type=int, default=161)
    parser.add_argument("--end-time", type=float, default=360.0)
    parser.add_argument("--output-interval", type=float, default=1.0)
    parser.add_argument("--rtol", type=float, default=1e-5)
    parser.add_argument("--atol", type=float, default=1e-7)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    presets = figure9_regimes()
    selected = presets if args.regime == "all" else {args.regime: presets[args.regime]}
    simulation = SimulationParameters(
        n_cells=args.n_cells,
        end_time=args.end_time,
        output_interval=args.output_interval,
        rtol=args.rtol,
        atol=args.atol,
    )
    results = {}
    for name, (mechanics, initial) in selected.items():
        started = time.perf_counter()
        result = simulate(mechanics=mechanics, initial=initial, simulation=simulation)
        elapsed = time.perf_counter()-started
        directory = args.output / name
        result.save(directory)
        plot_kymograph(result, directory / "active_kymograph.png",
                       title=name.replace("_", " ").title())
        results[name] = result
        print(
            f"{name}: {elapsed:.2f}s, final length={result.length[-1]:.3f} um, "
            f"max mass error={result.relative_mass_error.max():.3e}"
        )
    if len(results) == 4:
        plot_four_regimes(results, args.output / "four_regime_kymographs.png")


if __name__ == "__main__":
    main()
