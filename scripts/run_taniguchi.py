#!/usr/bin/env python3
"""Run published Taniguchi fixed-circle or Fig. 4F--I simulations.

This driver intentionally has no cell-size perturbation mode.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from taniguchi import (  # noqa: E402
    PHASE_PARAMETER_SETS,
    FixedParameters,
    simulate_fixed,
    simulate_phase_field,
)
from taniguchi.metrics import fixed_diagnostics  # noqa: E402


def save_fixed(path: Path, kk: float, end_time: float, seed: int) -> None:
    result = simulate_fixed(end_time, parameters=FixedParameters(K_k=kk), seed=seed,
                            output_interval=0.05)
    path.parent.mkdir(parents=True, exist_ok=True)
    summary = fixed_diagnostics(result)
    summary.update({"K_k": kk, "dt": 5e-4, "dx": 5.0 / 31.0,
                    "seed": seed, "end_time": end_time})
    np.savez_compressed(path, times=result.times, U=result.U, V=result.V,
                        x=result.grid.x, y=result.grid.y, mask=result.grid.mask,
                        summary_json=json.dumps(summary))
    print(json.dumps(summary, indent=2))


def save_phase(path: Path, parameter_set: str, end_time: float, seed: int) -> None:
    result = simulate_phase_field(end_time, parameter_set=parameter_set, seed=seed,
                                  output_interval=0.02)
    path.parent.mkdir(parents=True, exist_ok=True)
    area = np.sum(result.phi, axis=(1, 2)) * 0.1**2
    summary = {"parameter_set": parameter_set, "dt": 8e-5, "dx": 0.1,
               "seed": seed, "end_time": end_time,
               "event_count": result.event_count,
               "area_min": float(area.min()), "area_max": float(area.max()),
               "peak_V": float(result.V.max())}
    np.savez_compressed(path, times=result.times, phi=result.phi, U=result.U,
                        V=result.V, x=result.x, y=result.y,
                        summary_json=json.dumps(summary))
    print(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, choices=("fixed", "phase"))
    parser.add_argument("--condition", required=True,
                        help="K_K (fixed) or one of 4F, 4G, 4H, 4I (phase)")
    parser.add_argument("--end-time", type=float, required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "reports/output/taniguchi")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.model == "fixed":
        kk = float(args.condition)
        if kk not in {6.0, 5.7, 4.2, 2.5}:
            raise SystemExit("published fixed conditions are 6.0, 5.7, 4.2, 2.5")
        save_fixed(args.output / f"fixed_KK_{kk:.1f}.npz", kk,
                   args.end_time, args.seed)
    else:
        condition = args.condition.upper()
        if condition not in PHASE_PARAMETER_SETS:
            raise SystemExit("published phase conditions are 4F, 4G, 4H, 4I")
        save_phase(args.output / f"phase_{condition}.npz", condition,
                   args.end_time, args.seed)


if __name__ == "__main__":
    main()
