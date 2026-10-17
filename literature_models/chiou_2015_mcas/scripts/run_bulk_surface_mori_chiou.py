#!/usr/bin/env python3
"""Run shared-chemistry Mori-like and Chiou-like bulk--surface conditions."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "giese_2015_bulk_surface_wave_pinning" / "src"))

from bulk_surface import SimulationConfig  # noqa: E402
from bulk_surface.render import render_comparison, render_run  # noqa: E402
from bulk_surface.shared_conditions import (  # noqa: E402
    SharedChemistry,
    baseline_geometry,
    chiou_condition,
    mori_condition,
)


DEFAULT_OUTPUT = ROOT / "dev" / "bulk_surface" / "mori_chiou_shared" / "v1"
SNAPSHOTS = (0.0, 50.0, 60.0, 70.0, 100.0, 200.0, 300.0)


def write_json(path: Path, values: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(values, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--duration", type=float, default=300.0)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--output-interval", type=float, default=1.0)
    parser.add_argument("--dx", type=float, default=0.1)
    parser.add_argument("--n-surface", type=int, default=256)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    if args.smoke:
        args.duration = 1.0
        args.dt = 0.02
        args.output_interval = 0.2
        args.dx = 0.2
        args.n_surface = 64

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "README.md").write_text(
        """# Shared-chemistry Mori/Chiou bulk--surface runs

These are Mori-like stimulation and Chiou-like competition conditions using
one common bulk--surface response and parameter set. They are not literal
reproductions of the native one-dimensional published models.

- `configs/` records shared and condition-specific inputs.
- `runs/<condition>/volume_1x/trajectory.npz` contains the full saved fields.
- `profiles.npz` contains surface fields and bulk-to-membrane traces.
- `snapshots/` contains bulk, membrane, and profile panels at selected times.
- `figures/` contains profiles, kymographs, and conservation diagnostics.
- `comparison/` contains cross-condition plots.

Regenerate from the repository root with:

```bash
conda run -n scaling_sandbox python scripts/run_bulk_surface_mori_chiou.py
```
""",
        encoding="utf-8",
    )
    geometry = baseline_geometry(dx=args.dx, n_surface=args.n_surface)
    chemistry = SharedChemistry()
    cases = {
        "mori": mori_condition(geometry, chemistry),
        "chiou": chiou_condition(geometry, chemistry),
    }
    config = SimulationConfig(
        end_time=args.duration, dt=args.dt, output_interval=args.output_interval
    )

    write_json(output / "configs" / "common.json", {
        "chemistry": chemistry.as_dict(),
        "geometry": {
            "volume_label": "1x",
            "perimeter": geometry.perimeter,
            "area": geometry.area,
            "radius": 10.0 / (2.0 * 3.141592653589793),
            "dx": args.dx,
            "n_surface": args.n_surface,
        },
        "simulation": asdict(config),
        "interpretation": (
            "Shared bulk-surface Mori chemistry; Mori-like transient cue and "
            "Chiou-like unforced 60:40 two-peak competition."
        ),
    })
    write_json(output / "configs" / "mori_1x.json", {
        "condition": "mori", "volume": "1x",
        "stimulus": {"kind": "wrapped_top_hat", "amplitude": 10.0,
                     "center_fraction": 0.0, "width_fraction": 0.04,
                     "start": 50.0, "stop": 70.0},
    })
    write_json(output / "configs" / "chiou_1x.json", {
        "condition": "chiou", "volume": "1x", "stimulus": None,
        "initial_seeds": {"centers": [0.25, 0.75], "mass_weights": [0.6, 0.4],
                          "width_fraction": 0.04, "peak_amplitude": 0.8},
    })

    results = {}
    for name, case in cases.items():
        run_dir = output / "runs" / name / "volume_1x"
        print(f"Running {name} condition...", flush=True)
        result = case.simulator.run(case.initial, config)
        result.save(run_dir)
        selected = tuple(t for t in SNAPSHOTS if t <= args.duration)
        if not selected or selected[-1] != args.duration:
            selected = (*selected, args.duration)
        diagnostics = render_run(
            result, run_dir, condition=name, snapshot_times=selected
        )
        print(f"{name}: {diagnostics}", flush=True)
        results[name] = result

    render_comparison(results, output / "comparison")
    print(f"Outputs written to {output}", flush=True)


if __name__ == "__main__":
    main()
