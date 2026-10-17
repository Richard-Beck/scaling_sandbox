"""Prepare matched 0.5x/1x/2x area runs from the unmodified M7682 XML."""

from __future__ import annotations

import hashlib
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "model/model.xml"
OUT = ROOT / "results/area_scaling_80k"
AREAS = (0.5, 1.0, 2.0)
SEEDS = tuple(range(19, 29))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "binary_sha256": hashlib.sha256((ROOT / "software/morpheus-2.3.9").read_bytes()).hexdigest(),
        "area_multipliers": AREAS,
        "seeds": SEEDS,
        "stop_time": 80000,
        "model_time_note": "Morpheus StopTime, with CPM MCSDuration=0.5",
        "scaling": {
            "target_area": "1200 * area_multiplier",
            "initial_radius": "20 * sqrt(area_multiplier)",
            "target_perimeter": "200 * sqrt(area_multiplier)",
            "area_constraint_strength": "0.0004 / area_multiplier",
            "fixed": "lattice, node length, membrane resolution, biochemical concentrations/rates/diffusion, surface strength, CPM kinetics, plotting and tracking intervals",
        },
        "runs": [],
    }
    for multiplier in AREAS:
        label = f"area_{multiplier:g}x"
        scale = math.sqrt(multiplier)
        for seed in SEEDS:
            run = OUT / label / f"seed_{seed}"
            run.mkdir(parents=True, exist_ok=True)
            tree = ET.parse(SOURCE)
            model = tree.getroot()
            model.find("./Time/StopTime").set("value", "80000")
            model.find("./Time/RandomSeed").set("value", str(seed))
            cell = model.find("./CellTypes/CellType")
            cell.find("VolumeConstraint").set("target", f"{1200 * multiplier:.12g}")
            cell.find("VolumeConstraint").set("strength", f"{0.0004 / multiplier:.12g}")
            cell.find("SurfaceConstraint").set("target", f"{200 * scale:.12g}")
            model.find(".//InitCellObjects/Arrangement/Sphere").set("radius", f"{20 * scale:.12g}")
            tree.write(run / "model.xml", encoding="utf-8", xml_declaration=True)
            manifest["runs"].append({
                "task_id": len(manifest["runs"]),
                "area_multiplier": multiplier,
                "seed": seed,
                "relative_dir": str(run.relative_to(OUT)),
            })
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {len(manifest['runs'])} runs in {OUT}")


if __name__ == "__main__":
    main()
