"""Render a compact directional-behaviour atlas for focused length sweeps.

Each trajectory is reduced to cumulative leftward and rightward centre travel.
The dominant direction defines a three-class hue (left, right, or both), while
the unnormalised total travel in micrometres defines colour intensity.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch, Rectangle
import numpy as np


PROTOCOLS = ("single", "dual_equal", "dual_biased", "reversal")
MECHANICAL_CASES = ("soft", "stiff")
BASE_COLOURS = {
    "left": "#2166AC",
    "right": "#B2182B",
    "both": "#762A83",
}


@dataclass(frozen=True)
class Condition:
    chemistry: str
    tension_mode: str
    relative_path: str
    title: str


CONDITIONS = (
    Condition(
        "single_species", "global",
        "single_species/testing/focused_motile_length",
        "Single species\nGlobal tension",
    ),
    Condition(
        "single_species", "spatial",
        "single_species/testing/focused_motile_length_spatial_tension",
        "Single species\nSpatial tension",
    ),
    Condition(
        "three_species", "global",
        "three_species/testing/focused_motile_length",
        "Three species\nGlobal tension",
    ),
    Condition(
        "three_species", "spatial",
        "three_species/testing/focused_motile_length_spatial_tension",
        "Three species\nSpatial tension",
    ),
)


def directional_travel(times: np.ndarray, velocity: np.ndarray) -> tuple[float, float]:
    """Return cumulative leftward and rightward centre travel in micrometres."""
    left = float(np.trapezoid(np.maximum(-velocity, 0.0), times))
    right = float(np.trapezoid(np.maximum(velocity, 0.0), times))
    return left, right


def classify_direction(left: float, right: float, dominance: float) -> str:
    total = left + right
    if total <= np.finfo(float).eps:
        return "both"
    if left / total >= dominance:
        return "left"
    if right / total >= dominance:
        return "right"
    return "both"


def load_records(results_root: Path, dominance: float) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for condition_index, condition in enumerate(CONDITIONS):
        root = results_root / condition.relative_path
        task_dir = root / "tasks"
        json_paths = sorted(task_dir.glob("task_*.json"))
        if not json_paths:
            raise FileNotFoundError(f"no task JSON files found in {task_dir}")
        for json_path in json_paths:
            payload = json.loads(json_path.read_text())
            chemistry, length, mechanics, protocol, alpha = payload["task_key"]
            if chemistry != condition.chemistry:
                raise ValueError(f"unexpected chemistry in {json_path}: {chemistry}")
            trajectory_path = task_dir / payload["trajectory"]
            with np.load(trajectory_path) as trajectory:
                left, right = directional_travel(
                    trajectory["times"], trajectory["center_velocity"]
                )
            total = left + right
            direction = classify_direction(left, right, dominance)
            records.append({
                "condition_index": condition_index,
                "condition": condition.title.replace("\n", " / "),
                "chemistry_class": chemistry,
                "tension_mode": condition.tension_mode,
                "task_index": payload["task_index"],
                "reference_length_um": float(length),
                "mechanical_case": mechanics,
                "protocol": protocol,
                "tension_feedback_alpha": float(alpha),
                "leftward_travel_um": left,
                "rightward_travel_um": right,
                "total_travel_um": total,
                "leftward_fraction": left / total if total > 0.0 else 0.5,
                "rightward_fraction": right / total if total > 0.0 else 0.5,
                "direction_class": direction,
            })
    return records


def intensity_colour(direction: str, total: float, cap: float) -> tuple[float, float, float]:
    fraction = np.clip(total / cap, 0.0, 1.0) if cap > 0.0 else 0.0
    base = np.asarray(to_rgb(BASE_COLOURS[direction]))
    return tuple(1.0 - fraction * (1.0 - base))


def write_records(path: Path, records: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(records[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def render_atlas(
    output: Path,
    records: list[dict[str, object]],
    dominance: float,
    cap_quantile: float,
) -> float:
    lengths = sorted({float(record["reference_length_um"]) for record in records})
    alphas = sorted({float(record["tension_feedback_alpha"]) for record in records})
    length_index = {value: index for index, value in enumerate(lengths)}
    alpha_index = {value: index for index, value in enumerate(alphas)}
    totals = np.asarray([float(record["total_travel_um"]) for record in records])
    cap = float(np.quantile(totals, cap_quantile))

    lookup = {
        (
            int(record["condition_index"]), record["protocol"],
            float(record["reference_length_um"]),
            float(record["tension_feedback_alpha"]), record["mechanical_case"],
        ): record
        for record in records
    }

    fig, axes = plt.subplots(
        len(PROTOCOLS), len(CONDITIONS), figsize=(20, 13.2),
        sharex=True, sharey=True,
    )
    for row_index, protocol in enumerate(PROTOCOLS):
        for column_index, condition in enumerate(CONDITIONS):
            ax = axes[row_index, column_index]
            for alpha in alphas:
                y = alpha_index[alpha]
                for length in lengths:
                    x = length_index[length]
                    for mechanics, y_offset in (("stiff", -0.5), ("soft", 0.0)):
                        record = lookup[
                            (column_index, protocol, length, alpha, mechanics)
                        ]
                        colour = intensity_colour(
                            str(record["direction_class"]),
                            float(record["total_travel_um"]), cap,
                        )
                        ax.add_patch(Rectangle(
                            (x - 0.5, y + y_offset), 1.0, 0.5,
                            facecolor=colour, edgecolor="white", linewidth=0.45,
                        ))

            ax.set_xlim(-0.5, len(lengths) - 0.5)
            ax.set_ylim(-0.5, len(alphas) - 0.5)
            ax.set_xticks(range(len(lengths)))
            ax.set_yticks(range(len(alphas)))
            ax.grid(False)
            for spine in ax.spines.values():
                spine.set_color("#555555")
                spine.set_linewidth(0.7)
            if row_index == 0:
                ax.set_title(condition.title, fontsize=13, pad=10)
            if column_index == 0:
                ax.set_yticklabels([f"{alpha:g}" for alpha in alphas])
                ax.set_ylabel(
                    protocol.replace("_", " ").title() + "\n"
                    + r"Tension feedback $\alpha_t$",
                    fontsize=11,
                )
            if row_index == len(PROTOCOLS) - 1:
                ax.set_xticklabels(
                    [f"{length:g}" for length in lengths], rotation=45,
                    ha="right", rotation_mode="anchor",
                )
                ax.set_xlabel(r"Reference length $L_0$ ($\mu$m)", fontsize=11)

    fig.suptitle(
        "Directional motility atlas",
        fontsize=18, fontweight="bold", y=0.985,
    )
    fig.text(
        0.5, 0.952,
        "Each tile: soft mechanics on top, stiff mechanics below; "
        "intensity is cumulative centre travel (not length-normalised)",
        ha="center", va="center", fontsize=11,
    )

    direction_handles = [
        Patch(facecolor=BASE_COLOURS[label], edgecolor="none", label=label.title())
        for label in ("left", "both", "right")
    ]
    fig.legend(
        handles=direction_handles, loc="lower left", bbox_to_anchor=(0.19, 0.012),
        ncol=3, frameon=False, title=f"Direction class ({dominance:.0%} rule)",
    )

    strength_ax = fig.add_axes((0.54, 0.028, 0.25, 0.018))
    gradient = np.linspace(0.0, 1.0, 256)[None, :]
    both_base = np.asarray(to_rgb(BASE_COLOURS["both"]))
    rgb = 1.0 - gradient[..., None] * (1.0 - both_base)
    strength_ax.imshow(rgb, aspect="auto", extent=(0.0, cap, 0.0, 1.0))
    strength_ax.set_yticks([])
    strength_ax.set_xlabel(
        rf"Cumulative travel ($\mu$m); colours saturate at the "
        rf"{cap_quantile:.0%} quantile ({cap:.2f} $\mu$m)",
        fontsize=9,
    )
    strength_ax.tick_params(axis="x", labelsize=8)
    for spine in strength_ax.spines.values():
        spine.set_color("#555555")
        spine.set_linewidth(0.6)

    fig.subplots_adjust(
        left=0.08, right=0.985, top=0.91, bottom=0.115,
        wspace=0.10, hspace=0.16,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, facecolor="white")
    fig.savefig(output.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)
    return cap


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root", type=Path, default=Path("results"),
        help="Root containing single_species and three_species results",
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/focused_motile_length_directional_atlas.png"),
    )
    parser.add_argument(
        "--dominance", type=float, default=0.8,
        help="Travel fraction required for a left-only or right-only class",
    )
    parser.add_argument(
        "--cap-quantile", type=float, default=0.95,
        help="Quantile at which colour intensity reaches full saturation",
    )
    args = parser.parse_args()
    if not 0.5 < args.dominance <= 1.0:
        parser.error("--dominance must be in (0.5, 1]")
    if not 0.0 < args.cap_quantile <= 1.0:
        parser.error("--cap-quantile must be in (0, 1]")

    records = load_records(args.results_root, args.dominance)
    csv_path = args.output.with_suffix(".csv")
    write_records(csv_path, records)
    cap = render_atlas(args.output, records, args.dominance, args.cap_quantile)
    print(f"wrote {args.output}")
    print(f"wrote {args.output.with_suffix('.pdf')}")
    print(f"wrote {csv_path}")
    print(f"classified {len(records)} trajectories; intensity cap={cap:.6g} um")


if __name__ == "__main__":
    main()
