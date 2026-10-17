#!/usr/bin/env python3
"""Generate the complete manuscript figure set from the migration1d results.

The script uses the merged condition-level CSV files for population summaries,
the retained NPZ files for moving-cell trajectories, and reproducibly reruns a
small set of fixed-domain representative trajectories whose fields were not
retained by the original campaigns.  All publication outputs are PNG files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mpl_colors
from matplotlib import patches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
MIGRATION = PROJECT / "migration1d"
RESULTS = MIGRATION / "results"
OUTPUT = HERE / "figures"
CACHE = HERE / ".figure_cache"

if str(MIGRATION / "src") not in sys.path:
    sys.path.insert(0, str(MIGRATION / "src"))


# Okabe-Ito-derived colors, reused for the same concepts throughout.
INK = "#17212B"
MUTED = "#5F6B73"
GRID = "#DDE3E7"
PAPER = "#FFFFFF"
PANEL_BG = "#F7F9FA"
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
MAGENTA = "#CC79A7"
PURPLE = "#7A5195"
YELLOW = "#F0E442"
GRAY = "#A6ADB2"

VOLUME_COLORS = {0.5: SKY, 1.0: BLUE, 2.0: VERMILLION}
MECH_COLORS = {"soft": VERMILLION, "stiff": BLUE}
MODEL_COLORS = {"single_species": BLUE, "three_species": MAGENTA}
FEEDBACK_COLORS = {"none": INK, "global": ORANGE, "spatial": PURPLE}

PHENOTYPE_ORDER = ["homogeneous", "single", "dual", "other"]
PHENOTYPE_COLORS = {
    "homogeneous": "#D9DEE2",
    "single": BLUE,
    "dual": MAGENTA,
    "other": ORANGE,
}
REVERSAL_ORDER = [
    "reversed",
    "retained_original_pole",
    "multiple_polarized_domains",
    "lost_polarity",
    "ambiguous",
    "not_initially_left_polarized",
]
REVERSAL_LABELS = {
    "reversed": "clean reversal",
    "retained_original_pole": "old front retained",
    "multiple_polarized_domains": "dual / multiple",
    "lost_polarity": "depolarized",
    "ambiguous": "other",
    "not_initially_left_polarized": "not eligible",
}
REVERSAL_COLORS = {
    "reversed": GREEN,
    "retained_original_pole": BLUE,
    "multiple_polarized_domains": MAGENTA,
    "lost_polarity": GRAY,
    "ambiguous": ORANGE,
    "not_initially_left_polarized": "#E5E8EA",
}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.titlesize": 10.0,
        "axes.labelsize": 9.0,
        "axes.titleweight": "semibold",
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "axes.linewidth": 0.75,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "legend.frameon": False,
        "figure.facecolor": PAPER,
        "savefig.facecolor": PAPER,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.08,
        "mathtext.fontset": "dejavusans",
    }
)


@dataclass(frozen=True)
class DataBundle:
    one_single: pd.DataFrame
    one_dual: pd.DataFrame
    one_reversal: pd.DataFrame
    three_single: pd.DataFrame
    three_dual: pd.DataFrame
    three_reversal: pd.DataFrame
    moving: dict[tuple[str, str], pd.DataFrame]
    directional: pd.DataFrame
    one_roots: pd.DataFrame
    one_regime: pd.DataFrame
    one_dual_calibration: pd.DataFrame
    one_reversal_calibration: pd.DataFrame
    three_regime: pd.DataFrame


def read_data() -> DataBundle:
    one = RESULTS / "single_species"
    three = RESULTS / "three_species"
    moving = {
        ("single_species", "global"): pd.read_csv(
            one / "testing/focused_motile_length/focused_motile_length.csv"
        ),
        ("single_species", "spatial"): pd.read_csv(
            one
            / "testing/focused_motile_length_spatial_tension/focused_motile_length.csv"
        ),
        ("three_species", "global"): pd.read_csv(
            three / "testing/focused_motile_length/focused_motile_length.csv"
        ),
        ("three_species", "spatial"): pd.read_csv(
            three
            / "testing/focused_motile_length_spatial_tension/focused_motile_length.csv"
        ),
    }
    return DataBundle(
        one_single=pd.read_csv(one / "testing/wgd_null/wgd_null_screen.csv"),
        one_dual=pd.read_csv(
            one / "testing/wgd_null/double_stimulus/two_seeded_domains.csv"
        ),
        one_reversal=pd.read_csv(
            one / "testing/wgd_null/polarity_reversal/polarity_reversal.csv"
        ),
        three_single=pd.read_csv(three / "testing/wgd_null/wgd_null_screen.csv"),
        three_dual=pd.read_csv(
            three / "testing/wgd_null/double_stimulus/two_seeded_domains.csv"
        ),
        three_reversal=pd.read_csv(
            three / "testing/wgd_null/polarity_reversal/polarity_reversal.csv"
        ),
        moving=moving,
        directional=pd.read_csv(RESULTS / "focused_motile_length_directional_atlas.csv"),
        one_roots=pd.read_csv(one / "calibration/local_bistability/roots.csv"),
        one_regime=pd.read_csv(one / "calibration/regime_map/regime_map.csv"),
        one_dual_calibration=pd.read_csv(one / "testing/dual_pole_map/dual_pole_map.csv"),
        one_reversal_calibration=pd.read_csv(
            one / "testing/polarity_reversal/polarity_reversal_map.csv"
        ),
        three_regime=pd.read_csv(
            three / "calibration/tier1_regime_map/tier1_regime_map.csv"
        ),
    )


def clean_axis(ax: plt.Axes, grid: bool = False) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="y", color=GRID, linewidth=0.65, zorder=0)
        ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str, x: float = -0.10, y: float = 1.08) -> None:
    if ax.figure.get_constrained_layout() and y == 1.08:
        y = 0.94
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        color=INK,
        va="top",
        ha="left",
        clip_on=False,
    )


def figure_title(fig: plt.Figure, number: str, title: str) -> None:
    fig.suptitle(
        f"{number}  |  {title}",
        x=0.015,
        y=1.035 if fig.get_constrained_layout() else 0.995,
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
        color=INK,
    )


def conclusion_band(fig: plt.Figure, text: str, y: float = 0.006) -> None:
    fig.text(
        0.5,
        y,
        text,
        ha="center",
        va="bottom",
        fontsize=10,
        color=INK,
        fontweight="semibold",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#EEF4F6", "edgecolor": "none"},
    )


def save_figure(fig: plt.Figure, filename: str, title: str) -> Path:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT / filename
    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.12,
        pil_kwargs={"compress_level": 6},
        metadata={"Title": title, "Author": "migration1d figure pipeline"},
    )
    plt.close(fig)
    return path


def phenotype_from_status(value: str) -> str:
    value = str(value)
    if value == "single_polarized":
        return "single"
    if value in {"dual_polarized", "multi_peak"}:
        return "dual"
    if value.startswith("homogeneous"):
        return "homogeneous"
    return "other"


def outcome_from_count(value: float) -> str:
    if not np.isfinite(value):
        return "other"
    count = int(value)
    return "homogeneous" if count == 0 else "single" if count == 1 else "dual" if count == 2 else "other"


def endpoint_state_table(df: pd.DataFrame, state_col: str, keys: Sequence[str]) -> pd.DataFrame:
    pivot = df.pivot(index=list(keys), columns="volume_multiplier", values=state_col)
    return pd.crosstab(pivot[1.0], pivot[2.0], dropna=False)


def draw_transition_matrix(ax: plt.Axes, table: pd.DataFrame, highlight: str) -> None:
    all_states = sorted(set(table.index).union(table.columns))
    matrix = table.reindex(index=all_states, columns=all_states, fill_value=0).to_numpy()
    image = ax.imshow(matrix, cmap="Blues", norm=mpl_colors.PowerNorm(0.45), aspect="equal")
    threshold = max(1.0, matrix.max() * 0.38)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = int(matrix[row, col])
            ax.text(
                col,
                row,
                str(value),
                ha="center",
                va="center",
                fontsize=7.2,
                color="white" if value > threshold else INK,
                fontweight="semibold" if row == col else "normal",
            )
    labels = [state.replace("_", "\n") for state in all_states]
    ax.set_xticks(range(len(all_states)), labels, rotation=40, ha="right")
    ax.set_yticks(range(len(all_states)), labels)
    ax.set(xlabel="2× endpoint", ylabel="1× endpoint", title="Matched endpoint states")
    ax.text(
        0.5,
        -0.34,
        highlight,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8.2,
        color=INK,
        fontweight="bold",
    )
    for spine in ax.spines.values():
        spine.set_visible(False)


def grid_from_frame(
    df: pd.DataFrame,
    x: str,
    y: str,
    value: str,
    x_values: Sequence[float] | None = None,
    y_values: Sequence[float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs = np.asarray(sorted(df[x].unique()) if x_values is None else x_values, dtype=float)
    ys = np.asarray(sorted(df[y].unique()) if y_values is None else y_values, dtype=float)
    lookup = {(float(row[x]), float(row[y])): row[value] for _, row in df.iterrows()}
    grid = np.asarray([[lookup.get((float(xv), float(yv)), np.nan) for xv in xs] for yv in ys])
    return xs, ys, grid


def draw_categorical_map(
    ax: plt.Axes,
    df: pd.DataFrame,
    x: str,
    y: str,
    category: str,
    order: Sequence[str],
    colors: dict[str, str],
    title: str,
    xlabel: str,
    ylabel: str,
    marker: tuple[float, float] | None = None,
) -> None:
    work = df.copy()
    work["_code"] = work[category].map({item: index for index, item in enumerate(order)})
    xs, ys, grid = grid_from_frame(work, x, y, "_code")
    cmap = mpl_colors.ListedColormap([colors[item] for item in order])
    ax.imshow(
        grid,
        origin="lower",
        interpolation="nearest",
        aspect="auto",
        cmap=cmap,
        vmin=-0.5,
        vmax=len(order) - 0.5,
        extent=(xs.min(), xs.max(), ys.min(), ys.max()),
    )
    if marker is not None:
        ax.scatter(
            [marker[0]],
            [marker[1]],
            s=68,
            marker="*",
            facecolor="white",
            edgecolor=INK,
            linewidth=0.8,
            zorder=5,
        )
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(length=2.5)
    unique = [item for item in work[category].dropna().unique()]
    if len(unique) == 1:
        label = (labels := {"homogeneous": "homogeneous throughout", "not_initially_left_polarized": "not eligible throughout"}).get(str(unique[0]), f"all {str(unique[0]).replace('_', ' ')}")
        ax.text(
            0.5,
            0.5,
            label,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=7.5,
            color=INK,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": (1, 1, 1, 0.86), "edgecolor": GRID},
        )


def legend_handles(order: Sequence[str], colors: dict[str, str], labels: dict[str, str] | None = None):
    return [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markersize=7,
            color=colors[item],
            label=(labels or {}).get(item, item.replace("_", " ")),
        )
        for item in order
    ]


def draw_matched_strip(
    ax: plt.Axes,
    df: pd.DataFrame,
    keys: Sequence[str],
    title: str,
) -> None:
    pivot = df.pivot(index=list(keys), columns="volume_multiplier", values="reversal_outcome")
    volumes = [0.5, 1.0, 2.0]
    eligible = pivot.apply(
        lambda row: any(value != "not_initially_left_polarized" for value in row), axis=1
    )
    pivot = pivot.loc[eligible, volumes]
    rank = {state: index for index, state in enumerate(REVERSAL_ORDER)}
    order = sorted(
        range(len(pivot)),
        key=lambda index: tuple(rank.get(value, 99) for value in pivot.iloc[index].values),
    )
    codes = np.asarray(
        [[rank.get(value, rank["ambiguous"]) for value in row] for row in pivot.iloc[order].values]
    )
    cmap = mpl_colors.ListedColormap([REVERSAL_COLORS[state] for state in REVERSAL_ORDER])
    ax.imshow(
        codes,
        aspect="auto",
        interpolation="nearest",
        origin="upper",
        cmap=cmap,
        vmin=-0.5,
        vmax=len(REVERSAL_ORDER) - 0.5,
    )
    ax.set_xticks(range(3), ["0.5×", "1×", "2×"])
    ax.set_yticks([])
    ax.set(xlabel="matched volume", title=title)
    ax.text(
        -0.03,
        0.5,
        f"{len(pivot)} conditions\neligible at ≥1 size",
        transform=ax.transAxes,
        rotation=90,
        ha="right",
        va="center",
        color=MUTED,
        fontsize=7.5,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)


def reversal_prevalence(df: pd.DataFrame) -> pd.DataFrame:
    eligible = df[df.reversal_outcome != "not_initially_left_polarized"].copy()
    rows = []
    for volume, subset in eligible.groupby("volume_multiplier"):
        rows.append(
            {
                "volume": float(volume),
                "n": len(subset),
                "clean": 100.0 * subset.reversal_outcome.eq("reversed").mean(),
                "multiple": 100.0
                * subset.reversal_outcome.eq("multiple_polarized_domains").mean(),
            }
        )
    return pd.DataFrame(rows).sort_values("volume")


def competence_prevalence(df: pd.DataFrame, state_col: str) -> pd.Series:
    # The writer's definition counts all organized/dynamic non-homogeneous endpoints.
    organized = ~df[state_col].astype(str).str.startswith("homogeneous")
    return organized.groupby(df.volume_multiplier).mean() * 100.0


def dual_prevalence(df: pd.DataFrame) -> pd.Series:
    return (df.polarized_domain_count_t300 >= 2).groupby(df.volume_multiplier).mean() * 100.0


def fixed_key_columns(species: str) -> list[str]:
    return ["baseline_h", "baseline_mu", "baseline_Lambda"] if species == "single_species" else ["n", "kappa_R", "baseline_Lambda"]


def task_index(directory: Path) -> dict[tuple, Path]:
    result: dict[tuple, Path] = {}
    for path in sorted((directory / "tasks").glob("task_*.json")):
        payload = json.loads(path.read_text())
        key = tuple(payload["task_key"])
        result[key] = directory / "tasks" / payload["trajectory"]
    return result


_TASK_INDICES: dict[tuple[str, str], dict[tuple, Path]] = {}


def moving_trajectory(
    species: str,
    mode: str,
    length: float,
    mechanical_case: str,
    protocol: str,
    alpha: float,
) -> dict[str, np.ndarray]:
    key = (species, mode)
    if key not in _TASK_INDICES:
        suffix = "focused_motile_length" if mode == "global" else "focused_motile_length_spatial_tension"
        directory = RESULTS / species / "testing" / suffix
        _TASK_INDICES[key] = task_index(directory)
    task_key = (species, float(length), mechanical_case, protocol, float(alpha))
    path = _TASK_INDICES[key][task_key]
    with np.load(path) as loaded:
        return {name: loaded[name].copy() for name in loaded.files}


def fixed_cache_path(species: str, kind: str, volume: float, nodes: int) -> Path:
    return CACHE / f"{species}_{kind}_v{volume:g}_n{nodes}.npz"


def simulate_fixed(species: str, kind: str, volume: float, nodes: int = 81) -> dict[str, np.ndarray]:
    """Run or load one selected fixed-domain representative trajectory."""
    path = fixed_cache_path(species, kind, volume, nodes)
    if path.exists():
        with np.load(path) as loaded:
            return {name: loaded[name].copy() for name in loaded.files}

    CACHE.mkdir(parents=True, exist_ok=True)
    if species == "single_species":
        from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
        from migration1d.single_species.model import (
            ActivationPulse,
            DualPulseInitialCondition,
            PulseTrain,
            StripWavePinningChemistry,
        )
        from migration1d.single_species.wgd import parser as wgd_parser, scaled_parameters

        args = wgd_parser().parse_args(["--mode", "shard"])
        args.nodes = nodes
        representatives = {
            "single": (1.0, 2.46875, 49.375),
            "dual": (1.0, 2.375, 25.0),
            "reversal": (1.0, 2.1875, 71.875),
        }
        h, mu, baseline_lambda = representatives[kind]
        parameters, _, _, _ = scaled_parameters(args, mu, baseline_lambda, h, volume)
        if kind == "single":
            chemistry = StripWavePinningChemistry(
                parameters,
                activation=ActivationPulse(1.0, 0.1, 0.05, 10.0, 30.0),
            )
        elif kind == "dual":
            chemistry = StripWavePinningChemistry(
                parameters,
                initial_condition=DualPulseInitialCondition(
                    amplitude=2.0, width_fraction=0.04, asymmetry=0.01
                ),
            )
        else:
            chemistry = StripWavePinningChemistry(
                parameters,
                activation=PulseTrain(
                    (
                        ActivationPulse(1.0, 0.1, 0.05, 10.0, 30.0),
                        ActivationPulse(1.0, 0.9, 0.05, 150.0, 170.0),
                    )
                ),
            )
        result = run_simulation(
            MigrationModel(
                mechanics=MechanicsParameters(reference_length=parameters.domain_length),
                chemistry=chemistry,
            ),
            SimulationConfig(
                nodes=nodes,
                duration=300.0,
                output_interval=2.0,
                max_step=0.5,
                noise_amplitude=0.0,
                mechanics_activation_time=None,
            ),
        )
    else:
        from migration1d.holmes import HolmesChemistry
        from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
        from migration1d.three_species.model import GaussianRacPulse, RacPulseTrain
        from migration1d.three_species.wgd_extended import (
            SeededRacHolmesChemistry,
            parameters as scaled_holmes_parameters,
        )

        representatives = {
            "single": (3.0, 0.05333333333333334, 68.51601597035),
            "dual": (3.0, 0.08, 84.3274042712),
            "reversal": (4.0, 0.04, 63.2455532034),
        }
        hill, kappa, baseline_lambda = representatives[kind]
        parameters, _, _, _, length = scaled_holmes_parameters(
            kappa, hill, baseline_lambda, volume
        )
        if kind == "single":
            chemistry = HolmesChemistry(
                parameters,
                cue=GaussianRacPulse(1.0, 0.1, 0.05, 10.0, 30.0),
            )
        elif kind == "dual":
            chemistry = SeededRacHolmesChemistry(parameters, 2.0, 0.04, 0.01)
        else:
            chemistry = HolmesChemistry(
                parameters,
                cue=RacPulseTrain(
                    (
                        GaussianRacPulse(1.0, 0.1, 0.05, 10.0, 30.0),
                        GaussianRacPulse(1.0, 0.9, 0.05, 150.0, 170.0),
                    )
                ),
            )
        result = run_simulation(
            MigrationModel(MechanicsParameters(reference_length=length), chemistry),
            SimulationConfig(
                nodes=nodes,
                duration=300.0,
                output_interval=2.0,
                max_step=0.5,
                noise_amplitude=0.0,
                mechanics_activation_time=None,
                jacobian_sparsity=True,
            ),
        )

    arrays = {
        "times": result.times,
        "position": result.position,
        "membrane": result.membrane,
        "cytoplasm": result.cytoplasm,
    }
    np.savez_compressed(path, **arrays)
    return arrays


def material_fraction(position: np.ndarray) -> np.ndarray:
    return (position - position[:, [0]]) / (position[:, [-1]] - position[:, [0]])


def add_fixed_cues(ax: plt.Axes, kind: str) -> None:
    if kind == "single":
        ax.add_patch(
            patches.Rectangle((10, 0.035), 20, 0.13, facecolor=SKY, edgecolor="white", linewidth=0.5, alpha=0.9)
        )
    elif kind == "reversal":
        ax.add_patch(
            patches.Rectangle((10, 0.035), 20, 0.13, facecolor=SKY, edgecolor="white", linewidth=0.5, alpha=0.9)
        )
        ax.add_patch(
            patches.Rectangle((150, 0.835), 20, 0.13, facecolor=ORANGE, edgecolor="white", linewidth=0.5, alpha=0.9)
        )


def draw_fixed_kymograph(
    ax: plt.Axes,
    data: dict[str, np.ndarray],
    species_index: int,
    title: str,
    kind: str,
    vmax: float | None = None,
    ylabel: bool = True,
) -> matplotlib.collections.QuadMesh:
    times = data["times"]
    fraction = material_fraction(data["position"])
    profile = data["membrane"][:, species_index]
    time_grid = np.broadcast_to(times[:, None], fraction.shape)
    image = ax.pcolormesh(
        time_grid,
        fraction,
        profile,
        shading="auto",
        cmap="magma",
        vmin=0.0,
        vmax=vmax,
        rasterized=True,
    )
    add_fixed_cues(ax, kind)
    ax.set(title=title, xlabel="time (s)", ylabel="cell position" if ylabel else "")
    ax.set_yticks([0, 0.5, 1], ["left", "mid", "right"] if ylabel else [])
    ax.set_xlim(0, 300)
    return image


def draw_moving_kymograph(
    ax: plt.Axes,
    data: dict[str, np.ndarray],
    species_index: int,
    title: str,
    vmax: float | None = None,
    ylabel: bool = True,
) -> matplotlib.image.AxesImage:
    times = data["times"]
    position = data["position"]
    profile = data["membrane"][:, species_index]
    laboratory_grid = np.linspace(float(position.min()), float(position.max()), 360)
    raster = np.full((len(times), len(laboratory_grid)), np.nan)
    for index in range(len(times)):
        inside = (laboratory_grid >= position[index, 0]) & (laboratory_grid <= position[index, -1])
        raster[index, inside] = np.interp(
            laboratory_grid[inside], position[index], profile[index]
        )
    image = ax.imshow(
        raster.T,
        origin="lower",
        interpolation="bilinear",
        aspect="auto",
        extent=(times[0], times[-1], laboratory_grid[0], laboratory_grid[-1]),
        cmap="magma",
        vmin=0.0,
        vmax=vmax,
    )
    ax.plot(times, position[:, 0], color="white", linewidth=0.7)
    ax.plot(times, position[:, -1], color="white", linewidth=0.7)
    ax.set(title=title, xlabel="time (s)", ylabel="laboratory position (µm)" if ylabel else "")
    return image


def draw_cell_scaling(ax: plt.Axes, network: bool = False) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if network:
        nodes = {"Cdc42": (0.21, 0.68, SKY), "Rac": (0.50, 0.68, MAGENTA), "Rho": (0.79, 0.68, ORANGE)}
        for label, (x, y, color) in nodes.items():
            ax.add_patch(patches.Circle((x, y), 0.09, facecolor=color, edgecolor="white", linewidth=1.4))
            ax.text(x, y, label, ha="center", va="center", fontsize=8, fontweight="bold", color="white" if label != "Rho" else INK)
        ax.annotate("", (0.42, 0.68), (0.30, 0.68), arrowprops={"arrowstyle": "->", "color": INK, "lw": 1.5})
        ax.annotate("", (0.70, 0.68), (0.58, 0.68), arrowprops={"arrowstyle": "-[", "color": INK, "lw": 1.5})
        ax.annotate("", (0.27, 0.61), (0.71, 0.61), arrowprops={"arrowstyle": "-[", "color": INK, "lw": 1.5, "connectionstyle": "arc3,rad=-0.35"})
        ax.text(0.5, 0.49, "active ↔ inactive pools  •  Rac is the common readout", ha="center", va="center", fontsize=7.5, color=MUTED)
    else:
        ax.add_patch(patches.FancyBboxPatch((0.08, 0.59), 0.34, 0.18, boxstyle="round,pad=0.03,rounding_size=0.08", facecolor="#E7F2F8", edgecolor=BLUE, linewidth=1.2))
        ax.text(0.25, 0.68, "active u\n(membrane)", ha="center", va="center", fontsize=8, color=INK)
        ax.add_patch(patches.FancyBboxPatch((0.58, 0.59), 0.34, 0.18, boxstyle="round,pad=0.03,rounding_size=0.08", facecolor="#F1EDF6", edgecolor=PURPLE, linewidth=1.2))
        ax.text(0.75, 0.68, "inactive v\n(fast reservoir)", ha="center", va="center", fontsize=8, color=INK)
        ax.annotate("activation", (0.56, 0.71), (0.44, 0.71), ha="center", fontsize=7, arrowprops={"arrowstyle": "->", "color": GREEN, "lw": 1.4})
        ax.annotate("inactivation", (0.44, 0.62), (0.56, 0.62), ha="center", fontsize=7, arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.4})
        ax.text(0.5, 0.49, "mass conserved  •  localized activation cue", ha="center", va="center", fontsize=7.5, color=MUTED)
    lengths = [0.18, 0.23, 0.29]
    centers = [0.18, 0.50, 0.82]
    labels = ["0.5× V", "1× V", "2× V"]
    for length, center, label, color in zip(lengths, centers, labels, VOLUME_COLORS.values()):
        x0 = center - length / 2
        ax.add_patch(patches.FancyBboxPatch((x0, 0.22), length, 0.10, boxstyle="round,pad=0.015,rounding_size=0.05", facecolor=color, edgecolor="none", alpha=0.9))
        ax.text(center, 0.14, label, ha="center", va="center", fontsize=8, fontweight="semibold")
    ax.annotate("geometry only", (0.70, 0.27), (0.31, 0.27), ha="center", va="center", fontsize=7.5, color=MUTED, arrowprops={"arrowstyle": "->", "color": MUTED, "lw": 1.0})
    ax.text(0.5, 0.035, "concentrations + kinetics fixed  •  L ∝ V¹ᐟ³", ha="center", va="bottom", fontsize=7.7, color=INK)


def endpoint_means(data: dict[str, np.ndarray], species_index: int) -> tuple[np.ndarray, np.ndarray]:
    profile = data["membrane"][:, species_index]
    width = max(3, int(round(profile.shape[1] * 0.1)))
    return profile[:, :width].mean(axis=1), profile[:, -width:].mean(axis=1)


def domain_counts(data: dict[str, np.ndarray], species_index: int) -> np.ndarray:
    from migration1d.single_species.observables import polarized_domain_count

    return np.asarray([polarized_domain_count(profile) for profile in data["membrane"][:, species_index]])


def figure_1(data: DataBundle) -> Path:
    """Minimal wave-pinning size response."""
    representatives = {
        kind: {volume: simulate_fixed("single_species", kind, volume) for volume in (0.5, 1.0, 2.0)}
        for kind in ("single", "dual", "reversal")
    }
    fig, axes = plt.subplot_mosaic(
        [
            ["A", "A", "B05", "B1", "B2", "BM"],
            ["C05", "C1", "C2", "CT", "D1", "D2"],
            ["E", "E", "E", "F", "F", "F"],
        ],
        figsize=(19.0, 12.2),
        gridspec_kw={"height_ratios": [1.0, 1.0, 0.95], "hspace": 0.48, "wspace": 0.38},
    )
    figure_title(
        fig,
        "Figure 1",
        "Cell enlargement preserves polarization but impairs polarity replacement",
    )

    draw_cell_scaling(axes["A"], network=False)
    axes["A"].set_title("Minimal conserved polarity circuit", pad=4)
    panel_label(axes["A"], "A", x=-0.04)

    single_vmax = max(item["membrane"][:, 0].max() for item in representatives["single"].values())
    image = None
    for index, volume in enumerate((0.5, 1.0, 2.0)):
        ax = axes[{0.5: "B05", 1.0: "B1", 2.0: "B2"}[volume]]
        image = draw_fixed_kymograph(
            ax,
            representatives["single"][volume],
            0,
            f"{volume:g}× volume",
            "single",
            vmax=single_vmax,
            ylabel=index == 0,
        )
        if index == 0:
            panel_label(ax, "B")
    cbar = fig.colorbar(image, ax=[axes["B05"], axes["B1"], axes["B2"]], location="bottom", pad=0.15, shrink=0.86)
    cbar.set_label("active u (common scale)")
    state_table = endpoint_state_table(
        data.one_single,
        "terminal_outcome",
        ["baseline_h", "baseline_mu", "baseline_Lambda"],
    )
    unchanged = int(np.trace(state_table.reindex(index=sorted(set(state_table.index).union(state_table.columns)), columns=sorted(set(state_table.index).union(state_table.columns)), fill_value=0).to_numpy()))
    draw_transition_matrix(
        axes["BM"],
        state_table,
        f"{unchanged}/867 unchanged ({100 * unchanged / 867:.1f}%)",
    )

    dual_subset = data.one_dual[
        data.one_dual.baseline_h.eq(1.0) & data.one_dual.epsilon.eq(0.01)
    ].copy()
    dual_subset["phenotype"] = dual_subset.polarized_domain_count_t300.map(outcome_from_count)
    for index, volume in enumerate((0.5, 1.0, 2.0)):
        ax = axes[{0.5: "C05", 1.0: "C1", 2.0: "C2"}[volume]]
        draw_categorical_map(
            ax,
            dual_subset[dual_subset.volume_multiplier.eq(volume)],
            "baseline_Lambda",
            "baseline_mu",
            "phenotype",
            PHENOTYPE_ORDER,
            PHENOTYPE_COLORS,
            f"{volume:g}× volume",
            "baseline Λ",
            "baseline μ" if index == 0 else "",
            marker=(25.0, 2.375),
        )
        if index == 0:
            panel_label(ax, "C")
    axes["C05"].legend(
        handles=legend_handles(PHENOTYPE_ORDER, PHENOTYPE_COLORS),
        loc="lower left",
        ncol=2,
        fontsize=7.0,
    )
    prevalence = dual_prevalence(data.one_dual)
    ax = axes["CT"]
    ax.plot(prevalence.index, prevalence.values, color=MAGENTA, marker="o", linewidth=2.0)
    for volume, value in prevalence.items():
        ax.annotate(f"{value:.2f}%", (volume, value), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8, fontweight="semibold")
    ax.set(
        title="Persistent two-domain incidence",
        xlabel="volume multiplier",
        ylabel="all matched conditions (%)",
        xticks=[0.5, 1, 2],
    )
    ax.set_ylim(0, max(5.2, prevalence.max() * 1.3))
    clean_axis(ax, grid=True)

    reversal_vmax = max(
        representatives["reversal"][volume]["membrane"][:, 0].max() for volume in (1.0, 2.0)
    )
    for index, volume in enumerate((1.0, 2.0)):
        ax = axes[{1.0: "D1", 2.0: "D2"}[volume]]
        im = draw_fixed_kymograph(
            ax,
            representatives["reversal"][volume],
            0,
            "1×: clean reversal" if volume == 1.0 else "2×: competing fronts",
            "reversal",
            vmax=reversal_vmax,
            ylabel=index == 0,
        )
        if index == 0:
            panel_label(ax, "D")
    fig.colorbar(im, ax=[axes["D1"], axes["D2"]], location="bottom", pad=0.15, shrink=0.83, label="active u")

    draw_matched_strip(
        axes["E"],
        data.one_reversal,
        ["baseline_h", "baseline_mu", "baseline_Lambda"],
        "Matched reversal outcomes",
    )
    panel_label(axes["E"], "E", x=-0.03)
    prev = reversal_prevalence(data.one_reversal)
    axes["E"].text(
        0.985,
        0.96,
        "\n".join(
            [
                f"clean reversal  {row.volume:g}×: {row.clean:.1f}%" for row in prev.itertuples()
            ]
            + [
                f"multiple domains  1×→2×: {prev.loc[prev.volume.eq(1), 'multiple'].iloc[0]:.1f}% → {prev.loc[prev.volume.eq(2), 'multiple'].iloc[0]:.1f}%"
            ]
        ),
        transform=axes["E"].transAxes,
        ha="right",
        va="top",
        fontsize=8.2,
        linespacing=1.55,
        color=INK,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": (1, 1, 1, 0.88), "edgecolor": GRID},
    )
    axes["E"].legend(
        handles=legend_handles(REVERSAL_ORDER[:5], REVERSAL_COLORS, REVERSAL_LABELS),
        loc="lower center",
        bbox_to_anchor=(0.50, -0.29),
        ncol=5,
    )

    ax = axes["F"]
    for volume, linestyle in ((1.0, "-"), (2.0, "--")):
        trajectory = representatives["reversal"][volume]
        left, right = endpoint_means(trajectory, 0)
        ax.plot(trajectory["times"], left, color=BLUE, linestyle=linestyle, linewidth=1.8, label=f"old left front · {volume:g}×")
        ax.plot(trajectory["times"], right, color=VERMILLION, linestyle=linestyle, linewidth=1.8, label=f"new right front · {volume:g}×")
    ax.axvspan(150, 170, color=ORANGE, alpha=0.13, linewidth=0)
    ax.set(title="Competitive elimination after the opposing cue", xlabel="time (s)", ylabel="mean endpoint active u")
    clean_axis(ax, grid=True)
    ax2 = ax.twinx()
    for volume, linestyle in ((1.0, "-"), (2.0, "--")):
        trajectory = representatives["reversal"][volume]
        ax2.step(trajectory["times"], domain_counts(trajectory, 0), color=PURPLE, linestyle=linestyle, linewidth=1.0, alpha=0.75)
    ax2.set_ylabel("detected domains", color=PURPLE)
    ax2.set_yticks([0, 1, 2])
    ax.legend(loc="upper left", ncol=2, fontsize=7.4)
    panel_label(ax, "F", x=-0.04)

    conclusion_band(
        fig,
        "Front generation is retained after enlargement, whereas winner-take-all front replacement is impaired.",
    )
    return save_figure(
        fig,
        "figure_1_minimal_polarity.png",
        "Figure 1. Minimal polarity size response",
    )


def add_three_species_profile_inset(ax: plt.Axes, trajectory: dict[str, np.ndarray]) -> None:
    inset = ax.inset_axes([0.55, 0.56, 0.41, 0.38])
    fraction = material_fraction(trajectory["position"])[-1]
    terminal = trajectory["membrane"][-1]
    for profile, color, label in zip(terminal, (SKY, MAGENTA, ORANGE), ("Cdc42", "Rac", "Rho")):
        scale = max(float(np.max(profile)), 1.0e-12)
        inset.plot(fraction, profile / scale, color=color, linewidth=1.0, label=label)
    inset.set(xticks=[], yticks=[], xlim=(0, 1), ylim=(0, 1.05))
    inset.set_facecolor((1, 1, 1, 0.88))
    for spine in inset.spines.values():
        spine.set_color("white")


def figure_2(data: DataBundle) -> Path:
    """Three-species Rho-family size response."""
    representatives = {
        kind: {volume: simulate_fixed("three_species", kind, volume) for volume in (0.5, 1.0, 2.0)}
        for kind in ("single", "dual", "reversal")
    }
    fig, axes = plt.subplot_mosaic(
        [
            ["A", "A", "B05", "B1", "B2", "BM"],
            ["C05", "C1", "C2", "CT", "D1", "D2"],
            ["E1", "E2", "F", "F", "F", "F"],
        ],
        figsize=(19.0, 12.4),
        gridspec_kw={"height_ratios": [1.0, 1.0, 0.95], "hspace": 0.50, "wspace": 0.38},
    )
    figure_title(
        fig,
        "Figure 2",
        "A Cdc42–Rac–Rho network reveals stronger size-dependent multipolarity",
    )
    draw_cell_scaling(axes["A"], network=True)
    axes["A"].set_title("Rho-family polarity network", pad=4)
    panel_label(axes["A"], "A", x=-0.04)

    vmax = max(item["membrane"][:, 1].max() for item in representatives["single"].values())
    image = None
    for index, volume in enumerate((0.5, 1.0, 2.0)):
        ax = axes[{0.5: "B05", 1.0: "B1", 2.0: "B2"}[volume]]
        image = draw_fixed_kymograph(
            ax,
            representatives["single"][volume],
            1,
            f"{volume:g}× volume",
            "single",
            vmax=vmax,
            ylabel=index == 0,
        )
        add_three_species_profile_inset(ax, representatives["single"][volume])
        if index == 0:
            panel_label(ax, "B")
    fig.colorbar(image, ax=[axes["B05"], axes["B1"], axes["B2"]], location="bottom", pad=0.15, shrink=0.86, label="active Rac (common scale)")
    table = endpoint_state_table(
        data.three_single,
        "terminal_polarity_status",
        ["n", "kappa_R", "baseline_Lambda"],
    )
    states = sorted(set(table.index).union(table.columns))
    unchanged = int(np.trace(table.reindex(index=states, columns=states, fill_value=0).to_numpy()))
    draw_transition_matrix(axes["BM"], table, f"{unchanged}/507 unchanged ({100 * unchanged / 507:.1f}%)")

    dual_subset = data.three_dual[
        data.three_dual.n.eq(3.0) & data.three_dual.epsilon.eq(0.01)
    ].copy()
    dual_subset["phenotype"] = dual_subset.polarized_domain_count_t300.map(outcome_from_count)
    for index, volume in enumerate((0.5, 1.0, 2.0)):
        ax = axes[{0.5: "C05", 1.0: "C1", 2.0: "C2"}[volume]]
        draw_categorical_map(
            ax,
            dual_subset[dual_subset.volume_multiplier.eq(volume)],
            "baseline_Lambda",
            "kappa_R",
            "phenotype",
            PHENOTYPE_ORDER,
            PHENOTYPE_COLORS,
            f"{volume:g}× volume",
            "baseline Λ",
            "κR" if index == 0 else "",
            marker=(84.3274042712, 0.08),
        )
        if index == 0:
            panel_label(ax, "C")
    axes["C05"].legend(
        handles=legend_handles(PHENOTYPE_ORDER, PHENOTYPE_COLORS),
        loc="lower left",
        ncol=2,
        fontsize=7.0,
    )
    prevalence = dual_prevalence(data.three_dual)
    ax = axes["CT"]
    ax.plot(prevalence.index, prevalence.values, color=MAGENTA, marker="o", linewidth=2.0)
    for volume, value in prevalence.items():
        ax.annotate(f"{value:.1f}%", (volume, value), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8, fontweight="semibold")
    ax.set(title="Persistent two-domain incidence", xlabel="volume multiplier", ylabel="all matched conditions (%)", xticks=[0.5, 1, 2])
    ax.set_ylim(0, prevalence.max() * 1.27)
    clean_axis(ax, grid=True)

    large_dual = representatives["dual"][2.0]
    draw_fixed_kymograph(axes["D1"], large_dual, 1, "Large-cell active Rac", "dual", ylabel=True)
    panel_label(axes["D1"], "D")
    ax = axes["D2"]
    for volume, linestyle in ((0.5, "-"), (2.0, "--")):
        trajectory = representatives["dual"][volume]
        fraction = material_fraction(trajectory["position"])[-1]
        terminal = trajectory["membrane"][-1]
        for species_index, (color, label) in enumerate(zip((SKY, MAGENTA, ORANGE), ("Cdc42", "Rac", "Rho"))):
            ax.plot(fraction, terminal[species_index], color=color, linestyle=linestyle, linewidth=1.5, label=f"{label} · {volume:g}×")
    ax.set(title="Coherent late signaling profiles", xlabel="cell position", ylabel="active concentration")
    ax.set_xticks([0, 0.5, 1], ["left", "mid", "right"])
    clean_axis(ax, grid=True)
    ax.legend(ncol=2, fontsize=7.0)

    reversal_vmax = max(representatives["reversal"][volume]["membrane"][:, 1].max() for volume in (1.0, 2.0))
    for index, volume in enumerate((1.0, 2.0)):
        ax = axes[{1.0: "E1", 2.0: "E2"}[volume]]
        image = draw_fixed_kymograph(
            ax,
            representatives["reversal"][volume],
            1,
            "1×: Rac front replaced" if volume == 1.0 else "2×: incomplete replacement",
            "reversal",
            vmax=reversal_vmax,
            ylabel=index == 0,
        )
        add_three_species_profile_inset(ax, representatives["reversal"][volume])
        if index == 0:
            panel_label(ax, "E")
    fig.colorbar(image, ax=[axes["E1"], axes["E2"]], location="bottom", pad=0.15, shrink=0.83, label="active Rac")

    draw_matched_strip(
        axes["F"],
        data.three_reversal,
        ["n", "kappa_R", "baseline_Lambda"],
        "Population reversal outcomes",
    )
    panel_label(axes["F"], "F", x=-0.03)
    prev = reversal_prevalence(data.three_reversal)
    annotation = "\n".join(
        [f"clean reversal  {row.volume:g}×: {row.clean:.1f}%" for row in prev.itertuples()]
        + [f"multiple domains  0.5×→2×: {prev.iloc[0].multiple:.1f}% → {prev.iloc[-1].multiple:.1f}%"]
    )
    axes["F"].text(0.985, 0.96, annotation, transform=axes["F"].transAxes, ha="right", va="top", fontsize=8.2, linespacing=1.55, bbox={"boxstyle": "round,pad=0.3", "facecolor": (1, 1, 1, 0.88), "edgecolor": GRID})
    axes["F"].legend(
        handles=legend_handles(REVERSAL_ORDER[:5], REVERSAL_COLORS, REVERSAL_LABELS),
        loc="lower center",
        bbox_to_anchor=(0.50, -0.29),
        ncol=5,
    )

    conclusion_band(
        fig,
        "Network crosstalk preserves the reversal defect and enlarges the parameter region that sustains competing fronts.",
    )
    return save_figure(
        fig,
        "figure_2_rho_network.png",
        "Figure 2. Rho-family network size response",
    )


def draw_model_ladder(ax: plt.Axes, include_feedback: bool = False) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    count = 4 if include_feedback else 3
    centers = np.linspace(0.13, 0.87, count)
    labels = ["1-species\npolarity", "Cdc42–Rac–Rho\npolarity", "polarity +\nmechanics"]
    colors = [BLUE, MAGENTA, ORANGE]
    if include_feedback:
        labels.append("polarity + mechanics\n+ tension feedback")
        colors.append(PURPLE)
    for index, (center, label, color) in enumerate(zip(centers, labels, colors)):
        ax.add_patch(
            patches.FancyBboxPatch(
                (center - 0.105, 0.47),
                0.21,
                0.34,
                boxstyle="round,pad=0.02,rounding_size=0.04",
                facecolor=mpl_colors.to_rgba(color, 0.12),
                edgecolor=color,
                linewidth=1.4,
            )
        )
        ax.text(center, 0.64, label, ha="center", va="center", color=INK, fontsize=9, fontweight="semibold")
        if index < len(centers) - 1:
            ax.annotate("", (centers[index + 1] - 0.12, 0.64), (center + 0.12, 0.64), arrowprops={"arrowstyle": "->", "color": MUTED, "lw": 1.4})
    challenges = ["single cue", "competing cues", "cue reversal", "migration"]
    for x, challenge in zip(np.linspace(0.19, 0.81, 4), challenges):
        ax.text(x, 0.18, challenge, ha="center", va="center", fontsize=8, color=INK, bbox={"boxstyle": "round,pad=0.24", "facecolor": PANEL_BG, "edgecolor": GRID})
    ax.text(0.5, 0.02, "common phenotype challenges", ha="center", va="bottom", fontsize=7.5, color=MUTED)


def raw_fixed_phenotypes(data: DataBundle, species: str) -> dict[str, pd.Series]:
    if species == "single_species":
        single, dual, reversal, state = data.one_single, data.one_dual, data.one_reversal, "terminal_outcome"
    else:
        single, dual, reversal, state = data.three_single, data.three_dual, data.three_reversal, "terminal_polarity_status"
    rev = reversal_prevalence(reversal).set_index("volume")["clean"]
    return {
        "single-cue competence": competence_prevalence(single, state),
        "persistent dual fronts": dual_prevalence(dual),
        "clean reversal": rev,
    }


def figure_3(data: DataBundle) -> Path:
    fig = plt.figure(figsize=(17.5, 12.0))
    figure_title(fig, "Figure 3", "Increasing model complexity preserves a common size-response fingerprint")
    outer = GridSpec(3, 1, figure=fig, height_ratios=[0.65, 1.75, 1.15], hspace=0.43, top=0.94, bottom=0.075)

    ax_a = fig.add_subplot(outer[0])
    draw_model_ladder(ax_a)
    panel_label(ax_a, "A", x=-0.015, y=1.03)

    middle = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], width_ratios=[1.65, 1.0], wspace=0.25)
    bgrid = GridSpecFromSubplotSpec(3, 2, subplot_spec=middle[0], hspace=0.38, wspace=0.28)
    phenotypes = ["single-cue competence", "persistent dual fronts", "clean reversal"]
    raw = {
        "single_species": raw_fixed_phenotypes(data, "single_species"),
        "three_species": raw_fixed_phenotypes(data, "three_species"),
    }
    b_axes = []
    for row, phenotype in enumerate(phenotypes):
        for column, species in enumerate(("single_species", "three_species")):
            ax = fig.add_subplot(bgrid[row, column])
            series = raw[species][phenotype]
            ax.plot(series.index, series.values, marker="o", linewidth=2.0, color=MODEL_COLORS[species])
            for x, value in series.items():
                ax.annotate(f"{value:.1f}", (x, value), xytext=(0, 5), textcoords="offset points", ha="center", fontsize=7)
            if row == 0:
                ax.set_title("minimal 1-species" if species == "single_species" else "Cdc42–Rac–Rho")
            ax.set_ylabel(f"{phenotype}\n(%)" if column == 0 else "")
            ax.set_xlabel("volume multiplier" if row == 2 else "")
            ax.set_xticks([0.5, 1, 2], ["0.5×", "1×", "2×"] if row == 2 else ["", "", ""])
            padding = max(4.0, 0.15 * (series.max() - series.min() + 1))
            ax.set_ylim(max(0, series.min() - padding), min(100, series.max() + padding))
            clean_axis(ax, grid=True)
            b_axes.append(ax)
    panel_label(b_axes[0], "B", x=-0.26, y=1.16)

    ax_c = fig.add_subplot(middle[1])
    effects = []
    for phenotype in phenotypes:
        for species in ("single_species", "three_species"):
            series = raw[species][phenotype]
            effects.append((phenotype, species, float(series.loc[2.0] - series.loc[1.0])))
    ybase = np.arange(len(phenotypes))[::-1]
    offsets = {"single_species": -0.12, "three_species": 0.12}
    for species in ("single_species", "three_species"):
        values = [next(value for phenotype2, species2, value in effects if phenotype2 == phenotype and species2 == species) for phenotype in phenotypes]
        ax_c.scatter(values, ybase + offsets[species], s=62, color=MODEL_COLORS[species], edgecolor="white", linewidth=0.8, zorder=3, label="1-species" if species == "single_species" else "3-species")
        for value, y in zip(values, ybase + offsets[species]):
            ax_c.text(value + (0.8 if value >= 0 else -0.8), y, f"{value:+.1f} pp", ha="left" if value >= 0 else "right", va="center", fontsize=7.6, color=MODEL_COLORS[species])
    ax_c.axvline(0, color=INK, linewidth=1.0)
    ax_c.set_yticks(ybase, ["single-cue\ncompetence", "persistent\ndual fronts", "clean\nreversal"])
    ax_c.set(xlabel="paired prevalence change, 1× → 2× (percentage points)", title="Size-effect fingerprint")
    ax_c.set_xlim(-28, 10)
    clean_axis(ax_c, grid=True)
    ax_c.legend(loc="upper center", ncol=2)
    panel_label(ax_c, "C", x=-0.16)

    bottom = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[2], width_ratios=[1.0, 1.0, 1.15], wspace=0.35)
    ax_d1 = fig.add_subplot(bottom[0])
    draw_matched_strip(ax_d1, data.one_reversal, ["baseline_h", "baseline_mu", "baseline_Lambda"], "Minimal model")
    panel_label(ax_d1, "D", x=-0.12)
    ax_d1.text(0.5, -0.20, "68.3% of 1× reversing conditions\nstill reverse at 2×", transform=ax_d1.transAxes, ha="center", va="top", fontsize=8, fontweight="semibold")
    ax_d2 = fig.add_subplot(bottom[1])
    draw_matched_strip(ax_d2, data.three_reversal, ["n", "kappa_R", "baseline_Lambda"], "Rho-family network")
    ax_d2.text(0.5, -0.20, "63.0% of 1× reversing conditions\nstill reverse at 2×", transform=ax_d2.transAxes, ha="center", va="top", fontsize=8, fontweight="semibold")
    fig.legend(
        handles=legend_handles(REVERSAL_ORDER[:5], REVERSAL_COLORS, REVERSAL_LABELS),
        loc="lower center",
        bbox_to_anchor=(0.38, 0.045),
        ncol=5,
        fontsize=7.3,
    )

    ax_e = fig.add_subplot(bottom[2])
    ax_e.set_axis_off()
    panel_label(ax_e, "E", x=-0.08)
    ax_e.set_xlim(0, 1)
    ax_e.set_ylim(0, 1)
    ax_e.text(0.5, 0.94, "Mechanistic synthesis", ha="center", va="top", fontsize=10, fontweight="semibold")
    for y, label, color in ((0.72, "front generation", BLUE), (0.47, "front selection", MAGENTA), (0.22, "front replacement", VERMILLION)):
        ax_e.add_patch(patches.FancyBboxPatch((0.08, y - 0.08), 0.35, 0.16, boxstyle="round,pad=0.02", facecolor=mpl_colors.to_rgba(color, 0.13), edgecolor=color, linewidth=1.2))
        ax_e.text(0.255, y, label, ha="center", va="center", fontsize=8.5, fontweight="semibold")
    ax_e.text(0.58, 0.72, "largely retained", ha="left", va="center", fontsize=8.5, color=BLUE)
    ax_e.text(0.58, 0.47, "size-sensitive", ha="left", va="center", fontsize=8.5, color=MAGENTA)
    ax_e.text(0.58, 0.22, "consistently impaired", ha="left", va="center", fontsize=8.5, color=VERMILLION)
    ax_e.text(0.52, 0.84, "larger domain separates\ncandidate fronts", ha="center", va="center", fontsize=7.5, color=MUTED)
    ax_e.annotate("", (0.52, 0.53), (0.52, 0.77), arrowprops={"arrowstyle": "->", "color": MUTED})

    conclusion_band(fig, "Cell enlargement reduces polarity decisiveness rather than polarity competence.")
    return save_figure(fig, "figure_3_size_response_fingerprint.png", "Figure 3. Cross-model size-response fingerprint")


def draw_mechanochemical_schematic(ax: plt.Axes, feedback: bool = False) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.add_patch(
        patches.FancyBboxPatch(
            (0.08, 0.55),
            0.84,
            0.18,
            boxstyle="round,pad=0.025,rounding_size=0.09",
            facecolor="#EDF3F6",
            edgecolor=INK,
            linewidth=1.2,
        )
    )
    ax.add_patch(patches.Circle((0.78, 0.64), 0.065, facecolor=MAGENTA, edgecolor="white"))
    ax.text(0.78, 0.64, "front", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
    ax.annotate("active force", (0.96, 0.64), (0.84, 0.64), fontsize=7.5, ha="left", va="center", arrowprops={"arrowstyle": "->", "color": VERMILLION, "lw": 1.7})
    ax.annotate("rear / inward force", (0.04, 0.64), (0.15, 0.64), fontsize=7.2, ha="right", va="center", arrowprops={"arrowstyle": "<-", "color": BLUE, "lw": 1.5})
    ax.annotate("", (0.63, 0.39), (0.37, 0.39), arrowprops={"arrowstyle": "|-|", "color": ORANGE, "lw": 2})
    ax.text(0.5, 0.34, "Kelvin–Voigt spring + dashpot", ha="center", va="top", fontsize=8, color=INK)
    ax.plot([0.5, 0.5], [0.52, 0.78], color=MUTED, linestyle=":", linewidth=1)
    ax.text(0.5, 0.82, "cell center", ha="center", va="bottom", fontsize=7.5, color=MUTED)
    ax.text(0.5, 0.19 if feedback else 0.08, "reference length L₀ sets the explicit mechanical size", ha="center", va="bottom", fontsize=8.3, fontweight="semibold")
    if feedback:
        ax.annotate("positive tension", (0.77, 0.25), (0.52, 0.25), fontsize=7.5, ha="center", arrowprops={"arrowstyle": "->", "color": PURPLE, "lw": 1.4})
        ax.annotate("front inhibition", (0.78, 0.54), (0.78, 0.28), fontsize=7.5, ha="center", arrowprops={"arrowstyle": "-[", "color": PURPLE, "lw": 1.4})


def moving_success(row: pd.Series, phenotype: str) -> float:
    if phenotype == "single-cue competence":
        return float(row.terminal_pole_count == 1)
    if phenotype == "persistent dual fronts":
        return float(row.terminal_pole_count >= 2)
    if phenotype == "clean reversal":
        return float(row.terminal_pole_count == 1 and row.terminal_signed_polarity > 0)
    raise KeyError(phenotype)


def moving_protocol_for_phenotype(phenotype: str) -> str:
    return {
        "single-cue competence": "single",
        "persistent dual fronts": "dual_equal",
        "clean reversal": "reversal",
    }[phenotype]


def plot_moving_phenotype_curve(
    ax: plt.Axes,
    frame: pd.DataFrame,
    phenotype: str,
    title: str,
    show_ylabel: bool,
    show_xlabel: bool,
) -> None:
    protocol = moving_protocol_for_phenotype(phenotype)
    subset = frame[
        frame.protocol.eq(protocol) & frame.tension_feedback_alpha.eq(0.0)
    ].copy()
    subset["success"] = subset.apply(lambda row: moving_success(row, phenotype), axis=1)
    for mechanical_case in ("soft", "stiff"):
        selected = subset[subset.mechanical_case.eq(mechanical_case)].sort_values("reference_length_um")
        ax.plot(
            selected.reference_length_um,
            selected.success,
            marker="o",
            markersize=3.6,
            linewidth=1.6,
            color=MECH_COLORS[mechanical_case],
            label=mechanical_case,
        )
    ax.set(
        title=title,
        xlabel="reference length L₀ (µm)" if show_xlabel else "",
        ylabel="outcome (0 / 1)" if show_ylabel else "",
        ylim=(-0.08, 1.08),
        yticks=[0, 1],
        yticklabels=["no", "yes"] if show_ylabel else ["", ""],
    )
    clean_axis(ax, grid=True)


def moving_size_effect(frame: pd.DataFrame, phenotype: str) -> float:
    protocol = moving_protocol_for_phenotype(phenotype)
    subset = frame[
        frame.protocol.eq(protocol)
        & frame.tension_feedback_alpha.eq(0.0)
        & frame.reference_length_um.isin([10.0, 20.0])
    ].copy()
    subset["success"] = subset.apply(lambda row: moving_success(row, phenotype), axis=1)
    means = subset.groupby("reference_length_um").success.mean() * 100.0
    return float(means.loc[20.0] - means.loc[10.0])


def figure_4(data: DataBundle) -> Path:
    fig = plt.figure(figsize=(19.0, 20.0))
    figure_title(fig, "Figure 4", "Coupling polarity to mechanics connects signaling state to migration")
    outer = GridSpec(4, 1, figure=fig, height_ratios=[1.55, 1.30, 0.85, 1.0], hspace=0.48, top=0.955, bottom=0.055)

    top = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[0.72, 2.65], wspace=0.22)
    ax_a = fig.add_subplot(top[0])
    draw_mechanochemical_schematic(ax_a)
    ax_a.set_title("Moving-cell model", pad=4)
    panel_label(ax_a, "A", x=-0.08)

    bgrid = GridSpecFromSubplotSpec(3, 4, subplot_spec=top[1], height_ratios=[1.45, 0.75, 0.60], hspace=0.18, wspace=0.23)
    b_trajectories: list[tuple[str, float, dict[str, np.ndarray]]] = []
    for species in ("single_species", "three_species"):
        for length in (10.0, 40.0):
            b_trajectories.append((species, length, moving_trajectory(species, "global", length, "soft", "single", 0.0)))
    vmax_by_species = {
        species: max(traj["membrane"][:, 0 if species == "single_species" else 1].max() for sp, _, traj in b_trajectories if sp == species)
        for species in ("single_species", "three_species")
    }
    for column, (species, length, trajectory) in enumerate(b_trajectories):
        field = 0 if species == "single_species" else 1
        ax_k = fig.add_subplot(bgrid[0, column])
        draw_moving_kymograph(
            ax_k,
            trajectory,
            field,
            f"{'1-species' if species == 'single_species' else '3-species'} · L₀={length:g} µm",
            vmax=vmax_by_species[species],
            ylabel=column == 0,
        )
        if column == 0:
            panel_label(ax_k, "B", x=-0.20)
        ax_t = fig.add_subplot(bgrid[1, column])
        times = trajectory["times"]
        pos = trajectory["position"]
        center = trajectory["center"]
        ax_t.plot(times, pos[:, 0], color=BLUE, linewidth=1.0, label="left edge")
        ax_t.plot(times, pos[:, -1], color=VERMILLION, linewidth=1.0, label="right edge")
        ax_t.plot(times, center, color=INK, linewidth=1.4, label="center")
        ax_t.set(xlabel="", ylabel="position" if column == 0 else "")
        ax_t.tick_params(labelbottom=False)
        clean_axis(ax_t, grid=True)
        ax_l = fig.add_subplot(bgrid[2, column])
        ax_l.plot(times, trajectory["length"], color=ORANGE, linewidth=1.4)
        ax_l.axhline(length, color=MUTED, linestyle=":", linewidth=0.8)
        ax_l.set(xlabel="time (s)", ylabel="length" if column == 0 else "")
        clean_axis(ax_l, grid=True)

    cgrid = GridSpecFromSubplotSpec(3, 2, subplot_spec=outer[1], hspace=0.42, wspace=0.25)
    c_axes = []
    phenotypes = ["single-cue competence", "persistent dual fronts", "clean reversal"]
    for row, phenotype in enumerate(phenotypes):
        for column, species in enumerate(("single_species", "three_species")):
            ax = fig.add_subplot(cgrid[row, column])
            plot_moving_phenotype_curve(
                ax,
                data.moving[(species, "global")],
                phenotype,
                f"{phenotype} · {'moving 1-species' if species == 'single_species' else 'moving 3-species'}",
                column == 0,
                row == 2,
            )
            if row == 0 and column == 0:
                ax.legend(loc="lower right")
            c_axes.append(ax)
    panel_label(c_axes[0], "C", x=-0.10)

    dgrid = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[2], wspace=0.28)
    measures = [
        ("net_displacement_um", "|net displacement| (µm)"),
        ("total_travel_um", "cumulative travel (µm)"),
        ("mean_polarity_projected_speed_established", "polarity-projected speed (µm s⁻¹)"),
    ]
    for column, (field, label) in enumerate(measures):
        ax = fig.add_subplot(dgrid[column])
        for species in ("single_species", "three_species"):
            if field == "total_travel_um":
                frame = data.directional[
                    data.directional.chemistry_class.eq(species)
                    & data.directional.tension_mode.eq("global")
                    & data.directional.protocol.eq("single")
                    & data.directional.tension_feedback_alpha.eq(0.0)
                ]
            else:
                frame = data.moving[(species, "global")]
                frame = frame[frame.protocol.eq("single") & frame.tension_feedback_alpha.eq(0.0)]
            for mechanical_case, linestyle in (("soft", "-"), ("stiff", "--")):
                selected = frame[frame.mechanical_case.eq(mechanical_case)].sort_values("reference_length_um")
                values = selected[field].abs() if field == "net_displacement_um" else selected[field]
                ax.plot(selected.reference_length_um, values, color=MODEL_COLORS[species], linestyle=linestyle, marker="o", markersize=3.2, linewidth=1.4, label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {mechanical_case}")
        ax.set(title=label, xlabel="reference length L₀ (µm)", ylabel=label)
        clean_axis(ax, grid=True)
        if column == 0:
            panel_label(ax, "D", x=-0.10)
        if column == 2:
            ax.legend(loc="best", ncol=2, fontsize=7.1)

    bottom = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[3], width_ratios=[1.2, 1.0], wspace=0.30)
    ax_e = fig.add_subplot(bottom[0])
    categories = ["single pole", "dual / multiple", "failed reversal", "depolarized / other"]
    values: dict[str, list[float]] = {category: [] for category in categories}
    for species in ("single_species", "three_species"):
        frame = data.moving[(species, "global")]
        frame = frame[frame.tension_feedback_alpha.eq(0.0)]
        for row in frame.itertuples():
            success = row.protocol == "reversal" and row.terminal_pole_count == 1 and row.terminal_signed_polarity > 0
            if row.protocol == "reversal" and not success:
                category = "failed reversal"
            elif row.terminal_pole_count >= 2:
                category = "dual / multiple"
            elif row.terminal_pole_count == 1:
                category = "single pole"
            else:
                category = "depolarized / other"
            values[category].append(abs(float(row.net_displacement_um)))
    rng = np.random.default_rng(7)
    for index, category in enumerate(categories):
        arr = np.asarray(values[category])
        if len(arr) == 0:
            continue
        x = index + rng.normal(0, 0.055, len(arr))
        ax_e.scatter(x, arr, s=13, color=[BLUE, MAGENTA, VERMILLION, GRAY][index], alpha=0.38, edgecolor="none")
        ax_e.boxplot(arr, positions=[index], widths=0.45, showfliers=False, patch_artist=True, boxprops={"facecolor": mpl_colors.to_rgba([BLUE, MAGENTA, VERMILLION, GRAY][index], 0.16), "edgecolor": INK}, medianprops={"color": INK, "linewidth": 1.5}, whiskerprops={"color": INK}, capprops={"color": INK})
    ax_e.set_xticks(range(len(categories)), categories, rotation=18, ha="right")
    ax_e.set(title="Polarity state explains migration output", ylabel="|net displacement| (µm)")
    clean_axis(ax_e, grid=True)
    panel_label(ax_e, "E", x=-0.08)

    ax_f = fig.add_subplot(bottom[1])
    rows = ["single-cue\ncompetence", "persistent\ndual fronts", "clean\nreversal", "directional\ndisplacement"]
    columns = ["fixed\n1-species", "fixed\n3-species", "moving\n1-species", "moving\n3-species"]
    fixed_raw = {species: raw_fixed_phenotypes(data, species) for species in ("single_species", "three_species")}
    matrix = np.full((4, 4), np.nan)
    for col, species in enumerate(("single_species", "three_species")):
        for row, phenotype in enumerate(phenotypes):
            series = fixed_raw[species][phenotype]
            matrix[row, col] = series.loc[2.0] - series.loc[1.0]
    for offset, species in enumerate(("single_species", "three_species"), start=2):
        for row, phenotype in enumerate(phenotypes):
            matrix[row, offset] = moving_size_effect(data.moving[(species, "global")], phenotype)
        frame = data.moving[(species, "global")]
        frame = frame[frame.protocol.eq("single") & frame.tension_feedback_alpha.eq(0.0) & frame.reference_length_um.isin([10.0, 20.0])]
        grouped = frame.assign(abs_disp=frame.net_displacement_um.abs()).groupby("reference_length_um").abs_disp.mean()
        matrix[3, offset] = 100.0 * (grouped.loc[20.0] - grouped.loc[10.0]) / max(grouped.loc[10.0], 1.0e-12)
    cmap = mpl_colors.LinearSegmentedColormap.from_list("effect", [VERMILLION, "#F7F7F7", BLUE])
    image = ax_f.imshow(matrix, cmap=cmap, norm=mpl_colors.TwoSlopeNorm(vmin=-100, vcenter=0, vmax=100), aspect="auto")
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            if np.isfinite(matrix[row, col]):
                unit = "%" if row == 3 else " pp"
                ax_f.text(col, row, f"{matrix[row, col]:+.1f}{unit}", ha="center", va="center", fontsize=7.5, color=INK, fontweight="semibold")
            else:
                ax_f.text(col, row, "—", ha="center", va="center", color=MUTED)
    ax_f.set_xticks(range(4), columns)
    ax_f.set_yticks(range(4), rows)
    ax_f.set_title("Cross-model size sensitivity\n(fixed: 1×→2×; moving: L₀ 10→20 µm)")
    ax_f.tick_params(length=0)
    for spine in ax_f.spines.values():
        spine.set_visible(False)
    panel_label(ax_f, "F", x=-0.14)

    conclusion_band(fig, "Mechanical coupling makes the polarity phenotype consequential for directional displacement.")
    return save_figure(fig, "figure_4_mechanochemical_coupling.png", "Figure 4. Polarity–mechanics coupling")


def feedback_map_matrix(frame: pd.DataFrame, observable: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    protocol = {
        "terminal polarized": "single",
        "clean reversal": "reversal",
        "persistent dual fronts": "dual_equal",
        "directional displacement": "single",
    }[observable]
    subset = frame[frame.protocol.eq(protocol)].copy()
    if observable == "terminal polarized":
        subset["value"] = (subset.terminal_pole_count == 1).astype(float)
    elif observable == "clean reversal":
        subset["value"] = ((subset.terminal_pole_count == 1) & (subset.terminal_signed_polarity > 0)).astype(float)
    elif observable == "persistent dual fronts":
        subset["value"] = (subset.terminal_pole_count >= 2).astype(float)
    else:
        subset["value"] = subset.net_displacement_um.abs()
    grouped = subset.groupby(["tension_feedback_alpha", "reference_length_um"]).value.mean().reset_index()
    xs = np.asarray(sorted(grouped.reference_length_um.unique()))
    ys = np.asarray(sorted(grouped.tension_feedback_alpha.unique()))
    lookup = {(row.tension_feedback_alpha, row.reference_length_um): row.value for row in grouped.itertuples()}
    matrix = np.asarray([[lookup[(y, x)] for x in xs] for y in ys], dtype=float)
    return xs, ys, matrix


def draw_feedback_modes(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    draw_mechanochemical_schematic(ax, feedback=True)
    ax.add_patch(patches.FancyBboxPatch((0.05, 0.005), 0.40, 0.105, boxstyle="round,pad=0.012", facecolor=mpl_colors.to_rgba(ORANGE, 0.12), edgecolor=ORANGE))
    ax.text(0.25, 0.057, "global: one tension value acts everywhere", ha="center", va="center", fontsize=6.5, wrap=True)
    ax.add_patch(patches.FancyBboxPatch((0.55, 0.005), 0.40, 0.105, boxstyle="round,pad=0.012", facecolor=mpl_colors.to_rgba(PURPLE, 0.12), edgecolor=PURPLE))
    ax.text(0.75, 0.057, "spatial: propagated local tension signal", ha="center", va="center", fontsize=6.5, wrap=True)


def feedback_condition(frame: pd.DataFrame, mode: str, alpha: float) -> pd.DataFrame:
    if mode == "none":
        return frame[frame.tension_feedback_alpha.eq(0.0)]
    return frame[frame.tension_feedback_alpha.eq(alpha)]


def figure_5(data: DataBundle) -> Path:
    fig = plt.figure(figsize=(19.0, 20.5))
    figure_title(fig, "Figure 5", "Tension feedback creates a mechanochemical constraint on polarity and migration")
    outer = GridSpec(4, 1, figure=fig, height_ratios=[1.45, 1.85, 0.90, 1.0], hspace=0.48, top=0.955, bottom=0.055)

    top = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[0.82, 2.4], wspace=0.26)
    ax_a = fig.add_subplot(top[0])
    draw_feedback_modes(ax_a)
    ax_a.set_title("Tension-to-polarity feedback", pad=3)
    panel_label(ax_a, "A", x=-0.08)

    bgrid = GridSpecFromSubplotSpec(3, 3, subplot_spec=top[1], height_ratios=[1.45, 0.72, 0.72], hspace=0.20, wspace=0.25)
    b_trajectories = [
        moving_trajectory("single_species", "global", 20.0, "soft", "single", alpha)
        for alpha in (0.0, 5.0, 20.0)
    ]
    b_vmax = max(trajectory["membrane"][:, 0].max() for trajectory in b_trajectories)
    for column, (alpha, trajectory) in enumerate(zip((0.0, 5.0, 20.0), b_trajectories)):
        ax_k = fig.add_subplot(bgrid[0, column])
        draw_moving_kymograph(ax_k, trajectory, 0, f"global αₜ={alpha:g}", vmax=b_vmax, ylabel=column == 0)
        if column == 0:
            panel_label(ax_k, "B", x=-0.18)
        ax_l = fig.add_subplot(bgrid[1, column])
        ax_l.plot(trajectory["times"], trajectory["length"], color=ORANGE, linewidth=1.3, label="length")
        tension = trajectory.get("elastic_tension", np.zeros_like(trajectory["times"]))
        ax_l2 = ax_l.twinx()
        ax_l2.plot(trajectory["times"], tension, color=PURPLE, linewidth=1.0, alpha=0.85, label="tension")
        ax_l.set(ylabel="length" if column == 0 else "")
        ax_l2.set_ylabel("tension" if column == 2 else "")
        ax_l.tick_params(labelbottom=False)
        clean_axis(ax_l, grid=True)
        ax_c = fig.add_subplot(bgrid[2, column])
        ax_c.plot(trajectory["times"], trajectory["center"] - trajectory["center"][0], color=INK, linewidth=1.35)
        ax_c.set(xlabel="time (s)", ylabel="center Δx" if column == 0 else "")
        clean_axis(ax_c, grid=True)

    cgrid = GridSpecFromSubplotSpec(4, 4, subplot_spec=outer[1], hspace=0.70, wspace=0.26)
    columns = [
        ("single_species", "global", "1-species · global"),
        ("single_species", "spatial", "1-species · spatial"),
        ("three_species", "global", "3-species · global"),
        ("three_species", "spatial", "3-species · spatial"),
    ]
    observables = ["terminal polarized", "clean reversal", "persistent dual fronts", "directional displacement"]
    map_cache = {(species, mode, obs): feedback_map_matrix(data.moving[(species, mode)], obs) for species, mode, _ in columns for obs in observables}
    displacement_max = max(matrix.max() for (species, mode, obs), (_, _, matrix) in map_cache.items() if obs == "directional displacement")
    c_axes = []
    for row, observable in enumerate(observables):
        for column, (species, mode, title) in enumerate(columns):
            ax = fig.add_subplot(cgrid[row, column])
            xs, ys, matrix = map_cache[(species, mode, observable)]
            if observable == "directional displacement":
                image = ax.imshow(matrix, origin="lower", aspect="auto", cmap="viridis", vmin=0, vmax=displacement_max, extent=(xs.min(), xs.max(), ys.min(), ys.max()))
            else:
                image = ax.imshow(matrix, origin="lower", aspect="auto", cmap=mpl_colors.LinearSegmentedColormap.from_list("binary", ["#ECEFF1", GREEN]), vmin=0, vmax=1, extent=(xs.min(), xs.max(), ys.min(), ys.max()))
            if row == 0:
                ax.set_title(title)
            ax.set_xlabel("L₀ (µm)" if row == 3 else "")
            ax.set_ylabel("feedback αₜ" if column == 0 else "")
            if column == 0:
                ax.text(-0.34, 0.5, observable, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=7.5, fontweight="semibold")
            ax.set_yticks(ys)
            ax.tick_params(length=2)
            c_axes.append(ax)
    panel_label(c_axes[0], "C", x=-0.23)
    fig.colorbar(c_axes[0].images[0], ax=c_axes[:12], location="right", shrink=0.62, pad=0.01, label="fraction across soft/stiff mechanics")
    fig.colorbar(c_axes[-1].images[0], ax=c_axes[-4:], location="right", shrink=0.75, pad=0.01, label="mean |net displacement| (µm)")

    dgrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[2], wspace=0.27)
    for column, species in enumerate(("single_species", "three_species")):
        ax = fig.add_subplot(dgrid[column])
        for mode, alpha, label, color in (
            ("global", 0.0, "no feedback", FEEDBACK_COLORS["none"]),
            ("global", 20.0, "global αₜ=20", FEEDBACK_COLORS["global"]),
            ("spatial", 20.0, "spatial αₜ=20", FEEDBACK_COLORS["spatial"]),
        ):
            frame = data.moving[(species, mode)]
            subset = frame[frame.protocol.eq("single") & frame.tension_feedback_alpha.eq(alpha)].copy()
            grouped = subset.assign(abs_disp=subset.net_displacement_um.abs()).groupby("reference_length_um").abs_disp.mean()
            ax.plot(grouped.index, grouped.values, marker="o", linewidth=1.8, markersize=4, color=color, label=label)
        ax.set(title="moving 1-species" if species == "single_species" else "moving 3-species", xlabel="reference length L₀ (µm)", ylabel="mean |net displacement| (µm)" if column == 0 else "")
        clean_axis(ax, grid=True)
        if column == 0:
            panel_label(ax, "D", x=-0.10)
            ax.legend()

    bottom = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[3], width_ratios=[1.15, 1.0], wspace=0.32)
    ax_e = fig.add_subplot(bottom[0])
    metric_rows = ["polarization", "reversal", "dual fronts", "displacement"]
    condition_columns = [
        ("single_species", "global", 0.0, "1-sp\nnone"),
        ("single_species", "global", 20.0, "1-sp\nglobal"),
        ("single_species", "spatial", 20.0, "1-sp\nspatial"),
        ("three_species", "global", 0.0, "3-sp\nnone"),
        ("three_species", "global", 20.0, "3-sp\nglobal"),
        ("three_species", "spatial", 20.0, "3-sp\nspatial"),
    ]
    matrix = np.zeros((4, 6))
    for col, (species, mode, alpha, _) in enumerate(condition_columns):
        frame = data.moving[(species, mode)]
        for row, (metric, protocol) in enumerate(zip(metric_rows[:3], ("single", "reversal", "dual_equal"))):
            subset = frame[frame.protocol.eq(protocol) & frame.tension_feedback_alpha.eq(alpha) & frame.reference_length_um.isin([10.0, 20.0])].copy()
            if metric == "polarization":
                subset["value"] = (subset.terminal_pole_count == 1).astype(float)
            elif metric == "reversal":
                subset["value"] = ((subset.terminal_pole_count == 1) & (subset.terminal_signed_polarity > 0)).astype(float)
            else:
                subset["value"] = (subset.terminal_pole_count >= 2).astype(float)
            grouped = subset.groupby("reference_length_um").value.mean()
            matrix[row, col] = 100.0 * (grouped.loc[20.0] - grouped.loc[10.0])
        subset = frame[frame.protocol.eq("single") & frame.tension_feedback_alpha.eq(alpha) & frame.reference_length_um.isin([10.0, 20.0])].copy()
        grouped = subset.assign(value=subset.net_displacement_um.abs()).groupby("reference_length_um").value.mean()
        matrix[3, col] = 100.0 * (grouped.loc[20.0] - grouped.loc[10.0]) / max(grouped.loc[10.0], 1.0e-12)
    image = ax_e.imshow(matrix, cmap=mpl_colors.LinearSegmentedColormap.from_list("effect", [VERMILLION, "#F7F7F7", BLUE]), norm=mpl_colors.TwoSlopeNorm(vmin=-100, vcenter=0, vmax=100), aspect="auto")
    for row in range(4):
        for col in range(6):
            ax_e.text(col, row, f"{matrix[row, col]:+.0f}{'%' if row == 3 else ' pp'}", ha="center", va="center", fontsize=7.2, fontweight="semibold")
    ax_e.set_xticks(range(6), [item[3] for item in condition_columns])
    ax_e.set_yticks(range(4), metric_rows)
    ax_e.set_title("Does tension amplify the L₀ 10→20 µm size effect?")
    ax_e.tick_params(length=0)
    for spine in ax_e.spines.values():
        spine.set_visible(False)
    panel_label(ax_e, "E", x=-0.09)

    ax_f = fig.add_subplot(bottom[1])
    trajectory = moving_trajectory("single_species", "spatial", 40.0, "soft", "single", 5.0)
    times = trajectory["times"]
    abs_polarity = np.abs(trajectory["polarity"])
    signal = trajectory["tension_feedback_signal"]
    signal_mean = signal.mean(axis=1) if signal.ndim == 2 else signal
    center = np.abs(trajectory["center"] - trajectory["center"][0])
    strain = np.abs(trajectory["length"] / 40.0 - 1.0)
    series = [
        (abs_polarity / max(abs_polarity.max(), 1e-12), MAGENTA, "front activity"),
        (signal_mean / max(signal_mean.max(), 1e-12), PURPLE, "tension signal"),
        (strain / max(strain.max(), 1e-12), ORANGE, "strain"),
        (center / max(center.max(), 1e-12), INK, "displacement"),
    ]
    for values, color, label in series:
        ax_f.plot(times, values, color=color, linewidth=1.7, label=label)
    ax_f.set(title="Matched mechanistic chain · large cell, spatial αₜ=5", xlabel="time (s)", ylabel="normalized trajectory")
    clean_axis(ax_f, grid=True)
    ax_f.legend(ncol=2)
    ax_f.text(0.02, 0.08, "size → polarity / strain → tension feedback → curtailed displacement", transform=ax_f.transAxes, fontsize=8.0, fontweight="semibold", bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": GRID})
    panel_label(ax_f, "F", x=-0.08)

    conclusion_band(fig, "Tension feedback suppresses sustained front activity; spatial feedback is the stronger mechanochemical constraint.")
    return save_figure(fig, "figure_5_tension_feedback.png", "Figure 5. Tension-feedback constraint")


def moving_effect_at_alpha(frame: pd.DataFrame, phenotype: str, alpha: float) -> float:
    protocol = moving_protocol_for_phenotype(phenotype)
    subset = frame[
        frame.protocol.eq(protocol)
        & frame.tension_feedback_alpha.eq(alpha)
        & frame.reference_length_um.isin([10.0, 20.0])
    ].copy()
    subset["success"] = subset.apply(lambda row: moving_success(row, phenotype), axis=1)
    means = subset.groupby("reference_length_um").success.mean() * 100.0
    return float(means.loc[20.0] - means.loc[10.0])


def displacement_effect(frame: pd.DataFrame, alpha: float) -> float:
    subset = frame[
        frame.protocol.eq("single")
        & frame.tension_feedback_alpha.eq(alpha)
        & frame.reference_length_um.isin([10.0, 20.0])
    ].copy()
    means = subset.assign(value=subset.net_displacement_um.abs()).groupby("reference_length_um").value.mean()
    return float(100.0 * (means.loc[20.0] - means.loc[10.0]) / max(means.loc[10.0], 1.0e-12))


def figure_6(data: DataBundle) -> Path:
    fig = plt.figure(figsize=(16.5, 7.8))
    figure_title(fig, "Figure 6", "Integrated model-richness summary")
    grid = GridSpec(1, 2, figure=fig, width_ratios=[1.35, 1.0], wspace=0.28, left=0.06, right=0.98, top=0.90, bottom=0.18)

    ax_a = fig.add_subplot(grid[0])
    rows = ["establish one front", "suppress competing fronts", "reverse polarity", "sustain directional migration"]
    columns = ["minimal\npolarity", "Rho-family\npolarity", "passive\nmechanics", "spatial tension\nfeedback"]
    fixed = {species: raw_fixed_phenotypes(data, species) for species in ("single_species", "three_species")}
    matrix = np.full((4, 4), np.nan)
    for col, species in enumerate(("single_species", "three_species")):
        matrix[0, col] = fixed[species]["single-cue competence"].loc[2.0] - fixed[species]["single-cue competence"].loc[1.0]
        matrix[1, col] = -(fixed[species]["persistent dual fronts"].loc[2.0] - fixed[species]["persistent dual fronts"].loc[1.0])
        matrix[2, col] = fixed[species]["clean reversal"].loc[2.0] - fixed[species]["clean reversal"].loc[1.0]
    passive_frames = [data.moving[(species, "global")] for species in ("single_species", "three_species")]
    tension_frames = [data.moving[(species, "spatial")] for species in ("single_species", "three_species")]
    matrix[0, 2] = np.mean([moving_effect_at_alpha(frame, "single-cue competence", 0.0) for frame in passive_frames])
    matrix[1, 2] = -np.mean([moving_effect_at_alpha(frame, "persistent dual fronts", 0.0) for frame in passive_frames])
    matrix[2, 2] = np.mean([moving_effect_at_alpha(frame, "clean reversal", 0.0) for frame in passive_frames])
    matrix[3, 2] = np.mean([displacement_effect(frame, 0.0) for frame in passive_frames])
    matrix[0, 3] = np.mean([moving_effect_at_alpha(frame, "single-cue competence", 20.0) for frame in tension_frames])
    matrix[1, 3] = -np.mean([moving_effect_at_alpha(frame, "persistent dual fronts", 20.0) for frame in tension_frames])
    matrix[2, 3] = np.mean([moving_effect_at_alpha(frame, "clean reversal", 20.0) for frame in tension_frames])
    matrix[3, 3] = np.mean([displacement_effect(frame, 20.0) for frame in tension_frames])
    matrix[np.abs(matrix) < 0.05] = 0.0
    image = ax_a.imshow(matrix, cmap=mpl_colors.LinearSegmentedColormap.from_list("summary", [VERMILLION, "#FAFAFA", BLUE]), norm=mpl_colors.TwoSlopeNorm(vmin=-100, vcenter=0, vmax=100), aspect="auto")
    for row in range(4):
        for col in range(4):
            value = matrix[row, col]
            if not np.isfinite(value):
                text_value = "not modeled"
            else:
                symbol = "≈" if abs(value) < 5 else "↑" if value > 0 else "↓"
                text_value = f"{symbol}\n{value:+.1f}{'%' if row == 3 else ' pp'}"
            ax_a.text(col, row, text_value, ha="center", va="center", fontsize=8.4, fontweight="semibold", color=INK)
    ax_a.set_xticks(range(4), columns)
    ax_a.set_yticks(range(4), rows)
    ax_a.tick_params(length=0)
    ax_a.set_title("Direction of the enlargement effect across model richness")
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    panel_label(ax_a, "A", x=-0.10)
    ax_a.text(0.5, -0.075, "fixed chemistry: 1×→2× volume  •  moving cells: L₀ 10→20 µm", transform=ax_a.transAxes, ha="center", va="top", fontsize=7.6, color=MUTED)
    ax_a.text(0.5, -0.115, "color encodes effect direction; intensity is clipped at |100|", transform=ax_a.transAxes, ha="center", va="top", fontsize=7.0, color=MUTED)

    ax_b = fig.add_subplot(grid[1])
    ax_b.set_axis_off()
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    panel_label(ax_b, "B", x=-0.06)
    ax_b.text(0.25, 0.94, "stable external cue", ha="center", va="top", fontsize=10, fontweight="semibold", color=BLUE)
    ax_b.text(0.75, 0.94, "changing / competing cues", ha="center", va="top", fontsize=10, fontweight="semibold", color=VERMILLION)
    for center in (0.25, 0.75):
        ax_b.add_patch(patches.FancyBboxPatch((center - 0.20, 0.53), 0.40, 0.17, boxstyle="round,pad=0.025,rounding_size=0.08", facecolor="#EDF3F6", edgecolor=INK, linewidth=1.1))
    ax_b.add_patch(patches.Ellipse((0.10, 0.615), 0.065, 0.11, facecolor=MAGENTA, edgecolor="none"))
    ax_b.annotate("cue", (0.06, 0.615), (0.01, 0.615), fontsize=7, arrowprops={"arrowstyle": "->", "color": BLUE})
    ax_b.text(0.25, 0.43, "one stable front forms\neven in the large cell", ha="center", va="top", fontsize=8.5, color=INK)
    ax_b.add_patch(patches.Ellipse((0.60, 0.615), 0.065, 0.11, facecolor=MAGENTA, edgecolor="none"))
    ax_b.add_patch(patches.Ellipse((0.90, 0.615), 0.065, 0.11, facecolor=MAGENTA, edgecolor="none"))
    ax_b.annotate("old", (0.56, 0.615), (0.51, 0.615), fontsize=7, arrowprops={"arrowstyle": "->", "color": BLUE})
    ax_b.annotate("new", (0.94, 0.615), (0.99, 0.615), fontsize=7, ha="right", arrowprops={"arrowstyle": "->", "color": ORANGE})
    ax_b.text(0.75, 0.43, "spatially separated fronts\ncoexist or replace slowly", ha="center", va="top", fontsize=8.5, color=INK)
    ax_b.annotate("robust polarization", (0.25, 0.28), (0.25, 0.36), ha="center", fontsize=8, color=BLUE, arrowprops={"arrowstyle": "->", "color": BLUE})
    ax_b.text(0.25, 0.18, "directed displacement sustained", ha="center", va="center", fontsize=7.8, color=BLUE, fontweight="semibold")
    ax_b.annotate("weak / indecisive migration", (0.75, 0.28), (0.75, 0.36), ha="center", fontsize=8, color=VERMILLION, arrowprops={"arrowstyle": "->", "color": VERMILLION})
    ax_b.text(0.75, 0.18, "directional output curtailed", ha="center", va="center", fontsize=7.8, color=VERMILLION, fontweight="semibold")
    ax_b.set_title("Biological interpretation")

    conclusion_band(fig, "Enlargement preserves polarity competence but weakens spatial winner-take-all selection and, with feedback, migration.")
    return save_figure(fig, "figure_6_integrated_summary.png", "Figure 6. Integrated model-richness summary")


def draw_numeric_map(
    ax: plt.Axes,
    frame: pd.DataFrame,
    x: str,
    y: str,
    value: str,
    title: str,
    xlabel: str,
    ylabel: str,
    cmap: str = "viridis",
    norm=None,
) -> matplotlib.image.AxesImage:
    xs, ys, matrix = grid_from_frame(frame, x, y, value)
    image = ax.imshow(
        matrix,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=(xs.min(), xs.max(), ys.min(), ys.max()),
        cmap=cmap,
        norm=norm,
    )
    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    return image


def supplementary_1(data: DataBundle) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(17.0, 9.8), constrained_layout=True)
    figure_title(fig, "Figure S1", "Single-species calibration and numerical validation")
    roots = data.one_roots
    ax = axes[0, 0]
    for stable, color, label in ((True, BLUE, "stable"), (False, VERMILLION, "unstable")):
        subset = roots[roots.stable.eq(stable)]
        ax.scatter(subset.inactive, subset.active, s=10, color=color, alpha=0.75, label=label)
    ax.set(title="Local reaction roots", xlabel="inactive v", ylabel="active u")
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "A")

    regime = data.one_regime.copy()
    regime["phenotype"] = regime.terminal_polarity_status.map(phenotype_from_status)
    for column, h in enumerate((1.0, 2.0, 4.0), start=0):
        ax = axes[0 if column < 2 else 1, 1 if column == 0 else 2 if column == 1 else 0]
        draw_categorical_map(
            ax,
            regime[regime.h.eq(h)],
            "Lambda",
            "mu",
            "phenotype",
            PHENOTYPE_ORDER,
            PHENOTYPE_COLORS,
            f"Single-cue regime map · h={h:g}",
            "Λ",
            "μ",
        )
        panel_label(ax, chr(ord("B") + column))

    ax = axes[1, 1]
    for epsilon, color in zip(sorted(data.one_dual_calibration.epsilon.unique()), plt.get_cmap("viridis")(np.linspace(0.15, 0.85, 4))):
        subset = data.one_dual_calibration[data.one_dual_calibration.epsilon.eq(epsilon)]
        prevalence = (subset.polarized_domain_count_t300 >= 2).groupby(subset.h).mean() * 100
        ax.plot(prevalence.index, prevalence.values, marker="o", color=color, label=f"ε={epsilon:g}")
    ax.set(title="Dual-domain sensitivity", xlabel="cytoplasm depth h", ylabel="persistent dual outcomes (%)")
    clean_axis(ax, grid=True)
    ax.legend(ncol=2)
    panel_label(ax, "E")

    ax = axes[1, 2]
    errors = np.concatenate(
        [
            data.one_regime.maximum_conservation_error.to_numpy(),
            data.one_dual_calibration.maximum_conservation_error.to_numpy(),
            data.one_reversal_calibration.maximum_conservation_error.to_numpy(),
            data.one_single.maximum_conservation_error.to_numpy(),
            data.one_dual.maximum_conservation_error.to_numpy(),
            data.one_reversal.maximum_conservation_error.to_numpy(),
        ]
    )
    errors = errors[np.isfinite(errors) & (errors > 0)]
    ax.hist(np.log10(errors), bins=36, color=BLUE, alpha=0.82)
    ax.axvline(-7, color=VERMILLION, linestyle="--", linewidth=1.2, label="10⁻⁷ acceptance limit")
    ax.set(title="Mass-conservation diagnostics", xlabel="log₁₀ relative conservation error", ylabel="simulations")
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "F")
    axes[1, 0].legend(handles=legend_handles(PHENOTYPE_ORDER, PHENOTYPE_COLORS), loc="upper right", fontsize=7.0)
    return save_figure(fig, "figure_s1_single_species_validation.png", "Figure S1. Single-species calibration and validation")


def supplementary_2(data: DataBundle) -> Path:
    fig, axes = plt.subplots(10, 3, figsize=(17.0, 27.0), constrained_layout=True)
    figure_title(fig, "Figure S2", "Full single-species size-response maps")
    volumes = (0.5, 1.0, 2.0)
    hs = (1.0, 2.0, 4.0)
    row = 0
    for h in hs:
        for column, volume in enumerate(volumes):
            subset = data.one_single[data.one_single.baseline_h.eq(h) & data.one_single.volume_multiplier.eq(volume)].copy()
            subset["phenotype"] = subset.terminal_outcome.map(phenotype_from_status)
            draw_categorical_map(
                axes[row, column], subset, "baseline_Lambda", "baseline_mu", "phenotype", PHENOTYPE_ORDER, PHENOTYPE_COLORS,
                f"single cue · h={h:g} · {volume:g}×", "baseline Λ", "baseline μ" if column == 0 else "",
            )
        row += 1
    dual_grouped = data.one_dual.groupby(["baseline_h", "baseline_mu", "baseline_Lambda", "volume_multiplier"], as_index=False).polarized_domain_count_t300.max()
    for h in hs:
        for column, volume in enumerate(volumes):
            subset = dual_grouped[dual_grouped.baseline_h.eq(h) & dual_grouped.volume_multiplier.eq(volume)].copy()
            subset["phenotype"] = subset.polarized_domain_count_t300.map(outcome_from_count)
            draw_categorical_map(
                axes[row, column], subset, "baseline_Lambda", "baseline_mu", "phenotype", PHENOTYPE_ORDER, PHENOTYPE_COLORS,
                f"competing seeds · h={h:g} · {volume:g}×", "baseline Λ", "baseline μ" if column == 0 else "",
            )
        row += 1
    for h in hs:
        for column, volume in enumerate(volumes):
            subset = data.one_reversal[data.one_reversal.baseline_h.eq(h) & data.one_reversal.volume_multiplier.eq(volume)]
            draw_categorical_map(
                axes[row, column], subset, "baseline_Lambda", "baseline_mu", "reversal_outcome", REVERSAL_ORDER, REVERSAL_COLORS,
                f"cue reversal · h={h:g} · {volume:g}×", "baseline Λ", "baseline μ" if column == 0 else "",
            )
        row += 1
    for column, volume in enumerate(volumes):
        ax = axes[row, column]
        subset = data.one_reversal[data.one_reversal.volume_multiplier.eq(volume)].reversal_time_after_challenge.dropna()
        ax.hist(subset, bins=np.arange(0, 165, 10), color=VOLUME_COLORS[volume], alpha=0.85)
        ax.set(title=f"successful reversal times · {volume:g}×", xlabel="time after challenge (s)", ylabel="conditions" if column == 0 else "")
        clean_axis(ax, grid=True)
    panel_label(axes[0, 0], "A", x=-0.14)
    panel_label(axes[3, 0], "B", x=-0.14)
    panel_label(axes[6, 0], "C", x=-0.14)
    panel_label(axes[9, 0], "D", x=-0.14)
    axes[0, 2].legend(handles=legend_handles(PHENOTYPE_ORDER, PHENOTYPE_COLORS), loc="upper left", bbox_to_anchor=(1.02, 1.0))
    axes[6, 2].legend(handles=legend_handles(REVERSAL_ORDER, REVERSAL_COLORS, REVERSAL_LABELS), loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return save_figure(fig, "figure_s2_single_species_size_maps.png", "Figure S2. Full single-species size-response maps")


def supplementary_3(data: DataBundle) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(17.0, 10.0), constrained_layout=True)
    figure_title(fig, "Figure S3", "Three-species calibration, signaling coherence, and conservation")
    regime = data.three_regime.copy()
    regime["phenotype"] = regime.terminal_polarity_status.map(phenotype_from_status)
    for column, hill in enumerate((2.0, 3.0, 4.0)):
        draw_categorical_map(
            axes[0, column], regime[regime.n.eq(hill)], "Lambda", "kappa_R", "phenotype", PHENOTYPE_ORDER, PHENOTYPE_COLORS,
            f"single Rac cue · n={hill:g}", "Λ", "κR" if column == 0 else "",
        )
        panel_label(axes[0, column], chr(ord("A") + column))
    ax = axes[1, 0]
    trajectory = simulate_fixed("three_species", "single", 1.0)
    fraction = material_fraction(trajectory["position"])[-1]
    for profile, color, label in zip(trajectory["membrane"][-1], (SKY, MAGENTA, ORANGE), ("Cdc42", "Rac", "Rho")):
        ax.plot(fraction, profile, color=color, linewidth=2.0, label=label)
    ax.set(title="Organized late signaling profile", xlabel="cell position", ylabel="active concentration")
    ax.set_xticks([0, 0.5, 1], ["left", "mid", "right"])
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "D")
    ax = axes[1, 1]
    scatter = ax.scatter(regime.cdc42_rac_correlation, regime.cdc42_rho_correlation, c=regime.rac_range_contrast, cmap="viridis", s=15, alpha=0.65)
    ax.axhline(0, color=GRID, linewidth=0.8)
    ax.axvline(0, color=GRID, linewidth=0.8)
    ax.set(title="Network coherence metrics", xlabel="Cdc42–Rac correlation", ylabel="Cdc42–Rho correlation")
    fig.colorbar(scatter, ax=ax, label="Rac range contrast")
    clean_axis(ax)
    panel_label(ax, "E")
    ax = axes[1, 2]
    error_columns = ["maximum_cdc42_conservation_error", "maximum_rac_conservation_error", "maximum_rho_conservation_error"]
    for field, color, label in zip(error_columns, (SKY, MAGENTA, ORANGE), ("Cdc42", "Rac", "Rho")):
        values = regime[field].to_numpy()
        values = values[np.isfinite(values) & (values > 0)]
        ax.hist(np.log10(values), bins=28, histtype="step", linewidth=1.6, color=color, label=label)
    ax.axvline(-6, color=VERMILLION, linestyle="--", linewidth=1.1, label="10⁻⁶ limit")
    ax.set(title="Species-wise conservation", xlabel="log₁₀ relative error", ylabel="simulations")
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "F")
    axes[0, 2].legend(handles=legend_handles(PHENOTYPE_ORDER, PHENOTYPE_COLORS), loc="lower right", fontsize=7.0)
    return save_figure(fig, "figure_s3_three_species_validation.png", "Figure S3. Three-species calibration and validation")


def supplementary_4(data: DataBundle) -> Path:
    fig, axes = plt.subplots(10, 3, figsize=(17.0, 27.0), constrained_layout=True)
    figure_title(fig, "Figure S4", "Full three-species size-response maps")
    volumes = (0.5, 1.0, 2.0)
    hills = (2.0, 3.0, 4.0)
    row = 0
    for hill in hills:
        for column, volume in enumerate(volumes):
            subset = data.three_single[data.three_single.n.eq(hill) & data.three_single.volume_multiplier.eq(volume)].copy()
            subset["phenotype"] = subset.terminal_polarity_status.map(phenotype_from_status)
            draw_categorical_map(
                axes[row, column], subset, "baseline_Lambda", "kappa_R", "phenotype", PHENOTYPE_ORDER, PHENOTYPE_COLORS,
                f"single cue · n={hill:g} · {volume:g}×", "baseline Λ", "κR" if column == 0 else "",
            )
        row += 1
    dual_grouped = data.three_dual.groupby(["n", "kappa_R", "baseline_Lambda", "volume_multiplier"], as_index=False).polarized_domain_count_t300.max()
    for hill in hills:
        for column, volume in enumerate(volumes):
            subset = dual_grouped[dual_grouped.n.eq(hill) & dual_grouped.volume_multiplier.eq(volume)].copy()
            subset["phenotype"] = subset.polarized_domain_count_t300.map(outcome_from_count)
            draw_categorical_map(
                axes[row, column], subset, "baseline_Lambda", "kappa_R", "phenotype", PHENOTYPE_ORDER, PHENOTYPE_COLORS,
                f"competing Rac seeds · n={hill:g} · {volume:g}×", "baseline Λ", "κR" if column == 0 else "",
            )
        row += 1
    for hill in hills:
        for column, volume in enumerate(volumes):
            subset = data.three_reversal[data.three_reversal.n.eq(hill) & data.three_reversal.volume_multiplier.eq(volume)]
            draw_categorical_map(
                axes[row, column], subset, "baseline_Lambda", "kappa_R", "reversal_outcome", REVERSAL_ORDER, REVERSAL_COLORS,
                f"Rac cue reversal · n={hill:g} · {volume:g}×", "baseline Λ", "κR" if column == 0 else "",
            )
        row += 1
    for column, volume in enumerate(volumes):
        ax = axes[row, column]
        subset = data.three_reversal[data.three_reversal.volume_multiplier.eq(volume)].reversal_time_after_challenge.dropna()
        ax.hist(subset, bins=np.arange(0, 165, 10), color=VOLUME_COLORS[volume], alpha=0.85)
        ax.set(title=f"successful reversal times · {volume:g}×", xlabel="time after challenge (s)", ylabel="conditions" if column == 0 else "")
        clean_axis(ax, grid=True)
    panel_label(axes[0, 0], "A", x=-0.14)
    panel_label(axes[3, 0], "B", x=-0.14)
    panel_label(axes[6, 0], "C", x=-0.14)
    panel_label(axes[9, 0], "D", x=-0.14)
    axes[0, 2].legend(handles=legend_handles(PHENOTYPE_ORDER, PHENOTYPE_COLORS), loc="upper left", bbox_to_anchor=(1.02, 1.0))
    axes[6, 2].legend(handles=legend_handles(REVERSAL_ORDER, REVERSAL_COLORS, REVERSAL_LABELS), loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return save_figure(fig, "figure_s4_three_species_size_maps.png", "Figure S4. Full three-species size-response maps")


def supplementary_5(data: DataBundle) -> Path:
    fig, axes = plt.subplots(3, 3, figsize=(17.0, 13.5), constrained_layout=True)
    figure_title(fig, "Figure S5", "Moving-domain mechanical validation")
    no_feedback = {
        species: data.moving[(species, "global")][data.moving[(species, "global")].tension_feedback_alpha.eq(0.0)].copy()
        for species in ("single_species", "three_species")
    }
    measures = [
        ("final_length_um", "final length (µm)", "Final versus reference length"),
        ("maximum_absolute_strain", "maximum |strain|", "Deformation range"),
        ("net_displacement_um", "|net displacement| (µm)", "No-feedback displacement"),
        ("mean_absolute_speed_established", "mean speed (µm s⁻¹)", "Established-front speed"),
        ("terminal_pole_count", "terminal pole count", "Terminal polarity domains"),
    ]
    positions = [(0, 0), (0, 1), (1, 0), (1, 1), (1, 2)]
    measure_labels = ["A", "B", "D", "E", "F"]
    for index, ((field, ylabel, title), position) in enumerate(zip(measures, positions)):
        ax = axes[position]
        for species in ("single_species", "three_species"):
            frame = no_feedback[species]
            frame = frame[frame.protocol.eq("single")]
            for mechanical_case, linestyle in (("soft", "-"), ("stiff", "--")):
                subset = frame[frame.mechanical_case.eq(mechanical_case)].sort_values("reference_length_um")
                values = subset[field].abs() if field == "net_displacement_um" else subset[field]
                ax.plot(subset.reference_length_um, values, color=MODEL_COLORS[species], linestyle=linestyle, marker="o", markersize=3.2, linewidth=1.4, label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {mechanical_case}")
        if field == "final_length_um":
            ax.plot([5, 50], [5, 50], color=GRAY, linestyle=":", linewidth=1.0, label="identity")
        ax.set(title=title, xlabel="reference length L₀ (µm)", ylabel=ylabel)
        clean_axis(ax, grid=True)
        panel_label(ax, measure_labels[index])
    axes[0, 0].legend(ncol=2, fontsize=7)

    ax = axes[0, 2]
    for mechanical_case, color in MECH_COLORS.items():
        frame = no_feedback["single_species"]
        subset = frame[frame.protocol.eq("single") & frame.mechanical_case.eq(mechanical_case)].sort_values("reference_length_um")
        ax.plot(subset.reference_length_um, subset.spring_stiffness * subset.reference_length_um, color=color, marker="o", label=f"kL₀ · {mechanical_case}")
        ax.plot(subset.reference_length_um, subset.dashpot_viscosity * subset.reference_length_um, color=color, marker="s", linestyle="--", label=f"ηL₀ · {mechanical_case}")
    ax.set(title="Mechanical normalization", xlabel="reference length L₀ (µm)", ylabel="length-normalized coefficient")
    clean_axis(ax, grid=True)
    ax.legend(ncol=2, fontsize=7)
    panel_label(ax, "C")

    representative_trajectories = {
        mechanical_case: moving_trajectory("single_species", "global", 20.0, mechanical_case, "single", 0.0)
        for mechanical_case in ("soft", "stiff")
    }
    representative_vmax = max(item["membrane"][:, 0].max() for item in representative_trajectories.values())
    representative_image = None
    for column, mechanical_case in enumerate(("soft", "stiff")):
        ax = axes[2, column]
        trajectory = representative_trajectories[mechanical_case]
        representative_image = draw_moving_kymograph(ax, trajectory, 0, f"No feedback · {mechanical_case} mechanics", vmax=representative_vmax, ylabel=column == 0)
        inset = ax.inset_axes([0.53, 0.07, 0.43, 0.28])
        inset.plot(trajectory["times"], trajectory["length"], color=ORANGE, linewidth=1.0)
        inset.axhline(20, color=MUTED, linestyle=":", linewidth=0.7)
        inset.set(xticks=[], yticks=[], facecolor=(1, 1, 1, 0.88))
        panel_label(ax, "G" if column == 0 else "H")
    fig.colorbar(representative_image, ax=[axes[2, 0], axes[2, 1]], location="bottom", shrink=0.78, pad=0.11, label="active u (common scale)")
    ax = axes[2, 2]
    error_fields = {
        "single_species": "maximum_active_gtpase_conservation_error",
        "three_species": "maximum_rac_conservation_error",
    }
    for species in ("single_species", "three_species"):
        frame = no_feedback[species]
        grouped = frame.groupby("reference_length_um")[error_fields[species]].max()
        ax.semilogy(grouped.index, grouped.values, color=MODEL_COLORS[species], marker="o", label="1-species" if species == "single_species" else "Rac in 3-species")
    ax.axhline(1e-7, color=VERMILLION, linestyle="--", linewidth=1.0, label="acceptance limit")
    ax.set(title="Numerical accuracy across moving meshes", xlabel="reference length L₀ (µm)", ylabel="maximum relative conservation error")
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "I")
    return save_figure(fig, "figure_s5_mechanical_validation.png", "Figure S5. Moving-domain mechanical validation")


def supplementary_6(data: DataBundle) -> Path:
    metrics = [
        ("terminal_absolute_polarity", "terminal |polarity|"),
        ("terminal_pole_count", "terminal domain count"),
        ("net_displacement_um", "net displacement (µm)"),
        ("mean_absolute_speed_late", "late speed (µm s⁻¹)"),
        ("final_length_um", "final length (µm)"),
        ("maximum_absolute_strain", "maximum |strain|"),
        ("mean_polarity_projected_speed_late", "late polarity-projected speed"),
    ]
    protocols = ["single", "dual_equal", "dual_biased", "reversal"]
    fig, axes = plt.subplots(len(metrics), len(protocols), figsize=(18.0, 23.0), constrained_layout=True)
    figure_title(fig, "Figure S6", "Full reference-length phenotype atlas without tension feedback")
    for row, (field, ylabel) in enumerate(metrics):
        for column, protocol in enumerate(protocols):
            ax = axes[row, column]
            for species in ("single_species", "three_species"):
                frame = data.moving[(species, "global")]
                frame = frame[frame.protocol.eq(protocol) & frame.tension_feedback_alpha.eq(0.0)]
                for mechanical_case, linestyle in (("soft", "-"), ("stiff", "--")):
                    subset = frame[frame.mechanical_case.eq(mechanical_case)].sort_values("reference_length_um")
                    ax.plot(subset.reference_length_um, subset[field], color=MODEL_COLORS[species], linestyle=linestyle, marker="o", markersize=2.5, linewidth=1.1, label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {mechanical_case}")
            if row == 0:
                ax.set_title(protocol.replace("_", " "))
            ax.set_xlabel("L₀ (µm)" if row == len(metrics) - 1 else "")
            ax.set_ylabel(ylabel if column == 0 else "")
            clean_axis(ax, grid=True)
    axes[0, 0].legend(ncol=2, fontsize=7)
    for row, label in enumerate("ABCDEFG"):
        panel_label(axes[row, 0], label, x=-0.20)
    return save_figure(fig, "figure_s6_reference_length_atlas.png", "Figure S6. Full reference-length phenotype atlas")


def supplementary_7(data: DataBundle) -> Path:
    fig, axes = plt.subplots(2, 4, figsize=(19.0, 10.5), constrained_layout=True)
    figure_title(fig, "Figure S7", "Tension-feedback diagnostics")
    for column, mode in enumerate(("global", "spatial")):
        ax = axes[0, column]
        for species in ("single_species", "three_species"):
            frame = data.moving[(species, mode)]
            frame = frame[frame.protocol.eq("single")]
            for alpha, linestyle in ((5.0, "-"), (20.0, "--")):
                grouped = frame[frame.tension_feedback_alpha.eq(alpha)].groupby("reference_length_um").integrated_tension_feedback_exposure.mean()
                ax.plot(grouped.index, grouped.values, color=MODEL_COLORS[species], linestyle=linestyle, marker="o", markersize=3, label=f"{'1-sp' if species == 'single_species' else '3-sp'} · αₜ={alpha:g}")
        ax.set(title=f"{mode.capitalize()} feedback exposure", xlabel="L₀ (µm)", ylabel="integrated added deactivation")
        clean_axis(ax, grid=True)
        panel_label(ax, "A" if column == 0 else "B")
        if column == 0:
            ax.legend(ncol=2, fontsize=7)

    ax = axes[0, 2]
    for species in ("single_species", "three_species"):
        for mode, marker in (("global", "o"), ("spatial", "s")):
            frame = data.moving[(species, mode)]
            subset = frame[frame.tension_feedback_alpha.gt(0.0)]
            ax.scatter(
                subset.maximum_elastic_tension,
                subset.maximum_tension_added_deactivation_rate,
                s=10 + 1.5 * subset.tension_feedback_alpha,
                color=MODEL_COLORS[species],
                marker=marker,
                alpha=0.50,
                edgecolor="none",
                label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {mode}",
            )
    ax.set(title="Tension-to-deactivation coupling", xlabel="maximum elastic tension", ylabel="maximum added deactivation rate")
    clean_axis(ax, grid=True)
    ax.legend(ncol=2, fontsize=6.8)
    panel_label(ax, "C")

    ax = axes[0, 3]
    for species in ("single_species", "three_species"):
        for mode, linestyle in (("global", "-"), ("spatial", "--")):
            frame = data.moving[(species, mode)]
            grouped = frame[frame.protocol.eq("single")].groupby("tension_feedback_alpha").fraction_time_under_tension.mean()
            ax.plot(grouped.index, grouped.values, color=MODEL_COLORS[species], linestyle=linestyle, marker="o", label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {mode}")
    ax.set(title="Time under positive elastic tension", xlabel="feedback αₜ", ylabel="trajectory fraction")
    clean_axis(ax, grid=True)
    ax.legend(ncol=2, fontsize=6.8)
    panel_label(ax, "D")

    ax = axes[1, 0]
    for length, color in ((10.0, BLUE), (40.0, VERMILLION)):
        trajectory = moving_trajectory("single_species", "spatial", length, "soft", "single", 5.0)
        signal = trajectory["tension_feedback_signal"]
        index = int(np.argmax(signal.mean(axis=1))) if signal.ndim == 2 else int(np.argmax(signal))
        fraction = material_fraction(trajectory["position"])[index]
        profile = signal[index] if signal.ndim == 2 else np.full_like(fraction, signal[index])
        ax.plot(fraction, profile, color=color, linewidth=2, label=f"L₀={length:g} µm")
    ax.set(title="Spatial fields at each trajectory's peak mean signal", xlabel="material position", ylabel="tension-associated signal")
    ax.set_xticks([0, 0.5, 1], ["left", "mid", "right"])
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "E")

    ax = axes[1, 1]
    for species in ("single_species", "three_species"):
        for mode, marker in (("global", "o"), ("spatial", "s")):
            frame = data.directional[
                data.directional.chemistry_class.eq(species)
                & data.directional.tension_mode.eq(mode)
            ]
            grouped = frame.groupby("tension_feedback_alpha").total_travel_um.mean()
            ax.plot(grouped.index, grouped.values, color=MODEL_COLORS[species], marker=marker, linewidth=1.7, label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {mode}")
    ax.set(title="Feedback suppresses total travel\n(all protocols and mechanics)", xlabel="feedback αₜ", ylabel="mean cumulative travel (µm)")
    clean_axis(ax, grid=True)
    ax.legend(ncol=2, fontsize=7)
    panel_label(ax, "F")

    ax = axes[1, 2]
    for species in ("single_species", "three_species"):
        for mode, linestyle in (("global", "-"), ("spatial", "--")):
            frame = data.moving[(species, mode)]
            subset = frame[frame.protocol.eq("single")].copy()
            subset["polarized"] = subset.terminal_pole_count.eq(1)
            grouped = subset.groupby("tension_feedback_alpha").polarized.mean() * 100
            ax.plot(grouped.index, grouped.values, color=MODEL_COLORS[species], linestyle=linestyle, marker="o", label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {mode}")
    ax.set(title="Terminal one-domain retention", xlabel="feedback αₜ", ylabel="conditions (%)")
    clean_axis(ax, grid=True)
    ax.legend(ncol=2, fontsize=7)
    panel_label(ax, "G")

    ax = axes[1, 3]
    records = []
    for (species, mode), frame in data.moving.items():
        error_fields = [field for field in frame.columns if field.startswith("maximum_") and field.endswith("_conservation_error")]
        for field in error_fields:
            for value in frame[field].dropna():
                if value > 0:
                    records.append((species, mode, math.log10(value)))
    for index, ((species, mode), group) in enumerate(pd.DataFrame(records, columns=["species", "mode", "log_error"]).groupby(["species", "mode"])):
        ax.hist(group.log_error, bins=28, histtype="step", linewidth=1.4, color=MODEL_COLORS[species], linestyle="-" if mode == "global" else "--", label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {mode}")
    ax.axvline(-7, color=VERMILLION, linestyle=":", label="10⁻⁷")
    ax.set(title="Conservation with feedback", xlabel="log₁₀ relative error", ylabel="trajectories")
    clean_axis(ax, grid=True)
    ax.legend(ncol=2, fontsize=6.8)
    panel_label(ax, "H")
    return save_figure(fig, "figure_s7_tension_diagnostics.png", "Figure S7. Tension-feedback diagnostics")


def conserved_error_from_arrays(trajectory: dict[str, np.ndarray], species: str) -> float:
    position = trajectory["position"]
    intervals = np.diff(position, axis=1)
    weights = np.empty_like(position)
    weights[:, 0] = 0.5 * intervals[:, 0]
    weights[:, -1] = 0.5 * intervals[:, -1]
    weights[:, 1:-1] = 0.5 * (intervals[:, :-1] + intervals[:, 1:])
    membrane = np.sum(trajectory["membrane"] * weights[:, None, :], axis=2)
    cytoplasm = np.sum(trajectory["cytoplasm"] * weights[:, None, :], axis=2)
    total = membrane + cytoplasm
    relative = np.abs(total - total[0]) / np.maximum(np.abs(total[0]), 1e-14)
    return float(relative.max())


def relative_threshold_domains(profile: np.ndarray, fraction: float) -> int:
    from scipy.signal import find_peaks

    threshold = float(profile.min() + fraction * (profile.max() - profile.min()))
    above = profile >= threshold
    # Count contiguous above-threshold regions; endpoints are legitimate poles.
    transitions = np.diff(np.r_[False, above, False].astype(int))
    return int(np.sum(transitions == 1))


def supplementary_8(data: DataBundle) -> Path:
    fig, axes = plt.subplots(2, 4, figsize=(19.0, 10.5), constrained_layout=True)
    figure_title(fig, "Figure S8", "Numerical and classification robustness")
    resolutions = [41, 81, 161]
    resolution_data = {
        species: {nodes: simulate_fixed(species, "single", 1.0, nodes=nodes) for nodes in resolutions}
        for species in ("single_species", "three_species")
    }
    for column, species in enumerate(("single_species", "three_species")):
        ax = axes[0, column]
        field = 0 if species == "single_species" else 1
        for nodes, color in zip(resolutions, plt.get_cmap("viridis")(np.linspace(0.15, 0.85, 3))):
            trajectory = resolution_data[species][nodes]
            fraction = material_fraction(trajectory["position"])[-1]
            ax.plot(fraction, trajectory["membrane"][-1, field], color=color, linewidth=1.6, label=f"{nodes} nodes")
        ax.set(title=f"{'1-species u' if species == 'single_species' else '3-species Rac'} resolution", xlabel="material position", ylabel="terminal active profile" if column == 0 else "")
        ax.set_xticks([0, 0.5, 1], ["left", "mid", "right"])
        clean_axis(ax, grid=True)
        ax.legend()
        panel_label(ax, "A" if column == 0 else "B")

    ax = axes[0, 2]
    for species, marker in (("single_species", "o"), ("three_species", "s")):
        field = 0 if species == "single_species" else 1
        amplitudes = []
        for nodes in resolutions:
            terminal = resolution_data[species][nodes]["membrane"][-1, field]
            amplitudes.append(float(terminal.max() - terminal.min()))
        ax.plot(resolutions, amplitudes, marker=marker, color=MODEL_COLORS[species], label="1-species" if species == "single_species" else "3-species Rac")
    ax.set(title="Resolution convergence", xlabel="spatial nodes", ylabel="terminal range amplitude")
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "C")

    ax = axes[0, 3]
    for species, marker in (("single_species", "o"), ("three_species", "s")):
        errors = [conserved_error_from_arrays(resolution_data[species][nodes], species) for nodes in resolutions]
        ax.semilogy(resolutions, errors, marker=marker, color=MODEL_COLORS[species], label="1-species" if species == "single_species" else "3-species")
    ax.axhline(1e-7, color=VERMILLION, linestyle="--", label="10⁻⁷")
    ax.set(title="Integration / conservation convergence", xlabel="spatial nodes", ylabel="maximum relative error")
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "D")

    ax = axes[1, 0]
    fractions = np.linspace(0.25, 0.75, 21)
    for species, color in MODEL_COLORS.items():
        field = 0 if species == "single_species" else 1
        for volume, linestyle in ((1.0, "-"), (2.0, "--")):
            trajectory = simulate_fixed(species, "reversal", volume)
            counts = [relative_threshold_domains(trajectory["membrane"][-1, field], fraction) for fraction in fractions]
            ax.plot(fractions, counts, color=color, linestyle=linestyle, linewidth=1.6, label=f"{'1-sp' if species == 'single_species' else '3-sp'} · {volume:g}×")
    ax.set(title="Domain-count threshold sensitivity", xlabel="relative high-state threshold", ylabel="terminal domains", yticks=[0, 1, 2, 3])
    clean_axis(ax, grid=True)
    ax.legend(ncol=2, fontsize=7)
    panel_label(ax, "E")

    ax = axes[1, 1]
    subset = data.one_single[data.one_single.volume_multiplier.eq(1.0)].copy()
    category_order = [
        category
        for category in ("homogeneous_low", "homogeneous_high", "single_polarized", "dual_polarized")
        if category in set(subset.terminal_outcome)
    ]
    mapping = {category: index for index, category in enumerate(category_order)}
    category_colors = {
        "homogeneous_low": GRAY,
        "homogeneous_high": MUTED,
        "single_polarized": BLUE,
        "dual_polarized": MAGENTA,
    }
    rng = np.random.default_rng(9)
    for category in category_order:
        group = subset[subset.terminal_outcome.eq(category)]
        jitter = rng.normal(0, 0.045, len(group))
        ax.scatter(
            mapping[category] + jitter,
            group.terminal_polarization_amplitude,
            s=10,
            alpha=0.38,
            color=category_colors[category],
            edgecolors="none",
        )
    amplitude_cutoff = 0.5
    ax.axhline(amplitude_cutoff, color=VERMILLION, linestyle="--", linewidth=1.1)
    ax.text(0.02, amplitude_cutoff + 0.04, "polarization cutoff = 0.5", transform=ax.get_yaxis_transform(), fontsize=7.2, color=VERMILLION, va="bottom")
    homogeneous = subset[subset.terminal_outcome.str.startswith("homogeneous")]
    polarized = subset[subset.terminal_outcome.isin(["single_polarized", "dual_polarized"])]
    boundary_examples = [
        homogeneous.loc[(homogeneous.terminal_polarization_amplitude - amplitude_cutoff).abs().idxmin()],
        polarized.loc[(polarized.terminal_polarization_amplitude - amplitude_cutoff).abs().idxmin()],
    ]
    annotation_offsets = [(10, 30), (-78, 18)]
    for row, offset in zip(boundary_examples, annotation_offsets):
        x_value = mapping[str(row.terminal_outcome)]
        y_value = float(row.terminal_polarization_amplitude)
        ax.scatter(x_value, y_value, s=55, marker="D", facecolor="white", edgecolor=INK, linewidth=1.1, zorder=5)
        ax.annotate(
            f"Λ={row.baseline_Lambda:g}, μ={row.baseline_mu:g}, h={row.baseline_h:g}",
            (x_value, y_value),
            xytext=offset,
            textcoords="offset points",
            fontsize=6.8,
            color=INK,
            arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.7},
        )
    ax.set_xticks(range(len(category_order)), [value.replace("_", "\n") for value in category_order], rotation=20, ha="right")
    ax.set(title="Classification boundary audit · 1× volume", ylabel="terminal polarization amplitude")
    clean_axis(ax, grid=True)
    panel_label(ax, "F")

    ax = axes[1, 2]
    all_errors = []
    labels = []
    for label, frame, fields in (
        ("1-sp fixed", data.one_single, ["maximum_conservation_error"]),
        ("3-sp fixed", data.three_single, ["maximum_cdc42_conservation_error", "maximum_rac_conservation_error", "maximum_rho_conservation_error"]),
        ("1-sp moving", data.moving[("single_species", "global")], ["maximum_active_gtpase_conservation_error"]),
        ("3-sp moving", data.moving[("three_species", "global")], ["maximum_cdc42_conservation_error", "maximum_rac_conservation_error", "maximum_rho_conservation_error"]),
    ):
        values = np.concatenate([frame[field].dropna().to_numpy() for field in fields])
        values = np.log10(values[values > 0])
        all_errors.append(values)
        labels.append(label)
    ax.boxplot(all_errors, tick_labels=labels, showfliers=False, patch_artist=True, boxprops={"facecolor": "#E7EEF2", "edgecolor": INK}, medianprops={"color": VERMILLION, "linewidth": 1.5})
    ax.axhline(-7, color=VERMILLION, linestyle="--", linewidth=1)
    ax.tick_params(axis="x", rotation=18)
    ax.set(title="Numerical error across model classes", ylabel="log₁₀ relative conservation error")
    clean_axis(ax, grid=True)
    panel_label(ax, "G")

    ax = axes[1, 3]
    for species, color in MODEL_COLORS.items():
        frame = data.moving[(species, "global")]
        error_field = "maximum_active_gtpase_conservation_error" if species == "single_species" else "maximum_rac_conservation_error"
        grouped = frame.groupby("nodes")[error_field].max()
        ax.semilogy(grouped.index, grouped.values, marker="o", color=color, label="1-species" if species == "single_species" else "3-species Rac")
    ax.axhline(1e-7, color=VERMILLION, linestyle="--")
    ax.set(title="Adaptive mesh-size audit", xlabel="nodes used by length-scaled mesh", ylabel="maximum relative error")
    clean_axis(ax, grid=True)
    ax.legend()
    panel_label(ax, "H")
    return save_figure(fig, "figure_s8_numerical_robustness.png", "Figure S8. Numerical and classification robustness")


def validate_inputs(data: DataBundle) -> dict[str, object]:
    expected_rows = {
        "one_single": 2601,
        "one_dual": 6084,
        "one_reversal": 2601,
        "three_single": 1521,
        "three_dual": 6084,
        "three_reversal": 2601,
        "directional": 960,
    }
    for name, expected in expected_rows.items():
        observed = len(getattr(data, name))
        if observed != expected:
            raise ValueError(f"{name} has {observed} rows; expected {expected}")
    for key, frame in data.moving.items():
        if len(frame) != 240:
            raise ValueError(f"moving campaign {key} has {len(frame)} rows; expected 240")

    one_table = endpoint_state_table(data.one_single, "terminal_outcome", fixed_key_columns("single_species"))
    one_states = sorted(set(one_table.index).union(one_table.columns))
    one_unchanged = int(np.trace(one_table.reindex(index=one_states, columns=one_states, fill_value=0).to_numpy()))
    three_table = endpoint_state_table(data.three_single, "terminal_polarity_status", fixed_key_columns("three_species"))
    three_states = sorted(set(three_table.index).union(three_table.columns))
    three_unchanged = int(np.trace(three_table.reindex(index=three_states, columns=three_states, fill_value=0).to_numpy()))
    if one_unchanged != 847 or three_unchanged != 487:
        raise ValueError(f"matched endpoint counts changed: one={one_unchanged}, three={three_unchanged}")

    maximum_errors = {
        "one_fixed": float(data.one_single.maximum_conservation_error.max()),
        "three_fixed": float(data.three_single[["maximum_cdc42_conservation_error", "maximum_rac_conservation_error", "maximum_rho_conservation_error"]].max().max()),
        "one_moving": max(float(frame.maximum_active_gtpase_conservation_error.max()) for (species, _), frame in data.moving.items() if species == "single_species"),
        "three_moving": max(float(frame[["maximum_cdc42_conservation_error", "maximum_rac_conservation_error", "maximum_rho_conservation_error"]].max().max()) for (species, _), frame in data.moving.items() if species == "three_species"),
    }
    if maximum_errors["one_fixed"] > 1e-7 or maximum_errors["three_fixed"] > 1e-6:
        raise ValueError(f"fixed-domain conservation limit exceeded: {maximum_errors}")
    if max(maximum_errors["one_moving"], maximum_errors["three_moving"]) > 1e-7:
        raise ValueError(f"moving-domain conservation limit exceeded: {maximum_errors}")
    return {
        "row_counts": expected_rows,
        "moving_rows_per_campaign": 240,
        "matched_endpoint_retention": {
            "single_species": {"unchanged": one_unchanged, "denominator": 867},
            "three_species": {"unchanged": three_unchanged, "denominator": 507},
        },
        "maximum_conservation_errors": maximum_errors,
    }


def validate_pngs(paths: Sequence[Path]) -> list[dict[str, object]]:
    from PIL import Image, ImageStat

    records = []
    for path in paths:
        with Image.open(path) as image:
            image.verify()
        with Image.open(path) as image:
            rgb = image.convert("RGB")
            width, height = rgb.size
            if width < 2400 or height < 1400:
                raise ValueError(f"{path.name} is unexpectedly small at {width}×{height}")
            thumbnail = rgb.copy()
            thumbnail.thumbnail((600, 600))
            extrema = ImageStat.Stat(thumbnail).extrema
            dynamic_range = max(high - low for low, high in extrema)
            if dynamic_range < 40:
                raise ValueError(f"{path.name} has insufficient visual dynamic range")
            records.append(
                {
                    "filename": path.name,
                    "width_px": width,
                    "height_px": height,
                    "size_bytes": path.stat().st_size,
                    "mode": rgb.mode,
                    "dynamic_range": int(dynamic_range),
                }
            )
    return records


def build_manifest(paths: Sequence[Path], validation: dict[str, object]) -> Path:
    records = validate_pngs(paths)
    payload = {
        "format_policy": "PNG only",
        "render_dpi": 300,
        "figures": records,
        "input_validation": validation,
        "representative_fixed_trajectories": {
            "source": "reproducibly rerun from recorded campaign definitions",
            "single_species": {
                "single": {"h": 1.0, "mu": 2.46875, "Lambda": 49.375},
                "dual": {"h": 1.0, "mu": 2.375, "Lambda": 25.0, "epsilon": 0.01},
                "reversal": {"h": 1.0, "mu": 2.1875, "Lambda": 71.875},
            },
            "three_species": {
                "single": {"n": 3.0, "kappa_R": 0.05333333333333334, "Lambda": 68.51601597035},
                "dual": {"n": 3.0, "kappa_R": 0.08, "Lambda": 84.3274042712, "epsilon": 0.01},
                "reversal": {"n": 4.0, "kappa_R": 0.04, "Lambda": 63.2455532034},
            },
        },
    }
    destination = HERE / "figure_manifest.json"
    destination.write_text(json.dumps(payload, indent=2) + "\n")
    return destination


def parse_selection(text: str) -> list[str]:
    available = ["1", "2", "3", "4", "5", "6", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]
    if text.lower() == "all":
        return available
    requested = [item.strip().upper() for item in text.split(",") if item.strip()]
    invalid = [item for item in requested if item not in available]
    if invalid:
        raise ValueError(f"unknown figures {invalid}; choose from {available}")
    return requested


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", default="all", help="comma-separated figure IDs (1..6,S1..S8) or all")
    parser.add_argument("--refresh-fixed-cache", action="store_true", help="rerun selected fixed-domain representatives")
    args = parser.parse_args(argv)
    if args.refresh_fixed_cache and CACHE.exists():
        for path in CACHE.glob("*.npz"):
            path.unlink()

    data = read_data()
    validation = validate_inputs(data)
    functions = {
        "1": figure_1,
        "2": figure_2,
        "3": figure_3,
        "4": figure_4,
        "5": figure_5,
        "6": figure_6,
        "S1": supplementary_1,
        "S2": supplementary_2,
        "S3": supplementary_3,
        "S4": supplementary_4,
        "S5": supplementary_5,
        "S6": supplementary_6,
        "S7": supplementary_7,
        "S8": supplementary_8,
    }
    selected = parse_selection(args.only)
    paths = []
    for identifier in selected:
        print(f"Rendering {identifier}...", flush=True)
        paths.append(functions[identifier](data))
        print(paths[-1], flush=True)
    manifest = build_manifest(paths, validation)
    print(manifest, flush=True)


if __name__ == "__main__":
    main()
