#!/usr/bin/env python3
"""Plot the focused Giese abundance/gamma pilot."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "reports/output/giese_wp/boundary_pilot"
FIGURES = OUTPUT / "figures"
MASSES = (0.92, 0.94, 0.96, 0.98, 1.00)
GAMMAS = (0.75, 1.00, 1.25)


def load(task: int) -> dict:
    with np.load(OUTPUT / "checkpoints" / f"task_{task:03d}.npz") as raw:
        return {key: raw[key] for key in raw.files}


def load_selected(task: int) -> dict:
    confirmation = (ROOT / "reports/output/giese_wp/boundary_confirmation/checkpoints"
                    / f"task_{task:03d}.npz")
    if confirmation.exists():
        with np.load(confirmation) as raw:
            return {key: raw[key] for key in raw.files}
    return load(task)


def summary(task: int) -> dict:
    return json.loads(str(load(task)["summary_json"]))


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, diameter_index, title in zip(
        axes, (0, 1), ("1×: 10 µm", "WGD: 12.6 µm")
    ):
        values = np.zeros((len(GAMMAS), len(MASSES)))
        labels = np.empty(values.shape, dtype=object)
        for gi, gamma in enumerate(GAMMAS):
            for mi, mass in enumerate(MASSES):
                row = summary(gi * 10 + mi * 2 + diameter_index)
                if np.isfinite(row["winner_selection_time_s"]):
                    values[gi, mi] = row["winner_selection_time_s"]
                    labels[gi, mi] = f'{row["winner_selection_time_s"]:.0f}s'
                elif row["winner_selection_right_censored"]:
                    values[gi, mi] = 500
                    labels[gi, mi] = ">500"
                else:
                    values[gi, mi] = 0
                    labels[gi, mi] = "no pol."
        artist = ax.imshow(values, origin="lower", aspect="auto", cmap="viridis",
                           vmin=0, vmax=500)
        ax.set_xticks(range(len(MASSES)), MASSES)
        ax.set_yticks(range(len(GAMMAS)), GAMMAS)
        ax.set_xlabel("Total-abundance scale")
        ax.set_ylabel(r"Feedback strength $\gamma$")
        ax.set_title(title)
        for gi in range(len(GAMMAS)):
            for mi in range(len(MASSES)):
                ax.text(mi, gi, labels[gi, mi], ha="center", va="center",
                        color="white", fontsize=9)
    fig.colorbar(artist, ax=axes, label="Winner-selection time / censoring bound (s)")
    fig.suptitle("Focused two-site competition pilot")
    fig.savefig(FIGURES / "01_boundary_outcomes.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for task, label, color in ((12, "1×: 10 µm", "#0072B2"),
                               (13, "WGD: 12.6 µm", "#D55E00")):
        data = load_selected(task)
        time = data["times"]
        axes[0].plot(time, data["competition_mass_bias"], label=label, color=color)
        axes[1].step(time, data["metric_domains"], where="post", label=label,
                     color=color)
    axes[0].set(xlabel="Time (s)", ylabel="Strong-vs-weak membrane mass bias",
                title=r"$\gamma=1$, abundance = 0.94")
    axes[1].set(xlabel="Time (s)", ylabel="Domains above mean", ylim=(-0.1, 2.2),
                title="Winner selection")
    axes[1].set_yticks((0, 1, 2))
    for ax in axes:
        ax.axvspan(0, 10, color="0.8", alpha=0.5, label="stimulus" if ax is axes[0] else None)
        ax.legend(frameon=False)
    fig.savefig(FIGURES / "02_selected_split_trajectories.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
