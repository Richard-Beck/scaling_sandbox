#!/usr/bin/env python3
"""Map the endpoint-cue asymmetry required to resolve Eroume Rac fronts."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from run_eroume_endpoint_competition import (
    CUE,
    CUE_STOP,
    OUTPUT,
    competition_series,
    run_condition,
)
from migration1d import EroumeTwoEdgeCue

RATIOS = (0.0, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025,
          0.05, 0.10, 0.20, 0.40, 0.79, 1.0)
DURATION = 150.0


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    left_increment = CUE["left_level"]-CUE["basal_level"]
    records = []
    trajectories = {}
    conditions = {
        "1x": {"length": 10.0, "nodes": 61, "linear_scale": 1.0},
        "wgd": {
            "length": 10.0*2.0**(1.0/3.0), "nodes": 77,
            "linear_scale": 2.0**(1.0/3.0),
        },
    }
    for scale, condition in conditions.items():
        for ratio in RATIOS:
            cue = EroumeTwoEdgeCue(
                left_level=CUE["left_level"],
                right_level=CUE["basal_level"]+ratio*left_increment,
                basal_level=CUE["basal_level"],
                width_fraction=CUE["width_fraction"],
                start=0.0, stop=CUE_STOP,
            )
            result = run_condition(
                **condition, cue=cue, duration=DURATION,
            )
            metrics = competition_series(result)
            trajectories[(scale, ratio)] = (result, metrics)
            records.append({
                "scale": scale,
                "right_left_increment_ratio": ratio,
                "right_level": cue.right_level,
                "final_left_mass_fraction": float(metrics["left_mass_fraction"][-1]),
                "final_loser_winner_mass_ratio": float(
                    metrics["loser_winner_mass_ratio"][-1]
                ),
                "final_loser_winner_peak_ratio": float(
                    metrics["loser_winner_peak_ratio"][-1]
                ),
                "final_loser_winner_endpoint_ratio": float(
                    metrics["loser_winner_endpoint_ratio"][-1]
                ),
                "final_signed_endpoint_contrast": float(
                    metrics["signed_endpoint_contrast"][-1]
                ),
                "net_displacement": float(result.center[-1]-result.center[0]),
            })

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"1x": "#0072B2", "wgd": "#D55E00"}
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 8.0), sharex=True,
                             constrained_layout=True)
    for scale in ("1x", "wgd"):
        subset = [item for item in records if item["scale"] == scale]
        ratios = np.asarray([item["right_left_increment_ratio"] for item in subset])
        axes[0].plot(
            ratios, [item["final_signed_endpoint_contrast"] for item in subset],
            "o-", color=colors[scale], label=scale.upper(),
        )
        axes[1].plot(
            ratios, [item["final_loser_winner_endpoint_ratio"] for item in subset],
            "o-", color=colors[scale], label=scale.upper(),
        )
    axes[0].axhline(0.0, color="0.4", ls=":")
    axes[0].set(ylabel="signed endpoint Rac contrast", ylim=(-1.02, 0.05))
    axes[0].legend(frameon=False)
    axes[1].axhline(0.2, color="0.4", ls=":", label="20% two-pole threshold")
    axes[1].set(
        xlabel="right/left endpoint Ic increment", ylabel="loser/winner endpoint Rac",
        ylim=(-0.02, 1.02),
    )
    axes[1].legend(frameon=False)
    fig.suptitle("Eroumé–Marée two-front competition threshold at t=150")
    fig.savefig(OUTPUT/"04_cue_asymmetry_screen.png", dpi=220)
    plt.close(fig)

    selected = (("1x", 0.01), ("wgd", 0.01), ("1x", 0.025), ("wgd", 0.025))
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.5), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, (scale, ratio) in zip(axes.flat, selected):
        result, metrics = trajectories[(scale, ratio)]
        time_grid = np.broadcast_to(result.times[:, None], result.position.shape)
        image = ax.pcolormesh(
            result.position, time_grid, result.membrane[:, 3], shading="auto",
            cmap="magma", rasterized=True,
        )
        ax.plot(result.rear, result.times, color="white", lw=0.7)
        ax.plot(result.front, result.times, color="white", lw=0.7)
        ax.axhline(CUE_STOP, color="#54d6ff", ls="--", lw=0.8)
        ax.set(
            title=f"{scale.upper()} · right/left increment={ratio:g}\n"
                  f"endpoint loser/winner={metrics['loser_winner_endpoint_ratio'][-1]:.3f}",
            xlabel="laboratory position", ylabel="time",
        )
        plt.colorbar(image, ax=ax, pad=0.01, label="active Rac")
    fig.suptitle("Endpoint competition across cue asymmetry")
    fig.savefig(OUTPUT/"05_cue_asymmetry_kymographs.png", dpi=220)
    plt.close(fig)

    (OUTPUT/"cue_asymmetry_screen.json").write_text(
        json.dumps({"duration": DURATION, "records": records}, indent=2)+"\n"
    )
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
