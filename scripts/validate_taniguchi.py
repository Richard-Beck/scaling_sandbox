#!/usr/bin/env python3
"""Run acceptance checks and render their diagnostic figures.

These are reproduction/validation cases only.  No cell-size experiment is
implemented here or in the production driver.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from taniguchi.fixed import FixedParameters, TriangularGrid, simulate_fixed  # noqa: E402
from taniguchi.metrics import directional_radius, fixed_diagnostics  # noqa: E402
from taniguchi.phase_field import simulate_phase_field  # noqa: E402

ASSETS = ROOT / "reports/assets/taniguchi/v1"


def _field(ax, result, index, field="V", vmax=None):
    values = getattr(result, field)[index]
    image = ax.pcolormesh(result.grid.x, result.grid.y,
                          np.where(result.grid.mask, values, np.nan),
                          shading="nearest", vmin=0, vmax=vmax, cmap="viridis")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    return image


def threshold_figure() -> dict:
    p = FixedParameters(K_k=5.7)
    low = simulate_fixed(30, parameters=p, stochastic=False, output_interval=.05,
                         initial_transfers=((0, 0, 0.5, 0.3),))
    high = simulate_fixed(30, parameters=p, stochastic=False, output_interval=.05,
                          initial_transfers=((0, 0, 1.0, 0.3),))
    mask = high.grid.mask
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    axes[0].plot(low.times, low.V[:, mask].max(1), label="subthreshold (0.5)")
    axes[0].plot(high.times, high.V[:, mask].max(1), label="suprathreshold (1.0)")
    axes[0].set(xlabel="time (10 s units)", ylabel="maximum PIP3",
                title="Excitable excursion and recovery")
    axes[0].legend(frameon=False)
    axes[1].plot(high.times, high.U[:, mask].mean(1), label="mean PIP2")
    axes[1].plot(high.times, high.V[:, mask].mean(1), label="mean PIP3")
    axes[1].axhline(1 / p.gamma, color="0.6", ls="--", lw=1, label="PIP2 rest")
    axes[1].set(xlabel="time (10 s units)", ylabel="domain mean",
                title="Refractory recovery")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(ASSETS / "acceptance_excitable_recovery.png", dpi=180)
    plt.close(fig)
    return {"subthreshold_peak_V": float(low.V.max()),
            "suprathreshold_peak_V": float(high.V.max()),
            "suprathreshold_final_mean_V": float(high.V[-1, mask].mean()),
            "suprathreshold_final_mean_U": float(high.U[-1, mask].mean())}


def stochastic_figure() -> dict:
    result = simulate_fixed(20, parameters=FixedParameters(K_k=5.7), seed=4,
                            output_interval=.1)
    mask = result.grid.mask
    active = (result.V[:, mask] > .5).mean(1)
    peak = int(np.argmax(active))
    candidates = np.unique(np.clip([peak - 25, peak - 12, peak, peak + 20],
                                   0, len(result.times) - 1))
    while len(candidates) < 4:
        candidates = np.unique(np.append(candidates, len(candidates) * 10))
    candidates = candidates[:4]
    fig, axes = plt.subplots(2, 4, figsize=(11, 5.5))
    for column, index in enumerate(candidates):
        _field(axes[0, column], result, index, "V", vmax=4)
        _field(axes[1, column], result, index, "U", vmax=1 / result.parameters.gamma)
        axes[0, column].set_title(f"t={result.times[index]:.1f}")
    axes[0, 0].set_ylabel("PIP3")
    axes[1, 0].set_ylabel("PIP2")
    fig.suptitle("Stochastic nucleation produces propagating, anticorrelated domains")
    fig.tight_layout()
    fig.savefig(ASSETS / "acceptance_stochastic_waves.png", dpi=180)
    plt.close(fig)
    return fixed_diagnostics(result)


def collision_spiral_figure() -> dict:
    p = FixedParameters(K_k=5.7)
    collision = simulate_fixed(10, parameters=p, stochastic=False,
                               output_interval=.1,
                               initial_transfers=((-2, 0, 1.0, .35),
                                                  (2, 0, 1.0, .35)))
    g = TriangularGrid.circle()
    U = np.zeros(g.mask.shape)
    U[g.mask] = 1 / p.gamma
    V = np.zeros_like(U)
    wave = g.mask & (g.x > -2.0) & (g.x < -1.2) & (g.y < .3)
    refractory = g.mask & (g.x < -1.2) & (g.y > .3)
    V[wave] = 2.0
    U[wave] -= 2.0
    U[refractory] = 1.2
    spiral = simulate_fixed(20, parameters=p, grid=g, stochastic=False,
                            output_interval=.1, initial_U=U, initial_V=V)
    fig, axes = plt.subplots(2, 4, figsize=(10, 5.2))
    for ax, time in zip(axes[0], (0, 2, 4, 8)):
        index = int(np.argmin(abs(collision.times - time)))
        _field(ax, collision, index, vmax=5)
        ax.set_title(f"t={collision.times[index]:.0f}")
    for ax, time in zip(axes[1], (3, 8, 12, 20)):
        index = int(np.argmin(abs(spiral.times - time)))
        _field(ax, spiral, index, vmax=5)
        ax.set_title(f"t={spiral.times[index]:.0f}")
    axes[0, 0].set_ylabel("collision")
    axes[1, 0].set_ylabel("rotating wave")
    fig.suptitle("Wave collision/annihilation and a rotating broken front")
    fig.tight_layout()
    fig.savefig(ASSETS / "acceptance_collision_rotating_wave.png", dpi=180)
    plt.close(fig)
    return {"collision_final_mean_V": float(collision.V[-1, g.mask].mean()),
            "rotating_wave_V_at_t20": float(spiral.V[-1, g.mask].max())}


def phase_and_convergence_figure() -> dict:
    phase = simulate_phase_field(.30, dx=.1, dt=8e-5, box_half_width=7,
                                 stochastic=False, output_interval=.025,
                                 initial_transfers=((4.5, 0, 1.0, .5),))
    right = np.asarray([directional_radius(z, phase.x, phase.y, "right")
                        for z in phase.phi])
    left = np.asarray([directional_radius(z, phase.x, phase.y, "left")
                       for z in phase.phi])

    cases = {
        "published": simulate_fixed(8, stochastic=False, output_interval=.05,
                                    initial_transfers=((0, 0, 1.0, .3),)),
        "half dt": simulate_fixed(8, dt=2.5e-4, stochastic=False,
                                  output_interval=.05,
                                  initial_transfers=((0, 0, 1.0, .3),)),
        "finer dx": simulate_fixed(8, grid=TriangularGrid.circle(dx=.125),
                                   stochastic=False, output_interval=.05,
                                   initial_transfers=((0, 0, 1.0, .3),)),
    }
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7))
    axes[0].imshow(phase.V[-1], origin="lower",
                   extent=[phase.x[0], phase.x[-1], phase.y[0], phase.y[-1]],
                   vmin=0, vmax=1.5, cmap="viridis")
    axes[0].contour(phase.x, phase.y, phase.phi[0], [.5], colors="white",
                    linestyles="--", linewidths=1)
    axes[0].contour(phase.x, phase.y, phase.phi[-1], [.5], colors="black",
                    linewidths=1.5)
    axes[0].set(title="PIP3-driven protrusion", aspect="equal",
                xlabel="x", ylabel="y")
    axes[1].plot(phase.times, right, label="PIP3-rich right edge")
    axes[1].plot(phase.times, -left, label="opposite edge magnitude")
    axes[1].set(title="Directional boundary response", xlabel="time",
                ylabel="axis radius")
    axes[1].legend(frameon=False, fontsize=8)
    for label, result in cases.items():
        axes[2].plot(result.times, result.V[:, result.grid.mask].mean(1), label=label)
    axes[2].set(title="Time/space refinement", xlabel="time",
                ylabel="mean PIP3")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(ASSETS / "acceptance_phase_and_refinement.png", dpi=180)
    plt.close(fig)
    peaks = {label: float(result.V[:, result.grid.mask].mean(1).max())
             for label, result in cases.items()}
    return {"initial_right_radius": float(right[0]),
            "final_right_radius": float(right[-1]),
            "initial_left_radius_magnitude": float(-left[0]),
            "final_left_radius_magnitude": float(-left[-1]),
            "mean_V_peaks_by_resolution": peaks}


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    results = {
        "excitable_recovery": threshold_figure(),
        "stochastic_waves": stochastic_figure(),
        "collision_and_rotation": collision_spiral_figure(),
        "phase_and_refinement": phase_and_convergence_figure(),
    }
    (ASSETS / "acceptance_metrics.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
