#!/usr/bin/env python3
"""Run a motile 1x/WGD endpoint-competition assay with Eroume chemistry."""

from __future__ import annotations

import json

import numpy as np

from ._paths import output_dir

from migration1d import (  # noqa: E402
    EroumeChemistry,
    EroumeTwoEdgeCue,
    MechanicsParameters,
    MigrationModel,
    RacRhoForce,
    SimulationConfig,
    run_simulation,
)
from migration1d.model import control_volumes  # noqa: E402

OUTPUT = output_dir("three_species/testing/competing_stimuli")
WGD_SCALE = 2.0**(1.0/3.0)
CUE_STOP = 10.0
DURATION = 300.0
CUE = {
    "left_level": 5.10,
    # A 1% right/left increment lies on opposite sides of the 1x and WGD
    # second-pole nucleation thresholds in the local asymmetry screen.
    "right_level": 2.9715,
    "basal_level": 2.95,
    "width_fraction": 0.064,
}


def make_model(
    length: float,
    linear_scale: float,
    cue: EroumeTwoEdgeCue | None = None,
) -> MigrationModel:
    cue = cue or EroumeTwoEdgeCue(**CUE, start=0.0, stop=CUE_STOP)
    # reference_length=10 fixes the physical diffusion coefficients under the
    # WGD length change. Using reference_length=length would normalize away
    # the chemical size perturbation.
    chemistry = EroumeChemistry(reference_length=10.0, cue=cue)
    force = RacRhoForce(
        rac_field="Rac_act", rho_field="Rho_act",
        rac_gain=0.08, rho_gain=0.08,
        rac_half_level=3.0, rho_half_level=1.25, hill_exponent=4.0,
    )
    return MigrationModel(
        mechanics=MechanicsParameters(
            reference_length=length,
            spring_stiffness=1.0/linear_scale,
            dashpot_viscosity=2.0/linear_scale,
            left_friction=1.0, right_friction=1.0,
        ),
        chemistry=chemistry, left_force=force, right_force=force,
    )


def run_condition(
    length: float,
    nodes: int,
    linear_scale: float,
    cue: EroumeTwoEdgeCue | None = None,
    duration: float = DURATION,
):
    return run_simulation(
        make_model(length, linear_scale, cue),
        SimulationConfig(
            nodes=nodes, duration=duration, output_interval=1.0, max_step=0.2,
            noise_amplitude=0.0, mechanics_activation_time=CUE_STOP,
        ),
    )


def competition_series(result, field: int = 3) -> dict[str, np.ndarray]:
    count = len(result.times)
    left_mass = np.empty(count)
    right_mass = np.empty(count)
    left_peak = np.empty(count)
    right_peak = np.empty(count)
    signed_endpoint_contrast = np.empty(count)
    loser_winner_endpoint_ratio = np.empty(count)
    for index, (position, raw_profile) in enumerate(
        zip(result.position, result.membrane[:, field])
    ):
        profile = np.maximum(raw_profile, 0.0)
        weights = control_volumes(position)
        midpoint = len(weights)//2
        left_weights = weights.copy()
        right_weights = weights.copy()
        left_weights[midpoint+1:] = 0.0
        right_weights[:midpoint] = 0.0
        left_weights[midpoint] *= 0.5
        right_weights[midpoint] *= 0.5
        excess = np.maximum(profile-profile.min(), 0.0)
        left_mass[index] = excess@left_weights
        right_mass[index] = excess@right_weights
        band = max(2, len(profile)//5)
        left_peak[index] = profile[:band].max()
        right_peak[index] = profile[-band:].max()
        signed_endpoint_contrast[index] = (
            (profile[-1]-profile[0])
            / max(abs(profile[-1])+abs(profile[0]), 1.0e-12)
        )
        loser_winner_endpoint_ratio[index] = (
            min(profile[0], profile[-1])/max(profile[0], profile[-1], 1.0e-12)
        )
    loser = np.minimum(left_mass, right_mass)
    winner = np.maximum(left_mass, right_mass)
    return {
        "left_mass": left_mass,
        "right_mass": right_mass,
        "left_peak": left_peak,
        "right_peak": right_peak,
        "loser_winner_mass_ratio": np.divide(
            loser, winner, out=np.zeros_like(loser), where=winner > 1.0e-14,
        ),
        "loser_winner_peak_ratio": np.divide(
            np.minimum(left_peak, right_peak), np.maximum(left_peak, right_peak),
            out=np.zeros_like(left_peak),
            where=np.maximum(left_peak, right_peak) > 1.0e-14,
        ),
        "left_mass_fraction": np.divide(
            left_mass, left_mass+right_mass, out=np.full_like(left_mass, 0.5),
            where=(left_mass+right_mass) > 1.0e-14,
        ),
        "signed_endpoint_contrast": signed_endpoint_contrast,
        "loser_winner_endpoint_ratio": loser_winner_endpoint_ratio,
    }


def resolution_time(result, series, threshold: float = 0.2) -> float | None:
    """First post-cue time after which the endpoint loser stays resolved."""
    ratio = series["loser_winner_endpoint_ratio"]
    for index in np.flatnonzero((result.times >= CUE_STOP) & (ratio <= threshold)):
        if np.all(ratio[index:] <= threshold):
            return float(result.times[index])
    return None


def render(results, metrics) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"1x": "#0072B2", "wgd": "#D55E00"}
    titles = {"1x": "1× motile", "wgd": "WGD-null motile"}
    vmax = max(float(result.membrane[:, 3].max()) for result in results.values())
    xmin = min(float(result.rear.min()) for result in results.values())
    xmax = max(float(result.front.max()) for result in results.values())
    margin = 0.025*(xmax-xmin)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.5), sharex=True, sharey=True,
                             constrained_layout=True)
    image = None
    for ax, key in zip(axes, ("1x", "wgd")):
        result = results[key]
        time_grid = np.broadcast_to(result.times[:, None], result.position.shape)
        image = ax.pcolormesh(
            result.position, time_grid, result.membrane[:, 3], shading="auto",
            cmap="magma", vmin=0.0, vmax=vmax, rasterized=True,
        )
        ax.plot(result.rear, result.times, color="white", lw=0.7)
        ax.plot(result.front, result.times, color="white", lw=0.7)
        ax.axhline(CUE_STOP, color="#54d6ff", ls="--", lw=0.9)
        endpoint_ratio = metrics[key]["loser_winner_endpoint_ratio"][-1]
        if endpoint_ratio > 0.2:
            outcome = "two persistent poles"
        else:
            outcome = (
                "left wins" if metrics[key]["left_peak"][-1]
                > metrics[key]["right_peak"][-1] else "right wins"
            )
        resolved = resolution_time(result, metrics[key])
        timing = "unresolved" if resolved is None else f"resolved at t={resolved:g}"
        ax.text(
            0.98, 0.98,
            f"{outcome}\n"
            f"endpoint loser/winner={endpoint_ratio:.3g} · {timing}\n"
            f"Δx={result.center[-1]-result.center[0]:+.2f}",
            transform=ax.transAxes, ha="right", va="top", color="white",
            fontsize=8,
            bbox={"facecolor": "black", "alpha": 0.42, "edgecolor": "none"},
        )
        ax.set(title=titles[key], xlabel="laboratory position",
               xlim=(xmin-margin, xmax+margin))
    axes[0].set_ylabel("time")
    fig.colorbar(image, ax=axes, pad=0.02, label="active Rac")
    fig.suptitle(
        "Eroumé–Marée endpoint competition with motility\n"
        f"left Ic={CUE['left_level']:.4g} · right Ic={CUE['right_level']:.4g} · "
        "cyan dash: cues removed and mechanics released"
    )
    fig.savefig(OUTPUT/"00_eroume_1x_vs_wgd_motile_rac_kymographs.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(16.0, 10.0), sharey=True,
                             constrained_layout=True)
    for row, key in enumerate(("1x", "wgd")):
        result = results[key]
        time_grid = np.broadcast_to(result.times[:, None], result.position.shape)
        for ax, field, label in zip(
            axes[row], (1, 3, 5), ("Cdc42", "Rac", "Rho")
        ):
            values = result.membrane[:, field]
            image = ax.pcolormesh(result.position, time_grid, values,
                                  shading="auto", cmap="magma", rasterized=True)
            ax.plot(result.rear, result.times, color="white", lw=0.7)
            ax.plot(result.front, result.times, color="white", lw=0.7)
            ax.axhline(CUE_STOP, color="#54d6ff", ls="--", lw=0.8)
            ax.set(title=f"{titles[key]} · active {label}",
                   xlabel="laboratory position")
            plt.colorbar(image, ax=ax, pad=0.01)
        axes[row, 0].set_ylabel("time")
    fig.suptitle("Competing Cdc42 inputs propagate through Rac–Rho polarity")
    fig.savefig(OUTPUT/"01_active_gtpase_kymographs.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(4, 1, figsize=(9.5, 11.0), sharex=True,
                             constrained_layout=True)
    for key in ("1x", "wgd"):
        result = results[key]
        series = metrics[key]
        axes[0].plot(result.times, result.center-result.center[0],
                     color=colors[key], label=key.upper())
        axes[1].plot(result.times, 100*(result.length/result.length[0]-1),
                     color=colors[key])
        axes[2].plot(result.times, series["left_mass_fraction"], color=colors[key])
        axes[3].plot(result.times, series["loser_winner_endpoint_ratio"],
                     color=colors[key])
    axes[0].set_ylabel("center displacement")
    axes[0].legend(frameon=False)
    axes[1].set_ylabel("length change (%)")
    axes[2].axhline(0.5, color="0.4", ls=":")
    axes[2].set_ylabel("left Rac mass fraction")
    axes[3].axhline(0.2, color="0.4", ls=":")
    axes[3].set(xlabel="time", ylabel="loser/winner endpoint Rac")
    for ax in axes:
        ax.axvspan(0.0, CUE_STOP, color="#54d6ff", alpha=0.12)
    fig.suptitle("Competition, migration, and deformation")
    fig.savefig(OUTPUT/"02_competition_and_migration_metrics.png", dpi=220)
    plt.close(fig)

    selected_times = (CUE_STOP, 30.0, 100.0, DURATION)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, time in zip(axes.flat, selected_times):
        for key in ("1x", "wgd"):
            result = results[key]
            index = int(np.argmin(abs(result.times-time)))
            relative = ((result.position[index]-result.rear[index])
                        / result.length[index])
            ax.plot(relative, result.membrane[index, 3], color=colors[key],
                    label=key.upper())
        ax.set(title=f"t={time:g}", xlabel="relative position", ylabel="active Rac")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Rac competition in relative cell coordinates")
    fig.savefig(OUTPUT/"03_relative_coordinate_profiles.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    conditions = {
        "1x": {"length": 10.0, "nodes": 101, "linear_scale": 1.0},
        "wgd": {
            "length": 10.0*WGD_SCALE, "nodes": 127,
            "linear_scale": WGD_SCALE,
        },
    }
    results = {key: run_condition(**values) for key, values in conditions.items()}
    metrics = {key: competition_series(result) for key, result in results.items()}
    for key, result in results.items():
        np.savez_compressed(
            OUTPUT/f"{key}_motile_trajectory.npz",
            times=result.times, position=result.position, velocity=result.velocity,
            membrane=result.membrane, cytoplasm=result.cytoplasm,
            applied_edge_force=result.applied_edge_force,
            spring_force=result.spring_force, dashpot_force=result.dashpot_force,
            species=np.asarray(result.model.chemistry.species), **metrics[key],
        )
    render(results, metrics)

    summary = {
        "experiment": "Eroume nine-field endpoint competition, 1x versus WGD",
        "cue": {**CUE, "start": 0.0, "stop": CUE_STOP,
                "interpretation": "competing endpoint Ic levels; uniform Ic=2.95 afterward"},
        "mechanical_coupling": {
            "rac_gain": 0.08, "rho_gain": 0.08,
            "mechanics_release_time": CUE_STOP,
        },
        "wgd_null": {
            "length_scale": WGD_SCALE,
            "concentrations": "preserved",
            "diffusion_coefficients": "preserved in normalized physical units",
            "cue_width": "preserved as a fraction of cell length",
            "spring_and_dashpot": "scaled inversely with length",
        },
        "conditions": {},
    }
    for key, result in results.items():
        series = metrics[key]
        totals, _ = result.compartment_totals()
        pool_drifts = {}
        for label, pair in (("cdc42", (0, 1)), ("rac", (2, 3)), ("rho", (4, 5))):
            pool = totals[:, list(pair)].sum(axis=1)
            pool_drifts[label] = float(np.max(abs(pool/pool[0]-1.0)))
        summary["conditions"][key] = {
            **conditions[key],
            "terminal_outcome": (
                "two_persistent_poles"
                if series["loser_winner_endpoint_ratio"][-1] > 0.2
                else (
                    "left_wins" if series["left_peak"][-1] > series["right_peak"][-1]
                    else "right_wins"
                )
            ),
            "final_left_rac_mass_fraction": float(series["left_mass_fraction"][-1]),
            "final_loser_winner_rac_mass_ratio": float(
                series["loser_winner_mass_ratio"][-1]
            ),
            "final_loser_winner_rac_peak_ratio": float(
                series["loser_winner_peak_ratio"][-1]
            ),
            "final_loser_winner_rac_endpoint_ratio": float(
                series["loser_winner_endpoint_ratio"][-1]
            ),
            "final_signed_rac_endpoint_contrast": float(
                series["signed_endpoint_contrast"][-1]
            ),
            "resolution_time": resolution_time(result, series),
            "net_displacement": float(result.center[-1]-result.center[0]),
            "minimum_length": float(result.length.min()),
            "maximum_length": float(result.length.max()),
            "gtpase_pair_relative_drift": pool_drifts,
        }
    (OUTPUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
