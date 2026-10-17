#!/usr/bin/env python3
"""Run matched 1x/WGD six-field Holmes endpoint-migration experiments."""

from __future__ import annotations

from dataclasses import replace
import json

import numpy as np
from scipy.integrate import solve_ivp

from ._paths import output_dir

from migration1d import (  # noqa: E402
    HolmesChemistry,
    HolmesParameters,
    MechanicsParameters,
    MigrationModel,
    RacRhoForce,
    SimulationConfig,
    run_simulation,
)
from wgd_v1.polarity_models import (  # noqa: E402
    holmes_accessibility,
    holmes_inactive_diffusion,
)

OUTPUT = output_dir("three_species/testing/single_stimulus")
LINEAR_WGD_SCALE = 2.0 ** (1.0 / 3.0)
BASELINE_INACTIVE_FRACTION = 0.5
RELEASE_TIME = 70.0
RAC_GAIN = 0.08
RHO_GAIN = 0.08


class RelativeGradientPulse:
    """Holmes Rac-GEF cue fixed in relative cell coordinates."""

    def __init__(self, start=50.0, stop=70.0, mean=0.3, half_range=0.1):
        self.start = start
        self.stop = stop
        self.mean = mean
        self.half_range = half_range

    def __call__(self, time: float, position: np.ndarray) -> np.ndarray:
        if not self.start <= time < self.stop:
            return np.zeros_like(position)
        relative = (position - position[0]) / (position[-1] - position[0])
        return self.mean + self.half_range * (2.0 * relative - 1.0)


def prestimulus_equilibrium(parameters: HolmesParameters) -> tuple[float, ...]:
    """Relax the homogeneous circuit at zero external Rac-GEF signal."""
    totals = parameters.totals
    initial = np.concatenate((np.full(3, 0.1), totals - 0.1))

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        cdc42, rac, rho, cdc42_i, rac_i, rho_i = state
        accessibility = parameters.inactive_accessibility
        rates = np.asarray((
            parameters.cdc42_activation
            / (1.0 + (rho/parameters.rho_inhibits_cdc42_half_level)
               ** parameters.hill_exponent)
            * accessibility*cdc42_i/parameters.cdc42_total - cdc42,
            (parameters.rac_basal_activation
             + parameters.cdc42_to_rac_feedback*cdc42)
            * accessibility*rac_i/parameters.rac_total - rac,
            parameters.rho_activation
            / (1.0 + (rac/parameters.rac_inhibits_rho_half_level)
               ** parameters.hill_exponent)
            * accessibility*rho_i/parameters.rho_total - rho,
        ))
        return np.concatenate((rates, -rates))

    result = solve_ivp(
        rhs, (0.0, 300.0), initial, method="BDF", rtol=1.0e-10, atol=1.0e-12
    )
    if not result.success:
        raise RuntimeError(f"Holmes homogeneous relaxation failed: {result.message}")
    return tuple(float(value) for value in result.y[:3, -1])


def make_model(volume_ratio: float) -> MigrationModel:
    linear_scale = volume_ratio ** (1.0 / 3.0)
    parameters = replace(
        HolmesParameters(),
        # The complete signal is supplied only during the pulse.
        signal_baseline=0.0,
        inactive_accessibility=holmes_accessibility(
            volume_ratio, BASELINE_INACTIVE_FRACTION
        ),
        inactive_diffusivity=holmes_inactive_diffusion(
            volume_ratio, BASELINE_INACTIVE_FRACTION
        ),
    )
    parameters = replace(
        parameters, initial_active=prestimulus_equilibrium(parameters)
    )
    chemistry = HolmesChemistry(parameters, cue=RelativeGradientPulse())
    force = RacRhoForce(
        rac_gain=RAC_GAIN,
        rho_gain=RHO_GAIN,
        # Balance the two signed responses in the prestimulus equilibrium.
        rac_half_level=parameters.initial_active[1],
        rho_half_level=parameters.initial_active[2],
        hill_exponent=4.0,
    )
    return MigrationModel(
        mechanics=MechanicsParameters(
            reference_length=10.0 * linear_scale,
            spring_stiffness=1.0 / linear_scale,
            dashpot_viscosity=2.0 / linear_scale,
            left_friction=1.0,
            right_friction=1.0,
        ),
        chemistry=chemistry,
        left_force=force,
        right_force=force,
    )


def run_condition(volume_ratio: float, nodes: int):
    return run_simulation(
        make_model(volume_ratio),
        SimulationConfig(
            nodes=nodes,
            duration=300.0,
            output_interval=1.0,
            max_step=0.25,
            noise_amplitude=0.0,
            mechanics_activation_time=RELEASE_TIME,
        ),
    )


def hill(values: np.ndarray, half: float, exponent: float = 4.0) -> np.ndarray:
    positive = np.maximum(values, 0.0)
    power = positive**exponent
    return power / (half**exponent + power)


def signed_drive(result) -> np.ndarray:
    p = result.model.chemistry.parameters
    rac = hill(result.membrane[:, 1], p.initial_active[1])
    rho = hill(result.membrane[:, 2], p.initial_active[2])
    return RAC_GAIN * rac - RHO_GAIN * rho


def contrast(values: np.ndarray) -> np.ndarray:
    left, right = values[:, 0], values[:, -1]
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1.0e-12)


def render_kymographs(results: dict, species: int, name: str, label: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vmax = max(float(result.membrane[:, species].max()) for result in results.values())
    xmin = min(float(result.rear.min()) for result in results.values())
    xmax = max(float(result.front.max()) for result in results.values())
    margin = 0.03 * (xmax - xmin)
    fig, axes = plt.subplots(
        1, 2, figsize=(12.5, 6.5), sharex=True, sharey=True,
        constrained_layout=True,
    )
    image = None
    for ax, key, title in zip(axes, ("1x", "wgd"), ("1× motile", "WGD motile")):
        result = results[key]
        time_grid = np.broadcast_to(result.times[:, None], result.position.shape)
        image = ax.pcolormesh(
            result.position, time_grid, result.membrane[:, species],
            shading="auto", cmap="magma", vmin=0.0, vmax=vmax,
            rasterized=True,
        )
        ax.plot(result.rear, result.times, color="white", linewidth=0.75)
        ax.plot(result.front, result.times, color="white", linewidth=0.75)
        ax.axhline(50.0, color="#54d6ff", linestyle="--", linewidth=0.8)
        ax.axhline(70.0, color="#54d6ff", linestyle="--", linewidth=0.8)
        ax.text(
            0.98, 0.98,
            f"Δx = {result.center[-1] - result.center[0]:+.3f}\n"
            f"L(300) = {result.length[-1]:.3f}",
            transform=ax.transAxes, ha="right", va="top", color="white",
            fontsize=9,
            bbox={"facecolor": "black", "alpha": 0.4, "edgecolor": "none"},
        )
        ax.set(title=title, xlabel="laboratory position",
               xlim=(xmin - margin, xmax + margin))
    axes[0].set_ylabel("time")
    colorbar = fig.colorbar(image, ax=axes, pad=0.02)
    colorbar.set_label(f"active {label}")
    fig.suptitle(
        f"Holmes six-field endpoint migration: active {label}\n"
        "white curves: material boundaries · cyan dashes: Rac-GEF cue window",
        fontsize=14,
    )
    fig.savefig(OUTPUT / f"{name}_holmes_1x_vs_wgd_{label.lower()}_kymographs.png",
                dpi=220)
    plt.close(fig)


def render_diagnostics(results: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"1x": "#0072B2", "wgd": "#D55E00"}
    fig, axes = plt.subplots(4, 1, figsize=(9.5, 11.0), sharex=True,
                             constrained_layout=True)
    for key, result in results.items():
        axes[0].plot(result.times, result.center - result.center[0],
                     color=colors[key], label=key.upper())
        axes[1].plot(result.times, 100.0 * (result.length/result.length[0] - 1.0),
                     color=colors[key])
        axes[2].plot(result.times, contrast(result.membrane[:, 1]),
                     color=colors[key])
        drive = signed_drive(result)
        axes[3].plot(result.times, drive[:, -1] - drive[:, 0],
                     color=colors[key])
    for ax in axes:
        ax.axvspan(50.0, 70.0, color="#54d6ff", alpha=0.12)
        ax.axvline(70.0, color="0.45", linestyle="--", linewidth=0.8)
    axes[0].set_ylabel("center displacement")
    axes[0].legend(frameon=False)
    axes[1].set_ylabel("length change (%)")
    axes[2].set_ylabel("Rac front/back contrast")
    axes[3].set(xlabel="time", ylabel="front − rear signed drive")
    fig.suptitle("Holmes motility, deformation, and endpoint polarity")
    fig.savefig(OUTPUT / "02_motility_and_polarity_metrics.png", dpi=220)
    plt.close(fig)

    selected_times = (50.0, 70.0, 150.0, 300.0)
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), sharex=True,
                             constrained_layout=True)
    styles = {"1x": "-", "wgd": "--"}
    species = ((0, "Cdc42", "#009E73"), (1, "Rac", "#CC79A7"),
               (2, "Rho", "#E69F00"))
    for ax, time in zip(axes.flat, selected_times):
        for key, result in results.items():
            index = int(np.argmin(np.abs(result.times - time)))
            relative = ((result.position[index] - result.rear[index])
                        / result.length[index])
            for field, label, color in species:
                ax.plot(relative, result.membrane[index, field],
                        color=color, linestyle=styles[key],
                        label=f"{label} {key.upper()}")
        ax.set(title=f"t = {time:g}", xlabel="relative position",
               ylabel="active concentration")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=3, frameon=False)
    fig.suptitle("Active GTPase profiles in relative cell coordinates")
    fig.savefig(OUTPUT / "03_relative_coordinate_profiles.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    conditions = {
        "1x": {"volume_ratio": 1.0, "nodes": 101},
        "wgd": {"volume_ratio": 2.0, "nodes": 127},
    }
    results = {key: run_condition(**values) for key, values in conditions.items()}
    for key, result in results.items():
        np.savez_compressed(
            OUTPUT / f"{key}_motile_trajectory.npz",
            times=result.times,
            position=result.position,
            velocity=result.velocity,
            membrane=result.membrane,
            cytoplasm=result.cytoplasm,
            signed_drive=signed_drive(result),
            applied_edge_force=result.applied_edge_force,
            spring_force=result.spring_force,
            dashpot_force=result.dashpot_force,
        )
    render_kymographs(results, 1, "00", "Rac")
    render_kymographs(results, 2, "01", "Rho")
    render_diagnostics(results)

    summary = {
        "experiment": "Holmes six-field 1x/WGD endpoint migration; cue 50-70 s; mechanics released at 70 s",
        "chemical_cue": {
            "target": "Rac activation",
            "coordinates": "relative material coordinate",
            "mean": 0.3,
            "half_range": 0.1,
            "start": 50.0,
            "stop": 70.0,
        },
        "mechanical_coupling": {
            "signed_outward_drive": "0.08*h(Rac)-0.08*h(Rho)",
            "hill_exponent": 4.0,
            "half_levels": "condition-specific prestimulus active equilibrium",
            "mechanics_release_time": RELEASE_TIME,
            "baseline_spring_stiffness": 1.0,
            "baseline_dashpot_viscosity": 2.0,
            "endpoint_friction": 1.0,
        },
        "wgd_null": {
            "length": "x 2^(1/3)",
            "concentrations_and_gtpase_totals": "preserved",
            "native_line_amount": "x 2^(1/3)",
            "inactive_accessibility": "Holmes volume-to-membrane scaling",
            "effective_inactive_diffusion": "Holmes fast-exchange scaling",
            "spring_and_dashpot": "x 2^(-1/3)",
            "endpoint_friction_and_force_mapping": "preserved",
        },
        "conditions": {},
    }
    for key, result in results.items():
        membrane_total, inactive_total = result.compartment_totals()
        total = membrane_total + inactive_total
        p = result.model.chemistry.parameters
        summary["conditions"][key] = {
            **conditions[key],
            "reference_length": result.model.mechanics.reference_length,
            "inactive_accessibility": p.inactive_accessibility,
            "inactive_diffusivity": p.inactive_diffusivity,
            "prestimulus_active": list(p.initial_active),
            "net_displacement": float(result.center[-1] - result.center[0]),
            "minimum_length": float(result.length.min()),
            "maximum_length": float(result.length.max()),
            "final_length": float(result.length[-1]),
            "final_rac_contrast": float(contrast(result.membrane[:, 1])[-1]),
            "final_rho_contrast": float(contrast(result.membrane[:, 2])[-1]),
            "maximum_relative_material_drift": float(
                np.max(np.abs(total / total[0] - 1.0))
            ),
        }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
