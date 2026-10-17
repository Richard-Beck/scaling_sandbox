#!/usr/bin/env python3
"""Challenge Holmes and Eroume/Maree polarity with transient and persistent cues."""

from __future__ import annotations

from dataclasses import replace
import json

import numpy as np
from scipy.integrate import solve_ivp

from ._paths import output_dir

from migration1d import (  # noqa: E402
    EroumeChemistry,
    HolmesChemistry,
    HolmesParameters,
    MechanicsParameters,
    MigrationModel,
    RacRhoForce,
    SimulationConfig,
    run_simulation,
)

OUTPUT = output_dir("three_species/testing/stimulus_panel")
START, STOP, DURATION = 50.0, 70.0, 300.0
EROUME_CUE_END = 10.0


class RelativeGradient:
    def __init__(self, *, persistent: bool):
        self.persistent = persistent

    def __call__(self, time, position):
        active = time >= START and (self.persistent or time < STOP)
        if not active:
            return np.zeros_like(position)
        fraction = (position-position[0])/(position[-1]-position[0])
        return 0.3 + 0.1*(2.0*fraction-1.0)


class LeftTopHatPulse:
    def __call__(self, time, position):
        if not START <= time < STOP:
            return np.zeros_like(position)
        fraction = (position-position[0])/(position[-1]-position[0])
        return 10.0*(fraction <= 0.04)


class PolarizedHolmesChemistry(HolmesChemistry):
    """Strong right-polarized, pointwise mass-balanced initial condition."""

    def initial_profiles(self, position):
        fraction = (position-position[0])/(position[-1]-position[0])
        polarity = 2.0*fraction-1.0
        active = np.repeat(self.initial_membrane[:, None], len(position), axis=1)
        active += np.asarray((0.25, 0.30, -0.80))[:, None]*polarity
        inactive = self.parameters.totals[:, None]-active
        return active, inactive


def holmes_model(cue, *, polarized=False):
    parameters = replace(HolmesParameters(), signal_baseline=0.0)
    chemistry_class = PolarizedHolmesChemistry if polarized else HolmesChemistry
    chemistry = chemistry_class(parameters, cue=cue)
    force = RacRhoForce(
        rac_gain=0.08, rho_gain=0.08,
        rac_half_level=parameters.initial_active[1],
        rho_half_level=parameters.initial_active[2], hill_exponent=4.0,
    )
    return MigrationModel(
        mechanics=MechanicsParameters(
            spring_stiffness=1.0, dashpot_viscosity=2.0,
            left_friction=1.0, right_friction=1.0,
        ),
        chemistry=chemistry, left_force=force, right_force=force,
    )


def eroume_model():
    # The supported adapter preserves the dimensionless diffusion length when
    # the native 70.54-unit WGD axis is represented by this framework's L=10.
    chemistry = EroumeChemistry(reference_length=10.0)
    rac0 = float(chemistry.initial_membrane[3])
    rho0 = float(chemistry.initial_membrane[5])
    force = RacRhoForce(
        rac_field="Rac_act", rho_field="Rho_act",
        rac_gain=0.08, rho_gain=0.08,
        rac_half_level=rac0, rho_half_level=rho0, hill_exponent=4.0,
    )
    return MigrationModel(
        mechanics=MechanicsParameters(
            spring_stiffness=1.0, dashpot_viscosity=2.0,
            left_friction=1.0, right_friction=1.0,
        ),
        chemistry=chemistry, left_force=force, right_force=force,
    )


def simulate(model, *, release_time, rtol=2.0e-5, atol=1.0e-8):
    return run_simulation(
        model,
        SimulationConfig(
            nodes=101, duration=DURATION, output_interval=1.0,
            max_step=0.25, noise_amplitude=0.0,
            mechanics_activation_time=release_time,
            rtol=rtol, atol=atol,
        ),
    )


def contrast(values):
    left, right = values[:, 0], values[:, -1]
    return (right-left)/np.maximum(np.abs(right)+np.abs(left), 1.0e-12)


def save_result(name, result):
    np.savez_compressed(
        OUTPUT/f"{name}_trajectory.npz",
        times=result.times, position=result.position, velocity=result.velocity,
        membrane=result.membrane, cytoplasm=result.cytoplasm,
        applied_edge_force=result.applied_edge_force,
        spring_force=result.spring_force, dashpot_force=result.dashpot_force,
        species=np.asarray(result.model.chemistry.species),
    )


def draw_kymograph(ax, result, field, title):
    import matplotlib.pyplot as plt

    values = result.membrane[:, field]
    times = np.broadcast_to(result.times[:, None], result.position.shape)
    image = ax.pcolormesh(
        result.position, times, values, shading="auto", cmap="magma",
        vmin=float(values.min()), vmax=float(values.max()), rasterized=True,
    )
    ax.plot(result.rear, result.times, color="white", lw=0.7)
    ax.plot(result.front, result.times, color="white", lw=0.7)
    ax.axhline(START, color="#54d6ff", ls="--", lw=0.8)
    ax.axhline(STOP, color="#54d6ff", ls="--", lw=0.8)
    ax.text(
        0.98, 0.98, f"Δx={result.center[-1]-result.center[0]:+.3f}\n"
        f"final contrast={contrast(values)[-1]:+.3g}",
        transform=ax.transAxes, ha="right", va="top", color="white", fontsize=8,
        bbox={"facecolor": "black", "alpha": 0.4, "edgecolor": "none"},
    )
    ax.set(title=title, xlabel="laboratory position", ylabel="time")
    plt.colorbar(image, ax=ax, pad=0.01, label="active Rac")


def render(results):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    holmes_order = (
        ("holmes_transient", "Native transient gradient"),
        ("holmes_persistent", "Native persistent gradient"),
        ("holmes_top_hat", "Strong 4% top-hat pulse"),
        ("holmes_polarized_ic", "Polarized initial state; no cue"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0),
                             constrained_layout=True)
    for ax, (key, title) in zip(axes.flat, holmes_order):
        draw_kymograph(ax, results[key], 1, title)
    fig.suptitle(
        "Holmes six-field Rac persistence challenges\n"
        "white curves: moving boundaries · cyan dashes: standard cue window"
    )
    fig.savefig(OUTPUT/"00_holmes_rac_persistence_kymographs.png", dpi=220)
    plt.close(fig)

    colors = {
        "holmes_transient": "#0072B2", "holmes_persistent": "#009E73",
        "holmes_top_hat": "#D55E00", "holmes_polarized_ic": "#CC79A7",
        "eroume_transient": "#7A3E9D",
    }
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 8.0), sharex=True,
                             constrained_layout=True)
    for key, result in results.items():
        rac_index = 3 if key == "eroume_transient" else 1
        axes[0].plot(result.times, contrast(result.membrane[:, rac_index]),
                     color=colors[key], label=key.replace("_", " "))
        axes[1].plot(result.times, result.center-result.center[0],
                     color=colors[key])
    axes[0].set(ylabel="active Rac front/back contrast")
    axes[0].legend(frameon=False, ncol=2)
    axes[1].set(xlabel="time", ylabel="center displacement")
    for ax in axes:
        ax.axvspan(START, STOP, color="#54d6ff", alpha=0.12)
        ax.axvline(EROUME_CUE_END, color="#7A3E9D", ls="--", lw=0.8)
    fig.suptitle("Persistence and migration across chemical challenges")
    fig.savefig(OUTPUT/"01_persistence_and_displacement.png", dpi=220)
    plt.close(fig)

    eroume = results["eroume_transient"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.5), sharey=True,
                             constrained_layout=True)
    for ax, field, label in zip(axes, (1, 3, 5), ("Cdc42", "Rac", "Rho")):
        values = eroume.membrane[:, field]
        time_grid = np.broadcast_to(eroume.times[:, None], eroume.position.shape)
        image = ax.pcolormesh(eroume.position, time_grid, values, shading="auto",
                              cmap="magma", rasterized=True)
        ax.plot(eroume.rear, eroume.times, color="white", lw=0.7)
        ax.plot(eroume.front, eroume.times, color="white", lw=0.7)
        ax.axhline(EROUME_CUE_END, color="#54d6ff", ls="--", lw=0.8)
        ax.set(title=f"active {label}", xlabel="laboratory position")
        plt.colorbar(image, ax=ax, pad=0.01)
    axes[0].set_ylabel("time")
    fig.suptitle(
        "Nine-field Eroumé/Marée chemistry with endpoint mechanics\n"
        "cyan dash: directional cue removed; mechanics released at t=10"
    )
    fig.savefig(OUTPUT/"02_eroume_nine_field_kymographs.png", dpi=220)
    plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = {
        "holmes_transient": simulate(
            holmes_model(RelativeGradient(persistent=False)), release_time=STOP
        ),
        "holmes_persistent": simulate(
            holmes_model(RelativeGradient(persistent=True)), release_time=STOP
        ),
        "holmes_top_hat": simulate(
            holmes_model(LeftTopHatPulse()), release_time=STOP
        ),
        "holmes_polarized_ic": simulate(
            holmes_model(lambda time, position: np.zeros_like(position), polarized=True),
            release_time=0.0,
        ),
        "eroume_transient": simulate(
            eroume_model(), release_time=EROUME_CUE_END,
            rtol=1.0e-7, atol=1.0e-10
        ),
    }
    for name, result in results.items():
        save_result(name, result)
    render(results)

    summary = {"experiment": "GTPase polarity persistence campaign", "conditions": {}}
    for name, result in results.items():
        rac_index = 3 if name == "eroume_transient" else 1
        rho_index = 5 if name == "eroume_transient" else 2
        rac_contrast = contrast(result.membrane[:, rac_index])
        rho_contrast = contrast(result.membrane[:, rho_index])
        membrane_total, cytoplasm_total = result.compartment_totals()
        initial_total = membrane_total[0]+cytoplasm_total[0]
        final_total = membrane_total[-1]+cytoplasm_total[-1]
        if name == "eroume_transient":
            drift = max(
                abs(
                    final_total[list(pair)].sum()
                    / initial_total[list(pair)].sum() - 1.0
                )
                for pair in ((0, 1), (2, 3), (4, 5))
            )
        else:
            drift = max(
                abs(final_total[index]/initial_total[index]-1.0)
                for index in (0, 1, 2)
            )
        summary["conditions"][name] = {
            "final_rac_contrast": float(rac_contrast[-1]),
            "maximum_absolute_rac_contrast": float(np.max(np.abs(rac_contrast))),
            "final_rho_contrast": float(rho_contrast[-1]),
            "net_displacement": float(result.center[-1]-result.center[0]),
            "minimum_length": float(result.length.min()),
            "maximum_length": float(result.length.max()),
            "maximum_gtpase_relative_drift": float(drift),
        }
    summary["eroume_storage_note"] = (
        "Nine exact dynamic PDE fields occupy migration1d.membrane; the paired "
        "cytoplasm array is an unused zero-valued storage companion."
    )
    (OUTPUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
