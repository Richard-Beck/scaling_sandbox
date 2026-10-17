#!/usr/bin/env python3
"""Validate and render persistent nine-field Eroume polarity in migration1d."""

from __future__ import annotations

import json

import numpy as np

from ._paths import output_dir

from migration1d import (  # noqa: E402
    EroumeChemistry,
    EroumeGeometry,
    MechanicsParameters,
    MigrationModel,
    RacRhoForce,
    SimulationConfig,
    run_simulation,
)

OUTPUT = output_dir("three_species/calibration/geometry_validation")
NATIVE_LENGTH = EroumeGeometry().native_axis_length


def make_model(length: float, motile: bool):
    chemistry = EroumeChemistry(reference_length=length)
    force = RacRhoForce(
        rac_field="Rac_act", rho_field="Rho_act",
        rac_gain=0.08, rho_gain=0.08,
        rac_half_level=3.0, rho_half_level=1.25, hill_exponent=4.0,
    )
    return MigrationModel(
        mechanics=MechanicsParameters(
            reference_length=length, spring_stiffness=1.0,
            dashpot_viscosity=2.0, left_friction=1.0, right_friction=1.0,
        ),
        chemistry=chemistry,
        left_force=force if motile else (lambda state: 0.0),
        right_force=force if motile else (lambda state: 0.0),
    )


def run_condition(length: float, motile: bool):
    return run_simulation(
        make_model(length, motile),
        SimulationConfig(
            nodes=101, duration=200.0, output_interval=1.0, max_step=0.2,
            noise_amplitude=0.0,
            mechanics_activation_time=10.0 if motile else None,
        ),
    )


def contrast(values):
    maximum, minimum = values.max(axis=1), values.min(axis=1)
    return (maximum-minimum)/np.maximum(maximum+minimum, 1.0e-12)


def save(name, result):
    np.savez_compressed(
        OUTPUT/f"{name}_trajectory.npz",
        times=result.times, position=result.position, velocity=result.velocity,
        membrane=result.membrane, cytoplasm=result.cytoplasm,
        applied_edge_force=result.applied_edge_force,
        spring_force=result.spring_force, dashpot_force=result.dashpot_force,
        species=np.asarray(result.model.chemistry.species),
    )


def render(results):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = (
        ("native_fixed", "Native 70.54-µm coordinate · fixed"),
        ("normalized_fixed", "Normalized 10-unit coordinate · fixed"),
        ("normalized_motile", "Normalized coordinate · mechanics released at 10 s"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.8), sharey=True,
                             constrained_layout=True)
    for ax, (key, title) in zip(axes, order):
        result = results[key]
        values = result.membrane[:, 3]
        time_grid = np.broadcast_to(result.times[:, None], result.position.shape)
        image = ax.pcolormesh(
            result.position, time_grid, values, shading="auto", cmap="magma",
            vmin=float(values.min()), vmax=float(values.max()), rasterized=True,
        )
        ax.plot(result.rear, result.times, color="white", lw=0.75)
        ax.plot(result.front, result.times, color="white", lw=0.75)
        ax.axhline(10.0, color="#54d6ff", ls="--", lw=0.9)
        ax.text(
            0.98, 0.98,
            f"final contrast={contrast(values)[-1]:.3f}\n"
            f"Δx={result.center[-1]-result.center[0]:+.2f}",
            transform=ax.transAxes, ha="right", va="top", color="white",
            fontsize=8,
            bbox={"facecolor": "black", "alpha": 0.42, "edgecolor": "none"},
        )
        ax.set(title=title, xlabel="laboratory position")
        plt.colorbar(image, ax=ax, pad=0.01, label="active Rac")
    axes[0].set_ylabel("time")
    fig.suptitle(
        "Eroumé–Marée nine-field polarity after the directional cue ends\n"
        "cyan dash: cue becomes spatially uniform"
    )
    fig.savefig(OUTPUT/"00_rac_persistence_validation_kymographs.png", dpi=220)
    plt.close(fig)

    motile = results["normalized_motile"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.5), sharey=True,
                             constrained_layout=True)
    for ax, field, label in zip(axes, (1, 3, 5), ("Cdc42", "Rac", "Rho")):
        values = motile.membrane[:, field]
        time_grid = np.broadcast_to(motile.times[:, None], motile.position.shape)
        image = ax.pcolormesh(motile.position, time_grid, values, shading="auto",
                              cmap="magma", rasterized=True)
        ax.plot(motile.rear, motile.times, color="white", lw=0.75)
        ax.plot(motile.front, motile.times, color="white", lw=0.75)
        ax.axhline(10.0, color="#54d6ff", ls="--", lw=0.9)
        ax.set(title=f"active {label}", xlabel="laboratory position")
        plt.colorbar(image, ax=ax, pad=0.01)
    axes[0].set_ylabel("time")
    fig.suptitle("Persistent Cdc42–Rac front and Rho rear on the moving interval")
    fig.savefig(OUTPUT/"01_motile_active_gtpase_kymographs.png", dpi=220)
    plt.close(fig)

    colors = {"Cdc42": "#009E73", "Rac": "#CC79A7", "Rho": "#E69F00"}
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 9.5), sharex=True,
                             constrained_layout=True)
    for field, label in ((1, "Cdc42"), (3, "Rac"), (5, "Rho")):
        signed = contrast(motile.membrane[:, field])
        if label == "Rho":
            signed *= np.sign(
                motile.membrane[:, field, -1]-motile.membrane[:, field, 0]
            )
        axes[0].plot(motile.times, signed, color=colors[label], label=label)
    axes[0].set_ylabel("spatial contrast")
    axes[0].legend(frameon=False, ncol=3)
    axes[1].plot(motile.times, motile.center-motile.center[0], color="#0072B2")
    axes[1].set_ylabel("center displacement")
    axes[2].plot(motile.times, 100*(motile.length/motile.length[0]-1),
                 color="#D55E00")
    axes[2].set(xlabel="time", ylabel="length change (%)")
    for ax in axes:
        ax.axvspan(0.0, 10.0, color="#54d6ff", alpha=0.12)
    fig.suptitle("Chemical persistence and endpoint migration")
    fig.savefig(OUTPUT/"02_polarity_and_migration_metrics.png", dpi=220)
    plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = {
        "native_fixed": run_condition(NATIVE_LENGTH, False),
        "normalized_fixed": run_condition(10.0, False),
        "normalized_motile": run_condition(10.0, True),
    }
    for name, result in results.items():
        save(name, result)
    render(results)

    native = results["native_fixed"].membrane[-1]
    normalized = results["normalized_fixed"].membrane[-1]
    summary = {
        "experiment": "native-protocol Eroume nine-field persistence in migration1d",
        "adaptation": {
            "native_axis_length": NATIVE_LENGTH,
            "normalized_axis_length": 10.0,
            "diffusion_scale": (10.0/NATIVE_LENGTH)**2,
            "initial_state": "native public VCML concentrations; no pre-relaxation",
            "directional_cue": "Ic=2.6+0.05*x_source through t=10 s",
            "postcue_input": "spatially uniform Ic=2.95",
            "mechanics_release_time": 10.0,
            "storage": "nine dynamic rows plus unused zero companion rows",
        },
        "terminal_native_normalized_rms_difference": float(
            np.sqrt(np.mean((native-normalized)**2))
        ),
        "conditions": {},
    }
    for name, result in results.items():
        totals, _ = result.compartment_totals()
        pair_drift = {}
        for label, pair in (("cdc42", (0, 1)), ("rac", (2, 3)), ("rho", (4, 5))):
            pool = totals[:, list(pair)].sum(axis=1)
            pair_drift[label] = float(np.max(np.abs(pool/pool[0]-1.0)))
        summary["conditions"][name] = {
            "final_cdc42_contrast": float(contrast(result.membrane[:, 1])[-1]),
            "final_rac_contrast": float(contrast(result.membrane[:, 3])[-1]),
            "final_rho_contrast": float(contrast(result.membrane[:, 5])[-1]),
            "net_displacement": float(result.center[-1]-result.center[0]),
            "minimum_length": float(result.length.min()),
            "maximum_length": float(result.length.max()),
            "gtpase_pair_relative_drift": pair_drift,
        }
    (OUTPUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
