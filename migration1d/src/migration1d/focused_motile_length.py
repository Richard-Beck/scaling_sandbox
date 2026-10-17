"""Focused motile length and global-tension sweep for both chemistry models.

The production design holds chemistry fixed, varies reference length and a
100-fold mechanical-stiffness contrast, applies matched single, dual, and
reversal cues, and compares at most three global tension-feedback strengths.
One SLURM array task writes one JSON record and one compressed trajectory; a
dependent merge validates the complete design and renders representative
kymographs with feedback strength as the row variable.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
import json
import math
import os
from pathlib import Path

import numpy as np

from migration1d.forces import EndpointState, RacRhoForce
from migration1d.holmes import HolmesChemistry
from migration1d.model import (
    MechanicsParameters,
    MigrationModel,
    SimulationConfig,
    run_simulation,
)
from migration1d.single_species.model import (
    ActivationPulse,
    PulseTrain,
    SingleSpeciesParameters,
    StripWavePinningChemistry,
)
from migration1d.single_species.observables import classify_trajectory
from migration1d.tension import GlobalTensionFeedback, SpatialTensionFeedback
from migration1d.three_species.model import (
    GaussianRacPulse,
    RacPulseTrain,
    canonical_holmes_parameters,
    homogeneous_equilibrium,
)
from migration1d.three_species.observables import profile_metrics


CHEMISTRY_CLASSES = ("single_species", "three_species")
MECHANICAL_CASES = ("soft", "stiff")
PROTOCOLS = ("single", "dual_equal", "dual_biased", "reversal")
DEFAULT_LENGTHS = (5.0, 7.5, 10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0)
DEFAULT_TENSION_STRENGTHS = (0.0, 5.0, 20.0)


def parse_values(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def parse_names(text: str) -> list[str]:
    return [value.strip() for value in text.split(",") if value.strip()]


def output_root(args: argparse.Namespace) -> Path:
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    return Path(configured).resolve() if configured else args.output.resolve()


def all_tasks(args: argparse.Namespace) -> list[tuple[str, float, str, str, float]]:
    lengths = parse_values(args.lengths)
    mechanics = parse_names(args.mechanical_cases)
    protocols = parse_names(args.protocols)
    tension_strengths = parse_values(args.tension_strengths)
    if not lengths or any(length <= 0.0 for length in lengths):
        raise ValueError("lengths must contain positive values")
    if not tension_strengths or len(tension_strengths) > 3:
        raise ValueError("tension-strengths must contain between one and three values")
    if len(set(tension_strengths)) != len(tension_strengths):
        raise ValueError("tension-strengths must be unique")
    if any(not math.isfinite(value) or value < 0.0 for value in tension_strengths):
        raise ValueError("tension-strengths must be finite and non-negative")
    unknown_mechanics = set(mechanics) - set(MECHANICAL_CASES)
    unknown_protocols = set(protocols) - set(PROTOCOLS)
    if unknown_mechanics:
        raise ValueError(f"unknown mechanical cases: {sorted(unknown_mechanics)}")
    if unknown_protocols:
        raise ValueError(f"unknown protocols: {sorted(unknown_protocols)}")
    return [
        (args.chemistry_class, float(length), mechanical_case, protocol, float(alpha_t))
        for length in lengths
        for mechanical_case in mechanics
        for protocol in protocols
        for alpha_t in tension_strengths
    ]


def mechanics_parameters(length: float, mechanical_case: str) -> MechanicsParameters:
    """Return length-normalized mechanics with a 100-fold stiffness contrast."""
    multiplier = {"soft": 0.1, "stiff": 10.0}[mechanical_case]
    return MechanicsParameters(
        reference_length=length,
        spring_stiffness=multiplier * 10.0 / length,
        dashpot_viscosity=20.0 / length,
        left_friction=1.0,
        right_friction=1.0,
    )


def node_count(length: float, minimum: int, maximum_spacing: float) -> int:
    return max(minimum, int(math.ceil(length / maximum_spacing)) + 1)


@dataclass(frozen=True)
class CenteredHillForce:
    """Outward active-GTPase drive centered on the homogeneous resting state."""

    baseline: float
    gain: float = 0.08
    half_level: float = 1.0
    exponent: float = 4.0

    def hill(self, value: float) -> float:
        positive = max(float(value), 0.0)
        power = positive**self.exponent
        return power / (self.half_level**self.exponent + power)

    def __call__(self, state: EndpointState) -> float:
        centered = self.gain * (
            self.hill(state.sample("active_gtpase", "membrane"))
            - self.hill(self.baseline)
        )
        outward_sign = -1.0 if state.side == "left" else 1.0
        return outward_sign * centered


def activation_protocol(args: argparse.Namespace, protocol: str):
    left_amplitude = args.stimulus_amplitude
    right_amplitude = args.stimulus_amplitude
    if protocol == "dual_biased":
        left_amplitude *= 1.0 + args.dual_bias
        right_amplitude *= 1.0 - args.dual_bias
    pulses = []
    if protocol in {"single", "dual_equal", "dual_biased", "reversal"}:
        pulses.append(ActivationPulse(
            amplitude=left_amplitude,
            center_fraction=0.1,
            width_fraction=args.stimulus_width,
            start=args.stimulus_start,
            stop=args.stimulus_stop,
        ))
    if protocol in {"dual_equal", "dual_biased"}:
        pulses.append(ActivationPulse(
            amplitude=right_amplitude,
            center_fraction=0.9,
            width_fraction=args.stimulus_width,
            start=args.stimulus_start,
            stop=args.stimulus_stop,
        ))
    if protocol == "reversal":
        pulses.append(ActivationPulse(
            amplitude=args.reversal_amplitude,
            center_fraction=0.9,
            width_fraction=args.stimulus_width,
            start=args.reversal_start,
            stop=args.reversal_stop,
        ))
    return PulseTrain(tuple(pulses))


def rac_protocol(args: argparse.Namespace, protocol: str):
    left_amplitude = args.stimulus_amplitude
    right_amplitude = args.stimulus_amplitude
    if protocol == "dual_biased":
        left_amplitude *= 1.0 + args.dual_bias
        right_amplitude *= 1.0 - args.dual_bias
    pulses = []
    if protocol in {"single", "dual_equal", "dual_biased", "reversal"}:
        pulses.append(GaussianRacPulse(
            left_amplitude, 0.1, args.stimulus_width,
            args.stimulus_start, args.stimulus_stop,
        ))
    if protocol in {"dual_equal", "dual_biased"}:
        pulses.append(GaussianRacPulse(
            right_amplitude, 0.9, args.stimulus_width,
            args.stimulus_start, args.stimulus_stop,
        ))
    if protocol == "reversal":
        pulses.append(GaussianRacPulse(
            args.reversal_amplitude, 0.9, args.stimulus_width,
            args.reversal_start, args.reversal_stop,
        ))
    return RacPulseTrain(tuple(pulses))


def make_model(
    args: argparse.Namespace,
    length: float,
    mechanical_case: str,
    protocol: str,
    alpha_t: float = 0.0,
) -> MigrationModel:
    mechanics = mechanics_parameters(length, mechanical_case)
    if args.chemistry_class == "single_species":
        parameters = SingleSpeciesParameters(
            total_abundance=args.single_mu,
            dimensionless_length=length * math.sqrt(1.0 / args.single_active_diffusivity),
            cytoplasm_depth=args.single_h,
            diffusion_ratio=args.single_diffusion_ratio,
            basal_ratio=args.single_kappa,
            feedback_ratio=args.single_gamma,
            hill_coefficient=args.single_hill,
            inactivation_rate=1.0,
            domain_length=length,
        )
        chemistry = StripWavePinningChemistry(
            parameters, activation=activation_protocol(args, protocol)
        )
        force = CenteredHillForce(
            baseline=float(chemistry.initial_membrane[0]),
            gain=args.force_gain,
            half_level=parameters.half_saturation,
        )
    elif args.chemistry_class == "three_species":
        parameters = canonical_holmes_parameters(
            kappa_r=args.three_kappa_r, hill_exponent=args.three_hill
        )
        rest = homogeneous_equilibrium(parameters)
        parameters = replace(
            parameters, initial_active=tuple(float(value) for value in rest)
        )
        chemistry = HolmesChemistry(parameters, cue=rac_protocol(args, protocol))
        force = RacRhoForce(
            rac_gain=args.force_gain,
            rho_gain=args.force_gain,
            rac_half_level=float(rest[1]),
            rho_half_level=float(rest[2]),
            hill_exponent=4.0,
        )
    else:
        raise ValueError(f"unknown chemistry class {args.chemistry_class!r}")
    if args.tension_mode == "global":
        tension_feedback = GlobalTensionFeedback(alpha_t=alpha_t)
    else:
        tension_feedback = SpatialTensionFeedback(
            diffusivity=args.spatial_tension_diffusivity,
            relaxation_time=args.spatial_tension_relaxation_time,
            generation_rate=args.spatial_tension_generation_rate,
            alpha_t=alpha_t,
        )
    return MigrationModel(
        mechanics=mechanics, chemistry=chemistry,
        left_force=force, right_force=force,
        tension_feedback=tension_feedback,
    )


def endpoint_contrast(profiles: np.ndarray) -> np.ndarray:
    return (profiles[:, -1] - profiles[:, 0]) / np.maximum(
        np.abs(profiles[:, -1]) + np.abs(profiles[:, 0]), 1.0e-12
    )


def polarity_series(result, chemistry_class: str) -> np.ndarray:
    if chemistry_class == "single_species":
        return endpoint_contrast(result.membrane[:, 0])
    return np.asarray([
        (
            float(metrics["cdc42_moment"])
            + float(metrics["rac_moment"])
            - float(metrics["rho_moment"])
        ) / 3.0
        for metrics in (profile_metrics(profiles) for profiles in result.membrane)
    ])


def interval_mean(times: np.ndarray, values: np.ndarray, start: float, stop: float) -> float | None:
    mask = (times >= start) & (times <= stop)
    return float(np.mean(values[mask])) if np.any(mask) else None


def conservation_errors(result, chemistry_class: str, h: float) -> dict[str, float]:
    membrane, cytoplasm = result.compartment_totals()
    total = membrane + (h * cytoplasm if chemistry_class == "single_species" else cytoplasm)
    relative = np.max(
        np.abs(total - total[0]) / np.maximum(np.abs(total[0]), 1.0e-14), axis=0
    )
    names = ("active_gtpase",) if chemistry_class == "single_species" else (
        "cdc42", "rac", "rho"
    )
    return {
        f"maximum_{name}_conservation_error": float(value)
        for name, value in zip(names, relative)
    }


def tension_feedback_signal(result) -> np.ndarray:
    """Return the non-negative signal multiplied by ``alpha_t`` in the model."""
    if isinstance(result.model.tension_feedback, SpatialTensionFeedback):
        if result.tension_excess is None:
            raise ValueError("spatial tension result is missing its tension field")
        return np.maximum(result.tension_excess, 0.0)
    return result.elastic_tension


def result_row(args, task, result, polarity: np.ndarray) -> dict[str, object]:
    chemistry_class, reference_length, mechanical_case, protocol, alpha_t = task
    center_velocity = 0.5 * (result.edge_velocity[:, 0] + result.edge_velocity[:, 1])
    absolute_speed = np.abs(center_velocity)
    directed_speed = center_velocity * np.sign(polarity)
    intervals = {
        "cue": (args.stimulus_start, args.stimulus_stop),
        "early": (args.stimulus_stop, min(100.0, args.duration)),
        "established": (100.0, min(args.reversal_start, args.duration)),
        "reversal_cue": (args.reversal_start, min(args.reversal_stop, args.duration)),
        "reorientation": (args.reversal_stop, min(350.0, args.duration)),
        "late": (350.0, args.duration),
    }
    row: dict[str, object] = {
        "chemistry_class": chemistry_class,
        "reference_length_um": reference_length,
        "mechanical_case": mechanical_case,
        "protocol": protocol,
        "tension_feedback_mode": args.tension_mode,
        "tension_feedback_alpha": alpha_t,
        "nodes": result.config.nodes,
        "spring_stiffness": result.model.mechanics.spring_stiffness,
        "dashpot_viscosity": result.model.mechanics.dashpot_viscosity,
        "duration_seconds": args.duration,
        "net_displacement_um": float(result.center[-1] - result.center[0]),
        "minimum_length_um": float(np.min(result.length)),
        "maximum_length_um": float(np.max(result.length)),
        "final_length_um": float(result.length[-1]),
        "maximum_absolute_strain": float(np.max(np.abs(result.strain))),
        "terminal_signed_polarity": float(polarity[-1]),
        "terminal_absolute_polarity": float(abs(polarity[-1])),
    }
    elastic_tension = result.elastic_tension
    feedback_signal = tension_feedback_signal(result)
    added_deactivation = alpha_t * feedback_signal
    mean_added_deactivation = (
        added_deactivation
        if added_deactivation.ndim == 1
        else np.mean(added_deactivation, axis=1)
    )
    row.update({
        "median_elastic_tension": float(np.median(elastic_tension)),
        "p90_elastic_tension": float(np.quantile(elastic_tension, 0.9)),
        "maximum_elastic_tension": float(np.max(elastic_tension)),
        "fraction_time_under_tension": float(np.mean(elastic_tension > 1.0e-12)),
        "median_tension_feedback_signal": float(np.median(feedback_signal)),
        "p90_tension_feedback_signal": float(np.quantile(feedback_signal, 0.9)),
        "maximum_tension_feedback_signal": float(np.max(feedback_signal)),
        "p90_tension_added_deactivation_rate": float(
            np.quantile(added_deactivation, 0.9)
        ),
        "maximum_tension_added_deactivation_rate": float(
            np.max(added_deactivation)
        ),
        "integrated_tension_feedback_exposure": float(
            np.trapezoid(mean_added_deactivation, result.times)
        ),
    })
    if isinstance(result.model.tension_feedback, SpatialTensionFeedback):
        feedback = result.model.tension_feedback
        row.update({
            "spatial_tension_diffusivity": feedback.diffusivity,
            "spatial_tension_relaxation_time": feedback.relaxation_time,
            "spatial_tension_generation_rate": feedback.generation_rate,
            "spatial_tension_characteristic_length": feedback.characteristic_length,
        })
    active_profile = result.membrane[:, 0 if chemistry_class == "single_species" else 1]
    classified = classify_trajectory(result.times, active_profile)
    row.update({
        "terminal_polarity_status": classified["state"],
        "terminal_polarization_amplitude": classified["terminal_polarization_amplitude"],
        "terminal_pole_count": classified["terminal_pole_count"],
        "time_to_polarization": classified["time_to_polarization"],
    })
    if chemistry_class == "three_species":
        terminal = profile_metrics(result.membrane[-1])
        row.update({f"terminal_{key}": value for key, value in terminal.items()})
    for label, (start, stop) in intervals.items():
        row[f"mean_signed_velocity_{label}"] = interval_mean(
            result.times, center_velocity, start, stop
        )
        row[f"mean_absolute_speed_{label}"] = interval_mean(
            result.times, absolute_speed, start, stop
        )
        row[f"mean_polarity_projected_speed_{label}"] = interval_mean(
            result.times, directed_speed, start, stop
        )
    row.update(conservation_errors(result, chemistry_class, args.single_h))
    return row


def simulate_task(args: argparse.Namespace, task: tuple[str, float, str, str, float]):
    _, length, mechanical_case, protocol, alpha_t = task
    model = make_model(args, length, mechanical_case, protocol, alpha_t)
    nodes = node_count(length, args.minimum_nodes, args.maximum_spacing)
    result = run_simulation(
        model,
        SimulationConfig(
            nodes=nodes,
            duration=args.duration,
            output_interval=args.output_interval,
            max_step=args.max_step,
            noise_amplitude=0.0,
            mechanics_activation_time=0.0,
            jacobian_sparsity=args.chemistry_class == "three_species",
        ),
    )
    polarity = polarity_series(result, args.chemistry_class)
    return result, polarity, result_row(args, task, result, polarity)


def task_paths(root: Path, index: int) -> tuple[Path, Path]:
    base = root / "tasks" / f"task_{index:03d}"
    return base.with_suffix(".json"), base.with_suffix(".npz")


def run_task(args: argparse.Namespace) -> tuple[Path, Path]:
    tasks = all_tasks(args)
    index = int(os.environ.get("SLURM_ARRAY_TASK_ID", args.task_index))
    if not 0 <= index < len(tasks):
        raise ValueError(f"task index {index} outside 0..{len(tasks) - 1}")
    task = tasks[index]
    alpha_t = task[4]
    json_path, trajectory_path = task_paths(output_root(args), index)
    result, polarity, row = simulate_task(args, task)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        trajectory_path,
        times=result.times,
        position=result.position,
        membrane=result.membrane,
        center=result.center,
        length=result.length,
        center_velocity=0.5 * (result.edge_velocity[:, 0] + result.edge_velocity[:, 1]),
        polarity=polarity,
        elastic_tension=result.elastic_tension,
        tension_feedback_signal=tension_feedback_signal(result),
        tension_added_deactivation_rate=alpha_t * tension_feedback_signal(result),
    )
    json_path.write_text(json.dumps({
        "task_index": index,
        "task_key": list(task),
        "row": row,
        "trajectory": trajectory_path.name,
    }, separators=(",", ":")) + "\n")
    print(json_path)
    return json_path, trajectory_path


def load_complete_results(args: argparse.Namespace):
    root = output_root(args)
    payloads = []
    tasks = all_tasks(args)
    for index, expected in enumerate(tasks):
        json_path, trajectory_path = task_paths(root, index)
        if not json_path.exists() or not trajectory_path.exists():
            raise FileNotFoundError(f"missing task artifacts for index {index}")
        payload = json.loads(json_path.read_text())
        if payload.get("task_index") != index or tuple(payload.get("task_key", ())) != expected:
            raise ValueError(f"task {index} does not match the current design")
        payloads.append(payload)
    return payloads


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def render_length_summary(root: Path, rows: list[dict[str, object]]) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"soft": "#D55E00", "stiff": "#0072B2"}
    figures = []
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, protocol in zip(axes.flat, PROTOCOLS):
        for mechanical_case in MECHANICAL_CASES:
            subset = sorted(
                (row for row in rows if row["protocol"] == protocol
                 and row["mechanical_case"] == mechanical_case),
                key=lambda row: float(row["reference_length_um"]),
            )
            ax.plot(
                [row["reference_length_um"] for row in subset],
                [row["terminal_absolute_polarity"] for row in subset],
                marker="o", color=colors[mechanical_case], label=mechanical_case,
            )
        ax.set(title=protocol.replace("_", " "), xlabel="reference length (µm)",
               ylabel="terminal |signed polarity|")
    axes[0, 0].legend(frameon=False)
    path = root / "length_vs_polarity.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path.name)

    speed_intervals = ("cue", "early", "established", "reversal_cue", "reorientation", "late")
    interval_colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, len(speed_intervals)))
    fig, axes = plt.subplots(4, 2, figsize=(13, 16), sharex=True, constrained_layout=True)
    for row_index, protocol in enumerate(PROTOCOLS):
        for column_index, mechanical_case in enumerate(MECHANICAL_CASES):
            ax = axes[row_index, column_index]
            subset = sorted(
                (row for row in rows if row["protocol"] == protocol
                 and row["mechanical_case"] == mechanical_case),
                key=lambda row: float(row["reference_length_um"]),
            )
            for label, color in zip(speed_intervals, interval_colors):
                values = [row[f"mean_absolute_speed_{label}"] for row in subset]
                if any(value is not None for value in values):
                    ax.plot(
                        [row["reference_length_um"] for row in subset], values,
                        marker="o", label=label.replace("_", " "), color=color,
                    )
            ax.set(title=f"{protocol.replace('_', ' ')} · {mechanical_case}",
                   xlabel="reference length (µm)", ylabel="mean |center velocity|")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=3, frameon=False)
    path = root / "length_vs_interval_speed.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path.name)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    measures = (
        ("net_displacement_um", "net displacement (µm)"),
        ("maximum_absolute_strain", "maximum |strain|"),
        ("final_length_um", "final length (µm)"),
    )
    for ax, (field, label) in zip(axes, measures):
        for protocol, linestyle in zip(PROTOCOLS, ("-", "--", "-.", ":")):
            for mechanical_case in MECHANICAL_CASES:
                subset = sorted(
                    (row for row in rows if row["protocol"] == protocol
                     and row["mechanical_case"] == mechanical_case),
                    key=lambda row: float(row["reference_length_um"]),
                )
                ax.plot(
                    [row["reference_length_um"] for row in subset],
                    [row[field] for row in subset],
                    color=colors[mechanical_case], linestyle=linestyle,
                    label=f"{protocol} · {mechanical_case}",
                )
        ax.set(xlabel="reference length (µm)", ylabel=label)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=4, frameon=False)
    path = root / "length_vs_motility_and_deformation.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path.name)
    return figures


def render_kymograph_montages(args, root: Path, payloads: list[dict]) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    available = np.asarray(sorted(set(parse_values(args.lengths))))
    requested = (5.0, 10.0, 20.0, 35.0, 50.0)
    selected = []
    for value in requested:
        nearest = float(available[np.argmin(np.abs(available - value))])
        if nearest not in selected:
            selected.append(nearest)
    lookup = {tuple(payload["task_key"]): payload for payload in payloads}
    tension_strengths = parse_values(args.tension_strengths)
    figures = []
    for protocol in parse_names(args.protocols):
        for mechanical_case in parse_names(args.mechanical_cases):
            fig, axes = plt.subplots(
                len(tension_strengths), len(selected),
                figsize=(3.2 * len(selected), 3.6 * len(tension_strengths)),
                squeeze=False, constrained_layout=True,
            )
            trajectories = {}
            for alpha_t in tension_strengths:
                for length in selected:
                    payload = lookup[
                        (args.chemistry_class, length, mechanical_case, protocol, alpha_t)
                    ]
                    trajectories[(alpha_t, length)] = np.load(
                        root / "tasks" / payload["trajectory"]
                    )
            field = 0 if args.chemistry_class == "single_species" else 1
            vmax = max(
                float(data["membrane"][:, field].max())
                for data in trajectories.values()
            )
            image = None
            for row, alpha_t in enumerate(tension_strengths):
                for column, length in enumerate(selected):
                    data = trajectories[(alpha_t, length)]
                    times = data["times"]
                    position = data["position"]
                    profile = data["membrane"][:, field]
                    time_grid = np.broadcast_to(times[:, None], position.shape)
                    image = axes[row, column].pcolormesh(
                        position, time_grid, profile, shading="auto", cmap="magma",
                        vmin=0.0, vmax=vmax, rasterized=True,
                    )
                    axes[row, column].plot(
                        position[:, 0], times, color="white", linewidth=0.6
                    )
                    axes[row, column].plot(
                        position[:, -1], times, color="white", linewidth=0.6
                    )
                    if row == 0:
                        axes[row, column].set_title(f"L₀={length:g} µm")
                    if row == len(tension_strengths) - 1:
                        axes[row, column].set_xlabel("laboratory x (µm)")
                axes[row, 0].set_ylabel(f"αₜ={alpha_t:g}\ntime (s)")
            fig.colorbar(image, ax=axes, label="active GTPase" if field == 0 else "active Rac")
            fig.suptitle(f"{args.chemistry_class.replace('_', ' ')} · {protocol.replace('_', ' ')} · {mechanical_case}")
            path = root / f"kymographs_{protocol}_{mechanical_case}.png"
            fig.savefig(path, dpi=180)
            plt.close(fig)
            figures.append(path.name)
    return figures


def merge(args: argparse.Namespace) -> None:
    root = output_root(args)
    payloads = load_complete_results(args)
    rows = [payload["row"] for payload in payloads]
    maximum_error = max(
        float(value)
        for row in rows
        for key, value in row.items()
        if key.startswith("maximum_") and key.endswith("_conservation_error")
    )
    if maximum_error > args.maximum_conservation_error:
        raise ValueError(f"maximum conservation error {maximum_error:g} exceeds threshold")
    if min(float(row["minimum_length_um"]) for row in rows) <= args.minimum_safe_length:
        raise ValueError("at least one trajectory crossed the minimum safe cell length")
    write_csv(root / "focused_motile_length.csv", rows)
    figures = render_kymograph_montages(args, root, payloads)
    summary = {
        "analysis": f"focused_motile_length_{args.tension_mode}_tension_sweep",
        "chemistry_class": args.chemistry_class,
        "tasks": len(rows),
        "lengths_um": parse_values(args.lengths),
        "mechanical_cases": parse_names(args.mechanical_cases),
        "protocols": parse_names(args.protocols),
        "tension_feedback_strengths": parse_values(args.tension_strengths),
        "duration_seconds": args.duration,
        "maximum_conservation_error": maximum_error,
        "minimum_realized_length_um": min(float(row["minimum_length_um"]) for row in rows),
        "figures": figures,
    }
    if args.tension_mode == "spatial":
        summary["spatial_tension_parameters"] = {
            "diffusivity": args.spatial_tension_diffusivity,
            "relaxation_time": args.spatial_tension_relaxation_time,
            "generation_rate": args.spatial_tension_generation_rate,
            "characteristic_length": math.sqrt(
                args.spatial_tension_diffusivity
                * args.spatial_tension_relaxation_time
            ),
        }
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--mode", choices=("task", "merge"), default="task")
    result.add_argument("--chemistry-class", choices=CHEMISTRY_CLASSES, required=True)
    result.add_argument("--output", type=Path, default=Path("results/focused_motile_length"))
    result.add_argument("--task-index", type=int, default=0)
    result.add_argument("--lengths", default=",".join(f"{value:g}" for value in DEFAULT_LENGTHS))
    result.add_argument("--mechanical-cases", default=",".join(MECHANICAL_CASES))
    result.add_argument("--protocols", default=",".join(PROTOCOLS))
    result.add_argument(
        "--tension-mode", choices=("global", "spatial"), default="global"
    )
    result.add_argument(
        "--tension-strengths",
        default=",".join(f"{value:g}" for value in DEFAULT_TENSION_STRENGTHS),
        help="one to three feedback alpha_t values (default: 0,5,20)",
    )
    result.add_argument("--spatial-tension-diffusivity", type=float, default=10.0)
    result.add_argument("--spatial-tension-relaxation-time", type=float, default=10.0)
    result.add_argument("--spatial-tension-generation-rate", type=float, default=0.01)
    result.add_argument("--duration", type=float, default=500.0)
    result.add_argument("--output-interval", type=float, default=1.0)
    result.add_argument("--max-step", type=float, default=0.25)
    result.add_argument("--minimum-nodes", type=int, default=81)
    result.add_argument("--maximum-spacing", type=float, default=0.25)
    result.add_argument("--stimulus-amplitude", type=float, default=1.0)
    result.add_argument("--stimulus-width", type=float, default=0.05)
    result.add_argument("--stimulus-start", type=float, default=10.0)
    result.add_argument("--stimulus-stop", type=float, default=30.0)
    result.add_argument("--dual-bias", type=float, default=0.01)
    result.add_argument("--reversal-amplitude", type=float, default=1.0)
    result.add_argument("--reversal-start", type=float, default=250.0)
    result.add_argument("--reversal-stop", type=float, default=270.0)
    result.add_argument("--force-gain", type=float, default=0.08)
    result.add_argument("--single-mu", type=float, default=2.1875)
    result.add_argument("--single-h", type=float, default=1.0)
    result.add_argument("--single-active-diffusivity", type=float, default=0.05)
    result.add_argument("--single-diffusion-ratio", type=float, default=100.0)
    result.add_argument("--single-kappa", type=float, default=0.067)
    result.add_argument("--single-gamma", type=float, default=1.0)
    result.add_argument("--single-hill", type=float, default=2.0)
    result.add_argument("--three-kappa-r", type=float, default=0.04)
    result.add_argument("--three-hill", type=float, default=3.0)
    result.add_argument("--maximum-conservation-error", type=float, default=1.0e-7)
    result.add_argument("--minimum-safe-length", type=float, default=0.25)
    return result


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    if args.mode == "task":
        run_task(args)
    else:
        merge(args)


if __name__ == "__main__":
    main()
