"""Executable adapters for self-contained WGD v1 model implementations."""

from __future__ import annotations

from dataclasses import replace
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Mapping

import numpy as np
from scipy.optimize import brentq
import scipy.io

from giese_wp import (
    GieseParameters,
    GieseSolver,
    Stimulus,
    make_circular_system,
    membrane_metrics,
    summarize_trajectory,
)
from taniguchi.phase_field import (
    PHASE_PARAMETER_SETS,
    PhaseFieldResult,
    simulate_phase_field,
)
from taniguchi.wgd_metrics import (
    extract_phase_wgd_observables,
    summarize_phase_wgd,
    summarize_phase_wgd_observables,
)
from .contract import geometry_scales, make_summary, write_summary
from .polarity_models import (
    simulate_buttenschoen as simulate_buttenschoen_model,
    simulate_holmes as simulate_holmes_model,
    simulate_mori as simulate_mori_model,
)


class AdapterError(ValueError):
    """A model adapter cannot execute the requested configuration."""


CHIOU_INITIAL_STATE_SHA256 = (
    "5db8e4157e64f56adf5072851b75dd72b3b8dd5b438bd62c61c3c2121c786d80"
)


def _model_config(config: Mapping[str, Any]) -> dict[str, Any]:
    return dict(config.get("model_config", {}))


def _positive_float(values: Mapping[str, Any], key: str, default: float) -> float:
    value = values.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AdapterError(f"model_config.{key} must be numeric")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise AdapterError(f"model_config.{key} must be finite and positive")
    return value


def _choice(values: Mapping[str, Any], key: str, default: str,
            choices: set[str]) -> str:
    value = values.get(key, default)
    if value not in choices:
        raise AdapterError(
            f"model_config.{key} must be one of {sorted(choices)}"
        )
    return str(value)


def _boolean(values: Mapping[str, Any], key: str, default: bool) -> bool:
    value = values.get(key, default)
    if not isinstance(value, bool):
        raise AdapterError(f"model_config.{key} must be boolean")
    return value


def _write_resolved(output: Path, resolved: Mapping[str, Any]) -> None:
    output.mkdir(parents=True, exist_ok=True)
    (output / "resolved_config.json").write_text(
        json.dumps(resolved, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


ONE_D_MODELS = {
    "mori_wave_pinning": (10.0, 300.0, 1.0, 81),
    "holmes_rac_cdc42_rho": (20.0, 300.0, 5.0, 41),
    "buttenschoen_moving_domain": (55.0, 1200.0, 5.0, 81),
}

def _csv_value(value: str | None) -> float | str | None:
    if value is None or value.strip() in {"", "NA", "NaN", "nan"}:
        return None
    try:
        number = float(value)
    except ValueError:
        return value
    return number if math.isfinite(number) else None


def analyze_eroume_export(config: Mapping[str, Any], output: Path) -> dict[str, Any]:
    """Translate exported Eroume output into common WGD statistics."""

    path = output / "raw_metrics.csv"
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    except OSError as error:
        raise AdapterError(f"cannot read raw R output {path}: {error}") from error
    if len(rows) != 1:
        raise AdapterError(f"{path} must contain exactly one data row")
    raw = {key: _csv_value(value) for key, value in rows[0].items()}
    state = raw.pop("state_class", None)
    if not isinstance(state, str):
        raise AdapterError("raw_metrics.csv must contain a state_class string")

    metric_names = (
        "spatial_contrast", "spatial_cv", "response_time",
        "conservation_error", "front_back_contrast", "native_integral",
    )
    model_metrics = {name: raw.get(name) for name in metric_names
                     if raw.get(name) is not None}
    contrast = raw.get("spatial_contrast")
    polarity_strength = float(contrast) if isinstance(contrast, float) else None
    n_domains = None
    domain_value = raw.get("n_domains")
    if isinstance(domain_value, float) and domain_value >= 0:
        n_domains = int(domain_value)
    numerical_failure = state == "numerical_failure"
    status = "numerical_failure" if numerical_failure else "complete"
    response = raw.get("response_time")
    polarity_time = float(response) if isinstance(response, float) else None
    summary = make_summary(
        config, run_status=status,
        status_message=("R solver did not converge" if numerical_failure else ""),
        state_class=state, polarity_strength=polarity_strength,
        n_domains=n_domains,
        time_to_polarity=polarity_time,
        model_metrics=model_metrics,
        provenance={
            "manifest": f"experiments/wgd_v1/manifests/{config['model']}.md",
            "runner": "scripts/generate_wgd_spatial_results.R",
            "raw_metrics": "raw_metrics.csv",
        },
    )
    write_summary(summary, output)
    return summary


def resolve_eroume(config: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve one of the two predeclared Eroume-Maree geometry pairs."""

    if config["scaling_condition"] == "fixed_total":
        raise AdapterError(
            "eroume_maree fixed_total is undefined because the 2D model has "
            "concentration fields but no physical extrusion thickness"
        )
    if config["volume_ratio"] not in {1, 2}:
        raise AdapterError("eroume_maree supports the paired 1x and 2x modes")
    values = _model_config(config)
    aspect_ratio = values.get("aspect_ratio", 1.0)
    if aspect_ratio not in {1, 1.0, 1.29}:
        raise AdapterError("model_config.aspect_ratio must be 1.0 or 1.29")
    environment = values.get("r_environment", "scaling_sandbox")
    if not isinstance(environment, str) or not environment:
        raise AdapterError("model_config.r_environment must be a non-empty string")
    linear_scale, area_scale = geometry_scales(config["volume_ratio"])
    pair_offset = 0 if float(aspect_ratio) == 1.0 else 3
    worker = pair_offset + (2 if config["volume_ratio"] == 1 else 3)
    return {
        "model": "eroume_maree",
        "volume_ratio": config["volume_ratio"],
        "linear_scale": linear_scale,
        "area_scale": area_scale,
        "aspect_ratio": float(aspect_ratio),
        "worker": worker,
        "r_environment": environment,
    }


def simulate_eroume(
    config: Mapping[str, Any], repo_root: Path, output: Path
) -> Path:
    """Run the production Eroume solver and retain its native checkpoint."""

    resolved = resolve_eroume(config)
    _write_resolved(output, resolved)
    checkpoint = output / "trajectory.rds"
    command = [
        "conda", "run", "-n", resolved["r_environment"], "Rscript",
        str(repo_root / "scripts" / "generate_wgd_spatial_results.R"),
        f"--worker={resolved['worker']}", f"--checkpoint={checkpoint}",
    ]
    try:
        subprocess.run(
            command, cwd=repo_root, check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        detail = getattr(error, "stderr", None) or str(error)
        raise AdapterError(f"Eroume runner failed: {detail.strip()}") from error
    if not checkpoint.is_file():
        raise AdapterError(f"Eroume runner did not create {checkpoint}")
    return checkpoint


def analyze_eroume(
    config: Mapping[str, Any], repo_root: Path, output: Path
) -> dict[str, Any]:
    """Derive common statistics from an existing Eroume checkpoint."""

    resolved = resolve_eroume(config)
    checkpoint = output / "trajectory.rds"
    command = [
        "conda", "run", "-n", resolved["r_environment"], "Rscript",
        str(repo_root / "scripts" / "export_wgd_eroume_raw.R"),
        str(checkpoint), str(output),
    ]
    try:
        subprocess.run(
            command, cwd=repo_root, check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        detail = getattr(error, "stderr", None) or str(error)
        raise AdapterError(f"Eroume output adapter failed: {detail.strip()}") from error
    return analyze_eroume_export(config, output)


def resolve_one_d(config: Mapping[str, Any]) -> dict[str, Any]:
    model = config["model"]
    baseline, default_end, default_interval, default_n = ONE_D_MODELS[model]
    values = _model_config(config)
    linear_scale, area_scale = geometry_scales(config["volume_ratio"])
    baseline_length = _positive_float(values, "baseline_length_um", baseline)
    end_time = _positive_float(values, "end_time_s", default_end)
    interval = _positive_float(values, "output_interval_s", default_interval)
    n = values.get("n", default_n)
    if isinstance(n, bool) or not isinstance(n, int) or n < 11:
        raise AdapterError("model_config.n must be an integer of at least 11")
    concentration_scale = (1.0 / linear_scale
                           if config["scaling_condition"] == "fixed_total" else 1.0)
    resolved = {
        "model": model, "volume_ratio": config["volume_ratio"],
        "linear_scale": linear_scale, "area_scale": area_scale,
        "baseline_length_um": baseline_length,
        "length_um": baseline_length * linear_scale,
        "concentration_scale": concentration_scale,
        "n": n, "end_time_s": end_time, "output_interval_s": interval,
    }
    if model == "holmes_rac_cdc42_rho":
        fraction = values.get("baseline_membrane_fraction", 0.5)
        if (isinstance(fraction, bool) or not isinstance(fraction, (int, float))
                or not 0 < float(fraction) < 1):
            raise AdapterError("baseline_membrane_fraction must be between 0 and 1")
        resolved["baseline_membrane_fraction"] = float(fraction)
    return resolved


def simulate_one_d(config: Mapping[str, Any], output: Path) -> Path:
    resolved = resolve_one_d(config)
    _write_resolved(output, resolved)
    common = dict(length=resolved["length_um"],
                  concentration_scale=resolved["concentration_scale"],
                  n=resolved["n"], end_time=resolved["end_time_s"],
                  output_interval=resolved["output_interval_s"])
    if config["model"] == "mori_wave_pinning":
        raw = simulate_mori_model(**common)
    elif config["model"] == "holmes_rac_cdc42_rho":
        raw = simulate_holmes_model(
            **common, volume_ratio=config["volume_ratio"],
            baseline_fraction=resolved["baseline_membrane_fraction"])
    else:
        raw = simulate_buttenschoen_model(**common)
    path = output / "trajectory.npz"
    np.savez_compressed(path, **raw)
    return path


def _profile_domain_count(profile: np.ndarray) -> int:
    amplitude = float(np.max(profile) - np.min(profile))
    if amplitude <= 1e-8 * max(abs(float(np.mean(profile))), 1.0):
        return 0
    threshold = float(np.min(profile) + 0.1 * amplitude)
    peaks = ((profile >= np.roll(profile, 1)) &
             (profile > np.roll(profile, -1)) & (profile >= threshold))
    # These models are no-flux, not periodic: remove wrap-around comparisons.
    peaks[0] = profile[0] > profile[1] and profile[0] >= threshold
    peaks[-1] = profile[-1] > profile[-2] and profile[-1] >= threshold
    return max(1, int(np.count_nonzero(peaks)))


def analyze_one_d(config: Mapping[str, Any], output: Path) -> dict[str, Any]:
    resolved = resolve_one_d(config)
    with np.load(output / "trajectory.npz", allow_pickle=False) as raw:
        times = raw["times"]
        x = raw["x"]
        active = raw["active"]
        inactive = raw["inactive"]
        extras = {key: raw[key] for key in raw.files
                  if key not in {"times", "x", "active", "inactive",
                                 "support_length"}}
    maxima, minima = np.max(active, axis=1), np.min(active, axis=1)
    contrasts = (maxima-minima) / np.maximum(maxima+minima, 1e-12)
    means = np.mean(active, axis=1)
    cvs = np.std(active, axis=1) / np.maximum(np.abs(means), 1e-12)
    final_domains = _profile_domain_count(active[-1])
    target = 0.9 * float(np.max(contrasts))
    reached = np.flatnonzero(contrasts >= target)
    response_time = float(times[reached[0]]) if reached.size and target > 0 else None
    masses = np.mean(active+inactive, axis=1) * resolved["length_um"]
    conservation_error = float(np.max(np.abs(masses/masses[0]-1.0)))
    final_contrast = float(contrasts[-1])
    relative_x = (x - np.min(x)) / max(float(np.max(x) - np.min(x)), 1e-12)
    front = float(np.mean(active[-1, relative_x >= 0.8]))
    back = float(np.mean(active[-1, relative_x <= 0.2]))
    state = ("homogeneous" if final_domains == 0 or final_contrast < 0.05 else
             ("multi_domain" if final_domains > 1 else "polarized"))
    model_metrics = {
        "spatial_contrast": final_contrast,
        "spatial_cv": float(cvs[-1]),
        "max_spatial_contrast": float(np.max(contrasts)),
        "front_back_contrast": abs(front-back) / max(abs(front)+abs(back), 1e-12),
        "conservation_error": conservation_error,
        "initial_native_mass": float(masses[0]),
        "final_native_mass": float(masses[-1]),
    }
    for key, value in extras.items():
        if np.asarray(value).ndim == 0:
            model_metrics[key] = float(value)
    if config["model"] == "holmes_rac_cdc42_rho":
        for species in ("cdc42", "rho"):
            field = np.asarray(extras[species])[-1]
            maximum, minimum = float(np.max(field)), float(np.min(field))
            species_front = float(np.mean(field[relative_x >= 0.8]))
            species_back = float(np.mean(field[relative_x <= 0.2]))
            model_metrics[f"{species}_spatial_contrast"] = (
                (maximum-minimum) / max(maximum+minimum, 1e-12)
            )
            model_metrics[f"{species}_front_back_contrast"] = (
                abs(species_front-species_back) /
                max(abs(species_front)+abs(species_back), 1e-12)
            )
    status = "complete" if conservation_error <= 1e-5 else "quality_failure"
    summary = make_summary(
        config, run_status=status,
        status_message=("" if status == "complete" else
                        f"relative native-mass error {conservation_error:.3g}"),
        state_class=state, polarity_strength=final_contrast,
        n_domains=final_domains, time_to_polarity=response_time,
        model_metrics=model_metrics,
        provenance={
            "manifest": f"experiments/wgd_v1/manifests/{config['model']}.md",
            "solver": "python/wgd_v1/polarity_models.py",
        })
    write_summary(summary, output)
    return summary


def resolve_debelly_one_d(
    config: Mapping[str, Any], repo_root: Path
) -> dict[str, Any]:
    """Resolve the concentration-preserving 1D De Belly WGD experiment."""

    if config["scaling_condition"] not in {
        "wgd_homeostatic", "native_concentration"
    }:
        raise AdapterError(
            "debelly_one_d currently supports only homeostatic MCA density; "
            "the released nodal MCA sum is not a physical whole-cell total"
        )
    if config["stimulus_condition"] != "relative_default":
        raise AdapterError(
            "debelly_one_d currently supports stimulus_condition=relative_default"
        )
    if config["seed"] is not None:
        raise AdapterError("debelly_one_d is deterministic and requires seed=null")
    if config["volume_ratio"] not in {1, 2}:
        raise AdapterError("debelly_one_d supports the predeclared 1x/2x pair")

    values = _model_config(config)
    linear_scale, area_scale = geometry_scales(config["volume_ratio"])
    baseline_length = _positive_float(values, "baseline_length_um", 20.0)
    baseline_partition = values.get("baseline_partition", 100)
    if (isinstance(baseline_partition, bool)
            or not isinstance(baseline_partition, int)
            or baseline_partition < 20):
        raise AdapterError(
            "model_config.baseline_partition must be an integer of at least 20"
        )
    default_partition = round(baseline_partition * linear_scale)
    partition = values.get("partition", default_partition)
    if isinstance(partition, bool) or not isinstance(partition, int) or partition < 20:
        raise AdapterError("model_config.partition must be an integer of at least 20")

    end_time = _positive_float(values, "end_time_s", 350.0)
    dt = _positive_float(values, "dt_s", 1.0)
    opto_start = _positive_float(values, "opto_start_s", 50.0)
    if abs(end_time / dt - round(end_time / dt)) > 1e-10:
        raise AdapterError("model_config.end_time_s must be an integer multiple of dt_s")
    if abs(opto_start / dt - round(opto_start / dt)) > 1e-10:
        raise AdapterError("model_config.opto_start_s must be an integer multiple of dt_s")
    if end_time <= 3 * opto_start:
        raise AdapterError("end_time_s must extend beyond the optogenetic pulse")

    initial_mca_density = _positive_float(
        values, "initial_mca_density", 10.0 / 100.0
    )
    length = baseline_length * linear_scale
    opto_width = length / 5.0
    julia = values.get("julia_executable")
    if julia is not None and (not isinstance(julia, str) or not julia):
        raise AdapterError("model_config.julia_executable must be a non-empty string")
    if julia is None:
        julia = shutil.which("julia")
    if julia is None:
        module_julia = Path(
            "/app/eb/software/Julia/1.8.2-linux-x86_64/bin/julia"
        )
        julia = str(module_julia) if module_julia.is_file() else "julia"

    return {
        "model": "debelly_one_d",
        "volume_ratio": config["volume_ratio"],
        "linear_scale": linear_scale,
        "area_scale": area_scale,
        "baseline_length_um": baseline_length,
        "length_um": length,
        "partition": partition,
        "spatial_step_um": length / partition,
        "initial_mca_density": initial_mca_density,
        "mca_nodal_total": initial_mca_density * partition,
        "opto_width_um": opto_width,
        "opto_width_fraction": opto_width / length,
        "opto_start_s": opto_start,
        "opto_stop_s": 3 * opto_start,
        "end_time_s": end_time,
        "dt_s": dt,
        "julia_executable": julia,
        "solver": str(
            repo_root / "julia" / "debelly" / "one_d" / "protocols" /
            "mechanochemical_opto_front2back.jl"
        ),
    }


def simulate_debelly_one_d(
    config: Mapping[str, Any], repo_root: Path, output: Path
) -> Path:
    """Run the released Julia 1D model and retain tidy native trajectories."""

    resolved = resolve_debelly_one_d(config, repo_root)
    _write_resolved(output, resolved)
    environment = os.environ.copy()
    environment.update({
        "DEBELLY_WRITE_OUTPUTS": "0",
        "DEBELLY_QUIET": "1",
        "DEBELLY_TRAJECTORY_OUTPUT": str(output),
        "DEBELLY_LENGTH": str(resolved["length_um"]),
        "DEBELLY_PARTITION": str(resolved["partition"]),
        "DEBELLY_MCA_DENSITY": str(resolved["initial_mca_density"]),
        "DEBELLY_MCA_NODAL_TOTAL": str(resolved["mca_nodal_total"]),
        "DEBELLY_OPTO_WIDTH": str(resolved["opto_width_um"]),
        "DEBELLY_OPTO_TIME": str(resolved["opto_start_s"]),
        "DEBELLY_DT": str(resolved["dt_s"]),
        "DEBELLY_T": str(resolved["end_time_s"]),
        "DEBELLY_ALPHA_OPTO": "1",
        "DEBELLY_BETA_OPTO": "0",
    })
    project = repo_root / "julia" / "debelly" / "one_d"
    command = [
        resolved["julia_executable"], f"--project={project}", resolved["solver"]
    ]
    try:
        completed = subprocess.run(
            command, cwd=project / "protocols", env=environment, check=True,
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        detail = getattr(error, "stderr", None) or str(error)
        raise AdapterError(f"De Belly Julia runner failed: {detail.strip()}") from error
    (output / "runner_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output / "runner_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    trajectory = output / "trajectory.csv"
    required = (trajectory, output / "mechanics.csv", output / "tension.csv")
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise AdapterError(f"De Belly runner did not create: {', '.join(missing)}")
    return trajectory


def _read_debelly_field(
    path: Path, field_names: tuple[str, ...]
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    try:
        table = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    except OSError as error:
        raise AdapterError(f"cannot read De Belly output {path}: {error}") from error
    if table.size == 0:
        raise AdapterError(f"empty De Belly output {path}")
    times = np.unique(table["time_s"])
    positions = np.unique(table["position_um"])
    expected = times.size * positions.size
    if table.size != expected:
        raise AdapterError(f"nonrectangular De Belly output {path}")
    fields = {
        name: np.asarray(table[name], dtype=float).reshape(times.size, positions.size)
        for name in field_names
    }
    return times, positions, fields


def _debelly_region_metrics(
    times: np.ndarray, positions: np.ndarray, field: np.ndarray,
    prefix: str, opto_start: float
) -> dict[str, float]:
    relative = positions / float(positions[-1])
    proximal = np.mean(field[:, relative <= 0.2], axis=1)
    distal = np.mean(field[:, relative >= 0.8], axis=1)
    proximal_delta = proximal - proximal[0]
    distal_delta = distal - distal[0]
    after = times >= opto_start

    def peak(values: np.ndarray) -> tuple[float, float]:
        indices = np.flatnonzero(after)
        index = indices[int(np.argmax(np.abs(values[after])))]
        return float(values[index]), float(times[index])

    proximal_peak, proximal_time = peak(proximal_delta)
    distal_peak, distal_time = peak(distal_delta)
    return {
        f"{prefix}_proximal_peak_delta": proximal_peak,
        f"{prefix}_distal_peak_delta": distal_peak,
        f"{prefix}_distal_to_proximal_peak_ratio": (
            abs(distal_peak) / max(abs(proximal_peak), 1e-12)
        ),
        f"{prefix}_proximal_peak_time_s": proximal_time,
        f"{prefix}_distal_peak_time_s": distal_time,
        f"{prefix}_distal_peak_lag_s": distal_time - proximal_time,
        f"{prefix}_proximal_absolute_response_auc": float(
            np.trapezoid(np.abs(proximal_delta), times)
        ),
        f"{prefix}_distal_absolute_response_auc": float(
            np.trapezoid(np.abs(distal_delta), times)
        ),
        f"{prefix}_terminal_proximal": float(proximal[-1]),
        f"{prefix}_terminal_distal": float(distal[-1]),
        f"{prefix}_terminal_front_back_contrast": (
            abs(float(proximal[-1] - distal[-1])) /
            max(abs(float(proximal[-1])) + abs(float(distal[-1])), 1e-12)
        ),
    }


def analyze_debelly_one_d(
    config: Mapping[str, Any], repo_root: Path, output: Path
) -> dict[str, Any]:
    """Calculate predeclared biochemical and mechanical communication metrics."""

    resolved = resolve_debelly_one_d(config, repo_root)
    times, positions, chemistry = _read_debelly_field(
        output / "trajectory.csv", ("rac", "rho", "mca_bound", "mca_unbound")
    )
    tension_times, tension_x, mechanics = _read_debelly_field(
        output / "tension.csv", ("tension",)
    )
    if not np.array_equal(times, tension_times):
        raise AdapterError("De Belly chemistry and mechanics time grids differ")
    all_fields = [*chemistry.values(), *mechanics.values()]
    if not all(np.all(np.isfinite(field)) for field in all_fields):
        raise AdapterError("De Belly output contains non-finite values")

    rac = chemistry["rac"]
    maxima, minima = np.max(rac, axis=1), np.min(rac, axis=1)
    contrasts = (maxima - minima) / np.maximum(maxima + minima, 1e-12)
    final_domains = _profile_domain_count(rac[-1])
    target = 0.9 * float(np.max(contrasts))
    reached = np.flatnonzero((times >= resolved["opto_start_s"]) & (contrasts >= target))
    response_time = float(times[reached[0]]) if reached.size and target > 0 else None

    metrics: dict[str, Any] = {
        "spatial_contrast": float(contrasts[-1]),
        "max_spatial_contrast": float(np.max(contrasts)),
        "spatial_cv": float(np.std(rac[-1]) / max(abs(float(np.mean(rac[-1]))), 1e-12)),
        "native_integral": float(np.trapezoid(rac[-1], positions)),
    }
    for field_name in ("rac", "rho", "mca_bound"):
        metrics.update(_debelly_region_metrics(
            times, positions, chemistry[field_name], field_name,
            resolved["opto_start_s"],
        ))
    metrics.update(_debelly_region_metrics(
        times, tension_x, mechanics["tension"], "tension",
        resolved["opto_start_s"],
    ))

    mca = chemistry["mca_bound"] + chemistry["mca_unbound"]
    nodal_totals = np.sum(mca, axis=1)
    conservation_error = float(np.max(np.abs(nodal_totals / nodal_totals[0] - 1)))
    line_integrals = np.trapezoid(mca, positions, axis=1)
    metrics.update({
        "conservation_error": conservation_error,
        "initial_mca_mean_density": float(np.mean(mca[0])),
        "initial_mca_native_line_integral": float(line_integrals[0]),
        "final_mca_native_line_integral": float(line_integrals[-1]),
        "maximum_absolute_displacement": float(np.max(np.abs(
            _read_debelly_field(
                output / "mechanics.csv", ("velocity", "displacement")
            )[2]["displacement"]
        ))),
    })
    initial_density_error = abs(
        metrics["initial_mca_mean_density"] - resolved["initial_mca_density"]
    ) / resolved["initial_mca_density"]
    quality_messages = []
    if conservation_error > 5e-3:
        quality_messages.append(f"native MCA conservation error {conservation_error:.3g}")
    if initial_density_error > 1e-8:
        quality_messages.append(f"initial MCA density error {initial_density_error:.3g}")
    status = "quality_failure" if quality_messages else "complete"
    final_contrast = float(contrasts[-1])
    state = ("homogeneous" if final_domains == 0 or final_contrast < 0.05 else
             ("multi_domain" if final_domains > 1 else "polarized"))
    summary = make_summary(
        config, run_status=status, status_message="; ".join(quality_messages),
        state_class=state, polarity_strength=final_contrast,
        n_domains=final_domains, time_to_polarity=response_time,
        model_metrics=metrics,
        provenance={
            "manifest": "experiments/wgd_v1/manifests/debelly_one_d.md",
            "solver": "julia/debelly/one_d/protocols/mechanochemical_opto_front2back.jl",
            "output_adapter": "python/wgd_v1/adapters.py",
        },
    )
    write_summary(summary, output)
    return summary


def resolve_chiou(config: Mapping[str, Any], repo_root: Path) -> dict[str, Any]:
    """Resolve Chiou geometry/amount scaling and check its archived state."""

    values = _model_config(config)
    baseline_length = _positive_float(values, "baseline_length", 10.0)
    baseline_mass = _positive_float(values, "baseline_mass", 20.0)
    linear_scale, area_scale = geometry_scales(config["volume_ratio"])
    if config["scaling_condition"] in {
        "wgd_homeostatic", "native_concentration"
    }:
        total_mass = baseline_mass * linear_scale
    else:
        total_mass = baseline_mass
    template = Path(values.get(
        "template_path",
        repo_root / "model_resources" / "chiou" / "initial_state_1D.mat",
    ))
    if not template.is_absolute():
        template = repo_root / template
    if "n" in values and (
        isinstance(values["n"], bool) or not isinstance(values["n"], int)
        or values["n"] <= 0 or values["n"] % 2
    ):
        raise AdapterError("model_config.n must be a positive even integer")
    n = int(values.get("n", 500))
    equilibration_time = _positive_float(values, "equilibration_time", 2000.0)
    competition_limit = _positive_float(values, "competition_limit", 1e5)
    output_interval = _positive_float(values, "output_interval", 1.0)
    winner_fraction = _positive_float(values, "winner_fraction", 0.99)
    if winner_fraction <= 0.5 or winner_fraction > 1.0:
        raise AdapterError("model_config.winner_fraction must be in (0.5, 1]")
    digest = None
    if template.is_file():
        digest = hashlib.sha256(template.read_bytes()).hexdigest()
        if digest != CHIOU_INITIAL_STATE_SHA256:
            raise AdapterError(
                "Chiou initial-state checksum does not match the tracked published state"
            )
    return {
        "model": "chiou_mcas",
        "volume_ratio": config["volume_ratio"],
        "linear_scale": linear_scale,
        "area_scale": area_scale,
        "length": baseline_length * linear_scale,
        "total_mass": total_mass,
        "mean_concentration": total_mass / (baseline_length * linear_scale),
        "template_path": str(template),
        "template_available": template.is_file(),
        "template_sha256": digest,
        "n": n,
        "equilibration_time": equilibration_time,
        "competition_limit": competition_limit,
        "output_interval": output_interval,
        "winner_fraction": winner_fraction,
    }


def _load_chiou_kernel(repo_root: Path):
    path = repo_root / "scripts" / "chiou_solver_fast.py"
    specification = importlib.util.spec_from_file_location("chiou_solver_fast", path)
    if specification is None or specification.loader is None:
        raise AdapterError(f"cannot load Chiou solver kernel: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def _periodic_resample(values: np.ndarray, n: int) -> np.ndarray:
    old_n = values.size
    old_x = np.arange(old_n, dtype=float) / old_n
    new_x = np.arange(n, dtype=float) / n
    return np.interp(new_x, np.append(old_x, 1.0), np.append(values, values[0]))


def _chiou_mass(u: np.ndarray, v: np.ndarray, length: float) -> float:
    return float(np.mean(u + v) * length)


def _chiou_peak_count(u: np.ndarray) -> int:
    amplitude = float(np.max(u) - np.min(u))
    threshold = float(np.min(u) + 0.1 * amplitude)
    peaks = ((u >= np.roll(u, 1)) & (u > np.roll(u, -1)) & (u >= threshold))
    return max(1, int(np.count_nonzero(peaks)))


def simulate_chiou(config: Mapping[str, Any], repo_root: Path,
                    output: Path) -> Path:
    resolved = resolve_chiou(config, repo_root)
    _write_resolved(output, resolved)
    if not resolved["template_available"]:
        raise AdapterError(
            "Chiou archived initial state is missing: " + resolved["template_path"]
        )
    state = scipy.io.loadmat(resolved["template_path"], squeeze_me=True)
    if "u" not in state or "v" not in state:
        raise AdapterError("Chiou MAT resource must contain u and v arrays")
    template_u = np.ravel(state["u"]).astype(float)
    template_v = np.ravel(state["v"]).astype(float)
    if template_u.size != template_v.size or template_u.size % 2:
        raise AdapterError("Chiou MAT u and v arrays must have equal, even lengths")

    n = resolved["n"]
    half = n // 2
    old_half = template_u.size // 2
    u_half = _periodic_resample(template_u[:old_half], half)
    v_half = _periodic_resample(template_v[:old_half], half)
    half_length = resolved["length"] / 2.0
    current_half_mass = _chiou_mass(u_half, v_half, half_length)
    target_halves = np.asarray([0.6, 0.4]) * resolved["total_mass"]
    initial_u = np.concatenate([
        u_half * target_halves[0] / current_half_mass,
        u_half * target_halves[1] / current_half_mass,
    ])
    initial_v = np.concatenate([
        v_half * target_halves[0] / current_half_mass,
        v_half * target_halves[1] / current_half_mass,
    ])
    parameters = {
        "D_u": 0.01, "D_v": 1.0, "a": 1.0, "b": 1.0, "k": 0.0,
        "reaction_error": 0.05, "steady_error": 2e-4,
        "initial_dt": 1e-4, "dt_growth": 1.5,
        "output_interval": resolved["output_interval"],
    }
    kernel = _load_chiou_kernel(repo_root)
    equilibrated = kernel.integrate(
        initial_u, initial_v, resolved["length"], resolved["equilibration_time"],
        parameters, separate=True,
    )
    mass_before = _chiou_mass(equilibrated["u"], equilibrated["v"],
                               resolved["length"])
    competition = kernel.integrate(
        equilibrated["u"], equilibrated["v"], resolved["length"],
        resolved["competition_limit"], parameters, separate=False,
        record_masses=True, stop_fraction=resolved["winner_fraction"],
    )
    mass_after = _chiou_mass(competition["u"], competition["v"],
                              resolved["length"])
    fractions = np.column_stack((competition["mass1"], competition["mass2"]))
    denominators = np.sum(fractions, axis=1)
    fractions = np.divide(
        fractions, denominators[:, None], out=np.full_like(fractions, 0.5),
        where=denominators[:, None] > 0,
    )
    path = output / "trajectory.npz"
    np.savez_compressed(
        path,
        times=competition["times"],
        competitor_mass_1=competition["mass1"],
        competitor_mass_2=competition["mass2"],
        initial_u=equilibrated["u"], initial_v=equilibrated["v"],
        final_u=competition["u"], final_v=competition["v"],
        equilibration_steps=equilibrated["accepted_steps"],
        competition_steps=competition["accepted_steps"],
        rejected_steps=competition["rejected_steps"],
    )
    return path


def analyze_chiou(config: Mapping[str, Any], repo_root: Path,
                   output: Path) -> dict[str, Any]:
    resolved = resolve_chiou(config, repo_root)
    with np.load(output / "trajectory.npz", allow_pickle=False) as raw:
        times = raw["times"]
        mass1 = raw["competitor_mass_1"]
        mass2 = raw["competitor_mass_2"]
        initial_u, initial_v = raw["initial_u"], raw["initial_v"]
        final_u, final_v = raw["final_u"], raw["final_v"]
        equilibration_steps = int(raw["equilibration_steps"])
        competition_steps = int(raw["competition_steps"])
        rejected_steps = int(raw["rejected_steps"])
    masses = np.column_stack((mass1, mass2))
    denominators = np.sum(masses, axis=1)
    fractions = np.divide(
        masses, denominators[:, None], out=np.full_like(masses, 0.5),
        where=denominators[:, None] > 0,
    )
    final_winner_fraction = float(np.max(fractions[-1]))
    competition_imbalance = max(0.0, 2.0 * final_winner_fraction - 1.0)
    reached = final_winner_fraction >= resolved["winner_fraction"]
    competition_time = float(times[-1]) if reached else None
    final_domains = _chiou_peak_count(final_u)
    mass_before = _chiou_mass(initial_u, initial_v, resolved["length"])
    mass_after = _chiou_mass(final_u, final_v, resolved["length"])
    conservation_error = abs(mass_after - mass_before) / max(abs(mass_before), 1e-12)
    status = "complete" if conservation_error <= 1e-9 else "quality_failure"
    summary = make_summary(
        config,
        run_status=status,
        status_message=("" if status == "complete"
                        else f"relative native-mass error {conservation_error:.3g}"),
        state_class=("polarized" if reached else
                     ("multi_domain" if final_domains > 1 else "indeterminate")),
        polarity_strength=competition_imbalance,
        n_domains=final_domains,
        time_to_polarity=competition_time,
        model_metrics={
            "competition_time": competition_time,
            "censored_at": (None if reached else float(times[-1])),
            "winner_fraction": final_winner_fraction,
            "competition_imbalance": competition_imbalance,
            "initial_native_mass": mass_before,
            "final_native_mass": mass_after,
            "maximum_relative_mass_error": conservation_error,
            "equilibration_steps": equilibration_steps,
            "competition_steps": competition_steps,
            "rejected_steps": rejected_steps,
            "spatial_contrast": float(
                (np.max(final_u)-np.min(final_u)) /
                max(np.max(final_u)+np.min(final_u), 1e-12)
            ),
            "spatial_cv": float(np.std(final_u) / max(abs(np.mean(final_u)), 1e-12)),
            "conservation_error": conservation_error,
        },
        provenance={
            "manifest": "experiments/wgd_v1/manifests/chiou_mcas.md",
            "solver": "scripts/chiou_solver_fast.py",
            "initial_state": "model_resources/chiou/initial_state_1D.mat",
            "initial_state_sha256": resolved["template_sha256"],
        },
    )
    write_summary(summary, output)
    return summary


def run_chiou(config: Mapping[str, Any], repo_root: Path,
               output: Path) -> dict[str, Any]:
    simulate_chiou(config, repo_root, output)
    return analyze_chiou(config, repo_root, output)


def _giese_equilibrium_roots(parameters: GieseParameters, perimeter: float,
                             area: float, target_mass: float) -> list[tuple[float, float]]:
    p = parameters

    def residual(u: float) -> float:
        v = (target_mass - perimeter * u) / area
        activation = p.k0 + p.gamma * u * u / (p.K * p.K + u * u)
        return v * activation - p.delta * u

    upper = target_mass / perimeter
    samples = np.linspace(0.0, upper, 2001)
    values = np.asarray([residual(value) for value in samples])
    roots: list[float] = []
    for left, right, f_left, f_right in zip(
        samples[:-1], samples[1:], values[:-1], values[1:]
    ):
        if f_left == 0.0:
            roots.append(float(left))
        elif f_left * f_right < 0.0:
            roots.append(float(brentq(residual, left, right)))
    if values[-1] == 0.0:
        roots.append(float(samples[-1]))
    unique = []
    for u in roots:
        v = (target_mass - perimeter * u) / area
        if u >= 0.0 and v >= 0.0 and not any(abs(u - old[0]) < 1e-9 for old in unique):
            unique.append((u, float(v)))
    return unique


def _giese_fixed_total_state(baseline_diameter: float, diameter: float,
                             parameters: GieseParameters) -> tuple[float, float, float]:
    baseline_perimeter = math.pi * baseline_diameter
    baseline_area = math.pi * (baseline_diameter / 2.0) ** 2
    target_mass = (baseline_perimeter * parameters.u0
                   + baseline_area * parameters.v0)

    # Continue the equilibrium nearest the published state through geometry.
    prior = (parameters.u0, parameters.v0)
    steps = max(2, int(math.ceil(abs(diameter / baseline_diameter - 1.0) / 0.02)) + 1)
    for current in np.linspace(baseline_diameter, diameter, steps):
        perimeter = math.pi * current
        area = math.pi * (current / 2.0) ** 2
        roots = _giese_equilibrium_roots(parameters, perimeter, area, target_mass)
        if not roots:
            raise AdapterError("no positive homogeneous Giese equilibrium at fixed total")
        distances = np.asarray([
            ((u - prior[0]) / max(parameters.u0, 1e-12)) ** 2
            + ((v - prior[1]) / max(parameters.v0, 1e-12)) ** 2
            for u, v in roots
        ])
        order = np.argsort(distances)
        if len(order) > 1 and math.isclose(
            float(distances[order[0]]), float(distances[order[1]]), rel_tol=1e-8
        ):
            raise AdapterError("Giese fixed-total equilibrium branch is not unique")
        prior = roots[int(order[0])]
    return prior[0], prior[1], target_mass


def resolve_giese(config: Mapping[str, Any]) -> dict[str, Any]:
    values = _model_config(config)
    baseline_diameter = _positive_float(values, "baseline_diameter_um", 10.0)
    target_h = _positive_float(values, "target_h_um", 0.5)
    dt = _positive_float(values, "dt_s", 0.05)
    end_time = _positive_float(values, "end_time_s", 300.0)
    output_interval = _positive_float(values, "output_interval_s", 1.0)
    mode = _choice(values, "stimulus_mode", "two_site", {"none", "single", "two_site"})
    scheme = _choice(values, "scheme", "imex", {"imex", "giese_semi_implicit"})
    linear_scale, area_scale = geometry_scales(config["volume_ratio"])
    diameter = baseline_diameter * linear_scale
    parameters = GieseParameters()
    baseline_mass = (math.pi * baseline_diameter * parameters.u0
                     + math.pi * (baseline_diameter / 2.0) ** 2 * parameters.v0)
    if config["scaling_condition"] == "fixed_total":
        u0, v0, target_mass = _giese_fixed_total_state(
            baseline_diameter, diameter, parameters
        )
    else:
        u0, v0 = parameters.u0, parameters.v0
        target_mass = (math.pi * diameter * u0
                       + math.pi * (diameter / 2.0) ** 2 * v0)
    return {
        "model": "giese_wp",
        "volume_ratio": config["volume_ratio"],
        "linear_scale": linear_scale,
        "area_scale": area_scale,
        "baseline_diameter_um": baseline_diameter,
        "diameter_um": diameter,
        "target_h_um": target_h,
        "dt_s": dt,
        "end_time_s": end_time,
        "output_interval_s": output_interval,
        "stimulus_mode": mode,
        "scheme": scheme,
        "u0": u0,
        "v0": v0,
        "baseline_native_mass": baseline_mass,
        "target_native_mass": target_mass,
    }


def _giese_stimuli(mode: str) -> tuple[Stimulus, ...]:
    if mode == "none":
        return ()
    strong = Stimulus(0.44, 0.0, membrane_fraction=0.05, end_time=10.0)
    if mode == "single":
        return (strong,)
    return (strong, Stimulus(0.40, math.pi, membrane_fraction=0.05, end_time=10.0))


def simulate_giese(config: Mapping[str, Any], output: Path) -> Path:
    resolved = resolve_giese(config)
    _write_resolved(output, resolved)
    system = make_circular_system(
        resolved["diameter_um"], target_h=resolved["target_h_um"]
    )
    parameters = replace(
        GieseParameters(), u0=resolved["u0"], v0=resolved["v0"]
    )
    solver = GieseSolver(system, parameters)
    result = solver.run(
        end_time=resolved["end_time_s"],
        dt=resolved["dt_s"],
        stimuli=_giese_stimuli(resolved["stimulus_mode"]),
        output_interval=resolved["output_interval_s"],
        record_bulk=True,
        scheme=resolved["scheme"],
    )
    path = output / "trajectory.npz"
    np.savez_compressed(
        path,
        times=result.times,
        membrane=result.membrane,
        bulk=result.bulk,
        total_mass=result.total_mass,
        boundary_angles=system.boundary_angles,
        edge_lengths=system.edge_lengths,
        mesh_points=system.mesh.p,
        mesh_triangles=system.mesh.t,
    )
    return path


def analyze_giese(config: Mapping[str, Any], output: Path) -> dict[str, Any]:
    with np.load(output / "trajectory.npz", allow_pickle=False) as raw:
        times = raw["times"]
        membrane = raw["membrane"]
        bulk = raw["bulk"]
        total_mass = raw["total_mass"]
        edge_lengths = raw["edge_lengths"]
    metrics = [membrane_metrics(profile, edge_lengths) for profile in membrane]
    trajectory = summarize_trajectory(times, metrics)
    relative_mass_error = float(np.max(np.abs(
        total_mass / total_mass[0] - 1.0
    )))
    final = metrics[-1]
    membrane_mean = float(np.average(membrane[-1], weights=edge_lengths))
    membrane_variance = float(np.average(
        (membrane[-1]-membrane_mean)**2, weights=edge_lengths
    ))
    if trajectory["polarized_final"]:
        state = "multi_domain" if final["n_above_mean_clusters"] > 1 else "polarized"
    else:
        state = "homogeneous"
    status = "complete" if relative_mass_error <= 1e-9 else "quality_failure"
    summary = make_summary(
        config,
        run_status=status,
        status_message=("" if status == "complete"
                        else f"relative native-mass error {relative_mass_error:.3g}"),
        state_class=state,
        polarity_strength=float(final["relative_amplitude"]),
        n_domains=int(final["n_above_mean_clusters"]),
        time_to_polarity=float(trajectory["polarization_time_90_max_PF_s"]),
        model_metrics={
            **trajectory,
            "maximum_relative_mass_error": relative_mass_error,
            "initial_native_mass": float(total_mass[0]),
            "final_native_mass": float(total_mass[-1]),
            "final_membrane_mean": float(final["membrane_mean"]),
            "final_bulk_mean": float(np.mean(bulk[-1])),
            "spatial_contrast": float(final["relative_amplitude"]),
            "spatial_cv": math.sqrt(membrane_variance) / max(abs(membrane_mean), 1e-12),
            "conservation_error": relative_mass_error,
        },
        provenance={
            "manifest": "experiments/wgd_v1/manifests/giese_wp.md",
            "solver": "python/giese_wp",
        },
    )
    write_summary(summary, output)
    return summary


def run_giese(config: Mapping[str, Any], output: Path) -> dict[str, Any]:
    simulate_giese(config, output)
    return analyze_giese(config, output)


def resolve_taniguchi(config: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve the Fig. 4F deformable-cell homeostatic size experiment."""

    if config["scaling_condition"] == "fixed_total":
        raise AdapterError(
            "taniguchi_phase_field fixed_total is undefined because U and V "
            "have production, degradation, and interface-loss terms and do "
            "not form an independently conserved molecular pool"
        )
    if config["stimulus_condition"] != "relative_default":
        raise AdapterError(
            "taniguchi_phase_field currently supports relative_default only; "
            "its firing events are an intrinsic per-area stochastic process"
        )
    values = _model_config(config)
    parameter_set = _choice(values, "parameter_set", "4F", {"4F"})
    backend = _choice(values, "backend", "fast", {"fast", "reference"})
    stochastic = _boolean(values, "stochastic", True)
    if stochastic and config["seed"] is None:
        raise AdapterError("stochastic Taniguchi runs require an explicit seed")
    end_time = _positive_float(values, "end_time", 20.0)
    output_interval = _positive_float(values, "output_interval", 0.05)
    dx = _positive_float(values, "dx", 0.1)
    dt = _positive_float(values, "dt", 5e-4)
    baseline_box_half_width = _positive_float(
        values, "baseline_box_half_width", 10.0
    )
    activity_threshold = _positive_float(values, "activity_threshold", 0.5)
    minimum_domain_area = _positive_float(values, "minimum_domain_area", 0.2)
    checkpoint_interval = _positive_float(values, "checkpoint_interval", 2.0)
    linear_scale, area_scale = geometry_scales(config["volume_ratio"])
    baseline = PHASE_PARAMETER_SETS[parameter_set]
    baseline_radius = math.sqrt(baseline.area0 / math.pi)
    if baseline_box_half_width <= baseline_radius + 2.0 * baseline.epsilon:
        raise AdapterError(
            "model_config.baseline_box_half_width must leave at least two "
            "interface widths around the baseline cell"
        )
    return {
        "model": "taniguchi_phase_field",
        "volume_ratio": config["volume_ratio"],
        "linear_scale": linear_scale,
        "area_scale": area_scale,
        "parameter_set": parameter_set,
        "backend": backend,
        "stochastic": stochastic,
        "seed": config["seed"],
        "end_time": end_time,
        "output_interval": output_interval,
        "dx": dx,
        "dt": dt,
        "baseline_area0": baseline.area0,
        "area0": baseline.area0 * area_scale,
        "baseline_radius": baseline_radius,
        "initial_radius": baseline_radius * linear_scale,
        "baseline_box_half_width": baseline_box_half_width,
        "box_half_width": baseline_box_half_width * linear_scale,
        "activity_threshold": activity_threshold,
        "minimum_domain_area": minimum_domain_area,
        "checkpoint_interval": checkpoint_interval,
        "noise_rate_per_area": baseline.theta,
        "noise_width": baseline.noise_d,
        "noise_amplitude_mean": baseline.noise_sigma,
    }


def simulate_taniguchi(config: Mapping[str, Any], output: Path) -> Path:
    resolved = resolve_taniguchi(config)
    _write_resolved(output, resolved)
    parameters = replace(
        PHASE_PARAMETER_SETS[resolved["parameter_set"]],
        area0=resolved["area0"],
    )
    result = simulate_phase_field(
        resolved["end_time"],
        parameters=parameters,
        dx=resolved["dx"],
        dt=resolved["dt"],
        box_half_width=resolved["box_half_width"],
        output_interval=resolved["output_interval"],
        seed=(int(resolved["seed"]) if resolved["seed"] is not None else 1),
        stochastic=resolved["stochastic"],
        backend=resolved["backend"],
    )
    observables = extract_phase_wgd_observables(
        result,
        activity_threshold=resolved["activity_threshold"],
        minimum_domain_area=resolved["minimum_domain_area"],
    )
    path = output / "trajectory.npz"
    np.savez_compressed(
        path,
        format_version=np.asarray("taniguchi_wgd_compact_v1"),
        **observables,
        event_count=np.asarray(result.event_count, dtype=np.int64),
        area0=np.asarray(parameters.area0, dtype=float),
        parameter_set=np.asarray(resolved["parameter_set"]),
        backend=np.asarray(resolved["backend"]),
        dx=np.asarray(resolved["dx"], dtype=float),
        dt=np.asarray(resolved["dt"], dtype=float),
        activity_threshold=np.asarray(resolved["activity_threshold"], dtype=float),
        minimum_domain_area=np.asarray(resolved["minimum_domain_area"], dtype=float),
    )
    checkpoint_targets = np.arange(
        result.times[0], result.times[-1] + resolved["checkpoint_interval"],
        resolved["checkpoint_interval"],
    )
    checkpoint_indices = np.searchsorted(result.times, checkpoint_targets)
    checkpoint_indices = np.clip(checkpoint_indices, 0, result.times.size - 1)
    checkpoint_indices = np.unique(np.concatenate((
        np.asarray([0], dtype=int), checkpoint_indices,
        np.asarray([result.times.size - 1], dtype=int),
    )))
    np.savez_compressed(
        output / "checkpoints.npz",
        times=result.times[checkpoint_indices],
        phi=result.phi[checkpoint_indices],
        U=result.U[checkpoint_indices],
        V=result.V[checkpoint_indices],
        x=result.x,
        y=result.y,
        source_indices=checkpoint_indices,
        source_output_interval=np.asarray(resolved["output_interval"], dtype=float),
        checkpoint_interval=np.asarray(resolved["checkpoint_interval"], dtype=float),
    )
    return path


def analyze_taniguchi(config: Mapping[str, Any], output: Path) -> dict[str, Any]:
    resolved = resolve_taniguchi(config)
    path = output / "trajectory.npz"
    try:
        with np.load(path, allow_pickle=False) as raw:
            parameter_set = str(raw["parameter_set"])
            backend = str(raw["backend"])
            raw_area0 = float(raw["area0"])
            raw_dx = float(raw["dx"])
            raw_dt = float(raw["dt"])
            event_count = int(raw["event_count"])
            compact = (
                "format_version" in raw.files
                and str(raw["format_version"]) == "taniguchi_wgd_compact_v1"
            )
            if compact:
                raw_activity_threshold = float(raw["activity_threshold"])
                raw_minimum_domain_area = float(raw["minimum_domain_area"])
                observable_names = {
                    "times", "areas", "active_fractions", "spatial_contrasts",
                    "spatial_cvs", "domain_counts", "centroid_x", "centroid_y",
                    "shape_anisotropy", "directional_moments",
                    "directional_angles", "pip2_native_integrals",
                    "pip3_native_integrals", "box_edge_max_phi",
                }
                observables = {
                    key: np.asarray(raw[key]) for key in observable_names
                }
                result = None
            else:
                raw_activity_threshold = resolved["activity_threshold"]
                raw_minimum_domain_area = resolved["minimum_domain_area"]
                observables = None
                result = PhaseFieldResult(
                    x=np.asarray(raw["x"], dtype=float),
                    y=np.asarray(raw["y"], dtype=float),
                    parameters=replace(
                        PHASE_PARAMETER_SETS[parameter_set], area0=raw_area0
                    ),
                    times=np.asarray(raw["times"], dtype=float),
                    phi=np.asarray(raw["phi"]),
                    U=np.asarray(raw["U"]),
                    V=np.asarray(raw["V"]),
                    event_count=event_count,
                )
    except (OSError, KeyError, ValueError) as error:
        raise AdapterError(f"cannot read Taniguchi raw output {path}: {error}") from error
    if parameter_set != resolved["parameter_set"] or backend != resolved["backend"]:
        raise AdapterError("Taniguchi raw output does not match the resolved configuration")
    if not math.isclose(raw_area0, resolved["area0"], rel_tol=1e-12):
        raise AdapterError("Taniguchi raw target area does not match the requested geometry")
    if not math.isclose(raw_dx, resolved["dx"], rel_tol=1e-12) or not math.isclose(
        raw_dt, resolved["dt"], rel_tol=1e-12
    ):
        raise AdapterError("Taniguchi raw discretization does not match the configuration")
    if not math.isclose(
        raw_activity_threshold, resolved["activity_threshold"], rel_tol=1e-12
    ) or not math.isclose(
        raw_minimum_domain_area, resolved["minimum_domain_area"], rel_tol=1e-12
    ):
        raise AdapterError("Taniguchi raw activity definition does not match the configuration")

    if observables is not None:
        metrics = summarize_phase_wgd_observables(
            observables,
            area0=raw_area0,
            event_count=event_count,
            activity_threshold=resolved["activity_threshold"],
            minimum_domain_area=resolved["minimum_domain_area"],
        )
        raw_format = "taniguchi_wgd_compact_v1"
    else:
        assert result is not None
        metrics = summarize_phase_wgd(
            result,
            activity_threshold=resolved["activity_threshold"],
            minimum_domain_area=resolved["minimum_domain_area"],
        )
        raw_format = "legacy_full_fields"
    edge_phi = float(metrics["box_edge_max_phi"])
    status = "complete" if edge_phi <= 0.05 else "quality_failure"
    status_message = (
        "" if status == "complete" else
        f"diffuse cell reached the computational box (edge phi {edge_phi:.3g})"
    )
    model_metrics = {
        key: value for key, value in metrics.items()
        if key not in {"state_class", "terminal_n_domains", "time_to_directional_peak"}
    }
    model_metrics["spatial_contrast"] = metrics["terminal_spatial_contrast"]
    model_metrics["spatial_cv"] = metrics["terminal_spatial_cv"]
    model_metrics["native_integral"] = metrics["final_pip3_native_integral"]
    summary = make_summary(
        config,
        run_status=status,
        status_message=status_message,
        state_class=str(metrics["state_class"]),
        polarity_strength=float(metrics["terminal_spatial_contrast"]),
        n_domains=int(metrics["terminal_n_domains"]),
        time_to_polarity=(
            float(metrics["time_to_directional_peak"])
            if metrics["time_to_directional_peak"] is not None else None
        ),
        switching_or_oscillation_rate=float(metrics["domain_count_transition_rate"]),
        model_metrics=model_metrics,
        provenance={
            "manifest": "experiments/wgd_v1/manifests/taniguchi_phase_field.md",
            "solver": "python/taniguchi/phase_field.py",
            "metrics": "python/taniguchi/wgd_metrics.py",
            "backend": backend,
            "parameter_set": parameter_set,
            "raw_format": raw_format,
        },
    )
    write_summary(summary, output)
    return summary


def run_taniguchi(config: Mapping[str, Any], output: Path) -> dict[str, Any]:
    simulate_taniguchi(config, output)
    return analyze_taniguchi(config, output)


def resolve_config(config: Mapping[str, Any], repo_root: Path) -> dict[str, Any]:
    if config["model"] in ONE_D_MODELS:
        return resolve_one_d(config)
    if config["model"] == "debelly_one_d":
        return resolve_debelly_one_d(config, repo_root)
    if config["model"] == "chiou_mcas":
        return resolve_chiou(config, repo_root)
    if config["model"] == "giese_wp":
        return resolve_giese(config)
    if config["model"] == "eroume_maree":
        return resolve_eroume(config)
    if config["model"] == "taniguchi_phase_field":
        return resolve_taniguchi(config)
    raise AdapterError(f"no executable WGD v1 adapter for model {config['model']!r}")


def run_config(config: Mapping[str, Any], repo_root: Path) -> dict[str, Any]:
    output = Path(config["output_path"])
    if not output.is_absolute():
        output = repo_root / output
    if config["model"] in ONE_D_MODELS:
        simulate_one_d(config, output)
        return analyze_one_d(config, output)
    if config["model"] == "debelly_one_d":
        simulate_debelly_one_d(config, repo_root, output)
        return analyze_debelly_one_d(config, repo_root, output)
    if config["model"] == "chiou_mcas":
        return run_chiou(config, repo_root, output)
    if config["model"] == "giese_wp":
        return run_giese(config, output)
    if config["model"] == "eroume_maree":
        simulate_eroume(config, repo_root, output)
        return analyze_eroume(config, repo_root, output)
    if config["model"] == "taniguchi_phase_field":
        return run_taniguchi(config, output)
    raise AdapterError(f"no executable WGD v1 adapter for model {config['model']!r}")


def simulate_config(config: Mapping[str, Any], repo_root: Path) -> Path:
    output = Path(config["output_path"])
    if not output.is_absolute():
        output = repo_root / output
    if config["model"] in ONE_D_MODELS:
        return simulate_one_d(config, output)
    if config["model"] == "debelly_one_d":
        return simulate_debelly_one_d(config, repo_root, output)
    if config["model"] == "chiou_mcas":
        return simulate_chiou(config, repo_root, output)
    if config["model"] == "giese_wp":
        return simulate_giese(config, output)
    if config["model"] == "eroume_maree":
        return simulate_eroume(config, repo_root, output)
    if config["model"] == "taniguchi_phase_field":
        return simulate_taniguchi(config, output)
    raise AdapterError(
        f"separate simulation adapter is not available for model {config['model']!r}"
    )


def analyze_config(config: Mapping[str, Any], repo_root: Path) -> dict[str, Any]:
    output = Path(config["output_path"])
    if not output.is_absolute():
        output = repo_root / output
    if config["model"] in ONE_D_MODELS:
        return analyze_one_d(config, output)
    if config["model"] == "debelly_one_d":
        return analyze_debelly_one_d(config, repo_root, output)
    if config["model"] == "chiou_mcas":
        return analyze_chiou(config, repo_root, output)
    if config["model"] == "giese_wp":
        return analyze_giese(config, output)
    if config["model"] == "eroume_maree":
        return analyze_eroume(config, repo_root, output)
    if config["model"] == "taniguchi_phase_field":
        return analyze_taniguchi(config, output)
    raise AdapterError(
        f"separate output adapter is not available for model {config['model']!r}"
    )
