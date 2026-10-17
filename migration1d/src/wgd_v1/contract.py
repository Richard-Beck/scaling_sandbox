"""Minimal runtime checks for the WGD v1 JSON adapter contract.

The JSON schemas remain the serialized specification.  These dependency-free
checks cover the invariants adapters need at runtime on the HPC system.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any, Mapping


PROTOCOL_VERSION = "wgd_v1"
SCALING_CONDITIONS = {
    "wgd_homeostatic",
    "fixed_total",
    "native_concentration",
}
RUN_STATUSES = {
    "complete",
    "invalid_config",
    "numerical_failure",
    "quality_failure",
}
CONFIG_REQUIRED = {
    "protocol_version",
    "model",
    "volume_ratio",
    "scaling_condition",
    "stimulus_condition",
    "seed",
    "output_path",
}
CONFIG_ALLOWED = CONFIG_REQUIRED | {"model_config"}
SUMMARY_REQUIRED = {
    "protocol_version",
    "model_id",
    "volume_ratio",
    "linear_scale",
    "area_scale",
    "scaling_condition",
    "stimulus_condition",
    "seed",
    "state_class",
    "polarity_strength",
    "n_domains",
    "time_to_polarity",
    "switching_or_oscillation_rate",
    "run_status",
    "status_message",
    "model_metrics",
    "provenance",
}


class ContractError(ValueError):
    """A configuration or summary violates the WGD v1 contract."""


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def geometry_scales(volume_ratio: float) -> tuple[float, float]:
    """Return the linear and area scale for a positive volume ratio."""

    _require(_is_number(volume_ratio), "volume_ratio must be numeric")
    _require(math.isfinite(volume_ratio) and volume_ratio > 0,
             "volume_ratio must be finite and positive")
    linear_scale = float(volume_ratio) ** (1.0 / 3.0)
    return linear_scale, linear_scale**2


def validate_config(config: Mapping[str, Any], repo_root: Path | None = None) -> None:
    """Validate the common configuration fields and model manifest."""

    _require(isinstance(config, Mapping), "config must be a JSON object")
    missing = CONFIG_REQUIRED - set(config)
    extra = set(config) - CONFIG_ALLOWED
    _require(not missing, f"config is missing fields: {sorted(missing)}")
    _require(not extra, f"config has unknown fields: {sorted(extra)}")
    _require(config["protocol_version"] == PROTOCOL_VERSION,
             f"protocol_version must be {PROTOCOL_VERSION!r}")
    _require(isinstance(config["model"], str) and config["model"],
             "model must be a non-empty string")
    geometry_scales(config["volume_ratio"])
    _require(config["scaling_condition"] in SCALING_CONDITIONS,
             "unknown scaling_condition")
    _require(isinstance(config["stimulus_condition"], str)
             and bool(config["stimulus_condition"]),
             "stimulus_condition must be a non-empty string")
    seed = config["seed"]
    _require(seed is None or (isinstance(seed, int) and not isinstance(seed, bool)),
             "seed must be an integer or null")
    _require(isinstance(config["output_path"], str) and config["output_path"],
             "output_path must be a non-empty string")
    _require("model_config" not in config or isinstance(config["model_config"], Mapping),
             "model_config must be an object")

    root = Path(repo_root) if repo_root is not None else _repo_root()
    manifest = root / "experiments" / "wgd_v1" / "manifests" / f"{config['model']}.md"
    _require(manifest.is_file(), f"no WGD v1 manifest for model {config['model']!r}")


def make_summary(
    config: Mapping[str, Any],
    *,
    run_status: str,
    status_message: str = "",
    state_class: str | None = None,
    polarity_strength: float | None = None,
    n_domains: int | None = None,
    time_to_polarity: float | None = None,
    switching_or_oscillation_rate: float | None = None,
    model_metrics: Mapping[str, Any] | None = None,
    provenance: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Construct and validate a terminal summary for a valid config."""

    validate_config(config)
    linear_scale, area_scale = geometry_scales(config["volume_ratio"])
    manifest = f"experiments/wgd_v1/manifests/{config['model']}.md"
    provenance_out = dict(provenance or {})
    provenance_out.setdefault("manifest", manifest)
    summary = {
        "protocol_version": PROTOCOL_VERSION,
        "model_id": config["model"],
        "volume_ratio": config["volume_ratio"],
        "linear_scale": linear_scale,
        "area_scale": area_scale,
        "scaling_condition": config["scaling_condition"],
        "stimulus_condition": config["stimulus_condition"],
        "seed": config["seed"],
        "state_class": state_class,
        "polarity_strength": polarity_strength,
        "n_domains": n_domains,
        "time_to_polarity": time_to_polarity,
        "switching_or_oscillation_rate": switching_or_oscillation_rate,
        "run_status": run_status,
        "status_message": status_message,
        "model_metrics": dict(model_metrics or {}),
        "provenance": provenance_out,
    }
    validate_summary(summary, config=config)
    return summary


def _optional_nonnegative_number(summary: Mapping[str, Any], field: str) -> None:
    value = summary[field]
    _require(value is None or (_is_number(value) and math.isfinite(value) and value >= 0),
             f"{field} must be a finite nonnegative number or null")


def validate_summary(
    summary: Mapping[str, Any], config: Mapping[str, Any] | None = None
) -> None:
    """Validate required summary fields and geometry/config consistency."""

    _require(isinstance(summary, Mapping), "summary must be a JSON object")
    missing = SUMMARY_REQUIRED - set(summary)
    extra = set(summary) - SUMMARY_REQUIRED
    _require(not missing, f"summary is missing fields: {sorted(missing)}")
    _require(not extra, f"summary has unknown fields: {sorted(extra)}")
    _require(summary["protocol_version"] == PROTOCOL_VERSION,
             f"protocol_version must be {PROTOCOL_VERSION!r}")
    _require(isinstance(summary["model_id"], str) and summary["model_id"],
             "model_id must be a non-empty string")
    linear_scale, area_scale = geometry_scales(summary["volume_ratio"])
    _require(_is_number(summary["linear_scale"])
             and math.isclose(summary["linear_scale"], linear_scale, rel_tol=1e-12),
             "linear_scale is inconsistent with volume_ratio")
    _require(_is_number(summary["area_scale"])
             and math.isclose(summary["area_scale"], area_scale, rel_tol=1e-12),
             "area_scale is inconsistent with volume_ratio")
    _require(summary["scaling_condition"] in SCALING_CONDITIONS,
             "unknown scaling_condition")
    _require(isinstance(summary["stimulus_condition"], str)
             and bool(summary["stimulus_condition"]),
             "stimulus_condition must be a non-empty string")
    _require(summary["seed"] is None or (
        isinstance(summary["seed"], int) and not isinstance(summary["seed"], bool)
    ), "seed must be an integer or null")
    _require(summary["state_class"] is None or isinstance(summary["state_class"], str),
             "state_class must be a string or null")
    _optional_nonnegative_number(summary, "polarity_strength")
    n_domains = summary["n_domains"]
    _require(n_domains is None or (
        isinstance(n_domains, int) and not isinstance(n_domains, bool) and n_domains >= 0
    ), "n_domains must be a nonnegative integer or null")
    _optional_nonnegative_number(summary, "time_to_polarity")
    _optional_nonnegative_number(summary, "switching_or_oscillation_rate")
    _require(summary["run_status"] in RUN_STATUSES, "unknown run_status")
    _require(isinstance(summary["status_message"], str),
             "status_message must be a string")
    _require(isinstance(summary["model_metrics"], Mapping),
             "model_metrics must be an object")
    provenance = summary["provenance"]
    _require(isinstance(provenance, Mapping), "provenance must be an object")
    _require(isinstance(provenance.get("manifest"), str) and provenance["manifest"],
             "provenance.manifest must be a non-empty string")

    if config is not None:
        validate_config(config)
        pairs = {
            "model_id": "model",
            "volume_ratio": "volume_ratio",
            "scaling_condition": "scaling_condition",
            "stimulus_condition": "stimulus_condition",
            "seed": "seed",
        }
        for summary_key, config_key in pairs.items():
            _require(summary[summary_key] == config[config_key],
                     f"summary {summary_key} does not match config {config_key}")


def write_summary(summary: Mapping[str, Any], output_path: str | Path) -> Path:
    """Validate and atomically write ``summary.json`` in an output directory."""

    validate_summary(summary)
    directory = Path(output_path)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "summary.json"
    temporary = directory / f".summary.{os.getpid()}.tmp"
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    temporary.replace(target)
    return target
