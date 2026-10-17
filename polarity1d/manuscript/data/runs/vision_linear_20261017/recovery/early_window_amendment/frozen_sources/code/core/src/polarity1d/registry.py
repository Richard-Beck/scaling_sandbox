"""Stable names used by configs and command-line runs."""

from __future__ import annotations

from dataclasses import fields
from typing import Any

from .debelly import DeBelly, DeBellyParameters
from .holmes import HolmesModel3
from .jilkine import Goryachev, LEGI, Otsuji, WavePinning


MODEL_REGISTRY = {
    "wave_pinning": WavePinning,
    "otsuji": Otsuji,
    "goryachev": Goryachev,
    "legi": LEGI,
    "holmes_model3": HolmesModel3,
    "debelly": DeBelly,
    "debelly_2026": DeBelly,
}


def make_model(name: str, parameters: dict[str, Any] | None = None):
    parameters = {} if parameters is None else dict(parameters)
    try:
        constructor = MODEL_REGISTRY[name]
    except KeyError as error:
        raise ValueError(f"unknown model {name!r}; choose from {sorted(MODEL_REGISTRY)}") from error
    if constructor is DeBelly:
        valid = {item.name for item in fields(DeBellyParameters)}
        unknown = set(parameters) - valid
        if unknown:
            raise ValueError(f"unknown De Belly parameters: {sorted(unknown)}")
        return DeBelly(DeBellyParameters(**parameters))
    return constructor(**parameters)
