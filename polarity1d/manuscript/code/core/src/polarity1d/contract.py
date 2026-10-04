"""Small, model-independent contracts used by every WP1 simulation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Protocol, runtime_checkable

import numpy as np


@dataclass(frozen=True)
class Grid1D:
    """A conservative cell-centred finite-volume grid on ``[0, length]``."""

    length: float = 10.0
    cells: int = 101

    def __post_init__(self) -> None:
        if not np.isfinite(self.length) or self.length <= 0:
            raise ValueError("length must be finite and positive")
        if isinstance(self.cells, bool) or self.cells < 5:
            raise ValueError("cells must be an integer of at least five")

    @property
    def dx(self) -> float:
        return self.length / self.cells

    @property
    def x(self) -> np.ndarray:
        return (np.arange(self.cells, dtype=float) + 0.5) * self.dx


@runtime_checkable
class Stimulus(Protocol):
    """A physical activation-rate profile shared by all model adapters."""

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        ...

    @property
    def event_times(self) -> tuple[float, ...]:
        ...


@runtime_checkable
class ReactionDiffusionModel(Protocol):
    """Contract for models solved by the common reaction--diffusion engines."""

    name: str
    state_names: tuple[str, ...]
    diffusivities: tuple[float, ...]
    observable_name: str
    preferred_step: float
    nonnegative: bool

    def initial_state(self, grid: Grid1D) -> np.ndarray:
        ...

    def reaction(
        self,
        time: float,
        x: np.ndarray,
        state: np.ndarray,
        stimulus: np.ndarray,
    ) -> np.ndarray:
        ...

    def metadata(self) -> Mapping[str, Any]:
        ...


@dataclass
class Trajectory:
    """Tidy in-memory result returned by both biochemical and mechanical models."""

    model: str
    time: np.ndarray
    x: np.ndarray
    fields: dict[str, np.ndarray]
    observable_name: str
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.time = np.asarray(self.time, dtype=float)
        self.x = np.asarray(self.x, dtype=float)
        if self.time.ndim != 1 or self.x.ndim != 1:
            raise ValueError("time and x must be one-dimensional")
        if len(self.time) < 1 or np.any(np.diff(self.time) <= 0):
            raise ValueError("time must be nonempty and strictly increasing")
        if len(self.x) < 2 or np.any(np.diff(self.x) <= 0):
            raise ValueError("x must contain at least two increasing coordinates")
        expected = (len(self.time), len(self.x))
        for name, values in tuple(self.fields.items()):
            array = np.asarray(values, dtype=float)
            if array.shape != expected:
                raise ValueError(f"field {name!r} has shape {array.shape}, expected {expected}")
            if not np.all(np.isfinite(array)):
                raise FloatingPointError(f"field {name!r} contains non-finite values")
            self.fields[name] = array
        if self.observable_name not in self.fields:
            raise ValueError(f"observable field {self.observable_name!r} is absent")

    @property
    def observable(self) -> np.ndarray:
        return self.fields[self.observable_name]

    def save(self, path: str | Path) -> Path:
        """Save arrays without pickles; metadata is stored as a JSON string."""
        import json

        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        payload: dict[str, Any] = {
            "model": np.asarray(self.model),
            "time": self.time,
            "x": self.x,
            "observable_name": np.asarray(self.observable_name),
            "metadata_json": np.asarray(json.dumps(self.metadata, sort_keys=True)),
        }
        payload.update({f"field__{name}": value for name, value in self.fields.items()})
        np.savez_compressed(destination, **payload)
        return destination

    @classmethod
    def load(cls, path: str | Path) -> "Trajectory":
        import json

        with np.load(path, allow_pickle=False) as raw:
            fields = {
                key.removeprefix("field__"): np.asarray(raw[key])
                for key in raw.files
                if key.startswith("field__")
            }
            return cls(
                model=str(raw["model"]),
                time=np.asarray(raw["time"]),
                x=np.asarray(raw["x"]),
                fields=fields,
                observable_name=str(raw["observable_name"]),
                metadata=json.loads(str(raw["metadata_json"])),
            )


def output_times(end_time: float, interval: float, *, start: float = 0.0) -> np.ndarray:
    """Return an increasing output grid that always includes ``end_time``."""
    if end_time <= start or interval <= 0:
        raise ValueError("require end_time > start and interval > 0")
    times = np.arange(start, end_time + 0.5 * interval, interval, dtype=float)
    times = times[times <= end_time]
    if not np.isclose(times[-1], end_time):
        times = np.append(times, end_time)
    return times

