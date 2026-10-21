"""Declared shape/timing controls for local-native signature assays.

Area means spatially integrated commanded activation input during ON time.
It is neither an optical calibration nor an equal accumulated molecular dose.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math
import numpy as np


def integration_weights(x, length, grid_kind="cell"):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError("A finite one-dimensional native grid is required")
    if not math.isfinite(length) or length <= 0:
        raise ValueError("A finite positive domain length is required")
    if grid_kind == "cell":
        expected = (np.arange(len(x)) + .5) * length / len(x)
        weights = np.full(len(x), length / len(x))
    elif grid_kind == "node":
        expected = np.linspace(0., length, len(x))
        weights = np.full(len(x), length / (len(x) - 1))
        weights[[0, -1]] *= .5
    else:
        raise ValueError("grid_kind must be cell or node")
    if not np.allclose(x, expected, rtol=1e-9, atol=1e-10):
        raise ValueError("Input coordinates do not match the declared native grid")
    return weights


class AreaCue:
    """Native raised cosine, normalized to exact native quadrature area.

    For Spring, normalization uses the instantaneous physical domain. The
    commanded area therefore remains fixed as its native boundaries move.
    Domain truncation is reported and is never hidden by a nominal width label.
    """
    event_times = ()

    def __init__(self, area, width, grid_kind="cell", side="right", background=0.):
        self.area, self.width = float(area), float(width)
        self.grid_kind, self.side = grid_kind, side
        self.background = float(background)
        if (not all(math.isfinite(v) for v in
                    (self.area, self.width, self.background))
                or self.area < 0 or self.width <= 0 or self.background < 0
                or side not in ("left", "right")):
            raise ValueError("Invalid area-cue parameters")

    def profile(self, x, length):
        x = np.asarray(x, dtype=float)
        weights = integration_weights(x, length, self.grid_kind)
        if self.area == 0:
            return np.zeros_like(x), weights, 0.
        distance = length - x if self.side == "right" else x
        scaled = distance / self.width
        raw = np.where((scaled >= 0) & (scaled <= 1),
                       .5 * (1 + np.cos(np.pi * np.clip(scaled, 0, 1))), 0.)
        raw_area = float(raw @ weights)
        if raw_area <= 0:
            raise ValueError("Native grid does not resolve the cue footprint")
        multiplier = self.area / raw_area
        return raw * multiplier, weights, multiplier

    def __call__(self, time, x, length):
        del time
        return self.background + self.profile(x, length)[0]

    def profile_diagnostics(self, x, length):
        profile, weights, multiplier = self.profile(x, length)
        return dict(requested_width=self.width, requested_area=self.area,
                    sampled_area=float(profile @ weights),
                    sampled_peak=float(profile.max()),
                    reconstructed_boundary_peak=multiplier,
                    nominal_untruncated_continuum_peak=2 * self.area / self.width,
                    support_truncated=bool(self.width > length),
                    support_cells=int(np.count_nonzero(profile)),
                    grid_kind=self.grid_kind, actual_domain_length=float(length),
                    area_definition="Spatially integrated commanded input per ON time; native quadrature")


@dataclass(frozen=True)
class TrainSchedule:
    pulse_count: int
    requested_total_on: float
    actual_total_on: float
    pulse_duration: float
    requested_gap_multiple: float
    actual_gap_duration: float
    actual_gap_multiple: float
    native_step: float | None
    segments: tuple[tuple[str, float, float], ...]

    @property
    def end_time(self):
        return self.segments[-1][2]

    @property
    def total_on(self):
        return self.actual_total_on

    @property
    def on_intervals(self):
        return tuple((start, stop) for kind, start, stop in self.segments if kind == "on")

    @property
    def event_times(self):
        return tuple(sorted({t for _, start, stop in self.segments for t in (start, stop)}))

    def as_dict(self):
        return asdict(self)


def make_train(calibration_duration, total_on_factor, pulse_count, gap_multiple,
               native_step=None, pulse_counts=(2, 4, 8)):
    """Preserve equal realized total ON across pulse counts on native lattices.

    deBelly's total ON is rounded UP to a common multiple allowing every pulse
    count to use its half-unit step. Gaps are rounded to the nearest step (ties
    upward). Both requested and realized exposure/timing are always retained.
    There are n-1 inter-pulse gaps; the final pulse is followed by native recovery.
    """
    if (pulse_count not in pulse_counts or pulse_count < 1
            or not all(math.isfinite(v) for v in
                       (calibration_duration, total_on_factor, gap_multiple))
            or calibration_duration <= 0 or total_on_factor <= 0 or gap_multiple < 0):
        raise ValueError("Invalid cue-train factors")
    requested = float(calibration_duration * total_on_factor)
    actual = requested
    if native_step is not None:
        if not math.isfinite(native_step) or native_step <= 0:
            raise ValueError("Invalid native lattice spacing")
        quantum = math.lcm(*pulse_counts) * native_step
        actual = math.ceil(requested / quantum - 1e-12) * quantum
    pulse = actual / pulse_count
    gap = pulse * gap_multiple
    if native_step is not None:
        gap = math.floor(gap / native_step + .5 + 1e-12) * native_step
    now = 0.
    segments = []
    for i in range(pulse_count):
        segments.append(("on", now, now + pulse))
        now += pulse
        if i < pulse_count - 1 and gap > 0:
            segments.append(("gap", now, now + gap))
            now += gap
    return TrainSchedule(
        pulse_count, requested, actual, pulse, gap_multiple, gap,
        gap / pulse, native_step, tuple(segments))
