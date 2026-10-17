#!/usr/bin/env python3
"""Build A/C/E-only L=15 kymographs for review before a full Figure 2 sweep."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import (  # noqa: E402
    DeBelly, DeBellyParameters, FinitePulseSequence, LocalizedPulse, simulate,
)
from polarity1d.metrics import profile_domain_count  # noqa: E402


LENGTH = 15.0
CELLS = 75
STEP = 0.5
ALPHA0 = 1.8
BETA0 = 1.8
FIELDS = (
    "cue", "rac", "rho", "tension", "mca_bound", "mca_unbound",
    "delta_beta_tension", "delta_alpha_mca",
)
LABELS = {
    "cue": "Cue", "rac": "Rac", "rho": "Rho", "tension": "Tension",
    "mca_bound": "MCA bound", "mca_unbound": "MCA unbound",
    "delta_beta_tension": r"$\Delta\beta_T$",
    "delta_alpha_mca": r"$\Delta\alpha_{MCA}$",
}


def _model():
    return DeBelly(DeBellyParameters(
        rac_basal_activation=ALPHA0,
        rho_basal_activation=BETA0,
        rac_mechanical_feedback=0.3 * ALPHA0,
        rho_mechanical_feedback=0.3 * BETA0,
        assay_mode="reflection_symmetric_extension",
        initialization_mode="algebraic_fixed_point",
        input_sampling_mode="cell_centred",
    ))


def _pulse(side, start, duration, ratio, width):
    return {
        "amplitude": ratio * ALPHA0,
        "start": start,
        "stop": start + duration,
        "taper": duration / 2.0,
        "width_fraction": width,
        "sides": [side],
    }


def _orientation(values, x):
    left = np.mean(values[:, x <= 0.2], axis=1)
    right = np.mean(values[:, x >= 0.8], axis=1)
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)


def _contrast(values):
    high, low = np.max(values, axis=1), np.min(values, axis=1)
    return (high - low) / np.maximum(np.abs(high) + np.abs(low), 1e-12)


def _first_sustained(time, condition, *, after, before=None, duration=50.0):
    start = int(np.searchsorted(time, after, side="left"))
    limit = len(time) if before is None else int(np.searchsorted(time, before, side="left"))
    for index in range(start, limit):
        stop = int(np.searchsorted(time, time[index] + duration, side="left"))
        if stop < limit and np.all(condition[index : stop + 1]):
            return float(time[index])
    return None


def _simulate(protocol, end_time):
    return simulate(
        _model(), protocol=protocol, length=LENGTH, cells=CELLS,
        end_time=end_time, output_interval=2.0, step=STEP,
    )


def _payload(trajectory, protocol):
    fields = {
        "cue": np.asarray([
            protocol(float(time), trajectory.x, LENGTH)
            for time in trajectory.time
        ])
    }
    fields.update({field: trajectory.fields[field] for field in FIELDS[1:]})
    return fields


def _competition():
    protocol = LocalizedPulse(
        amplitude=2.0 * ALPHA0, start=0.0, stop=80.0, taper=40.0,
        width_fraction=0.2, sides=("left", "right"), asymmetry=0.15,
    )
    trajectory = _simulate(protocol, 1_500.0)
    x = trajectory.x / LENGTH
    rac = trajectory.fields["rac"]
    domains = np.asarray([profile_domain_count(profile) for profile in rac])
    orientation = _orientation(rac, x)
    contrast = _contrast(rac)
    left = np.mean(rac[:, x <= 0.15], axis=1)
    right = np.mean(rac[:, x >= 0.85], axis=1)
    centre = np.mean(rac[:, (x >= 0.4) & (x <= 0.6)], axis=1)
    two_pole_strength = np.minimum(
        (left - centre) / np.maximum(np.abs(left) + np.abs(centre), 1e-12),
        (right - centre) / np.maximum(np.abs(right) + np.abs(centre), 1e-12),
    )
    two = (domains >= 2) & (two_pole_strength >= 0.05)
    resolved = (domains == 1) & (contrast >= 0.3) & (np.abs(orientation) >= 0.25)
    two_time = _first_sustained(
        trajectory.time, two, after=0.0, before=250.0, duration=20.0
    )
    resolution = _first_sustained(
        trajectory.time, resolved, after=80.0, duration=100.0
    )
    final = rac[-1]
    low, amplitude = float(np.min(final)), float(np.max(final) - np.min(final))
    half_width = float(np.trapezoid(
        (final >= low + 0.5 * amplitude).astype(float), x
    ))
    core_width = float(np.trapezoid(
        (final >= low + 0.8 * amplitude).astype(float), x
    ))
    metrics = {
        "two_pole_establishment_s": two_time,
        "peak_two_pole_strength": float(np.max(two_pole_strength[trajectory.time <= 230.0])),
        "strong_single_front_s": resolution,
        "final_rac_contrast": float(contrast[-1]),
        "final_rac_orientation": float(orientation[-1]),
        "final_rac_domain_count": int(domains[-1]),
        "final_front_half_height_width_L": half_width,
        "final_front_80pct_width_L": core_width,
    }
    return protocol, trajectory, metrics


def _reversal():
    pulses = (
        _pulse("left", 0.0, 160.0, 2.0, 0.2),
        _pulse("right", 600.0, 320.0, 4.0, 0.4),
    )
    protocol = FinitePulseSequence(pulses=pulses)
    trajectory = _simulate(protocol, 1_800.0)
    x = trajectory.x / LENGTH
    rac = _orientation(trajectory.fields["rac"], x)
    rho = _orientation(trajectory.fields["rho"], x)
    contrast = _contrast(trajectory.fields["rac"])
    before = int(np.argmin(np.abs(trajectory.time - 600.0)))
    target = (rac >= 0.25) & (rho <= -0.25) & (contrast >= 0.3)
    reversal = _first_sustained(
        trajectory.time, target, after=920.0, duration=100.0
    )
    metrics = {
        "pre_opposite_rac_directed": float(-rac[before]),
        "pre_opposite_rho_directed": float(rho[before]),
        "pre_opposite_rac_contrast": float(contrast[before]),
        "strong_reversal_s": reversal,
        "latency_after_opposite_start_s": None if reversal is None else reversal - 600.0,
        "latency_after_withdrawal_s": None if reversal is None else reversal - 920.0,
        "final_rac_directed": float(rac[-1]),
        "final_rho_directed": float(-rho[-1]),
        "final_rac_contrast": float(contrast[-1]),
    }
    return protocol, trajectory, metrics


def _train(period):
    half_period = period / 2.0
    starts = [0.0] + [600.0 + index * half_period for index in range(10)]
    sides = ["left"] + ["right" if index % 2 == 0 else "left" for index in range(10)]
    pulses = [_pulse("left", 0.0, 160.0, 2.0, 0.2)]
    pulses.extend(
        {
            "amplitude": 4.0 * ALPHA0,
            "start": start,
            "stop": start + half_period,
            "taper": 0.0,
            "width_fraction": 0.4,
            "sides": [side],
        }
        for side, start in zip(sides[1:], starts[1:])
    )
    protocol = FinitePulseSequence(pulses=tuple(pulses))
    end_time = starts[-1] + half_period
    trajectory = _simulate(protocol, end_time)
    x = trajectory.x / LENGTH
    rac = _orientation(trajectory.fields["rac"], x)
    rho = _orientation(trajectory.fields["rho"], x)
    contrast = _contrast(trajectory.fields["rac"])
    event_rows = []
    for index, (side, start) in enumerate(zip(sides[1:], starts[1:]), start=1):
        next_start = starts[index + 1] if index + 1 < len(starts) else end_time
        sign = 1.0 if side == "right" else -1.0
        target = (
            (sign * rac >= 0.25)
            & (-sign * rho >= 0.25)
            & (contrast >= 0.3)
        )
        establishment = _first_sustained(
            trajectory.time, target, after=start, before=next_start,
            duration=min(20.0, half_period / 4.0),
        )
        event_rows.append({
            "period_s": period,
            "event": index,
            "side": side,
            "establishment_s": establishment,
            "latency_s": None if establishment is None else establishment - start,
            "tracked": establishment is not None,
            "peak_directed_rac": float(np.max(sign * rac[
                (trajectory.time >= start) & (trajectory.time < next_start)
            ])),
        })
    return protocol, trajectory, event_rows


def main():
    output = ROOT / "results/tmp/trial/debelly_l15_protocol_design"
    output.mkdir(parents=True, exist_ok=True)
    competition = _competition()
    reversal = _reversal()
    periods = (100.0, 200.0, 400.0, 800.0, 1_600.0)
    train_results = [_train(period) for period in periods]
    event_rows = [row for _, _, rows in train_results for row in rows]
    with (output / "frequency_screen.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(event_rows[0]))
        writer.writeheader()
        writer.writerows(event_rows)
    summaries = []
    for period, result in zip(periods, train_results):
        rows = result[2]
        summaries.append({
            "period_s": period,
            "tracked_events": sum(row["tracked"] for row in rows),
            "events": len(rows),
            "mean_latency_s": float(np.mean([
                row["latency_s"] for row in rows if row["latency_s"] is not None
            ])) if any(row["tracked"] for row in rows) else None,
        })
    # Prefer the shortest train that tracks every event. If none does, show the
    # train with the most tracked events and then the shortest period.
    fully_tracked = [item for item in summaries if item["tracked_events"] == item["events"]]
    selected_summary = min(
        fully_tracked or summaries,
        key=lambda item: (
            0 if item["tracked_events"] == item["events"] else -item["tracked_events"],
            item["period_s"],
        ),
    )
    selected_index = [item["period_s"] for item in summaries].index(selected_summary["period_s"])
    train = train_results[selected_index]

    entries = [
        ("A", "Two competing end cues", *competition[:2], competition[2]),
        ("C", "Induce → stabilize → opposite pulse", *reversal[:2], reversal[2]),
        (
            "E", f"Alternating left/right train (cycle {selected_summary['period_s']:g} s)",
            train[0], train[1], selected_summary,
        ),
    ]
    payloads = [_payload(trajectory, protocol) for _, _, protocol, trajectory, _ in entries]

    cmaps = {
        "cue": "Greys", "rac": "viridis", "rho": "magma",
        "tension": "coolwarm", "mca_bound": "cividis",
        "mca_unbound": "plasma", "delta_beta_tension": "Reds",
        "delta_alpha_mca": "Blues",
    }
    norms = {}
    for field in FIELDS[1:]:
        values = np.concatenate([payload[field].ravel() for payload in payloads])
        if field == "tension":
            limit = max(abs(float(np.min(values))), abs(float(np.max(values))), 1e-12)
            norms[field] = mpl.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
        else:
            low, high = float(np.min(values)), float(np.max(values))
            norms[field] = mpl.colors.Normalize(low, high if high > low else low + 1e-12)

    mpl.rcParams.update({"font.size": 7, "axes.linewidth": 0.7})
    figure, axes = plt.subplots(3, len(FIELDS), figsize=(23, 8.2), constrained_layout=True)
    images = {}
    for row_index, ((label, title, protocol, trajectory, metrics), payload) in enumerate(zip(entries, payloads)):
        for column, field in enumerate(FIELDS):
            axis = axes[row_index, column]
            norm = (
                mpl.colors.Normalize(0.0, max(float(np.max(payload["cue"])), 1e-12))
                if field == "cue" else norms[field]
            )
            images[field] = axis.imshow(
                payload[field].T, origin="lower", aspect="auto", interpolation="nearest",
                extent=(trajectory.time[0], trajectory.time[-1], 0.0, 1.0),
                cmap=cmaps[field], norm=norm,
            )
            axis.set_xlabel("Time (s)")
            if column == 0:
                axis.set_ylabel("Relative position")
                axis.text(
                    -0.38, 1.13, f"{label}  {title}", transform=axis.transAxes,
                    fontsize=9, fontweight="bold", va="top",
                )
            else:
                axis.set_yticklabels([])
            if row_index == 0:
                axis.set_title(LABELS[field], fontsize=8)
            axis.tick_params(length=2, pad=1)
        axes[row_index, -1].text(
            1.02, 0.5, "\n".join(
                f"{key.replace('_', ' ')}: {value:.3g}" if isinstance(value, float)
                else f"{key.replace('_', ' ')}: {value}"
                for key, value in metrics.items()
            ),
            transform=axes[row_index, -1].transAxes, va="center", fontsize=6,
        )
    for column, field in enumerate(FIELDS):
        color_axis = axes[0, column].inset_axes((0.03, -0.24, 0.94, 0.035))
        colorbar = figure.colorbar(images[field], cax=color_axis, orientation="horizontal")
        colorbar.ax.tick_params(labelsize=5, length=1, pad=1)
    figure.suptitle(
        r"L=15 µm protocol preview — reflection-symmetric 2026 de Belly extension "
        r"($\alpha_0=\beta_0=1.8$, $\alpha=\beta=0.3\alpha_0$)",
        fontsize=12, fontweight="bold",
    )
    destination = output / "Figure2_ACE_L15_protocol_preview.png"
    figure.savefig(destination, dpi=300)
    figure.savefig(destination.with_suffix(".pdf"))
    plt.close(figure)

    np.savez_compressed(
        output / "Figure2_ACE_L15_preview_trajectories.npz",
        **{
            f"{label}__{key}": value
            for (label, _, _, trajectory, _), payload in zip(entries, payloads)
            for key, value in {"time": trajectory.time, "x": trajectory.x, **payload}.items()
        },
    )
    report = {
        "model": "reflection-symmetric extension of the 2026 de Belly mechanochemical model",
        "length_um": LENGTH,
        "parameters": {
            "alpha0": ALPHA0, "beta0": BETA0,
            "alpha_mechanical": 0.3 * ALPHA0,
            "beta_mechanical": 0.3 * BETA0,
        },
        "competition": competition[2],
        "reversal": reversal[2],
        "frequency_screen": summaries,
        "selected_frequency_train": selected_summary,
    }
    (output / "Figure2_ACE_L15_preview_metrics.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))
    print(destination)


if __name__ == "__main__":
    main()
