#!/usr/bin/env python3
"""Quick De Belly 1D length and two-site competition screen."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "wgd_v1" / "debelly_regime_screen"
JULIA = Path("/app/eb/software/Julia/1.8.2-linux-x86_64/bin/julia")
PROJECT = ROOT / "src" / "debelly" / "one_d"
SOLVER = PROJECT / "protocols" / "mechanochemical_opto_front2back.jl"

SINGLE_LENGTHS = (12.0, 16.0, 20.0, 2 ** (1 / 3) * 20, 32.0, 36.0, 40.0, 50.0, 64.0)
TWO_LENGTHS = (20.0, 2 ** (1 / 3) * 20, 32.0, 40.0, 50.0)
ASYMMETRIES = (0.0, 0.02, 0.10)
FRICTIONS = (100.0, 150.0, 200.0, 250.0, 300.0, 400.0)


def tasks():
    jobs = [dict(mode="single", length=length, asymmetry=0.0, friction=100.0)
            for length in SINGLE_LENGTHS]
    jobs.extend(dict(mode="two_site", length=length, asymmetry=asymmetry,
                     friction=100.0)
                for length in TWO_LENGTHS for asymmetry in ASYMMETRIES)
    return jobs


def friction_tasks():
    return [dict(mode="two_site", length=length, asymmetry=0.02,
                 friction=friction)
            for length in TWO_LENGTHS[:2] for friction in FRICTIONS[1:]]


def diagnostic_tasks():
    jobs = []
    for length in TWO_LENGTHS[:2]:
        jobs.append(dict(mode="two_site", length=length, asymmetry=0.02,
                         friction=100.0, variant="absolute_width",
                         opto_width=1.6, feedback=True))
        jobs.append(dict(mode="two_site", length=length, asymmetry=0.02,
                         friction=100.0, variant="no_mechanics",
                         feedback=False))
    return jobs


def label(job):
    length = f"{job['length']:.3f}".rstrip("0").rstrip(".")
    asymmetry = f"{100*job['asymmetry']:.0f}pct"
    suffix = f"_{job['variant']}" if "variant" in job else ""
    return f"{job['mode']}_L{length}_asym{asymmetry}_chi{job['friction']:.0f}{suffix}"


def run_job(job):
    output = OUT / label(job)
    output.mkdir(parents=True, exist_ok=True)
    partition = max(30, round(job["length"] / 0.4))
    width = job.get(
        "opto_width", job["length"] * (0.2 if job["mode"] == "single" else 0.08)
    )
    feedback = job.get("feedback", True)
    environment = os.environ.copy()
    environment.update({
        "DEBELLY_WRITE_OUTPUTS": "0",
        "DEBELLY_QUIET": "1",
        "DEBELLY_TRAJECTORY_OUTPUT": str(output),
        "DEBELLY_LENGTH": str(job["length"]),
        "DEBELLY_PARTITION": str(partition),
        "DEBELLY_MCA_DENSITY": "0.1",
        "DEBELLY_MCA_NODAL_TOTAL": str(0.1 * partition),
        "DEBELLY_OPTO_WIDTH": str(width),
        "DEBELLY_OPTO_TIME": "50",
        "DEBELLY_DT": "1",
        "DEBELLY_T": "700",
        "DEBELLY_ALPHA_OPTO": "1",
        "DEBELLY_BETA_OPTO": "0",
        "DEBELLY_FRICTION": str(job["friction"]),
        "DEBELLY_STIMULUS_MODE": job["mode"],
        "DEBELLY_STIMULUS_ASYMMETRY": str(job["asymmetry"]),
        "DEBELLY_PROTRUSION_SPEED": "-0.15" if job["mode"] == "single" else "0",
        "DEBELLY_ACTIVE_STRESS": "100" if feedback else "0",
        "DEBELLY_RAC_MECHANICAL_FEEDBACK": "0.2" if feedback else "0",
        "DEBELLY_RHO_MECHANICAL_FEEDBACK": "0.2" if feedback else "0",
    })
    completed = subprocess.run(
        [str(JULIA), f"--project={PROJECT}", str(SOLVER)],
        cwd=PROJECT / "protocols", env=environment, check=True,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    metadata = {**job, "partition": partition, "dx": job["length"] / partition,
                "opto_width": width, "end_time": 700.0}
    (output / "task.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (output / "runner_stdout.txt").write_text(completed.stdout)
    (output / "runner_stderr.txt").write_text(completed.stderr)
    print(output)


def run_task(index):
    run_job(tasks()[index])


def load(path, filename, fields):
    raw = np.genfromtxt(path / filename, delimiter=",", names=True)
    times = np.unique(raw["time_s"])
    x = np.unique(raw["position_um"])
    arrays = {field: np.asarray(raw[field]).reshape(times.size, x.size)
              for field in fields}
    return times, x, arrays


def peak(values, times, start=50):
    indices = np.flatnonzero(times >= start)
    index = indices[np.argmax(np.abs(values[indices]))]
    return float(values[index]), float(times[index])


def sustained_selection(counts, times, start=150, duration=50):
    first = np.searchsorted(times, start)
    width = max(1, round(duration / np.median(np.diff(times))))
    for index in range(first, len(times) - width + 1):
        if np.all(counts[index:index + width] == 1):
            return float(times[index])
    return np.nan


def analyze_task(path):
    job = json.loads((path / "task.json").read_text())
    times, x, fields = load(
        path, "trajectory.csv", ("rac", "rho", "mca_bound", "mca_unbound")
    )
    _, tension_x, tension = load(path, "tension.csv", ("tension",))
    _, _, mechanics = load(path, "mechanics.csv", ("velocity", "displacement"))
    relative = x / x[-1]
    initial = float(np.mean(fields["rac"][0]))
    row = dict(job, task=path.name, initial_rac=initial)
    for field in ("rac", "rho"):
        proximal = np.mean(fields[field][:, relative <= 0.2], axis=1)
        distal = np.mean(fields[field][:, relative >= 0.8], axis=1)
        ppeak, ptime = peak(proximal - proximal[0], times)
        dpeak, dtime = peak(distal - distal[0], times)
        row.update({f"{field}_proximal_peak": ppeak,
                    f"{field}_distal_peak": dpeak,
                    f"{field}_distal_peak_time": dtime,
                    f"{field}_distal_to_proximal": abs(dpeak) / max(abs(ppeak), 1e-12)})
    for at in (350, 700):
        index = int(np.argmin(np.abs(times - at)))
        profile = fields["rac"][index]
        row[f"rac_contrast_{at}"] = float(
            (np.max(profile) - np.min(profile)) /
            max(np.max(profile) + np.min(profile), 1e-12)
        )
    row["maximum_displacement"] = float(np.max(np.abs(mechanics["displacement"])))
    total_mca = fields["mca_bound"] + fields["mca_unbound"]
    sums = np.sum(total_mca, axis=1)
    row["mca_conservation_error"] = float(np.max(np.abs(sums / sums[0] - 1)))

    tension_relative = tension_x / tension_x[-1]
    distal_tension = np.mean(tension["tension"][:, tension_relative >= 0.8], axis=1)
    row["tension_distal_peak"], row["tension_distal_peak_time"] = peak(
        distal_tension - distal_tension[0], times
    )

    if job["mode"] == "two_site":
        distance = max(1, round(0.25 * x.size))
        counts = np.zeros(times.size, dtype=int)
        for index, profile in enumerate(fields["rac"]):
            peaks, _ = find_peaks(profile, height=initial + 0.05,
                                  prominence=0.03, distance=distance)
            counts[index] = len(peaks)
        after = times >= 150
        row["maximum_domains"] = int(np.max(counts))
        row["domains_350"] = int(counts[np.argmin(np.abs(times - 350))])
        row["domains_700"] = int(counts[-1])
        row["multipolar_lifetime"] = float(np.sum((counts >= 2) & after))
        multi_times = times[(counts >= 2) & after]
        row["last_multipolar_time"] = float(multi_times[-1]) if multi_times.size else np.nan
        row["winner_selection_time"] = sustained_selection(counts, times)
        left = np.max(fields["rac"][:, (relative >= 0.1) & (relative <= 0.35)], axis=1)
        right = np.max(fields["rac"][:, (relative >= 0.65) & (relative <= 0.9)], axis=1)
        row["final_site_imbalance"] = float(
            abs(left[-1] - right[-1]) / max(abs(left[-1]) + abs(right[-1]), 1e-12)
        )
    return row


def analyze():
    paths = [OUT / label(job) for job in tasks()]
    missing = [path for path in paths if not (path / "trajectory.csv").is_file()]
    if missing:
        raise SystemExit(f"missing {len(missing)} tasks, first: {missing[0]}")
    rows = [analyze_task(path) for path in paths]
    keys = sorted({key for row in rows for key in row})
    with (OUT / "screen_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

    single = [row for row in rows if row["mode"] == "single"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    axes[0, 0].plot([r["length"] for r in single],
                    [abs(r["rac_distal_peak"]) for r in single], "o-")
    axes[0, 0].set_ylabel("|Distal Rac peak|")
    axes[0, 1].plot([r["length"] for r in single],
                    [r["rac_distal_peak_time"] for r in single], "o-")
    axes[0, 1].set_ylabel("Distal Rac peak time (s)")
    axes[1, 0].plot([r["length"] for r in single],
                    [r["rac_contrast_350"] for r in single], "o-", label="350 s")
    axes[1, 0].plot([r["length"] for r in single],
                    [r["rac_contrast_700"] for r in single], "s--", label="700 s")
    axes[1, 0].set_ylabel("Rac contrast")
    axes[1, 0].legend(frameon=False)
    axes[1, 1].plot([r["length"] for r in single],
                    [r["maximum_displacement"] for r in single], "o-")
    axes[1, 1].set_ylabel("Maximum displacement")
    for axis in axes[1]:
        axis.set_xlabel("Length (µm)")
    fig.suptitle("Single-cue length continuation")
    fig.tight_layout()
    fig.savefig(OUT / "single_length_screen.png", dpi=180)
    plt.close(fig)

    single_runs = []
    for row in single:
        path = OUT / row["task"]
        times, x, fields = load(path, "trajectory.csv", ("rac", "rho"))
        single_runs.append((row, times, x, fields))

    selected_lengths = (20.0, TWO_LENGTHS[1], 40.0, 64.0)
    selected_runs = [run for run in single_runs
                     if any(abs(run[0]["length"] - value) < 1e-6
                            for value in selected_lengths)]
    fig, axes = plt.subplots(len(selected_runs), 2, figsize=(11, 11),
                             sharex=True, sharey=True)
    for column, field in enumerate(("rac", "rho")):
        limit = max(np.max(np.abs(fields[field] - fields[field][0]))
                    for _, _, _, fields in selected_runs)
        for row_index, (row, times, x, fields) in enumerate(selected_runs):
            change = fields[field] - fields[field][0]
            image = axes[row_index, column].imshow(
                change, origin="lower", aspect="auto", cmap="RdBu_r",
                vmin=-limit, vmax=limit,
                extent=(0, 1, times[0], times[-1]),
            )
            axes[row_index, column].axhline(50, color="black", ls=":", lw=0.7)
            axes[row_index, column].axhline(150, color="black", ls=":", lw=0.7)
            axes[row_index, column].set_title(
                f"L={row['length']:.1f} µm: Δ{field.capitalize()}"
            )
            fig.colorbar(image, ax=axes[row_index, column], shrink=0.78)
    for axis in axes[:, 0]:
        axis.set_ylabel("Time (s)")
    for axis in axes[-1]:
        axis.set_xlabel("Relative position (stimulated → distal)")
    fig.suptitle("Where and when the single-cue response persists\n"
                 "Dotted lines mark stimulus onset and withdrawal")
    fig.tight_layout()
    fig.savefig(OUT / "single_length_response_maps.png", dpi=180)
    plt.close(fig)

    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, len(single_runs)))
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for color, (row, times, x, fields) in zip(colors, single_runs):
        relative = x / x[-1]
        for axis, field, region in (
            (axes[0], "rac", relative <= 0.2),
            (axes[1], "rac", relative >= 0.8),
            (axes[2], "rho", relative >= 0.8),
        ):
            mean = np.mean(fields[field][:, region], axis=1)
            axis.plot(times, mean - mean[0], color=color,
                      label=f"{row['length']:.1f} µm")
    axes[0].set_ylabel("Proximal ΔRac")
    axes[1].set_ylabel("Distal ΔRac")
    axes[2].set_ylabel("Distal ΔRho")
    for axis in axes:
        axis.axvspan(50, 150, color="#BBBBBB", alpha=0.22)
        axis.axhline(0, color="black", lw=0.6)
    axes[0].legend(ncol=3, fontsize=8, frameon=False)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("Region-mean responses; proximal/distal are the first/last 20%")
    fig.tight_layout()
    fig.savefig(OUT / "single_length_region_timecourses.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    profile_times = (150, 250, 350, 500, 700)
    profile_colors = plt.get_cmap("plasma")(np.linspace(0.05, 0.9, len(profile_times)))
    for axis, (row, times, x, fields) in zip(axes.flat, selected_runs):
        relative = x / x[-1]
        for time, color in zip(profile_times, profile_colors):
            index = int(np.argmin(np.abs(times - time)))
            axis.plot(relative, fields["rac"][index] - fields["rac"][0],
                      color=color, label=f"{time} s")
        axis.axhline(0, color="black", lw=0.6)
        axis.set_title(f"L={row['length']:.1f} µm")
        axis.set_ylabel("ΔRac from pre-stimulus state")
    for axis in axes[-1]:
        axis.set_xlabel("Relative position (stimulated → distal)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Rac profiles show delayed distal inhibition and recovery")
    fig.tight_layout()
    fig.savefig(OUT / "single_length_rac_profiles.png", dpi=180)
    plt.close(fig)

    fixed_times = (250, 350, 500, 700)
    fixed_colors = plt.get_cmap("magma")(np.linspace(0.15, 0.85, len(fixed_times)))
    lengths = []
    suppression = {time: [] for time in fixed_times}
    strongest_times = []
    half_recovery_times = []
    strongest_censored = []
    recovery_censored = []
    remaining_350 = []
    remaining_700 = []
    for row, times, x, fields in single_runs:
        lengths.append(row["length"])
        relative = x / x[-1]
        rac_change = fields["rac"] - fields["rac"][0]
        distal = np.mean(rac_change[:, relative >= 0.8], axis=1)
        for time in fixed_times:
            suppression[time].append(-distal[np.argmin(np.abs(times - time))])
        minimum_index = np.argmin(distal[times >= 50]) + np.flatnonzero(times >= 50)[0]
        strongest_times.append(times[minimum_index])
        strongest_censored.append(minimum_index == len(times) - 1)
        recovered = np.flatnonzero(
            (np.arange(len(times)) > minimum_index) &
            (np.abs(distal) <= 0.5 * abs(distal[minimum_index]))
        )
        half_recovery_times.append(times[recovered[0]] if recovered.size else times[-1])
        recovery_censored.append(not recovered.size)
        rms = np.sqrt(np.mean(rac_change**2, axis=1))
        peak_rms = np.max(rms[times >= 50])
        remaining_350.append(rms[np.argmin(np.abs(times - 350))] / peak_rms)
        remaining_700.append(rms[-1] / peak_rms)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for time, color in zip(fixed_times, fixed_colors):
        axes[0].plot(lengths, suppression[time], "o-", color=color,
                     label=f"{time} s")
    axes[0].set_ylabel("Distal Rac suppression (−ΔRac)")
    axes[0].legend(frameon=False)
    axes[1].plot(lengths, strongest_times, "o-", label="Strongest suppression")
    axes[1].plot(lengths, half_recovery_times, "s--", label="Half recovery")
    for x_value, y_value, censored in zip(lengths, strongest_times, strongest_censored):
        if censored:
            axes[1].scatter(x_value, y_value, marker="^", s=90, color="#0072B2")
    for x_value, y_value, censored in zip(lengths, half_recovery_times, recovery_censored):
        if censored:
            axes[1].scatter(x_value, y_value, marker="^", s=90, color="#D55E00")
    axes[1].set_ylabel("Time (s); triangles are right-censored")
    axes[1].legend(frameon=False)
    axes[2].plot(lengths, remaining_350, "o-", label="350 s")
    axes[2].plot(lengths, remaining_700, "s--", label="700 s")
    axes[2].set_ylabel("Whole-cell Rac response remaining\n(fraction of peak RMS response)")
    axes[2].legend(frameon=False)
    for axis in axes:
        axis.set_xlabel("Length (µm)")
    fig.suptitle("Length changes persistence and recovery, not simply 'polarity strength'")
    fig.tight_layout()
    fig.savefig(OUT / "single_length_recovery_summary.png", dpi=180)
    plt.close(fig)

    two = [row for row in rows if row["mode"] == "two_site"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for asymmetry in ASYMMETRIES:
        subset = [row for row in two if row["asymmetry"] == asymmetry]
        axes[0].plot([r["length"] for r in subset],
                     [r["multipolar_lifetime"] for r in subset], "o-",
                     label=f"{100*asymmetry:.0f}% asymmetry")
        axes[1].plot([r["length"] for r in subset],
                     [r["final_site_imbalance"] for r in subset], "o-")
    axes[0].set_ylabel("Multipolar lifetime after pulse (s)")
    axes[0].legend(frameon=False)
    axes[1].set_ylabel("Final two-site imbalance")
    for axis in axes:
        axis.set_xlabel("Length (µm)")
    fig.suptitle("Two-site Rac competition")
    fig.tight_layout()
    fig.savefig(OUT / "two_site_screen.png", dpi=180)
    plt.close(fig)
    print(OUT / "screen_summary.csv")


def analyze_friction():
    jobs = [dict(mode="two_site", length=length, asymmetry=0.02,
                 friction=friction)
            for length in TWO_LENGTHS[:2] for friction in FRICTIONS]
    paths = [OUT / label(job) for job in jobs]
    missing = [path for path in paths if not (path / "trajectory.csv").is_file()]
    if missing:
        raise SystemExit(f"missing {len(missing)} friction tasks, first: {missing[0]}")
    rows = [analyze_task(path) for path in paths]
    keys = sorted({key for row in rows for key in row})
    with (OUT / "friction_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for length in TWO_LENGTHS[:2]:
        subset = [row for row in rows if row["length"] == length]
        tag = "1x" if length == TWO_LENGTHS[0] else "WGD"
        axes[0].plot([r["friction"] for r in subset],
                     [r["multipolar_lifetime"] for r in subset], "o-", label=tag)
        axes[1].plot([r["friction"] for r in subset],
                     [r["domains_350"] for r in subset], "o-", label=tag)
    axes[0].set_ylabel("Multipolar lifetime after pulse (s)")
    axes[1].set_ylabel("Rac domains at 350 s")
    for axis in axes:
        axis.set_xlabel("Friction coefficient χ")
        axis.legend(frameon=False)
    fig.suptitle("Friction tunes the 1x/WGD mechanical crossover")
    fig.tight_layout()
    fig.savefig(OUT / "friction_crossover_screen.png", dpi=180)
    plt.close(fig)
    print(OUT / "friction_summary.csv")


def analyze_diagnostics():
    reference = [dict(mode="two_site", length=length, asymmetry=0.02,
                      friction=100.0, variant="relative_full")
                 for length in TWO_LENGTHS[:2]]
    paths = []
    labels = []
    for job in reference:
        base = {key: value for key, value in job.items() if key != "variant"}
        paths.append(OUT / label(base))
        labels.append(job["variant"])
    for job in diagnostic_tasks():
        paths.append(OUT / label(job))
        labels.append(job["variant"])
    rows = []
    for path, variant in zip(paths, labels):
        row = analyze_task(path)
        row["variant"] = variant
        rows.append(row)
    keys = sorted({key for row in rows for key in row})
    with (OUT / "diagnostic_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    print(OUT / "diagnostic_summary.csv")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=int)
    parser.add_argument("--friction-task", type=int)
    parser.add_argument("--diagnostic-task", type=int)
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--analyze-friction", action="store_true")
    parser.add_argument("--analyze-diagnostics", action="store_true")
    parser.add_argument("--count", action="store_true")
    parser.add_argument("--friction-count", action="store_true")
    parser.add_argument("--diagnostic-count", action="store_true")
    args = parser.parse_args()
    if args.count:
        print(len(tasks()))
    elif args.friction_count:
        print(len(friction_tasks()))
    elif args.diagnostic_count:
        print(len(diagnostic_tasks()))
    elif args.analyze:
        analyze()
    elif args.analyze_friction:
        analyze_friction()
    elif args.analyze_diagnostics:
        analyze_diagnostics()
    elif args.task is not None:
        run_task(args.task)
    elif args.friction_task is not None:
        run_job(friction_tasks()[args.friction_task])
    elif args.diagnostic_task is not None:
        run_job(diagnostic_tasks()[args.diagnostic_task])
    else:
        parser.error("choose --task, --analyze, or --count")


if __name__ == "__main__":
    main()
