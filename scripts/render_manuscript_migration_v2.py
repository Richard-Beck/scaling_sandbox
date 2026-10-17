"""Reassemble Figure 5 with truthful MSD points and matched-window slopes.

This is an outline-owned editorial renderer: it performs no fit or simulation
and never writes to the canonical 3 October campaign. Run with python -B.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter, NullLocator
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "wgd_langevin/refined_ou_20261003"
CANONICAL = CAMPAIGN / "figures/main_03_migration_OU_provenance.json"
sys.path.insert(0, str(CAMPAIGN / "scripts"))
import figure_panels as panels
import figure_tracks as tracks

NAMES = ["Xuan2018", "Mallin2023", "Jemaa2017", "Jemaa2023", "Hosny_baseline"]
SCENARIOS = ["primary", "prior5", "prior20", "strong", "weak"]
GROUPS = ["WGD−", "WGD+"]
COLORS = [tracks.MINUS, tracks.PLUS]
NAVY = tracks.NAVY
SURFACE = CAMPAIGN / "results/main_relative_surface.csv"
VISIBLE_CAVEAT = (
    "Points sample digitized published curves; original experimental lag grid unknown."
)
B_CAPTION = (
    "B, Upper: twelve chosen display samples from each digitized Xuan population "
    "MSD curve (small points) and the unchanged fitted analytical OU expectation "
    "(solid lines), on logarithmic lag and MSD axes. These are samples of the "
    "published raster curves, not recovered experimental lag positions or "
    "independent cell observations; the original experimental lag grid is unknown. "
    "The open circle identifies the crossover of the full digitized traces at "
    "2.32014 h, computed by linear interpolation in native lag and MSD. The open "
    "diamond identifies the different fitted OU crossover at 1.11406 h. Lower: "
    "descriptive slopes of log(MSD) versus log(lag) fitted in six equal "
    "log-lag windows of the full digitized curves (points), and matching-window "
    "moving OU slopes (solid lines). Each window spans approximately 0.304 "
    "decades; both observed and model slopes use the same full window support. "
    "Observed raster columns remain correlated and give no biological sample "
    "size or confidence interval. These slopes are newly calculated curve-shape "
    "summaries, not author-reported local exponents or fitted OU memory times. "
    "A constant power-law exponent would appear horizontal; reference values "
    "1 and 2 indicate diffusive and straight-line scaling. The model does not "
    "reproduce the complete observed curvature or the observed crossover. "
    "The analytical fit, all fitting targets and numerical inputs to A and "
    "C–F are unchanged."
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def log_slope(log_t: np.ndarray, log_y: np.ndarray) -> float:
    x = np.asarray(log_t) - np.mean(log_t)
    return float(np.dot(x, np.asarray(log_y) - np.mean(log_y)) / np.dot(x, x))


def letter(ax, label: str) -> None:
    ax.text(-.055, 1.04, label, transform=ax.transAxes, color=NAVY,
            fontsize=14, fontweight="bold", va="bottom")


def prepare_data(original: dict, out: Path) -> dict:
    """Read exact retained sources, validate analytical predictions, export display data."""
    scenarios = {}
    dynamic = panels.load_fit_scenarios()
    for name in SCENARIOS:
        source = original["scenario_fit_files"][name]
        path = ROOT / source["path"]
        assert digest(path) == source["sha256"], (name, "changed scenario input")
        assert dynamic[name]["path"] == path, (name, "canonical scenario selection changed")
        record = json.loads(path.read_text())
        scenarios[name] = dict(path=path, record=record, z=record["z"])
    primary = scenarios["primary"]
    assert primary["path"].name == original["selected_primary_fit"]
    pars = np.exp(np.asarray(primary["z"][:20], dtype=float).reshape(5, 2, 2))
    observation_rows, _ = panels.build_observation_rows(primary)
    assert len(observation_rows) == 9
    # Read the original table independently before accepting the reconstructed values.
    with (CAMPAIGN / "figures/main_03_observation_folds.csv").open() as stream:
        saved_rows = list(csv.DictReader(stream))
    for current, saved in zip(observation_rows, saved_rows):
        assert (current["study"], current["metric"]) == (saved["study"], saved["metric"])
        for key in ("observed_WGDminus", "observed_WGDplus", "predicted_WGDminus",
                    "predicted_WGDplus", "observed_fold_WGDplus_over_WGDminus",
                    "predicted_fold_WGDplus_over_WGDminus"):
            assert current[key] == float(saved[key]), (current["study"], key)
    write_csv(out / "panel_a_observed_predicted.csv", observation_rows)

    curve_rows = primary["record"]["rows"]
    curves = []
    prediction_checks = []
    for group in (0, 1):
        rows = sorted((r for r in curve_rows if r["study"] == "Xuan2018"
                       and r["metric"] == "msd_curve" and int(r["group"]) == group),
                      key=lambda row: row["lag_h"])
        t = np.array([r["lag_h"] for r in rows], dtype=float)
        y = np.array([r["obs"] for r in rows], dtype=float)
        expected = tracks.msd_expected(t, *pars[0, group])
        retained = np.array([r["predicted"] for r in rows], dtype=float)
        np.testing.assert_allclose(expected, retained, rtol=1e-12, atol=1e-10)
        prediction_checks.append(dict(group=group, native_fitting_lags=len(rows),
                                      max_absolute_difference=float(np.max(np.abs(expected-retained))),
                                      passed=True))
        curves.append((t, y))
    lower = max(t[0] for t, _ in curves)
    upper = min(t[-1] for t, _ in curves)
    edges = np.linspace(np.log(lower), np.log(upper), 7)
    centers = (edges[:-1] + edges[1:]) / 2
    halfwidth = (edges[1] - edges[0]) / 2
    model_centers = np.linspace(centers[0], centers[-1], 300)

    def model_window(center, group):
        x = center + np.linspace(-halfwidth, halfwidth, 151)
        return log_slope(x, np.log(tracks.msd_expected(np.exp(x), *pars[0, group])))

    synthetic = []
    for exponent in (.75, 1., 1.5, 2., 2.4):
        for x in (np.linspace(np.log(lower), np.log(upper), 151),
                  np.array([-.9, -.5, -.15, .21, .8, 1.3])):
            recovered = log_slope(x, np.log(7.) + exponent*x)
            assert abs(recovered-exponent) < 1e-12
            synthetic.append(dict(expected=exponent, recovered=recovered, passed=True))
    display, slopes = [], []
    for group, (t, y) in enumerate(curves):
        indices = np.unique(np.rint(np.linspace(0, len(t)-1, 12)).astype(int))
        assert len(indices) == 12
        for index in indices:
            display.append(dict(group=group, group_label=GROUPS[group],
                                native_raster_index=int(index), lag_h=float(t[index]),
                                population_msd_um2=float(y[index]),
                                evidence="chosen display sample of digitized published line"))
        x = np.log(t)
        for i, center in enumerate(centers):
            mask = (x >= edges[i]) & (x < edges[i+1] if i < 5 else x <= edges[i+1])
            assert mask.sum() >= 20
            slopes.append(dict(group=group, group_label=GROUPS[group], window=i+1,
                               lag_low_h=float(np.exp(edges[i])),
                               lag_center_h=float(np.exp(center)),
                               lag_high_h=float(np.exp(edges[i+1])),
                               correlated_raster_samples=int(mask.sum()),
                               observed_window_slope=log_slope(x[mask], np.log(y[mask])),
                               model_window_slope=model_window(center, group)))
    model_slopes = np.array([[model_window(c, group) for c in model_centers] for group in (0, 1)])
    assert np.all((model_slopes >= 1) & (model_slopes <= 2))
    assert np.all(np.diff(model_slopes, axis=1) < 0)
    native = np.unique(np.concatenate([t[(t >= lower) & (t <= upper)] for t, _ in curves]))
    difference = lambda t: np.interp(t, *curves[0]) - np.interp(t, *curves[1])
    differences = difference(native)
    crossings = [float(brentq(difference, native[i], native[i+1]))
                 for i in np.flatnonzero(differences[:-1]*differences[1:] < 0)]
    grid = np.geomspace(lower, upper, 600)
    model_difference = lambda t: float(tracks.msd_expected(t, *pars[0, 0]) -
                                       tracks.msd_expected(t, *pars[0, 1]))
    values = np.array([model_difference(t) for t in grid])
    fitted_crossings = [float(brentq(model_difference, grid[i], grid[i+1]))
                       for i in np.flatnonzero(values[:-1]*values[1:] < 0)]
    np.testing.assert_allclose(crossings, original["population_curve_crossovers_h"], rtol=0, atol=1e-11)
    np.testing.assert_allclose(fitted_crossings, original["fitted_ou_crossovers_h"], rtol=0, atol=1e-11)
    write_csv(out / "panel_b_raster_display_samples.csv", display)
    write_csv(out / "panel_b_window_slopes.csv", slopes)
    write_csv(out / "panel_b_model_window_slopes.csv", [
        dict(lag_center_h=float(np.exp(c)), reference_slope=float(model_slopes[0, i]),
             enlarged_slope=float(model_slopes[1, i])) for i, c in enumerate(model_centers)])

    chosen, metrics, selections = [], [], {}
    with np.load(CAMPAIGN / "results/tracks_primary.npz") as archive:
        for group in (0, 1):
            path, metric, selection = tracks.choose_representative(archive[f"Xuan2018_{group}"])
            assert selection["index"] == original["track_selection_indices"][GROUPS[group]]
            chosen.append(path.copy())
            metrics.append(metric)
            selections[GROUPS[group]] = selection["index"]
    write_csv(out / "panel_c_track_metrics.csv", [
        dict(group=group, group_label=GROUPS[group], selected_track_index=selections[GROUPS[group]],
             **metric) for group, metric in enumerate(metrics)])
    write_csv(out / "panel_c_selected_tracks.csv", [
        dict(group=group, sample=i, time_h=i*tracks.CAD, x_um=float(x), y_um=float(y))
        for group, path in enumerate(chosen) for i, (x, y) in enumerate(path)])
    parameter_rows = [dict(study=name, group=group, group_label=GROUPS[group],
                           tau_h=float(pars[i, group, 0]), tau_min=float(pars[i, group, 0]*60),
                           sigma_um_h_1p5=float(pars[i, group, 1]))
                      for i, name in enumerate(NAMES) for group in (0, 1)]
    write_csv(out / "panel_e_parameters.csv", parameter_rows)
    scenario_folds = []
    for name in SCENARIOS:
        sp = np.exp(np.asarray(scenarios[name]["z"][:20], dtype=float).reshape(5, 2, 2))
        scenario_folds.append(sp[:, 1, 0] / sp[:, 0, 0])
    folds = np.array(scenario_folds).T
    display_folds = np.vstack([folds, np.exp(np.mean(np.log(folds), axis=0))])
    fold_rows = [dict(study=name, **{s: float(display_folds[i, j]) for j, s in enumerate(SCENARIOS)})
                 for i, name in enumerate(NAMES + ["Pooled geometric mean"])]
    with (CAMPAIGN / "figures/figure3_persistence.csv").open() as stream:
        saved_folds = list(csv.DictReader(stream))
    for current, saved in zip(fold_rows, saved_folds):
        for name in SCENARIOS:
            assert current[name] == float(saved[name])
    for i, row in enumerate(saved_folds[:5]):
        assert float(row["reference_minutes"]) == pars[i, 0, 0]*60
        assert float(row["polyploid_minutes"]) == pars[i, 1, 0]*60
    write_csv(out / "panel_f_scenario_folds.csv", fold_rows)
    return dict(pars=pars, observations=observation_rows, curves=curves, lower=lower,
                upper=upper, display=display, slopes=slopes, centers=model_centers,
                model_slopes=model_slopes, crossings=crossings, fitted_crossings=fitted_crossings,
                chosen=chosen, metrics=metrics, selections=selections, display_folds=display_folds,
                validation=dict(original_fitting_lag_predictions=prediction_checks,
                                synthetic_power_laws=synthetic,
                                panel_a_exact_table_match=True, panel_c_exact_indices_match=True,
                                panel_e_exact_parameter_table_match=True, panel_f_exact_folds_match=True,
                                observed_and_model_crossovers_unchanged=True,
                                model_window_slopes_between_one_and_two=True))


def draw_b(fig, slot, data: dict, size=8, label=True):
    sub = slot.subgridspec(2, 1, height_ratios=[2.15, 1.25], hspace=.17)
    upper = fig.add_subplot(sub[0])
    lower = fig.add_subplot(sub[1], sharex=upper)
    if label:
        letter(upper, "B")
    grid = np.geomspace(data["lower"], data["upper"], 600)
    for group, color in enumerate(COLORS):
        selected = [r for r in data["display"] if r["group"] == group]
        upper.plot([r["lag_h"] for r in selected], [r["population_msd_um2"] for r in selected],
                   ls="none", marker="o", ms=3, color=color, zorder=4)
        upper.plot(grid, tracks.msd_expected(grid, *data["pars"][0, group]), color=color, lw=1.7)
        group_slopes = [r for r in data["slopes"] if r["group"] == group]
        lower.plot([r["lag_center_h"] for r in group_slopes],
                   [r["observed_window_slope"] for r in group_slopes],
                   ls="none", marker="o", ms=4.3, color=color, mec="white", mew=.35, zorder=4)
        lower.plot(np.exp(data["centers"]), data["model_slopes"][group], color=color, lw=1.7)
    observed_cross = data["crossings"][0]
    model_cross = data["fitted_crossings"][0]
    upper.scatter(observed_cross, np.interp(observed_cross, *data["curves"][0]),
                  s=35, marker="o", facecolor="white", edgecolor=NAVY, lw=1., zorder=5)
    upper.scatter(model_cross, tracks.msd_expected(model_cross, *data["pars"][0, 0]),
                  s=30, marker="D", facecolor="white", edgecolor="#526876", lw=1., zorder=5)
    upper.text(.98, .14, f"Curve crossover {observed_cross:.2f} h ○",
               transform=upper.transAxes, ha="right", fontsize=size, color=NAVY)
    upper.text(.98, .035, f"OU crossover {model_cross:.2f} h ◇",
               transform=upper.transAxes, ha="right", fontsize=size, color="#526876")
    upper.set(xscale="log", yscale="log", ylabel="MSD (µm²)",
              xlim=(data["lower"]*.94, data["upper"]*1.07))
    upper.tick_params(axis="x", which="both", labelbottom=False)
    upper.set_title("Xuan 2018 · population curve shape", loc="left", fontsize=size+1, color=NAVY)
    upper.legend(handles=[
        Line2D([], [], color=COLORS[0], marker="o", ls="none", ms=3, label="WGD− curve samples"),
        Line2D([], [], color=COLORS[1], marker="o", ls="none", ms=3, label="WGD+ curve samples"),
        Line2D([], [], color="#526876", lw=1.7, label="OU model lines")],
        loc="upper left", frameon=False, fontsize=size-1, ncol=1)
    lower.set(ylabel="Window slope α", xlabel="Lag time (h; log scale)",
              ylim=(.7, 2.1), yticks=[1, 1.5, 2])
    for level in (1, 2):
        lower.axhline(level, ls=":", color="#cad4da", lw=.65, zorder=0)
    lower.xaxis.set_major_locator(FixedLocator([.2, .5, 1, 2, 5, 10]))
    lower.xaxis.set_major_formatter(FixedFormatter(["0.2", "0.5", "1", "2", "5", "10"]))
    lower.xaxis.set_minor_formatter(NullFormatter())
    for ax in (upper, lower):
        ax.tick_params(labelsize=size)
        ax.xaxis.label.set_size(size)
        ax.yaxis.label.set_size(size)
        ax.grid(axis="y", which="major", color="#e9edf0", lw=.6)
    return upper, lower


def integrated(data: dict, out: Path):
    fig = plt.figure(figsize=(16, 13.4))
    gs = fig.add_gridspec(3, 2, left=.075, right=.965, bottom=.115, top=.92,
                          height_ratios=[1.35, 1.15, .78], width_ratios=[1.08, 1.],
                          hspace=.42, wspace=.24)
    ax_a = fig.add_subplot(gs[0, 0])
    panels.draw_observed_predicted(ax_a, data["observations"])
    letter(ax_a, "A")
    ax_a.get_legend().remove()
    ax_a.legend(handles=[
        Line2D([], [], marker="o", color="none", mfc="#5b7380", mec="white", label="Observed fold", ms=6),
        Line2D([], [], marker="D", color="none", mfc="#d36c32", mec="white", label="Predicted fold", ms=5.5),
        Line2D([], [], color="#33454e", ls="--", label="Equal groups")],
        frameon=False, ncol=3, fontsize=7, loc="lower left", bbox_to_anchor=(0, 1.02))
    ax_a.xaxis.set_major_locator(FixedLocator([.5, 1, 2, 3]))
    ax_a.xaxis.set_major_formatter(FixedFormatter(["0.5", "1", "2", "3"]))
    ax_a.xaxis.set_minor_locator(NullLocator())
    ax_a.xaxis.set_minor_formatter(NullFormatter())
    ax_a.tick_params(labelsize=8)
    draw_b(fig, gs[0, 1], data, size=8)

    sub = gs[1, 0].subgridspec(1, 2, wspace=.24)
    lim = float(max(np.abs(t).max() for t in data["chosen"])*1.13)
    for group, (name, color) in enumerate(zip(GROUPS, COLORS)):
        ax = fig.add_subplot(sub[0, group])
        if group == 0:
            letter(ax, "C")
        path = data["chosen"][group]
        metric = data["metrics"][group]
        tracks.axes_track(ax, path, color, linewidth=1.5)
        tracks.style_axis(ax, lim)
        tracks.arrow(ax, path[0], path[-1], NAVY, 1.)
        for j in (18, 36, 54):
            tracks.arrow(ax, path[j], path[j+1], color, .65)
        for a, b in ((0, 12), (60, 72)):
            ax.plot(path[[a, b], 0], path[[a, b], 1], ls="--", color=tracks.GRAY, lw=.85)
        ax.set_title(f"{name} · {metric['mean_sampled_speed_um_h']:.1f} µm/h\n"
                     f"{metric['path_length_um']:.0f} µm path · directness {metric['directness']:.2f}",
                     loc="left", color=color, fontweight="bold", fontsize=9, pad=6)
        ax.set_xlabel("x (µm)", fontsize=8)
        ax.set_ylabel("y (µm)", fontsize=8)
        ax.tick_params(labelsize=7)

    surface_grid = gs[1, 1].subgridspec(1, 3, width_ratios=[1, 1, .055], wspace=.38)
    a1, a2 = fig.add_subplot(surface_grid[0, 0]), fig.add_subplot(surface_grid[0, 1])
    cax = fig.add_subplot(surface_grid[0, 2])
    with SURFACE.open() as stream:
        source_rows = list(csv.DictReader(stream))
    scores = [float(r[k]) for r in source_rows for k in ("fixed_drift_delta_data", "fixed_drift_delta_total")]
    norm = matplotlib.colors.AsinhNorm(linear_width=1., vmin=min(-.01, min(scores)), vmax=max(scores))
    mesh = panels.draw_relative_surface(a1, csv_path=SURFACE, score_column="fixed_drift_delta_data",
                                       norm=norm, title="Data misfit", cmap="magma")
    panels.draw_relative_surface(a2, csv_path=SURFACE, score_column="fixed_drift_delta_total",
                                 norm=norm, title="Total objective", cmap="magma")
    a2.set_ylabel("")
    for ax in (a1, a2):
        ax.set_xlabel("WGD+ / WGD− persistence fold", fontsize=8)
        ax.set_ylabel("WGD+ / WGD− noise fold", fontsize=8)
        ax.tick_params(labelsize=7)
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label("Change in score from fitted point", fontsize=8)
    cb.set_ticks([0, 1, 10, 100, 1000, 10000])
    cb.ax.tick_params(labelsize=6)
    letter(a1, "D")
    a1.text(0, 1.13, "WGD− fixed; WGD+ persistence and noise shifted", transform=a1.transAxes,
            ha="left", va="bottom", color=NAVY, fontsize=9, fontweight="bold")
    a2.legend(handles=[
        Line2D([], [], marker="o", color="none", mfc="white", mec="#242f35", label="Equal group means", ms=5),
        Line2D([], [], marker="*", color="none", mfc="#111111", mec="white", label="Selected fit", ms=8)],
        frameon=False, fontsize=6.5, loc="upper right")

    sub_e = gs[2, 0].subgridspec(1, 2, wspace=.35)
    e1, e2 = fig.add_subplot(sub_e[0, 0]), fig.add_subplot(sub_e[0, 1])
    for k, (ax, title, unit) in enumerate(((e1, "Persistence", "τ (min)"),
                                         (e2, "Velocity noise", "σ (µm h⁻³ᐟ²)"))):
        offsets = [-.12, .12]
        for j, name in enumerate(NAMES):
            values = data["pars"][j, :, k]*(60 if k == 0 else 1)
            ax.plot(values, [j+v for v in offsets], color="#b8c2c9", lw=1, zorder=1)
            for group, (color, off) in enumerate(zip(COLORS, offsets)):
                ax.scatter(values[group], j+off, s=25, color=color, ec="white", lw=.45, zorder=3)
        ax.set_yticks(np.arange(5), NAMES)
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlabel(unit, fontsize=8)
        ax.set_title(title, loc="left", color=NAVY, fontsize=9)
        ax.grid(axis="x", alpha=.15)
        ax.tick_params(labelsize=7)
    e1.legend(handles=[Line2D([], [], marker="o", color="none", mfc=c, label=l, ms=5)
                       for l, c in zip(GROUPS, COLORS)], frameon=False, fontsize=7, loc="lower right")
    letter(e1, "E")

    ax_f = fig.add_subplot(gs[2, 1])
    folds = data["display_folds"]
    logfold = np.log2(folds)
    im = ax_f.imshow(logfold, aspect="auto", cmap="YlOrRd", vmin=0, vmax=max(1., float(np.max(logfold))))
    ax_f.set_xticks(range(5), ["Primary", "5-min prior", "20-min prior", "Strong pooling", "Weak pooling"],
                    rotation=22, ha="right")
    ax_f.set_yticks(range(6), NAMES + ["Pooled geometric mean"])
    ax_f.set_title("WGD+ / WGD− persistence fold", loc="left", color=NAVY, fontsize=9)
    for i in range(6):
        for j in range(5):
            ax_f.text(j, i, f"{folds[i, j]:.1f}", ha="center", va="center", fontsize=6.5,
                      color="white" if logfold[i, j] > 4.6 else "#40230f")
    cb = fig.colorbar(im, ax=ax_f, fraction=.046, pad=.025)
    cb.set_label("log₂ persistence fold", fontsize=7)
    cb.ax.tick_params(labelsize=7)
    letter(ax_f, "F")
    fig.suptitle("Figure 5 | Migration summaries and conditional OU dynamics", fontsize=15,
                 color=NAVY, fontweight="bold", y=.98)
    fig.text(.075, .052, VISIBLE_CAVEAT, fontsize=9, color=NAVY)
    fig.text(.075, .032, "B slopes describe curve shape, not independent measurements or OU memory times. "
             "C tracks are simulations. D/F show conditional objective/prior sensitivity.", fontsize=8, color=tracks.GRAY)
    fig.text(.075, .015, "No new fit or simulation; A and C–F retain their original numerical evidence. "
             "See the enlarged B detail for readable curve-shape comparison.", fontsize=8, color=tracks.GRAY)
    save(fig, out, "figure5_migration_v2")


def detail(data: dict, out: Path):
    fig = plt.figure(figsize=(7.6, 6.2))
    slot = fig.add_gridspec(1, 1, left=.14, right=.96, bottom=.19, top=.88)[0]
    draw_b(fig, slot, data, size=10, label=False)
    fig.suptitle("Figure 5B | Digitized population MSD and curve shape", x=.14, ha="left",
                 fontsize=12, color=NAVY, y=.97)
    fig.text(.14, .095, "Points sample digitized published curves;", fontsize=9, color=NAVY)
    fig.text(.14, .069, "original experimental lag grid unknown.", fontsize=9, color=NAVY)
    fig.text(.14, .035, "Lower points: newly estimated window slopes. Correlated raster samples; no biological error bars.",
             fontsize=7.4, color=tracks.GRAY)
    save(fig, out, "figure5b_msd_detail")


def save(fig, out: Path, stem: str):
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(out / f"{stem}.{suffix}", dpi=230, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "manuscript_outline/assets/migration_v2")
    args = parser.parse_args()
    out = args.out.resolve()
    assert not out.is_relative_to(CAMPAIGN), "Outputs must remain outside canonical campaign"
    original = json.loads(CANONICAL.read_text())
    protected = {ROOT / p: h for p, h in original["source_sha256"].items()}
    for item in original["scenario_fit_files"].values():
        protected[ROOT / item["path"]] = item["sha256"]
    for name in ("main_03_migration_OU.png", "main_03_migration_OU.pdf", "main_03_migration_OU.svg",
                 "main_03_migration_OU_caption.md", "main_03_migration_OU_provenance.json",
                 "figure3_persistence.csv", "figure3_primary_noise.csv"):
        path = CAMPAIGN / "figures" / name
        protected.setdefault(path, digest(path))
    for path, expected in protected.items():
        assert digest(path) == expected, (str(path), "original input changed before rendering")
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.titlesize": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    data = prepare_data(original, out)
    integrated(data, out)
    detail(data, out)
    for path, expected in protected.items():
        assert digest(path) == expected, (str(path), "canonical file changed during rendering")
    validation = dict(passed=True, fit_or_simulation_performed=False,
                      recorded_original_input_count=len(original["source_sha256"]),
                      protected_file_count=len(protected), protected_hashes_unchanged=True,
                      **data["validation"])
    write_json(out / "validation.json", validation)
    integration = dict(
        figure_id="figure5", title="Figure 5. Published migration summaries and conditional OU interpretation",
        assets=[dict(path=str((out / "figure5_migration_v2.png").relative_to(ROOT)),
                     caption="Reassembled Figure 5: original A and C–F numerical evidence, revised MSD/curve-shape panel B."),
                dict(path=str((out / "figure5b_msd_detail.png").relative_to(ROOT)),
                     caption="Figure 5B continued: enlarged detail with explicit raster-sample caveat.")],
        replace_legend_index=1, replacement_legend=B_CAPTION,
        additional_narrative=[
            "The observed Xuan curves cross near 2.32 h while the unchanged fitted OU curves cross near 1.11 h. "
            "Matched-window slopes make the remaining curve-shape discrepancy visible; fitted persistence is a conditional summary.",
            "A slowing cell can produce greater late displacement, but the enlarged-group migration phenotype "
            "and its causal connection to polarity mechanisms require direct experimental measurement."],
        resolved_limitation="The 13 October point-versus-model and exponent-display request is implemented, "
                            "with explicit raster-sample provenance; the unknown original experimental lag grid remains unmeasured.",
        preserved_other_legend_entries=True,
        user_feedback_evidence=[
            dict(session="/home/4473331/.codex/sessions/2026/10/13/rollout-2026-10-13T13-27-58-01a1364e-870d-7dd3-b0b8-0021c62eb57a.jsonl",
                 line=148, timestamp="2026-10-13T17:52:30.419Z", role="user",
                 excerpt="reserve lines only for the model",
                 meaning="Use observed curve-sample points and expose exponent/curvature behaviour."),
            dict(session="/home/4473331/.codex/sessions/2026/10/13/rollout-2026-10-13T13-27-58-01a1364e-870d-7dd3-b0b8-0021c62eb57a.jsonl",
                 line=214, timestamp="2026-10-13T18:09:00.378Z", role="user",
                 excerpt="Did Xuan really calculate values of MSD at as many lags",
                 meaning="Do not equate raster display samples with original experimental lag measurements.")])
    write_json(out / "outline_integration.json", integration)
    write_json(out / "caption.json", dict(figure_number=5, panel_b=B_CAPTION,
                                          visible_caveat=VISIBLE_CAVEAT))
    provenance = dict(
        generated_at=datetime.now(timezone.utc).isoformat(), producer="finished",
        status="editorial_reassembly_no_refit", renderer=str(Path(__file__).relative_to(ROOT)),
        renderer_sha256=digest(Path(__file__)), canonical_source=str(CANONICAL.relative_to(ROOT)),
        selected_primary_fit=original["selected_primary_fit"],
        protected_source_sha256={str(p.relative_to(ROOT)): h for p, h in protected.items()},
        original_recorded_input_hashes_all_verified=True, canonical_files_unmodified=True,
        population_curve_crossovers_h=data["crossings"], fitted_ou_crossovers_h=data["fitted_crossings"],
        observed_display_count_per_group=12, original_experimental_lag_grid="unknown",
        display_rule="Twelve evenly spaced native-raster indices per curve including endpoints; layout samples, not author lag locations.",
        slope_windows=6, slope_window_log10_width=float(np.log10(data["upper"]/data["lower"])/6),
        slope_method="OLS log(MSD) versus log(lag); observed native raster columns and matching-width moving model windows.",
        model_slope_window_points=151, slope_uncertainty="No biological error bars or sample size inferred from correlated raster samples.",
        preserved_panel_numerics=["A observed/predicted folds", "C selected native simulated tracks and metrics",
                                 "D both saved surfaces and shared asinh normalization",
                                 "E primary parameters", "F scenario persistence folds"],
        validation=validation,
        output_sha256={p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "provenance.json"})
    write_json(out / "provenance.json", provenance)
    print(json.dumps(dict(out=str(out), status=provenance["status"], validation_passed=True,
                          original_inputs_verified=len(original["source_sha256"]),
                          observed_crossover_h=data["crossings"][0],
                          model_crossover_h=data["fitted_crossings"][0])))


if __name__ == "__main__":
    main()
