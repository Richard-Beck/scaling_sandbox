"""Render the complete, audited ring functional comparison; never integrate.

The primary figures retain the full prospective size/placement matrix. Selected
trajectories are prespecified supporting anchors, not outcome-selected examples.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import os
import re
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm, SymLogNorm, ListedColormap
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
import length_dev_output_layout as layout

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/runs/chemical_ring_functional_v2_20261017"
TOPIC = ROOT / "figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017"
MODELS = ("wave_pinning", "goryachev", "holmes_model3")
NAMES = {"wave_pinning": "Wave-pinning", "goryachev": "Goryachev", "holmes_model3": "Holmes 3"}
PEAKS = {"wave_pinning": .1, "goryachev": 6.4, "holmes_model3": 1.}
PRIMARY = .005
SITE_LABELS = ("Single peak near old, unresolved", "Cue-associated peak", "Distinct old + cue peaks",
               "One region spans old + cue", "Shifted / merged, unresolved", "No resolved peak", "Ambiguous multiple peaks")
SITE_COLORS = ("#4678a3", "#d77c35", "#8b65a9", "#69a7a0", "#bcb069", "#e8e8e8", "#6e777c")
ARTIFACTS = []
SOURCES = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def source(path):
    path = Path(path).resolve()
    SOURCES[str(path)] = sha(path)
    return path


def rows(path):
    with source(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def f(row, key, default=np.nan):
    try:
        value = float(row.get(key, ""))
        return value if np.isfinite(value) else default
    except (TypeError, ValueError):
        return default


def truth(row, key):
    return str(row.get(key, "")).lower() in ("true", "1")


def pick(data, model=None, stage=None, epsilon=None, circumference=None, fraction=None):
    return [r for r in data if (model is None or r["model"] == model)
            and (stage is None or np.isclose(f(r, "time_from_onset"), stage))
            and (epsilon is None or np.isclose(f(r, "epsilon"), epsilon))
            and (circumference is None or np.isclose(f(r, "circumference"), circumference))
            and (fraction is None or np.isclose(f(r, "fraction"), fraction))]


def site_index(text):
    if text in ("one peak associated with old core", "one old-core-associated peak; incumbent identity unresolved"): return 0
    if text == "one peak associated with cue patch": return 1
    if text == "distinct peaks associated with old core and cue patch": return 2
    if text.startswith("one connected region spanning"): return 3
    if text.startswith("one shifted"): return 4
    if text.startswith("no peak"): return 5
    return 6


def old_associated(text):
    # A sole old-core-associated peak can be translated or merged. Only the
    # explicit distinct old + cue label supports this restricted association.
    return text == "distinct peaks associated with old core and cue patch"


def save(fig, stem, caption):
    for suffix in (".png", ".pdf", ".svg"):
        path = layout.dev_output_path(TOPIC, stem + suffix)
        fig.savefig(path, dpi=190 if suffix == ".png" else None, bbox_inches="tight")
        ARTIFACTS.append(path)
    plt.close(fig)
    path = layout.dev_output_path(TOPIC, stem + "_caption.txt")
    path.write_text(caption + "\n")
    ARTIFACTS.append(path)


def norm_values(values, signed=False, categorical=False):
    values = np.array([x for x in values if np.isfinite(x)])
    if categorical: return Normalize(-.5, 6.5), ListedColormap(SITE_COLORS)
    if not len(values): return Normalize(0, 1), plt.get_cmap("viridis")
    low, high = float(values.min()), float(values.max())
    if signed and low < 0:
        bound = max(abs(low), abs(high))
        # Linear around zero; logarithmic tails retain both decline and growth.
        return SymLogNorm(max(bound * 1e-3, 1e-12), vmin=-bound, vmax=bound), plt.get_cmap("RdBu")
    if np.isclose(low, high, rtol=1e-10, atol=1e-14):
        pad = max(abs(low) * .01, 1e-6)
        low, high = low - pad, high + pad
    return Normalize(low, high), plt.get_cmap("viridis")


def matrix(ax, base, model, key, title, sizes, values=None, signed=False,
           categorical=False, overlay=False, colorbar=True):
    """Row-specific rectangles use realized fraction, preserving center rounding."""
    subset = pick(base, model=model)
    lookup = values if values is not None else {(f(r, "circumference"), f(r, "fraction")): r for r in subset}
    observed = []
    for row in subset:
        record = lookup.get((f(row, "circumference"), f(row, "fraction")))
        if record is not None and row["status"] == "complete":
            val = site_index(record.get(key, "")) if categorical else f(record, key)
            if key == "net_repair_per_damage" and not truth(record, "normalized_gain_available"):
                val = np.nan
            observed.append(val)
    norm, cmap = norm_values(observed, signed=signed, categorical=categorical)
    for iy, circumference in enumerate(sizes):
        rr = sorted(pick(subset, circumference=circumference), key=lambda r: f(r, "fraction"))
        centers = np.array([f(r, "realized_fraction", f(r, "fraction")) for r in rr])
        edges = np.r_[centers[0] - .0125, (centers[:-1] + centers[1:]) / 2, centers[-1] + .0125]
        for j, row in enumerate(rr):
            record = lookup.get((circumference, f(row, "fraction")))
            available = row["status"] == "complete" and record is not None
            val = site_index(record.get(key, "")) if available and categorical else f(record, key) if available else np.nan
            if available and key == "net_repair_per_damage" and not truth(record, "normalized_gain_available"):
                val = np.nan
            face = cmap(norm(val)) if np.isfinite(val) else "white" if available else "#d5d5d5"
            corners = [(edges[j], iy - .5), (edges[j + 1], iy - .5),
                       (edges[j + 1], iy + .5), (edges[j], iy + .5)]
            ax.add_patch(Polygon(corners, facecolor=face, edgecolor="white", linewidth=.25))
            if not available:
                ax.plot(centers[j], iy, "x", color="#737373", ms=2.1, mew=.45)
            elif not np.isfinite(val):
                ax.plot(centers[j], iy, "+", color="#737373", ms=2.6, mew=.55)
            elif overlay and not old_associated(record.get("location_association", "")):
                ax.plot(centers[j], iy, ".", color="black", ms=2.5)
    ax.set(xlim=(.087, .513), ylim=(-.5, len(sizes) - .5), title=title,
           xlabel="Realized shorter-route fraction", ylabel="Circumference (native units)")
    ax.set_xticks((.1, .2, .3, .4, .5))
    ax.set_yticks(range(len(sizes)), [str(int(s)) for s in sizes])
    ax.tick_params(axis="y", labelsize=7)
    if colorbar and not categorical:
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=.048, pad=.025)
        cbar.ax.tick_params(labelsize=7)
    return norm, cmap


def format_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def figure1(data, basal, sizes):
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 9.1), gridspec_kw={"height_ratios": [1.8, 1.]})
    for col, model in enumerate(MODELS):
        matrix(axes[0, col], data, model, "initial_cue_weighted_gain",
               f"{chr(65 + col)}  {NAMES[model]}: initial cue gain", sizes)
        ax = axes[1, col]
        for epsilon, style, color in ((.001, "--", "#969696"), (.01, ":", "#969696"), (.005, "-", "#1b1b1b")):
            rr = sorted(pick(basal, model=model, epsilon=epsilon), key=lambda r: f(r, "circumference"))
            yy = [f(r, "net_repair_per_damage") if truth(r, "normalized_gain_available") else np.nan for r in rr]
            ax.plot([f(r, "circumference") for r in rr], yy, style, color=color,
                    marker="o" if epsilon == .005 else None, ms=3.3, lw=1.5,
                    label=f"{100 * epsilon:g}% removed")
        ax.set(xlabel="Circumference (native units)", ylabel="Net-rate excess / active-material deficit\n(native time⁻¹)",
               title=f"{chr(68 + col)}  Basal conditional repair")
        ax.set_xticks((10, 19, 30, 39, 50, 65, 80))
        ax.legend(frameon=False, fontsize=7)
        format_axes(ax)
    fig.suptitle("The same ring basals supply input gain and a defined conditional repair response", fontsize=14, y=1.01)
    fig.text(.01, -.015, "Every planned size and placement is retained. Gray ×: preparation/overlap exclusion. Gain is cue-weighted inactive substrate\n(Wave-pinning/Goryachev), or inactive Rac / 7.5 (Holmes). Basal repair has no future-placement dependence. Axes retain native model units.", fontsize=9)
    fig.tight_layout(h_pad=2.1, w_pad=2.)
    save(fig, "ring_functional_size_location_summary", "Full prospective matrix of exact onset cue-RHS gains, with one basal conditional inactive-pool removal/return measurement per model/circumference. Top axes use actual sampled shorter-route fractions, not interpolated requested locations. Gain equals the cue-weighted inactive concentration for Wave-pinning/Goryachev and cue-weighted inactive Rac divided by 7.5 for Holmes; it is not total activation propensity. Each bottom curve uses a 0.1-unit whole-inactive-profile removal and exact pointwise return at 0.1%, 0.5% (primary) or 1%. The 0.5% pairs have tenfold-tighter matched references; signed ratios require both damage and net-rate numerator to exceed 100 times their numerical discrepancies. Gray crosses preserve unavailable matrix cells. Native time units and cue amplitudes differ across models, so numerical rates are not cross-calibrated. Across circumference, material, incumbent architecture and both routes change together.")


def stage_figure(data, stages, sizes, stage, headline=False):
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 13.1))
    for col, model in enumerate(MODELS):
        lookup = {(f(r, "circumference"), f(r, "fraction")): r for r in pick(stages, model=model, stage=stage, epsilon=.005)}
        matrix(axes[0, col], data, model, "net_repair_per_damage", f"{chr(65 + col)}  {NAMES[model]}: relative repair", sizes,
               values=lookup, overlay=True, signed=True)
        matrix(axes[1, col], data, model, "sham_absolute_net_rate", f"{chr(68 + col)}  Fixed original core: sham net rate", sizes,
               values=lookup, overlay=True, signed=True)
        matrix(axes[2, col], data, model, "final_location_association", f"{chr(71 + col)}  Final site association, t = 14000", sizes,
               categorical=True)
    handles = [Patch(facecolor=c, label=l) for c, l in zip(SITE_COLORS, SITE_LABELS)]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(.5, -.003))
    cue = "cue on" if stage < 1000 else "cue off"
    rate_label = "loss" if stage == 10 else "rate"
    fig.suptitle(f"Conditional repair and absolute fixed-core {rate_label} at t = {stage:g} ({cue})", fontsize=15, y=1.005)
    fig.text(.01, .045, "All rates use the fixed original basal core, not a moving site. Repair is read at t + 0.1 after return.\nBlack dot: no resolved distinct old site; the fixed-core response is regional. Middle row: signed amount per native time (red < 0).\nGray ×: preparation/overlap exclusion; white +: gated ratio, raw rate retained. Final flags do not establish stability.", fontsize=9)
    fig.tight_layout(rect=(0, .10, 1, .99), h_pad=1.5, w_pad=1.5)
    stem = "ring_functional_challenge_summary" if headline else f"ring_functional_stage{stage:g}_support"
    save(fig, stem, f"Complete prespecified matrix of primary 0.5% inactive-profile removal/return probes launched at t={stage:g}, read at t={stage + .1:g}, with the ongoing cue {'on' if stage < 1000 else 'off'} in both sham and probe. The top row is excess net core amount-rate divided by sham-relative active-material deficit; gated ratios are withheld while raw values remain in stage_probes.csv. The middle row shows the absolute unperturbed sham net core rate, so positive relative repair can be compared with concurrent decline. All rates use the fixed original basal half-contrast core. Black dots mark conditions without a resolved distinct old site; any single peak is treated as lacking that identity, including a sole peak near the original core. Those probes remain regional responses. The bottom row records final circular-peak location association at 14000 native units using the frozen 0.05 prominence reference and association regions; merged/unresolved states remain explicit. The complete three-stage histories and all perturbation amplitudes are retained. Cross-size and cross-model causal stability rankings are not implied.")


def support_figures(data, stages, basal, sizes):
    # Existing amount-rate expressed relative to the same frozen basal core amount.
    denominators = {(r["model"], f(r, "circumference")): f(r, "basal_core_amount") for r in basal}
    derived = []
    for row in pick(stages, epsilon=.005):
        den = denominators[(row["model"], f(row, "circumference"))]
        derived.append(dict(model=row["model"], circumference=f(row, "circumference"), fraction=f(row, "fraction"),
                            time_from_onset=f(row, "time_from_onset"), location_association=row["location_association"],
                            sham_absolute_net_rate=f(row, "sham_absolute_net_rate"), basal_core_amount=den,
                            fractional_basal_core_net_rate=f(row, "sham_absolute_net_rate") / den))
    path = layout.dev_output_path(TOPIC, "ring_functional_fractional_core_rates.csv")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(derived[0])); writer.writeheader(); writer.writerows(derived)
    ARTIFACTS.append(path)
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 13))
    for iy, stage in enumerate((10, 1000, 2000)):
        for col, model in enumerate(MODELS):
            lookup = {(f(r, "circumference"), f(r, "fraction")): r for r in pick(derived, model=model, stage=stage)}
            matrix(axes[iy, col], data, model, "fractional_basal_core_net_rate",
                   f"{NAMES[model]}, t = {stage}", sizes, values=lookup, signed=True, overlay=True)
    fig.suptitle("Secondary display: core net amount-rate / original basal core amount", y=1.005, fontsize=14)
    fig.tight_layout(h_pad=1.5)
    save(fig, "ring_functional_fractional_core_rates", "Secondary normalization of the retained primary sham net active-material amount-rate by the fixed original basal core amount, at all three prespecified stages. This changes units to a fractional basal-core amount per native time; it introduces neither a fitted rate nor a new stability metric. Black dots mark a regional response without a resolved distinct old site. All measurements use the fixed original core; negative rates need not mean loss from an entire moving or broadening site. The absolute rates remain primary and appear in the stage summary figures. The exact numerator, denominator and derived values are exported in ring_functional_fractional_core_rates.csv.")
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.))
    for col, model in enumerate(MODELS):
        rr = sorted(pick(basal, model=model, epsilon=.005), key=lambda r: f(r, "circumference"))
        for iy, key in enumerate(("basal_core_amount", "basal_width")):
            ax = axes[iy, col]
            ax.plot([f(r, "circumference") for r in rr], [f(r, key) for r in rr], "o-", color="#426b80", ms=3.4)
            ax.set(title=NAMES[model] if iy == 0 else "", xlabel="Circumference (native units)",
                   ylabel="Basal active amount in original core" if iy == 0 else "Basal half-contrast core width")
            format_axes(ax)
    fig.suptitle("Starting organization: amount and width change with circumference", y=1.01)
    fig.tight_layout(h_pad=1.7)
    save(fig, "ring_functional_basal_organization", "Basal primary active material within the fixed half-contrast incumbent core and that core's width, for every qualified circumference in every model. These descriptors use the same exact full-state checkpoints as the onset gains and all placements. Increasing circumference changes conserved total material, core geometry and both routes together; this plot is not an isolated width intervention.")
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 12.))
    for iy, stage in enumerate((10, 1000, 2000)):
        for col, model in enumerate(MODELS):
            lookup = {(f(r, "circumference"), f(r, "fraction")): r for r in pick(stages, model=model, stage=stage, epsilon=.005)}
            matrix(axes[iy, col], data, model, "location_association", f"{NAMES[model]}, t = {stage}",
                   sizes, values=lookup, categorical=True)
    fig.legend(handles=[Patch(facecolor=c, label=l) for c, l in zip(SITE_COLORS, SITE_LABELS)],
               loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle("Spatial identity is checked separately at every probe stage", y=1.005)
    fig.tight_layout(rect=(0, .055, 1, 1))
    save(fig, "ring_functional_stage_identity", "All three prespecified probe-stage circular-peak location associations, using exactly the same fixed prominence and association rules as the final matrices. A single shifted peak remains unresolved; the flags are descriptive regional identities rather than stable-branch labels. Full primary-field kymographs are retained in the referenced cue and recovery arrays.")


def actual_budget(data, stages, sizes):
    """Unperturbed instantaneous native RHS at the exact t10 snapshot."""
    keys = ("snapshot_absolute_reaction_rate", "snapshot_absolute_diffusion_rate", "snapshot_absolute_net_rate")
    assert all(key in stages[0] for key in keys), "Collect exact snapshot allocations before rendering"
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 12.7))
    for col, model in enumerate(MODELS):
        lookup = {(f(r, "circumference"), f(r, "fraction")): r for r in pick(stages, model=model, stage=10, epsilon=.005)}
        for iy, (key, name) in enumerate(zip(keys, ("Reaction", "Diffusion", "Net = reaction + diffusion"))):
            matrix(axes[iy, col], data, model, key, f"{NAMES[model]}: {name}", sizes, values=lookup, signed=True, overlay=True)
    fig.suptitle("Actual old-core budget at exact t = 10: native reaction and diffusion", y=1.005, fontsize=14)
    fig.text(.01, .008, "Unperturbed snapshot at t = 10, before the 0.1-unit probe; amount per native time. Blue: positive; red: negative.\nThe fixed basal core is an open region. Negative diffusion is net outward diffusion across its edges, not a disappearing molecule.\nBlack dot: no resolved distinct old site. The reaction and diffusion maps use their own stated native color scales.", fontsize=9)
    fig.tight_layout(rect=(0, .07, 1, .99), h_pad=1.4)
    save(fig, "ring_functional_actual_t10_budget", "Actual instantaneous primary active-material reaction, diffusion and net amount-rates integrated over the frozen original basal core at the exact unperturbed t10 snapshot, before the removal/return probe is launched. The same complete size/placement grid is retained for every model. These are allocations of the native RHS on coupled trajectories, not independently evolved interventions. Negative diffusion denotes net outward redistribution across the core boundary. Per-model colorbars retain the signed magnitudes; differing native units do not provide a cross-model rate calibration. Black dots identify fixed-core regional measurements without a resolved distinct old site.")


def read_phase(path, offset):
    with np.load(source(path)) as z:
        time = z["time"].copy() + offset
        primary = z["states"][:, 0].copy() if "states" in z else z["primary"].copy()
        x = z["x"].copy()
    return time, x, primary


def kymograph(ax, row, title):
    tc, x, pc = read_phase(row["cue_path"], 0)
    tr, xr, pr = read_phase(row["recovery_path"], 1000)
    np.testing.assert_array_equal(x, xr)
    t, p = np.r_[tc, tr[1:]], np.r_[pc, pr[1:]]
    # Display a split time coordinate so the three prespecified early probes
    # remain visible while the full retained 14000-unit window is still shown.
    transform = lambda a: np.arcsinh(np.asarray(a) / 100.)
    mesh = ax.pcolormesh(x, transform(t), p, shading="nearest", cmap="magma", rasterized=True)
    for stage in (10, 1000, 2000):
        ax.axhline(transform(stage), color="cyan", lw=.65, ls="--")
    ax.axvline(x[0], color="white", lw=.7)
    ax.axvline(x[int(f(row, "cue_center_cell"))], color="#7ed5ed", lw=.7)
    times = (0, 10, 100, 1000, 2000, 5000, 14000)
    ax.set(yticks=transform(times), yticklabels=[str(t) for t in times], xlabel="Ring position (native units)",
           ylabel="Time from cue onset (native units)\nasinh display, scale = 100", title=title)
    plt.colorbar(mesh, ax=ax, fraction=.048, pad=.025, label="Primary active field")


def anchors(data, stages):
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 12.), gridspec_kw={"height_ratios": [2., 1., 1.]})
    for col, model in enumerate(MODELS):
        row = pick(data, model=model, circumference=40, fraction=.4)[0]
        kymograph(axes[0, col], row, f"{NAMES[model]}: C40, d16")
        rr = sorted(pick(stages, model=model, circumference=40, fraction=.4, epsilon=.005), key=lambda r: f(r, "time_from_onset"))
        xx = np.arange(3)
        yy = [f(r, "net_repair_per_damage") if truth(r, "normalized_gain_available") else np.nan for r in rr]
        axes[1, col].plot(xx, yy, "o-", color="#4d7091")
        for j, r in enumerate(rr):
            if not old_associated(r["location_association"]):
                axes[1, col].plot(j, yy[j], "o", mfc="none", mec="black", ms=10)
        axes[1, col].set(ylabel="Excess net rate / deficit\n(native time⁻¹)", title="Primary relative repair (open ring = regional)")
        axes[2, col].plot(xx, [f(r, "sham_absolute_net_rate") for r in rr], "o-", label="Sham", color="#202020")
        axes[2, col].plot(xx, [f(r, "probe_absolute_net_rate") for r in rr], "s--", label="Probe after return", color="#b7683c")
        axes[2, col].axhline(0, color="gray", lw=.7)
        axes[2, col].set(ylabel="Absolute core net amount-rate", title="Both absolute rates retained")
        for iy in (1, 2):
            axes[iy, col].set_xticks(xx, ["10\ncue on", "1000\ncue off", "2000\ncue off"])
            axes[iy, col].set_xlabel("Probe launch time; return at +0.1")
            format_axes(axes[iy, col])
        axes[2, col].legend(frameon=False, fontsize=8)
    fig.suptitle("Prespecified same-geometry anchor: native trajectories and their matching conditional probes", fontsize=13, y=1.005)
    fig.tight_layout(h_pad=1.5)
    save(fig, "ring_functional_C40_d16_anchors", "Prespecified C40 shorter-route d16 (fraction 0.4) primary-field histories for all three models, with full 14000-unit windows shown on an explicitly asinh-transformed time coordinate. Dashed lines mark the t=10,1000,2000 probe launches; the cyan vertical line marks the challenger center and the left edge marks the original peak origin. The old site crosses the periodic cut: activity bands at x=0 and x=C belong to the same site. Each lower curve uses that exact trajectory snapshot; site identity comes from that native snapshot, not the probe. Open rings flag relative rates without a resolved distinct old site, so they are regional responses. The absolute sham and probe rates reveal cases in which positive deficit decay coexists with fixed-original-core material loss. Native model units differ; colors and rates are not calibrated across models.")
    fig, axes = plt.subplots(1, 2, figsize=(10, 6))
    for ax, fraction in zip(axes, (.3, .4)):
        row = pick(data, model="holmes_model3", circumference=40, fraction=fraction)[0]
        kymograph(ax, row, f"Holmes 3: C40, d{f(row, 'short_route'):g}")
    fig.suptitle("Prespecified Holmes location control: exact shared full basal state", y=1.01)
    fig.tight_layout()
    save(fig, "ring_functional_Holmes_C40_location_control", "Prespecified Holmes C40 d12 versus d16 location control. Both challenges use the exact same complete native Rac/Rho basal state, native kinetics, width-3 cue of peak 1, 1000-unit cue and 13000-unit recovery; only placement changes. The full primary histories and three probe stages are shown. Activity at x=0 and x=C is one site crossing the periodic cut. This fixed-circumference comparison isolates cue location within this ring protocol; it does not isolate one route, since moving the cue changes both routes, nor establish a universal distance-delay mechanism.")


def goryachev_history(data):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9), gridspec_kw={"height_ratios": [2., 1.]})
    for col, circumference in enumerate((40, 80)):
        row = pick(data, model="goryachev", circumference=circumference, fraction=.5)[0]
        kymograph(axes[0, col], row, f"Goryachev: C{circumference}, antipodal d{circumference//2}")
        tc, x, pc = read_phase(row["cue_path"], 0)
        tr, xr, pr = read_phase(row["recovery_path"], 1000)
        t, p = np.r_[tc, tr[1:]], np.r_[pc, pr[1:]]
        cue_cell = int(f(row, "cue_center_cell"))
        ax = axes[1, col]
        ax.plot(np.arcsinh(t/100), p[:, 0], color="#4678a3", label="Original peak coordinate")
        ax.plot(np.arcsinh(t/100), p[:, cue_cell], color="#d77c35", label="Cue center")
        ax.axvline(np.arcsinh(1000/100), color="gray", ls="--", lw=.8)
        times = (0, 10, 100, 1000, 2000, 5000, 14000)
        ax.set(xticks=np.arcsinh(np.array(times)/100), xticklabels=[str(t) for t in times],
               xlabel="Time from cue onset (native units; asinh display)", ylabel="Active field at fixed coordinate")
        ax.legend(frameon=False, fontsize=8)
        format_axes(ax)
    fig.suptitle("Goryachev histories resolve two geographically different single-peak endpoints", y=1.005, fontsize=13)
    fig.tight_layout(h_pad=1.7)
    save(fig, "ring_functional_Goryachev_C40_C80_antipodal_history", "Supporting history comparison chosen after the full expanded matrix exposed geographically different single-peak outcomes. Both Goryachev conditions use the prespecified antipodal fraction 0.5, unchanged native kinetics and cue protocol, but circumference, total material, architecture and both absolute routes differ together. At C40 activity relocates toward the challenged site; at C80 activity remains near the original site throughout, while the challenged-site transient declines after withdrawal. The lower traces are primary-field values at two fixed coordinates, not tracked peak amplitudes or whole-site material. Full kymographs provide spatial context; the original region crosses the periodic cut at 0/C. This pair supports a retained old-region interpretation for this C80 history; the complete matrix remains labelled by geography rather than automatically assigning tracked fates to every sole peak. It does not isolate a transport route or a causal size effect.")


def mdtable(headers, values):
    return "| " + " | ".join(headers) + " |\n| " + " | ".join("---" for _ in headers) + " |\n" + "\n".join("| " + " | ".join(str(v) for v in row) + " |" for row in values)


def fmt(value):
    return f"{value:.5g}" if np.isfinite(value) else "unavailable"


def write_report(data, basal, stages, audit, run):
    generated = now()
    counts = audit["counts"]
    size_rows, gain_rows, stage_rows, fold_rows, budget_rows, anchor_rows = [], [], [], [], [], []
    for model in MODELS:
        bb = sorted(pick(basal, model=model, epsilon=.005), key=lambda r: f(r, "circumference"))
        for circumference in (19, 39):
            b = pick(bb, circumference=circumference)
            rr = pick(data, model=model, circumference=circumference, fraction=.4)
            if not b or not rr or rr[0]["status"] != "complete": continue
            b, r = b[0], rr[0]
            size_rows.append((NAMES[model], circumference, fmt(f(b, "basal_core_amount")), fmt(f(b, "basal_width")),
                              fmt(f(r, "initial_cue_weighted_gain")), fmt(f(b, "net_repair_per_damage"))))
            q = pick(stages, model=model, circumference=circumference, fraction=.5, stage=10, epsilon=.005)[0]
            budget_rows.append((NAMES[model], circumference, fmt(f(q,"snapshot_absolute_reaction_rate")),
                                fmt(f(q,"snapshot_absolute_diffusion_rate")), fmt(f(q,"snapshot_absolute_net_rate")),
                                fmt(f(q,"snapshot_absolute_net_rate")/f(q,"basal_core_amount"))))
        b19, b39 = (pick(bb, circumference=c)[0] for c in (19, 39))
        r19, r39 = (pick(data, model=model, circumference=c, fraction=.4)[0] for c in (19, 39))
        fold_rows.append((NAMES[model], fmt(f(r39,"initial_cue_weighted_gain")/f(r19,"initial_cue_weighted_gain")),
                          fmt(f(b39,"net_repair_per_damage")/f(b19,"net_repair_per_damage")),
                          fmt(f(b39,"basal_core_amount")/f(b19,"basal_core_amount"))))
        anchor = pick(data, model=model, circumference=40, fraction=.4)[0]
        anchor_rows.append((NAMES[model], fmt(f(anchor,"initial_cue_weighted_gain")),
                            fmt(f(anchor,"observed_patch_mean_extra_at_1")), fmt(f(anchor,"observed_patch_mean_extra_at_10")),
                            anchor["final_location_association"]))
        valid = [r for r in pick(data, model=model) if r["status"] == "complete"]
        for c in sorted({f(r, "circumference") for r in valid}):
            values = [f(r, "initial_cue_weighted_gain") for r in valid if f(r, "circumference") == c]
            gain_rows.append(dict(model=model, circumference=c, gain_min=min(values), gain_max=max(values),
                                  gain_placement_range=max(values)-min(values)))
        rr = pick(stages, model=model, stage=10, epsilon=.005)
        pos = [r for r in rr if truth(r, "normalized_gain_available") and f(r, "net_repair_per_damage") > 0]
        conflict = [r for r in pos if f(r, "sham_absolute_net_rate") < 0 and f(r, "probe_absolute_net_rate") < 0]
        stage_rows.append((NAMES[model], len(rr), len(pos), len(conflict),
                           sum(not old_associated(r["location_association"]) for r in rr)))
    # A small direct endpoint table makes the complete curves reviewable without
    # promoting a monotonic fit or inventing a model-wide stability ranking.
    range_rows = []
    for model in MODELS:
        bb = sorted(pick(basal, model=model, epsilon=.005), key=lambda r: f(r, "circumference"))
        lo, hi = bb[0], bb[-1]
        gainlo = pick(data, model=model, circumference=f(lo, "circumference"), fraction=.5)[0]
        gainhi = pick(data, model=model, circumference=f(hi, "circumference"), fraction=.5)[0]
        range_rows.append((NAMES[model], f"{f(lo,'circumference'):g} → {f(hi,'circumference'):g}",
                           f"{fmt(f(gainlo,'initial_cue_weighted_gain'))} → {fmt(f(gainhi,'initial_cue_weighted_gain'))}",
                           f"{fmt(f(lo,'net_repair_per_damage'))} → {fmt(f(hi,'net_repair_per_damage'))}",
                           f"{fmt(f(lo,'basal_core_amount'))} → {fmt(f(hi,'basal_core_amount'))}"))
    primary = pick(stages, epsilon=.005)
    missing = sum(not truth(r, "normalized_gain_available") for r in primary)
    text = f"""# Native ring functional comparison: positive conditional repair during fixed-core loss

Generated {generated}. Producer: completed; collection audit passed at {audit['generated_at']}. This development report preserves the main user-reviewed integrated figure. The [v1 plan](PLAN.md), [expanded plan](PLAN_v2.md) and [prospective design critique](supporting/metadata/DESIGN_CRITIQUE.md) precede these outcomes.

At the prespecified early time 10, **all 652 challenged conditions have a positive reliable conditional repair coefficient while the unperturbed fixed original core is losing active material**. In 621 conditions both the sham and the probe still decline at 10.1; in 31 the perturbed return has a positive absolute rate. All amount-rates integrate the fixed original basal half-contrast core; they do not follow an entire translating or broadening site. The positive coefficient measures a restoring contribution relative to the sham along this inactive-pool perturbation. It does not imply net core maintenance, survival of a distinct incumbent, or explain the final spatial state by itself.

The complete ring matrix connects three quantities measured from the **same native ring states**: the immediate challenger input gain, the response to a specified inactive-pool disturbance, and subsequent spatial competition. All three models are necessary: Goryachev couples a substantial reduction in incremental cue gain to a strong increase in basal conditional repair, whereas Wave-pinning and Holmes have nearly unchanged onset gain over C19–39 and much smaller coefficient changes. Their spatial responses still change with circumference and placement. Fixed-circumference location comparisons share the exact full basal state. Across circumference, total material, incumbent architecture, both physical routes and the fraction of the ring occupied by the width-3 cue change together.

![Full initial gain maps and one basal conditional-repair curve per size](ring_functional_size_location_summary.png)

Figure 1. Full initial-gain maps and basal conditional-repair curves use every eligible planned placement. The onset gain is exactly cue-weighted inactive substrate for Wave-pinning/Goryachev, or inactive Rac divided by 7.5 for Holmes. It measures the incremental cue input at the starting state, not total activation propensity. Across eligible placements within any one size, its maximum relative range is only 0.0392%, 0.0891% and 0.00208%, respectively. The full maps therefore expose nearly placement-independent starting substrate despite the later placement dependence. Narrow model-specific color ranges and y axes must not be read as comparable enhancement magnitudes. The 0.1%, 0.5% and 1% curves probe one declared disturbance direction. Basal repair is measured once per model/size and has no distance dependence.

Endpoint values below summarize the displayed full curves; they are not fitted scaling laws. The antipodal fraction 0.5 is used for the gain endpoint comparison to avoid near-core exclusions. Rates retain native model units and the cue peaks differ (0.1, 6.4 and 1), so their numerical magnitudes do not rank models.

{mdtable(('Model','Qualified C range','Onset gain at fraction 0.5','Primary basal repair','Basal core active amount'),range_rows)}

The 19-to-39 comparison below is recalculated in the ring geometry. These are not imported interval values or the same full checkpoints as the separate integrated figure. It provides a transparent size anchor while the complete curves expose the remainder of the range.

{mdtable(('Model','C','Core amount','Core width','Onset gain, fraction 0.4','Primary basal repair'),size_rows)}

For the ring-native C19 → C39 comparison, the following ratios quantify the different trends directly. Basal core amount roughly doubles in all three, whereas the conditional coefficient increases 4.27-fold only in Goryachev. Wave-pinning increases 1.16-fold and Holmes 1.015-fold. The assay therefore supports a model-specific conditional restoring contribution; it does not supply a common intrinsic stability law.

{mdtable(('Model','Onset gain ratio, fraction 0.4','Basal repair ratio','Basal core amount ratio'),fold_rows)}

![Complete challenged-stage repair, absolute core rate and final site maps](ring_functional_challenge_summary.png)

Figure 2. Every active condition is shown at the prespecified early time 10 with cue on. The repair readout is taken at 10.1, immediately after exact pointwise return. The middle row retains the absolute sham core net amount-rate; all these values are negative. The final map describes the spatial state at 14000; a location flag does not establish a tracked fate or stable branch. A black dot in the rate maps marks the absence of a resolved distinct old site, so the fixed-core probe is a regional response. This rule applies to **every single-peak state**, including a sole peak near the original core and one connected region spanning the original and challenged sites.

{mdtable(('Model','t10 primary probes','Positive gated repair','Positive repair, both absolute rates negative','No resolved distinct old site at t10'),stage_rows)}

The [complete exact-t10 native budget maps](ring_functional_actual_t10_budget.png) explain the sign of actual core loss without equating it to the conditional probe difference. At the unperturbed t10 snapshot, diffusion makes a negative contribution in all 652 conditions. Reaction is positive in 500 conditions but smaller than outward diffusion; in 152, reaction also removes active material. The original core is an open region: negative diffusion is redistribution across its boundary, not whole-ring loss. These terms allocate the coupled native RHS; they are not independently evolved pathway ablations.

For the prespecified antipodal placement (fraction 0.5), the same allocation is shown at C19 and C39. At C39 reaction contributes positively in all three models, while the larger outward diffusive term leaves a negative net rate. Absolute net decline and decline relative to starting core amount both become smaller in these six size anchors. This is compatible with a larger material capacity and a conditional restoring contribution, while leaving their separate causal effects unresolved.

{mdtable(('Model','C','Actual reaction at t10','Actual diffusion at t10','Actual net at t10','Net / basal core amount'),budget_rows)}

The displayed disagreement between relative repair and absolute loss is part of the result. It prevents interpreting a positive conditional coefficient as generic resistance to competition, or treating that coefficient alone as an explanation of final location. [Absolute and fractional core-rate support](ring_functional_fractional_core_rates.png) expresses the retained sham amount-rate per original basal core amount; [its CSV](supporting/data/ring_functional_fractional_core_rates.csv) contains the numerator and denominator. This secondary unit change can describe how much of the starting core is being lost per native time; it is not a new stability metric. [Basal amount and core-width curves](ring_functional_basal_organization.png) expose the changing starting architecture.

The final spatial map retains contrasts that a universal maintenance claim would obscure. Wave-pinning has cue-associated single peaks at all eligible C10–39 placements, but also distinct old- and cue-associated peaks for many C50–80 placements. Holmes has cue-associated peaks at C10–25. Its C30 row is transitional: seven cue-associated, two spanning and one shifted/unresolved state; distinct old/cue peaks or merged/bridged regions occur from C35 through C75. At C40 its d10–13 conditions form one region spanning both sites, whereas d14–20 retain distinct associated peaks. Goryachev has cue-associated single peaks through C50; at C55–65 placement separates cue-associated peaks from sole old-core-associated peaks, and at C70–80 every tested placement has a sole old-core-associated peak whose incumbent identity remains unresolved by that geographic flag alone. These are finite-window geographic classifications, not a common succession of stable states.

The [Goryachev C40 versus C80 antipodal histories](ring_functional_Goryachev_C40_C80_antipodal_history.png) inspect this matrix contrast directly. At C40 activity relocates toward the cue; at C80 the original region remains active throughout while the challenged-site transient disappears after withdrawal. The full C80 history therefore supports retained original-region activity beyond the endpoint's geographic label. This supporting pair was selected after the expanded map exposed the contrast; the prespecified C40 d16 anchors below remain unchanged. The comparison jointly changes size, material, incumbent architecture and both absolute routes and does not isolate which change causes the different competition result.

The same systematic maps remain available at [1000 immediately after withdrawal](ring_functional_stage1000_support.png) and [2000, after 1000 cue-free units](ring_functional_stage2000_support.png). Their probes return at 1000.1 and 2000.1. [All three stage-identity maps](ring_functional_stage_identity.png) keep disappearance, merger and ambiguous association separate from the numerical repair readout. Neither a regional probe after site loss nor a sole shifted maximum establishes incumbent stability.

![Prespecified exact C40 d16 trajectories and matching probe rates for all three models](ring_functional_C40_d16_anchors.png)

The prespecified C40 d16 anchor connects full primary-field histories to its own three probes for every model. The full 14000-unit histories use a labelled asinh time axis so early stages remain legible. The original site crosses the periodic cut, so bands at positions 0 and C are the same site. These examples complement the complete maps; they do not supply the matrix selection rule. The separate [Holmes C40 d12 versus d16 location control](ring_functional_Holmes_C40_location_control.png) starts from an exact shared complete Rac/Rho state. It demonstrates the consequence of moving the cue in this fixed ring protocol, while changing both routes together.

The table links each exact C40 d16 onset gain to its native 1- and 10-unit active-field increments in the fixed width-3 cue patch (challenged trajectory minus the same basal state evolved without cue), and then to its final spatial association. A separate [finite-response audit]({os.path.relpath(run/'finite_response_audit/audit.json',TOPIC)}) retains 50 exact matched no-cue controls and independently reproduces all 1304 worker-recorded increments exactly. Finite responses evolve all coupled native fields; dividing by the different cue peaks would not calibrate them across models.

{mdtable(('Model','Exact onset gain','Patch mean increment at 1','Patch mean increment at 10','Final location association'),anchor_rows)}

## Assay and numerical boundary

The challenger is a width-3 raised cosine, on for 1000 then off for 13000 native units, with unchanged native kinetics, diffusion and pool densities. Cue centers use the nearest sampled cell at mesh spacing 0.05. Maps use realized route fractions; the CSV retains requested fractions, both actual routes and both edge-to-cue gaps. The probe withdraws 0.5% of the whole current inactive-primary profile, evolves the full native model for 0.1, and returns exactly the removed profile pointwise. A same-start sham evolves under the same ongoing cue condition. Holmes evolves all four Rac/Rho fields with the native Eq. (5)–(6) factor-two implementation and basal Rac activation 3.5. Native time units and input strengths are model-specific.

The active-material deficit D is sham minus probe in the fixed original basal half-contrast core, after return. Conditional repair is the probe-minus-sham **net** core amount-rate divided by D. Reaction and diffusion allocations, both absolute rates, raw ratio estimates and all amplitude rows are retained. Every primary pair has a tenfold-tighter matched reference. Damage and numerator must each exceed 100 times the corresponding numerical discrepancy; damage also exceeds 10⁻¹¹ of the basal core amount. {missing} of {len(primary)} primary stage ratios are withheld, all in Goryachev (83 at 1000 and 130 at 2000); raw measurements remain available. All 652 primary t10 ratios pass. These gates qualify subtraction accuracy, not continuum convergence or generic stability.

Sparse primary-field histories and complete predefined native snapshots are retained for new conditions. They do not constitute continuously sampled four-field Holmes mechanistic histories. This experiment tests neither freely migrating cells nor biological persistence, and does not isolate width, a single transport route or a causal partition of replacement time.

## Coverage and reproducibility

The grid contains {audit['planned_matrix_cells']} planned cells from {audit['model_size_workers']} model/size workers: {', '.join(str(v)+' '+k for k,v in counts.items())}. The unavailable 17-cell row is Holmes C80, whose fixed preparation produces multiple sites and does not qualify as the declared single-incumbent starting state. Preparation failures and cue/core-overlap exclusions remain explicit rows; neither is a biological replacement failure. The 36 exact matching v1 conditions and shared basals are reused by immutable source reference.

- [Complete matrix CSV]({os.path.relpath(run/'results/matrix.csv',TOPIC)}) — geometry, onset and short-time responses, final identity, exact trajectory paths and stage summaries.
- [Basal probe CSV]({os.path.relpath(run/'results/basal.csv',TOPIC)}) — one measurement per size at each of three removal amplitudes.
- [All stage probes]({os.path.relpath(run/'results/stage_probes.csv',TOPIC)}) — all stages and amplitudes, raw rates, damage, error bounds and unavailable reasons.
- [Independent reconstruction audit]({os.path.relpath(run/'results/audit.json',TOPIC)}) and [retained source hashes]({os.path.relpath(run/'results/source_sha256.json',TOPIC)}).
- [Vector exports](supporting/vectors/) and [report artifact manifest](supporting/metadata/ring_functional_report.manifest.json).

Reproduce the exports after collection with the existing environment:

```bash
/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/chemical_ring_functional_report_20261017.py
```

The reconstruction audit verifies full planned-cell coverage, exact starting checkpoints and cue/recovery continuity, the native cue-RHS identity, conserved-pool accounting and each recorded probe directly from retained native arrays. The separate matched-control audit verifies all finite 1- and 10-unit increments. Their maximum errors and scopes are in the linked audits. The [independent scientific result critique](supporting/metadata/RESULT_CRITIQUE.md) checks the metrics, numerical gates, native model choices, full spatial comparison and interpretation limits; these remain development figures rather than a promoted manuscript conclusion.
"""
    path = TOPIC / "RESULTS.md"; path.write_text(text); ARTIFACTS.append(path)
    path = layout.dev_output_path(TOPIC, "ring_functional_gain_placement_ranges.csv")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(gain_rows[0])); writer.writeheader(); writer.writerows(gain_rows)
    ARTIFACTS.append(path)


def main(run, support_refresh=False):
    if support_refresh:
        previous = json.loads((TOPIC / "supporting/metadata/ring_functional_report.manifest.json").read_text())
        for path, digest in previous["source_sha256"].items():
            if Path(path).resolve() != Path(__file__).resolve():
                assert sha(source(path)) == digest, path
        ARTIFACTS.extend(TOPIC / name for name in previous["sha256"])
    audit = json.loads(source(run / "results/audit.json").read_text())
    assert audit["passed"] and audit["full_867_cell_schema_coverage"]
    finite_path = source(run / "finite_response_audit/audit.json")
    assert sha(finite_path) == audit["finite_time_response_validation"]["sha256"]
    finite = json.loads(finite_path.read_text())
    assert finite["passed"] and finite["finite_time_readouts"] == 1304
    for path, digest in finite["source_sha256"].items():
        assert sha(source(path)) == digest
    data, basal, stages = (rows(run / "results" / f"{name}.csv") for name in ("matrix", "basal", "stage_probes"))
    assert len(data) == 867
    source(__file__)
    sizes = sorted({f(r, "circumference") for r in data})
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                         "svg.fonttype": "none", "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
    if not support_refresh:
        figure1(data, basal, sizes)
        stage_figure(data, stages, sizes, 10, headline=True)
        for stage in (1000, 2000): stage_figure(data, stages, sizes, stage)
        support_figures(data, stages, basal, sizes)
        anchors(data, stages)
    actual_budget(data, stages, sizes)
    goryachev_history(data)
    write_report(data, basal, stages, audit, run)
    ARTIFACTS[:] = list(dict.fromkeys(ARTIFACTS))
    links = 0
    for target in re.findall(r'\]\(([^)]+)\)', (TOPIC / "RESULTS.md").read_text()):
        if target != "supporting/metadata/ring_functional_report.manifest.json":
            assert (TOPIC / target).exists(), target
        links += 1
    record = dict(generated_at=now(), artifact_root="../..", completed_data_audit=audit,
                  sha256={str(p.relative_to(TOPIC)): sha(p) for p in ARTIFACTS},
                  source_sha256=SOURCES, local_links_checked=links,
                  simulation_and_frozen_sources_unchanged=True)
    path = layout.dev_output_path(TOPIC, "ring_functional_report.manifest.json")
    path.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"generated_at":record["generated_at"], "artifacts":len(ARTIFACTS), "local_links":links,
                      "report":str(TOPIC / "RESULTS.md")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, default=DEFAULT)
    parser.add_argument("--support-refresh", action="store_true", help="Retain verified existing headline exports, refresh budget/history support and prose")
    args = parser.parse_args()
    main(args.run.resolve(), support_refresh=args.support_refresh)
