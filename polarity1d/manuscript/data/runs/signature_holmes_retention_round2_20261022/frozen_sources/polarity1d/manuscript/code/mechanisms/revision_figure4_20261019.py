"""Read-only reanalysis and compact manuscript Figure 4 from frozen noisy rings.

No native integration, changed forcing, fitted event classifier or site ancestry.
The common-circumference projection acts on the actual endpoint distribution.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[4]
RUN = ROOT / "polarity1d/manuscript/data/runs/vision_noise_main_20261017"
OLD = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/vision_noise_20261017"
OUT = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4"
DATA = OUT / "supporting/data"
META = OUT / "supporting/metadata"
VECTORS = OUT / "supporting/vectors"
MODELS = ("wave_pinning", "holmes_model3", "goryachev")
NAMES = dict(zip(MODELS, ("Wave-pinning", "Holmes 3", "Goryachev")))
COLORS = {20: "#397bab", 40: "#40916c", 60: "#cf7930"}
MODEL_COLORS = {"wave_pinning": "#66549a", "holmes_model3": "#bc5374", "goryachev": "#398c8c"}
LAGS = (5, 40, 100, 500, 2000, 8000)
DRAW_COUNT = 2000
BOOTSTRAP_SEED = 2026101904
C_REFERENCE = 60.0


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def read_rows(path, predicate=None):
    with Path(path).open(newline="") as f:
        return [r for r in csv.DictReader(f) if predicate is None or predicate(r)]


def write_rows(name, rows):
    p = DATA / (name + ".csv")
    with p.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)
    return p


def select(rows, **query):
    return [r for r in rows if all(r[k] == v if isinstance(v, str) else float(r[k]) == float(v)
                                  for k, v in query.items())]


def safe_ratio(a, b):
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    return np.divide(a, b, out=np.full(a.shape, np.nan), where=b > 0)


def interval(draws):
    d = np.asarray(draws)
    return np.nanpercentile(d, (2.5, 97.5), axis=0) if np.isfinite(d).any() else (np.nan, np.nan)


def weights(n, rng):
    if n < 2:
        return np.empty((0, n))
    return rng.multinomial(n, np.full(n, 1 / n), DRAW_COUNT) / n


def summary(values, w, denominators=None):
    x = np.asarray(values, float)
    if denominators is None:
        estimate = np.mean(x, axis=0)
        draws = w @ x
    else:
        d = np.asarray(denominators, float)
        estimate = safe_ratio(np.sum(x, axis=0), np.sum(d, axis=0))
        draws = safe_ratio(w @ x, w @ d)
    lo, hi = interval(draws) if len(w) else (np.nan, np.nan)
    return dict(value=float(estimate), ci_lower=float(lo), ci_upper=float(hi),
                bootstrap_defined_fraction=float(np.isfinite(draws).mean()) if len(w) else np.nan)


def save(fig, stem):
    paths = []
    for suffix in ("png", "pdf", "svg"):
        p = (OUT if suffix == "png" else VECTORS) / f"{stem}.{suffix}"
        fig.savefig(p, dpi=200 if suffix == "png" else None, bbox_inches="tight")
        paths.append(p)
    plt.close(fig)
    return paths


def selfcheck():
    # Exact circular translation: geometry remapping changes no displacement.
    C, d = 20., np.array([0., 1., 3., -8., 9.])
    angle = 2 * np.pi * d / C
    np.testing.assert_allclose(C * np.angle(np.exp(1j * angle)) / (2 * np.pi), d)
    np.testing.assert_allclose(np.cos(2 * np.pi * d / 60), np.cos(angle / 3))
    # A dominant-site label can jump while a field moment is unchanged.
    theta = 2 * np.pi * (np.arange(600) + .5) / 600
    antipodal = 2 + np.cos(2 * theta)
    assert abs(np.sum(antipodal * np.exp(1j * theta))) < 1e-10
    assert np.sum(antipodal) > 0
    # Bootstrap operates on whole-seed sums/counts, preserving missing pairs.
    w = np.array([[1., 0.], [0., 1.], [.5, .5]])
    np.testing.assert_allclose(safe_ratio(w @ np.array([2., 0.]), w @ np.array([4., 0.]))[[0, 2]], .5)
    assert np.isnan(safe_ratio(0, 0))
    return {"circular_conversion": "pass", "cancellation_fixture": "pass", "paired_seed_resampling": "pass"}


def analyze():
    manifest = json.loads((RUN / "manifest.json").read_text())
    assert len(manifest["jobs"]) == 243
    old_seed = {r["case_id"]: r for r in read_rows(OLD / "supporting/data/vision_noise_seed_summary.csv")}
    old_audit = json.loads((OLD / "supporting/metadata/vision_noise_analysis_audit.json").read_text())
    assert old_audit["cases"] == 243
    assert old_audit["index_sha256"] == sha(RUN / "manifest.json")
    samples, pairs, increments, cdf_parts, sources = [], [], [], [], []
    trace_cache = {}
    cdf_thresholds = np.r_[0., np.geomspace(1e-6, np.pi, 181)]
    for job in manifest["jobs"]:
        folder = RUN / "workers" / f"{job['slot']:04d}"
        status = json.loads((folder / "status.json").read_text())
        assert status["complete"] and status["status"] == "complete" and status["native_time"] == 32000
        trace_path = folder / "trace.npz"
        trace_hash = sha(trace_path)
        assert trace_hash == status["output_sha256"]["trace.npz"]
        assert trace_hash == old_seed[f"slot_{job['slot']:04d}"]["archive_sha256"]
        sources.append(dict(slot=job["slot"], condition=job["kind"], model=job["model"],
                            circumference=job["circumference"], fraction=job["mean_fraction"], seed=job["seed"],
                            path=str(trace_path), sha256=trace_hash,
                            status_path=str(folder / "status.json"), status_sha256=sha(folder / "status.json")))
        with np.load(trace_path) as z:
            t, P = z["time"], z["P_real"] + 1j * z["P_imag"]
            R, reg, peaks = abs(P), z["current_half_contrast_regions"], z["peaks_prominence_ge_005"]
            contrast, pos = z["contrast"], z["dominant_region_position"]
            q1, q2 = z["q1"], z["q2"]
        np.testing.assert_array_equal(t, np.arange(0, 32005, 5))
        base = dict(model=job["model"], circumference=job["circumference"], fraction=job["mean_fraction"],
                    condition=job["kind"], seed=job["seed"], slot=job["slot"], observations=len(t))
        good_direction = R >= .1
        one_region = (reg == 1) & (contrast >= .05) & good_direction
        strict = one_region & (peaks <= 1)
        samples.append(dict(**base, weak_R_fraction=float(np.mean(R < .1)),
                            multiregion_fraction=float(np.mean(reg > 1)),
                            strict_site_fraction=float(np.mean(strict)),
                            one_region_fraction=float(np.mean(one_region)),
                            mean_R=float(np.mean(R)), final_R=float(R[-1]),
                            mean_q1=float(np.mean(q1)), mean_q2=float(np.mean(q2)),
                            minimum_R=float(np.min(R)),
                            mean_q2_when_weak=float(np.mean(q2[R < .1])) if (R < .1).any() else np.nan))
        C = job["circumference"]
        for lag in LAGS:
            k = lag // 5
            dtheta = np.angle(P[k:] * P[:-k].conj())
            m = np.minimum(R[k:], R[:-k])
            for gate in (.05, .1, .2):
                v = m >= gate
                angles = abs(dtheta[v])
                row = dict(**base, lag=lag, gate=gate, pairs_available=len(v), pairs_eligible=int(v.sum()),
                           eligible_fraction=float(np.mean(v)),
                           angular_second_moment_sum=float(np.sum(angles**2)),
                           absolute_angle_sum=float(angles.sum()),
                           angle_q50=float(np.quantile(angles, .5)) if len(angles) else np.nan,
                           angle_q95=float(np.quantile(angles, .95)) if len(angles) else np.nan,
                           angle_q99=float(np.quantile(angles, .99)) if len(angles) else np.nan,
                           max_absolute_angle=float(angles.max()) if len(angles) else np.nan,
                           diagnostic_above_pi_over_4=int(np.sum(angles > np.pi / 4)),
                           interpretation="endpoint-conditioned increments; no event class")
                increments.append(row)
                if job["kind"] == "noise" and lag in (5, 40, 100, 500):
                    cdf_parts.append((base, lag, gate, np.searchsorted(np.sort(angles), cdf_thresholds, side="right"), len(angles), len(v)))
            ds = (pos[k:] - pos[:-k] + C / 2) % C - C / 2
            for gate_name, gate in (("strict_resolved_site", strict), ("one_halfcontrast_region", one_region)):
                v = gate[k:] & gate[:-k]
                native = np.cos(2 * np.pi * ds[v] / C)
                reference = np.cos(2 * np.pi * ds[v] / C_REFERENCE)
                pairs.append(dict(**base, lag=lag, gate=gate_name, pairs_available=len(v),
                                  pairs_eligible=int(v.sum()), eligible_fraction=float(np.mean(v)),
                                  displacement_squared_sum=float(np.sum(ds[v]**2)),
                                  native_cosine_sum=float(native.sum()), reference_cosine_sum=float(reference.sum()),
                                  geometry_remapping_sum=float((reference - native).sum()),
                                  reference_circumference=C_REFERENCE,
                                  same_site_ancestry_assumed=False))
        trace_cache[job["slot"]] = dict(time=t, P=P, R=R, regions=reg, strict=strict, pos=pos)
    write_rows("seed_observation_summary", samples)
    write_rows("seed_endpoint_statistics", pairs)
    write_rows("seed_angular_increments", increments)
    write_rows("native_sources", sources)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    obs_group, endpoint_group, increment_group, cdf_group = [], [], [], []
    for model in MODELS:
        for C in COLORS:
            for fraction in (.001, .01):
                query = dict(model=model, circumference=C, fraction=fraction, condition="noise")
                group = select(samples, **query)
                assert len(group) == 12
                w = weights(12, rng)
                groupbase = dict(**query, n_seeds=12, observations_per_seed=6401,
                                 bootstrap_draws=DRAW_COUNT, uncertainty_unit="whole independent seed")
                for metric in ("weak_R_fraction", "multiregion_fraction", "strict_site_fraction", "one_region_fraction", "mean_R", "final_R", "mean_q1", "mean_q2"):
                    obs_group.append(dict(**groupbase, metric=metric, **summary([r[metric] for r in group], w)))
                for lag in LAGS:
                    for gate in ("strict_resolved_site", "one_halfcontrast_region"):
                        pp = select(pairs, **query, lag=lag, gate=gate)
                        count = np.array([r["pairs_eligible"] for r in pp])
                        common = dict(**groupbase, lag=lag, gate=gate, reference_circumference=C_REFERENCE,
                                      eligible_pairs=int(count.sum()), pairs_available=sum(r["pairs_available"] for r in pp),
                                      eligible_pair_fraction=float(sum(count) / sum(r["pairs_available"] for r in pp)),
                                      independent_seeds_with_pairs=int(np.sum(count > 0)))
                        for metric in ("displacement_squared", "native_cosine", "reference_cosine", "geometry_remapping"):
                            endpoint_group.append(dict(**common, metric=metric,
                                                       **summary([r[metric + "_sum"] for r in pp], w, count)))
                    for gate in (.05, .1, .2):
                        pp = select(increments, **query, lag=lag, gate=gate)
                        count = np.array([r["pairs_eligible"] for r in pp])
                        increment_group.append(dict(**groupbase, lag=lag, gate=gate,
                                                    eligible_pair_fraction=float(sum(count) / sum(r["pairs_available"] for r in pp)),
                                                    independent_seeds_with_pairs=int(np.sum(count > 0)),
                                                    metric="angular_second_moment",
                                                    **summary([r["angular_second_moment_sum"] for r in pp], w, count)))
                        if lag in (5, 40, 100, 500):
                            cp = [(a, n, total) for b, l, g, a, n, total in cdf_parts
                                  if all(b[k] == v for k, v in query.items()) and l == lag and g == gate]
                            counts = np.array([n for a, n, total in cp])
                            numerators = np.stack([a for a, n, total in cp])
                            estimate = safe_ratio(numerators.sum(axis=0), counts.sum())
                            draws = safe_ratio(w @ numerators, (w @ counts)[:, None])
                            lo, hi = np.nanpercentile(draws, (2.5, 97.5), axis=0)
                            for j, threshold in enumerate(cdf_thresholds):
                                cdf_group.append(dict(**groupbase, lag=lag, gate=gate, threshold_radians=threshold,
                                                      cdf=estimate[j], ci_lower=lo[j], ci_upper=hi[j],
                                                      eligible_pair_fraction=sum(counts) / sum(total for a, n, total in cp)))
    write_rows("observation_statistics", obs_group)
    write_rows("endpoint_statistics", endpoint_group)
    write_rows("angular_increment_statistics", increment_group)
    write_rows("angular_increment_cdf", cdf_group)
    return manifest, sources, samples, obs_group, endpoint_group, increments, cdf_group, trace_cache


def history(ax, row, trace_cache, number, detail=False):
    with np.load(row["profile_archive"]) as z:
        t, u, x = z["time"], z["primary"], z["x"]
    C = float(row["circumference"])
    position = (x / C + .5) % 1 - .5
    order = np.argsort(position)
    mesh = ax.pcolormesh(position[order], t / 1000, u[:, order], shading="nearest", cmap="viridis", vmin=0, vmax=float(u.max()), rasterized=True)
    ax.set_ylim(32, 0)
    ax.set_xlim(-.5, .5)
    ax.set_xlabel("Boundary position / C")
    ax.set_ylabel("Time (10³ native units)")
    title = f"{number}  {NAMES[row['model']]}, C{C:g}"
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    cb = plt.colorbar(mesh, ax=ax, fraction=.042, pad=.02)
    cb.set_label("Active field", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.text(.02, .04, f"seed {row['seed']}", transform=ax.transAxes, color="white", fontsize=7,
            bbox=dict(facecolor="black", alpha=.35, edgecolor="none", pad=2))
    if detail:
        ax.set_title(title + f", f={float(row['fraction']):g}", loc="left", fontsize=8)
    return dict(profile_path=row["profile_archive"], profile_sha256=sha(row["profile_archive"]),
                seed=row["seed"], model=row["model"], circumference=C,
                fraction=float(row["fraction"]), role=row["role"])


def plot_main(obs, endpoints, trace_cache):
    representatives = read_rows(OLD / "supporting/data/vision_noise_representative_histories.csv",
                                lambda r: r["role"] == "within_condition_descriptor_medoid")
    direction = read_rows(OLD / "supporting/data/vision_noise_direction_correlations.csv",
                          lambda r: r["condition"] == "noise" and float(r["fraction"]) == .01 and
                          float(r["gate"]) == .1 and r["window"] == "full" and float(r["lag"]) % 40 == 0)
    fig = plt.figure(figsize=(12.1, 10.9))
    gs = fig.add_gridspec(3, 3, height_ratios=(1.35, 1, 1.05), hspace=.64, wspace=.40)
    histories = []
    for j, (model, C) in enumerate((("wave_pinning", 20), ("holmes_model3", 60), ("goryachev", 20))):
        row, = select(representatives, model=model, circumference=C, fraction=.01)
        histories.append(history(fig.add_subplot(gs[0, j]), row, trace_cache, "A" + str(j + 1)))
    for j, model in enumerate(MODELS):
        ax = fig.add_subplot(gs[1, j])
        for C, color in COLORS.items():
            r = select(direction, model=model, circumference=C)
            x = np.array([float(a["lag"]) for a in r]) / 1000
            y, lo, hi = [np.array([float(a[k]) for a in r]) for k in ("direction_correlation", "ci_lower", "ci_upper")]
            ax.plot(x, y, color=color)
            ax.fill_between(x, lo, hi, color=color, alpha=.13, lw=0)
        ax.set_ylim(-.08, 1.04)
        ax.set_xlim(0, 8)
        ax.axhline(0, color=".8", lw=.6)
        ax.set_xlabel("Lag (10³ native units)")
        ax.set_ylabel("Activity-direction cosine")
        ax.set_title(f"B{j + 1}  {NAMES[model]}", loc="left", fontweight="bold", fontsize=10)
        if model == "goryachev":
            ax.text(.03, .09, "C20: mean input also removes polarity", fontsize=7.3, transform=ax.transAxes,
                    bbox=dict(facecolor="white", alpha=.9, edgecolor="none", pad=1))
    ax = fig.add_subplot(gs[2, 0])
    for model, offset, marker in (("wave_pinning", -1., "o"), ("holmes_model3", 1., "s")):
        r = select(endpoints, model=model, fraction=.01, gate="strict_resolved_site", lag=8000, metric="displacement_squared")
        for q in r:
            C = int(q["circumference"])
            x, y = C + offset, np.sqrt(q["value"])
            lo, hi = np.sqrt(q["ci_lower"]), np.sqrt(q["ci_upper"])
            ax.plot([x, x], [lo, hi], color=MODEL_COLORS[model], lw=1.1)
            ax.scatter([x], [y], color=MODEL_COLORS[model], marker=marker, s=27)
            ax.annotate(f"{q['eligible_pair_fraction']:.0%}", (x, y), xytext=(3, 7), textcoords="offset points", fontsize=7)
    ax.set_xticks((20, 40, 60))
    ax.set_xlabel("Circumference C")
    ax.set_ylabel("Site RMS displacement\n(native length units), lag 8000")
    ax.set_ylim(0, 7.8)
    ax.text(.04, .025, "Goryachev unavailable: 1 / 0 / 0 pairs", transform=ax.transAxes, fontsize=6.8)
    ax.set_title("C  Physical site displacement", loc="left", fontweight="bold", fontsize=10)
    ax.legend(handles=[Line2D([], [], color=MODEL_COLORS[m], marker=k, ls="", label=NAMES[m])
                       for m, k in (("wave_pinning", "o"), ("holmes_model3", "s"))], frameon=False, fontsize=7.5, loc="upper left")
    ax = fig.add_subplot(gs[2, 1])
    for model, marker in (("wave_pinning", "o"), ("holmes_model3", "s")):
        for metric, style in (("native_cosine", "-"), ("reference_cosine", "--")):
            r = select(endpoints, model=model, fraction=.01, gate="strict_resolved_site", lag=8000, metric=metric)
            x, y, lo, hi = [np.array([q[k] for q in r]) for k in ("circumference", "value", "ci_lower", "ci_upper")]
            ax.plot(x, y, style, marker=marker, color=MODEL_COLORS[model], ms=3.5, lw=1.1)
            ax.fill_between(x, lo, hi, color=MODEL_COLORS[model], alpha=.07, lw=0)
    ax.set_xticks((20, 40, 60))
    ax.set_xlabel("Circumference C")
    ax.set_ylabel("Eligible site-direction cosine\nlag 8000")
    ax.set_ylim(.05, 1.04)
    ax.set_title("D  Same endpoints, common geometry", loc="left", fontweight="bold", fontsize=9.5)
    ax.legend(handles=[Line2D([], [], color=".3", label="Native circumference"),
                       Line2D([], [], color=".3", ls="--", label="Remapped to C* = 60")], frameon=False, fontsize=7.5, loc="lower right")
    ax = fig.add_subplot(gs[2, 2])
    for model, marker in (("wave_pinning", "o"), ("holmes_model3", "s"), ("goryachev", "^")):
        for metric, style in (("multiregion_fraction", "-"), ("weak_R_fraction", "--")):
            r = select(obs, model=model, fraction=.01, metric=metric)
            x, y, lo, hi = [np.array([q[k] for q in r]) for k in ("circumference", "value", "ci_lower", "ci_upper")]
            ax.plot(x, y, style, marker=marker, color=MODEL_COLORS[model], ms=3.5, lw=1.1)
            ax.fill_between(x, lo, hi, color=MODEL_COLORS[model], alpha=.07, lw=0)
    ax.set_xticks((20, 40, 60))
    ax.set_ylim(-.015, .61)
    ax.set_xlabel("Circumference C")
    ax.set_ylabel("Fraction of observed sample times")
    ax.set_title("E  Local pattern and magnitude", loc="left", fontweight="bold", fontsize=10)
    ax.legend(handles=[Line2D([], [], color=".3", label="Multiple half-contrast regions"),
                       Line2D([], [], color=".3", ls="--", label="R < 0.1")], frameon=False, fontsize=7.2, loc="upper right")
    fig.suptitle("Figure 4 · Polarity memory combines local dynamics and boundary geometry", fontsize=14, y=.995)
    fig.text(.5, .956, "Higher declared input f = 0.01 · 12 independent seeds per condition · pointwise 95% whole-seed intervals", ha="center", fontsize=9)
    fig.legend(handles=[Line2D([], [], color=col, label=f"C = {C}") for C, col in COLORS.items()] +
                       [Line2D([], [], color=MODEL_COLORS[m], label=NAMES[m]) for m in MODELS],
               loc="lower center", bbox_to_anchor=(.5, -.006), ncol=6, frameon=False, fontsize=8.5)
    paths = save(fig, "figure4_stochastic_revision")
    return paths, histories


def plot_supplement(obs, endpoints, cdf, increments, cache, manifest):
    artifacts = []
    fig, axes = plt.subplots(3, 2, figsize=(10, 9), sharex=True, sharey=True)
    for i, model in enumerate(MODELS):
        for j, fraction in enumerate((.001, .01)):
            ax = axes[i, j]
            for C, color in COLORS.items():
                r = select(cdf, model=model, circumference=C, fraction=fraction, lag=5, gate=.1)
                x, y, lo, hi = [np.array([float(q[k]) for q in r]) for k in ("threshold_radians", "cdf", "ci_lower", "ci_upper")]
                ax.plot(x, 1-y, color=color)
                ax.fill_between(x, 1-hi, 1-lo, color=color, alpha=.12, lw=0)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(.00001, 1.)
            ax.set_ylim(1e-5, 1.05)
            ax.set_title(f"{NAMES[model]} · f = {fraction:g}", loc="left")
            if j == 0: ax.set_ylabel("Fraction above |Δ activity angle|")
            if i == 2: ax.set_xlabel("Absolute 5-unit angle increment (radians)")
    fig.suptitle("Descriptive increment tails; both endpoint R ≥ 0.1", fontsize=13)
    fig.text(.5, .015, "No abrupt-event classification · temporal origins remain clustered within seeds · exact maxima and gates .05/.2 in tables", ha="center", fontsize=8)
    fig.legend(handles=[Line2D([], [], color=col, label=f"C{C}") for C, col in COLORS.items()], loc="lower center", bbox_to_anchor=(.5, -.02), ncol=3, frameon=False)
    fig.tight_layout(rect=(0, .045, 1, .97))
    artifacts += save(fig, "figure4_increment_tails")
    fig, axes = plt.subplots(3, 3, figsize=(12, 9), sharex=True)
    for i, model in enumerate(MODELS):
        for j, C in enumerate(COLORS):
            ax = axes[i, j]
            for fraction, style in ((.001, "--"), (.01, "-")):
                jobs = [job for job in manifest["jobs"] if job["model"] == model and job["circumference"] == C and job["mean_fraction"] == fraction]
                noisy = [cache[job["slot"]]["R"] for job in jobs if job["kind"] == "noise"]
                t = cache[jobs[0]["slot"]]["time"][::8] / 1000
                value = np.stack(noisy)[:, ::8]
                w = weights(12, np.random.default_rng(BOOTSTRAP_SEED + i * 3 + j))
                lo, hi = np.percentile(w @ value, (2.5, 97.5), axis=0)
                ax.plot(t, value.mean(axis=0), color=COLORS[C], ls=style)
                ax.fill_between(t, lo, hi, color=COLORS[C], alpha=.12, lw=0)
                mean_job, = [job for job in jobs if job["kind"] == "mean"]
                ax.plot(t, cache[mean_job["slot"]]["R"][::8], color=".2", ls=style, lw=.9)
            basal, = [job for job in manifest["jobs"] if job["kind"] == "basal" and job["model"] == model and job["circumference"] == C]
            ax.plot(t, cache[basal["slot"]]["R"][::8], color=".65", lw=.7, ls=":")
            ax.axhline(.1, color=".65", ls=":", lw=.5)
            ax.set_ylim(0, 1.3)
            ax.set_title(f"{NAMES[model]} · C{C}", loc="left", fontsize=9)
            if j == 0: ax.set_ylabel("Relative activity magnitude R")
            if i == 2: ax.set_xlabel("Time (10³ native units)")
    fig.suptitle("Magnitude and mean/no-input controls distinguish forced-state changes", fontsize=13)
    fig.legend(handles=[Line2D([], [], color=".2", label="Mean input"), Line2D([], [], color=".65", ls=":", label="No input"),
                        Line2D([], [], color="#397bab", label="Noisy mean ± 95% seed interval"),
                        Line2D([], [], color=".4", ls="--", label="f = .001"), Line2D([], [], color=".4", label="f = .01")],
               loc="lower center", bbox_to_anchor=(.5, -.015), ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, .055, 1, .97))
    artifacts += save(fig, "figure4_magnitude_controls")
    # Both strengths of direction, native/site projection and coverage remain explicit.
    direction = read_rows(OLD / "supporting/data/vision_noise_direction_correlations.csv",
                          lambda r: r["condition"] == "noise" and float(r["gate"]) == .1 and r["window"] == "full" and float(r["lag"]) % 40 == 0)
    fig, axes = plt.subplots(3, 2, figsize=(10, 8), sharex=True, sharey=True)
    for i, model in enumerate(MODELS):
        for j, f in enumerate((.001, .01)):
            ax = axes[i, j]
            for C, color in COLORS.items():
                r = select(direction, model=model, circumference=C, fraction=f)
                x, y, lo, hi = [np.array([float(q[k]) for q in r]) for k in ("lag", "direction_correlation", "ci_lower", "ci_upper")]
                ax.plot(x / 1000, y, color=color)
                ax.fill_between(x / 1000, lo, hi, color=color, alpha=.13, lw=0)
            ax.set_ylim(-.08, 1.04)
            ax.set_title(f"{NAMES[model]} · f={f:g}", loc="left")
            if j == 0: ax.set_ylabel("Activity-direction cosine")
            if i == 2: ax.set_xlabel("Lag (10³ native units)")
    fig.tight_layout()
    artifacts += save(fig, "figure4_both_strengths_memory")
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.0), sharex=True)
    for i, f in enumerate((.001, .01)):
        for j, metric in enumerate(("displacement_squared", "native_cosine", "reference_cosine")):
            ax = axes[i, j]
            for model in MODELS:
                for gate, style in (("strict_resolved_site", "-"), ("one_halfcontrast_region", "--")):
                    r = select(endpoints, model=model, fraction=f, lag=8000, gate=gate, metric=metric)
                    allowed = [q for q in r if q["independent_seeds_with_pairs"] >= 6 and q["eligible_pair_fraction"] >= .01]
                    x = np.array([q["circumference"] for q in allowed])
                    y = np.array([q["value"] for q in allowed])
                    ax.plot(x, np.sqrt(y) if metric == "displacement_squared" else y, color=MODEL_COLORS[model], ls=style, marker="o", ms=3)
            metric_label = {"displacement_squared": "Site RMS displacement", "native_cosine": "Native site-direction cosine", "reference_cosine": "Site cosine remapped to C*=60"}[metric]
            ax.set_title(f"f={f:g} · {metric_label}", loc="left", fontsize=8)
            ax.set_ylabel("Native length units" if metric == "displacement_squared" else "Cosine")
            ax.set_xticks((20, 40, 60))
            ax.set_xlabel("C")
    fig.suptitle("Site-gate sensitivity at lag 8000; full eligibility and uncertainty remain in tables", fontsize=12)
    fig.legend(handles=[Line2D([], [], color=MODEL_COLORS[m], label=NAMES[m]) for m in MODELS] +
                       [Line2D([], [], color=".3", label="Original resolved-site gate"),
                        Line2D([], [], color=".3", ls="--", label="One half-contrast region, no peak veto")],
               loc="lower center", bbox_to_anchor=(.5, -.015), ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, .08, 1, .96))
    artifacts += save(fig, "figure4_site_gate_sensitivity")
    return artifacts


def validate_selected(histories):
    results = []
    for h in histories:
        p = Path(h["profile_path"])
        status = json.loads((p.parent / "status.json").read_text())
        assert sha(p) == status["output_sha256"]["primary.npz"]
        with np.load(p) as z:
            t, u, x = z["time"], z["primary"], z["x"]
        with np.load(p.parent / "trace.npz") as z:
            indices = (t / 5).astype(int)
            saved = z["M1_real"][indices] + 1j * z["M1_imag"][indices]
            m0 = z["M0"][indices]
        C = h["circumference"]
        actual = C / len(x) * (u @ np.exp(2j * np.pi * x / C))
        error = float(np.max(abs(actual - saved)))
        error0 = float(np.max(abs(C / len(x) * u.sum(axis=1) - m0)))
        np.testing.assert_allclose(actual, saved, rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(C / len(x) * u.sum(axis=1), m0, rtol=1e-11, atol=1e-11)
        results.append(dict(**h, primary_frames=len(t), moment_max_error=error, total_activity_max_error=error0))
    return results


def write_docs(obs, endpoints, increments, artifacts, histories, sources, selected_audit):
    def ob(model, C, metric):
        return select(obs, model=model, circumference=C, fraction=.01, metric=metric)[0]
    def ep(model, C, metric):
        return select(endpoints, model=model, circumference=C, fraction=.01, lag=8000, gate="strict_resolved_site", metric=metric)[0]
    highinc = [r for r in increments if r["condition"] == "noise" and r["lag"] == 5 and r["gate"] == .1]
    maxangle = max(r["max_absolute_angle"] for r in highinc)
    pmax = sum(r["diagnostic_above_pi_over_4"] for r in highinc)
    main_claims = [
        "At the higher declared input, larger Wave-pinning and Holmes rings retain more activity-direction memory; their physical site-displacement trends differ.",
        "Remapping the actual eligible displacements to a shared C*=60 markedly reduces the size contrast in site-angle memory. This is a geometric coordinate comparison conditional on the observed endpoint pairs.",
        "Passive local-pattern and magnitude readouts distinguish predominantly single-region wandering from multiregion organization and weak-resultant intervals; region counts do not establish functional competition.",
        "High-input Goryachev C20 is a different forced-state regime: its deterministic mean input also removes polarized asymmetry."
    ]
    caption = """# Figure 4. Polarity memory combines local dynamics and boundary geometry

The compact main sheet uses the higher of the two predeclared input strengths,
f = 0.01. Every quantitative noisy group retains 12 independent realizations;
bands and bars are pointwise 95% whole-seed percentile bootstrap intervals.
The lower strength, original full-vector memory and window diagnostics, all
18 medoids and mean/no-input controls remain linked supporting evidence.

**A1–A3.** Native active-field histories from the predeclared within-condition
descriptor medoids: Wave-pinning C20 seed 202610180018, Holmes C60 seed
202610180205 and Goryachev C20 seed 202610180086. Each image uses its own
zero-based native activity scale; identical relative boundary positions span different
physical lengths. Saved profiles are 40 native units apart. These histories
illustrate observed local organization and are not event-class exemplars.

**B1–B3.** Whole-window conditional activity-direction cosine at lag τ, using
both endpoint relative first-moment magnitudes R ≥ 0.1. The vector is the
first angular moment of native activity divided by its initial complex moment;
its angle is an activity axis, not a movement direction. Larger Wave-pinning
and Holmes rings decorrelate less over the measured lag horizon. The high-input
Goryachev C20 mean control also loses polarity, precluding a common-intact-state
wandering interpretation. Supplementary controls display this loss directly.

**C.** RMS shortest physical displacement of the dominant half-contrast region
at lag 8000, conditional on the original resolved-site endpoint gate. Numbers
are the eligible endpoint-pair percentages. Wave-pinning C20→C60 decreases
from 3.25 to 2.25 native length units; Holmes is 4.57 versus 5.54, with overlapping
seed intervals and coverage falling from 92% to 37%. No physical Goryachev
population curve is inferred from its sparse gate coverage.

**D.** The solid curves use each ring's native site-angle cosine. Dashed curves
remap every same observed eligible physical displacement to a common C*=60
before averaging. Thus both statistics use exactly the same pairs and weights;
no trajectories or kinetics are altered. C and D concern site centers, whereas
B concerns the global activity moment. Neither endpoint gating nor remapping
establishes uninterrupted site identity or a causal fraction of the memory effect.

**E.** Time fractions with multiple current half-contrast regions (solid) and
weak relative resultant R < 0.1 (dashed), using all 6401 scalar sample times per
seed. Line colours identify Wave-pinning (purple), Holmes (rose) and Goryachev
(teal). Larger Holmes rings spend more observed time in multiple-region
configurations without commensurate weak-resultant occupancy. Goryachev C20
shows both pattern multiplicity and weak-resultant intervals. Regions are
measured local activity domains, not assigned protrusions or competitors.

Native kinetics and local stochastic-input statistics are unchanged across
sizes within a model. The positive lognormal input has Gaussian local variance
one, physical spatial correlation parameter 1.5, temporal knot correlation time
50 and knot spacing 5; analytically variance-normalized linear interpolation
does not produce an exact continuous-time OU path. Fixed local statistics do
not equate total input or angular forcing across circumferences. Native model
units are not calibrated across models or to biological times.
"""
    (META / "CAPTION.md").write_text(caption)
    methods = r"""# Figure 4: methods, evidence and limits

This is a read-only reanalysis of the original 243 completed runs: 216 noisy
histories (3 models × 3 circumferences × 2 strengths × 12 independent seeds),
18 deterministic mean-input controls and nine no-input controls. All runs span
0–32000 native time units. No new simulation, burn-in selection, fitted
exponential time, event classifier or model parameter change is introduced.
Every native trace hash is checked against both producer status and the
previous audited seed table. Scalar cadence is 5; profile cadence is 40.

For circumference $C$, the native field and normalized activity moment are

$$
M_k(t)=\int_0^C u(s,t)e^{2\pi i k s/C}\,\mathrm ds,
\qquad P(t)=\frac{M_1(t)}{M_1(0)},\qquad R(t)=|P(t)|.\tag{1}
$$

The original direction table supplies B, retaining its frozen pooled-endpoint
estimator and 1000 whole-seed bootstrap draws. The full-vector raw correlation,
first/last-half windows and signed restricted areas remain in the unchanged
source report. New metrics use 2000 complete-seed bootstrap resamples, RNG seed
2026101904, and percentile 2.5–97.5 bounds. Time origins within a seed are never
independent resampling units; no confidence bands are supplied for one-run controls.

For the original resolved-site gate, both endpoints must have one current
half-contrast region, contrast ≥0.05 native units, at most one ≥0.05-prominence
peak and R≥0.1. A distinct sensitivity table removes only the peak-count veto.
Its less restrictive estimand does not repair or replace the primary gate.
The original full-profile independent descriptor audit remains authoritative;
all profiles displayed in A additionally reproduce their stored raw moments
within the numerical tolerance recorded in validation.json.

For shortest physical endpoint displacement, native site-angle memory and the
common-reference projection, use

$$
d_C(t,\tau)=\left[(s(t+\tau)-s(t)+C/2)\bmod C\right]-C/2,\tag{2}
$$

$$
K_C(\tau)=\left\langle\cos\left(\frac{2\pi d_C}{C}\right)\right\rangle_E,
\qquad
K_{60}(\tau)=\left\langle\cos\left(\frac{2\pi d_C}{60}\right)\right\rangle_E.\tag{3}
$$

The expectation pools the eligible endpoint pairs E across complete seeds.
Bootstrap resamples recompute the ratio of pooled seed numerator sums and
eligible-pair counts. Every native and remapped value uses identical E. The
reference circumference 60 is the largest tested, so every observed shortest
displacement maps to an angle at most π. This compares coordinates using
observed displacement distributions; it changes neither forcing nor dynamics.
It is not an additive causal attribution, a counterfactual evolved trajectory,
a diffusion estimate or evidence of uninterrupted site lineage. Original
shortest endpoint displacements already fold any travel beyond half the native
ring; remapping cannot recover winding or path length. Lags 5–8000, including
500, are exported so short-lag sensitivity to that folding is inspectable. Conditioning
on sparse single-region intervals can alter comparisons across sizes. CSVs
retain all sparse/undefined values, eligible-pair fractions, independent seeds
with pairs and bootstrap-defined fractions. Supplementary sensitivity curves
are shown only where at least six seeds contribute and coverage is ≥1%; this
display rule does not censor the tables. Goryachev's main physical curve remains
unavailable because the original peak gate is too restrictive.

Signed activity-direction increments use

$$
\Delta_\tau\phi(t)=\arg\left[P(t+\tau)P(t)^*\right],\qquad
\min(R(t),R(t+\tau))\ge r_0.\tag{4}
$$

Tables retain lags 5,40,100,500,2000,8000 and gates r₀=0.05,0.1,0.2. The
supplementary 5-unit survival curves average the pooled observed increment
distribution, with whole-seed uncertainty. Quantiles, maxima and the descriptive
count above π/4 are retained per seed. That threshold is a tail diagnostic,
not an abrupt-reorientation event definition. Small sampled increments do not
rule out subcadence fast changes; large wrapped increments alone would not
identify a physical site switch. The full native histories and R controls
must accompany interpretation.

Multiple-region and weak-R fractions use all 6401 sample times, with no chosen
adaptation exclusion. Half-contrast regions are defined relative to each
current profile's own minimum and range. Multiple regions need not be equal,
opposed, persistent, independent or functional competitors. Weak global
asymmetry can coexist with local activity. R is initial-moment normalized and
may exceed one; it is not a concentration, coherence ratio or migration speed.

The noise law is a declared model input; biological noise and migration are not
calibrated. These bounded histories do not establish stationary persistence,
universal enlargement effects or a force-to-motion relationship. The native
source campaign's fixed-forcing ODE and interpolation sensitivity checks remain
linked in its report; they are not long-horizon continuum convergence proofs.
"""
    (META / "METHODS_EVIDENCE_LIMITS.md").write_text(methods.replace(chr(92) * 2, chr(92)))
    geom = []
    for model in ("wave_pinning", "holmes_model3"):
        for C in (20, 60):
            a, b = ep(model, C, "native_cosine"), ep(model, C, "reference_cosine")
            geom.append(f"- {NAMES[model]} C{C}: native eligible site cosine {a['value']:.4f}; remapped C*=60 {b['value']:.4f} [{b['ci_lower']:.4f}, {b['ci_upper']:.4f}], endpoint coverage {a['eligible_pair_fraction']:.1%}.")
    result = "# Figure 4 chemical-ring interim evidence\n\nThe current four-model result, including Spring, is [RESULTS_WITH_SPRING.md](RESULTS_WITH_SPRING.md). This preserved interim note concerns the original three chemical models only.\n\n" + "\n\n".join(main_claims) + "\n\n" + "\n".join(geom)
    result += f"\n\nAt scalar cadence 5, the largest observed |activity-angle increment| among all chemical-model noisy runs with both endpoint R≥0.1 is {maxangle:.6f} rad ({np.degrees(maxangle):.3f}°). The π/4 tail diagnostic has {pmax} observations. This conditional sampled result supplies no evidence for a separate jump class and cannot exclude faster unresolved dynamics. Spring's cadence-.5 diagnostics include rare larger increments and are reported separately in the current result.\n"
    for model in MODELS:
        for C in (20, 60):
            a, b = ob(model, C, "multiregion_fraction"), ob(model, C, "weak_R_fraction")
            result += f"\n{NAMES[model]} C{C}: multiregion fraction {a['value']:.3%} [{a['ci_lower']:.3%}, {a['ci_upper']:.3%}]; weak-R fraction {b['value']:.3%} [{b['ci_lower']:.3%}, {b['ci_upper']:.3%}].\n"
    result += "\nGenerated analysis changes authoritative state by adding paired actual-displacement remapping, independent-seed pattern occupancy and increment distributions. All original simulations, tables, controls, representative selection and source outputs are preserved.\n"
    (OUT / "RESULTS.md").write_text(result)
    content = dict(title="Figure 4 · Polarity memory combines local dynamics and boundary geometry", role="stochastic dynamics",
                   claim="Polarity memory combines local dynamics and boundary geometry.",
                   main_asset=str(OUT / "figure4_stochastic_revision.png"),
                   assets=[dict(kind="main", path=str(OUT / "figure4_stochastic_revision.png"),
                                pdf=str(VECTORS / "figure4_stochastic_revision.pdf"), svg=str(VECTORS / "figure4_stochastic_revision.svg"),
                                caption="Native field histories, activity-direction memory, physical site displacement, common-circumference remapping and passive local-pattern readouts.")],
                   legend=[re.sub(r"\s+", " ", re.sub(r"[*#]", "", paragraph)).strip() for paragraph in caption.split("\n\n") if paragraph and not paragraph.startswith("# ")], narrative=main_claims,
                   measurements=["Conditional activity-direction memory", "Physical endpoint RMS displacement and pair coverage", "Native versus common-C60 site-angle projection on identical pairs", "Multiregion and weak-R time fractions", "Gated sampled increment distributions without event classes"],
                   limitations=["Activity angle does not assign cell movement.", "Endpoint pairs do not establish site identity.", "Projection is descriptive geometry, not a causal partition.", "Goryachev C20 high input changes mean-control state.", "No sampled increment class implies no absence of faster events.", "Finite observation window and prescribed noise law remain conditional."],
                   sources=[dict(label=label, path=str(path.relative_to(ROOT))) for label, path in
                            (("Frozen native campaign", RUN / "manifest.json"), ("Original results", OLD / "RESULTS.md"),
                             ("Declared noise design", OLD / "DESIGN.md"), ("Methods and limits", META / "METHODS_EVIDENCE_LIMITS.md"),
                             ("Native source hashes", DATA / "native_sources.csv"))],
                   supplementary=[dict(title=title, path=str(OUT / (stem + ".png")), caption=title) for title, stem in
                                (("Both declared strengths: memory", "figure4_both_strengths_memory"),
                                 ("Magnitude with all controls", "figure4_magnitude_controls"),
                                 ("Increment tails", "figure4_increment_tails"),
                                 ("Site-gate sensitivity", "figure4_site_gate_sensitivity"))],
                   inherited_support=[str(OLD / "vision_noise_selected_histories.png"), str(OLD / "vision_noise_correlations.png"), str(OLD / "vision_noise_restricted_integrals.png"), str(OLD / "vision_noise_lowest_resultant_diagnostic.png")])
    (META / "content.json").write_text(json.dumps(content, indent=2))
    validation = dict(generated_at=datetime.now(timezone.utc).isoformat(), analysis_only=True, new_native_integrations=0,
                      original_case_count=len(sources), independent_noisy_cases=sum(r["condition"] == "noise" for r in sources),
                      all_trace_hashes_match_producer_and_prior_audit=True, scalar_grid_verified="6401 points,0–32000,step5",
                      selfcheck=selfcheck(), representative_profile_checks=selected_audit,
                      uncertainty=dict(new_metrics_draws=DRAW_COUNT, seed=BOOTSTRAP_SEED, unit="whole independent seed",
                                       inherited_memory_draws=1000),
                      largest_eligible_5unit_angle_radians=maxangle, above_pi_over_4_observations=pmax,
                      projection_reference_circumference=C_REFERENCE,
                      source_script=str(Path(__file__).resolve()), source_script_sha256=sha(__file__),
                      outputs=[dict(path=str(p), sha256=sha(p)) for p in artifacts],
                      rendered_inspection="pending separate image inspection")
    (META / "validation.json").write_text(json.dumps(validation, indent=2))


def main():
    for p in (OUT, DATA, META, VECTORS):
        p.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": .65, "lines.linewidth": 1.4,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "svg.fonttype": "none", "pdf.fonttype": 42})
    manifest, sources, seeds, obs, endpoints, increments, cdf, cache = analyze()
    paths, histories = plot_main(obs, endpoints, cache)
    paths += plot_supplement(obs, endpoints, cdf, increments, cache, manifest)
    selected_audit = validate_selected(histories)
    write_docs(obs, endpoints, increments, paths, histories, sources, selected_audit)
    print(json.dumps(dict(output=str(OUT), original_cases=len(sources), artifacts=len(paths), new_integrations=0)))


if __name__ == "__main__":
    main()
