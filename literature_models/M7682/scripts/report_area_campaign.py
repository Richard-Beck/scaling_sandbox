"""Validate M7682 area campaign and make track, MSD, and GIF report."""

from __future__ import annotations

import csv
import datetime as dt
import html
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/area_scaling_80k"
OUT = BASE / "report"
NODE_LENGTH = 0.1
LATTICE_WIDTH = 250.0
TRACK_INTERVAL = 5
STOP_TIME = 80000
LABELS = {0.5: "0.5× area", 1.0: "1× area", 2.0: "2× area"}
COLORS = {0.5: "#177e89", 1.0: "#d17121", 2.0: "#684f9e"}
LAGS = np.unique(np.rint(np.geomspace(1, 4000, 105)).astype(int))
KEY_LAGS = (100, 1000, 5000, 10000, 20000)
LAGS = np.unique(np.r_[LAGS, np.array(KEY_LAGS) // TRACK_INTERVAL])


def get_track(run: Path) -> np.ndarray:
    root = ET.parse(run / "celltracks.xml").getroot()
    particles = root.findall(".//particle")
    if len(particles) != 1:
        raise ValueError(f"Expected one tracked cell in {run}, found {len(particles)}")
    detections = particles[0].findall("detection")
    times = np.array([int(d.get("t")) for d in detections])
    if not np.array_equal(times, np.arange(STOP_TIME // TRACK_INTERVAL + 1)):
        raise ValueError(f"Incomplete track in {run}: {len(times)} detections")
    wrapped = np.array([[float(d.get("x")), float(d.get("y"))] for d in detections])
    if not np.isfinite(wrapped).all():
        raise ValueError(f"Nonfinite track in {run}")
    delta = np.diff(wrapped, axis=0)
    delta -= LATTICE_WIDTH * np.round(delta / LATTICE_WIDTH)
    if np.max(np.abs(delta)) > LATTICE_WIDTH / 4:
        raise ValueError(f"Implausible displacement after unwrap in {run}")
    return np.vstack((np.zeros(2), np.cumsum(delta, axis=0))) * NODE_LENGTH


def msd(track: np.ndarray) -> np.ndarray:
    return np.array([
        np.mean(np.sum((track[lag:] - track[:-lag]) ** 2, axis=1))
        for lag in LAGS
    ])


def style_axes(ax):
    ax.grid(alpha=0.18)
    ax.spines[["top", "right"]].set_visible(False)


def make_plots(data: dict[float, list[dict]]) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), sharex=True, sharey=True, constrained_layout=True)
    max_extent = max(np.abs(run["track"]).max() for group in data.values() for run in group)
    limit = max(5, math.ceil(max_extent / 5) * 5)
    for ax, (size, group) in zip(axes, data.items()):
        for run in group:
            track = run["track"]
            ax.plot(track[:, 0], track[:, 1], color=COLORS[size], alpha=0.56, linewidth=0.8)
            ax.plot(track[0, 0], track[0, 1], ".", color="#111827", markersize=5)
        ax.set(title=LABELS[size], xlim=(-limit, limit), ylim=(-limit, limit), aspect="equal", xlabel="X displacement (model length units)")
        style_axes(ax)
    axes[0].set_ylabel("Y displacement (model length units)")
    fig.suptitle("Cell centroid tracks · 10 matched seeds per area · 80,000 model time units")
    fig.savefig(OUT / "tracks.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
    for size, group in data.items():
        values = np.vstack([run["msd"] for run in group])
        for row in values:
            ax.loglog(LAGS * TRACK_INTERVAL, row, color=COLORS[size], alpha=0.16, linewidth=0.8)
        mean = values.mean(axis=0)
        ax.loglog(LAGS * TRACK_INTERVAL, mean, color=COLORS[size], linewidth=2.5, label=LABELS[size])
    ax.set(xlabel="Lag (model time units)", ylabel="Time-origin-averaged MSD (model length units²)",
           title="Mean squared displacement across lags")
    style_axes(ax)
    ax.legend(frameon=False)
    fig.savefig(OUT / "msd.png", dpi=170)
    plt.close(fig)


def make_gif(run: Path, size: float, seed: int) -> str:
    # Native Morpheus four-panel snapshots, sampled every 2,000 time units.
    frames = []
    for time in range(0, STOP_TIME + 1, 2000):
        path = run / f"plot_{time:05d}.png"
        if not path.is_file():
            raise ValueError(f"Missing native snapshot: {path}")
        frame = Image.open(path).convert("RGB")
        canvas = Image.new("RGB", (frame.width, frame.height + 42), "white")
        canvas.paste(frame, (0, 42))
        draw = ImageDraw.Draw(canvas)
        draw.text((14, 14), f"M7682 · {LABELS[size]} · seed {seed} · time {time:,}", fill="black")
        frames.append(canvas.convert("P", palette=Image.Palette.ADAPTIVE, colors=128))
    name = f"representative_{size:g}x.gif"
    frames[0].save(OUT / name, save_all=True, append_images=frames[1:], duration=[130] * (len(frames) - 1) + [950], loop=0, optimize=True)
    return name


def main() -> None:
    manifest = json.loads((BASE / "manifest.json").read_text())
    OUT.mkdir(exist_ok=True)
    data = {size: [] for size in LABELS}
    per_run = []
    for task in manifest["runs"]:
        run = BASE / task["relative_dir"]
        if not (run / "COMPLETE").is_file():
            raise RuntimeError(f"Run incomplete: {run}")
        if not (run / "plot_80000.png").is_file():
            raise RuntimeError(f"Final plot missing: {run}")
        track = get_track(run)
        values = msd(track)
        if not np.isfinite(values).all():
            raise ValueError(f"Nonfinite MSD in {run}")
        size = task["area_multiplier"]
        entry = {**task, "track": track, "msd": values, "run": run}
        data[size].append(entry)
        per_run.append({"area_multiplier": size, "seed": task["seed"],
                        "net_displacement": float(np.linalg.norm(track[-1])),
                        "path_length": float(np.linalg.norm(np.diff(track, axis=0), axis=1).sum()),
                        **{f"msd_lag_{lag}": float(values[np.where(LAGS * TRACK_INTERVAL == lag)[0][0]]) for lag in KEY_LAGS}})
    for size, group in data.items():
        if len(group) != 10 or [x["seed"] for x in group] != list(range(19, 29)):
            raise ValueError(f"Unexpected seed set for {size}")
    with (OUT / "per_run.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=per_run[0].keys())
        writer.writeheader()
        writer.writerows(per_run)
    with (OUT / "msd_by_lag.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["area_multiplier", "seed", "lag_model_time", "msd_model_length_squared"])
        for size, group in data.items():
            for run in group:
                writer.writerows((size, run["seed"], int(lag * TRACK_INTERVAL), float(value)) for lag, value in zip(LAGS, run["msd"]))
    make_plots(data)
    representative = {}
    reference_lag = int(np.where(LAGS * TRACK_INTERVAL == 5000)[0][0])
    for size, group in data.items():
        median = np.median([run["msd"][reference_lag] for run in group])
        chosen = min(group, key=lambda run: abs(run["msd"][reference_lag] - median))
        representative[size] = {"seed": chosen["seed"], "gif": make_gif(chosen["run"], size, chosen["seed"])}
    generated = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    rows = []
    for lag in KEY_LAGS:
        idx = int(np.where(LAGS * TRACK_INTERVAL == lag)[0][0])
        means = {size: np.mean([run["msd"][idx] for run in group]) for size, group in data.items()}
        rows.append(f"<tr><td>{lag:,}</td><td>{means[0.5]:.3g}</td><td>{means[1.0]:.3g}</td><td>{means[2.0]:.3g}</td><td>{means[0.5]/means[1.0]:.2f}</td><td>{means[2.0]/means[1.0]:.2f}</td></tr>")
    gifs = "".join(f'<figure><img src="{html.escape(representative[size]["gif"])}" alt="Native model animation for {LABELS[size]}"><figcaption>{LABELS[size]} · seed {representative[size]["seed"]}</figcaption></figure>' for size in LABELS)
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>M7682 area scaling · 80,000 time units</title><style>
body{{font:16px/1.55 system-ui,sans-serif;color:#202937;background:#f7f8fa;margin:0}}main{{max-width:1160px;margin:auto;padding:36px 24px 70px}}h1,h2{{line-height:1.2}}.sub{{color:#5b6573}}.card{{background:white;border:1px solid #dde2e8;border-radius:12px;padding:22px;margin:22px 0}}img{{max-width:100%;height:auto}}.wide{{width:100%}}.gifs{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}}figure{{margin:0}}figcaption{{text-align:center;color:#4b5563}}table{{border-collapse:collapse;width:100%}}th,td{{padding:8px 11px;border-bottom:1px solid #e1e4e8;text-align:right}}th:first-child,td:first-child{{text-align:left}}a{{color:#126384}}@media(max-width:800px){{.gifs{{grid-template-columns:1fr}}}}
</style></head><body><main><h1>M7682 · area scaling experiment</h1><p class="sub">30 completed simulations · 10 matched seeds per area · 80,000 model time units · generated {generated}</p>
<div class="card"><h2>Design</h2><p>Areas are 0.5×, 1×, and 2× the source model target of 1,200 lattice sites. Initial radius and target perimeter scale with the square root of the area multiplier; area constraint strength scales inversely with it. Local biochemical initial values, kinetics, diffusion, membrane resolution, lattice, and CPM timing remain fixed. Seeds 19–28 are matched across groups. One cell is simulated per run.</p><p>The source file <code>model/model.xml</code> remains unchanged. “80,000” refers to Morpheus StopTime; CPM MCSDuration is 0.5. Coordinates use the model node length of 0.1. Centroid positions were unwrapped across the periodic 250 × 250 lattice.</p></div>
<div class="card"><h2>Centroid tracks</h2><img class="wide" src="tracks.png" alt="Thirty centroid tracks in three area groups"><p>All tracks start at the origin and use the same physical axis limits. Each line is one seed; this plot does not imply cells share the same stochastic history across sizes.</p></div>
<div class="card"><h2>Mean squared displacement</h2><img class="wide" src="msd.png" alt="MSD curves over lags from 5 to 20,000 model time units"><p>Thin lines show each seed; thick lines show the equal-seed mean. MSD at each lag averages squared displacement over all available time origins in a track. Long lags have fewer time origins. Ratios below compare group means and are descriptive, not confidence intervals.</p><table><thead><tr><th>Lag</th><th>0.5× MSD</th><th>1× MSD</th><th>2× MSD</th><th>0.5× / 1×</th><th>2× / 1×</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
<div class="card"><h2>Representative native simulations</h2><div class="gifs">{gifs}</div><p>Each GIF samples the native Morpheus four-panel output every 2,000 model time units. The representative seed is the run nearest its group's median MSD at lag 5,000. Native color limits are fixed by the source model and can saturate.</p></div>
<div class="card"><h2>Data and interpretation</h2><p><a href="per_run.csv">Per-run summaries</a> · <a href="msd_by_lag.csv">Seed-level MSD by lag</a> · <a href="../manifest.json">Run manifest and exact scaling</a>. This is an area-scaling model experiment; it does not include gene dosage, biochemical remodeling, or measured WGD cell geometry. The native model's dynamics and parameterization are retained.</p><p>Report data generated {generated}. Reload this page manually after rerunning the analysis.</p></div></main></body></html>"""
    (OUT / "index.html").write_text(page)
    summary = {"generated_utc": generated, "runs": len(per_run), "representative": representative,
               "lag_range": [int(LAGS[0] * TRACK_INTERVAL), int(LAGS[-1] * TRACK_INTERVAL)]}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
