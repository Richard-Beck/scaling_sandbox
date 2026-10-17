"""Assemble native simulation plots without changing the plotted data."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "results/default"
OUT = ROOT / "results/snapshots"
OUT.mkdir(exist_ok=True)
times = list(range(0, 5001, 1000))
paths = [RUN / f"plot_{t:05d}.png" for t in times]
images = [np.asarray(Image.open(p).convert("RGB")) for p in paths]
assert all(im.shape == (636, 720, 3) for im in images)
assert "=== Simulation finished ===" in (RUN / "simulation.log").read_text()
data = np.loadtxt(RUN / "cells_displacement.dat")
assert data[-1, 0] == 5000 and np.isfinite(data).all()

fig, axes = plt.subplots(2, 3, figsize=(18, 12.8))
for ax, t, im in zip(axes.flat, times, images):
    ax.imshow(im)
    ax.set_title(f"Time = {t:,}", fontsize=15)
    ax.axis("off")
fig.suptitle("M7682 · default conditions · native Morpheus snapshots", fontsize=21)
fig.text(.5, .015, "Unmodified model.xml | LK = 5.77 | seed = 19 | each snapshot: cell shape, R, B, kb", ha="center", fontsize=13)
fig.tight_layout(rect=(0, .035, 1, .965))
fig.savefig(OUT / "default_snapshots.png", dpi=160)
fig.savefig(OUT / "default_snapshots.pdf")
plt.close(fig)

# Native Cell ID plot interior. Use one shared crop across all times so
# translation, orientation, and size remain directly comparable.
masks = [im[24:299, 9:285] for im in images]
union = np.any(np.stack([np.all(im > 220, axis=2) for im in masks]), axis=0)
yy, xx = np.where(union)
x0, x1 = max(0, xx.min()-18), min(276, xx.max()+19)
y0, y1 = max(0, yy.min()-18), min(275, yy.max()+19)
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
for ax, t, im in zip(axes.flat, times, masks):
    ax.imshow(im[y0:y1, x0:x1], interpolation="nearest")
    ax.set_title(f"Time = {t:,}")
    ax.axis("off")
fig.suptitle("M7682 · cell shape and movement", fontsize=20)
fig.text(.5, .018, "White: cell | Black: medium | Same fixed window and scale in every panel", ha="center")
fig.subplots_adjust(left=.03, right=.97, bottom=.09, top=.87, wspace=.07, hspace=.22)
fig.savefig(OUT / "cell_shapes.png", dpi=180)
plt.close(fig)

model = ROOT / "model/model.xml"
xml = ET.parse(model).getroot()
manifest = {
    "source": "https://morpheus.gitlab.io/model/m7682/",
    "model_url": "https://morpheus.gitlab.io/media/model/m7682/model.xml",
    "model_sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
    "model_modified": False,
    "morpheus_version": "2.3.9 (revision 25.02.2025)",
    "binary_url": "https://gitlab.com/morpheus.lab/morpheus/-/package_files/178666119/download",
    "binary_sha256": hashlib.sha256((ROOT / "software/morpheus-2.3.9").read_bytes()).hexdigest(),
    "gnuplot_version": "5.4.2",
    "threads": 1,
    "slurm_job": 20510538,
    "random_seed": int(xml.find("Time/RandomSeed").get("value")),
    "stop_time": float(xml.find("Time/StopTime").get("value")),
    "snapshot_times": times,
    "final_displacement_lattice_units": float(data[-1, 1]),
    "validation": "Simulation finished; six readable PNGs; displacement output finite and reaches time 5000",
}
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(manifest, indent=2))
