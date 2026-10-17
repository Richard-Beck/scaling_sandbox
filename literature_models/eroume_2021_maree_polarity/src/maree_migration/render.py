"""Cell-movie, montage, and diagnostic rendering."""

from __future__ import annotations

from pathlib import Path
import math
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
import numpy as np
from PIL import Image

from .metrics import trajectory_metrics
from .model import SimulationResult


def _limit(stack: np.ndarray, phi: np.ndarray) -> float:
    values = stack[(phi >= 0.5) & np.isfinite(stack)]
    return max(float(np.percentile(values, 99.5)), 1e-12) if values.size else 1.0


class CellRenderer:
    def __init__(self, result: SimulationResult, *, pixels: int = 640,
                 diagnostic: bool = False, trajectory: bool = False):
        self.result = result
        self.diagnostic = diagnostic
        self.show_trajectory = trajectory
        self.metrics = trajectory_metrics(result)
        self.rac_max = _limit(result.concentrations[:, 3], result.phi)
        self.rho_max = _limit(result.concentrations[:, 5], result.phi)
        self.p3_max = _limit(result.concentrations[:, 8], result.phi)
        width = pixels * (3 if diagnostic else 1)
        self.figure = Figure(figsize=(width / 100, pixels / 100), dpi=100,
                             facecolor="black")
        self.canvas = FigureCanvasAgg(self.figure)

    def frame(self, index: int) -> np.ndarray:
        self.figure.clear()
        n = 3 if self.diagnostic else 1
        axes = self.figure.subplots(1, n, squeeze=False)[0]
        fields = ((3, "Rac", self.rac_max), (5, "Rho", self.rho_max),
                  (8, "PIP3", self.p3_max)) if self.diagnostic else ((3, "Rac", self.rac_max),)
        phi = self.result.phi[index]
        extent = (self.result.x[0], self.result.x[-1],
                  self.result.y[0], self.result.y[-1])
        for ax, (species, label, maximum) in zip(axes, fields):
            ax.set_facecolor("black")
            support = np.clip(phi, 0.0, 1.0)
            image = np.sqrt(np.clip(self.result.concentrations[index, species] / maximum,
                                    0.0, 1.0)) * support
            ax.imshow(image, origin="lower", extent=extent, cmap="magma",
                      vmin=0.0, vmax=1.0, interpolation="bilinear")
            ax.contour(self.result.x, self.result.y, phi, (0.5,), colors="white",
                       linewidths=1.1)
            if self.show_trajectory:
                ax.plot(self.metrics[:index + 1, 1], self.metrics[:index + 1, 2],
                        color="#5ee7ff", lw=1.0, alpha=0.85)
            angle = (self.result.config.stimulus.angle_at(float(self.result.times[index]))
                     if self.result.config.stimulus else None)
            if angle is not None:
                radians = math.radians(angle)
                ax.annotate("", xy=(0.16 + 0.08 * math.cos(radians),
                                     0.88 + 0.08 * math.sin(radians)),
                            xytext=(0.16, 0.88), xycoords="axes fraction",
                            arrowprops={"arrowstyle": "->", "color": "white", "lw": 2})
            ax.text(0.03, 0.97, label, color="white", transform=ax.transAxes,
                    ha="left", va="top")
            ax.text(0.97, 0.97, f"t={self.result.times[index]:.1f} s", color="white",
                    transform=ax.transAxes, ha="right", va="top", family="monospace")
            ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:]); ax.set_aspect("equal")
            ax.set_axis_off()
        self.figure.subplots_adjust(0, 0, 1, 1, wspace=0.01)
        self.canvas.draw()
        return np.ascontiguousarray(np.asarray(self.canvas.buffer_rgba())[..., :3])


def _write_mp4(path: Path, frames, fps: float) -> None:
    executable = shutil.which("ffmpeg")
    if executable is None:
        raise RuntimeError("ffmpeg not found; load the HPC FFmpeg module before --render")
    iterator = iter(frames)
    first = next(iterator)
    height, width = first.shape[:2]
    path.parent.mkdir(parents=True, exist_ok=True)
    command = [executable, "-y", "-loglevel", "error", "-f", "rawvideo",
               "-pix_fmt", "rgb24", "-s:v", f"{width}x{height}", "-r", str(fps),
               "-i", "-", "-an", "-c:v", "libx264", "-crf", "19",
               "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(path)]
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        assert process.stdin is not None
        process.stdin.write(first.tobytes())
        for frame in iterator:
            process.stdin.write(frame.tobytes())
        process.stdin.close()
        stderr = process.stderr.read().decode() if process.stderr else ""
        status = process.wait()
    except Exception:
        process.kill(); process.wait()
        raise
    if status:
        raise RuntimeError(f"ffmpeg failed: {stderr.strip()}")


def render_outputs(result: SimulationResult, output_dir: Path, *, fps: float = 15.0,
                   montage_count: int = 7, pixels: int = 640,
                   trajectory: bool = False, diagnostic: bool = True) -> list[Path]:
    """Write a main MP4, 6--8 frame montage, PNG sequence, and diagnostic MP4."""
    output_dir.mkdir(parents=True, exist_ok=True)
    indices = np.unique(np.rint(np.linspace(0, len(result.times) - 1,
                                           min(montage_count, len(result.times)))).astype(int))
    renderer = CellRenderer(result, pixels=pixels, trajectory=trajectory)
    movie = output_dir / "cell_rac.mp4"
    _write_mp4(movie, (renderer.frame(i) for i in range(len(result.times))), fps)
    still_dir = output_dir / "frames"
    still_dir.mkdir(exist_ok=True)
    selected_frames = []
    for i in indices:
        selected_frames.append(renderer.frame(int(i)))
        Image.fromarray(selected_frames[-1]).save(still_dir / f"frame_{i:04d}.png")
    thumb_width = max(1, pixels // 3)
    thumbs = [Image.fromarray(frame).resize((thumb_width, thumb_width))
              for frame in selected_frames]
    montage = Image.new("RGB", (thumb_width * len(thumbs), thumb_width), "black")
    for column, thumb in enumerate(thumbs):
        montage.paste(thumb, (column * thumb_width, 0))
    montage_path = output_dir / "montage.png"
    montage.save(montage_path)
    created = [movie, montage_path, *sorted(still_dir.glob("*.png"))]
    if diagnostic:
        diagnostic_renderer = CellRenderer(result, pixels=max(240, pixels // 2),
                                           diagnostic=True, trajectory=trajectory)
        diagnostic_path = output_dir / "diagnostic_rac_rho_pip3.mp4"
        _write_mp4(diagnostic_path,
                   (diagnostic_renderer.frame(i) for i in range(len(result.times))), fps)
        created.append(diagnostic_path)
    return created
