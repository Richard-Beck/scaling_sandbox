"""Microscopy-style rendering for deformable Taniguchi trajectories."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
from matplotlib import colormaps  # noqa: E402
from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
import matplotlib.patheffects as path_effects  # noqa: E402
import numpy as np


@dataclass(frozen=True)
class RenderScale:
    """Movie-wide upper limits used by the false-colour channels."""

    pip2_max: float
    pip3_max: float


def _positive_percentile(field: np.ndarray, phi: np.ndarray,
                         percentile: float = 99.5) -> float:
    """Estimate a stable movie-wide limit without copying the full movie."""
    if not 0 < percentile <= 100:
        raise ValueError("percentile must be in (0, 100]")
    frame_step = max(1, field.shape[0] // 100)
    pixel_step = max(1, field.shape[1] // 128)
    sampled_field = field[::frame_step, ::pixel_step, ::pixel_step]
    sampled_phi = phi[::frame_step, ::pixel_step, ::pixel_step]
    values = sampled_field[(sampled_phi >= 0.5) & np.isfinite(sampled_field)
                           & (sampled_field > 0)]
    if values.size == 0:
        return 1.0
    limit = float(np.percentile(values, percentile))
    return limit if limit > np.finfo(float).eps else 1.0


def determine_scale(U: np.ndarray, V: np.ndarray, phi: np.ndarray,
                    percentile: float = 99.5) -> RenderScale:
    """Determine robust intensity limits shared by every rendered frame."""
    return RenderScale(
        pip2_max=_positive_percentile(U, phi, percentile),
        pip3_max=_positive_percentile(V, phi, percentile),
    )


def microscopy_rgb(phi: np.ndarray, U: np.ndarray, V: np.ndarray, *,
                   mode: str, scale: RenderScale) -> np.ndarray:
    """Convert one state to an RGB microscopy-like image on black."""
    if mode not in {"pip3", "dual"}:
        raise ValueError("mode must be 'pip3' or 'dual'")
    support = np.clip((phi - 0.06) / 0.88, 0.0, 1.0)
    support = support * support * (3.0 - 2.0 * support)
    v = np.clip(V / scale.pip3_max, 0.0, 1.0)
    # Gamma brightens propagating fronts without changing the movie-wide LUT.
    v = np.sqrt(v)
    if mode == "pip3":
        rgb = colormaps["magma"](v)[..., :3]
        rgb *= support[..., None]
    else:
        u = np.sqrt(np.clip(U / scale.pip2_max, 0.0, 1.0))
        rgb = np.empty((*phi.shape, 3), dtype=float)
        # Additive cyan/magenta fluorescence composite: PIP2=(G,B), PIP3=(R,B).
        rgb[..., 0] = v
        rgb[..., 1] = 0.82 * u
        rgb[..., 2] = np.clip(0.68 * u + 0.72 * v, 0.0, 1.0)
        rgb *= support[..., None]
    return np.asarray(np.clip(rgb, 0.0, 1.0), dtype=np.float32)


class MicroscopyRenderer:
    """Reusable, axis-free canvas for movie frames and stills."""

    def __init__(self, x: np.ndarray, y: np.ndarray, *, condition: str,
                 mode: str, scale: RenderScale, pixels: int = 768,
                 dpi: int = 120):
        if pixels <= 0 or pixels % 2:
            raise ValueError("pixels must be a positive even number")
        self.x = x
        self.y = y
        self.condition = condition
        self.mode = mode
        self.scale = scale
        self.figure = Figure(figsize=(pixels / dpi, pixels / dpi), dpi=dpi,
                             facecolor="black")
        self.canvas = FigureCanvasAgg(self.figure)
        self.axes = self.figure.add_axes((0, 0, 1, 1), facecolor="black")
        self.axes.set_xlim(float(x[0]), float(x[-1]))
        self.axes.set_ylim(float(y[0]), float(y[-1]))
        self.axes.set_aspect("equal")
        self.axes.set_axis_off()
        self._pixels = pixels

    def render(self, phi: np.ndarray, U: np.ndarray, V: np.ndarray,
               time: float) -> np.ndarray:
        self.axes.clear()
        self.axes.set_facecolor("black")
        self.axes.set_xlim(float(self.x[0]), float(self.x[-1]))
        self.axes.set_ylim(float(self.y[0]), float(self.y[-1]))
        self.axes.set_aspect("equal")
        self.axes.set_axis_off()
        rgb = microscopy_rgb(phi, U, V, mode=self.mode, scale=self.scale)
        self.axes.imshow(rgb, origin="lower", interpolation="bilinear",
                         extent=(self.x[0], self.x[-1], self.y[0], self.y[-1]))
        contour = self.axes.contour(self.x, self.y, phi, levels=(0.5,),
                                    colors="#f4f4f4", linewidths=1.35,
                                    antialiased=True)
        contour.set_path_effects([
            path_effects.Stroke(linewidth=2.5, foreground=(0, 0, 0, 0.55)),
            path_effects.Normal(),
        ])
        label = "PIP3" if self.mode == "pip3" else "PIP2 / PIP3"
        self.axes.text(0.035, 0.955, f"Taniguchi {self.condition}  |  {label}",
                       transform=self.axes.transAxes, va="top", color="white",
                       fontsize=11)
        self.axes.text(0.965, 0.955, f"t = {time:05.2f}",
                       transform=self.axes.transAxes, ha="right", va="top",
                       color="white", fontsize=11, family="monospace")
        # A two-unit scale bar is short enough for the 20-unit published box.
        x0, x1 = float(self.x[-1]) - 3.0, float(self.x[-1]) - 1.0
        y0 = float(self.y[0]) + 1.0
        self.axes.plot((x0, x1), (y0, y0), color="white", lw=2.2,
                       solid_capstyle="butt")
        self.axes.text((x0 + x1) / 2, y0 + 0.28, "2 units", color="white",
                       fontsize=9, ha="center", va="bottom")
        self.canvas.draw()
        rgba = np.asarray(self.canvas.buffer_rgba())
        return np.ascontiguousarray(rgba[..., :3])


def video_indices(n_frames: int, fps: float, duration: float) -> np.ndarray:
    """Select trajectory frames for a fixed-duration time-lapse."""
    if n_frames < 1 or fps <= 0 or duration <= 0:
        raise ValueError("n_frames, fps, and duration must be positive")
    count = min(n_frames, max(2, int(round(fps * duration))))
    return np.unique(np.rint(np.linspace(0, n_frames - 1, count)).astype(int))


def still_indices(n_frames: int, count: int) -> np.ndarray:
    if n_frames < 1 or count < 1:
        raise ValueError("n_frames and count must be positive")
    return np.unique(np.rint(np.linspace(0, n_frames - 1,
                                        min(count, n_frames))).astype(int))


def representative_indices(V: np.ndarray, phi: np.ndarray, count: int) -> np.ndarray:
    """Choose endpoints and well-separated peaks of intracellular activity."""
    if V.shape != phi.shape or V.ndim != 3:
        raise ValueError("V and phi must be equally shaped frame stacks")
    n_frames = V.shape[0]
    if count <= 2 or n_frames <= 2:
        return still_indices(n_frames, count)
    inside = phi >= 0.5
    pixels = np.maximum(inside.sum(axis=(1, 2)), 1)
    masked = np.where(inside, np.maximum(V, 0), 0)
    peak = masked.max(axis=(1, 2))
    mean = masked.sum(axis=(1, 2)) / pixels
    peak_scale = peak.max()
    mean_scale = mean.max()
    if peak_scale <= np.finfo(float).eps:
        return still_indices(n_frames, count)
    activity = peak / peak_scale
    if mean_scale > np.finfo(float).eps:
        activity += mean / mean_scale
    candidates = np.flatnonzero(
        (activity[1:-1] >= activity[:-2])
        & (activity[1:-1] >= activity[2:])
        & (activity[1:-1] >= 0.01 * activity.max())
    ) + 1
    minimum_gap = max(1, n_frames // (4 * count))
    chosen = [0, n_frames - 1]
    for index in candidates[np.argsort(activity[candidates])[::-1]]:
        if all(abs(int(index) - old) >= minimum_gap for old in chosen):
            chosen.append(int(index))
            if len(chosen) == min(count, n_frames):
                break
    # Very quiet trajectories may not supply enough distinct activity peaks.
    for index in still_indices(n_frames, count):
        if len(chosen) == min(count, n_frames):
            break
        if int(index) not in chosen:
            chosen.append(int(index))
    return np.asarray(sorted(chosen), dtype=int)


def write_png(path: Path, rgb: np.ndarray) -> None:
    from PIL import Image

    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgb, mode="RGB").save(path, optimize=True)


def write_mp4(path: Path, frames, *, fps: float,
              ffmpeg: str | None = None, crf: int = 18) -> None:
    """Stream RGB frames directly to H.264, avoiding a temporary PNG set."""
    executable = ffmpeg or shutil.which("ffmpeg")
    if executable is None:
        raise RuntimeError(
            "ffmpeg was not found; on the HPC run "
            "'module load FFmpeg/7.1.2-GCCcore-14.3.0' first"
        )
    iterator = iter(frames)
    try:
        first = next(iterator)
    except StopIteration as exc:
        raise ValueError("cannot write an empty movie") from exc
    height, width, channels = first.shape
    if channels != 3:
        raise ValueError("movie frames must be RGB")
    path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        executable, "-y", "-loglevel", "error", "-f", "rawvideo",
        "-pix_fmt", "rgb24", "-s:v", f"{width}x{height}",
        "-r", str(fps), "-i", "-", "-an", "-c:v", "libx264",
        "-preset", "medium", "-crf", str(crf), "-pix_fmt", "yuv420p",
        "-movflags", "+faststart", str(path),
    ]
    process = subprocess.Popen(command, stdin=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    try:
        assert process.stdin is not None
        process.stdin.write(first.tobytes())
        for frame in iterator:
            if frame.shape != first.shape:
                raise ValueError("all movie frames must have the same shape")
            process.stdin.write(frame.tobytes())
        process.stdin.close()
        stderr = process.stderr.read().decode(errors="replace") if process.stderr else ""
        return_code = process.wait()
    except Exception:
        process.kill()
        process.wait()
        raise
    if return_code:
        raise RuntimeError(f"ffmpeg failed with status {return_code}: {stderr.strip()}")
