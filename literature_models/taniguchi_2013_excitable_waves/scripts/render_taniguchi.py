#!/usr/bin/env python3
"""Render Taniguchi 4F--4I trajectories as microscopy-style time lapses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from taniguchi.render import (  # noqa: E402
    MicroscopyRenderer,
    RenderScale,
    determine_scale,
    representative_indices,
    video_indices,
    write_mp4,
    write_png,
)


def _condition(data: np.lib.npyio.NpzFile, path: Path) -> str:
    if "summary_json" in data.files:
        summary = json.loads(str(data["summary_json"]))
        if "parameter_set" in summary:
            return str(summary["parameter_set"]).upper()
    for value in ("4F", "4G", "4H", "4I"):
        if value.lower() in path.stem.lower():
            return value
    raise ValueError(f"cannot infer 4F--4I condition from {path}")


def render_trajectory(path: Path, output_dir: Path, *, modes: tuple[str, ...],
                      fps: float, duration: float, still_count: int,
                      pixels: int, percentile: float,
                      pip2_max: float | None, pip3_max: float | None,
                      ffmpeg: str | None) -> list[Path]:
    if pip2_max is not None and pip2_max <= 0:
        raise ValueError("pip2_max must be positive")
    if pip3_max is not None and pip3_max <= 0:
        raise ValueError("pip3_max must be positive")
    created: list[Path] = []
    with np.load(path, allow_pickle=False) as data:
        required = {"times", "phi", "U", "V", "x", "y"}
        missing = required.difference(data.files)
        if missing:
            raise ValueError(f"{path} is missing arrays: {sorted(missing)}")
        times, phi = data["times"], data["phi"]
        U, V, x, y = data["U"], data["V"], data["x"], data["y"]
        if not (len(times) == len(phi) == len(U) == len(V)):
            raise ValueError(f"inconsistent trajectory lengths in {path}")
        condition = _condition(data, path)
        auto_scale = determine_scale(U, V, phi, percentile)
        scale = RenderScale(
            auto_scale.pip2_max if pip2_max is None else pip2_max,
            auto_scale.pip3_max if pip3_max is None else pip3_max,
        )
        movie_frames = video_indices(len(times), fps, duration)
        representative = representative_indices(V, phi, still_count)

        for mode in modes:
            tag = "pip3" if mode == "pip3" else "pip2_pip3"
            renderer = MicroscopyRenderer(x, y, condition=condition, mode=mode,
                                          scale=scale, pixels=pixels)
            movie_path = output_dir / f"phase_{condition}_{tag}.mp4"
            frames = (renderer.render(phi[i], U[i], V[i], float(times[i]))
                      for i in movie_frames)
            write_mp4(movie_path, frames, fps=fps, ffmpeg=ffmpeg)
            created.append(movie_path)

            still_dir = output_dir / f"phase_{condition}_{tag}_stills"
            for i in representative:
                still_path = still_dir / f"phase_{condition}_{tag}_t{times[i]:06.2f}.png"
                write_png(still_path, renderer.render(phi[i], U[i], V[i],
                                                      float(times[i])))
                created.append(still_path)
        metadata = {
            "source": str(path), "condition": condition, "modes": list(modes),
            "fps": fps, "duration_seconds": duration,
            "movie_frame_count": int(len(movie_frames)),
            "still_count": int(len(representative)), "pixels": pixels,
            "still_times": [float(times[i]) for i in representative],
            "intensity_percentile": percentile,
            "pip2_max": scale.pip2_max, "pip3_max": scale.pip3_max,
            "outline_level": 0.5,
        }
        metadata_path = output_dir / f"phase_{condition}_render.json"
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
        created.append(metadata_path)
    return created


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="+", help="phase_4F.npz etc.")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "reports/output/taniguchi/rendered")
    parser.add_argument("--mode", choices=("pip3", "dual", "both"),
                        default="both")
    parser.add_argument("--fps", type=float, default=20.0)
    parser.add_argument("--duration", type=float, default=20.0,
                        help="target movie duration in seconds")
    parser.add_argument("--still-count", type=int, default=5)
    parser.add_argument("--pixels", type=int, default=768)
    parser.add_argument("--percentile", type=float, default=99.5)
    parser.add_argument("--pip2-max", type=float)
    parser.add_argument("--pip3-max", type=float)
    parser.add_argument("--ffmpeg", help="path to ffmpeg executable")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    modes = ("pip3", "dual") if args.mode == "both" else (args.mode,)
    for source in args.inputs:
        created = render_trajectory(
            source, args.output_dir, modes=modes, fps=args.fps,
            duration=args.duration, still_count=args.still_count,
            pixels=args.pixels, percentile=args.percentile,
            pip2_max=args.pip2_max, pip3_max=args.pip3_max,
            ffmpeg=args.ffmpeg,
        )
        for path in created:
            print(path)


if __name__ == "__main__":
    main()
