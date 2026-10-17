"""Two-dimensional ungated polarity/velocity SDE; reproducible illustrative paths."""
import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


@dataclass(frozen=True)
class Parameters:
    gp: float = 3.0
    gv: float = 3.0
    r: float = 1.0


@dataclass(frozen=True)
class DimensionalParameters:
    tau_p: float = 1.0
    tau_v: float = 10.0
    alpha: float = 1.0
    beta: float = 1.0
    vs: float = 1.0
    ps: float = 1.0
    a: float = 0.0
    K: float = 0.06
    n: float = 2.0

    def __post_init__(self):
        if not all(np.isfinite(value) for value in asdict(self).values()):
            raise ValueError('Dimensional parameters must be finite')
        if self.a < 0 or self.K <= 0 or self.n < 1:
            raise ValueError('Require a>=0, K>0 and n>=1')
        if min(self.tau_p, self.tau_v, self.vs, self.ps) <= 0:
            raise ValueError('tau_p, tau_v, vs, ps must be positive')


def drift(p, v, params):
    if isinstance(params, DimensionalParameters):
        dp = -p / params.tau_p + params.beta * v / (1 + np.linalg.norm(v, axis=-1, keepdims=True) / params.vs)
        magnitude = np.linalg.norm(p, axis=-1, keepdims=True)
        dp += params.a * p * magnitude**(params.n-1) / (params.K**params.n + magnitude**params.n)
        dv = -v / params.tau_v + params.alpha * p / (1 + np.linalg.norm(p, axis=-1, keepdims=True) / params.ps)
        return dp, dv
    dp = -p + params.gp * v / (1 + np.linalg.norm(v, axis=-1, keepdims=True))
    dv = params.r * (-v + params.gv * p / (1 + np.linalg.norm(p, axis=-1, keepdims=True)))
    return dp, dv


def simulate(params, p0, v0, sigma_p, sigma_v, duration, dt, seed):
    """Euler-Maruyama, independent Wiener components; left-point position update."""
    if duration <= 0 or dt <= 0 or min(sigma_p, sigma_v) < 0:
        raise ValueError("Require duration, dt > 0 and noise strengths >= 0")
    steps = round(duration / dt)
    if steps < 1 or not np.isclose(steps * dt, duration, rtol=1e-10, atol=1e-12):
        raise ValueError("duration must be an integer multiple of dt")
    p0, v0 = np.asarray(p0, dtype=float), np.asarray(v0, dtype=float)
    if p0.shape != v0.shape or p0.ndim != 2 or p0.shape[1] != 2:
        raise ValueError("p0 and v0 must have matching (n_tracks, 2) shapes")
    rng = np.random.default_rng(seed)
    p = np.empty((steps + 1,) + p0.shape)
    v = np.empty_like(p)
    x = np.zeros_like(p)
    p[0], v[0] = p0, v0
    for i in range(steps):
        dp, dv = drift(p[i], v[i], params)
        p[i + 1] = p[i] + dt * dp + sigma_p * np.sqrt(dt) * rng.standard_normal(p0.shape)
        v[i + 1] = v[i] + dt * dv + sigma_v * np.sqrt(dt) * rng.standard_normal(v0.shape)
        x[i + 1] = x[i] + dt * v[i]
    if not all(np.isfinite(a).all() for a in (p, v, x)):
        raise FloatingPointError("Nonfinite trajectory; reduce timestep/check parameters")
    return np.arange(steps + 1) * dt, p, v, x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimensional', action='store_true', help='Integrate the dimensional equations directly')
    for name, default in asdict(Parameters()).items():
        parser.add_argument('--' + name.lower(), type=float, default=default, dest=name)
    for name, default in asdict(DimensionalParameters()).items():
        flag = '--' + name.replace('_', '-')
        parser.add_argument(flag, type=float, default=default, dest=name)
    parser.add_argument('--noise', nargs='+', type=float, default=[0.15, 0.45, 0.9], help='Equal sigma_p and sigma_v in the selected units')
    parser.add_argument('--duration', type=float, default=100.0)
    parser.add_argument('--dt', type=float, default=0.005)
    parser.add_argument('--sample-interval', type=float, default=0.05)
    parser.add_argument('--seed', type=int, default=20260924)
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'results_ungated')
    args = parser.parse_args()
    params = Parameters(**{k: getattr(args, k) for k in asdict(Parameters())})
    if args.dimensional:
        params = DimensionalParameters(**{k: getattr(args, k) for k in asdict(DimensionalParameters())})
    stride = round(args.sample_interval / args.dt)
    samples = round(args.duration / args.sample_interval)
    if stride < 1 or not np.isclose(stride * args.dt, args.sample_interval) or samples < 1 or not np.isclose(samples * args.sample_interval, args.duration):
        parser.error('sample-interval must be a positive multiple of dt and divide duration')
    args.out.mkdir(parents=True, exist_ok=True)
    angles = np.arange(3) * 2 * np.pi / 3
    initial = np.vstack([np.zeros((3, 2)), 2 * np.column_stack([np.cos(angles), np.sin(angles)])])
    title = ', '.join(f'{key} = {value:g}' for key, value in asdict(params).items())
    time_label = 'Time t' if args.dimensional else 'Nondimensional time T'
    p_label, v_label, x_label = ('p', 'v', 'x') if args.dimensional else ('P', 'V', 'X')
    overview, axes = plt.subplots(1, len(args.noise), figsize=(6 * len(args.noise), 5), squeeze=False)
    diagnostics, dax = plt.subplots(2, len(args.noise), figsize=(6 * len(args.noise), 7), squeeze=False)
    scenarios = []
    for j, sigma in enumerate(args.noise):
        seed = args.seed + j
        t, p, v, x = simulate(params, initial, initial, sigma, sigma, args.duration, args.dt, seed)
        keep = np.unique(np.r_[np.arange(0, len(t), stride), len(t) - 1])
        t, p, v, x = t[keep], p[keep], v[keep], x[keep]
        stem = f'noise_{sigma:g}'.replace('.', 'p')
        np.savez_compressed(args.out / f'{stem}.npz', **dict(zip(('t', 'p', 'v', 'x') if args.dimensional else ('T', 'P', 'V', 'X'), (t, p, v, x))))
        fig, grid = plt.subplots(2, 3, figsize=(12, 8), layout='constrained')
        for track, ax in enumerate(grid.flat):
            points = x[:, track]
            segments = np.stack([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap='viridis', norm=plt.Normalize(0, args.duration), linewidth=1.5)
            lc.set_array(t[:-1])
            ax.add_collection(lc)
            ax.autoscale()
            ax.scatter(*points[0], c='black', marker='o', s=28, label='Start', zorder=3)
            ax.scatter(*points[-1], c='crimson', marker='X', s=40, label='End', zorder=3)
            ax.set(title=f'Track {track + 1}: ' + ('rest start' if track < 3 else 'aligned start'), xlabel=x_label+'1', ylabel=x_label+'2')
            ax.set_aspect('equal', adjustable='datalim')
            ax.grid(alpha=0.2)
            axes[0, j].plot(*points.T, lw=1.2, label=f'{track + 1}: ' + ('rest' if track < 3 else 'aligned'))
            dax[0, j].plot(t, np.linalg.norm(v[:, track], axis=-1), lw=0.7, alpha=0.8)
            dax[1, j].plot(t, np.linalg.norm(p[:, track], axis=-1), lw=0.7, alpha=0.8)
        grid[0, 0].legend(fontsize=8)
        fig.colorbar(lc, ax=grid.ravel().tolist(), label=time_label, shrink=0.8)
        fig.suptitle(f'Ungated feedback trajectories | sigma_{p_label} = sigma_{v_label} = {sigma:g}\n{title}\ndt = {args.dt:g}, sampling = {args.sample_interval:g}, seed = {seed}', fontsize=11)
        fig.savefig(args.out / f'trajectories_{stem}.png', dpi=180)
        plt.close(fig)
        axes[0, j].set(title=f'sigma_{p_label} = sigma_{v_label} = {sigma:g}', xlabel=x_label+'1', ylabel=x_label+'2')
        axes[0, j].set_aspect('equal', adjustable='datalim')
        axes[0, j].scatter(0, 0, color='black', s=25, zorder=3)
        axes[0, j].legend(fontsize=7)
        axes[0, j].grid(alpha=0.2)
        dax[0, j].set(title=f'sigma_{p_label} = sigma_{v_label} = {sigma:g}', ylabel=f'Speed |{v_label}|')
        dax[1, j].set(xlabel=time_label, ylabel=f'Polarity magnitude |{p_label}|')
        scenarios.append({'sigma_p': sigma, 'sigma_v': sigma, 'seed': seed, 'data': f'{stem}.npz'})
    for fig, name in [(overview, 'trajectories_overview.png'), (diagnostics, 'speed_and_polarity.png')]:
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(args.out / name, dpi=180)
        plt.close(fig)
    metadata = {'coupling': 'ungated', 'parameters': asdict(params), 'mode': 'dimensional' if args.dimensional else 'nondimensional',
                'duration': args.duration, 'dt': args.dt, 'sample_interval': args.sample_interval,
                p_label+'0': initial.tolist(), v_label+'0': initial.tolist(),
                x_label+'0': 'zero', 'position_unit': 'user dimensional units' if args.dimensional else 'v_s * tau_p', 'method': 'Euler-Maruyama',
                'scenarios': scenarios, 'generated_utc': datetime.now(timezone.utc).isoformat(),
                'numpy': np.__version__, 'matplotlib': matplotlib.__version__,
                'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'Saved {6 * len(scenarios)} trajectories across {len(scenarios)} scenarios to {args.out.resolve()}')


if __name__ == '__main__':
    main()
