"""Illustrative isotropic 2D extension of the ungated De Belly reference fit."""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from simulate import DimensionalParameters, drift, simulate

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=['ungated', 'self_activation'], default='self_activation')
    parser.add_argument('--noise-fractions', nargs='+', type=float, default=[.1, .3, .6, .9, 1.2, 1.8])
    args = parser.parse_args()
    levels = args.noise_fractions
    if any(f <= 0 or not np.isfinite(f) for f in levels) or len(set(levels)) != len(levels):
        parser.error('Noise fractions must be finite, positive and distinct')
    model_label = 'ungated' if args.model == 'ungated' else 'n=2 self-activation'
    source = HERE/f'debelly_matching/results_{args.model}/best_parameters.json'
    fit = json.loads(source.read_text())
    params = DimensionalParameters(**{key: fit.get(key, default) for key, default in asdict(DimensionalParameters()).items()})
    from scipy.optimize import brentq
    def velocity(p):
        return params.alpha*params.tau_v*p/(1+p/params.ps)
    def residual(p):
        v = velocity(p)
        return float(drift(np.array([[p, 0.]]), np.array([[v, 0.]]), params)[0][0, 0])
    grid = np.geomspace(1e-9, 10., 2000)
    roots = [brentq(residual, l, r, xtol=1e-15) for l, r in zip(grid[:-1], grid[1:])
             if residual(l)>0 and residual(r)<0]
    if not roots:
        raise ValueError('No attracting positive polarized amplitude found')
    pstar = min(roots, key=lambda p: abs(p-.06))
    vstar = velocity(pstar)
    angles = np.arange(6)*2*np.pi/6
    directions = np.column_stack([np.cos(angles), np.sin(angles)])
    p0, v0 = pstar*directions, vstar*directions
    dp, dv = drift(p0, v0, params)
    np.testing.assert_allclose(dp, 0, atol=1e-14)
    np.testing.assert_allclose(dv, 0, atol=1e-14)
    out = HERE/('results_fitted_2d' if args.model == 'ungated' else 'results_fitted_2d_self_activation')
    out.mkdir(exist_ok=True)
    duration, dt, stride = 10000., .1, 10
    ncols = min(3, len(levels))
    nrows = (len(levels)+ncols-1)//ncols
    overview, axes = plt.subplots(nrows, ncols, figsize=(16, 5*nrows), squeeze=False, layout='constrained')
    axes = axes.ravel()
    for unused in axes[len(levels):]:
        unused.axis('off')
    diagnostics, dax = plt.subplots(2, len(levels), figsize=(5*len(levels), 7), squeeze=False, sharex=True, layout='constrained')
    records = []
    for j, fraction in enumerate(levels):
        # Each uncoupled OU component has variance sigma^2*tau/2;
        # in 2D, total RMS magnitude equals sigma*sqrt(tau).
        sp = fraction*pstar/np.sqrt(params.tau_p)
        sv = fraction*vstar/np.sqrt(params.tau_v)
        seed = 20260925+j
        t, p, v, x = simulate(params, p0, v0, sp, sv, duration, dt, seed)
        t, p, v, x = t[::stride], p[::stride], v[::stride], x[::stride]
        np.testing.assert_array_equal(t, np.arange(10001))
        stem = 'noise_' + f'{fraction*100:02g}'.replace('.', 'p') + 'pct'
        np.savez_compressed(out/(stem+'.npz'), t=t, p=p, v=v, x=x)
        fig, grid = plt.subplots(2, 3, figsize=(12, 8), layout='constrained')
        for track, ax in enumerate(grid.flat):
            points = x[:, track]
            lines = LineCollection(np.stack([points[:-1], points[1:]], axis=1),
                                   cmap='viridis', norm=plt.Normalize(0, duration), linewidth=1.5)
            lines.set_array(t[:-1])
            ax.add_collection(lines)
            ax.autoscale()
            ax.scatter(*points[0], c='black', s=25, zorder=3, label='Start')
            ax.scatter(*points[-1], c='crimson', marker='X', s=35, zorder=3, label='End')
            ax.set(title=f'Track {track+1}', xlabel='Relative displacement x1', ylabel='Relative displacement x2')
            ax.set_aspect('equal', adjustable='datalim')
            ax.grid(alpha=.2)
            axes[j].plot(*points.T, lw=1.1, label=str(track+1))
            dax[0, j].plot(t, np.linalg.norm(p[:, track], axis=-1)/pstar, lw=.5, alpha=.7)
            dax[1, j].plot(t, np.linalg.norm(v[:, track], axis=-1)/vstar, lw=.5, alpha=.7)
        grid[0, 0].legend(fontsize=8)
        fig.colorbar(lines, ax=grid.ravel().tolist(), label='Time (s)', shrink=.8)
        fig.suptitle(f'Fitted {model_label} model: relative displacement | noise fraction {fraction:g}\n'
                     f'sigma_p={sp:.4g}, sigma_v={sv:.4g}; duration 10000 s; seed {seed}')
        fig.savefig(out/f'trajectories_{stem}.png', dpi=180)
        plt.close(fig)
        axes[j].scatter(0, 0, c='black', s=25, zorder=3)
        axes[j].set(title=f'Noise fraction {fraction:g}\nsigma_p={sp:.3g}, sigma_v={sv:.3g}',
                    xlabel='Relative displacement x1', ylabel='Relative displacement x2')
        axes[j].set_aspect('equal', adjustable='datalim')
        axes[j].grid(alpha=.2)
        dax[0, j].set(title=f'Noise fraction {fraction:g}', ylabel='|p| / polarized |p|')
        dax[1, j].set(xlabel='Time (s)', ylabel='|v| / polarized |v|')
        for ax in dax[:, j]:
            ax.axhline(1, color='black', ls='--', lw=.8)
        records.append(dict(noise_fraction=fraction, sigma_p=sp, sigma_v=sv, seed=seed,
                            data=stem+'.npz', max_p_norm=float(np.linalg.norm(p, axis=-1).max()),
                            max_v_norm=float(np.linalg.norm(v, axis=-1).max())))
    overview.suptitle(f'Isotropic 2D extension: {model_label}\n{6*len(levels)} illustrative tracks; 10000 s; displacement has no calibrated distance unit')
    overview.savefig(out/'trajectories_overview.png', dpi=180)
    plt.close(overview)
    diagnostics.suptitle('State amplitudes relative to the deterministic polarized state')
    diagnostics.savefig(out/'state_amplitudes.png', dpi=180)
    plt.close(diagnostics)
    meta = dict(generated_utc=datetime.now(timezone.utc).isoformat(), parameters=asdict(params),
                model=args.model, coupling=f'isotropic Euclidean-norm saturation; {model_label}', cue='none',
                pstar=pstar, vstar=vstar, p0=p0.tolist(), v0=v0.tolist(), x0='zero',
                duration_s=duration, dt_s=dt, sample_interval_s=1.,
                noise_convention='sigma_p=f*pstar/sqrt(tau_p); sigma_v=f*vstar/sqrt(tau_v); independent isotropic additive noise',
                noise_interpretation='f is the 2D RMS relative to polarized amplitude for uncoupled OU relaxation, not the coupled stationary variation',
                position='integral of relative velocity over seconds; arbitrary relative-displacement units',
                scenarios=records, numpy=np.__version__,
                source_hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [source, Path(__file__), HERE/'simulate.py']})
    (out/'metadata.json').write_text(json.dumps(meta, indent=2)+'\n')
    rows = '\n'.join(f'| {r["noise_fraction"]:g} | {r["sigma_p"]:.6g} | {r["sigma_v"]:.6g} |' for r in records)
    (out/'README.md').write_text(f'''# Illustrative 2D tracks from fitted parameters

[Overview](trajectories_overview.png), [state amplitudes](state_amplitudes.png).

Uses the {model_label} 29-um De Belly calibration without parameter changes. The
external calibration cue is absent. Scalar states become 2D vectors with
Euclidean saturation; this isotropic extension and its noise are assumptions,
not calibrated turning dynamics or experimentally measured migration.

The positive deterministic fixed point has |p|={pstar:.6g} and |v|={vstar:.6g}.
Each scenario starts six tracks there, aligned at evenly spaced angles; x=0.
No burn-in or path selection is applied. Integration is Euler–Maruyama,
dt=0.1 s (about tau_p/{params.tau_p/dt:.0f}), duration 10000 s, sampled every 1 s.

| Noise fraction | sigma_p | sigma_v |
| --- | --- | --- |
{rows}

Independent isotropic noise uses sigma=f*state_amplitude/sqrt(relaxation_time).
For an uncoupled 2D OU process this gives RMS magnitude f times the reference
amplitude. Coupling changes the actual variance; these are noise conventions,
not measured percentages of variation. The p and v noise amplitudes differ
because their fitted amplitudes and relaxation times differ. No clipping is used.

Displacement is the integral of relative v over seconds, with no calibrated
distance unit. Equal coordinate scaling preserves path geometry; individual
panels have independent limits. NPZ arrays contain t,p,v,x with state shape
10001 x 6 x 2. Metadata records seeds, initial conditions, assumptions and hashes.

Reproduce from repository root using the existing environment:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 conda run --no-capture-output -n scaling_sandbox python wgd_langevin/logistic_feedback/fitted_tracks.py --model {args.model} --noise-fractions {' '.join(map(str, levels))}
''')
    print(json.dumps(dict(pstar=pstar, vstar=vstar, scenarios=records, output=str(out)), indent=2))


if __name__ == '__main__':
    main()
