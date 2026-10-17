"""New, frozen spatial probes; existing retained results are not changed."""
import argparse
import sys
import time
from datetime import datetime, timezone
import numpy as np
from design import *
sys.path.insert(0, str(ROOT/'polarity1d/src'))
from polarity1d.registry import make_model
from polarity1d.simulation import simulate
from polarity1d.metrics import spatial_contrast, profile_domain_count

def halfmass(profiles):
    n = profiles.shape[1]
    result = []
    for row in profiles:
        mass = np.maximum(row, 0)
        if mass.sum() < 1e-12:
            result.append(0.)
        else:
            result.append(np.interp(.5*mass.sum(), np.r_[0, np.cumsum(mass)], np.linspace(0, 1, n+1))-.5)
    return np.array(result)

def run(index, refine=False):
    model, L = jobs()[index]
    folder = OUT/'spatial'/f'{model}_{L:g}'
    if refine: folder = OUT/'refinement'/f'{model}_{L:g}'
    folder.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    for assay in ASSAYS:
        target = folder/f'{assay}.npz'
        protocol = Protocol(model, assay)
        if target.exists():
            with np.load(target) as previous:
                if np.isclose(previous['time'][-1],protocol.end): continue
        unit = protocol.unit
        # Coarse conditioning tail, dense perturbation/early response sampling.
        t = np.unique(np.r_[np.arange(0, min(protocol.condition, 10*unit)+.001, unit/20),
                            np.arange(0, protocol.end+.001, unit/2),
                            np.arange(protocol.condition, protocol.end+.001, unit/20),
                            protocol.event_times, protocol.end])
        t = t[(t >= 0) & (t <= protocol.end)]
        options = dict(length=L, cells=max(25, int(np.ceil(L/(.1 if refine else .2)))), times=t)
        if model == 'debelly': options['step'] = .25 if refine else .5
        else: options.update(backend='bdf', max_step=1. if refine else 2.5,
                            rtol=5e-7 if refine else 2e-6)
        trajectory = simulate(make_model(model, DEBELLY if model == 'debelly' else {}),
                              protocol=protocol, **options)
        p = halfmass(trajectory.observable)
        q = -halfmass(trajectory.fields['mca_bound']) if model == 'debelly' else np.full_like(p, np.nan)
        contrast = np.array([spatial_contrast(row) for row in trajectory.observable])
        fronts = np.array([profile_domain_count(row) for row in trajectory.observable])
        np.savez_compressed(target, time=t, p=p, q=q, cue=np.array([protocol.projected(x) for x in t]),
                            contrast=contrast, fronts=fronts, x=trajectory.x,
                            profiles=trajectory.observable.astype(np.float32),
                            zero_mass=np.sum(np.abs(trajectory.observable), axis=1)<1e-12)
        save_json(folder/f'{assay}.json', dict(model=model, length=L, assay=assay,
                   metadata=trajectory.metadata, generated=datetime.now(timezone.utc).isoformat(),
                   condition_seconds=protocol.condition, end_seconds=protocol.end,
                   source_hashes=digest([Path(__file__), HERE/'design.py']),
                   weak_conditioning_drift=float(np.ptp(p[(t >= protocol.condition-unit) & (t <= protocol.condition)]))))
        print(model, L, assay, 'seconds', round(time.monotonic()-start, 2), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('index', type=int)
    parser.add_argument('--refine', action='store_true')
    args = parser.parse_args()
    run(args.index, args.refine)
