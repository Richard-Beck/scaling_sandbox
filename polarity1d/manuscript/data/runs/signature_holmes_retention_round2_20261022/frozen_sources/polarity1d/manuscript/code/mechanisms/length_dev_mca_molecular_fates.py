"""First-exit MCA residence and splitting-probability companion to Figure X.

Passive continuous-time, two-state tracer interpretation; never applies the
native staggered exchange residual or uniform concentration correction to labels.
Workers use frozen inputs and a frozen copy of this script. No native dynamics
are integrated or modified. Run `launch`, then `render` after all workers finish.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import solve_banded

ROOT = next(p for p in Path(__file__).resolve().parents if p.name == 'manuscript')
RUN = ROOT / 'data/runs/debelly_mca_molecular_fates_20261010_v2'
DEV = ROOT / 'figures/length_mechanisms_20261008/dev'
OUT = DEV / 'mca_molecular_fates'
X = DEV / 'figure_x_reversal_overview/supporting'
STEM = 'MCA_molecular_fates'
LENGTHS = [23, 31, 39]
PYTHON = '/home/4473331/.conda/envs/scaling_sandbox/bin/python'
CONFIGS = [('main', 1601, .5), ('time_refined', 1601, .25),
           ('space_coarse', 801, .5), ('space_refined', 3201, .5), ('native_mesh', 201, .5)]


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def read(path):
    return json.loads(path.read_text())


def output(name):
    category = {'.pdf': 'vectors', '.svg': 'vectors', '.csv': 'data',
                '.npz': 'data', '.json': 'metadata', '.txt': 'metadata',
                '.html': 'viewers'}.get(Path(name).suffix)
    p = OUT / ('supporting/' + category if category else '') / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def inputs():
    manifest = read(X / 'metadata/Figure_X_reversal_overview.manifest.json')
    audit = read(X / 'metadata/Figure_X_reversal_overview_audit.json')
    assert manifest['passed'] and audit['passed']
    sources = {str(p.relative_to(ROOT)): sha(p) for p in [
        X / 'metadata/Figure_X_reversal_overview.manifest.json',
        X / 'metadata/Figure_X_reversal_overview_audit.json',
        X / 'data/Figure_X_reversal_overview_fields.npz']}
    native = ROOT / 'data/runs/length_narrative_memory_20261009/frozen_sources/code/core/src/polarity1d/debelly.py'
    sources[str(native.relative_to(ROOT))] = sha(native)
    tree = ast.parse(native.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'DeBellyParameters')
    parameters = {n.target.id: ast.literal_eval(n.value) for n in cls.body if isinstance(n, ast.AnnAssign)}
    endpoints = DEV / 'deBelly_redraft/supporting/data/Figure_alpha_steady_state_mechanism_endpoint_metrics.csv'
    sources[str(endpoints.relative_to(ROOT))] = sha(endpoints)
    rows = {int(r['length_um']): r for r in csv.DictReader(endpoints.open())}
    jobs = []
    for L in LENGTHS:
        a = next(a for a in audit['checks'] if a['length_um'] == L)
        archive = ROOT / a['archive']
        meta = read(archive.with_suffix('.json'))
        cp = next(ROOT / name for name, digest in manifest['source_sha256'].items()
                  if digest == a['checkpoint_sha256'])
        assert sha(archive) == meta['archive_sha256'] and sha(cp) == a['checkpoint_sha256']
        assert meta['metrics']['full_polarity_success']
        for p in [archive, archive.with_suffix('.json'), cp]:
            sources[str(p.relative_to(ROOT))] = sha(p)
        jobs.append(dict(length=L, archive=a['archive'], checkpoint=str(cp.relative_to(ROOT)),
                         field_tail_pass=a['native_field_tail_pass'],
                         threshold=float(rows[L]['rac_switch_threshold']),
                         steady_left_command=float(rows[L]['derived_left_command_signed_um_per_s']),
                         steady_right_command=float(rows[L]['derived_right_command_signed_um_per_s'])))
    return jobs, parameters, sources


class Background:
    def __init__(self, job, parameters):
        self.L, self.p = job['length'], parameters
        with np.load(ROOT / job['archive']) as z:
            self.x = z['x'].copy()
            self.duration = float(z['pulse_time'][-1])
            self.t = np.r_[z['pulse_time'], self.duration + z['recovery_time'][1:]]
            self.f = {key: np.r_[z['pulse_field_' + key], z['recovery_field_' + key][1:]]
                      for key in ['rac', 'tension', 'velocity', 'displacement',
                                  'mca_bound', 'mca_unbound', 'left_boundary_displacement',
                                  'right_boundary_displacement']}
            with np.load(ROOT / job['checkpoint']) as cp:
                for key, values in self.f.items():
                    np.testing.assert_array_equal(values[0], cp['state_' + key])
                for key in self.f:
                    np.testing.assert_array_equal(z['pulse_field_' + key][-1], z['recovery_field_' + key][0])
        assert np.all(np.diff(self.t) > 0) and np.isfinite(np.concatenate(list(self.f.values()))).all()
        switch = .5 * (1 + np.tanh((self.f['rac'][:, [0, -1]] - job['threshold']) / parameters['protrusion_switch_scale']))
        self.edges = parameters['protrusion_speed'] * switch / (1 + self.f['tension'][:, [0, -1]]**2 / parameters['protrusion_tension_scale'])
        self.edges[:, 1] *= -1
        np.testing.assert_allclose(self.edges[0], [job['steady_left_command'], job['steady_right_command']], atol=1e-13)
        self.mv = np.diff(self.f['displacement'], axis=0) / np.diff(self.t)[:, None]
        self.me = np.diff(np.column_stack([self.f[k][:, 0] for k in ['left_boundary_displacement', 'right_boundary_displacement']]), axis=0) / np.diff(self.t)[:, None]
        with np.load(X / 'data/Figure_X_reversal_overview_fields.npz') as z:
            tx = z[f'L{self.L}_time']; use = (tx >= 0) & (tx <= 1200)
            self.display = tx[use].copy()
            for key in ['rac', 'tension', 'velocity', 'mca_bound']:
                ids = np.searchsorted(self.t, self.display)
                np.testing.assert_array_equal(self.f[key][ids], z[f'L{self.L}_{key}'][use])

    def grid(self, cells):
        return (np.arange(cells) + .5) * self.L / cells

    def rates(self, when, cells, stationary=None):
        faces = np.linspace(0, self.L, cells + 1)
        if stationary:
            i = 0 if stationary == 'initial' else -1
            v = self.f['velocity'][i]
            e = self.edges[i]
            if stationary == 'reflected_initial':
                v, e = -self.f['velocity'][0, ::-1], -self.edges[0, ::-1]
            bface = np.interp(faces, np.r_[0, self.x, self.L], np.r_[e[0], v, e[1]])
            uface = np.zeros_like(bface)
        else:
            i = min(len(self.t) - 2, max(0, int(np.searchsorted(self.t, when, side='right') - 1)))
            w = (when - self.t[i]) / (self.t[i + 1] - self.t[i])
            v = (1 - w) * self.f['velocity'][i] + w * self.f['velocity'][i + 1]
            e = (1 - w) * self.edges[i] + w * self.edges[i + 1]
            bface = np.interp(faces, np.r_[0, self.x, self.L], np.r_[e[0], v, e[1]])
            uface = np.interp(faces, np.r_[0, self.x, self.L], np.r_[self.me[i, 0], self.mv[i], self.me[i, 1]])
        return generator(bface, uface, self.L / cells, self.p)

    def weights(self, times, cells):
        x = self.grid(cells)
        out = []
        for t in times:
            i = min(len(self.t) - 2, max(0, int(np.searchsorted(self.t, t, side='right') - 1)))
            w = (t - self.t[i]) / (self.t[i + 1] - self.t[i])
            b, u = [np.interp(x, self.x, (1 - w) * self.f[k][i] + w * self.f[k][i + 1]) for k in ['mca_bound', 'mca_unbound']]
            out.append(b / (b + u))
        return np.array(out)


def generator(bface, uface, dx, p):
    """Backward subgenerator, interleaved [bound_i, unbound_i].

    Interior diffusion is reflecting at endpoints; only outward advection kills
    labels. Inward advective concentration input contains no cohort labels.
    """
    n = len(bface) - 1
    q = np.zeros((5, 2 * n))
    left, right = np.zeros(2 * n), np.zeros(2 * n)
    for state, face, diffusion in [(0, bface, p['mca_bound_diffusion']), (1, uface, p['mca_unbound_diffusion'])]:
        if diffusion == 0:
            rightward = np.maximum(face[1:-1], 0) / dx
            leftward = np.maximum(-face[1:-1], 0) / dx
        else:
            # Exponentially fitted flux avoids the extra diffusion introduced
            # by first-order upwinding in a molecular tracer interpretation.
            pe = face[1:-1] * dx / diffusion
            def bernoulli(z):
                out = np.empty_like(z)
                small, high, low = abs(z) < 1e-5, z > 50, z < -50
                middle = ~(small | high | low)
                out[small] = 1-z[small]/2+z[small]**2/12
                out[high] = z[high]*np.exp(-z[high])
                out[low] = -z[low]
                out[middle] = z[middle]/np.expm1(z[middle])
                return out
            rightward = diffusion / dx**2 * bernoulli(-pe)
            leftward = diffusion / dx**2 * bernoulli(pe)
        # Row i is current state; column j is next state.
        q[0, state + 2::2] = rightward
        q[4, state:-2:2] = leftward
        q[2, state:-2:2] -= rightward
        q[2, state + 2::2] -= leftward
        left[state] = max(-face[0], 0) / dx
        right[2 * n - 2 + state] = max(face[-1], 0) / dx
    q[1, 1::2] = p['mca_off']
    q[3, 0::2] = p['mca_on']
    q[2, 0::2] -= p['mca_off']
    q[2, 1::2] -= p['mca_on']
    q[2] -= left + right
    return q, left, right


def apply_banded(q, y):
    out = q[2, :, None] * y
    for offset in [-2, -1, 1, 2]:
        if offset > 0:
            out[:-offset] += q[2 - offset, offset:, None] * y[offset:]
        else:
            out[-offset:] += q[2 - offset, :offset, None] * y[:offset]
    return out


def stationary(q, left, right):
    rhs = np.column_stack([np.ones(len(left)), left, right])
    result = solve_banded((2, 2), -q, rhs)
    residual = float(np.max(np.abs(apply_banded(q, result) + rhs)))
    assert residual < 1e-6 and result[:, 0].min() > 0
    np.testing.assert_allclose(result[:, 1:].sum(axis=1), 1, atol=2e-8)
    assert result[:, 1:].min() > -2e-8 and result[:, 1:].max() < 1 + 2e-8
    return result, residual


def backward_step(q, forcing, y, h):
    """Second-order L-stable SDIRK, with midpoint background coefficients."""
    gamma = 1 - 1 / np.sqrt(2)
    matrix = -gamma*h*q
    matrix[2] += 1
    first = solve_banded((2,2), matrix, y+gamma*h*forcing, check_finite=False)
    return solve_banded((2,2), matrix, y+(1-gamma)/gamma*(first-y)+gamma*h*forcing, check_finite=False)


def sample_states(bg, values, cells):
    values = values.reshape(cells,2,-1)
    return np.stack([np.column_stack([np.interp(bg.x,bg.grid(cells),values[:,s,k]) for k in range(values.shape[-1])])
                     for s in range(2)],axis=1).reshape(2*len(bg.x),-1)


def propagate(bg, cells, step, closures=('terminal', 'reflected_initial'), end=None):
    end = bg.t[-1] if end is None else float(end)
    outputs = bg.display[bg.display <= end]
    # Preserve native sample/withdrawal events; subdivide each interval.
    knots = np.unique(np.r_[bg.t[bg.t <= end], outputs, end])
    segments = [np.linspace(a, b, max(1, int(np.ceil((b - a) / step - 1e-10))) + 1)[:-1] for a, b in zip(knots[:-1], knots[1:])]
    grid = np.r_[np.concatenate(segments), end]
    saved = np.searchsorted(grid, outputs)
    assert np.max(abs(grid[saved] - outputs)) < 1e-10
    terms = []
    terminal_residual = []
    for closure in closures:
        y, residual = stationary(*bg.rates(end, cells, stationary=closure))
        terms.append(y)
        terminal_residual.append(residual)
    y = np.column_stack(terms + [np.ones((2 * cells, 1))])
    result = np.empty((len(outputs), 2 * len(bg.x), y.shape[1]))
    lookup = {int(j): i for i, j in enumerate(saved)}
    max_probability_error = 0.
    minimum_probability, maximum_probability = 1., 0.
    for j in range(len(grid) - 2, -1, -1):
        h = grid[j + 1] - grid[j]
        q, left, right = bg.rates(.5 * (grid[j + 1] + grid[j]), cells)
        forcing = np.zeros_like(y)
        for c in range(len(closures)):
            forcing[:, 3*c] = 1
            forcing[:, 3*c+1] = left
            forcing[:, 3*c+2] = right
        y = backward_step(q,forcing,y,h)
        max_probability_error = max(max_probability_error, float(abs(y[:, 1] + y[:, 2] - 1).max()))
        minimum_probability = min(minimum_probability, float(y[:, 1:3].min()))
        maximum_probability = max(maximum_probability, float(y[:, 1:3].max()))
        if j in lookup:
            result[lookup[j]] = sample_states(bg,y,cells)
    assert minimum_probability > -2e-8 and maximum_probability < 1 + 2e-8
    assert max_probability_error < 2e-8 and result[:, :, 0].min() > 0
    return outputs, result, dict(steps=len(grid)-1, maximum_step_s=float(np.diff(grid).max()),
                                terminal_stationary_equation_residuals=terminal_residual,
                                probability_sum_max_error=max_probability_error,
                                probability_minimum=minimum_probability, probability_maximum=maximum_probability)


def mix(values, weights):
    return weights[..., None] * values[:, 0::2] + (1 - weights[..., None]) * values[:, 1::2]


def selfcheck():
    """Analytic transport, exchange, reflection and forward/backward tests."""
    p = dict(mca_on=5., mca_off=1.4, mca_bound_diffusion=0., mca_unbound_diffusion=0.)
    n, dx, speed = 9, .2, .3
    q, l, r = generator(np.full(n+1, speed), np.full(n+1, speed), dx, p)
    y, _ = stationary(q, l, r)
    expected = np.repeat(np.arange(n, 0, -1) * dx / speed, 2)
    np.testing.assert_allclose(y[:, 0], expected, atol=1e-12)
    np.testing.assert_allclose(y[:, 2], 1, atol=1e-12)
    p.update(mca_bound_diffusion=.3, mca_unbound_diffusion=.3)
    faces = np.zeros(n+1); faces[0], faces[-1] = -.07, .07
    q, l, r = generator(faces, faces, dx, p)
    y, _ = stationary(q, l, r)
    i = np.arange(n)
    exact_t = n * dx / (2 * .07) + i * (n-1-i) * dx**2 / (2*.3)
    exact_p = (i*dx + .3/.07) / ((n-1)*dx + 2*.3/.07)
    np.testing.assert_allclose(y[:, 0], np.repeat(exact_t, 2), atol=1e-10)
    np.testing.assert_allclose(y[:, 2], np.repeat(exact_p, 2), atol=1e-12)
    np.testing.assert_allclose(apply_banded(q, np.ones((2*n, 1)))[:, 0], -l-r, atol=1e-12)
    # Independently assemble a dense generator from jump events.
    dense = np.zeros((2*n, 2*n))
    for j in range(n):
        for state in [0, 1]:
            a = 2*j+state
            for k in [j-1, j+1]:
                if 0 <= k < n:
                    dense[a, 2*k+state] += .3/dx**2
            dense[a, 2*j+1-state] += p['mca_off'] if state == 0 else p['mca_on']
            dense[a, a] = -dense[a].sum()-l[a]-r[a]
    np.testing.assert_allclose(apply_banded(q, np.eye(2*n)), dense, atol=1e-12)
    np.testing.assert_allclose(np.linalg.solve(-dense, np.column_stack([np.ones(2*n), l, r])), y, atol=1e-10)
    # A time-varying protocol, integrated independently forward and backward.
    h = .04; schedule = []
    terminal = y.copy()
    backward = terminal.copy()
    for j in range(100):
        qq = dense.copy() * (1 + .3*np.sin(j*.17))
        ll, rr = l*(1+.3*np.sin(j*.17)), r*(1+.3*np.sin(j*.17))
        schedule.append((qq, ll, rr))
    gamma = 1 - 1/np.sqrt(2)
    for qq, ll, rr in schedule[::-1]:
        band = np.zeros_like(q)
        for offset in range(-2,3):
            for col in range(2*n):
                row = col-offset
                if 0 <= row < 2*n:
                    band[2-offset,col] = qq[row,col]
        backward = backward_step(band,np.column_stack([np.ones(2*n),ll,rr]),backward,h)
    cohort = np.zeros(2*n); cohort[5] = 1
    budget = np.zeros(3)
    for qq, ll, rr in schedule:
        inverse = np.linalg.inv(np.eye(2*n)-gamma*h*qq)
        transition = (1-(1-gamma)/gamma)*inverse+(1-gamma)/gamma*inverse@inverse
        reward = (1-gamma)*h*inverse@inverse+gamma*h*inverse
        occupied = reward.T@cohort
        budget += np.array([occupied.sum(),ll@occupied,rr@occupied])
        cohort = transition.T@cohort
    budget += cohort@terminal
    np.testing.assert_allclose(budget, backward[5], atol=1e-11)
    return dict(passed=True, analytic_uniform_advection=True, analytic_symmetric_diffusion_and_edge_loss=True,
                independently_assembled_dense_generator=True, zero_diffusive_boundary_loss=True,
                time_varying_forward_backward_agreement=True, forward_backward_max_error=float(abs(budget-backward[5]).max()))


def launch():
    assert not RUN.exists(), 'Campaign exists; do not overwrite frozen producer inputs'
    tests = selfcheck()
    jobs, parameters, sources = inputs()
    for folder in ['workers', 'logs', 'frozen_sources']:
        (RUN / folder).mkdir(parents=True, exist_ok=True)
    script = RUN / 'frozen_sources' / Path(__file__).name
    shutil.copy2(Path(__file__), script)
    sources[str(script.relative_to(ROOT))] = sha(script)
    write(RUN/'manifest.json', dict(generated_at=now(), jobs=jobs, parameters=parameters,
                                   source_sha256=sources, configurations=CONFIGS, selfcheck=tests))
    write(RUN/'status.json', dict(generated_at=now(), producer='running', complete=False))
    command = f'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; {PYTHON} {script} worker --slot "$SLURM_ARRAY_TASK_ID"'
    args = ['sbatch', '--parsable', '--qos=small', '--ntasks=1', '--cpus-per-task=1', '--mem=3G', '--time=00:20:00',
            '--array=0-2', '--job-name=mca-fates', '--output='+str(RUN/'logs/%A_%a.out'),
            '--error='+str(RUN/'logs/%A_%a.err'), '--wrap', command]
    job = subprocess.check_output(args, text=True).strip()
    write(RUN/'launch.json', dict(generated_at=now(), array=job, resources=dict(qos='small', cpus=1, memory='3G'), configs=CONFIGS))
    print(job)


def worker(slot):
    start = time.monotonic()
    m = read(RUN/'manifest.json')
    for name, digest in m['source_sha256'].items():
        assert sha(ROOT/name) == digest, name
    bg = Background(m['jobs'][slot], m['parameters'])
    audits, arrays = [], {}
    for name, cells, step in m['configurations']:
        print(bg.L, name, 'start', flush=True)
        t, values, audit = propagate(bg, cells, step)
        weights = bg.weights(t, len(bg.x))
        average = mix(values, weights)
        arrays.update({name+'_time': t, name+'_x': bg.x, name+'_conditioned': values,
                       name+'_weights': weights, name+'_mixed': average})
        steady, residual = stationary(*bg.rates(0, cells, stationary='initial'))
        steady = sample_states(bg,steady,cells)
        arrays[name+'_steady_conditioned'] = steady
        arrays[name+'_steady_mixed'] = mix(steady[None], bg.weights([0], len(bg.x)))[0]
        audit.update(name=name, cells=cells, requested_step_s=step, steady_residual=residual)
        audits.append(audit)
        print(bg.L, name, 'finished', round(time.monotonic()-start, 2), flush=True)
    # Alternative tail closure at 5000 s measures dependence on saved future extent.
    t, values, audit = propagate(bg, 1601, .5, closures=('terminal',), end=bg.duration+5000)
    arrays['short_tail_mixed'] = mix(values, bg.weights(t, len(bg.x)))
    t, values, audit_late = propagate(bg, 1601, .5, closures=('terminal',), end=bg.duration+9000)
    arrays['late_tail_mixed'] = mix(values, bg.weights(t, len(bg.x)))
    arrays['duration'] = np.array(bg.duration)
    p = RUN/'workers'/f'{slot:03d}.npz'
    np.savez_compressed(p, **arrays)
    write(p.with_suffix('.json'), dict(generated_at=now(), length=bg.L, archive_sha256=sha(p),
                                     audits=audits, short_tail_audit=audit, late_tail_audit=audit_late, wall_seconds=time.monotonic()-start,
                                     native_background_last_time_s=float(bg.t[-1]), native_field_tail_pass=m['jobs'][slot]['field_tail_pass'],
                                     maximum_inward_membrane_edge_speed_um_per_s=float(max(bg.me[:,0].max(), -bg.me[:,1].min(), 0))))
    print(bg.L, 'complete', time.monotonic()-start, flush=True)


def comparison(a, b):
    return dict(mean_time_max_relative_difference=float((abs(a[..., 0]-b[..., 0])/a[..., 0]).max()),
                mean_time_max_absolute_difference_s=float(abs(a[...,0]-b[...,0]).max()),
                mean_time_max_difference_relative_to_spatial_max=float((abs(a[...,0]-b[...,0])/a[...,0].max(axis=1)[:,None]).max()),
                right_probability_max_absolute_difference=float(abs(a[..., 2]-b[..., 2]).max()))


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    m = read(RUN/'manifest.json')
    for name, digest in m['source_sha256'].items():
        assert sha(ROOT/name) == digest, name
    tests = selfcheck()
    arrays, traces, audits, rows = {}, {}, [], []
    for slot, job in enumerate(m['jobs']):
        path = RUN/'workers'/f'{slot:03d}.npz'
        meta = read(path.with_suffix('.json'))
        assert sha(path) == meta['archive_sha256']
        with np.load(path) as z:
            a = {k: z[k].copy() for k in z.files}
        L, t, x = job['length'], a['main_time'], a['main_x']
        main = a['main_mixed']
        checks = dict(length_um=L, worker=meta, time_refinement=comparison(main, a['time_refined_mixed']),
                      terminal_vs_reflected_basal_tail=comparison(main, main[..., [3,4,5]]),
                      recovery_5000_vs_10000_s=comparison(main, a['short_tail_mixed']),
                      recovery_9000_vs_10000_s=comparison(main, a['late_tail_mixed']),
                      maximum_probability_surviving_to_saved_history_end=float(main[..., -1].max()))
        for name in ['space_coarse', 'space_refined', 'native_mesh']:
            other = np.stack([np.array([np.interp(x, a[name+'_x'], row) for row in a[name+'_mixed'][:, :, k]]) for k in range(3)], axis=-1)
            checks[name+'_vs_1601'] = comparison(main, other)
        # Pointwise relative errors overemphasize a few-second residence beside
        # an open end. Require both relative and absolute time-step agreement.
        assert checks['time_refinement']['mean_time_max_relative_difference'] < .02
        assert checks['time_refinement']['mean_time_max_absolute_difference_s'] < 3
        assert checks['time_refinement']['right_probability_max_absolute_difference'] < .0001
        assert checks['terminal_vs_reflected_basal_tail']['mean_time_max_relative_difference'] < .01
        assert checks['terminal_vs_reflected_basal_tail']['right_probability_max_absolute_difference'] < .005
        assert checks['recovery_9000_vs_10000_s']['mean_time_max_relative_difference'] < .01
        assert checks['space_refined_vs_1601']['mean_time_max_relative_difference'] < .025
        assert checks['space_refined_vs_1601']['right_probability_max_absolute_difference'] < .005
        audits.append(checks)
        pre = np.arange(-60., 0., .5)
        v = np.concatenate([np.broadcast_to(a['main_steady_mixed'], (len(pre), len(x), 3)), main[..., :3]])
        time_axis = np.r_[pre, t]
        traces[L] = dict(time=time_axis, x=x, values=v, duration=float(a['duration']))
        # Full refinement arrays remain in immutable worker evidence; the
        # companion export contains main fields and both initial-state results.
        arrays.update({f'L{L}_{k}': value for k, value in a.items() if k.startswith('main_') or k == 'duration'})
        arrays[f'L{L}_display_time'] = time_axis
        arrays[f'L{L}_display_mean_residence_s'] = v[...,0]
        arrays[f'L{L}_display_p_left'] = v[...,1]
        arrays[f'L{L}_display_p_right'] = v[...,2]
        for when in [None, 0., 60., float(a['duration']), 1200.]:
            i = 0 if when is None else int(np.argmin(abs(t-when)))
            vals = a['main_steady_mixed'] if when is None else main[i]
            for fraction in [.05, .25, .5, .75, .95]:
                rows.append(dict(length_um=L, phase='stationary_basal' if when is None else 'reversal',
                                 time_after_onset_s='' if when is None else float(t[i]), position_fraction=fraction,
                                 mean_remaining_residence_s=float(np.interp(L*fraction,x,vals[:,0])),
                                 p_left=float(np.interp(L*fraction,x,vals[:,1])),
                                 p_right=float(np.interp(L*fraction,x,vals[:,2])),
                                 cue_duration_s=float(a['duration']), native_field_tail_pass=job['field_tail_pass']))
    np.savez_compressed(output(STEM+'_fields.npz'), **arrays)
    with output(STEM+'_measurements.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig = plt.figure(figsize=(12, 7.8))
    gs = fig.add_gridspec(2, 4, width_ratios=[23,31,39,2.2], left=.105, right=.91, bottom=.17, top=.865, wspace=.18, hspace=.30)
    maximum = max(tr['values'][...,0].max() for tr in traces.values())
    upper = float(np.ceil(maximum/500)*500)
    for row in range(2):
        norm = Normalize(0, upper if row == 0 else 1)
        for col, L in enumerate(LENGTHS):
            tr = traces[L]; ax = fig.add_subplot(gs[row,col])
            mesh = ax.pcolormesh(tr['x'],tr['time'],tr['values'][...,0 if row == 0 else 2],shading='nearest',
                                 norm=norm,cmap='viridis' if row == 0 else 'RdBu_r',rasterized=True)
            ax.set_yscale('symlog',linthresh=120,linscale=.65)
            ax.set(xlim=(0,L),ylim=(1200,-60),xticks=[0,L/2,L],yticks=[-60,0,60,120,300,600,1200])
            ax.set_yticklabels(['−60','0','60','120','300','600','1200'])
            ax.tick_params(labelsize=9,length=3)
            if col:
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(('Mean remaining residence' if row == 0 else 'First-exit side')+'\nTime from cue onset (s)',labelpad=8)
            ax.axhline(0,color='white',lw=2); ax.axhline(0,color='#222222',lw=.65)
            ax.axhline(60,color='white',ls=':',lw=1.2)
            ax.axhline(tr['duration'],color='white',ls='--',lw=2)
            ax.axhline(tr['duration'],color='#222222',ls='--',lw=.7)
            ax.text(.025,.975,chr(65+3*row+col),transform=ax.transAxes,va='top',weight='bold',fontsize=11,
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.6))
            if row == 0:
                ax.set_title(f'L = {L} µm',fontsize=12,pad=12)
            else:
                ax.set_xlabel('Reference position x (µm)')
        cb = fig.colorbar(mesh,cax=fig.add_subplot(gs[row,3]))
        cb.ax.tick_params(labelsize=9)
        if row == 0:
            cb.set_label('Mean time to first exit (s)',labelpad=8)
        else:
            cb.set_ticks([0,.25,.5,.75,1]); cb.set_ticklabels(['0 · LHS','.25','.50','.75','1 · RHS'])
            cb.set_label('Probability of RHS first exit',labelpad=8)
    fig.suptitle('MCA molecular residence and exit fate',fontsize=15,y=.975)
    fig.text(.105,.93,'Local bound/unbound mixture · passive first-exit predictions on the Figure X histories',fontsize=10,color='#444444')
    fig.text(.105,.078,'−60–0 s: stationary basal reference, with no future cue.  From 0 s: predictions include the future cue and recovery.',fontsize=9)
    fig.text(.105,.046,'Solid: cue onset     Dotted: 60 s     Dashed: cue withdrawal     Zero diffusive flux at both ends',fontsize=9)
    for ext in ['png','pdf','svg']:
        fig.savefig(output(STEM+'.'+ext),dpi=220,facecolor='white')
    plt.close(fig)
    caption = (
        'MCA molecular residence and first-exit fate, companion to Figure X. Columns are L=23,31,39 µm, with the exact same '
        'corrected relative-20% histories and shortest tested canonical scalar-strength-passing pulses (232.244898,285.280612,809.081633 s). '
        'A–C: expected remaining time until first domain exit, including all intervening binding cycles. D–F: probability of first exit '
        'through the RHS; LHS probability is its complement. Each location/time represents an MCA molecule whose initial bound/unbound '
        'state is drawn from the local native pool fractions. Conditioned results for each initial state are also exported. These are future '
        'outcome predictions, not instantaneous fluxes, molecule age, or uninterrupted attachment lifetime. Positive-time predictions use '
        'the subsequent recorded cue and complete 10000-s recovery, followed by a stationary extrapolation of the final cortical velocity '
        'and endpoint commands with zero membrane velocity. The alternative reflected initial basal tail tests the unrecorded future. '
        'Premature switching to that same final-recorded terminal background at 5000 or 9000 s of recovery tests how much of the '
        'recorded recovery must be retained; no response tail is relabeled as fully converged. The −60–0-s band repeats an independently calculated '
        'stationary basal first-exit diagnostic with no future cue, analogous to Figure X\'s repeated basal field strip. It is not a '
        'pre-onset prediction of the later imposed cue. The change in diagnostic at zero is explicit. Time increases downward, linear through '
        '120 s then logarithmic, through 1200 s; column widths scale with length. The background is interpolated from recorded native fields '
        'without changing their dynamics. Cortex face velocities use piecewise linear space/time interpolation; membrane velocities are '
        'interval-average displacement derivatives, including endpoint motion. Bound/unbound diffusivities are 0.001/0.3 µm²/s, '
        'with reflecting diffusion at both endpoints; exchange is a conservative continuous-time two-state process with off/on rates '
        '1.4/5 s^-1. Only outward boundary advection removes labels. Inward advective input carries no cohort labels. Uniform native '
        'replenishment and the staggered-exchange numerical residual do not create labeled molecules; no inference is made about identity '
        'after first exit. This explicit passive tracer interpretation is distinct from exactly decomposing the native concentration '
        'mass-closure step. The backward subgenerator equations use second-order L-stable SDIRK with midpoint background coefficients, '
        'nominal step 0.5 s and 1601 tracer cells with conservative exponentially fitted interior fluxes '
        'on the fixed 201-cell native background. All displays and comparisons are evaluated at the original 201 Figure X cell centers. '
        'Time-step halving, 801/3201/201-cell tracer calculations, two terminal continuations and '
        'analytic/independent forward-backward checks accompany the figure. Native field-tail flags remain separate (false,true,false). '
        'Solid/dotted/dashed lines mark onset,60 s,and withdrawal. All three lengths share each row\'s absolute color scale without clipping. '
        'Raw histories are immutable; no new native integration or main report change.')
    output(STEM+'_caption.txt').write_text(caption+'\n')
    audit = dict(generated_at=now(), passed=True, selfcheck=tests, checks=audits,
                 background_identity_with_Figure_X=True, raw_sources_unchanged=True,
                 new_native_integrations=0, tracer_interpretation='conservative two-state first exit; unlabeled replenishment and inward input')
    write(output(STEM+'_audit.json'), audit)
    methods = ('At each point, remaining residence and first-exit probabilities condition on the entire future protocol. '
               'The upper basal band instead holds the verified initial background fixed forever; it has only LHS loss. '
               'The tagged backward generator has conservative transport/exchange and nonnegative outward absorption rates. '
               'It does not propagate the native split-exchange numerical residual or mass correction as tagged matter. '
               'Inward membrane boundary input during relaxation is unlabeled; it neither absorbs nor creates tags. '
               'This makes first-exit identity explicit while leaving native mechanics untouched. '
               'The final background extrapolation is tested rather than asserted to be a verified steady response.')
    html = ('<!doctype html><html lang="en"><meta charset="utf-8"><title>MCA molecular fates</title>'
            '<style>body{max-width:1250px;margin:32px auto;padding:0 24px;font:16px/1.55 system-ui}img{width:100%;height:auto}p{max-width:1100px}td,th{padding:6px 12px;text-align:left;border-bottom:1px solid #ddd}table{border-collapse:collapse}</style>'
            '<h1>MCA residence and first-exit fate</h1><p>Generated '+now()+' · Producer finished</p>'
            '<a href="../vectors/'+STEM+'.pdf"><img src="../../'+STEM+'.png" alt="MCA residence and exit-fate kymographs"></a>'
            '<p><a href="../vectors/'+STEM+'.pdf">PDF</a> · <a href="../vectors/'+STEM+'.svg">SVG</a> · '
            '<a href="../data/'+STEM+'_fields.npz">Figure fields and state-conditioned results</a> · <a href="../data/'+STEM+'_measurements.csv">Selected measurements</a> · '
            '<a href="../metadata/'+STEM+'_audit.json">Numerical audit</a></p><h2>Interpretation</h2><p>'+methods+'</p>'
            '<h2>Numerical sensitivity</h2><table><tr><th>Length</th><th>Time halving, max residence change</th><th>1601→3201 cells</th><th>Alternative tail</th></tr>'+
            ''.join(f'<tr><td>{a["length_um"]} µm</td><td>{100*a["time_refinement"]["mean_time_max_relative_difference"]:.3f}%</td><td>{100*a["space_refined_vs_1601"]["mean_time_max_relative_difference"]:.3f}%</td><td>{100*a["terminal_vs_reflected_basal_tail"]["mean_time_max_relative_difference"]:.3f}%</td></tr>' for a in audits)+
            '</table><h2>Caption</h2><p>'+caption+'</p></html>')
    output(STEM+'.html').write_text(html)
    artifacts = [p for p in OUT.rglob('*') if p.is_file() and p.name != STEM+'.manifest.json' and p.name != 'README.md']
    sources = dict(m['source_sha256'])
    sources[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    for slot in range(3):
        for suffix in ['.npz','.json']:
            p = RUN/'workers'/f'{slot:03d}{suffix}'; sources[str(p.relative_to(ROOT))] = sha(p)
    write(output(STEM+'.manifest.json'), dict(artifact_root='../..',generated_at=now(),producer='finished',passed=True,
                                            legend=caption,description=caption,sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts},
                                            source_sha256=sources,audit=audit,
                                            reproduction=[PYTHON+' '+str(Path(__file__))+' render'], configurations=CONFIGS,
                                            native_background=dict(assay_mode='reflection_symmetric_extension',center_boundary_motion=False,fields='Exact Figure X histories'),
                                            tracer_parameters={k:m['parameters'][k] for k in ['mca_on','mca_off','mca_bound_diffusion','mca_unbound_diffusion']},
                                            array_definitions=dict(main_conditioned='time × interleaved [bound_i,unbound_i] × channels',
                                                main_mixed='time × original Figure X cell center × channels; local pool-weighted initial state',
                                                channels=['mean_remaining_residence_s_terminal_tail','p_left_terminal_tail','p_right_terminal_tail',
                                                          'mean_remaining_residence_s_reflected_basal_tail','p_left_reflected_basal_tail',
                                                          'p_right_reflected_basal_tail','probability_surviving_to_recorded_history_end'],
                                                main_weights='Local bound fraction at each time and original Figure X cell center',
                                                steady='Stationary initial-background diagnostic; columns mean_remaining_residence_s,p_left,p_right',
                                                display='Stationary basal band followed by future-protocol predictions; see caption for different conditioning')))
    write(RUN/'status.json', dict(generated_at=now(),producer='finished',complete=True,audit_passed=True,figure=str(OUT.relative_to(ROOT))))
    print(json.dumps(dict(output=str(output(STEM+'.png')),checks=audits,passed=True),indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['check','launch','worker','render'])
    parser.add_argument('--slot', type=int)
    args = parser.parse_args()
    if args.command == 'check':
        print(json.dumps(selfcheck(),indent=2))
    elif args.command == 'launch':
        launch()
    elif args.command == 'worker':
        worker(args.slot)
    else:
        render()
