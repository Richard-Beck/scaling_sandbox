"""Fixed-physical-slope finite countergradients from corrected basal checkpoints.

Reuse the corrected runner's trials, recovery, gates, arrays, refinement and
independent audit. A bounded pilot chooses one slope per model, never per length.
"""
from pathlib import Path
import argparse
import copy
import json
import math
import shutil
import subprocess
import sys
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE / 'core'))
import corrected_reversal as r
import native as n
import polarity_establishment as e

DEFAULT = n.DATA / 'runs/fixed_slope_reversal_20261006'
PILOT = n.DATA / 'runs/fixed_slope_pilot_20261006'
BASE = n.DATA / 'runs/corrected_reversal_20261003/manifest.json'
FACTORS = [1., 4., 16.]
base_initialize = r.initialize
base_sources = r.source_hashes


class LinearGradient:
    """Added stimulus g*x, with x in the solver's native physical coordinates."""
    def __init__(self, slope, background=0., stop=None):
        if not np.isfinite(slope) or slope <= 0:
            raise ValueError('Gradient slope must be finite and positive')
        self.slope, self.background, self.stop = float(slope), float(background), stop
        self.event_times = () if stop is None else (float(stop),)

    def __call__(self, time, x, length):
        x = np.asarray(x, dtype=float)
        return np.full_like(x, self.background) if self.stop is not None and time >= self.stop else self.background + self.slope*x


def response_cue(case):
    return LinearGradient(case['gradient_slope'], case.get('background_stimulus', 0.))


def sources():
    result = base_sources()
    for path in [Path(__file__), CODE/'gradient_reversal_plot.py', CODE/'figure_builder.py', CODE/'build.py', BASE]:
        result[str(path.relative_to(n.ROOT))] = n.sha(path)
    return result


def freeze(out, manifest):
    r.write(out/'manifest.json', manifest)
    hashes = sources()
    for rel in hashes:
        dest = out/'frozen_sources'/rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(n.ROOT/rel, dest)
    r.write(out/'provenance.json', dict(generated_at=n.now(), sources=hashes,
                                      manifest_sha256=n.sha(out/'manifest.json')))


def initialize(out):
    selection = r.read(PILOT/'results/slope_selection.json')
    if not selection['complete'] or not selection['audit_passed']:
        raise RuntimeError('Slope pilot must complete and pass audit first')
    manifest = base_initialize(out)
    manifest.update(experiment='fixed_physical_slope_countergradient',
                    cue='Added right linear ramp g*x for each independent duration, then basal withdrawal. One positive physical slope per model, fixed across lengths and preparation widths. No crossfade or gradient reconditioning.',
                    slope_by_model=selection['slope_by_model'],
                    slope_selection=selection,
                    preparation_conditions=manifest['conditions'])
    for case in manifest['cases']:
        case['gradient_slope'] = selection['slope_by_model'][case['model']]
        case['response_peak_at_nominal_length'] = case['gradient_slope']*case['length']
        case['response_shape'] = 'linear_fixed_physical_slope'
    freeze(out, manifest)
    r.publish(out, 'initialized')
    return manifest


def pilot_init(out):
    if out.exists():
        raise ValueError('Pilot requires a fresh destination')
    cases = r.read(BASE)['cases']; selected = []
    for model in r.MODELS:
        eligible = sorted((c for c in cases if c['model']==model and c['condition']=='relative_20pct' and c['eligible']), key=lambda c:c['length'])
        for i in [0, len(eligible)//2, len(eligible)-1]:
            for factor in FACTORS:
                case = copy.deepcopy(eligible[i])
                case.update(case_index=len(selected), source_case_index=case['case_index'],
                            gradient_factor=factor, gradient_slope=factor*case['amplitude']/e.REFERENCES[model])
                maximum = max(case['durations'])
                case['pilot_durations'] = [0.] + [maximum/q for q in [1000.,100.,10.,1.]]
                selected.append(case)
    for folder in ['work/chunks', 'work/logs', 'results', 'frozen_sources']:
        (out/folder).mkdir(parents=True, exist_ok=True)
    manifest = dict(created_at=n.now(), experiment='fixed_slope_calibration', cases=selected,
                    factors=FACTORS, pilot_preparation='relative_20pct; shortest, middle, longest qualifying length',
                    selection_rule='Prefer the lowest factor achieving both gates at all three lengths. Otherwise maximize number of full-strength successes, then orientation successes, then prefer the lower factor; retain all failures. Each success means at least one of the four tested positive durations passes.',
                    full_campaign=str(DEFAULT), recovery_horizon=r.RECOVERY)
    freeze(out, manifest)
    r.write(out/'results/status.json', dict(generated_at=n.now(), complete=False, producer='not started', cases=len(selected)))


def pilot_submit(out, action, array=None, dependency=None):
    command = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '
    import shlex
    command += shlex.join([sys.executable, str(Path(__file__)), action, '--out', str(out)])
    if array is not None:
        command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    opts = ['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=02:00:00',
            '--job-name=fixed-slope-'+action,'--chdir='+str(n.ROOT),
            '--output='+str(out/'work/logs'/f'{action}_%A_%a.out'),
            '--error='+str(out/'work/logs'/f'{action}_%A_%a.err')]
    if array is not None: opts.append('--array='+array)
    if dependency: opts.append('--dependency=afterany:'+dependency)
    return subprocess.check_output(opts+['--wrap',command], text=True).strip().split(';')[0]


def pilot_launch(out):
    manifest = r.verify(out)
    if (out/'launch.json').exists():raise ValueError('Pilot already launched')
    array = pilot_submit(out, 'pilot-worker', f'0-{len(manifest["cases"])-1}')
    controller = pilot_submit(out, 'pilot-finish', dependency=array)
    r.write(out/'launch.json', dict(generated_at=n.now(), array=array, controller=controller, stage='running'))
    r.write(out/'results/status.json', dict(generated_at=n.now(), complete=False, producer='running', cases=len(manifest['cases'])))


def pilot_worker(out, slot):
    case = r.verify(out)['cases'][slot]; rows = []; arrays = {}
    for i,duration in enumerate(case['pilot_durations']):
        row, values = r.trial(case, duration)
        row['trace_prefix'] = f'trial{i}_'
        arrays.update({row['trace_prefix']+key:value for key,value in values.items()})
        rows.append(row)
    path = out/'work/chunks'/f'pilot_{slot:03d}.npz'
    np.savez_compressed(path, **arrays)
    r.write(path.with_suffix('.json'), dict(case=case, tests=rows, trace_sha256=n.sha(path)))


def pilot_finish(out):
    manifest = r.verify(out)
    sys.path.insert(0, str(CODE))
    import audit_reversal as audit
    groups = {}; table = []
    for case in manifest['cases']:
        path = out/'work/chunks'/f'pilot_{case["case_index"]:03d}.npz'
        result = r.read(path.with_suffix('.json'))
        if result['case'] != case or n.sha(path) != result['trace_sha256']:
            raise RuntimeError('Pilot archive identity failed')
        if [t['duration'] for t in result['tests']] != case['pilot_durations']:
            raise RuntimeError('Pilot duration coverage failed')
        checkpoint = n.ROOT/case['checkpoint']
        if n.sha(checkpoint) != case['checkpoint_sha256']:raise RuntimeError('Pilot checkpoint changed')
        with np.load(checkpoint) as z: initial = audit.checkpoint_fields(z, case)
        with np.load(path) as z:
            for row in result['tests']:
                prefix = row['trace_prefix']
                audit.audit_trial(row, {k[len(prefix):]:z[k] for k in z.files if k.startswith(prefix)}, case, initial)
        if not r.preparation_valid(case, result['tests']):raise RuntimeError('Pilot no-cue control failed')
        summary = dict(model=case['model'], length=case['length'], factor=case['gradient_factor'], slope=case['gradient_slope'])
        for gate in r.DEFINITIONS:
            successes = [t['duration'] for t in result['tests'] if t['duration']>0 and t[gate+'_success']]
            summary[gate+'_duration'] = min(successes) if successes else None
        table.append(summary)
        groups.setdefault((case['model'],case['gradient_factor']), []).append(summary)
    slopes = {}; selected = []
    for model in r.MODELS:
        candidates = []
        for factor in FACTORS:
            rows = groups[model,factor]
            counts = [sum(row[gate+'_duration'] is not None for row in rows) for gate in r.DEFINITIONS]
            candidates.append(dict(model=model, factor=factor, orientation_successes=counts[0], full_strength_successes=counts[1], lengths=[row['length'] for row in rows], slope=rows[0]['slope']))
        chosen = max(candidates, key=lambda c:(c['full_strength_successes'],c['orientation_successes'],-c['factor']))
        slopes[model] = chosen['slope']; selected.append(chosen)
    n.csvwrite(out/'results/pilot_readouts.csv', table)
    selection = dict(generated_at=n.now(), complete=True, audit_passed=True,
                     slope_by_model=slopes, selected=selected, selection_rule=manifest['selection_rule'],
                     pilot_manifest_sha256=n.sha(out/'manifest.json'), tested_factors=FACTORS)
    r.write(out/'results/slope_selection.json', selection)
    r.write(out/'results/status.json', dict(generated_at=n.now(), complete=True, producer='finished', trials=sum(len(c['pilot_durations']) for c in manifest['cases']), audit_passed=True))
    (out/'results/index.html').write_text(f'<h1>Fixed-slope calibration</h1><p>Generated {n.now()}; producer finished. All trial arrays passed independent gate audit.</p><a href="pilot_readouts.csv">All pilot readouts</a> · <a href="slope_selection.json">Selected slopes and rule</a>')
    if out.resolve()!=PILOT.resolve():raise RuntimeError('Automatic full launch requires canonical pilot destination')
    initialize(DEFAULT)
    r.launch_round(DEFAULT, 0)
    r.submit(DEFAULT, 'monitor', 0)


r.RUNNER = Path(__file__).resolve()
r.RENDERER = CODE/'gradient_reversal_plot.py'
r.DEFAULT = DEFAULT
r.response_cue = response_cue
r.source_hashes = sources
r.initialize = initialize


if __name__ == '__main__':
    if len(sys.argv)>1 and sys.argv[1].startswith('pilot-'):
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('action', choices=['pilot-init','pilot-launch','pilot-worker','pilot-finish'])
        parser.add_argument('--out', type=Path, default=PILOT)
        parser.add_argument('--slot', type=int)
        args = parser.parse_args(); out = args.out.resolve()
        try:
            if args.action=='pilot-init':pilot_init(out)
            elif args.action=='pilot-launch':pilot_launch(out)
            elif args.action=='pilot-worker':pilot_worker(out,args.slot)
            else:pilot_finish(out)
        except Exception:
            if args.action in ['pilot-worker','pilot-finish']:
                import traceback
                r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='stopped with error',error=traceback.format_exc()))
            raise
    else:r.main()
