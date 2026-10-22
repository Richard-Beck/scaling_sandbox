"""Three response geometries from one exact basal checkpoint per model/length.

A small, fixed-factor Goryachev pilot precedes the full duration comparison.
Unchanged trials can be reused only with identical checkpoint and stimulus.
"""
from pathlib import Path
import argparse, copy, gzip, json, os, shlex, shutil, subprocess, sys
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE / 'core'))
import corrected_reversal as r
import native as n
import polarity_establishment as e

DEFAULT = n.DATA / 'runs/three_test_reversal_20261006'
PILOT = n.DATA / 'runs/goryachev_fixed_amplitude_pilot_20261006'
PATCH = n.DATA / 'runs/corrected_reversal_20261003'
GRADIENT = n.DATA / 'runs/fixed_slope_reversal_20261006'
FACTORS = [4., 16., 64.]
TESTS = ['fixed_physical', 'fixed_relative', 'fixed_slope']


class Gradient:
    event_times = ()
    def __init__(self, slope): self.slope = float(slope)
    def __call__(self, time, x, length): return self.slope * np.asarray(x)


def response_cue(case):
    if case['condition'] == 'fixed_slope': return Gradient(case['gradient_slope'])
    return e.Patch(case['response_amplitude'], case['width'], side='right',
                   width_mode=case['width_mode'], fraction=.2,
                   background=case.get('background_stimulus', 0.))


def sources():
    paths = [Path(__file__), CODE/'three_test_audit.py', CODE/'audit_reversal.py',
             PATCH/'manifest.json', GRADIENT/'manifest.json',
             GRADIENT/'results/independent_audit.json', PATCH/'results/independent_audit.json', r.GRID]
    result = n.sources()
    result.update({str(p.relative_to(n.ROOT)): n.sha(p) for p in paths})
    return result


def freeze(out, manifest):
    for folder in ['work/chunks', 'work/jobs', 'work/logs', 'results', 'frozen_sources']:
        (out/folder).mkdir(parents=True, exist_ok=True)
    r.write(out/'manifest.json', manifest)
    hashes = sources()
    for name in hashes:
        dest = out/'frozen_sources'/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(n.ROOT/name, dest)
    r.write(out/'provenance.json', dict(generated_at=n.now(), sources=hashes,
                                       manifest_sha256=n.sha(out/'manifest.json')))


def base_cases():
    return [c for c in r.read(PATCH/'manifest.json')['cases'] if c['condition']=='relative_20pct']


def make_case(base, geometry, factor=1.):
    c = copy.deepcopy(base)
    c.update(preparation_condition=base['condition'], condition=geometry,
             source_case_index=base['case_index'], response_amplitude=base['amplitude']*factor,
             response_factor=factor)
    if geometry != 'fixed_slope':
        c.update(width_mode='physical' if geometry=='fixed_physical' else 'relative',
                 width=(.2 if c['model']=='spring' else 3.) if geometry=='fixed_physical' else .2*c['length'],
                 width_fraction_at_reference=.2)
    else:
        c['gradient_slope'] = r.read(GRADIENT/'manifest.json')['slope_by_model'][c['model']]
    return c


def pilot_init(out):
    if out.exists(): raise ValueError('Use a fresh pilot destination')
    eligible = sorted([c for c in base_cases() if c['model']=='goryachev' and c['eligible']], key=lambda c:c['length'])
    cases=[]
    for i in [0,len(eligible)//2,len(eligible)-1]:
        for geometry in TESTS[:2]:
            for factor in FACTORS:
                c=make_case(eligible[i],geometry,factor)
                c.update(case_index=len(cases),pilot_durations=[0.,5.12,51.2,512.,5120.])
                cases.append(c)
    freeze(out,dict(created_at=n.now(),experiment='constant_amplitude_goryachev_pilot',cases=cases,
                    factors=FACTORS,selection_rule='Lowest factor with both gates at all six geometry-length pairs; otherwise maximize strength successes, then orientation successes, then prefer the smaller factor.'))
    r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='not started'))


def pilot_worker(out,slot):
    c=r.verify(out)['cases'][slot]; rows=[];arrays={}
    for i,d in enumerate(c['pilot_durations']):
        row,z=r.trial(c,d);row['trace_prefix']=f'trial{i}_';rows.append(row)
        arrays.update({row['trace_prefix']+k:v for k,v in z.items()})
    path=out/'work/chunks'/f'pilot_{slot:03d}.npz'
    np.savez_compressed(path,**arrays)
    r.write(path.with_suffix('.json'),dict(case=c,tests=rows,trace_sha256=n.sha(path)))


def pilot_finish(out):
    import audit_reversal as a
    m=r.verify(out); summaries=[]
    for c in m['cases']:
        p=out/'work/chunks'/f'pilot_{c["case_index"]:03d}.npz';result=r.read(p.with_suffix('.json'))
        if result['case']!=c or n.sha(p)!=result['trace_sha256']: raise RuntimeError('Pilot identity changed')
        with np.load(n.ROOT/c['checkpoint']) as z: initial=a.checkpoint_fields(z,c)
        with np.load(p) as z:
            for row in result['tests']:
                prefix=row['trace_prefix']
                a.audit_trial(row,{k[len(prefix):]:z[k] for k in z.files if k.startswith(prefix)},c,initial)
        if not r.preparation_valid(c,result['tests']): raise RuntimeError('Invalid no-cue control')
        summaries.append(dict(length=c['length'],geometry=c['condition'],factor=c['response_factor'],
                              **{gate: min([t['duration'] for t in result['tests'] if t['duration']>0 and t[gate+'_success']],default=None) for gate in r.DEFINITIONS}))
    counts=[dict(factor=f,**{gate:sum(t[gate] is not None for t in summaries if t['factor']==f) for gate in r.DEFINITIONS}) for f in FACTORS]
    chosen=max(counts,key=lambda t:(t['full_polarity'],t['orientation'],-t['factor']))
    selection=dict(generated_at=n.now(),audit_passed=True,selected_factor=chosen['factor'],
                   amplitude=.1*chosen['factor'],counts=counts,selection_rule=m['selection_rule'])
    n.csvwrite(out/'results/pilot_readouts.csv',summaries)
    r.write(out/'results/selection.json',selection)
    r.write(out/'results/status.json',dict(generated_at=n.now(),complete=True,producer='finished',audit_passed=True,trials=90,selection=chosen))
    print(json.dumps(selection,indent=2))
    if out.resolve()==PILOT.resolve():
        initialize(DEFAULT)
        r.launch_round(DEFAULT,0)
        r.submit(DEFAULT,'monitor',0)


def import_trials(out, case, source, old, chosen, jobs):
    """Keep native archives unchanged; each imported chunk has an explicit source."""
    groups={}
    for t in chosen: groups.setdefault(t['trace'],[]).append(t)
    for trace,ts in groups.items():
        slot=len(jobs); dest=out/'work/chunks'/f'r0_{slot:05d}.npz'; origin=source/trace
        try: os.link(origin,dest)
        except OSError: shutil.copy2(origin,dest)
        rows=copy.deepcopy(ts)
        for row in rows:row['trace']=str(dest.relative_to(out))
        jobs.append(dict(slot=slot,case_index=case['case_index'],durations=[t['duration'] for t in rows]))
        r.write(dest.with_suffix('.json'),dict(generated_at=n.now(),round=0,slot=slot,case_index=case['case_index'],
                checkpoint_sha256=case['checkpoint_sha256'],trials=rows,trace_sha256=n.sha(dest),
                reuse_source=dict(campaign=str(source.relative_to(n.ROOT)),case_index=old['case_index'],
                                  trace=trace,archive_sha256=n.sha(origin))))


def initialize(out):
    if out.exists():raise ValueError('Use a fresh campaign destination')
    selection=r.read(PILOT/'results/selection.json')
    if not selection['audit_passed']:raise RuntimeError('Pilot not audited')
    cases=[]
    for base in base_cases():
        for geometry in TESTS:
            factor=selection['selected_factor'] if base['model']=='goryachev' and geometry!='fixed_slope' else 1.
            c=make_case(base,geometry,factor);c['case_index']=len(cases);cases.append(c)
    manifest=dict(created_at=n.now(),models=r.MODELS,excluded_models=['otsuji','legi'],cases=cases,
                  conditions=TESTS,experiment='three_geometry_shared_preparation',amplitude_factor=1.,
                  preparation='Exact nominal relative-20% corrected basal-release checkpoint shared by all three response geometries.',
                  goryachev_selection=selection,recovery_horizon=r.RECOVERY,sample_interval=r.STEP,gate_window=r.TAIL,
                  refinement=dict(relative_bracket_tolerance=.01,maximum_rounds=r.MAX_ROUNDS),
                  cue='Fixed width 3 (Spring native 0.2); relative 20%; finite g*x at retained constant slope. Patch response amplitude is constant across lengths; only Goryachev uses the audited fixed factor.')
    freeze(out,manifest)
    records={}
    for source in [PATCH,GRADIENT]:
        with gzip.open(source/'results/records.json.gz','rt') as f:
            records[source]={ (c['model'],c['length'],c['condition']):c for c in json.load(f)['cases'] }
    jobs=[]
    for c in cases:
        if not c['eligible']:continue
        source=GRADIENT if c['condition']=='fixed_slope' else PATCH
        old=records[source][c['model'],c['length'],'relative_20pct']
        if old['checkpoint_sha256']!=c['checkpoint_sha256']:raise RuntimeError('Reuse checkpoint differs')
        reuse=c['condition']=='fixed_slope' or (c['condition']=='fixed_relative' and c['response_factor']==1.)
        c['reused_positive_trials']=reuse
        import_trials(out,c,source,old,old['tests'] if reuse else [t for t in old['tests'] if t['duration']==0.],jobs)
    # Metadata only: scientific cases are refrozen before any worker launches.
    freeze(out,manifest)
    r.save_jobs(out,0,jobs)
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='initialized',round=0,arrays=[]))
    r.publish(out,'initialized')


def advance(out,round):
    manifest,groups,pending,errors=r.collect(out,partial=True)
    if pending or errors:raise RuntimeError(f'Incomplete execution: {len(pending)} pending, {len(errors)} errors')
    new=[]
    for c in manifest['cases']:
        if not r.preparation_valid(c,groups[c['case_index']]):continue
        if round==0:
            existing={t['duration'] for t in groups[c['case_index']]}
            ds=[d for d in c['durations'] if d not in existing]
        else: ds=r.refinement(groups[c['case_index']])
        for start in range(0,len(ds),10):new.append(dict(slot=len(new),case_index=c['case_index'],durations=ds[start:start+10]))
    if new and round<=r.MAX_ROUNDS:
        r.save_jobs(out,round+1,new);return r.launch_round(out,round+1)
    if new:raise RuntimeError('Refinement cap exhausted')
    rows,cases=r.rows_and_records(manifest,groups,True)
    unresolved=sum(q['duration_boundary_unresolved'] for c in cases for q in c['result'].values())
    if unresolved:raise RuntimeError('Unresolved selected boundaries')
    audit=out/'results/independent_audit.json'
    subprocess.run([sys.executable,str(CODE/'three_test_audit.py'),'--campaign',str(out),'--before-publication','--output',str(audit)],check=True)
    result=r.read(audit)
    if not result['passed'] or not result['complete']:raise RuntimeError('Independent audit failed')
    lookup={(t['case_index'],t['definition']):t for t in result['selected_readouts']}
    for row in rows:
        expected=lookup[row['case_index'],row['definition']]
        for key in ['status','endpoint','selected_tail_field_pass_0_1_percent']:
            if row[key]!=expected[key]:raise RuntimeError('Independent selection differs')
    r.write(out/'results/validation.json',dict(generated_at=n.now(),passed=True,complete_execution=True,
            independent_audit_passed=True,independent_audit_sha256=n.sha(audit),unresolved_readouts=0,
            refinement_rounds=round,shared_exact_checkpoint=True))
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='complete',round=round,arrays=[]))
    r.publish(out,'complete',complete=True)
    subprocess.run([sys.executable,str(CODE/'three_test_figures.py'),'--campaign',str(out)],check=True)
    subprocess.run([sys.executable,str(CODE/'mechanisms/three_test_mechanisms.py'),'init','--campaign',str(out)],check=True)


r.RUNNER=Path(__file__).resolve()
r.DEFAULT=DEFAULT
r.response_cue=response_cue
r.source_hashes=sources
r.initialize=initialize
r.advance=advance
original_launch_round=r.launch_round
def launch_round(out,round):
    if round!=0:return original_launch_round(out,round)
    r.verify(out)
    controller=r.submit(out,'advance',0)
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='no_cue_controls_running',round=0,arrays=[],controller=controller))
    r.publish(out,'no_cue_controls_running')
r.launch_round=launch_round


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1].startswith('pilot-'):
        p=argparse.ArgumentParser();p.add_argument('action');p.add_argument('--out',type=Path,default=PILOT);p.add_argument('--slot',type=int)
        p.add_argument('--round',type=int,default=0);p.add_argument('--offset',type=int,default=0)
        args=p.parse_args();out=args.out.resolve()
        if args.action=='pilot-init':pilot_init(out)
        elif args.action=='pilot-worker':pilot_worker(out,args.slot)
        elif args.action=='pilot-finish':pilot_finish(out)
        elif args.action=='pilot-launch':
            job=r.submit(out,'pilot-worker',0,array=f'0-{len(r.verify(out)["cases"])-1}')
            final=r.submit(out,'pilot-finish',0,dependency=job)
            r.write(out/'launch.json',dict(generated_at=n.now(),array=job,controller=final))
            r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='running',array=job,controller=final))
            print(json.dumps(dict(array=job,controller=final)))
    else:r.main()
