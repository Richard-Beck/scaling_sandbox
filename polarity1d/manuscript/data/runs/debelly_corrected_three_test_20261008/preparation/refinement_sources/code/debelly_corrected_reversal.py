"""Corrected independent-endpoint preparation and three-test replacement campaign."""
from pathlib import Path
from dataclasses import replace
import argparse,copy,gzip,json,os,shutil,subprocess,sys
import numpy as np
CODE=Path(__file__).resolve().parent
sys.path.insert(0,str(CODE/'core'))
import native as n
import corrected_reversal as r
import polarity_establishment as e
DEFAULT=n.DATA/'runs/debelly_corrected_three_test_20261008'
BASE=n.DATA/'runs/three_test_reversal_20261006'
TESTS=['fixed_physical','fixed_relative','fixed_slope']
_original_model=n.model_for

def model_for(key):
    m=_original_model(key)
    return n.ds.DeBelly(replace(m.parameters,center_boundary_motion=False)) if key=='debelly' else m
n.model_for=model_for

class Gradient:
    event_times=()
    def __init__(self,g):self.g=float(g)
    def __call__(self,t,x,L):return self.g*np.asarray(x)

def response_cue(c):
    if c['condition']=='fixed_slope':return Gradient(c['gradient_slope'])
    return e.Patch(c['response_amplitude'],c['width'],side='right',width_mode=c['width_mode'],fraction=.2)

def sources():
    result=n.sources()
    for p in [Path(__file__),CODE/'debelly_corrected_reversal_audit.py',CODE/'audit_reversal.py',BASE/'manifest.json',BASE/'results/independent_audit.json',r.GRID,n.DATA/'establishment_basal/manifest.json']:
        result[str(p.relative_to(n.ROOT))]=n.sha(p)
    return result

def freeze(out):
    hashes=sources()
    for name in hashes:
        dest=out/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(n.ROOT/name,dest)
    return hashes

def prep_init(out):
    if out.exists():raise ValueError('Fresh destination required')
    for f in ['preparation/work/cases','work/chunks','work/jobs','work/logs','results','frozen_sources']:(out/f).mkdir(parents=True,exist_ok=True)
    m=e.make_manifest();cases=[c for c in m['cases'] if c['model']=='debelly' and c['condition']=='relative_20pct' and c['amplitude_factor']==1]
    for i,c in enumerate(cases):c['original_establishment_case_index']=c['case_index'];c['case_index']=i
    m.update(cases=cases,boundary_method='independent_endpoints',center_boundary_motion=False)
    r.write(out/'preparation/manifest.json',m)
    r.write(out/'preparation/provenance.json',dict(sources=freeze(out),manifest_sha256=n.sha(out/'preparation/manifest.json')))
    r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='not started',stage='fresh_corrected_preparation',lengths=len(cases)))

def prep_verify(out):
    p=r.read(out/'preparation/provenance.json')
    if sources()!=p['sources'] or n.sha(out/'preparation/manifest.json')!=p['manifest_sha256']:raise RuntimeError('Frozen preparation inputs changed')
    return r.read(out/'preparation/manifest.json')

def prep_worker(out,slot):
    prep_verify(out)
    e.sources=sources
    e.worker(out/'preparation',slot)

def initialize(out):
    plan=prep_verify(out);original=r.read(BASE/'manifest.json');cases=copy.deepcopy(original['cases'])
    bylength={c['length']:c for c in plan['cases']}
    for c in cases:
        if c['model']!='debelly':continue
        pc=bylength[c['length']];p=out/'preparation/work/cases'/f'case_{pc["case_index"]:04d}/result.json'
        if not p.exists():raise RuntimeError('Missing corrected preparation '+str(p))
        result=r.read(p);s=result['summary']
        if s['status']=='unsettled' or s['duration_boundary_unresolved'] or s['upper_optimum_unresolved']:raise RuntimeError('Unresolved preparation '+str(s))
        c.update(eligible=False,checkpoint=None,checkpoint_sha256=None,initial_orientation=None,initial_contrast=None,preparation_status=s['status'],boundary_method='independent_endpoints',center_boundary_motion=False,corrected_preparation=str(p.relative_to(n.ROOT)))
        if s['status']=='polarized':
            cp=out/'preparation'/s['selected_checkpoint'];state,tr=n.read_checkpoint(cp,'debelly')
            c.update(checkpoint=str(cp.relative_to(n.ROOT)),checkpoint_sha256=n.sha(cp),initial_orientation=float(tr.orientation[-1]),initial_contrast=float(tr.contrast[-1]),eligible=True)
        c.pop('reused_positive_trials',None)
    manifest=copy.deepcopy(original);manifest.update(created_at=n.now(),cases=cases,experiment='three_geometry_corrected_debelly',debelly_boundary_method='independent_endpoints',center_boundary_motion=False,replaced_model='debelly',previous_campaign=str(BASE.relative_to(n.ROOT)),preparation='All three tests share own-model nominal relative20 basal checkpoint. deBelly independently re-established from rest with local endpoint commands; unchanged four models retain exact audited evidence.')
    r.write(out/'manifest.json',manifest)
    r.write(out/'provenance.json',dict(generated_at=n.now(),sources=freeze(out),manifest_sha256=n.sha(out/'manifest.json')))
    with gzip.open(BASE/'results/records.json.gz','rt') as f:old={c['case_index']:c for c in json.load(f)['cases']}
    from three_test_reversal import import_trials
    hooks()
    # Import only immutable unchanged four-model trials; never any deBelly state/outcome.
    jobs=[]
    for c in cases:
        if not c['eligible']:continue
        if c['model']=='debelly':jobs.append(dict(slot=len(jobs),case_index=c['case_index'],durations=[0.]))
        else:import_trials(out,c,BASE,old[c['case_index']],old[c['case_index']]['tests'],jobs)
    r.save_jobs(out,0,jobs)
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='initialized',round=0,arrays=[]))
    r.publish(out,'initialized')

_original_launch=r.launch_round

def launch_round(out,round):
    # Cached chunks are already present; submit only missing slots.
    r.verify(out);planned=r.jobs(out,round)
    slots=[j['slot'] for j in planned if not (out/'work/chunks'/f'r{round}_{j["slot"]:05d}.json').exists()]
    if not slots:return advance(out,round)
    arrays=[r.submit(out,'worker',round,array=e.array_spec(slots))]
    controller=r.submit(out,'advance',round,dependency=':'.join(arrays))
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='duration_grid_running' if round else 'no_cue_controls_running',round=round,arrays=arrays,controller=controller))
    r.publish(out,'duration_grid_running' if round else 'no_cue_controls_running')

def advance(out,round):
    m,groups,pending,errors=r.collect(out,partial=True)
    if pending or errors:raise RuntimeError(f'Incomplete: {len(pending)} pending, {len(errors)} errors')
    new=[]
    for c in m['cases']:
        if c['model']!='debelly' or not r.preparation_valid(c,groups[c['case_index']]):continue
        ds=[d for d in c['durations'] if d not in {t['duration'] for t in groups[c['case_index']]}] if round==0 else r.refinement(groups[c['case_index']])
        for start in range(0,len(ds),10):new.append(dict(slot=len(new),case_index=c['case_index'],durations=ds[start:start+10]))
    if new:
        if round>=r.MAX_ROUNDS:raise RuntimeError('Refinement cap exhausted')
        r.save_jobs(out,round+1,new);return launch_round(out,round+1)
    rows,cases=r.rows_and_records(m,groups,True)
    if any(q['duration_boundary_unresolved'] or q['unresolved_transition_count'] for c in cases for q in c['result'].values()):raise RuntimeError('Unresolved sampled boundaries')
    audit=out/'results/independent_audit.json'
    subprocess.run([sys.executable,str(CODE/'debelly_corrected_reversal_audit.py'),'--campaign',str(out),'--before-publication','--output',str(audit)],check=True)
    result=r.read(audit)
    if not result['passed'] or not result['complete']:raise RuntimeError('Independent audit failed')
    lookup={(a['case_index'],a['definition']):a for a in result['selected_readouts']}
    for row in rows:
        for key in ['status','endpoint','selected_tail_field_pass_0_1_percent']:
            if row[key]!=lookup[row['case_index'],row['definition']][key]:raise RuntimeError('Independent selection differs')
    r.write(out/'results/validation.json',dict(generated_at=n.now(),passed=True,complete_execution=True,independent_audit_passed=True,independent_audit_sha256=n.sha(audit),refinement_rounds=round))
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='complete',round=round,arrays=[]));r.publish(out,'complete',True)
    subprocess.run([sys.executable,str(CODE/'three_test_figures.py'),'--campaign',str(out)],check=True)

# Importing three_test_reversal above resets hooks; explicitly reset before dispatch.
def hooks():
    r.RUNNER=Path(__file__).resolve();r.DEFAULT=DEFAULT;r.response_cue=response_cue;r.source_hashes=sources;r.advance=advance;r.launch_round=launch_round;r.initialize=initialize
hooks()
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1].startswith('prep-'):
        p=argparse.ArgumentParser();p.add_argument('action');p.add_argument('--out',type=Path,default=DEFAULT);p.add_argument('--round',type=int,default=0);p.add_argument('--offset',type=int,default=0);p.add_argument('--slot',type=int)
        a=p.parse_args();out=a.out.resolve()
        if a.action=='prep-init':prep_init(out)
        elif a.action=='prep-worker':prep_worker(out,a.slot)
        elif a.action=='prep-finish':initialize(out);hooks();launch_round(out,0)
        elif a.action=='prep-launch':
            plan=prep_verify(out);array=r.submit(out,'prep-worker',0,array=f'0-{len(plan["cases"])-1}');final=r.submit(out,'prep-finish',0,dependency=array)
            r.write(out/'preparation_launch.json',dict(array=array,finalizer=final,generated_at=n.now()));print(array,final)
    else:r.main()
