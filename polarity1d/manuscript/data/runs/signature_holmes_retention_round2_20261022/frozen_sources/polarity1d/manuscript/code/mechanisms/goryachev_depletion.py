"""Factorial substrate-availability diagnostics from exact native checkpoints.

During the pulse only, equalize cue and/or challenger feedback substrate to
the initial right-end concentration at L=9. The incumbent retains native
feedback. Opposite reaction transfers preserve mass; a regularized kinetic
clamp vanishes at zero actual substrate. Native recovery is always restored.
"""
from pathlib import Path
from dataclasses import dataclass
import argparse, copy, gzip, json, shlex, shutil, subprocess, sys, traceback
import numpy as np

CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE))
import three_test_reversal as t
import audit_reversal as independent
sys.path.insert(0, str(CODE/'mechanisms'))
import three_test_mechanisms as spatial
import three_test_analysis as analysis
n, r, e = t.n, t.r, t.e
from models.goryachev import Goryachev

OUT = n.DATA/'runs/goryachev_depletion_20261008_v2'
LENGTHS = [9., 15., 21., 27., 33., 39.]
KINDS = ['native', 'cue_clamp', 'feedback_clamp', 'both_clamp', 'both_initial']


@dataclass(frozen=True)
class AvailabilityControl(Goryachev):
    kind: str = 'native'
    length: float = 9.
    reference_v: float = .1
    initial_v: float = .1
    epsilon_fraction: float = .001

    def effective(self, v):
        # Equals reference_v at v=reference_v, approaches reference_v+epsilon,
        # and vanishes at v=0. Positivity and local mass transfer are retained.
        eps = self.reference_v*self.epsilon_fraction
        return (self.reference_v+eps)*v/(v+eps)

    def challenger_weight(self, x):
        # Native left half through q=.4; smooth transition; full right half.
        q = np.clip((np.asarray(x)/self.length-.4)/.1, 0., 1.)
        return q*q*(3.-2.*q)

    def reaction(self, time, x, state, stimulus):
        del time
        u, v = state
        cv, fv = v, v
        if self.kind in ['cue_clamp', 'both_clamp']:
            cv = self.effective(v)
        if self.kind in ['feedback_clamp', 'both_clamp']:
            fv = v+self.challenger_weight(x)*(self.effective(v)-v)
        if self.kind == 'both_initial':
            factor = self.reference_v/self.initial_v
            cv = factor*v
            fv = v*(1.+self.challenger_weight(x)*(factor-1.))
        transfer = self.a*self.ec*u*u*fv+self.b*self.ec*u*fv-self.c*u
        transfer += stimulus*cv
        return np.vstack((transfer, -transfer))


def submit(out, action, array=None, dependency=None):
    cmd = 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join(
        [sys.executable, str(Path(__file__)), action, '--out', str(out)])
    if array: cmd += ' --slot "$SLURM_ARRAY_TASK_ID"'
    opts = ['sbatch', '--parsable', '--qos=small', '--ntasks=1', '--cpus-per-task=1',
            '--mem=3G', '--time=02:00:00', '--job-name=gory-depletion-'+action,
            '--chdir='+str(n.ROOT), '--output='+str(out/'logs'/f'{action}_%A_%a.out'),
            '--error='+str(out/'logs'/f'{action}_%A_%a.err')]
    if array: opts += ['--array='+array]
    if dependency: opts += ['--dependency=afterany:'+dependency]
    return subprocess.check_output(opts+['--wrap', cmd], text=True).strip().split(';')[0]


def verify(out):
    p = r.read(out/'provenance.json')
    assert n.sha(out/'manifest.json') == p['manifest_sha256']
    for name, digest in p['sources'].items():
        assert n.sha(n.ROOT/name) == digest, 'Source changed: '+name
        assert n.sha(out/'frozen_sources'/name) == digest
    return r.read(out/'manifest.json')


def initialize(out):
    if out.exists(): raise ValueError('Use a fresh output directory')
    campaign = t.DEFAULT
    assert r.read(campaign/'results/status.json')['complete']
    with gzip.open(campaign/'results/records.json.gz', 'rt') as stream:
        records = json.load(stream)['cases']
    cs = [c for c in records if c['model']=='goryachev' and c['length'] in LENGTHS]
    ref = next(c for c in cs if c['length']==9. and c['condition']=='fixed_relative')
    state, _ = n.read_checkpoint(n.ROOT/ref['checkpoint'], 'goryachev')
    vref = float(state[1,-1])
    jobs = []
    for c0 in cs:
        assert c0['eligible']
        c = {k:copy.deepcopy(v) for k,v in c0.items() if k not in ['tests','result']}
        s, _ = n.read_checkpoint(n.ROOT/c['checkpoint'], 'goryachev')
        for kind in KINDS:
            jobs.append(dict(slot=len(jobs),case=c,kind=kind,reference_v=vref,
                             initial_v=float(s[1,-1]),epsilon_fraction=.001,
                             native_result=c0['result']))
    for folder in ['logs','workers','results','frozen_sources']:
        (out/folder).mkdir(parents=True, exist_ok=True)
    r.write(out/'manifest.json', dict(created_at=n.now(),jobs=jobs,lengths=LENGTHS,
        reference_length=9.,reference_inactive=vref,kinds=KINDS,geometries=t.TESTS,
        intervention='Cue only; native recovery. Kinetic clamp for cue term everywhere cue acts; feedback clamp only in challenger half (smooth q=.4-.5 transition). Incumbent feedback stays native. Transfers remain equal and opposite. This does not add actual inactive molecules.',
        both_initial='Constant coefficient multiplier reference initial right-end v / own initial right-end v for cue and challenger feedback; retains all subsequent substrate variation.',
        search='17 logarithmic doses .05-5120 plus native endpoints and their halves; refine every sampled fail-to-pass boundary for both gates to 1%, at most 16 rounds. No monotonicity assumption; unsampled narrow intervals remain possible.',
        controls='Matched no-cue intervention at cap and each selected endpoint. Native zero control and native endpoint replays. Same exact initial state for every dose.',
        recovery=10000.,sample=2.5,gate_tail=500.))
    paths = [Path(__file__), CODE/'three_test_reversal.py', CODE/'audit_reversal.py',
             CODE/'mechanisms/three_test_mechanisms.py', CODE/'mechanisms/three_test_analysis.py',
             campaign/'manifest.json', campaign/'results/endpoints.csv']
    hashes = n.sources(); hashes.update({str(p.relative_to(n.ROOT)):n.sha(p) for p in paths})
    for name in hashes:
        dest = out/'frozen_sources'/name; dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(n.ROOT/name,dest)
    r.write(out/'provenance.json',dict(sources=hashes,manifest_sha256=n.sha(out/'manifest.json')))
    r.write(out/'results/status.json',dict(generated_at=n.now(),producer='running',complete=False,jobs=len(jobs)))
    array = submit(out,'worker',f'0-{len(jobs)-1}')
    final = submit(out,'finish',dependency=array)
    r.write(out/'launch.json',dict(generated_at=n.now(),array=array,finalizer=final))
    print(json.dumps(dict(jobs=len(jobs),array=array,finalizer=final,reference_v=vref)))


def run_trial(job, D, no_cue=False, full=False):
    c = job['case']; L = c['length']; cp = n.ROOT/c['checkpoint']
    assert n.sha(cp) == c['checkpoint_sha256']
    state, prep = n.read_checkpoint(cp,'goryachev')
    original = state.copy()
    model = AvailabilityControl(kind=job['kind'],length=L,reference_v=job['reference_v'],
                                initial_v=job['initial_v'],epsilon_fraction=job['epsilon_fraction'])
    pulse = None
    if D>0:
        tr = n._sim(model,state,L,c['cells'],n.ConstantInput(0.) if no_cue else t.response_cue(c),
                    D,interval=min(2.5,D/20))
        pulse = tr
        state = np.vstack([tr.fields['active'][-1],tr.fields['inactive'][-1]])
    recovery = n.run('goryachev',L,state,e.withdrawal(c),10000.,interval=2.5,n=c['cells'])
    # Use standard pulse adapter for recorded scalar/field contract.
    if pulse:
        pulse.active = pulse.fields['active']; pulse.orientation = n._orientation(pulse)
        pulse.contrast = np.ptp(pulse.active,axis=1)
    row = spatial.gates(recovery,c)
    row.update(duration=float(D),no_cue=bool(no_cue),no_cue_maximum_scaled_excursion=None,
               no_cue_left_tail_valid=None,no_cue_preparation_valid=None)
    z = dict(x=recovery.x,originalstate=original,
             pulse_time=pulse.time if pulse else np.array([0.]),
             pulse_orientation=pulse.orientation if pulse else np.array([c['initial_orientation']]),
             pulse_contrast=pulse.contrast if pulse else np.array([c['initial_contrast']]),
             recovery_time=recovery.time,recovery_orientation=recovery.orientation,
             recovery_contrast=recovery.contrast,finalstate=recovery.state)
    tail = recovery.time>=9500.
    row['field_tail_metrics'] = {}
    for name,y in recovery.fields.items():
        a = y[tail]
        scale = c['fixed_scales'][name]
        row['field_tail_metrics'][name] = dict(range_scaled=float(np.ptp(a,axis=0).max())/scale,
                 projected_drift_500=float(abs(a[-1]-a[-2]).max())/2.5*500/scale)
        for suffix,value in [('min',a.min(axis=0)),('max',a.max(axis=0)),('penultimate',a[-2]),('last',a[-1])]:
            z['fieldtail_'+name+'_'+suffix] = value
    if D==0:
        excursion=0.
        for name,y in recovery.fields.items():
            initial=prep.fields[name][-1]
            excursion=max(excursion,float(np.max(abs(y-initial)))/c['fixed_scales'][name])
            for suffix,value in [('min',y.min(axis=0)),('max',y.max(axis=0)),('initial',initial)]:
                z['no_cue_'+name+'_'+suffix] = value
        o=recovery.orientation[tail];q=recovery.contrast[tail]
        left=bool(o.max()<-.05 and np.ptp(o)<.005 and q.min()>=.05)
        row.update(no_cue_maximum_scaled_excursion=excursion,no_cue_left_tail_valid=left,
                   no_cue_preparation_valid=bool(excursion<=.01 and left and not row['orientation_success']))
    times=[];left=[];right=[];vright=[];vmin=[];vmean=[]
    initialmass=None;masserror=0.;minimum=0.
    for phase,tr in [('pulse',pulse),('recovery',recovery)]:
        if tr is None:continue
        u=tr.fields['active'];v=tr.fields['inactive'];q=(np.arange(u.shape[1])+.5)/u.shape[1]
        background=np.median(u[:,(q>=.475)&(q<=.525)],axis=1)
        excess=np.maximum(u-background[:,None],0.)
        lm=excess[:,q<.5].sum(axis=1)*L/u.shape[1]
        rm=excess[:,q>.5].sum(axis=1)*L/u.shape[1]
        tm=tr.time if phase=='pulse' else D+tr.time
        times.extend(tm);left.extend(lm);right.extend(rm);vright.extend(v[:,-1]);vmin.extend(v.min(axis=1));vmean.extend(v.mean(axis=1))
        mass=(u+v).sum(axis=1)
        masserror=max(masserror,float(np.max(abs(mass/original.sum()-1.))))
        minimum=min(minimum,float(min(u.min(),v.min())))
        z[phase+'_mass']=mass
        if full:
            z[phase+'_active']=u;z[phase+'_field_active']=u;z[phase+'_field_inactive']=v
            z[phase+'_x']=tr.x
        if initialmass is None:initialmass=lm[0]
    times=np.asarray(times);left=np.asarray(left)
    half=np.flatnonzero(left<=.5*initialmass)
    row.update(old_pole_half_removal_time=float(times[half[0]]) if len(half) else None,
               conservation_error=masserror,minimum_concentration=minimum,
               no_cue_stays_left=bool(recovery.orientation[tail].max()<-.05))
    z.update(history_time=times,old_pole_excess=left,new_pole_excess=np.asarray(right),
             inactive_right=np.asarray(vright),inactive_min=np.asarray(vmin),inactive_mean=np.asarray(vmean))
    independent.audit_trial(row,z,c,{name:prep.fields[name][-1] for name in ['active','inactive']})
    assert masserror<1e-6 and minimum>=-1e-7
    return row,z


def worker(out,slot):
    m=verify(out);job=m['jobs'][slot];directory=out/'workers'/f'{slot:03d}'
    directory.mkdir(exist_ok=True)
    rows=[]
    def trial(D,no_cue=False,full=False,label=None):
        name=label or ('nocue_' if no_cue else 'dose_')+format(D,'.12g').replace('.','p')
        path=directory/(name+'.npz')
        if path.with_suffix('.json').exists():return r.read(path.with_suffix('.json'))['row']
        row,z=run_trial(job,D,no_cue,full)
        np.savez_compressed(path,**z)
        row['archive']=str(path.relative_to(out));row['archive_sha256']=n.sha(path)
        r.write(path.with_suffix('.json'),dict(job=job,row=row,generated_at=n.now()))
        return row
    zero=trial(0.);assert zero['no_cue_preparation_valid']
    if job['kind']=='native':
        doses=sorted(set(v['endpoint'] for v in job['native_result'].values() if v['endpoint'] is not None))
        rows=[zero]+[trial(float(d),full=True) for d in doses]
        for gate,res in job['native_result'].items():
            if res['endpoint'] is not None:assert next(x for x in rows if x['duration']==res['endpoint'])[gate+'_success']
    else:
        doses=set(np.geomspace(.05,5120.,17).tolist())
        for res in job['native_result'].values():
            if res['endpoint'] is not None:doses.update([res['endpoint'],res['endpoint']/2.])
        rows=[zero]+[trial(float(d)) for d in sorted(doses)]
        for iteration in range(16):
            rows.sort(key=lambda x:x['duration']);new=set()
            for gate in r.DEFINITIONS:
                for a,b in zip(rows,rows[1:]):
                    if not a[gate+'_success'] and b[gate+'_success'] and (b['duration']-a['duration'])/b['duration']>.01:
                        new.add((a['duration']+b['duration'])/2.)
            if not new:break
            rows += [trial(float(d)) for d in sorted(new)]
        else:raise RuntimeError('Unresolved sampled boundaries')
    selected={}
    controls=[]
    for gate in r.DEFINITIONS:
        passing=[x for x in rows if x['duration']>0 and x[gate+'_success']]
        if passing:
            best=min(passing,key=lambda x:x['duration'])
            lower=max(x['duration'] for x in rows if x['duration']<best['duration'] and not x[gate+'_success'])
            selected[gate]=dict(status='success',endpoint=best['duration'],lower=lower,
                                 relative_bracket=(best['duration']-lower)/best['duration'],
                                 field_tail_pass=best['field_tail_pass'],
                                 final_orientation=best['final_orientation'],
                                 old_pole_half_removal_time=best['old_pole_half_removal_time'])
            if job['kind']!='native':
                trial(best['duration'],full=True,label='selected_'+gate)
                controls.append(trial(best['duration'],no_cue=True,full=True))
        else:selected[gate]=dict(status='censored',endpoint=None,lower=None,relative_bracket=None,
                                  field_tail_pass=None,final_orientation=None,old_pole_half_removal_time=None)
    if job['kind']!='native':controls.append(trial(5120.,no_cue=True,full=True))
    r.write(directory/'summary.json',dict(generated_at=n.now(),job=job,trials=rows,selected=selected,
              no_cue_controls=controls,no_cue_valid=all(x['no_cue_stays_left'] and not x['orientation_success'] for x in controls)))


def finish(out):
    m=verify(out);summaries=[];audit=[];endpoints=[]
    for job in m['jobs']:
        directory=out/'workers'/f'{job["slot"]:03d}'
        summary=r.read(directory/'summary.json');assert summary['job']==job
        assert summary['no_cue_valid'], 'Intervention alone reverses'
        c=job['case'];cp=n.ROOT/c['checkpoint'];assert n.sha(cp)==c['checkpoint_sha256']
        _,prep=n.read_checkpoint(cp,'goryachev')
        for p in sorted(directory.glob('*.json')):
            if p.name=='summary.json':continue
            record=r.read(p);row=record['row'];assert record['job']==job
            zpath=p.with_suffix('.npz');assert n.sha(zpath)==row['archive_sha256']
            with np.load(zpath) as z:
                independent.audit_trial(row,z,c,{name:prep.fields[name][-1] for name in ['active','inactive']})
                np.testing.assert_array_equal(z['originalstate'],np.vstack([prep.fields['active'][-1],prep.fields['inactive'][-1]]))
                for phase in ['pulse','recovery']:
                    if phase+'_mass' in z:assert np.max(abs(z[phase+'_mass']/z['originalstate'].sum()-1.))<1e-6
                if 'pulse_active' in z:
                    u=z['pulse_active'];v=z['pulse_field_inactive']
                    np.testing.assert_allclose(n._orientation(type('Trace',(),dict(fields={'active':u},x=z['pulse_x']))()),z['pulse_orientation'],atol=1e-10)
                    assert min(u.min(),v.min())>=-1e-7
            audit.append(dict(slot=job['slot'],archive=str(zpath.relative_to(out)),passed=True))
        for gate,sel in summary['selected'].items():
            # Independent selection and bracket check from audited trial rows.
            passing=sorted(x['duration'] for x in summary['trials'] if x['duration']>0 and x[gate+'_success'])
            assert sel['endpoint']==(passing[0] if passing else None)
            if job['kind']!='native' and passing:assert sel['relative_bracket']<=.01
            endpoints.append(dict(length=c['length'],geometry=c['condition'],kind=job['kind'],gate=gate,
                                  reference_v=job['reference_v'],initial_v=job['initial_v'],**sel))
        summaries.append(summary)
    n.csvwrite(out/'results/endpoints.csv',endpoints)
    r.write(out/'results/independent_audit.json',dict(generated_at=n.now(),passed=True,archives=len(audit),checks=audit))
    r.write(out/'results/summaries.json',summaries)
    render(out,m,endpoints)
    r.write(out/'results/status.json',dict(generated_at=n.now(),producer='finished',complete=True,
             jobs=len(summaries),audited_archives=len(audit),audit_passed=True,no_cue_controls_passed=True))
    print(json.dumps(dict(jobs=len(summaries),archives=len(audit),complete=True)))


def render(out,m,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True)
    colors=['#222222','#0072B2','#D55E00','#009E73','#CC79A7']
    for i,gate in enumerate(r.DEFINITIONS):
        for j,geometry in enumerate(t.TESTS):
            ax=axes[i,j]
            for kind,color in zip(KINDS,colors):
                rs=sorted([x for x in rows if x['gate']==gate and x['geometry']==geometry and x['kind']==kind],key=lambda x:x['length'])
                ax.plot([x['length'] for x in rs],[x['endpoint'] if x['endpoint'] else np.nan for x in rs],'-',color=color,label=kind)
                for x in rs:
                    if x['endpoint']:ax.plot(x['length'],x['endpoint'],'o',color=color,mfc=color if x['field_tail_pass'] else 'white')
                    else:ax.plot(x['length'],5120.,'^',color=color,mfc='white')
            ax.set(title=geometry+' / '+gate,yscale='log',xlabel='Length',ylabel='Shortest passing tested cue duration')
            ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);fig.suptitle('Goryachev: cue versus challenger-feedback substrate compensation')
    fig.tight_layout();fig.savefig(out/'results/duration_comparison.png',dpi=170)
    fig.savefig(out/'results/duration_comparison.pdf');plt.close(fig)
    text = f'''# Goryachev depletion controls\n\nGenerated {n.now()}; producer finished.\n\nAll controls start exact native basal checkpoints. During the pulse only, cue and/or challenger feedback use a regularized constant effective substrate equal to initial right-end inactive at L=9 ({m['reference_inactive']:.8g}). Native feedback remains at the old pole; the transition occupies normalized positions .4-.5. Native recovery lasts 10,000 units. Equal/opposite reaction transfers preserve the actual total pool; this is kinetic compensation, not replenishment of actual inactive material. Both initial uses fixed initial concentration ratios instead, retaining dynamic depletion.\n\nThe two scalar gates and separate field-tail flag are unchanged. Search refines all sampled fail-to-pass transitions to 1%; it cannot exclude narrower unsampled success intervals. Native curves use audited retained thresholds with endpoint replays, rather than a fresh coarse search. Matched no-cue intervention controls at selected doses and the cap must remain left oriented. Open markers retain failed field flags; triangles denote censoring at 5,120.\n\n[Duration figure](duration_comparison.png), [endpoints](endpoints.csv), [independent audit](independent_audit.json).\n\nInterpretation limits: the feedback control targets the challenger, not incumbent sequestration. Effective-substrate clamping is regularized to vanish at zero actual inactive concentration and can cease to compensate when transport drives substrate close to zero. A residual length effect does not establish independence from substrate limitation unless clamp effectiveness is checked.\n'''
    (out/'results/README.md').write_text(text)
    (out/'results/index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Goryachev depletion tests</title></head><body><h1>Goryachev depletion tests</h1><p>Generated '+n.now()+'; producer finished.</p><img style="width:100%;max-width:1400px" src="duration_comparison.png"><p><a href="README.md">Design and limits</a> · <a href="endpoints.csv">Endpoints</a> · <a href="independent_audit.json">Audit</a></p></body></html>')


def self_check():
    L=39.;x=(np.arange(201)+.5)*L/201
    s=np.vstack([np.linspace(.001,100,201),np.linspace(.02,.4,201)])
    signal=np.linspace(0,6.4,201)
    native=Goryachev()
    model=AvailabilityControl(length=L,reference_v=.3,initial_v=.09)
    np.testing.assert_array_equal(model.reaction(0,x,s,signal),native.reaction(0,x,s,signal))
    for kind in KINDS:
        from dataclasses import replace
        a=replace(model,kind=kind)
        flux=a.reaction(0,x,s,signal)
        np.testing.assert_array_equal(flux.sum(axis=0),np.zeros(201))
        if kind=='feedback_clamp':np.testing.assert_array_equal(flux[:,:80],native.reaction(0,x,s,signal)[:,:80])
        depleted=np.vstack([s[0],np.zeros(201)])
        np.testing.assert_allclose(a.reaction(0,x,depleted,signal)[0],-a.c*depleted[0])
    print('Native identity, term separation, incumbent preservation, zero-substrate boundary and conservation checks passed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['init','worker','finish','self-check'])
    parser.add_argument('--out',type=Path,default=OUT);parser.add_argument('--slot',type=int)
    args=parser.parse_args();out=args.out.resolve()
    try:
        if args.action=='self-check':self_check()
        elif args.action=='init':initialize(out)
        elif args.action=='worker':worker(out,args.slot)
        else:finish(out)
    except Exception:
        if out.exists():
            path=out/'results'/('worker_'+str(args.slot)+'.error.json' if args.action=='worker' else 'error.json')
            r.write(path,dict(generated_at=n.now(),error=traceback.format_exc()))
        raise
