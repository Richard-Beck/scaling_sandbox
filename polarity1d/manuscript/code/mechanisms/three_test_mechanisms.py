"""Native spatial replays and bounded causal controls for all 15 comparison cells."""
from pathlib import Path
from dataclasses import replace
import argparse, copy, gzip, json, shlex, shutil, subprocess, sys, traceback
import numpy as np

CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE))
import three_test_reversal as t
n,r,e=t.n,t.r,t.e
OUT=n.DATA/'runs/three_test_mechanisms_20261006'


def submit(out,action,array=None,dependency=None):
    command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__)),action,'--out',str(out)])
    if array:command+=' --slot "$SLURM_ARRAY_TASK_ID"'
    opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=4G','--time=02:00:00',
          '--job-name=three-test-mechanism-'+action,'--chdir='+str(n.ROOT),'--output='+str(out/'logs'/f'{action}_%A_%a.out'),'--error='+str(out/'logs'/f'{action}_%A_%a.err')]
    if array:opts.append('--array='+array)
    if dependency:opts.append('--dependency=afterany:'+dependency)
    return subprocess.check_output(opts+['--wrap',command],text=True).strip().split(';')[0]


def initialize(out,campaign):
    if out.exists():raise ValueError('Use a fresh mechanism destination')
    status=r.read(campaign/'results/status.json');assert status['complete']
    with gzip.open(campaign/'results/records.json.gz','rt') as f:cases=json.load(f)['cases']
    jobs=[];coverage=[]
    for model in r.MODELS:
        for geometry in t.TESTS:
            cs=sorted([c for c in cases if c['model']==model and c['condition']==geometry and c['eligible']],key=lambda c:c['length'])
            lengths={cs[0]['length'],cs[len(cs)//2]['length'],cs[-1]['length']}
            successes=[c for c in cs if c['result']['orientation']['status']=='success']
            # Add an interior extremum when it is not among the three anchors.
            if len(successes)>3:
                minimum=min(successes,key=lambda c:c['result']['orientation']['endpoint'])
                maximum=max(successes,key=lambda c:c['result']['orientation']['endpoint'])
                if minimum['length'] not in [cs[0]['length'],cs[-1]['length']]:lengths.add(minimum['length'])
                if maximum['length'] not in [cs[0]['length'],cs[-1]['length']]:lengths.add(maximum['length'])
            selected=[c for c in cs if c['length'] in lengths]
            coverage.append(dict(model=model,geometry=geometry,lengths=sorted(lengths),eligible_lengths=len(cs),scalar_successes=len(successes)))
            for old in selected:
                c={k:v for k,v in old.items() if k not in ['tests','result']}
                replays={}
                for gate in r.DEFINITIONS:
                    passing=[row for row in old['tests'] if row['duration']>0 and row[gate+'_success']]
                    if passing:
                        best=min(passing,key=lambda row:row['duration']);failed=max([row for row in old['tests'] if row['duration']<best['duration'] and not row[gate+'_success']],key=lambda row:row['duration'])
                        for label,row in [('pass',best),('fail',failed)]:replays.setdefault(row['duration'],dict(row=row,labels=[]))['labels'].append(gate+'_'+label)
                    else:
                        row=max(old['tests'],key=lambda row:row['duration']);replays.setdefault(row['duration'],dict(row=row,labels=[]))['labels'].append(gate+'_censored')
                for duration,item in replays.items():jobs.append(dict(slot=len(jobs),case=c,duration=duration,kind='native',expected=item['row'],labels=item['labels']))
                duration=old['result']['orientation']['endpoint'] or max(old['durations'])
                kinds=['mobility_4x']
                if geometry=='fixed_slope':kinds.append('matched_gradient_peak')
                if model=='debelly':kinds.append('mca_frozen_during_pulse')
                if model=='spring':kinds.append('feedback_off_during_pulse')
                if old['length']==cs[-1]['length']:kinds.append('mesh_401')
                for kind in kinds:jobs.append(dict(slot=len(jobs),case=c,duration=duration,kind=kind,labels=['orientation_selected_diagnostic']))
    for folder in ['logs','trajectories','results','frozen_sources']:(out/folder).mkdir(parents=True,exist_ok=True)
    r.write(out/'jobs.json',jobs)
    paths=[Path(__file__),CODE/'three_test_reversal.py',CODE/'audit_reversal.py',campaign/'manifest.json',campaign/'results/endpoints.csv']
    hashes=n.sources();hashes.update({str(p.relative_to(n.ROOT)):n.sha(p) for p in paths})
    for name in hashes:
        dest=out/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(n.ROOT/name,dest)
    r.write(out/'manifest.json',dict(generated_at=n.now(),campaign=str(campaign.relative_to(n.ROOT)),coverage=coverage,jobs=len(jobs),
            design='Short, middle, long eligible lengths plus interior duration extrema in every model-test cell. Adjacent native passing/failing doses for both gates, deduplicated. Controls alter declared physics only and never populate comparison curves.',
            control_definitions=dict(mobility_4x='Fourfold inactive diffusivity (both Rac/Rho in Holmes); deBelly fourfold unbound-MCA diffusion; Spring fourfold inactive diffusion. Pulse and recovery; exact original starting checkpoint.',
            matched_gradient_peak='Diagnostic slope times reference length/L, so ramp peak matches its reference-length peak; this length-varying input is a control only.',
            mca_frozen_during_pulse='MCA exchange/advection/diffusion suppressed only during pulse; native recovery restored.',
            feedback_off_during_pulse='Spring extension-dependent global deactivation coefficient d1=0 during pulse only; native recovery restored.',
            mesh_401='Interpolate all fields to 401 samples, preserve spatial means (trapezoidal for Spring); modified discretization diagnostic, not exact preparation or new threshold.')))
    r.write(out/'provenance.json',dict(generated_at=n.now(),sources=hashes,jobs_sha256=n.sha(out/'jobs.json'),manifest_sha256=n.sha(out/'manifest.json')))
    r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='running',jobs=len(jobs),coverage_cells=15))
    array=submit(out,'worker',f'0-{len(jobs)-1}');final=submit(out,'finish',dependency=array)
    r.write(out/'launch.json',dict(generated_at=n.now(),array=array,controller=final))
    print(json.dumps(dict(jobs=len(jobs),array=array,controller=final,coverage=coverage),indent=2))


def verify(out):
    p=r.read(out/'provenance.json')
    if n.sha(out/'jobs.json')!=p['jobs_sha256'] or n.sha(out/'manifest.json')!=p['manifest_sha256']:raise RuntimeError('Mechanism inputs changed')
    for name,digest in p['sources'].items():
        if n.sha(n.ROOT/name)!=digest:raise RuntimeError('Mechanism source changed: '+name)


def gates(tr,case):
    tail=tr.time>=9500.;o=tr.orientation[tail];q=tr.contrast[tail]
    orientation=bool(o.min()>.05 and q.min()>=.05 and np.ptp(o)<.005)
    error=float(np.max(abs(o/abs(case['initial_orientation'])-1.)))
    fields=r.settling.monitored(tr,case['model'],case['length']);maximum=0.
    for name,y in fields.items():
        a=y[tail];scale=case['fixed_scales'][name]
        maximum=max(maximum,float(np.ptp(a,axis=0).max())/scale,float(abs(a[-1]-a[-2]).max())/2.5*500/scale)
    return dict(orientation_success=orientation,full_polarity_success=bool(orientation and error<=.01),
                final_orientation=float(o[-1]),minimum_contrast_last_500=float(q.min()),orientation_range_last_500=float(np.ptp(o)),
                maximum_relative_polarity_error_last_500=error,field_tail_pass=bool(maximum<=.001),field_tail_max_scaled=maximum)


def interpolate_state(state,model,L,cells):
    old=len(state['rac']) if isinstance(state,dict) else (len(state)-2)//2 if model=='spring' else state.shape[1]
    a=np.linspace(0,1,old) if model=='spring' else (np.arange(old)+.5)/old
    b=np.linspace(0,1,cells) if model=='spring' else (np.arange(cells)+.5)/cells
    def remap(y):
        z=np.interp(b,a,y)
        oldmean=np.trapezoid(y,a) if model=='spring' else np.mean(y)
        newmean=np.trapezoid(z,b) if model=='spring' else np.mean(z)
        if np.min(y)>=0 and newmean>0:z*=oldmean/newmean
        else:z+=oldmean-newmean
        return z
    if isinstance(state,dict):return {key:remap(value) for key,value in state.items()}
    if model=='spring':
        new=np.empty(2*cells+2);new[:-2:2]=remap(state[:-2:2]);new[1:-2:2]=remap(state[1:-2:2]);new[-2:]=state[-2:];return new
    return np.vstack([remap(y) for y in state])


def worker(out,slot):
    verify(out);job=r.read(out/'jobs.json')[slot];c=copy.deepcopy(job['case']);kind=job['kind'];D=job['duration'];model=c['model'];L=c['length']
    checkpoint=n.ROOT/c['checkpoint']
    if n.sha(checkpoint)!=c['checkpoint_sha256']:raise RuntimeError('Checkpoint changed')
    state,prep=n.read_checkpoint(checkpoint,model)
    initial_state=copy.deepcopy(state);original_model=n.model_for
    original_transport=n.ds._transport;original_rhs=n.sp.rhs
    if kind=='mesh_401':state=interpolate_state(state,model,L,401);c['cells']=401
    if kind=='matched_gradient_peak':c['gradient_slope']*= (1. if model=='spring' else 15.)/L
    if kind=='mobility_4x':
        m=original_model(model)
        if model in ['wave_pinning','goryachev']:m=replace(m,dv=m.dv*4)
        elif model=='holmes_model3':m=replace(m,inactive_diffusion=m.inactive_diffusion*4)
        elif model=='debelly':m=replace(m,parameters=replace(m.parameters,mca_unbound_diffusion=m.parameters.mca_unbound_diffusion*4))
        else:
            def mobility_rhs(*args,**kwargs):
                y=args[1];domain=y[-1]-y[-2];v=y[1:-2:2];dx=1/(len(v)-1)
                result=original_rhs(*args,**kwargs)
                lap=np.r_[2*(v[1]-v[0]),np.diff(v,2),2*(v[-2]-v[-1])]/dx**2
                result[1:-2:2]+=30/domain**2*lap
                return result
            n.sp.rhs=mobility_rhs
        if model!='spring':n.model_for=lambda key:m if key==model else original_model(key)
    elif kind=='mca_frozen_during_pulse':n.ds._transport=lambda values,*args,**kwargs:values.copy()
    elif kind=='feedback_off_during_pulse':
        def no_feedback(*args,**kwargs):kwargs['d1']=0.;return original_rhs(*args,**kwargs)
        n.sp.rhs=no_feedback
    pulse=n.run(model,L,state,t.response_cue(c),D,interval=min(2.5,D/20),n=c['cells']) if D>0 else None
    if kind=='mca_frozen_during_pulse':n.ds._transport=original_transport
    if kind=='feedback_off_during_pulse':n.sp.rhs=original_rhs
    recovery=n.run(model,L,pulse.state if pulse else state,e.withdrawal(c),10000.,interval=2.5,n=c['cells'])
    metrics=gates(recovery,c)
    if kind=='native':
        for key in ['orientation_success','full_polarity_success','field_tail_pass']:
            if metrics[key]!=job['expected'][key]:raise RuntimeError('Native replay gate differs: '+key)
        if abs(metrics['final_orientation']-job['expected']['final_orientation'])>1e-9:raise RuntimeError('Native replay endpoint differs')
    arrays={}
    for phase,tr in [('pulse',pulse),('recovery',recovery)]:
        if tr is None:continue
        for key in ['time','x','active','orientation','contrast']:arrays[phase+'_'+key]=getattr(tr,key)
        arrays.update({phase+'_field_'+key:value for key,value in tr.fields.items()})
    if isinstance(initial_state,dict):arrays.update({'originalstate_'+k:v for k,v in initial_state.items()})
    else:arrays['originalstate']=initial_state
    path=out/'trajectories'/f'trajectory_{slot:04d}.npz';np.savez_compressed(path,**arrays)
    r.write(path.with_suffix('.json'),dict(generated_at=n.now(),job=job,metrics=metrics,archive_sha256=n.sha(path),checkpoint_sha256=c['checkpoint_sha256']))


def finish(out):
    verify(out);jobs=r.read(out/'jobs.json');missing=[j['slot'] for j in jobs if not (out/'trajectories'/f'trajectory_{j["slot"]:04d}.json').exists()]
    if missing:raise RuntimeError('Missing mechanism trajectories: '+str(missing))
    subprocess.run([sys.executable,str(CODE/'mechanisms/three_test_analysis.py'),'--out',str(out)],check=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','worker','finish']);p.add_argument('--out',type=Path,default=OUT);p.add_argument('--campaign',type=Path,default=t.DEFAULT);p.add_argument('--slot',type=int)
    args=p.parse_args();out=args.out.resolve()
    try:
        if args.action=='init':initialize(out,args.campaign.resolve())
        elif args.action=='worker':worker(out,args.slot)
        else:finish(out)
    except Exception:
        if out.exists():r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='stopped with error',error=traceback.format_exc()))
        raise
