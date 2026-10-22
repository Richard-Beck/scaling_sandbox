"""Development beta=0 reversal campaign from one settled start per length.

Zero only delta_beta_tension, retaining the native homogeneous reference and
protrusion threshold, as in the preceding global tension-to-Rho control.
"""
from pathlib import Path
from dataclasses import replace
import argparse, copy, json, os, shlex, shutil, subprocess, sys
import numpy as np
CODE=Path(__file__).resolve().parent
sys.path[:0]=[str(CODE/'core'),str(CODE/'mechanisms')]
import native as n
import corrected_reversal as r
import polarity_establishment as e
import settling
import length_dev_output_layout as layout
DEFAULT=n.DATA/'runs/debelly_beta0_reversal_20261010'
BASE=n.DATA/'runs/debelly_corrected_three_test_20261008'
OUT=n.ROOT/'figures/length_mechanisms_20261008/dev'
TESTS=['fixed_physical','fixed_relative','fixed_slope']
RECOVERY=r.RECOVERY;STEP=r.STEP;TAIL=r.TAIL
FEEDBACK_AUDITS=[]
_original_feedback=n.ds._feedback_terms
_original_run=n.run
_original_model=n.model_for

def feedback(p,mca,tension):
    alpha,beta=_original_feedback(p,mca,tension)
    return alpha,np.zeros_like(beta)
n.ds._feedback_terms=feedback

def model_for(key):
    assert key=='debelly'
    return n.ds.DeBelly(replace(n.db._parameters(),center_boundary_motion=False))
n.model_for=model_for

def checked_run(key,L,state,cue,end,**kwargs):
    tr=_original_run(key,L,state,cue,end,**kwargs)
    p=model_for(key).parameters
    if state is not None:
        for name in n.DB_NAMES:np.testing.assert_array_equal(tr.fields[name][0],state[name])
    assert all(np.isfinite(v).all() for v in tr.fields.values())
    np.testing.assert_array_equal(tr.fields['delta_beta_tension'],np.zeros_like(tr.fields['tension']))
    da,_=_original_feedback(p,tr.fields['mca_bound'],tr.fields['tension'])
    np.testing.assert_array_equal(tr.fields['delta_alpha_mca'],da)
    np.testing.assert_array_equal(tr.fields['rho_activation'],np.ones_like(da)*p.rho_basal_activation)
    mass=(tr.fields['mca_bound']+tr.fields['mca_unbound']).sum(axis=1)
    err=float(np.max(abs(mass/mass[0]-1)));assert err<1e-10
    FEEDBACK_AUDITS.append(dict(initial_full_state_exact=True,delta_beta_zero_exact=True,native_delta_alpha_exact=True,mca_mass_error=err,end_time=float(end)))
    return tr
n.run=checked_run

class Gradient:
    event_times=()
    def __init__(self,g):self.g=g
    def __call__(self,t,x,L):return self.g*np.asarray(x)
def response_cue(c):
    if c['condition']=='fixed_slope':return Gradient(c['gradient_slope'])
    return e.Patch(c['response_amplitude'],c['width'],side='right',width_mode=c['width_mode'],fraction=.2)

def sources():
    paths=[Path(__file__),CODE/'debelly_beta0_reversal_audit.py',CODE/'audit_reversal.py',CODE/'mechanisms/length_dev_output_layout.py',BASE/'manifest.json',BASE/'results/independent_audit.json',n.DATA/'establishment_basal/manifest.json',r.GRID]
    return dict(n.sources(),**{str(p.relative_to(n.ROOT)):n.sha(p) for p in paths})
def freeze(out):
    hashes=sources()
    for name,digest in hashes.items():
        target=out/'frozen_sources'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(n.ROOT/name,target)
        assert n.sha(target)==digest
    return hashes

def verify_preparation(out):
    p=r.read(out/'preparation_provenance.json')
    assert sources()==p['sources'] and n.sha(out/'preparation_plan.json')==p['plan_sha256']
    return r.read(out/'preparation_plan.json')

def prepare_init(out,policy):
    if out.exists():raise ValueError('Fresh destination required')
    for name in ['work/preparations','work/chunks','work/jobs','work/logs','results','frozen_sources']:(out/name).mkdir(parents=True,exist_ok=True)
    assert r.read(BASE/'results/independent_audit.json')['passed']
    cases=[copy.deepcopy(c) for c in r.read(BASE/'manifest.json')['cases'] if c['model']=='debelly' and c['condition']=='fixed_relative']
    for slot,c in enumerate(cases):c['preparation_slot']=slot
    rule=r.read(n.DATA/'establishment_basal/manifest.json')['steady_state']
    plan=dict(created_at=n.now(),cases=cases,steady_state=rule,eligibility_policy=policy,
        intervention='delta_beta_tension exactly zero throughout settling, cue and recovery; original parameters, homogeneous reference, and Rac protrusion threshold retained.',
        branch='Continue the twelve existing native polarized branches; ten original preparation exclusions remain explicit, without a new de novo establishment search.')
    r.write(out/'preparation_plan.json',plan)
    r.write(out/'preparation_provenance.json',dict(sources=freeze(out),plan_sha256=n.sha(out/'preparation_plan.json')))
    r.write(out/'results/status.json',dict(generated_at=n.now(),producer='not started',complete=False,stage='beta0_preparation',lengths=len(cases)))
    write_index(out)

def prepare_worker(out,slot):
    plan=verify_preparation(out);c=plan['cases'][slot]
    record=dict(length=c['length'],original_checkpoint=c['checkpoint'],original_checkpoint_sha256=c['checkpoint_sha256'],eligible=False,verified=False,original_eligible=c['eligible'])
    if c['eligible']:
        assert n.sha(n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
        state,_=n.read_checkpoint(n.ROOT/c['checkpoint'],'debelly')
        tr,check=e.release(c,state,plan);record.update(check)
        if not check['verified']:raise RuntimeError('Unresolved beta=0 steady state')
        path=out/'work/preparations'/f'L{c["length"]:g}.npz';n.save_checkpoint(path,'debelly',tr.state,tr)
        ori=float(tr.orientation[-1]);contrast=float(tr.contrast[-1])
        original_gate=ori<-.05 and contrast>=.05
        asymmetric=ori<-.005 and contrast>=.05
        record.update(checkpoint=str(path.relative_to(n.ROOT)),checkpoint_sha256=n.sha(path),initial_orientation=ori,initial_contrast=contrast,
            original_cutoff_eligible=bool(original_gate),asymmetric_eligible=bool(asymmetric),eligible=bool(original_gate if plan['eligibility_policy']=='original' else asymmetric),feedback_audits=list(FEEDBACK_AUDITS))
    r.write(out/'work/preparations'/f'{slot:03d}.json',record)

def initialize(out):
    plan=verify_preparation(out);cases=[]
    slope=next(c['gradient_slope'] for c in r.read(BASE/'manifest.json')['cases'] if c['model']=='debelly' and c['condition']=='fixed_slope')
    for base in plan['cases']:
        prep_path=out/'work/preparations'/f'{base["preparation_slot"]:03d}.json';prep=r.read(prep_path)
        for geometry in TESTS:
            c=copy.deepcopy(base);c.update(case_index=len(cases),condition=geometry,eligible=prep['eligible'],variant_preparation=str(prep_path.relative_to(n.ROOT)),original_checkpoint=prep['original_checkpoint'],original_checkpoint_sha256=prep['original_checkpoint_sha256'],original_cutoff_eligible=prep.get('original_cutoff_eligible',False),orientation_cutoff=.05 if plan['eligibility_policy']=='original' else .005)
            if base['eligible']:
                for k in ['checkpoint','checkpoint_sha256','initial_orientation','initial_contrast']:c[k]=prep[k]
            c.update(preparation_status='polarized' if c['eligible'] else 'original_or_beta0_ineligible',response_amplitude=3.6,response_factor=1.)
            if geometry=='fixed_physical':c.update(width=3.,width_mode='physical')
            elif geometry=='fixed_relative':c.update(width=.2*c['length'],width_mode='relative')
            else:c['gradient_slope']=slope
            cases.append(c)
    m=dict(created_at=n.now(),models=['debelly'],cases=cases,conditions=TESTS,experiment='debelly_beta0_shared_steady_reversal',intervention=plan['intervention'],branch=plan['branch'],eligibility_policy=plan['eligibility_policy'],center_boundary_motion=False,preparation_plan_sha256=n.sha(out/'preparation_plan.json'),recovery_horizon=RECOVERY,sample_interval=STEP,gate_window=TAIL,refinement=dict(relative_bracket_tolerance=.01,maximum_rounds=r.MAX_ROUNDS))
    r.write(out/'manifest.json',m);r.write(out/'provenance.json',dict(sources=freeze(out),manifest_sha256=n.sha(out/'manifest.json')))
    jobs=[dict(slot=i,case_index=c['case_index'],durations=[0.]) for i,c in enumerate(c for c in cases if c['eligible'])]
    r.save_jobs(out,0,jobs);r.write(out/'launch.json',dict(generated_at=n.now(),stage='initialized',round=0,arrays=[]));publish(out,'initialized')

def trial(case,duration):
    checkpoint=n.ROOT/case['checkpoint']
    if n.sha(checkpoint)!=case['checkpoint_sha256']:raise RuntimeError('Checkpoint changed')
    state,prep=n.read_checkpoint(checkpoint,case['model'])
    pulse=None
    if duration>0:
        # Integrate through the left limit, then restart under basal withdrawal.
        cue=response_cue(case)
        pulse=n.run(case['model'],case['length'],state,cue,duration,interval=min(STEP,duration/20),n=case['cells'])
        state=pulse.state
    recovery=n.run(case['model'],case['length'],state,e.withdrawal(case),RECOVERY,interval=STEP,n=case['cells'])
    tail=recovery.time>=RECOVERY-TAIL;o=recovery.orientation[tail];contrast=recovery.contrast[tail]
    if len(o)!=201 or not all(np.isfinite(v).all() for v in [o,contrast]):raise RuntimeError('Invalid tail samples')
    orientation=bool(np.min(o)>case['orientation_cutoff'] and np.min(contrast)>=.05 and np.ptp(o)<.005)
    error=float(np.max(abs(o/abs(case['initial_orientation'])-1)))
    full=bool(orientation and error<=.01)
    fields=settling.monitored(recovery,case['model'],case['length']);initial_fields=settling.monitored(prep,case['model'],case['length'])
    arrays=dict(x=recovery.x,pulse_time=pulse.time if pulse else np.array([0.]),
                pulse_orientation=pulse.orientation if pulse else np.array([case['initial_orientation']]),
                pulse_contrast=pulse.contrast if pulse else np.array([case['initial_contrast']]),
                recovery_time=recovery.time,recovery_orientation=recovery.orientation,recovery_contrast=recovery.contrast)
    maximum=0.;no_cue_excursion=0.;metrics={}
    for name,y in fields.items():
        if not np.isfinite(y).all():raise RuntimeError('Nonfinite full native field')
        scale=case['fixed_scales'][name];a=y[tail]
        minimum=a.min(axis=0);maximum_field=a.max(axis=0)
        change=float(np.max(maximum_field-minimum))/scale
        drift=float(np.max(abs(a[-1]-a[-2])))/STEP*TAIL/scale
        maximum=max(maximum,change,drift);metrics[name]=dict(range_scaled=change,projected_drift_500=drift)
        for suffix,v in [('min',minimum),('max',maximum_field),('penultimate',a[-2]),('last',a[-1])]:arrays['fieldtail_'+name+'_'+suffix]=v
        if duration==0:
            no_cue_excursion=max(no_cue_excursion,float(np.max(abs(y-initial_fields[name][-1])))/scale)
            arrays['no_cue_'+name+'_min']=y.min(axis=0)
            arrays['no_cue_'+name+'_max']=y.max(axis=0)
            arrays['no_cue_'+name+'_initial']=initial_fields[name][-1]
    if isinstance(recovery.state,dict):arrays.update({'finalstate_'+k:v for k,v in recovery.state.items()})
    else:arrays['finalstate']=recovery.state
    row=dict(duration=float(duration),orientation_success=orientation,full_polarity_success=full,
             final_orientation=float(o[-1]),minimum_contrast_last_500=float(np.min(contrast)),orientation_range_last_500=float(np.ptp(o)),
             maximum_relative_polarity_error_last_500=error,field_tail_pass=bool(maximum<=.001),
             field_tail_max_scaled=maximum,field_tail_metrics=metrics,
             no_cue_maximum_scaled_excursion=no_cue_excursion if duration==0 else None,
             no_cue_left_tail_valid=bool(np.max(o)<-case['orientation_cutoff'] and np.ptp(o)<.005 and np.min(contrast)>=.05) if duration==0 else None,
             no_cue_preparation_valid=bool(no_cue_excursion<=.01 and not orientation and np.max(o)<-case['orientation_cutoff'] and np.ptp(o)<.005 and np.min(contrast)>=.05) if duration==0 else None)
    row['original_orientation_success']=bool(np.min(o)>.05 and np.min(contrast)>=.05 and np.ptp(o)<.005)
    row['feedback_audits']=list(FEEDBACK_AUDITS)
    FEEDBACK_AUDITS.clear()
    return row,arrays


_original_publish=r.publish

def write_index(out):
    s=r.read(out/'results/status.json');stamp=s['generated_at']
    p=layout.dev_output_path(OUT,'beta0_reversal_status.json');r.write(p,s)
    viewer=layout.dev_output_path(OUT,'beta0_reversal_status.html')
    link=os.path.relpath(out/'results',viewer.parent)
    viewer.write_text(f'<!doctype html><meta charset="utf-8"><title>Beta=0 reversal status</title><h1>Tension→Rho disabled: three reversal cues</h1><p>Data generated {stamp}; producer {s["producer"]}; stage {s["stage"]}.</p><p>One newly settled beta=0 checkpoint per length is shared exactly across all three cues. Pulse duration excludes 10,000-s basal recovery. Native reference and direct tension inhibition are retained.</p><p><a href="{link}/status.json">Status</a> · <a href="{link}/endpoints.csv">Current endpoints</a> · <a href="{link}/independent_audit.json">Audit</a></p>'+ ('<a href="../vectors/beta0_reversal.pdf"><img src="../../beta0_reversal.png" style="max-width:100%"></a>' if s['complete'] else '<p>Final figure appears after execution, refinement and independent audit.</p>'))

def publish(out,stage,complete=False):
    s=_original_publish(out,stage,complete);write_index(out);return s

def advance(out,round):
    m,groups,pending,errors=r.collect(out,partial=True)
    if pending or errors:raise RuntimeError(f'Incomplete execution: {pending}; {errors}')
    new=[]
    for c in m['cases']:
        if not c['eligible']:continue
        if not r.preparation_valid(c,groups[c['case_index']]):raise RuntimeError('Beta0 no-cue stability failed')
        ds=[d for d in c['durations'] if d not in {t['duration'] for t in groups[c['case_index']]}] if round==0 else r.refinement(groups[c['case_index']])
        for start in range(0,len(ds),5):new.append(dict(slot=len(new),case_index=c['case_index'],durations=ds[start:start+5]))
    if new:
        if round>=r.MAX_ROUNDS:raise RuntimeError('Refinement cap exhausted')
        r.save_jobs(out,round+1,new);return r.launch_round(out,round+1)
    rows,cases=r.rows_and_records(m,groups,True)
    assert not any(q['duration_boundary_unresolved'] or q['unresolved_transition_count'] for c in cases for q in c['result'].values())
    audit=out/'results/independent_audit.json'
    subprocess.run([sys.executable,str(CODE/'debelly_beta0_reversal_audit.py'),'--campaign',str(out),'--before-publication','--output',str(audit)],check=True)
    result=r.read(audit);assert result['passed'] and result['complete']
    lookup={(t['case_index'],t['definition']):t for t in result['selected_readouts']}
    for row in rows:
        for k in ['status','endpoint','selected_tail_field_pass_0_1_percent']:assert row[k]==lookup[row['case_index'],row['definition']][k]
    r.write(out/'results/validation.json',dict(generated_at=n.now(),passed=True,complete_execution=True,refinement_rounds=round,independent_audit_sha256=n.sha(audit)))
    render(out,rows)
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='complete',round=round,arrays=[]));publish(out,'complete',True)

def render(out,rows=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    assert r.read(out/'results/validation.json')['passed']
    m,groups,_,_=r.collect(out)
    if rows is None:rows,_=r.rows_and_records(m,groups,True)
    import csv
    native=list(csv.DictReader((BASE/'results/endpoints.csv').open()))
    fig,axs=plt.subplots(2,3,figsize=(12.7,8.4),layout='constrained',sharex=True)
    colors={'orientation':'#246a9a','full_polarity':'#b14f36'}
    for col,(geometry,title) in enumerate(zip(TESTS,['Fixed patch · 3 µm','Relative patch · 20%','Constant slope · 0.24'])):
        ax=axs[0,col];ax.set_title(title);ax.set_ylabel('Shortest passing cue (s)');ax.set_xlabel('Cell length (µm)');ax.grid(alpha=.15)
        for gate,color in colors.items():
            rs=sorted([q for q in rows if q['condition']==geometry and q['definition']==gate],key=lambda q:q['length'])
            ax.plot([q['length'] for q in rs],[q['endpoint'] if q['status']=='success' else np.nan for q in rs],color=color,lw=1.6,ls='-' if gate=='orientation' else '--')
            for q in rs:
                if q['status']=='success':ax.plot(q['length'],q['endpoint'],'o' if gate=='orientation' else 's',color=color,mfc=color if q['selected_tail_field_pass_0_1_percent'] else 'white',ms=5)
                elif q['status']=='censored':ax.plot(q['length'],q['maximum_tested_duration'],'^',color=color,mfc='white')
                elif gate=='orientation':ax.plot(q['length'],-.06,'x',color='#777',transform=ax.get_xaxis_transform(),clip_on=False)
            ns=sorted([q for q in native if q['model']=='debelly' and q['condition']==geometry and q['definition']==gate and q['status']=='success'],key=lambda q:float(q['length']))
            ax.plot([float(q['length']) for q in ns],[float(q['endpoint']) for q in ns],color=color,alpha=.35,lw=1,ls=':')
        ax.set_yscale('log')
    for ax,(field,label) in zip(axs[1],[('initial_orientation','Basal Rac polarity magnitude'),('initial_contrast','Basal Rac contrast'),('rac','Verified beta=0 Rac profiles')]):
        ax.set_title(label);ax.grid(alpha=.15);ax.set_xlabel('Cell length (µm)' if field!='rac' else 'Position x/L')
        bs=[c for c in m['cases'] if c['condition']=='fixed_relative' and c['checkpoint']]
        if field!='rac':
            vals=[abs(c[field]) for c in bs];ax.plot([c['length'] for c in bs],vals,'o-',color='#333',ms=4)
            if field=='initial_orientation':ax.axhline(.05,ls=':',color='#777',label='Original cutoff 0.05');ax.legend(frameon=False,fontsize=8)
        else:
            selected=[bs[0],bs[len(bs)//2],bs[-1]]
            for c in selected:
                _,tr=n.read_checkpoint(n.ROOT/c['checkpoint'],'debelly');ax.plot(tr.x/c['length'],tr.fields['rac'][-1],label=f'{c["length"]:g} µm')
            ax.legend(frameon=False,fontsize=8)
    handles=[Line2D([],[],color=c,marker='o' if g=='orientation' else 's',label='β=0 orientation' if g=='orientation' else 'β=0 1% strength') for g,c in colors.items()]
    handles+=[Line2D([],[],color='#777',ls=':',label='Native feedback baseline'),Line2D([],[],color='#777',marker='o',mfc='white',ls='',label='Field stability flag fails'),Line2D([],[],color='#777',marker='^',mfc='white',ls='',label='Censored')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False,fontsize=9)
    fig.suptitle('Reversal with tension→Rho feedback disabled',fontsize=15)
    artifacts={}
    for ext in ['png','pdf','svg']:
        path=layout.dev_output_path(OUT,'beta0_reversal.'+ext);fig.savefig(path,dpi=200);artifacts[layout.dev_artifact_key(OUT,path)]=n.sha(path)
    plt.close(fig)
    caption='Tension-dependent Rho activation is zero throughout new steady-state preparation, finite countercue and recovery. The original homogeneous reference and protrusion threshold, direct tension inhibition and MCA→Rac feedback remain. Each length shares its exact new checkpoint across all three cues. Top: shortest tested successful pulse, with all sampled fail→pass boundaries refined to ≤1%; recovery of 10,000 s is excluded. Blue: orientation; red: 1% replacement of the own-start polarity magnitude. Pale dotted lines show the native manuscript baseline. Open symbols fail the separate field-tail flag. Bottom: new basal asymmetry and selected full Rac profiles. Eligibility policy: '+m['eligibility_policy']+'. Asymmetric policy requires left moment <−0.005 and contrast ≥0.05; reversal requires right moment >0.005, contrast ≥0.05 and final-500-s moment range <0.005, followed by the separate 1% strength requirement. Original 0.05 eligibility and reversal flags are exported. Ten native preparation exclusions are retained; no de novo establishment search was performed.'
    cp=layout.dev_output_path(OUT,'beta0_reversal_caption.txt');cp.write_text(caption+'\n')
    for name in ['endpoints.csv','independent_audit.json','validation.json']:
        dest=layout.dev_output_path(OUT,'beta0_reversal_'+name);shutil.copy2(out/'results'/name,dest);artifacts[layout.dev_artifact_key(OUT,dest)]=n.sha(dest)
    r.write(layout.dev_output_path(OUT,'beta0_reversal.manifest.json'),dict(generated_at=n.now(),producer='finished',artifact_root='../..',caption=caption,campaign=str(out.relative_to(n.ROOT)),sha256=artifacts,checkpoint_sha256={str(c['length']):c['checkpoint_sha256'] for c in m['cases'] if c['condition']=='fixed_relative' and c['eligible']},sources=sources(),reproduction=f'{sys.executable} {__file__} render --out {out}'))

r.RUNNER=Path(__file__).resolve();r.DEFAULT=DEFAULT;r.response_cue=response_cue;r.source_hashes=sources;r.trial=trial;r.advance=advance;r.publish=publish

def main():
    if len(sys.argv)>1 and sys.argv[1] in ['prepare-init','prepare-worker','prepare-finish','prepare-launch','render']:
        ap=argparse.ArgumentParser();ap.add_argument('action');ap.add_argument('--out',type=Path,default=DEFAULT);ap.add_argument('--slot',type=int);ap.add_argument('--offset',type=int,default=0);ap.add_argument('--round',type=int,default=0);ap.add_argument('--policy',choices=['asymmetric','original'],default='asymmetric');a=ap.parse_args();out=a.out.resolve()
        if a.action=='prepare-init':prepare_init(out,a.policy)
        elif a.action=='prepare-worker':prepare_worker(out,a.slot+a.offset)
        elif a.action=='prepare-finish':initialize(out);r.launch_round(out,0);r.submit(out,'monitor',0)
        elif a.action=='render':render(out);write_index(out)
        else:
            plan=verify_preparation(out);array=r.submit(out,'prepare-worker',0,array=f'0-{len(plan["cases"])-1}');final=r.submit(out,'prepare-finish',0,dependency=array)
            r.write(out/'preparation_launch.json',dict(generated_at=n.now(),array=array,finalizer=final));r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='running',stage='beta0_preparation',lengths=len(plan['cases'])));write_index(out);print(array,final)
    else:r.main()
if __name__=='__main__':main()
