"""Isolated tenfold-reduced active/inactive diffusion, three finite-cue tests.

Restart the retained nominal relative-20% basal state, settle the variant under
basal input, then share that exact checkpoint across all three response tests.
Original manuscript sources and evidence are never modified.
"""
from pathlib import Path
import argparse, copy, gzip, json, os, shlex, shutil, subprocess, sys
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE/'core'))
import native as n
import corrected_reversal as r
import polarity_establishment as e
from models.wave_pinning import WavePinning

DEFAULT = n.DATA/'runs/wave_pinning_slow_diffusion_20261008'
BASE = n.DATA/'runs/three_test_reversal_20261006'
TESTS = ['fixed_physical', 'fixed_relative', 'fixed_slope']
PARAMETERS = dict(du=.01, dv=1.)
_native_model = n.model_for

def model_for(key):
    return WavePinning(**PARAMETERS) if key == 'wave_pinning' else _native_model(key)
n.model_for = model_for

class Gradient:
    event_times = ()
    def __init__(self, slope): self.slope = float(slope)
    def __call__(self, time, x, length): return self.slope*np.asarray(x)

def response_cue(c):
    if c['condition']=='fixed_slope': return Gradient(c['gradient_slope'])
    return e.Patch(c['response_amplitude'],c['width'],side='right',width_mode=c['width_mode'],fraction=.2)

def sources():
    result=n.sources()
    for p in [Path(__file__),CODE/'wave_pinning_slow_diffusion_audit.py',CODE/'audit_reversal.py',BASE/'manifest.json',r.GRID,n.DATA/'establishment_basal/manifest.json']:
        result[str(p.relative_to(n.ROOT))]=n.sha(p)
    return result

def freeze_sources(out):
    hashes=sources()
    for name in hashes:
        dest=out/'frozen_sources'/name
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(n.ROOT/name,dest)
    return hashes

def verify_preparation(out):
    provenance=r.read(out/'preparation_provenance.json')
    if sources()!=provenance['sources'] or n.sha(out/'preparation_plan.json')!=provenance['plan_sha256']:
        raise RuntimeError('Preparation scientific inputs changed')
    return r.read(out/'preparation_plan.json')

def prepare_init(out,mode):
    n.protect_retained_output(out)
    if out.exists():raise ValueError('Use a fresh output destination')
    for folder in ['work/preparations','work/chunks','work/jobs','work/logs','results','frozen_sources']:
        (out/folder).mkdir(parents=True,exist_ok=True)
    cases=[copy.deepcopy(c) for c in r.read(BASE/'manifest.json')['cases'] if c['model']=='wave_pinning' and c['condition']=='fixed_relative']
    for slot,c in enumerate(cases):c['preparation_slot']=slot
    plan=dict(created_at=n.now(),parameters=PARAMETERS,preparation_mode=mode,cases=cases,
              steady_state=r.read(n.DATA/'establishment_basal/manifest.json')['steady_state'])
    r.write(out/'preparation_plan.json',plan)
    r.write(out/'preparation_provenance.json',dict(sources=freeze_sources(out),plan_sha256=n.sha(out/'preparation_plan.json')))
    r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='not started',stage='variant_preparation',lengths=len(cases)))
    write_index(out)

def prepare_worker(out,slot):
    plan=verify_preparation(out);c=plan['cases'][slot]
    record=dict(original_checkpoint=c['checkpoint'],original_checkpoint_sha256=c['checkpoint_sha256'],length=c['length'],eligible=False,verified=False)
    if c['eligible']:
        if n.sha(n.ROOT/c['checkpoint'])!=c['checkpoint_sha256']:raise RuntimeError('Original checkpoint changed')
        state,tr=n.read_checkpoint(n.ROOT/c['checkpoint'],c['model'])
        if plan['preparation_mode']=='settle':
            tr,check=e.release(c,state,plan)
            record.update(check)
            path=out/'work/preparations'/f'L{c["length"]:g}.npz'
            n.save_checkpoint(path,c['model'],tr.state,tr)
        else:
            record.update(verified=True,elapsed=0.,checks=[])
            path=n.ROOT/c['checkpoint']
        record.update(checkpoint=str(path.relative_to(n.ROOT)),checkpoint_sha256=n.sha(path),
                      initial_orientation=float(tr.orientation[-1]),initial_contrast=float(tr.contrast[-1]))
        record['eligible']=bool(record['verified'] and tr.orientation[-1]<-.05 and tr.contrast[-1]>=.05)
        if not record['verified']:raise RuntimeError('Variant basal settling unresolved; do not treat as biological exclusion')
    r.write(out/'work/preparations'/f'{slot:03d}.json',record)

def initialize(out):
    plan=verify_preparation(out);cases=[]
    slope=next(c['gradient_slope'] for c in r.read(BASE/'manifest.json')['cases'] if c['model']=='wave_pinning' and c['condition']=='fixed_slope')
    for base in plan['cases']:
        p=out/'work/preparations'/f'{base["preparation_slot"]:03d}.json'
        if not p.exists():raise RuntimeError('Missing variant preparation: '+str(p))
        prep=r.read(p)
        for geometry in TESTS:
            c=copy.deepcopy(base)
            c.update(case_index=len(cases),condition=geometry,variant_parameters=PARAMETERS,
                     original_checkpoint=prep['original_checkpoint'],original_checkpoint_sha256=prep['original_checkpoint_sha256'],
                     eligible=prep['eligible'],variant_preparation=str(p.relative_to(n.ROOT)))
            if base['eligible']:
                for k in ['checkpoint','checkpoint_sha256','initial_orientation','initial_contrast']:c[k]=prep[k]
            c.update(preparation_status='polarized' if c['eligible'] else 'variant_or_original_unpolarized',
                     response_amplitude=.1,response_factor=1.,preparation_condition='relative_20pct')
            if geometry=='fixed_physical':c.update(width=3.,width_mode='physical',width_fraction_at_reference=.2)
            elif geometry=='fixed_relative':c.update(width=.2*c['length'],width_mode='relative',width_fraction_at_reference=.2)
            else:c['gradient_slope']=slope
            cases.append(c)
    manifest=dict(created_at=n.now(),models=['wave_pinning'],cases=cases,conditions=TESTS,
                  experiment='wave_pinning_tenfold_reduced_diffusion',parameters=PARAMETERS,baseline_parameters=dict(du=.1,dv=10.),
                  preparation_mode=plan['preparation_mode'],preparation_plan_sha256=n.sha(out/'preparation_plan.json'),
                  recovery_horizon=r.RECOVERY,sample_interval=r.STEP,gate_window=r.TAIL,
                  refinement=dict(relative_bracket_tolerance=.01,maximum_rounds=r.MAX_ROUNDS),
                  preparation='Retained nominal relative-20% polarized branch continued under reduced diffusion to verified basal state, shared exactly across all three cues.' if plan['preparation_mode']=='settle' else 'Exact retained relative-20% baseline state; diffusion reduced at reversal onset.',
                  cue='Fixed width 3 and relative 20% cosine patches at peak 0.1; finite constant physical slope inherited from the manuscript.',
                  readout='Shortest tested finite cue duration giving sustained right orientation; separate 1% scalar-strength replacement gate; excludes basal recovery.',
                  no_cue_policy='Native original no-cue full-profile stability gate; variant-only onset can fail this diagnostic.',
                  retained_length_membership=True)
    r.write(out/'manifest.json',manifest)
    hashes=freeze_sources(out)
    r.write(out/'provenance.json',dict(generated_at=n.now(),sources=hashes,manifest_sha256=n.sha(out/'manifest.json')))
    jobs=[dict(slot=i,case_index=c['case_index'],durations=[0.]) for i,c in enumerate(c for c in cases if c['eligible'])]
    r.save_jobs(out,0,jobs)
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='initialized',round=0,arrays=[]))
    publish(out,'initialized')

_original_publish=r.publish

def write_index(out):
    s=r.read(out/'results/status.json');stamp=s['generated_at']
    images=''.join(f'<h2>{title}</h2><a href="{g}.pdf"><img src="{g}.png" style="max-width:900px;width:100%"></a>' for g,title in zip(TESTS,['Fixed support 3','Relative support 20%','Constant physical slope'])) if s['complete'] else '<p>Final plots will appear after complete execution and independent audit.</p>'
    (out/'results/index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>Wave-pinning: reduced diffusion</title><h1>Wave-pinning: tenfold reduced diffusion</h1><p>Active diffusion 0.01; inactive diffusion 1. Data generated {stamp}; producer {s["producer"]}; stage {s["stage"]}.</p><p>Readout: shortest passing finite cue duration, excluding 10,000-unit basal recovery. Filled markers pass the separate native-field stability flag; open markers fail it. Triangles mark tested censor limits; crosses mark ineligible preparations.</p><p><a href="status.json">Status</a> · <a href="endpoints.csv">Endpoint data</a> · <a href="../manifest.json">Protocol</a> · <a href="independent_audit.json">Audit</a></p>{images}')

def publish(out,stage,complete=False):
    result=_original_publish(out,stage,complete)
    write_index(out)
    return result

def advance(out,round):
    manifest,groups,pending,errors=r.collect(out,partial=True)
    if pending or errors:raise RuntimeError(f'Incomplete execution: {len(pending)} pending, {len(errors)} errors')
    planned=[]
    for c in manifest['cases']:
        if not r.preparation_valid(c,groups[c['case_index']]):continue
        ds=[d for d in c['durations'] if d not in {t['duration'] for t in groups[c['case_index']]}] if round==0 else r.refinement(groups[c['case_index']])
        for start in range(0,len(ds),10):planned.append(dict(slot=len(planned),case_index=c['case_index'],durations=ds[start:start+10]))
    if planned:
        if round>=r.MAX_ROUNDS:raise RuntimeError('Refinement cap exhausted')
        r.save_jobs(out,round+1,planned)
        return r.launch_round(out,round+1)
    rows,cases=r.rows_and_records(manifest,groups,True)
    if any(q['duration_boundary_unresolved'] or q['unresolved_transition_count'] for c in cases for q in c['result'].values()):raise RuntimeError('Unresolved sampled boundaries')
    audit=out/'results/independent_audit.json'
    subprocess.run([sys.executable,str(CODE/'wave_pinning_slow_diffusion_audit.py'),'--campaign',str(out),'--before-publication','--output',str(audit)],check=True)
    result=r.read(audit)
    if not result['passed'] or not result['complete']:raise RuntimeError('Independent audit failed')
    lookup={(t['case_index'],t['definition']):t for t in result['selected_readouts']}
    for row in rows:
        expected=lookup[row['case_index'],row['definition']]
        for key in ['status','endpoint','selected_tail_field_pass_0_1_percent']:
            if row[key]!=expected[key]:raise RuntimeError('Independent selection differs')
    r.write(out/'results/validation.json',dict(generated_at=n.now(),passed=True,independent_audit_passed=True,
                 independent_audit_sha256=n.sha(audit),complete_execution=True,refinement_rounds=round))
    render(out,rows)
    r.write(out/'launch.json',dict(generated_at=n.now(),stage='complete',round=round,arrays=[]))
    publish(out,'complete',True)

def render(out,rows=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    if rows is None:
        if not r.read(out/'results/validation.json')['passed']:raise RuntimeError('Audit required')
        m,groups,_,_=r.collect(out);rows,_=r.rows_and_records(m,groups,True)
    panels=[]
    for geometry,title in zip(TESTS,['Fixed support 3','Relative support 20%','Constant physical slope']):
        fig,ax=plt.subplots(figsize=(7.8,5.1))
        for gate,color,marker,label in [('orientation','#246a9a','o','Stable orientation reversal'),('full_polarity','#b14f36','s','1% scalar strength replacement')]:
            rs=sorted([q for q in rows if q['condition']==geometry and q['definition']==gate],key=lambda q:q['length'])
            ax.plot([q['length'] for q in rs],[q['endpoint'] if q['status']=='success' else np.nan for q in rs],color=color,lw=1.4,ls='-' if gate=='orientation' else '--')
            for q in rs:
                if q['status']=='success':ax.plot(q['length'],q['endpoint'],marker,color=color,mfc=color if q['selected_field_tail_pass'] else 'white',ms=5)
                elif q['status']=='censored':ax.plot(q['length'],q['maximum_tested_duration'],'^',color=color,mfc='white',ms=6)
                elif gate=='orientation':ax.plot(q['length'],-.065,'x' if q['status']=='ineligible' else '+',transform=ax.get_xaxis_transform(),color='#777',clip_on=False)
        ax.set(xlabel='Cell length (native units)',ylabel='Shortest passing finite cue duration (native time)',title=f'Wave-pinning · {title}\nActive diffusion 0.01; inactive diffusion 1')
        ax.set_yscale('log');ax.grid(axis='y',alpha=.2)
        handles=[Line2D([],[],color=c,marker=m,label=label) for c,m,label in [('#246a9a','o','Orientation reversal'),('#b14f36','s','1% strength replacement')]]
        handles += [Line2D([],[],color='#777',marker='^',mfc='white',ls='',label='Censored at tested limit'),Line2D([],[],color='#777',marker='o',mfc='white',ls='',label='Field stability flag fails')]
        ax.legend(handles=handles,fontsize=8,frameon=False)
        fig.text(.12,.03,'Independent finite pulses; 10,000-unit recovery excluded.\nSame verified variant basal checkpoint for all three cues at each length.',fontsize=8)
        fig.subplots_adjust(bottom=.21,top=.83,left=.12,right=.97)
        for ext in ['png','pdf','svg']:fig.savefig(out/'results'/f'{geometry}.{ext}',dpi=180)
        panels.append(fig)
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out/'results/length_vs_reversal_time.pdf') as pdf:
        for fig in panels:pdf.savefig(fig)
    for fig in panels:plt.close(fig)
    r.write(out/'results/figure_provenance.json',dict(generated_at=n.now(),producer='finished',audit_sha256=n.sha(out/'results/independent_audit.json'),output_sha256={p.name:n.sha(p) for p in (out/'results').iterdir() if p.suffix in ['.png','.pdf','.svg']}))

r.RUNNER=Path(__file__).resolve();r.DEFAULT=DEFAULT
r.response_cue=response_cue;r.source_hashes=sources;r.advance=advance;r.publish=publish

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1] in ['prepare-init','prepare-worker','prepare-finish','prepare-launch','render']:
        parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action')
        parser.add_argument('--out',type=Path,default=DEFAULT);parser.add_argument('--slot',type=int)
        parser.add_argument('--round',type=int,default=0);parser.add_argument('--offset',type=int,default=0)
        parser.add_argument('--mode',choices=['settle','onset'],default='settle')
        args=parser.parse_args();out=args.out.resolve()
        if args.action=='prepare-init':prepare_init(out,args.mode)
        elif args.action=='prepare-worker':
            try:prepare_worker(out,args.slot)
            except Exception:
                import traceback
                r.write(out/'work/preparations'/f'{args.slot:03d}.error.json',dict(generated_at=n.now(),error=traceback.format_exc()))
                raise
        elif args.action=='prepare-finish':
            try:
                initialize(out);r.launch_round(out,0);r.submit(out,'monitor',0)
            except Exception:
                import traceback
                r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='stopped with error',stage='variant_preparation_failed',error=traceback.format_exc()))
                write_index(out)
                raise
        elif args.action=='render':render(out);write_index(out)
        else:
            plan=verify_preparation(out)
            array=r.submit(out,'prepare-worker',0,array=f'0-{len(plan["cases"])-1}')
            finalizer=r.submit(out,'prepare-finish',0,dependency=array)
            r.write(out/'preparation_launch.json',dict(array=array,finalizer=finalizer,generated_at=n.now()))
            r.write(out/'results/status.json',dict(generated_at=n.now(),complete=False,producer='running',stage='variant_preparation',lengths=len(plan['cases'])))
            write_index(out)
    else:r.main()
