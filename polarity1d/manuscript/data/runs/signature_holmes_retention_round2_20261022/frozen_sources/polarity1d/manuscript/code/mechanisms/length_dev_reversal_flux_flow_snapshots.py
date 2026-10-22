"""Alpha C cortical-flow components during native relative-20% reversals.

Short restarts recover the exact last native step at each retained snapshot.
No full reversal integrations or changes to retained simulation evidence.
"""
from pathlib import Path
from datetime import datetime, timezone
from dataclasses import replace, asdict
import csv
import gzip
import hashlib
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from length_dev_basal_front_tension_decomposition import solve, gradient
from length_dev_output_layout import dev_output_path

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/'figures/length_mechanisms_20261008'
OUT = REPORT/'dev/reversal_flux_flow'
PARENT = ROOT/'data/runs/length_narrative_memory_20261009'
BASE = ROOT/'data/runs/debelly_corrected_three_test_20261008'
STEM = 'reversal_MCA_flux_cortical_flow_snapshots'
LENGTHS = [21, 25, 31, 35, 39]
PHASES = ['Basal', 'Early cue', 'Later cue', 'Withdrawal', 'Relaxation']
SOURCES = {}
FLOW = [('endpoint','Endpoint forcing','#D55E00','-'),
        ('rho','Rho stress gradient','#0072B2','-'),
        ('membrane','Membrane motion','#929292','-'),
        ('total','Total flow','#222222','--')]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def source(p):
    p = Path(p)
    SOURCES[str(p.relative_to(ROOT))] = sha(p)
    return p


def read(p):
    return json.loads(source(p).read_text())


def output(name):
    return dev_output_path(OUT, name)


def write(name, value):
    output(name).write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def csvout(name, rows):
    with output(name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def decompose(tr, length, p, threshold):
    """Match the native final cortical solve and alpha's displayed-state flux."""
    f = {key: a[-1] for key, a in tr.fields.items()}
    old = {key: a[-2] for key, a in tr.fields.items()}
    dt = float(tr.time[-1]-tr.time[-2])
    dx = length/len(tr.x)
    switch = .5*(1+np.tanh((old['rac'][[0,-1]]-threshold)/p.protrusion_switch_scale))
    boundary = p.rac_to_rho_coupling*p.protrusion_speed*switch/(1+old['tension'][[0,-1]]**2/p.protrusion_tension_scale)
    boundary[1] *= -1
    membrane = (f['displacement']-old['displacement'])/dt
    membrane_edge = np.array([(f[k][0]-old[k][0])/dt for k in ['left_boundary_displacement','right_boundary_displacement']])
    friction = p.linker_friction*np.maximum(f['mca_bound'],1e-10)
    components = dict(endpoint=solve(friction,p.cortical_viscosity,np.zeros_like(tr.x),dx,*boundary),
                      rho=solve(friction,p.cortical_viscosity,p.active_stress*p.rho_to_rac_coupling*np.gradient(old['rho'],dx),dx),
                      membrane=solve(friction,p.cortical_viscosity,friction*membrane,dx),
                      total=f['velocity'].copy())
    flow_error = float(np.max(abs(sum(components[k] for k in ['endpoint','rho','membrane'])-components['total'])))
    assert flow_error < 1e-11, flow_error
    flux = {}
    for label, pool, velocity, edge, diffusion in [
            ('bound','mca_bound',f['velocity'],boundary,p.mca_bound_diffusion),
            ('unbound','mca_unbound',membrane,membrane_edge,p.mca_unbound_diffusion)]:
        y = f[pool]
        face = .5*(velocity[:-1]+velocity[1:])
        donor = np.where(face>=0,y[:-1],y[1:])
        flux[label+'_advection'] = np.r_[edge[0]*y[0],face*donor,edge[1]*y[-1]]
        flux[label+'_diffusion'] = np.r_[0.,-diffusion*np.diff(y)/dx,0.]
    flux['total_flux'] = sum(flux.values())
    divergence = -np.diff(flux['total_flux'])/dx
    conservation = float(abs(divergence.sum()*dx-(flux['total_flux'][0]-flux['total_flux'][-1])))
    assert conservation < 1e-12
    tension_error = float(np.max(abs(p.membrane_stiffness*gradient(f['displacement'],dx,f['left_boundary_displacement'][0],f['right_boundary_displacement'][0])-f['tension'])))
    assert tension_error < 1e-9
    return dict(x=tr.x,face_x=np.linspace(0,length,len(tr.x)+1),fields=f,flow=components,flux=flux,
                membrane_velocity=membrane,membrane_edge_velocity=membrane_edge,cortical_edge_velocity=boundary,
                audit=dict(last_native_step_s=dt,flow_superposition_max_error_um_per_s=flow_error,
                           flux_integral_closure_error=conservation,tension_identity_max_error_pN=tension_error))


def collect():
    manifest = read(PARENT/'manifest.json')
    assert manifest['center_boundary_motion'] is False
    assert read(PARENT/'status.json')['audit_passed']
    assert read(BASE/'results/independent_audit.json')['passed']
    with gzip.open(source(BASE/'results/records.json.gz'),'rt') as f:
        records = {r['case_index']:r for r in json.load(f)['cases']}
    endpoint_path = REPORT/'dev/deBelly_redraft/supporting/data/Figure_alpha_steady_state_mechanism_endpoint_metrics.csv'
    endpoint_meta = read(endpoint_path.parent.parent/'metadata/Figure_alpha_steady_state_mechanism.manifest.json')
    assert sha(endpoint_path) == endpoint_meta['sha256']['supporting/data/'+endpoint_path.name]
    with source(endpoint_path).open() as f:
        thresholds = {int(r['length_um']):float(r['rac_switch_threshold']) for r in csv.DictReader(f)}
    # The immutable native sources that produced the exact retained histories.
    sys.path[:0] = [str(PARENT/'frozen_sources/code'),str(PARENT/'frozen_sources/code/core')]
    import native as n
    import polarity_establishment as establishment
    p = replace(n.db._parameters(),center_boundary_motion=False)
    n.model_for = lambda key:n.ds.DeBelly(p)
    assert p.center_boundary_motion is False
    snapshots, audits, jobs = {}, [], []
    for length in LENGTHS:
        job = next(j for j in manifest['jobs'] if j['case']['length']==length and j['case']['condition']=='fixed_relative')
        case,D = job['case'],float(job['duration'])
        r = records[case['case_index']]
        assert D == r['result']['full_polarity']['endpoint']
        assert D == min(t['duration'] for t in r['tests'] if t['full_polarity_success'])
        assert D > 180
        archive = ROOT/job['source'] if job['source'] else PARENT/f"workers/{manifest['missing_jobs'].index(job):03d}.npz"
        meta = read(archive.with_suffix('.json'))
        assert sha(source(archive)) == meta['archive_sha256']
        checkpoint = source(ROOT/case['checkpoint'])
        assert sha(checkpoint) == case['checkpoint_sha256']
        initial,_ = n.read_checkpoint(checkpoint,'debelly')
        cue = establishment.Patch(case['response_amplitude'],case['width'],side='right',width_mode=case['width_mode'],fraction=.2)
        assert case['response_amplitude']==3.6 and case['width_mode']=='relative'
        with np.load(archive) as z:
            for key in n.DB_NAMES:
                np.testing.assert_array_equal(z['pulse_field_'+key][0],initial[key])
                np.testing.assert_array_equal(z['pulse_field_'+key][-1],z['recovery_field_'+key][0])
            tail = z['recovery_time']>=9500
            rac = z['recovery_field_rac'][tail]
            orientation = rac@(2*z['x']/length-1)/rac.sum(axis=1)
            strength_error = float(np.max(abs(orientation/abs(case['initial_orientation'])-1)))
            assert orientation.min()>.05 and np.ptp(rac,axis=1).min()>=.05 and np.ptp(orientation)<.005 and strength_error<=.01
            assert meta['metrics']['full_polarity_success']
            for row,(phase,prefix,target) in enumerate(zip(PHASES,['pulse']*4+['recovery'],[0.,60.,180.,D,120.])):
                if row==0:
                    tr = n.run('debelly',length,initial,n.ConstantInput(0),.5,interval=.5,n=case['cells'])
                    replay_error = max(float(abs(tr.state[k]-initial[k]).max()) for k in n.DB_NAMES)
                    assert replay_error < 1e-5
                    actual,offset = .5,0.
                else:
                    times = z[prefix+'_time']
                    i = int(np.argmin(abs(times-target)))
                    assert i>0 and abs(times[i]-target)<1e-8
                    state = {key:z[prefix+'_field_'+key][i-1].copy() for key in n.DB_NAMES}
                    dt = float(times[i]-times[i-1])
                    tr = n.run('debelly',length,state,cue if prefix=='pulse' else n.ConstantInput(0),dt,interval=.5,n=case['cells'])
                    replay_error = max(float(abs(tr.state[k]-z[prefix+'_field_'+k][i]).max()) for k in n.DB_NAMES)
                    assert replay_error < 1e-10,(length,phase,replay_error)
                    actual,offset = float(times[i]),D if prefix=='recovery' else 0.
                s = decompose(tr,length,p,thresholds[length])
                s.update(length=length,row=row,phase=phase,time_after_onset_s=actual+offset,
                         time_since_withdrawal_s=actual+offset-D,duration=D)
                s['audit'].update(length_um=length,phase=phase,time_after_onset_s=actual+offset,
                                  maximum_native_replay_field_error=replay_error,basal_continuation=row==0)
                snapshots[row,length] = s
                audits.append(s['audit'])
        jobs.append(dict(length_um=length,duration_s=D,archive=str(archive.relative_to(ROOT)),checkpoint=case['checkpoint'],
                         checkpoint_sha256=case['checkpoint_sha256'],exact_checkpoint_and_phase_continuity=True,
                         shortest_tested_scalar_strength_selection_verified=True,independent_scalar_gate_pass=True,
                         maximum_strength_error=strength_error,native_field_tail_pass=meta['metrics']['field_tail_pass']))
        print(length,'all five snapshots verified',flush=True)
    # Register all imported frozen Python dependencies and check their original hashes.
    for module in list(sys.modules.values()):
        path = getattr(module,'__file__',None)
        if path and Path(path).suffix=='.py' and PARENT/'frozen_sources' in Path(path).parents:
            path = source(Path(path))
            relative = str(path.relative_to(PARENT/'frozen_sources'))
            if relative in manifest['source_sha256']:
                assert sha(path)==manifest['source_sha256'][relative]
    return snapshots,audits,jobs,asdict(p)


def render(snapshots):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none',
                         'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6})
    fig = plt.figure(figsize=(16.5,11.8))
    gs = fig.add_gridspec(5,5,left=.095,right=.988,bottom=.085,top=.865,hspace=.48,wspace=.23)
    values = [v for s in snapshots.values() for v in s['flow'].values()]
    values += [s['cortical_edge_velocity'] for s in snapshots.values()]
    lo,hi = min(float(v.min()) for v in values),max(float(v.max()) for v in values)
    pad = .08*max(hi-lo,1e-12)
    limits = {'flow':[lo-pad,hi+pad]}
    for row,phase in enumerate(PHASES):
        for col,length in enumerate(LENGTHS):
            s = snapshots[row,length]
            ax = fig.add_subplot(gs[row,col])
            for key,label,color,ls in FLOW:
                edge = s['cortical_edge_velocity'] if key in ['endpoint','total'] else [0.,0.]
                ax.plot(np.r_[0,s['x'],length],np.r_[edge[0],s['flow'][key],edge[1]],
                        color=color,ls=ls,lw=1.6 if key=='total' else 1.4)
            ax.set_xlim(0,length)
            ax.set_ylim(limits['flow'])
            ax.axhline(0,color='#BBBBBB',lw=.55,zorder=0)
            ax.set_xticks([0,length/2,length],[f'{v:g}' for v in [0,length/2,length]])
            ax.yaxis.set_major_locator(MaxNLocator(3))
            ax.tick_params(labelsize=8,pad=2,length=2.5)
            ax.grid(alpha=.1,lw=.5)
            if col:
                ax.tick_params(labelleft=False)
            if col==0:
                ax.set_ylabel('Cortical velocity\n(µm/s)',fontsize=9)
            ax.set_xlabel('x (µm)',fontsize=9,labelpad=1)
            time = s['time_after_onset_s']
            title = 'Basal' if row==0 else f't = {time:.2f} s' if row>=3 else f't = {time:g} s'
            ax.set_title(title,fontsize=10,pad=5)
            if row==0:
                pos = ax.get_position()
                fig.text((pos.x0+pos.x1)/2,.902,f'L = {length} µm',ha='center',fontsize=13)
            if col==0:
                pos = gs[row,0].get_position(fig)
                text = ['Basal','Cue: 60 s','Cue: 180 s','Cue withdrawal','Withdrawal + 120 s'][row]
                fig.text(.016,(pos.y0+pos.y1)/2,text,rotation=90,ha='center',va='center',fontsize=12)
    fig.suptitle('Cortical-flow components during reversal',fontsize=19,y=.995)
    fig.text(.54,.957,'Relative-20% cues · positive velocity points rightward · shared vertical scale throughout',ha='center',fontsize=10)
    handles = [Line2D([],[],color=color,ls=ls,lw=1.6,label=label) for key,label,color,ls in FLOW]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.944),ncol=len(handles),frameon=False,fontsize=10)
    fig.text(.095,.04,'Left = old front; right = challenged front. Withdrawal rows follow each cell’s own cue duration; printed times are from cue onset.',fontsize=9)
    fig.text(.095,.02,'Endpoint forcing + Rho stress gradient + membrane motion = total cortical flow. Basal profiles follow one 0.5-s no-cue step.',fontsize=9)
    for ext in ['.png','.pdf','.svg']:
        fig.savefig(output(STEM+ext),dpi=220,facecolor='white')
    plt.close(fig)
    return limits


def main():
    generated = datetime.now(timezone.utc).isoformat()
    for path in [Path(__file__),Path(__file__).with_name('length_dev_basal_front_tension_decomposition.py'),Path(__file__).with_name('length_dev_output_layout.py')]:
        source(path)
    report_hashes = {name:sha(REPORT/name) for name in ['index.html','length_mechanisms.pdf']}
    snapshots,audits,jobs,parameters = collect()
    arrays,flux_rows,flow_rows,state_rows = {},[],[],[]
    for (row,length),s in snapshots.items():
        ident = f'L{length}_snapshot{row}'
        arrays.update({ident+'_x':s['x'],ident+'_face_x':s['face_x'],ident+'_time_after_onset_s':np.array(s['time_after_onset_s'])})
        arrays.update({ident+'_flow_'+k:v for k,v in s['flow'].items()})
        arrays.update({ident+'_flux_'+k:v for k,v in s['flux'].items()})
        arrays.update({ident+'_field_'+k:v for k,v in s['fields'].items()})
        arrays[ident+'_membrane_velocity'] = s['membrane_velocity']
        arrays[ident+'_membrane_edge_velocity'] = s['membrane_edge_velocity']
        arrays[ident+'_cortical_edge_velocity'] = s['cortical_edge_velocity']
        common = dict(length_um=length,snapshot_index=row,phase=s['phase'],time_after_onset_s=s['time_after_onset_s'],time_since_withdrawal_s=s['time_since_withdrawal_s'])
        for i,x in enumerate(s['face_x']):
            flux_rows.append(dict(**common,face_index=i,x_um=float(x),**{k:float(v[i]) for k,v in s['flux'].items()}))
        for i,x in enumerate(s['x']):
            flow_rows.append(dict(**common,cell_index=i,x_um=float(x),**{k:float(v[i]) for k,v in s['flow'].items()}))
            state_rows.append(dict(**common,cell_index=i,x_um=float(x),**{k:float(s['fields'][k][i]) for k in ['rac','rho','mca_bound','mca_unbound','tension','displacement']},membrane_velocity=float(s['membrane_velocity'][i])))
    np.savez_compressed(output(STEM+'_fields.npz'),**arrays)
    csvout(STEM+'_fluxes.csv',flux_rows)
    csvout(STEM+'_flows.csv',flow_rows)
    csvout(STEM+'_profiles.csv',state_rows)
    limits = render(snapshots)
    caption = ('Spatial cortical-flow components during corrected native deBelly relative-20% reversals. '
        'Columns are L=21,25,31,35,39 µm. The five rows show basal state, 60 s and 180 s under the cue, each length’s own cue withdrawal, '
        'and 120 s of basal recovery after withdrawal. Each panel shows endpoint forcing, direct Rho stress-gradient forcing, membrane-motion drag and total cortical velocity, following Figure alpha C. '
        'All lengths remain under cue at 60 and 180 s. Printed withdrawal/recovery times are absolute time from cue onset, '
        'not a common elapsed time across lengths. Cue durations are '+', '.join(f"{j['length_um']} µm: {j['duration_s']:.6f} s" for j in jobs)+'. '
        'Every trial starts at its exact qualifying basal checkpoint and uses the shortest tested canonical scalar-1%-strength-passing relative cue, peak 3.6 and support 20% of nominal length. '
        'No state is reset at withdrawal. Basal transport diagnostics follow one native 0.5-s no-cue step, as in alpha B–C; basal continuation changes are audited. '
        'Other displayed fields reproduce the original retained snapshots. Short restarts over only the preceding saved output interval recover every native step and are independently checked against all eleven retained fields. '
        'Cortical velocity is decomposed by independently assembled Dirichlet solves into prescribed endpoint forcing, direct Rho stress-gradient forcing and membrane-motion drag. '
        'The final-step inputs respect the native staggering: previous Rho and endpoint Rac/tension, updated bound-MCA friction and actual final-step membrane velocity. Their sum reconstructs the saved total flow. '
        'These are instantaneous source superpositions, not independently evolved causal ablations. Positive velocity means rightward; left is the old front, right the challenged front. '
        'Curves show the 201 native cell centres plus prescribed endpoint values. Endpoint forcing carries the native endpoint commands; Rho and membrane contributions have zero endpoint velocities. '
        'One common vertical scale is used across all 25 panels. Supporting data also retain native MCA flux components and bound/unbound profiles. '
        'Positions are nominal reference coordinates in µm. Independent scalar gates and native-field tail flags remain distinct and are exported. Original simulation evidence and the main report are unchanged.\n')
    output(STEM+'_caption.txt').write_text(caption)
    audit = dict(generated_at=generated,passed=True,snapshots=len(audits),lengths_um=LENGTHS,jobs=jobs,checks=audits,
                 parameters=parameters,common_axis_limits=limits,
                 maximum_flow_superposition_error_um_per_s=max(a['flow_superposition_max_error_um_per_s'] for a in audits),
                 maximum_nonbasal_replay_field_error=max(a['maximum_native_replay_field_error'] for a in audits if not a['basal_continuation']),
                 main_report_unchanged=all(sha(REPORT/k)==v for k,v in report_hashes.items()))
    assert audit['main_report_unchanged']
    write(STEM+'_audit.json',audit)
    viewer = ('<!doctype html><html><head><meta charset="utf-8"><title>Cortical-flow components during reversal</title>'
        '<style>body{font:16px system-ui;max-width:1700px;margin:24px auto;padding:0 20px}img{width:100%}p{line-height:1.5}</style></head><body>'
        '<h1>Cortical-flow components during reversal</h1><p>Generated '+generated+' · producer finished · all 25 snapshots verified.</p>'
        '<p><a href="../../'+STEM+'.png">Full PNG</a> · <a href="../vectors/'+STEM+'.pdf">Vector PDF</a> · '
        '<a href="../data/'+STEM+'_fields.npz">Fields</a> · <a href="../metadata/'+STEM+'_audit.json">Audit</a></p>'
        '<img src="../../'+STEM+'.png" alt="Endpoint, Rho and membrane contributions to cortical velocity at five lengths and five snapshots"><p>'+caption+'</p></body></html>')
    output(STEM+'.html').write_text(viewer)
    artifacts = [output(STEM+ext) for ext in ['.png','.pdf','.svg','.html','_fields.npz','_fluxes.csv','_flows.csv','_profiles.csv','_caption.txt','_audit.json']]
    write(STEM+'.manifest.json',dict(generated_at=generated,passed=True,artifact_root='../..',lengths_um=LENGTHS,
          source_sha256=SOURCES,sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts},
          reproduction='python polarity1d/manuscript/code/mechanisms/length_dev_reversal_flux_flow_snapshots.py',
          native_histories_unchanged=True,main_report_sha256=report_hashes))
    print(output(STEM+'.png'),flush=True)


if __name__=='__main__':
    main()
