"""Native right-end protrusion drive and cortical flow; retained fields only."""
from pathlib import Path
from datetime import datetime, timezone
from dataclasses import replace
import csv, gzip, hashlib, importlib.util, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_output_layout as layout

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'data/runs/length_narrative_memory_20261009'
BASE = ROOT/'data/runs/debelly_corrected_three_test_20261008'
OUT = ROOT/'figures/length_mechanisms_20261008/dev/stimulated_endpoint'
LENGTHS = [18,19,21,23,25,27,29,31,33,35,37,39]
GEOS = ['fixed_physical','fixed_relative','fixed_slope']
TITLES = ['Fixed patch: 3 µm','Relative patch: 20%','Constant slope: 0.24']
STEM = 'stimulated_endpoint'
NAMES = ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension',
         'raw_left_displacement','raw_right_displacement',
         'left_boundary_displacement','right_boundary_displacement']
SOURCES = {}
ARTIFACTS = {}

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def register(p):
    SOURCES[str(p.relative_to(ROOT))] = sha(p)
    return p

def read(p): return json.loads(register(p).read_text())
def output(name): return layout.dev_output_path(OUT,name)
def keep(p): ARTIFACTS[str(p.relative_to(OUT))] = sha(p)
def write(name,value):
    p=output(name);p.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n');keep(p)
def csvwrite(name,rows):
    p=output(name)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    keep(p)

def save(fig,name):
    for ext in ['png','pdf','svg']:
        p=output(name+'.'+ext);fig.savefig(p,dpi=210,facecolor='white');keep(p)
    plt.close(fig)

def main():
    register(Path(__file__));register(Path(layout.__file__))
    m=read(RUN/'manifest.json')
    assert m['center_boundary_motion'] is False
    assert read(RUN/'status.json')['audit_passed']
    assert read(BASE/'results/independent_audit.json')['passed']
    assert read(ROOT/'data/runs/length_mechanisms_20261008/narrative_20261009/memory/independent_audit.json')['passed']
    # Use the actual frozen model law, independent of subsequent live development.
    frozen=register(RUN/'frozen_sources/code/core/models/debelly.py')
    assert sha(frozen)==m['source_sha256']['code/core/models/debelly.py']
    sys.path.insert(0,str(RUN/'frozen_sources/code/core/src'))
    spec=importlib.util.spec_from_file_location('endpoint_frozen_debelly',frozen)
    ds=importlib.util.module_from_spec(spec);sys.modules[spec.name]=ds;spec.loader.exec_module(ds)
    p=replace(ds._parameters(),center_boundary_motion=False)
    rac0,rho0=ds._homogeneous_rac_rho(p,.5)
    threshold=1.3*rac0
    records_path=register(BASE/'results/records.json.gz')
    with gzip.open(records_path,'rt') as f:
        records={r['case_index']:r for r in json.load(f)['cases']}
    traces={};samples=[];summary=[];audits=[];arrays={}
    for j in m['jobs']:
        c=j['case'];L=int(c['length']);geo=c['condition'];D=float(j['duration'])
        assert L in LENGTHS and geo in GEOS and not c['center_boundary_motion']
        r=records[c['case_index']]
        assert D==r['result']['full_polarity']['endpoint']
        assert D==min(t['duration'] for t in r['tests'] if t['full_polarity_success'])
        path=register(ROOT/j['source'] if j['source'] else RUN/f"workers/{m['missing_jobs'].index(j):03d}.npz")
        meta=read(path.with_suffix('.json'));assert sha(path)==meta['archive_sha256']
        cp=register(ROOT/c['checkpoint']);assert sha(cp)==c['checkpoint_sha256']
        phases={};tail_max=0.;recurrence_checks=0;recurrence_error=0.
        with np.load(cp) as start, np.load(path) as z:
            x=z['x'];np.testing.assert_allclose(x,(np.arange(c['cells'])+.5)*L/c['cells'],rtol=0,atol=1e-12)
            tp=z['pulse_time'];tr=z['recovery_time'];tail=tr>=9500
            assert tp[0]==tr[0]==0 and abs(tp[-1]-D)<1e-9 and tr[-1]==10000
            assert np.all(np.diff(tp)>0) and np.all(np.diff(tr)>0)
            for name in NAMES:
                a=z['pulse_field_'+name];b=z['recovery_field_'+name]
                assert np.isfinite(a).all() and np.isfinite(b).all()
                np.testing.assert_array_equal(a[0],start['state_'+name])
                np.testing.assert_array_equal(b[0],a[-1])
                if name in c['fixed_scales']:
                    y=b[tail]
                    if name=='displacement': y=y-y.mean(axis=1)[:,None]
                    scale=c['fixed_scales'][name]
                    tail_max=max(tail_max,float(np.ptp(y,axis=0).max()/scale),
                                 float(abs(y[-1]-y[-2]).max()/np.diff(tr)[-1]*500/scale))
            separation=z['recovery_field_boundary_separation'][tail]
            tail_max=max(tail_max,float(np.ptp(separation,axis=0).max()/L),
                         float(abs(separation[-1]-separation[-2]).max()/np.diff(tr)[-1]*500/L))
            rac=z['recovery_field_rac'][tail]
            orientation=rac@(2*x/L-1)/rac.sum(axis=1)
            strength_error=float(abs(orientation/abs(c['initial_orientation'])-1).max())
            assert orientation.min()>.05 and np.ptp(rac,axis=1).min()>=.05
            assert np.ptp(orientation)<.005 and strength_error<=.01
            flag=bool(tail_max<=.001)
            assert flag==meta['metrics']['field_tail_pass']==j['expected']['field_tail_pass']
            for phase,t in [('pulse',tp),('recovery',tr)]:
                R=z[phase+'_field_rac'][:,-1];T=z[phase+'_field_tension'][:,-1]
                v=z[phase+'_field_velocity'][:,-1]
                u=z[phase+'_field_right_boundary_displacement'][:,-1]
                switch=.5*(1+np.tanh((R-threshold)/p.protrusion_switch_scale))
                inhibition=1/(1+T*T/p.protrusion_tension_scale)
                drive=-p.rac_to_rho_coupling*p.protrusion_speed*switch*inhibition
                native=-p.rac_to_rho_coupling*p.protrusion_speed*ds.smooth_threshold(R,p.protrusion_switch_scale,threshold)/(1+T*T/p.protrusion_tension_scale)
                np.testing.assert_allclose(drive,native,rtol=1e-14,atol=1e-16)
                assert drive.min()>=0 and drive.max()<=abs(p.protrusion_speed)
                # Exactly consecutive half-second states verify the actual command
                # timing against the independent stored endpoint recurrence.
                adjacent=np.isclose(np.diff(t),.5,rtol=0,atol=1e-10)
                if adjacent.any():
                    inferred=(u[1:]-p.boundary_relaxation*u[:-1])/np.diff(t)
                    err=float(abs(inferred[adjacent]-drive[:-1][adjacent]).max())
                    assert err<1e-12
                    recurrence_error=max(recurrence_error,err);recurrence_checks+=int(adjacent.sum())
                if geo=='fixed_slope': stimulus=float(c['gradient_slope']*x[-1])
                else:
                    width=.2*L if c['width_mode']=='relative' else c['width']
                    stimulus=float(c['response_amplitude']*.5*(1+np.cos(np.pi*(L-x[-1])/width)))
                applied=np.full(len(t),stimulus if phase=='pulse' else 0.)
                phases[phase]=dict(time=t.copy(),drive=drive,velocity=v.copy(),stimulus=applied,
                                   rac=R.copy(),tension=T.copy(),switch=switch,inhibition=inhibition)
        traces[geo,L]=dict(phases=phases,duration=D,flag=flag)
        prefix=f'{geo}_L{L}'
        arrays[prefix+'_cue_duration_s']=np.array(D)
        arrays[prefix+'_field_tail_pass']=np.array(flag)
        arrays[prefix+'_nominal_distance_from_right_endpoint_um']=np.array(L/(2*c['cells']))
        for phase,a in phases.items():
            for key,y in a.items():arrays[prefix+'_'+phase+'_'+key]=y
            for i,t in enumerate(a['time']):
                samples.append(dict(condition=geo,length_um=L,phase='cue' if phase=='pulse' else 'basal_recovery',
                    time_from_onset_s=float(t+(D if phase=='recovery' else 0)),
                    time_from_withdrawal_s=float(t if phase=='recovery' else t-D),
                    applied_extra_Rac_stimulus=float(a['stimulus'][i]),right_protrusion_drive_um_s=float(a['drive'][i]),
                    right_adjacent_cortical_velocity_um_s=float(a['velocity'][i]),right_rac=float(a['rac'][i]),
                    right_tension_pN=float(a['tension'][i]),rac_switch=float(a['switch'][i]),tension_inhibition=float(a['inhibition'][i])))
        a=phases['pulse'];b=phases['recovery']
        row=dict(condition=geo,length_um=L,cue_duration_s=D,nominal_distance_from_right_endpoint_um=L/(2*c['cells']),
                 initial_drive_um_s=float(a['drive'][0]),withdrawal_drive_um_s=float(a['drive'][-1]),
                 final_drive_um_s=float(b['drive'][-1]),initial_velocity_um_s=float(a['velocity'][0]),
                 withdrawal_velocity_um_s=float(a['velocity'][-1]),final_velocity_um_s=float(b['velocity'][-1]),
                 peak_cue_drive_um_s=float(a['drive'].max()),peak_cue_velocity_um_s=float(a['velocity'].max()),
                 native_field_tail_pass=flag,archive=str(path.relative_to(ROOT)))
        summary.append(row)
        audits.append(dict(condition=geo,length_um=L,passed=True,archive_sha256=sha(path),checkpoint_sha256=c['checkpoint_sha256'],
            exact_initial_fields=len(NAMES),exact_phase_continuity_fields=len(NAMES),shortest_tested_strength_pulse_verified=True,
            independent_scalar_strength_relative_error=strength_error,independent_field_tail_max_scaled=tail_max,
            field_tail_pass=flag,half_second_command_recurrence_checks=recurrence_checks,maximum_command_recurrence_error_um_s=recurrence_error))
        print(geo,L,'audited',flush=True)
    assert set(traces)=={(g,L) for g in GEOS for L in LENGTHS}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    colors={L:plt.get_cmap('turbo')(.06+.88*i/11) for i,L in enumerate(LENGTHS)}
    def legend(fig,phase=True):
        handles=[Line2D([],[],color=colors[L],lw=2,label=str(L)) for L in LENGTHS]
        if phase:
            handles += [Line2D([],[],color='#444',lw=1.5,label='Cue on'),Line2D([],[],color='#444',lw=1.5,ls='--',label='Basal relaxation'),
                Line2D([],[],color='#444',marker='o',ls='',mfc='#444',label='Final fields stable'),
                Line2D([],[],color='#444',marker='o',ls='',mfc='white',label='Final fields changing')]
        fig.legend(handles=handles,loc='outside lower center',ncol=8 if phase else 12,frameon=False,fontsize=8.5,
                   title='Nominal cell length (µm)'+(' · circles mark cue withdrawal' if phase else ''))
    fig,axs=plt.subplots(3,3,figsize=(14,9),sharex=True,sharey='row',layout='constrained',gridspec_kw={'height_ratios':[.65,1,1]})
    keys=['stimulus','drive','velocity']
    labels=['Applied extra Rac stimulus\n(native input units)','Right-end protrusion drive\n(µm/s)','Cortical velocity beneath\nright endpoint (µm/s)']
    for col,(geo,title) in enumerate(zip(GEOS,TITLES)):
        axs[0,col].set_title(title,pad=12)
        for row,key in enumerate(keys):
            ax=axs[row,col];ax.text(-.12,1.04,'ABCDEFGHI'[row*3+col],transform=ax.transAxes,fontweight='bold',fontsize=12)
            for L in LENGTHS:
                case=traces[geo,L];a=case['phases']['pulse'];b=case['phases']['recovery'];D=case['duration']
                if key=='stimulus':
                    ax.plot([0,D,D],[a[key][0],a[key][0],0],color=colors[L],lw=1.15,alpha=.8)
                else:
                    ax.plot(a['time'],a[key],color=colors[L],lw=1.25)
                    ax.plot(D+b['time'],b[key],color=colors[L],lw=1.25,ls='--')
                    ax.plot(D,a[key][-1],'o',ms=4.5,mec=colors[L],mfc=colors[L] if case['flag'] else 'white',mew=1,zorder=4)
            ax.set_xlim(0,1200);ax.set_xticks([0,300,600,900,1200]);ax.grid(alpha=.13);ax.axhline(0,color='#999',lw=.6)
            if col==0:ax.set_ylabel(labels[row])
            if row==2:ax.set_xlabel('Time from reversal-cue onset (s)')
    legend(fig);save(fig,STEM+'_cue_and_relaxation')
    fig,axs=plt.subplots(2,3,figsize=(14,7.2),sharex=True,sharey='row',layout='constrained')
    for col,(geo,title) in enumerate(zip(GEOS,TITLES)):
        axs[0,col].set_title(title,pad=12)
        for row,key in enumerate(['drive','velocity']):
            ax=axs[row,col];ax.text(-.12,1.04,'ABCDEF'[row*3+col],transform=ax.transAxes,fontweight='bold',fontsize=12)
            for L in LENGTHS:
                case=traces[geo,L];b=case['phases']['recovery']
                ax.plot(b['time'],b[key],color=colors[L],lw=1.25)
            ax.set_xscale('symlog',linthresh=10,linscale=1);ax.set_xlim(0,10000)
            ax.set_xticks([0,10,100,1000,10000],['0','10','100','1,000','10,000']);ax.grid(alpha=.13);ax.axhline(0,color='#999',lw=.6)
            if col==0:ax.set_ylabel(labels[row+1])
            if row==1:ax.set_xlabel('Time after cue withdrawal (s)')
    legend(fig,False);save(fig,STEM+'_full_relaxation')
    csvwrite(STEM+'_samples.csv',samples);csvwrite(STEM+'_summary.csv',summary)
    path=output(STEM+'_traces.npz');np.savez_compressed(path,**arrays);keep(path)
    caption=('Corrected native deBelly reversal at all twelve eligible lengths, using each geometry and length’s shortest tested canonical 1%-strength-passing finite pulse from its exact basal checkpoint. Columns: fixed support 3 µm, relative support 20%, constant physical slope 0.24. The cue-and-relaxation figure shows the applied extra Rac stimulus at the last native cell, the right-end state-evaluated protrusion drive, and recorded cortical velocity in that cell. Solid lines: cue on; dashed: basal relaxation. Circles mark withdrawal; filled/open markers report the independently recomputed final native-field tail flag, not stability at withdrawal. Positive velocity points toward the stimulated right endpoint. The native sampling point is x=L−L/(2×201), 0.045–0.097 µm beneath the nominal endpoint. The protrusion drive is +0.15 × Rac switch / (1+T²/10), with Rac switch = 0.5 × (1+tanh((Rac−1.3×resting Rac)/0.01)), using the frozen native parameters. It is calculated from saved rightmost Rac and tension; the staggered solver uses this drive in the NEXT mechanical step. Recorded cortical velocity belongs to the completed step, so the two curves have a native-step timing offset (up to 0.5 s) as well as a different spatial location. The exact boundary cortical velocity is prescribed by the applied protrusion command; the plotted cell-centred flow is not an independent membrane-endpoint speed. Independent endpoint commands and protrusion remain native throughout, including after external cue withdrawal. Relaxation means basal cue-free evolution, not forced zero protrusion. The full-relaxation companion aligns all curves to withdrawal and retains the entire 10000-s recovery on an axis linear from 0–10 s and logarithmic thereafter. All native samples, Rac/tension factors, summaries, hashes and audits are exported. No new integration or change to the main report.')
    for name in [STEM+'_cue_and_relaxation',STEM+'_full_relaxation']:
        path=output(name+'_caption.txt');path.write_text(caption+'\n');keep(path)
    write(STEM+'_audit.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,trajectories=36,
        new_integrations=0,checks=audits,resting_rac=rac0,resting_rho=rho0,rac_threshold=threshold,
        command_timing='State-evaluated next-step drive; recorded velocity is from the completed staggered step.'))
    write(STEM+'.manifest.json',dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',
        independent_audit_passed=True,trajectories=36,new_integrations=0,source_sha256=SOURCES,sha256=dict(ARTIFACTS),
        legend=caption,reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))]))
    print(output(STEM+'_cue_and_relaxation.png'))

if __name__=='__main__':main()
