"""Audited absolute and baseline-subtracted local Rac budgets and input factors."""
from pathlib import Path
from datetime import datetime, timezone
import csv, gzip, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_output_layout as layout

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'data/runs/debelly_old_pole_rac_budget_20261010'
BASE=ROOT/'data/runs/debelly_corrected_three_test_20261008'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/native_reversal'
STEM='old_pole_Rac_budget'
SOURCES={}
TERMS=[('basal_activation','#0072B2','Rho-suppressed basal activation'),
       ('mca_activation','#D55E00','Rho-suppressed MCA activation'),
       ('deactivation','#7B3294','Deactivation'),('diffusion','#009E73','Diffusion'),
       ('net_rate','#222222','Net Rac rate')]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def register(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return p
def read(p):return json.loads(register(p).read_text())
def output(name):return layout.dev_output_path(OUT,name)
def write(name,value):output(name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
def csvwrite(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def measure():
    m=read(RUN/'manifest.json');prov=read(RUN/'provenance.json')
    assert sha(RUN/'manifest.json')==prov['manifest_sha256'] and m['center_boundary_motion'] is False
    original=read(BASE/'provenance.json')
    for name,digest in prov['sources'].items():
        assert sha(RUN/'frozen_sources'/name)==digest
        if name.startswith('code/core/'):assert original['sources'][name]==digest
    path=register(BASE/'results/records.json.gz')
    with gzip.open(path,'rt') as f:canonical={r['case_index']:r for r in json.load(f)['cases']}
    traces=[];checks=[];summaries=[];samples=[];arrays={}
    for job in m['jobs']:
        L=job['length_um'];D=job['duration'];c=job['case'];rec=canonical[c['case_index']]
        assert D==rec['result']['full_polarity']['endpoint']==min(r['duration'] for r in rec['tests'] if r['full_polarity_success'])
        prefix=RUN/'workers'/f"{job['slot']:03d}";meta=read(prefix.with_suffix('.json'))
        archive=register(prefix.with_suffix('.npz'));assert sha(archive)==meta['archive_sha256'] and meta['passed']
        source=register(Path(job['source']));assert sha(source)==job['source_sha256']
        assert sha(register(Path(job['checkpoint'])))==job['checkpoint_sha256']
        with np.load(archive) as z:data={k:z[k].copy() for k in z.files}
        p=meta['parameters'];rho=data['old_rho'];b=data['updated_bound_mca'];r=data['rac_after']
        q=1/(1+rho*rho);switch=.5*(1-np.tanh((b-p['mca_rac_threshold'])/p['mca_switch_scale']))
        basal=p['rac_deactivation']*p['rac_basal_activation']*q
        mca=p['rac_deactivation']*p['rac_mechanical_feedback']*p['rho_to_rac_coupling']*switch*q
        loss=-p['rac_deactivation']*r
        for name,expected in [('rho_suppression_factor',q),('mca_switch',switch),('basal_activation',basal),('mca_activation',mca),('deactivation',loss)]:
            np.testing.assert_allclose(data[name],expected,rtol=0,atol=1e-15)
        assert np.max(abs(data['local_cue_activation']))==0
        total=basal+mca+loss+data['diffusion']
        actual=(r-data['rac_before'])/data['step_s']
        np.testing.assert_allclose(actual,data['net_rate'],rtol=0,atol=1e-15)
        residual=float(abs(actual-total).max());assert residual<1e-12
        np.testing.assert_allclose(data['rac_before'][1:],data['rac_after'][:-1],rtol=0,atol=1e-15)
        np.testing.assert_allclose(np.cumsum(data['step_s']),data['time_s'],rtol=0,atol=1e-8)
        assert data['rac_before'][0]==data['initial_old_pole_rac']
        integrated=float(np.sum(total*data['step_s']));change=float(r[-1]-data['initial_old_pole_rac'])
        integral_error=abs(integrated-change);assert integral_error<1e-10
        # Reconstruct first-cell Neumann diffusion directly from the original,
        # independently retained Rac profiles at every common sample time.
        differror=0.;count=0;dx=L/c['cells']
        with np.load(source) as z:
            for phase in ['pulse','recovery']:
                t=z[phase+'_time']+(D if phase=='recovery' else 0.)
                mask=(t>0)&(t<=job['limit_s']+1e-9);times=t[mask]
                nearest=np.argmin(abs(data['time_s'][:,None]-times[None,:]),axis=0)
                assert np.max(abs(data['time_s'][nearest]-times))<1e-8
                rac=z[phase+'_field_rac'][mask]
                expected=p['rac_deactivation']*p['rac_diffusion']*(rac[:,1]-rac[:,0])/dx**2
                err=float(abs(data['diffusion'][nearest]-expected).max())
                differror=max(differror,err);count+=len(times)
        assert differror<1e-12
        check=dict(length_um=L,passed=True,maximum_budget_residual=residual,maximum_original_profile_diffusion_residual=differror,
                   original_diffusion_samples_checked=count,integrated_budget_error=integral_error,zero_local_cue_verified=True,
                   exact_native_replay_field_error=meta['maximum_native_replay_field_error'],
                   native_replay_samples_checked=meta['original_sample_times_compared'],
                   canonical_shortest_strength_pulse_verified=True,actual_Rac_steps=len(r))
        checks.append(check);traces.append((job,data))
        summary=dict(length_um=L,cue_duration_s=D,initial_old_pole_rac=float(data['initial_old_pole_rac']),
                     old_pole_rac_at_1200_s=float(r[-1]),total_Rac_change=change,
                     initial_rate_reference_time_s=float(data['time_s'][0]),initial_rate_reference_dt_s=float(data['step_s'][0]),
                     maximum_budget_residual=residual,source_archive=str(source.relative_to(ROOT)))
        for term,_,_ in TERMS:
            summary[term+'_initial']=float(data[term][0]);summary[term+'_final']=float(data[term][-1])
            summary[term+'_time_integral']=float(np.sum(data[term]*data['step_s']))
        summaries.append(summary)
        for i in range(len(r)):
            row=dict(length_um=L,time_from_onset_s=float(data['time_s'][i]),phase=str(data['phase'][i]),step_s=float(data['step_s'][i]))
            for name in ['basal_activation','mca_activation','deactivation','diffusion','net_rate','component_sum',
                         'rho_suppression_factor','old_rho','updated_bound_mca','mca_switch','local_cue_activation','rac_before','rac_after']:
                row[name]=float(data[name][i])
            samples.append(row)
        arrays.update({f'L{L}_{k}':v for k,v in data.items()})
    assert [j['length_um'] for j,d in traces]==[23,31,39]
    return traces,checks,summaries,samples,arrays

def save(fig,stem,caption,hashes):
    for ext in ['png','pdf','svg']:
        p=output(stem+'.'+ext);fig.savefig(p,dpi=220,facecolor='white');hashes[str(p.relative_to(OUT))]=sha(p)
    plt.close(fig)
    p=output(stem+'_caption.txt');p.write_text(caption+'\n');hashes[str(p.relative_to(OUT))]=sha(p)

def main():
    register(Path(__file__));register(Path(layout.__file__))
    traces,checks,summaries,samples,arrays=measure()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axs=plt.subplots(2,3,figsize=(13.8,7.4),layout='constrained',sharex=True,sharey='row')
    for col,(job,d) in enumerate(traces):
        t=d['time_s'];D=job['duration'];L=job['length_um']
        for row in range(2):
            ax=axs[row,col]
            for term,color,label in TERMS:
                y=d[term] if row==0 else d[term]-d[term][0]
                ax.plot(t,1000*y,color=color,lw=1.7 if term!='net_rate' else 2.1,ls='-' if term!='net_rate' else '--',label=label,zorder=4 if term=='net_rate' else 2)
            ax.axhline(0,color='#999999',lw=.8,zorder=1)
            ax.axvline(D,color='#777777',ls=':',lw=1.2)
            ax.set(xlim=(0,1200),xticks=[0,300,600,900,1200]);ax.grid(alpha=.12)
            letter=chr(65+3*row+col)
            title=f'{letter}   L = {L} µm' if row==0 else f'{letter}   Change from initial rate'
            ax.set_title(title,loc='left',fontweight='bold')
            if row==0:ax.text(D+12,.95,'Cue off',rotation=90,transform=ax.get_xaxis_transform(),va='top',fontsize=8,color='#555555')
            else:ax.set_xlabel('Time from reversal-cue onset (s)')
    axs[0,0].set_ylabel('Rac rate contribution\n(10⁻³ native concentration/s)')
    axs[1,0].set_ylabel('Change in rate contribution\n(10⁻³ native concentration/s)')
    fig.suptitle('Contributions to Rac at the old pole · native relative-20% reversal',fontsize=15,fontweight='bold')
    handles=[Line2D([],[],color=color,lw=2,ls='--' if term=='net_rate' else '-',label=label) for term,color,label in TERMS]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False,fontsize=9,
               title='Both activation terms include local Rho suppression · old-pole cue input = 0')
    caption=('Contributions to the local Rac rate at the old left pole during native corrected deBelly relative-width-20% reversals at L=23,31,39 µm. Each trial uses the shortest tested canonical scalar-1%-strength-passing finite cue from its exact basal checkpoint; vertical dotted lines mark withdrawal at '+', '.join(f"{j['duration']:.6f}" for j,d in traces)+' s, respectively. The old pole is the first native finite-volume cell, x=L/(2×201), followed throughout reversal. A–C show actual signed rate contributions: basal activation, MCA-dependent activation, linear deactivation, diffusion, and their net sum, which equals the measured Rac increment divided by the native step duration. The two activation terms each include the multiplicative factor 1/(1+Rho_old²). Therefore variation in the basal activation curve directly reports changing Rho suppression, while variation in the MCA activation curve reflects both the bound-MCA switch and Rho suppression. Deactivation is -0.04 times newly solved Rac; diffusion is 0.02 times the discrete Neumann Laplacian of newly solved Rac. D–F subtract each term’s own value in the first recorded 0.5-s cue step; they expose the shifts hidden by cancellation between large positive activation and negative deactivation. These lower panels are changes in rates, not cumulative Rac changes; the sum of component changes equals the change in net rate. The first old-pole update is locally at its native resting balance to numerical precision. The cue is applied at the opposite end and is exactly zero at the old-pole cell throughout, so no cue contribution is plotted. Rac has no advective term. Endogenous mechanical effects enter these local terms through Rho, bound MCA and evolving neighboring Rac. This instantaneous additive budget does not measure independently evolved causal pathway removals. The native update is staggered: activation uses previous-step Rho and already updated MCA; decay and diffusion use the newly solved Rac. Observational wrappers record every actual Rac solve, including the fractional terminal cue step, without altering the dynamics. All eleven native state fields replay the retained original histories exactly at all common sample times. Independent audits recompute activation, decay, component sums, measured increments, integrated closure, the first-cell diffusion from original retained spatial profiles, canonical cue selection and zero local cue input. Supporting data also export the Rho suppression factor, bound-MCA switch and local Rac. Display and replay extend through 1200 s from cue onset. Native protrusion remains enabled.')
    hashes={};save(fig,STEM,caption,hashes)
    fig,axs=plt.subplots(1,3,figsize=(13.8,4.3),layout='constrained',sharex=True,sharey=True)
    for ax,(job,d),letter in zip(axs,traces,'ABC'):
        ax.plot(d['time_s'],d['rho_suppression_factor'],color='#0072B2',lw=2,label='Rho suppression factor')
        ax.plot(d['time_s'],d['mca_switch'],color='#D55E00',lw=2,label='MCA activation switch')
        ax.axvline(job['duration'],color='#777777',ls=':',lw=1.2)
        ax.set(title=f"{letter}   L = {job['length_um']} µm",xlim=(0,1200),ylim=(0,1.04),xlabel='Time from reversal-cue onset (s)',xticks=[0,300,600,900,1200]);ax.grid(alpha=.12)
    axs[0].set_ylabel('Dimensionless activation factor')
    fig.suptitle('Factors controlling old-pole Rac activation · relative-20% reversal',fontsize=14,fontweight='bold')
    fig.legend(handles=axs[0].get_legend_handles_labels()[0],labels=axs[0].get_legend_handles_labels()[1],loc='outside lower center',ncol=2,frameon=False)
    factor_caption='Companion to the exact old-pole Rac rate budget. Blue: 1/(1+Rho_old²), the fraction of all Rac activation remaining after Rho inhibition; a smaller factor means stronger suppression. Orange: 0.5×[1−tanh((bound_MCA_new−0.068)/0.01)], the MCA activation switch; a smaller switch means less MCA-dependent activation. Basal activation is 0.04 times the blue factor; MCA-dependent activation is 0.008 times the product of the blue and orange factors. These are the actual staggered inputs used by the Rac solve, at the same old-pole cell, lengths and canonical pulses as the contribution figure. Vertical dotted lines mark withdrawal. The local cue is identically zero. Curves describe coevolving inputs and are not independent pathway ablations.'
    save(fig,STEM+'_factors',factor_caption,hashes)
    for name,rows in [(STEM+'_summary.csv',summaries),(STEM+'_samples.csv',samples)]:
        p=output(name);csvwrite(p,rows);hashes[str(p.relative_to(OUT))]=sha(p)
    p=output(STEM+'_traces.npz');np.savez_compressed(p,**arrays);hashes[str(p.relative_to(OUT))]=sha(p)
    write(STEM+'_audit.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,checks=checks,
        native_core_hashes_match_canonical_campaign=True,local_cue_exactly_zero=True,native_replays=3))
    p=output(STEM+'_audit.json');hashes[str(p.relative_to(OUT))]=sha(p)
    write(STEM+'.manifest.json',dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',
        legend=caption,sha256=hashes,source_sha256=SOURCES,independent_audit_passed=True,
        native_replays=3,model='native corrected independent-endpoint deBelly',selection='shortest tested scalar1%-strength-passing relative20 cue',
        reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))]))
    print(output(STEM+'.png'));print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
