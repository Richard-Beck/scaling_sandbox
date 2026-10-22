"""Independent checks and standalone figure; does not edit canonical reports."""
from pathlib import Path
from dataclasses import replace
import json, sys
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=next(p for p in Path(__file__).resolve().parents if p.name=='manuscript')
OUT=ROOT/'data/runs/debelly_short_length_diagnostic_20261010_v2'
sys.path.insert(0,str(OUT/'frozen_sources/code/mechanisms'))
import debelly_short_length_diagnostic_20261010 as s

def stationary_balance(n,L,p):
    """Continuum stationary first cosine mode, independently of native stepping.

    Eliminate quasistatic membrane displacement and cortex velocity, then
    Rac/Rho and unbound MCA at zero growth. Compare cortical convergence
    replenishing bound-MCA asymmetry with diffusion/exchange removing it.
    These are linear loop coefficients per unit bound-MCA perturbation,
    not full nonlinear eigenvalues or fractions of a finite-cue outcome.
    """
    R,H=n.ds._homogeneous_rac_rho(p,.5)
    B=p.mca_density*p.mca_on/(p.mca_on+p.mca_off)
    alpha,beta=n.ds._derived_rates(p,np.array(B),np.array(0.),np.array(0.))
    fH=-p.rac_deactivation*alpha*2*H/(1+H*H)**2
    fB=-p.rac_deactivation*p.rac_mechanical_feedback/(2*p.mca_switch_scale)/np.cosh((B-p.mca_rac_threshold)/p.mca_switch_scale)**2/(1+H*H)
    gR=-p.rho_deactivation*beta*2*R/(1+R*R)**2
    gT=p.rho_deactivation*p.rho_mechanical_feedback/(2*p.tension_switch_scale)/np.cosh(p.tension_rho_threshold/p.tension_switch_scale)**2/(1+R*R)
    q2=(np.pi/L)**2;friction=p.linker_friction*B
    den=p.cortical_viscosity*q2+friction
    tension_per_rho=friction*p.active_stress/den
    ar=p.rac_deactivation*(1+p.rac_diffusion*q2)
    ah=p.rho_deactivation*(1+p.rho_diffusion*q2)
    rho_per_bound=gR*fB/(ar*(ah+gT*tension_per_rho)-gR*fH)
    gain=B*p.active_stress*q2/den*rho_per_bound
    loss=p.mca_bound_diffusion*q2+p.mca_off*p.mca_unbound_diffusion*q2/(p.mca_on+p.mca_unbound_diffusion*q2)
    return float(gain),float(loss)

def main():
    n,e,m=s.configure();checks=[];rows=[];traces={}
    for path in sorted((OUT/'workers').glob('*.json')):
        r=json.loads(path.read_text());p=path.with_suffix('.npz')
        assert n.sha(p)==r['trace_sha256']
        state,tr=n.read_checkpoint(p,'debelly');L=r['job']['length'];v=r['job']['variant']
        active=tr.fields['rac'];weights=2*(np.arange(201)+.5)/201-1
        moment=active@weights/active.sum(axis=1);contrast=np.ptp(active,axis=1)
        np.testing.assert_allclose(moment,tr.orientation,atol=2e-15)
        np.testing.assert_allclose(contrast,tr.contrast,atol=2e-15)
        assert abs(moment[-1]-r['final_orientation'])<2e-15
        total=(tr.fields['mca_bound']+tr.fields['mca_unbound']).sum(axis=1)
        mass_error=float(np.max(abs(total/total[0]-1)))
        assert mass_error<1e-11
        assert all(np.isfinite(f).all() for f in tr.fields.values())
        assert min(tr.fields[k].min() for k in ['rac','rho','mca_bound','mca_unbound'])>=0
        # Independently recompute the two pre-holdout windows, holdout
        # excursion and final window using retained full native arrays.
        rule=m['steady_state'];scales=m['case_template']['fixed_scales']
        end=float(tr.time[-1]);before=end-rule['holdout']; metrics=[]
        fields={k:tr.fields[k] for k in ['rac','rho','mca_bound','mca_unbound','velocity','tension']}
        fields['displacement']=tr.fields['displacement']-tr.fields['displacement'].mean(axis=1)[:,None]
        fields['physical_length']=L+tr.fields['boundary_separation']
        for stop in [before-rule['chunk'],before,end]:
            sel=(tr.time>=stop-rule['chunk'])&(tr.time<=stop);t=tr.time[sel]
            vals=[]
            for k,f in fields.items():
                y=f[sel]
                vals.extend([float(np.ptp(y,axis=0).max()/scales[k]),
                    float(abs(y[-1]-y[-2]).max()/(t[-1]-t[-2])*rule['chunk']/scales[k])])
            metrics.append(max(vals))
        i=int(np.argmin(abs(tr.time-before)));hold=tr.time>=before
        excursion=max(float(abs(f[hold]-f[i]).max()/scales[k]) for k,f in fields.items())
        assert r['convergence']['verified'] and max(metrics)<=rule['tolerance']
        assert excursion<=rule['holdout_tolerance']
        row=dict(length_um=L,variant=v,abs_rac_moment=float(abs(moment[-1])),rac_contrast=float(contrast[-1]),
            verified=True,elapsed_s=end,maximum_window_scaled_drift=max(metrics),holdout_excursion=excursion,
            mca_mass_relative_error=mass_error,trace_sha256=n.sha(p))
        rows.append(row);traces[L,v]=(tr.time,abs(moment));checks.append(row)
    assert len(rows)==9
    variants=[('native',s.parameters(n,'native')),
              ('half_cortical_viscosity',s.parameters(n,'half_cortical_viscosity')),
              ('half_unbound_mca_diffusion',replace(s.parameters(n,'native'),mca_unbound_diffusion=.15))]
    native_roots=json.loads((OUT/'results/linear_additional.json').read_text())['roots']
    roots=[]
    for v,p in variants:
        root=brentq(lambda L:np.subtract(*stationary_balance(n,L,p)),10,25)
        numeric=next(r['linear_onset_um'] for r in native_roots if r['variant']==v)
        assert abs(root-numeric)<.001
        roots.append(dict(variant=v,continuum_onset_um=root,native_tangent_onset_um=numeric))
    n.csvwrite(OUT/'results/continuation_summary.csv',rows)
    balances=[dict(length_um=L,**dict(zip(['flow_feedback_gain','mca_smoothing_loss'],stationary_balance(n,L,variants[0][1])))) for L in [15.,17.,18.,19.5344,23.]]
    n.write(OUT/'results/stationary_balance.json',dict(roots=roots,rows=balances,definition=stationary_balance.__doc__))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axs=plt.subplots(1,3,figsize=(14,4.6));colors=['#444444','#0072B2','#D55E00']
    Ls=np.linspace(12,25,220)
    ax=axs[0]
    gains,losses=np.array([stationary_balance(n,L,variants[0][1]) for L in Ls]).T
    ax.plot(Ls,gains*1000,label='Flow / MCA / Rac feedback',color=colors[1])
    ax.plot(Ls,losses*1000,label='MCA diffusion / exchange smoothing',color=colors[2])
    ax.axvline(roots[0]['continuum_onset_um'],ls=':',color='#555555')
    ax.set(title='A  Native stationary linear balance',xlabel='Nominal length (µm)',ylabel='Loop coefficient (10⁻³ s⁻¹)');ax.legend(frameon=False,fontsize=8)
    ax=axs[1]
    for L in [16.,17.,17.5,17.75,18.]:
        t,p=traces[L,'native'];ax.plot(t,p,label=f'{L:g} µm')
    ax.set_yscale('log');ax.axhline(.05,color='#888888',ls=':',label='Report eligibility = 0.05')
    ax.set(title='B  Continue the polarized L = 18 state',xlabel='Basal continuation time (s)',ylabel='Absolute Rac first moment',ylim=(1e-7,.12));ax.legend(frameon=False,fontsize=8)
    ax=axs[2]
    selected=[r for r in rows if r['length_um']==17]
    names={'native':'Native','half_gtpase_diffusion':'½ Rac/Rho\ndiffusion','half_cortical_viscosity':'½ cortical\nviscosity','half_unbound_mca_diffusion':'½ unbound-MCA\ndiffusion'}
    order=list(names)
    selected=sorted(selected,key=lambda r:order.index(r['variant']))
    ax.bar(np.arange(4),[r['abs_rac_moment'] for r in selected],color=['#888888','#999999',colors[1],colors[2]])
    ax.set_xticks(np.arange(4),[names[r['variant']] for r in selected],fontsize=8)
    ax.axhline(.05,color='#888888',ls=':');ax.set(title='C  Rescue at L = 17 µm',ylabel='Verified basal Rac first moment',ylim=(0,.085))
    fig.suptitle('Short-cell polarity: attachment smoothing outruns mechanical reinforcement',fontsize=14)
    fig.tight_layout(rect=(0,0.04,1,.93))
    fig.text(.02,.015,'Corrected independent endpoints. Native basal input, no cue. All nine continuations pass two convergence windows and a 10,000-s holdout.',fontsize=9)
    for ext in ['png','pdf','svg']:fig.savefig(OUT/'results'/f'short_length_diagnostic.{ext}',dpi=180)
    plt.close(fig)
    audit=dict(generated_at=n.now(),producer='finished',passed=True,independent_continuation_checks=checks,
               independently_derived_continuum_onsets=roots,script_sha256=n.sha(Path(__file__)),
               figures={ext:n.sha(OUT/'results'/f'short_length_diagnostic.{ext}') for ext in ['png','pdf','svg']})
    n.write(OUT/'results/audit.json',audit)
    print(json.dumps(dict(rows=rows,roots=roots),indent=2))

if __name__=='__main__':main()
