"""Render and independently audit one deBelly threshold figure."""
from pathlib import Path
import json,csv,textwrap,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from debelly_zero_crossing import OUT,LENGTHS,n,r,pol
FIG=n.ROOT/'figures/debelly_zero_crossing'
COLORS=['#0072B2','#E69F00','#009E73']
def main():
 FIG.mkdir(parents=True,exist_ok=True)
 jobs=r.read(OUT/'jobs.json');assert all((OUT/f'done_{i}.json').exists() for i,j in enumerate(jobs) if j['kind']!='linear')
 linear=r.read(OUT/'linear.json');ll=np.array([a['length'] for a in linear]);rates=np.array([a['rate'] for a in linear]);i=np.flatnonzero(rates[:-1]*rates[1:]<0)[0]
 bracket=[float(ll[i]),float(ll[i+1])];Lc=float(np.interp(0,rates[i:i+2],ll[i:i+2]));rows=[];seedcurves={};audits=[]
 for L in LENGTHS:
  for eps in [1e-5,1e-6]:
   path=OUT/f'seed_L{L:g}_e{eps:g}.npz'
   with np.load(path) as z:
    t=z['time'];x=z['x'];b=z['field_mca_bound'];rac=z['field_rac'];q=np.cos(np.pi*x/L);amp=b@q/(q@q)
    np.testing.assert_allclose(amp[0],eps,rtol=1e-8)
    np.testing.assert_allclose(b[0].sum(),201*.1*5/6.4,atol=1e-12)
    pm=pol(b,x,L);pr=pol(rac,x,L);mask=(t>=1000)&(t<=4000);rate=float(np.polyfit(t[mask],np.log(abs(amp[mask])),1)[0]);pred=float(np.interp(L,ll,rates))
    assert abs(rate-pred)<2e-6,(L,rate,pred)
    assert np.all(np.isfinite(b)) and np.min(b)>=0 and np.min(rac)>=0
    rows.append(dict(length=L,seed_amplitude=eps,measured_growth_rate=rate,tangent_rate=pred,final_abs_mca_polarity=abs(pm[-1]),final_abs_rac_polarity=abs(pr[-1]),horizon=40000,file_sha256=n.sha(path)))
    seedcurves[L,eps]=(t,np.abs(pr),np.abs(pm))
 n.csvwrite(OUT/'spontaneous_polarizability.csv',rows)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,axes=plt.subplots(3,3,figsize=(14.5,11.6),gridspec_kw={'height_ratios':[1,1,.98]})
 summaries=[];exports=[]
 for col,L in enumerate([18.,23.,39.]):
  ax=axes[0,col];ar=axes[1,col]
  for color,delta in zip(COLORS,[-100.,-50.,0.]):
   path=OUT/f'withdrawal_L{L:g}_offset{delta:g}.npz'
   with np.load(path) as z:
    D=float(z['duration']);x=z['x'];tp=z['pulse_time']-D;tr=z['recovery_time'];t=np.r_[tp,tr[1:]]
    b=np.r_[z['pulse_field_mca_bound'],z['recovery_field_mca_bound'][1:]];rac=np.r_[z['pulse_field_rac'],z['recovery_field_rac'][1:]];rho=np.r_[z['pulse_field_rho'],z['recovery_field_rho'][1:]]
    pm=pol(b,x,L);p=pol(rac,x,L);peak=rac.max(1);hpeak=rho.max(1);rec=tr>0
    post=pol(z['recovery_field_mca_bound'],x,L);crosses=np.flatnonzero(post[:-1]*post[1:]<0);cross=None
    if len(crosses):
     j=crosses[0];cross=float(np.interp(0,post[j:j+2][::-1],tr[j:j+2][::-1]))
    label=f'{abs(delta):g} s before zero' if delta else 'At zero'
    label+=f' (D={D:g} s)'
    ax.plot(t,pm,color=color,lw=1.9,label=label)
    if cross is not None:ax.scatter(cross,0,color=color,s=28,zorder=5,edgecolors='white',linewidths=.6)
    recovery_peak=z['recovery_field_rac'].max(1);threshold=recovery_peak[-1]+.1*(recovery_peak[0]-recovery_peak[-1]);relax=float(tr[np.flatnonzero(recovery_peak<=threshold)[0]])
    if delta==0:
     ar.axvspan(0,relax,color='#D55E00',alpha=.12);ar.text(.46,.94,f'90% excess decay\n{relax:g} s',transform=ar.transAxes,fontsize=8,va='top')
    ax.plot(10000,pm[-1],marker='o',color=color,ms=4)
    # Peak concentration relative to each branch's late basal reference; Rac solid, Rho dotted.
    ar.plot(t,peak/peak[-1],color=color,lw=1.7);ar.plot(t,hpeak/hpeak[-1],color=color,lw=1.1,ls=':')
    summaries.append(dict(length=L,offset=delta,duration=D,mca_at_withdrawal=float(post[0]),first_postwithdrawal_zero=cross,rac_90pct_excess_decay=relax,final_mca=float(pm[-1]),final_rac=float(p[-1]),rac_peak_at_withdrawal=float(peak[len(tp)-1]),rac_late_peak=float(peak[-1]),archive_sha256=n.sha(path)))
    for j in range(len(t)):exports.append(dict(length=L,offset=delta,cue_duration=D,time_from_withdrawal=float(t[j]),mca_polarity=float(pm[j]),rac_polarity=float(p[j]),rac_peak=float(peak[j]),rho_peak=float(hpeak[j])))
  for a in [ax,ar]:
   a.axvspan(-500,0,color='#9EA7AB',alpha=.12);a.axvline(0,color='#333333',ls='--',lw=1);a.set_xscale('symlog',linthresh=100);a.set_xlim(-600,10000);a.set_xticks([-300,0,100,1000,10000]);a.set_xticklabels(['−300','0','100','1,000','10,000']);a.grid(alpha=.15);a.set_xlabel('Time from cue withdrawal (s)')
  ax.axhline(0,color='#333333',lw=.8);ax.set_ylim(-.205,.205);ax.set_title(f'{chr(65+col)}   L = {L:g} µm',loc='left',weight='bold');ax.legend(frameon=False,fontsize=8,loc='lower left');ax.text(.98,.95,'cue off →',transform=ax.transAxes,ha='right',va='top',fontsize=9)
  ar.axhline(1,color='#777777',lw=.8);ar.set_ylim(0,13);ar.set_title('Rac relaxation after withdrawal',fontsize=11,loc='left');ar.text(.98,.94,'Rac: solid\nRho: dotted',transform=ar.transAxes,ha='right',va='top',fontsize=9)
 axes[0,0].text(.04,.88,'Zero crossed; polarity fades',transform=axes[0,0].transAxes,fontsize=9,color='#43525B')
 axes[0,1].text(.04,.88,'Reversed bias retained',transform=axes[0,1].transAxes,fontsize=9,color='#43525B')
 axes[0,2].text(.04,.88,'Reversed bias retained',transform=axes[0,2].transAxes,fontsize=9,color='#43525B')
 axes[0,0].set_ylabel('Bound-MCA polarity\npositive = old attachment bias')
 axes[1,0].set_ylabel('Peak concentration / late basal peak')
 a=axes[2,0];a.plot(ll,rates*1000,color='#65418A',lw=2);a.axhline(0,color='#555555',lw=.8);a.axvline(Lc,color='#555555',ls='--',lw=1)
 a.scatter([q['length'] for q in rows if q['seed_amplitude']==1e-5],[q['measured_growth_rate']*1000 for q in rows if q['seed_amplitude']==1e-5],s=22,color='#65418A',zorder=3)
 a.set_title('D   Growth of small spatial perturbations',loc='left',weight='bold');a.set_xlabel('Cell length (µm)');a.set_ylabel('Growth rate (10⁻³ s⁻¹)');a.text(.97,.06,f'Critical length ≈ {Lc:.2f} µm\nDecay below; growth above',transform=a.transAxes,ha='right',fontsize=9)
 a=axes[2,1];a.axvspan(Lc,23,color='#E69F00',alpha=.14);a.text(.04,.95,'Slow-onset regime: ≈19.5–23 µm',transform=a.transAxes,va='top',fontsize=8)
 for eps,ls,m in [(1e-5,'-','o'),(1e-6,'--','s')]:
  rr=[q for q in rows if q['seed_amplitude']==eps];a.plot([q['length'] for q in rr],[q['final_abs_rac_polarity'] for q in rr],ls=ls,marker=m,color='#0072B2',ms=4,label=f'MCA seed {eps:g}')
 a.axvline(Lc,color='#555555',ls='--',lw=1);a.axhline(.05,color='#999999',ls=':',lw=1);a.set_title('E   Spontaneous polarity from basal rest',loc='left',weight='bold');a.set_xlabel('Cell length (µm)');a.set_ylabel('|Rac polarity| at 40,000 s');a.legend(frameon=False,fontsize=9);a.set_ylim(bottom=0)
 a=axes[2,2]
 for L,color in zip([18.,21.,23.,39.],['#999999','#E69F00','#009E73','#0072B2']):
  t,p,m=seedcurves[L,1e-5];a.plot(t[1:],p[1:],label=f'L={L:g}',color=color,lw=1.8)
 a.set_yscale('log');a.set_ylim(1e-9,.3);a.axhline(.05,color='#999999',ls=':',lw=1);a.set_title('F   Slow growth near the transition',loc='left',weight='bold');a.set_xlabel('Time after small perturbation (s)');a.set_ylabel('|Rac polarity|');a.legend(frameon=False,fontsize=9,ncol=2)
 for a in axes[2]:a.grid(alpha=.15)
 fig.suptitle('deBelly: MCA zero crossing tracks memory reorientation',fontsize=19,weight='bold',y=.975)
 fig.text(.07,.94,'Residual Rac/Rho activity shifts withdrawal timing; short cells resist spontaneous polarization.',fontsize=12,color='#43525B')
 fig.subplots_adjust(left=.075,right=.98,top=.90,bottom=.19,hspace=.62,wspace=.29)
 caption=('A–C: Exact native basal-polarized checkpoints receive the same rightward gradient (slope 0.24). Withdraw 100 or 50 s before, or at, the uninterrupted-cue MCA zero crossing; D is pulse duration. Gray: cue on; dots: post-withdrawal MCA zero crossing. Orange: Rac excess decays by 90% (at-zero trial). Time after withdrawal uses a compressed scale beyond ±100 s. MCA polarity is its normalized first spatial moment; negative marks the new attachment bias. Residual Rac activity can carry MCA across zero after cue removal. Zero crossing is a memory-coordinate landmark, not a guarantee of retained reversal. D–F: Relaxed homogeneous basal state, no cue, zero-mean cosine bound-MCA perturbation (0.013% or 0.0013% of basal bound MCA). D compares native-step linear growth (line) with measured seed growth (dots); E shows finite-horizon polarity, not a probability. In the sampled ≈19.5–23 µm transition regime, spontaneous growth is slow and seed-dependent at 40,000 s. Short cells can retain a cued polarized branch while small perturbations decay. Native parameter regime: 201 cells, 0.5 s step; critical length is a discretized-model estimate.')
 fig.text(.075,.135,'\n'.join(textwrap.wrap(caption,170)),fontsize=9,va='top',linespacing=1.4)
 for ext in ['png','pdf','svg']:fig.savefig(FIG/f'debelly_zero_crossing.{ext}',dpi=220)
 plt.close(fig)
 n.csvwrite(OUT/'withdrawal_trajectories.csv',exports);r.write(OUT/'withdrawal_summary.json',summaries)
 audit=dict(generated_at=n.now(),producer='finished',all_jobs_complete=True,positive_finite_seed_fields=True,zero_mean_seed_and_amplitude_verified=True,measured_growth_matches_tangent_within_2e6=True,critical_length_bracket=bracket,interpolated_critical_length=Lc,script_sha256=n.sha(Path(__file__)),native_solver_sha256=n.sha(Path(n.ds.__file__)),figure_sha256={ext:n.sha(FIG/f'debelly_zero_crossing.{ext}') for ext in ['png','pdf','svg']})
 r.write(OUT/'figure_audit.json',audit)
 (FIG/'caption.txt').write_text(caption+'\nGenerated '+audit['generated_at']+'; simulation producer finished.\n')
 print(json.dumps(audit,indent=2));print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
