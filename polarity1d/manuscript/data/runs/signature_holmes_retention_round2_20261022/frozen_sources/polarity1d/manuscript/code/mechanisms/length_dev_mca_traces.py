"""Development-only long-time views of the exact Figure 8 native pulse histories.

No integrations; writes only figures/length_mechanisms_20261008/dev/.
"""
from pathlib import Path
import length_dev_output_layout as devlayout
from datetime import datetime,timezone
import csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'data/runs/length_quality_debelly_traces_20261009'
OLD=ROOT/'data/runs/length_debelly_review_20261009'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/native_reversal'
def output(name):return devlayout.dev_output_path(OUT,name)
def artifact_key(path):return devlayout.dev_artifact_key(OUT,path)
LENGTHS=[18,19,21,23,25,27,29,31,33,35,37,39]
GEOS=['fixed_physical','fixed_relative','fixed_slope']
TITLES=['Fixed width: 3 µm','Relative width: 20%','Constant slope: 0.24']
COLORS={L:plt.get_cmap('turbo')(.06+.88*i/11) for i,L in enumerate(LENGTHS)}
SOURCES={str(Path(devlayout.__file__).relative_to(ROOT)):hashlib.sha256(Path(devlayout.__file__).read_bytes()).hexdigest()}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def load(p,keys):
 SOURCES[str(p.relative_to(ROOT))]=sha(p)
 with np.load(p) as z:return {k:z[k].copy() for k in keys}
def polarity(y,x,L):return y@(2*x/L-1)/y.sum(axis=-1)
def crossing(t,y):
 ix=np.flatnonzero(y>=0)
 if not len(ix):return None
 i=ix[0]
 return float(t[i]) if i==0 else float(np.interp(0.,y[i-1:i+1],t[i-1:i+1]))
def persistent(t,valid):
 good=np.logical_and.accumulate(np.asarray(valid)[::-1])[::-1];ix=np.flatnonzero(good)
 return float(t[ix[0]]) if len(ix) else None

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 SOURCES[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 helper=ROOT/'code/mechanisms/length_quality_debelly_traces_20261009.py';SOURCES[str(helper.relative_to(ROOT))]=sha(helper)
 assert read(ROOT/'data/runs/length_mechanisms_20261008/quality_20261009_debelly/independent_audit.json')['passed']
 manifest=read(RUN/'manifest.json');assert manifest['center_boundary_motion'] is False
 traces={};rows=[]
 for slot,j in enumerate(manifest['jobs']):
  c=j['case'];L=int(c['length']);geo=c['condition'];D=j['duration']
  p=RUN/f'workers/{slot:03d}.npz' if j['existing_pulse_slot'] is None else OLD/f"workers/{j['existing_pulse_slot']:03d}.npz"
  meta=read(p.with_suffix('.json'));assert sha(p)==meta['archive_sha256']
  z=load(p,['x']+[phase+'_'+key for phase in ['pulse','recovery'] for key in ['time','field_rac','field_mca_bound']])
  cp=ROOT/c['checkpoint'];assert sha(cp)==c['checkpoint_sha256'];checkpoint=load(cp,['state_rac','state_mca_bound'])
  for field in ['rac','mca_bound']:
   np.testing.assert_allclose(z['pulse_field_'+field][0],checkpoint['state_'+field],atol=1e-12)
   np.testing.assert_allclose(z['recovery_field_'+field][0],z['pulse_field_'+field][-1],atol=1e-12)
  x=z['x'];tp=z['pulse_time'];tr=z['recovery_time'];op=polarity(z['pulse_field_rac'],x,L);o=polarity(z['recovery_field_rac'],x,L);q=np.ptp(z['recovery_field_rac'],axis=1)
  assert abs(tp[-1]-D)<1e-9 and tr[-1]==10000.
  tail=tr>=9500.;gated=bool(o[tail].min()>.05 and q[tail].min()>=.05 and np.ptp(o[tail])<.005);assert gated==meta['metrics']['orientation_success'] and gated
  target=abs(c['initial_orientation']);strength=bool(gated and np.max(abs(o[tail]/target-1))<=.01);assert strength==meta['metrics']['full_polarity_success']
  threshold=persistent(tr,(o>.05)&(q>=.05))
  suffixmin=np.minimum.accumulate(o[::-1])[::-1];suffixmax=np.maximum.accumulate(o[::-1])[::-1];suffixcontrast=np.minimum.accumulate(q[::-1])[::-1]
  stable=persistent(tr,(suffixmin>.05)&(suffixcontrast>=.05)&((suffixmax-suffixmin)<.005))
  within=persistent(tr,(o>.05)&(q>=.05)&(abs(o/target-1)<=.01))
  reflected=checkpoint['state_rac'][::-1]
  profile_error=np.max(abs(z['recovery_field_rac']-reflected)/np.maximum(abs(reflected),1e-12),axis=1)
  profile_time=persistent(tr,profile_error<=.01)
  pm=-polarity(z['pulse_field_mca_bound'],x,L);rm=-polarity(z['recovery_field_mca_bound'],x,L)
  t=np.r_[tp,D+tr[1:]];m=np.r_[pm,rm[1:]]
  rows.append(dict(condition=geo,length=L,cue_duration_s=D,recovery_horizon_s=10000.,initial_rac_polarity=c['initial_orientation'],first_rac_zero_from_onset_s=crossing(tp,op),first_mca_zero_from_onset_s=crossing(t,m),persistent_mca_new_rear_sign_from_onset_s=persistent(t,m>=0),persistent_rac_orientation_and_contrast_from_withdrawal_s=threshold,persistent_rac_orientation_and_contrast_from_onset_s=D+threshold if threshold is not None else None,persistent_rac_orientation_contrast_range_from_withdrawal_s=stable,persistent_rac_orientation_contrast_range_from_onset_s=D+stable if stable is not None else None,persistent_reflected_strength_1pct_from_withdrawal_s=within,persistent_reflected_strength_1pct_from_onset_s=D+within if within is not None else None,persistent_reflected_rac_profile_1pct_from_withdrawal_s=profile_time,persistent_reflected_rac_profile_1pct_from_onset_s=D+profile_time if profile_time is not None else None,final_reflected_rac_profile_max_relative_error=float(profile_error[-1]),reflected_rac_profile_final500_1pct_pass=bool(np.max(profile_error[tail])<=.01),orientation_final500_pass=gated,reflected_strength_final500_pass=strength,independent_full_field_tail_pass=meta['metrics']['field_tail_pass'],final_rac_polarity=float(o[-1]),source_archive=str(p.relative_to(ROOT))))
  traces[geo,L]=(z,op,o,pm,rm,meta['metrics']['field_tail_pass'])
 table=output('MCA_traces_timings.csv')
 with table.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 fig,axs=plt.subplots(2,3,figsize=(12.5,8.2),layout='constrained',sharex=True,sharey='row')
 for col,geo in enumerate(GEOS):
  for row in [0,1]:
   ax=axs[row,col];ax.text(-.13,1.04,chr(65+3*row+col),transform=ax.transAxes,fontsize=12,fontweight='bold',fontfamily='DejaVu Sans');ax.grid(alpha=.13);ax.axhline(0,color='#8c8c8c',lw=.7)
   if row==0:ax.set_title(TITLES[col],pad=10);ax.axhspan(.99,1.01,color='#dddddd',alpha=.9,zorder=0);ax.axhline(1.,color='#777777',lw=.8,ls=':')
   for L in LENGTHS:
    z,op,o,pm,rm,flag=traces[geo,L];D=float(z['pulse_time'][-1]);pulse,rec=(op/abs(op[0]),o/abs(op[0])) if row==0 else (pm,rm)
    ax.plot(z['pulse_time'],pulse,color=COLORS[L],lw=1.35);ax.plot(D+z['recovery_time'],rec,color=COLORS[L],lw=1.3,ls='--');ax.plot([D],[pulse[-1]],'o',ms=4,mec=COLORS[L],mfc=COLORS[L] if flag else 'white',mew=1,zorder=4)
   ax.set_xscale('symlog',linthresh=10,linscale=.7);ax.set_xlim(0,11200);ax.set_xticks([0,10,100,1000,10000],['0','10','100','1,000','10,000']);ax.set_ylim((-1.3,7.8) if row==0 else (-.17,.17))
   if col==0:ax.set_ylabel('Rac polarity / initial magnitude' if row==0 else '−MCA polarity toward new rear')
   if row==1:ax.set_xlabel('Time from cue onset (s; log above 10 s)')
 handles=[Line2D([],[],color=COLORS[L],lw=2,label=str(L)) for L in LENGTHS]+[Line2D([],[],color='#333333',lw=1.5,label='Cue'),Line2D([],[],color='#333333',ls='--',lw=1.5,label='Recovery'),Line2D([],[],color='#333333',marker='o',ls='',mfc='#333333',label='Stable fields'),Line2D([],[],color='#333333',marker='o',ls='',mfc='white',label='Changing fields')]
 fig.legend(handles=handles,loc='outside lower center',ncol=8,frameon=False,title='Cell length (µm) · circles mark withdrawal; fill gives final native-field tail flag')
 hashes={}
 for ext in ['png','pdf','svg']:
  p=output('MCA_traces.'+ext);fig.savefig(p,dpi=220,facecolor='white');hashes[artifact_key(p)]=sha(p)
 plt.close(fig)
 caption='Development view of the exact Figure 8 selected pulses, extended through all 10000 s of basal recovery. A–C show Rac polarity divided by its own pre-cue magnitude; D–F show sign-reversed bound-MCA polarity, positive for the new rear. Colors identify all 12 eligible lengths; columns compare the fixed 3 µm patch, 20% patch and slope 0.24 gradient. Solid lines are cue-on; dashed lines are recovery; circles mark each selected cue withdrawal. Filled circles pass the independent native-field tail flag; open circles retain changing fields. Time is measured from cue onset, linear from 0 to 10 s and logarithmic thereafter. The Rac target band spans 0.99–1.01, with its center at 1 marking the reflected initial scalar magnitude. The exact selected pulses satisfy the official final 500-s 1%-strength criterion in 15 of 36 cases (7 fixed-width, 6 relative-width, 2 gradient); 20 reach a persistent 1%strength band before the final recorded sample, but some have fewer than 500 s remaining. All 36 satisfy the weaker orientation gate. Early Rac sign change is followed by relaxation and much later scalar-strength regrowth in many cells. Timing measurements use recorded samples only; they are conditional on the observed 10000-s recovery and do not extrapolate beyond it.'
 (output('MCA_traces_caption.txt')).write_text(caption+'\n')
 methods=dict(persistent_orientation='First recorded recovery sample after which Rac polarity>0.05 and contrast>=0.05 at every remaining sample through 10000 s. Range criterion is reported separately.',persistent_orientation_with_range='First recorded recovery sample after which all remaining samples have orientation>0.05, contrast>=0.05 and their entire remaining orientation range<0.005. This is stricter than the original final 500-s gate and is a descriptive arrival readout.',persistent_reflected_strength='First recorded recovery sample after which all remaining samples have orientation>0.05, contrast>=0.05 and Rac polarity within 1% of absolute initial polarity. Empty CSV cells mean not reached within recorded recovery.',reflected_rac_profile='Pointwise relative difference from the exact pre-cue Rac profile reflected on nominal coordinates: max_x abs(Rac(x,t)-Rac_initial(L-x))/max(abs(Rac_initial(L-x)),1e-12)<=0.01 at every remaining recorded sample. This checks Rac shape only, not all mechanical/chemical state fields.',first_zero='First linearly interpolated sign crossing between recorded samples; MCA first crossing can precede a later recrossing. Persistent MCA sign uses recorded samples and is exported separately.',time_resolution='Pulse sampling 0.5 or 1 s as retained; recovery 2.5 s. Persistent arrivals are sampled bounds, not interpolated.',unchanged='Exact Figure8 histories, cue durations, colors and corrected independent-endpoint physics. No new integrations; current report files untouched.')
 stats={}
 for geo in GEOS:
  rr=[r for r in rows if r['condition']==geo];stats[geo]={}
  for key in ['first_rac_zero_from_onset_s','persistent_rac_orientation_and_contrast_from_withdrawal_s','persistent_rac_orientation_contrast_range_from_withdrawal_s','persistent_reflected_strength_1pct_from_withdrawal_s','persistent_reflected_rac_profile_1pct_from_withdrawal_s']:
   values=[r[key] for r in rr if r[key] is not None];stats[geo][key]=dict(reached=len(values),not_reached_lengths=[r['length'] for r in rr if r[key] is None],minimum=min(values) if values else None,maximum=max(values) if values else None,median=float(np.median(values)) if values else None)
 stats['official_final500_strength_passes']={geo:sum(r['reflected_strength_final500_pass'] for r in rows if r['condition']==geo) for geo in GEOS}
 (output('MCA_traces_summary.json')).write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),methods=methods,summary=stats),indent=2)+'\n')
 hashes[artifact_key(table)]=sha(table)
 record=dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',description=caption,legend=caption,methods=methods,sha256=hashes,renderer_sha256=sha(Path(__file__)),source_sha256=SOURCES,reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))],native_histories=36,new_integrations=0,current_report_changed=False)
 (output('MCA_traces.manifest.json')).write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(stats,indent=2));print(output('MCA_traces.png'))
if __name__=='__main__':main()
