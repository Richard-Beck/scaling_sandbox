"""Development-only Figure 9F native MCA budget at two cue-relative times.

Reuses corrected continuous-cue histories because both evaluation points precede
withdrawal in the exact selected finite pulses. Performs only the existing
half-second diagnostic budget steps; no new long trajectories.
"""
from pathlib import Path
import length_dev_output_layout as devlayout
from datetime import datetime,timezone
from types import SimpleNamespace
import csv,hashlib,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'code/core'),str(ROOT/'code/mechanisms')]
import debelly_corrected_reversal as corrected
import debelly_corrected_advection as ad
n=corrected.n
RUN=ROOT/'data/runs/length_quality_debelly_traces_20261009'
OLD=ROOT/'data/runs/length_debelly_review_20261009'
FLOW=ROOT/'data/runs/debelly_corrected_advection_20261008'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/native_reversal'
def output(name):return devlayout.dev_output_path(OUT,name)
def artifact_key(path):return devlayout.dev_artifact_key(OUT,path)
LENGTHS=[18,19,21,23,25,27,29,31,33,35,37,39]
SOURCES={str(Path(devlayout.__file__).relative_to(ROOT)):hashlib.sha256(Path(devlayout.__file__).read_bytes()).hexdigest()}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':9,'pdf.fonttype':42,'svg.fonttype':'none'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def csvwrite(p,rs):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 for p in [Path(__file__),Path(ad.__file__),Path(corrected.__file__),Path(n.__file__),Path(n.ds.__file__)]:SOURCES[str(p.relative_to(ROOT))]=sha(p)
 assert n.model_for('debelly').parameters.center_boundary_motion is False
 assert read(FLOW/'independent_audit.json')['passed']
 assert read(OLD/'independent_audit.json')['passed']
 old=read(OLD/'manifest.json');manifest=read(RUN/'manifest.json')
 jobs={int(j['case']['length']):j for j in manifest['jobs'] if j['case']['condition']=='fixed_slope'}
 assert sorted(jobs)==LENGTHS
 budgets=[];audits=[]
 for L in LENGTHS:
  j=jobs[L];c=j['case'];D=float(j['duration']);assert c['condition']=='fixed_slope' and np.isclose(c['gradient_slope'],.24,rtol=0,atol=1e-15)
  cp=ROOT/c['checkpoint'];assert sha(cp)==c['checkpoint_sha256'];SOURCES[str(cp.relative_to(ROOT))]=sha(cp)
  if L<=21:
   slot=next(i for i,q in enumerate(old['jobs']) if q['case']['length']==L and q['case']['condition']=='fixed_slope' and q['purpose']=='continuous_and_zero_withdrawal')
   p=OLD/f'workers/{slot:03d}.npz';meta=read(p.with_suffix('.json'));assert sha(p)==meta['archive_sha256'];prefix='pulse_'
  else:
   p=FLOW/f'L{L:g}_native.npz';meta=read(FLOW/f'L{L:g}_summary.json');assert sha(p)==meta['native_sha256'];prefix=''
  SOURCES[str(p.relative_to(ROOT))]=sha(p)
  with np.load(p) as z:
   t=z[prefix+'time'];x=z['x'];assert np.allclose(np.diff(t),.5,atol=1e-10)
   # Nearest retained native sample preserves the original exact mechanical
   # identity; interpolating mechanical fields would not do so.
   requests=[('A',60.),('B',D-60.)];indices=[]
   for panel,target in requests:
    if target<0 or target>D or target>t[-1]:
     budgets.append(dict(panel=panel,length=L,condition='fixed_slope',cue_duration_s=D,requested_time_from_onset_s=target,actual_time_from_onset_s=None,offset_s=None,time_before_withdrawal_s=None,status='unavailable',phase='outside_selected_cue',rate_boundary=None,rate_rho=None,rate_other=None,reorientation_rate=None));indices.append(None)
    else:
     i=int(np.argmin(abs(t-target)));assert i>0;assert abs(t[i]-target)<=.25+1e-9;indices.append(i)
   needed=sorted({i for i in indices if i is not None}|{i-1 for i in indices if i is not None});lookup={v:i for i,v in enumerate(needed)}
   fields={k:z[prefix+'field_'+k][needed].copy() for k in n.DB_NAMES}
   with np.load(cp) as checkpoint:
    for k in n.DB_NAMES:np.testing.assert_allclose(z[prefix+'field_'+k][0],checkpoint['state_'+k],atol=1e-12)
  for (label,target),i in zip(requests,indices):
   if i is None:continue
   pair=[lookup[i-1],lookup[i]];tr=SimpleNamespace(x=x,time=t[[i-1,i]],fields={k:f[pair].copy() for k,f in fields.items()})
   r,_=ad.decompose(tr,1,L,corrected.response_cue(c))
   others=sum(r['rate_'+k] for k in ['membrane','diffusion','exchange','replenishment'])
   r.update(panel=label,condition='fixed_slope',cue_duration_s=D,requested_time_from_onset_s=target,actual_time_from_onset_s=float(t[i]),offset_s=float(t[i]-target),time_before_withdrawal_s=float(D-t[i]),status='available',phase='cue_on',rate_other=others,source_archive=str(p.relative_to(ROOT)))
   assert abs(r['rate_boundary']+r['rate_rho']+others-r['reorientation_rate'])<1e-11
   budgets.append(r);audits.append(dict(panel=label,length=L,exact_corrected_checkpoint=True,source_archive_sha256=sha(p),mechanical_residual=r['mechanical_residual'],source_budget_residual=r['budget_residual'],timing_offset_s=r['offset_s']))
 # All 12 selected gradient cues exceed60 s; no requested point falls in recovery.
 assert len(budgets)==24 and all(r['status']=='available' and r['phase']=='cue_on' for r in budgets)
 keys=list(dict.fromkeys(k for r in budgets for k in r));csvwrite(output('Fig9F_cue_timing.csv'),[{k:r.get(k) for k in keys} for r in budgets])
 fig,axs=plt.subplots(1,2,figsize=(10.7,4.6),layout='constrained',sharex=True,sharey=True)
 curves=[('rate_boundary','#D55E00','-','Endpoint advection'),('rate_rho','#0072B2','-','Rho advection'),('rate_other','#999999','-','Other terms'),('reorientation_rate','#222222','--','Net rate')]
 for ax,label,title in zip(axs,'AB',['60 s after cue initiation','60 s before cue withdrawal']):
  ax.text(-.13,1.04,label,transform=ax.transAxes,fontsize=12,fontweight='bold',fontfamily='DejaVu Sans');ax.set_title(title,pad=10);ax.grid(alpha=.13);ax.axhline(0,color='#777777',lw=.8)
  rs=sorted([r for r in budgets if r['panel']==label],key=lambda r:r['length'])
  for key,color,ls,name in curves:ax.plot([r['length'] for r in rs],np.array([r[key] for r in rs])*1e4,color=color,ls=ls,lw=1.8,label=name)
  ax.set_xlim(17,40);ax.set_xticks([18,23,27,31,35,39]);ax.set_xlabel('Cell length (µm)')
 axs[0].set_ylabel('−dP(MCA)/dt (10⁻⁴ s⁻¹)')
 handles=[Line2D([],[],color=color,ls=ls,lw=1.8,label=name) for key,color,ls,name in curves]
 fig.legend(handles=handles,loc='outside lower center',ncol=4,frameon=False)
 hashes={}
 for ext in ['png','pdf','svg']:
  p=output('Fig9F_cue_timing.'+ext);fig.savefig(p,dpi=220,facecolor='white');hashes[artifact_key(p)]=sha(p)
 plt.close(fig)
 maximum_offset=max(abs(r['offset_s']) for r in budgets)
 caption='Development comparison of the Figure 9F bound-MCA polarity-rate budget at two cue-relative times. A: 60 s after initiation. B: 60 s before removal of each length’s selected shortest orientation-reversing gradient pulse. All 12 eligible lengths (18–39 µm) retain the exact corrected basal starts, constant slope 0.24 and selected cue durations used for Figure 8/9. Endpoint-driven advection (orange), Rho-driven advection (blue), and other terms (gray: membrane motion, diffusion, exchange and uniform restoration) sum to the net rate (black dashed). Positive rates erase the original attachment bias. The rate includes the changing MCA normalization and is the same half-second native-step identity as Figure 9F, evaluated at the requested states rather than at MCA zero. Panel A uses the exact recorded time 60 s. Panel B uses the nearest retained 0.5-s native sample to duration minus 60 s; maximum absolute timing offset is '+f'{maximum_offset:.6f}'+' s, and every requested/actual time is exported. Both panels share the same axes. Every selected cue exceeds 60 s, so all 24 states are cue-on and none is unavailable or in recovery. This is instantaneous source bookkeeping at saved states, not independently evolved pathway interventions.'
 (output('Fig9F_cue_timing_caption.txt')).write_text(caption+'\n')
 methods=dict(state_selection='Corrected continuous-cue native histories are identical to the selected finite pulse before withdrawal; both requested times precede withdrawal at every length. No clamping, extrapolation or state interpolation.',native_budget='Unchanged debelly_corrected_advection.decompose: source superposition at fixed state with full-flow upwind donors frozen; native 0.5-s diagnostic step and exact changing normalization. Mechanical velocity superposition and all terms independently checked for residual<1e-11.',selection_gate='Same shortest tested ORIENTATION-reversing pulses used in Figure 8 and its Figure 9 duration comparison. These are not newly selected all-field or full-strength pulses.',timing='A=60 s exactly; B=nearest0.5-s sample to selected duration−60 s. Actual times and offsets recorded per length.',new_compute='24 brief native half-second budget steps; no new long integrations or parameter/cue changes.',unavailable='Out-of-cue or unrecorded requested times are explicitly unavailable; no such cases occur in the current selected constant-slope series.')
 hashes[artifact_key(output('Fig9F_cue_timing.csv'))]=sha(output('Fig9F_cue_timing.csv'))
 summary={}
 for label in 'AB':
  rs=[r for r in budgets if r['panel']==label];summary[label]={key:dict(minimum=min(r[key] for r in rs),maximum=max(r[key] for r in rs),at_shortest=rs[0][key],at_longest=rs[-1][key]) for key in ['rate_boundary','rate_rho','rate_other','reorientation_rate']}
 record=dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',description=caption,legend=caption,methods=methods,sha256=hashes,source_sha256=SOURCES,renderer_sha256=sha(Path(__file__)),reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))],audit=audits,maximum_timing_offset_s=maximum_offset,summary=summary,current_report_changed=False)
 (output('Fig9F_cue_timing.manifest.json')).write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(summary,indent=2));print(output('Fig9F_cue_timing.png'))
if __name__=='__main__':main()
