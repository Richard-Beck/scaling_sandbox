"""Three audited, vector-exported landscape pages for deBelly mechanism report."""
from pathlib import Path
import csv,gzip,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
CODE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE));import gradient_reversal as g
r,n,e=g.r,g.n,g.e
OUT=n.DATA/'runs/fixed_slope_mechanisms_20261006/debelly'
FIG=n.ROOT/'figures/fixed_slope_mechanisms/debelly'
RUN=n.DATA/'runs/fixed_slope_reversal_20261006'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':10,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.6})
COL={'pass':'#007f78','fail':'#b54d28','clamp':'#785aa6'}

def load(L,label='pass'):return np.load(OUT/f'L{L}_{label}.npz')
def moment(z,phase):
 b=z[phase+'_field_mca_bound']+z[phase+'_field_mca_unbound'];q=2*z['x']/(z['x'][0]+z['x'][-1])-1
 return b@q/b.sum(1)

def save(fig,name):
 fig.savefig(FIG/(name+'.pdf'),facecolor='white');fig.savefig(FIG/(name+'.png'),dpi=220,facecolor='white');plt.close(fig)

def audit():
 with gzip.open(RUN/'results/records.json.gz','rt') as f:record=json.load(f)
 reports=[]
 for job in r.read(OUT/'jobs.json'):
  case=next(c for c in record['cases'] if c['case_index']==job['case']['case_index']);trial=next(t for t in case['tests'] if t['duration']==job['duration'])
  z=load(int(case['length']),job['label']);errors={}
  with np.load(RUN/trial['trace']) as old:
   prefix=trial['trace_prefix']
   for key in ['pulse_time','pulse_orientation','pulse_contrast','recovery_time','recovery_orientation','recovery_contrast']:
    assert z[key].shape==old[prefix+key].shape
    errors[key]=float(np.max(abs(z[key]-old[prefix+key])))
   # Full native selected tail information, independently regenerated.
   tail=z['recovery_time']>=9500
   for key in case['fixed_scales']:
    if key=='physical_length':y=case['length']+z['recovery_field_boundary_separation']
    else:y=z['recovery_field_'+key]
    if key=='displacement':y=y-y.mean(1)[:,None]
    a=y[tail]
    for suffix,value in [('min',a.min(0)),('max',a.max(0)),('penultimate',a[-2]),('last',a[-1])]:
     errors[key+'_'+suffix]=float(np.max(abs(value-old[prefix+'fieldtail_'+key+'_'+suffix])))
  summary=r.read(OUT/f'L{case["length"]:g}_{job["label"]}.json')
  assert summary['orientation_success']==trial['orientation_success'];assert summary['full_polarity_success']==trial['full_polarity_success'];assert summary['field_tail_pass']==trial['field_tail_pass']
  assert max(errors.values())<1e-10
  reports.append(dict(length=case['length'],label=job['label'],duration=job['duration'],original_trace=trial['trace'],trace_prefix=trial['trace_prefix'],maximum_absolute_array_difference=max(errors.values()),checkpoint_hash_pass=n.sha(n.ROOT/case['checkpoint'])==case['checkpoint_sha256'],gates_match=True))
 r.write(OUT/'independent_reconstruction_audit.json',dict(generated_at=n.now(),passed=True,trials=reports))

def snapshots():
 colors=['#777777','#3999bf','#d68c28','#007f78']
 fig,axes=plt.subplots(4,3,figsize=(11,7),layout='constrained')
 fields=[('rac','Rac (native concentration)'),('rho','Rho (native concentration)'),('mca_bound','Bound MCA'),('velocity','Cortical velocity (µm / time)')]
 for j,L in enumerate([18,23,39]):
  z=load(L);t=z['pulse_time'];d=float(z['pulse_duration']);x=z['x']
  for i,(field,ylabel) in enumerate(fields):
   ax=axes[i,j]
   for tt,col in zip([0,25,100,d],colors):
    k=int(np.argmin(abs(t-tt)));ax.plot(x,z['pulse_field_'+field][k],color=col,label=f't = {t[k]:.0f}'+(' (cutoff)' if tt==d else ''))
   if field=='mca_bound':ax.axhline(.068,color='#777777',linestyle=':',linewidth=.8)
   ax.set_xlim(0,L);ax.grid(alpha=.12)
   if j==0:ax.set_ylabel(ylabel)
   if i==0:
    ax.set_title(f'L = {L}; selected pulse = {d:.1f}\nRamp peak = {L*.24:.2f}')
    ax.legend(frameon=False,ncol=2,loc='upper left',fontsize=6.8)
   if i==3:ax.set_xlabel('Nominal position (µm)')
 fig.suptitle('deBelly: the cue turns Rac quickly, while adhesion memory changes slowly',fontweight='bold',fontsize=14)
 save(fig,'debelly_01_on_cue_snapshots')

def recovery():
 colors=['#5a9dc1','#9dc8da','#db933e','#b45344','#222222']
 fig,axes=plt.subplots(4,3,figsize=(11,7),layout='constrained')
 fields=[('rac','Rac'),('mca_bound','Bound MCA'),('velocity','Cortical velocity (µm / time)'),('tension','Membrane tension (native force)')]
 for j,L in enumerate([18,23,39]):
  z=load(L);cont=np.load(OUT/f'L{L}_pass_continuation.npz');x=z['x'];t=z['recovery_time']
  for i,(field,ylabel) in enumerate(fields):
   ax=axes[i,j]
   for tt,col in zip([250,1000,5000,10000,30000],colors):
    if tt<=10000:y=z['recovery_field_'+field][int(np.argmin(abs(t-tt)))]
    else:y=cont['field_'+field][-1]
    ax.plot(x,y,color=col,label=f'+{tt:,}')
   if field=='mca_bound':ax.axhline(.068,color='#777777',linestyle=':',linewidth=.8)
   ax.set_xlim(0,L);ax.grid(alpha=.12)
   if j==0:ax.set_ylabel(ylabel)
   if i==0:
    ax.set_title(f'L = {L}\n10,000-unit field-tail '+('passes' if L==18 else 'fails'))
    ax.legend(frameon=False,ncol=3,loc='upper left',fontsize=6.5)
   if i==3:ax.set_xlabel('Nominal position (µm)')
 fig.suptitle('deBelly: persistent right Rac can precede settled mechanics and MCA',fontweight='bold',fontsize=14)
 save(fig,'debelly_02_basal_recovery_snapshots')

def diagnostics():
 fig,axes=plt.subplots(2,3,figsize=(11,7),layout='constrained')
 rows=[x for x in csv.DictReader(open(RUN/'results/endpoints.csv')) if x['model']=='debelly' and x['condition']=='relative_20pct' and x['status']=='success']
 ax=axes[0,0]
 for definition,label,col in [('orientation','Ordinary orientation','#007f78'),('full_polarity','Strength replacement','#777777')]:
  a=[x for x in rows if x['definition']==definition];xx=[float(x['length']) for x in a];yy=[float(x['endpoint']) for x in a];ax.plot(xx,yy,color=col,label=label)
  for x,y,row in zip(xx,yy,a):ax.plot(x,y,'o',markerfacecolor=col if row['selected_field_tail_pass']=='True' else 'white',markeredgecolor=col,markersize=4)
 for L in [18,23,39]:
  y=float(next(x['endpoint'] for x in rows if x['definition']=='orientation' and float(x['length'])==L));ax.annotate(str(L),(L,y),xytext=(4,5),textcoords='offset points',fontweight='bold')
 ax.set_title('A  The measured duration has a minimum');ax.set_xlabel('Nominal length (µm)');ax.set_ylabel('Shortest passing sampled pulse');ax.legend(frameon=False)
 ax.text(.04,.42,'Open markers: native field tail fails\nDuration excludes basal recovery',transform=ax.transAxes,fontsize=7,color='#555555')
 ax=axes[0,1]
 for L,col in [(23,'#007f78'),(39,'#785aa6')]:
  z=load(L);ax.plot(z['pulse_time'],moment(z,'pulse'),color=col,label=f'L = {L}')
  ax.plot(float(z['pulse_duration']),moment(z,'pulse')[-1],'o',color=col,markersize=4)
 control=np.load(OUT/'L39_reduced_mca_memory_pulse_300.npz');b=control['pulse_field_mca_bound']+control['pulse_field_mca_unbound'];q=2*control['pulse_x']/39-1;mc=b@q/b.sum(1)
 ax.plot(control['pulse_time'],mc,color='#d58a2b',linestyle='--',label='L39 reduced initial MCA bias')
 ax.plot(control['pulse_time'][-1],mc[-1],'o',color='#d58a2b',markersize=4)
 ax.axhline(0,color='#888888',linestyle=':',linewidth=1)
 ax.set_title('B  Larger cells start with more old MCA bias');ax.set_xlabel('Time after cue onset');ax.set_ylabel('Total MCA first moment (right positive)');ax.legend(frameon=False,fontsize=6.3)
 ax.text(.04,.08,'By t = 170, native losses are both ≈ 0.085.\nReduced initial MCA bias: 300 passes, 173 / 250 fail.\nOriginal L39: the same 300-unit pulse fails.',transform=ax.transAxes,fontsize=6.3)
 ax=axes[0,2]
 L=18
 for label,ls in [('pass','-'),('fail','--')]:
  z=load(L,label);cont=np.load(OUT/f'L{L}_{label}_continuation.npz');t=np.r_[z['recovery_time'],10000+cont['time'][1:]];o=np.r_[z['recovery_orientation'],cont['orientation'][1:]]
  k=t>=9000;ax.plot(t[k],o[k],color=COL[label],linestyle=ls,label=label.title())
 ax.axhline(.05,color='#333333',linewidth=1,linestyle=':');ax.axvspan(9500,10000,color='#bbbbbb',alpha=.25)
 ax.set_title('C  L18: failure is just below the +0.05 gate');ax.set_xlabel('Basal recovery time');ax.set_ylabel('Rac orientation');ax.ticklabel_format(axis='y',style='plain',useOffset=False);ax.legend(frameon=False)
 ax.text(.04,.08,'Both end on the right by 30,000.\nThe weak endpoint makes the fixed gate demanding.',transform=ax.transAxes,fontsize=7)
 ax=axes[1,0]
 for L,col in [(23,'#007f78'),(39,'#785aa6')]:
  for label,ls in [('pass','-'),('fail','--')]:
   z=load(L,label);cont=np.load(OUT/f'L{L}_{label}_continuation.npz');t=np.r_[z['recovery_time'],10000+cont['time'][1:]];o=np.r_[z['recovery_orientation'],cont['orientation'][1:]]
   ax.plot(t,o,color=col,linestyle=ls,label=f'L{L} {label}')
 ax.axvspan(9500,10000,color='#bbbbbb',alpha=.25);ax.axhline(.05,color='#888888',linestyle=':',linewidth=.8);ax.set_ylim(-.085,.085)
 ax.set_title('D  L23: relaxation gate; L39: different fates');ax.set_xlabel('Basal recovery time');ax.set_ylabel('Rac orientation');ax.legend(frameon=True,facecolor='white',framealpha=.92,ncol=2,fontsize=6.3,loc='upper right')
 ax.text(.03,.10,'L23 fail: all tail samples > 0.05,\nbut range 0.00771 exceeds 0.005.\nL39 neighboring pulses choose opposite sides.',transform=ax.transAxes,fontsize=6.6,bbox=dict(facecolor='white',alpha=.88,edgecolor='none',pad=2))
 ax=axes[1,1]
 for L,col in [(23,'#007f78'),(39,'#785aa6')]:
  z=load(L);clamp=np.load(OUT/f'L{L}_mca_clamp.npz');ax.plot(z['recovery_time'],z['recovery_orientation'],color=col,label=f'L{L} native')
  ax.plot(clamp['recovery_time'],clamp['recovery_orientation'],color=col,linestyle='--',label=f'L{L} MCA fixed during cue')
 ax.set_ylim(-.08,.3);ax.axhline(0,color='#888888',linestyle=':',linewidth=.8);ax.set_title('E  Fast Rac alone cannot replace the memory');ax.set_xlabel('Basal recovery time');ax.set_ylabel('Rac orientation');ax.legend(frameon=False,fontsize=6.5)
 ax.text(.03,.51,'Diagnostic: suppress all MCA changes only\nduring the pulse; restore native recovery.\nRac still points right under the cue.',transform=ax.transAxes,fontsize=6.7)
 ax=axes[1,2]
 mode=r.read(OUT/'linear_mode_growth.json');ls=[x['length'] for x in mode['rows']];rates=[x['modes'][0]['dominant_growth_rate']*1000 for x in mode['rows']];ax.plot(ls,rates,color='#555555');ax.axhline(0,color='#888888',linestyle=':',linewidth=1)
 for L in [18,23,39]:
  y=next(x['modes'][0]['dominant_growth_rate']*1000 for x in mode['rows'] if x['length']==L);ax.plot(L,y,'o',color=COL['pass']);ax.annotate(str(L),(L,y),xytext=(4,5),textcoords='offset points')
 ax.set_title('F  Late amplification gets faster at larger L');ax.set_xlabel('Nominal length (µm)');ax.set_ylabel('Small basal polarity growth per 1,000 units')
 ax.text(.03,.05,'Local diagnostic around homogeneous rest.\nSeeded native runs reproduce these rates.\nLong-arm delay is not slower final amplification.',transform=ax.transAxes,fontsize=6.7)
 for ax in axes.flat:ax.grid(alpha=.12)
 fig.suptitle('deBelly: three different limits shape the apparent U-shaped response',fontweight='bold',fontsize=14)
 save(fig,'debelly_03_mechanism_diagnostics')

def main():
 audit();snapshots();recovery();diagnostics()
 pages=[dict(pdf='debelly_01_on_cue_snapshots.pdf',png='debelly_01_on_cue_snapshots.png',main_claim='The fixed-slope cue turns Rac within 5–7.5 time units, but rewriting the MCA/mechanical memory requires a far longer pulse.',caption='Native spatial snapshots from the exact relative-20%-prepared checkpoints at L18, L23 (sampled duration minimum) and L39. Each column uses its shortest ordinary-orientation-passing pulse. All chemical, mechanical and solver parameters retain the manuscript protocol. The dotted bound-MCA line is the Rac-feedback threshold 0.068: low MCA promotes Rac. A rightward cue makes Rac high and Rho low at the right, changes cortical flow, and gradually redistributes MCA. At L23 and L39, MCA remains biased toward the right at cue cutoff, so successful reversal need not mean adhesion memory has already completely reversed.'),dict(pdf='debelly_02_basal_recovery_snapshots.pdf',png='debelly_02_basal_recovery_snapshots.png',main_claim='The scalar assay can call right reversal while coupled MCA, velocity and tension still evolve; extended basal recovery approaches the reflected preparation.',caption='Spatial profiles at 250, 1,000, 5,000, 10,000 and 30,000 units after withdrawal of each selected pulse. The first four times are within the unchanged manuscript basal recovery. The last is an additional diagnostic continuation with the original native basal dynamics. The ordinary selected scalar endpoint passes at all three lengths, but the independent 0.1% field-tail flag fails at L23 and L39. The 30,000-unit continuation is not a new preparation-style convergence/holdout certification.'),dict(pdf='debelly_03_mechanism_diagnostics.pdf',png='debelly_03_mechanism_diagnostics.png',main_claim='The falling arm is strongly influenced by weak-polarity/finite-recovery scalar gates; the rising arm requires rewriting a spatial adhesion/mechanical state; greater old MCA imbalance contributes but does not fully predict the threshold.',caption='A: Original duration curves and field flags. B: Unmodified total-MCA first moments during the pulse; at a common time of 170 both lengths have lost about 0.085, but L39 starts at 0.1798 versus 0.1225 at L23. A labeled altered-initial-state diagnostic reduces only the odd reflection components of bound/unbound MCA at L39 to match the initial L23 moment, retaining total mass, both even components and every other original field. That altered state reverses under a 300-unit pulse while the original same pulse fails; altered 173- and 250-unit pulses still fail. This establishes a contribution of the larger adhesion imprint but rules out first moment alone as a sufficient predictor. C: Adjacent L18 pulses both ultimately recover right; the original failing pulse is only about 0.000014 below +0.05 over the assay tail. D: At L23 both pulses recover right, but the failing pulse still changes by 0.00771 in the final 500 units (allowed <0.005); at L39 the two neighboring pulses approach opposite polarized branches. E: Explicit causal perturbation, not an assay modification: all MCA exchange, diffusion and advection are suppressed only during the cue; the exact original Rac/Rho/mechanical state, cue and unchanged native basal recovery remain. Strong right Rac under the cue no longer persists at L23/L39, supporting a necessary role for cue-induced adhesion-memory change in these examples. F: Validated local linearization about homogeneous basal rest, with measured native seeded trajectories agreeing; amplification is faster at L39 than at L23, so slower amplification cannot explain the rising duration arm. No universal length scaling law or isolated-MCA sufficiency claim is made.')]
 r.write(FIG/'assets.json',dict(generated_at=n.now(),model='debelly',ordered_pages=pages,ordered_pdf=[p['pdf'] for p in pages],ordered_png=[p['png'] for p in pages],input_audit=str((OUT/'independent_reconstruction_audit.json').relative_to(n.ROOT))))
 print('Saved',FIG)

if __name__=='__main__':main()
