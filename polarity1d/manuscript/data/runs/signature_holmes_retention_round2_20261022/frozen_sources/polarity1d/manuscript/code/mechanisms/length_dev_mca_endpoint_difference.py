"""Exact Figure8 pulses: bound MCA at the first minus last native cell."""
from pathlib import Path
from datetime import datetime,timezone
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_output_layout as layout
import length_dev_mca_traces as source
ROOT=Path(__file__).resolve().parents[2];OUT=source.OUT
def output(name):return layout.dev_output_path(OUT,name)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
SOURCES={str(Path(__file__).relative_to(ROOT)):sha(__file__),str(Path(layout.__file__).relative_to(ROOT)):sha(layout.__file__),str(Path(source.__file__).relative_to(ROOT)):sha(source.__file__)}
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def main():
 m=read(source.RUN/'manifest.json');assert m['center_boundary_motion'] is False
 assert read(ROOT/'data/runs/length_mechanisms_20261008/quality_20261009_debelly/independent_audit.json')['passed']
 traces={};rows=[];arrays={}
 for slot,j in enumerate(m['jobs']):
  c=j['case'];L=int(c['length']);geo=c['condition'];D=float(j['duration'])
  path=source.RUN/f'workers/{slot:03d}.npz' if j['existing_pulse_slot'] is None else source.OLD/f"workers/{j['existing_pulse_slot']:03d}.npz"
  meta=read(path.with_suffix('.json'));assert sha(path)==meta['archive_sha256'];SOURCES[str(path.relative_to(ROOT))]=sha(path)
  cp=ROOT/c['checkpoint'];assert sha(cp)==c['checkpoint_sha256'];SOURCES[str(cp.relative_to(ROOT))]=sha(cp)
  with np.load(path) as z:
   x=z['x'];tp=z['pulse_time'];tr=z['recovery_time'];bp=z['pulse_field_mca_bound'];br=z['recovery_field_mca_bound']
  with np.load(cp) as z:np.testing.assert_array_equal(bp[0],z['state_mca_bound'])
  np.testing.assert_array_equal(br[0],bp[-1]);assert tr[-1]==10000 and abs(tp[-1]-D)<1e-9
  ep=bp[:,0]-bp[:,-1];er=br[:,0]-br[:,-1]
  t=np.r_[tp,D+tr[1:]];e=np.r_[ep,er[1:]]
  moment=np.r_[-source.polarity(bp,x,L),-source.polarity(br,x,L)[1:]]
  rows.append(dict(length_um=L,condition=geo,cue_duration_s=D,first_sample_x_um=float(x[0]),last_sample_x_um=float(x[-1]),initial_bound_mca_left_minus_right=float(ep[0]),at_withdrawal=float(ep[-1]),final=float(er[-1]),first_endpoint_difference_zero_from_onset_s=source.crossing(t,e),persistent_positive_endpoint_difference_from_onset_s=source.persistent(t,e>=0),first_sign_reversed_mca_moment_zero_from_onset_s=source.crossing(t,moment),persistent_positive_sign_reversed_mca_moment_from_onset_s=source.persistent(t,moment>=0),source_archive=str(path.relative_to(ROOT))))
  traces[geo,L]=(tp,D+tr,ep,er,D)
  for key,value in dict(time=t,left_minus_right=e,minus_mca_polarity=moment).items():arrays[f'{geo}_L{L}_{key}']=value
 fig,axs=plt.subplots(2,3,figsize=(12.7,8.1),layout='constrained',sharey=True)
 for col,geo in enumerate(source.GEOS):
  for row in [0,1]:
   ax=axs[row,col];ax.text(-.14,1.04,chr(65+3*row+col),transform=ax.transAxes,fontsize=12,fontweight='bold',fontfamily='DejaVu Sans');ax.grid(alpha=.13);ax.axhline(0,color='#888888',lw=.8)
   if row==0:ax.set_title(source.TITLES[col],pad=10)
   for L in source.LENGTHS:
    tp,tr,ep,er,D=traces[geo,L];color=source.COLORS[L]
    ax.plot(tp,ep,color=color,lw=1.4);ax.plot(tr,er,color=color,lw=1.3,ls='--');ax.plot(D,ep[-1],'o',ms=4,color=color,zorder=4)
   if col==0:ax.set_ylabel('Bound MCA: left − right (native concentration)')
   if row==0:ax.set_xlim(0,1800);ax.set_xticks([0,600,1200,1800]);ax.set_xlabel('Time from cue onset (s)')
   else:ax.set_xscale('symlog',linthresh=10,linscale=.7);ax.set_xlim(0,max(r['cue_duration_s'] for r in rows)+10000);ax.set_xticks([0,10,100,1000,10000],['0','10','100','1,000','10,000']);ax.set_xlabel('Time from cue onset (s; log above 10 s)')
  axs[0,col].margins(y=.06)
 handles=[Line2D([],[],color=source.COLORS[L],lw=2,label=str(L)) for L in source.LENGTHS]+[Line2D([],[],color='#333333',lw=1.5,label='Cue'),Line2D([],[],color='#333333',ls='--',lw=1.5,label='Recovery'),Line2D([],[],color='#333333',marker='o',ls='',label='Withdrawal')]
 fig.legend(handles=handles,loc='outside lower center',ncol=8,frameon=False,title='Cell length (µm) · exact Figure 8 selected pulses')
 hashes={}
 for ext in ['png','pdf','svg']:
  p=output('MCA_endpoint_difference_traces.'+ext);fig.savefig(p,dpi=220,facecolor='white');hashes[layout.dev_artifact_key(OUT,p)]=sha(p)
 plt.close(fig)
 p=output('MCA_endpoint_difference_traces.npz');np.savez_compressed(p,**arrays);hashes[layout.dev_artifact_key(OUT,p)]=sha(p)
 p=output('MCA_endpoint_difference_traces_summary.csv')
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 hashes[layout.dev_artifact_key(OUT,p)]=sha(p)
 caption='Bound MCA at the old left end minus the opposite right end for the exact 36 Figure 8 selected cue/recovery histories: all 12 eligible lengths, with fixed 3 µm patch, relative 20% patch, and slope 0.24 gradient in columns. End values are the first and last native finite-volume cell centres, x=L/402 and x=L−L/402, rather than interpolated boundary concentrations. Positive difference indicates more bound MCA at the new rear (left end). Colors consistently identify length; solid curves are cue-on, dashed curves are recovery, and plain circles mark each actual withdrawal. A–C resolve the first 1800 s on a linear axis; D–F retain every recorded sample through 10000 s after withdrawal, using a linear axis from 0 to 10 s and logarithmic thereafter. The selected cues are the previously retained shortest tested pulses satisfying the orientation selection gate; this plot adds no success or stability encoding. The endpoint difference is a concentration difference, distinct from the normalized spatial MCA polarity moment. No integrations, interpolation, new cue selection or report changes.'
 p=output('MCA_endpoint_difference_traces_caption.txt');p.write_text(caption+'\n');hashes[layout.dev_artifact_key(OUT,p)]=sha(p)
 result=dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',description=caption,legend=caption,sha256=hashes,source_sha256=SOURCES,renderer_sha256=sha(__file__),reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))],native_histories=36,new_integrations=0,current_report_changed=False,time_resolution='Exact retained samples: pulse 0.5 or 1 s, recovery 2.5 s; zero-crossing summary interpolates only adjacent samples; persistent sign is first recorded sample positive through final sample.')
 output('MCA_endpoint_difference_traces.manifest.json').write_text(json.dumps(result,indent=2)+'\n')
 print(output('MCA_endpoint_difference_traces.png'))
 for geo in source.GEOS:
  rr=[r for r in rows if r['condition']==geo];delta=[r['first_endpoint_difference_zero_from_onset_s']-r['first_sign_reversed_mca_moment_zero_from_onset_s'] for r in rr if r['first_endpoint_difference_zero_from_onset_s'] is not None and r['first_sign_reversed_mca_moment_zero_from_onset_s'] is not None]
  print(geo,'endpoint-minus-moment first-zero timing:',min(delta),np.median(delta),max(delta))
if __name__=='__main__':main()
