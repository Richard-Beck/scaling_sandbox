"""Development-only inherited basal endpoint state across qualified lengths.

Reads exact corrected shared pre-countercue checkpoints; no integration.
"""
from pathlib import Path
import length_dev_output_layout as devlayout
from datetime import datetime,timezone
from dataclasses import asdict
import csv,hashlib,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'code/core')]
import debelly_corrected_reversal as corrected
n=corrected.n
BASE=corrected.DEFAULT
OUT=ROOT/'figures/length_mechanisms_20261008/dev/resting_state'
def output(name):return devlayout.dev_output_path(OUT,name)
def artifact_key(path):return devlayout.dev_artifact_key(OUT,path)
LENGTHS=[18,19,21,23,25,27,29,31,33,35,37,39]
GEOS=['fixed_physical','fixed_relative','fixed_slope']
SOURCES={str(Path(devlayout.__file__).relative_to(ROOT)):hashlib.sha256(Path(devlayout.__file__).read_bytes()).hexdigest()}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 p=n.model_for('debelly').parameters
 assert p.center_boundary_motion is False and p.assay_mode=='reflection_symmetric_extension'
 assert p.rac_to_rho_coupling==1. and p.protrusion_speed==-.15 and p.protrusion_switch_scale==.01 and p.protrusion_tension_scale==10.
 frozen=read(BASE/'provenance.json');verified={}
 for name in ['code/core/native.py','code/core/models/debelly.py','code/core/src/polarity1d/debelly.py','code/debelly_corrected_reversal.py']:
  digest=frozen['sources'][name];assert sha(BASE/'frozen_sources'/name)==digest and sha(ROOT/name)==digest,name;SOURCES[name]=digest;verified[name]=digest
 SOURCES[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 manifest=read(BASE/'manifest.json');assert manifest['center_boundary_motion'] is False
 assert read(BASE/'results/independent_audit.json')['passed']
 eligible={int(c['length']) for c in manifest['cases'] if c['model']=='debelly' and c['eligible']};assert eligible==set(LENGTHS)
 R,H=n.ds._homogeneous_rac_rho(p,.5);threshold=1.3*R
 assert abs(threshold-.899061648345644)<1e-13
 rows=[];audits=[]
 for L in LENGTHS:
  cases=[next(c for c in manifest['cases'] if c['model']=='debelly' and c['length']==L and c['condition']==geo) for geo in GEOS]
  assert all(c['eligible'] and c['center_boundary_motion'] is False for c in cases)
  assert len({c['checkpoint'] for c in cases})==1 and len({c['checkpoint_sha256'] for c in cases})==1
  c=cases[0];cp=ROOT/c['checkpoint'];assert sha(cp)==c['checkpoint_sha256'];SOURCES[str(cp.relative_to(ROOT))]=sha(cp)
  with np.load(cp) as z:fields={k:z['state_'+k].copy() for k in n.DB_NAMES}
  for key,value in fields.items():assert np.isfinite(value).all()
  for key in ['left_boundary_displacement','right_boundary_displacement','raw_left_displacement','raw_right_displacement']:assert np.ptp(fields[key])==0
  np.testing.assert_array_equal(fields['left_boundary_displacement'],fields['raw_left_displacement']);np.testing.assert_array_equal(fields['right_boundary_displacement'],fields['raw_right_displacement'])
  left_rac=float(fields['rac'][0]);right_rac=float(fields['rac'][-1]);left_tension=float(fields['tension'][0]);right_tension=float(fields['tension'][-1])
  sl=float(n.ds.smooth_threshold(left_rac,p.protrusion_switch_scale,threshold));sr=float(n.ds.smooth_threshold(right_rac,p.protrusion_switch_scale,threshold))
  np.testing.assert_allclose([sl,sr],.5*(1+np.tanh((np.array([left_rac,right_rac])-threshold)/.01)),rtol=0,atol=1e-15)
  il=1/(1+left_tension**2/p.protrusion_tension_scale);ir=1/(1+right_tension**2/p.protrusion_tension_scale)
  vl=p.rac_to_rho_coupling*p.protrusion_speed*sl*il;vr=-p.rac_to_rho_coupling*p.protrusion_speed*sr*ir
  rows.append(dict(length_um=L,homogeneous_rac=R,homogeneous_rho=H,rac_switch_threshold=threshold,left_rac=left_rac,right_rac=right_rac,left_protrusion_switch=sl,right_protrusion_switch=sr,left_tension_pN=left_tension,right_tension_pN=right_tension,left_tension_inhibition_factor=il,right_tension_inhibition_factor=ir,derived_left_command_signed_um_per_s=vl,derived_left_command_outward_um_per_s=-vl,derived_right_command_signed_um_per_s=vr,derived_right_command_outward_um_per_s=vr,left_command_before_tension_inhibition_outward_um_per_s=-p.protrusion_speed*sl,left_endpoint_displacement_um=float(fields['left_boundary_displacement'][0]),right_endpoint_displacement_um=float(fields['right_boundary_displacement'][0]),left_equilibrium_displacement_from_current_command_um=.5*vl/(1-p.boundary_relaxation),left_boundary_one_step_displacement_change_um=(p.boundary_relaxation-1)*float(fields['left_boundary_displacement'][0])+.5*vl,left_first_cell_displacement_um=float(fields['displacement'][0]),right_last_cell_displacement_um=float(fields['displacement'][-1]),left_first_cell_cortical_velocity_um_per_s=float(fields['velocity'][0]),right_last_cell_cortical_velocity_um_per_s=float(fields['velocity'][-1]),checkpoint=str(cp.relative_to(ROOT)),checkpoint_sha256=sha(cp)))
  audits.append(dict(length=L,all_three_assays_share_exact_checkpoint=True,checkpoint_hash_pass=True,corrected_independent_endpoints=True,native_global_boundary_fields_verified=True))
 csvpath=output('basal_endpoint_metrics_vs_length.csv')
 with csvpath.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 fig,axs=plt.subplots(1,5,figsize=(15.5,4.2),layout='constrained',sharex=True)
 x=[r['length_um'] for r in rows];blue='#0072B2';gray='#777777'
 for ax,letter in zip(axs,'ABCDE'):
  ax.text(-.15,1.04,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',fontfamily='DejaVu Sans');ax.grid(alpha=.13);ax.set_xlim(17,40);ax.set_xticks([18,23,29,35,39]);ax.set_xlabel('Cell length (µm)')
 axs[0].plot(x,[r['left_rac'] for r in rows],color=blue,marker='o',ms=4,label='Old left pole');axs[0].axhline(threshold,color=gray,ls=':',lw=1.1,label='Protrusion threshold');axs[0].set_ylabel('Left-end Rac (native concentration)');axs[0].legend(loc='lower right',frameon=False,fontsize=7.5)
 axs[1].plot(x,[r['left_protrusion_switch'] for r in rows],color=blue,marker='o',ms=4);axs[1].set_ylabel('Left protrusion switch, S');axs[1].set_ylim(0,1)
 axs[2].plot(x,[r['derived_left_command_outward_um_per_s'] for r in rows],color=blue,marker='o',ms=4,label='With native tension inhibition');axs[2].plot(x,[r['left_command_before_tension_inhibition_outward_um_per_s'] for r in rows],color=gray,ls='--',lw=1.4,label='Before tension inhibition');axs[2].set_ylabel('Derived outward left command (µm/s)');axs[2].set_ylim(bottom=0);axs[2].legend(loc='upper left',frameon=False,fontsize=7.3)
 axs[3].plot(x,[r['left_endpoint_displacement_um'] for r in rows],color=blue,marker='o',ms=4);axs[3].axhline(0,color=gray,lw=.7);axs[3].set_ylabel('Native left endpoint displacement (µm)')
 axs[4].plot(x,[r['right_tension_pN'] for r in rows],color=blue,marker='o',ms=4,label='Right end');axs[4].plot(x,[r['left_tension_pN'] for r in rows],color=gray,ls='--',marker='s',ms=3,label='Left end');axs[4].set_ylabel('End tension (pN)');axs[4].axhline(0,color=gray,lw=.7);axs[4].legend(loc='upper left',frameon=False,fontsize=8)
 outputs={}
 for ext in ['png','pdf','svg']:
  path=output('basal_endpoint_metrics_vs_length.'+ext);fig.savefig(path,dpi=220,facecolor='white');outputs[artifact_key(path)]=sha(path)
 plt.close(fig)
 caption='Inherited basal endpoint state before the countercue, at all 12 qualifying corrected lengths (18–39 µm). All three response assays share the exact same pre-cue checkpoint at each length; unqualified preparations are excluded. A: Rac at the first native cell beside the old left pole; the dotted line is the native protrusion threshold, 1.3 times homogeneous basal Rac = '+f'{threshold:.15g}'+'. B: left protrusion switch S = 0.5 times [1 + tanh((Rac_left − threshold)/0.01)]. C: derived instantaneous outward left command, 0.15 times S divided by [1 + T_left squared/10]; the dashed comparison removes only that algebraic tension divisor. Outward magnitude is positive, whereas the signed left command is negative. D: exact stored global left endpoint displacement, with its native negative sign; this is distinct from displacement at the first interior cell. E: right-end tension (blue) and left-end tension (gray dashed), sampled at the last and first native cells. The model uses independent endpoint commands without centering. Derived commands evaluate the current endpoint law at the checkpoint; they are not stored previous-step velocity or basal endpoint drift. The native half-second update is u_next = 0.98 times u + 0.5 times v_command, so stationary endpoints can retain nonzero command and displacement. No new integration or interpolation.'
 (output('basal_endpoint_metrics_vs_length_caption.txt')).write_text(caption+'\n')
 outputs[artifact_key(csvpath)]=sha(csvpath)
 rec=dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',description=caption,legend=caption,sha256=outputs,source_sha256=SOURCES,renderer_sha256=sha(Path(__file__)),parameters=asdict(p),parameter_identity=dict(current_native_sources_match_corrected_campaign_frozen_sources=verified,center_boundary_motion=False,homogeneous_step_s=.5,homogeneous_rac=R,homogeneous_rho=H,protrusion_threshold=threshold),independent_audit=audits,reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))],new_integrations=0,current_report_changed=False)
 (output('basal_endpoint_metrics_vs_length.manifest.json')).write_text(json.dumps(rec,indent=2)+'\n')
 for r in rows:print(r['length_um'],r['left_rac'],r['left_protrusion_switch'],r['derived_left_command_signed_um_per_s'],r['left_endpoint_displacement_um'],r['right_tension_pN'])
 print(output('basal_endpoint_metrics_vs_length.png'))
if __name__=='__main__':main()
