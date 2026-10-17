"""Publication figures for native resting states and global tension→Rho ablation.

Read-only use of completed, independently audited corrected native campaigns.
No integration, reconditioning, or changes to source evidence.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT/'figures/length_mechanisms_20261008'
OUT = FIG/'narrative_20261009/resting_ablation'
DATA = ROOT/'data/runs/length_mechanisms_20261008/narrative_20261009/resting_ablation'
DEV = FIG/'dev/tension_rho_ablation/supporting'
RUN = ROOT/'data/runs/debelly_global_tension_rho_ablation_20261009'
BASE = ROOT/'data/runs/debelly_corrected_three_test_20261008'
COLORS = {23:'#0072B2', 31:'#762a83', 39:'#D55E00'}
FRONT, REAR = '#0072B2', '#777777'
NATIVE, ABLATION = '#333333', '#009E73'
SOURCES = {}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
    'axes.titlesize':12,'axes.labelsize':11,'axes.spines.top':False,
    'axes.spines.right':False,'xtick.labelsize':10,'ytick.labelsize':10,
    'legend.fontsize':10,'pdf.fonttype':42,'svg.fonttype':'none'})

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def source(path):
    path=Path(path); SOURCES[str(path.relative_to(ROOT))]=sha(path)
    return path

def read(path):
    return json.loads(source(path).read_text())

def rows(path):
    with source(path).open(newline='') as handle:
        return list(csv.DictReader(handle))

def writecsv(name, records):
    with (DATA/name).open('w',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=list(records[0])); w.writeheader();w.writerows(records)

def panel(ax, letter, title, ylabel, spatial=False):
    ax.set_title(title,loc='left',pad=9)
    ax.text(-.18,1.045,letter,transform=ax.transAxes,fontweight='bold',fontsize=13)
    ax.set_ylabel(ylabel);ax.grid(alpha=.14,lw=.6)
    if spatial:
        ax.set_xlim(0,1);ax.set_xticks([0,.5,1]);ax.set_xlabel('Position x/L')
    else:
        ax.set_xlim(17,40);ax.set_xticks([18,23,31,39]);ax.set_xlabel('Cell length (µm)')

def output(fig, stem, caption, audit):
    hashes={}
    for ext in ['pdf','svg','png']:
        path=OUT/f'{stem}.{ext}';fig.savefig(path,dpi=200,facecolor='white')
        hashes[str(path.relative_to(FIG))]=sha(path)
    plt.close(fig)
    (DATA/f'{stem}_caption.txt').write_text(caption+'\n')
    (OUT/f'{stem}_caption.txt').write_text(caption+'\n')
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',
        description=caption,legend=caption,sha256=hashes,source_sha256=SOURCES.copy(),
        renderer_sha256=sha(__file__),new_integrations=0,audit=audit,
        data_sha256={str(p.relative_to(ROOT)):sha(p) for p in DATA.iterdir() if p.suffix in ['.csv','.npz']},
        reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))])
    (DATA/f'{stem}.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/f'{stem}.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def resting():
    assert read(BASE/'results/independent_audit.json')['passed']
    manifest=read(BASE/'manifest.json')
    assert manifest['center_boundary_motion'] is False
    metrics=rows(FIG/'dev/resting_state/supporting/data/basal_endpoint_metrics_vs_length.csv')
    numerical=[{k:float(v) for k,v in r.items() if k not in ['checkpoint','checkpoint_sha256']} for r in metrics]
    threshold=numerical[0]['rac_switch_threshold']
    eligible=sorted({int(c['length']) for c in manifest['cases'] if c['model']=='debelly' and c['eligible']})
    assert eligible==[int(r['length_um']) for r in numerical]
    profiles={}; audits=[];spatial=[]
    for raw,r in zip(metrics,numerical):
        L=int(r['length_um']); cp=source(ROOT/raw['checkpoint']);assert sha(cp)==raw['checkpoint_sha256']
        with np.load(cp) as z:
            f={k[6:]:z[k].copy() for k in z.files if k.startswith('state_')};x=z['x'].copy()
            assert all(np.isfinite(y).all() for y in f.values())
            np.testing.assert_equal([f['rac'][0],f['rac'][-1]],[r['left_rac'],r['right_rac']])
            np.testing.assert_equal([f['tension'][0],f['tension'][-1]],[r['left_tension_pN'],r['right_tension_pN']])
            switch=.5*(1+np.tanh((f['rac'][[0,-1]]-threshold)/.01))
            command=.15*switch/(1+f['tension'][[0,-1]]**2/10)
            np.testing.assert_allclose(switch,[r['left_protrusion_switch'],r['right_protrusion_switch']],atol=1e-15)
            np.testing.assert_allclose(command,[r['derived_left_command_outward_um_per_s'],r['derived_right_command_outward_um_per_s']],atol=1e-15)
            # Global endpoint values, not first/last interior displacements.
            np.testing.assert_equal([f['left_boundary_displacement'][0],f['right_boundary_displacement'][0]],[r['left_endpoint_displacement_um'],r['right_endpoint_displacement_um']])
            if L in COLORS:
                profiles[L]=(x/L,f)
                spatial.extend(dict(length_um=L,x_over_L=float(xx),**{k:float(f[k][i]) for k in ['rac','rho','mca_bound','tension']}) for i,xx in enumerate(x/L))
        audits.append(dict(length_um=L,checkpoint_hash_verified=True,endpoint_fields_exact=True,command_and_switch_recomputed=True))
    writecsv('resting_endpoint_metrics.csv',metrics);writecsv('resting_selected_profiles.csv',spatial)
    fig,axs=plt.subplots(3,3,figsize=(12,9.5),layout='constrained')
    x=np.array(eligible)
    specs=[('Endpoint Rac','Rac (native conc.)','left_rac','right_rac'),
        ('Protrusion switch','Activation, S','left_protrusion_switch','right_protrusion_switch'),
        ('Endpoint command','Outward command (µm/s)','derived_left_command_outward_um_per_s','derived_right_command_outward_um_per_s'),
        ('Endpoint displacement','Outward displacement (µm)','left_endpoint_displacement_um','right_endpoint_displacement_um'),
        ('Endpoint tension','Tension (pN)','left_tension_pN','right_tension_pN')]
    for ax,(title,ylabel,left,right),letter in zip(axs.flat,specs,'ABCDE'):
        panel(ax,letter,title,ylabel)
        for field,color,ls,side in [(left,FRONT,'-','Front / left'),(right,REAR,'--','Rear / right')]:
            yy=np.array([r[field] for r in numerical]);yy=-yy if field=='left_endpoint_displacement_um' else yy
            ax.plot(x,yy,color=color,ls=ls,lw=2,marker='o',ms=3.5,label=side)
        if letter=='A':
            ax.axhline(threshold,color='#9c2525',ls=':',lw=1.4,label='Switch threshold')
            ax.legend(frameon=False,fontsize=8.5,loc='lower left',handlelength=2)
        if letter=='B':ax.set_ylim(-.015,1)
        if letter in 'CD':ax.set_ylim(bottom=-.005)
        if letter=='E':ax.axhline(0,color='#aaaaaa',lw=.8)
    for ax,(name,title,ylabel),letter in zip(list(axs.flat)[5:],[('rac','Rac profile','Rac (native conc.)'),('rho','Rho profile','Rho (native conc.)'),('mca_bound','Attachment profile','Bound MCA (native density)'),('tension','Tension profile','Tension (pN)')],'FGHI'):
        panel(ax,letter,title,ylabel,spatial=True)
        for L,color in COLORS.items():xx,f=profiles[L];ax.plot(xx,f[name],color=color,lw=2,label=f'{L} µm')
        if name=='rac':ax.axhline(threshold,color='#9c2525',ls=':',lw=1)
        if name=='tension':ax.axhline(0,color='#aaaaaa',lw=.8)
        if name=='rho':ax.legend(frameon=False,fontsize=10)
    caption=('Native resting states before the countercue. A–E show both endpoints at all 12 eligible corrected preparations (L = 18–39 µm); front means the old left pole and rear the right endpoint. A: Rac and its protrusion threshold, 0.8990616483 (= 1.3 times homogeneous basal Rac). B: switch S = 0.5 × [1 + tanh((Rac − threshold)/0.01)]. C: instantaneous outward endpoint command = 0.15 × S/[1 + T²/10], evaluated at the checkpoint. D: actual stored global endpoint displacement, expressed outward (left sign reversed for display); this is distinct from first/last interior-cell displacement. E: first/last-cell membrane tension. F–I: native Rac, Rho, bound-MCA and tension profiles for L = 23, 31 and 39 µm (blue, purple, orange); only position is divided by L. The dotted line in F repeats the Rac switch threshold. Independent endpoint commands are retained. The half-second update is u_next = 0.98 × u + 0.5 × v_command, with relaxation time −0.5/log(0.98) = 24.75 s; nonzero resting command balances endpoint relaxation. Each length starts from the exact qualified basal checkpoint shared by all three reversal assays. No new integration.')
    audit=dict(passed=True,eligible_lengths=eligible,threshold=threshold,checks=audits)
    output(fig,'resting_length_state',caption,audit)
    return dict(image='narrative_20261009/resting_ablation/resting_length_state.png',
        title='Length changes the mechanical starting state',
        paragraphs=['The deBelly cells enter the reversal assay with different mechanical starting conditions. As length increases, Rac at the old front approaches and crosses the steep protrusion switch. A small endpoint concentration difference therefore produces a large change in activation and outward command. The resting endpoint displacement grows with that command, while the rear remains constrained by its inactive switch.',
        'The spatial profiles show the chemical and mechanical state that accompanies this change: longer cells carry a stronger Rac–Rho asymmetry, a redistributed attachment field and a larger tension difference between the ends. These are the inherited conditions that a rightward cue must overturn. The next comparison separates how the mechanics responds to a prescribed command from how the native feedback selects that command.'],
        caption=caption,sources=[[os.path.relpath(DATA/'resting_endpoint_metrics.csv',FIG),'Endpoint measurements'],[os.path.relpath(DATA/'resting_selected_profiles.csv',FIG),'Selected profiles'],[os.path.relpath(DATA/'resting_length_state.manifest.json',FIG),'Figure provenance']],renderer=Path(__file__).name)

def ablation():
    campaign=read(RUN/'manifest.json'); provenance=read(RUN/'provenance.json')
    assert sha(RUN/'manifest.json')==provenance['manifest_sha256']
    assert read(RUN/'audits/summary.json')['passed']
    assert read(DEV/'metadata/global_tension_rho_ablation_independent_audit.json')['passed']
    for name,digest in provenance['frozen_sources'].items():assert sha(source(RUN/'frozen_sources'/name))==digest
    endpoint=rows(DEV/'data/global_tension_rho_ablation_endpoints.csv')
    by={(r['condition'],int(r['length_um'])):r for r in endpoint}
    final={};audits=[];outprofiles=[]
    for job,control in zip(campaign['jobs'],campaign['native_controls_reused']):
        L=int(job['length']);cp=source(RUN/'checkpoints'/f'L{L}.npz');assert sha(cp)==job['checkpoint_sha256']
        for cond,trace,record_path,commands_path in [('native',ROOT/control['trace'],ROOT/control['result'],ROOT/control['commands']),('ablation',RUN/'workers'/f"{job['slot']:02d}.npz",RUN/'workers'/f"{job['slot']:02d}.json",RUN/'workers'/f"{job['slot']:02d}.commands.npz")]:
            record=read(record_path); assert record['verified_steady_state'] and record['initial_state_exact']
            assert sha(source(trace))==record['trace_sha256']
            if cond=='native': assert sha(trace)==control['trace_sha256']
            with np.load(cp) as initial,np.load(trace) as z,np.load(source(commands_path)) as command:
                for name in initial.files:
                    if name.startswith('state_'):np.testing.assert_array_equal(z['field_'+name[6:]][0],initial[name])
                f={k[6:]:z[k][-1].copy() for k in z.files if k.startswith('field_')};xx=z['x']/L
                moment=float(f['rac']@(2*xx-1)/f['rac'].sum());contrast=float(np.ptp(f['rac']))
                np.testing.assert_allclose([moment,contrast,float(f['left_boundary_displacement'][0]),float(command['left_command'][-1])],[float(by[cond,L][k]) for k in ['rac_polarity','rac_contrast','left_boundary_displacement','left_command']],atol=1e-14)
                expected=-.15*.5*(1+np.tanh((f['rac'][0]-campaign['original_rac_threshold'])/.01))/(1+f['tension'][0]**2/10)
                np.testing.assert_allclose(command['left_command'][-1],expected,atol=1e-15)
                if cond=='ablation':
                    np.testing.assert_array_equal(z['field_delta_beta_tension'],0.)
                    np.testing.assert_array_equal(z['field_rho_activation'],1.)
                final[cond,L]=(xx,f)
                outprofiles.extend(dict(length_um=L,condition=cond,x_over_L=float(x),**{k:float(f[k][i]) for k in f}) for i,x in enumerate(xx))
                audits.append(dict(length_um=L,condition=cond,exact_restart_verified=True,native_trace_hash_verified=True,command_recomputed=True,orientation_recomputed=moment,contrast_recomputed=contrast,verified_steady_state=True))
    u={cond:np.array([-float(by[cond,L]['left_boundary_displacement']) for L in COLORS]) for cond in ['native','ablation']}
    spread={cond:float(np.ptp(a)) for cond,a in u.items()};reduction=1-spread['ablation']/spread['native']
    assert abs(reduction-.5038045341972212)<1e-12
    writecsv('ablation_endpoint_metrics.csv',endpoint);writecsv('ablation_complete_final_profiles.csv',outprofiles)
    # Retain the complete native global histories and all-field exports with provenance.
    for name in ['global_tension_rho_ablation_endpoint_timecourses.csv','global_tension_rho_ablation_initial_final_all_fields.csv']:
        shutil.copyfile(source(DEV/'data'/name),DATA/name)
    fig=plt.figure(figsize=(12,7.25),layout='constrained')
    grid=fig.add_gridspec(2,12,height_ratios=[1,1.22])
    top=[fig.add_subplot(grid[0,3*i:3*(i+1)]) for i in range(4)]
    specs=[('Front displacement','Outward displacement (µm)','left_boundary_displacement',-1),('Front command','Outward command (µm/s)','outward_left_command',1),('Rac orientation','Leftward moment magnitude','rac_polarity',-1),('Rac contrast','max(Rac) − min(Rac)','rac_contrast',1)]
    for ax,(title,ylabel,field,sign),letter in zip(top,specs,'ABCD'):
        panel(ax,letter,title,ylabel);ax.set_xlim(21.5,40.5);ax.set_xticks([23,31,39])
        for cond,color,style in [('native',NATIVE,'--'),('ablation',ABLATION,'-')]:
            values=[sign*float(by[cond,L][field]) for L in COLORS]
            ax.plot(list(COLORS),values,color=color,ls=style,lw=2,marker='o',ms=5,label='Native' if cond=='native' else 'Tension→Rho removed')
        if letter=='A':ax.text(.04,.94,f'{100*reduction:.1f}% smaller spread',transform=ax.transAxes,ha='left',va='top',color=ABLATION,fontsize=10,fontweight='bold')
        if letter=='C':ax.axhline(.05,color='#999999',ls=':',lw=1)
        if letter=='D':ax.set_ylim(0,.36)
    for i,(L,color) in enumerate(COLORS.items()):
        ax=fig.add_subplot(grid[1,4*i:4*(i+1)])
        panel(ax,chr(69+i),f'Rac profile · L = {L} µm','Rac (native conc.)',True)
        for cond,ls in [('native','--'),('ablation','-')]:
            xx,f=final[cond,L];ax.plot(xx,f['rac'],color=color,ls=ls,lw=2.2,label='Native' if cond=='native' else 'Tension→Rho removed')
        ax.set_ylim(.53,.93)
    fig.legend(handles=[Line2D([],[],color=NATIVE,ls='--',lw=2,label='Native'),Line2D([],[],color=ABLATION,lw=2,label='Tension→Rho removed')],loc='outside lower center',ncol=2,frameon=False)
    caption=('Removing tension→Rho feedback reduces resting mechanical length dependence while preserving asymmetric Rac profiles. Native controls (dashed) and globally ablated states (solid) at L = 23, 31 and 39 µm start from identical corrected basal checkpoints, without a countercue. A: actual global front displacement, expressed outward; its max-minus-min spread falls from 1.968526 to 0.976774 µm (50.4% reduction). B: instantaneous outward front command; the native direct tension inhibition of protrusion remains. C: magnitude of the leftward normalized Rac first moment; dotted line is the separate 0.05 orientation cutoff. D: Rac spatial contrast, max minus min. E–G: complete final Rac profiles on x/L, with a common concentration axis and length colors blue/purple/orange. The L = 31 and 39 ablations fall below the orientation cutoff but retain clear spatial contrast; all three are left biased, with contrast 0.205–0.210. Only the tension-dependent increment to Rho activation is removed globally. Original parameters, homogeneous reference, Rac threshold, MCA→Rac feedback, Rac–Rho mutual inhibition, Rho active stress, direct protrusion tension inhibition and independent endpoint relaxation are retained. All three ablations pass two 1,000-s settling windows and a 10,000-s holdout; final times are 15,010, 16,010 and 18,010 s. This independently evolved perturbation demonstrates a contribution of tension feedback to resting displacement differences; it does not measure a contribution to reversal duration.')
    audit=dict(passed=True,displacement_spread_um=spread,spread_reduction_fraction=reduction,checks=audits)
    (DATA/'ablation_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    output(fig,'tension_rho_ablation',caption,audit)
    return dict(image='narrative_20261009/resting_ablation/tension_rho_ablation.png',
        title='Tension feedback amplifies the resting displacement difference',
        paragraphs=['The feedback ablation tests whether tension-dependent Rho activation helps generate these different mechanical starting states. Removing this increment globally and allowing the same initial cells to relax reduces the displacement spread across L = 23, 31 and 39 µm by 50.4%, from 1.969 to 0.977 µm. The front commands also become more similar, supporting a substantial role for tension feedback in amplifying the resting mechanical difference.',
        'Rac remains visibly asymmetric at every length. Its spatial contrast settles near 0.21, and the length dependence of its orientation moment is greatly reduced. The profiles show a weaker but retained chemical asymmetry alongside reduced mechanical separation. Together with the equal-command comparison, this gives a concrete account of part of the resting length dependence: the mechanical response varies with length, and native tension feedback amplifies the commands that drive it.'],
        caption=caption,sources=[[os.path.relpath(DATA/'ablation_endpoint_metrics.csv',FIG),'Ablation endpoint measurements'],[os.path.relpath(DATA/'ablation_complete_final_profiles.csv',FIG),'Complete final profiles'],[os.path.relpath(DATA/'ablation_validation.json',FIG),'Independent validation'],[os.path.relpath(DATA/'tension_rho_ablation.manifest.json',FIG),'Figure provenance']],renderer=Path(__file__).name)

def main():
    OUT.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
    source(Path(__file__))
    chapters={'resting':resting(),'ablation':ablation()}
    (DATA/'chapter_overrides.json').write_text(json.dumps(chapters,indent=2)+'\n')
    print(json.dumps(dict(figures=str(OUT),chapters=str(DATA/'chapter_overrides.json')),indent=2))

if __name__=='__main__':main()
