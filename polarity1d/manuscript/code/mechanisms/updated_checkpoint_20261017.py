"""Assemble one updated review checkpoint from checked component exports.

No model integration. Linear gain includes all 97 retained basals; recovery
contains only the 15 measured active-deactivation pilot basals. Ring recovery
baselines/maps remain explicitly unmeasured in this checkpoint.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv,json,hashlib,os,sys
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import length_dev_output_layout as layout
import expanded_ring_components_20261017 as ring

ROOT=Path(__file__).resolve().parents[2];REPO=ROOT.parents[1]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/updated_checkpoint_20261017'
OLD=ROOT/'figures/length_mechanisms_20261008/dev/integrated_framework_candidate_20261017'
GAIN=ROOT/'figures/length_mechanisms_20261008/dev/all_five_linear_function_20261017'
RING=ROOT/'figures/length_mechanisms_20261008/dev/expanded_ring_components_20261017'
GLOBAL=ROOT/'figures/length_mechanisms_20261008/dev/global_polarity_vector_20261017'
PILOT=ROOT/'data/runs/all_five_incumbent_recovery_pilot_20261017/results'
MODELS=('wave_pinning','goryachev','debelly','holmes_model3','spring')
NAMES={'wave_pinning':'Wave-pinning','goryachev':'Goryachev','debelly':'deBelly','holmes_model3':'Holmes 3','spring':'Spring'}
COLORS=dict(zip(MODELS,('#3976a8','#c96926','#2a9288','#8560a4','#a34f7f')))
REFS={'wave_pinning':19.,'goryachev':19.,'debelly':19.,'holmes_model3':20.,'spring':1.}
SOURCES={};ARTIFACTS=[]

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(2**20),b''):h.update(block)
    return h.hexdigest()

def source(path):
    path=Path(path).resolve();SOURCES[str(path.relative_to(REPO))]=sha(path);return path

def rows(path):
    with source(path).open() as stream:return list(csv.DictReader(stream))

def save(fig,name):
    for ext in ('.png','.pdf','.svg'):
        p=layout.dev_output_path(OUT,name+ext);fig.savefig(p,dpi=200 if ext=='.png' else None,bbox_inches='tight');ARTIFACTS.append(p)
    plt.close(fig)

def export(name,data):
    p=layout.dev_output_path(OUT,name)
    with p.open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    ARTIFACTS.append(p)

def figure1(folds,gains,profiles,clocks):
    fig=plt.figure(figsize=(9.4,8.4))
    a=fig.add_axes([.13,.65,.28,.25]);b=fig.add_axes([.56,.65,.39,.25])
    for key,marker,label,offset in [('excess_active_fold','o','Excess active amount',-.19),('width_fold','s','Core width',0),('peak_concentration_fold','^','Peak concentration',.19)]:
        a.scatter([float(r[key]) for r in folds],np.arange(5)+offset,marker=marker,s=25,facecolor='#333333' if marker=='o' else 'white',edgecolor='#333333',lw=.8,label=label)
    labels=[r['model_label']+'\n'+('L₀ 1→2' if r['model']=='spring' else f"{r['nominal_size_small']}→{r['nominal_size_large']}") for r in folds]
    a.set(yticks=np.arange(5),yticklabels=labels,xlim=(.8,3.2),xticks=[1,2,3],ylim=(4.6,-.6),xlabel='Larger / smaller value')
    a.tick_params(axis='y',length=0,labelsize=7.3);a.axvline(1,color='#aaaaaa',ls='--',lw=.7)
    a.set_title('A  Incumbent organization',loc='left',fontsize=9,pad=9)
    fig.legend(*a.get_legend_handles_labels(),loc='upper left',bbox_to_anchor=(.12,.992),frameon=False,fontsize=7.2,labelspacing=.25)
    for model,marker in zip(MODELS,('o','s','D','^','v')):
        rr=[r for r in gains if r['model']==model]
        b.plot([r['relative_nominal_size'] for r in rr],[r['gain_over_own_reference'] for r in rr],color=COLORS[model],marker=marker,ms=2.8,lw=1.1,label=NAMES[model],mec='white',mew=.2)
    b.axhline(1,color='#bdbdbd',ls='--',lw=.7);b.set_yscale('log',base=2)
    b.set(xlim=(.1,2.75),ylim=(.5,8),xticks=[.5,1,1.5,2,2.5],yticks=[.5,1,2,4,8],yticklabels=['½','1','2','4','8'],
        xlabel='Nominal size / own reference size',ylabel='Cue-weighted gain / own reference\n(log₂ scale)')
    b.set_title('B  Initial opposite-end input gain: all 97 basals',loc='left',fontsize=9,pad=9)
    fig.legend(*b.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.75,.992),frameon=False,ncol=3,fontsize=7.2,columnspacing=1.1,handlelength=1.6)
    b.text(.98,.94,'Refs: L19; Holmes L20; Spring resting L₀ = 1\nRelative nominal size is not a physical calibration.',transform=b.transAxes,ha='right',va='top',fontsize=6.5,linespacing=1.5)
    fig.text(.08,.546,'C  Restored local deficit over the model-specific early window',fontsize=9,ha='left')
    fig.text(.08,.521,'Exactly three measured pilot sizes per model; common deactivation rule, different fixed native clocks.',fontsize=7.2,ha='left')
    for left,model in zip((.09,.27,.45,.63,.81),MODELS):
        ax=fig.add_axes([left,.355,.13,.115]);block=next(r for r in clocks['models'] if r['model']==model)
        assert block['status']=='available';tau=float(block['selected_horizon_native']);rr=block['primary_results']
        xx=[float(r['size_native']) for r in rr];yy=[float(r['primary_F']) for r in rr];ee=[3*float(r['primary_E']) for r in rr]
        ax.errorbar(xx,yy,yerr=ee,fmt='o',color=COLORS[model],ms=4,lw=.7,capsize=2)
        pad=(max(xx)-min(xx))*.12
        ax.set(xlim=(min(xx)-pad,max(xx)+pad),ylim=(-.01,.21),xticks=xx,xticklabels=[f'{x:g}' for x in xx],yticks=[0,.1,.2])
        ax.set_title(f'{NAMES[model]}\nτ = {tau:.4g} native',fontsize=7.5,pad=7)
        ax.tick_params(labelsize=6.7);ax.axhline(0,color='#dddddd',lw=.6)
        if model=='wave_pinning':ax.set_ylabel('Initial deficit restored\n(fraction)',fontsize=7)
    fig.text(.51,.306,'Native pilot sizes (Spring: resting L₀). Points are unjoined; compare size dependence within each model.\nBars show empirical numerical sensitivity (3E), not biological confidence intervals.',ha='center',fontsize=7,linespacing=1.4)
    fig.text(.08,.254,'D  Native basal profiles: the same small/large pairs as A',ha='left',fontsize=9)
    fig.legend(handles=[Line2D([],[],color='#969696',ls='--',lw=1.2,label='Smaller'),Line2D([],[],color='#222222',lw=1.3,label='Larger')],loc='center',bbox_to_anchor=(.8,.257),frameon=False,ncol=2,fontsize=7)
    for left,fold in zip((.09,.27,.45,.63,.81),folds):
        ax=fig.add_axes([left,.09,.13,.10]);model=fold['model'];gg=[r for r in profiles if r['model']==model]
        for size,col,style in ((float(fold['nominal_size_small']),'#969696','--'),(float(fold['nominal_size_large']),'#222222','-')):
            rr=sorted([r for r in gg if float(r['nominal_size'])==size],key=lambda r:float(r['physical_position']))
            ax.plot([float(r['physical_position']) for r in rr],[float(r['active_concentration']) for r in rr],color=col,ls=style,lw=1.1)
        L=float(fold['physical_size_large']);ax.set(xlim=(0,L),ylim=(0,1.1*max(float(r['active_concentration']) for r in gg)),xticks=[0,L],xticklabels=['0',f'{L:.2f}' if model=='spring' else f'{L:g}'])
        ax.set_title(NAMES[model],fontsize=7.5,pad=5);ax.tick_params(labelsize=6.5)
        if model=='wave_pinning':ax.set_ylabel('Active concentration\n(native)',fontsize=7)
    fig.text(.51,.028,'Physical interval coordinate (native units); profile concentration scales differ across models.\nAll panels use linear 1D states. Spring profiles use measured physical domains 1.049→2.049.',ha='center',fontsize=7,linespacing=1.4)
    save(fig,'updated_checkpoint_figure1_linear')

def figure2_response(summary,gainmap):
    fig=plt.figure(figsize=(9.4,8.2));fig.text(.51,.985,'Figure 2 · Ring correspondence and local response',ha='center',va='top',fontsize=11,weight='bold')
    fig.text(.51,.95,'Continued review sheet 1 of 2 · retained measurements plus an explicit recovery-map gap',ha='center',fontsize=8)
    ax=fig.add_axes([.06,.70,.90,.17]);ax.axis('off')
    observables=list(dict.fromkeys(r['observable'] for r in summary));observables=[x for x in observables if 'Opposing-point' not in x]
    table=[]
    for name in observables:
        block=[r for r in summary if r['observable']==name];first=block[0]
        cells=[name,f"{float(first['expected_ring_over_line']):g}"]
        for model in ring.MODELS:
            r=next(r for r in block if r['model']==model);cells.append(f"{float(r['measured_ratio_min']):.4f}–{float(r['measured_ratio_max']):.4f}")
        table.append(cells)
    tab=ax.table(cellText=table,colLabels=['A  Ring / interval observable','Expected','Wave-pinning','Goryachev','Holmes 3'],cellLoc='center',colWidths=[.34,.09,.19,.19,.19],bbox=[0,0,1,1])
    tab.auto_set_font_size(False);tab.set_fontsize(6.7)
    for (i,j),cell in tab.get_celld().items():
        cell.set_edgecolor('#d4d9dd');cell.set_linewidth(.5)
        if i==0:cell.set_facecolor('#e8eef2');cell.set_text_props(weight='bold')
        elif j==0:cell.set_text_props(ha='left')
    fig.text(.06,.646,'13 independently prepared pairs at C = 2L: width/amount double; peak and matched-cue gain agree.\nBaseline input matching: endpoint support 3 ↔ full ring support 6. The retained finite challenges below use full support 3.',fontsize=7,linespacing=1.4)
    sizes=sorted({float(r['circumference']) for r in gainmap})
    for left,model in zip((.07,.385,.70),ring.MODELS):
        ax=fig.add_axes([left,.285,.235,.285]);ring.map_cells(ax,gainmap,model,sizes,'gain_over_own_C30_antipodal',Normalize(0,3),plt.get_cmap('viridis'),highlights=True)
        ax.set_title('B  '+ring.NAMES[model],fontsize=8,pad=7);ax.set_xlabel('Shorter-route fraction d/C',fontsize=7);ax.set_ylabel('Sampled circumference C',fontsize=7);ax.tick_params(labelsize=6.4)
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,3),cmap='viridis'),cax=fig.add_axes([.95,.285,.012,.285]));cb.set_ticks([0,1,2,3]);cb.ax.tick_params(labelsize=6.5);cb.set_label('Gain / own C30 antipodal gain',fontsize=7)
    fig.text(.06,.213,'One fixed 0–3 scale: within-size gain varies by <0.09% across eligible placements. Numbered cells link to sheet 2.\nRows are discrete sampled circumferences; × = 198 cue/core overlaps; dark gray = 17 unavailable Holmes C80 basals.\nNative input/time units are not calibrated across models.',fontsize=7,linespacing=1.4,va='top')
    gap=fig.add_axes([.06,.055,.90,.105]);gap.set_facecolor('#f5f5f5');gap.set_xticks([]);gap.set_yticks([])
    for s in gap.spines.values():s.set_edgecolor('#bdbdbd')
    gap.text(.025,.76,'C  Common active-deactivation recovery baseline and separation maps: not measured yet',transform=gap.transAxes,fontsize=8.5,weight='bold',va='top')
    gap.text(.025,.39,'This checkpoint contains the 15-case linear clock pilot only. No full linear or ring recovery grid is inferred.\nThe basal ring coefficient has no future separation; a separation map needs the same kick on declared challenged states.',transform=gap.transAxes,fontsize=7,va='top',linespacing=1.5)
    save(fig,'updated_checkpoint_figure2_response')

def vector(rows,C,distance):
    rr=[r for r in rows if r['model']=='holmes_model3' and float(r['circumference'])==C and np.isclose(float(r['fraction']),distance/C)]
    return {key:np.array([float(r[key]) for r in rr]) for key in ('time','vector_x','vector_y','relative_magnitude')}

def figure2_fates(outcomes,display,vectors):
    fig=plt.figure(figsize=(9.4,10.6));fig.text(.51,.985,'Figure 2 · From local patterns to global activity polarity',ha='center',va='top',fontsize=11,weight='bold')
    fig.text(.51,.954,'Continued review sheet 2 of 2 · source-map cells 1/2/3 match the histories and vector curves below',ha='center',fontsize=8)
    sizes=sorted({float(r['circumference']) for r in outcomes})
    for left,model in zip((.07,.385,.70),ring.MODELS):
        ax=fig.add_axes([left,.715,.235,.19]);ring.map_cells(ax,outcomes,model,sizes,'final_location_association',Normalize(-.5,6.5),ListedColormap(ring.COLORS),categorical=True,highlights=True)
        ax.set_title('D  '+ring.NAMES[model],fontsize=8,pad=7);ax.set_xlabel('Shorter-route fraction d/C',fontsize=7);ax.set_ylabel('Sampled C',fontsize=7);ax.tick_params(labelsize=6.1)
    handles=[Patch(facecolor=ring.COLORS[j],label=ring.LABELS[j]) for j in range(5)]
    handles += [Line2D([],[],marker='x',ls='',color='#888888',label='Cue/core overlap excluded'),Patch(facecolor='#777777',label='Holmes C80 unavailable')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.51,.672),ncol=3,frameon=False,fontsize=6.4,handlelength=1.3,columnspacing=1.0,labelspacing=.5)
    vmax=max(float(display[f'example{i}_primary'].max()) for i in (1,2,3))
    for left,i,C,distance,label in zip((.07,.385,.70),(1,2,3),(30,40,40),(12,12,16),('Cue-associated single peak','Merged / bridged region','Two distinct sites')):
        ax=fig.add_axes([left,.355,.235,.225]);x=display[f'example{i}_x'];p=display[f'example{i}_primary'];t=display[f'example{i}_time']
        mesh=ax.pcolormesh(np.r_[x,x[0]+C],t,np.c_[p,p[:,0]],cmap='magma',vmin=0,vmax=vmax,shading='nearest',rasterized=True)
        ax.axvline(0,color='white',ls=':',lw=.6);ax.axvline(distance,color='#83d6e2',ls=':',lw=.6);ax.hlines(1000,-10,C-10,color='#83d6e2',ls='--',lw=.8)
        if C==30:ax.axvspan(20,30,color='#d0d0d0');ax.text(25,2400,'Outside C30',rotation=90,ha='center',va='center',fontsize=6.5,color='#555555')
        ax.set(xlim=(-10,30),ylim=(4000,0),xticks=[-10,0,10,20,30],yticks=[0,1000,2000,4000],xlabel='Position from original peak',ylabel='Native time')
        ax.set_title(f'E{i}  Holmes C{C} d{distance}\n{label}',fontsize=7.7,pad=7);ax.tick_params(labelsize=6.4);ax.xaxis.label.set_size(7);ax.yaxis.label.set_size(7)
    cb=fig.colorbar(mesh,cax=fig.add_axes([.95,.355,.012,.225]));cb.ax.tick_params(labelsize=6.5);cb.set_label('Native active Rac',fontsize=7)
    fig.text(.06,.303,'1→2 changes size/material/architecture/routes. Only 2→3 moves the cue in the identical full basal state.\nCue off at 1000; histories show 0–4000, outcomes at 14,000. Gray denotes space outside C30; periodic cut keeps incumbent at 0 intact.',fontsize=7,linespacing=1.4)
    f=fig.add_axes([.10,.065,.35,.17]);g=fig.add_axes([.59,.065,.35,.17]);theta=np.linspace(0,2*np.pi,361)
    f.plot(np.cos(theta),np.sin(theta),color='#d1d1d1',lw=.7);f.axhline(0,color='#dddddd',lw=.5);f.axvline(0,color='#dddddd',lw=.5)
    for i,C,distance,col,style in ((1,30,12,ring.COLORS[1],'--'),(2,40,12,ring.COLORS[3],'-'),(3,40,16,ring.COLORS[2],'-')):
        z=vector(vectors,C,distance);keep=z['time']<=4000
        f.plot(z['vector_x'][keep],z['vector_y'][keep],color=col,ls=style,lw=1.2,label=str(i))
        for moment,marker in ((1000,'s'),(4000,'^')):
            j=np.flatnonzero(z['time']==moment)[0];f.plot(z['vector_x'][j],z['vector_y'][j],marker=marker,ms=4,color=col,mec='white',mew=.3)
        end=np.flatnonzero(z['time']==4000)[0]
        f.text(z['vector_x'][end]+(.08 if i==3 else -.12),z['vector_y'][end]+.06,str(i),color=col,fontsize=7,fontweight='bold')
        g.plot(z['time'][keep],z['relative_magnitude'][keep],color=col,ls=style,lw=1.2,label=f'{i}: cue {360*distance/C:g}°')
    f.plot(1,0,'o',mfc='white',mec='black',ms=4,zorder=8)
    f.set(xlim=(-1.15,1.15),ylim=(-1.15,1.15),xticks=[-1,0,1],yticks=[-1,0,1],xlabel='Re P',ylabel='Im P');f.set_aspect('equal',adjustable='box')
    f.set_title('F  Activity-vector space (not cell position)',loc='left',fontsize=8,pad=7);f.tick_params(labelsize=6.5)
    z=vector(vectors,40,20);keep=z['time']<=4000;g.plot(z['time'][keep],z['relative_magnitude'][keep],color='#707070',ls=':',lw=1.2,label='Cancellation: cue 180°')
    g.axvline(1000,color='#888888',ls='--',lw=.6);g.axhline(1,color='#bdbdbd',ls=':',lw=.6)
    g.set(xlim=(0,4000),ylim=(0,1.15),xticks=[0,1000,2000,4000],yticks=[0,.5,1],xlabel='Native time from cue onset',ylabel='|P| relative to own initial moment')
    g.set_title('G  Relative oriented activity',loc='left',fontsize=8,pad=7);g.tick_params(labelsize=6.5);g.xaxis.label.set_size(7);g.yaxis.label.set_size(7)
    g.legend(loc='upper right',bbox_to_anchor=(1,.28),frameon=False,fontsize=6.1,ncol=2,handlelength=1.5,columnspacing=.7)
    fig.text(.50,.013,'P = M₁(t)/M₁(0): initial (1,0); magnitude can exceed 1. Squares: cue off; triangles: t4000.\nA weak resultant can reflect cancellation: antipodal Holmes retains two sites and nearly unchanged total Rac. See linked full-vector component.',ha='center',fontsize=7,linespacing=1.4)
    save(fig,'updated_checkpoint_figure2_fates_vectors')

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    OUT.mkdir(parents=True,exist_ok=True);source(__file__);source(Path(ring.__file__))
    folds=rows(OLD/'supporting/data/all_five_incumbent_size_folds.csv');profiles=rows(OLD/'supporting/data/all_five_native_basal_profiles.csv')
    rawgain=rows(GAIN/'supporting/data/native_gain_design_evidence.csv');assert len(rawgain)==97
    gain=[]
    for r in rawgain:
        ref=next(q for q in rawgain if q['model']==r['model'] and float(q['size_native'])==REFS[r['model']]);gain.append(dict(**r,reference_nominal_size=REFS[r['model']],reference_gain_native=float(ref['cue_weighted_gain_native']),relative_nominal_size=float(r['size_native'])/REFS[r['model']],gain_over_own_reference=float(r['cue_weighted_gain_native'])/float(ref['cue_weighted_gain_native'])))
    clocks=json.loads(source(PILOT/'clock_decision.json').read_text());assert clocks['complete'] and len(clocks['models'])==5
    verification=json.loads(source(PILOT/'verification.json').read_text());assert verification['passed'] and verification['completed_cases']==15
    short=rows(PILOT/'short_curves.csv');pilot=[]
    for block in clocks['models']:
        assert block['status']=='available' and len(block['primary_results'])==3
        for r in block['primary_results']:
            assert r['passes'];tau=float(r['horizon_native']);np.testing.assert_allclose(r['primary_F'],tau*r['primary_K'],rtol=2e-14)
            q,=[q for q in short if q['model']==r['model'] and float(q['size_native'])==float(r['size_native']) and float(q['epsilon'])==.005 and np.isclose(float(q['time_native']),tau,rtol=1e-13,atol=1e-16)]
            np.testing.assert_allclose(float(q['restored_fraction']),r['primary_F'],rtol=2e-14)
            pilot.append(dict(model=r['model'],size_native=r['size_native'],horizon_native=tau,restored_fraction=r['primary_F'],conditional_K_per_native_time=r['primary_K'],empirical_E=r['primary_E']))
    summary=rows(RING/'supporting/data/expanded_ring_baseline_ratio_summary.csv');gainmap=rows(RING/'supporting/data/expanded_ring_initial_gain_map.csv');outcomes=rows(RING/'supporting/data/expanded_ring_outcome_map.csv')
    with np.load(source(RING/'supporting/data/expanded_ring_Holmes_native_display_samples.npz')) as z:display={k:z[k].copy() for k in z.files}
    vectors=rows(GLOBAL/'supporting/data/global_polarity_traces.csv')
    figure1(folds,gain,profiles,clocks);figure2_response(summary,gainmap);figure2_fates(outcomes,display,vectors)
    export('updated_checkpoint_all_five_gain.csv',gain);export('updated_checkpoint_recovery_pilot.csv',pilot)
    checks=dict(linear_gain_cases=len(gain),linear_recovery_pilot_cases=len(pilot),pilot_sizes_per_model=3,pilot_points_unjoined=True,
        recovery_assay='new tapered active-excess deactivation; no inactive hold/return coefficients',ring_common_recovery_map='not measured; explicit pending box',
        ring_planned_map_cells=len(gainmap),new_native_integrations=0,old_figures_modified=False,figure2_review_sheets=2)
    p=layout.dev_output_path(OUT,'updated_checkpoint_verification.json');p.write_text(json.dumps(checks,indent=2)+'\n');ARTIFACTS.append(p)
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root='../..',source_sha256=SOURCES,sha256={str(p.relative_to(OUT)):sha(p) for p in ARTIFACTS},checks=checks)
    p=layout.dev_output_path(OUT,'updated_checkpoint.manifest.json');p.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),artifacts=len(ARTIFACTS),checks=checks),indent=2))

if __name__=='__main__':main()
