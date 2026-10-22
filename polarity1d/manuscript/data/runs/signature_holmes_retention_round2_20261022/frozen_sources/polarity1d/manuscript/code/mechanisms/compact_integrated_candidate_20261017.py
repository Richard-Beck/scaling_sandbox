"""Compact integration of retained incumbent and native-ring evidence.

Rendering and exact data export only. No native integrations or new metrics.
The reviewed integrated figure and ring report are read-only inputs.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import shutil
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import length_dev_output_layout as layout

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
OUT = ROOT / 'figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017'
OLD = ROOT / 'figures/length_mechanisms_20261008/dev/integrated_framework_candidate_20261017'
RING = ROOT / 'figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017'
RUN = ROOT / 'data/runs/chemical_ring_functional_v2_20261017'
MODELS = ('wave_pinning', 'goryachev', 'holmes_model3')
NAMES = {'wave_pinning':'Wave-pinning', 'goryachev':'Goryachev', 'holmes_model3':'Holmes 3'}
COLORS = {'wave_pinning':'#3976a8', 'goryachev':'#c96926', 'holmes_model3':'#8560a4'}
MARKERS = {'wave_pinning':'o', 'goryachev':'s', 'holmes_model3':'^'}
SITE_COLORS = ['#4678a3', '#d77c35', '#8b65a9', '#69a7a0', '#bcb069']
SITE_LABELS = ['Single peak near old site*', 'Single peak near cue', 'Distinct old + cue peaks',
               'One region spans both sites', 'Shifted / merged, unresolved']
SOURCES, ARTIFACTS = {}, []


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def source(path):
    path = Path(path).resolve()
    SOURCES[str(path.relative_to(REPO))] = sha(path)
    return path


def read_csv(path):
    with source(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def export_csv(name, rows):
    path = layout.dev_output_path(OUT, name)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    ARTIFACTS.append(path)


def pick(rows, model, C=None, fraction=None, epsilon=None):
    return [r for r in rows if r['model']==model
            and (C is None or float(r['circumference'])==C)
            and (fraction is None or np.isclose(float(r['fraction']),fraction))
            and (epsilon is None or np.isclose(float(r['epsilon']),epsilon))]


def save(fig, name):
    for ext in ('.png','.pdf','.svg'):
        path = layout.dev_output_path(OUT, name+ext)
        fig.savefig(path, dpi=200 if ext=='.png' else None, bbox_inches='tight')
        ARTIFACTS.append(path)
    plt.close(fig)


def site_index(text):
    if text=='one old-core-associated peak; incumbent identity unresolved': return 0
    if text=='one peak associated with cue patch': return 1
    if text=='distinct peaks associated with old core and cue patch': return 2
    if text.startswith('one connected region spanning'): return 3
    if text.startswith('one shifted'): return 4
    raise ValueError(text)


def size_curves(matrix, basal):
    result=[]
    for model in MODELS:
        reference_basal=pick(basal,model,19,epsilon=.005)[0]
        reference_gain=pick(matrix,model,19,.5)[0]
        for b in sorted(pick(basal,model,epsilon=.005),key=lambda r:float(r['circumference'])):
            C=float(b['circumference'])
            m=pick(matrix,model,C,.5)[0]
            assert m['status']=='complete' and b['normalized_gain_available']=='True'
            valid=[r for r in pick(matrix,model,C) if r['status']=='complete']
            gains=np.array([float(r['initial_cue_weighted_gain']) for r in valid])
            result.append(dict(model=model,circumference=C,initial_gain_native=float(m['initial_cue_weighted_gain']),
                initial_gain_C19_reference_native=float(reference_gain['initial_cue_weighted_gain']),
                initial_gain_relative_to_C19=float(m['initial_cue_weighted_gain'])/float(reference_gain['initial_cue_weighted_gain']),
                basal_conditional_repair_native=float(b['net_repair_per_damage']),
                basal_conditional_repair_C19_reference_native=float(reference_basal['net_repair_per_damage']),
                basal_conditional_repair_relative_to_C19=float(b['net_repair_per_damage'])/float(reference_basal['net_repair_per_damage']),
                initial_gain_placement_min_native=float(gains.min()),initial_gain_placement_max_native=float(gains.max()),
                initial_gain_placement_relative_range=float(np.ptp(gains)/gains.mean()),
                ring_basal_core_amount=float(b['basal_core_amount']),ring_basal_width=float(b['basal_width']),
                ring_basal_peak=float(b['basal_peak']),checkpoint_sha256=b['checkpoint_sha256']))
    return result


def figure1(folds, curves):
    fig=plt.figure(figsize=(8.5,3.9))
    gs=fig.add_gridspec(1,3,left=.135,right=.98,bottom=.24,top=.78,wspace=.56,width_ratios=[1.15,1,1])
    axes=[fig.add_subplot(gs[0,i]) for i in range(3)]
    ax=axes[0]
    for key,marker,label,offset in [('excess_active_fold','o','Excess active amount',-.19),
            ('width_fold','s','Core width',0),('peak_concentration_fold','^','Peak concentration',.19)]:
        ax.scatter([float(r[key]) for r in folds],np.arange(5)+offset,marker=marker,s=25,
            facecolor='#333333' if marker=='o' else 'white',edgecolor='#333333',lw=.8,label=label,zorder=4)
    labels=[r['model_label']+'\n'+('L₀ 1→2' if r['model']=='spring' else f"{r['nominal_size_small']}→{r['nominal_size_large']}") for r in folds]
    ax.set(yticks=np.arange(5),yticklabels=labels,xlim=(.8,3.2),xticks=[1,2,3],ylim=(4.6,-.6),xlabel='Larger / smaller value')
    ax.tick_params(axis='y',length=0,labelsize=7.2)
    ax.axvline(1,color='#aaaaaa',ls='--',lw=.7)
    ax.set_title('A  Incumbent organization',loc='left',fontsize=9,pad=10)
    ax.legend(loc='lower left',bbox_to_anchor=(-.14,1.15),ncol=1,frameon=False,fontsize=7,handlelength=1.4,labelspacing=.25)
    titles=['B  Initial input gain','C  Basal conditional repair']
    keys=['initial_gain_relative_to_C19','basal_conditional_repair_relative_to_C19']
    for ax,key,title in zip(axes[1:],keys,titles):
        for model in MODELS:
            rr=[r for r in curves if r['model']==model]
            ax.plot([r['circumference'] for r in rr],[r[key] for r in rr],color=COLORS[model],
                marker=MARKERS[model],ms=3.1,lw=1.2,label=NAMES[model],mec='white',mew=.25)
            for C in (19,39):
                row=next(r for r in rr if r['circumference']==C)
                ax.scatter([C],[row[key]],s=33,facecolors='none',edgecolors=COLORS[model],lw=.8,zorder=5)
        ax.axhline(1,color='#aaaaaa',lw=.75,ls='--')
        ax.set(xlim=(7,83),xticks=[10,40,80],xlabel='Ring circumference C',ylabel='Relative to own C = 19 value')
        ax.set_title(title,loc='left',fontsize=9,pad=10)
        ax.grid(axis='y',color='#eeeeee',lw=.6,zorder=0)
    axes[1].set(ylim=(0,2.08),yticks=[0,.5,1,1.5,2])
    axes[2].set_yscale('log',base=2)
    axes[2].set(ylim=(.2,17),yticks=[.25,.5,1,2,4,8,16],yticklabels=['¼','½','1','2','4','8','16'])
    axes[2].set_ylabel('Relative to own C = 19 value\n(log₂ scale)')
    handles=[Line2D([],[],color=COLORS[m],marker=MARKERS[m],lw=1.2,ms=4,label=NAMES[m]) for m in MODELS]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.70,.995),ncol=3,frameon=False,fontsize=7.4,
        handlelength=1.4,columnspacing=1.2)
    axes[0].set_xlabel('Larger / smaller value\n(original interval preparations)')
    fig.text(.70,.925,'B–C: independently prepared native rings',ha='center',fontsize=7.4,color='#555555')
    fig.text(.56,.035,'B: antipodal width-3 input. C: 0.5% inactive-pool hold/return, hold 0.1.\nOwn C19 references; native rates are not calibrated across models.\nPositive conditional repair does not imply net core gain or incumbent survival.',
        ha='center',va='center',fontsize=7.1,linespacing=1.5)
    save(fig,'compact_incumbent_and_function')


def load_holmes(matrix,fraction):
    row=pick(matrix,'holmes_model3',40,fraction)[0]
    with np.load(source(row['cue_path'])) as z:
        tc=z['time'].copy();x=z['x'].copy();pc=z['states'][:,0].copy();basal=z['states'][0].copy()
    with np.load(source(row['recovery_path'])) as z:
        tr=z['time'].copy()+1000;pr=z['states'][:,0].copy()
        np.testing.assert_array_equal(z['states'][0,0],pc[-1])
    t=np.r_[tc,tr[1:]];p=np.r_[pc,pr[1:]]
    position=(x-x[0]+10)%40-10;order=np.argsort(position)
    keep=t<=4000
    out=dict(time=t[keep],position=position[order],active_rac=p[keep][:,order],
        complete_basal_state=basal,full_time=t,full_primary=p)
    path=layout.dev_output_path(OUT,f'Holmes_C40_d{float(row["short_route"]):g}_exact_native.npz')
    np.savez_compressed(path,**out);ARTIFACTS.append(path)
    return row,out


def draw_map(ax,matrix,model,sizes):
    cy=np.array(sizes)
    yedges=np.r_[cy[0]-(cy[1]-cy[0])/2,(cy[:-1]+cy[1:])/2,cy[-1]+(cy[-1]-cy[-2])/2]
    for iy,C in enumerate(sizes):
        rr=sorted(pick(matrix,model,C),key=lambda r:float(r['fraction']))
        cc=np.array([float(r.get('realized_fraction') or r['fraction']) for r in rr])
        xe=np.r_[cc[0]-.0125,(cc[:-1]+cc[1:])/2,cc[-1]+.0125]
        for j,row in enumerate(rr):
            st=row['status']
            color=SITE_COLORS[site_index(row['final_location_association'])] if st=='complete' else '#eeeeee'
            rect=Rectangle((xe[j],yedges[iy]),xe[j+1]-xe[j],yedges[iy+1]-yedges[iy],
                facecolor=color,edgecolor='white',linewidth=.2,hatch='////' if st=='unqualified_preparation' else None)
            ax.add_patch(rect)
            if st=='geometric_overlap_exclusion': ax.plot(cc[j],C,'x',color='#999999',ms=2.1,mew=.4)
    ax.set(xlim=(.0875,.5125),ylim=(yedges[0],yedges[-1]),xticks=[.1,.2,.3,.4,.5],
        yticks=[10,20,30,40,50,60,70,80],xlabel='Shorter-route fraction d/C',ylabel='Circumference C')
    if model=='holmes_model3':
        for fraction,tag in ((.3,'D'),(.4,'E')):
            row=pick(matrix,model,40,fraction)[0]
            ax.scatter(float(row['realized_fraction']),40,marker='o',s=98,facecolor='white',edgecolor='#111111',lw=.9,zorder=8)
            ax.text(float(row['realized_fraction']),40,tag,ha='center',va='center',fontsize=7.2,weight='bold',zorder=9)


def figure2(matrix):
    selected=[load_holmes(matrix,fraction) for fraction in (.3,.4)]
    np.testing.assert_array_equal(selected[0][1]['complete_basal_state'],selected[1][1]['complete_basal_state'])
    assert selected[0][0]['checkpoint_sha256']==selected[1][0]['checkpoint_sha256']
    fig=plt.figure(figsize=(8.5,7.6))
    maps=[fig.add_axes([left,.595,.245,.335]) for left in (.075,.395,.715)]
    sizes=sorted({float(r['circumference']) for r in matrix})
    for i,(ax,model) in enumerate(zip(maps,MODELS)):
        draw_map(ax,matrix,model,sizes)
        ax.set_title(f'{chr(65+i)}  {NAMES[model]}',loc='left',fontsize=9,pad=8)
    fig.text(.515,.977,'Size and placement alter the resulting activity pattern',ha='center',va='top',fontsize=10,weight='bold')
    handles=[Patch(facecolor=c,edgecolor='none',label=l) for c,l in zip(SITE_COLORS,SITE_LABELS)]
    handles += [Line2D([],[],color='#999999',marker='x',linestyle='none',ms=4,label='Cue/core overlap: excluded'),
                Patch(facecolor='#eeeeee',edgecolor='#999999',hatch='////',label='Preparation unavailable')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.555),ncol=3,frameon=False,fontsize=7,
        columnspacing=1.5,handlelength=1.4,labelspacing=.55)
    fig.text(.52,.455,'D–E  Holmes C40: same complete starting state; only cue location changes',ha='center',fontsize=8.3)
    axes=[fig.add_axes([.12,.11,.33,.285]),fig.add_axes([.56,.11,.33,.285])]
    vmax=2.4
    for ax,(row,z),tag in zip(axes,selected,('D','E')):
        mesh=ax.pcolormesh(z['position'],z['time'],z['active_rac'],shading='nearest',cmap='magma',norm=Normalize(0,vmax),rasterized=True)
        ax.axhline(1000,color='white',lw=.85,ls='--')
        ax.axvline(0,color='white',lw=.6,ls=':')
        d=float(row['short_route']);ax.axvline(d,color='#74d8ee',lw=.7,ls=':')
        ax.set(xlim=(-10,30),ylim=(4000,0),xticks=[-10,0,10,20,30],yticks=[0,1000,2000,3000,4000],
            xlabel='Ring coordinate relative to original peak',ylabel='Time from cue onset (native units)')
        result='Merger' if tag=='D' else 'Two retained sites'
        ax.set_title(f'{tag}  d = {d:g}: {result}',loc='left',fontsize=9,pad=8)
        ax.text(29,1050,'Cue off',color='white',ha='right',va='top',fontsize=7)
    cax=fig.add_axes([.91,.11,.015,.285])
    fig.colorbar(mesh,cax=cax,ticks=[0,.8,1.6,2.4],label='Active Rac (native concentration)')
    fig.text(.52,.025,'Maps: 1000 cue-on + 13,000 cue-free native units. Kymographs: first 4000 units.\n*Old-associated endpoint geography alone does not establish incumbent survival. No stable-branch classification.',
        ha='center',va='center',fontsize=7,linespacing=1.4)
    save(fig,'compact_size_location_and_histories')
    return selected


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,
        'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    OUT.mkdir(parents=True,exist_ok=True)
    source(__file__)
    for p in (RING/'RESULTS.md',RING/'supporting/metadata/RESULT_CRITIQUE.md',
              RUN/'results/audit.json',OLD/'supporting/metadata/FROZEN_HASHES.json'):
        source(p)
    audit=json.loads((RUN/'results/audit.json').read_text())
    assert audit['passed'] and audit['full_867_cell_schema_coverage']
    matrix=read_csv(RUN/'results/matrix.csv');basal=read_csv(RUN/'results/basal.csv')
    folds=read_csv(OLD/'supporting/data/all_five_incumbent_size_folds.csv')
    assert len(matrix)==867 and sum(r['status']=='complete' for r in matrix)==652
    target=layout.dev_output_path(OUT,'all_five_original_interval_folds.csv')
    shutil.copy2(OLD/'supporting/data/all_five_incumbent_size_folds.csv',target);ARTIFACTS.append(target)
    curves=size_curves(matrix,basal)
    assert len(curves)==50
    export_csv('ring_native_size_curves.csv',curves)
    export_csv('complete_endpoint_map.csv',[{k:r[k] for k in ('model','circumference','fraction','realized_fraction','short_route','long_route','status','final_location_association','checkpoint_sha256')} for r in matrix])
    figure1(folds,curves)
    selected=figure2(matrix)
    old_manifest=json.loads(source(OLD/'supporting/metadata/manifest.json').read_text())
    original_folds=OLD/'supporting/data/all_five_incumbent_size_folds.csv'
    assert sha(original_folds)==old_manifest['sha256'][str(original_folds.relative_to(OLD))]
    ring_manifest=json.loads(source(RING/'supporting/metadata/ring_functional_report.manifest.json').read_text())
    native_inputs=[RUN/'results/matrix.csv',RUN/'results/basal.csv']
    native_inputs += [Path(row[key]) for row,z in selected for key in ('cue_path','recovery_path')]
    for path in native_inputs:
        assert sha(path)==ring_manifest['source_sha256'][str(path)]
    checks=dict(planned_matrix_rows=len(matrix),complete_conditions=652,excluded_overlap=198,unavailable_preparation=17,
        size_curves=len(curves),models=list(MODELS),normalization='each native-ring model at C19',
        retained_input_hashes_match_prior_reports=True,
        maximum_onset_gain_relative_placement_range={m:max(r['initial_gain_placement_relative_range'] for r in curves if r['model']==m) for m in MODELS},
        selected_Holmes_exact_full_initial_state_equality=True,selected_matrix_mesh=.05,
        selected_Holmes_checkpoint_sha256=selected[0][0]['checkpoint_sha256'],
        display_windows=dict(endpoint_maps=[0,14000],kymographs=[0,4000]),native_integrations_performed=0,
        scientific_scope='A uses original interval preparations; functional curves and maps use separately prepared native rings')
    path=layout.dev_output_path(OUT,'verification.json');path.write_text(json.dumps(checks,indent=2)+'\n');ARTIFACTS.append(path)
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root='../..',source_sha256=SOURCES,
        artifacts={str(p.relative_to(OUT)):sha(p) for p in ARTIFACTS},checks=checks)
    path=layout.dev_output_path(OUT,'manifest.json');path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),artifacts=len(ARTIFACTS),checks=checks),indent=2))


if __name__=='__main__': main()
