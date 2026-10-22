"""Render one checkpoint with a strict linear-then-ring geometry split.

Figure 1 uses the original all-five interval states and interval functional
assays. Figure 2 reuses the reviewed native-ring spatial comparison unchanged.
No native integration, new metric, or change to previous candidates.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import os
import shutil
import sys

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import compact_integrated_candidate_20261017 as base

OUT=base.ROOT/'figures/length_mechanisms_20261008/dev/linear_then_ring_checkpoint_20261017'
PREVIOUS=base.ROOT/'figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017'
FUNCTION_MODELS=('wave_pinning','goryachev')


def copy_file(path,name=None):
    path=base.source(path)
    target=base.layout.dev_output_path(OUT,name or path.name)
    shutil.copy2(path,target);base.ARTIFACTS.append(target)
    return target


def figure1(folds,functions,profiles):
    fig=plt.figure(figsize=(8.5,6.1))
    axes=[fig.add_axes(bounds) for bounds in ((.135,.55,.225,.34),(.475,.55,.195,.34),(.785,.55,.195,.34))]
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
    ax.legend(loc='upper left',bbox_to_anchor=(-.14,1.40),frameon=False,fontsize=7,handlelength=1.4,labelspacing=.25)
    for ax,key,title in zip(axes[1:],('cue_gain_relative_to_L19','net_repair_relative_to_L19'),
            ('B  Opposite-end input gain','C  Conditional core repair')):
        for model in FUNCTION_MODELS:
            rr=sorted([r for r in functions if r['model']==model],key=lambda r:float(r['length_native']))
            ax.plot([float(r['length_native']) for r in rr],[float(r[key]) for r in rr],
                color=base.COLORS[model],marker=base.MARKERS[model],ms=3,lw=1.2,mec='white',mew=.25)
            for L in (19,39):
                r=next(r for r in rr if float(r['length_native'])==L)
                ax.scatter([L],[float(r[key])],s=31,facecolors='none',edgecolors=base.COLORS[model],lw=.8,zorder=5)
        ax.axhline(1,color='#aaaaaa',lw=.75,ls='--')
        ax.set(xlim=(8,41),xticks=[9,19,29,39],xlabel='Interval length L',ylabel='Relative to own L = 19 value')
        ax.set_title(title,loc='left',fontsize=9,pad=10)
        ax.grid(axis='y',color='#eeeeee',lw=.6,zorder=0)
    axes[1].set(ylim=(0,2.08),yticks=[0,.5,1,1.5,2])
    axes[2].set(ylim=(0,4.2),yticks=[0,1,2,3,4])
    axes[2].text(.04,.95,'G conversion: 3.56× (19→39)',transform=axes[2].transAxes,
        ha='left',va='top',fontsize=6.5,color=base.COLORS['goryachev'])
    handles=[Line2D([],[],color=base.COLORS[m],marker=base.MARKERS[m],lw=1.2,ms=4,label=base.NAMES[m]) for m in FUNCTION_MODELS]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.70,.995),ncol=2,frameon=False,fontsize=7.5,
        handlelength=1.5,columnspacing=1.5)
    fig.text(.70,.94,'B–C: tall versus broad interval exemplars',ha='center',fontsize=7.4,color='#555555')
    fig.text(.55,.462,'All panels: retained linear 1D interval preparations; no ring-derived measurements.',ha='center',fontsize=7.1)
    fig.text(.55,.438,'C: 0.5% whole inactive profile withheld for 0.1, then returned exactly.\nConditional response to this disturbance; not an all-five stability measurement.',ha='center',va='center',fontsize=7,linespacing=1.4)
    fig.text(.075,.366,'D  Native basal profiles: the same small/large pairs as A',ha='left',fontsize=8.8)
    fig.legend(handles=[Line2D([],[],color='#969696',ls='--',lw=1.2,label='Smaller'),
            Line2D([],[],color='#222222',lw=1.4,label='Larger')],loc='center',bbox_to_anchor=(.79,.37),
        ncol=2,frameon=False,fontsize=7.2,handlelength=2,columnspacing=1.3)
    for left,fold in zip((.075,.265,.455,.645,.835),folds):
        ax=fig.add_axes([left,.135,.135,.165])
        model=fold['model'];group=[r for r in profiles if r['model']==model]
        for size,color,style in ((float(fold['nominal_size_small']),'#969696','--'),
                                (float(fold['nominal_size_large']),'#222222','-')):
            rr=sorted([r for r in group if float(r['nominal_size'])==size],key=lambda r:float(r['physical_position']))
            ax.plot([float(r['physical_position']) for r in rr],[float(r['active_concentration']) for r in rr],color=color,ls=style,lw=1.25)
        xmax=float(fold['physical_size_large'])
        ymax=max(float(r['active_concentration']) for r in group)
        ax.set(xlim=(0,xmax),ylim=(0,ymax*1.10),title=fold['model_label'])
        ax.set_xticks([0,xmax]);ax.set_xticklabels(['0',f'{xmax:.2f}' if model=='spring' else f'{xmax:g}'])
        ax.tick_params(labelsize=6.7)
        ax.set_title(fold['model_label'],fontsize=7.5,pad=5)
        if model=='wave_pinning':ax.set_ylabel('Active concentration\n(native units)',fontsize=7)
    fig.text(.52,.065,'Physical interval coordinate (native units); concentration scales differ between models.\nSpring uses its measured domains 1.049→2.049, not resting lengths as spatial coordinates.',
        ha='center',va='center',fontsize=7,linespacing=1.4)
    base.save(fig,'linear_incumbent_organization_and_mechanisms')


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,
        'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    OUT.mkdir(parents=True,exist_ok=True);base.OUT=OUT
    base.source(__file__);base.source(base.__file__)
    prior=json.loads(base.source(PREVIOUS/'supporting/metadata/manifest.json').read_text())
    linear=json.loads(base.source(base.OLD/'supporting/metadata/manifest.json').read_text())
    protected={str(PREVIOUS/rel):digest for rel,digest in prior['artifacts'].items()}
    pdf=base.ROOT/'figures/length_mechanisms_20261008/length_mechanisms.pdf'
    pdf_hash=base.sha(base.source(pdf));protected[str(pdf)]=pdf_hash
    original_verification=json.loads(base.source(base.OLD/'supporting/metadata/verification.json').read_text())
    assert pdf_hash==original_verification['pdf_sha256']
    linear_names=('all_five_incumbent_size_folds.csv','conditional_functional_readouts.csv','all_five_native_basal_profiles.csv')
    data=[]
    for name in linear_names:
        path=base.OLD/'supporting/data'/name
        assert base.sha(path)==linear['sha256'][str(path.relative_to(base.OLD))]
        data.append(base.read_csv(path));copy_file(path)
    folds,functions,profiles=data
    assert len(folds)==5 and len(functions)==36
    assert set(r['model'] for r in functions)==set(FUNCTION_MODELS)
    assert set(r['model'] for r in profiles)==set(r['model'] for r in folds)
    for r in functions:
        np.testing.assert_allclose(float(r['cue_gain_relative_to_L19']),float(r['initial_local_cue_gain_native'])/float(r['cue_gain_L19_reference_native']),rtol=2e-15)
        np.testing.assert_allclose(float(r['net_repair_relative_to_L19']),float(r['immediate_net_core_repair_per_damage_native'])/float(r['net_repair_L19_reference_native']),rtol=2e-15)
    figure1(folds,functions,profiles)
    # The spatial comparison is already reviewed. Preserve its actual pixels
    # and exact retained arrays rather than silently altering the second figure.
    spatial_files=['compact_size_location_and_histories.png',
        'supporting/vectors/compact_size_location_and_histories.pdf',
        'supporting/vectors/compact_size_location_and_histories.svg',
        'supporting/data/complete_endpoint_map.csv',
        'supporting/data/Holmes_C40_d12_exact_native.npz','supporting/data/Holmes_C40_d16_exact_native.npz']
    for rel in spatial_files:
        path=PREVIOUS/rel
        assert base.sha(path)==prior['artifacts'][rel]
        copy_file(path)
    with np.load(OUT/'supporting/data/Holmes_C40_d12_exact_native.npz') as a, np.load(OUT/'supporting/data/Holmes_C40_d16_exact_native.npz') as b:
        np.testing.assert_array_equal(a['complete_basal_state'],b['complete_basal_state'])
    for path,digest in protected.items():assert base.sha(path)==digest,path
    endpoints={model:{key:float(next(r for r in functions if r['model']==model and float(r['length_native'])==39)[key])
        for key in ('cue_gain_relative_to_L19','net_repair_relative_to_L19')} for model in FUNCTION_MODELS}
    checks=dict(figure1_geometry='linear 1D intervals only',figure1_models=[r['model'] for r in folds],
        linear_functional_exemplars=list(FUNCTION_MODELS),linear_functional_cases=len(functions),linear_L19_to_L39=endpoints,
        all_linear_input_hashes_match_reviewed_integrated_manifest=True,figure2_geometry='native periodic chemical rings only',
        figure2_pixels_identical_to_reviewed_spatial_comparison=True,Holmes_complete_basal_equality=True,
        original_pdf_sha256=pdf_hash,previous_compact_artifacts_unchanged=True,
        native_integrations_performed=0,new_scientific_metrics=0,debelly_only_in_linear_organization_and_profiles=True)
    path=base.layout.dev_output_path(OUT,'verification.json');path.write_text(json.dumps(checks,indent=2)+'\n');base.ARTIFACTS.append(path)
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root='../..',source_sha256=base.SOURCES,
        artifacts={str(p.relative_to(OUT)):base.sha(p) for p in base.ARTIFACTS},checks=checks)
    path=base.layout.dev_output_path(OUT,'manifest.json');path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),artifacts=len(base.ARTIFACTS),checks=checks),indent=2))


if __name__=='__main__':main()
