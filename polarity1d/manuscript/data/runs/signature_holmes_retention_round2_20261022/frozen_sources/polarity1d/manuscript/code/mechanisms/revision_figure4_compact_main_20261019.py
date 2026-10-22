"""One-page F4 composition from the completed, guarded signature assembly.

Reuses exact plotted points/intervals and declared native examples; no estimator
or native integration is introduced. The original two-sheet package is intact.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
FIG=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4'
SOURCE=FIG/'signature_assembly_20261019'
OUT=FIG/'main_plate_20261019'
MODELS=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes 3',spring_ring='Spring ring')
LEVELS=(.001,.003,.01,.03,.1)
COLORS=('#4055a8','#0d8a93','#66923d','#d2882c','#ba445b')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--out',type=Path,default=OUT);a=p.parse_args()
    source=a.source;out=a.out
    for sub in ('data','metadata','vectors'):(out/'supporting'/sub).mkdir(parents=True,exist_ok=True)
    complete=json.loads((source/'supporting/metadata/FINAL_POSTPROCESS_RECEIPT.json').read_text())
    assert complete['status']=='complete'
    sources=[source/'supporting/data'/n for n in ('main_signature_points.csv','local_example_context.csv','local_example_profile_cells.csv','local_example_kymograph_data.npz')]
    points=pd.read_csv(sources[0]);context=pd.read_csv(sources[1]);profiles=pd.read_csv(sources[2]);fields=np.load(sources[3])
    assert len(points)==240 and (points.seeds==12).all()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'svg.fonttype':'none','pdf.fonttype':42})
    fig=plt.figure(figsize=(7.2,9.6));grid=fig.add_gridspec(4,4,height_ratios=(1.2,1.35,.25,1.35),left=.105,right=.94,bottom=.16,top=.805,hspace=.58,wspace=.4)
    fig.suptitle('Figure 4 · Input regime shapes size effects on activity memory',x=.075,y=.976,ha='left',fontsize=12.5)
    fig.text(.075,.945,'Native local organization and second-half summaries across four ring models',fontsize=9)
    handles=[Line2D([],[],color=c,marker='o',lw=1.25,ms=3,label=f'f={f:g}') for f,c in zip(LEVELS,COLORS)]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.064,.928),ncol=5,frameon=False,fontsize=8,columnspacing=1.3,handlelength=1.5)
    fig.text(.075,.885,'B–C: 12 independent seeds per condition; whiskers show 95% whole-seed intervals.',fontsize=7.8)
    fig.text(.075,.858,'A: native activity (color); largest size, f=0.1; fixed window and declared medoid.',fontsize=7.8)
    axes=[]
    for col,(model,name) in enumerate(MODELS.items()):
        row=context[context.model==model].iloc[0]
        ax=fig.add_subplot(grid[0,col]);axes.append(ax)
        t=fields[model+'_time'];q=fields[model+'_coordinate'];u=fields[model+'_primary']
        original=profiles[profiles.model==model]
        maximum=max(float(u.max()),float(original.basal.max()))
        length=1. if model=='spring_ring' else float(row['size'])
        image=ax.imshow(u.T,aspect='auto',origin='lower',interpolation='nearest',extent=[t[0],t[-1],0,length],vmin=0,vmax=maximum,cmap='viridis')
        ax.axvline(row.actual_snapshot_time,color='white',lw=.7,ls='--')
        ax.set_title(f'A{col+1}  {name}',loc='left',fontsize=8.7,pad=17)
        ax.text(0,1.025,'C₀=3 · material q' if model=='spring_ring' else 'C=60 · physical arc',transform=ax.transAxes,fontsize=7)
        ax.set_xticks([t[0],t[-1]],[f'{v:g}' for v in (t[0],t[-1])]);ax.set_xlabel('Native time',fontsize=7.4)
        ax.set_yticks([0,length]);ax.tick_params(labelsize=7)
        if col==0:ax.set_ylabel('Position',fontsize=8)
        elif model!='spring_ring':ax.tick_params(labelleft=False)
        colorbar=fig.colorbar(image,ax=ax,fraction=.055,pad=.035,ticks=[0,maximum]);colorbar.ax.tick_params(labelsize=6.1,pad=1);colorbar.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2g'))
        ax.text(.02,.04,'low contrast' if model=='goryachev' else '',transform=ax.transAxes,color='white',fontsize=7.3)
        bm=fig.add_subplot(grid[1,col]);cov=fig.add_subplot(grid[2,col]);org=fig.add_subplot(grid[3,col])
        axes.extend((bm,cov,org))
        for strength,color in zip(LEVELS,COLORS):
            for axis,metric,style in ((bm,'direction_cosine_reported','-'),(cov,'direction_coverage','-'),(org,'qualified_multiregion_fraction','-'),(org,'low_contrast_fraction','--')):
                d=points[(points.model==model)&(points.fraction==strength)&(points.metric==metric)].sort_values('size');assert len(d)==3
                x=d.basal_circumference.to_numpy();y=d.estimate.to_numpy()
                axis.plot(x,y,color=color,lw=1.05,ls=style,marker='o' if style=='-' else None,ms=2.8)
                axis.vlines(x,d.ci_low,d.ci_high,color=color,lw=.65,alpha=.62,linestyles=style)
        ticks=points[(points.model==model)&(points.fraction==.001)&(points.metric=='direction_coverage')].sort_values('size').basal_circumference.to_numpy()
        for axis in (bm,cov,org):
            axis.set_xlim(ticks[0]-.1*np.ptp(ticks),ticks[-1]+.1*np.ptp(ticks));axis.set_xticks(ticks,[f'{x:.2f}' if model=='spring_ring' else f'{x:g}' for x in ticks])
            axis.spines[['top','right']].set_visible(False);axis.tick_params(labelsize=7,pad=2);axis.grid(axis='y',color='#e7ebed',lw=.55)
        bm.set_ylim(-1.05,1.05);bm.set_yticks([-1,0,1]);bm.axhline(0,color='#9aa5aa',lw=.55)
        bm.set_title(f'B{col+1}  lag {800 if model=="spring_ring" else 8000}',loc='left',fontsize=8.6,pad=7)
        bm.tick_params(labelbottom=False)
        cov.set_ylim(-.04,1.04);cov.set_yticks([0,1]);cov.tick_params(labelbottom=False,labelsize=6.7)
        org.set_ylim(-.035,1.04);org.set_yticks([0,.5,1]);org.set_title(f'C{col+1}',loc='left',fontsize=8.6,pad=7)
        org.set_xlabel('Basal physical C' if model=='spring_ring' else 'Circumference',fontsize=7.5)
        if col==0:
            bm.set_ylabel('Heading cosine',fontsize=8);cov.set_ylabel('Usable\npairs',fontsize=7);org.set_ylabel('Fraction of sampled times',fontsize=8)
        else:
            for axis in (bm,cov,org):axis.tick_params(labelleft=False)
    styles=[Line2D([],[],color='#40515c',lw=1.2,label='Contrast-qualified multiple regions'),Line2D([],[],color='#40515c',lw=1.2,ls='--',label='Low contrast')]
    fig.legend(handles=styles,loc='lower left',bbox_to_anchor=(.065,.102),ncol=2,frameon=False,fontsize=7.7,columnspacing=1.4)
    fig.text(.075,.085,'Direction requires both endpoint R≥0.1; the thin strip retains all-pair coverage.',fontsize=7.5)
    fig.text(.075,.064,'Contrast is qualified at ≥10% of basal. Full profiles, matched mean controls and windows: supplement.',fontsize=7.5)
    fig.text(.075,.043,'Native clocks/laws differ. Whole positive input scales with f. Activity heading is distinct from migration.',fontsize=7.5)
    fig.text(.075,.022,'First half discarded without assuming stationarity. Spring sampled run-duration inference remains excluded.',fontsize=7.5)
    target=out/'figure4_main.png';fig.savefig(target,dpi=220)
    for extension in ('pdf','svg'):fig.savefig(out/'supporting/vectors'/('figure4_main.'+extension))
    plt.close(fig)
    points.to_csv(out/'supporting/data/main_signature_points.csv',index=False);context.to_csv(out/'supporting/data/native_example_context.csv',index=False)
    receipt=dict(generated_at=datetime.now(timezone.utc).isoformat(),source=str(Path(__file__).relative_to(ROOT)),source_sha256=sha(__file__),native_integrations=0,
        estimators_changed=False,source_assembly_receipt_sha256=sha(source/'supporting/metadata/FINAL_POSTPROCESS_RECEIPT.json'),
        source_data_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},layout_inches=[7.2,9.6],main_points=len(points),
        sources_preserved=True,independent_render_review='pending',outputs={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.suffix in ('.png','.pdf','.svg','.csv')})
    (out/'supporting/metadata/COMPOSITION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
