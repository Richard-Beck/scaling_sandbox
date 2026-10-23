"""Slim Figure4 with honest partially identified birth/takeover observables."""
from __future__ import annotations
import csv, hashlib, json, math, shutil
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import SymLogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import numpy as np
import figure4_round3_render_20261022 as R
import figure4_round5_analysis_20261023 as A

OUT=A.OUT;OLD=A.OLD;DATA=OUT/'supporting/data';META=OUT/'supporting/metadata'
MODELS=R.MODELS;NAMES=R.NAMES
COMP_LABELS=('Motion','Amount','Membership','Residual')
COMP_SHORT=('M','A','B','R')
read=R.read_csv;n=R.n

def copy_baseline():
    DATA.mkdir(parents=True,exist_ok=True);META.mkdir(parents=True,exist_ok=True)
    # Hard links preserve bytes without duplicating retained compact evidence.
    for src in (OLD/'supporting/data').iterdir():
        if src.is_file() and not src.name.startswith('plotted_'):
            dest=DATA/src.name
            if not dest.exists():dest.hardlink_to(src)

def main_axes(ax,model=None,xlabel=True):
    R.style_axis(ax)
    if model:R.facet(ax,NAMES[model])
    if xlabel:R.log_f(ax)

def section(fig,spec,label,text=''):
    pos=spec.get_position(fig)
    fig.text(.018,pos.y1+.012,label,ha='left',va='bottom',fontweight='bold',fontsize=10,color='black')
    if text:fig.text(.082,pos.y1+.013,text,fontsize=6.5,ha='left',va='bottom')

def nulls(fig,spec):
    examples=R.load_kymographs(DATA,'null_A');acf=read(DATA/'null_A_acf.csv')
    grid=spec.subgridspec(1,3,wspace=.33)
    labels=('Drifting peak','Same-rate jumps','Perimeter-rate jumps')
    for j,kind in enumerate(R.NULLS):
        sub=grid[j].subgridspec(2,1,height_ratios=[1.3,1],hspace=.32)
        pair=sub[0].subgridspec(1,2,width_ratios=[1,3],wspace=.22)
        for k,size in enumerate((20,60)):
            record=next(r for r in examples if R.null_kind(r)==kind and R.equivalent_size(r)==size)
            ax=fig.add_subplot(pair[k]);R.kymo(ax,record,1.,ylabel=k==0,label=f'C{size}')
        ax=fig.add_subplot(sub[1]);R.correlation(ax,acf,kind,null=True,legend=j==0)
        if j==0:ax.set_ylabel(r'$C(\ell)$',labelpad=2)
        ax.text(.96,.88,r'$\alpha='+('2','0','-1')[j]+'$',transform=ax.transAxes,ha='right',fontsize=6)
        pos=grid[j].get_position(fig);fig.text((pos.x0+pos.x1)/2,pos.y1+.004,labels[j],fontsize=6.6,ha='center')

def alpha(fig,spec):
    grid=spec.subgridspec(1,6,width_ratios=[1,1,1,1,.15,1.18],wspace=.43)
    rows=R.primary(read(DATA/'local_alpha.csv'));folds=R.primary(read(DATA/'noise_tolerance_folds.csv'))
    for j,model in enumerate(MODELS):
        ax=fig.add_subplot(grid[j]);R.alpha_panel(ax,rows,model,ylabel=j==0)
        if j>0:ax.tick_params(labelleft=False)
        bounds=grid[j].get_position(fig)
        R.ramp_key(fig,[model],bounds.y1+.012,x=bounds.x0+.1*bounds.width,width=.8*bounds.width)
    ax=fig.add_subplot(grid[5]);R.fold_panel(ax,folds)
    ax.set_ylabel(r'$f_{1/2}/f_{1/2}(C_{\min})$',labelpad=2)
    ax.set_ylim(.95,5.3);ax.set_yticks([1,2,4],['1','2','4'])
    ax.set_xlabel('C (Spring ×20)',fontsize=5.6)

def centred_examples(fig,spec):
    records=R.load_kymographs(DATA,'model_D');selection=read(DATA/'model_D_selection.csv')
    acf=R.primary(read(DATA/'condition_lag_metrics.csv'));centering=[]
    grid=spec.subgridspec(1,4,wspace=.39)
    for j,model in enumerate(MODELS):
        sub=grid[j].subgridspec(2,1,height_ratios=[1.45,1],hspace=.40)
        selected=[r for r in records if r['model']==model];chosen=float(next(r['f'] for r in selection if r['model']==model))
        pair=sub[0].subgridspec(1,2,width_ratios=[1,3],wspace=.28)
        commonmax=max(r['field'].max() for r in selected)
        for k,size in enumerate((20,60)):
            row=next(r for r in selected if R.equivalent_size(r)==size);field=row['field'];nb=field.shape[1]
            theta=2*np.pi*(np.arange(nb)+.5)/nb
            phase=float(np.angle(field[0]@np.exp(1j*theta))%(2*np.pi))
            shift=int(np.rint(nb*(.5-phase/(2*np.pi))))
            centred=dict(row,field=np.roll(field,shift,axis=1))
            vmax=(9.2 if k==0 else 83.) if model=='goryachev' else commonmax
            ax=fig.add_subplot(pair[k]);native=size/20 if model=='spring_ring' else size
            switches=json.loads(row.get('switch_times','[]'))
            R.kymo(ax,centred,vmax,ylabel=k==0,label=f'C{native:g}',switches=switches,elapsed=True)
            if model=='goryachev':ax.text(.5,1.02,f'0–{vmax:g}',ha='center',transform=ax.transAxes,fontsize=5.1)
            centering.append(dict(model=model,C_label=native,f=chosen,initial_heading=phase,roll_bins=shift,
                residual_centering_radians=float(A.wrapped(phase+shift*2*np.pi/nb-np.pi)),activity_max=vmax,
                activity_scale='per_size' if model=='goryachev' else 'shared_within_model'))
        ax=fig.add_subplot(sub[1]);rr=[r for r in acf if r['model']==model and math.isclose(n(r,'f'),chosen)]
        R.correlation(ax,rr,model,legend=j==0)
        if j==0:ax.set_ylabel(r'$C(\ell)$',labelpad=2)
        pos=grid[j].get_position(fig);fig.text((pos.x0+pos.x1)/2,pos.y1+.005,NAMES[model]+f', f = {chosen:g}',ha='center',fontsize=6.2)
    A.write_csv(DATA/'kymograph_centering.csv',centering)

def drift(fig,spec):
    rows=[r for r in R.primary(read(DATA/'drift_decomposition.csv')) if n(r,'flicker_frames')==5 and n(r,'MSD_range_fraction')==.25]
    grid=spec.subgridspec(2,4,wspace=.35,hspace=.38);plotrows=[];fits=[]
    for j,model in enumerate(MODELS):
        for i,key in enumerate(('Ds_over_smallest','kpred_over_kobs')):
            ax=fig.add_subplot(grid[i,j])
            ax.set_yscale('log')
            if i==0:
                ax.set_ylim((1e-7,2.) if model=='goryachev' else (.2,8.))
                ax.set_yticks([1e-6,1e-3,1] if model=='goryachev' else [.25,1,4],
                    ['$10^{-6}$','$10^{-3}$','1'] if model=='goryachev' else ['0.25','1','4'])
            else:
                ax.set_ylim((.003,3.3) if model=='goryachev' else (.3,3.3))
                ax.set_yticks([.01,.1,1] if model=='goryachev' else [.5,1,2],
                    ['0.01','0.1','1'] if model=='goryachev' else ['0.5','1','2'])
            ax.yaxis.set_minor_formatter(NullFormatter())
            for size in sorted({n(r,'C_label') for r in rows if r['model']==model}):
                rr=sorted([r for r in rows if r['model']==model and n(r,'C_label')==size],key=lambda r:n(r,'f'))
                if i==1:rr=[r for r in rr if n(r,'multiple_front_fraction')<=.05 and n(r,'f')<=.01]
                color=R.size_color(model,size*20 if model=='spring_ring' else size)
                ax.plot([n(r,'f') for r in rr],[n(r,key) for r in rr],color=color,lw=.65)
                for r in rr:
                    value=n(r,key)
                    if np.isfinite(value) and value>0:
                        R.bounded_point(ax,n(r,'f'),value,n(r,key+'_ci_lower'),n(r,key+'_ci_upper'),color,
                            open_marker=R.yes(r.get('forced_combined')),ms=1.55,lw=.35)
                        plotrows.append(dict(r,quantity=key,value=value,cutoff='f<=.01 and multi<=.05' if i else 'all inputs'))
                if i==0 and model!='goryachev':
                    weak=[r for r in rr if n(r,'f')<=.01 and n(r,'Ds')>0]
                    exponent=np.polyfit(np.log([n(r,'f') for r in weak]),np.log([n(r,'Ds') for r in weak]),1)[0]
                    fits.append(dict(model=model,C_label=size,exponent=exponent,n_points=len(weak),input_max=.01))
            ax.axhline(1,color='.65',ls='--',lw=.5)
            R.log_f(ax)
            if i==1:ax.set_xlim(8e-5,.013);ax.set_xticks([1e-4,1e-3,.01],['$10^{-4}$','$10^{-3}$','$10^{-2}$'])
            if j==0:ax.set_ylabel((r'$D_s/D_s(C_{\min})$' if i==0 else r'$k_{\rm pred}/k_{\rm obs}$'),labelpad=2)
            elif model!='goryachev':ax.tick_params(labelleft=False)
            if i==0:R.facet(ax,NAMES[model])
            if i==1 and not any(r['model']==model and r['quantity']==key for r in plotrows):
                ax.text(.5,.5,'No eligible\nclosure',transform=ax.transAxes,ha='center',va='center',fontsize=6)
            R.style_axis(ax)
    A.write_csv(DATA/'plotted_D.csv',plotrows);A.write_csv(DATA/'weak_input_drift_exponents.csv',fits)

def birth_panels(fig,spec):
    rows=read(DATA/'birth_heading_condition_metrics.csv');takeovers=read(DATA/'takeover_by_size_distance.csv')
    grid=spec.subgridspec(3,4,height_ratios=[1,1,2.35],wspace=.36,hspace=.38)
    symbols=('o','s','^');linestyles=('-','--',':');labels=(r'$d\leq0$',r'$0<d\leq C/6$',r'$d>C/6$')
    maps=[]
    for j,model in enumerate(MODELS):
        sizes=sorted({n(r,'C_label') for r in rows if r['model']==model})
        ax=fig.add_subplot(grid[0,j]);R.facet(ax,NAMES[model])
        for size in sizes:
            rr=sorted([r for r in rows if r['model']==model and n(r,'C_label')==size],key=lambda r:n(r,'f'))
            color=R.size_color(model,size*20 if model=='spring_ring' else size)
            ax.plot([n(r,'f') for r in rr],[n(r,'certified_birth_rate_lower_bound') for r in rr],'-o',ms=1.6,color=color,lw=.6)
        ax.set_yscale('symlog',linthresh=1e-8,linscale=.6);ax.set_ylim(0,.1)
        ax.set_yticks([0,1e-7,1e-5,1e-3],['0','$10^{-7}$','$10^{-5}$','$10^{-3}$'])
        R.log_f(ax);R.style_axis(ax)
        if j==0:ax.set_ylabel('i   Birth rate\n(certified lower bound)',labelpad=2,fontsize=6)
        else:ax.tick_params(labelleft=False)
        ax=fig.add_subplot(grid[1,j])
        for b in range(3):
            rr=sorted([r for r in takeovers if r['model']==model and n(r,'distance_bin')==b],key=lambda r:n(r,'C_label'))
            x=[n(r,'C_label') for r in rr];y=[n(r,'probability') for r in rr]
            ax.plot(x,y,linestyles[b],color='.25',marker=symbols[b],ms=2.2,lw=.6,mfc='white',mew=.6,label=labels[b])
            for r in rr:
                if np.isfinite(n(r,'probability')):
                    ax.plot([n(r,'C_label')]*2,[n(r,'ci_lower'),n(r,'ci_upper')],color='.4',lw=.4,zorder=0)
        ax.set_ylim(-.03,1.03);ax.set_xticks([sizes[0],sizes[2],sizes[-1]])
        ax.set_xlabel('Circumference C',labelpad=1.5)
        if j==0:
            ax.set_ylabel('ii   Takeover within 20τ\n(certified subset)',labelpad=2,fontsize=6)
            ax.legend(frameon=False,fontsize=4.9,handlelength=1.1,labelspacing=.1,loc='upper left')
        else:ax.tick_params(labelleft=False)
        R.style_axis(ax)
        ax=fig.add_subplot(grid[2,j]);fs=sorted({n(r,'f') for r in rows if r['model']==model})
        matrix=np.full((len(sizes)*5-1,len(fs)),np.nan)
        for si,size in enumerate(sizes):
            for fi,f in enumerate(fs):
                r=next(r for r in rows if r['model']==model and n(r,'C_label')==size and n(r,'f')==f)
                for ci,c in enumerate(A.COMPONENTS):matrix[si*5+ci,fi]=n(r,c+'_fraction')
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#eeeeee')
        # Explicit clipping preserves source signed fractions (which may exceed
        # one); black dot denotes values beyond this displayed scale.
        im=ax.imshow(matrix,aspect='auto',cmap=cmap,norm=SymLogNorm(linthresh=.1,vmin=-10,vmax=10),interpolation='nearest',origin='upper')
        for yi,xi in zip(*np.where(np.isfinite(matrix)&(abs(matrix)>10))):ax.plot(xi,yi,'.',color='black',ms=1.6)
        for si,size in enumerate(sizes):
            if si:ax.axhline(si*5-.9,color='white',lw=2)
            ax.text(-.085,1-(si*5+1.5)/matrix.shape[0],f'{size:g}',ha='right',va='center',transform=ax.transAxes,fontsize=5.5)
        ax.set_yticks([si*5+ci for si in range(5) for ci in range(4)],list(COMP_SHORT)*5,fontsize=5.1)
        selected=[fi for fi,f in enumerate(fs) if f in (.0001,.001,.01,.1,1)]
        ax.set_xticks(selected,[f'{fs[fi]:g}' for fi in selected],rotation=50,ha='right',fontsize=5.1)
        ax.set_xlabel('Input fraction f',labelpad=1.5)
        ax.tick_params(length=0,pad=1.3)
        for spine in ax.spines.values():spine.set_visible(False)
        if j==0:ax.set_ylabel('iii   Signed heading-loss fraction       C',labelpad=16,fontsize=6)
        maps.append(im)
    pos=grid[2,3].get_position(fig)
    cax=fig.add_axes([pos.x1+.012,pos.y0,.009,pos.height]);bar=fig.colorbar(maps[-1],cax=cax,ticks=[-10,-1,-.1,0,.1,1,10],extend='both')
    bar.set_ticklabels(['−10','−1','−0.1','0','0.1','1','10'])
    bar.ax.tick_params(labelsize=5,length=2,pad=1);bar.outline.set_linewidth(.35)
    bounds=grid[2,0].get_position(fig)
    fig.text(.082,bounds.y1+.005,'M motion   A amount   B membership   R residual; gray unavailable; • |fraction| > 10',fontsize=5.5)
    bounds=grid[1,0].get_position(fig)
    fig.text(.082,bounds.y1+.005,r'Pooled $f > 0.01$; $\tau=50$ (Spring 5); $d=\mathrm{arc}-A_{\rm inc}/(2h)$, $h$ = global contrast',fontsize=5.5)

def supplement():
    temp=OUT/'figure4_round5_new_supplement.pdf'
    with PdfPages(temp) as pdf:
        def save(fig,title,caption=''):
            fig.text(.06,.97,title,ha='left',va='top',fontsize=11)
            if caption:fig.text(.06,.025,caption,ha='left',va='bottom',fontsize=6.7,wrap=True)
            pdf.savefig(fig);plt.close(fig)
        fig=plt.figure(figsize=(8.27,7));grid=fig.add_gridspec(1,4,left=.075,right=.97,bottom=.16,top=.85,wspace=.48)
        for j,kind in enumerate(R.NULLS):R.alpha_panel(fig.add_subplot(grid[j]),read(DATA/'null_B_local_alpha.csv'),kind,null=True,ylabel=j==0)
        R.fold_panel(fig.add_subplot(grid[3]),read(DATA/'null_B_noise_tolerance_folds.csv'),null=True)
        save(fig,'S4.1   Displaced null-mixture exponents and tolerance folds','Historical simulated nulls; unchanged source records. These diagnostic curves were removed from the main figure.')
        rows=[r for r in R.primary(read(DATA/'drift_decomposition.csv')) if n(r,'flicker_frames')==5 and n(r,'MSD_range_fraction')==.25]
        free=read(DATA/'detector_free_decomposition.csv')
        fig,axes=plt.subplots(3,4,figsize=(8.27,9.2));fig.subplots_adjust(left=.085,right=.97,bottom=.11,top=.91,hspace=.48,wspace=.35)
        for j,model in enumerate(MODELS):
            for i,key in enumerate(('lambda_eff','kpred_over_kobs','multiple_front_fraction')):
                ax=axes[i,j]
                for size in sorted({n(r,'C_label') for r in rows if r['model']==model}):
                    rr=sorted([r for r in rows if r['model']==model and n(r,'C_label')==size],key=lambda r:n(r,'f'))
                    color=R.size_color(model,size*20 if model=='spring_ring' else size)
                    ax.plot([n(r,'f') for r in rr],[n(r,key) for r in rr],color=color,lw=.7)
                    if i<2:
                        rr=sorted([r for r in free if r['model']==model and n(r,'C_label')==size],key=lambda r:n(r,'f'))
                        ax.plot([n(r,'f') for r in rr],[n(r,key) for r in rr],color=color,lw=.65,ls='--')
                R.log_f(ax);R.style_axis(ax)
                if i==0:R.facet(ax,NAMES[model]);ax.set_yscale('symlog',linthresh=1e-8)
                if i==1:ax.set_yscale('log');ax.axhline(1,color='.6',ls=':',lw=.5)
                if i==2:ax.set_ylim(0,1)
                if j==0:ax.set_ylabel((r'$\lambda_{\rm eff}$',r'$k_{\rm pred}/k_{\rm obs}$','Multiple-front fraction')[i])
        save(fig,'S4.2   Displaced drift/relocation diagnostics','Solid: tracked; dashed: detector-free fit. All inputs and failures retained, with no rare-multiple filter. Native clocks differ across models.')
        birth=read(DATA/'birth_heading_condition_metrics.csv')
        fig,axes=plt.subplots(3,4,figsize=(8.27,9.2));fig.subplots_adjust(left=.09,right=.97,bottom=.11,top=.91,hspace=.5,wspace=.35)
        for j,model in enumerate(MODELS):
            for size in sorted({n(r,'C_label') for r in birth if r['model']==model}):
                rr=sorted([r for r in birth if r['model']==model and n(r,'C_label')==size],key=lambda r:n(r,'f'))
                color=R.size_color(model,size*20 if model=='spring_ring' else size)
                for i,key in enumerate(('certification_fraction','pair_coverage','persistent_onset_rate_upper_bound')):
                    axes[i,j].plot([n(r,'f') for r in rr],[n(r,key) for r in rr],color=color,lw=.7,marker='o',ms=2)
            for i in range(3):
                ax=axes[i,j];R.log_f(ax);R.style_axis(ax)
                if i<2:ax.set_ylim(-.03,1.03)
                else:ax.set_yscale('symlog',linthresh=1e-8)
                if i==0:R.facet(ax,NAMES[model])
                if j==0:ax.set_ylabel(('Certified/onsets','Lag-pair coverage','All secondary onset rate')[i])
        save(fig,'S4.3   Identification and endpoint coverage','Certified births are lower bounds. All persistent secondary onsets give upper bounds including splits. Different certification fractions prevent population-level comparisons of takeover probability.')
        fig=plt.figure(figsize=(8.27,6.3));ax=fig.add_axes([.075,.17,.86,.67]);ax.axis('off')
        text=('Requested split counts are unavailable for every one of the 260 conditions.\n\n'
              'The compact producer stored maximum-overlap identities, excess-weighted contour centroids,\n'
              'full-basin amounts and contour-cell widths. It discarded contour masks/endpoints, parent IDs,\n'
              'per-region plateau heights and event lists. A centroid and width do not determine the component\n'
              'interval, because the centroid is weighted by the unknown within-component profile.\n\n'
              'A certified no-parent birth is a new identity that persists for at least five frames and whose\n'
              'conservative angular support enclosure misses every previous enclosure. This is a subset of\n'
              'all no-overlap births. No uncertified track is called a split or a true birth.\n\n'
              'Incumbent half-width is estimated as full-basin amount / (2 × current global contrast),\n'
              'using the prior frame. It is an explicit plateau approximation, not a recovered local plateau.\n'
              'Takeovers require five consecutive raw dominant-ID frames, within a common follow-up of\n'
              '20 noise blocks (1,000 chemical units; 100 Spring). Shorter follow-up is excluded. Distances\n'
              'are grouped as d ≤ 0, 0 < d ≤ C/6 and d > C/6; designed inputs f > 0.01 are pooled.\n\n'
              'The historical aggregate region audits used a two-frame pairwise-dominance rule and include\n'
              'split onsets. They cannot replace the new five-frame event definition. No noise runs were added.')
        ax.text(0,1,text,va='top',fontsize=9,linespacing=1.45)
        save(fig,'S4.4   Missing information and exact scope of the birth observable')
        for model in MODELS:
            sizes=sorted({n(r,'C_label') for r in birth if r['model']==model})
            fig,axes=plt.subplots(5,1,figsize=(8.27,10.5));fig.subplots_adjust(left=.11,right=.95,bottom=.08,top=.91,hspace=.48)
            for ax,size in zip(axes,sizes):
                rr=sorted([r for r in birth if r['model']==model and n(r,'C_label')==size],key=lambda r:n(r,'f'))
                for ci,c in enumerate(A.COMPONENTS):
                    xx=np.array([n(r,'f') for r in rr]);yy=np.array([n(r,c+'_fraction') for r in rr])
                    ax.plot(xx,yy,label=COMP_LABELS[ci],lw=.8,color=plt.get_cmap('tab10')(ci))
                    ax.fill_between(xx,[n(r,c+'_fraction_ci_lower') for r in rr],[n(r,c+'_fraction_ci_upper') for r in rr],color=plt.get_cmap('tab10')(ci),alpha=.12)
                ax.axhline(0,color='.7',lw=.4);ax.axhline(1,color='.7',lw=.4);R.log_f(ax);R.style_axis(ax)
                ax.set_ylabel(f'C{size:g}\nSigned fraction');ax.set_ylim(-2,3)
            axes[0].legend(frameon=False,ncol=4,fontsize=7,loc='upper left')
            save(fig,'S4.5   '+NAMES[model]+': signed heading-loss fractions','Pointwise whole-seed bootstrap intervals (2,000 draws). Values outside displayed limits remain in CSV. Fractions sum to one; negative values represent cancellation. This allocation is descriptive and convention dependent.')
    # Every preceding historical supplementary page is appended unchanged.
    from pypdf import PdfReader, PdfWriter
    writer=PdfWriter()
    newreader=PdfReader(temp);oldreader=PdfReader(OLD/'figure4_supplement.pdf')
    for source in (newreader,oldreader):
        for page in source.pages:writer.add_page(page)
    with (OUT/'figure4_supplement.pdf').open('wb') as stream:writer.write(stream)
    META.joinpath('supplement_preservation.json').write_text(json.dumps(dict(new_pages=len(newreader.pages),historical_pages=len(oldreader.pages),
        historical_source=str((OLD/'figure4_supplement.pdf').relative_to(A.REPO)),historical_sha256=R.sha(OLD/'figure4_supplement.pdf'),
        total_pages=len(writer.pages),historical_pages_appended_without_rendering=True),indent=2))

def render():
    copy_baseline()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6.2,'axes.labelsize':6.3,'axes.linewidth':.5,
        'xtick.labelsize':5.4,'ytick.labelsize':5.4,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
    fig=plt.figure(figsize=(180/25.4,420/25.4))
    grid=fig.add_gridspec(5,1,left=.082,right=.945,bottom=.045,top=.971,
        height_ratios=[1.9,1.18,2.,1.9,4.3],hspace=.26)
    nulls(fig,grid[0]);alpha(fig,grid[1]);centred_examples(fig,grid[2]);drift(fig,grid[3]);birth_panels(fig,grid[4])
    for label,spec in zip('ABCDE',grid):section(fig,spec,label)
    for ext in ('png','pdf','svg'):fig.savefig(OUT/('figure4_main.'+ext),dpi=400 if ext=='png' else 300)
    plt.close(fig);supplement()
    META.joinpath('round5_render_receipt.json').write_text(json.dumps(dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        width_mm=180,height_mm=420,main_titles=False,black_panel_labels=True,
        rare_multiple_cutoff=.05,closure_input_max=.01,initial_heading_centered=True,
        source_sha256={p.name:R.sha(p) for p in Path(__file__).parent.glob('figure4_round5*.py')},
        output_sha256={p.name:R.sha(p) for p in OUT.glob('figure4_*.*') if p.is_file()}),indent=2))

if __name__=='__main__':render()
