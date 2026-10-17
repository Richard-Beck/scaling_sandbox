"""Build review candidates for annotated WGD−/WGD+ representative track plots.

The tracks are existing model-generated draws under the primary Gaussian OU
hierarchical fit. They are not observed cell tracks. Run from the repository
root with the existing scaling_sandbox environment.
"""
from pathlib import Path
import csv, hashlib, json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parents[2]
FIT_DIR = ROOT / 'wgd_langevin/gaussian_hierarchy/results'
OUT = ROOT / 'reports/figure_revision/figures/figure3_review'
OUT.mkdir(parents=True, exist_ok=True)
NAVY = '#19394d'
MINUS = '#168b87'
PLUS = '#d87529'
GRAY = '#6d7c86'
CAD, DURATION, LAG = 1/6, 12., 11.1111111111
STUDY_INDEX = 0  # Xuan2018, whose fitted measures include sampled speed and late-lag MSD.

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'pdf.fonttype':42})

# Use the same exact integrated-OU displacement function as the fit pipeline.
sys.path.insert(0, str(ROOT / 'wgd_langevin/persistence_only'))
from campaign import displacement_variance

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def selected_fit():
    fs=[]
    for p in FIT_DIR.glob('fit_*.json'):
        d=json.loads(p.read_text())
        if d['scenario']=='primary': fs.append((d['total_score'],p,d))
    if not fs: raise FileNotFoundError('No completed primary hierarchy fit was found')
    return min(fs,key=lambda x:(x[0],x[1].name))

def track_metrics(track):
    step=np.diff(track,axis=0)
    path=float(np.linalg.norm(step,axis=1).sum())
    endpoint=float(np.linalg.norm(track[-1]-track[0]))
    return {'path_length_um':path,'mean_sampled_speed_um_h':path/DURATION,
            'endpoint_displacement_um':endpoint,'directness':endpoint/path}

def choose_representative(tracks):
    """Closest robustly standardized candidate to group medians.

    Score speed, directness, and endpoint displacement; path length is omitted
    because it is speed multiplied by the fixed 12 h duration.
    """
    metrics=[track_metrics(t) for t in tracks]
    keys=['mean_sampled_speed_um_h','directness','endpoint_displacement_um']
    a=np.array([[m[k] for k in keys] for m in metrics])
    med=np.median(a,axis=0)
    q25,q75=np.quantile(a,[.25,.75],axis=0)
    scale=np.maximum((q75-q25)/1.349,np.maximum(np.abs(med)*1e-8,1e-10))
    score=np.sum(((a-med)/scale)**2,axis=1)
    i=int(np.argmin(score))
    return tracks[i],metrics[i],{'index':i,'score':float(score[i]),
        'group_medians':dict(zip(keys,map(float,med))),
        'robust_scales_iqr_over_1.349':dict(zip(keys,map(float,scale))),
        'scored_metrics':keys,'path_length_omitted_as_redundant_with_speed':True}

def msd_expected(t,tau,eta):
    # Exact two-dimensional expectation; eta is the Wiener coefficient.
    return 2*eta**2*displacement_variance(np.asarray(t),np.array([1/tau]),np.array([tau/2]))

def arrow(ax,a,b,color=GRAY,lw=1.5,style='-|>'):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=11,
                                 color=color,lw=lw,shrinkA=0,shrinkB=0))

def axes_track(ax,track,color,linewidth=1.55,markers=True):
    ax.plot(track[:,0],track[:,1],color=color,lw=linewidth,zorder=2)
    ax.scatter(*track[0],s=24,marker='o',facecolor='white',edgecolor=color,lw=1.2,zorder=4)
    ax.scatter(*track[-1],s=38,marker='*',facecolor=color,edgecolor='white',lw=.4,zorder=5)

def fmt_metrics(m):
    return (f"Path length: {m['path_length_um']:.0f} µm\n"
            f"Mean sampled speed: {m['mean_sampled_speed_um_h']:.1f} µm/h\n"
            f"End displacement: {m['endpoint_displacement_um']:.0f} µm\n"
            f"Directness: {m['directness']:.2f}")

def style_axis(ax,lim):
    ax.set_aspect('equal',adjustable='box');ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim)
    ax.set_xlabel('Position x (µm)');ax.set_ylabel('Position y (µm)')
    ax.grid(alpha=.15,lw=.5)

def option1(tracks,metrics,lim,msd):
    fig=plt.figure(figsize=(12.4,7.0),layout='constrained')
    gs=fig.add_gridspec(2,2,width_ratios=[1,1.12])
    for g,(label,c) in enumerate([('WGD−',MINUS),('WGD+',PLUS)]):
        ax=fig.add_subplot(gs[g,0]);t=tracks[g];axes_track(ax,t,c);style_axis(ax,lim)
        arrow(ax,t[0],t[-1],NAVY,1.25)
        arrow(ax,t[18],t[19],GRAY,.9)
        ax.annotate('one 10-min step',xy=t[19],xytext=(.98,.72),textcoords='axes fraction',
                    ha='right',va='top',fontsize=7.5,color=GRAY,
                    arrowprops=dict(arrowstyle='->',color=GRAY,lw=.7))
        ax.set_title(label+' · one median-near simulated track',loc='left',color=NAVY,weight='bold')
        ax.text(.02,.98,'○ start     ★ end     arrow = net displacement',transform=ax.transAxes,
                va='top',fontsize=8,color=GRAY)
        ax=fig.add_subplot(gs[g,1]);ax.axis('off')
        ax.text(.04,.85,fmt_metrics(metrics[g]),va='top',linespacing=1.65,fontsize=11,
                color=NAVY,bbox=dict(boxstyle='round,pad=.6',fc='#f3f6f7',ec='#d8e0e4'))
        ax.text(.04,.35,f"MSD at {LAG:.2f}-h lag: {msd[g]:,.0f} µm²\n(mean squared displacement across fitted cells)",
                va='top',fontsize=9.2,color=NAVY,linespacing=1.5)
        ax.text(.04,.12,'Speed = path length / 12 h\nDirectness = end displacement / path length',
                va='top',fontsize=9,color=GRAY,linespacing=1.55)
    fig.suptitle('Option 1 · Direct labels beside each track\nXuan2018 model simulations · 12 h · 10-min sampling',fontsize=14,weight='bold',color=NAVY)
    return fig

def option2(tracks,metrics,lim,msd):
    fig=plt.figure(figsize=(10.4,8.0),layout='constrained')
    gs=fig.add_gridspec(1,2,width_ratios=[1.28,.72])
    ax=fig.add_subplot(gs[0,0])
    for g,(label,c) in enumerate([('WGD−',MINUS),('WGD+',PLUS)]):
        t=tracks[g];axes_track(ax,t,c,1.65)
        # Sparse time dots make temporal ordering visible without covering path.
        ids=np.array([0,18,36,54,72])
        ax.scatter(t[ids,0],t[ids,1],s=16,color=c,edgecolor='white',lw=.4,zorder=3)
        arrow(ax,t[0],t[-1],NAVY,1.15)
        arrow(ax,t[18],t[19],c,.8)
    style_axis(ax,lim);ax.set_title('Same origin and spatial scale',loc='left',color=NAVY,weight='bold')
    ax.legend(handles=[Line2D([],[],color=MINUS,lw=2,label='WGD−'),Line2D([],[],color=PLUS,lw=2,label='WGD+')],frameon=False,loc='lower left')
    ax=fig.add_subplot(gs[0,1]);ax.axis('off')
    ax.text(0,1,'Read the track',fontsize=13,weight='bold',color=NAVY,va='top')
    notes=[('1','Path length','sum of all sampled step lengths'),('2','Speed','path length ÷ elapsed time'),('3','Directness','end displacement ÷ path length'),('4','MSD','mean squared displacement across simulated cells at a fixed lag')]
    y=.88
    for n,title,body in notes:
        ax.text(.02,y,n,ha='center',va='center',color='white',weight='bold',
                bbox=dict(boxstyle='circle,pad=.35',fc=NAVY,ec='none'))
        ax.text(.14,y+.015,title,va='center',weight='bold',color=NAVY,fontsize=10)
        ax.text(.14,y-.035,body,va='top',color=GRAY,fontsize=8.5,wrap=True)
        y-=.18
    ax.text(.02,.14,f'Model-expected MSD ({LAG:.2f} h lag)\nWGD−  {msd[0]:,.0f} µm²\nWGD+  {msd[1]:,.0f} µm²',
            va='top',color=NAVY,fontsize=9.5,linespacing=1.5,
            bbox=dict(boxstyle='round,pad=.5',fc='#f3f6f7',ec='#d8e0e4'))
    fig.suptitle('Option 2 · Overlay with numbered measurement key\nXuan2018 model simulations · 12 h · 10-min sampling',fontsize=14,weight='bold',color=NAVY)
    ax.text(.02,.01,'Circle: start  ·  star: end  ·  dark arrow: net displacement  ·  colored arrow: one 10-min step',transform=ax.transAxes,ha='left',va='bottom',fontsize=7.2,color=GRAY)
    return fig

def option3(tracks,metrics,lim,msd,tau_eta,review=True):
    fig=plt.figure(figsize=(12.0,8.2),layout='constrained')
    gs=fig.add_gridspec(2,2,height_ratios=[1.1,.9],hspace=.2)
    for g,(label,c) in enumerate([('WGD−',MINUS),('WGD+',PLUS)]):
        ax=fig.add_subplot(gs[0,g]);t=tracks[g];axes_track(ax,t,c,1.5);style_axis(ax,lim)
        arrow(ax,t[0],t[-1],NAVY,1.05)
        for j in (18,36,54):arrow(ax,t[j],t[j+1],c,.75)
        # Equal-duration displacement chords illustrate the fixed-lag samples
        # that are squared and averaged to form an empirical track MSD.
        for a,b in ((0,12),(60,72)):
            ax.plot(t[[a,b],0],t[[a,b],1],ls='--',color=GRAY,lw=1.15,zorder=4)
            ax.scatter(t[[a,b],0],t[[a,b],1],s=12,facecolor='white',edgecolor=GRAY,lw=.8,zorder=5)
        ax.annotate('same Δt = 2 h',xy=t[12],xytext=(.97,.78),textcoords='axes fraction',
                    ha='right',va='top',fontsize=7.5,color=GRAY,
                    arrowprops=dict(arrowstyle='->',color=GRAY,lw=.65))
        ax.set_title(f"{label} · speed {metrics[g]['mean_sampled_speed_um_h']:.1f} µm/h · path {metrics[g]['path_length_um']:.0f} µm",loc='left',weight='bold',color=NAVY)
        ax.text(.03,.96,f"Directness = {metrics[g]['directness']:.2f}\n(net displacement / path length)",transform=ax.transAxes,va='top',fontsize=8,color=NAVY)
    ax=fig.add_subplot(gs[1,:])
    tt=np.linspace(.03,DURATION,300)
    for g,(label,c) in enumerate([('WGD−',MINUS),('WGD+',PLUS)]):
        tau,eta=tau_eta[g]
        ax.plot(tt,msd_expected(tt,tau,eta),color=c,lw=2.3,label=f'{label} fitted OU expectation')
        ax.scatter([LAG],[msd[g]],color=c,s=38,zorder=3)
    ax.axvline(LAG,color=GRAY,ls=':',lw=1)
    ax.set(xlabel='Lag time (h)',ylabel='Mean squared displacement (µm²)',xlim=(0,DURATION),ylim=(0,None))
    ax.grid(alpha=.17);ax.legend(frameon=False,ncol=2,loc='upper left')
    ax.text(.99,.04,f'At fitted lag {LAG:.2f} h:  WGD− {msd[0]:,.0f} µm²  |  WGD+ {msd[1]:,.0f} µm²',
            transform=ax.transAxes,ha='right',va='bottom',fontsize=9,color=NAVY,
            bbox=dict(boxstyle='round,pad=.35',fc='white',ec='#d8e0e4'))
    heading='Option 3 · Tracks above, population MSD below' if review else 'Representative WGD− and WGD+ tracks with fitted MSD'
    fig.suptitle(heading+'\nXuan2018 model simulations · 12 h · 10-min sampling',fontsize=14,weight='bold',color=NAVY)
    ax.text(.99,.96,'Dashed chords: equal 2 h windows. Track MSD averages squared\ndisplacements across tracks and start times at a fixed lag; curve = fitted model expectation.',
            transform=ax.transAxes,ha='right',va='top',fontsize=7.5,color=GRAY,
            bbox=dict(boxstyle='round,pad=.35',fc='white',ec='#d8e0e4',alpha=.94))
    return fig

def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=240,bbox_inches='tight')
    plt.close(fig)

def main():
    score,fit_path,fit=selected_fit()
    tr_path=FIT_DIR/'tracks_primary.npz'
    data=np.load(tr_path)
    pars=np.exp(np.asarray(fit['z'][:20]).reshape(5,2,2))
    chosen=[];metrics=[];selection={}
    for g in (0,1):
        t,m,s=choose_representative(data[f'Xuan2018_{g}'])
        chosen.append(t);metrics.append(m);selection['WGD−' if g==0 else 'WGD+']=s
    tau_eta=pars[STUDY_INDEX]
    msd=[float(msd_expected(LAG,*tau_eta[g])) for g in (0,1)]
    fit_msd={int(r['group']):float(r['predicted']) for r in fit['rows']
             if r['study']=='Xuan2018' and r['metric']=='msd'}
    if any(not np.isclose(msd[g],fit_msd[g],rtol=1e-6) for g in (0,1)):
        raise ValueError(f"Displayed MSD does not reproduce selected fit: {msd} vs {fit_msd}")
    lim=float(max(np.abs(t).max() for t in chosen)*1.13)
    fs=[option1(chosen,metrics,lim,msd),option2(chosen,metrics,lim,msd),option3(chosen,metrics,lim,msd,tau_eta)]
    for i,f in enumerate(fs,1):save(f,f'option{i}')
    save(option3(chosen,metrics,lim,msd,tau_eta,review=False),'figure3_annotated_tracks')
    # Review sheet: preserve labels and sufficient resolution for readable text.
    from PIL import Image,ImageOps,ImageDraw,ImageFont
    ims=[Image.open(OUT/f'option{i}.png').convert('RGB') for i in range(1,4)]
    target_w=1600; gaps=24; resized=[]
    for im in ims:
        h=round(im.height*target_w/im.width);resized.append(im.resize((target_w,h),Image.Resampling.LANCZOS))
    pad=42; sheet=Image.new('RGB',(target_w+2*pad,sum(x.height for x in resized)+gaps*2+2*pad),(245,247,248))
    y=pad
    for im in resized:
        sheet.paste(im,(pad,y));y+=im.height+gaps
    sheet.save(OUT/'comparison.png',quality=95)
    for g,label in enumerate(['WGD_minus','WGD_plus']):
        np.savetxt(OUT/f'{label}_representative_track.csv',
                   np.column_stack((np.arange(73)*CAD,chosen[g])),delimiter=',',
                   header='time_h,x_um,y_um',comments='')
    metrics_doc={label:dict(metrics[g],fit_parameters={'tau_h':float(tau_eta[g,0]),'eta_um_h_3_2':float(tau_eta[g,1])},
                            msd_expected_at_11_111h_um2=msd[g],selection=selection[label])
                 for g,label in enumerate(['WGD−','WGD+'])}
    (OUT/'figure3_tracks_provenance.json').write_text(json.dumps({
        'created_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),
        'inputs':{str(fit_path.relative_to(ROOT)):digest(fit_path),str(tr_path.relative_to(ROOT)):digest(tr_path)},
        'fit_scenario':'primary','fit_start_file':fit_path.name,'fit_total_score':score,
        'user_selected_design':'Option 3, saved without the review-option label as figure3_annotated_tracks.png',
        'study':'Xuan2018','group_order':['WGD−','WGD+'],'tracks_source':'Existing simulated primary-fit OU track ensemble, tracks_primary.npz; these are model-generated and are not observed cell tracks.',
        'schedule':{'duration_h':DURATION,'cadence_h':CAD,'frames':72},
        'selection_rule':'From 32 existing simulated tracks per group, choose track minimizing sum of squared robust standardized deviations from group medians in sampled speed, directness, and endpoint displacement. Path length excluded from score because it is exactly speed × fixed duration. The same two tracks appear in all options.',
        'metric_definitions':{'sampled_speed':'sum of sampled step lengths / elapsed time','path_length':'sum of sampled step lengths','directness':'single displayed track endpoint displacement / its path length','MSD':'exact fitted OU ensemble expected squared displacement using the campaign displacement_variance helper at lag 11.111 h; distinct from a single track endpoint squared displacement'},
        'parameters_and_track_metrics':metrics_doc,
        'limitations':'Primary hierarchy fit to approximate literature summary extractions. WGD+ is the study-specific giant polyploid contrast for Xuan2018, not evidence that every contrast is an isogenic WGD experiment. Selected single paths illustrate fitted groups and do not represent observed individuals or experimental variability.'},indent=2,ensure_ascii=False)+'\n')
    with (OUT/'track_metrics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['group','path_length_um','mean_sampled_speed_um_h','endpoint_displacement_um','directness','msd_expected_at_11_111h_um2','tau_h','eta_um_h_3_2','representative_index_of_32'])
        w.writeheader()
        for g,label in enumerate(['WGD−','WGD+']):
            w.writerow(dict(group=label,**metrics[g],msd_expected_at_11_111h_um2=msd[g],tau_h=tau_eta[g,0],eta_um_h_3_2=tau_eta[g,1],representative_index_of_32=selection[label]['index']))
    print(f'Built options and comparison sheet in {OUT}')

if __name__=='__main__': main()
